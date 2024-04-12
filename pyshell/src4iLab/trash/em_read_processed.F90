!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module em_read_processed

  ! --- modules ------------------------------

  use GO, only : gol, goErr, goPr
  use TIPP_Base,              only : TEmis
  use File_HDF

  implicit none

  ! --- in/out -----------------------------

  private

  public :: Em_Read_Processed_Init, Em_Read_Processed_Done

  ! --- const ------------------------------

  character(len=*), parameter     ::  mname = 'Em_Read_Processed'


contains

  !---------------------------------------------------------
  ! Open HDF emission file and initialise em structure
  ! for reading (read attributes and allocate fields)
  !---------------------------------------------------------

  subroutine Em_Read_Processed_Init( em, dir, species, category, source, &
                                      year, hybrid_name, region_name, &
                                      status, actual_year )

    ! --- modules ------------------------------

    use GO                  , only : pathsep
    use GO                  , only : NewDate, IncrDate, Get_End_Of, wrtgol, operator(+)
    use GO                  , only : goReadFromLine
    use tipp_base,            only : Emis_Sparse_Init
    use tipp_base,            only : Init
    use os_specs,             only : MAX_FILENAME_LEN, LONG_STR_LEN

    ! --- in/out ---------------------------------

    type(TEmis), intent(out)         :: em
    character(len=*), intent(in)     :: dir   ! directory in which file should be
    character(len=*), intent(in)     :: species
    character(len=*), intent(in)     :: category
    character(len=*), intent(in)     :: source
    integer, intent(in)              :: year
    character(len=*), intent(in)     :: hybrid_name
    character(len=*), intent(in)     :: region_name
    integer, intent(out)             :: status
    integer, intent(in), optional    :: actual_year

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//'/Em_Read_Processed_Init'

    ! --- local ------------------------------

    integer                          :: ayear
    character(len=MAX_FILENAME_LEN)  :: fname, fname2
    character(len=LONG_STR_LEN)      :: line
    logical                          :: exist
    type(THdfFile)                   :: hdf
    type(TSds)                       :: sds, sds_valid
    integer(1), allocatable          :: bvalid(:,:,:,:)
    integer, dimension(4)            :: help
    character(len=80)                :: title, date_time
    integer                          :: attr_index
    logical                          :: no_lat, no_lon
    integer, allocatable             :: times(:,:)
    integer                          :: ip
    integer                          :: lenp
    integer                          :: jd

    ! --- begin ------------------------------

    ! actual year ...
    ayear = year
    if ( present(actual_year) ) ayear = actual_year

    ! Construct file name
    if ( dir == '.' ) then
      write( fname, '(2a,a,"-",a,"-",a,"-",i4.4,"-",a,"-",a,".hdf")' ) &
           trim(dir), pathsep, &
           trim(species), trim(category), trim(source), year, &
           trim(hybrid_name), trim(region_name)
    else
      write( fname, '(6a,a,"-",a,"-",a,"-",i4.4,"-",a,"-",a,".hdf")' ) &
           trim(dir), pathsep, trim(hybrid_name), pathsep, trim(region_name), pathsep, &
           trim(species), trim(category), trim(source), year, &
           trim(hybrid_name), trim(region_name)
    end if

    ! file exists?
    inquire( file=trim(fname), exist=exist )

    if ( .not. exist ) then
       ! If 3D does not exist, then try 2D (i.e. 'sfc' instead of 'ml..')
      if ( dir == '.' ) then
         write( fname2, '(2a,a,"-",a,"-",a,"-",i4.4,"-sfc-",a,".hdf")' ) &
              trim(dir), pathsep, &
              trim(species), trim(category), trim(source), year, &
              trim(region_name)
      else
         write( fname2, '(6a,a,"-",a,"-",a,"-",i4.4,"-sfc-",a,".hdf")' ) &
              trim(dir), pathsep, 'sfc', pathsep, trim(region_name), pathsep, &
              trim(species), trim(category), trim(source), year, &
              trim(region_name)
      end if

       ! file exists?
       inquire( file=trim(fname2), exist=exist )

       if ( exist ) then
          ! copy name:
          fname = trim(fname2)
          ! store full path:
          em%filename = trim(fname2)
       else
          write (gol,'("cannot find one of these emission files:")'); call goErr
          write (gol,'("  3d file    : ",a)') trim(fname); call goErr
          write (gol,'("  2d file    : ",a)') trim(fname2); call goErr
          write (gol,'("parameters:")'); call goErr
          write (gol,'("  directory  : ",a)') trim(dir); call goErr
          write (gol,'("  species    : ",a)') trim(species); call goErr
          write (gol,'("  category   : ",a)') trim(category); call goErr
          write (gol,'("  source     : ",a)') trim(source); call goErr
          write (gol,'("  year       : ",i4.4)') year; call goErr
          write (gol,'("  region key : ",a)') trim(region_name); call goErr
          write (gol,'("  hybrid key : ",a," or sfc")') trim(hybrid_name); call goErr
          TRACEBACK; status=1; return
       end if

    end if

    !! info ...
    !write (gol,'("opening ",a)') trim(fname); call goPr

    ! Init HDF
    call Init( hdf, fname, 'read', status )

    ! file format, only present for format 2 and newer;
    status = 1  ! silent
    call FindAttribute( hdf, 'format', attr_index, status )
    IF_ERROR_RETURN(status=1)
    if ( status < 0 ) then
      em%format_nr = 1
    else
      call ReadAttribute( hdf, 'format', em%format_nr, status )
      IF_ERROR_RETURN(status=1)
    end if

    ! how stored ?
    if ( em%format_nr == 1 ) then
      em%storage = 'full'
    else
      call ReadAttribute( hdf, 'storage', em%storage, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! Read global attributes
    call ReadAttribute( hdf, 'title', title, status )
    IF_NOTOK_RETURN(status=1)

    call ReadAttribute( hdf, 'date_time', date_time, status )
    IF_NOTOK_RETURN(status=1)

    ! how stored ?
    select case ( trim(em%storage) )

      !~~ full 4D storage:
      case ( 'full' )

        ! Init SDS
        call Init( sds, hdf, 'field', status )
        IF_NOTOK_RETURN(status=1)

      !~~ sparse storage:
      case ( 'sparse' )

        ! open dataset with indices:
        call Init( sds, hdf, 'sparse_index', status )
        IF_NOTOK_RETURN(status=1)
        ! get shape:
        call ReadAttribute( sds, 'sparse_n', em%sparse_n, status )
        IF_NOTOK_RETURN(status=1)
        ! storage:
        allocate( em%sparse_index(4,em%sparse_n) )
        ! read:
        call ReadData( sds, em%sparse_index, status )
        IF_NOTOK_RETURN(status=1)
        ! close:
        call Done( sds, status )
        IF_NOTOK_RETURN(status=1)

        ! open dataset with values:
        call Init( sds, hdf, 'sparse_field', status )
        IF_NOTOK_RETURN(status=1)

    end select

    ! Read attributes
    call ReadAttribute( sds, 'name', em%name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'original_filename', em%original_filename, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'original_dimensions', help, status )
    IF_NOTOK_RETURN(status=1)
    em%nx_org = help(1)
    em%ny_org = help(2)
    em%nz_org = help(3)
    em%np_org = help(4)
    call ReadAttribute( sds, 'filetype', em%filetype, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'source', em%source, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'species', em%species, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'category', em%category, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'category_long', em%category_long, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'original_unit', em%original_unit, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'unit', em%unit, status )
    IF_NOTOK_RETURN(status=1)
    !>>> only scaling at runtime given factor in emission list
    !call ReadAttribute( sds, 'scaling_factor', em%scaling_factor, status )
    !IF_NOTOK_RETURN(status=1)
    !<<<
    call ReadAttribute( sds, 'year', em%year, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'region_name', em%region_name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'hybrid_name', em%hybrid_name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'dimensions', help, status )
    IF_NOTOK_RETURN(status=1)
    em%nx = help(1)
    em%ny = help(2)
    em%nz = help(3)
    em%np = help(4)

    call ReadAttribute( sds, 'mole_mass_full_molecule', em%xm_fm, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'mole_mass_element', em%xm_el, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'mole_mass_unit', em%xm_unit, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'mole_mass_undefined', em%xm_undefined, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'total_field', em%total, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'total_field_unit', em%total_unit, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'additional_comments', em%additional_comments, status)
    IF_NOTOK_RETURN(status=1)

    !! info ...
    !write (gol,'("    ",a," - ",a)') trim(em%category), trim(em%category_long); call goPr
    !write (gol,'("     global total = ",e10.4," ",a)') em%total, trim(em%total_unit); call goPr

    ! Allocate arrays
    allocate( em%x_grid(em%nx) )
    allocate( em%y_grid(em%ny) )
    if ( em%nz == 1 ) then
       allocate( em%z_a_grid(1) )
       allocate( em%z_b_grid(1) )
    else
       allocate( em%z_a_grid(em%nz+1) )
       allocate( em%z_b_grid(em%nz+1) )
    end if
    allocate( em%p_grid(em%np) )
    allocate( em%t1(em%np) )
    allocate( em%t2(em%np) )

    ! no 4D-field read yet:
    nullify( em%field )

    ! init sparse storage:
    call EMIS_Sparse_Init( em, status )
    IF_NOTOK_RETURN(status=0)

    ! how stored ?
    select case ( trim(em%storage) )

      !~~ full 4D storage:
      case ( 'full' )

        ! storage:
        allocate( em%field(em%nx,em%ny,em%nz,em%np) )
        ! Read data
        call ReadData( sds, em%field, status )
        IF_NOTOK_RETURN(status=1)

        ! storage:
        allocate( em%valid(em%nx,em%ny,em%nz,em%np) )
        ! data set with valid flags available not available for first format:
        if ( em%format_nr == 1 ) then
          ! assume all are valid ...
          em%valid = .true.
        else
          ! storage:
          allocate( bvalid(em%nx,em%ny,em%nz,em%np) )
          ! open data set:
          call Init( sds_valid, hdf, 'valid', status )
          IF_NOTOK_RETURN(status=1)
          ! read byte flags:
          call ReadData( sds_valid, bvalid, status )
          IF_NOTOK_RETURN(status=1)
          ! close:
          call Done( sds_valid, status )
          IF_NOTOK_RETURN(status=1)
          ! convert:  1 = true (valid), 0 = false (invalid)
          em%valid = bvalid == 1
          ! clear:
          deallocate( bvalid )
        end if

      !~~ sparse storage:
      case ( 'sparse' )

        ! get shape:
        call ReadAttribute( sds, 'sparse_n', em%sparse_n, status )
        IF_NOTOK_RETURN(status=1)

        ! storage:
        allocate( em%sparse_field(em%sparse_n) )

        ! read:
        call ReadData( sds, em%sparse_field, status )
        IF_NOTOK_RETURN(status=1)

    end select

    ! axes stored as sds attributes:
    if ( em%format_nr == 1 ) then
      ! axes definition:
      call ReadAttribute( sds, 'x_def', em%x_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'y_def', em%y_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'z_def', em%z_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'p_def', em%p_def, status )
      IF_NOTOK_RETURN(status=1)
      ! Read attributes
      call ReadAttribute( sds, 'x_grid', em%x_grid, status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'y_grid', em%y_grid, status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'z_a_grid', em%z_a_grid, status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'z_b_grid', em%z_b_grid, status )
      IF_NOTOK_RETURN(status=1)

      ! character line:
      call ReadAttribute( sds, 'p_grid', line, status )
      IF_NOTOK_RETURN(status=1)
      ! time axes
      do ip = 1, em%np
        ! extract part:
        lenp = ceiling(len_trim(line)/real(em%np))
        em%p_grid(ip) = trim(line((ip-1)*lenp+1:ip*lenp))
        ! switch ...
        if ( em%p_grid(ip)(1:4) == 'Year' ) then
          em%t1(ip) = NewDate( year=ayear, month=1, day=1 )
          em%t2(ip) = Get_End_Of( em%t1(ip), 'year' )
        else if ( any( em%p_grid(ip)(1:3) == (/'Jan','Feb','Mar','Apr','May','Jun',&
                                               'Jul','Aug','Sep','Oct','Nov','Dec'/) ) ) then
          em%t1(ip) = NewDate( year=ayear, month=ip, day=1 )
          em%t2(ip) = Get_End_Of( em%t1(ip), 'month' )
        else if ( em%p_grid(ip)(1:2) == 'JD' ) then
          read (em%p_grid(ip)(3:5),'(i3)') jd
          em%t1(ip) = NewDate( year=ayear, month=1, day=1 ) + IncrDate( day=jd-1 )
          if ( ip > 1 ) em%t2(ip-1) = em%t1(ip)
          if ( ip == em%np ) em%t2(ip) = Get_End_Of( em%t1(ip), 'year' )
        else
          write (gol,'("unsupported period name to set t1/t2 : ",a)') trim(em%p_grid(ip)); call goErr
          write (gol,'("filename : ",a)') trim(em%filename); call goErr
          write (gol,'("period   : ",i6)') ip; call goErr
          TRACEBACK; status=1; return
        end if
      end do

    end if

    ! Close SDS
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! axes stored as data sets:
    if ( em%format_nr > 1 ) then
      ! read x-axis:
      call Init( sds, hdf, 'x_grid', status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'x_def', em%x_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadData( sds, em%x_grid, status )
      IF_NOTOK_RETURN(status=1)
      call Done( sds, status )
      IF_NOTOK_RETURN(status=1)
      ! read y-axis:
      call Init( sds, hdf, 'y_grid', status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'y_def', em%y_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadData( sds, em%y_grid, status )
      IF_NOTOK_RETURN(status=1)
      call Done( sds, status )
      IF_NOTOK_RETURN(status=1)
      ! read z-axis, a-coefficients:
      call Init( sds, hdf, 'z_a_grid', status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'z_def', em%z_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadData( sds, em%z_a_grid, status )
      IF_NOTOK_RETURN(status=1)
      call Done( sds, status )
      IF_NOTOK_RETURN(status=1)
      ! read z-axis, b-coefficients:
      call Init( sds, hdf, 'z_b_grid', status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'z_def', em%z_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadData( sds, em%z_b_grid, status )
      IF_NOTOK_RETURN(status=1)
      call Done( sds, status )
      IF_NOTOK_RETURN(status=1)
      ! read p-axis:
      call Init( sds, hdf, 'p_grid', status )
      IF_NOTOK_RETURN(status=1)
      call ReadAttribute( sds, 'p_def', em%p_def, status )
      IF_NOTOK_RETURN(status=1)
      call ReadData( sds, em%p_grid, status )
      IF_NOTOK_RETURN(status=1)
      call Done( sds, status )
      IF_NOTOK_RETURN(status=1)
      ! read time values:
      !~ storage:
      allocate( times(6,em%np) )
      !~ read start times:
      call Init( sds, hdf, 'time1', status )
      IF_NOTOK_RETURN(status=1)
      call ReadData( sds, times, status )
      IF_NOTOK_RETURN(status=1)
      call Done( sds, status )
      IF_NOTOK_RETURN(status=1)
      !~ actual year:
      times(1,:) = ayear
      !~ convert:
      do ip = 1, em%np
        em%t1(ip) = NewDate( time6=times(:,ip) )
      end do
      !~ read start times:
      call Init( sds, hdf, 'time2', status )
      IF_NOTOK_RETURN(status=1)
      call ReadData( sds, times, status )
      IF_NOTOK_RETURN(status=1)
      call Done( sds, status )
      IF_NOTOK_RETURN(status=1)
      !~ actual year:
      times(1,:) = ayear
      times(1,em%np) = ayear+1
      !~ convert:
      do ip = 1, em%np
        em%t2(ip) = NewDate( time6=times(:,ip) )
      end do
      ! clear:
      deallocate( times )
    end if

    !! debug ...
    !write (gol,'("aaa1 ",a)') trim(em%name); call goPr
    !do ip = 1, em%np
    !  call wrtgol( '  a2 "'//trim(em%p_grid(ip))//'" ',  em%t1(ip), ' - ', em%t2(ip) ); call goPr
    !end do

    ! Close HDF file
    call Done( hdf, status )
    IF_NOTOK_RETURN(status=1)

    ! Initialize ll grid
    no_lon = em%nx == 1   ! zonal field
    no_lat = em%ny == 1   ! meridional field
    call Init( em%region_name, no_lon, no_lat, em%lli, status )
    IF_NOTOK_RETURN(status=1)

    ! Initialize hybrid grid
    call Init( em%hybrid_name, em%levi, status )
    IF_NOTOK_RETURN(status=1)

    ! ok:
    status = 0

  end subroutine Em_Read_Processed_Init


  !---------------------------------------------------------
  ! Close HDF emission file and finalise em structure
  ! (deallocate fields)
  !---------------------------------------------------------

  subroutine Em_Read_Processed_Done( em, status )

    ! --- modules --------------------------------

    use grid,                    only : Done
    use tipp_base,               only : Emis_Sparse_Done

    ! --- in/out ---------------------------------

    type(TEmis), intent(inout)       :: em
    integer, intent(out)             :: status

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//'/Em_Read_Processed_Done'

    !--- begin ------------------------------------------------

    ! o full storage:
    if ( associated(em%field) ) deallocate( em%field )
    if ( associated(em%valid) ) deallocate( em%valid )
    ! o sparse storage:
    call Emis_Sparse_Done( em, status )
    IF_NOTOK_RETURN(status=1)

    ! grid:
    deallocate( em%x_grid )
    deallocate( em%y_grid )
    deallocate( em%z_a_grid )
    deallocate( em%z_b_grid )
    deallocate( em%p_grid )
    deallocate( em%t1 )
    deallocate( em%t2 )

    call Done( em%lli, status )
    IF_NOTOK_RETURN(status=1)

    call Done( em%levi, status )
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine Em_Read_Processed_Done

end module em_read_processed
