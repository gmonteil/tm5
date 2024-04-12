!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module bc_read_processed

  ! --- modules ------------------------------

  use go,                     only : gol, goErr
  use tipp_base,              only : TBouCon
  use file_hdf

  implicit none

  ! --- in/out -----------------------------

  private

  public :: Init, ReadField, Done

  ! --- const ------------------------------

  character(len=*), parameter     ::  mname = 'bc_read_processed'

  ! --- interfaces -------------------------------------

  interface Init
     module procedure bc_read_processed_Init
  end interface

  interface ReadField
     module procedure bc_read_processed_ReadField
  end interface

  interface Done
     module procedure bc_read_processed_Done
  end interface

contains

  subroutine bc_read_processed_Init( hdf, dir, description, source, year, &
       hybrid_name, region_name, bc, status )
    !---------------------------------------------------------
    ! Open HDF boundary condition file and initialise bc structure
    ! for reading (read attributes and allocate fields)
    !---------------------------------------------------------

    ! --- modules ------------------------------
    use os_specs, only : MAX_FILENAME_LEN, SHORT_STR_LEN

    ! --- in/out ---------------------------------

    type(THdfFile), intent(out)      :: hdf
    character(len=*), intent(in)     :: dir   ! directory in which file should be
    character(len=*), intent(in)     :: description
    character(len=*), intent(in)     :: source
    integer, intent(in)              :: year
    character(len=*), intent(in)     :: hybrid_name
    character(len=*), intent(in)     :: region_name
    type(TBouCon), intent(out)       :: bc
    integer, intent(out)             :: status

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//'/bc_read_processed_Init'

    ! --- local ------------------------------

    character(len=MAX_FILENAME_LEN)  :: fname
    logical                          :: exist
    type(TSds)                       :: sds
    integer, dimension(4)            :: help
    character(len=SHORT_STR_LEN)     :: title, date_time

    ! --- begin ------------------------------

    ! Construct file name
    write( fname, '(a,"-",a,"-",i4.4,"-",a,"-",a,".hdf")' ) &
         trim(description), trim(source), year, &
         trim(hybrid_name), trim(region_name)
    fname = trim(dir)//'/'//trim(fname)

    ! file exists?
    inquire( file=fname, exist=exist )

    if ( .not. exist ) then
       ! If 3D does not exist, then try 2D (i.e. 'sfc' instead of 'ml..')
       write( fname, '(a,"-",a,"-",i4.4,"-sfc-",a,".hdf")' ) &
            trim(description), trim(source), year, trim(region_name)
       fname = trim(dir)//'/'//trim(fname)

       ! file exists?
       inquire( file=fname, exist=exist )

       if ( .not. exist ) then
          write (*,'("ERROR - cannot find boundary condition file")')
          write (*,'("      - directory:  ",a)') trim(dir)
          write (*,'("      - description: ",a)') trim(description)
          write (*,'("      - source:      ",a)') trim(source)
          write (*,'("      - year:        ",i4.4)') year
          write (*,'("      - hybrid key:  ",a," or sfc")') trim(hybrid_name)
          write (*,'("      - region key:  ",a)') trim(region_name)
          write (*,'("ERROR in ",a)') rname; status=1; return
       end if

    end if

    write (*,'(a,": opening ",a)') rname, trim(fname)

    ! Init HDF
    call Init( hdf, fname, 'read', status )
    IF_NOTOK_RETURN(status=1)

    ! Read global attributes
    call ReadAttribute( hdf, 'title', title, status )
    IF_NOTOK_RETURN(status=1)

    call ReadAttribute( hdf, 'date_time', date_time, status )
    IF_NOTOK_RETURN(status=1)

    ! Init SDS
    call Init( sds, hdf, 'field', status )
    IF_NOTOK_RETURN(status=1)

    ! Read attributes
    call ReadAttribute( sds, 'name', bc%name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'original_filename', bc%original_filename, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'original_dimensions', help, status )
    IF_NOTOK_RETURN(status=1)
    bc%nx_org = help(1)
    bc%ny_org = help(2)
    bc%nz_org = help(3)
    bc%np_org = help(4)
    call ReadAttribute( sds, 'description', bc%description, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'description_long', bc%description_long, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'source', bc%source, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'unit', bc%unit, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'scaling_factor', bc%scaling_factor, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'year', bc%year, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'region_name', bc%region_name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'hybrid_name', bc%hybrid_name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'dimensions', help, status )
    IF_NOTOK_RETURN(status=1)
    bc%nx = help(1)
    bc%ny = help(2)
    bc%nz = help(3)
    bc%np = help(4)
    call ReadAttribute( sds, 'x_def', bc%x_def, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'y_def', bc%y_def, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'z_def', bc%z_def, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'p_def', bc%p_def, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'additional_comments', bc%additional_comments, status )
    IF_NOTOK_RETURN(status=1)

    ! Close SDS
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! Allocate arrays
    allocate( bc%field(bc%nx,bc%ny,bc%nz,bc%np) )
    allocate( bc%x_grid(bc%nx) )
    allocate( bc%y_grid(bc%ny) )
    if ( bc%nz == 1 ) then
       allocate( bc%z_a_grid(1) )
       allocate( bc%z_b_grid(1) )
    else
       allocate( bc%z_a_grid(bc%nz+1) )
       allocate( bc%z_b_grid(bc%nz+1) )
    end if
    allocate( bc%p_grid(bc%np) )

    ! ok
    status = 0

  end subroutine bc_read_processed_Init


  subroutine bc_read_processed_ReadField( hdf, bc, status )
    !---------------------------------------------------------
    ! Read general HDF boundary condition file
    !---------------------------------------------------------

    ! --- modules ------------------------------

    use tipp_base,               only : Init

    ! --- in/out ---------------------------------

    type(THdfFile), intent(inout)    :: hdf
    type(TBouCon), intent(inout)     :: bc
    integer, intent(out)             :: status

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//'/bc_read_processed_ReadField'

    !--- local ------------------------------------------------

    type(TSds)                       :: sds
    logical                          :: no_lon, no_lat

    !--- begin ------------------------------------------------

    ! Init SDS
    call Init( sds, hdf, 'field', status )
    IF_NOTOK_RETURN(status=1)

    ! Read attributes
    call ReadAttribute( sds, 'x_grid', bc%x_grid, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'y_grid', bc%y_grid, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'z_a_grid', bc%z_a_grid, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'z_b_grid', bc%z_b_grid, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'p_grid', bc%p_grid, status )
    IF_NOTOK_RETURN(status=1)

    ! Read data
    call ReadData( sds, bc%field, status )
    IF_NOTOK_RETURN(status=1)

    ! Initialize ll grid
    no_lon = .false.
    if ( bc%nx == 1 ) no_lon = .true.   ! zonal field
    no_lat = .false.
    if ( bc%ny == 1 ) no_lat = .true.   ! meridional field

    call Init( bc%region_name, no_lon, no_lat, bc%lli, status )
    IF_NOTOK_RETURN(status=1)

    ! Initialize (possibly dummy) hybrid grid
    call Init( bc%hybrid_name, bc%levi, status )
    IF_NOTOK_RETURN(status=1)

    ! Close SDS
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine bc_read_processed_ReadField


  subroutine bc_read_processed_Done( hdf, bc, status )
    !---------------------------------------------------------
    ! Close HDF boundary condition file and finalise bc structure
    ! (deallocate fields)
    !---------------------------------------------------------

    ! --- modules ------------------------------

    use grid,                    only : Done

    ! --- in/out ---------------------------------

    type(THdfFile), intent(inout)    :: hdf
    type(TBouCon), intent(inout)     :: bc
    integer, intent(out)             :: status

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//'/bc_read_processed_Done'

    !--- begin ------------------------------------------------

    deallocate( bc%field )
    deallocate( bc%x_grid )
    deallocate( bc%y_grid )
    deallocate( bc%z_a_grid )
    deallocate( bc%z_b_grid )
    deallocate( bc%p_grid )

    call Done( bc%lli, status )
    IF_NOTOK_RETURN(status=1)

    call Done( bc%levi, status )
    IF_NOTOK_RETURN(status=1)

    ! Close HDF file
    call Done( hdf, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine bc_read_processed_Done

end module bc_read_processed
