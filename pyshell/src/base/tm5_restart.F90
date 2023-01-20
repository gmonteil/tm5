#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_NOTOK_MDF(action) if (status/=0) then; TRACEBACK; action; call MDF_CLose(fid,status); status=1; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"
!----------------------------------------------------------------------------
!                  TM5                                                      !
!----------------------------------------------------------------------------
!BOP
!
! !MODULE:  Restart
!
! !DESCRIPTION: Write and read restart files.
!\\
!\\
! !INTERFACE:
!
module Restart
  !
  ! !USES:
  !
  use GO     , only : gol, goPr, goErr
  use dims   , only : nregions
  use os_specs, only : MAX_RCKEY_LEN, MAX_FILENAME_LEN

  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public  ::  Restart_Init   ! read restart keys in rc file
  public  ::  Restart_Done   ! nothing yet
  public  ::  Restart_Save   ! wrapper around Restart_Write
  public  ::  Restart_Write  ! write a restart file
  public  ::  Restart_Read   ! read a restart file
  public  ::  rs_write        ! model must write restart
  !
  ! !PRIVATE DATA MEMBERS:
  !
  character(len=*), parameter  ::  mname = 'Restart'
  character(len=MAX_FILENAME_LEN)           ::  rs_write_dir
  logical                      ::  rs_write
  logical                      ::  rs_write_extra
  integer                      ::  rs_write_extra_dhour, rs_write_extra_hour
  integer                      ::  fid ! file id for IF_NOTOK_MDF macro
  !
  ! !REVISION HISTORY:
  !    8 Apr 2011 - P. Le Sager - Close MDF file if error occurs. This is
  !                 needed for mpi_abort not to hang. See TM5_MPI_Abort in
  !                 partools, and remarks below. Made IF_NOTOK_MDF macro for
  !                 that purpose.
  !   28 Apr 2011 - P. Le Sager - Read method : handle restart file with extra
  !                 tracers.
  !
  ! !REMARKS:
  ! (1) when an error occurs when accessing MDF files, you should first close
  !   the file before returning. The IF_NOTOK_MDF macro takes care of that.
  !   The only thing you need is to call it like that :
  !
  !      IF_NOTOK_MDF(fid=xxxx)
  !
  !   where you replace xxxx with the integer id (file handler) of the file
  !   you are accessing. Note that this does not solve all problems (but
  !   probably most of them): it is still possible that MDF_Close hangs...
  !
  !EOP
  !------------------------------------------------------------------------

contains


  ! ================================================================


  subroutine Restart_Init( status )

    use GO,           only : Init, Done, ReadRc
    use global_data,  only : rcF
    use global_data,  only : outdir
    use TM5_Geometry, only : lli

    implicit none

    ! --- in/out -------------------------------

    integer, intent(out)    ::  status

    ! --- const --------------------------------

    character(len=*), parameter  ::  rname = 'Restart_Init'

    ! --- local --------------------------------
    character(len=MAX_RCKEY_LEN)      :: default_restart_dir

    ! --- begin --------------------------------

    !
    ! read settings from rcfile
    !

    ! open:
    !call Init( rcF, rcfile, status )
    !IF_NOTOK_RETURN(status=1)

    ! write restart files at all ?
    call ReadRc( rcF, 'restart.write', rs_write, status, default=.false. )
    IF_ERROR_RETURN(status=1)

    ! further settings ...
    if ( rs_write ) then

      ! output directory:
      write(default_restart_dir, '(a,"/restart")') trim(outdir)
      call ReadRc( rcF, 'restart.write.dir', rs_write_dir, status, default=default_restart_dir )
      IF_ERROR_RETURN(status=1)

      ! extra restart files ?
      call ReadRc( rcF, 'restart.write.extra', rs_write_extra, status, default=.false. )
      IF_ERROR_RETURN(status=1)
      if ( rs_write_extra ) then
        call ReadRc( rcF, 'restart.write.extra.hour', rs_write_extra_hour, status, default=0 )
        IF_ERROR_RETURN(status=1)
        call ReadRc( rcF, 'restart.write.extra.dhour', rs_write_extra_dhour, status, default=24 )
        IF_ERROR_RETURN(status=1)
      end if

    end if  ! write restart files

    ! close:
    !call Done( rcF, status )
    !IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0

  end subroutine Restart_Init


  ! ***


  subroutine Restart_Done( status )

    implicit none

    ! --- in/out -------------------------------

    integer, intent(out)    ::  status

    ! --- const --------------------------------

    character(len=*), parameter  ::  rname = 'Restart_Done'

    ! --- begin --------------------------------

    ! nothing to be done ...

    ! ok
    status = 0

  end subroutine Restart_Done


  ! ***


  subroutine Restart_Save( status, extra, isfirst )

    use dims, only : idate

    implicit none

    ! --- in/out -------------------------------

    integer, intent(out)            ::  status
    logical, intent(in), optional   ::  extra
    logical, intent(in), optional   ::  isfirst

    ! --- const --------------------------------

    character(len=*), parameter  ::  rname = 'Restart_Save'

    ! --- local --------------------------------

    logical  ::  is_extra

    ! --- begin --------------------------------

    ! options ...
    is_extra = .false.
    if ( present(extra) ) is_extra = extra

    ! write restart files at all ?
    if ( rs_write ) then

      ! end or extra ?
      if ( is_extra ) then

        ! save extra restart files ?
        if ( rs_write_extra ) then

          ! every hour+n*dhour only :
          if (  modulo( idate(4) - rs_write_extra_hour, rs_write_extra_dhour ) == 0 .and. &
                all( idate(5:6) == 0 ) ) then

            ! write restart file for this time:
            call Restart_Write( status, isfirst=isfirst )
            IF_NOTOK_RETURN(status=1)

          end if  ! for this hour

        end if   ! extra restart files ?

      else

        ! write restart file :
        call Restart_Write( status, isfirst=isfirst )
        IF_NOTOK_RETURN(status=1)

      end if  ! not extra

    end if  ! write at all

    ! ok
    status = 0

  end subroutine Restart_Save


  ! ***

  subroutine Restart_FileName( fname, status, key, dir, isfirst )

    use dims,           only : idate
    use global_data,    only : outdir

    implicit none

    ! --- in/out -------------------------------

    character(len=*), intent(out)   ::  fname
    integer, intent(out)            ::  status

    character(len=*), intent(in), optional  ::  dir
    character(len=*), intent(in), optional  ::  key
    logical, intent(in), optional           ::  isfirst

    ! --- local -------------------------------

    character(len=*), parameter  ::  rname = 'Restart_FileName'
    character(len=MAX_FILENAME_LEN)           ::  adir
    character(len=32)            ::  akey

    ! --- begin --------------------------------

    ! destination directory:
    adir = trim(outdir)
    if ( present(dir) ) adir = trim(dir)

    ! extra key, for example '_x' to denote that
    ! a restart file was dumped after process 'x':
    akey = ''
    if ( present(key) )  akey = trim(key)

    ! if this is the initial time, add an extra key to avoid
    ! that the restart file for this hour from the previous
    ! run is overwritten:
    if ( present(isfirst) ) then
      if ( isfirst ) akey = trim(akey)//'_initial'
    end if

    ! write filename:
    write (fname,'(a,"/TM5_restart_",i4.4,4i2.2,a,".nc")') &
                trim(adir), idate(1:5), trim(akey)

    ! ok
    status = 0

  end subroutine Restart_FileName

  !--------------------------------------------------------------------------
  !                    TM5                                                  !
  !--------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE:  Restart_Write
  !
  ! !DESCRIPTION:
  !\\
  !\\
  ! !INTERFACE:
  !
  subroutine Restart_Write( status, key, isfirst )

    use GO          , only : Get
    use dims        , only : nregions, at, bt, region_name
    use chem_param  , only : ntracet, ntrace_chem, ntrace, names, tracer_name_len
    use global_data , only : mass_dat
!#ifdef with_tendencies
    !use tm5_tendency, only : plc_ntr, plc_trname
    !use tm5_tendency, only : plc_npr, plc_prname
    !use tracer_data , only : plc_dat
!#endif
    use TM5_Geometry, only : lli, levi
    use meteodata   , only : sp_dat, phlb_dat, m_dat
    use misctools   , only : check_dir
    use file_netcdf
    use netcdf,       only : nf90_put_att, NF90_GLOBAL

    implicit none
    !
    ! !OUTPUT PARAMETERS:
    !
    integer,          intent(out)           ::  status
    !
    ! !INPUT PARAMETERS:
    !
    character(len=*), intent(in), optional  ::  key
    logical,          intent(in), optional  ::  isfirst
    !
    ! !REVISION HISTORY:
    !    8 Apr 2011 - P. Le Sager - use IF_NOTOK_MDF macro
    !      Aug 2015 - S. Basu - adapt for 4DVAR, use file_netcdf
    !
    ! !REMARKS:
    !
    !EOP
    !------------------------------------------------------------------------
    !BOC

    ! --- local --------------------------------

    character(len=*), parameter  ::  rname = 'Restart_Write'

    integer               ::  imr, jmr, lmr
    integer               ::  region
    character(len=MAX_FILENAME_LEN)    ::  fname
    integer               ::  ftype
    integer               ::  ncid, grpid, io_status, itr
    integer               ::  dimid_lon, dimid_lat, dimid_lev, dimid_hlev
    integer               ::  dimid_lon_sfc, dimid_lat_sfc
    integer               ::  dimid_trace, dimid_trace_transp, dimid_trace_chem
    integer               ::  dimid_name
    integer               ::  varid
    integer               ::  varid_sp, varid_ph, varid_m, varid_at, varid_bt
    integer               ::  varid_names
    integer               ::  varid_rm
#ifdef slopes
    integer               ::  varid_rxm, varid_rym, varid_rzm
#endif
    integer               ::  varid_rmc
#ifdef with_tendencies
    integer               ::  varid_plc(plc_ntr,plc_npr)
    integer               ::  itr, ipr
    integer               ::  time6(6)
#endif
    integer               ::  rtype
    character(len=ntrace*(tracer_name_len+1))    :: attr_value

    ! --- begin --------------------------------

    write (gol,'("write restart file(s) ...")'); call goPr

    ! name of restart file
    call Restart_FileName( fname, status, key=key, dir=rs_write_dir, isfirst=isfirst )
    IF_NOTOK_RETURN(status=1)
    call check_dir(fname)

    write (gol,'("  destination : ",a)') trim(fname); call goPr

    ! o open netcdf file
    ncid = nc_open(trim(fname), 'c', status)
    IF_NOTOK_RETURN(status=1)

    nc_variables_deflate = .true.
    nc_deflate_level = 9

    call nc_create_dim(ncid, 'lev', levi%nlev)
    call nc_create_dim(ncid, 'hlev', levi%nlev+1)
    call nc_create_dim(ncid, 'trace_transp', ntracet)
    if (ntrace_chem > 0) call nc_create_dim(ncid, 'trace_chem', ntrace_chem)
    call nc_create_dim(ncid, 'trace', ntrace)

    lmr = levi%nlev

    ! at and bt coefficients
    call nc_dump_var(ncid, 'at', (/'hlev'/), at(1:lmr+1), (/'long_name'/), (/'hybrid grid a_t coefficient'/))
    call nc_dump_var(ncid, 'bt', (/'hlev'/), bt(1:lmr+1), (/'long+name'/), (/'hybrid grid b_t coefficient'/))

    ! tracer names
    attr_value = names(1)
    do itr = 2, ntrace
        write(attr_value, '(a,",",a)') trim(attr_value), names(itr)
    end do
    io_status = nf90_put_att(ncid, NF90_GLOBAL, 'tracer_names', trim(attr_value))

    ! loop over regions:
    do region = 1, nregions

      grpid = nc_create_group(ncid, region_name(region))

      ! grid size
      imr = lli(region)%nlon
      jmr = lli(region)%nlat

      ! o define dimensions
      call nc_create_dim(grpid, 'lon', imr+4)
      call nc_create_dim(grpid, 'lat', jmr+4)

      ! o write variables
      ! surface pressure
      call nc_dump_var(grpid, 'sp', (/'lon','lat'/), sp_dat(region)%data(-1:imr+2, -1:jmr+2, 1), (/'long_name', 'unit     '/), &
        (/'surface pressure', 'Pascals         '/))
      ! half level pressure
      call nc_dump_var(grpid, 'ph', (/'lon ','lat ','hlev'/), phlb_dat(region)%data(-1:imr+2, -1:jmr+2, 1:lmr+1), &
        (/'long_name', 'unit     '/), (/'half level pressure', 'Pascals            '/))
      ! air mass
      call nc_dump_var(grpid, 'm', (/'lon','lat','lev'/), m_dat(region)%data(-1:imr+2, -1:jmr+2, 1:lmr), (/'long_name', 'unit     '/), &
        (/'air mass', 'Kg      '/))
      ! tracer mass
      call nc_dump_var(grpid, 'rm', (/'lon         ','lat         ','lev         ','trace_transp'/), &
        mass_dat(region)%rm_t(-1:imr+2, -1:jmr+2, 1:lmr, 1:ntracet), (/'long_name', 'unit     '/), &
        (/'transported tracer mass', 'Kg                     '/))
      ! tracer mass slopes
#ifdef slopes
      call nc_dump_var(grpid, 'rxm', (/'lon         ','lat         ','lev         ','trace_transp'/), &
        mass_dat(region)%rxm_t(-1:imr+2, -1:jmr+2, 1:lmr, 1:ntracet), (/'long_name', 'unit     '/), &
        (/'tracer mass slope along x', 'Kg/(half cell)           '/))
      call nc_dump_var(grpid, 'rym', (/'lon         ','lat         ','lev         ','trace_transp'/), &
        mass_dat(region)%rym_t(-1:imr+2, -1:jmr+2, 1:lmr, 1:ntracet), (/'long_name', 'unit     '/), &
        (/'tracer mass slope along y', 'Kg/(half cell)           '/))
      call nc_dump_var(grpid, 'rzm', (/'lon         ','lat         ','lev         ','trace_transp'/), &
        mass_dat(region)%rzm_t(-1:imr+2, -1:jmr+2, 1:lmr, 1:ntracet), (/'long_name', 'unit     '/), &
        (/'tracer mass slope along z', 'Kg/(half cell)           '/))
#endif
      ! non-transported tracers:
      if (ntrace_chem > 0) then
        call nc_dump_var(grpid, 'rmc', (/'lon       ','lat       ','lev       ','trace_chem'/), &
            mass_dat(region)%rm_k(-1:imr+2, -1:jmr+2, 1:lmr, ntracet+1:ntracet+ntrace_chem), (/'long_name', 'unit     '/), &
            (/'non-transported tracer mass', 'Kg                         '/))
      end if

!#ifdef with_tendencies
      !! production, loss, and concentration:
      !do itr = 1, plc_ntr
        !do ipr = 1, plc_npr
          !! define netcdf variable:
          !call MDF_Def_Var( ncid, trim(plc_trname(itr))//'_'//trim(plc_prname(ipr)), rtype, &
                                   !(/dimid_lon,dimid_lat,dimid_lev/), varid, status )
          !IF_NOTOK_MDF(fid=ncid)
          !call MDF_Put_Att( ncid, varid, 'long_name', 'chemical tendency', status )
          !IF_NOTOK_MDF(fid=ncid)
          !call MDF_Put_Att( ncid, varid, 'unit', trim(plc_dat(region,itr,ipr)%unit), status )
          !IF_NOTOK_MDF(fid=ncid)
          !! extract time as 6 integers:
          !call Get( plc_dat(region,itr,ipr)%t, time6=time6 )
          !! add time attribute:
          !call MDF_Put_Att( ncid, varid, 'time', time6, status )
          !IF_NOTOK_MDF(fid=ncid)
          !! store variable id:
          !varid_plc(itr,ipr) = varid
        !end do
      !end do
!#endif

!#ifdef with_tendencies
!#if defined MPI && defined with_netcdf4_par
      !! set independent data mode (not all processes are writing):
      !do itr = 1, plc_ntr
        !do ipr = 1, plc_npr
          !call MDF_Var_Par_Access( ncid, varid_plc(itr,ipr) , MDF_INDEPENDENT, status )
          !IF_NOTOK_MDF(fid=ncid)
        !end do
      !end do
!#endif

      !! write production/loss/concentration levels on this pe:
      !if ( lmloc > 0 ) then

        !do itr = 1, plc_ntr
          !do ipr = 1, plc_npr
            !call MDF_Put_Var( ncid, varid_plc(itr,ipr), &
                            !plc_dat(region,itr,ipr)%rm_k(1:imr,1:jmr,1:lmloc), status, &
                            !start=(/1,1,offsetl+1/), count=(/imr,jmr,lmloc/) )
            !IF_NOTOK_MDF(fid=ncid)
          !end do
        !end do

      !end if
!#endif

      !! o close file
      !call MDF_Close( ncid, status )
      !IF_NOTOK_RETURN(status=1)

    end do   ! regions

    nc_variables_deflate = .false.

    call nc_close(ncid)

    !write (gol,'("  ok")'); call goPr

    ! ok
    status = 0

  end subroutine Restart_Write
  !EOC
  !--------------------------------------------------------------------------
  !                    TM5                                                  !
  !--------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE:  Restart_Read
  !
  ! !DESCRIPTION:
  !\\
  !\\
  ! !INTERFACE:
  !
  subroutine Restart_Read( status, region, surface_pressure, pressure, air_mass, tracer_mass, override_istart)
    !
    ! !USES:
    !
    use GO          , only : Init, Done, ReadRc
    use GO          , only : Set
    use GO          , only : goMatchValue
    use dims        , only : nregions, region_name
    use dims        , only : istart, idate, idatei
    use chem_param  , only : ntracet, ntrace_chem, ntrace
    use chem_param  , only : names, tracer_name_len
    use global_data , only : rcF, outdir
    use global_data , only : mass_dat
!#ifdef with_tendencies
    !use tm5_tendency, only : plc_ntr, plc_trname
    !use tm5_tendency, only : plc_npr, plc_prname
    !use tracer_data , only : plc_dat
!#endif
    use TM5_Geometry, only : lli, levi
    use meteodata   , only : sp_dat, phlb_dat, m_dat
    use file_netcdf

    implicit none
    !
    ! !OUTPUT PARAMETERS:
    !
    integer, intent(out)             ::  status
    !
    ! !INPUT PARAMETERS:
    !
    integer, intent(in), optional   ::  region
    logical, intent(in), optional   ::  surface_pressure, pressure, air_mass, tracer_mass
    logical, intent(in), optional   ::  override_istart
    !
    ! !REVISION HISTORY:
    !    8 Apr 2011 - P. Le Sager - use IF_NOTOK_MDF macro
    !   28 Apr 2011 - P. Le Sager - Check on tracer availability in restart file.
    !                             - Allows for more tracers in restart file than needed
    !   10 May 2011 - P. Le Sager - Added deallocate statement to work with zoom regions
    !      Aug 2015 - S. Basu - Adapted to 4DVAR, use NC instead of HDF
    !
    ! !REMARKS:
    !
    !EOP
    !------------------------------------------------------------------------
    !BOC

    ! --- const --------------------------------

    character(len=*), parameter  ::  rname = mname//'/Restart_Read'

    ! --- local --------------------------------

    character(len=MAX_RCKEY_LEN)              ::  rs_read_dir, default_restart_dir
    logical                         ::  exist
    logical                         ::  do_sp, do_ph, do_m, do_sflux, do_rm, do_plc, do_megan, do_pulse
    integer                         ::  imr, jmr, lmr
    integer                         ::  n
    character(len=MAX_FILENAME_LEN)              ::  fname
    integer                         ::  ncid, grpid
    integer                         ::  varid_sp
    integer                         ::  varid_ph
    integer                         ::  varid_m
    !integer                         ::  varid_slhf, varid_sshf
    integer                         ::  varid_names
    integer                         ::  iname
    character(len=tracer_name_len), allocatable ::  values_names(:)
    integer                         ::  itr, itr_loc, itr_file
    integer                         ::  varid_rm
    integer                         ::  ntracet_restart, dimid
    integer                         ::  shp(2)
#ifdef slopes
    integer                         ::  varid_rxm, varid_rym, varid_rzm
#endif
    integer                         ::  varid_rmc
#ifdef with_tendencies
    integer                         ::  varid_plc(plc_ntr,plc_npr)
    integer                         ::  itr, ipr
    integer                         ::  time6(6)
#endif
    real, allocatable               ::  dummy_1d(:), dummy_2d(:,:), dummy_3d(:,:,:), dummy_4d(:,:,:,:)
    integer                         ::  im_lo, im_hi, jm_lo, jm_hi, file_lat, file_lon

    ! --- begin --------------------------------

    write (gol,'("read restart file ...")'); call goPr

    ! correct istart ?
    if ( istart /= 33 .and. .not. present(override_istart) ) then
      write (gol,'("  skip; istart not 33 but ",i3)') istart; call goPr
      status=0; return
    end if

    ! init time ?
    if ( any( idate /= idatei ) ) then
      write (gol,'("  skip; idate not idatei but ",i4,5i2.2)') idate; call goPr
      status=0; return
    end if

    ! input directory:
    write(default_restart_dir, '(a,"/restart")') trim(outdir)
    call ReadRc( rcF, 'restart.read.dir', rs_read_dir, status, default=default_restart_dir )
    IF_ERROR_RETURN(status=1)

    ! data sets:
    do_sp      = .false.  ;  if ( present(surface_pressure ) ) do_sp      = surface_pressure
    do_ph      = .false.  ;  if ( present(pressure         ) ) do_ph      = pressure
    do_m       = .false.  ;  if ( present(air_mass         ) ) do_m       = air_mass
    do_rm      = .false.  ;  if ( present(tracer_mass      ) ) do_rm      = tracer_mass
    !do_plc     = .false.  ;  if ( present(tendencies       ) ) do_plc     = tendencies

    ! name of restart file
    call Restart_FileName( fname, status, dir=trim(rs_read_dir) )
    IF_NOTOK_RETURN(status=1)
    write (gol,'("  restore from ",a)') trim(fname); call goPr

    ! test ...
    inquire( file=fname, exist=exist )
    if ( .not. exist ) then
        write (gol,'("restart file not found : ",a)') trim(fname); call goErr
        TRACEBACK; status=1; return
    end if

    ncid = nc_open(trim(fname), 'r', status)
    IF_NOTOK_RETURN(status=1)

    if ( do_rm ) then
        ! check the number of tracers
        if (ntrace /= nc_get_dim(ncid, 'trace')) then
            write(gol, '("Number of tracers ", i2, " in model but ", i2, " in restart file")') ntrace, nc_get_dim(ncid, 'trace')
            call goErr
            status = 1 ; IF_NOTOK_RETURN(status=1)
        end if
        if (ntracet /= nc_get_dim(ncid, 'trace_transp')) then
            write(gol, '("Number of transported tracers ", i2, " in model but ", i2, " in restart file")') ntracet, nc_get_dim(ncid, 'trace_transp')
            call goErr
            status = 1 ; IF_NOTOK_RETURN(status=1)
        end if

    end if

    ! loop over ns
    do n = 1, nregions

      if (present(region)) then
        if (n /= region) cycle
      end if

      grpid = nc_get_group(ncid, region_name(n))

      ! grid size
      imr = lli(n)%nlon
      jmr = lli(n)%nlat
      lmr = levi%nlev

      ! o read variables
      ! For old restart files, padded masses and pressures were not stored. If we're reading one of those files, we need
      ! to read in just the 'core' of each region. For newer restart files, we need to read in padded masses.
      file_lat = nc_get_dim(grpid, 'lat')
      ! Stupid fortran. I can't use a 'select case (file_lat)' because 'case(jmr)' is not allowed. Expressions within
      ! 'case' statements must be constants.
      if (file_lat == jmr) then
        jm_lo = 1
        jm_hi = jmr
        write(gol, '(a, ": WARNING - expecting padded masses in restart file, got core masses for region ", a)') &
            rname, region_name(n) ; call goPr
        write(gol, '(a, ": WARNING - latitude dimension ", i3, " instead of ", i3)') &
            rname, file_lat, jmr+4 ; call goPr
      else if (file_lat == jmr+4) then
        jm_lo = -1
        jm_hi = jmr+2
      else
        write(gol, '(a, ": dimension lat for region ", a, " in file ", a, " is ", i3, ", which is neither ", i3, " nor ", i3)') &
            rname, region_name(n), trim(fname), file_lat, jmr, jmr+4 ; call goErr
            status = 1 ; IF_NOTOK_RETURN(status=1)
      end if

      file_lon = nc_get_dim(grpid, 'lon')
      if (file_lon == imr) then
        im_lo = 1
        im_hi = imr
        write(gol, '(a, ": WARNING - expecting padded masses in restart file, got core masses for region ", a)') &
            rname, region_name(n) ; call goPr
        write(gol, '(a, ": WARNING - longitude dimension ", i3, " instead of ", i3)') &
            rname, file_lon, imr+4 ; call goPr
      else if (file_lon == imr+4) then
        im_lo = -1
        im_hi = imr+2
      else
        write(gol, '(a, ": dimension lon for region ", a, " in file ", a, " is ", i3, ", which is neither ", i3, " nor ", i3)') &
            rname, region_name(n), trim(fname), file_lon, imr, imr+4 ; call goErr
            status = 1 ; IF_NOTOK_RETURN(status=1)
      end if

      ! surface pressure
      if ( do_sp ) then
        sp_dat(n)%data = 0.0
        dummy_2d = nc_read_var(grpid, 'sp', status)
        IF_NOTOK_RETURN(status=1)
        sp_dat(n)%data(im_lo:im_hi, jm_lo:jm_hi, 1) = dummy_2d
        deallocate(dummy_2d)
      end if

      ! half level pressure
      if ( do_ph ) then
        dummy_3d = nc_read_var(grpid, 'ph', status)
        IF_NOTOK_RETURN(status=1)
        phlb_dat(n)%data(im_lo:im_hi, jm_lo:jm_hi, 1:lmr+1) = dummy_3d
        deallocate(dummy_3d)
      end if

      ! air mass
      if ( do_m ) then
        dummy_3d = nc_read_var(grpid, 'm', status)
        IF_NOTOK_RETURN(status=1)
        m_dat(n)%data(im_lo:im_hi, jm_lo:jm_hi, 1:lmr) = dummy_3d
        deallocate(dummy_3d)
      end if

      ! tracer mass
      if ( do_rm ) then
        dummy_4d = nc_read_var(grpid, 'rm', status)
        IF_NOTOK_RETURN(status=1)
        mass_dat(n)%rm_t(im_lo:im_hi, jm_lo:jm_hi, 1:lmr, 1:ntracet) = dummy_4d
        deallocate(dummy_4d)

#ifdef slopes
        dummy_4d = nc_read_var(grpid, 'rxm', status)
        IF_NOTOK_RETURN(status=1)
        mass_dat(n)%rxm_t(im_lo:im_hi, jm_lo:jm_hi, 1:lmr, 1:ntracet) = dummy_4d
        deallocate(dummy_4d)

        dummy_4d = nc_read_var(grpid, 'rym', status)
        IF_NOTOK_RETURN(status=1)
        mass_dat(n)%rym_t(im_lo:im_hi, jm_lo:jm_hi, 1:lmr, 1:ntracet) = dummy_4d
        deallocate(dummy_4d)

        dummy_4d = nc_read_var(grpid, 'rzm', status)
        IF_NOTOK_RETURN(status=1)
        mass_dat(n)%rzm_t(im_lo:im_hi, jm_lo:jm_hi, 1:lmr, 1:ntracet) = dummy_4d
        deallocate(dummy_4d)
#endif
        if (ntrace_chem > 0) then
            dummy_4d = nc_read_var(grpid, 'rmc', status)
            IF_NOTOK_RETURN(status=1)
            mass_dat(n)%rm_k(im_lo:im_hi, jm_lo:jm_hi, 1:lmr, ntracet+1:ntracet+ntrace_chem) = dummy_4d
            deallocate(dummy_4d)
        end if
      end if ! tracer mass

    end do   ! regions

    call nc_close(ncid)

    ! ok
    status = 0

  end subroutine Restart_Read



end module Restart
