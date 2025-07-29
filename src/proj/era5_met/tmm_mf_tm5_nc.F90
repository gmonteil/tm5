!###############################################################################
!
! Input/output of meteofiles produced by TM5 : NetCDF version.
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module tmm_mf_tm5_nc

  use GO, only : gol, goErr, goPr
  use GO , only : TDate
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN, DUMMY_STR_LEN, LONG_STR_LEN

  !use MDF, only : MDF_NETCDF
  use MDF, only : MDF_NETCDF4
  use Grid, only : TllGridInfo, TLevelInfo

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  TMeteoFile_tm5_nc
  public  ::  TMM_MF_TM5_NC_Init, TMM_MF_TM5_NC_Done
  public  ::  Init, Done
  public  ::  Get
  public  ::  ReadRecord
  public  ::  WriteRecord


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf_tm5_nc'

  ! ~~~ output keys and defaults

  ! current format version
  character(len=*), parameter  ::  output_format = 'tm5-nc'

  ! extension and type:
  !integer, parameter           ::  output_type = MDF_NETCDF
  integer, parameter           ::  output_type = MDF_NETCDF4
  character(len=*), parameter  ::  output_ext  = '.nc'

  ! standard timevalue is seconds since ...
  integer, parameter  ::  since_time6(6) = (/1900,01,01,00,00,00/)


  !--- type --------------------------------------

  ! single field, all time records:
  type T_Cache_Field
    ! parameter name:
    character(len=256)     ::  paramkey
    ! data:
    real, pointer          ::  field(:,:,:,:)  ! (nx,ny,nz,nt)
    real, pointer          ::  sp   (:,:  ,:)  ! (nx,ny   ,nt)
  end type T_Cache_Field
  
  ! *

  type TMeteoFile_tm5_nc
    ! input/output ?
    character(len=1)       ::  io
    !
    ! field collection
    !
    character(len=MAX_FILENAME_LEN)     ::  fname
    character(len=MAX_RCKEY_LEN)     ::  paramkeys      ! -aa-bb-cc-
    integer                ::  nparam
    character(len=64)      ::  tres
    type(TDate)            ::  trange(2)
    logical                ::  is_aver
    integer                ::  rnk
    !
    ! file
    !
    integer                ::  hid
    integer                ::  dimid_nv
    integer                ::  dimid_lon, varid_lon, varid_lon_bounds
    integer                ::  dimid_lat, varid_lat, varid_lat_bounds
    integer                ::  dimid_lonb, varid_lonb
    integer                ::  dimid_latb, varid_latb
    integer                ::  dimid_lev, varid_lev, varid_ap, varid_b, varid_ap_bounds, varid_b_bounds
    integer                ::  dimid_lev_u, varid_lev_u
    integer                ::  dimid_lev_v, varid_lev_v
    integer                ::  dimid_time, varid_time, varid_time_bounds, varid_reftime
    integer                ::  dimid_timeval, varid_timevalues, varid_timevalues_bounds, varid_reftimevalues
    integer                ::  varid_cell_area
    integer                ::  varid_ps, varid_ps_u, varid_ps_v
    !integer                ::  dimid_tm5_lm, dimid_tm5_lmb, varid_tm5_at, varid_tm5_bt
    integer, allocatable              ::  varid_param(:)
    character(len=16), allocatable    ::  varname_param(:)      ! (/'aa','bb','cc'/)
    character(len=MAX_RCKEY_LEN), allocatable   ::  cfname_param(:)       ! from CF table
    character(len=64), allocatable    ::  cfunit_param(:)       ! from CF table, following UDUnits
    integer, allocatable              ::  itrec_param(:)
    !
    ! input chache
    !
    ! current cached file:
    character(len=1024)               ::  cache_fname
    ! grid:
    type(TllGridInfo)                 ::  cache_lli
    ! levels:
    type(TLevelInfo)                  ::  cache_levi
    ! time:
    integer, pointer                  ::  cache_timevalues_bounds(:,:,:)  ! (2,6,ntime)
    ! cached params:
    type(T_Cache_Field), allocatable  ::  cache_field(:)  ! (nparam)
    !
    ! output
    !
    logical                ::  output_initialised
    integer                ::  output_nrec
    integer                ::  output_ntrec
    !
    ! adhoc ...
    integer                ::  fixyear
  end type TMeteoFile_tm5_nc


  ! --- interfaces -------------------------------

  interface Init
    module procedure mf_Init
  end interface

  interface Done
    module procedure mf_Done
  end interface

  interface Get
    module procedure mf_Get
  end interface

  interface ReadRecord
    module procedure mf_ReadRecord
  end interface

  interface WriteRecord
    module procedure mf_WriteRecord_2d
    module procedure mf_WriteRecord_3d
  end interface


  ! --- var --------------------------------------

  ! timer id's:
  integer           ::  itim_readrecord


contains


  ! ==============================================================


  subroutine TMM_MF_TM5_NC_Init( rcf, status )

    use GO, only : TrcFile
    use GO, only : GO_Timer_Def
    use TMM_CF, only : TMM_CF_Init

    ! --- in/out ---------------------------------

    type(TRcFile), intent(in)           ::  rcf
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TMM_MF_TM5_NC_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! init cf table and udunits package:
    call TMM_CF_Init( rcf, status )
    IF_NOTOK_RETURN(status=1)

    ! define timers:
    call GO_Timer_Def( itim_readrecord, 'tmm tm5-nc readrecord', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TMM_MF_TM5_NC_Init


  ! ***


  subroutine TMM_MF_TM5_NC_Done( status )

    use TMM_CF, only : TMM_CF_Done

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TMM_MF_TM5_NC_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! done with standard names table and udunits package:
    call TMM_CF_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TMM_MF_TM5_NC_Done


  ! ==============================================================


  subroutine mf_Init( mf, io, dir, archivekeys, paramkey, &
                               tref, t1, t2, status )

    use GO, only : TDate, Set, Get, NewDate, AnyDate, IsAnyDate
    use GO, only : rTotal, operator(-), operator(>=)
    use GO, only : goVarValue, goWriteKeyNum, goReplace
    use PArray, only : pa_Init

    ! --- in/out ----------------------------

    type(TMeteoFile_tm5_nc), intent(out)   ::  mf
    character(len=1), intent(in)        ::  io
    character(len=*), intent(in)        ::  dir
    character(len=*), intent(in)        ::  archivekeys
    character(len=*), intent(in)        ::  paramkey
    type(TDate), intent(in)             ::  tref, t1, t2
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Init'

    ! --- local --------------------------------

    character(len=64)     ::  mf_mdir
    character(len=1)      ::  mf_psep, mf_nsep

    character(len=64)     ::  mf_filekey

    character(len=4)      ::  mf_fckey

    type(TDate)           ::  tfile
    integer               ::  ccyy, mm, dd, dh
    type(TDate)           ::  tc

    integer               ::  iparam

    ! --- begin --------------------------------
    

    ! store i/o :
    mf%io = io

    ! default flags:
    mf%is_aver = .false.
    mf%rnk     = 2         ! most 2D fields

    !
    ! extract fields from archivekey :
    !   mdir=ec-fg_3h-ml60-glb3x2;tres=_21p06
    !
    mf%fixyear  = -1
    call goVarValue( archivekeys, ';', 'fixyear', '=', mf%fixyear, status )
    if (status>0) then; TRACEBACK; status=1; return; end if

    mf_mdir  = 'no_mdir'
    call goVarValue( archivekeys, ';', 'mdir', '=', mf_mdir, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    mf%tres  = 'no_tres'
    call goVarValue( archivekeys, ';', 'tres', '=', mf%tres, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    ! path seperation character:
    mf_psep  = '/'
    call goVarValue( archivekeys, ';', 'pathsep', '=', mf_psep, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    ! name seperation character:
    mf_nsep  = '-'
    call goVarValue( archivekeys, ';', 'namesep', '=', mf_nsep, status )
    if (status>0) then; TRACEBACK; status=1; return; end if


    !
    ! main file
    !

    ! * set mf_filekey (uvsp,t,etc) and parmeters:
    select case ( paramkey )
      case ( 'sp', 'tsp' )
        mf_filekey = paramkey
        mf%paramkeys = '-'//trim(paramkey)//'-'
        mf%nparam    = 1
      case ( 'mfu', 'mfv' )
        mf_filekey = 'mfuv'
        mf%paramkeys = '-mfu-mfv-'
        mf%nparam    = 2
        mf%rnk       = 3
      case ( 'mfw' )
        mf_filekey = 'mfw'
        mf%paramkeys = '-mfw-'
        mf%nparam    = 1
        mf%rnk       = 3
      case ( 'T' )
        mf_filekey   = 't'
        mf%paramkeys = '-T-'
        mf%nparam    = 1
        mf%rnk       = 3
      case ( 'Q' )
        mf_filekey   = 'q'
        mf%paramkeys = '-Q-'
        mf%nparam    = 1
        mf%rnk       = 3
      case ( 'CLWC', 'CIWC', 'CC', 'CCO', 'CCU', &
             'clwc', 'ciwc', 'cc', 'cco', 'ccu'  )
        mf_filekey   = 'cld'
        mf%paramkeys = '-CLWC-CIWC-CC-CCO-CCU-'
        mf%nparam    = 5
        mf%rnk       = 3
      case ( 'eu', 'ed', 'du', 'dd' ) ! computed online:  'cloud_base', 'cloud_top', 'cloud_lfs'
        mf_filekey   = 'convec'
        !mf%paramkeys = '-eu-ed-du-dd-cloud_base-cloud_top-cloud_lfs-'
        mf%paramkeys = '-eu-ed-du-dd-'
        mf%nparam    = 4
        mf%rnk       = 3
      ! o constant fields
      case ( 'oro', 'lsm' )
        mf%tres      = 'constant'
        mf_filekey   = trim(paramkey)
        mf%paramkeys = '-'//trim(paramkey)//'-'
        mf%nparam    = 1
      ! o monthly fields:
      case ( 'srols' )
        mf%tres      = 'month'
        mf_filekey   = trim(paramkey)
        mf%paramkeys = '-'//trim(paramkey)//'-'
        mf%nparam    = 1
      ! o vegetation fields
      case ( 'cvl', 'cvh', &
             'tv01', 'tv02', 'tv03', 'tv04', 'tv05', &
             'tv06', 'tv07',         'tv09', 'tv10', &
             'tv11',         'tv13',                 &
             'tv16', 'tv17', 'tv18', 'tv19'         )
        mf_filekey = 'veg'
        mf%paramkeys = '-'
        mf%paramkeys = trim(mf%paramkeys)//'cvl-cvh-'
        mf%paramkeys = trim(mf%paramkeys)//'tv01-tv02-tv03-tv04-tv05-'
        mf%paramkeys = trim(mf%paramkeys)//'tv06-tv07-tv09-tv10-'
        mf%paramkeys = trim(mf%paramkeys)//'tv11-tv13-'
        mf%paramkeys = trim(mf%paramkeys)//'tv16-tv17-tv18-tv19-'
        mf%nparam    = 17
      ! o each surface file in a seperate file:
      case ( 'sr', 'srmer', &
             'swvl1', &
             'albedo', 'lsrh', 'ci', 'g10m', 'u10m', 'v10m', 'sd', &
             'blh', &
             't2m', 'd2m', &
             'sstr', 'src', 'raero', 'ustar', &
             'sst', 'sps', 'skt' )
        mf_filekey = trim(paramkey)
        mf%paramkeys = '-'//trim(paramkey)//'-'
        mf%nparam    = 1
      case ( 'sshf', 'slhf', 'ewss', 'nsss', 'lsp', 'cp', 'sf', &
             'ssr', 'ssrd', 'str', 'strd' )
          mf_filekey   = trim(paramkey)
          mf%paramkeys = '-'//trim(paramkey)//'-'
          mf%nparam    = 1
          mf%is_aver = .true.
      case default
        write (gol,'("unsupported paramkey `",a,"`")') paramkey; call goErr
        TRACEBACK; status=1; return
    end select

    ! eventually replace:
    call goVarValue( archivekeys, ';', 'filekey', '=', mf_filekey, status )
    if (status>0) then; TRACEBACK; status=1; return; end if

    ! storage:
    allocate( mf%varid_param  (mf%nparam) )
    allocate( mf%varname_param(mf%nparam) )
    allocate( mf%cfname_param (mf%nparam) )
    allocate( mf%cfunit_param (mf%nparam) )
    allocate( mf%itrec_param  (mf%nparam) )

    ! convert input times to file name times:
    call GetTime( mf_filekey, mf%tres, tref, t1, t2, status, &
                       tfile=tfile, trange=mf%trange )
    IF_NOTOK_RETURN(status=1)

    ! adhoc: fixed year ?
    if ( mf%fixyear > 0 ) call Set( tfile, year=mf%fixyear )

    ! extract time values:
    call Get( tfile, year=ccyy, month=mm, day=dd )

    ! replace time values in directory name:
    call goReplace( mf_mdir, '<yyyy>', '(i4.4)', ccyy, status )
    IF_NOTOK_RETURN(status=1)
    call goReplace( mf_mdir, '<mm>', '(i2.2)', mm, status )
    IF_NOTOK_RETURN(status=1)

    ! special data set: trap change from fg to fc data:
    if ( mf%tres == '_fg006up4tr3' ) then
      tc = NewDate( 2000, 09, 12 )
      if ( tfile >= tc ) mf%tres = '_fc012up2tr3'
    end if

    ! forecast key: '', 'f1', .., 'f10' ;
    ! no key for constant files (t1 and t2 are anydate)
    ! or month files (t1 is begin of month thus probably < tref)
    mf_fckey = ''
    if ( (.not. IsAnyDate(t1)) .and. (t1 >= tref) ) then
      dh = floor( rTotal( t1 - tref, 'day' ) )
      if ( dh > 0 ) call goWriteKeyNum( mf_fckey, 'f', dh )
    end if

    ! create file name:
    !   dir / ec_od-ml60-T159 - oro.hdf
    !   dir / ec_od-ml60-T159 - T_20000101_fg006up4tr3.hdf
    select case ( mf%tres )
      case ( 'constant' )
        ! filename without date:
        write (mf%fname,'(6a)') trim(dir), mf_psep, trim(mf_mdir), mf_nsep, trim(mf_filekey), trim(output_ext)
      case ( 'month' )
        ! filename without day and forecast key:
        write (mf%fname,'(5a,"_",i4.4,i2.2,a)') &
                trim(dir), mf_psep, trim(mf_mdir), mf_nsep, trim(mf_filekey), ccyy, mm, trim(output_ext)
      case default
        ! filename including date:
        write (mf%fname,'(5a,"_",i4.4,2i2.2,3a)') &
                trim(dir), mf_psep, trim(mf_mdir), mf_nsep, &
                trim(mf_filekey), ccyy, mm, dd, trim(mf_fckey), trim(mf%tres), trim(output_ext)
    end select

    ! in case of output, not initialised yet ...
    mf%output_initialised = .false.

    ! number of expected time records in a file:
    call GetTime( mf_filekey, mf%tres, tref, t1, t2, status, nrec=mf%output_ntrec )
    IF_NOTOK_RETURN(status=1)
    
    ! no chached file content yet:
    mf%cache_fname = ''
    ! no time info yet:
    call pa_Init( mf%cache_timevalues_bounds )
    ! storage:
    allocate( mf%cache_field(mf%nparam), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! loop:
    do iparam = 1, mf%nparam
      ! no data yet:
      call pa_Init( mf%cache_field(iparam)%field )
      call pa_Init( mf%cache_field(iparam)%sp    )
    end do

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    use MDF, only : MDF_Close
    use Grid, only : Done
    use PArray, only : pa_Done

    ! --- in/out ------------------------------------

    type(TMeteoFile_tm5_nc), intent(inout) ::  mf
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Done'

    ! --- local -------------------------------------
    
    integer          ::  iparam

    ! --- begin -------------------------------------
    
    !! info ...
    !write (gol,'("<",a,">")') trim(rname); call goPr

    !! info ..
    !write (gol,'(a,": done for ",a)') rname, trim(mf%fname); call goPr    
    ! clear time info:
    call pa_Done( mf%cache_timevalues_bounds )
    ! clear cache:
    do iparam = 1, mf%nparam
      ! done:
      call pa_Done( mf%cache_field(iparam)%field )
      call pa_Done( mf%cache_field(iparam)%sp    )
    end do ! param
    ! clear:
    deallocate( mf%cache_field, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! grid defined?
    if ( len_trim(mf%cache_fname) > 0 ) then
      ! done with levels:
      call Done( mf%cache_levi, status )
      IF_NOTOK_RETURN(status=1)
      ! done with grid:
      call Done( mf%cache_lli, status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! reset:
    mf%cache_fname = ''

    ! file opend for output ?
    if ( mf%output_initialised ) then
      ! close file:
      call MDF_Close( mf%hid, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! clear:
    deallocate( mf%varid_param   )
    deallocate( mf%varname_param )
    deallocate( mf%cfname_param  )
    deallocate( mf%cfunit_param  )
    deallocate( mf%itrec_param   )

    ! ok
    status = 0

  end subroutine mf_Done



  ! ***


  subroutine mf_Get( mf, status, trange1, trange2, paramkeys, filename )

    use GO, only : TDate

    ! --- in/out ----------------------------

    type(TMeteoFile_tm5_nc), intent(in)          ::  mf
    integer, intent(out)                      ::  status

    type(TDate), intent(out), optional        ::  trange1, trange2
    character(len=*), intent(out), optional   ::  paramkeys
    character(len=*), intent(out), optional   ::  filename

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Get'

    ! --- local --------------------------------

    ! --- begin --------------------------------

    ! time range:
    if ( present(trange1) ) trange1 = mf%trange(1)
    if ( present(trange2) ) trange2 = mf%trange(2)

    ! parameter names:
    if ( present(paramkeys) ) paramkeys = trim(mf%paramkeys)

    ! file name:
    if ( present(filename) ) filename = trim(mf%fname)

    ! ok
    status = 0

  end subroutine mf_Get



  ! ******************************************************************
  ! ***
  ! *** time range, parameters, file names
  ! ***
  ! ******************************************************************

  !
  ! Return time parameters:
  !  o tfile   :  date in filename
  !  o trange  :  time interval covered by fields in file
  !  o nrec    :  number of time records in completed file
  !

  subroutine GetTime( filekey, tres, tref, t1, t2, status, &
                         tfile, trange, nrec )

    use GO, only : TDate, NewDate, AnyDate, Get, Set, wrtgol, IncrDate, IsAnyDate
    use GO,     only : operator(<), operator(+), operator(-), operator(/)
    use Go,     only : Pretty, rTotal
    use dims,   only : okdebug_tmm

    ! --- in/out --------------------------------

    character(len=*), intent(in)            ::  filekey
    character(len=*), intent(in)            ::  tres
    type(TDate), intent(in)                 ::  tref, t1, t2
    integer, intent(out)                    ::  status

    type(TDate), intent(out), optional      ::  tfile
    type(TDate), intent(out), optional      ::  trange(2)
    integer, intent(out), optional          ::  nrec

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GetTime'

    ! --- local --------------------------------

    integer          ::  year, month
    integer          ::  hour1, time6(6)
    integer          ::  dd, hh, step
    logical          ::  interval
    real             ::  dhr

    ! --- begin --------------------------------

    ! set day shift, start hour, and step
    select case ( tres )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tmpp [21,21]
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_21p06', '_21p03', '_av21' )

        ! routine is called with  tref,t1,t2:
        !    t1,t1,t2
        !    t1,any,any    (oro and other constant fields)
        ! thus use tref to construct the file times
        if ( IsAnyDate(t1) ) then
          call Get( tref, hour=hour1 )
          interval = .false.
        else
          call Get( t1, hour=hour1 )
          interval = t1 < t2
        end if

        ! file ccyymmdd contains fields for (21,21];
        ! only uvsp is valid for [21,21] since it contains surface pressure for 21:00
        if ( present(tfile) ) then
          tfile = tref
          call Set( tfile, hour=0, min=0, sec=0 )
          if ( (hour1 > 21) .or. ((interval .or. filekey=='uvsp') .and. hour1==21) ) then
            tfile = tfile + IncrDate(day=1)
          end if
       end if

        ! fields by default valid for (21,21];
        ! only uvsp is valid for [21,21] since it contains surface pressure for 21:00
        if ( present(trange) ) then
          trange(1) = tref
          call Set( trange(1), hour=0, min=0, sec=0 )  !  00:00 today
          if ( (hour1 > 21) .or. ((interval .or. filekey=='uvsp') .and. hour1==21) ) then
            trange(1) = trange(1) + IncrDate(day=1)
          end if
          trange(1) = trange(1) - IncrDate(hour=3)     !  previous 21:00
          trange(2) = trange(1) + IncrDate(day=1)      !  next     21:00
          ! boundary not included in most cases:
          if ( filekey /= 'uvsp' ) trange(1) = trange(1) + IncrDate(mili=1)
        end if

        ! number of records in file:
        if ( present(nrec) ) then
          select case ( tres )
            case ( '_21p06' ) ; nrec = 24/6
            case ( '_21p03' ) ; nrec = 24/3
            case ( '_av21'  ) ; nrec = 24/24
            case default
              write (gol,'("unsupported tres for setting nrec : ",a)') tres; call goErr
              TRACEBACK; status=1; return
          end select
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tm5 constant
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'constant' )

        ! no date in filename ...
        if ( present(tfile) ) tfile = AnyDate()

        ! fields always valid ...
        if ( present(trange) ) then
          trange(1) = AnyDate()
          trange(2) = AnyDate()
        end if

        ! only one output record in constant file:
        if ( present(nrec) ) nrec = 1

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tm5 monthly file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'month' )

        ! file ccyymmdd contains fields for this month:
        if ( present(tfile) ) then
          call Get( t1, year=year, month=month )
          tfile = NewDate( year=year, month=month, day=1 )
        end if

        ! field valid from begin to end of month:
        if ( present(trange) ) then
          call Get( t1, year=year, month=month )
          trange(1) = NewDate( year=year, month=month, day=1, hour=00  )
          month = month + 1
          if ( month > 12 ) then
            year = year + 1
            month = 1
          end if
          trange(2) = NewDate( year=year, month=month, day=1, hour=00  )
        end if

        ! only one output record in month file:
        if ( present(nrec) ) nrec = 1

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tm5 [00,24]
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_00p06', '_00p03', '_an0tr6', '_fg006up4tr3', '_fc012up2tr3', '_00p01' )
        ! file ccyymmdd contains fields for [00,24) :
        if ( present(tfile) ) then
          tfile = t1 + (t2-t1)/2
          call Set( tfile, hour=0, min=0, sec=0 )
        end if

        ! fields valid for [00,24) :
        if ( present(trange) ) then
          trange(1) = t1 + (t2-t1)/2
          call Set( trange(1), hour=0, min=0, sec=0 )  !  00:00 today
          trange(2) = trange(1) + IncrDate(hour=24) - IncrDate(mili=1)
        end if

        ! number of records in file:
        if ( present(nrec) ) then
          select case ( tres )
            case ( '_00p06'  ) ; nrec = 24/6
            case ( '_an0tr6' ) ; nrec = 24/6
            case ( '_00p03', '_fg006up4tr3', '_fc012up2tr3' )
              ! by default: 3 hourly files
              ! for forecasts after 12+72, only 6 hourly available:
              !  f0  [ 00, 24)  : 00 03 06 09 12 15 18 21   :  nrec=8
              !  f1  [ 24, 48)  : 00 03 06 09 12 15 18 21   :  nrec=8
              !  f2  [ 48, 72)  : 00 03 06 09 12 15 18 21   :  nrec=8
              !  f3  [ 72, 96)  : 00 03 06 09 12    18      :  nrec=6
              !  f4  [ 96,120)  : 00    06    12    18      :  nrec=4
              !   :
              !  f9  [192,216)  : 00    06    12    18      :  nrec=4
              !  f10 [216,240)  : 00    06    12            :  nrec=3
              dhr = rTotal( t1 - tref, 'hour' )
              if ( dhr < 72.0 ) then
                nrec = 8
              else if ( dhr < 96 ) then
                nrec = 6
              else if ( dhr < 216 ) then
                nrec = 4
              else
                nrec = 3
              end if
            case ( '_00p01'       ) ; nrec = 24/1 ! TODO: older versions of this file have "24/24". Doesn't make much sense to me, but we need to check the impact
            case default
              write (gol,'("unsupported tres for setting nrec : ",a)') tres; call goErr
              TRACEBACK; status=1; return
          end select
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ???
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported time resolution key:")'); call goErr
        write (gol,'("  ",a)') trim(tres); call goErr
        TRACEBACK; status=1; return

    end select

    ! Debug
    if (present(tfile) .and. okdebug_tmm) then
        write(gol,'(a,": for filekey ", a, " and tres ", a, ", t1 = ", a, ", t2 = ", a, ", calculated tfile = ", a)') &
            rname, trim(filekey), trim(tres), trim(Pretty(t1)), trim(Pretty(t2)), trim(Pretty(tfile))
            call goPr
    end if
    ! End debug

    ! ok
    status = 0

  end subroutine GetTime


  ! ******************************************************************
  ! ***
  ! *** input
  ! ***
  ! ******************************************************************


  !
  ! initialiase grid info from sds
  !

  subroutine lli_Init_mf( lli, nuv, mf, status )

    use Grid, only : TllGridInfo, Init
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var

    ! --- in/out ----------------------------------

    type(TllGridInfo), intent(out)          ::  lli
    character(len=1), intent(in)            ::  nuv
    type(TMeteoFile_tm5_nc), intent(in)     ::  mf
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/lli_Init_mf'

    ! --- local -----------------------------------

    integer             ::  dimid, varid
    real, allocatable   ::  values(:)
    real                ::  lon_deg, dlon_deg
    integer             ::  nlon
    real                ::  lat_deg, dlat_deg
    integer             ::  nlat

    ! --- begin ------------------------------------

    ! number of longitudes:
    call MDF_Inq_DimID( mf%hid, 'lon', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( mf%hid, dimid, status, length=nlon )
    IF_NOTOK_RETURN(status=1)
    ! lon axis:
    allocate( values(nlon) )
    ! read:
    call MDF_Inq_VarID( mf%hid, 'lon', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( mf%hid, varid, values, status )
    IF_NOTOK_RETURN(status=1)
    ! extract:
    lon_deg  = values(1)
    dlon_deg = values(2) - values(1)
    ! clear:
    deallocate( values )

    ! number of latgitudes:
    call MDF_Inq_DimID( mf%hid, 'lat', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( mf%hid, dimid, status, length=nlat )
    IF_NOTOK_RETURN(status=1)
    ! lat axis:
    allocate( values(nlat) )
    ! read:
    call MDF_Inq_VarID( mf%hid, 'lat', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( mf%hid, varid, values, status )
    IF_NOTOK_RETURN(status=1)
    ! extract:
    lat_deg  = values(1)
    dlat_deg = values(2) - values(1)
    ! clear:
    deallocate( values )

    ! define grid:
    call Init( lli, lon_deg, dlon_deg, nlon, &
                    lat_deg, dlat_deg, nlat, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine lli_Init_mf


  ! ***


  !
  ! initialiase level info from sds
  !

  subroutine levi_Init_mf( levi, mf, status )

    use Grid, only : TLevelInfo, Init
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var
    use MDF, only : MDF_Get_Att, MDF_GLOBAL

    ! --- in/out ----------------------------------

    type(TLevelInfo), intent(out)           ::  levi
    type(TMeteoFile_tm5_nc), intent(in)     ::  mf
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Init_mf'

    ! --- local -----------------------------------

    integer                ::  dimid, varid
    integer                ::  lev
    character(len=1)       ::  nw
    integer                ::  lm
    real, allocatable      ::  at(:), bt(:)
    real, allocatable      ::  ap_bounds(:,:), b_bounds(:,:)

    ! --- begin ------------------------------------

    ! 2D or 3D ?
    select case ( mf%rnk )

      case ( 2 )

        ! set dummy values ...
        call Init( levi, 1, (/0.0,0.0/), (/0.0,0.0/), status )
        IF_NOTOK_RETURN(status=1)

      case ( 3 )

        ! read number of levels from global attribute:
        call MDF_Inq_DimID( mf%hid, 'lev', dimid, status )
        IF_NOTOK_RETURN(status=1)
        call MDF_Inquire_Dimension( mf%hid, dimid, status, length=lev )
        IF_NOTOK_RETURN(status=1)

        ! layers or levels ?
        call MDF_Get_Att( mf%hid, MDF_GLOBAL, 'nw', nw, status )
        IF_NOTOK_RETURN(status=1)

        ! extract:
        select case ( nw )
          !~ layers
          case ( 'n' )!, '*' )
            ! layers is lev ...
            lm = lev
            ! storage:
            allocate( at(lm+1), bt(lm+1) )
            allocate( ap_bounds(2,lev), b_bounds(2,lev) )
            ! extract hybride coeff
            call MDF_Inq_VarID( mf%hid, 'ap_bounds', varid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Get_Var( mf%hid, varid, ap_bounds, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Inq_VarID( mf%hid, 'b_bounds', varid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Get_Var( mf%hid, varid, b_bounds, status )
            IF_NOTOK_RETURN(status=1)
            ! extract values:
            at = (/ ap_bounds(1,1), ap_bounds(2,:) /)
            bt = (/  b_bounds(1,1),  b_bounds(2,:) /)
            ! clear:
            deallocate( ap_bounds, b_bounds )
          !~ half levels
          case ( 'w' )
            ! layers is one less than lev ...
            lm = lev-1
            ! storage:
            allocate( at(lm+1), bt(lm+1) )
            ! extract hybride coeff
            call MDF_Inq_VarID( mf%hid, 'ap', varid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Get_Var( mf%hid, varid, at, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Inq_VarID( mf%hid, 'b', varid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Get_Var( mf%hid, varid, bt, status )
            IF_NOTOK_RETURN(status=1)
          !~ unknown ...
          case default
            write (gol,'("found unsupported value of nw attribute : ",a)') trim(nw); call goErr
            write (gol,'("  file : ",a)') trim(mf%fname); call goErr
            TRACEBACK; status=1; return
        end select

        ! fill ...
        call Init( levi, lm, at, bt, status )
        IF_NOTOK_RETURN(status=1)

        ! clear:
        deallocate( at, bt )

      case default

        write (gol,'("unsupported rank : ",i1)') mf%rnk; call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine levi_Init_mf


  ! ***

  subroutine GetTimeRecordNumber( mf, t1, t2, irec, status )

    use GO, only : TDate, Get
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Inquire_Variable, MDF_Get_Var

    ! --- in/out -------------------------------

    type(TMeteoFile_tm5_nc), intent(inout)  ::  mf
    type(TDate), intent(in)                 ::  t1, t2
    integer, intent(out)                    ::  irec
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GetTimeRecordNumber'

    ! --- local -------------------------------

    integer               ::  timevalues1(6), timevalues2(6)
    integer               ::  dimid, varid
    integer               ::  ntime, itime
    integer, allocatable  ::  timevalues_bounds(:,:,:)

    ! --- begin ---------------------------------

    ! target times:
    call Get( t1, time6=timevalues1 )
    call Get( t2, time6=timevalues2 )

    ! number of time records:
    call MDF_Inq_DimID( mf%hid, 'time', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( mf%hid, dimid, status, length=ntime )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    allocate( timevalues_bounds(2,6,ntime) )

    ! read time boundaries:
    call MDF_Inq_VarID( mf%hid, 'timevalues_bounds', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( mf%hid, varid, timevalues_bounds, status )
    IF_NOTOK_RETURN(status=1)

    ! search ...
    irec = -1
    do itime = 1, ntime
      ! compare:
      if ( all(timevalues_bounds(1,:,itime) == timevalues1) .and. &
           all(timevalues_bounds(2,:,itime) == timevalues2) ) then
        irec = itime
        exit
      end if
    end do

    ! check ...
    if ( irec < 0 ) then
      write (gol,'("could not find time record for:")'); call goErr
      write (gol,'("  t1, t2  : ",i4,2("-",i2.2)," ",i2.2,2(":",i2.2)," , ",'// &
                               & 'i4,2("-",i2.2)," ",i2.2,2(":",i2.2))') timevalues1, timevalues2; call goErr
      write (gol,'("in time bounds:")'); call goErr
      do itime = 1, ntime
        write (gol,'("  ",i6,"  : ",i4,2("-",i2.2)," ",i2.2,2(":",i2.2)," , ",'// &
                                 & 'i4,2("-",i2.2)," ",i2.2,2(":",i2.2))') &
                       itime, timevalues_bounds(1,:,itime), timevalues_bounds(2,:,itime); call goErr
      end do
      write (gol,'("in file:")'); call goErr
      write (gol,'("  ",a)') trim(mf%fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! clear:
    deallocate( timevalues_bounds )

    ! ok
    status = 0

  end subroutine GetTimeRecordNumber


  ! ***


  subroutine mf_ReadRecord( mf, paramkey, unit, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                status )

    use PArray, only : pa_Done, pa_SetShape
    use GO, only : TDate, Get
    use GO, only : GO_Timer_Start, GO_Timer_End
    use GO, only : goReadFromLine
    use Grid, only : TllGridInfo, TLevelInfo
    use Grid, only : Init, Done
    use MDF, only : MDF_Open, MDF_Close
    use MDF, only : MDF_READ
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Inquire_Variable, MDF_Get_Var
    use MDF, only : MDF_Get_Att

    use TMM_CF, only : TMM_CF_Convert_Units

    ! --- in/out -------------------------------

    type(TMeteoFile_tm5_nc), intent(inout)  ::  mf
    character(len=*), intent(in)         ::  paramkey
    character(len=*), intent(in)         ::  unit
    type(TDate), intent(in)              ::  t1, t2
    character(len=1), intent(in)         ::  nuv
    character(len=2), intent(out)        ::  gridtype
    type(TLevelInfo), intent(out)        ::  levi
    character(len=1), intent(in)         ::  nw
    type(TllGridInfo), intent(inout)     ::  lli
    real, pointer                        ::  ll(:,:,:)
    real, pointer                        ::  sp_ll(:,:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord'

    ! --- local -------------------------------

    integer               ::  varid
    integer               ::  irec
    integer               ::  shp(7), ndims
    character(len=64)     ::  cfunits
    real                  ::  ufac

    integer               ::  l, n
    character(len=1024)   ::  line
    character(len=32)     ::  vname
    integer               ::  iparam

    integer               ::  timevalues1(6), timevalues2(6)
    integer               ::  dimid_time, varid_time
    integer               ::  ntime, itime

    real, allocatable     :: dummy(:, :, :)

    ! --- begin ---------------------------------

    ! start timing:
    call GO_Timer_Start( itim_readrecord, status )
    IF_NOTOK_RETURN(status=1)

    ! input ?
    if ( mf%io /= 'i' ) then
      write (gol,'("file should have been opened for input, but io=",a)') mf%io; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! *** chache
    
    ! input file not cached?
    if ( trim(mf%fname) /= trim(mf%cache_fname) ) then
    
      ! ! info ...
      ! write (gol,'("fill cache with ",a)') trim(mf%fname); call goPr
      
      ! store name:
      mf%cache_fname = trim(mf%fname)

      ! open for reading:
      call MDF_Open( trim(mf%fname), output_type, MDF_READ, mf%hid, status )
      IF_NOTOK_RETURN(status=1)

      ! setup grid definition:
      call lli_Init_mf( mf%cache_lli, nuv, mf, status )
      IF_NOTOK_RETURN(status=1)

      ! setup level definition:
      call levi_Init_mf( mf%cache_levi, mf, status )
      IF_NOTOK_RETURN(status=1)

      ! copy the paramkeys ('-aa-bb-cc-') except for the first '-' :
      l = len_trim(mf%paramkeys)
      line = mf%paramkeys(2:l)
      ! loop over all parameters:
      do iparam = 1, mf%nparam
    
        ! split at '-', read first part:
        call goReadFromLine( line, vname, status, sep='-' )
        IF_NOTOK_RETURN(status=1)
        
        !! info ...
        !write (gol,'("  param ",a," ...")') trim(vname); call goPr

        ! store name for searching:
        mf%cache_field(iparam)%paramkey = trim(vname)

        ! search variable:
        call MDF_Inq_VarID( mf%hid, trim(vname), varid, status )
        IF_NOTOK_RETURN(status=1)

        ! number of dimensions:
        call MDF_Inquire_Variable( mf%hid, varid, status, ndims=ndims )
        IF_NOTOK_RETURN(status=1)
        ! init shape:
        shp = 1
        ! reset first dimensions to actual shape:
        call MDF_Inquire_Variable( mf%hid, varid, status, shp=shp(1:ndims) )
        IF_NOTOK_RETURN(status=1)
        
        ! set shape:
        !if ( associated(mf%cache_field(iparam)%field) ) then
        !  print *, 'xxx1 ', shape(mf%cache_field(iparam)%field)
        !else
        !  print *, 'xxx1 cached field not allocated yet'
        !end if
        !print *, '  x1 ', shp(1:4)
        call pa_SetShape( mf%cache_field(iparam)%field, shp(1:4) )

        ! constant field ?
        if ( trim(mf%tres) == 'constant' ) then

          ! extract data array:
          select case ( ndims )
            !* 2D fields
            case ( 2 )
              ! 2D, single record
              call MDF_Get_Var( mf%hid, varid, mf%cache_field(iparam)%field(:,:,1,1), status )
              IF_NOTOK_RETURN(status=1)
            !* 3d fields
            case ( 3 )
              ! 3D, single record
              call MDF_Get_Var( mf%hid, varid, mf%cache_field(iparam)%field(:,:,:,1), status )
              IF_NOTOK_RETURN(status=1)
            !* unknown
            case default
              write (gol,'("unsupported data rank:",i6)') ndims; call goErr
              TRACEBACK; status=1; return
          end select

        else

          ! number of time records:
          call MDF_Inq_DimID( mf%hid, 'time', dimid_time, status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Inquire_Dimension( mf%hid, dimid_time, status, length=ntime )
          IF_NOTOK_RETURN(status=1)
          ! storage:
          call pa_SetShape( mf%cache_timevalues_bounds, (/2,6,ntime/) )
          ! read time boundaries:
          call MDF_Inq_VarID( mf%hid, 'timevalues_bounds', varid_time, status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Get_Var( mf%hid, varid_time, mf%cache_timevalues_bounds, status )
          IF_NOTOK_RETURN(status=1)

          ! extract data array:
          select case ( ndims )
            !* 2D temporal fields
            case ( 2+1 )
              ! complete record:
              call MDF_Get_Var( mf%hid, varid, mf%cache_field(iparam)%field(:,:,:,1), status )
              IF_NOTOK_RETURN(status=1)
            !* 3d temporal fields
            case ( 3+1 )
              ! complete record:
              call MDF_Get_Var( mf%hid, varid, mf%cache_field(iparam)%field(:,:,:,:), status )
              IF_NOTOK_RETURN(status=1)
            !* unknown
            case default
              write (gol,'("unsupported data rank:",i6)') ndims; call goErr
              TRACEBACK; status=1; return
          end select

        end if

        ! get unit from field in file:
        call MDF_Get_Att( mf%hid, varid, 'units', cfunits, status )
        IF_NOTOK_RETURN(status=1)
        ! conversion factor:
        call TMM_CF_Convert_Units( trim(cfunits), trim(unit), ufac, status )
        IF_NOTOK_RETURN(status=1)
        ! apply ?
        if ( ufac /= 1.0 ) then
          ! convert:
          mf%cache_field(iparam)%field = mf%cache_field(iparam)%field * ufac
          !! info ...
          !write (gol,'("      convert `",a,"` from `",a,"` to `",a,"` with factor ",f8.2," ; new range ",2f8.2)') &
          !               trim(vname), trim(cfunits), trim(unit), ufac, minval(ll), maxval(ll); call goPr
        end if
        
        ! ~ 

        ! surface pressure required ?
        if ( mf%rnk == 3 ) then
          ! search variable:
          select case ( trim(vname) )
            case ( 'mfu' )
              call MDF_Inq_VarID( mf%hid, 'ps_u', varid, status )
              IF_NOTOK_RETURN(status=1)
            case ( 'mfv' )
              call MDF_Inq_VarID( mf%hid, 'ps_v', varid, status )
              IF_NOTOK_RETURN(status=1)
            case default
              call MDF_Inq_VarID( mf%hid, 'ps'  , varid, status )
              IF_NOTOK_RETURN(status=1)
          end select
          ! number of dimensions:
          call MDF_Inquire_Variable( mf%hid, varid, status, ndims=ndims )
          IF_NOTOK_RETURN(status=1)
          ! init shape:
          shp = 1
          ! reset first dimensions to actual shape:
          call MDF_Inquire_Variable( mf%hid, varid, status, shp=shp(1:ndims) )
          IF_NOTOK_RETURN(status=1)
          ! storage:
          call pa_SetShape( mf%cache_field(iparam)%sp, shp(1:3) )
          ! complete record:
          call MDF_Get_Var( mf%hid, varid, mf%cache_field(iparam)%sp, status )
          IF_NOTOK_RETURN(status=1)
        else
          ! for safety ...
          call pa_Done( mf%cache_field(iparam)%sp )
        end if

      end do  ! params

      ! close
      call MDF_Close( mf%hid, status )
      IF_NOTOK_RETURN(status=1)
      
    end if  ! fill cache


    ! *** variable id

    
    ! init varid (index of chached field):
    varid = -999
    ! loop:
    do iparam = 1, mf%nparam
      ! match?
      if ( trim(paramkey) == trim(mf%cache_field(iparam)%paramkey) ) then
        varid = iparam
        exit
      end if
    end do  ! params
    ! check ..
    if ( varid < 0 ) then
      write (gol,'("param `",a,"` not chached")') trim(paramkey); call goErr
      TRACEBACK; status=1; return
    end if
    !! info ...
    !write (gol,'("selected cached param ",i0," `",a,"` ...")') varid, trim(paramkey); call goPr

    ! *** grid definition

    ! always regular lat/lon grid ..
    gridtype = 'll'

    
    ! setup grid definition as copy:
    call Init( lli, mf%cache_lli%lon_deg(1), mf%cache_lli%dlon_deg, mf%cache_lli%im, &
                    mf%cache_lli%lat_deg(1), mf%cache_lli%dlat_deg, mf%cache_lli%jm, status  )
    IF_NOTOK_RETURN(status=1)


    ! setup level definition as copy:
    call Init( levi, mf%cache_levi, status )
    IF_NOTOK_RETURN(status=1)

    
    ! get dimensions:
    shp = 1
    shp(1:4) = shape(mf%cache_field(varid)%field)

    ! *** data

    ! constant field ?
    if ( trim(mf%tres) == 'constant' ) then

      ! extract data array:
      select case ( mf%rnk )
        !* 2D fields
        case ( 2 )
          ! storage:
          shp(3) = 1
          call pa_SetShape( ll, shp(1:3) )
          !! complete record:
          !call MDF_Get_Var( mf%hid, varid, ll(:,:,1), status )
          !IF_NOTOK_RETURN(status=1)
          ! single record:
          ll(:,:,1) = mf%cache_field(varid)%field(:,:,1,1)
          !! info ...
          !write (gol,'("  copied `",a,"` 2D constant field from cache ...")') trim(paramkey); call goPr
        !* 3d fields
        case ( 3 )
          ! storage:
          call pa_SetShape( ll, shp(1:3) )
          !! complete record:
          !call MDF_Get_Var( mf%hid, varid, ll, status )
          !IF_NOTOK_RETURN(status=1)
          ! single record:
          ll = mf%cache_field(varid)%field(:,:,:,1)
          !! info ...
          !write (gol,'("  copied `",a,"` 3D constant field from cache ...")') trim(paramkey); call goPr
        !* unknown
        case default
          write (gol,'("unsupported data rank:",i6)') mf%rnk; call goErr
          TRACEBACK; status=1; return
      end select

      ! for safety ...
      call pa_Done( sp_ll )

    else

      !! which time record ?
      !call GetTimeRecordNumber( mf, t1, t2, irec, status )
      !IF_NOTOK_RETURN(status=1)

      ! target times:
      call Get( t1, time6=timevalues1 )
      call Get( t2, time6=timevalues2 )
      ! count:
      ntime = size(mf%cache_timevalues_bounds,3)
      ! search ...
      irec = -1
      do itime = 1, ntime
        ! compare:
        if ( all(mf%cache_timevalues_bounds(1,:,itime) == timevalues1) .and. &
             all(mf%cache_timevalues_bounds(2,:,itime) == timevalues2) ) then
          irec = itime
          exit
        end if
      end do
      ! check ...
      if ( irec < 0 ) then
        write (gol,'("could not find time record for:")'); call goErr
        write (gol,'("  t1, t2  : ",i4,2("-",i2.2)," ",i2.2,2(":",i2.2)," , ",'// &
                                 & 'i4,2("-",i2.2)," ",i2.2,2(":",i2.2))') timevalues1, timevalues2; call goErr
        write (gol,'("in time bounds:")'); call goErr
        do itime = 1, ntime
          write (gol,'("  ",i6,"  : ",i4,2("-",i2.2)," ",i2.2,2(":",i2.2)," , ",'// &
                                   & 'i4,2("-",i2.2)," ",i2.2,2(":",i2.2))') &
                         itime, mf%cache_timevalues_bounds(1,:,itime), mf%cache_timevalues_bounds(2,:,itime); call goErr
        end do
        write (gol,'("in file:")'); call goErr
        write (gol,'("  ",a)') trim(mf%fname); call goErr
        TRACEBACK; status=1; return
      end if

      ! extract data array:
      select case ( mf%rnk )
        !* 2D temporal fields
        case ( 2 )
          ! storage:
          shp(3) = 1
          call pa_SetShape( ll, shp(1:3) )
          !! complete record:
          !call MDF_Get_Var( mf%hid, varid, ll(:,:,1), status, &
          !                      start=(/1,1,irec/), count=(/shp(1),shp(2),1/) )
          !IF_NOTOK_RETURN(status=1)
          ! copy record:
          ll(:,:,1) = mf%cache_field(varid)%field(:,:,irec,1)
          !! info ...
          !write (gol,'("  copied `",a,"` 2D record from cache record ",i0," ...")') trim(paramkey), irec; call goPr
          ! for safety ...
          call pa_Done( sp_ll )
        !* 3d temporal fields
        case ( 3 )
          ! storage:
          call pa_SetShape( ll, shp(1:3) )
          !! complete record:
          !call MDF_Get_Var( mf%hid, varid, ll, status, &
          !                      start=(/1,1,1,irec/), count=(/shp(1),shp(2),shp(3),1/) )
          !IF_NOTOK_RETURN(status=1)
          ! copy record:
          allocate(dummy(shp(1), shp(2), shp(3)))
          dummy = mf%cache_field(varid)%field(:,:,:,irec)
          ll = dummy ! mf%cache_field(varid)%field(:,:,:,irec)
          deallocate(dummy)
          !! info ...
          !write (gol,'("  copied `",a,"` 3D record from cache record ",i0," ...")') trim(paramkey), irec; call goPr
          ! storage:
          call pa_SetShape( sp_ll, shp(1:2) )
          !! complete record:
          !call MDF_Get_Var( mf%hid, varid, sp_ll, status, &
          !                      start=(/1,1,irec/), count=(/shp(1),shp(2),1/) )
          !IF_NOTOK_RETURN(status=1)
          ! copy record:
          sp_ll = mf%cache_field(varid)%sp(:,:,irec)
          !! info ...
          !write (gol,'("  copied surface pressure from cache record ",i0," ...")') irec; call goPr
        !* unknown
        case default
          write (gol,'("unsupported data rank:",i6)') mf%rnk; call goErr
          TRACEBACK; status=1; return
      end select

    end if  ! 2D or 3D

    ! *** unit conversion

    !! get unit from field in file:
    !call MDF_Get_Att( mf%hid, varid, 'units', cfunits, status )
    !IF_NOTOK_RETURN(status=1)
    !! conversion factor:
    !call TMM_CF_Convert_Units( trim(cfunits), trim(unit), ufac, status )
    !IF_NOTOK_RETURN(status=1)
    !! apply ?
    !if ( ufac /= 1.0 ) then
    !  ! convert:
    !  ll = ll * ufac
    !  !! info ...
    !  !write (gol,'("      convert `",a,"` from `",a,"` to `",a,"` with factor ",f8.2," ; new range ",2f8.2)') &
    !  !               trim(paramkey), trim(cfunits), trim(unit), ufac, minval(ll), maxval(ll); call goPr
    !end if
    !
    !! *** surface pressure
    !
    !! surface pressure required ?
    !if ( mf%rnk == 3 ) then
    !  ! search variable:
    !  select case ( nuv )
    !    case ( 'n' )
    !      call MDF_Inq_VarID( mf%hid, 'ps', varid, status )
    !      IF_NOTOK_RETURN(status=1)
    !    case ( 'u' )
    !      call MDF_Inq_VarID( mf%hid, 'ps_u', varid, status )
    !      IF_NOTOK_RETURN(status=1)
    !    case ( 'v' )
    !      call MDF_Inq_VarID( mf%hid, 'ps_v', varid, status )
    !      IF_NOTOK_RETURN(status=1)
    !    case default
    !      write (gol,'("unsupported nuv :",a)') nuv; call goErr
    !      TRACEBACK; status=1; return
    !  end select
    !  ! storage:
    !  call pa_SetShape( sp_ll, shp(1:2) )
    !  ! complete record:
    !  call MDF_Get_Var( mf%hid, varid, sp_ll, status, &
    !                        start=(/1,1,irec/), count=(/shp(1),shp(2),1/) )
    !  IF_NOTOK_RETURN(status=1)
    !else
    !  ! for safety ...
    !  call pa_Done( sp_ll )
    !end if

    ! ***

    !! close
    !call MDF_Close( mf%hid, status )
    !IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End( itim_readrecord, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine mf_ReadRecord


  ! ******************************************************************
  ! ***
  ! *** output
  ! ***
  ! ******************************************************************


  subroutine CF_Put_Standard_Atts( hid, varid, cf_standard_name, cf_units, status )

    use MDF, only : MDF_Put_Att

    ! --- in/out ---------------------------------

    integer, intent(in)                           ::  hid
    integer, intent(in)                           ::  varid
    character(len=*), intent(in)                  ::  cf_standard_name
    character(len=*), intent(in)                  ::  cf_units
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/CF_Put_Standard_Atts'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! add standard name attribute:
    call MDF_Put_Att( hid, varid, 'standard_name', trim(cf_standard_name), status )
    IF_NOTOK_RETURN(status=1)

    ! add units attribute:
    call MDF_Put_Att( hid, varid, 'units', trim(cf_units), status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine CF_Put_Standard_Atts


  ! ***


  subroutine Define_File( mf, lli, nuv, status, levi, nw )!, nlev )

    use GO    , only : goReadFromLine
    use Binas , only : grav, ae
    use Grid  , only : TllGridInfo, TLevelInfo, AreaOper
    use MDF   , only : MDF_Put_Att, MDF_GLOBAL
    use MDF   , only : MDF_Def_Var, MDF_Put_Var, MDF_FLOAT, MDF_DOUBLE, MDF_INT
    use MDF   , only : MDF_Def_Dim, MDF_UNLIMITED
    use MDF   , only : MDF_EndDef
    use TMM_CF, only : TMM_CF_Standard_Units, TMM_CF_Convert_Name

    ! --- in/out -------------------------------

    type(TMeteoFile_tm5_nc), intent(inout)   ::  mf
    type(TllGridInfo), intent(in)            ::  lli
    character(len=1), intent(in)             ::  nuv
    integer, intent(out)                     ::  status

    type(TLevelInfo), intent(in), optional   ::  levi
    character(len=1), intent(in), optional   ::  nw
    !integer, intent(in), optional            ::  nlev

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Define_File'

    ! --- local ----------------------------------

    character(len=DUMMY_STR_LEN)    ::  units
    character(len=DUMMY_STR_LEN)    ::  cell_measure
    character(len=DUMMY_STR_LEN)    ::  cell_methods
    character(len=DUMMY_STR_LEN)    ::  coordinates
    real, allocatable               ::  pat(:,:)

    integer                         ::  l, n
    character(len=LONG_STR_LEN)     ::  line
    character(len=32)               ::  paramkey
    integer                         ::  iparam
    character(len=LONG_STR_LEN)     ::  cf_standard_name
    character(len=DUMMY_STR_LEN)    ::  cf_units
    integer                         ::  varid

    ! --- begin ----------------------------------

    !
    ! CF standard global attribues
    !

    ! convention id:
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'Conventions', 'CF-1.4', status )
    IF_NOTOK_RETURN(status=1)

    ! descriptive title:
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'title', 'TM meteo file', status )
    IF_NOTOK_RETURN(status=1)

    ! who?
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'institution', 'TM community', status )
    IF_NOTOK_RETURN(status=1)

    ! from ?
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'source', 'TM produced meteo file', status )
    IF_NOTOK_RETURN(status=1)

    ! how ?
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'history', 'None', status )
    IF_NOTOK_RETURN(status=1)

    ! published material:
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'references', 'None', status )
    IF_NOTOK_RETURN(status=1)

    ! other:
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'comment', 'None', status )
    IF_NOTOK_RETURN(status=1)

    !
    ! TMM specific global attributes
    !

    ! file format
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'tmm_format', trim(output_format), status )
    IF_NOTOK_RETURN(status=1)
    ! grid type
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'tmm_gridtype', 'll', status )
    IF_NOTOK_RETURN(status=1)
    ! gravity constant:
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'grav', grav, status )
    IF_NOTOK_RETURN(status=1)
    ! earth radius:
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'ae', ae, status )
    IF_NOTOK_RETURN(status=1)

    ! save first and last lon/lat (center) for use with HIPHOP
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'lonmin', lli%lon_deg(1)       , status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'lonmax', lli%lon_deg(lli%nlon), status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'latmin', lli%lat_deg(1)       , status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'latmax', lli%lat_deg(lli%nlat), status )
    IF_NOTOK_RETURN(status=1)

    !
    ! axis
    !

    ! extra coordinate axes (with other names than the dimensions):
    coordinates = ''

    ! auxilary cell info:
    cell_measure = ''

    ! how data was formed:
    cell_methods = ''

    ! vertices:
    call MDF_Def_Dim( mf%hid, 'nv', 2, mf%dimid_nv, status )
    IF_NOTOK_RETURN(status=1)

    ! * longitudes

    ! longitude dimension:
    call MDF_Def_Dim( mf%hid, 'lon', lli%nlon, mf%dimid_lon, status )
    IF_NOTOK_RETURN(status=1)

    ! variable:
    call MDF_Def_Var( mf%hid, 'lon', MDF_DOUBLE, (/mf%dimid_lon/), mf%varid_lon, status )
    IF_NOTOK_RETURN(status=1)

    ! convert to CF name and unit:
    cf_standard_name = 'longitude'
    call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
    IF_NOTOK_RETURN(status=1)
    ! add standard_name and units attributes:
    call CF_Put_Standard_Atts( mf%hid, mf%varid_lon, trim(cf_standard_name), trim(cf_units), status )
    IF_NOTOK_RETURN(status=1)

    ! generic axis name:
    call MDF_Put_Att( mf%hid, mf%varid_lon, 'axis', 'X', status )
    IF_NOTOK_RETURN(status=1)

    ! add attribute with name of variable with boundaries:
    call MDF_Put_Att( mf%hid, mf%varid_lon, 'bounds', 'lon_bounds', status )
    IF_NOTOK_RETURN(status=1)
    ! add variable for boundaries:
    call MDF_Def_Var( mf%hid, 'lon_bounds', MDF_DOUBLE, (/mf%dimid_nv,mf%dimid_lon/), mf%varid_lon_bounds, status )
    IF_NOTOK_RETURN(status=1)

    ! extra:
    if ( nuv /= 'n' ) then
      ! dimension:
      call MDF_Def_Dim( mf%hid, 'lonb', lli%nlon+1, mf%dimid_lonb, status )
      IF_NOTOK_RETURN(status=1)
      ! variable:
      call MDF_Def_Var( mf%hid, 'lonb', MDF_DOUBLE, (/mf%dimid_lonb/), mf%varid_lonb, status )
      IF_NOTOK_RETURN(status=1)
      ! add standard_name and units attributes:
      call CF_Put_Standard_Atts( mf%hid, mf%varid_lonb, trim(cf_standard_name), trim(cf_units), status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! * latitudes

    ! latitude dimension:
    call MDF_Def_Dim( mf%hid, 'lat', lli%nlat, mf%dimid_lat, status )
    IF_NOTOK_RETURN(status=1)

    ! variable:
    call MDF_Def_Var( mf%hid, 'lat', MDF_DOUBLE, (/mf%dimid_lat/), mf%varid_lat, status )
    IF_NOTOK_RETURN(status=1)

    ! convert to CF name and unit:
    cf_standard_name = 'latitude'
    call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
    IF_NOTOK_RETURN(status=1)
    ! add standard_name and units attributes:
    call CF_Put_Standard_Atts( mf%hid, mf%varid_lat, trim(cf_standard_name), trim(cf_units), status )
    IF_NOTOK_RETURN(status=1)

    ! generic axis name:
    call MDF_Put_Att( mf%hid, mf%varid_lat, 'axis', 'Y', status )
    IF_NOTOK_RETURN(status=1)

    ! add attribute with name of variable with boundaries:
    call MDF_Put_Att( mf%hid, mf%varid_lat, 'bounds', 'lat_bounds', status )
    IF_NOTOK_RETURN(status=1)
    ! add variable for boundaries:
    call MDF_Def_Var( mf%hid, 'lat_bounds', MDF_DOUBLE, (/mf%dimid_nv,mf%dimid_lat/), mf%varid_lat_bounds, status )
    IF_NOTOK_RETURN(status=1)

    ! extra:
    if ( nuv /= 'n' ) then
      ! dimension:
      call MDF_Def_Dim( mf%hid, 'latb', lli%nlat+1, mf%dimid_latb, status )
      IF_NOTOK_RETURN(status=1)
      ! variable:
      call MDF_Def_Var( mf%hid, 'latb', MDF_DOUBLE, (/mf%dimid_latb/), mf%varid_latb, status )
      IF_NOTOK_RETURN(status=1)
      ! add standard_name and units attributes:
      call CF_Put_Standard_Atts( mf%hid, mf%varid_latb, trim(cf_standard_name), trim(cf_units), status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! * area

    ! cell variable ?
    if ( nuv == 'n' ) then
      ! also provide the area:
      cell_measure = trim(cell_measure)//' area: cell_area'
      ! cell area:
      call MDF_Def_Var( mf%hid, 'cell_area', MDF_FLOAT, (/mf%dimid_lon,mf%dimid_lat/), mf%varid_cell_area, status )
      IF_NOTOK_RETURN(status=1)
      ! convert to CF name and unit:
      cf_standard_name = 'cell_area'
      call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
      IF_NOTOK_RETURN(status=1)
      ! add standard_name and units attributes:
      call CF_Put_Standard_Atts( mf%hid, mf%varid_cell_area, trim(cf_standard_name), trim(cf_units), status )
      IF_NOTOK_RETURN(status=1)
      ! add description:
      call MDF_Put_Att( mf%hid, mf%varid_cell_area, 'long_name', 'area of grid cell', status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! * time

    if ( trim(mf%tres) /= 'constant' ) then

      ! time dimension:
      call MDF_Def_Dim( mf%hid, 'time', MDF_UNLIMITED, mf%dimid_time, status )
      IF_NOTOK_RETURN(status=1)

      ! standard units is single value only ...
      write (units,'("seconds since ",i4.4,2("-",i2.2)," ",i2.2,2(":",i2.2))') since_time6

      ! variable:
      call MDF_Def_Var( mf%hid, 'time', MDF_DOUBLE, (/mf%dimid_time/), mf%varid_time, status )
      IF_NOTOK_RETURN(status=1)
      ! add units attribute:
      call MDF_Put_Att( mf%hid, mf%varid_time, 'standard_name', 'time', status )
      IF_NOTOK_RETURN(status=1)
      ! add units attribute:
      call MDF_Put_Att( mf%hid, mf%varid_time, 'units', trim(units), status )
      IF_NOTOK_RETURN(status=1)

      ! name of variable with interval bounds:
      call MDF_Put_Att( mf%hid, mf%varid_time, 'bounds', 'time_bounds', status )
      IF_NOTOK_RETURN(status=1)

      ! variable:
      call MDF_Def_Var( mf%hid, 'time_bounds', MDF_DOUBLE, (/mf%dimid_nv,mf%dimid_time/), mf%varid_time_bounds, status )
      IF_NOTOK_RETURN(status=1)
      ! add units attribute:
      call MDF_Put_Att( mf%hid, mf%varid_time, 'units', trim(units), status )
      IF_NOTOK_RETURN(status=1)

      ! time averages ?
      if ( mf%is_aver ) then
        ! value is the mean over a time interval:
        cell_methods = trim(cell_methods)//' time: mean'
      else
        ! instant fields:
        cell_methods = trim(cell_methods)//' time: point'
      end if

      ! variable:
      call MDF_Def_Var( mf%hid, 'reftime', MDF_DOUBLE, (/mf%dimid_time/), mf%varid_reftime, status )
      IF_NOTOK_RETURN(status=1)
      ! add units attribute:
      call MDF_Put_Att( mf%hid, mf%varid_reftime, 'units', trim(units), status )
      IF_NOTOK_RETURN(status=1)
      ! description:
      call MDF_Put_Att( mf%hid, mf%varid_reftime, 'long_name', 'reference time', status )
      IF_NOTOK_RETURN(status=1)
      ! extra time coordinates:
      coordinates = trim(coordinates)//' reftime'

    end if

    ! * time values

    if ( trim(mf%tres) /= 'constant' ) then

      ! time values dimension:
      call MDF_Def_Dim( mf%hid, 'timeval', 6, mf%dimid_timeval, status )
      IF_NOTOK_RETURN(status=1)

      ! variable:
      call MDF_Def_Var( mf%hid, 'timevalues', MDF_INT, (/mf%dimid_timeval,mf%dimid_time/), mf%varid_timevalues, status )
      IF_NOTOK_RETURN(status=1)
      ! description:
      call MDF_Put_Att( mf%hid, mf%varid_timevalues, 'long_name', 'year month day hour minute second', status )
      IF_NOTOK_RETURN(status=1)
      ! add units attribute:
      call MDF_Put_Att( mf%hid, mf%varid_timevalues, 'units', '1', status )
      IF_NOTOK_RETURN(status=1)

      ! name of variable with interval bounds:
      call MDF_Put_Att( mf%hid, mf%varid_timevalues, 'bounds', 'timevalues_bounds', status )
      IF_NOTOK_RETURN(status=1)

      ! variable:
      call MDF_Def_Var( mf%hid, 'timevalues_bounds', MDF_INT, &
                           (/mf%dimid_nv,mf%dimid_timeval,mf%dimid_time/), &
                           mf%varid_timevalues_bounds, status )
      IF_NOTOK_RETURN(status=1)
      ! description:
      call MDF_Put_Att( mf%hid, mf%varid_timevalues_bounds, 'long_name', 'year month day hour minute second', status )
      IF_NOTOK_RETURN(status=1)
      ! add units attribute:
      call MDF_Put_Att( mf%hid, mf%varid_timevalues_bounds, 'units', '1', status )
      IF_NOTOK_RETURN(status=1)

      ! variable:
      call MDF_Def_Var( mf%hid, 'reftimevalues', MDF_INT, (/mf%dimid_timeval,mf%dimid_time/), mf%varid_reftimevalues, status )
      IF_NOTOK_RETURN(status=1)
      ! description:
      call MDF_Put_Att( mf%hid, mf%varid_reftimevalues, 'long_name', 'reference time values', status )
      IF_NOTOK_RETURN(status=1)
      ! add units attribute:
      call MDF_Put_Att( mf%hid, mf%varid_reftimevalues, 'units', '1', status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! * levels

    ! levels ?
    if ( present(levi) ) then

      ! check ...
      if ( .not. present(nw) ) then
        write (gol,'("optional argument levi without nw ...")'); call goErr
        TRACEBACK; status=1; return
      end if

      !! tm5 specials ...
      !call MDF_Def_Dim( mf%hid, 'tm5_lm', levi%nlev, mf%dimid_tm5_lm, status )
      !IF_NOTOK_RETURN(status=1)
      !call MDF_Def_Dim( mf%hid, 'tm5_lmb', levi%nlev+1, mf%dimid_tm5_lmb, status )
      !IF_NOTOK_RETURN(status=1)
      !call MDF_Def_Var( mf%hid, 'tm5_at', MDF_DOUBLE, (/mf%dimid_tm5_lmb/), mf%varid_tm5_at, status )
      !IF_NOTOK_RETURN(status=1)
      !call MDF_Def_Var( mf%hid, 'tm5_bt', MDF_DOUBLE, (/mf%dimid_tm5_lmb/), mf%varid_tm5_bt, status )
      !IF_NOTOK_RETURN(status=1)

      ! save nw to facilitate reading:
      call MDF_Put_Att( mf%hid, MDF_GLOBAL, 'nw', trim(nw), status )
      IF_NOTOK_RETURN(status=1)

      ! where defined ?
      select case ( nw )

        ! layer mid, or selection
        case ( 'n' )!, '*' )

          ! hybride layers:
          !if ( present(nlev) ) then
          !  ! lowest layers only:
          !  call MDF_Def_Dim( mf%hid, 'lev', nlev, mf%dimid_lev, status )
          !  IF_NOTOK_RETURN(status=1)
          !else
            ! full dimension:
            call MDF_Def_Dim( mf%hid, 'lev', levi%nlev, mf%dimid_lev, status )
            IF_NOTOK_RETURN(status=1)
          !end if

          ! variable:
          call MDF_Def_Var( mf%hid, 'lev', MDF_DOUBLE, (/mf%dimid_lev/), mf%varid_lev, status )
          IF_NOTOK_RETURN(status=1)
          ! convert to CF name and unit:
          cf_standard_name = 'atmosphere_hybrid_sigma_pressure_coordinate'
          call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
          IF_NOTOK_RETURN(status=1)
          ! add standard_name and units attributes:
          call CF_Put_Standard_Atts( mf%hid, mf%varid_lev, trim(cf_standard_name), trim(cf_units), status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'long_name', 'pressure at layer midpoints', status )
          IF_NOTOK_RETURN(status=1)
          ! direction of increasing pressure:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'positive', 'down', status )
          IF_NOTOK_RETURN(status=1)
          ! how to compute pressure:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'forumula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)', status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'forumula_terms', 'ap: ap b: b ps: ps', status )
          IF_NOTOK_RETURN(status=1)
          ! generic axis name:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'axis', 'Z', status )
          IF_NOTOK_RETURN(status=1)

          ! extra
          if ( nuv /= 'n' ) then
            ! hybride layers:
            call MDF_Def_Dim( mf%hid, 'lev_u', levi%nlev, mf%dimid_lev_u, status )
            IF_NOTOK_RETURN(status=1)
            ! variable:
            call MDF_Def_Var( mf%hid, 'lev_u', MDF_DOUBLE, (/mf%dimid_lev_u/), mf%varid_lev_u, status )
            IF_NOTOK_RETURN(status=1)
            ! convert to CF name and unit:
            cf_standard_name = 'atmosphere_hybrid_sigma_pressure_coordinate'
            call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
            IF_NOTOK_RETURN(status=1)
            ! add standard_name and units attributes:
            call CF_Put_Standard_Atts( mf%hid, mf%varid_lev_u, trim(cf_standard_name), trim(cf_units), status )
            IF_NOTOK_RETURN(status=1)
            ! description:
            call MDF_Put_Att( mf%hid, mf%varid_lev_u, 'long_name', 'pressure at layer midpoints', status )
            IF_NOTOK_RETURN(status=1)
            ! direction of increasing pressure:
            call MDF_Put_Att( mf%hid, mf%varid_lev_u, 'positive', 'down', status )
            IF_NOTOK_RETURN(status=1)
            ! how to compute pressure:
            call MDF_Put_Att( mf%hid, mf%varid_lev_u, 'forumula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)', status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Put_Att( mf%hid, mf%varid_lev_u, 'forumula_terms', 'ap: ap b: b ps: ps_u', status )
            IF_NOTOK_RETURN(status=1)

            ! hybride layers:
            call MDF_Def_Dim( mf%hid, 'lev_v', levi%nlev, mf%dimid_lev_v, status )
            IF_NOTOK_RETURN(status=1)
            ! variable:
            call MDF_Def_Var( mf%hid, 'lev_v', MDF_DOUBLE, (/mf%dimid_lev_v/), mf%varid_lev_v, status )
            IF_NOTOK_RETURN(status=1)
            ! convert to CF name and unit:
            cf_standard_name = 'atmosphere_hybrid_sigma_pressure_coordinate'
            call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
            IF_NOTOK_RETURN(status=1)
            ! add standard_name and units attributes:
            call CF_Put_Standard_Atts( mf%hid, mf%varid_lev_v, trim(cf_standard_name), trim(cf_units), status )
            IF_NOTOK_RETURN(status=1)
            ! description:
            call MDF_Put_Att( mf%hid, mf%varid_lev_v, 'long_name', 'pressure at layer midpoints', status )
            IF_NOTOK_RETURN(status=1)
            ! direction of increasing pressure:
            call MDF_Put_Att( mf%hid, mf%varid_lev_v, 'positive', 'down', status )
            IF_NOTOK_RETURN(status=1)
            ! how to compute pressure:
            call MDF_Put_Att( mf%hid, mf%varid_lev_v, 'forumula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)', status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Put_Att( mf%hid, mf%varid_lev_v, 'forumula_terms', 'ap: ap b: b ps: ps_v', status )
            IF_NOTOK_RETURN(status=1)
          end if

          ! variable:
          call MDF_Def_Var( mf%hid, 'ap', MDF_DOUBLE, (/mf%dimid_lev/), mf%varid_ap, status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_ap, 'long_name', 'atmosphere hybrid sigma pressure coefficient ap at layer midpoints', status )
          IF_NOTOK_RETURN(status=1)
          ! units:
          call MDF_Put_Att( mf%hid, mf%varid_ap, 'units', 'Pa', status )
          IF_NOTOK_RETURN(status=1)
          ! name of variable with boundary values:
          call MDF_Put_Att( mf%hid, mf%varid_ap, 'bounds', 'ap_bounds', status )
          IF_NOTOK_RETURN(status=1)

          ! variable:
          call MDF_Def_Var( mf%hid, 'b', MDF_DOUBLE, (/mf%dimid_lev/), mf%varid_b, status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_b, 'long_name', 'atmosphere hybrid sigma pressure coefficient b at layer midpoints', status )
          IF_NOTOK_RETURN(status=1)
          ! units:
          call MDF_Put_Att( mf%hid, mf%varid_b, 'units', '1', status )
          IF_NOTOK_RETURN(status=1)
          ! name of variable with boundary values:
          call MDF_Put_Att( mf%hid, mf%varid_b, 'bounds', 'b_bounds', status )
          IF_NOTOK_RETURN(status=1)

          ! variable:
          call MDF_Def_Var( mf%hid, 'ap_bounds', MDF_DOUBLE, (/mf%dimid_nv,mf%dimid_lev/), mf%varid_ap_bounds, status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_ap_bounds, 'long_name', 'atmosphere hybrid sigma pressure coefficient ap at layer interfaces', status )
          IF_NOTOK_RETURN(status=1)
          ! units:
          call MDF_Put_Att( mf%hid, mf%varid_ap_bounds, 'units', 'Pa', status )
          IF_NOTOK_RETURN(status=1)

          ! variable:
          call MDF_Def_Var( mf%hid, 'b_bounds', MDF_DOUBLE, (/mf%dimid_nv,mf%dimid_lev/), mf%varid_b_bounds, status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_b_bounds, 'long_name', 'atmosphere hybrid sigma pressure coefficient b at layer interfaces', status )
          IF_NOTOK_RETURN(status=1)
          ! units:
          call MDF_Put_Att( mf%hid, mf%varid_b_bounds, 'units', '1', status )
          IF_NOTOK_RETURN(status=1)

        ! layer interfaces
        case ( 'w' )

          ! hybride layers:
          call MDF_Def_Dim( mf%hid, 'lev', levi%nlev+1, mf%dimid_lev, status )
          IF_NOTOK_RETURN(status=1)

          ! variable:
          call MDF_Def_Var( mf%hid, 'lev', MDF_DOUBLE, (/mf%dimid_lev/), mf%varid_lev, status )
          IF_NOTOK_RETURN(status=1)
          ! convert to CF name and unit:
          cf_standard_name = 'atmosphere_hybrid_sigma_pressure_coordinate'
          call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
          IF_NOTOK_RETURN(status=1)
          ! add standard_name and units attributes:
          call CF_Put_Standard_Atts( mf%hid, mf%varid_lev, trim(cf_standard_name), trim(cf_units), status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'long_name', 'pressure at layer interfaces', status )
          IF_NOTOK_RETURN(status=1)
          ! direction of increasing pressure:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'positive', 'down', status )
          IF_NOTOK_RETURN(status=1)
          ! how to compute pressure:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'forumula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)', status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'forumula_terms', 'ap: ap b: b ps: ps', status )
          IF_NOTOK_RETURN(status=1)
          ! generic axis name:
          call MDF_Put_Att( mf%hid, mf%varid_lev, 'axis', 'Z', status )
          IF_NOTOK_RETURN(status=1)

          ! variable:
          call MDF_Def_Var( mf%hid, 'ap', MDF_DOUBLE, (/mf%dimid_lev/), mf%varid_ap, status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_ap, 'long_name', 'atmosphere hybrid sigma pressure coefficient ap at layer interfaces', status )
          IF_NOTOK_RETURN(status=1)
          ! units:
          call MDF_Put_Att( mf%hid, mf%varid_ap, 'units', 'Pa', status )
          IF_NOTOK_RETURN(status=1)

          ! variable:
          call MDF_Def_Var( mf%hid, 'b', MDF_DOUBLE, (/mf%dimid_lev/), mf%varid_b, status )
          IF_NOTOK_RETURN(status=1)
          ! description:
          call MDF_Put_Att( mf%hid, mf%varid_b, 'long_name', 'atmosphere hybrid sigma pressure coefficient b at layer interfaces', status )
          IF_NOTOK_RETURN(status=1)
          ! units:
          call MDF_Put_Att( mf%hid, mf%varid_b, 'units', '1', status )
          IF_NOTOK_RETURN(status=1)

        ! other ...
        case default

          write (gol,'("unsupported nw : ",a)') trim(nw); call goErr
          TRACEBACK; status=1; return

      end select

      ! surface pressure:
      call MDF_Def_Var( mf%hid, 'ps', MDF_FLOAT, (/mf%dimid_lon ,mf%dimid_lat ,mf%dimid_time/), mf%varid_ps, status )
      IF_NOTOK_RETURN(status=1)
      ! convert to CF name and unit:
      cf_standard_name = 'surface_air_pressure'
      call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
      IF_NOTOK_RETURN(status=1)
      ! add standard_name and units attributes:
      call CF_Put_Standard_Atts( mf%hid, mf%varid_ps, trim(cf_standard_name), trim(cf_units), status )
      IF_NOTOK_RETURN(status=1)
      ! for post processing ...
      l = len_trim(cell_measure)
      if ( l > 0 ) then
        call MDF_Put_Att( mf%hid, mf%varid_ps, 'cell_measure', cell_measure(2:l), status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! specify how computed:
      l = len_trim(cell_methods)
      if ( l > 0 ) then
        call MDF_Put_Att( mf%hid, mf%varid_ps, 'cell_methods', cell_methods(2:l), status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! extra coordinates that apply to this variable:
      l = len_trim(coordinates)
      if ( l > 0 ) then
        call MDF_Put_Att( mf%hid, mf%varid_ps, 'coordinates', coordinates(2:l), status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! extra:
      if ( nuv /= 'n' ) then
        ! surface pressure:
        call MDF_Def_Var( mf%hid, 'ps_u', MDF_FLOAT, (/mf%dimid_lonb,mf%dimid_lat ,mf%dimid_time/), mf%varid_ps_u, status )
        IF_NOTOK_RETURN(status=1)
        ! convert to CF name and unit:
        cf_standard_name = 'surface_air_pressure'
        call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
        IF_NOTOK_RETURN(status=1)
        ! add standard_name and units attributes:
        call CF_Put_Standard_Atts( mf%hid, mf%varid_ps_u, trim(cf_standard_name), trim(cf_units), status )
        IF_NOTOK_RETURN(status=1)

        ! surface pressure:
        call MDF_Def_Var( mf%hid, 'ps_v', MDF_FLOAT, (/mf%dimid_lon ,mf%dimid_latb,mf%dimid_time/), mf%varid_ps_v, status )
        IF_NOTOK_RETURN(status=1)
        ! convert to CF name and unit:
        cf_standard_name = 'surface_air_pressure'
        call TMM_CF_Standard_Units( trim(cf_standard_name), cf_units, status )
        IF_NOTOK_RETURN(status=1)
        ! add standard_name and units attributes:
        call CF_Put_Standard_Atts( mf%hid, mf%varid_ps_v, trim(cf_standard_name), trim(cf_units), status )
        IF_NOTOK_RETURN(status=1)
      end if

    end if  ! levels

    ! * meteo variables

    ! copy the paramkeys ('-aa-bb-cc-') except for the first '-' :
    l = len_trim(mf%paramkeys)
    line = mf%paramkeys(2:l)
    ! loop over all parameters:
    do iparam = 1, mf%nparam

      ! split at '-', read first part:
      call goReadFromLine( line, paramkey, status, sep='-' )
      IF_NOTOK_RETURN(status=1)

      ! define variable:
      if ( present(levi) ) then
        ! 3D field:
        if ( trim(paramkey) == 'mfu' ) then
          call MDF_Def_Var( mf%hid, trim(paramkey), MDF_FLOAT, (/mf%dimid_lonb,mf%dimid_lat ,mf%dimid_lev_u,mf%dimid_time/), varid, status )
          IF_NOTOK_RETURN(status=1)
        else if ( trim(paramkey) == 'mfv' ) then
          call MDF_Def_Var( mf%hid, trim(paramkey), MDF_FLOAT, (/mf%dimid_lon ,mf%dimid_latb,mf%dimid_lev_v,mf%dimid_time/), varid, status )
          IF_NOTOK_RETURN(status=1)
        else
          call MDF_Def_Var( mf%hid, trim(paramkey), MDF_FLOAT, (/mf%dimid_lon ,mf%dimid_lat ,mf%dimid_lev,mf%dimid_time/), varid, status )
          IF_NOTOK_RETURN(status=1)
        end if
      else
        ! 2D field:
        if ( trim(mf%tres) == 'constant' ) then
          call MDF_Def_Var( mf%hid, trim(paramkey), MDF_FLOAT, (/mf%dimid_lon,mf%dimid_lat/), varid, status )
          IF_NOTOK_RETURN(status=1)
        else
          call MDF_Def_Var( mf%hid, trim(paramkey), MDF_FLOAT, (/mf%dimid_lon,mf%dimid_lat,mf%dimid_time/), varid, status )
          IF_NOTOK_RETURN(status=1)
        end if
      end if

      ! convert from tm5 name to standard name:
      call TMM_CF_Convert_Name( trim(paramkey), mf%cfname_param(iparam), status )
      IF_NOTOK_RETURN(status=1)
      ! get standard units:
      call TMM_CF_Standard_Units( trim(mf%cfname_param(iparam)), mf%cfunit_param(iparam), status )
      IF_NOTOK_RETURN(status=1)
      ! add standard_name and units attributes:
      call CF_Put_Standard_Atts( mf%hid, varid, trim(mf%cfname_param(iparam)), trim(mf%cfunit_param(iparam)), status )
      IF_NOTOK_RETURN(status=1)
      ! for post processing ...
      l = len_trim(cell_measure)
      if ( l > 0 ) then
        call MDF_Put_Att( mf%hid, varid, 'cell_measure', cell_measure(2:l), status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! specify how computed:
      l = len_trim(cell_methods)
      if ( l > 0 ) then
        call MDF_Put_Att( mf%hid, varid, 'cell_methods', cell_methods(2:l), status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! extra coordinates that apply to this variable:
      l = len_trim(coordinates)
      if ( l > 0 ) then
        call MDF_Put_Att( mf%hid, varid, 'coordinates', coordinates(2:l), status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! store variable id:
      mf%varid_param(iparam) = varid
      ! store variable name:
      mf%varname_param(iparam) = trim(paramkey)
      ! no records written yet:
      mf%itrec_param(iparam) = 0

    end do

    !
    ! end of defintion phase:
    !

    call MDF_EndDef( mf%hid, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! fill time independent variables
    !

    ! axis:
    call MDF_Put_Var( mf%hid, mf%varid_lon, lli%lon_deg, status )
    IF_NOTOK_RETURN(status=1)
    ! axis bounds:
    call MDF_Put_Var( mf%hid, mf%varid_lon_bounds, lli%lon_bounds_deg, status )
    IF_NOTOK_RETURN(status=1)
    ! extra:
    if ( nuv /= 'n' ) then
      call MDF_Put_Var( mf%hid, mf%varid_lonb, lli%blon_deg, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! axis:
    call MDF_Put_Var( mf%hid, mf%varid_lat, lli%lat_deg, status )
    IF_NOTOK_RETURN(status=1)
    ! axis bounds:
    call MDF_Put_Var( mf%hid, mf%varid_lat_bounds, lli%lat_bounds_deg, status )
    IF_NOTOK_RETURN(status=1)
    ! extra:
    if ( nuv /= 'n' ) then
      call MDF_Put_Var( mf%hid, mf%varid_latb, lli%blat_deg, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! cell area ?
    if ( nuv == 'n' ) then
      ! storage:
      allocate( pat(lli%nlon,lli%nlat) )
      ! fill:
      call AreaOper( lli, pat, '=', 'm2', status )
      IF_NOTOK_RETURN(status=1)
      ! write:
      call MDF_Put_Var( mf%hid, mf%varid_cell_area, pat, status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( pat )
    end if

    ! levels ?
    if ( present(levi) ) then
      !! tm5 variables:
      !call MDF_Put_Var( mf%hid, mf%varid_tm5_at, levi%a, status )
      !IF_NOTOK_RETURN(status=1)
      !call MDF_Put_Var( mf%hid, mf%varid_tm5_bt, levi%b, status )
      !IF_NOTOK_RETURN(status=1)
      ! where defined ?
      select case ( nw )
        !~ layer mid points:
        case ( 'n' )!, '*' )
          ! number of layers:
          !if ( present(nlev) ) then
          !  n = nlev
          !else
            n = levi%nlev
          !end if
          ! standard pressures (full levels):
          call MDF_Put_Var( mf%hid, mf%varid_lev, levi%fp0(1:n), status )
          IF_NOTOK_RETURN(status=1)
          ! full level coefficients:
          call MDF_Put_Var( mf%hid, mf%varid_ap, levi%fa(1:n), status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Put_Var( mf%hid, mf%varid_b , levi%fb(1:n), status )
          IF_NOTOK_RETURN(status=1)
          ! half level coefficients:
          call MDF_Put_Var( mf%hid, mf%varid_ap_bounds, levi%fa_bounds(:,1:n), status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Put_Var( mf%hid, mf%varid_b_bounds , levi%fb_bounds(:,1:n), status )
          IF_NOTOK_RETURN(status=1)
          ! extra ?
          if ( nuv /= 'n' ) then
            ! standard pressures (full levels):
            call MDF_Put_Var( mf%hid, mf%varid_lev_u, levi%fp0(1:n), status )
            IF_NOTOK_RETURN(status=1)
            ! standard pressures (full levels):
            call MDF_Put_Var( mf%hid, mf%varid_lev_v, levi%fp0(1:n), status )
            IF_NOTOK_RETURN(status=1)
          end if
        !~ layer interfaces:
        case ( 'w' )
          ! standard pressures (half levels):
          call MDF_Put_Var( mf%hid, mf%varid_lev, levi%p0, status )
          IF_NOTOK_RETURN(status=1)
          ! half level coefficients:
          call MDF_Put_Var( mf%hid, mf%varid_ap, levi%a, status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Put_Var( mf%hid, mf%varid_b , levi%b, status )
          IF_NOTOK_RETURN(status=1)
        !~ other ...
        case default
          write (gol,'("unsupported nw : ",a)') trim(nw); call goErr
          TRACEBACK; status=1; return
      end select
    end if

    !
    ! done
    !

    ! ok
    status = 0

  end subroutine Define_File


  ! ***


  subroutine WriteTimes( mf, itrec, tref, t1, t2, status )

    use GO      , only : TDate
    use GO      , only : operator(-), operator(+), operator(/), rTotal, Get, NewDate, IsAnyDate
    use MDF     , only : MDF_Put_Var

    ! --- in/out -------------------------------

    type(TMeteoFile_tm5_nc), intent(inout)  ::  mf
    integer, intent(in)                  ::  itrec
    type(TDate), intent(in)              ::  tref, t1, t2
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/WriteTimes'

    ! --- local ------------------------------

    type(TDate)   ::  tsince
    type(TDate)   ::  tmid
    real(8)       ::  tsec
    integer       ::  time6(6)

    ! --- begin ---------------------------------

    ! not for constant fields ...
    if ( trim(mf%tres) /= 'constant' ) then

      ! base time:
      tsince = NewDate( time6=since_time6 )

      ! write reference times:
      if ( IsAnyDate(tref) ) then
        tsec = 0.0
        time6 = 0
      else
        tsec = rTotal( tref - tsince, 'sec' )
        call Get( tref, time6=time6 )
      end if
      call MDF_Put_Var( mf%hid, mf%varid_reftime, (/tsec/), status, start=(/itrec/), count=(/1/) )
      IF_NOTOK_RETURN(status=1)
      call MDF_Put_Var( mf%hid, mf%varid_reftimevalues, time6, status, start=(/1,itrec/), count=(/6,1/) )
      IF_NOTOK_RETURN(status=1)

      ! write mid time:
      if ( IsAnyDate(t1) .or. IsAnyDate(t2) ) then
        tsec = 0.0
        time6 = 0
      else
        tmid = t1 + (t2-t1)/2
        tsec = rTotal( tmid - tsince, 'sec' )
        call Get( tmid, time6=time6 )
      end if
      call MDF_Put_Var( mf%hid, mf%varid_time, (/tsec/), status, start=(/itrec/), count=(/1/) )
      IF_NOTOK_RETURN(status=1)
      call MDF_Put_Var( mf%hid, mf%varid_timevalues, time6, status, start=(/1,itrec/), count=(/6,1/) )
      IF_NOTOK_RETURN(status=1)

      ! start time:
      if ( IsAnyDate(t1) ) then
        tsec = 0.0
        time6 = 0
      else
        tsec = rTotal( t1 - tsince, 'sec' )
        call Get( t1, time6=time6 )
      end if
      call MDF_Put_Var( mf%hid, mf%varid_time_bounds, (/tsec/), status, start=(/1,itrec/), count=(/1,1/) )
      IF_NOTOK_RETURN(status=1)
      call MDF_Put_Var( mf%hid, mf%varid_timevalues_bounds, time6, status, start=(/1,1,itrec/), count=(/1,6,1/) )
      IF_NOTOK_RETURN(status=1)

      ! end time:
      if ( IsAnyDate(t2) ) then
        tsec = 0.0
        time6 = 0
      else
        tsec = rTotal( t2 - tsince, 'sec' )
        call Get( t2, time6=time6 )
      end if
      call MDF_Put_Var( mf%hid, mf%varid_time_bounds, (/tsec/), status, start=(/2,itrec/), count=(/1,1/) )
      IF_NOTOK_RETURN(status=1)
      call MDF_Put_Var( mf%hid, mf%varid_timevalues_bounds, time6, status, start=(/2,1,itrec/), count=(/1,6,1/) )
      IF_NOTOK_RETURN(status=1)

    end if   ! not const

    ! ok
    status = 0

  end subroutine WriteTimes


  ! ***


  subroutine mf_WriteRecord_2d( mf, tmi, paramkey, unit, tref, t1, t2, &
                                lli, nuv, ll, status )

    use GO      , only : TDate
    use Grid    , only : TllGridInfo
    use MDF     , only : MDF_Create, MDF_Close, MDF_REPLACE
    use MDF     , only : MDF_Put_Var
    use tmm_info, only : TMeteoInfo
    use TMM_CF  , only : TMM_CF_Convert_Units
    use misctools,  only : check_dir

    ! --- in/out -------------------------------

    type(TMeteoFile_tm5_nc), intent(inout)  ::  mf
    type(TMeteoInfo), intent(in)         ::  tmi
    character(len=*), intent(in)         ::  paramkey, unit
    type(TDate), intent(in)              ::  tref, t1, t2
    type(TllGridInfo), intent(in)        ::  lli
    character(len=1), intent(in)         ::  nuv
    real, intent(in)                     ::  ll(:,:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_WriteRecord_2d'

    ! --- local ------------------------------

    integer       ::  iparam
    integer       ::  i
    integer       ::  itrec
    type(TDate)   ::  tsince
    type(TDate)   ::  tmid
    real(8)       ::  tsec
    integer       ::  time6(6)
    real          ::  ufac

    ! --- begin ---------------------------------

    ! output ?
    if ( mf%io /= 'o' ) then
      write (gol,'("file should have been opened for output, but io=",a)') mf%io; call goErr
      TRACEBACK; status=1; return
    end if

    ! new or existing ?
    if ( .not. mf%output_initialised ) then
      ! if directory does not exist, create it
      call check_dir(mf%fname)
      ! open new file, destroy old:
      call MDF_Create( trim(mf%fname), output_type, MDF_REPLACE, mf%hid, status )
      IF_NOTOK_RETURN(status=1)
!      write(gol,'(a," created file ",a)') rname, trim(mf%fname) ; call goPr
      ! write file header:
      call Define_File( mf, lli, nuv, status )
      IF_NOTOK_RETURN(status=1)
      ! status new
      call WriteStatus( mf, 'in-progress', status )
      IF_NOTOK_RETURN(status=1)
      ! no records written yet:
      mf%output_nrec = 0
      ! now the file is initialised
      mf%output_initialised = .true.
    endif

    ! search parameter:
    iparam = -1
    do i = 1, mf%nparam
      if ( trim(paramkey) == trim(mf%varname_param(i)) ) then
        iparam = i
        exit
      end if
    end do
    if ( iparam < 0 ) then
      write (gol,'("could not find parameter `",a,"` in list : ",a)') trim(paramkey), trim(mf%paramkeys); call goErr
      TRACEBACK; status=1; return
    end if

    ! next record:
    mf%itrec_param(iparam) = mf%itrec_param(iparam) + 1
    ! short:
    itrec = mf%itrec_param(iparam)

    ! update time fields:
    call WriteTimes( mf, itrec, tref, t1, t2, status )
    IF_NOTOK_RETURN(status=1)

    ! conversion factor:
    call TMM_CF_Convert_Units( trim(unit), trim(mf%cfunit_param(iparam)), ufac, status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    if ( ufac /= 1.0 ) then
      write (gol,'("      convert `",a,"` from `",a,"` to `",a,"` with factor ",f8.2," ; write range ",2f8.2)') &
               trim(paramkey), trim(unit), trim(mf%cfunit_param(iparam)), ufac, minval(ll*ufac), maxval(ll*ufac); call goPr
      !if ( trim(paramkey) == 'tv01' ) stop 'break after tv01'
    end if

    ! write data:
    if ( trim(mf%tres) == 'constant' ) then
      call MDF_Put_Var( mf%hid, mf%varid_param(iparam), ll*ufac, status )
      IF_NOTOK_RETURN(status=1)
    else
      call MDF_Put_Var( mf%hid, mf%varid_param(iparam), ll*ufac, status, &
                            start=(/1,1,itrec/), count=(/size(ll,1),size(ll,2),1/) )
      IF_NOTOK_RETURN(status=1)
    end if

    ! next record has been written:
    mf%output_nrec = mf%output_nrec + 1

    ! completed ? then re-write status file:
    if ( mf%output_nrec == mf%output_ntrec*mf%nparam ) then
      ! close file:
      call MDF_Close( mf%hid, status )
      IF_NOTOK_RETURN(status=1)
      ! rewrite status file:
      call WriteStatus( mf, 'completed', status )
      IF_NOTOK_RETURN(status=1)
      ! reset flag:
      mf%output_initialised = .false.
!      write(gol,'(a," for file ",a,", output un-initialized")') rname, trim(mf%fname) ; call goPr
!      write(gol,'(a," output_nrec = ",i6, ", output_ntrec = ", i6, ", nparam = ", i6)') rname, mf%output_nrec, mf%output_ntrec, mf%nparam ; call goPr
!      write(gol,'(a)') "  " ; call goPr
    end if

    ! ok
    status = 0

  end subroutine mf_WriteRecord_2d


  ! ***


  subroutine mf_WriteRecord_3d( mf, tmi, spname, paramkey, unit, tref, t1, t2, &
                                lli, nuv, levi, nw, ps, ll, status )!, &
                                !nlev )

    use GO      , only : TDate
    use Grid    , only : TllGridInfo, TLevelInfo
    use MDF     , only : MDF_Create, MDF_Close, MDF_REPLACE
    use MDF     , only : MDF_Put_Var
    use tmm_info, only : TMeteoInfo
    use TMM_CF  , only : TMM_CF_Convert_Units
    use misctools,  only : check_dir

    ! --- in/out -------------------------------

    type(TMeteoFile_tm5_nc), intent(inout)  ::  mf
    type(TMeteoInfo), intent(in)         ::  tmi
    character(len=*), intent(in)         ::  spname, paramkey, unit
    type(TDate), intent(in)              ::  tref, t1, t2
    type(TllGridInfo), intent(in)        ::  lli
    character(len=1), intent(in)         ::  nuv
    type(TLevelInfo), intent(in)         ::  levi
    character(len=1), intent(in)         ::  nw
    real, intent(in)                     ::  ps(:,:)
    real, intent(in)                     ::  ll(:,:,:)
    integer, intent(out)                 ::  status

    !integer, intent(in), optional        ::  nlev

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_WriteRecord_3d'

    ! --- local ------------------------------

    integer             ::  iparam
    integer             ::  i
    integer             ::  itrec
    type(TDate)         ::  tsince
    type(TDate)         ::  tmid
    real(8)             ::  tsec
    integer             ::  time6(6)
    real                ::  ufac
    real, allocatable   ::  llX(:,:,:)
    character(len=1)    ::  nwX

    ! --- begin ---------------------------------

    ! output ?
    if ( mf%io /= 'o' ) then
      write (gol,'("file should have been opened for output, but io=",a)') mf%io; call goErr
      TRACEBACK; status=1; return
    end if

    ! convection fields provided on lowest layers only ...
    if ( nw == '*' ) then
      ! extend to all layers:
      allocate( llX(size(ll,1),size(ll,2),levi%nlev) )
      ! pad with zeros:
      llX = 0.0
      llX(:,:,1:size(ll,3)) = ll
      ! new description:
      nwX = 'n'
    else
      ! just copy:
      allocate( llX(size(ll,1),size(ll,2),size(ll,3)) )
      llX = ll
      nwX = nw
    end if

    ! new or existing ?
    if ( .not. mf%output_initialised ) then
      ! if directory does not exist, create it
      call check_dir(mf%fname)
      ! open new file, destroy old:
      call MDF_Create( trim(mf%fname), output_type, MDF_REPLACE, mf%hid, status )
      IF_NOTOK_RETURN(status=1)
      ! write file header:
      call Define_File( mf, lli, nuv, status, levi=levi, nw=nwX )!, nlev=nlev )
      IF_NOTOK_RETURN(status=1)
      ! status new
      call WriteStatus( mf, 'in-progress', status )
      IF_NOTOK_RETURN(status=1)
      ! no records written yet:
      mf%output_nrec = 0
      ! now the file is initialised
      mf%output_initialised = .true.
    endif

    ! search parameter:
    iparam = -1
    do i = 1, mf%nparam
      if ( trim(paramkey) == trim(mf%varname_param(i)) ) then
        iparam = i
        exit
      end if
    end do
    if ( iparam < 0 ) then
      write (gol,'("could not find parameter `",a,"` in list : ",a)') trim(paramkey), trim(mf%paramkeys); call goErr
      TRACEBACK; status=1; return
    end if

    ! next record:
    mf%itrec_param(iparam) = mf%itrec_param(iparam) + 1
    ! short:
    itrec = mf%itrec_param(iparam)

    ! update time fields:
    call WriteTimes( mf, itrec, tref, t1, t2, status )
    IF_NOTOK_RETURN(status=1)

    ! conversion factor:
    call TMM_CF_Convert_Units( trim(unit), trim(mf%cfunit_param(iparam)), ufac, status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    if ( ufac /= 1.0 ) then
      write (gol,'("      convert `",a,"` from `",a,"` to `",a,"` with factor ",f8.2," ; write range ",2f8.2)') &
               trim(paramkey), trim(unit), trim(mf%cfunit_param(iparam)), ufac, minval(llX*ufac), maxval(llX*ufac); call goPr
    end if

    ! write data:
    call MDF_Put_Var( mf%hid, mf%varid_param(iparam), llX*ufac, status, &
                          start=(/1,1,1,itrec/), count=(/size(llX,1),size(llX,2),size(llX,3),1/) )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( llX )

    ! write surface pressure:
    if ( nuv == 'u' ) then
      call MDF_Put_Var( mf%hid, mf%varid_ps_u, ps, status, &
                            start=(/1,1,itrec/), count=(/size(ps,1),size(ps,2),1/) )
      IF_NOTOK_RETURN(status=1)
    else if ( nuv == 'v' ) then
      call MDF_Put_Var( mf%hid, mf%varid_ps_v, ps, status, &
                            start=(/1,1,itrec/), count=(/size(ps,1),size(ps,2),1/) )
      IF_NOTOK_RETURN(status=1)
    else
      call MDF_Put_Var( mf%hid, mf%varid_ps, ps, status, &
                            start=(/1,1,itrec/), count=(/size(ps,1),size(ps,2),1/) )
      IF_NOTOK_RETURN(status=1)
    end if

    ! next record has been written:
    mf%output_nrec = mf%output_nrec + 1

    ! completed ? then re-write status file:
    if ( mf%output_nrec == mf%output_ntrec*mf%nparam ) then
      ! close file:
      call MDF_Close( mf%hid, status )
      IF_NOTOK_RETURN(status=1)
      ! rewrite status file:
      call WriteStatus( mf, 'completed', status )
      IF_NOTOK_RETURN(status=1)
      ! reset flag:
      mf%output_initialised = .false.
    end if

    ! ok
    status = 0

  end subroutine mf_WriteRecord_3d


  ! ***


  subroutine WriteStatus( mf, msg, status )

    ! --- in/out -------------------------------

    type(TMeteoFile_tm5_nc), intent(inout)  ::  mf
    character(len=*), intent(in)         ::  msg
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/WriteStatus'

    ! --- local ------------------------------

    integer       ::  fu
    logical       ::  opened

    ! --- begin ---------------------------------

    ! select unused file unit:
    fu = 1234
    do
      inquire( unit=fu, opened=opened )
      if ( .not. opened ) exit
      fu = fu + 1
    end do

    ! open:
    open( fu, file=trim(mf%fname)//'.status', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening status file:")'); call goErr
      write (gol,'("  file   : ",a)') trim(mf%fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! write message:
    write (fu,'(a)',iostat=status) msg
    if (status/=0) then
      write (gol,'("writing status:")'); call goErr
      write (gol,'("  file   : ",a)') trim(mf%fname); call goErr
      write (gol,'("  msg    : ",a)') msg; call goErr
      TRACEBACK; status=1; return
    end if

    ! done:
    close( fu, iostat=status )
    if (status/=0) then
      write (gol,'("closing status file:")'); call goErr
      write (gol,'("  file   : ",a)') trim(mf%fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine WriteStatus


end module tmm_mf_tm5_nc
