!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module tmm_mf

  use GO  , only : gol, goErr, goPr, goBug
  use GO  , only : TDate
  use Grid, only : TllGridInfo, TggGridInfo, TshGridInfo, TshGrid, TLevelInfo
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN, DUMMY_STR_LEN

#ifdef with_tmm_tmpp
  use tmm_mf_tmpp      , only : TMeteoFile_tmpp
#endif
#ifdef with_tmm_tm5
  use tmm_mf_tm5_hdf    , only : TMeteoFile_tm5_hdf
  use tmm_mf_tm5_nc     , only : TMeteoFile_tm5_nc
#endif
#ifdef with_tmm_ecmwf
  use tmm_mf_ecmwf_tmpp, only : TMeteoFile_ecmwf_tmpp
  use tmm_mf_ecmwf_tm5 , only : TMeteoFile_ecmwf_tm5
#endif
#ifdef with_tmm_ncep
  use tmm_mf_ncep_cdc  , only : TMeteoFile_ncep_cdc
  use tmm_mf_ncep_gfs  , only : TMeteoFile_ncep_gfs
#endif
#ifdef with_prism
  use tmm_mf_prism     , only : TMeteoFile_prism
#endif
#ifdef with_tmm_msc
  use tmm_mf_msc       , only : TMeteoFile_msc
#endif

  implicit none

  ! --- in/out -------------------------------------

  private

  public   ::  TMeteoFile
  public   ::  Init, Done
  public   ::  Opened, CheckTime, CheckParam

  public   ::  SetupInput
  public   ::  ReadRecord

  public   ::  SetupOutput
  public   ::  WriteRecord


  ! --- const ---------------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf'

  ! --- types --------------------------------------

  type TMeteoFile
    ! opened yet ?
    logical                    ::  opened
    character(len=1)           ::  io             = ''
    ! meteo archive keys:
    character(len=MAX_FILENAME_LEN)         ::  dir            = ''
    character(len=MAX_RCKEY_LEN)         ::  archivekey     = ''
    ! parameter keys for fields in this file
    character(len=MAX_RCKEY_LEN)         ::  paramkeys      = ''
    ! time range for which file is valid
    type(TDate)                ::  t1, t2
    !
    ! access to current meteo file
    !
    character(len=10)          ::  filetype       = ''
    character(len=MAX_FILENAME_LEN)         ::  filename       = ''
    character(len=MAX_FILENAME_LEN)         ::  spm_filename   = ''
#ifdef with_tmm_tmpp
    type(TMeteoFile_tmpp)          ::  mf_tmpp         ! tmpp written hdf file
#endif
#ifdef with_tmm_tm5
    type(TMeteoFile_tm5_hdf)       ::  mf_tm5_hdf      ! tm5 written hdf file
    type(TMeteoFile_tm5_nc)        ::  mf_tm5_nc       ! tm5 written netcdf file
#endif
#ifdef with_tmm_ecmwf
    type(TMeteoFile_ecmwf_tmpp)    ::  mf_ecmwf_tmpp   ! grib file retrieved with tmpp
    type(TMeteoFile_ecmwf_tm5)     ::  mf_ecmwf_tm5    ! grib file retrieved with tm5
#endif
#ifdef with_tmm_ncep
    type(TMeteoFile_ncep_cdc)      ::  mf_ncep_cdc     ! ncep file from cdc archive
    type(TMeteoFile_ncep_gfs)      ::  mf_ncep_gfs     ! ncep gfs file
#endif
#ifdef with_prism
    type(TMeteoFile_prism)         ::  mf_prism        ! prism file
#endif
#ifdef with_tmm_msc
    type(TMeteoFile_msc)           ::  mf_msc          ! msc file
#endif
  end type TMeteoFile


  ! --- interfaces -------------------------------

  interface Init
    module procedure mf_Init
  end interface

  interface Done
    module procedure mf_Done
  end interface

  interface Opened
    module procedure mf_Opened
  end interface

  interface CheckTime
    module procedure mf_CheckTime
  end interface

  interface CheckParam
    module procedure mf_CheckParam
  end interface

  interface SetupInput
    module procedure mf_SetupInput
  end interface

  interface ReadRecord
    module procedure mf_ReadRecord
  end interface

!  interface ReadEqvLatStuff
!    module procedure mf_ReadEqvLatStuff
!  end interface

  interface SetupOutput
    module procedure mf_SetupOutput
  end interface

  interface WriteRecord
    module procedure mf_WriteRecord_2d
    module procedure mf_WriteRecord_3d
  end interface


contains


  ! ===========================================================
  !
  !   init/done
  !
  ! ===========================================================


  subroutine mf_Init( mf, io, status )

    ! --- begin -------------------------------------------

    type(TMeteoFile), intent(out)           ::  mf
    character(len=1), intent(in)            ::  io
    integer, intent(out)                    ::  status

    ! --- const -------------------------------------------

    character(len=*), parameter ::  rname = mname//'/mf_Init'

    ! --- begin -------------------------------------------

    ! input or output ?
    mf%io = io

    ! file not opened yet
    mf%opened = .false.

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    use GO, only : goSystem
    use Grid, only : Done

#ifdef with_tmm_tmpp
    use tmm_mf_tmpp      , only : Done
#endif
#ifdef with_tmm_tm5
    use tmm_mf_tm5_hdf    , only : Done
    use tmm_mf_tm5_nc     , only : Done
#endif
#ifdef with_tmm_ecmwf
    use tmm_mf_ecmwf_tmpp, only : Done
    use tmm_mf_ecmwf_tm5 , only : Done
#endif
#ifdef with_tmm_ncep
    use tmm_mf_ncep_cdc  , only : Done
    use tmm_mf_ncep_gfs  , only : Done
#endif
#ifdef with_prism
    use tmm_mf_prism     , only : Done
#endif
#ifdef with_tmm_msc
    use tmm_mf_msc       , only : Done
#endif

    ! --- begin -------------------------------------------

    type(TMeteoFile), intent(inout)         ::  mf
    integer, intent(out)                    ::  status

    ! --- const -------------------------------------------

    character(len=*), parameter ::  rname = mname//'/mf_Done'

    ! --- begin -------------------------------------------

    ! close file if necessary
    if ( mf%opened ) then
      select case ( mf%filetype )
#ifdef with_tmm_tmpp
        case ( 'tmpp' )
          call Done( mf%mf_tmpp, status )
          IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_tmm_tm5
        case ( 'hdf', 'tm5-hdf' )
          call Done( mf%mf_tm5_hdf, status )
          IF_NOTOK_RETURN(status=1)
        case ( 'tm5-nc' )
          call Done( mf%mf_tm5_nc, status )
          IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_tmm_ecmwf
        case ( 'ecmwf-tmpp' )
          call Done( mf%mf_ecmwf_tmpp, status )
          IF_NOTOK_RETURN(status=1)
        case ( 'ecmwf-tm5' )
          call Done( mf%mf_ecmwf_tm5, status )
          IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_tmm_ncep
        case ( 'ncep-cdc' )
          call Done( mf%mf_ncep_cdc, status )
          IF_NOTOK_RETURN(status=1)
        case ( 'ncep-gfs' )
          call Done( mf%mf_ncep_gfs, status )
          IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_prism
        case ( 'prism' )
          call Done( mf%mf_prism, status )
          IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_tmm_msc
        case ( 'msc-data' )
          call Done( mf%mf_msc, status )
          IF_NOTOK_RETURN(status=1)
#endif
        case default
          write (gol,'("unsupported filetype `",a,"`")') trim(mf%filetype); call goErr
          TRACEBACK; status=1; return
      end select
      mf%opened = .false.
    end if

    ! ok
    status = 0

  end subroutine mf_Done


  ! ***


  logical function mf_Opened( mf )

    ! --- begin -------------------------------------------

    type(TMeteoFile), intent(in)         ::  mf

    ! --- begin -------------------------------------------

    mf_Opened = mf%opened

  end function mf_Opened


  ! ===========================================================
  !
  !   check contents of open meteo file
  !
  ! ===========================================================


  ! Check time in meteo file;
  ! status:
  !   <0    : mf does not include [t1,t2]
  !    0    : mf includes [t,t2]
  !   >0    : error; mf not open ?
  !

  subroutine mf_CheckTime( mf, t1, t2, status )

    use GO, only : TDate, IncrDate, wrtgol, IsAnyDate
    use GO, only : operator(+), operator(-), operator(==), operator(<), operator(<=)

    ! --- begin -------------------------------------------

    type(TMeteoFile), intent(in)           ::  mf
    type(TDate), intent(in)                ::  t1, t2
    integer, intent(out)                   ::  status

    ! --- const -------------------------------------------

    character(len=*), parameter ::  rname = mname//'/mf_CheckTime'

    ! --- begin -------------------------------------------

    ! not open ?
    if ( .not. Opened(mf) ) then
      write (gol,'("meteo file not opened")'); call goErr
      TRACEBACK; status = 1; return
    end if

    ! trap any date:
    if ( IsAnyDate(t1) .and. IsAnyDate(t2) ) then
      status =  0; return
    end if

    ! [t1,t2] is either:
    !   covered by mf   ->  status =  0
    !   older than mf   ->  status = -2
    !   newer than mf   ->  status = -1
    !   error ... (half in, half outside mf)
    ! seperate tests for intervals and instant time:
    if ( t1 == t2 ) then

      ! instant time
      if ( ( (mf%t1 <= t1) .and. (t1 <= mf%t2) ) ) then
        status =  0; return
      else if ( t1 < mf%t1 ) then
        status = -2; return
      else if ( mf%t2 <= t1 ) then
        status = -1; return
      else
        write (gol,'("requested instant time t1 (=t2) overlaps part of mf time:")'); call goErr
        call wrtgol( '  t1     : ', t1 ); call goErr
        call wrtgol( '  t2     : ', t2 ); call goErr
        call wrtgol( '  mf%t1  : ', mf%t1 ); call goErr
        call wrtgol( '  mf%t2  : ', mf%t2 ); call goErr
        write (gol,'("  params : ",a)') trim(mf%paramkeys); call goErr
        TRACEBACK; status = 1; return
      end if

    else if ( t1 < t2 ) then

      ! interval
      ! extra: [t1,t2] is covered by mf%(t1,t2) ...
      if ( ( (mf%t1                  <= t1) .and. (t2 <= mf%t2                 ) ) .or. &
           ( (mf%t1-IncrDate(mili=1) <= t1) .and. (t2 <= mf%t2+IncrDate(mili=1)) )      ) then
        status =  0; return
      else if ( t2 <= mf%t1 ) then
        ! request for field older than those in file
        status = -2; return
      else if ( mf%t2 <= t1 ) then
        ! request for field newer than those in file
        status = -1; return
      else
        write (gol,'("requested interval [t1,t2] overlaps part of mf time:")'); call goErr
        call wrtgol( '  t1     : ', t1 ); call goErr
        call wrtgol( '  t2     : ', t2 ); call goErr
        call wrtgol( '  mf%t1  : ', mf%t1 ); call goErr
        call wrtgol( '  mf%t2  : ', mf%t2 ); call goErr
        write (gol,'("  params : ",a)') trim(mf%paramkeys); call goErr
        TRACEBACK; status = 1; return
      end if

    else

      write (gol,'("arguments should specify an instant time or valid interval :")'); call goErr
      call wrtgol( '  t1 : ', t1 ); call goErr
      call wrtgol( '  t2 : ', t2 ); call goErr
      TRACEBACK; status = 1; return

    end if

    ! something wrong if this point is reached ...
    status = 1

  end subroutine mf_CheckTime


  ! ***


  ! Check if param is included in meteo file;
  ! status:
  !   <0    : mf does not include param
  !    0    : mf includes param
  !   >0    : error; mf not open ?
  !

  subroutine mf_CheckParam( mf, io, archivekey, paramkey, status )

    use GO, only : goLoCase

    ! --- begin -------------------------------------------

    type(TMeteoFile), intent(in)           ::  mf
    character(len=*), intent(in)           ::  io
    character(len=*), intent(in)           ::  archivekey
    character(len=*), intent(in)           ::  paramkey
    integer, intent(out)                   ::  status

    ! --- const -------------------------------------------

    character(len=*), parameter ::  rname = mname//'/mf_CheckParam'

    ! --- local --------------------------------------------

    integer             ::  pos

    ! --- begin -------------------------------------------

    ! not open ?
    if ( .not. Opened(mf) ) then
      write (gol,'("meteo file not opened")'); call goErr
      TRACEBACK; status = 1; return
    end if

    ! by default not found ..
    status = -1

    ! wrong input/output ? then leave:
    if ( io /= mf%io ) return

    ! wrong grid ? then leave
    if ( archivekey /= mf%archivekey ) return

    ! param list is for example: '-ps-pu-pv-',
    ! thus search for example for '-pu-' ...
    ! convert all to lowercase
    pos = index( goLoCase(trim(mf%paramkeys)), '-'//goLoCase(trim(paramkey))//'-' )
    if ( pos < 1 ) return

    ! ok
    status = 0

  end subroutine mf_CheckParam



  ! ===========================================================
  !
  !   open meteo file for input
  !
  ! ===========================================================


  !
  ! Open the meteo file that contains the field specified by
  ! archivekey, parameter key, time,
  ! or do nothing if the requested file has been opened already.
  !
  ! <archivekey> = <archivetype>:<archivename>
  !
  !    tmpp:od-fc-ml60-glb3x2
  !   tmppS:od-fc-ml60-glb3x2
  !    grib:od-fc-ml60-glb3x2
  !   prism:
  !

  subroutine mf_SetupInput( mf, archivekey, paramkey, tday, t1, t2, &
                                rcfilename, dir, status )

    use GO, only : goSplitLine, goReadFromLine
    use GO, only : goSystem
    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : TDate, IncrDate, Get, NewDate, wrtgol, &
                   Operator(+), Operator(-), Operator(/)
#ifdef with_tmm_tmpp
    use tmm_mf_tmpp      , only : Init, Get
#endif
#ifdef with_tmm_tm5
    use tmm_mf_tm5_hdf    , only : Init, Get
    use tmm_mf_tm5_nc     , only : Init, Get
#endif
#ifdef with_tmm_ecmwf
    use tmm_mf_ecmwf_tmpp, only : Init, Get
    use tmm_mf_ecmwf_tm5 , only : Init, Get
#endif
#ifdef with_tmm_ncep
    use tmm_mf_ncep_cdc  , only : Init, Get
    use tmm_mf_ncep_gfs  , only : Init, Get
#endif
#ifdef with_prism
    use tmm_mf_prism     , only : Init
#endif
#ifdef with_tmm_msc
    use tmm_mf_msc       , only : Init, Get
#endif

    ! --- in/out -------------------------------------

    type(TMeteoFile), intent(inout)        ::  mf
    character(len=*), intent(in)           ::  archivekey
    character(len=*), intent(in)           ::  paramkey
    type(TDate), intent(in)                ::  tday, t1, t2
    character(len=*), intent(in)           ::  rcfilename
    character(len=*), intent(in)           ::  dir
    integer, intent(inout)                 ::  status

    ! --- const -------------------------------------

    character(len=*), parameter ::  rname = mname//'/mf_SetupInput'

    ! name of info file:
    character(len=*), parameter  ::  infofilename = 'tmm_info.rc'

    ! --- local -------------------------------------

    character(len=10)        ::  archivetype
    character(len=MAX_FILENAME_LEN) ::  archivename

    character(len=DUMMY_STR_LEN)       ::  command
    integer                  ::  year1, month1, day1, hour1
    integer                  ::  year2, month2, day2, hour2
    integer                  ::  dth
    type(TrcFile)            ::  infofile

    character(len=MAX_FILENAME_LEN)       ::  archivename2
    character(len=10)        ::  mclass
    character(len=10)        ::  mtype
    character(len=10)        ::  mlevs
    character(len=10)        ::  mgrid
    character(len=10)        ::  filekey
    character(len=16)        ::  treskey
    logical                  ::  with_spm
    logical                  ::  constant

    ! --- begin -------------------------------------

    ! store archive key:
    mf%archivekey = trim(archivekey)
    mf%dir        = trim(dir)

    ! split archive key in type and name:
    call goSplitLine( archivekey, archivetype, ':', archivename, status )
    IF_NOTOK_RETURN(status=1)

    ! usually, meteo is storred in file;
    ! for PRISM project, meteo is in the memory ...
    select case ( archivetype )

#ifdef with_tmm_tmpp

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! hdf files written by tmpp
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tmpp' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'tmpp'

        ! setup file:
        call Init( mf%mf_tmpp, 'i', dir, archivename, paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_tmpp%fname

        ! extract time range:
        call Get( mf%mf_tmpp, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract paramkeys for fields in file:
        call Get( mf%mf_tmpp, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

#endif


#ifdef with_tmm_tm5

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! hdf or netcdf files written by tm5
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tm5-hdf' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'tm5-hdf'

        ! setup file:
        call Init( mf%mf_tm5_hdf, 'i', dir, archivename, paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_tm5_hdf%fname

        ! extract time range:
        call Get( mf%mf_tm5_hdf, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract paramkeys for fields in file:
        call Get( mf%mf_tm5_hdf, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

      case ( 'tm5-nc' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'tm5-nc'

        ! setup file:
        call Init( mf%mf_tm5_nc, 'i', dir, archivename//';form='//trim(archivetype), paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_tm5_nc%fname

        ! extract time range:
        call Get( mf%mf_tm5_nc, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract paramkeys for fields in file:
        call Get( mf%mf_tm5_nc, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

#endif


#ifdef with_tmm_ecmwf

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ecmwf grib files
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tmpp' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'ecmwf-tmpp'

        ! setup file:
        call Init( mf%mf_ecmwf_tmpp, dir, archivename, paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_ecmwf_tmpp%fname

        ! extract time range:
        call Get( mf%mf_ecmwf_tmpp, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract list of parameters in files:
        call Get( mf%mf_ecmwf_tmpp, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

      case ( 'ecmwf-tm5' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'ecmwf-tm5'

        ! setup file:
        call Init( mf%mf_ecmwf_tm5, dir, trim(archivename), paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_ecmwf_tm5%fname

        ! extract time range:
        call Get( mf%mf_ecmwf_tm5, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract list of parameters in files:
        call Get( mf%mf_ecmwf_tm5, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

#endif


#ifdef with_tmm_ncep

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ncep files
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ncep-cdc' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'ncep-cdc'

        ! setup file:
        call Init( mf%mf_ncep_cdc, dir, archivename, paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_ncep_cdc%fname

        ! extract time range:
        call Get( mf%mf_ncep_cdc, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract list of parameters in files:
        call Get( mf%mf_ncep_cdc, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

      case ( 'ncep-gfs' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'ncep-gfs'

        ! setup file:
        call Init( mf%mf_ncep_gfs, dir, archivename, paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_ncep_gfs%fname

        ! extract time range:
        call Get( mf%mf_ncep_gfs, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract list of parameters in files:
        call Get( mf%mf_ncep_gfs, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

#endif


#ifdef with_tmm_msc

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! msc-data text files
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'msc-data' )

        ! wich of the 'mf%mf_???' is used ?
        mf%filetype = 'msc-data'

        ! setup file:
        call Init( mf%mf_msc, dir, archivename, paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        mf%filename = mf%mf_msc%fname

        ! extract time range:
        call Get( mf%mf_msc, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract which fields are stored in the file:
        call Get( mf%mf_msc, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

#endif

#ifdef with_prism

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! prism meteo in memory
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'prism' )

        ! only the requested parameter is provided by this prism 'file' ...
        mf%paramkeys = '-'//paramkey//'-'

        ! infinite time range ...
        mf%t1 = NewDate( year=1900, month=1, day=1, hour=1 )
        mf%t2 = NewDate( year=9999, month=9, day=9, hour=9 )
        !call wrtgol( '     fields valid from : ', mf%t1 )
        !call wrtgol( '                  to   : ', mf%t2 )

        ! set file type and file name:
        mf%filetype = 'prism'
        mf%filename = 'dummy'

        ! setup prism access;
        ! tday is used for orography date
        ! (adhoc solution; at the moment only [t1,t2] is provided to ReadRecord
        ! but this should become tday, [t1,t2] )
        call Init( mf%mf_prism, tday, status )
        IF_NOTOK_RETURN(status=1)

#endif

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error ...
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default

        write (gol,'("unsupported archivetype `",a,"` for parameter ", a)') trim(archivetype), trim(paramkey); call goErr
        TRACEBACK; status=1; return

    end select

    ! file is opened (or, at least file name is known)
    mf%opened = .true.

    ! ok
    status = 0

  end subroutine mf_SetupInput



  ! ===========================================================
  !
  !   read fields, grid definition, etc
  !
  ! ===========================================================


  subroutine mf_ReadRecord( mf, paramkey, unit, tday, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use GO             , only : TDate, operator(+), operator(-), operator(/)
    use Grid           , only : TllGridInfo, TggGridInfo, TshGridInfo, TLevelInfo
    use tmm_info       , only : TMeteoInfo, Init, AddHistory
#ifdef with_tmm_tmpp
    use tmm_mf_tmpp      , only : ReadRecord
#endif
#ifdef with_tmm_tm5
    use tmm_mf_tm5_hdf    , only : ReadRecord
    use tmm_mf_tm5_nc     , only : ReadRecord
#endif
#ifdef with_tmm_ecmwf
    use tmm_mf_ecmwf_tmpp, only : ReadRecord
    use tmm_mf_ecmwf_tm5 , only : ReadRecord
#endif
#ifdef with_tmm_ncep
    use tmm_mf_ncep_cdc  , only : ReadRecord
    use tmm_mf_ncep_gfs  , only : ReadRecord
#endif
#ifdef with_prism
    use tmm_mf_prism     , only : ReadRecord
#endif
#ifdef with_tmm_msc
    use tmm_mf_msc       , only : ReadRecord
#endif

    ! --- in/out -------------------------------

    type(TMeteoFile), intent(inout)      ::  mf
    character(len=*), intent(in)         ::  paramkey
    character(len=*), intent(in)         ::  unit
    type(TDate), intent(in)              ::  tday, t1, t2
    character(len=1), intent(in)         ::  nuv
    character(len=1), intent(in)         ::  nw

    character(len=2), intent(out)        ::  gridtype
    type(TLevelInfo), intent(out)        ::  levi
    type(TllGridInfo), intent(inout)     ::  lli
    real, pointer                        ::  ll(:,:,:)
    real, pointer                        ::  sp_ll(:,:)
    type(TggGridInfo), intent(inout)     ::  ggi
    real, pointer                        ::  gg(:,:)
    real, pointer                        ::  sp_gg(:)
    type(TshGridInfo), intent(inout)     ::  shi
    complex, pointer                     ::  sh(:,:)
    complex, pointer                     ::  lnsp_sh(:)
    type(TMeteoInfo), intent(out)        ::  tmi
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord'

    ! --- local --------------------------------------

    type(TDate)         ::  tmid

    ! --- begin ---------------------------------

    !write (*,'(a,": begin")') name

    select case ( mf%filetype )
#ifdef with_tmm_tmpp
      case ( 'tmpp' )
        ! read from hdf file:
        call ReadRecord( mf%mf_tmpp, paramkey, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           status )
        IF_NOTOK_RETURN(status=1)
        ! fill some info values:
        call Init( tmi, paramkey, 'unkown', status )
        call AddHistory( tmi, 'archivekey=='//trim(mf%archivekey), status )
#endif
#ifdef with_tmm_tm5
      case ( 'hdf', 'tm5-hdf' )
        ! read from hdf file:
        call ReadRecord( mf%mf_tm5_hdf, paramkey, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           status )
        IF_NOTOK_RETURN(status=1)
        ! fill some info values:
        call Init( tmi, paramkey, 'unkown', status )
        call AddHistory( tmi, 'archivekey=='//trim(mf%archivekey), status )
      case ( 'tm5-nc' )
        ! read from hdf file:
        call ReadRecord( mf%mf_tm5_nc, paramkey, unit, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           status )
        IF_NOTOK_RETURN(status=1)
        ! fill some info values:
        call Init( tmi, paramkey, 'unkown', status )
        call AddHistory( tmi, 'archivekey=='//trim(mf%archivekey), status )
#endif
#ifdef with_tmm_ecmwf
      case ( 'ecmwf-tmpp' )
        ! read from grib file:
        call ReadRecord( mf%mf_ecmwf_tmpp, paramkey, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           ggi, gg, sp_gg, &
                           shi, sh, lnsp_sh, &
                           tmi, status )
        IF_NOTOK_RETURN(status=1)
      case ( 'ecmwf-tm5' )
        ! read from grib file:
        call ReadRecord( mf%mf_ecmwf_tm5, paramkey, tday, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           ggi, gg, sp_gg, &
                           shi, sh, lnsp_sh, &
                           tmi, status )
        IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_tmm_ncep
      case ( 'ncep-cdc' )
        ! read from ncep file:
        call ReadRecord( mf%mf_ncep_cdc, paramkey, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           ggi, gg, sp_gg, &
                           shi, sh, lnsp_sh, &
                           tmi, status )
        IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_tmm_ncep
      case ( 'ncep-gfs' )
        ! read from ncep file:
        call ReadRecord( mf%mf_ncep_gfs, paramkey, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           ggi, gg, sp_gg, &
                           shi, sh, lnsp_sh, &
                           tmi, status )
        IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_prism
      case ( 'prism' )
        ! receive from oasis coupler:
        call ReadRecord( mf%mf_prism, paramkey, t1, t2, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           ggi, gg, sp_gg, &
                           shi, sh, lnsp_sh, &
                           tmi, status )
        IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_tmm_msc
      case ( 'msc-data' )
        ! read from grib file:
        tmid = t1 + (t2-t1)/2
        call ReadRecord( mf%mf_msc, paramkey, tmid, tmid, nuv, nw, &
                           gridtype, levi, &
                           lli, ll, sp_ll, &
                           ggi, gg, sp_gg, &
                           shi, sh, lnsp_sh, &
                           tmi, status )
        IF_NOTOK_RETURN(status=1)
#endif
      case default
        write (gol,'("unsupported filetype `",a,"`")') trim(mf%filetype); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

    !write (*,'(a,": end")') name

  end subroutine mf_ReadRecord


  ! ***


  ! ===========================================================
  !
  !   open meteo file for output
  !
  ! ===========================================================


  !
  ! Open the meteo file that should contain the field specified by
  ! archivekey, parameter key, time,
  ! or do nothing if the requested file has been opened already.
  !
  ! <archivekey> = <archivetype>:<archivename>
  !
  !    tmpp:od-fc-ml60-glb3x2
  !

  subroutine mf_SetupOutput( mf, archivekey, paramkey, tday, t1, t2, &
                                 rcfilename, dir, status )

    use GO, only : goSplitLine
    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : TDate

#ifdef with_tmm_tm5
    use tmm_mf_tm5_hdf , only : Init, Get
    use tmm_mf_tm5_nc  , only : Init, Get
#endif

    ! --- in/out -------------------------------------

    type(TMeteoFile), intent(inout)        ::  mf
    character(len=*), intent(in)           ::  archivekey
    character(len=*), intent(in)           ::  paramkey
    type(TDate), intent(in)                ::  tday, t1, t2
    character(len=*), intent(in)           ::  rcfilename
    character(len=*), intent(in)           ::  dir
    integer, intent(inout)                 ::  status

    ! --- const -------------------------------------

    character(len=*), parameter ::  rname = mname//'/mf_SetupOutput'

    ! --- local -------------------------------------

    character(len=10)        ::  archivetype
    character(len=MAX_FILENAME_LEN)       ::  archivename

!    character(len=DUMMY_STR_LEN)       ::  command
!    integer                  ::  year1, month1, day1, hour1
!    integer                  ::  year2, month2, day2, hour2
!    integer                  ::  dth
    type(TrcFile)            ::  infofile

    character(len=MAX_FILENAME_LEN)       ::  archivename2
    character(len=10)        ::  mclass
    character(len=10)        ::  mtype
    character(len=10)        ::  mlevs
    character(len=10)        ::  mgrid
    character(len=10)        ::  filekey
    character(len=16)        ::  treskey
    logical                  ::  with_spm

    ! --- begin -------------------------------------

    ! store archive key:
    mf%archivekey = trim(archivekey)

    ! split archive key in type and name:
    call goSplitLine( archivekey, archivetype, ':', archivename, status )
    IF_NOTOK_RETURN(status=1)

    ! deceide on archive type:
    select case ( archivetype )

#ifdef with_tmm_tm5
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! daily hdf day files
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tm5-hdf' )

        ! always hdf files:
        mf%filetype = 'tm5-hdf'

        ! setup file:
        call Init( mf%mf_tm5_hdf, 'o', dir, trim(archivename), paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        call Get( mf%mf_tm5_hdf, status, filename=mf%filename )
        IF_NOTOK_RETURN(status=1)

        ! extract time range:
        call Get( mf%mf_tm5_hdf, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract paramkeys for fields in file:
        call Get( mf%mf_tm5_hdf, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

      case ( 'tm5-nc' )

        ! always netcdf files:
        mf%filetype = 'tm5-nc'

        ! setup file:
        call Init( mf%mf_tm5_nc, 'o', dir, trim(archivename), paramkey, &
                               tday, t1, t2, status )
        IF_NOTOK_RETURN(status=1)

        ! store filename:
        call Get( mf%mf_tm5_nc, status, filename=mf%filename )
        IF_NOTOK_RETURN(status=1)

        ! extract time range:
        call Get( mf%mf_tm5_nc, status, trange1=mf%t1, trange2=mf%t2 )
        IF_NOTOK_RETURN(status=1)

        ! extract paramkeys for fields in file:
        call Get( mf%mf_tm5_nc, status, paramkeys=mf%paramkeys )
        IF_NOTOK_RETURN(status=1)

#endif

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error ...
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default

        write (gol,'("unsupported archivetype `",a,"`")') trim(archivetype); call goErr
        TRACEBACK; status=1; return

    end select

    ! file is opened (or, at least file name is known)
    mf%opened = .true.

    ! ok
    status = 0

  end subroutine mf_SetupOutput


  ! ***


  subroutine mf_WriteRecord_2d( mf, tmi, paramkey, unit, tday, t1, t2, &
                                lli, nuv, ll, status )

    use GO        , only : TDate
    use Grid      , only : TllGridInfo
    use tmm_info  , only : TMeteoInfo
#ifdef with_tmm_tm5
    use tmm_mf_tm5_hdf, only : WriteRecord
    use tmm_mf_tm5_nc , only : WriteRecord
#endif

    ! --- in/out -------------------------------

    type(TMeteoFile), intent(inout)      ::  mf
    type(TMeteoInfo), intent(in)         ::  tmi
    character(len=*), intent(in)         ::  paramkey, unit
    type(TDate), intent(in)              ::  tday, t1, t2
    type(TllGridInfo), intent(in)        ::  lli
    character(len=1), intent(in)         ::  nuv
    real, intent(in)                     ::  ll(:,:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_WriteRecord_2d'

    ! --- begin ---------------------------------

    select case ( mf%filetype )
#ifdef with_tmm_tm5
      case ( 'hdf', 'tm5-hdf' )
        call WriteRecord( mf%mf_tm5_hdf, tmi, paramkey, unit, tday, t1, t2, &
                            lli, nuv, ll, status )
        IF_NOTOK_RETURN(status=1)
      case ( 'tm5-nc' )
        call WriteRecord( mf%mf_tm5_nc, tmi, paramkey, unit, tday, t1, t2, &
                            lli, nuv, ll, status )
        IF_NOTOK_RETURN(status=1)
#endif
      case default
        write (gol,'("unsupported filetype `",a,"`")') trim(mf%filetype); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine mf_WriteRecord_2d


  ! ***


  subroutine mf_WriteRecord_3d( mf, tmi, spname, paramkey, unit, tday, t1, t2, &
                                lli, nuv, levi, nw, ps, ll, status )!, &
                                !nlev )

    use GO        , only : TDate
    use Grid      , only : TllGridInfo, TLevelInfo
    use tmm_info  , only : TMeteoInfo
#ifdef with_tmm_tm5
    use tmm_mf_tm5_hdf, only : WriteRecord
    use tmm_mf_tm5_nc , only : WriteRecord
#endif

    ! --- in/out -------------------------------

    type(TMeteoFile), intent(inout)      ::  mf
    type(TMeteoInfo), intent(in)         ::  tmi
    character(len=*), intent(in)         ::  spname, paramkey, unit
    type(TDate), intent(in)              ::  tday, t1, t2
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

    ! --- begin ---------------------------------

    select case ( mf%filetype )
#ifdef with_tmm_tm5
      case ( 'tm5-hdf' )
        call WriteRecord( mf%mf_tm5_hdf, tmi, spname, paramkey, unit, tday, t1, t2, &
                            lli, nuv, levi, nw, ps, ll, status )!, &
                            !nlev=nlev )
        IF_NOTOK_RETURN(status=1)
      case ( 'tm5-nc' )
        call WriteRecord( mf%mf_tm5_nc, tmi, spname, paramkey, unit, tday, t1, t2, &
                            lli, nuv, levi, nw, ps, ll, status )!, &
                            !nlev=nlev )
        IF_NOTOK_RETURN(status=1)
#endif
      case default
        write (gol,'("unsupported filetype `",a,"`")') trim(mf%filetype); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine mf_WriteRecord_3d


end module tmm_mf
