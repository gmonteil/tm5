!###############################################################################
!
! Input/output of meteofiles : grib version.
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

module tmm_mf_ecmwf_tmpp

  use GO       , only : gol, goErr, goPr
  use GO       , only : TDate
  use file_grib, only : TGribFile
  use os_specs,  only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

  implicit none

  ! --- in/out ----------------------------

  private

  public  ::  TMeteoFile_ecmwf_tmpp
  public  ::  Init, Done
  public  ::  Get
  public  ::  ReadRecord

  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf_ecmwf_tmpp'


  !--- type ---------------------------------

  type TMeteoFile_ecmwf_tmpp
    ! file name:
    character(len=MAX_FILENAME_LEN)         ::  fname
    ! time range covered by file:
    type(TDate)                ::  trange(2)
    logical                    ::  constant
    ! other time keys for this file:
    character(len=16)          ::  treskey
    type(TDate)                ::  tday
    ! current time range covered by grib record:
    type(TDate)                ::  tref, t1, t2
    !
    ! file description
    !
    character(len=16)          ::  ec_class, ec_type
    character(len=MAX_RCKEY_LEN)         ::  paramkeys
    !
    ! 3d field
    !
    type(TGribFile)            ::  grib
    !
    ! surface pressure field
    !
    logical                    ::  do_spm
    character(len=MAX_FILENAME_LEN)         ::  spm_fname
    type(TGribFile)            ::  spm_grib
    logical                    ::  spm_lnsp2sp
    !
  end type TMeteoFile_ecmwf_tmpp


  ! --- interfaces -------------------------

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


contains


  ! ==============================================================
  ! ====
  ! ==== original
  ! ====
  ! ==============================================================


  subroutine mf_Init( mf, dir, archivekeys, paramkey, &
                               tday, t1, t2, status )

    use GO, only : TDate, Get, NewDate, operator(>)
    use GO, only : goVarValue, goWriteKeyNum

    ! --- in/out ----------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(out)  ::  mf
    character(len=*), intent(in)        ::  dir
    character(len=*), intent(in)        ::  archivekeys
    character(len=*), intent(in)        ::  paramkey
    type(TDate), intent(in)             ::  tday, t1, t2
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Init'

    ! --- local --------------------------------

    character(len=16)     ::  form
    character(len=16)     ::  ec_class, ec_type, ec_tres, ec_levs, mfile
    integer               ::  ec_sh, ec_gg
    character(len=16)     ::  grid, gridN, gridT
    character(len=16)     ::  levs
    character(len=16)     ::  spm_tres

    type(TDate)           ::  tfile
    integer               ::  ccyy, mm, dd
    type(TDate)           ::  tc

    logical               ::  exist

    ! --- begin --------------------------------

    !
    ! extract fields from archivekey :
    !   form=tmpp;class=od;type=fg;nlev=60;sh=159;gg=80;tres=_fg006up4tr3
    !
    form   = 'tmpp'
    call goVarValue( archivekeys, ';', 'form', '=', form, status )
    IF_ERROR_RETURN(status=1)
    !
    ec_class   = 'od'
    call goVarValue( archivekeys, ';', 'class', '=', ec_class, status )
    IF_ERROR_RETURN(status=1)
    !
    ec_type    = 'fc'
    call goVarValue( archivekeys, ';', 'type', '=', ec_type, status )
    IF_ERROR_RETURN(status=1)
    !
    ec_levs    = 'ml60'
    call goVarValue( archivekeys, ';', 'levs', '=', ec_levs, status )
    IF_ERROR_RETURN(status=1)
    !
    ec_sh      = 159
    call goVarValue( archivekeys, ';', 'sh', '=', ec_sh, status )
    IF_ERROR_RETURN(status=1)
    !
    ec_gg      = 80
    call goVarValue( archivekeys, ';', 'gg', '=', ec_gg, status )
    IF_ERROR_RETURN(status=1)
    !
    ec_tres    = ''
    call goVarValue( archivekeys, ';', 'tres', '=', ec_tres, status )
    IF_ERROR_RETURN(status=1)

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! fg data available up to 2000-09-12 00:00 ; after, use fc data:
    if ( ec_tres == '_fg006up4tr3' ) then
      tc = NewDate( 2000, 09, 12, 00, 00, 00 )
      if ( t1 > tc ) then
        write (gol,'("      WARNING - using fc data after 2000-09-12 00:00")'); call goPr
        ec_type = 'fc'
        ec_tres = '_fc012up2tr3'
      end if
    end if
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! store some values ...
    mf%ec_class = ec_class
    mf%ec_type  = ec_type

    !
    ! meteo file
    !

    ! grid : T159, N80, etc
    call goWriteKeyNum( gridT, 'T', ec_sh )
    call goWriteKeyNum( gridN, 'N', ec_gg )

    ! defaults
    levs           = ec_levs
    grid           = 'Xxx'
    mfile          = trim(paramkey)
    mf%paramkeys   = '-'//trim(paramkey)//'-'
    mf%constant    = .false.
    mf%do_spm      = .false.
    mf%spm_fname   = 'none'
    mf%spm_lnsp2sp = .false.
    spm_tres       = ec_tres

    ! specials
    select case ( paramkey )
      case ( 'LNSP' )
        levs           = 'ml1'
        grid           = gridT
        if ( form /= '3hr' ) ec_tres = '_03p03'
      case ( 'VO', 'D' )
        mfile = 'VOD'
        mf%paramkeys   = '-VO-D-'
        grid           = gridT
        mf%do_spm      = .true.
        spm_tres       = '_03p03'
        if ( form == '3hr' ) then
          mfile        = paramkey
          mf%paramkeys = '-'//trim(paramkey)//'-'
          spm_tres     = ec_tres
        end if
      case ( 'T', 'W' )
        grid           = gridT
        mf%do_spm      = .true.
        spm_tres       = '_03p03'
        if ( form == '3hr' ) then
          spm_tres     = ec_tres
        end if
      case ( 'Q' )
        grid           = gridN
        mf%do_spm      = .true.
        spm_tres       = '_03p03'
        mf%spm_lnsp2sp = .true.
        if ( form == '3hr' ) then
          spm_tres     = ec_tres
        end if
      case ( 'CLWC', 'CIWC', 'CC', 'CCO', 'CCU' )
        grid           = gridN
        mfile          = 'CLD'
        mf%paramkeys   = '-CLWC-CIWC-CC-CCO-CCU-'
        mf%do_spm      = .true.
        spm_tres       = '_03p03'
        mf%spm_lnsp2sp = .true.
        if ( form == '3hr' ) then
          spm_tres     = ec_tres
        end if
      case ( 'oro' )
        levs           = 'sfc'
        grid           = gridN
        mf%constant    = .true.
        ec_tres        = 'const'
      case ( 'lsm' )
        levs           = 'sfc'
        grid           = gridN
        mf%constant    = .true.
      case ( 'cvl', 'cvh', 'tvl', 'tvh', 'sr', 'al', 'lsrh' )
        ec_type        = 'an'
        levs           = 'sfc'
        grid           = gridN
        mfile          = 'S0'
        mf%paramkeys   = '-cvl-cvh-tvl-tvh-sr-al-lsrh-'
        ec_tres        = '_06p06'
      case ( 'ci', 'swvl1', 'swvl2', 'swvl3', 'swvl4', '10fg', 'sd', 'lsp', &
             'cp', 'sf', 'sshf', 'slhf', 'blh', 'u10m', 'v10m', 't2m', 'd2m', &
             'ssr', 'ewss', 'nsss', 'sstr' ,'src' )
        levs           = 'sfc'
        grid           = gridN
        mfile          = 'S1'
        mf%paramkeys   = '-'
        mf%paramkeys   = trim(mf%paramkeys)//'ci-swvl1-swvl2-swvl3-swvl4-10fg-sd-lsp-'
        mf%paramkeys   = trim(mf%paramkeys)//'cp-sf-sshf-slhf-blh-10u-10v-2t-2d-'
        mf%paramkeys   = trim(mf%paramkeys)//'ssr-ewss-nsss-sstr-src-'
        ec_tres        = '_03p03'
    end select

    ! convert input times to file name times:
    call GetGribTime( ec_tres, tday, t1, t2, status, tfile=tfile, trange=mf%trange )
    IF_NOTOK_RETURN(status=1)

    ! store some time params
    mf%treskey = ec_tres
    if ( ec_type == 'an' ) mf%treskey = 'an'//ec_tres
    mf%tday    = tday

    ! create file name:
    !   dir/od-fc-2000-01-ml60-T159-T_20000101_fg006up4tr3.gb
    select case ( paramkey )
      case ( 'oro' )
        write (mf%fname,'(a,"/",a,"-",a,"-",i4.4,"-",i2.2,"-",a,"-",a,"-",a,a)') &
                  trim(dir), &
                  trim(ec_class), trim(ec_type), 0000, 00, trim(levs), trim(grid), &
                  'Z', '.gb'
      case ( 'lsm' )
        write (mf%fname,'(a,"/",a,"-",a,"-",i4.4,"-",i2.2,"-",a,"-",a,"-",a,a)') &
                  trim(dir), &
                  trim(ec_class), trim(ec_type), 0000, 00, trim(levs), trim(grid), &
                  'LSM', '.gb'
      case default
        ! extract time values:
        call Get( tfile, year=ccyy, month=mm, day=dd )
        ! path and file include date:
        write (mf%fname,'(a,"/",a,"-",a,"-",i4.4,"-",i2.2,"-",a,"-",a,"-",a,"_",i4.4,2i2.2,a,a)') &
                  trim(dir), &
                  trim(ec_class), trim(ec_type), ccyy, mm, trim(levs), trim(grid), &
                  trim(mfile), ccyy, mm, dd, trim(ec_tres), '.gb'
    end select

    ! exist ?
    inquire( file=mf%fname, exist=exist )
    if ( .not. exist ) then
      write (gol,'("grib file does not exist:")'); call goErr
      write (gol,'("  ",a)') trim(mf%fname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if


    !
    ! sp or lnsp ?
    !

    if ( mf%do_spm ) then

      ! grid : T159
      grid = gridT

      ! levels: ml1
      levs = 'ml1'

      ! create file name:
      !   dir/od-fc-2000-01-ml60-T159-T_20000101_fg006up4tr3.gb
      write (mf%spm_fname,'(a,"/",a,"-",a,"-",i4.4,"-",i2.2,"-",a,"-",a,"-",a,"_",i4.4,2i2.2,a,a)') &
                trim(dir), &
                trim(ec_class), trim(ec_type), ccyy, mm, trim(levs), trim(grid), &
                'LNSP', ccyy, mm, dd, trim(spm_tres), '.gb'

      ! exist ?
      inquire( file=mf%spm_fname, exist=exist )
      if ( .not. exist ) then
        write (gol,'("grib file does not exist:")'); call goErr
        write (gol,'("  ",a)') trim(mf%spm_fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

    end if

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    use file_grib, only : Done

    ! --- in/out ------------------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(inout)  ::  mf
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Done'

    ! --- begin -------------------------------------

    ! nothing to be done ...

    ! ok
    status = 0

  end subroutine mf_Done


  ! ***


  subroutine mf_Get( mf, status, trange1, trange2, paramkeys )

    use GO, only : TDate

    ! --- in/out ----------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(in)   ::  mf
    integer, intent(out)                ::  status

    type(TDate), intent(out), optional        ::  trange1, trange2
    character(len=*), intent(out), optional   ::  paramkeys

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Get'

    ! --- local --------------------------------

    ! --- begin --------------------------------

    ! time range:
    if ( present(trange1) ) trange1 = mf%trange(1)
    if ( present(trange2) ) trange2 = mf%trange(2)

    ! contents:
    if ( present(paramkeys) ) paramkeys = trim(mf%paramkeys)

    ! ok
    status = 0

  end subroutine mf_Get


  ! ***


  subroutine mf_ReadRecord( mf, paramkey, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use parray    , only : pa_Init, pa_Done
    use GO        , only : TDate
    use Grid      , only : TLevelInfo
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use tmm_info  , only : TMeteoInfo

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(inout)  ::  mf
    character(len=*), intent(in)          ::  paramkey
    type(TDate), intent(in)               ::  t1, t2
    character(len=1), intent(in)          ::  nuv
    character(len=1), intent(in)          ::  nw
    character(len=2), intent(out)         ::  gridtype
    type(TLevelInfo), intent(out)         ::  levi
    type(TllGridInfo), intent(inout)      ::  lli
    real, pointer                         ::  ll(:,:,:)
    real, pointer                         ::  sp_ll(:,:)
    type(TggGridInfo), intent(inout)      ::  ggi
    real, pointer                         ::  gg(:,:)
    real, pointer                         ::  sp_gg(:)
    type(TshGridInfo), intent(inout)      ::  shi
    complex, pointer                      ::  sh(:,:)
    complex, pointer                      ::  lnsp_sh(:)
    type(TMeteoInfo), intent(out)         ::  tmi
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord'

    ! --- local -------------------------------

    real, pointer         ::  ll2(:,:,:)
    real, pointer         ::  gg2(:,:)
    complex, pointer      ::  sh2(:,:)

    ! --- begin ---------------------------------

    ! combined field ?
    select case ( paramkey )

      case ( 'sstr' )

        ! read first field:
        call mf_ReadRecord_1( mf, 'ewss', t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! init pointer:
        call pa_Init( ll2 )
        call pa_Init( gg2 )
        call pa_Init( sh2 )

        ! read second field:
        call mf_ReadRecord_1( mf, 'nsss', t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll2, sp_ll, &
                                ggi, gg2, sp_gg, &
                                shi, sh2, lnsp_sh, &
                                tmi, status )
         IF_NOTOK_RETURN(status=1)

        ! process:
        select case ( gridtype )
          case ( 'll' ) ; ll = sqrt( ll**2 + ll2**2 )
          case ( 'gg' ) ; gg = sqrt( gg**2 + gg2**2 )
          case ( 'sh' ) ; sh = sqrt( sh**2 + sh2**2 )
          case default
            write (gol,'("unsupported gridtype for substract :",a)') gridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! clear pointers:
        call pa_Done( ll2 )
        call pa_Done( gg2 )
        call pa_Done( sh2 )

      case default

        call mf_ReadRecord_1( mf, paramkey, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )
         IF_NOTOK_RETURN(status=1)

    end select

    ! ok
    status = 0

  end subroutine mf_ReadRecord



  ! ***


  subroutine mf_ReadRecord_1( mf, paramkey, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use parray    , only : pa_Init, pa_Done
    use GO        , only : TDate, operator(<), operator(-), rTotal, wrtgol
    use Grid      , only : TLevelInfo
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use tmm_info  , only : TMeteoInfo

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(inout)  ::  mf
    character(len=*), intent(in)          ::  paramkey
    type(TDate), intent(in)               ::  t1, t2
    character(len=1), intent(in)          ::  nuv
    character(len=1), intent(in)          ::  nw
    character(len=2), intent(out)         ::  gridtype
    type(TLevelInfo), intent(out)         ::  levi
    type(TllGridInfo), intent(inout)      ::  lli
    real, pointer                         ::  ll(:,:,:)
    real, pointer                         ::  sp_ll(:,:)
    type(TggGridInfo), intent(inout)      ::  ggi
    real, pointer                         ::  gg(:,:)
    real, pointer                         ::  sp_gg(:)
    type(TshGridInfo), intent(inout)      ::  shi
    complex, pointer                      ::  sh(:,:)
    complex, pointer                      ::  lnsp_sh(:)
    type(TMeteoInfo), intent(out)         ::  tmi
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord_1'

    ! --- local -------------------------------

    type(TDate)           ::  tref
    real, pointer         ::  ll1(:,:,:)
    real, pointer         ::  gg1(:,:)
    complex, pointer      ::  sh1(:,:)
    real                  ::  dt_sec

    ! --- begin ---------------------------------

    ! accumulated field ?
    select case ( paramkey )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! accumulated fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'lsp', 'cp', 'sf', 'sshf', 'slhf', 'ssr', 'ewss', 'nsss' )

        ! get reference time for requested time interval:
        call GetGribTime( mf%treskey, mf%tday, t1, t2, status, tref=tref )
        IF_NOTOK_RETURN(status=1)

        ! should be a time interval ...
        if ( .not. (t1 < t2) ) then
          write (gol,'("accumulated fields requires time interval:")'); call goErr
          write (gol,'("  paramkey : ",a)') paramkey; call goErr
          call wrtgol( '  t1       : ', t1 ); call goErr
          call wrtgol( '  t2       : ', t2 ); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if

        ! read field accumulated over [tref,t2] :
        call mf_ReadRecord_2( mf, paramkey, tref, t2, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! substract [tref,t1] if necessary:
        if ( tref < t1 ) then

          ! init pointer:
          call pa_Init( ll1 )
          call pa_Init( gg1 )
          call pa_Init( sh1 )

          ! read field accumulated over [tref,t1] :
          call mf_ReadRecord_2( mf, paramkey, tref, t1, t1, nuv, nw, &
                                  gridtype, levi, &
                                  lli, ll1, sp_ll, &
                                  ggi, gg1, sp_gg, &
                                  shi, sh1, lnsp_sh, &
                                  tmi, status )
           IF_NOTOK_RETURN(status=1)

          ! substract:
          select case ( gridtype )
            case ( 'll' ) ; ll = ll - ll1
            case ( 'gg' ) ; gg = gg - gg1
            case ( 'sh' ) ; sh = sh - sh1
            case default
              write (gol,'("unsupported gridtype for substract :",a)') gridtype; call goErr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select

          ! clear pointers:
          call pa_Done( ll1 )
          call pa_Done( gg1 )
          call pa_Done( sh1 )

        end if

        ! return time averages only:
        dt_sec = rTotal( t2 - t1, 'sec' )
        select case ( gridtype )
          case ( 'll' ) ; ll = ll / dt_sec
          case ( 'gg' ) ; gg = gg / dt_sec
          case ( 'sh' ) ; sh = sh / dt_sec
          case default
            write (gol,'("unsupported gridtype for time average :",a)') gridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! instantaneous fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        ! get reference time for requested time interval:
        call GetGribTime( mf%treskey, mf%tday, t1, t2, status, tref=tref )
        IF_NOTOK_RETURN(status=1)

        ! just read ..
        call mf_ReadRecord_2( mf, paramkey, tref, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )
         IF_NOTOK_RETURN(status=1)

    end select

    ! ok
    status = 0

  end subroutine mf_ReadRecord_1


  ! ***


  subroutine mf_ReadRecord_2( mf, paramkey, tref, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use GO        , only : TDate, wrtgol
    use Grid      , only : Init, Done
    use Grid      , only : TLevelInfo
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use Grid      , only : Interpol
    use file_grib , only : Init, Done, ReadRecord, Get, Check
    use file_grib , only : levtype_sfc, levtype_hyb
    use file_grib , only : gridtype_ll, gridtype_gg, gridtype_sh
    use grib_table, only : GetPid
    use PArray    , only : pa_Init, pa_Done, pa_SetShape
    use tmm_info  , only : TMeteoInfo, Init, AddHistory

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(inout)  ::  mf
    character(len=*), intent(in)          ::  paramkey
    type(TDate), intent(in)               ::  tref, t1, t2
    character(len=1), intent(in)          ::  nuv
    character(len=1), intent(in)          ::  nw
    character(len=2), intent(out)         ::  gridtype
    type(TLevelInfo), intent(out)         ::  levi
    type(TllGridInfo), intent(inout)      ::  lli
    real, pointer                         ::  ll(:,:,:)
    real, pointer                         ::  sp_ll(:,:)
    type(TggGridInfo), intent(inout)      ::  ggi
    real, pointer                         ::  gg(:,:)
    real, pointer                         ::  sp_gg(:)
    type(TshGridInfo), intent(inout)      ::  shi
    complex, pointer                      ::  sh(:,:)
    complex, pointer                      ::  lnsp_sh(:)
    type(TMeteoInfo), intent(out)         ::  tmi
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord_2'

    ! --- local -------------------------------

    logical                ::  isfirst
    logical                ::  reopened
    integer                ::  pid
    integer                ::  nlev, glevtype, glevel

    integer                ::  level

    integer                ::  ggridtype
    real                   ::  lon_first, lon_inc
    integer                ::  lon_n
    real                   ::  lat_first, lat_inc
    integer                ::  lat_n
    integer                ::  ggN
    integer                ::  shT

    integer                ::  greftime(5), gtimerange(4)
    character(len=64)      ::  key

    integer                ::  ilat
    real, pointer          ::  pat(:,:)

    type(TshGridInfo)      ::  tmp_shi
    complex, pointer       ::  tmp_sh(:)

    ! --- begin ---------------------------------

    ! no fluxes through boundaries ...
    if ( nuv /= 'n' ) then
      write (gol,'("unsupported nuv key : ",a)') nuv; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! init pointer arrays:
    call pa_Init( pat )

    !
    ! ~~~ 3d field or 2d stored in 3d array
    !

    ! open grib file
    call Init( mf%grib, mf%fname, 'r', status )
    IF_NOTOK_RETURN(status=1)

    ! arrays and grids not defined yet
    isfirst = .true.
    reopened = .false.

    ! loop over records
    level = 0
    do

      !
      ! read gribsection in file buffer
      !

      call ReadRecord( mf%grib, status )
      select case ( status )
        case ( 0 )
          ! no error
        case ( 1 )
          ! eof
          if ( .not. reopened ) then
            !write (*,'("grib read record: re-open ...")')
            ! close:
            call Done( mf%grib, status )
            IF_NOTOK_RETURN(status=1)
            ! reopen:
            call Init( mf%grib, mf%fname, 'r', status )
            IF_NOTOK_RETURN(status=1)
            reopened = .true.
            cycle
          else
            write (gol,'("reached eof before requested record was found")'); call goErr
            write (gol,'("  file     : ",a)') trim(mf%fname); call goErr
            write (gol,'("  paramkey : ",a)') trim(paramkey); call goErr
            call wrtgol( '  tref     : ', tref ); call goErr
            call wrtgol( '  t1       : ', t1 ); call goErr
            call wrtgol( '  t2       : ', t2 ); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
          end if
        case default
          write (gol,'("error from grib ReadRecord; status=",i6)') status; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select


      !
      ! checks ...
      !

      ! get param id for the requested field from grib table:
      select case ( paramkey )
        case ( 'spm' ) ; call GetPid( 'ec', 'SP'    , pid, status )
        case default   ; call GetPid( 'ec', paramkey, pid, status )
      end select
      IF_NOTOK_RETURN(status=1)

      ! check parameter; continue if not ok:
      call Check( mf%grib, pid=pid, debug=0, status=status )
      if (status/=0) cycle

      ! fill times ?
      if ( .not. mf%constant ) then

        ! extract time fields from grib, store in mf%tref/mf%t1/mf%t2
        call SetTime( mf, status )
        IF_NOTOK_RETURN(status=1)

        ! check time:
        call CheckTime( mf, tref, t1, t2, status )
        if (status/=0) then
          !write (*,'("grib read record: wrong time; skip ...")')
          cycle
          !write (gol,'("found unexpected times in grib file:")'); call goErr
          !write (gol,'("  paramkey     : ",a)') paramkey; call goErr
          !call wrtgol( '  req. t1 : ', t1 ); call goErr
          !call wrtgol( '       t2 : ', t2 ); call goErr
          !call wrtgol( '  grib t1 : ', mf%t1 ); call goErr
          !call wrtgol( '       t2 : ', mf%t2 ); call goErr
          !write (gol,'("  grib file    : ",a)') trim(mf%fname); call goErr
          !write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if

      end if  ! time checking

      ! extract level stuff:
      call Get( mf%grib, nlev=nlev, levtype=glevtype, level=glevel, status=status )
      IF_NOTOK_RETURN(status=1)
      ! check level type:
      select case ( glevtype )
        case ( levtype_sfc )
          ! surface field
          nlev = 1
          glevel = 1
        case ( levtype_hyb )
          select case ( paramkey )
            case ( 'LNSP' )
              nlev = 1
            case default
              ! level in 3d field
          end select
        case default
          write (gol,'("found unexpected level type: ")'); call goErr
          write (gol,'("  leveltype : ",i3)') glevtype; call goErr
          write (gol,'("  paramkey  : ",a)') paramkey; call goErr
          write (gol,'("  grib file : ",a)') trim(mf%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select
      ! check level:
      level = level + 1
      if ( glevel /= level ) then
        write (gol,'("found unexpected level: ")'); call goErr
        write (gol,'("  paramkey   : ",a)') paramkey; call goErr
        write (gol,'("  req. level : ",i6)') level; call goErr
        write (gol,'("  grib level : ",i6)') glevel; call goErr
        write (gol,'("  grib nlev  : ",i6)') nlev; call goErr
        write (gol,'("  grib file  : ",a)') trim(mf%fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

      !
      ! define grids and arrays
      !
      if ( isfirst ) then
        !
        ! * info
        !
        ! example of history:
        !   model=ecmwf;class=od;type=fc;tref=2000,12,31,12,00; ...
        !     trange=001,234,240,001;sh=159;nlev=60
        !
        call Init( tmi, paramkey, 'unknown', status )
        call AddHistory( tmi, 'model==ecmwf', status )
        call AddHistory( tmi, 'class=='//trim(mf%ec_class), status )
        call AddHistory( tmi, 'type=='//trim(mf%ec_type)  , status )
        !
        call Get( mf%grib, status, reftime=greftime, timerange=gtimerange )
        IF_NOTOK_RETURN(status=1)
        write (key,'("tref==",i4.4,4(",",i2.2))') greftime
        call AddHistory( tmi, trim(key), status )
        write (key,'("trange==",i3.3,3(",",i3.3))') gtimerange
        call AddHistory( tmi, trim(key), status )
        !
        write (key,'("nlev==",i3.3)') nlev
        call AddHistory( tmi, trim(key), status )

        !
        ! * define horizontal grid:
        !

        ! extract grid type:
        call Get( mf%grib, status, gridtype=ggridtype )
        IF_NOTOK_RETURN(status=1)

        ! setup:
        select case ( ggridtype )

          ! o lat/lon
          case ( gridtype_ll )

            ! routine returns lat/lon grid:
            gridtype = 'll'

            ! grib storage is north pole to south pole:
            call Get( mf%grib, status, &
                               lon_first=lon_first, lon_inc=lon_inc, lon_n=lon_n, &
                               lat_last =lat_first, lat_inc=lat_inc, lat_n=lat_n     )
            IF_NOTOK_RETURN(status=1)

            ! define grid structure:
            call Init( lli, lon_first, lon_inc, lon_n, &
                            lat_first, lat_inc, lat_n, status     )
            IF_NOTOK_RETURN(status=1)

            ! init array to store 2d field from grib file (north-south order):
            call pa_SetShape( pat, lon_n, lat_n )

            ! allocate output:
            call pa_SetShape( ll, lon_n, lat_n, nlev )

            ! add to history:
            write (key,'("longrid==",f7.2,",",f6.2,",",i4)') lon_first, lon_inc, lon_n
            call AddHistory( tmi, trim(key), status )
            write (key,'("latgrid==",f7.2,",",f6.2,",",i4)') lat_first, lat_inc, lat_n
            call AddHistory( tmi, trim(key), status )

          ! o gaussian grid
          case ( gridtype_gg )

            ! routine returns gg grid:
            gridtype = 'gg'

            ! extract grid number:
            call Get( mf%grib, status, N=ggN )
            IF_NOTOK_RETURN(status=1)

            ! define grid structure:
            call Init( ggi, ggN, .true., status )
            IF_NOTOK_RETURN(status=1)

            ! allocate output:
            call pa_SetShape( gg, ggi%np, nlev )

            ! add to history:
            write (key,'("gg==",i4.4)') ggN
            call AddHistory( tmi, trim(key), status )

          ! o spectral field:
          case ( gridtype_sh )

            ! routine returns sh grid:
            gridtype = 'sh'

            ! extract spectral truncation:
            call Get( mf%grib, status, T=shT )
            IF_NOTOK_RETURN(status=1)

            ! intialize spherical harmonic field info:
            call Init( shi, shT, status )
            IF_NOTOK_RETURN(status=1)

            ! allocate output:
            call pa_SetShape( sh, shi%np, nlev )

            ! add to history:
            write (key,'("sh==",i4.4)') shT
            call AddHistory( tmi, trim(key), status )

          case default
            write (gol,'("unsupported gridtype for setup : ",i6)') ggridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        !
        ! * levels
        !
        select case ( nlev )
          case ( 1 )
            call Init( levi, nlev, (/0.0,0.0/), (/0.0,0.0/), status )
            IF_NOTOK_RETURN(status=1)
          case ( 60 )
            call Init( levi, 'ec60', status )
            IF_NOTOK_RETURN(status=1)
        end select

        ! not again ...
        isfirst = .false.
      end if  ! isfirst (grid definition and allocation)

      !
      ! store
      !

      select case ( ggridtype )

        case ( gridtype_ll )

          ! read 2d pat from grib; storred from north to south
          call Get( mf%grib, status, ll=pat )
          IF_NOTOK_RETURN(status=1)

          ! store from south to north:
          do ilat = 1, lat_n
            ll(:,ilat,level) = pat(:,lat_n+1-ilat)
          end do

        case ( gridtype_gg )

          ! read 2d pat from grib:
          call Get( mf%grib, status, gg=gg(:,level) )
          IF_NOTOK_RETURN(status=1)

        case ( gridtype_sh )

          ! read 2d pat from grib:
          call Get( mf%grib, status, sh=sh(:,level) )
          IF_NOTOK_RETURN(status=1)

        case default
          write (gol,'("unsupported gridtype for 2d pat : ",i6)') gridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! last record for this field ?
      if ( glevel == nlev ) exit

    end do   ! records

    ! close grib file
    call Done( mf%grib, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ~~~ surface pressure
    !

    if ( mf%do_spm ) then

      ! open grib file
      call Init( mf%spm_grib, mf%spm_fname, 'r', status )
      IF_NOTOK_RETURN(status=1)

      ! loop over time records
      do

        ! read gribsection in file buffer
        call ReadRecord( mf%spm_grib, status )
        IF_NOTOK_RETURN(status=1)

        ! fill times
        call SetTime( mf, status )
        IF_NOTOK_RETURN(status=1)

        ! check time:
        call CheckTime( mf, tref, t1, t2, status )
        if (status/=0) then
          !write (*,'("grib read record: spm wrong time; skip ...")')
          cycle
          !write (gol,'("found unexpected times in grib file:")'); call goErr
          !write (gol,'("  paramkey     : ",a)') paramkey; call goErr
          !call wrtgol( '  req. t1 : ', t1 ); call goErr
          !call wrtgol( '       t2 : ', t2 ); call goErr
          !call wrtgol( '  grib t1 : ', mf%t1 ); call goErr
          !call wrtgol( '       t2 : ', mf%t2 ); call goErr
          !write (gol,'("  grib file    : ",a)') trim(mf%spm_fname); call goErr
          !write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if

        ! time ok
        exit

      end do  ! time loop

      ! set param id:
      select case ( ggridtype )
        case ( gridtype_ll )
          call GetPid( 'ec', 'SP', pid, status )
          IF_NOTOK_RETURN(status=1)
        case ( gridtype_gg )
          if ( mf%spm_lnsp2sp ) then
            call GetPid( 'ec', 'LNSP', pid, status )
            IF_NOTOK_RETURN(status=1)
          else
            call GetPid( 'ec', 'SP', pid, status )
            IF_NOTOK_RETURN(status=1)
          end if
        case ( gridtype_sh )
          call GetPid( 'ec', 'LNSP', pid, status )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("unsupported gridtype for setup sp/lnsp : ",i6)') ggridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select
      ! check parameter:
      call Check( mf%spm_grib, pid=pid, debug=1, status=status )
      IF_NOTOK_RETURN(status=1)

      ! check level:
      call Get( mf%spm_grib, levtype=glevtype, level=glevel, status=status )
      IF_NOTOK_RETURN(status=1)
      select case ( ggridtype )
        case ( gridtype_ll )
          if ( glevtype /= levtype_sfc ) then
            write (gol,'("found unexpected level type ")'); call goErr
            write (gol,'("  paramkey        : ",a)') paramkey; call goErr
            write (gol,'("  sfc level type  : ",i6)') levtype_sfc; call goErr
            write (gol,'("  grib level type : ",i6)') glevtype; call goErr
            write (gol,'("  grib file       : ",a)') trim(mf%spm_fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
          end if
        case ( gridtype_gg )
          if ( mf%spm_lnsp2sp ) then
            if ( (glevtype /= levtype_hyb) .or. (glevel /= 1) ) then
              write (gol,'("found unexpected level type (lnsp for 3d gg)")'); call goErr
              write (gol,'("  paramkey        : ",a)') paramkey; call goErr
              write (gol,'("  hyb level type  : ",i6)') levtype_hyb; call goErr
              write (gol,'("  grib level type : ",i6)') glevtype; call goErr
              write (gol,'("  grib level      : ",i6)') glevel; call goErr
              write (gol,'("  grib file       : ",a)') trim(mf%spm_fname); call goErr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
            end if
          else
            if ( glevtype /= levtype_sfc ) then
              write (gol,'("found unexpected level type ")'); call goErr
              write (gol,'("  paramkey        : ",a)') paramkey; call goErr
              write (gol,'("  sfc level type  : ",i6)') levtype_sfc; call goErr
              write (gol,'("  grib level type : ",i6)') glevtype; call goErr
              write (gol,'("  grib file       : ",a)') trim(mf%spm_fname); call goErr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
            end if
          end if
        case ( gridtype_sh )
          if ( (glevtype /= levtype_hyb) .or. (glevel /= 1) ) then
            write (gol,'("found unexpected level type ")'); call goErr
            write (gol,'("  paramkey        : ",a)') paramkey; call goErr
            write (gol,'("  hyb level type  : ",i6)') levtype_hyb; call goErr
            write (gol,'("  grib level type : ",i6)') glevtype; call goErr
            write (gol,'("  grib level      : ",i6)') glevel; call goErr
            write (gol,'("  grib file       : ",a)') trim(mf%spm_fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
          end if
        case default
          write (gol,'("unsupported gridtype for sp/lnsp levs : ",i6)') ggridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! read and store surface pressure field:
      select case ( ggridtype )

        case ( gridtype_ll )
          ! allocate storage
          call pa_SetShape( sp_ll, lon_n, lat_n )
          ! read 2d pat from grib; storred from north to south
          call Get( mf%spm_grib, status, ll=pat )
          IF_NOTOK_RETURN(status=1)
          ! store from south to north:
          do ilat = 1, lat_n
            sp_ll(:,ilat) = pat(:,lat_n+1-ilat)
          end do

        case ( gridtype_gg )
          ! allocate storage
          call pa_SetShape( sp_gg, ggi%np )
          ! convert from sh lnsp or read directly:
          if ( mf%spm_lnsp2sp ) then
            ! allocate output:
            call pa_SetShape( sh, shi%np, nlev )
            ! extract spectral truncation:
            call Get( mf%spm_grib, status, T=shT )
            IF_NOTOK_RETURN(status=1)
            ! intialize spherical harmonic field info:
            call Init( tmp_shi, shT, status )
            IF_NOTOK_RETURN(status=1)
            ! allocate storage for lnsp:
            call pa_Init( tmp_sh )
            call pa_SetShape( tmp_sh, tmp_shi%np )
            ! read lnsp on sh grid:
            call Get( mf%spm_grib, status, sh=tmp_sh )
            IF_NOTOK_RETURN(status=1)
            ! interpolate from sh to gg:
            call Interpol( tmp_shi, tmp_sh, ggi, sp_gg, status )  ! ln(Pa)
            IF_NOTOK_RETURN(status=1)
            ! clear
            call Done( tmp_shi )
            call pa_Done( tmp_sh )
            ! convert from lnsp to sp :
            sp_gg = exp(sp_gg)   ! Pa
          else
            ! read 2d pat from grib; storred from north to south
            call Get( mf%spm_grib, status, gg=sp_gg )
            IF_NOTOK_RETURN(status=1)
          end if

        case ( gridtype_sh )
          ! allocate storage
          call pa_SetShape( lnsp_sh, shi%np )
          ! read spectral coeff:
          call Get( mf%spm_grib, status, sh=lnsp_sh )
          IF_NOTOK_RETURN(status=1)

        case default
          write (gol,'("unsupported gridtype for reading sp/lnsp : ",i6)') ggridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! close grib file
      call Done( mf%spm_grib, status )
      IF_NOTOK_RETURN(status=1)

    end if

    !
    ! ~~~ end
    !

    ! deallocate arrays
    call pa_Done( pat )

    ! ok
    status = 0

  end subroutine mf_ReadRecord_2


  ! ****************************************************************************


  !
  ! In gribfile:
  !    reftime : for example time at which forecast is made
  !    timerange: increment or interval
  !
  ! arguments:          ok if:
  !
  !   time1 == time2      time1==time2 == reftime+timerange
  !
  !   time1 == 0          time2 == reftime+timerange
  !
  !   time2 == 0          time1 == reftime
  !
  !   time1 < time2       time1 == refitme, time2 == reftime+timerange
  !
  !
  !   grib [t1-----------t2]
  !                                                              status
  !    o                         time1/2        record too old      1
  !    o                time1/2                 ok                  0
  !    o       time1/2                          record too new      2
  !
  !
  ! SetTime( mf, status )
  !
  !   Extracts time values from current grib record,
  !   store in mf%t1, mf%t2
  !
  !   return status:
  !     0       : ok
  !     other   : some error
  !
  ! CheckTime( mf, time1, time2, status )
  !
  !   return status:
  !     0   : times match
  !     1   : times do not match, try next record
  !     2   : current record is newer than requested (reopen ?)
  !     3   : some error
  !


  ! ***


  !
  ! Return time parameters for grib files:
  !  o tfile   :  date in filename
  !  o trange  :  time interval covered by fields in file
  !  o tref    :  reference time (forecast start?) for tday,[t1,t2]
  !

  subroutine GetGribTime( ec_treskey, tday, t1, t2, status, tfile, trange, tref )

    use GO, only : TDate, Get, Set, wrtgol, NewDate, IncrDate, operator(+), operator(-)

    ! --- in/out --------------------------------

    character(len=*), intent(in)            ::  ec_treskey
    type(TDate), intent(in)                 ::  tday, t1, t2
    integer, intent(out)                    ::  status

    type(TDate), intent(out), optional      ::  tfile
    type(TDate), intent(out), optional      ::  trange(2)
    type(TDate), intent(out), optional      ::  tref

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GetGribTime'

    ! --- local --------------------------------

    integer          ::  hour2, time6(6)
    integer          ::  dd, hh, step

    ! --- begin --------------------------------

    ! set day shift, start hour, and step
    select case ( ec_treskey )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! constant field
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'const' )

        ! only tday is usefull ...
        if ( present(tfile ) ) tfile  = tday
        if ( present(trange) ) trange = (/t1,t2/)  ! any, any
        if ( present(tref  ) ) tref   = tday       ! dummy ...

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! fg, 3 hourly
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_fg006up4tr3' )

        ! extract end hour:
        call Get( t2, hour=hour2 )

        ! set forecast start time and step given hour:
        select case ( hour2 )
          case ( 00 ) ; dd =  0  ;  hh = 00  ;  step = 6
          case ( 03 ) ; dd =  0  ;  hh = 03  ;  step = 3
          case ( 06 ) ; dd =  0  ;  hh = 06  ;  step = 6
          case ( 09 ) ; dd =  0  ;  hh = 09  ;  step = 3
          case ( 12 ) ; dd =  0  ;  hh = 12  ;  step = 6
          case ( 15 ) ; dd =  0  ;  hh = 15  ;  step = 3
          case ( 18 ) ; dd =  0  ;  hh = 18  ;  step = 6
          case ( 21 ) ; dd =  0  ;  hh = 21  ;  step = 3
          case default
            write (gol,'("unsupported hour :")'); call goErr
            write (gol,'("  hour2     : ",i2)') hour2; call goErr
            write (gol,'("  timesteps : ",a )') ec_treskey; call goErr
            call wrtgol( '  time1     : ', t1 ); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! file ccyymmdd contains fields for [00,24) :
        if ( present(tfile) ) then
          tfile = t1
          call Set( tfile, hour=00 )
        end if

        ! fields valid for [00,24) :
        if ( present(trange) ) then
          trange(1) = t1
          call Set( trange(1), hour=0, min=0, sec=0, mili=0 )
          trange(2) = trange(1) + IncrDate(hour=23,min=59,sec=59,mili=999)
          ! trap last fg time 2000-09-12 00:00 ; data valid forr 00:00 only:
          call Get( tfile, time6=time6 )
          if ( all(time6(1:3)==(/2000,09,12/)) ) trange(2) = trange(1)
        end if

        ! reference time = start of forecast
        if ( present(tref) ) then
          call Get( t1, time6=time6 )
          time6(4:6) = 0
          tref = NewDate( time6=time6 ) + IncrDate( day=dd, hour=hh )
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! fc, 3 hourly
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_fc012up2tr3', '_03p03' )

        ! extract end hour:
        call Get( t2, hour=hour2 )

        ! set forecast start time and step given hour:
        select case ( hour2 )
          case ( 00 ) ; dd = -1  ;  hh = 12  ;  step = 12
          case ( 03 ) ; dd =  0  ;  hh = 00  ;  step =  3
          case ( 06 ) ; dd =  0  ;  hh = 00  ;  step =  6
          case ( 09 ) ; dd =  0  ;  hh = 00  ;  step =  9
          case ( 12 ) ; dd =  0  ;  hh = 00  ;  step = 12
          case ( 15 ) ; dd =  0  ;  hh = 12  ;  step =  3
          case ( 18 ) ; dd =  0  ;  hh = 12  ;  step =  6
          case ( 21 ) ; dd =  0  ;  hh = 12  ;  step =  9
          case default
            write (gol,'("unsupported hour :")'); call goErr
            write (gol,'("  hour2     : ",i2)') hour2; call goErr
            write (gol,'("  timesteps : ",a )') ec_treskey; call goErr
            call wrtgol( '  time1     : ', t1 ); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! file ccyymmdd contains fields for (00,24] :
        if ( present(tfile) ) then
          ! current day by default:
          tfile = t2
          call Set( tfile, hour=0, min=0, sec=0, mili=0 ) ! 00:00
          ! trap (..,00:00], this should be previous day:
          call Get( t2, time6=time6 )
          if ( all(time6(4:6)==0) ) tfile = tfile - IncrDate(day=1)
        end if

        ! fields valid for (00,24] :
        if ( present(trange) ) then
          trange(1) = t2
          call Set( trange(1), hour=0, min=0, sec=0, mili=0 ) ! 00:00
          ! trap 00:00, this should be previous day:
          call Get( t2, time6=time6 )
          if ( all(time6(4:6)==0) ) trange(1) = trange(1) - IncrDate(day=1)
          ! complete (00,24]
          trange(2) = trange(1) + IncrDate(day=1)  ! 24:00
          call Set( trange(1), mili=1 ) ! > 00:00
        end if

        ! reference time = start of forecast
        if ( present(tref) ) then
          call Get( t2, time6=time6 )
          time6(4:6) = 0
          tref = NewDate( time6=time6 ) + IncrDate( day=dd, hour=hh )
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! fc, 6 hourly
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_06p06' )

        ! extract end hour:
        call Get( t2, hour=hour2 )

        ! set forecast start time and step given hour:
        select case ( hour2 )
          case ( 00 ) ; dd = -1  ;  hh = 12  ;  step = 12
          case ( 06 ) ; dd =  0  ;  hh = 00  ;  step =  6
          case ( 12 ) ; dd =  0  ;  hh = 00  ;  step = 12
          case ( 18 ) ; dd =  0  ;  hh = 12  ;  step =  6
          case default
            write (gol,'("unsupported hour :")'); call goErr
            write (gol,'("  hour2     : ",i2)') hour2; call goErr
            write (gol,'("  timesteps : ",a )') ec_treskey; call goErr
            call wrtgol( '  time1     : ', t1 ); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! file ccyymmdd contains fields for (00,24] :
        if ( present(tfile) ) then
          ! current day by default:
          tfile = t1
          call Set( tfile, hour=00 )
          ! trap 00:00, this should be previous day:
          call Get( t1, time6=time6 )
          if ( all(time6(4:6)==0) ) tfile = tfile - IncrDate(day=1)
        end if

        ! fields valid for (00,24] :
        if ( present(trange) ) then
          trange(1) = t1
          call Set( trange(1), hour=0, min=0, sec=0, mili=0 ) ! 00:00
          ! trap 00:00, this should be previous day:
          call Get( t1, time6=time6 )
          if ( all(time6(4:6)==0) ) trange(1) = trange(1) - IncrDate(day=1)
          ! complete (00,24]
          trange(2) = trange(1) + IncrDate(day=1)  ! 24:00
          call Set( trange(1), mili=1 ) ! > 00:00
        end if

        ! reference time = start of forecast
        if ( present(tref) ) then
          call Get( t1, time6=time6 )
          time6(4:6) = 0
          tref = NewDate( time6=time6 ) + IncrDate( day=dd, hour=hh )
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! analysis, 6 hourly
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'an_06p06' )

        ! extract end hour:
        call Get( t2, hour=hour2 )

        ! set forecast start time and step given hour:
        select case ( hour2 )
          case ( 00 ) ; dd =  0  ;  hh = 00  ;  step = 00
          case ( 06 ) ; dd =  0  ;  hh = 06  ;  step = 00
          case ( 12 ) ; dd =  0  ;  hh = 12  ;  step = 00
          case ( 18 ) ; dd =  0  ;  hh = 18  ;  step = 00
          case default
            write (gol,'("unsupported hour :")'); call goErr
            write (gol,'("  hour2     : ",i2)') hour2; call goErr
            write (gol,'("  timesteps : ",a )') ec_treskey; call goErr
            call wrtgol( '  time1     : ', t1 ); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! file ccyymmdd contains fields for (00,24] :
        if ( present(tfile) ) then
          ! current day by default:
          tfile = t1
          call Set( tfile, hour=00 )
          ! trap 00:00, this should be previous day:
          call Get( t1, time6=time6 )
          if ( all(time6(4:6)==0) ) tfile = tfile - IncrDate(day=1)
        end if

        ! fields valid for (00,24] :
        if ( present(trange) ) then
          trange(1) = t1
          call Set( trange(1), hour=0, min=0, sec=0, mili=0 ) ! 00:00
          ! trap 00:00, this should be previous day:
          call Get( t1, time6=time6 )
          if ( all(time6(4:6)==0) ) trange(1) = trange(1) - IncrDate(day=1)
          ! complete (00,24]
          trange(2) = trange(1) + IncrDate(day=1)  ! 24:00
          call Set( trange(1), mili=1 ) ! > 00:00
        end if

        ! reference time = start of forecast
        if ( present(tref) ) then
          call Get( t1, time6=time6 )
          time6(4:6) = 0
          tref = NewDate( time6=time6 ) + IncrDate( day=dd, hour=hh )
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ???
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported time resolution key:")'); call goErr
        write (gol,'("  ",a)') trim(ec_treskey); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return

    end select

    ! ok
    status = 0

  end subroutine GetGribTime


  ! ***


  !
  ! Extract time fields of current grib record,
  ! store in mf%tref, mf%t1, mf%t2
  !


  subroutine SetTime( mf, status )

    use GO, only : TDate, NewDate, IncrDate, operator(+), wrtgol
    use file_grib, only : TGribFile, Check, Get

    ! --- const -------------------------------------

    character(len=*), parameter ::  rname = mname//'/SetTime'

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(inout)   ::  mf
    integer, intent(out)                   ::  status

    ! --- local -------------------------------

    integer          ::  reftime(5), timerange(4)

    ! --- begin -------------------------------

    ! extract time fields from grib record:
    call Get( mf%grib, status, reftime=reftime, timerange=timerange )
    IF_NOTOK_RETURN(status=1)

    ! Fill t1 and t2 with the time information; might be equal.
    ! Check time range indicator (WMO code table 5):
    select case ( timerange(4) )

      case ( 0, 1 )
        !
        ! 0 = Forecast product valid for reference time + P1 (P1>0),
        !     or uninitialized analysis product for reference time (P1=0)
        !
        ! 1 = Initialized analysis product for reference time (P1=0).
        !

        ! fill reference time:
        mf%tref = NewDate( time5=reftime )

        ! fill t1 with reftime+timerange;
        ! add P1 in hours; check time unit (WMO code table 4)
        mf%t1 = NewDate( time5=reftime )
        select case ( timerange(1) )
          case ( 1 )  ! hours
            mf%t1 = mf%t1 + IncrDate( hour=timerange(2) )
          case default
            write (gol,'("grib timerange units other than hours not supported yet")'); call goErr
            write (gol,'("  reftime   : ",i4,4i3)') reftime; call goErr
            write (gol,'("  timerange : ",4i3)') timerange; call goErr
            write (gol,'("  file      : ",a)') trim(mf%grib%fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! instant time:
        mf%t2 = mf%t1

      case ( 2 )
        !
        ! 2 = Product with a valid time ranging between
        !      reference time + P1 and reference time + P2
        !

        ! fill t1 with reftime+P1;
        ! add P1 in hours; check time unit (WMO code table 4)
        mf%t1 = NewDate( time5=reftime )
        select case ( timerange(1) )
          case ( 1 )  ! hours
            mf%t1 = mf%t1 + IncrDate( hour=timerange(2) )
          case default
            write (gol,'("grib timerange units other than hours not supported yet")'); call goErr
            write (gol,'("  file : ",a)') trim(mf%grib%fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! fill t2 with reftime+P1;
        ! add P1 in hours; check time unit (WMO code table 4)
        mf%t2 = NewDate( time5=reftime )
        select case ( timerange(1) )
          case ( 1 )  ! hours
            mf%t2 = mf%t2 + IncrDate( hour=timerange(3) )
          case default
            write (gol,'("grib timerange units other than hours not supported yet")'); call goErr
            write (gol,'("  file : ",a)') trim(mf%grib%fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

      case default
        write (gol,'("unsupported time range indicator:")'); call goErr
        write (gol,'("  indicator : ",i6)') timerange(4); call goErr
        write (gol,'("  file      : ",a)') trim(mf%grib%fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return

    end select

    ! ok
    status = 0

  end subroutine SetTime


  ! ***


  subroutine CheckTime( mf, tref, t1, t2, status )

    use GO
    use file_grib, only : Check, Get

    ! --- const -------------------------------------

    character(len=*), parameter ::  rname = mname//'/CheckTime'

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tmpp), intent(in)      ::  mf
    type(TDate), intent(in)                ::  tref, t1, t2
    integer, intent(out)                   ::  status

    ! --- local -------------------------------

    integer          ::  year1, year2

    ! --- begin -------------------------------

    ! requested year zero ? always ok
    call Get( t1, year=year1 )
    call Get( t2, year=year2 )
    if ( (year1 == 0) .and. (year2 == 0) ) then
      ! requested constant field, always ok
      status = 0; return    ! ok
    end if

    if ( year1 == 0 ) then
      ! do not test t1, only t2
      if ( t2 == mf%t2 ) then
        status = 0; return    ! ok
      else
        status = 1; return    ! not ok, try next
      end if
    else if ( year2 == 0 ) then
      ! do not test t2, only t1
      if ( t1 == mf%t1 ) then
        status = 0; return    ! ok
      else
        status = 1; return    ! not ok, try next
      end if
    end if

!    ! interval or instant time
!    if ( t1 < t2 ) then
!
!      !! time interval: [t1,t2] should be inside [t1,t2]
!      !if ( (t1 >= mf%t1) .and. (t2 <= mf%t2) ) then
!      !  status = 0; return  ! ok
!      !else
!      !  status = 1; return  ! try next
!      !end if
!
!      ! time interval: [t1,t2] should be equal to [t1,t2]
!      if ( (t1 == mf%t1) .and. (t2 == mf%t2) ) then
!        status = 0; return  ! ok
!      else
!        status = 1; return  ! try next
!      end if
!
!    else if ( t1 == t2 ) then
!
!      ! instant time: t2 should match t2
!      if ( t2 == mf%t2 ) then
!        status = 0; return  ! ok
!      else
!        status = 1; return  ! try next
!      end if
!
!    else
!
!      write (gol,'("t1 should not exceed t2:")'); call goErr
!      call wrtgol( '  t1 : ', t1 ); call goErr
!      call wrtgol( '  t2 : ', t2 ); call goErr
!      write (*,'("ERROR in ",a)') rname; status=3; return
!
!    end if


    ! compare all:
    if ( (tref == mf%tref) .and. (t1 == mf%t1) .and. (t2 == mf%t2) ) then
      status = 0; return  ! ok
    else
      status = 1; return  ! try next
    end if

    ! some error ...
    status = 1

  end subroutine CheckTime


end module tmm_mf_ecmwf_tmpp
