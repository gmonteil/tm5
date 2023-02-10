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

module tmm_mf_ecmwf_tm5

  use GO       , only : gol, goErr, goPr
  use GO       , only : TDate
  use os_specs , only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

  implicit none

  ! --- in/out ----------------------------

  private

  public  ::  TMeteoFile_ecmwf_tm5
  public  ::  Init, Done
  public  ::  Get
  public  ::  ReadRecord

  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf_ecmwf_tm5'


  !--- type ---------------------------------

  type TMeteoFile_ecmwf_tm5
    ! file name:
    character(len=MAX_FILENAME_LEN)         ::  dir
    character(len=MAX_FILENAME_LEN)         ::  fname
    ! time range covered by file:
    type(TDate)                ::  trange(2)
    ! other time keys for this file:
    character(len=16)          ::  treskey
    ! current time range covered by grib record:
    type(TDate)                ::  tref, t1, t2
    !
    ! file description
    !
    character(len=16)          ::  ec_class, ec_type, ec_levs
    integer                    ::  ec_sh, ec_gg
    character(len=MAX_RCKEY_LEN)         ::  paramkeys
    !
  end type TMeteoFile_ecmwf_tm5


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


  subroutine mf_Init( mf, dir, archivekeys, paramkey, &
                               tday, t1, t2, status )

    use GO, only : TDate
    use GO, only : goVarValue

    ! --- in/out ----------------------------

    type(TMeteoFile_ecmwf_tm5), intent(out)  ::  mf
    character(len=*), intent(in)        ::  dir
    character(len=*), intent(in)        ::  archivekeys
    character(len=*), intent(in)        ::  paramkey
    type(TDate), intent(in)             ::  tday, t1, t2
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Init'

    ! --- local --------------------------------

    ! --- begin --------------------------------

    ! store
    mf%dir  = dir

    !
    ! extract fields from archivekey :
    !   form=tmpp;class=od;type=fg;levs=ml60;sh=159;gg=80;tres=_fg006up4tr3
    !
    mf%ec_class   = 'od'
    call goVarValue( archivekeys, ';', 'class', '=', mf%ec_class, status )
    IF_ERROR_RETURN(status=1)
    !
    mf%ec_type    = 'fc'
    call goVarValue( archivekeys, ';', 'type', '=', mf%ec_type, status )
    IF_ERROR_RETURN(status=1)
    !
    mf%ec_levs    = 'ml60'
    call goVarValue( archivekeys, ';', 'levs', '=', mf%ec_levs, status )
    IF_ERROR_RETURN(status=1)
    !
    mf%ec_sh      = 159
    call goVarValue( archivekeys, ';', 'sh', '=', mf%ec_sh, status )
    IF_ERROR_RETURN(status=1)
    !
    mf%ec_gg      = 80
    call goVarValue( archivekeys, ';', 'gg', '=', mf%ec_gg, status )
    IF_ERROR_RETURN(status=1)
    !
    mf%treskey    = '_fc012up2tr3'
    call goVarValue( archivekeys, ';', 'tres', '=', mf%treskey, status )
    IF_ERROR_RETURN(status=1)

    ! specials
    select case ( paramkey )
      case ( 'oro', 'lsm' )
        ! overwrite timeresolutionkey, used later on to set trange
        mf%treskey        = 'const'
        ! tmm_convec tries to read oro using the default sourcekey,
        ! which probably contains type=fc ; force to use an for oro ...
        mf%ec_type = 'an'
    end select

    ! single parameter in a file:
    mf%paramkeys   = '-'//trim(paramkey)//'-'

    ! extract time range:
    call GetGribTime( mf%treskey, tday, t1, t2, status, trange=mf%trange )
    IF_NOTOK_RETURN(status=1)

    ! dummy filename, might be used in error diagnose
    write (mf%fname,'("ecmwf mars grib file for param ",a)') trim(paramkey)

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    ! --- in/out ------------------------------------

    type(TMeteoFile_ecmwf_tm5), intent(inout)  ::  mf
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

    type(TMeteoFile_ecmwf_tm5), intent(in)   ::  mf
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


  subroutine mf_ReadRecord( mf, paramkey, tday, t1, t2, nuv, nw, &
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

    type(TMeteoFile_ecmwf_tm5), intent(inout)  ::  mf
    character(len=*), intent(in)          ::  paramkey
    type(TDate), intent(in)               ::  tday, t1, t2
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

    integer               ::  iveg

    logical               ::  unit_change
    real                  ::  unit_fac

    ! --- begin ---------------------------------

    ! combined field ?
    select case ( paramkey )

      ! *** surface stress

      case ( 'sstr' )

        ! read first field:
        call mf_ReadRecord_1( mf, 'ewss', tday, t1, t2, nuv, nw, &
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
        call mf_ReadRecord_1( mf, 'nsss', tday, t1, t2, nuv, nw, &
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
            write (gol,'("unsupported gridtype for surface stress :",a)') gridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! clear pointers:
        call pa_Done( ll2 )
        call pa_Done( gg2 )
        call pa_Done( sh2 )

      ! *** vegetation types

      case ( 'tv01', 'tv02', 'tv03', 'tv04', 'tv05', 'tv06', 'tv07', 'tv08', 'tv09', 'tv10', &
             'tv11', 'tv12', 'tv13', 'tv14', 'tv15', 'tv16', 'tv17', 'tv18', 'tv19', 'tv20'  )

        ! extract number from name
        read (paramkey(3:4),'(i2)') iveg

        ! low vegetation types:
        call mf_ReadRecord_1( mf, 'tvl', tday, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! set elements that match requested vegetation type to 100%, zero elsewhere
        select case ( gridtype )
          case ( 'gg' )
            where ( nint(gg(:,1)) == iveg )
              gg(:,1) = 100.0  ! %
            elsewhere
              gg(:,1) = 0.0
            end where
          case default
            write (gol,'("unsupported gridtype for vegetation fractions :",a)') gridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! init pointer:
        call pa_Init( ll2 )
        call pa_Init( gg2 )
        call pa_Init( sh2 )

        ! high vegetation types:
        call mf_ReadRecord_1( mf, 'tvh', tday, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll2, sp_ll, &
                                ggi, gg2, sp_gg, &
                                shi, sh2, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! set elements that match requested vegetation type to 100%:
        select case ( gridtype )
          case ( 'gg' )
            where ( nint(gg2(:,1)) == iveg )
              gg(:,1) = 100.0  ! %
            end where
          case default
            write (gol,'("unsupported gridtype for vegetation fractions :",a)') gridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! clear pointers:
        call pa_Done( ll2 )
        call pa_Done( gg2 )
        call pa_Done( sh2 )

      ! *** default

      case default

        call mf_ReadRecord_1( mf, paramkey, tday, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )
         IF_NOTOK_RETURN(status=1)

    end select

    ! unit change ?
    unit_change = .true.
    select case ( paramkey )
      case ( 'lsm' ) ; unit_fac = 100.0    ! 0-1 -> 0-100%
      case default   ; unit_change = .false.
    end select
    ! apply ?
    if ( unit_change ) then
      select case ( gridtype )
        case ( 'll' ) ; ll = ll * unit_fac
        case ( 'gg' ) ; gg = gg * unit_fac
        !case ( 'sh' ) ; sh = sh * unit_fac
        case default
          write (gol,'("unsupported gridtype for unit change :",a)') gridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select
    end if

    ! ok
    status = 0

  end subroutine mf_ReadRecord



  ! ***


  subroutine mf_ReadRecord_1( mf, paramkey, tday, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use parray    , only : pa_Init, pa_Done
    use GO        , only : TDate, TIncrDate, operator(<), operator(-), rTotal, wrtgol
    use Grid      , only : TLevelInfo
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use tmm_info  , only : TMeteoInfo

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tm5), intent(inout)  ::  mf
    character(len=*), intent(in)          ::  paramkey
    type(TDate), intent(in)               ::  tday, t1, t2
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
    type(TIncrDate)       ::  tshift
    type(TDate)           ::  trefs, t1s, t2s
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

      case ( 'lsp', 'cp', 'sf', 'sshf', 'slhf', &
             'ssr', 'ssrd', 'str', 'strd', &
             'ewss', 'nsss', &
             'UDMF', 'UDDR', 'DDMF', 'DDDR' )

        ! get reference time for requested time interval:
        call GetGribTime( mf%treskey, tday, t1, t2, status, tref=tref )
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
        call wrtgol( '    accum ', tref, ' - ', t2 ); call goPr
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
          call wrtgol( '    accum ', tref, ' - ', t1 ); call goPr
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

        ! get reference time for requested time interval;
        ! eventually shift for analysed fields in case of forecasts:
        call GetGribTime( mf%treskey, tday, t1, t2, status, tref=tref, tshift=tshift )
        IF_NOTOK_RETURN(status=1)

        ! shift times (might be zero):
        trefs = tref - tshift
        t1s   = t1 - tshift
        t2s   = t2 - tshift

        ! just read ..
        call mf_ReadRecord_2( mf, paramkey, trefs, t1s, t2s, nuv, nw, &
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

    use GO        , only : TDate, wrtgol, Get, NewDate, operator(>)
    use GO        , only : goWriteKeyNum
    use Grid      , only : Init, Done
    use Grid      , only : TLevelInfo
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use Grid      , only : Interpol
    use file_grib , only : TGribFile
    use file_grib , only : Init, Done, ReadRecord, Get, Check
    use file_grib , only : levtype_sfc, levtype_hyb, levtype_land
    use file_grib , only : gridtype_ll, gridtype_gg, gridtype_sh
    use grib_table, only : GetPid
    use PArray    , only : pa_Init, pa_Done, pa_SetShape
    use tmm_info  , only : TMeteoInfo, Init, AddHistory

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tm5), intent(inout)  ::  mf
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

    character(len=16)      ::  ec_class, ec_type
    character(len=16)      ::  ec_grid, gridN, gridT
    character(len=16)      ::  levs
    character(len=16)      ::  treskey
    logical                ::  constant

    type(TGribFile)        ::  grib
    logical                ::  do_spm
    character(len=MAX_FILENAME_LEN)     ::  spm_fname
    type(TGribFile)        ::  spm_grib
    logical                ::  spm_lnsp
    logical                ::  spm_lnsp2sp

    integer                ::  pid
    character(len=7)       ::  gribcode
    character(len=4)       ::  timekey
    character(len=16)      ::  spm_levs, spm_paramkey, ec_paramkey

    type(TDate)            ::  tfile
    integer                ::  ccyy, mm, dd, hh
    type(TDate)            ::  tc

    logical                ::  exist

    logical                ::  isfirst
    logical                ::  reopened
    integer                ::  nlev, glevtype, glevel

    integer                ::  level

    integer                ::  nlev_out, level_out

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

          !write (gol,'("mf_ReadRecord_2:  paramkey : ",a)') trim(paramkey); call goPr
          !call wrtgol( 'mf_ReadRecord_2:  tref     : ', tref ); call goPr
          !call wrtgol( 'mf_ReadRecord_2:  t1       : ', t1 ); call goPr
          !call wrtgol( 'mf_ReadRecord_2:  t2       : ', t2 ); call goPr

    ! no fluxes through boundaries ...
    if ( nuv /= 'n' ) then
      write (gol,'("unsupported nuv key : ",a)') nuv; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! limitted support of fluxes ..
    !if ( nw /= 'n' ) then
    !  write (gol,'("unsupported nw key : ",a)') nw; call goErr
    !  write (gol,'("in ",a)') rname; call goErr; status=1; return
    !end if

    ! init pointer arrays:
    call pa_Init( pat )

    !
    ! ~~~ 3d field or 2d stored in 3d array
    !

    ! grid : T159, N80, etc
    call goWriteKeyNum( gridT, 'T', mf%ec_sh )
    call goWriteKeyNum( gridN, 'N', mf%ec_gg )

    ! defaults
    ec_paramkey    = paramkey
    ec_class       = mf%ec_class
    ec_type        = mf%ec_type
    levs           = mf%ec_levs
    ec_grid        = gridN
    treskey        = mf%treskey
    constant       = .false.
    do_spm         = .false.
    spm_lnsp       = .false.
    spm_lnsp2sp    = .false.

    ! specials
    select case ( paramkey )
      case ( 'LNSP' )
        ec_grid        = gridT
      case ( 'VO', 'D' )
        ec_grid        = gridT
        do_spm         = .true.
        spm_lnsp       = .true.
      case ( 'T', 'W', 'Q', 'CLWC', 'CIWC', 'CC', 'UDMF', 'UDDR', 'DDMF', 'DDDR' )
        do_spm         = .true.
        spm_lnsp2sp    = (ec_class == 'e4') .or. (ec_class == 'ei')
      case ( 'oro', 'lsm' )
        levs           = 'sfc'
        constant       = .true.
        treskey        = 'const'
        ! tmm_convec tries to read oro using the default sourcekey,
        ! which probably contains type=fc ; force to use an for oro ...
        ec_type = 'an'
!      case ( 'cvl', 'cvh', 'tvl', 'tvh', 'sr', 'albedo', 'lsrh' )
!        ec_type        = 'an'
!        levs           = 'sfc'
      case ( 'ci', 'sst', 'swvl1', 'swvl2', 'swvl3', 'swvl4', '10fg', 'sd', 'lsp', &
             'cp', 'sf', 'sshf', 'slhf', 'blh', 'u10m', 'v10m', 't2m', 'd2m', &
             'ssr', 'ewss', 'nsss', 'sstr' ,'src', 'skt' )
        levs           = 'sfc'
      case ( 'sp' )
        spm_lnsp2sp    = (ec_class == 'e4') .or. (ec_class == 'ei')
        if ( spm_lnsp2sp ) then
          ec_paramkey = 'LNSP'
        else
          levs = 'sfc'
        end if
    end select

    ! write gribcode
    call GetPid( 'ec', ec_paramkey, pid, status )
    IF_NOTOK_RETURN(status=1)
    write (gribcode,'(i3,".128")') pid
    gribcode = adjustl(gribcode)

    ! convert input times to file name times:
    call GetGribTime( treskey, tref, t1, t2, status, tfile=tfile )
    IF_NOTOK_RETURN(status=1)

          !call wrtgol( 'mf_ReadRecord_2:  tfile    : ', tfile ); call goPr

    ! extract time values:
    call Get( tfile, year=ccyy, month=mm, day=dd, hour=hh )
    ! in file name: 0000, 0600, 1200, 1800
    write (timekey,'(i4.4)') hh*100
    timekey = adjustl(timekey)

    ! create file name:
    !   dir/od-fc-20000101-1200-ml60-138-T159.gb
    !
    ! filename includes date:
    write (mf%fname,'(a,"/",a,"-",a,"-",i4.4,2i2.2,"-",a,"-",a,"-",a,"-",a,".gb")') &
              trim(mf%dir), &
              trim(ec_class), trim(ec_type), ccyy, mm, dd, trim(timekey), &
              trim(levs), trim(gribcode), trim(ec_grid)

    ! exist ?
    inquire( file=mf%fname, exist=exist )
    if ( .not. exist ) then
      write (gol,'("grib file does not exist:")'); call goErr
      write (gol,'("  ",a)') trim(mf%fname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! open grib file
    call Init( grib, mf%fname, 'r', status )
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

      call ReadRecord( grib, status )
      select case ( status )
        case ( 0 )
          ! no error
        case ( 1 )
          ! eof
          if ( .not. reopened ) then
            !write (gol,'("grib read record: re-open ...")'); call goPr
            ! close:
            call Done( grib, status )
            IF_NOTOK_RETURN(status=1)
            ! reopen:
            call Init( grib, mf%fname, 'r', status )
            IF_NOTOK_RETURN(status=1)
            reopened = .true.
            cycle
          else
            write (gol,'("reached eof before requested record was found")'); call goErr
            write (gol,'("  file     : ",a)') trim(mf%fname); call goErr
            write (gol,'("  paramkey : ",a)') trim(paramkey); call goErr
            call wrtgol( '  tref     : ', tref ); call goErr
            call wrtgol( '  t1       : ', t1   ); call goErr
            call wrtgol( '  t2       : ', t2   ); call goErr
            write (gol,'("tips:")'); call goErr
            write (gol,'("  o grib file corrupted or zero ?")'); call goErr
            write (gol,'("  o if accumulatd field,")'); call goErr
            write (gol,'("    check list of accumulated fields in mf_ReadRecord_1")'); call goErr
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
        case ( 'spm' ) ; call GetPid( 'ec', 'SP'       , pid, status )
        case default   ; call GetPid( 'ec', ec_paramkey, pid, status )
      end select
      IF_NOTOK_RETURN(status=1)

      ! check parameter; continue if not ok:
      call Check( grib, pid=pid, debug=0, status=status )
      if (status/=0) cycle

      ! fill times ?
      if ( .not. constant ) then

        ! extract time fields from grib, store in mf%tref/mf%t1/mf%t2
        call SetTime( mf, grib, status )
        IF_NOTOK_RETURN(status=1)

        ! check time:
        call CheckTime( mf, tref, t1, t2, status )
        if (status/=0) then
          !write (gol,'("grib read record: wrong time; skip ...")'); call goPr
          !write (gol,'("  paramkey     : ",a)') paramkey; call goPr
          !call wrtgol( '  req. tref : ', tref ); call goPr
          !call wrtgol( '       t1   : ', t1   ); call goPr
          !call wrtgol( '       t2   : ', t2   ); call goPr
          !call wrtgol( '  grib tref : ', mf%tref ); call goPr
          !call wrtgol( '       t1   : ', mf%t1   ); call goPr
          !call wrtgol( '       t2   : ', mf%t2   ); call goPr
          !write (gol,'("  grib file    : ",a)') trim(mf%fname); call goPr
          cycle
        end if

      end if  ! time checking

      ! extract level stuff:
      call Get( grib, nlev=nlev, levtype=glevtype, level=glevel, status=status )
      IF_NOTOK_RETURN(status=1)
      ! check level type:
      select case ( glevtype )
        case ( levtype_sfc, levtype_land )
          ! surface field
          nlev = 1
          glevel = 1
        case ( levtype_hyb )
          select case ( paramkey )
            case ( 'LNSP', 'sp' )
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

      ! number of output levels:
      nlev_out = nlev
      if ( nw == 'w' ) nlev_out = nlev_out + 1

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
        call Get( grib, status, reftime=greftime, timerange=gtimerange )
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
        call Get( grib, status, gridtype=ggridtype )
        IF_NOTOK_RETURN(status=1)

        ! setup:
        select case ( ggridtype )

          ! o lat/lon
          case ( gridtype_ll )

            ! routine returns lat/lon grid:
            gridtype = 'll'

            ! grib storage is north pole to south pole:
            call Get( grib, status, &
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
            call pa_SetShape( ll, lon_n, lat_n, nlev_out )
            ll = 0.0

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
            call Get( grib, status, N=ggN )
            IF_NOTOK_RETURN(status=1)

            ! define grid structure:
            call Init( ggi, ggN, .true., status )
            IF_NOTOK_RETURN(status=1)

            ! allocate output:
            call pa_SetShape( gg, ggi%np, nlev_out )
            gg = 0.0

            ! add to history:
            write (key,'("gg==",i4.4)') ggN
            call AddHistory( tmi, trim(key), status )

          ! o spectral field:
          case ( gridtype_sh )

            ! routine returns sh grid:
            gridtype = 'sh'

            ! extract spectral truncation:
            call Get( grib, status, T=shT )
            IF_NOTOK_RETURN(status=1)

            ! intialize spherical harmonic field info:
            call Init( shi, shT, status )
            IF_NOTOK_RETURN(status=1)

            ! allocate output:
            call pa_SetShape( sh, shi%np, nlev_out )
            sh = cmplx(0.0,0.0)

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
          case ( 91 )
            call Init( levi, 'ec91', status )
            IF_NOTOK_RETURN(status=1)
          case default
            write (gol,'("do not how to init levi for nlev = ",i6)') nlev; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! not again ...
        isfirst = .false.
      end if  ! isfirst (grid definition and allocation)

      !
      ! store
      !

      ! layers numbered 1..nlev, half levels numberd 1..nlev+1
      ! top-down, thus 1 is space and nlev+1 is surface
      if ( nw == 'w' ) then
        ! store half levels
        select case ( paramkey )
          ! store in upper half level of a layer:
          case ( 'UDMF', 'DDMF' )
            level_out = level
          !! store in lower half level of a layer:
          !case ( 'dummy' )
          !  level_out = level + 1
          ! to be implemented ...
          case default
            write (gol,'("do not if data is on upper or lower half level ...")'); call goErr
            write (gol,'("  paramkey : ",a)') paramkey; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
      else
        ! store full levels
        level_out = level
      end if

      select case ( ggridtype )

        case ( gridtype_ll )

          ! read 2d pat from grib; storred from north to south
          call Get( grib, status, ll=pat )
          IF_NOTOK_RETURN(status=1)

          ! store from south to north:
          do ilat = 1, lat_n
            ll(:,ilat,level_out) = pat(:,lat_n+1-ilat)
          end do

        case ( gridtype_gg )

          ! read 2d pat from grib:
          call Get( grib, status, gg=gg(:,level_out) )
          IF_NOTOK_RETURN(status=1)

          ! convert from lnsp to sp ?
          if ( paramkey == 'sp' .and. spm_lnsp2sp ) gg(:,level_out) = exp(gg(:,level_out)) ! Pa

        case ( gridtype_sh )

          ! read 2d pat from grib:
          call Get( grib, status, sh=sh(:,level_out) )
          IF_NOTOK_RETURN(status=1)

        case default
          write (gol,'("unsupported gridtype for 2d pat : ",i6)') gridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! last record for this field ?
      if ( glevel == nlev ) exit

    end do   ! records

    ! close grib file
    call Done( grib, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ~~~ surface pressure
    !

    if ( do_spm ) then

      ! read lnsp and covert to sp ?
      if ( spm_lnsp .or. spm_lnsp2sp) then
        spm_levs     = levs
        spm_paramkey = 'LNSP'
      else
        spm_levs     = 'sfc'
        spm_paramkey = 'SP'
      end if

      ! write gribcode:
      call GetPid( 'ec', spm_paramkey, pid, status )
      IF_NOTOK_RETURN(status=1)
      write (gribcode,'(i3,".128")') pid
      gribcode = adjustl(gribcode)

      ! create file name:
      !   dir/od-fc-2000-01-ml60-T159-T_20000101_fg006up4tr3.gb
      write (spm_fname,'(a,"/",a,"-",a,"-",i4.4,2i2.2,"-",a,"-",a,"-",a,"-",a,".gb")') &
              trim(mf%dir), &
              trim(ec_class), trim(ec_type), ccyy, mm, dd, trim(timekey), &
              trim(spm_levs), trim(gribcode), trim(ec_grid)

      ! exist ?
      inquire( file=spm_fname, exist=exist )
      if ( .not. exist ) then
        write (gol,'("grib file does not exist:")'); call goErr
        write (gol,'("  ",a)') trim(spm_fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

      ! open grib file
      call Init( spm_grib, spm_fname, 'r', status )
      IF_NOTOK_RETURN(status=1)

      ! loop over time records
      do

        ! read gribsection in file buffer
        call ReadRecord( spm_grib, status )
        IF_NOTOK_RETURN(status=1)

        ! fill times
        call SetTime( mf, spm_grib, status )
        IF_NOTOK_RETURN(status=1)

        ! check time:
        call CheckTime( mf, tref, t1, t2, status )
        if (status/=0) then
          !write (gol,'("grib read record: spm wrong time; skip ...")'); call goPr
          cycle
          !write (gol,'("found unexpected times in grib file:")'); call goErr
          !write (gol,'("  paramkey     : ",a)') paramkey; call goErr
          !call wrtgol( '  req. t1 : ', t1 ); call goErr
          !call wrtgol( '       t2 : ', t2 ); call goErr
          !call wrtgol( '  grib t1 : ', mf%t1 ); call goErr
          !call wrtgol( '       t2 : ', mf%t2 ); call goErr
          !write (gol,'("  grib file    : ",a)') trim(spm_fname); call goErr
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
          if ( spm_lnsp2sp ) then
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
      call Check( spm_grib, pid=pid, debug=1, status=status )
      IF_NOTOK_RETURN(status=1)

      ! check level:
      call Get( spm_grib, levtype=glevtype, level=glevel, status=status )
      IF_NOTOK_RETURN(status=1)
      select case ( ggridtype )
        case ( gridtype_ll  )
          if ( glevtype /= levtype_sfc ) then
            write (gol,'("found unexpected level type ")'); call goErr
            write (gol,'("  paramkey        : ",a)') paramkey; call goErr
            write (gol,'("  sfc level type  : ",i6)') levtype_sfc; call goErr
            write (gol,'("  grib level type : ",i6)') glevtype; call goErr
            write (gol,'("  grib file       : ",a)') trim(spm_fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
          end if
        case ( gridtype_gg )
          if ( spm_lnsp2sp ) then
            if ( (glevtype /= levtype_hyb) .or. (glevel /= 1) ) then
              write (gol,'("found unexpected level type (lnsp for 3d gg)")'); call goErr
              write (gol,'("  paramkey        : ",a)') paramkey; call goErr
              write (gol,'("  hyb level type  : ",i6)') levtype_hyb; call goErr
              write (gol,'("  grib level type : ",i6)') glevtype; call goErr
              write (gol,'("  grib level      : ",i6)') glevel; call goErr
              write (gol,'("  grib file       : ",a)') trim(spm_fname); call goErr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
            end if
          else
            if ( glevtype /= levtype_sfc ) then
              write (gol,'("found unexpected level type ")'); call goErr
              write (gol,'("  paramkey        : ",a)') paramkey; call goErr
              write (gol,'("  sfc level type  : ",i6)') levtype_sfc; call goErr
              write (gol,'("  grib level type : ",i6)') glevtype; call goErr
              write (gol,'("  grib file       : ",a)') trim(spm_fname); call goErr
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
            write (gol,'("  grib file       : ",a)') trim(spm_fname); call goErr
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
          call Get( spm_grib, status, ll=pat )
          IF_NOTOK_RETURN(status=1)
          ! store from south to north:
          do ilat = 1, lat_n
            sp_ll(:,ilat) = pat(:,lat_n+1-ilat)
          end do

        case ( gridtype_gg )
          ! allocate storage
          call pa_SetShape( sp_gg, ggi%np )
          ! read gg field from grib:
          call Get( spm_grib, status, gg=sp_gg )
          IF_NOTOK_RETURN(status=1)
          ! convert from lnsp to sp ?
          if ( spm_lnsp2sp ) sp_gg = exp(sp_gg)   ! Pa

        case ( gridtype_sh )
          ! allocate storage
          call pa_SetShape( lnsp_sh, shi%np )
          ! read spectral coeff:
          call Get( spm_grib, status, sh=lnsp_sh )
          IF_NOTOK_RETURN(status=1)

        case default
          write (gol,'("unsupported gridtype for reading sp/lnsp : ",i6)') ggridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! close grib file
      call Done( spm_grib, status )
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
  ! Called as:
  !   call GetGribTime( treskey, tday, t1, t2, status, tref=tref )
  !   call GetGribTime( treskey, tref, t1, t2, status, tfile=tfile )
  !

  subroutine GetGribTime( treskey, t0, t1, t2, status, tfile, trange, tref, tshift )

    use GO, only : TDate, TIncrDate, Get, Set, wrtgol, NewDate, IncrDate
    use GO, only : operator(+), operator(-), operator(<), rTotal, iTotal
    use GO, only : AnyDate, IsAnyDate
    use GO, only : wrtgol

    ! --- in/out --------------------------------

    character(len=*), intent(in)            ::  treskey
    type(TDate), intent(in)                 ::  t0, t1, t2
    integer, intent(out)                    ::  status

    type(TDate), intent(out), optional      ::  tfile
    type(TDate), intent(out), optional      ::  trange(2)
    type(TDate), intent(out), optional      ::  tref
    type(TIncrDate), intent(out), optional  ::  tshift

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GetGribTime'

    ! --- local --------------------------------

    integer          ::  hour1, hour2, time6(6)
    integer          ::  dd, hh, step
    integer          ::  anhh
    real             ::  ddr

    ! --- begin --------------------------------

          !write (gol,'("      GetGribTime: treskey  : ",a)') trim(treskey); call goPr
          !call wrtgol( '      GetGribTime: t0       : ', t0 ); call goPr
          !call wrtgol( '      GetGribTime: t1       : ', t1 ); call goPr
          !call wrtgol( '      GetGribTime: t2       : ', t2 ); call goPr

    ! files opend upon reading, thus no particular time range for which file is valid:
    if ( present(trange) ) then
      trange(1) = AnyDate()
      trange(2) = AnyDate()
    end if

    ! zero shift by default
    if ( present(tshift) ) tshift = IncrDate(hour=0)

    ! set day shift, start hour, and step
    select case ( treskey )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! constant field
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'const' )

        ! only t0 is usefull ...
        if ( present(tfile ) ) tfile  = t0
        !if ( present(trange) ) trange = (/t1,t2/)  ! any, any
        if ( present(tref  ) ) tref   = t0       ! dummy ...

        ! take analysed fields always at least 24 hour old,
        ! since these are the only analysed fields available in forecast mode
        if ( present(tshift) ) then
          !tshift = IncrDate(hour=24)
          ! FIX for start of ml91 test suite; try if 12 is ok too ...
          tshift = IncrDate(hour=12)
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! fc, 3 hourly
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_fc012up2tr3' )

        ! end hour counted from t0 00:00
        hour2 = iTotal( t2 - t0, 'hour' )

        ! set forecast start time and step given hour:
        if ( hour2 == 0 ) then
          ! 00:00   fc day 0
          dd = -1  ;  hh = 12  ;  step = 12
        else if ( (hour2 <= 12) .and. (modulo(hour2,3) == 0) ) then
          ! (00,12]   fc day 0
          dd   =  0
          hh   = 00
          step = hour2
        else if ( ( (t0 < NewDate(year=2006,month=03,day=14)                  ) .and. &
                    ( ((hour2 <= 12+ 72) .and. (modulo(hour2,3) == 0)) .or. &
                      ((hour2 <= 12+240) .and. (modulo(hour2,6) == 0))        )         ) &
             .or. ( ( ((hour2 <= 12+ 96) .and. (modulo(hour2,3) == 0)) .or. &
                      ((hour2 <= 12+240) .and. (modulo(hour2,6) == 0))        )         )   ) then
          ! (12,240]   fc days 1-10
          dd   =  0
          hh   = 12
          step = hour2 - 12
        else
          write (gol,'("unsupported hour :")'); call goErr
          write (gol,'("  hour2     : ",i3)') hour2; call goErr
          write (gol,'("  treskey   : ",a )') treskey; call goErr
          call wrtgol( '  time1     : ', t1 ); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if

        !! fields valid for hh+(00,12] :
        !if ( present(trange) ) then
        !  trange(1) = t2
        !  call Set( trange(1), hour=hh, min=0, sec=0, mili=0 ) ! hh:00
        !  ! trap 00:00, this should be previous day:
        !  call Get( t2, time6=time6 )
        !  if ( all(time6(4:6)==0) ) trange(1) = trange(1) - IncrDate(day=1)
        !  ! complete (00,12]
        !  trange(2) = trange(1) + IncrDate(hour=12)  ! 24:00
        !  call Set( trange(1), mili=1 ) ! > 00:00
        !end if

        ! reference time = start of forecast
        if ( present(tref) ) then
          tref = t0 + IncrDate(day=dd,hour=hh)
        end if

        ! adhoc: if tfile is requested, probably the 'tref' returned before
        ! is now in input 't0' ...
        if ( present(tfile) ) then
          tfile = t0  ! .. is tref !
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! analysis, files for hours 0, 6, 12, and 18
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_an0tr6' )

        ! reference time = analysis time
        if ( present(tref) ) then
          tref = t1
        end if

        !                           t0       t1,t2  ->
        !  -+-----------------------+-----------------------+--->
        !   00    06    12    18    00    06    12    18    00
        !   a
        !         a     a
        !                     a     a(-----------]                00 analysis/forecast
        !                                 a     a(--------------> 12 forecast
        !
        !                          -24    -24  -24   -24   -48   shift

        ! take analysed fields always at least 24 hour old,
        ! since these are the only analysed fields available in forecast mode,
        ! and to obtain a contineous time line
        ! t0 is always 00:00
        if ( present(tshift) ) then
          ! difference between t1 and t0 00:00 in fraction of days:
          ddr = rTotal( t1 - t0, 'day' )
          ! set time shift in days:
          tshift = IncrDate(day=floor(ddr)+1)
          !! FIX for start of ml91 test suite; try if 12 is ok too ...
          !tshift = IncrDate(day=floor(ddr),hour=12)
        end if

        ! one file for each time:
        if ( present(tfile) ) then
          tfile = t1
        end if

        !! fields in file valid for instant time:
        !if ( present(trange) ) then
        !  trange(1) = t1
        !  trange(2) = t1
        !end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ???
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported time resolution key:")'); call goErr
        write (gol,'("  `",a,"`")') trim(treskey); call goErr
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


  subroutine SetTime( mf, grib, status )

    use GO, only : TDate, NewDate, IncrDate, operator(+), wrtgol
    use file_grib, only : TGribFile, Check, Get

    ! --- const -------------------------------------

    character(len=*), parameter ::  rname = mname//'/SetTime'

    ! --- in/out -------------------------------

    type(TMeteoFile_ecmwf_tm5), intent(inout)   ::  mf
    type(TGribFile), intent(inout)         ::  grib
    integer, intent(out)                   ::  status

    ! --- local -------------------------------

    integer          ::  reftime(5), timerange(4)

    ! --- begin -------------------------------

    ! extract time fields from grib record:
    call Get( grib, status, reftime=reftime, timerange=timerange )
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
            write (gol,'("  file      : ",a)') trim(grib%fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! instant time:
        mf%t2 = mf%t1

      case ( 2 )
        !
        ! 2 = Product with a valid time ranging between
        !      reference time + P1 and reference time + P2
        !

        ! fill reftime:
        mf%tref = NewDate( time5=reftime )

        ! fill t1 with reftime+P1;
        ! add P1 in hours; check time unit (WMO code table 4)
        mf%t1 = mf%tref
        select case ( timerange(1) )
          case ( 1 )  ! hours
            mf%t1 = mf%t1 + IncrDate( hour=timerange(2) )
          case default
            write (gol,'("grib timerange units other than hours not supported yet")'); call goErr
            write (gol,'("  file : ",a)') trim(grib%fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! fill t2 with reftime+P2;
        ! add P2 in hours; check time unit (WMO code table 4)
        mf%t2 = mf%tref
        select case ( timerange(1) )
          case ( 1 )  ! hours
            mf%t2 = mf%t2 + IncrDate( hour=timerange(3) )
          case default
            write (gol,'("grib timerange units other than hours not supported yet")'); call goErr
            write (gol,'("  file : ",a)') trim(grib%fname); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

      case default
        write (gol,'("unsupported time range indicator:")'); call goErr
        write (gol,'("  indicator : ",i6)') timerange(4); call goErr
        write (gol,'("  file      : ",a)') trim(grib%fname); call goErr
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

    type(TMeteoFile_ecmwf_tm5), intent(in)      ::  mf
    type(TDate), intent(in)                ::  tref, t1, t2
    integer, intent(out)                   ::  status

    ! --- local -------------------------------

    integer          ::  year1, year2

    ! --- begin -------------------------------

    !call wrtgol( 'CheckTime: (', tref, ') ', t1, ' - ', t2 ); call goPr
    !call wrtgol( 'CheckTime: (', mf%tref, ') ', mf%t1, ' - ', mf%t2 ); call goPr

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
!      write (gol,'("in ",a)') rname; call goErr; status=3; return
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


end module tmm_mf_ecmwf_tm5
