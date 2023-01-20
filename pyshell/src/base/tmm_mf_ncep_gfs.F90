!###############################################################################
!
!  Interface to NCEP Global Forecast System data
!
!  Database
!  --------
!
!    From Phil Rasch and Dani Bundi (NCAR)
!
!  NCEP Cray Binary files
!  ----------------------
!
!  See http://dss.ucar.edu/pub/reanalysis//CRAY_bin.html .
!
!  Grib files
!  ----------
!
!  For grib codes, see:
!     http://www.nco.ncep.noaa.gov/pmb/docs/on388/
!
!  Spectral fields
!  -----------
!
!  From the 'discription' in the NetCDF files in the CDC archive
!  it was concluded that  'ncep spectral' / sqrt(2) = 'ecmwf spectral' .
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module tmm_mf_ncep_gfs

  use GO      , only : gol, goErr, goPr
  use GO      , only : TDate
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

  implicit none

  ! --- in/out ----------------------------

  private

  public  ::  TMeteoFile_ncep_gfs
  public  ::  Init, Done
  public  ::  Get
  public  ::  ReadRecord

  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf_ncep_gfs'

  !--- type ---------------------------------

  type TMeteoFile_ncep_gfs
    ! field collection
    character(len=MAX_RCKEY_LEN)     ::  paramkeys
    type(TDate)            ::  trange(2)
    ! file names
    character(len=MAX_FILENAME_LEN)     ::  dir
    character(len=MAX_FILENAME_LEN)     ::  fname
    character(len=16)      ::  ext
    ! time resolution
    character(len=32)      ::  treskey
    type(TDate)            ::  tref
  end type TMeteoFile_ncep_gfs


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


       ! --- adhoc ----------------------------------------

       integer                          ::  adhoc_fu    = 10
       character(len=MAX_FILENAME_LEN)               ::  adhoc_fname = 'none'


contains


  ! ==============================================================


  subroutine mf_Init( mf, dir, archivekeys, paramkey, &
                               tref, t1, t2, status )

    use GO, only : TDate, IsAnyDate, IncrDate, operator(+)
    use GO, only : goVarValue

    ! --- in/out ----------------------------

    type(TMeteoFile_ncep_gfs), intent(out)  ::  mf
    character(len=*), intent(in)            ::  dir
    character(len=*), intent(in)            ::  archivekeys
    character(len=*), intent(in)            ::  paramkey
    type(TDate), intent(in)                 ::  tref, t1, t2
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Init'

    ! --- local --------------------------------------

    ! --- begin --------------------------------

    ! store directory:
    mf%dir = trim(dir)//'/'

    !
    ! extract fields from archivekey :
    !   tres=fc024up1tr3
    !
    !mf%mdir  = 'no_mdir'
    !call goVarValue( archivekeys, ';', 'mdir', '=', mf%mdir, status )
    !if (status>0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! one specific set only yet ...
    mf%treskey = 'fc024up1tr3'

    ! store reftime:
    mf%tref = tref

    ! set file specific stuff:
    !   o list of parameters in file
    !   o file extension
    select case ( paramkey )
      ! * ncep cray binary files
      case ( 'oro', 'LNSP', 'Tv', 'D', 'VO', 'Q' )
        mf%paramkeys = '-oro-LNSP-Tv-D-VO-Q-'
        mf%ext = 'SF'
      ! * ncep grib files
      case ( 'lsm', 'albedo', 'sr', 'sps', 'ci', 'skt', 'sst', &
             'u10m', 'v10m', 'slhf', 'sshf', 'sstr', 'lsp', 'cp', &
             'blh' )
        mf%paramkeys = '-'
        mf%paramkeys = trim(mf%paramkeys)//'lsm-albedo-sr-sps-ci-skt-sst-'
        mf%paramkeys = trim(mf%paramkeys)//'u10m-v10m-slhf-sshf-sstr-lsp-cp-'
        mf%paramkeys = trim(mf%paramkeys)//'blh-'
        mf%ext = 'SFLUXGrbF'
      ! * error ...
      case default
        write (gol,'("do not know file for param ",a)') paramkey; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! files are valid for instant times only;
    ! averaged fields are stored at the end of the interval (t2) thus valid for (t1,t2]
    if ( IsAnyDate(t1) .and. IsAnyDate(t2) ) then
      mf%trange(1) = tref
      mf%trange(2) = tref
    else
      mf%trange(1) = t1 + IncrDate( mili=1 )
      mf%trange(2) = t2
    end if

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    ! --- in/out ------------------------------------

    type(TMeteoFile_ncep_gfs), intent(inout)   ::  mf
    integer, intent(out)                       ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Done'

    ! --- begin -------------------------------------

    ! files have been closed in ReadRecord/WriteRecord

    ! ok
    status = 0

  end subroutine mf_Done



  ! ***


  subroutine mf_Get( mf, status, trange1, trange2, paramkeys )

    use GO, only : TDate

    ! --- in/out ----------------------------

    type(TMeteoFile_ncep_gfs), intent(in)     ::  mf
    integer, intent(out)                      ::  status

    type(TDate), intent(out), optional        ::  trange1, trange2
    character(len=*), intent(out), optional   ::  paramkeys

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Get'

    ! --- local --------------------------------

    ! --- begin --------------------------------

    ! time range:
    if ( present(trange1) ) trange1 = mf%trange(1)
    if ( present(trange2) ) trange2 = mf%trange(2)

    ! parameter names:
    if ( present(paramkeys) ) paramkeys = mf%paramkeys

    ! ok
    status = 0

  end subroutine mf_Get


  ! ***

  !
  ! time arguments:
  !    tref , any, any   : oro and other constant fields
  !    dummy, t1 , t1    : file name with t1
  !    dummy, t1 , t2    : file name with t2
  !

  subroutine GetTime( treskey, tref, t1, t2, status, tfile, fcst_hr, aver_hr, taver )

    use go, only : TDate, Get, NewDate, IncrDate, IsAnyDate, operator(>), operator(-), operator(+)

    ! --- in/out ----------------------------------

    character(len=*), intent(in)            ::  treskey
    type(TDate), intent(in)                 ::  tref, t1, t2
    integer, intent(out)                    ::  status

    type(TDate), intent(out), optional      ::  tfile
    integer,  intent(out), optional         ::  fcst_hr
    integer,  intent(out), optional         ::  aver_hr
    type(TDate), intent(out), optional      ::  taver

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GetTime'

    ! --- local -----------------------------------

    integer              ::  year, month, day, hour
    integer              ::  dd, h0, dh, dha
    type(TDate)          ::  tfc

    ! --- begin -----------------------------------

    !
    ! set time values
    !

    ! resolution depends on data set:
    select case ( treskey )

      !case ( 'rasch-bundy' )
      !
      !  ! set forecast hour and time step:
      !  select case ( hour )
      !    case ( 00 ) ; h0 = 00 ; dh = 00
      !    case ( 03 ) ; h0 = 00 ; dh = 03
      !    case ( 06 ) ; h0 = 06 ; dh = 00
      !    case ( 09 ) ; h0 = 06 ; dh = 03
      !    case ( 12 ) ; h0 = 12 ; dh = 00
      !    case ( 15 ) ; h0 = 12 ; dh = 03
      !    case ( 18 ) ; h0 = 12 ; dh = 06
      !    case ( 21 ) ; h0 = 12 ; dh = 09
      !    case default
      !      write (gol,'("unsupported hour : ",i3)') hour; call goErr
      !      write (gol,'("in ",a)') rname; call goErr; status=1; return
      !  end select

      !
      ! ** Rasch-Bundy set ; combination of 'forecast_0' and forecast_1 files
      !
      !    instant field    : analysis or forecast
      !    interval [t1,t2] : forecast file stamped at t2
      !

      case ( 'fc024up1tr3' )

        if ( IsAnyDate(t1) .and. IsAnyDate(t2) ) then
          ! constant field; extract time stuff from tref:
          call Get( tref, year=year, month=month, day=day, hour=hour )
        else
          ! extract time stuff from t2, which is either t1 or end of interval
          call Get( t2, year=year, month=month, day=day, hour=hour )
          ! accumulated field for 12:00 is in files 12.24 ...
          if ( (t2>t1) .and. (hour==12) ) hour = 12+24
        end if

        ! set forecast hour, time step:
        select case ( hour )
          case ( 12    ) ; dd = 0 ; h0 = 12 ; dh = 00 ; dha=00
          case ( 15    ) ; dd = 0 ; h0 = 12 ; dh = 03 ; dha=00
          case ( 18    ) ; dd = 0 ; h0 = 12 ; dh = 06 ; dha=00
          case ( 21    ) ; dd = 0 ; h0 = 12 ; dh = 09 ; dha=06
          case ( 00    ) ; dd = 1 ; h0 = 12 ; dh = 12 ; dha=06
          case ( 03    ) ; dd = 1 ; h0 = 12 ; dh = 15 ; dha=12
          case ( 06    ) ; dd = 1 ; h0 = 12 ; dh = 18 ; dha=12
          case ( 09    ) ; dd = 1 ; h0 = 12 ; dh = 21 ; dha=18
          case ( 12+24 ) ; dd = 1 ; h0 = 12 ; dh = 24 ; dha=18
          case default
            write (gol,'("unsupported hour : ",i3)') hour; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

      case default
        write (gol,'("unsupported time resolution key : ",a)') treskey; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! fill start of forecast:
    tfc = NewDate( year=year, month=month, day=day, hour=h0 ) - IncrDate(day=dd)

    !
    ! return output
    !

    ! return time for filename ?
    if ( present(tfile) ) tfile = tfc

    ! return time step ?
    if ( present(fcst_hr) ) fcst_hr = dh

    ! return start of averaging interval ?
    if ( present(aver_hr) ) aver_hr = dha

    ! start of aver interval ?
    if ( present(taver) ) taver = tfc + IncrDate(hour=dha)

    ! ok
    status = 0

  end subroutine GetTime


  ! ***


  ! Return a field given parameter name, time, etc.
  ! Only one of grid types is filled!

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

    type(TMeteoFile_ncep_gfs), intent(inout)  ::  mf
    character(len=*), intent(in)              ::  paramkey
    type(TDate), intent(in)                   ::  t1, t2
    character(len=1), intent(in)              ::  nuv
    character(len=1), intent(in)              ::  nw
    character(len=2), intent(out)             ::  gridtype
    type(TLevelInfo), intent(out)             ::  levi
    type(TllGridInfo), intent(inout)          ::  lli
    real, pointer                             ::  ll(:,:,:)
    real, pointer                             ::  sp_ll(:,:)
    type(TggGridInfo), intent(inout)          ::  ggi
    real, pointer                             ::  gg(:,:)
    real, pointer                             ::  sp_gg(:)
    type(TshGridInfo), intent(inout)          ::  shi
    complex, pointer                          ::  sh(:,:)
    complex, pointer                          ::  lnsp_sh(:)
    type(TMeteoInfo), intent(out)             ::  tmi
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord'

    ! --- local -------------------------------

    real, pointer         ::  ll2(:,:,:)
    real, pointer         ::  gg2(:,:)
    complex, pointer      ::  sh2(:,:)

    ! --- begin ---------------------------------

    ! combined or other special field ?
    select case ( paramkey )

      !
      ! total surface stress :  sstr^2 =  ewss^2 + nsss^2
      !
      case ( 'sstr' )

!        ! read first field:
!        call mf_ReadRecord_1( mf, 'ewss', t1, t2, nuv, nw, gridtype, levi, &
!                                lli, ll, sp_ll, ggi, gg, sp_gg, shi, sh, lnsp_sh, &
!                                tmi, status )
!        IF_NOTOK_RETURN(status=1)
!
!        ! init pointer:
!        call pa_Init( ll2 ) ; call pa_Init( gg2 ) ; call pa_Init( sh2 )
!
!        ! read second field:
!        call mf_ReadRecord_1( mf, 'nsss', t1, t2, nuv, nw, gridtype, levi, &
!                                lli, ll2, sp_ll, ggi, gg2, sp_gg, shi, sh2, lnsp_sh, &
!                                tmi, status )
!         IF_NOTOK_RETURN(status=1)
!
!        ! process:
!        select case ( gridtype )
!          case ( 'll' ) ; ll = sqrt( ll**2 + ll2**2 )
!          case ( 'gg' ) ; gg = sqrt( gg**2 + gg2**2 )
!          case default
!            write (gol,'("unsupported gridtype for sstr :",a)') gridtype; call goErr
!            write (gol,'("in ",a)') rname; call goErr; status=1; return
!        end select
!
!        ! clear pointers:
!        call pa_Done( ll2 ) ; call pa_Done( gg2 ) ; call pa_Done( sh2 )

                ! read other time average field to setup output grid:
                call mf_ReadRecord_1( mf, 'slhf', t1, t2, nuv, nw, gridtype, levi, &
                                        lli, ll, sp_ll, ggi, gg, sp_gg, shi, sh, lnsp_sh, &
                                        tmi, status )
                IF_NOTOK_RETURN(status=1)

                ! dummy ...
                write (*,'("          WARNING - no surface stress in rasch-bundy set; use constant value ...")')
                gg = 0.1  ! N/m2

      !
      ! sea surface temperature :  skin temperture *  sea-mask
      !
      case ( 'sst' )

        ! read skin temperature:
        call mf_ReadRecord_1( mf, 'skt', t1, t2, nuv, nw, gridtype, levi, &
                                lli, ll, sp_ll, ggi, gg, sp_gg, shi, sh, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! init pointers:
        call pa_Init( ll2 ) ; call pa_Init( gg2 ) ; call pa_Init( sh2 )

        ! read land-sea mask:
        call mf_ReadRecord_1( mf, 'lsm', t1, t2, nuv, nw, gridtype, levi, &
                                lli, ll2, sp_ll, ggi, gg2, sp_gg, shi, sh2, lnsp_sh, &
                                tmi, status )
         IF_NOTOK_RETURN(status=1)

        ! process:
        select case ( gridtype )
          case ( 'll' ) ; ll = ll * max( 0.0, 1.0 - ll2 )
          case ( 'gg' ) ; gg = gg * max( 0.0, 1.0 - gg2 )
          case default
            write (gol,'("unsupported gridtype for sst :",a)') gridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! clear pointers:
        call pa_Done( ll2 ) ; call pa_Done( gg2 ) ; call pa_Done( sh2 )

      !
      ! not available ...
      !
      case ( 'sr' )

        ! read land-sea mask to setup grids:
        call mf_ReadRecord_1( mf, 'lsm', t1, t2, nuv, nw, gridtype, levi, &
                                lli, ll, sp_ll, ggi, gg, sp_gg, shi, sh, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! dummy ...
        write (*,'("          WARNING - no surface roughness in rasch-bundy set; use constant value ...")')
        gg = 0.001   ! m

      !
      ! no specials ...
      !
      case default

        call mf_ReadRecord_1( mf, paramkey, t1, t2, nuv, nw, gridtype, levi, &
                                lli, ll, sp_ll, ggi, gg, sp_gg, shi, sh, lnsp_sh, &
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
    use GO        , only : TDate, IncrDate, operator(<), operator(+), operator(-), rTotal, wrtgol
    use Grid      , only : TLevelInfo
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use tmm_info  , only : TMeteoInfo

    ! --- in/out -------------------------------

    type(TMeteoFile_ncep_gfs), intent(inout)  ::  mf
    character(len=*), intent(in)              ::  paramkey
    type(TDate), intent(in)                   ::  t1, t2
    character(len=1), intent(in)              ::  nuv
    character(len=1), intent(in)              ::  nw
    character(len=2), intent(out)             ::  gridtype
    type(TLevelInfo), intent(out)             ::  levi
    type(TllGridInfo), intent(inout)          ::  lli
    real, pointer                             ::  ll(:,:,:)
    real, pointer                             ::  sp_ll(:,:)
    type(TggGridInfo), intent(inout)          ::  ggi
    real, pointer                             ::  gg(:,:)
    real, pointer                             ::  sp_gg(:)
    type(TshGridInfo), intent(inout)          ::  shi
    complex, pointer                          ::  sh(:,:)
    complex, pointer                          ::  lnsp_sh(:)
    type(TMeteoInfo), intent(out)             ::  tmi
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord_1'

    ! --- local -------------------------------

    type(TDate)           ::  ta
    real, pointer         ::  ll1(:,:,:)
    real, pointer         ::  gg1(:,:)
    complex, pointer      ::  sh1(:,:)
    real                  ::  dt, dt1
    integer               ::  aver_hr

    ! --- begin ---------------------------------

    ! accumulated field ?
    select case ( paramkey )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! accumulated fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'lsp', 'cp', 'sf', 'sshf', 'slhf', 'ssr', 'ewss', 'nsss' )

        ! should be a time interval ...
        if ( .not. (t1 < t2) ) then
          write (gol,'("accumulated fields requires time interval:")'); call goErr
          write (gol,'("  paramkey : ",a)') paramkey; call goErr
          call wrtgol( '  t1       : ', t1 ); call goErr
          call wrtgol( '  t2       : ', t2 ); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if

        ! fields are averaged over intervals of 3 or 6 hour; set begin of averaging interval:
        call GetTime( mf%treskey, mf%tref, t1, t2, status, taver=ta )

        ! read field accumulated over [ta,t2] :
        call mf_ReadRecord_2( mf, paramkey, ta, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! substract [tf1,t1] if necessary:
        if ( ta < t1 ) then

          ! init pointer:
          call pa_Init( ll1 )
          call pa_Init( gg1 )
          call pa_Init( sh1 )

          ! read field accumulated over [ta,t1] :
          call mf_ReadRecord_2( mf, paramkey, ta, t1, nuv, nw, &
                                  gridtype, levi, &
                                  lli, ll1, sp_ll, &
                                  ggi, gg1, sp_gg, &
                                  shi, sh1, lnsp_sh, &
                                  tmi, status )
          IF_NOTOK_RETURN(status=1)

          ! ll contains aver over [ta,t2], ll1 contains aver over [ta,t1]
          dt  = rTotal( t2 - ta, 'sec' )
          dt1 = rTotal( t1 - ta, 'sec' )

          ! substract accumulations, return average:
          select case ( gridtype )
            case ( 'll' ) ; ll = ( ll*dt - ll1*dt1 ) / rTotal( t2 - t1, 'sec' )
            case ( 'gg' ) ; gg = ( gg*dt - gg1*dt1 ) / rTotal( t2 - t1, 'sec' )
            case ( 'sh' ) ; sh = ( sh*dt - sh1*dt1 ) / rTotal( t2 - t1, 'sec' )
            case default
              write (gol,'("unsupported gridtype for substract :",a)') gridtype; call goErr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select

          ! clear pointers:
          call pa_Done( ll1 )
          call pa_Done( gg1 )
          call pa_Done( sh1 )

        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! instantaneous fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        call mf_ReadRecord_2( mf, paramkey, t1, t2, nuv, nw, &
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


  ! Return a field given parameter name, time, etc.
  ! Only one of grid types is filled!

  subroutine mf_ReadRecord_2( mf, paramkey, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use PArray    , only : pa_Init, pa_Done, pa_SetShape
    use binas     , only : grav
    use GO        , only : TDate, Get
    use Grid      , only : TLevelInfo, Init, Done
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use Grid      , only : Interpol
    use tmm_info  , only : TMeteoInfo, Init, AddHistory

    ! --- in/out -------------------------------

    type(TMeteoFile_ncep_gfs), intent(inout)  ::  mf
    character(len=*), intent(in)              ::  paramkey
    type(TDate), intent(in)                   ::  t1, t2
    character(len=1), intent(in)              ::  nuv
    character(len=1), intent(in)              ::  nw
    character(len=2), intent(out)             ::  gridtype
    type(TLevelInfo), intent(out)             ::  levi
    type(TllGridInfo), intent(inout)          ::  lli
    real, pointer                             ::  ll(:,:,:)
    real, pointer                             ::  sp_ll(:,:)
    type(TggGridInfo), intent(inout)          ::  ggi
    real, pointer                             ::  gg(:,:)
    real, pointer                             ::  sp_gg(:)
    type(TshGridInfo), intent(inout)          ::  shi
    complex, pointer                          ::  sh(:,:)
    complex, pointer                          ::  lnsp_sh(:)
    type(TMeteoInfo), intent(out)             ::  tmi
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord_2'

    ! --- local -------------------------------

    type(TDate)          ::  tfile
    integer              ::  year, month, day, hour
    integer              ::  fcst_hr

    complex, pointer      ::  sh1(:)
    real, pointer         ::  gg1(:)

    character(len=64)     ::  key

    ! --- begin ---------------------------------

    ! no fluxes through boundaries ...
    if ( nuv /= 'n' ) then
      write (gol,'("unsupported nuv key : ",a)') nuv; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! init info; example of history:  model==msc;sh==159;nlev==60
    call Init( tmi, paramkey, 'unknown', status )
    call AddHistory( tmi, 'model==ncep gfs', status )
    call AddHistory( tmi, 'archive==ncar', status )

    !
    ! write filename
    !

    ! extract times used in filename:
    call GetTime( mf%treskey, mf%tref, t1, t2, status, tfile=tfile, fcst_hr=fcst_hr )
    IF_NOTOK_RETURN(status=1)

    ! write filename:
    !   <dir>/20040614/20040614.00.00.SF
    !   <dir>/20040614/20040614.00.00.SFLUXGrbF
    !
    call Get( tfile, year=year, month=month, day=day, hour=hour )
    write (mf%fname,'(i4.4,2i2.2,"-",i4.4,2i2.2,".",i2.2,".",i2.2,".",a)') &
             year, month, day, year, month, day, hour, fcst_hr, trim(mf%ext)

    !
    ! read from cray binary or grib file
    !

    select case ( paramkey )

      !
      ! * 2d spectral fields
      !
      case ( 'oro', 'LNSP' )

        ! output is spectral field
        gridtype = 'sh'

        ! read 2d field
        call Read_Bin_2d( trim(mf%dir)//trim(mf%fname), paramkey, t1, shi, sh1, status )
        IF_NOTOK_RETURN(status=1)

        ! copy from rank1 array:
        call pa_SetShape( sh, shi%np, 1 )
        sh(:,1) = sh1
        call pa_Done( sh1 )

        ! unit conversion:
        select case ( paramkey )
          !
          ! factor for conversion from cbar to Pa :
          !   [Pa] = [cbar] * 1e-2 [bar/cbar] * 1e5 [Pa/bar] = [mbar] * 1e3
          ! add to first complex coeff:
          !   sp * fac  =  exp( lnsp + nlog(fac) )
          !             =  exp( {sum_i=1,n c_i p_i} + nlog(fac) )
          !             =  exp( c_1 + {sum_i=2,n c_i p_i} + nlog(fac) )
          case ( 'LNSP' )
            sh(1,1) = sh(1,1) + cmplx(log(1.0e3),0.0)
          !
          ! ncep oro in in [m], should be [m][g] = [m m/s2]
          case ( 'oro' )
            sh(:,1) = sh(:,1) * grav
        end select

        ! dummy levels
        call Init( levi, 1, (/0.0,0.0/), (/1.0,0.0/), status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        call AddHistory( tmi, 'file=='//trim(mf%fname), status )
        write (key,'("sh==",i4.4)') shi%T
        call AddHistory( tmi, trim(key), status )

      !
      ! * 3d spectral fields
      !
      case ( 'VO', 'D', 'Tv', 'Q' )

        ! output is spectral field
        gridtype = 'sh'

        ! read 3d field
        call Read_Bin_3d( trim(mf%dir)//trim(mf%fname), paramkey, t1, levi, shi, sh, status )
        IF_NOTOK_RETURN(status=1)

        ! unit conversion:
        select case ( paramkey )
          ! For some reason, the u/v/w from VO/D needs a factor -1 ...
          ! The minus is probably caused by the upwards coordinate system of ncep
          ! instead of the downward from ecmwf.
          case ( 'VO', 'D' )
            sh = - sh
        end select

        ! read lnsp:
        call Read_Bin_2d( trim(mf%dir)//trim(mf%fname), 'LNSP', t1, shi, lnsp_sh, status )
        IF_NOTOK_RETURN(status=1)
        ! unit conversion:
        lnsp_sh(1) = lnsp_sh(1) + cmplx(log(1.0e3),0.0)

        ! info ...
        call AddHistory( tmi, 'file=='//trim(mf%fname), status )
        write (key,'("sh==",i4.4)') shi%T
        call AddHistory( tmi, trim(key), status )

      !
      ! * 2d surface fields
      !
      case ( 'lsm', 'albedo', 'sr', 'sps', 'ci', 'skt', 'u10m', 'v10m', 'slhf', 'sshf', 'ewss', 'nsss', 'lsp', 'cp' )

        ! output is gaussian grid
        gridtype = 'gg'

        ! read 2d field
        call Read_Grib_2d( mf, paramkey, t1, t2, ggi, gg1, status )
        IF_NOTOK_RETURN(status=1)

        ! copy from rank1 array:
        call pa_SetShape( gg, ggi%np, 1 )
        gg(:,1) = gg1
        call pa_Done( gg1 )

        ! unit conversion:
        select case ( paramkey )
          !
          ! for some probably historical reaseon, TM expects land/sea mask in % ...
          case ( 'lsm' )
            gg(:,1) = gg(:,1) * 100.0     ! 0-1  ->  %
          !
          ! fluxes downward (ecmwf direction) instead of upward (ncep direction)
          case ( 'slhf', 'sshf' )
            gg(:,1) = - gg(:,1)
          !
          ! kg water / m2  / s   ->   m water / s
          ! With density of 998 kg water / m3  :    kg/m2/s / (kg/m3) = m/s
          case ( 'lsp', 'cp' )
            gg(:,1) = gg(:,1) / 998.0     ! m water / s
        end select

        ! dummy levels
        call Init( levi, 1, (/0.0,0.0/), (/1.0,0.0/), status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        call AddHistory( tmi, 'file=='//trim(mf%fname), status )
        write (key,'("gg==",i4.4)') ggi%N
        call AddHistory( tmi, trim(key), status )

      !
      ! * ????
      !
      case default
        write (gol,'("unsupported paramkey : ",a)') paramkey; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select


    !
    ! ~~~ end
    !

    ! ok
    status = 0

  end subroutine mf_ReadRecord_2


  ! ***


  subroutine Check_Bin_Time( ncb, t, status )

    use go      , only : TDate, NewDate, IncrDate, operator(+), operator(/=), IsAnyDate, wrtgol
    use file_ncb, only : TNcepCrayBin

    ! --- in/out -----------------------------------

    type(TNcepCrayBin), intent(in)            ::  ncb
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Check_Bin_Time'

    ! --- local -----------------------------------

    type(TDate)          ::  tf

    ! --- begin -----------------------------------

    ! trap any date ..
    if ( IsAnyDate(t) ) then
      ! ok
      status = 0; return
    end if

    ! time in bin file:
    tf = NewDate( time4=ncb%idate ) + IncrDate( hour=ncb%fcst_hr )

    ! check:
    if ( tf /= t ) then
      write (gol,'("wrong time in binary file:")'); call goErr
      call wrtgol( '  file : ', tf ); call goErr
      call wrtgol( '  t    : ', t  ); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine Check_Bin_Time


  ! ***


  subroutine Read_Bin_2d( fname, paramkey, t1, shi, sh, status )

    use parray  , only : pa_SetShape
    use go      , only : TDate
    use grid    , only : TshGridInfo, Init
    use file_ncb, only : TNcepCrayBin, Init, Done, ReadRecord

    ! --- in/out -----------------------------------

    character(len=*), intent(in)              ::  fname
    character(len=*), intent(in)              ::  paramkey
    type(TDate), intent(in)                   ::  t1
    type(TshGridInfo), intent(out)            ::  shi
    complex, pointer                          ::  sh(:)
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Read_Bin_2d'

    ! --- local -----------------------------------

    type(TNcepCrayBin)   ::  ncb

    ! --- begin -----------------------------------

    ! open file
    call Init( ncb, fname, status )
    IF_NOTOK_RETURN(status=1)

    ! check time:
    call Check_Bin_Time( ncb, t1, status )
    IF_NOTOK_RETURN(status=1)

    ! init output grid definition:
    call Init( shi, ncb%shT, status )
    IF_NOTOK_RETURN(status=1)

    ! allocate output grid:
    call pa_SetShape( sh, shi%np )

    ! read record:
    call ReadRecord( ncb, paramkey, sh, status )
    IF_NOTOK_RETURN(status=1)

    ! convert from NCEP spectral coeff to ECMWF spectral coeff:
    sh = sh / sqrt(2.0)

    ! close file:
    call Done( ncb, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Read_Bin_2d


  ! ***


  subroutine Read_Bin_3d( fname, paramkey, t1, levi, shi, sh, status )

    use parray  , only : pa_SetShape
    use go      , only : TDate
    use grid    , only : TshGridInfo, Init
    use grid    , only : TLevelInfo, Init
    use file_ncb, only : TNcepCrayBin, Init, Done, ReadRecord

    ! --- in/out -----------------------------------

    character(len=*), intent(in)              ::  fname
    character(len=*), intent(in)              ::  paramkey
    type(TDate), intent(in)                   ::  t1
    type(TLevelInfo), intent(out)             ::  levi
    type(TshGridInfo), intent(out)            ::  shi
    complex, pointer                          ::  sh(:,:)
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Read_Bin_2d'

    ! --- local -----------------------------------

    type(TNcepCrayBin)   ::  ncb

    ! --- begin -----------------------------------

    ! open file
    call Init( ncb, fname, status )
    IF_NOTOK_RETURN(status=1)

    ! check time:
    call Check_Bin_Time( ncb, t1, status )
    IF_NOTOK_RETURN(status=1)

    ! init levels
    call Init( levi, ncb%nlev, ncb%sigma_half*0.0, ncb%sigma_half, status )
    IF_NOTOK_RETURN(status=1)

    ! init output grid definition:
    call Init( shi, ncb%shT, status )
    IF_NOTOK_RETURN(status=1)

    ! allocate output grid:
    call pa_SetShape( sh, shi%np, ncb%nlev )

    ! read record:
    call ReadRecord( ncb, paramkey, sh, status )
    IF_NOTOK_RETURN(status=1)

    ! convert from NCEP spectral coeff to ECMWF spectral coeff:
    sh = sh / sqrt(2.0)

    ! close file:
    call Done( ncb, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Read_Bin_3d


  ! ***


  subroutine Read_Grib_2d( mf, paramkey, t1, t2, ggi, gg, status )

    use parray  , only : pa_SetShape
    use go      , only : TDate, iTotal, Get, operator(==), operator(-), IsAnyDate
    use grid    , only : TggGridInfo, Init
    use file_ncg, only : TNcepGrib, Init, Done, ReadRecord, Get
    use file_ncg, only : CheckRecord

    ! --- in/out -----------------------------------

    type(TMeteoFile_ncep_gfs), intent(inout)  ::  mf
    character(len=*), intent(in)              ::  paramkey
    type(TDate), intent(in)                   ::  t1, t2
    type(TggGridInfo), intent(out)            ::  ggi
    real, pointer                             ::  gg(:)
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Read_Grib_2d'

    ! --- local -----------------------------------

    type(TNcepGrib)   ::  ncg
    character(len=5)  ::  param
    type(TDate)       ::  tref
    integer           ::  reftime(5)
    integer           ::  timerange(4), dh, dha
    character(len=4)  ::  levtype
    integer           ::  level
    integer           ::  ggN

          integer ::  irec

    ! --- begin -----------------------------------

    !! open file
    !call Init( ncg, trim(mf%dir)//trim(mf%fname), status )
    !IF_NOTOK_RETURN(status=1)

         ! For some reason, opening a new with the same file unit
         ! (while the old one is closed) in the same program gives errors  ...
         ! Thus, new file unit for each new file name ...
         if ( trim(mf%dir)//trim(mf%fname) /= adhoc_fname ) then
           adhoc_fu    = adhoc_fu + 1
           adhoc_fname = trim(mf%dir)//trim(mf%fname)
         end if
         call Init( ncg, adhoc_fu, trim(mf%dir)//trim(mf%fname), status )
         IF_NOTOK_RETURN(status=1)

    ! set grib param key
    select case ( paramkey )
      case ( 'lsm'    ) ; param = 'LAND'   ; levtype = 'SFC'  ; level = 0    ; irec=32
      case ( 'sps'    ) ; param = 'PRES'   ; levtype = 'SFC'  ; level = 0    ; irec=38
      case ( 'ci'     ) ; param = 'ICEC'   ; levtype = 'SFC'  ; level = 0    ; irec=33
      case ( 'albedo' ) ; param = 'ALBDO'  ; levtype = 'SFC'  ; level = 0    ; irec=48
      case ( 'skt'    ) ; param = 'TMP'    ; levtype = 'SFC'  ; level = 0    ; irec=5
      case ( 'u10m'   ) ; param = 'UGRD'   ; levtype = 'HTGL' ; level = 10   ; irec=34
      case ( 'v10m'   ) ; param = 'VGRD'   ; levtype = 'HTGL' ; level = 10   ; irec=35
      case ( 'slhf'   ) ; param = 'LHTFL'  ; levtype = 'SFC'  ; level = 0    ; irec=4
      case ( 'sshf'   ) ; param = 'SHTFL'  ; levtype = 'SFC'  ; level = 0    ; irec=3
      !case ( 'ewss'   ) ; param = ''
      !case ( 'nsss'   ) ; param = ''
      case ( 'lsp'    ) ; param = 'PRATE'  ; levtype = 'SFC'  ; level = 0    ; irec=29
      case ( 'cp'     ) ; param = 'CPRAT'  ; levtype = 'SFC'  ; level = 0    ; irec=30
      !case ( 'blh'    ) ; param = 'HPBL'
      case default
        write (gol,'("unsupported paramkey ",a)') paramkey; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! no time at all, instant time or time average ?
    if ( IsAnyDate(t1) .and. IsAnyDate(t2) ) then
      ! no search ...
      reftime   = -1
      timerange = -1
    else if ( t1 == t2 ) then
      ! extract ref time and time step:
      call GetTime( mf%treskey, t1, t1, t2, status, tfile=tref, fcst_hr=dh )
      IF_NOTOK_RETURN(status=1)
      ! set reftime and timerange arrays:
      call Get( tref, time5=reftime )
      timerange = (/1,dh,0,10/)      ! hours, P1, P2, valid for reftime + P1
    else
      ! extract ref time and time step:
      call GetTime( mf%treskey, t1, t1, t2, status, tfile=tref, fcst_hr=dh, aver_hr=dha )
      IF_NOTOK_RETURN(status=1)
      ! set reftime and timerange arrays:
      call Get( tref, time5=reftime )
      timerange = (/1,dha,dh,3/)      ! hours, P1, P2, aver over reftime + [P1,P2]
    end if

    ! >>> for some reason, this failes ... problems with w3 library ?
    !! search requested record:
    !call ReadRecord( ncg, status, param=param, reftime=reftime, timerange=timerange, &
    !                   levtype=levtype, level=level )
    !IF_NOTOK_RETURN(status=1)
    !<<<

         ! >>> not so nice, but effective ...

         ! read record given specified message number:
         call ReadRecord( ncg, status, irec=irec )
         IF_NOTOK_RETURN(status=1)

         ! check contents:
         call CheckRecord( ncg, status, param=param, reftime=reftime, timerange=timerange, &
                            levtype=levtype, level=level )
         IF_NOTOK_RETURN(status=1)

         ! <<<

    ! extract Gaussian grid number :
    call Get( ncg, status, ggN=ggN )
    IF_NOTOK_RETURN(status=1)

    ! init output grid, no reduced grid:
    call Init( ggi, ggN, .false., status )
    IF_NOTOK_RETURN(status=1)

    ! allocate output grid:
    call pa_SetShape( gg, ggi%np )

    ! extract field:
    call Get( ncg, status, gg=gg )
    IF_NOTOK_RETURN(status=1)

    ! close file:
    call Done( ncg, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Read_Grib_2d


end module tmm_mf_ncep_gfs


