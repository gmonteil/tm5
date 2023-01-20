!###############################################################################
!
! Input/output of meteofiles : MSC meteo version.
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

module tmm_mf_msc

  use GO, only : gol, goErr, goPr
  use GO, only : TDate
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN, DUMMY_STR_LEN

  implicit none

  ! --- in/out ---------------------------

  private

  public   ::  TMeteoFile_msc
  public   ::  Init, Done
  public   ::  Get
  public   ::  ReadRecord


  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf_msc'


  ! --- types -------------------------------

  type TMeteoFile_msc
    ! file name:
    character(len=MAX_FILENAME_LEN)         ::  fname
    ! current time range:
    type(TDate)                ::  t1, t2
    ! params stored in file:
    character(len=MAX_RCKEY_LEN)         ::  paramkeys
    ! extra ..
    character(len=10)          ::  sp_unit
  end type TMeteoFile_msc


  ! --- interfaces --------------------------

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


  ! ==================================================================

  ! Initialise msc file:
  !  o setup filename for given:
  !     * input dir
  !     * archive keys (see tmm_mf_hdf for inspiration, or tmm_mf_grib/Init2)
  !     * parameter name ('LNSP', 'VO', etc)
  !     * times (now only t1 is used)
  !

  subroutine mf_Init( mf, dir, archivekeys, paramkey, &
                               tday, t1, t2, status )

    use GO, only : TDate, NewDate, Get
    use GO, only : goVarValue

    ! --- in/out -----------------------

    type(TMeteoFile_msc), intent(out)         ::  mf
    character(len=*), intent(in)        ::  dir
    character(len=*), intent(in)        ::  archivekeys
    character(len=*), intent(in)        ::  paramkey
    type(TDate), intent(in)             ::  tday, t1, t2
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Init'

    ! --- local ---------------------

    logical              ::  exist
    character(len=20)    ::  msc_mdir, msc_tres, msc_type
    integer               ::  ccyy, mm, dd, mm2, dd2
    integer                  ::  year1, month1, day1, hour1
    integer                  ::  year2, month2, day2

    ! --- begin ------------------------

    !
    ! extract fields from archivekey :
    !   nlev=71;sh=47;mdir=cmam;tres=_1dag_6hrly
    !
    msc_mdir   = 'cmam'
    call goVarValue( archivekeys, ';', 'mdir', '=', msc_mdir, status )
    if (status>0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
    !
    msc_type   = 'unknown'
    call goVarValue( archivekeys, ';', 'type', '=', msc_type, status )
    if (status>0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
    !
    msc_tres   = '_20000101.dat'
    call goVarValue( archivekeys, ';', 'tres', '=', msc_tres, status )
    if (status>0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
    !
    mf%sp_unit = 'unknown'
    call goVarValue( archivekeys, ';', 'sp_unit', '=', mf%sp_unit, status )
    if (status>0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! extract time values:
    call Get( t1, year=ccyy, month=mm, day=dd )

    ! fill filename (should be function of t1 etc):
    write (mf%fname,'(a,"/",a,"-",a,"-",a,"_",i4.4,i2.2,i2.2,a,".dat")') &
    trim(dir), trim(msc_mdir), trim(msc_type),'all', ccyy, mm, dd, trim(msc_tres)

    ! file exist ?
    inquire( file=mf%fname, exist=exist )
    if ( .not. exist ) then
      write (gol,'("msc file not found:")'); call goErr
      write (gol,'("  ",a)') trim(mf%fname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! file contains data for the following time range:
    mf%t1 = NewDate( 2000, 01, 01, 00, 00, 00 )
    mf%t2 = NewDate( 2000, 01, 01, 18, 00, 00 )


    ! params stored in file; specify a 'minus' seperated list:
    mf%paramkeys = '-VO-D-T-LNSP-'

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    ! --- in/out -----------------------

    type(TMeteoFile_msc), intent(inout)    ::  mf
    integer, intent(out)             ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Done'

    ! --- begin ------------------------

    ! deallocate temporary arrays etc

    ! ok
    status = 0

  end subroutine mf_Done


  ! ***


  subroutine mf_Get( mf, status, trange1, trange2, paramkeys )

    use GO, only : TDate

    ! --- in/out ----------------------------

    type(TMeteoFile_msc), intent(in)    ::  mf
    integer, intent(out)                ::  status

    type(TDate), intent(out), optional            ::  trange1, trange2
    character(len=*), intent(out), optional       ::  paramkeys

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Get'

    ! --- local --------------------------------

    ! --- begin --------------------------------

    ! time range:
    if ( present(trange1) ) trange1 = mf%t1
    if ( present(trange2) ) trange2 = mf%t2

    ! params:
    if ( present(paramkeys) ) paramkeys = mf%paramkeys

    ! ok
    status = 0

  end subroutine mf_Get


  ! ***


  ! Return a field given parameter name, time, etc.
  ! Only one of grid types is filled!

  subroutine mf_ReadRecord( mf, paramkey, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use GO        , only : TDate, wrtgol, Get
    use GO        , only : goGetFU
    use Grid      , only : TLevelInfo, Init
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use PArray    , only : pa_Init, pa_Done, pa_SetShape
    use tmm_info  , only : TMeteoInfo, Init, AddHistory

    ! --- in/out -------------------------------

    type(TMeteoFile_msc), intent(inout)   ::  mf
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

    ! --- const --------------------------------

    ! number of layers
    integer, parameter   ::  nlev3d = 71

    ! spectral resolution
    integer, parameter   ::  shT = 47

    ! number of lines in record (one layer)
    integer, parameter   ::  nline = 393

    ! --- local -------------------------------

    integer             ::  fu
    integer             ::  year1, month1, day1, hour1, itrec
    integer             ::  nskip
    integer             ::  iline, iline_tot
    character(len=DUMMY_STR_LEN)  ::  line
    character(len=10)   ::  htime
    character(len=4)   ::   hname
    integer             ::  ilev, nlev
    character(len=64)   ::  key
    logical             ::  read_lnsp

    ! --- begin ---------------------------------

    !write (*,'("grib read record: paramkey : ",a)') paramkey
    !call wrtgol( 'grib read record: t1       : ', t1 )
    !call wrtgol( 'grib read record: t2       : ', t2 )

    ! no fluxes through boundaries ...
    if ( nuv /= 'n' ) then
      write (gol,'("unsupported nuv key : ",a)') nuv; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! extract year, month, day and hour:
    call Get( t1, year=year1, month=month1, day=day1, hour=hour1 )

    ! time records for 00/06/12/18 only:
    select case ( hour1 )
       case ( 00 )
         itrec = 1
    write (htime,'(i4.4,i2.2,i2.2,i2.2)') year1,month1,day1,hour1
       case ( 06 )
         itrec = 2
    write (htime,'(i4.4,i2.2,i2.2,i2.2)') year1,month1,day1,hour1
       case ( 12 )
         itrec = 3
    write (htime,'(i4.4,i2.2,i2.2,i2.2)') year1,month1,day1,hour1
       case ( 18 )
         itrec = 4
    write (htime,'(i4.4,i2.2,i2.2,i2.2)') year1,month1,day1,hour1
    end select

    ! example of history:
    !   model==msc;sh==159;nlev==60
    !
    call Init( tmi, paramkey, 'unknown', status )
    call AddHistory( tmi, 'model==msc', status )
    !
    write (key,'("nlev==",i3.3)') nlev
    call AddHistory( tmi, trim(key), status )
    !
    write (key,'("sh==",i4.4)') shT
    call AddHistory( tmi, trim(key), status )

    ! routine returns sh grid:
    gridtype = 'sh'

    ! intialize spherical harmonic field info:
    call Init( shi, shT, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ~~~ 3d field
    !

    ! select free file unit
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open text file
    open( fu, file=trim(mf%fname), form='formatted', iostat=status )
    IF_NOTOK_RETURN(status=1)

    ! by default: 3d field, read lnsp
    nlev = nlev3d
    read_lnsp = .true.

    ! skip previous time records:
    nskip = nline * ( nlev3d * 3 + 1 ) * (itrec-1)

    ! number of lines to skip;
    ! name of field in header line:
    select case ( paramkey )
      case ( 'VO' )
        nskip = nskip + 0
        hname = 'VORT'
      case ( 'D' )
        nskip = nskip + nline * nlev3d
        hname = ' DIV'
      case ( 'T' )
        nskip = nskip + nline * nlev3d * 2
        hname = 'TEMP'
      case ( 'LNSP' )
        nskip = nskip + nline * nlev3d * 3
        hname = 'LNSP'
        nlev = 1
        read_lnsp = .false.
      case default
        write (gol,'("unsupported paramkey `",a,"` for nskip")') trim(paramkey); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! allocate output:
    call pa_SetShape( sh, shi%np, nlev )

    ! no lines read yet:
    iline_tot = 0

    ! skip first lines
    do iline = 1, nskip
      iline_tot = iline_tot + 1
      read (fu,'(a)',iostat=status) line
      if ( status /= 0 ) then
        write (gol,'("while reading line : ",i10)') iline_tot; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
    end do

    ! read spectral coeff
    do ilev = 1, nlev
      ! skip header line
      iline_tot = iline_tot + 1
      read (fu,'(a)',iostat=status) line
      if ( status /= 0 ) then
        write (gol,'("while reading line : ",i10)') iline_tot; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! check: header line should contain correct variable name:
      if ( index(line,hname) < 1 ) then
        write (gol,'("record variable not ok:")'); call goErr
        write (gol,'("  line nr  : ",i6)') iline_tot; call goErr
        write (gol,'("  header   : ",a)') trim(line); call goErr
        write (gol,'("  searched : ",a)') htime; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! check: header line should contain correct time :
      if ( index(line,htime) < 1 ) then
        write (gol,'("record time not ok:")'); call goErr
        write (gol,'("  line nr  : ",i10)') iline_tot; call goErr
        write (gol,'("  header   : ",a)') trim(line); call goErr
        write (gol,'("  searched : ",a)') htime; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! read coeff
      do iline = 1, nline-1
        iline_tot = iline_tot + 1
        read (fu,'(6e22.15)',iostat=status) sh(iline*3-2:iline*3,ilev)
        if ( status /= 0 ) then
          write (gol,'("while reading line : ",i10)') iline_tot; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if
      end do
    end do

    ! close
    close( fu, iostat=status )
    IF_NOTOK_RETURN(status=1)

    ! unit conversion:
    select case ( paramkey )
      case ( 'LNSP' )
        !   sp * fac  =  exp( lnsp + ln(fac) )
        !             =  exp( {sum_i=1,n c_i p_i} + ln(fac) )
        !             =  exp( c_1 + {sum_i=2,n c_i p_i} + ln(fac) )
        select case ( mf%sp_unit )
          case ( 'hPa' )
            ! add ln(fac) to first complex coeff (only level 1 is in use):
            do ilev = 1, nlev
              sh(1,ilev) = sh(1,ilev) + cmplx(log(100.0),0.0)
            end do
          case default
            write (gol,'("unsupported sp unit `",a,"`")') trim(mf%sp_unit); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
    end select

    ! levels
    select case ( nlev )
      case ( 1 )
        call Init( levi, nlev, (/0.0,0.0/), (/0.0,0.0/), status )
        IF_NOTOK_RETURN(status=1)
      case ( 71 )
        call Init( levi, 'msc71', status )
        IF_NOTOK_RETURN(status=1)
        case default
        write (gol,'("unsupported nlev `",i4,"` for levi")') nlev; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select


    !
    ! ~~~ surface pressure
    !

    if ( read_lnsp ) then

      ! allocate output:
      call pa_SetShape( lnsp_sh, shi%np )

      ! open text file
      open( fu, file=trim(mf%fname), form='formatted', iostat=status )
      IF_NOTOK_RETURN(status=1)

      ! skip previous time records:
      nskip = nline * ( nlev3d * 3 + 1 ) * (itrec-1)

      ! skip 3D fields
      nskip = nskip + nline * nlev3d * 3
      hname = 'LNSP'

      ! no lines read yet:
      iline_tot = 0

      ! skip first lines
      do iline = 1, nskip
        iline_tot = iline_tot + 1
        read (fu,'(a)',iostat=status) line
        if ( status /= 0 ) then
          write (gol,'("while reading line : ",i6)') iline_tot; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if
      end do

      ! skip header line
      iline_tot = iline_tot + 1
      read (fu,'(a)',iostat=status) line
      if ( status /= 0 ) then
        write (gol,'("while reading line : ",i6)') iline_tot; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! check: header line should contain correct variable name:
      if ( index(line,hname) < 1 ) then
        write (gol,'("record variable not ok:")'); call goErr
        write (gol,'("  line nr  : ",i6)') iline_tot; call goErr
        write (gol,'("  header   : ",a)') trim(line); call goErr
        write (gol,'("  searched : ",a)') htime; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! check: header line should contain correct time :
      if ( index(line,htime) < 1 ) then
        write (gol,'("record time not ok:")'); call goErr
        write (gol,'("  line nr  : ",i6)') iline_tot; call goErr
        write (gol,'("  header   : ",a)') trim(line); call goErr
        write (gol,'("  searched : ",a)') htime; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

      ! read spectral coeff
      do iline = 1, nline-1
        iline_tot = iline_tot + 1
        read (fu,'(6e22.15)',iostat=status) lnsp_sh(iline*3-2:iline*3)
        if ( status /= 0 ) then
          write (gol,'("while reading line : ",i6)') iline_tot; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if
      end do

      ! close
      close( fu, iostat=status )
      IF_NOTOK_RETURN(status=1)

      ! unit conversion:
      !   sp * fac  =  exp( lnsp + ln(fac) )
      !             =  exp( {sum_i=1,n c_i p_i} + ln(fac) )
      !             =  exp( c_1 + {sum_i=2,n c_i p_i} + ln(fac) )
      select case ( mf%sp_unit )
        case ( 'hPa' )
          ! add ln(fac) to first complex coeff:
          lnsp_sh(1) = lnsp_sh(1) + cmplx(log(100.0),0.0)
        case default
          write (gol,'("unsupported sp unit `",a,"`")') trim(mf%sp_unit); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

    end if

    !
    ! ~~~ end
    !

    ! ok
    status = 0

  end subroutine mf_ReadRecord



end module tmm_mf_msc
