!###############################################################################
!
!  Interface to NCEP re-analysis data.
!
!  Database
!  --------
!
!    http://www.cdc.noaa.gov/cdc/reanalysis/reanalysis.shtml
!
!  NetCDF files
!  ------------
!
!  global attributes:
!    Conventions                = "CDC Non-gridded" ;
!    title                      = "4x Daily Spectral Coefficients from the NMC Reanalysis for Natural Log of Pressure at the Surface" ;
!    base_date                  = 2000s, 1s, 1s ;
!    history                    = ... ;
!    description                = "Data is from NMC initialized reanalysis",
!    platform                   = "Model" ;
!    m_fastest                  = "F" ;
!    storage_type               = "Spectral" ;
!    compression_used           = "T" ;
!
!
!  dimensions:
!    num_values = 4032
!    level      = 28
!    num_mean   = 2 ;
!    time       = UNLIMITED ;
!
!  level data sets:
!
!    level             float    array[level]
!    units                      = "sigma_level"
!    long_name                  = "Sigma"
!    positive                   = "down"
!
!  spectral data sets:
!
!    mean             float     array[2*level,time]
!      long_name      string    "First Spectral Coefficients for Natural Log of Pressure at the Surface"
!      units          string    "-"
!      missing_value  float
!
!    add_offset       float     array[level,time]
!      long_name                "Add offset of Spectral Coefficients for Natural Log of Pressure"
!      units                    "-"
!      missing_value
!
!    scale_factor     float     array[level,time]
!      long_name                "Scale Factor of Spectral Coefficients for Natural Log of Pressure"
!      units                    "-"
!      missing_value
!
!    time            double     array[time]
!      units                    "hours since 1-1-1 00:00:0.0" ;
!      long_name                "Time" ;
!      delta_t                  "0000-00-00 06:00:00" ;
!
!    pres            int        array[num_values,level,time]
!      long_name                "Spectral Coefficients for Natural Log of Pressure at the Surface"
!      units                    "natural log of pressure in centibars" ;
!      missing_value
!      precision                = 0s ;
!      least_significant_digit  = 32767s ;
!      trunc_count              = 62s ;
!      trunc_type               = "Triangular" ;
!      var_desc                 = "Pressure",              "P" ;
!      dataset                  = "NMC Reanalysis",         "L" ;
!      level_desc               = "Surface",             "F" ;
!      statistic                = "Individual Obs\n",       "I" ;
!      parent_stat              = "Other\n",              "-" ;
!
!  gaussian or lon/lat data sets:
!
!    short pres(time, lat, lon) ;
!      long_name                = "4xDaily Surface Pressure" ;
!      valid_range              = 40000.f, 115000.f ;
!      actual_range             = 48540.f, 108990.f ;
!      units                    = "Pascals" ;
!      add_offset               = 367650.f ;
!      scale_factor             = 10.f ;
!      missing_value            = 32766s ;
!      precision                = -1s ;
!      least_significant_digit  = -1s ;
!      GRIB_id                  = 1s ;
!      GRIB_name                = "PRES" ;
!      var_desc                 = "Pressure\n", "CC" ;
!      dataset                  = "NMC Reanalysis\n", "L" ;
!
!  Date and time
!  -------------
!
!   Global attribute 'base_date' specifies the first date :
!     base_date = (/2000,01,01/)
!
!   Data set 'time' contains hours since year 1 or there about.
!   Substract the first element from the array to have an offset
!   in hours from midnight at base_date.
!
!  Spectral fields
!  -----------
!
!  From the 'discription':
!
!  For each latitude (going from north pole to south pole).
!  the associated legendre functions are defined:
!  Pbar(m,n,theta) =
!    sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))*sin(theta)**m/(2**n*factorial(n))
!    times the (n+m)th derivative of (x**2-1)**n with respect to x=cos(theta)
!
!    note: theta = 0.5*pi - phi, where phi is latitude and theta is colatitude.
!           Therefore, cos(theta) = sin(phi) and sin(theta) = cos(phi).
!
!  where n is degree (subscript), m is order (superscript),
!  and theta is colatitude in radians.
!  The functions are orthogonal polynomials on the surface of the sphere and
!  are normalized so that the integral of (PBAR(n,m,theta)**2)*sin(theta)
!  on the interval theta=0 to theta=pi equals 1.
!  Note that Pbar(0,0,theta) = sqrt(2)/2
!
!  Truncation:
!    number of data values : 4032 = 2016 complex numbers
!    number of spectral coeff = (T+1)*(T+2)/2 = 2016, thus T=62
!
!  The 'mean' are the first spectral coeff and should be added.
!
!  Packing
!  -------
!
!  From http://www.cdc.noaa.gov/PublicData/faq.html#12 :
!  "Most of the data in our netCDF files are packed. That is to say they have been
!  transformed by a scale factor and an add offset to reduce the storage needed to
!  two bytes per value. When you extract the short integers, you must unpack the
!  data to recover the correct floating point data values. Data files that contain
!  packed data will have a non-zero add offset and/or a scale factor not equal to 1.
!  The transformation is:
!    float_value = (short_value * scale_factor) + add_offset                 "
!
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

module tmm_mf_ncep_cdc

  use GO      , only : gol, goErr, goPr
  use GO      , only : TDate
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

  implicit none

  ! --- in/out ----------------------------

  private

  public  ::  TMeteoFile_ncep_cdc
  public  ::  Init, Done
  public  ::  Get
  public  ::  ReadRecord

  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf_ncep_cdc'

  !--- type ---------------------------------

  type TMeteoFile_ncep_cdc
    ! field collection
    character(len=MAX_RCKEY_LEN)     ::  paramkeys
    type(TDate)            ::  trange(2)
    ! file names
    character(len=MAX_FILENAME_LEN)     ::  dir
    integer                ::  ccyy
    character(len=1)       ::  pathsep, namesep
  end type TMeteoFile_ncep_cdc


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
                               tref, t1, t2, status )

    use GO, only : TDate, Get, NewDate, operator(>)
    use GO, only : goVarValue

    ! --- in/out ----------------------------

    type(TMeteoFile_ncep_cdc), intent(out)  ::  mf
    character(len=*), intent(in)            ::  dir
    character(len=*), intent(in)            ::  archivekeys
    character(len=*), intent(in)            ::  paramkey
    type(TDate), intent(in)                 ::  tref, t1, t2
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Init'

    ! --- local -------------------------------

    type(TDate)         ::  tend

    ! --- begin --------------------------------

    ! store directory:
    mf%dir = dir

    !
    ! extract fields from archivekey :
    !   mdir=spectral
    !
    mf%pathsep  = '/'
    call goVarValue( archivekeys, ';', 'pathsep', '=', mf%pathsep, status )
    if (status>0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
    !
    mf%namesep  = '-'
    call goVarValue( archivekeys, ';', 'namesep', '=', mf%namesep, status )
    if (status>0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! extract year in filename from requested time range:
    call Get( t2, year=mf%ccyy )
    tend = NewDate( time6=(/mf%ccyy,12,31,18,00,00/) )
    if ( t2 > tend ) mf%ccyy = mf%ccyy + 1

    ! files valid for complete year:
    mf%trange(1) = NewDate( mf%ccyy-1, 12, 31, 18, 00, 00, 001 )
    mf%trange(2) = NewDate( mf%ccyy  , 12, 31, 18, 00, 00, 000 )

    ! files contain one param only:
    mf%paramkeys   = '-'//trim(paramkey)//'-'

    ! dummy filename (might be used in error messages)
    mf%filename = 'ncep cdc file'

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    ! --- in/out ------------------------------------

    type(TMeteoFile_ncep_cdc), intent(inout)   ::  mf
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

    type(TMeteoFile_ncep_cdc), intent(in)     ::  mf
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

    type(TMeteoFile_ncep_cdc), intent(inout)  ::  mf
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

    ! combined field ?
    select case ( paramkey )

      !
      ! total surface stress :  sstr^2 =  ewss^2 + nsss^2
      !
      case ( 'sstr' )

        ! read first field:
        call mf_ReadRecord_1( mf, 'ewss', t1, t2, nuv, nw, gridtype, levi, &
                                lli, ll, sp_ll, ggi, gg, sp_gg, shi, sh, lnsp_sh, &
                                tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! init pointer:
        call pa_Init( ll2 ) ; call pa_Init( gg2 ) ; call pa_Init( sh2 )

        ! read second field:
        call mf_ReadRecord_1( mf, 'nsss', t1, t2, nuv, nw, gridtype, levi, &
                                lli, ll2, sp_ll, ggi, gg2, sp_gg, shi, sh2, lnsp_sh, &
                                tmi, status )
         IF_NOTOK_RETURN(status=1)

        ! process:
        select case ( gridtype )
          case ( 'll' ) ; ll = sqrt( ll**2 + ll2**2 )
          case ( 'gg' ) ; gg = sqrt( gg**2 + gg2**2 )
          case default
            write (gol,'("unsupported gridtype for sstr :",a)') gridtype; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! clear pointers:
        call pa_Done( ll2 ) ; call pa_Done( gg2 ) ; call pa_Done( sh2 )

                write (gol,'("WARNING : adhoc constant value for ncep surface stress ...")'); call goPr
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

  ! Return a field given parameter name, time, etc.
  ! Only one of grid types is filled!

  subroutine mf_ReadRecord_1( mf, paramkey, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                ggi, gg, sp_gg, &
                                shi, sh, lnsp_sh, &
                                tmi, status )

    use PArray    , only : pa_Init, pa_Done, pa_SetShape
    use file_hdf  , only : THdfFile, TSds, Init, Done, ReadAttribute, ReadData, GetInfo
    use binas     , only : grav
    use GO        , only : TDate, NewDate, Get, rTotal, IncrDate
    use GO        , only : operator(+), operator(-), operator(/=), operator(<)
    use Grid      , only : TLevelInfo, Init, Done
    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
    use Grid      , only : Interpol
    use tmm_info  , only : TMeteoInfo, Init, AddHistory

    ! --- in/out -------------------------------

    type(TMeteoFile_ncep_cdc), intent(inout)  ::  mf
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

    logical               ::  read_2d, read_3d

    logical               ::  constant
    type(TDate)           ::  tbase
    integer               ::  year, hours
    integer               ::  tstart
    real                  ::  tfactor

    integer               ::  nlev, ilev
    integer               ::  data_dims(1)
    real, allocatable     ::  levels(:)
    character(len=16)     ::  levtype

    integer               ::  shT
    real, pointer         ::  gg1(:)
    complex, pointer      ::  sh1(:)

    character(len=MAX_FILENAME_LEN)    ::  fname
    character(len=64)     ::  mdir
    character(len=16)     ::  fkey
    character(len=16)     ::  gridkey
    logical               ::  exist
    type(THdfFile)        ::  hdf
    type(TSds)            ::  sds

    character(len=16)     ::  sds_name

    character(len=64)     ::  key

    ! --- begin ---------------------------------

    ! no fluxes through boundaries ...
    if ( nuv /= 'n' ) then
      write (gol,'("unsupported nuv key : ",a)') nuv; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! init info; example of history:  model==msc;sh==159;nlev==60
    call Init( tmi, paramkey, 'unknown', status )
    call AddHistory( tmi, 'model==ncep/ncar reanalysis 1', status )
    call AddHistory( tmi, 'archive==CDC netcdf archive', status )


    !
    ! ~~~ setup
    !

    read_2d  = .false.
    read_3d  = .false.
    constant = .false.

    ! set mf_filekey to first part of file name:
    select case ( paramkey )
      !
      case ( 'LNSP'   ) ; read_2d=.true.; mdir='spectral'     ; fkey='pres.nlog.sfc'; sds_name='pres'
      case ( 'VO'     ) ; read_3d=.true.; mdir='spectral'     ; fkey='vort'         ; sds_name='vort'
      case ( 'D'      ) ; read_3d=.true.; mdir='spectral'     ; fkey='div'          ; sds_name='div'
      case ( 'Tv'     ) ; read_3d=.true.; mdir='spectral'     ; fkey='vair'         ; sds_name='vair'
      case ( 'Q'      ) ; read_3d=.true.; mdir='spectral'     ; fkey='shum'         ; sds_name='shum'
      case ( 'W'      ) ; read_3d=.true.; mdir='pressure'     ; fkey='omega'        ; sds_name='omega'
      !
      case ( 'oro'    ) ; read_2d=.true.; mdir='surface_gauss'; fkey='hgt.sfc'      ; sds_name='hgt'   ; constant=.true.
      case ( 'lsm'    ) ; read_2d=.true.; mdir='surface_gauss'; fkey='land.sfc'     ; sds_name='land'  ; constant=.true.
      case ( 'sr'     ) ; read_2d=.true.; mdir='surface_gauss'; fkey='sfcr.sfc'     ; sds_name='sfcr'
      case ( 'sps'    ) ; read_2d=.true.; mdir='surface_gauss'; fkey='pres.sfc'     ; sds_name='pres'
      case ( 'ci'     ) ; read_2d=.true.; mdir='surface_gauss'; fkey='icec.sfc'     ; sds_name='icec'
      case ( 'skt'    ) ; read_2d=.true.; mdir='surface_gauss'; fkey='skt.sfc'      ; sds_name='skt'
      case ( 'u10m'   ) ; read_2d=.true.; mdir='surface_gauss'; fkey='uwnd.10m'     ; sds_name='uwnd'
      case ( 'v10m'   ) ; read_2d=.true.; mdir='surface_gauss'; fkey='vwnd.10m'     ; sds_name='vwnd'
      case ( 'slhf'   ) ; read_2d=.true.; mdir='surface_gauss'; fkey='lhtfl.sfc'    ; sds_name='lhtfl'
      case ( 'sshf'   ) ; read_2d=.true.; mdir='surface_gauss'; fkey='shtfl.sfc'    ; sds_name='shtfl'
      case ( 'ewss'   ) ; read_2d=.true.; mdir='surface_gauss'; fkey='ugwd.sfc'     ; sds_name='ugwd'
      case ( 'nsss'   ) ; read_2d=.true.; mdir='surface_gauss'; fkey='vgwd.sfc'     ; sds_name='vgwd'
      case ( 'lsp'    ) ; read_2d=.true.; mdir='surface_gauss'; fkey='prate.sfc'    ; sds_name='prate'
      case ( 'cp'     ) ; read_2d=.true.; mdir='surface_gauss'; fkey='cprat.sfc'    ; sds_name='cprat'
      case ( 't2m'    ) ; read_2d=.true.; mdir='surface_gauss'; fkey='air.2m'       ; sds_name='air'
      case ( 'ssr'    ) ; read_2d=.true.; mdir='surface_gauss'; fkey='dswrf.sfc'    ; sds_name='dswrf'
      case ( 'sd'     ) ; read_2d=.true.; mdir='surface_gauss'; fkey='weasd.sfc'    ; sds_name='weasd'
      case ( 'swvl1'  ) ; read_2d=.true.; mdir='surface_gauss'; fkey='soilw.0-10cm' ; sds_name='soilw'
      !
      case default
        write (gol,'("unsupported paramkey `",a,"` for ncep 2d files")') paramkey; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! time index: number of 6 hour intervals from start of year.
    if ( constant ) then
      ! index 0 in data sets ...
      tstart = 0
    else
      ! base date:
      tbase = NewDate( year=mf%ccyy, month=1, day=1 )
      ! last 6 hours of previous year ?
      if ( t2 < tbase ) then
        ! first record covers (-6,0]
        tstart = 0
      else
        ! number of hours:
        hours = nint(rTotal( t2 - tbase, 'hour' ))
        ! index = 0, 1, 2, ...
        tstart = ceiling( real(hours)/real(6.0) )
      end if
    end if


    !
    ! ~~~ 2d field
    !

    if ( read_2d ) then

      ! example file names:
      !   pres.nlog.sfc.spec.2000.nc

      select case ( mdir )
        case ( 'surface_gauss' ) ; gridtype = 'gg' ; gridkey = '.gauss'
        case ( 'spectral'      ) ; gridtype = 'sh' ; gridkey = '.spec'
        case default
          write (gol,'("do not know gridtype for mdir : ",a)') trim(mdir); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! write filename
      if ( constant ) then
        write (fname,'(a,a,a,a,a,a,".",a)') &
                  trim(mf%dir), mf%pathsep, trim(mdir), mf%namesep, trim(fkey), trim(gridkey), 'nc'
      else
        write (fname,'(a,a,a,a,a,a,".",i4.4,".",a)') &
                  trim(mf%dir), mf%pathsep, trim(mdir), mf%namesep, trim(fkey), trim(gridkey), mf%ccyy, 'nc'
      end if

      ! file exist ?
      inquire( file=fname, exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found:")'); call goErr
        write (gol,'("  ",a)') trim(fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

      ! open file:
      call Init( hdf, trim(fname), 'read', status )
      IF_NOTOK_RETURN(status=1)

      ! check time
      if ( constant ) then
        tfactor = 1.0
      else
        call Check_Time( hdf, tstart, t1, t2, tfactor, status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! read 2d field: fill ggi/gg, or shi/sh :
      select case ( gridtype )
        case ( 'gg' )
          call pa_Init( gg1 )
          call Read_Gaussian_2d( hdf, sds_name, tstart, ggi, gg1, status )
          IF_NOTOK_RETURN(status=1)
          call pa_SetShape( gg, ggi%np, 1 )
          gg(:,1) = gg1
          call pa_Done( gg1 )
          ! apply time factor:
          gg = gg * tfactor
        case ( 'sh' )
          call pa_Init( sh1 )
          call Read_Spectral_2d( hdf, sds_name, tstart, shi, sh1, status )
          IF_NOTOK_RETURN(status=1)
          call pa_SetShape( sh, shi%np, 1 )
          sh(:,1) = sh1
          call pa_Done( sh1 )
          ! apply time factor:
          sh = sh * tfactor
        case default
          write (gol,'("unsupported grid type : ",a)') gridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! info ...
      write (key,'("sh==",i4.4)') shi%T
      call AddHistory( tmi, trim(key), status )

      ! close file:
      call Done( hdf, status )
      IF_NOTOK_RETURN(status=1)

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
          gg(:,1) = gg(:,1) * grav
        !
        ! for some probably historical reaseon, TM expects land/sea mask in % ...
        case ( 'lsm' )
          gg(:,1) = gg(:,1) * 100.0     ! 0-1  ->  %
        !
        ! adhoc surface roughness
        case ( 'sr' )
          write (gol,'("WARNING - adhoc constant value for surface roughness ...")'); call goPr
          gg(:,1) = 0.001   ! m
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

    end if   ! read 2d ?


    !
    ! ~~~ 3d field
    !

    if ( read_3d ) then

      ! example file names:
      !
      !           div.spec.2000.nc
      ! pres.nlog.sfc.spec.2000.nc
      !          vort.spec.2000.nc

      select case ( mdir )
        case ( 'spectral' ) ; gridtype = 'sh' ; gridkey = '.spec' ; levtype = 'sigma'
        case ( 'pressure' ) ; gridtype = 'll' ; gridkey = ''      ; levtype = 'pressure'
        case default
          write (gol,'("do not know grid- and levtype for mdir : ",a)') trim(mdir); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! write filename
      write (fname,'(a,a,a,a,a,a,".",i4.4,".",a)') &
                trim(mf%dir), mf%pathsep, trim(mdir), mf%namesep, trim(fkey), trim(gridkey), mf%ccyy, 'nc'

      ! file exist ?
      inquire( file=fname, exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found:")'); call goErr
        write (gol,'("  ",a)') trim(fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

      ! open file:
      call Init( hdf, trim(fname), 'read', status )
      IF_NOTOK_RETURN(status=1)

      ! check time
      if ( .not. constant ) then
        call Check_Time( hdf, tstart, t1, t2, tfactor, status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! extract level stuff:
      select case ( levtype )
        !
        case ( 'sigma' )
          ! open level data set:
          call Init( sds, hdf, 'level', status )
          IF_NOTOK_RETURN(status=1)
          ! extract dimensions:
          call GetInfo( sds, status, data_dims=data_dims )
          IF_NOTOK_RETURN(status=1)
          nlev = data_dims(1)
          !call ReadData( sds, sigma, status )
          !IF_NOTOK_RETURN(status=1)
          call Done( sds, status )
          IF_NOTOK_RETURN(status=1)
          ! level defintion
          select case ( nlev )
            case ( 28 )
              call Init( levi, 'nc28', status )
              IF_NOTOK_RETURN(status=1)
            case default
              write (gol,'("level definition not supported for nlev ",i4)') nlev; call goErr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select
        !
        case ( 'pressure' )
          ! open level data set:
          call Init( sds, hdf, 'level', status )
          IF_NOTOK_RETURN(status=1)
          ! extract dimensions:
          call GetInfo( sds, status, data_dims=data_dims )
          IF_NOTOK_RETURN(status=1)
          nlev = data_dims(1)
          ! read pressure levels:
          allocate( levels(nlev) )
          call ReadData( sds, levels, status )
          IF_NOTOK_RETURN(status=1)
          ! close data set:
          call Done( sds, status )
          IF_NOTOK_RETURN(status=1)
          ! level defintion
          write (gol,'("WARNING - adhoc implementation of pressure levels!")'); call goPr
          call Init( levi, nlev, (/levels,0.0/), (/levels*0.0,0.0/), status )
          IF_NOTOK_RETURN(status=1)
          ! clear
          deallocate( levels )
        !
        case default
          write (gol,'("level type not supported :",a)') levtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! info ...
      write (key,'("nlev==",i3.3)') nlev
      call AddHistory( tmi, trim(key), status )

      ! read 3d field: fill lli/ll, ggi/gg, or shi/sh
      select case ( gridtype )
        case ( 'll' )
          call Read_LonLat_3d( hdf, sds_name, tstart, nlev, nw, lli, ll, status )
          IF_NOTOK_RETURN(status=1)
          ! apply time factor:
          ll = ll * tfactor
        case ( 'sh' )
          call Read_Spectral_3d( hdf, sds_name, tstart, nlev, shi, sh, status )
          IF_NOTOK_RETURN(status=1)
          ! apply time factor:
          sh = sh * tfactor
        case default
          write (gol,'("unsupported grid type : ",a)') gridtype; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! info ...
      write (key,'("sh==",i4.4)') shT
      call AddHistory( tmi, trim(key), status )

      ! close file:
      call Done( hdf, status )
      IF_NOTOK_RETURN(status=1)

      ! unit conversion:
      select case ( paramkey )
        ! For some reason, the u/v/w from VO/D needs a factor -1 ...
        ! The minus is probably caused by the upwards coordinate system of ncep
        ! instead of the downward from ecmwf.
        case ( 'VO', 'D' )
          sh = - sh
      end select

      !
      ! ~~~ surface pressure
      !

      ! name of ncep file and data set:
      mdir = 'spectral'
      fkey = 'pres.nlog.sfc' ;  sds_name = 'pres'

      ! write surface pressure filename
      write (fname,'(a,"/",a,"-",a,a,".",i4.4,".",a)') &
                trim(mf%dir), trim(mdir), trim(fkey), '.spec', mf%ccyy, 'nc'

      ! file exist ?
      inquire( file=fname, exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found:")'); call goErr
        write (gol,'("  ",a)') trim(fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

      ! open file:
      call Init( hdf, trim(fname), 'read', status )
      IF_NOTOK_RETURN(status=1)

      ! check time
      call Check_Time( hdf, tstart, t1, t2, tfactor, status )
      IF_NOTOK_RETURN(status=1)

      ! read spectral field
      call Read_Spectral_2d( hdf, sds_name, tstart, shi, lnsp_sh, status )
      IF_NOTOK_RETURN(status=1)

      ! lnsp never accumulated, thus no time factor ...

      ! close file:
      call Done( hdf, status )
      IF_NOTOK_RETURN(status=1)

      ! unit conversion:
      !   sp * fac  =  exp( lnsp + nlog(fac) )
      !             =  exp( {sum_i=1,n c_i p_i} + nlog(fac) )
      !             =  exp( c_1 + {sum_i=2,n c_i p_i} + nlog(fac) )
      ! factor for conversion from cbar to Pa :
      !   [Pa] = [cbar] * 1e-2 [bar/cbar] * 1e5 [Pa/bar] = [mbar] * 1e3
      lnsp_sh(1) = lnsp_sh(1) + cmplx(log(1.0e3),0.0)

    end if   ! read 3d ?


    !
    ! ~~~ end
    !

    ! ok
    status = 0

  end subroutine mf_ReadRecord_1


  ! ***


  subroutine Check_Time( hdf, tstart, t1, t2, tfactor, status )

    use file_hdf  , only : THdfFile, TSds, Init, Done, ReadAttribute, ReadData
    use GO        , only : TDate, NewDate, IncrDate, wrtgol
    use GO        , only : operator(+), operator(-), rTotal
    use GO        , only :operator(/=), operator(==), operator(>=), operator(<=)

    ! --- in/out ---------------------------------------------

    type(THdfFile), intent(inout)    ::  hdf
    integer, intent(in)              ::  tstart
    type(TDate), intent(in)          ::  t1, t2
    real, intent(out)                ::  tfactor
    integer, intent(out)             ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Check_Time'

    ! --- local ----------------------------------------------

    type(TSds)            ::  sds

    integer               ::  base_date(3)
    real(8)               ::  time1(1), time(1), dhour
    integer               ::  idh
    type(TDate)           ::  tbase
    type(TDate)           ::  trec
    character(len=32)     ::  avg_period
    integer               ::  avg_dhour

    ! --- begin ----------------------------------------------

    ! o read start date:
    call ReadAttribute( hdf, 'base_date', base_date, status )
    IF_NOTOK_RETURN(status=1)
    tbase = NewDate( year=base_date(1), month=base_date(2), day=base_date(3) )

    ! init time dataset:
    call Init( sds, hdf, 'time', status )
    IF_NOTOK_RETURN(status=1)


    ! read hour offsets:
    call ReadData( sds, time1, status, start=(/0/) )
    if (status/=0) then
      write (gol,'("reading time value for start=0")'); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if
    call ReadData( sds, time, status, start=(/tstart/) )
    if (status/=0) then
      write (gol,'("reading time value for start=",i6)') tstart; call goErr
      call wrtgol( '  base_date   : ', tbase ); call goErr
      call wrtgol( '  t1          : ', t1 ); call goErr
      call wrtgol( '  t2          : ', t2 ); call goErr
      write (gol,'("  file name   : ",a)') trim(hdf%fname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! record time :
    dhour = time(1) - time1(1)
    trec = tbase + IncrDate( hour=floor(dhour), min=int((dhour-floor(dhour))*60.0) )

    ! by default, no time factor needs to applied:
    tfactor = 1.0

    ! instant field or time average ?
    if ( t1 == t2 ) then

      ! instant time; check:  trec  ==  t1==t2  ?
      if ( trec /= t1 ) then
        write (gol,'("time of supposed record does not match with requested time:")'); call goErr
        write (gol,'("  index        : ",i6)') tstart; call goErr
        write (gol,'("  base_time    : ",i4,2i3.2)') base_date; call goErr
        call wrtgol( '  record time  : ', trec ); call goErr
        call wrtgol( '  t1           : ', t1 ); call goErr
        call wrtgol( '  t2           : ', t2 ); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

    else if ( t1 <= t2 ) then

      ! interval; check :   [t1,t2]  in  [trec-6,trec]  ?

      ! identify averaging period:
      call ReadAttribute( sds, 'avg_period', avg_period, status )
      if (status/=0) then
        write (gol,'("reading avg_period; time inteval requested while not time average fields ?")'); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      select case ( avg_period )
        case ( '0000-00-00 06:00:00' ) ; avg_dhour = 6
        case default
          write (gol,'("unsupported avg_period : ",a)') trim(avg_period); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select

      ! check:
      if ( (t1 >= trec-IncrDate(hour=avg_dhour)) .and. (t2 <= trec) ) then

         ! ok, compute fraction:
         idh = nint(rTotal( t2 - t1, 'hour' ))
         tfactor = real(idh)/real(avg_dhour)

      else

        write (gol,'("time of supposed record does not match with requested time:")'); call goErr
        write (gol,'("  index        : ",i6)') tstart; call goErr
        write (gol,'("  base_time    : ",i4,2i3.2)') base_date; call goErr
        call wrtgol( '  record time  : ', trec ); call goErr
        write (gol,'("  avg_period   : ",a)') trim(avg_period); call goErr
        call wrtgol( '  t1           : ', t1 ); call goErr
        call wrtgol( '  t2           : ', t2 ); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return

      end if

    else

      write (gol,'("times should be the same or interval [t1,t2] : ")'); call goErr
      call wrtgol( '  t1  : ', t1 ); call goErr
      call wrtgol( '  t2  : ', t2 ); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return

    end if  ! instant or interval

    ! ok
    status = 0

  end subroutine Check_Time


  ! ***


  subroutine Read_Spectral_2d( hdf, sds_name, tstart, shi, sh, status )

    use parray  , only : pa_SetShape
    use file_hdf, only : THdfFile, TSds, Init, Done, ReadAttribute, ReadData
    use Grid    , only : TshGridInfo, Init

    ! --- in/out -------------------------------------

    type(THdfFile), intent(inout)     ::  hdf
    character(len=*), intent(in)      ::  sds_name
    integer, intent(in)               ::  tstart
    type(TShGridInfo),  intent(out)   ::  shi
    complex, pointer                  ::  sh(:)
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Read_Spectral_2d'

    ! --- local -------------------------------

    type(TSds)            ::  sds
    character(len=1)      ::  compression_used
    real                  ::  meanr(2,1)
    real                  ::  add_offset(1)
    real                  ::  scale_factor(1)
    integer               ::  trunc_count

    ! some cdc nc files contain a float for the packed data ...
    real, allocatable     ::  idata(:,:)

    ! --- begin ---------------------------------

    ! mean = first complex coeff
    call Init( sds, hdf, 'mean', status )
    IF_NOTOK_RETURN(status=1)
    call ReadData( sds, meanr, status, start=(/0,tstart/) )
    IF_NOTOK_RETURN(status=1)
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! packed ?
    call ReadAttribute( hdf, 'compression_used', compression_used, status )
    IF_NOTOK_RETURN(status=1)

    ! only packed yet ...
    if ( scan(compression_used,'Tt') == 0 ) then
      write (gol,'("only packed ncep data supported yet ...")'); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! read packing offset:
    call Init( sds, hdf, 'add_offset', status )
    IF_NOTOK_RETURN(status=1)
    call ReadData( sds, add_offset, status, start=(/tstart/) )
    IF_NOTOK_RETURN(status=1)
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! read packing factor:
    call Init( sds, hdf, 'scale_factor', status )
    IF_NOTOK_RETURN(status=1)
    call ReadData( sds, scale_factor, status, start=(/tstart/) )
    IF_NOTOK_RETURN(status=1)
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! open actual data set:
    call Init( sds, hdf, sds_name, status )
    IF_NOTOK_RETURN(status=1)

    ! read spectral truncation:
    call ReadAttribute( sds, 'trunc_count', trunc_count, status )
    IF_NOTOK_RETURN(status=1)

    ! setup output spectral definition:
    call Init( shi, trunc_count, status )
    IF_NOTOK_RETURN(status=1)

    ! data array:
    allocate( idata(shi%np*2,1) )

    ! read data:
    call ReadData( sds, idata, status, start=(/0,tstart/) )
    IF_NOTOK_RETURN(status=1)

    ! close data set:
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! setup output grid:
    call pa_SetShape( sh, shi%np )

    ! unpack and transform from realreal to complex:
    sh = transfer( ( idata(:,1) * scale_factor(1) ) + add_offset(1) , sh )

    ! add mean to first coeff:
    select case ( sds_name )
      case ( 'orog' )
        ! for some reason, adding the mean messes gives a bias in orography ..
      case default
        sh(1) = sh(1) + cmplx(meanr(1,1),meanr(2,1))
    end select

    ! convert from NCEP spectral coeff to ECMWF spectral coeff:
    sh = sh / sqrt(2.0)

    ! clear
    deallocate(  idata )

    ! ok
    status = 0

  end subroutine Read_Spectral_2d


  ! ***


  subroutine Read_Spectral_3d( hdf, sds_name, tstart, nlev, shi, sh, status )

    use parray  , only : pa_SetShape
    use file_hdf, only : THdfFile, TSds, Init, Done, ReadAttribute, ReadData
    use grid    , only : TshGridInfo, Init

    ! --- in/out -------------------------------------

    type(THdfFile), intent(inout)     ::  hdf
    character(len=*), intent(in)      ::  sds_name
    integer, intent(in)               ::  tstart
    integer, intent(in)               ::  nlev
    type(TShGridInfo),  intent(out)   ::  shi
    complex, pointer                  ::  sh(:,:)
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Read_Spectral_3d'

    ! --- local -------------------------------

    integer               ::  ilev
    type(TSds)            ::  sds
    character(len=1)      ::  compression_used
    real                  ::  meanr(nlev*2,1)
    real                  ::  add_offset(nlev,1)
    real                  ::  scale_factor(nlev,1)
    integer, allocatable  ::  idata(:,:,:)
    integer               ::  trunc_count

    ! --- begin ---------------------------------

    ! mean = first complex coeff
    call Init( sds, hdf, 'mean', status )
    IF_NOTOK_RETURN(status=1)
    call ReadData( sds, meanr, status, start=(/0,tstart/) )
    IF_NOTOK_RETURN(status=1)
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! packed ?
    call ReadAttribute( hdf, 'compression_used', compression_used, status )
    IF_NOTOK_RETURN(status=1)

    ! only packed yet ...
    if ( scan(compression_used,'Tt') == 0 ) then
      write (gol,'("only packed ncep data supported yet ...")'); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! read packing offset:
    call Init( sds, hdf, 'add_offset', status )
    IF_NOTOK_RETURN(status=1)
    call ReadData( sds, add_offset, status, start=(/0,tstart/) )
    IF_NOTOK_RETURN(status=1)
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! read packing factor:
    call Init( sds, hdf, 'scale_factor', status )
    IF_NOTOK_RETURN(status=1)
    call ReadData( sds, scale_factor, status, start=(/0,tstart/) )
    IF_NOTOK_RETURN(status=1)
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! open actual data set:
    call Init( sds, hdf, sds_name, status )
    IF_NOTOK_RETURN(status=1)

    ! read spectral truncation:
    call ReadAttribute( sds, 'trunc_count', trunc_count, status )
    IF_NOTOK_RETURN(status=1)

    ! setup output spectral definition:
    call Init( shi, trunc_count, status )
    IF_NOTOK_RETURN(status=1)

    ! data array:
    allocate( idata(shi%np*2,nlev,1) )

    ! read data:
    call ReadData( sds, idata, status, start=(/0,0,tstart/) )
    IF_NOTOK_RETURN(status=1)

    ! close data set:
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! setup output grid:
    call pa_SetShape( sh, shi%np, nlev )

    ! unpack and transform from realreal to complex:
    do ilev = 1, nlev
      sh(:,ilev) = transfer( ( idata(:,ilev,1) * scale_factor(ilev,1) ) + add_offset(ilev,1) , sh(:,ilev) )
    end do

    ! add mean to first coeff:
    do ilev = 1, nlev
      sh(1,ilev) = sh(1,ilev) + cmplx(meanr((ilev-1)*2+1,1),meanr((ilev-1)*2+2,1))
    end do

    ! convert from NCEP spectral coeff to ECMWF spectral coeff:
    sh = sh / sqrt(2.0)

    ! clear
    deallocate( idata        )

    ! ok
    status = 0

  end subroutine Read_Spectral_3d

  ! ***


  subroutine Read_Gaussian_2d( hdf, sds_name, tstart, ggi, gg, status )

    use parray  , only : pa_SetShape
    use file_hdf, only : THdfFile, TSds, Init, Done, ReadAttribute, ReadData, GetInfo
    use Grid    , only : TggGridInfo, Init

    ! --- in/out -------------------------------------

    type(THdfFile), intent(inout)     ::  hdf
    character(len=*), intent(in)      ::  sds_name
    integer, intent(in)               ::  tstart
    type(TggGridInfo),  intent(out)   ::  ggi
    real, pointer                     ::  gg(:)
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Read_Gaussian_2d'

    ! --- local -------------------------------

    type(TSds)            ::  sds
    integer               ::  data_dims(3)
    real                  ::  add_offset
    real                  ::  scale_factor
    integer, allocatable  ::  idata(:,:,:)

    ! --- begin ---------------------------------

    ! open actual data set:
    call Init( sds, hdf, sds_name, status )
    IF_NOTOK_RETURN(status=1)

    ! extract grid size:
    call GetInfo( sds, status, data_dims=data_dims )
    IF_NOTOK_RETURN(status=1)

    ! setup grid definition:
    call Init( ggi, data_dims(2)/2, .false., status )
    IF_NOTOK_RETURN(status=1)

    ! data array:
    allocate( idata(data_dims(1),data_dims(2),1) )

    ! read packing:
    call ReadAttribute( sds, 'add_offset', add_offset, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'scale_factor', scale_factor, status )
    IF_NOTOK_RETURN(status=1)

    ! read data:
    call ReadData( sds, idata, status, start=(/0,0,tstart/) )
    IF_NOTOK_RETURN(status=1)

    ! close data set:
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! setup output grid:
    call pa_SetShape( gg, ggi%np )

    ! unpack and transform from rank2 to rank1 :
    gg = transfer( ( idata(:,:,1) * scale_factor ) + add_offset , gg )

    ! clear
    deallocate(  idata )

    ! ok
    status = 0

  end subroutine Read_Gaussian_2d


  ! ***


  subroutine Read_LonLat_3d( hdf, sds_name, tstart, nlev, nw, lli, ll, status )

    use parray  , only : pa_SetShape
    use file_hdf, only : THdfFile, TSds, Init, Done, ReadAttribute, ReadData, GetInfo
    use grid    , only : TllGridInfo, Init

    ! --- in/out -------------------------------------

    type(THdfFile), intent(inout)     ::  hdf
    character(len=*), intent(in)      ::  sds_name
    integer, intent(in)               ::  tstart
    integer, intent(in)               ::  nlev
    character(len=*), intent(in)      ::  nw
    type(TllGridInfo),  intent(out)   ::  lli
    real, pointer                     ::  ll(:,:,:)
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Read_LonLat_3d'

    ! --- local -------------------------------

    integer               ::  j
    type(TSds)            ::  sds
    integer               ::  data_dims(4)
    real                  ::  add_offset
    real                  ::  scale_factor
    integer, allocatable  ::  idata(:,:,:,:)

    ! --- begin ---------------------------------

    ! open data set:
    call Init( sds, hdf, sds_name, status )
    IF_NOTOK_RETURN(status=1)

    ! extract grid size:
    call GetInfo( sds, status, data_dims=data_dims )
    IF_NOTOK_RETURN(status=1)

    ! setup grid definition:
    call Init( lli,   0.0, 360.0/ data_dims(1)   , data_dims(1), &
                    -90.0, 180.0/(data_dims(2)-1), data_dims(2), status )
    IF_NOTOK_RETURN(status=1)

    ! data array:
    allocate( idata(data_dims(1),data_dims(2),nlev,1) )

    ! read packing:
    call ReadAttribute( sds, 'add_offset', add_offset, status )
    IF_NOTOK_RETURN(status=1)
    call ReadAttribute( sds, 'scale_factor', scale_factor, status )
    IF_NOTOK_RETURN(status=1)

    ! read data:
    call ReadData( sds, idata, status, start=(/0,0,0,tstart/) )
    IF_NOTOK_RETURN(status=1)

    ! close data set:
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! setup output grid:
    select case ( nw )
      case ( 'n' )
        call pa_SetShape( ll, lli%nlon, lli%nlat, nlev )
      case ( 'w' )
        call pa_SetShape( ll, lli%nlon, lli%nlat, nlev+1 )
      case default
        write (gol,'("unsupported nw : ",a)') nw; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! initial zero:
    ll = 0.0

    ! unpack, and transform from north->south to south->north:
    do j = 1, lli%nlat
      ll(:,j,1:nlev) = real( ( idata(:,lli%nlat+1-j,:,1) * scale_factor ) + add_offset )
    end do

    ! clear
    deallocate(  idata )

    ! ok
    status = 0

  end subroutine Read_LonLat_3d



end module tmm_mf_ncep_cdc


!program test
!
!  use file_hdf
!  use grid
!
!  ! --- const ---------------------------
!
!  integer, parameter ::  T = 62
!
!  character(len=*), parameter  ::  fdir = '/misc/p71/co2/ncep.reanalysis/spectral/'
!
!  ! --- local ----------------------------------
!
!  type(THdfFile)  ::  nc, hdf
!  type(TSds)      ::  sds
!  integer         ::  status
!
!  real, allocatable     ::  rr(:)
!  real                  ::  mean(2,1)
!  real                  ::  add_offset(1)
!  real                  ::  scale_factor(1)
!  real, allocatable     ::  idata(:,:)
!
!  type(TshGridInfo)     ::  shi
!  complex, allocatable  ::  cc(:)
!  type(TshGrid)         ::  shc
!
!  integer  ::  tstart
!
!  type(TllGridInfo)  ::  lli
!  real   ::  ll(120,90)
!
!  ! --- begin ------------------------------------
!
!  print *, ''
!  print *, 'test: start'
!
!  call Init( shi, T, status )
!  if (status/=0) stop 'error'
!
!  allocate( cc(shi%np) )
!
!  allocate( idata(shi%np*2,1) )
!  allocate( rr(shi%np*2) )
!
!  !
!  ! time
!  !
!
!  tstart = 0
!
!  !
!  ! read
!  !
!
!  call Init( nc, fdir//'pres.nlog.sfc.spec.2000.nc', 'read', status )
!  if (status/=0) stop 'error'
!
!  call Init( sds, nc, 'mean', status )
!  if (status/=0) stop 'error'
!  call ReadData( sds, mean, status, start=(/0,tstart/) )
!  if (status/=0) stop 'error'
!  call Done( sds, status )
!  if (status/=0) stop 'error'
!
!  call Init( sds, nc, 'add_offset', status )
!  if (status/=0) stop 'error'
!  call ReadData( sds, add_offset, status, start=(/tstart/) )
!  if (status/=0) stop 'error'
!  call Done( sds, status )
!  if (status/=0) stop 'error'
!
!  call Init( sds, nc, 'scale_factor', status )
!  if (status/=0) stop 'error'
!  call ReadData( sds, scale_factor, status, start=(/tstart/) )
!  if (status/=0) stop 'error'
!  call Done( sds, status )
!  if (status/=0) stop 'error'
!
!  call Init( sds, nc, 'pres', status )
!  if (status/=0) stop 'error'
!  call ReadData( sds, idata, status, start=(/0,tstart/) )
!  if (status/=0) stop 'error'
!  call Done( sds, status )
!  if (status/=0) stop 'error'
!
!  call Done( nc, status )
!  if (status/=0) stop 'error'
!
!  !
!  ! sh grid
!  !
!
!  ! From http://www.cdc.noaa.gov/PublicData/faq.html#12 :
!  !
!  ! Most of the data in our netCDF files are packed. That is to say they have been
!  ! transformed by a scale factor and an add offset to reduce the storage needed to
!  ! two bytes per value. When you extract the short integers, you must unpack the
!  ! data to recover the correct floating point data values. Data files that contain
!  ! packed data will have a non-zero add offset and/or a scale factor not equal to 1.
!  !
!  ! The transformation is:
!  !  float_value = (short_value * scale_factor) + add_offset
!  !
!
!  print *, mean(1,1)
!  print *, add_offset
!  print *, scale_factor
!  print *, sum(idata)/size(idata)
!
!  rr = ( idata(:,1) * scale_factor(1) ) + add_offset(1)
!
!  print *, minval(rr), maxval(rr)
!  print *, rr(1), exp(rr(1))
!
!  ! convert to complex coeff:
!  cc = transfer( rr, cc )
!
!  ! convert from NOAA spectral coeff to ECMWF spectral coeff:
!  cc = cc / sqrt(2.0)
!
!  !
!  ! convert units from  nlog(bar)  to  nlog(Pa)
!  !  1 bar = 1e5 Pa
!  !  sp_Pa = exp( sp_nlog_cbar ) * 1e5
!  !        = exp( sp_nlog_cbar + nlog(1e5) )
!  !        = exp( sum cnm pnm  + nlog(1e5) )
!  !
!  ! add conversion offset to real part first complex coeff,
!  ! which represent the first global constant mode of lnsp :
!  cc(1) = cc(1) + cmplx(log(1.0e5),0.0)
!
!  call Init( shc )
!  call Set( shc, T, cc )
!
!  !
!  ! interpol
!  !
!
!  call Init( lli, -180.0+3.0/2, 3.0, 120, -90.0+2.0/2, 2.0, 90, status )
!  if (status/=0) stop 'error'
!
!  call Interpol( shc, lli, ll )
!
!  call Done( lli, status )
!  if (status/=0) stop 'error'
!
!  !
!  ! dump
!  !
!
!  call Init( hdf, 'lnsp.hdf', 'create', status )
!  if (status/=0) stop 'error'
!  call Init( sds, hdf, 'lnsp', shape(ll), 'real(4)', status )
!  if (status/=0) stop 'error'
!  call WriteData( sds, ll, status )
!  if (status/=0) stop 'error'
!  call Done( sds, status )
!  if (status/=0) stop 'error'
!  call Done( hdf, status )
!  if (status/=0) stop 'error'
!
!  !
!  ! done
!  !
!
!  deallocate( rr )
!  deallocate( idata )
!
!  deallocate( cc )
!  call Done( shc )
!  call Done( shi )
!
!  print *, 'test: end'
!  print *, ''
!
!end program test
