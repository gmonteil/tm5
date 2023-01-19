!
! NCEP Grib files.
!
! Grib tables: http://www.nco.ncep.noaa.gov/pmb/docs/on388/
!


module file_ncg

  use os_specs, only : MAX_FILENAME_LEN

  implicit none

  ! --- in/out ----------------------------------

  private

  public  ::  TNcepGrib
  public  ::  Init, Done
  public  ::  ReadRecord
  public  ::  CheckRecord
  public  ::  ListRecord
  public  ::  Get


  ! --- const -----------------------------------------------

  character(len=*), parameter ::  mname = 'file_ncg'


  ! --- types ----------------------------------

  type TNcepGrib
    ! file name:
    character(len=MAX_FILENAME_LEN) ::  fname
    ! unit for grib data file:
    integer                    ::  lugb
    ! unit for index file (?):
    integer                    ::  lugi
    ! message number
    integer                    ::  mnr
    ! product definition section:
    integer                    ::  pds(200)
    ! grid description section:
    integer                    ::  gds(200)
    ! bit map and binary data sections:
    integer                    ::  n
    logical, pointer           ::  bms(:)
    real, pointer              ::  bds(:)
  end type TNcepGrib


  ! --- interfaces ----------------------------------------

  interface Init
    module procedure ncg_Init
  end interface

  interface Done
    module procedure ncg_Done
  end interface

  interface ReadRecord
    module procedure ncg_ReadRecord
  end interface

  interface CheckRecord
    module procedure ncg_CheckRecord
  end interface

  interface ListRecord
    module procedure ncg_ListRecord
  end interface

  interface Get
    module procedure ncg_Get
  end interface


  ! --- local ------------------------------------------

  ! grib tables:
  character(len=4)    ::  type_of_level(0:255)
  character(len=5)    ::  parameter_abbrev(0:255)
  character(len=6)    ::  forecast_time_unit(0:4)
  character(len=20)   ::  time_range_indicator(0:10)
  character(len=2)    ::  data_representation_type(0:50)



contains


  ! =======================================================================


  subroutine ncg_Init( ncg, unit, fname, status )

    ! --- in/out ----------------------------------------

    type(TNcepGrib), intent(out)      ::  ncg
integer, intent(in)     ::  unit
    character(len=*), intent(in)      ::  fname
    integer, intent(out)              ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncg_Init'

    ! --- local ------------------------------------------

    logical               ::  exist, opened

    ! --- begin ------------------------------------------

    ! store file name:
    ncg%fname = fname

    ! file exist ?
    inquire( file=trim(ncg%fname), exist=exist )
    if ( .not. exist ) then
      write (*,'("ERROR - file not found :")')
      write (*,'("ERROR -   ",a)') trim(ncg%fname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

ncg%lugb = unit
!    ! select free file unit:
!    ! note: baopenr expects number within {1,..,99}
!    ncg%lugb = 10
!    do
!      inquire( unit=ncg%lugb, opened=opened )
!      if ( .not. opened ) exit
!      ncg%lugb = ncg%lugb + 1
!    end do

    ! open file:
#ifdef with_w3
    call baOpenR( ncg%lugb, trim(ncg%fname), status )
    if ( status /= 0 ) then
      write (*,'("ERROR - from baOpenR :")')
      write (*,'("ERROR -   file : ",a )') trim(ncg%fname)
      write (*,'("ERROR -   iret : ",i6)') status
      select case ( status )
        case ( 6 ) ; write (*,'("ERROR -   (file unit not in 1..999  : ",i6,")")') ncg%lugb
      end select
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
#else
    write (*,'("ERROR - not compiled with w3 library")')
    write (*,'("ERROR in ",a)') rname; status=1; return
#endif

    ! no index file yet ...
    ncg%lugi = 0

    ! no data sections read yet, but allocate some space already:
    ncg%n = 0
    allocate( ncg%bms(1000) )
    allocate( ncg%bds(1000) )

    ! fill param abbrevs (GRIB table 2)
    parameter_abbrev(1  ) = 'PRES '
    parameter_abbrev(11 ) = 'TMP  '
    parameter_abbrev(15 ) = 'TMAX '
    parameter_abbrev(16 ) = 'TMIN '
    parameter_abbrev(33 ) = 'UGRD '
    parameter_abbrev(34 ) = 'VGRD '
    parameter_abbrev(51 ) = 'SPFH '
    parameter_abbrev(54 ) = 'PWAT '
    parameter_abbrev(59 ) = 'PRATE'
    parameter_abbrev(65 ) = 'WEASD'
    parameter_abbrev(71 ) = 'TCDC '
    parameter_abbrev(81 ) = 'LAND '
    parameter_abbrev(84 ) = 'ALBDO'
    parameter_abbrev(90 ) = 'WATR '
    parameter_abbrev(91 ) = 'ICEC '
    parameter_abbrev(121) = 'LHTFL'
    parameter_abbrev(122) = 'SHTFL'
    parameter_abbrev(124) = 'UFLX '
    parameter_abbrev(125) = 'VFLX '
    parameter_abbrev(144) = 'SOILW'
    parameter_abbrev(145) = 'PEVPR'
    parameter_abbrev(146) = 'CWORK'
    parameter_abbrev(147) = 'U-GWD'
    parameter_abbrev(148) = 'V-GWD'
    parameter_abbrev(155) = 'GFLUX'
    parameter_abbrev(204) = 'DSWRF'
    parameter_abbrev(205) = 'DLWRF'
    parameter_abbrev(212) = 'ULWRF'
    parameter_abbrev(211) = 'USWRF'
    parameter_abbrev(214) = 'CPRAT'
    parameter_abbrev(221) = 'HPBL '

    ! type of level (GRIB table 3)                            level
    type_of_level(  1) = 'SFC'
    type_of_level(  8) = 'NTAT'
    type_of_level(105) = 'HTGL'  ! height above ground level  (meters)
    type_of_level(112) = 'DBLY'  ! layer between two depths   10=0-10cm, 2760=10-200cm
    type_of_level(200) = 'EATM'  ! entire atmosphere
    type_of_level(211) = 'BCY'   ! boundary layer cloud layer
    type_of_level(212) = 'LCBL'  ! low  cloud bottom level
    type_of_level(213) = 'LCTL'  ! low  cloud top    level
    type_of_level(214) = 'LCY'   ! low  cloud layer
    type_of_level(222) = 'MCBL'  ! mid  cloud bottom level
    type_of_level(223) = 'MCTL'  ! mid  cloud top    level
    type_of_level(224) = 'MCY'   ! mid  cloud layer
    type_of_level(232) = 'HCBL'  ! high cloud bottom level
    type_of_level(233) = 'HCTL'  ! high cloud top    level
    type_of_level(234) = 'HCY'   ! high cloud layer
    type_of_level(242) = 'CCBL'  ! conv cloud bottom level
    type_of_level(243) = 'CCTL'  ! conv cloud top    level
    type_of_level(244) = 'CCY'   ! conv cloud layer

    ! forecast time unit (GRIB table 4)
    forecast_time_unit(  0) = 'minute'
    forecast_time_unit(  1) = 'hour'
    forecast_time_unit(  2) = 'day'
    forecast_time_unit(  3) = 'month'
    forecast_time_unit(  4) = 'year'

    ! time range indicator
    time_range_indicator(  0) = 'tref+P1'
    time_range_indicator(  3) = 'tref+[P1,P2] aver'
    time_range_indicator(  4) = 'tref+[P1,P2] accum'
    time_range_indicator( 10) = 'tref+P1'

    ! data representation type (GRIB table 6)
    data_representation_type( 0) = 'll'
    data_representation_type( 4) = 'gg'
    data_representation_type(50) = 'sh'

    ! ok
    status = 0

  end subroutine ncg_Init



  ! ***


  subroutine ncg_Done( ncg, status )

    ! --- in/out ----------------------------------------

    type(TNcepGrib), intent(inout)      ::  ncg
    integer, intent(out)                ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncg_Done'

    ! --- begin ------------------------------------------

    ! clear
    deallocate( ncg%bms )
    deallocate( ncg%bds )

    ! close file:
#ifdef with_w3
    call baClose( ncg%lugb, status )
    if ( status /= 0 ) then
      write (*,'("ERROR - from baClose :")')
      write (*,'("ERROR -   file : ",a )') trim(ncg%fname)
      write (*,'("ERROR -   iret : ",i6)') status
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
#else
    write (*,'("ERROR - not compiled with w3 library")')
    write (*,'("ERROR in ",a)') rname; status=1; return
#endif

    ! ok
    status = 0

  end subroutine ncg_Done


  ! ***


  subroutine ncg_ReadRecord( ncg, status, irec, param, reftime, timerange, levtype, level )

    ! --- in/out ----------------------------------------

    type(TNcepGrib), intent(inout)          ::  ncg
    integer, intent(out)                    ::  status

    integer, intent(in), optional           ::  irec
    character(len=*), intent(in), optional  ::  param
    integer, intent(in), optional           ::  reftime(5)
    integer, intent(in), optional           ::  timerange(4)
    character(len=*), intent(in), optional  ::  levtype
    integer, intent(in), optional           ::  level

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncg_ReadRecord'

    ! --- local -----------------------------------------

    integer        ::  pds(200)
    integer        ::  gds(200)
    integer        ::  nskip
    integer        ::  i
    integer        ::  n

    ! --- begin ------------------------------------------

    ! wild cards by default
    pds = -1
    gds = -1

    ! number of grib messages to skip:
    nskip = 0
    if ( present(irec) ) nskip = irec - 1

    ! search param ?
    if ( present(param) ) then
      do i = 0, 255
        if ( param == parameter_abbrev(i) ) then
          pds(5) = i
          exit
        end if
      end do
    end if

    ! search reftime ?
    if ( present(reftime) ) then
      if ( all(reftime >= 0) ) then
        pds(21) = int(ceiling(reftime(1)/100.0))
        pds( 8) = modulo(reftime(1),100)
        pds( 9) = reftime(2)
        pds(10) = reftime(3)
        pds(11) = reftime(4)
        pds(12) = reftime(5)
      end if
    end if

    ! search time range ?
    if ( present(timerange) ) pds(13:16) = timerange

    ! search level type ?
    if ( present(levtype) ) then
      do i = 0, 255
        if ( levtype == type_of_level(i) ) then
          pds(6) = i
          exit
        end if
      end do
    end if

    ! search level ?
    if ( present(level) ) pds(7) = level

    ! loop until data section is large enough to hold the data:
    do
      n = size(ncg%bds)
#ifdef with_w3
      ! find and unpack grib message;
      ! no search ...
      call GetGB( ncg%lugb, ncg%lugi, n, nskip, pds, gds, &
                   ncg%n, ncg%mnr, ncg%pds, ncg%gds, ncg%bms, ncg%bds, status )
      select case ( status )
        case (  0 )
          ! all ok, thus leave loop:
          exit
        case ( 96 )
          write (*,'("ERROR - from GetGB : reading index file")')
        case ( 97 )
          write (*,'("ERROR - from GetGB : reading grib file")')
        case ( 98 )
          ! number of data points greater than size ncg%bds;
          ! check for excesive size first ...
          if ( n >= 1000000 ) then
            write (*,'("ERROR - excesive size of bds still not large enough:")')
            write (*,'("ERROR -   size bds : ",i8)') n
          else
            ! reallocate:
            deallocate( ncg%bms ) ; allocate( ncg%bms(n*10) )
            deallocate( ncg%bds ) ; allocate( ncg%bds(n*10) )
            ! try again ...
            cycle
          end if
        case ( 99 )
          write (*,'("ERROR - from GetGB : request not found")')
        case default
          write (*,'("ERROR - from GetGB : W3FI63 grib unpacker return code : ",i6)') status
          select case ( status )
            case (  1 ) ; write (*,'("ERROR -   = `GRIB` NOT FOUND IN FIRST 100 CHARS                ")')
            case (  2 ) ; write (*,'("ERROR -   = `7777` NOT IN CORRECT LOCATION                     ")')
            case (  3 ) ; write (*,'("ERROR -   = UNPACKED FIELD IS LARGER THAN 260000               ")')
            case (  4 ) ; write (*,'("ERROR -   = GDS/ GRID NOT ONE OF CURRENTLY ACCEPTED VALUES     ")')
            case (  5 ) ; write (*,'("ERROR -   = GRID NOT CURRENTLY AVAIL FOR CENTER INDICATED      ")')
            case (  8 ) ; write (*,'("ERROR -   = TEMP GDS INDICATED, BUT GDS FLAG IS OFF            ")')
            case (  9 ) ; write (*,'("ERROR -   = GDS INDICATES SIZE MISMATCH WITH STD GRID          ")')
            case ( 10 ) ; write (*,'("ERROR -   = INCORRECT CENTER INDICATOR                         ")')
            case ( 11 ) ; write (*,'("ERROR -   = BINARY DATA SECTION (BDS) NOT COMPLETELY PROCESSED.")') &
                        ; write (*,'("ERROR -     PROGRAM IS NOT SET TO PROCESS FLAG COMBINATIONS    ")') &
                        ; write (*,'("ERROR -     SHOWN IN OCTETS 4 AND 14.                          ")')
            case ( 12 ) ; write (*,'("ERROR -   = BINARY DATA SECTION (BDS) NOT COMPLETELY PROCESSED.")') &
                        ; write (*,'("ERROR -     PROGRAM IS NOT SET TO PROCESS FLAG COMBINATIONS    ")')
          end select
      end select
#else
      write (*,'("ERROR - not compiled with w3 library")')
      write (*,'("ERROR in ",a)') rname; status=1; return
#endif
      write (*,'("ERROR -   file      : ",a )') trim(ncg%fname)
      write (*,'("ERROR -   irec      : ",i6     )') nskip+1
      if ( present(param    ) ) write (*,'("ERROR -   param     : ",a,i4   )') param, pds(5)
      if ( present(reftime  ) ) write (*,'("ERROR -   reftime   : ",i4,4i3 )') reftime
      if ( present(reftime  ) ) write (*,'("ERROR -     (pds    : ",6i3,")")') pds(21), pds(8:12)
      if ( present(timerange) ) write (*,'("ERROR -   timerange : ",4i4    )') timerange
      if ( present(levtype  ) ) write (*,'("ERROR -   levtype   : ",a,i4   )') levtype, pds(6)
      if ( present(level    ) ) write (*,'("ERROR -   level     : ",i6     )') level
      write (*,'("ERROR in ",a)') rname; status=1; return
    end do

    ! ok
    status = 0

  end subroutine ncg_ReadRecord


  ! ***


  subroutine ncg_CheckRecord( ncg, status, param, reftime, timerange, levtype, level )

    ! --- in/out ----------------------------------------

    type(TNcepGrib), intent(inout)          ::  ncg
    integer, intent(out)                    ::  status

    character(len=*), intent(in), optional  ::  param
    integer, intent(in), optional           ::  reftime(5)
    integer, intent(in), optional           ::  timerange(4)
    character(len=*), intent(in), optional  ::  levtype
    integer, intent(in), optional           ::  level

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncg_CheckRecord'

    ! --- local -----------------------------------------

    integer        ::  pds(200)

    ! --- begin ------------------------------------------

    ! initially ok
    status = 0

    ! check param ?
    if ( present(param) ) then
      if ( param /= parameter_abbrev(ncg%pds(5)) ) then
        write (*,'("ERROR - parameters do not match : ")')
        write (*,'("ERROR -   requested   : ",a)') param
        write (*,'("ERROR -   grib record : ",a)') parameter_abbrev(ncg%pds(5))
        status = 1
      end if
    end if

    ! check reftime ?
    if ( present(reftime) ) then
      if ( all(reftime >= 0) ) then
        pds(21) = int(ceiling(reftime(1)/100.0))
        pds( 8) = modulo(reftime(1),100)
        pds( 9) = reftime(2)
        pds(10) = reftime(3)
        pds(11) = reftime(4)
        pds(12) = reftime(5)
        if ( any( (/pds(21),pds(8:12)/) /= (/ncg%pds(21),ncg%pds(8:12)/) ) ) then
          write (*,'("ERROR - reference times do not match : ")')
          write (*,'("ERROR -   requested   : ",6i3)') pds(21), pds(8:12)
          write (*,'("ERROR -   grib record : ",6i3)') ncg%pds(21), ncg%pds(8:12)
          status = 1
        end if
      end if
    end if

    ! check time range ?
    if ( present(timerange) ) then
      if ( all(timerange >= 0) .and. any( timerange /= ncg%pds(13:16) ) ) then
        write (*,'("ERROR - time ranges do not match : ")')
        write (*,'("ERROR -   requested   : ",6i3)') timerange
        write (*,'("ERROR -   grib record : ",6i3)') ncg%pds(13:16)
        status = 1
      end if
    end if

    ! check level type ?
    if ( present(levtype) ) then
      if ( levtype /= type_of_level(ncg%pds(6)) ) then
        write (*,'("ERROR - level types do not match : ")')
        write (*,'("ERROR -   requested   : ",a)') levtype
        write (*,'("ERROR -   grib record : ",a)') type_of_level(ncg%pds(6))
        status = 1
      end if
    end if

    ! check level ?
    if ( present(level) ) then
      if ( level /= ncg%pds(7) ) then
        write (*,'("ERROR - levels do not match : ")')
        write (*,'("ERROR -   requested   : ",6i3)') level
        write (*,'("ERROR -   grib record : ",6i3)') ncg%pds(7)
        status = 1
      end if
    end if

    ! any errors ?
    if ( status /= 0 ) then
      write (*,'("ERROR - some checks failed ...")')
      write (*,'("ERROR -   checks requrested for:")')
      if ( present(param    ) ) write (*,'("ERROR -     param     : ",a,i4   )') param, pds(5)
      if ( present(reftime  ) ) write (*,'("ERROR -     reftime   : ",i4,4i3 )') reftime
      if ( present(reftime  ) ) write (*,'("ERROR -       (pds    : ",6i3,")")') pds(21), pds(8:12)
      if ( present(timerange) ) write (*,'("ERROR -     timerange : ",4i4    )') timerange
      if ( present(levtype  ) ) write (*,'("ERROR -     levtype   : ",a,i4   )') levtype, pds(6)
      if ( present(level    ) ) write (*,'("ERROR -     level     : ",i6     )') level
      write (*,'("ERROR -   file      : ",a )') trim(ncg%fname)
      write (*,'("ERROR -   record nr : ",i6)') ncg%mnr
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine ncg_CheckRecord


  ! ***


  subroutine ncg_ListRecord( ncg, status )

    ! --- in/out ----------------------------------------

    type(TNcepGrib), intent(inout)      ::  ncg
    integer, intent(out)                ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncg_ListRecord'

    ! --- local ---------------------------------------

    integer       ::  i

    ! --- begin ------------------------------------------

    write (*,'("                                                              ")')
    write (*,'("Section 1 - Product Definition Section                        ")')
    write (*,'("                                                              ")')
    write (*,'(" ( 1)  - ID OF CENTER                                       : ",i6)') ncg%pds(1)
    write (*,'(" (23)  - SUBCENTER NUMBER                                   : ",i6)') ncg%pds(23)
    write (*,'(" ( 2)  - GENERATING PROCESS ID NUMBER                       : ",i6)') ncg%pds(2)
    write (*,'(" (18)  - VERSION NR OF GRIB SPECIFICATION                   : ",i6)') ncg%pds(18)
    write (*,'(" ( 3)  - GRID DEFINITION                                    : ",i6)') ncg%pds(3)
    write (*,'(" ( 4)  - GDS/BMS FLAG (RIGHT ADJ COPY OF OCTET 8)           : ",i6)') ncg%pds(4)
    write (*,'(" (19)  - VERSION NR OF PARAMETER TABLE                      : ",i6)') ncg%pds(19)
    write (*,'(" (17)  - NUMBER INCLUDED IN AVERAGE                         : ",i6)') ncg%pds(17)
    write (*,'(" (20)  - NR MISSING FROM AVERAGE/ACCUMULATION               : ",i6)') ncg%pds(20)
    write (*,'(" (22)  - UNITS DECIMAL SCALE FACTOR                         : ",i6)') ncg%pds(22)
    write (*,'(" ( 5)  - INDICATOR OF PARAMETER                             : ",i6,"  =  ",a)') &
                            ncg%pds(5), parameter_abbrev(ncg%pds(5))
    write (*,'(" ( 6)  - TYPE OF LEVEL                                      : ",i6,"  =  ",a)') &
                            ncg%pds(6), type_of_level(ncg%pds(6))
    write (*,'(" ( 7)  - HEIGHT/PRESSURE , ETC OF LEVEL                     : ",i6)') ncg%pds(7)
    write (*,'(" (21)  - CENTURY OF REFERENCE TIME OF DATA                  : ",i6)') ncg%pds(21)
    write (*,'(" ( 8)  - YEAR INCLUDING (CENTURY-1)                         : ",i6)') ncg%pds(8)
    write (*,'(" ( 9)  - MONTH OF YEAR                                      : ",i6)') ncg%pds(9)
    write (*,'(" (10)  - DAY OF MONTH                                       : ",i6)') ncg%pds(10)
    write (*,'(" (11)  - HOUR OF DAY                                        : ",i6)') ncg%pds(11)
    write (*,'(" (12)  - MINUTE OF HOUR                                     : ",i6)') ncg%pds(12)
    write (*,'(" (13)  - INDICATOR OF FORECAST TIME UNIT                    : ",i6,"  =  ",a)') &
                            ncg%pds(13), forecast_time_unit(ncg%pds(13))
    write (*,'(" (14)  - TIME RANGE 1                                       : ",i6)') ncg%pds(14)
    write (*,'(" (15)  - TIME RANGE 2                                       : ",i6)') ncg%pds(15)
    write (*,'(" (16)  - TIME RANGE FLAG                                    : ",i6,"  =  ",a)') &
                            ncg%pds(16), time_range_indicator(ncg%pds(16))
    write (*,'("                                                              ")')
    write (*,'("Section 2 - Grid Description Section                          ")')
    write (*,'("                                                              ")')
    write (*,'(" ( 1)  - DATA REPRESENTATION TYPE                           : ",i6,"  =  ",a)') &
                            ncg%gds( 1), data_representation_type(ncg%gds( 1))
    select case ( data_representation_type(ncg%gds( 1)) )
      case ( 'll' )
        write (*,'("                                                              ")')
        write (*,'("LATITUDE/LONGITUDE GRID                                       ")')
        write (*,'(" ( 2)  - N(I) NR POINTS ON LATITUDE CIRCLE                  : ",i6)') ncg%gds( 2)
        write (*,'(" ( 3)  - N(J) NR POINTS ON LONGITUDE MERIDIAN               : ",i6)') ncg%gds( 3)
        write (*,'(" ( 4)  - LA(1) LATITUDE OF ORIGIN                           : ",i6)') ncg%gds( 4)
        write (*,'(" ( 5)  - LO(1) LONGITUDE OF ORIGIN                          : ",i6)') ncg%gds( 5)
        write (*,'(" ( 6)  - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)       : ",i6)') ncg%gds( 6)
        write (*,'(" ( 7)  - LA(2) LATITUDE OF EXTREME POINT                    : ",i6)') ncg%gds( 7)
        write (*,'(" ( 8)  - LO(2) LONGITUDE OF EXTREME POINT                   : ",i6)') ncg%gds( 8)
        write (*,'(" ( 9)  - DI LONGITUDINAL DIRECTION OF INCREMENT             : ",i6)') ncg%gds( 9)
        write (*,'(" (10)  - DJ LATITUDINAL DIRECTION INCREMENT                 : ",i6)') ncg%gds(10)
        write (*,'(" (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)    : ",i6)') ncg%gds(11)
      case ( 'gg' )
        write (*,'("                                                              ")')
        write (*,'("GAUSSIAN  GRID                                                ")')
        write (*,'(" ( 2)  - N(I) NR POINTS ON LATITUDE CIRCLE                  : ",i6)') ncg%gds( 2)
        write (*,'(" ( 3)  - N(J) NR POINTS ON LONGITUDE MERIDIAN               : ",i6)') ncg%gds( 3)
        write (*,'(" ( 4)  - LA(1) LATITUDE OF ORIGIN                           : ",i6)') ncg%gds( 4)
        write (*,'(" ( 5)  - LO(1) LONGITUDE OF ORIGIN                          : ",i6)') ncg%gds( 5)
        write (*,'(" ( 6)  - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)      : ",i6)') ncg%gds( 6)
        write (*,'(" ( 7)  - LA(2) LATITUDE OF EXTREME POINT                    : ",i6)') ncg%gds( 7)
        write (*,'(" ( 8)  - LO(2) LONGITUDE OF EXTREME POINT                   : ",i6)') ncg%gds( 8)
        write (*,'(" ( 9)  - DI LONGITUDINAL DIRECTION OF INCREMENT             : ",i6)') ncg%gds( 9)
        write (*,'(" (10)  - N - NR OF CIRCLES POLE TO EQUATOR                  : ",i6)') ncg%gds(10)
        write (*,'(" (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)    : ",i6)') ncg%gds(11)
      case ( 'sh' )
        write (*,'("                                                              ")')
        write (*,'("SPHERICAL HARMONIC COEFFICIENTS                               ")')
        write (*,'(" ( 2)  - J PENTAGONAL RESOLUTION PARAMETER                  : ",i6)') ncg%gds( 2)
        write (*,'(" ( 3)  - K      x          x         x                      : ",i6)') ncg%gds( 3)
        write (*,'(" ( 4)  - M      x          x         x                      : ",i6)') ncg%gds( 4)
        write (*,'(" ( 5)  - REPRESENTATION TYPE                                : ",i6)') ncg%gds( 5)
        write (*,'(" ( 6)  - COEFFICIENT STORAGE MODE                           : ",i6)') ncg%gds( 6)
    end select
    if ( (ncg%gds(12) > 0) .or.(ncg%gds(19) > 0) ) then
      write (*,'("                                                                       ")')
      write (*,'("VERTICAL COORDINATE PARAMETERS                                         ")')
      write (*,'(" (12)  - NV - NR OF VERT COORD PARAMETERS                            : ",i6)') ncg%gds(12)
      write (*,'(" (13)  - PV - OCTET NR OF LIST OF VERT COORD PARAMETERS              : ",i6)') ncg%gds(13)
      write (*,'("         PL - LOCATION OF THE LIST OF NUMBERS OF POINTS IN EACH ROW    ")')
      write (*,'("              (IF NO VERT COORD PARAMETERS ARE PRESENT)                ")')
      write (*,'("         255 IF NEITHER ARE PRESENT                                    ")')
      write (*,'(" (19)  - NUMBER OF VERTICAL COORDINATE PARAMETERS                    : ",i6)') ncg%gds(19)
      write (*,'(" (20)  - OCTET NUMBER OF THE LIST OF VERTICAL COORDINATE PARAMETERS  : ",i6)') ncg%gds(20)
      write (*,'("         OCTET NUMBER OF THE LIST OF NUMBERS OF POINTS IN EACH ROW   : ")')
      write (*,'("         255 IF NEITHER ARE PRESENT                                  : ")')
      write (*,'(" (21)  - FOR GRIDS WITH PL, NUMBER OF POINTS IN GRID                 : ",i6)') ncg%gds(21)
      write (*,'(" (22)  - NUMBER OF WORDS IN EACH ROW                                 : ",i6)') ncg%gds(22)
    end if
    write (*,'("                                                                       ")')
    write (*,'("Section 4 - Binary Data Section                                        ")')
    write (*,'("                                                                       ")')
    write (*,'("  first 20 values:                                                     ")')
    do i = 1, 20
      write (*,'("    ",es16.6)') ncg%bds(i)
    end do
    write (*,'("                                                                       ")')

    ! ok
    status = 0

  end subroutine ncg_ListRecord


  ! ***


  subroutine ncg_Get( ncg, status, ggN, gg )

    ! --- in/out ----------------------------------------

    type(TNcepGrib), intent(inout)      ::  ncg
    integer, intent(out)                ::  status

    integer, intent(out), optional      ::  ggN
    real, intent(out), optional         ::  gg(:)

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncg_Get'

    ! --- begin ------------------------------------------

    ! return gg number:
    if ( present(ggN) ) then
      if ( data_representation_type(ncg%gds(1)) /= 'gg' ) then
        write (*,'("ERROR - ggN not for grid type : ",a)') data_representation_type(ncg%gds(1))
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
      ggN = ncg%gds(10)
    end if

    ! return gg field:
    if ( present(gg) ) then
      if ( data_representation_type(ncg%gds(1)) /= 'gg' ) then
        write (*,'("ERROR - ggN not for grid type : ",a)') data_representation_type(ncg%gds(1))
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
      if ( size(gg) /= ncg%n ) then
        write (*,'("ERROR - size of output grid does not match with storred data:")')
        write (*,'("ERROR -   storred    : ",i6)') ncg%n
        write (*,'("ERROR -   output gg  : ",i6)') size(gg)
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
      gg = ncg%bds(1:ncg%n)
    end if

    ! ok
    status = 0

  end subroutine ncg_Get

end module file_ncg




! #####################################################################################
! ###
! ### test
! ###
! #####################################################################################
!
! ifort -o test.x -r8 file_ncg.f90 /ahome/segers/opt/ifort-default-R64/lib/libw3.a  &&  ./test.x
!
! Problem: reading from a new grib file gives an error from GetGB .
!
! Tested:
!  o set all input/output variables to zero before second file
!  o multiple open/read/close from same file is not a problem
!  o copy files to same temporary file name does not help
!
!
!program test
!
!  !use file_ncg
!
!  !type(TNcepGrib)   ::  ncg
!  integer           ::  status
!  integer           ::  irec
!  logical           ::  opened
!
!  integer           ::  pds(200)
!  integer           ::  gds(200)
!
!  integer                    ::  lugb
!  integer                    ::  lugi
!  integer                    ::  ncg_n
!  integer                    ::  ncg_pds(200)
!  integer                    ::  ncg_gds(200)
!  integer                    ::  mnr
!
!  integer, parameter         ::  n = 400000
!  logical                    ::  bms(n)
!  real                       ::  bds(n)
!
!  ! --- begin ------------------------------------
!
!  print *, 'test: begin'
!
!  do irec = 31, 31
!
!  print *, 'test: 1', irec
!
!    !call system( '/bin/cp -f /p71/co2/ncep.gfs/20040614/20040614.06.00.SFLUXGrbF /p71/co2/segers/tmp/tmp.gb' )
!
!    lugb = 10
!    call baOpenR( lugb, '/p71/co2/ncep.gfs/20040614/20040614.06.00.SFLUXGrbF', status )
!    if (status/=0) stop 'ERROR in test'
!
!    lugi = 0
!    pds = -1
!    gds = -1
!    call GetGB( lugb, lugi, n, irec-1, pds, gds, &
!                   ncg_n, mnr, ncg_pds, ncg_gds, bms, bds, status )
!    print *, '    getgb : ', status
!    if (status/=0) stop 'ERROR in test'
!
!    !call baClose( lugb, status )
!    !if (status/=0) stop 'ERROR in test'
!
!    inquire( lugb, opened=opened )
!    print *, '    opened : ', opened
!
!    close(lugb,iostat=status)
!    if (status/=0) stop 'ERROR in test from close'
!
!  end do
!
!  print *, 'test: 2'
!
!    !call system( '/bin/cp -f /p71/co2/ncep.gfs/20040614/20040614.12.00.SFLUXGrbF /p71/co2/segers/tmp/tmp.gb' )
!
!    lugb = 10
!    call baOpenR( lugb, '/p71/co2/ncep.gfs/20040614/20040614.12.00.SFLUXGrbF', status )
!    if (status/=0) stop 'ERROR in test'
!
!    lugi = 0
!    irec = 31
!    pds = -1
!    gds = -1
!    ncg_n = 0
!    mnr   = 0
!    ncg_pds = 0
!    ncg_gds = 0
!    bms     = .false.
!    bds     = 0.0
!    status  = 0
!    call GetGB( lugb, lugi, n, irec-1, pds, gds, &
!                   ncg_n, mnr, ncg_pds, ncg_gds, bms, bds, status )
!    print *, '    getgb : ', status
!    if (status/=0) stop 'ERROR in test'
!
!    call baClose( lugb, status )
!    if (status/=0) stop 'ERROR in test'
!
!  print *, 'test: end'
!
!end program test
