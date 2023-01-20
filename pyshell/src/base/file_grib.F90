!
!ProTeX: 1.14-AJS
!
!BOI
!
! !TITLE:        File_GRIB - Interface to GRIB files
! !AUTHORS:      Arjo Segers
! !AFFILIATION:  KNMI
! !DATE:         \today
!
! !INTRODUCTION: Usage
!
!   \bv
!   Provides a module:
!
!     use GribFile, only : ...
!
!   The module provides a structure type containing all what should
!   be known to open, read, and close a grib file:
!
!     type(TGribFile)   ::  gribfile
!
!   The gribfile is opened either for reading or writing with
!   an 'Init' routine, which is in fact an interface to 'PbOpen'.
!   A type TgribFile contains buffers to store one record
!   of a grib file (a 'grib message').
!
!   1. To open a grib file for reading, use:
!
!       call Init( gribfile, filename, 'r', status )
!          type(TGribFile), intent(inout)   ::  gribfile
!          character(len=*), intent(in)     ::  filename
!          character(len=*), intent(in)     ::  mode
!          integer, intent(out)             ::  status
!
!      To read next record in buffer, use:
!
!        call ReadRecord( gribfile, status )
!
!      Contents of the message is extracted with the 'Get' routine,
!      for example:
!
!        call Get( gribfile, status, pid=thepid )
!
!      Complete syntaxis and description of arguments
!      (see also the 'GRIB SECTIONS' part later on):
!
!        call Get( gribfile, status, &
!                  model_id, pid, &
!                  levtype, level, hyb_a, hyb_b, &
!                  reftime, timerange, &
!                  gridtype, &
!                  ll, &
!                  lon_first, lon_last, lon_inc, lon_n, &
!                  lat_first, lat_last, lat_inc, lat_n, &
!                  T, sh, &
!                  N, gg  )
!
!          integer, intent(out), optional   ::  model_id
!          ! Element 3 of section 1.
!
!          integer, intent(out), optional   ::  pid
!          ! Parameter identification.
!          ! See the integer parameters 'pid_' in the code.
!          ! Element 6 of section 1.
!
!          integer, intent(out), optional   ::  levtype, level
!          ! Elements 7-9 of section 1.
!          ! --> Predefined values for level type:
!          !  integer, parameter  ::  levtype_sfc = 1    ! surface
!          !  integer, parameter  ::  levtype_hyb = 109  ! hybrid level
!
!          real, intent(out), optional      ::  hyb_a(:), hyb_b(:)
!          ! Vertical parameters in real section 2
!
!          integer, intent(out), optional   ::  reftime(5)
!          ! Reference time: a 5 element integer array with
!          !   year_inc_century, month, day, hour, minutes
!          ! Elements 10-14,21 of section 1.
!
!          integer, intent(out), optional   ::  timerange(4)
!          ! Time range: a 4 element integer array with
!          !   time range unit, value 1, value 2, indicator
!          ! Elements 15-18 of section 1.
!
!          integer, intent(out), optional   ::  gridtype
!          ! Element 1 of section 2.
!          ! --> Predefined values:
!          !  integer, parameter  ::  gridtype_ll = 0      ! lat/lon
!          !  integer, parameter  ::  gridtype_gg = 4      ! gaussian grid
!          !  integer, parameter  ::  gridtype_sh = 50     ! spectral
!
!          integer, intent(out), optional   ::  lon_first, lon_last, lon_inc, lon_n
!          integer, intent(out), optional   ::  lat_first, lon_last, lat_inc, lat_n
!          ! Elements of section 2
!          ! NOTE: values for lat/lon in mili degrees !
!
!          real, intent(out), optional      ::  ll(:,:)
!          ! Value of lat/lon grid (stored from west->east, north->south!)
!
!          complex, intent(out), optional   ::  sh(:)
!          integer, intent(out), optional   ::  T
!          ! spectral coefficients;
!          ! size should match with spectral truncation 'T'
!
!          real, intent(out), optional      ::  gg(:)
!          integer, intent(out), optional   ::  N
!          ! data values on Gaussian grid;
!          ! number of values should match grid size 'N'
!
!      Contents could be checked with the 'Check' routine,
!      for example:
!
!        call Check( gribfile, status, pid=154 )
!
!      Complete syntaxis:
!
!        call Check( gribfile, status, debug=1, &
!                    model_id, pid, &
!                    levtype, level, &
!                    reftime, timerange, &
!                    gridtype, &
!                    lon_first, lon_last, lon_inc, lon_n, &
!                    lat_first, lat_last, lat_inc, lat_n, &
!                    T  )
!
!        type(TGribFile), intent(in)      ::  gribfile
!        integer, intent(out)             ::  status
!
!        integer, intent(in), optional    ::  debug
!
!        integer, intent(in), optional    ::  model_id
!        integer, intent(in), optional    ::  pid
!        integer, intent(in), optional    ::  levtype, level
!        integer, intent(in), optional    ::  reftime(5)
!        integer, intent(in), optional    ::  timerange(4)
!        integer, intent(in), optional    ::  gridtype
!        integer, intent(in), optional    ::  lon_first, lon_last, lon_inc, lon_n
!        integer, intent(in), optional    ::  lat_first, lat_last, lat_inc, lat_n
!        integer, intent(in), optional    ::  T
!
!        integer, intent(in), optional    ::  status
!        integer, intent(out), optional   ::  status
!
!      See the 'Get' command above for a description.
!      In case of one or more fields not matching the grib field:
!       o return status is equal to number of tests failed;
!       o info about the failed test is printed if debug is present and >0 .
!
!      Note: hybride coefficients can not checked due to a
!        lack of precission ... Similar for grid values (ll, gg, and sh).
!
!      The grib file is closed with a 'Done' routine
!      (interface to 'PbClose'):
!
!        call Done( gribfile, status )
!
!
!   2. To open a grib file for writing, use:
!
!        call Init( gribfile, 'test.gb', 'w', status )
!
!      For safety, set all sections in the buffer to zero
!      or an other safe value:
!
!        call Clear( gribfile, status )
!
!      Fill some or all of the appropriate fields with
!      the 'Set' command:
!
!        call Set( gribfile, status,  &
!                  model_id, pid, &
!                  levtype, hyb_a, hyb_b, level, &
!                  reftime, timerange, &
!                  gridtype, &
!                  ll, lon_n, lat_n, &
!                  lon_first, lon_inc, lon_last, &
!                  lat_first, lat_inc, lat_last, &
!                  scanning_mode, nbits )
!
!          type(TGribFile), intent(inout)  ::  gribfile
!          integer, intent(out)            ::  status
!
!          integer, intent(in), optional   ::  debug
!
!          integer, intent(in), optional   ::  model_id
!          integer, intent(in), optional   ::  pid
!          integer, intent(in), optional   ::  levtype
!          integer, intent(in), optional   ::  level
!          real, intent(in), optional      ::  hyb_a(:), hyb_b(:)
!          integer, intent(in), optional   ::  reftime(5)   ! (/yy,mm,dd,hh,min/)
!          integer, intent(in), optional   ::  timerange(4)
!          integer, intent(in), optional   ::  gridtype
!          real, intent(in), optional      ::  ll(:,:)
!          integer, intent(in), optional   ::  lon_n, lat_n
!          integer, intent(in), optional   ::  lon_first, lon_inc, lon_last  ! mili degree
!          integer, intent(in), optional   ::  lat_first, lat_inc, lat_last  ! mili degree
!
!      See the 'Get' command for a description of most of the fields.
!      Only used as arguments to 'Set':
!
!          integer, intent(in), optional        ::  scanning_mode
!          ! Something with storage order; does not work proper yet.
!
!          integer, intent(in), optional        ::  nbits
!          ! Number of bits used to store a real value; default 24 .
!
!      To write the record to the file, use:
!
!        call WriteRecord( gribfile, status )
!
!   3. To check wether a file is assigned to a grib structure,
!      use the logical function 'Opened' :
!
!        if ( Opened(gribfile) ) ...
!          logical                        ::  Opened
!          type(TgribFile), intent(in)    ::  gribfile
!
!
!   4. To copy the buffer from one to another sturcture, use:
!
!        call CopySections( gribIn, gribOut, status )
!          type(TGribFile), intent(in)     ::  gribIn
!          type(TGribFile), intent(out)    ::  gribOut
!
!      to copy all sections except real section 4 (the actual data);
!      to copy the data too, use :
!
!        call CopyAllSections( gribIn, gribOut, status )
!          type(TGribFile), intent(in)     ::  gribIn
!          type(TGribFile), intent(out)    ::  gribOut
!
!   \ev
!
! !INTRODUCTION: Compilation
!
!   Compile together with the GribEx or Emos library from ECMWF:
!   \bv
!     f90 -o test.exe ... file_grib.o ... -l gribex
!   \ev
!   See also \htmladdnormallink{GRIBEX}{http://www.ecmwf.int/publications/manuals/libraries/gribex/gribexIntroduction.html}
!   at the ECMWF website.
!
! !INTRODUCTION: Grib sections
!
!   Overview of some useful entries in grib files.
!   For code tables, see:
!   \begin{itemize}
!     \item \htmladdnormallink{WMO GRIB code tables}{http://www.ecmwf.int/publications/manuals/libraries/gribex/wmoCodeTables.html}
!     \item \htmladdnormallink{ECMWF local table 2 versions}{http://www.ecmwf.int/publications/manuals/libraries/tables/tables_index.html}
!   \end{itemize}
!
!   \bv
!   Section 1 - Product Definition Section.
!   ---------------------------------------
!
!    1. Code Table 2 Version Number.               128
!    2. Originating centre identifier.              98
!    3. Model identification.                      199
!    4. Grid definition.                           255
!    5. Flag (Code Table 1)                   10000000
!
!    6. Parameter identifier (Code Table 2).
!
!    7. Type of level (Code Table 3).              109
!    8. Value 1 of level (Code Table 3).
!    9. Value 2 of level (Code Table 3).
!
!   10. Year of reference time of data.            100
!   11. Month of reference time of data.             8
!   12. Day of reference time of data.              10
!   13. Hour of reference time of data.             12
!   14. Minute of reference time of data.            0
!
!   15. Time unit (Code Table 4).                    1
!   16. Time range one.                              0
!   17. Time range two.                              0
!   18. Time range indicator (Code Table 5)          0
!
!   21. Century of reference time of data.
!
!   Section 2 - Grid Description Section.
!   -------------------------------------
!
!       Southern latitudes and Western longitudes are negative.
!       Unit is mili degrees.
!
!    1. Data represent type = lat/long     (Table 6)         0
!    2. Number of points along a parallel.                 360
!    3. Number of points along a meridian.                 181
!    4. Latitude of first grid point.                    90000
!    5. Longitude of first grid point.                 -180000
!    6. Resolution and components flag.               10000000
!    7. Latitude of last grid point.                    -90000
!    8. Longitude of last grid point.                   179000
!    9. i direction (East-West) increment.                1000
!   10. j direction (North-South) increment.              1000
!   11. Scanning mode flags (Code Table 8)            00000000
!
!    1. Data represent type = spectral     (Table 6)        50
!    2. J - Pentagonal resolution parameter.               106  <--- T
!    3. K - Pentagonal resolution parameter.               106
!    4. M - Pentagonal resolution parameter.               106
!    5. Representation type (Table 9)                        1
!    6. Representation mode (Table 10).                      2
!
!    1. Data represent type = gaussian lat/lon (Table 6)     4
!    2. Number of points along a parallel           <0=reduced, >0=regular
!    3. Number of points along a meridian.                 160
!    4. Latitude of first grid point.                    89141
!    5. Longitude of first grid point.                       0
!    6. Resolution and components flag.               00000000
!    7. Latitude of last grid point.                    -89141
!    8. Longitude of last grid point.                   358878
!    9. i direction (East-West) increment            Not given
!   10. Number of parallels between pole and equator.       80   <--- N
!   11. Scanning mode flags (Code Table 8)            00000000
!   23-22*2N. Number of points along a parallel for reduced grid
!
!   12. Number of vertical coordinate parameters.          122
!       (parameters stored in real part of section 2, starting form index 11)
!
!   Section 4 - Binary Data  Section.
!   -------------------------------------
!
!    1. Number of data values coded/decoded.             11556
!    2. Number of bits per data value.                      16
!
!    3. Type of data       (0=grid pt, 128=spectral).      128
!    4. Type of packing    (0=simple, 64=complex).          64
!    5. Type of data       (0=float, 32=integer).            0
!    6. Additional flags   (0=none, 16=present).             0
!    7. Reserved.                                            0
!    8. Number of values   (0=single, 64=matrix).            0
!    9. Secondary bit-maps (0=none, 32=present).             0
!   10. Values width       (0=constant, 16=variable).        0
!   11. Byte offset of start of packed data (N).          2214
!   12. Power (P * 1000).                                  500
!   13. Pentagonal resolution parameter J for subset.       20
!   14. Pentagonal resolution parameter K for subset.       20
!   15. Pentagonal resolution parameter M for subset.       20
!
!    3. Type of data       (0=grid pt, 128=spectral).        0
!    4. Type of packing    (0=simple, 64=complex).           0
!    5. Type of data       (0=float, 32=integer).            0
!    6. Additional flags   (0=none, 16=present).             0
!    7. Reserved.                                            0
!    8. Number of values   (0=single, 64=matrix).            0
!    9. Secondary bit-maps (0=none, 32=present).             0
!   10. Values width       (0=constant, 16=variable).        0
!
!   \ev
!
!EOI
!

module file_grib

  use GO, only : gol, goErr, goPr
  use os_specs, only : MAX_FILENAME_LEN

  implicit none

  ! --- in/out -------------------------------

  private

  public   ::  TGribFile
  public   ::  Init, Done, Opened

  public   ::  Clear
  public   ::  CopySections, CopyAllSections
  public   ::  Get
  public   ::  Check
  public   ::  Set

  public   ::  ReadRecord
  public   ::  WriteRecord

  public   ::  levtype_sfc, levtype_hyb, levtype_land

  public   ::  gridtype_ll, gridtype_gg, gridtype_sh

  public   ::  pidname, pidmax
  public   ::  pid_T
  public   ::  pid_Q
  public   ::  pid_U, pid_V, pid_W
  public   ::  pid_SP, pid_LNSP
  public   ::  pid_VO, pid_D
  public   ::  pid_Z, pid_LSM
  public   ::  pid_SR, pid_AL, pid_LSRH
  public   ::  pid_SLHF
  public   ::  pid_CLWC, pid_CIWC, pid_CC

  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'file_grib'


  ! *** Type of level (Code Table 3)

  integer, parameter  ::  levtype_sfc  = 1    ! ground or water surface
  integer, parameter  ::  levtype_hyb  = 109  ! hybrid level
  integer, parameter  ::  levtype_land = 112  ! layer below surface

  ! *** parameter identifier (Code Table 2)

  ! o numbers

  integer, parameter  ::  pid_Z    = 129  ! geopotential (orography)
  integer, parameter  ::  pid_T    = 130  ! temperature
  integer, parameter  ::  pid_U    = 131  ! u-velocity
  integer, parameter  ::  pid_V    = 132  ! v-velocity
  integer, parameter  ::  pid_Q    = 133  ! specific humidity
  integer, parameter  ::  pid_SP   = 134  ! surface pressure
  integer, parameter  ::  pid_W    = 135  ! vertical velocity

  integer, parameter  ::  pid_VO   = 138  ! vorticity

  integer, parameter  ::  pid_SLHF = 147  ! surface latent heat flux  (W m**-2 s)

  integer, parameter  ::  pid_LNSP = 152  ! ln surface pressure

  integer, parameter  ::  pid_D    = 155  ! divergence

  integer, parameter  ::  pid_LSM  = 172  ! land sea mask
  integer, parameter  ::  pid_SR   = 173  ! Surface roughness                  m
  integer, parameter  ::  pid_AL   = 174  ! Albedo                             (0-1)

  integer, parameter  ::  pid_LSRH = 234  ! Logarithm of SR length for heat

  integer, parameter  ::  pid_CLWC = 246  ! cloud liquid water content
  integer, parameter  ::  pid_CIWC = 247  ! cloud ice water content
  integer, parameter  ::  pid_CC   = 248  ! cloud cover

  ! o maximum number
  integer, parameter  ::  pidmax = 260

  ! o names

  character(len=4), parameter ::  pidname(pidmax) = (/ &
    'a   ','b   ','c   ','d   ','e   ','f   ','g   ','p008','p009','p010', &
    'p011','p012','p013','p014','p015','p016','p017','p018','p019','p020', &
    'p021','p022','p023','p024','p025','p026','p027','p028','p029','p030', &
    'p031','p032','p033','p034','p035','p036','p037','p038','p039','p040', &
    'p041','p042','p043','p044','p045','p046','p047','p048','p049','p050', &
    'p051','p052','p053','p054','p055','p056','p057','p058','p059','p060', &
    'p061','p062','p063','p064','p065','p066','p067','p068','p069','p070', &
    'p071','p072','p073','p074','p075','p076','p077','p078','p079','p080', &
    'p081','p082','p083','p084','p085','p086','p087','p088','p089','p090', &
    'p091','p092','p093','p094','p095','p096','p097','p098','p099','p100', &
    'p101','p102','p103','p104','p105','p106','p107','p108','p109','p110', &
    'p111','p112','p113','p114','p115','p116','p117','p118','p119','p120', &
    'p121','p122','p123','p124','p125','p126','p127','p128','Z   ','T   ', &
    'U   ','V   ','Q   ','SP  ','W   ','p136','p137','VO  ','p139','p140', &
    'SD  ','LSP ','CP  ','SF  ','p145','SSHF','SLHF','p148','p149','p150', &
    'p151','LNSP','p153','p154','D   ','zg  ','pw  ','pu  ','pv  ','sp  ', &
    'eu  ','du  ','dk  ','ed  ','dd  ','p166','T2M ','D2M ','p169','p170', &
    'p171','LSM ','SR  ','AL  ','p175','SSR ','p177','p178','p179','EWSS', &
    'NSSS','p182','p183','p184','p185','p186','p187','p188','p189','clb ', &
    'clt ','clfs','p193','p194','p195','p196','p197','SRC ','p199','p200', &
    'p201','p202','p203','p204','p205','p206','p207','p208','p209','p210', &
    'p211','p212','p213','p214','p215','p216','p217','p218','p219','p220', &
    'p221','p222','p223','p224','p225','p226','p227','p228','p229','p230', &
    'p231','p232','p233','LSRH','p235','p236','p237','p238','p239','p240', &
    'p241','p242','p243','p244','p245','CLWC','CIWC','CC  ','cco ','ccu ', &
    'p251','p252','p253','p254','p255','p256','p257','p258','p259','p260' /)

  ! *** data representation type (Code Table 6)

  integer, parameter  ::  gridtype_ll = 0
  integer, parameter  ::  gridtype_gg = 4
  integer, parameter  ::  gridtype_sh = 50

  ! --- types -------------------------------

  type TGribFile
    integer                ::  fu
    character(len=MAX_FILENAME_LEN) ::  fname

    ! *** buffer for grib message
    integer                ::  inbuff_bytes     ! length of input buffer
    ! Section 0 - Indicator Section.
    integer                ::  isec0(2)
    ! Section 1 - Product Definition Section.
    integer                ::  isec1(512)
    ! Section 2 - Grid Description Section.
    integer                ::  isec2(512)
    real(8)                ::  rsec2(512)
    ! Section 3 - ??
    integer                ::  isec3(2)
    real(8)                ::  rsec3(2)
    ! Section 4 - Binary Data  Section.
    integer                ::  isec4(512)
    real(8), pointer       ::  rsec4(:)
    !
    integer                ::  kword
  end type TGribFile


  ! --- interfaces ---------------------------

  interface Init
    module procedure grib_Init
  end interface

  interface Done
    module procedure grib_Done
  end interface

  interface Opened
    module procedure grib_Opened
  end interface

  interface ReadRecord
    module procedure grib_ReadRecord
  end interface

  interface WriteRecord
    module procedure grib_WriteRecord
  end interface

  interface Clear
    module procedure grib_Clear
  end interface

  interface CopySections
    module procedure grib_CopySections
  end interface

  interface CopyAllSections
    module procedure grib_CopyAllSections
  end interface

  interface Get
    module procedure grib_Get
  end interface

  interface Check
    module procedure grib_Check
  end interface

  interface Set
    module procedure grib_Set
  end interface





contains


  ! =================================================================

  ! Open gribfile.
  !
  ! USAGE
  !   call Init( gribfile, 'input.gb', 'r'|'w'  )
  !
  ! DESCRIPTION
  !   Interface around routine 'pbOpen'.
  !   In 'gribfile', space is allocated to store grib sections.
  !


  subroutine grib_Init( F, file, mode, status )

    ! --- in/out -----------------------

    type(TGribFile), intent(out)              ::  F
    character(len=*), intent(in)              ::  file
    character(len=*), intent(in)              ::  mode
    integer, intent(out)                      ::  status

    ! --- const ------------------------

    character(len=*), parameter  ::  rname = mname//'/grib_Init'

    ! --- local ------------------------

    ! dummy buffer of insufficient length:
    integer        ::  inbuff(10)

    ! --- begin ------------------------

    ! *** initialize grib sections

    ! Section 0 - Indicator Section.
    F%isec0 = 0

    ! Section 1 - Product Definition Section.
    F%isec1 = 0

    ! Section 2 - Grid Description Section.
    F%isec2 = 0
    F%rsec2 = 0.0
    ! Section 3 - ??
    F%isec3 = 0
    F%rsec3 = 0.0

    ! Section 4 - Binary Data  Section.
    F%isec4 = 0
    nullify( F%rsec4 )

    F%kword = 0

    ! *** Open specified file with requested name.

    ! A free file unit seems to be assigned in pbOpen.
    call pbOpen( F%fu, file, mode, status )
    if ( status /= 0 ) then
      write (gol,'("from pbOpen:")'); call goErr
      select case ( status )
        case (-1)
          write (gol,'("  could not open file")'); call goErr
        case (-2)
          write (gol,'("  invalid file name")') ; call goErr
        case (-3)
          write (gol,'("  invalid open mode")'); call goErr
        case default
          write (gol,'("  unknown return status: ",i6)') status; call goErr
      end select
      write (gol,'("  file: ",a)') trim(file); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! save filename
    F%fname = file

    ! *** special tasks

    select case ( mode )

      case ( 'r' )

        ! Determine length of records ('messages' in grib context);
        ! open with insufficient buffer, and check error messages:

        status = 1   ! force return even if error occures
        call pbGrib( F%fu, inbuff, size(inbuff)*kind(inbuff), F%inbuff_bytes, status )
        select case ( status )
          case ( 0)
            ! no error ?
            write (gol,'("tried to read record into buffer with insuficient length,")'); call goErr
            write (gol,'("but no error occured ...")'); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
          case (-1)
            write (gol,'("end-of-file is hit before a GRIB product is read")'); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
          case (-3)
            ! size of input buffer was not sufficient (as expected)
          case default
            write (gol,'("unknown return status from pbGrib : ",i6)') status; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
        ! adhoc increment (sometimes first record is very different form others)
        F%inbuff_bytes = F%inbuff_bytes*2
        ! Close file:
        call pbClose( F%fu, status )
        if ( status /= 0 ) then
          write (gol,'("error in call to pbClose; status=",i6)') status; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if

        ! Reopen file:
        call pbOpen( F%fu, file, mode, status )
        if ( status /= 0 ) then
          write (gol,'("from second open; status=",i6)') status; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if

      case ( 'w' )

        ! nothing special

      case default

        write (gol,'(" unknown mode `",a,"`")') mode; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return

    end select

  end subroutine grib_Init


  ! ===


  ! Close gribfile, clear buffers
  !
  ! USAGE
  !   call Done( gribfile )
  !
  ! DESCRIPTION
  !   Interface around routine 'pbClose'.
  !

  subroutine grib_Done( F, status )

    ! --- in/out -----------------------

    type(TGribFile), intent(inout)  ::  F
    integer, intent(inout)          ::  status

    ! --- const ------------------------

    character(len=*), parameter  ::  rname = mname//'/grib_Done'

    ! --- begin ------------------------

    ! close file:
    call pbClose( F%fu, status )
    if ( status /= 0 ) then
      write (gol,'("from closing grib file:")'); call goErr
      write (gol,'("  status : ",i6)') status; call goErr
      write (gol,'("  file   : ",a)') trim(F%fname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! clear section 4 :
    if (associated(F%rsec4)) deallocate( F%rsec4 )

  end subroutine grib_Done


  ! ===


  ! Grib file opened ?

  logical function grib_Opened( gribfile )

    ! --- in/out ------------------------------

    type(TGribFile), intent(in)     ::  gribfile

    ! --- begin -------------------------

    inquire( unit=gribfile%fu, opened=grib_Opened )

  end function grib_Opened


  ! =============================================================


  subroutine grib_CopySections( gribIn, gribOut, status )

    ! ---in/out ------------------------------

    type(TGribFile), intent(in)     ::  gribIn
    type(TGribFile), intent(out)    ::  gribOut
    integer, intent(inout)          ::  status

    ! --- begin -------------------------

    !gribOut%isec0 = gribIn%isec0

    gribOut%isec1 = gribIn%isec1

    gribOut%isec2 = gribIn%isec2
    gribOut%rsec2 = gribIn%rsec2

    gribOut%isec3 = gribIn%isec3

    gribOut%isec4 = gribIn%isec4

    ! ok
    status = 0

  end subroutine grib_CopySections


  ! ***


  subroutine grib_CopyAllSections( gribIn, gribOut, status )

    ! ---in/out ------------------------------

    type(TGribFile), intent(in)     ::  gribIn
    type(TGribFile), intent(out)    ::  gribOut
    integer, intent(inout)          ::  status

    ! --- begin -------------------------

    !gribOut%isec0 = gribIn%isec0

    gribOut%isec1 = gribIn%isec1

    gribOut%isec2 = gribIn%isec2
    gribOut%rsec2 = gribIn%rsec2

    gribOut%isec3 = gribIn%isec3

    gribOut%isec4 = gribIn%isec4
    if ( associated(gribOut%rsec4) ) then
      if ( size(gribOut%rsec4) /= size(gribIn%rsec4) ) then
        deallocate( gribOut%rsec4 )
        allocate( gribOut%rsec4(size(gribIn%rsec4)) )
      end if
    else
      allocate( gribOut%rsec4(size(gribIn%rsec4)) )
    end if
    gribOut%rsec4 = gribIn%rsec4

    ! ok
    status = 0

  end subroutine grib_CopyAllSections


  ! ***


  subroutine grib_Clear( F, status )

    ! --- in/out -----------------------

    type(TGribFile), intent(out)    ::  F
    integer, intent(inout)          ::  status

    ! --- begin ------------------------

    ! Section 0 - Indicator Section.
    F%isec0 = 0

    ! Section 1 - Product Definition Section.
    F%isec1 = 0

    ! Section 2 - Grid Description Section.
    F%isec2 = 0
    F%rsec2 = 0.0

    ! Section 3 - ??
    F%isec3 = 0
    F%rsec3 = 0.0

    ! Section 4 - Binary Data  Section.
    F%isec4 = 0
    if ( associated(F%rsec4) ) F%rsec4 = 0.0

    ! default number of bits
    F%isec4(2) = 24

    F%kword = 0

    ! ok
    status = 0

  end subroutine grib_Clear


  ! ==================================================================

  ! USAGE
  !   call ReadRecord( gribfile, status )
  !
  ! DESCRIPTION
  !   Reads next grib record in buffers.
  !
  ! RETURN STATUS
  !    0 : no error
  !    1 : eof
  !    2 : some error
  !   Execution might stop if other errors are encountered.
  !


  subroutine grib_ReadRecord( F, status )

    ! --- in/out -------------------------

    type(TGribFile), intent(inout)   ::  F
    integer, intent(inout)           ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/grib_ReadRecord'

    ! --- local --------------------------

    integer, allocatable   ::  inbuff(:)
    integer                ::  length, len4
    logical              ::  verbose

    ! --- begin --------------------------

    ! write error messages ?
    verbose = status == 0

    ! allocate array to store packed data:
    allocate( inbuff(ceiling(F%inbuff_bytes*1.0/kind(inbuff))) )

    ! read record:
    call pbGrib( F%fu, inbuff, size(inbuff)*kind(inbuff), length, status )
    select case ( status )
      case ( 0)
        ! no error !
      case (-1)
        if ( verbose ) then
          write (gol,'("end-of-file is hit before a GRIB product is read")'); call goErr
          write (gol,'("  grib file : ",a)') F%fname; call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=1; return
      case (-3)
        ! size of input buffer is not sufficient for the GRIB product !
        if ( verbose ) then
          write (gol,'("size of input buffer is not sufficient for the GRIB product:")'); call goErr
          write (gol,'("  current size : ",i8)') size(inbuff)*kind(inbuff); call goErr
          write (gol,'("  required     : ",i8)') length; call goErr
          write (gol,'("  grib file    : ",a)') F%fname; call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=2; return
      case default
        if ( verbose ) then
          write (gol,'("unknown return status from pbGrib : ",i6)') status; call goErr
          write (gol,'("  grib file : ",a)') F%fname; call goErr
          write (gol,'("in ",a)') rname; call goErr
        endif
        status=2; return
    end select

    ! Decode.
    ! The array to hold the unpacked data may need to be 4 times
    ! as long as that for the packed data (or more!).
    len4 = ceiling(4*F%inbuff_bytes*1.0/kind(F%rsec4))
    if ( associated(F%rsec4) ) then
      if ( size(F%rsec4) < len4 ) then
        deallocate( F%rsec4 )
        allocate( F%rsec4(len4) )
      end if
    else
      allocate( F%rsec4(len4) )
    end if
    call GribEx(F%isec0,F%isec1,F%isec2,F%rsec2,F%isec3,F%rsec3,&
                F%isec4, F%rsec4, size(F%rsec4)*kind(F%rsec4), &
                inbuff, size(inbuff), F%kword, 'D' , status )
    if ( status /= 0 ) then
      select case ( status )
        case (-6)
          if ( verbose ) then
            write (gol,'("in call to GribEx : ")'); call goErr
            write (gol,'("  found pseudo data in grib file ")'); call goErr
            write (*,'("ERROR in ",a)') rname; call goErr
          end if
          status=2; return
        case default
          if ( verbose ) then
            write (gol,'("unknown return status from GribEx : ",i6)') status; call goErr
            write (gol,'("in ",a)') rname
          end if
          status=2; return
      end select
    end if

    ! clear input buffer:
    deallocate( inbuff )

    ! ok
    status = 0

  end subroutine grib_ReadRecord


  ! ***


  ! encode and write grib sections


  subroutine grib_WriteRecord( F, status )

    ! --- in/out --------------------

    type(TGribFile), intent(in)        ::  F
    integer, intent(out)               ::  status

    ! --- const ---------------------

    character(len=*), parameter  ::  rname = mname//'/grib_WriteRecord'

    ! --- local ---------------------

    integer, allocatable   ::  outbuff(:)

    ! --- begin ---------------------

    ! *** encode

    status = 1
    allocate( outbuff(size(F%rsec4)) )
    call GribEx(F%isec0,F%isec1,F%isec2,F%rsec2,F%isec3,F%rsec3,&
                F%isec4, F%rsec4, size(F%rsec4)*kind(F%rsec4), &
                outbuff, size(outbuff), F%kword, 'C' , status )
    select case ( status )
      case (-2)
        write (gol,'("WARNING: bitmap was encountered with all bits set to 1")'); call goPr
      case (-3)
        write (gol,'("WARNING: predefined bitmap was encountered.")'); call goPr
      case (-4)
        write (gol,'("WARNING: bitmap was encountered.")'); call goPr
      case (-5)
        write (gol,'("WARNING: bitmap was encountered and the data has not been decoded.")'); call goPr
      case (-6)
        write (gol,'("WARNING: ECMWF pseudo-GRIB data (TIDE or BUDG) has been encountered.")'); call goPr
      case ( 0 )
        ! ok
      case default
        write (gol,'("unknown status from GribEx : ",i6)') status
        write (gol,'("in ",a)') rname; status=1; return
    end select

    ! *** write GRIB record

    ! PBWRITE: Writes a block of bytes to a file.
    !   status :    -1 = Could not write to file.
    !            >=  0 = Number of bytes written.
    !
    call pbWrite( F%fu, outbuff, F%isec0(1), status )
    if ( status == -1 ) then
      write (gol,'("could not write to file")'); call goErr
      write (gol,'("in ",a)') rname; status=1; return
    else if ( status >= 0 ) then
      ! bytes written correctly
    else
      write (gol,'("unknown status from GribEx : ",i6)') status
      write (gol,'("in ",a)') rname; status=1; return
    end if

    ! done
    deallocate( outbuff )

    ! ok
    status = 0

  end subroutine grib_WriteRecord


  ! ============================================================

  ! Extract data from grib sections.
  ! Optional arguments are passed to read parameter id, date, etc.
  ! If the same option is preceded by 'check_', an error
  ! message is written if the grib file does not match
  !
  !  pid   :  parameter identifier (integer number)
  !
  !  level :  index of hybride level
  !
  !  reftime   :  5 element integer array : yy, mm, dd, hh, min
  !  timerange :  4 element ingeger array : unit, val1, val2, indicator
  !

  subroutine grib_Get( F, status, &
                           model_id, pid, &
                           levtype, level, nlev, hyb_a, hyb_b, &
                           reftime, timerange, &
                           gridtype, &
                           ll, &
                           lon_first, lon_last, lon_inc, lon_n, &
                           lat_first, lat_last, lat_inc, lat_n, &
                           T, sh, &
                           N, gg  )

    ! --- in/out -------------------------

    type(TGribFile), intent(in)      ::  F
    integer, intent(out)             ::  status

    integer, intent(out), optional   ::  model_id
    integer, intent(out), optional   ::  pid
    integer, intent(out), optional   ::  levtype, level, nlev
    real, intent(out), optional      ::  hyb_a(:), hyb_b(:)
    integer, intent(out), optional   ::  reftime(5)
    integer, intent(out), optional   ::  timerange(4)
    integer, intent(out), optional   ::  gridtype
    real, intent(out), optional      ::  lon_first, lon_last, lon_inc
    integer, intent(out), optional   ::  lon_n
    real, intent(out), optional      ::  lat_first, lat_last, lat_inc
    integer, intent(out), optional   ::  lat_n
    real, intent(out), optional      ::  ll(:,:)
    integer, intent(out), optional   ::  T
    complex, intent(out), optional   ::  sh(:)
    integer, intent(out), optional   ::  N
    real, intent(out), optional      ::  gg(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/grib_Get'

    ! --- local --------------------------

    integer                ::  grib_reftime(5)
    integer                ::  grib_timerange(4)
    integer                ::  base
    integer                ::  grib_n

    real, allocatable      ::  rsec4_adhoc(:)

    ! --- begin --------------------------

    ! -- model id

    if ( present( model_id ) ) model_id = F%isec1(3)

    ! --  parameter

    if ( present( pid ) ) pid = F%isec1(6)

    ! --  level

    if ( present(levtype) ) levtype = F%isec1(7)
    if ( present(level  ) ) level   = F%isec1(8)

    ! Vertical Coordinate Parameters strored in rsec2, start in element 11
    if ( present(hyb_a) .or. present(hyb_b) ) then
      if ( .not. (present(hyb_a) .and. present(hyb_b)) ) then
        write (gol,'("extract both hybride params or none")'); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      if ( F%isec2(12) /= size(hyb_a) + size(hyb_b) ) then
        write (gol,'("numbers of hybride parameters do not match:")'); call goErr
        write (gol,'("  expected : ",i6)') size(hyb_a)+size(hyb_b); call goErr
        write (gol,'("  found    : ",i6)') F%isec2(12); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      base = 10            ;  hyb_a = F%rsec2(base+1:base+size(hyb_a))
      base = 10+size(hyb_a);  hyb_b = F%rsec2(base+1:base+size(hyb_b))
    end if

    ! number of levels
    if ( present(nlev) ) nlev = F%isec2(12)/2 - 1


    ! --  reference time

    ! extract reference time:
    grib_reftime(1) = (F%isec1(21)-1)*100 + F%isec1(10)   ! year incl. century
    grib_reftime(2) = F%isec1(11)            ! month
    grib_reftime(3) = F%isec1(12)            ! day
    grib_reftime(4) = F%isec1(13)            ! hour
    grib_reftime(5) = F%isec1(14)            ! minutes

    if ( present(reftime) ) reftime = grib_reftime

    ! --  time range

    grib_timerange(1:4) = F%isec1(15:18)

    if ( present(timerange) ) timerange = grib_timerange

    ! ---  grid

    if ( present(gridtype)  ) gridtype  = F%isec2( 1)
    if ( present(lon_n)     ) lon_n     = F%isec2( 2)
    if ( present(lat_n)     ) lat_n     = F%isec2( 3)
    if ( present(lat_first) ) lat_first = F%isec2( 4) / 1000.0
    if ( present(lon_first) ) lon_first = F%isec2( 5) / 1000.0
    if ( present(lat_last)  ) lat_last  = F%isec2( 7) / 1000.0
    if ( present(lon_last)  ) lon_last  = F%isec2( 8) / 1000.0
    if ( present(lon_inc)   ) lon_inc   = F%isec2( 9) / 1000.0
    if ( present(lat_inc)   ) lat_inc   = F%isec2(10) / 1000.0

    ! -- lat/lon grid

    if ( present(ll) ) then
      !call Check( F, gridtype=gridtype_ll )
      ! >>> no recursive call >>>>>>>>>>>>>>
      grib_n = F%isec2(1)
      if ( grib_n /= gridtype_ll ) then
        write (gol,'("grid types do not match:")'); call goErr
        write (gol,'("  expected : ",i4," (ll)")') gridtype_ll; call goErr
        write (gol,'("  found    : ",i4)') grib_n; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      if ( (size(ll,1)/=F%isec2(2)) .or. (size(ll,2)/=F%isec2(3)) ) then
        write (gol,'("grid sizes do not match:")'); call goErr
        write (gol,'("  expected : ",i4," x ",i4)') size(ll,1), size(ll,2); call goErr
        write (gol,'("  found    : ",i4," x ",i4)') F%isec2(2), F%isec2(3); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      if ( .not. associated(F%rsec4) ) then
        write (gol,'("no grib section 4 read")'); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ll = reshape( F%rsec4(1:size(ll)), shape(ll) )
    end if

    ! -- gaussian lat/lon grid

    if ( present(N) ) then
      !call Check( F, gridtype=gridtype_gg )
      ! >>> no recursive call >>>>>>>>>>>>>>
      grib_n = F%isec2(1)
      if ( grib_n /= gridtype_gg ) then
        write (gol,'("grid types do not match:")'); call goErr
        write (gol,'("  expected : ",i4," (gg)")') gridtype_gg; call goErr
        write (gol,'("  found    : ",i4)') grib_n; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      N = F%isec2(10)
    end if

    if ( present(gg) ) then
      !call Check( F, gridtype=gridtype_gg )
      ! >>> no recursive call >>>>>>>>>>>>>>
      grib_n = F%isec2(1)
      if ( grib_n /= gridtype_gg ) then
        write (gol,'("grid typess do not match:")'); call goErr
        write (gol,'("  expected : ",i4," (gg)")') gridtype_gg; call goErr
        write (gol,'("  found    : ",i4)') grib_n; call goErr
        write (gol,'("in ",a)') rname; call goErr
        status=1; return
      end if
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      if ( size(gg) /= F%isec4(1) ) then
        write (gol,'("gg grid sizes do not match:")'); call goErr
        write (gol,'("  expected : ",i4)') size(gg); call goErr
        write (gol,'("  found    : ",i4)') F%isec4(1); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! data present ?
      if ( .not. associated(F%rsec4) ) then
        write (gol,'("no grib section 4 read")'); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      gg = F%rsec4(1:size(gg))
    end if

    ! -- spectral coef

    if ( present(T) ) then
      !call Check( F, gridtype=gridtype_sh )
      ! >>> no recursive call >>>>>>>>>>>>>>
      grib_n = F%isec2(1)
      if ( grib_n /= gridtype_sh ) then
        write (gol,'("grid types do not match:")'); call goErr
        write (gol,'("  expected : ",i4," (sh)")') gridtype_sh; call goErr
        write (gol,'("  found    : ",i4)') grib_n; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      T = F%isec2(4)
    end if

    if ( present(sh) ) then
      !call Check( F, gridtype=gridtype_sh )
      ! >>> no recursive call >>>>>>>>>>>>>>
      grib_n = F%isec2(1)
      if ( grib_n /= gridtype_sh ) then
        write (gol,'("grid typess do not match:")'); call goErr
        write (gol,'("  expected : ",i4," (sh)")') gridtype_sh; call goErr
        write (gol,'("  found    : ",i4)') grib_n; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! check number of elements; complex is 2 elements !
      if ( size(sh) /= F%isec4(1)/2 ) then
        write (gol,'("numbers of spectral coefficients do not match:")'); call goErr
        write (gol,'("  expected : ",i6)') size(sh); call goErr
        write (gol,'("  found    : ",i6)') F%isec4(1)/2; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
      ! convert from real to complex:
      ! >>> NOTE: this call lets 'f90 -C' fail on SUN compiler !
      !sh = transfer( F%rsec4(1:F%isec4(1)), (/(0.0,0.0)/) )
      ! <<< fix >>>
      allocate( rsec4_adhoc(F%isec4(1)) )
      rsec4_adhoc = F%rsec4(1:F%isec4(1))
      sh = transfer( rsec4_adhoc, (/(0.0,0.0)/) )
      deallocate( rsec4_adhoc )
      ! <<<
    end if

    ! ok
    status = 0

  end subroutine grib_Get


  ! ===


  ! NOTE
  !   Do not check hybride coeff; different number of bits etc
  !   cause small differences ...

  subroutine grib_Check( F, status, debug, &
                             model_id, pid, &
                             levtype, level, &
                             reftime, timerange, &
                             gridtype, &
                             lon_first, lon_last, lon_inc, lon_n, &
                             lat_first, lat_last, lat_inc, lat_n, &
                             T  )

    ! --- in/out -------------------------

    type(TGribFile), intent(in)      ::  F
    integer, intent(out)             ::  status

    integer, intent(in), optional    ::  debug

    integer, intent(in), optional    ::  model_id
    integer, intent(in), optional    ::  pid
    integer, intent(in), optional    ::  levtype, level
    !real, intent(in), optional       ::  hyb_a(:), hyb_b(:)
    integer, intent(in), optional    ::  reftime(5)
    integer, intent(in), optional    ::  timerange(4)
    integer, intent(in), optional    ::  gridtype
    integer, intent(in), optional    ::  lon_first, lon_last, lon_inc, lon_n
    integer, intent(in), optional    ::  lat_first, lat_last, lat_inc, lat_n
    integer, intent(in), optional    ::  T

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/grib_Check'

    ! --- local --------------------------

    logical                ::  match

    integer                ::  grib_model_id
    integer                ::  grib_pid
    integer                ::  grib_levtype, grib_level
    integer                ::  grib_reftime(5)
    integer                ::  grib_timerange(4)
    integer                ::  grib_n
    integer                ::  grib_T

    logical                ::  verbose

    ! --- begin --------------------------

    ! write error messages ?
    verbose = .false.
    if ( present(debug) ) verbose = debug > 0

    ! ok by default
    status = 0

    ! -- check model id

    if ( present( model_id ) ) then
      call Get( F, status, model_id=grib_model_id )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( grib_model_id /= model_id ) then
        if ( verbose ) then
          write (gol,'("model id''s do not match:")'); call goErr
          write (gol,'("  expected : ",i4)') model_id; call goErr
          write (gol,'("  found    : ",i4)') grib_model_id; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if


    ! -- check parameter

    if ( present( pid ) ) then
      call Get( F, status, pid=grib_pid )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( grib_pid /= pid ) then
        if ( verbose ) then
          write (gol,'("parameter id''s do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') pid; call goErr
          write (gol,'("  found    : ",i6)') grib_pid; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if


    ! -- check level

    if ( present(levtype) ) then
      call Get( F, status, levtype=grib_levtype )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( grib_levtype /= levtype ) then
        if ( verbose ) then
          write (gol,'("level types do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') levtype; call goErr
          write (gol,'("  found    : ",i6)') grib_levtype; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    if ( present(level) ) then
      call Get( F, status, level=grib_level )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( grib_level /= level ) then
        if ( verbose ) then
          write (gol,'("levels do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') level; call goErr
          write (gol,'("  found    : ",i6)') grib_level; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    ! -- check reference time

    ! extract reference time:
    grib_reftime(1) = (F%isec1(21)-1)*100 + F%isec1(10)   ! year incl. century
    grib_reftime(2) = F%isec1(11)            ! month
    grib_reftime(3) = F%isec1(12)            ! day
    grib_reftime(4) = F%isec1(13)            ! hour
    grib_reftime(5) = F%isec1(14)            ! minutes

    if ( present(reftime) ) then
      if ( any( grib_reftime /= reftime ) ) then
        if ( verbose ) then
          write (gol,'("reference times do not match:")'); call goErr
          write (gol,'("  expected : ",i4,4i3)') reftime; call goErr
          write (gol,'("  found    : ",i4,4i3)') grib_reftime; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if


    ! -- check time range

    grib_timerange(1:4) = F%isec1(15:18)
    if ( present(timerange) ) then
      if ( any( grib_timerange /= timerange ) ) then
        if ( verbose ) then
          write (gol,'("time ranges do not match:")'); call goErr
          write (gol,'("  expected : ",4i3)') timerange; call goErr
          write (gol,'("  found    : ",4i3)') grib_timerange; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    ! --- check grid

    grib_n = F%isec2(1)
    if ( present(gridtype) ) then
      if ( grib_n /= gridtype ) then
        if ( verbose ) then
          write (gol,'("gridtypes do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') gridtype; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(2)
    if ( present(lon_n) ) then
      if ( grib_n /= lon_n ) then
        if ( verbose ) then
          write (gol,'("number of longitudes do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lon_n; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(3)
    if ( present(lat_n) ) then
      if ( grib_n /= lat_n ) then
        if ( verbose ) then
          write (gol,'("number of latitudes do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lat_n; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(4)
    if ( present(lat_first) ) then
      if ( grib_n /= lat_first ) then
        if ( verbose ) then
          write (gol,'("first latitudes do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lat_first; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(5)
    if ( present(lon_first) ) then
      if ( grib_n /= lon_first ) then
        if ( verbose ) then
          write (gol,'("first longitudes do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lon_first; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(7)
    if ( present(lat_last) ) then
      if ( grib_n /= lat_last ) then
        if ( verbose ) then
          write (gol,'("last latitudes do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lat_last; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(8)
    if ( present(lon_last) ) then
      if ( grib_n /= lon_last ) then
        if ( verbose ) then
          write (gol,'("last longitudes do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lon_last; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(9)
    if ( present(lon_inc) ) then
      if ( grib_n /= lon_inc ) then
        if ( verbose ) then
          write (gol,'("longitude increments do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lon_inc; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    grib_n = F%isec2(10)
    if ( present(lat_inc) ) then
      if ( grib_n /= lat_inc ) then
        if ( verbose ) then
          write (gol,'("latitude increments do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') lat_inc; call goErr
          write (gol,'("  found    : ",i6)') grib_n; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
        end if
        status=status+1
      end if
    end if

    ! --- check truncation

    if ( present(T) ) then
      grib_T = F%isec2(4)
      if ( grib_T /= T ) then
        if ( verbose ) then
          write (gol,'("spectral truncations do not match:")'); call goErr
          write (gol,'("  expected : ",i6)') T; call goErr
          write (gol,'("  found    : ",i6)') grib_T; call goErr
          write (gol,'("  in ",a)') trim(F%fname); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=status+1
      end if
    end if

    ! return with status equal to number of failures

  end subroutine grib_Check


  ! ===


  ! fill gribsections

  subroutine grib_Set( F, status, &
                           model_id, pid, &
                           levtype, hyb_a, hyb_b, level, &
                           reftime, timerange, &
                           gridtype, &
                           ll, &
                           lon_first, lon_inc, lon_last, lon_n, &
                           lat_first, lat_inc, lat_last, lat_n, &
                           scanning_mode, nbits )

    ! --- in/out ----------------------

    type(TGribFile), intent(inout)       ::  F
    integer, intent(out)                 ::  status

    integer, intent(in), optional        ::  model_id

    integer, intent(in), optional        ::  pid

    integer, intent(in), optional        ::  levtype
    integer, intent(in), optional        ::  level
    real, intent(in), optional           ::  hyb_a(:), hyb_b(:)

    integer, intent(in), optional        ::  reftime(5)   ! (/yy,mm,dd,hh,min/)
    integer, intent(in), optional        ::  timerange(4)

    integer, intent(in), optional        ::  gridtype

    real, intent(in), optional           ::  ll(:,:)
    integer, intent(in), optional        ::  lon_n, lat_n
    integer, intent(in), optional        ::  lon_first, lon_inc, lon_last  ! mili degree
    integer, intent(in), optional        ::  lat_first, lat_inc, lat_last  ! mili degree
    integer, intent(in), optional        ::  scanning_mode

    integer, intent(in), optional        ::  nbits

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/grib_Set'

    ! --- local ------------------------

    integer           :: base

    ! --- begin ------------------------

    ! Section 1 - Product Definition Section.
    ! ---------------------------------------

    ! Code Table 2 Version Number.
    F%isec1(1) = 128

    ! Originating centre identifier.
    F%isec1(2) = 98    ! ECMWF

    ! Model identification.
    if ( present(model_id) ) F%isec1(3) = model_id

    ! Grid definition.
    F%isec1(4) = 255   ! ???

    ! Flag (Code Table 1)                   10000000
    F%isec1(5) = 128   ! ???

    ! Parameter identifier (Code Table 2, NO LOCAL TABLE FOR KNMI)
    if ( present(pid) ) F%isec1(6) = pid

    ! Type of level (Code Table 3).              109   ; hybride
    ! Value 1 of level (Code Table 3).             1
    ! Value 2 of level (Code Table 3).             0
    if ( present(levtype) ) F%isec1(7) = levtype
    if ( present(level    ) ) F%isec1(8) = level
    F%isec1(9) = 0

    ! 10. Year of reference time of data.
    ! 11. Month of reference time of data.
    ! 12. Day of reference time of data.
    ! 13. Hour of reference time of data.
    ! 14. Minute of reference time of data.
    ! 21. Century of reference time of data.
    if ( present(reftime) ) then
      F%isec1(21)    = int(reftime(1)/100.0) + 1
      F%isec1(10)    = mod(reftime(1),100)
      F%isec1(11:14) = reftime(2:5)
    end if

    ! Time unit (Code Table 4).
    ! Time range one.
    ! Time range two.
    ! Time range indicator (Code Table 5)
    if ( present(timerange) ) F%isec1(15:18) = timerange

    ! Number averaged.
    F%isec1(19) = 0

    ! Number missing from average.
    F%isec1(20) = 0


    ! Section 2 - Grid Description Section.
    ! -------------------------------------
    ! (Southern latitudes and Western longitudes are negative.)

    !  1. Data represent type (Table 6)         0 = lat/long
    if ( present(gridtype) ) F%isec2(1) = gridtype

    !  2. Number of points along a parallel.                 144
    !  3. Number of points along a meridian.                  72
    !  4. Latitude of first grid point.                   -88750  mdeg
    !  5. Longitude of first grid point.                 -180000  mdeg
    !  6. Resolution and components flag.               10000000
    !  7. Latitude of last grid point.                     88750
    !  8. Longitude of last grid point.                   177500
    !  9. i direction (East-West) increment.                2500
    ! 10. j direction (North-South) increment.              2500

    if ( present(lon_n) .or. present(lon_first) .or. present(lon_last) .or. present(lon_inc) ) then

      if ( present(lon_n)     ) F%isec2(2) = lon_n
      if ( present(lon_first) ) F%isec2(5) = lon_first
      if ( present(lon_last)  ) F%isec2(8) = lon_last
      if ( present(lon_inc)   ) F%isec2(9) = lon_inc

      if ( present(lon_first) .and. present(lon_last) .and. present(lon_inc) ) then
        ! lon_n
        F%isec2(2) = (lon_last-lon_first)/lon_inc + 1
      else if ( present(lon_n) .and. present(lon_first) .and. present(lon_inc) ) then
        ! lon_last
        F%isec2(8) = lon_first + lon_inc*(lon_n-1)
      else if ( present(lon_n) .and. present(lon_first) .and. present(lon_last) ) then
        ! lon_inc
        F%isec2(9) = (lon_last-lon_first)/lon_n
      else
        write (gol,'("expected at least 3 of lon_n,lon_first,lon_last,lon_inc")'); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

    end if

    if ( present(lat_n) .or. present(lat_first) .or. present(lat_last) .or. present(lat_inc) ) then

      if ( present(lat_n)     ) F%isec2(3)  = lat_n
      if ( present(lat_first) ) F%isec2(4)  = lat_first
      if ( present(lat_last)  ) F%isec2(7)  = lat_last
      if ( present(lat_inc)   ) F%isec2(10) = lat_inc

      if ( present(lat_first) .and. present(lat_last) .and. present(lat_inc) ) then
        ! lat_n
        F%isec2(3) = (lat_last-lat_first)/lat_inc + 1
      else if ( present(lat_n) .and. present(lat_first) .and. present(lat_inc) ) then
        ! lat_last
        F%isec2(7) = lat_first + lat_inc*(lat_n-1)
      else if ( present(lat_n) .and. present(lat_first) .and. present(lat_last) ) then
        ! lat_inc
        F%isec2(10) = (lat_last-lat_first)/lat_n
      else
        write (gol,'("expected at least 3 of lon_n,lon_first,lon_last,lon_inc")'); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

    end if

    F%isec2(6) = 128

    ! 11. Scanning mode flags (Code Table 8)
    !       bit 1 : 0 = west to east   ,    1 = east to west
    !       bit 2 : 0 = north to south ,    1 = south to north
    !       bit 3 : 0 = store lon rows ,    1 = store lat columns
    ! default scanning mode 0
    !F%isec2(11) = 0
    !if ( F%isec2(8) < F%isec2(5) ) F%isec2(11) = F%isec2(11) + 4  ! east to west
    !if ( F%isec2(7) > F%isec2(4) ) F%isec2(11) = F%isec2(11) + 2  ! south to north
    ! PROBLEM >>> GribEx seems to accept scanning mode 0 only ..
    F%isec2(11) = 0
    ! <<<
    if ( present(scanning_mode) ) F%isec2(11) = scanning_mode

    ! 12. Number of vertical coordinate parameters.
    ! Vertical Coordinate Parameters strored in rsec2, start in element 11
    if ( present(hyb_a) .and. present(hyb_b) ) then
      F%isec2(12) = size(hyb_a) + size(hyb_b)
      base = 10            ;  F%rsec2(base+1:base+size(hyb_a)) = hyb_a
      base = 10+size(hyb_a);  F%rsec2(base+1:base+size(hyb_b)) = hyb_b
    end if

    ! set rest to zero for safety
    F%isec2(13:size(F%isec2)) = 0


    ! Section 4 - Binary Data  Section.
    ! -------------------------------------

    !  1. Number of data values coded/decoded.             11556
    !  2. Number of bits per data value.                      16
    !
    !  3. Type of data       (0=grid pt, 128=spectral).      128
    !  4. Type of packing    (0=simple, 64=complex).          64
    !  5. Type of data       (0=float, 32=integer).            0
    !  6. Additional flags   (0=none, 16=present).             0
    !  7. Reserved.                                            0
    !  8. Number of values   (0=single, 64=matrix).            0
    !  9. Secondary bit-maps (0=none, 32=present).             0
    ! 10. Values width       (0=constant, 16=variable).        0
    ! 11. Byte offset of start of packed data (N).          2214
    ! 12. Power (P * 1000).                                  500
    ! 13. Pentagonal resolution parameter J for subset.       20
    ! 14. Pentagonal resolution parameter K for subset.       20
    ! 15. Pentagonal resolution parameter M for subset.       20


    ! set lat/lon grid:
    !  3. Type of data       (0=grid pt, 128=spectral).        0
    !  4. Type of packing    (0=simple, 64=complex).           0
    !  5. Type of data       (0=float, 32=integer).            0
    !  6. Additional flags   (0=none, 16=present).             0
    !  7. Reserved.                                            0
    !  8. Number of values   (0=single, 64=matrix).            0
    !  9. Secondary bit-maps (0=none, 32=present).             0
    ! 10. Values width       (0=constant, 16=variable).        0
    if ( present(ll) ) then

      ! check grid size:
      if ( F%isec2(2) /= size(ll,1) .or. F%isec2(3) /= size(ll,2) ) then
        write (gol,'("grid sizes do not match:")'); call goErr
        write (gol,'("  expected : ",i4," x ",i4)') size(ll,1), size(ll,2); call goErr
        write (gol,'("  found    : ",i4," x ",i4)') F%isec2(2), F%isec2(3); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if

      !  1. Number of data values coded/decoded.             11556
      F%isec4(1) = size(ll)

      !  3. Type of data       (0=grid pt, 128=spectral).        0
      !  4. Type of packing    (0=simple, 64=complex).           0
      !  5. Type of data       (0=float, 32=integer).            0
      !  6. Additional flags   (0=none, 16=present).             0
      !  7. Reserved.                                            0
      !  8. Number of values   (0=single, 64=matrix).            0
      !  9. Secondary bit-maps (0=none, 32=present).             0
      ! 10. Values width       (0=constant, 16=variable).        0
      F%isec4(3:10) = 0

      ! allocate space to store field
      ! NOTE: allocate extra, sometimes errors if fit exately!
      if ( associated(F%rsec4) ) then
        if ( size(F%rsec4) < size(ll) ) then
          deallocate( F%rsec4 )
          allocate( F%rsec4(size(ll)*2) )
        end if
      else
        allocate( F%rsec4(size(ll)*2) )
      end if

      ! convert to 1D array:
      F%rsec4(1:size(ll)) = reshape( ll, (/ size(ll) /) )

    end if

    !  2. Number of bits per data value.                      16
    if ( present(nbits) ) F%isec4(2) = nbits

    ! ok
    status = 0

  end subroutine grib_Set



end module file_grib

