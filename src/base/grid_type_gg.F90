module grid_type_gg

  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  TggGridInfo
  
  public  ::  Init, Done
  public  ::  Check
  
  public  ::  GetLons
  public  ::  Longitudes, Latitudes
  
  public  ::  AreaOper
  
  !public  ::  Eval
  public  ::  gg_Eval_lat_lon

  public  ::  Divergence, Nabla
  

  ! --- const ----------------------------------

  character(len=*), parameter  ::  mname = 'grid_type_gg'


  ! --- types ---------------------------------
  
  ! *** location, size, etc
  
  type TggGridInfo
    ! * resolution
    integer            ::  N
    logical            ::  reduced
    ! * Gaussian latitudes
    integer            ::  nlat
    real, pointer      ::  lat(:) 
    real, pointer      ::  lat_deg(:)
    logical, pointer   ::  latflag(:)
    ! * longitudinal points for a Gaussian latitude
    integer, pointer   ::  nlon(:)
    integer            ::  nlon_reg
    real, pointer      ::  dlon(:)
    real, pointer      ::  dlon_deg(:)
    ! * total number of grid points
    integer            ::  np
    ! * start and end indices for each row
    integer, pointer   ::  i1(:)
    integer, pointer   ::  im(:)
!    ! * row number for each cell
!    integer, pointer   ::  j(:)
    ! * lat of cell boundaries for a row
    real, pointer      ::  blat(:)           ! rad
    real, pointer      ::  blat_deg(:)       ! deg
    real, pointer      ::  dlat(:)
    real, pointer      ::  dlat_deg(:)
!    ! * lon of cell and center for each cell
!    real, pointer      ::  clon(:), blon(:)
    ! * area of cell in a row
    real, pointer      ::  area(:)                    ! rad^2
    real, pointer      ::  area_m2(:)                 ! m^2
  end type TggGridInfo
  
  
  ! --- interfaces ----------------------------
  
  interface Init
    module procedure ggi_Init  
  end interface
  
  interface Done
    module procedure ggi_Done
  end interface
    
  interface Longitudes
    module procedure gg_Longitudes
  end interface
  
  interface Latitudes
    module procedure gg_Latitudes
  end interface
  
  interface Check
    module procedure gg_Check
  end interface
  
  interface AreaOper
    module procedure gg_AreaOper
  end interface
  
  !interface Eval
  !  module procedure gg_Eval_lat_lon
  !end interface
  
  interface Divergence
    module procedure Divergence_gg
  end interface
  
  interface Nabla
    module procedure Nabla_gg
  end interface
  

contains


  ! ========================================================
  
  
  subroutine ggi_Init( ggi, N, reduced, status )
  
    use Grid_Tools, only : ll_area
    use Binas, only : deg2rad, ae
    
    ! --- in/out ---------------------------------
    
    type(TggGridInfo), intent(out)     ::  ggi
    integer, intent(in)                ::  N
    logical, intent(in)                ::  reduced
    integer, intent(out)               ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ggi_Init'
    
    ! --- local ---------------------------------
    
    integer          ::  j
    
    ! --- begin ---------------------------------
    
    ! * resolution
    ggi%N = N
    ggi%reduced = reduced
    
    ! * Gaussian latitudes

    ggi%nlat = 2*ggi%N
    allocate( ggi%lat(ggi%nlat) )
    allocate( ggi%lat_deg(ggi%nlat) )
    allocate( ggi%latflag(ggi%nlat) )
    allocate( ggi%nlon(ggi%nlat) )
    allocate( ggi%dlon(ggi%nlat) )
    allocate( ggi%dlon_deg(ggi%nlat) )

    ! set Gaussian grid at northern hemisphere:
    select case ( ggi%N )

      case ( 1 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,    1 ,         1 ,         1,      30.0 )

      case ( 16 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,    1 ,        64 ,        64 ,  85.76059 )
        call SetGLat( ggi,    2 ,        64 ,        64 ,  80.26878 )
        call SetGLat( ggi,    3 ,        64 ,        64 ,  74.74454 )
        call SetGLat( ggi,    4 ,        64 ,        64 ,  69.21297 )
        call SetGLat( ggi,    5 ,        64 ,        64 ,  63.67863 )
        call SetGLat( ggi,    6 ,        64 ,        64 ,  58.14296 )
        call SetGLat( ggi,    7 ,        64 ,        64 ,  52.60653 )
        call SetGLat( ggi,    8 ,        64 ,        64 ,  47.06964 )
        call SetGLat( ggi,    9 ,        64 ,        64 ,  41.53246 )
        call SetGLat( ggi,   10 ,        64 ,        64 ,  35.99508 )
        call SetGLat( ggi,   11 ,        64 ,        64 ,  30.45755 )
        call SetGLat( ggi,   12 ,        64 ,        64 ,  24.91993 )
        call SetGLat( ggi,   13 ,        64 ,        64 ,  19.38223 )
        call SetGLat( ggi,   14 ,        64 ,        64 ,  13.84448 )
        call SetGLat( ggi,   15 ,        64 ,        64 ,  8.306703 )
        call SetGLat( ggi,   16 ,        64 ,        64 ,  2.768903 )

      case ( 24 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,    1 ,        96 ,        96 ,  87.1591   )
        call SetGLat( ggi,    2 ,        96 ,        96 ,  83.47894  )
        call SetGLat( ggi,    3 ,        96 ,        96 ,  79.77705  )
        call SetGLat( ggi,    4 ,        96 ,        96 ,  76.07024  )
        call SetGLat( ggi,    5 ,        96 ,        96 ,  72.36158  )
        call SetGLat( ggi,    6 ,        96 ,        96 ,  68.65202  )
        call SetGLat( ggi,    7 ,        96 ,        96 ,  64.94195  )
        call SetGLat( ggi,    8 ,        96 ,        96 ,  61.23158  )
        call SetGLat( ggi,    9 ,        96 ,        96 ,  57.52099  )
        call SetGLat( ggi,   10 ,        96 ,        96 ,  53.81027  )
        call SetGLat( ggi,   11 ,        96 ,        96 ,  50.09945  )
        call SetGLat( ggi,   12 ,        96 ,        96 ,  46.38856  )
        call SetGLat( ggi,   13 ,        96 ,        96 ,  42.6776   )
        call SetGLat( ggi,   14 ,        96 ,        96 ,  38.96661  )
        call SetGLat( ggi,   15 ,        96 ,        96 ,  35.25558  )
        call SetGLat( ggi,   16 ,        96 ,        96 ,  31.54452  )
        call SetGLat( ggi,   17 ,        96 ,        96 ,  27.83344  )
        call SetGLat( ggi,   18 ,        96 ,        96 ,  24.12235  )
        call SetGLat( ggi,   19 ,        96 ,        96 ,  20.41124  )
        call SetGLat( ggi,   20 ,        96 ,        96 ,  16.70012  )
        call SetGLat( ggi,   21 ,        96 ,        96 ,  12.98899  )
        call SetGLat( ggi,   22 ,        96 ,        96 ,   9.277853 )
        call SetGLat( ggi,   23 ,        96 ,        96 ,   5.566714 )
        call SetGLat( ggi,   24 ,        96 ,        96 ,   1.855572 )

      case ( 32 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,    1 ,          1,         1,     87.8638   )
        call SetGLat( ggi,    2 ,          1,         1,     85.09653  )
        call SetGLat( ggi,    3 ,          1,         1,     82.31291  )
        call SetGLat( ggi,    4 ,          1,         1,     79.5256   )
        call SetGLat( ggi,    5 ,          1,         1,     76.7369   )
        call SetGLat( ggi,    6 ,          1,         1,     73.94752  )
        call SetGLat( ggi,    7 ,          1,         1,     71.15775  )
        call SetGLat( ggi,    8 ,          1,         1,     68.36776  )
        call SetGLat( ggi,    9 ,          1,         1,     65.57761  )
        call SetGLat( ggi,   10 ,          1,         1,     62.78735  )
        call SetGLat( ggi,   11 ,          1,         1,     59.99702  )
        call SetGLat( ggi,   12 ,          1,         1,     57.20663  )
        call SetGLat( ggi,   13 ,          1,         1,     54.4162   )
        call SetGLat( ggi,   14 ,          1,         1,     51.62573  )
        call SetGLat( ggi,   15 ,          1,         1,     48.83524  )
        call SetGLat( ggi,   16 ,          1,         1,     46.04473  )
        call SetGLat( ggi,   17 ,          1,         1,     43.2542   )
        call SetGLat( ggi,   18 ,          1,         1,     40.46365  )
        call SetGLat( ggi,   19 ,          1,         1,     37.67309  )
        call SetGLat( ggi,   20 ,          1,         1,     34.88252  )
        call SetGLat( ggi,   21 ,          1,         1,     32.09195  )
        call SetGLat( ggi,   22 ,          1,         1,     29.30136  )
        call SetGLat( ggi,   23 ,          1,         1,     26.51077  )
        call SetGLat( ggi,   24 ,          1,         1,     23.72017  )
        call SetGLat( ggi,   25 ,          1,         1,     20.92957  )
        call SetGLat( ggi,   26 ,          1,         1,     18.13897  )
        call SetGLat( ggi,   27 ,          1,         1,     15.34836  )
        call SetGLat( ggi,   28 ,          1,         1,     12.55776  )
        call SetGLat( ggi,   29 ,          1,         1,      9.767145 )
        call SetGLat( ggi,   30 ,          1,         1,      6.976533 )
        call SetGLat( ggi,   31 ,          1,         1,      4.185921 )
        call SetGLat( ggi,   32 ,          1,         1,      1.395307 )

      case ( 47 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,    1 ,       192,       192,   88.5420  )
        call SetGLat( ggi,    2 ,       192,       192,   86.6531  )
        call SetGLat( ggi,    3 ,       192,       192,   84.7532  )
        call SetGLat( ggi,    4 ,       192,       192,   82.8508  )
        call SetGLat( ggi,    5 ,       192,       192,   80.9473  )
        call SetGLat( ggi,    6 ,       192,       192,   79.0435  )
        call SetGLat( ggi,    7 ,       192,       192,   77.1394  )
        call SetGLat( ggi,    8 ,       192,       192,   75.2351  )
        call SetGLat( ggi,    9 ,       192,       192,   73.3307  )
        call SetGLat( ggi,   10 ,       192,       192,   71.4262  )
        call SetGLat( ggi,   11 ,       192,       192,   69.5217  )
        call SetGLat( ggi,   12 ,       192,       192,   67.6171  )
        call SetGLat( ggi,   13 ,       192,       192,   65.7125  )
        call SetGLat( ggi,   14 ,       192,       192,   63.8079  )
        call SetGLat( ggi,   15 ,       192,       192,   61.9033  )
        call SetGLat( ggi,   16 ,       192,       192,   59.9986  )
        call SetGLat( ggi,   17 ,       192,       192,   58.0939  )
        call SetGLat( ggi,   18 ,       192,       192,   56.1893  )
        call SetGLat( ggi,   19 ,       192,       192,   54.2846  )
        call SetGLat( ggi,   20 ,       192,       192,   52.3799  )
        call SetGLat( ggi,   21 ,       192,       192,   50.4752  )
        call SetGLat( ggi,   22 ,       192,       192,   48.5705  )
        call SetGLat( ggi,   23 ,       192,       192,   46.6658  )
        call SetGLat( ggi,   24 ,       192,       192,   44.7611  )
        call SetGLat( ggi,   25 ,       192,       192,   42.8564  )
        call SetGLat( ggi,   26 ,       192,       192,   40.9517  )
        call SetGLat( ggi,   27 ,       192,       192,   39.0470  )
        call SetGLat( ggi,   28 ,       192,       192,   37.1422  )
        call SetGLat( ggi,   29 ,       192,       192,   35.2375  )
        call SetGLat( ggi,   30 ,       192,       192,   33.3328  )
        call SetGLat( ggi,   31 ,       192,       192,   31.4281  )
        call SetGLat( ggi,   32 ,       192,       192,   29.5234  )
        call SetGLat( ggi,   33 ,       192,       192,   27.6186  )
        call SetGLat( ggi,   34 ,       192,       192,   25.7139  )
        call SetGLat( ggi,   35 ,       192,       192,   23.8092  )
        call SetGLat( ggi,   36 ,       192,       192,   21.9044  )
        call SetGLat( ggi,   37 ,       192,       192,   19.9997  )
        call SetGLat( ggi,   38 ,       192,       192,   18.0950  )
        call SetGLat( ggi,   39 ,       192,       192,   16.1902  )
        call SetGLat( ggi,   40 ,       192,       192,   14.2855  )
        call SetGLat( ggi,   41 ,       192,       192,   12.3808  )
        call SetGLat( ggi,   42 ,       192,       192,   10.4760  )
        call SetGLat( ggi,   43 ,       192,       192,   8.57131  )
        call SetGLat( ggi,   44 ,       192,       192,   6.66657  )
        call SetGLat( ggi,   45 ,       192,       192,   4.76184  )
        call SetGLat( ggi,   46 ,       192,       192,   2.85710  )
        call SetGLat( ggi,   47 ,       192,       192,   0.952368 )

      case ( 80 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,    1 ,        18 ,       320,  89.14152 )
        call SetGLat( ggi,    2 ,        25 ,       320,  88.02943 )
        call SetGLat( ggi,    3 ,        36 ,       320,  86.91077 )
        call SetGLat( ggi,    4 ,        40 ,       320,  85.79063 )
        call SetGLat( ggi,    5 ,        45 ,       320,  84.66992 )
        call SetGLat( ggi,    6 ,        54 ,       320,  83.54895 )
        call SetGLat( ggi,    7 ,        60 ,       320,  82.42782 )
        call SetGLat( ggi,    8 ,        64 ,       320,  81.30659 )
        call SetGLat( ggi,    9 ,        72 ,       320,  80.18531 )
        call SetGLat( ggi,   10 ,        72 ,       320,  79.06398 )
        call SetGLat( ggi,   11 ,        80 ,       320,  77.94262 )
        call SetGLat( ggi,   12 ,        90 ,       320,  76.82124 )
        call SetGLat( ggi,   13 ,        96 ,       320,  75.69984 )
        call SetGLat( ggi,   14 ,       100 ,       320,  74.57843 )
        call SetGLat( ggi,   15 ,       108 ,       320,  73.45701 )
        call SetGLat( ggi,   16 ,       120 ,       320,  72.33558 )
        call SetGLat( ggi,   17 ,       120 ,       320,  71.21414 )
        call SetGLat( ggi,   18 ,       128 ,       320,  70.09269 )
        call SetGLat( ggi,   19 ,       135 ,       320,  68.97124 )
        call SetGLat( ggi,   20 ,       144 ,       320,  67.84978 )
        call SetGLat( ggi,   21 ,       144 ,       320,  66.72833 )
        call SetGLat( ggi,   22 ,       150 ,       320,  65.60686 )
        call SetGLat( ggi,   23 ,       160 ,       320,  64.48540 )
        call SetGLat( ggi,   24 ,       160 ,       320,  63.36393 )
        call SetGLat( ggi,   25 ,       180 ,       320,  62.24246 )
        call SetGLat( ggi,   26 ,       180 ,       320,  61.12099 )
        call SetGLat( ggi,   27 ,       180 ,       320,  59.99952 )
        call SetGLat( ggi,   28 ,       192 ,       320,  58.87804 )
        call SetGLat( ggi,   29 ,       192 ,       320,  57.75657 )
        call SetGLat( ggi,   30 ,       200 ,       320,  56.63509 )
        call SetGLat( ggi,   31 ,       200 ,       320,  55.51361 )
        call SetGLat( ggi,   32 ,       216 ,       320,  54.39214 )
        call SetGLat( ggi,   33 ,       216 ,       320,  53.27066 )
        call SetGLat( ggi,   34 ,       216 ,       320,  52.14917 )
        call SetGLat( ggi,   35 ,       225 ,       320,  51.02769 )
        call SetGLat( ggi,   36 ,       225 ,       320,  49.90621 )
        call SetGLat( ggi,   37 ,       240 ,       320,  48.78473 )
        call SetGLat( ggi,   38 ,       240 ,       320,  47.66325 )
        call SetGLat( ggi,   39 ,       240 ,       320,  46.54176 )
        call SetGLat( ggi,   40 ,       256 ,       320,  45.42028 )
        call SetGLat( ggi,   41 ,       256 ,       320,  44.29879 )
        call SetGLat( ggi,   42 ,       256 ,       320,  43.17731 )
        call SetGLat( ggi,   43 ,       256 ,       320,  42.05582 )
        call SetGLat( ggi,   44 ,       288 ,       320,  40.93434 )
        call SetGLat( ggi,   45 ,       288 ,       320,  39.81285 )
        call SetGLat( ggi,   46 ,       288 ,       320,  38.69137 )
        call SetGLat( ggi,   47 ,       288 ,       320,  37.56988 )
        call SetGLat( ggi,   48 ,       288 ,       320,  36.44839 )
        call SetGLat( ggi,   49 ,       288 ,       320,  35.32691 )
        call SetGLat( ggi,   50 ,       288 ,       320,  34.20542 )
        call SetGLat( ggi,   51 ,       288 ,       320,  33.08393 )
        call SetGLat( ggi,   52 ,       288 ,       320,  31.96244 )
        call SetGLat( ggi,   53 ,       300 ,       320,  30.84096 )
        call SetGLat( ggi,   54 ,       300 ,       320,  29.71947 )
        call SetGLat( ggi,   55 ,       300 ,       320,  28.59798 )
        call SetGLat( ggi,   56 ,       300 ,       320,  27.47649 )
        call SetGLat( ggi,   57 ,       320 ,       320,  26.35500 )
        call SetGLat( ggi,   58 ,       320 ,       320,  25.23351 )
        call SetGLat( ggi,   59 ,       320 ,       320,  24.11203 )
        call SetGLat( ggi,   60 ,       320 ,       320,  22.99054 )
        call SetGLat( ggi,   61 ,       320 ,       320,  21.86905 )
        call SetGLat( ggi,   62 ,       320 ,       320,  20.74756 )
        call SetGLat( ggi,   63 ,       320 ,       320,  19.62607 )
        call SetGLat( ggi,   64 ,       320 ,       320,  18.50458 )
        call SetGLat( ggi,   65 ,       320 ,       320,  17.38309 )
        call SetGLat( ggi,   66 ,       320 ,       320,  16.26160 )
        call SetGLat( ggi,   67 ,       320 ,       320,  15.14011 )
        call SetGLat( ggi,   68 ,       320 ,       320,  14.01862 )
        call SetGLat( ggi,   69 ,       320 ,       320,  12.89713 )
        call SetGLat( ggi,   70 ,       320 ,       320,  11.77564 )
        call SetGLat( ggi,   71 ,       320 ,       320,  10.65415 )
        call SetGLat( ggi,   72 ,       320 ,       320,   9.53266 )
        call SetGLat( ggi,   73 ,       320 ,       320,   8.41117 )
        call SetGLat( ggi,   74 ,       320 ,       320,   7.28968 )
        call SetGLat( ggi,   75 ,       320 ,       320,   6.16819 )
        call SetGLat( ggi,   76 ,       320 ,       320,   5.04670 )
        call SetGLat( ggi,   77 ,       320 ,       320,   3.92521 )
        call SetGLat( ggi,   78 ,       320 ,       320,   2.80372 )
        call SetGLat( ggi,   79 ,       320 ,       320,   1.68223 )
        call SetGLat( ggi,   80 ,       320 ,       320,   0.56074 )
        
      case ( 192 )
      
        !
        ! NCEP GFS grib files
        ! N192  
        ! Source: Doni Bundy (NCAR)
        !
 
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,     1,        768,        768,  89.46281   )
        call SetGLat( ggi,     2,        768,        768,  89.17743   )
        call SetGLat( ggi,     3,        768,        768,  88.71048   )
        call SetGLat( ggi,     4,        768,        768,  88.2429    )
        call SetGLat( ggi,     5,        768,        768,  87.77509   )
        call SetGLat( ggi,     6,        768,        768,  87.30717   )
        call SetGLat( ggi,     7,        768,        768,  86.83917   )
        call SetGLat( ggi,     8,        768,        768,  86.37115   )
        call SetGLat( ggi,     9,        768,        768,  85.9031    )
        call SetGLat( ggi,    10,        768,        768,  85.43503   )
        call SetGLat( ggi,    11,        768,        768,  84.96694   )
        call SetGLat( ggi,    12,        768,        768,  84.49885   )
        call SetGLat( ggi,    13,        768,        768,  84.03075   )
        call SetGLat( ggi,    14,        768,        768,  83.56264   )
        call SetGLat( ggi,    15,        768,        768,  83.09453   )
        call SetGLat( ggi,    16,        768,        768,  82.62641   )
        call SetGLat( ggi,    17,        768,        768,  82.15829   )
        call SetGLat( ggi,    18,        768,        768,  81.69018   )
        call SetGLat( ggi,    19,        768,        768,  81.22205   )
        call SetGLat( ggi,    20,        768,        768,  80.75393   )
        call SetGLat( ggi,    21,        768,        768,  80.2858    )
        call SetGLat( ggi,    22,        768,        768,  79.81767   )
        call SetGLat( ggi,    23,        768,        768,  79.34955   )
        call SetGLat( ggi,    24,        768,        768,  78.88142   )
        call SetGLat( ggi,    25,        768,        768,  78.41328   )
        call SetGLat( ggi,    26,        768,        768,  77.94516   )
        call SetGLat( ggi,    27,        768,        768,  77.47703   )
        call SetGLat( ggi,    28,        768,        768,  77.0089    )
        call SetGLat( ggi,    29,        768,        768,  76.54076   )
        call SetGLat( ggi,    30,        768,        768,  76.07262   )
        call SetGLat( ggi,    31,        768,        768,  75.60449   )
        call SetGLat( ggi,    32,        768,        768,  75.13636   )
        call SetGLat( ggi,    33,        768,        768,  74.66822   )
        call SetGLat( ggi,    34,        768,        768,  74.20009   )
        call SetGLat( ggi,    35,        768,        768,  73.73196   )
        call SetGLat( ggi,    36,        768,        768,  73.26382   )
        call SetGLat( ggi,    37,        768,        768,  72.79568   )
        call SetGLat( ggi,    38,        768,        768,  72.32755   )
        call SetGLat( ggi,    39,        768,        768,  71.85941   )
        call SetGLat( ggi,    40,        768,        768,  71.39127   )
        call SetGLat( ggi,    41,        768,        768,  70.92313   )
        call SetGLat( ggi,    42,        768,        768,  70.455     )
        call SetGLat( ggi,    43,        768,        768,  69.98686   )
        call SetGLat( ggi,    44,        768,        768,  69.51872   )
        call SetGLat( ggi,    45,        768,        768,  69.05059   )
        call SetGLat( ggi,    46,        768,        768,  68.58245   )
        call SetGLat( ggi,    47,        768,        768,  68.11431   )
        call SetGLat( ggi,    48,        768,        768,  67.64618   )
        call SetGLat( ggi,    49,        768,        768,  67.17804   )
        call SetGLat( ggi,    50,        768,        768,  66.7099    )
        call SetGLat( ggi,    51,        768,        768,  66.24176   )
        call SetGLat( ggi,    52,        768,        768,  65.77363   )
        call SetGLat( ggi,    53,        768,        768,  65.30549   )
        call SetGLat( ggi,    54,        768,        768,  64.83735   )
        call SetGLat( ggi,    55,        768,        768,  64.36921   )
        call SetGLat( ggi,    56,        768,        768,  63.90107   )
        call SetGLat( ggi,    57,        768,        768,  63.43293   )
        call SetGLat( ggi,    58,        768,        768,  62.96479   )
        call SetGLat( ggi,    59,        768,        768,  62.49665   )
        call SetGLat( ggi,    60,        768,        768,  62.02852   )
        call SetGLat( ggi,    61,        768,        768,  61.56038   )
        call SetGLat( ggi,    62,        768,        768,  61.09224   )
        call SetGLat( ggi,    63,        768,        768,  60.6241    )
        call SetGLat( ggi,    64,        768,        768,  60.15596   )
        call SetGLat( ggi,    65,        768,        768,  59.68782   )
        call SetGLat( ggi,    66,        768,        768,  59.21968   )
        call SetGLat( ggi,    67,        768,        768,  58.75154   )
        call SetGLat( ggi,    68,        768,        768,  58.28341   )
        call SetGLat( ggi,    69,        768,        768,  57.81527   )
        call SetGLat( ggi,    70,        768,        768,  57.34713   )
        call SetGLat( ggi,    71,        768,        768,  56.87899   )
        call SetGLat( ggi,    72,        768,        768,  56.41085   )
        call SetGLat( ggi,    73,        768,        768,  55.94271   )
        call SetGLat( ggi,    74,        768,        768,  55.47457   )
        call SetGLat( ggi,    75,        768,        768,  55.00643   )
        call SetGLat( ggi,    76,        768,        768,  54.53829   )
        call SetGLat( ggi,    77,        768,        768,  54.07016   )
        call SetGLat( ggi,    78,        768,        768,  53.60202   )
        call SetGLat( ggi,    79,        768,        768,  53.13388   )
        call SetGLat( ggi,    80,        768,        768,  52.66574   )
        call SetGLat( ggi,    81,        768,        768,  52.1976    )
        call SetGLat( ggi,    82,        768,        768,  51.72946   )
        call SetGLat( ggi,    83,        768,        768,  51.26132   )
        call SetGLat( ggi,    84,        768,        768,  50.79318   )
        call SetGLat( ggi,    85,        768,        768,  50.32504   )
        call SetGLat( ggi,    86,        768,        768,  49.8569    )
        call SetGLat( ggi,    87,        768,        768,  49.38876   )
        call SetGLat( ggi,    88,        768,        768,  48.92062   )
        call SetGLat( ggi,    89,        768,        768,  48.45248   )
        call SetGLat( ggi,    90,        768,        768,  47.98434   )
        call SetGLat( ggi,    91,        768,        768,  47.5162    )
        call SetGLat( ggi,    92,        768,        768,  47.04806   )
        call SetGLat( ggi,    93,        768,        768,  46.57992   )
        call SetGLat( ggi,    94,        768,        768,  46.11178   )
        call SetGLat( ggi,    95,        768,        768,  45.64364   )
        call SetGLat( ggi,    96,        768,        768,  45.1755    )
        call SetGLat( ggi,    97,        768,        768,  44.70736   )
        call SetGLat( ggi,    98,        768,        768,  44.23922   )
        call SetGLat( ggi,    99,        768,        768,  43.77108   )
        call SetGLat( ggi,   100,        768,        768,  43.30294   )
        call SetGLat( ggi,   101,        768,        768,  42.8348    )
        call SetGLat( ggi,   102,        768,        768,  42.36666   )
        call SetGLat( ggi,   103,        768,        768,  41.89853   )
        call SetGLat( ggi,   104,        768,        768,  41.43039   )
        call SetGLat( ggi,   105,        768,        768,  40.96225   )
        call SetGLat( ggi,   106,        768,        768,  40.49411   )
        call SetGLat( ggi,   107,        768,        768,  40.02597   )
        call SetGLat( ggi,   108,        768,        768,  39.55783   )
        call SetGLat( ggi,   109,        768,        768,  39.08969   )
        call SetGLat( ggi,   110,        768,        768,  38.62155   )
        call SetGLat( ggi,   111,        768,        768,  38.15341   )
        call SetGLat( ggi,   112,        768,        768,  37.68527   )
        call SetGLat( ggi,   113,        768,        768,  37.21713   )
        call SetGLat( ggi,   114,        768,        768,  36.74899   )
        call SetGLat( ggi,   115,        768,        768,  36.28085   )
        call SetGLat( ggi,   116,        768,        768,  35.81271   )
        call SetGLat( ggi,   117,        768,        768,  35.34457   )
        call SetGLat( ggi,   118,        768,        768,  34.87643   )
        call SetGLat( ggi,   119,        768,        768,  34.40829   )
        call SetGLat( ggi,   120,        768,        768,  33.94015   )
        call SetGLat( ggi,   121,        768,        768,  33.47201   )
        call SetGLat( ggi,   122,        768,        768,  33.00387   )
        call SetGLat( ggi,   123,        768,        768,  32.53573   )
        call SetGLat( ggi,   124,        768,        768,  32.06759   )
        call SetGLat( ggi,   125,        768,        768,  31.59945   )
        call SetGLat( ggi,   126,        768,        768,  31.13131   )
        call SetGLat( ggi,   127,        768,        768,  30.66317   )
        call SetGLat( ggi,   128,        768,        768,  30.19503   )
        call SetGLat( ggi,   129,        768,        768,  29.72689   )
        call SetGLat( ggi,   130,        768,        768,  29.25875   )
        call SetGLat( ggi,   131,        768,        768,  28.79061   )
        call SetGLat( ggi,   132,        768,        768,  28.32247   )
        call SetGLat( ggi,   133,        768,        768,  27.85433   )
        call SetGLat( ggi,   134,        768,        768,  27.38619   )
        call SetGLat( ggi,   135,        768,        768,  26.91805   )
        call SetGLat( ggi,   136,        768,        768,  26.44991   )
        call SetGLat( ggi,   137,        768,        768,  25.98177   )
        call SetGLat( ggi,   138,        768,        768,  25.51363   )
        call SetGLat( ggi,   139,        768,        768,  25.04549   )
        call SetGLat( ggi,   140,        768,        768,  24.57735   )
        call SetGLat( ggi,   141,        768,        768,  24.10921   )
        call SetGLat( ggi,   142,        768,        768,  23.64107   )
        call SetGLat( ggi,   143,        768,        768,  23.17293   )
        call SetGLat( ggi,   144,        768,        768,  22.70479   )
        call SetGLat( ggi,   145,        768,        768,  22.23665   )
        call SetGLat( ggi,   146,        768,        768,  21.76851   )
        call SetGLat( ggi,   147,        768,        768,  21.30037   )
        call SetGLat( ggi,   148,        768,        768,  20.83223   )
        call SetGLat( ggi,   149,        768,        768,  20.36409   )
        call SetGLat( ggi,   150,        768,        768,  19.89595   )
        call SetGLat( ggi,   151,        768,        768,  19.42781   )
        call SetGLat( ggi,   152,        768,        768,  18.95967   )
        call SetGLat( ggi,   153,        768,        768,  18.49153   )
        call SetGLat( ggi,   154,        768,        768,  18.02339   )
        call SetGLat( ggi,   155,        768,        768,  17.55525   )
        call SetGLat( ggi,   156,        768,        768,  17.08711   )
        call SetGLat( ggi,   157,        768,        768,  16.61897   )
        call SetGLat( ggi,   158,        768,        768,  16.15083   )
        call SetGLat( ggi,   159,        768,        768,  15.68269   )
        call SetGLat( ggi,   160,        768,        768,  15.21455   )
        call SetGLat( ggi,   161,        768,        768,  14.74641   )
        call SetGLat( ggi,   162,        768,        768,  14.27827   )
        call SetGLat( ggi,   163,        768,        768,  13.81013   )
        call SetGLat( ggi,   164,        768,        768,  13.34199   )
        call SetGLat( ggi,   165,        768,        768,  12.87385   )
        call SetGLat( ggi,   166,        768,        768,  12.40571   )
        call SetGLat( ggi,   167,        768,        768,  11.93757   )
        call SetGLat( ggi,   168,        768,        768,  11.46943   )
        call SetGLat( ggi,   169,        768,        768,  11.00129   )
        call SetGLat( ggi,   170,        768,        768,  10.53315   )
        call SetGLat( ggi,   171,        768,        768,  10.06501   )
        call SetGLat( ggi,   172,        768,        768,   9.59687   )
        call SetGLat( ggi,   173,        768,        768,   9.128731  )
        call SetGLat( ggi,   174,        768,        768,   8.660591  )
        call SetGLat( ggi,   175,        768,        768,   8.192451  )
        call SetGLat( ggi,   176,        768,        768,   7.724311  )
        call SetGLat( ggi,   177,        768,        768,   7.256171  )
        call SetGLat( ggi,   178,        768,        768,   6.788031  )
        call SetGLat( ggi,   179,        768,        768,   6.31989   )
        call SetGLat( ggi,   180,        768,        768,   5.85175   )
        call SetGLat( ggi,   181,        768,        768,   5.383611  )
        call SetGLat( ggi,   182,        768,        768,   4.915471  )
        call SetGLat( ggi,   183,        768,        768,   4.44733   )
        call SetGLat( ggi,   184,        768,        768,   3.97919   )
        call SetGLat( ggi,   185,        768,        768,   3.51105   )
        call SetGLat( ggi,   186,        768,        768,   3.04291   )
        call SetGLat( ggi,   187,        768,        768,   2.57477   )
        call SetGLat( ggi,   188,        768,        768,   2.10663   )
        call SetGLat( ggi,   189,        768,        768,   1.63849   )
        call SetGLat( ggi,   190,        768,        768,   1.17035   )
        call SetGLat( ggi,   191,        768,        768,   0.7022101 )
        call SetGLat( ggi,   192,        768,        768,   0.23407   )

      case ( 160 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,     1,         18,        640,  89.57009 )
        call SetGLat( ggi,     2,         25,        640,  89.01318 )
        call SetGLat( ggi,     3,         36,        640,  88.45297 )
        call SetGLat( ggi,     4,         40,        640,  87.89203 )
        call SetGLat( ggi,     5,         45,        640,  87.33080 )
        call SetGLat( ggi,     6,         50,        640,  86.76944 )
        call SetGLat( ggi,     7,         60,        640,  86.20800 )
        call SetGLat( ggi,     8,         64,        640,  85.64651 )
        call SetGLat( ggi,     9,         72,        640,  85.08499 )
        call SetGLat( ggi,    10,         72,        640,  84.52345 )
        call SetGLat( ggi,    11,         80,        640,  83.96190 )
        call SetGLat( ggi,    12,         90,        640,  83.40033 )
        call SetGLat( ggi,    13,         90,        640,  82.83876 )
        call SetGLat( ggi,    14,         96,        640,  82.27718 )
        call SetGLat( ggi,    15,        108,        640,  81.71559 )
        call SetGLat( ggi,    16,        120,        640,  81.15400 )
        call SetGLat( ggi,    17,        120,        640,  80.59240 )
        call SetGLat( ggi,    18,        125,        640,  80.03080 )
        call SetGLat( ggi,    19,        128,        640,  79.46920 )
        call SetGLat( ggi,    20,        135,        640,  78.90760 )
        call SetGLat( ggi,    21,        144,        640,  78.34600 )
        call SetGLat( ggi,    22,        150,        640,  77.78439 )
        call SetGLat( ggi,    23,        160,        640,  77.22278 )
        call SetGLat( ggi,    24,        160,        640,  76.66117 )
        call SetGLat( ggi,    25,        180,        640,  76.09956 )
        call SetGLat( ggi,    26,        180,        640,  75.53795 )
        call SetGLat( ggi,    27,        180,        640,  74.97634 )
        call SetGLat( ggi,    28,        192,        640,  74.41473 )
        call SetGLat( ggi,    29,        192,        640,  73.85311 )
        call SetGLat( ggi,    30,        200,        640,  73.29150 )
        call SetGLat( ggi,    31,        216,        640,  72.72988 )
        call SetGLat( ggi,    32,        216,        640,  72.16827 )
        call SetGLat( ggi,    33,        225,        640,  71.60665 )
        call SetGLat( ggi,    34,        225,        640,  71.04504 )
        call SetGLat( ggi,    35,        240,        640,  70.48342 )
        call SetGLat( ggi,    36,        240,        640,  69.92181 )
        call SetGLat( ggi,    37,        243,        640,  69.36019 )
        call SetGLat( ggi,    38,        250,        640,  68.79857 )
        call SetGLat( ggi,    39,        256,        640,  68.23695 )
        call SetGLat( ggi,    40,        270,        640,  67.67534 )
        call SetGLat( ggi,    41,        270,        640,  67.11372 )
        call SetGLat( ggi,    42,        288,        640,  66.55210 )
        call SetGLat( ggi,    43,        288,        640,  65.99048 )
        call SetGLat( ggi,    44,        288,        640,  65.42886 )
        call SetGLat( ggi,    45,        300,        640,  64.86725 )
        call SetGLat( ggi,    46,        300,        640,  64.30563 )
        call SetGLat( ggi,    47,        320,        640,  63.74401 )
        call SetGLat( ggi,    48,        320,        640,  63.18239 )
        call SetGLat( ggi,    49,        320,        640,  62.62077 )
        call SetGLat( ggi,    50,        320,        640,  62.05915 )
        call SetGLat( ggi,    51,        324,        640,  61.49753 )
        call SetGLat( ggi,    52,        360,        640,  60.93591 )
        call SetGLat( ggi,    53,        360,        640,  60.37429 )
        call SetGLat( ggi,    54,        360,        640,  59.81267 )
        call SetGLat( ggi,    55,        360,        640,  59.25105 )
        call SetGLat( ggi,    56,        360,        640,  58.68943 )
        call SetGLat( ggi,    57,        360,        640,  58.12781 )
        call SetGLat( ggi,    58,        375,        640,  57.56619 )
        call SetGLat( ggi,    59,        375,        640,  57.00457 )
        call SetGLat( ggi,    60,        375,        640,  56.44295 )
        call SetGLat( ggi,    61,        384,        640,  55.88133 )
        call SetGLat( ggi,    62,        384,        640,  55.31971 )
        call SetGLat( ggi,    63,        400,        640,  54.75809 )
        call SetGLat( ggi,    64,        400,        640,  54.19647 )
        call SetGLat( ggi,    65,        400,        640,  53.63485 )
        call SetGLat( ggi,    66,        405,        640,  53.07323 )
        call SetGLat( ggi,    67,        432,        640,  52.51161 )
        call SetGLat( ggi,    68,        432,        640,  51.94999 )
        call SetGLat( ggi,    69,        432,        640,  51.38837 )
        call SetGLat( ggi,    70,        432,        640,  50.82675 )
        call SetGLat( ggi,    71,        432,        640,  50.26513 )
        call SetGLat( ggi,    72,        450,        640,  49.70351 )
        call SetGLat( ggi,    73,        450,        640,  49.14189 )
        call SetGLat( ggi,    74,        450,        640,  48.58026 )
        call SetGLat( ggi,    75,        450,        640,  48.01864 )
        call SetGLat( ggi,    76,        480,        640,  47.45702 )
        call SetGLat( ggi,    77,        480,        640,  46.89540 )
        call SetGLat( ggi,    78,        480,        640,  46.33378 )
        call SetGLat( ggi,    79,        480,        640,  45.77216 )
        call SetGLat( ggi,    80,        480,        640,  45.21054 )
        call SetGLat( ggi,    81,        480,        640,  44.64892 )
        call SetGLat( ggi,    82,        480,        640,  44.08730 )
        call SetGLat( ggi,    83,        500,        640,  43.52567 )
        call SetGLat( ggi,    84,        500,        640,  42.96405 )
        call SetGLat( ggi,    85,        500,        640,  42.40243 )
        call SetGLat( ggi,    86,        500,        640,  41.84081 )
        call SetGLat( ggi,    87,        500,        640,  41.27919 )
        call SetGLat( ggi,    88,        512,        640,  40.71757 )
        call SetGLat( ggi,    89,        512,        640,  40.15595 )
        call SetGLat( ggi,    90,        540,        640,  39.59433 )
        call SetGLat( ggi,    91,        540,        640,  39.03270 )
        call SetGLat( ggi,    92,        540,        640,  38.47108 )
        call SetGLat( ggi,    93,        540,        640,  37.90946 )
        call SetGLat( ggi,    94,        540,        640,  37.34784 )
        call SetGLat( ggi,    95,        540,        640,  36.78622 )
        call SetGLat( ggi,    96,        540,        640,  36.22460 )
        call SetGLat( ggi,    97,        540,        640,  35.66298 )
        call SetGLat( ggi,    98,        576,        640,  35.10136 )
        call SetGLat( ggi,    99,        576,        640,  34.53973 )
        call SetGLat( ggi,   100,        576,        640,  33.97811 )
        call SetGLat( ggi,   101,        576,        640,  33.41649 )
        call SetGLat( ggi,   102,        576,        640,  32.85487 )
        call SetGLat( ggi,   103,        576,        640,  32.29325 )
        call SetGLat( ggi,   104,        576,        640,  31.73163 )
        call SetGLat( ggi,   105,        576,        640,  31.17000 )
        call SetGLat( ggi,   106,        576,        640,  30.60838 )
        call SetGLat( ggi,   107,        576,        640,  30.04676 )
        call SetGLat( ggi,   108,        600,        640,  29.48514 )
        call SetGLat( ggi,   109,        600,        640,  28.92352 )
        call SetGLat( ggi,   110,        600,        640,  28.36190 )
        call SetGLat( ggi,   111,        600,        640,  27.80028 )
        call SetGLat( ggi,   112,        600,        640,  27.23865 )
        call SetGLat( ggi,   113,        600,        640,  26.67703 )
        call SetGLat( ggi,   114,        600,        640,  26.11541 )
        call SetGLat( ggi,   115,        600,        640,  25.55379 )
        call SetGLat( ggi,   116,        600,        640,  24.99217 )
        call SetGLat( ggi,   117,        640,        640,  24.43055 )
        call SetGLat( ggi,   118,        640,        640,  23.86892 )
        call SetGLat( ggi,   119,        640,        640,  23.30730 )
        call SetGLat( ggi,   120,        640,        640,  22.74568 )
        call SetGLat( ggi,   121,        640,        640,  22.18406 )
        call SetGLat( ggi,   122,        640,        640,  21.62244 )
        call SetGLat( ggi,   123,        640,        640,  21.06082 )
        call SetGLat( ggi,   124,        640,        640,  20.49919 )
        call SetGLat( ggi,   125,        640,        640,  19.93757 )
        call SetGLat( ggi,   126,        640,        640,  19.37595 )
        call SetGLat( ggi,   127,        640,        640,  18.81433 )
        call SetGLat( ggi,   128,        640,        640,  18.25271 )
        call SetGLat( ggi,   129,        640,        640,  17.69109 )
        call SetGLat( ggi,   130,        640,        640,  17.12946 )
        call SetGLat( ggi,   131,        640,        640,  16.56784 )
        call SetGLat( ggi,   132,        640,        640,  16.00622 )
        call SetGLat( ggi,   133,        640,        640,  15.44460 )
        call SetGLat( ggi,   134,        640,        640,  14.88298 )
        call SetGLat( ggi,   135,        640,        640,  14.32136 )
        call SetGLat( ggi,   136,        640,        640,  13.75973 )
        call SetGLat( ggi,   137,        640,        640,  13.19811 )
        call SetGLat( ggi,   138,        640,        640,  12.63649 )
        call SetGLat( ggi,   139,        640,        640,  12.07487 )
        call SetGLat( ggi,   140,        640,        640,  11.51325 )
        call SetGLat( ggi,   141,        640,        640,  10.95162 )
        call SetGLat( ggi,   142,        640,        640,  10.39000 )
        call SetGLat( ggi,   143,        640,        640,   9.82838 )
        call SetGLat( ggi,   144,        640,        640,   9.26676 )
        call SetGLat( ggi,   145,        640,        640,   8.70514 )
        call SetGLat( ggi,   146,        640,        640,   8.14352 )
        call SetGLat( ggi,   147,        640,        640,   7.58189 )
        call SetGLat( ggi,   148,        640,        640,   7.02027 )
        call SetGLat( ggi,   149,        640,        640,   6.45865 )
        call SetGLat( ggi,   150,        640,        640,   5.89703 )
        call SetGLat( ggi,   151,        640,        640,   5.33541 )
        call SetGLat( ggi,   152,        640,        640,   4.77379 )
        call SetGLat( ggi,   153,        640,        640,   4.21216 )
        call SetGLat( ggi,   154,        640,        640,   3.65054 )
        call SetGLat( ggi,   155,        640,        640,   3.08892 )
        call SetGLat( ggi,   156,        640,        640,   2.52730 )
        call SetGLat( ggi,   157,        640,        640,   1.96568 )
        call SetGLat( ggi,   158,        640,        640,   1.40405 )
        call SetGLat( ggi,   159,        640,        640,   0.84243 )
        call SetGLat( ggi,   160,        640,        640,   0.28081 )

      case ( 256 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,     1,         18,       1024,   89.73115 )
        call SetGLat( ggi,     2,         25,       1024,   89.38287 )
        call SetGLat( ggi,     3,         32,       1024,   89.03254 )
        call SetGLat( ggi,     4,         40,       1024,   88.68175 )
        call SetGLat( ggi,     5,         45,       1024,   88.33077 )
        call SetGLat( ggi,     6,         50,       1024,   87.97972 )
        call SetGLat( ggi,     7,         60,       1024,   87.62861 )
        call SetGLat( ggi,     8,         64,       1024,   87.27748 )
        call SetGLat( ggi,     9,         72,       1024,   86.92632 )
        call SetGLat( ggi,    10,         72,       1024,   86.57515 )
        call SetGLat( ggi,    11,         75,       1024,   86.22398 )
        call SetGLat( ggi,    12,         81,       1024,   85.87279 )
        call SetGLat( ggi,    13,         90,       1024,   85.52160 )
        call SetGLat( ggi,    14,         96,       1024,   85.17041 )
        call SetGLat( ggi,    15,        100,       1024,   84.81921 )
        call SetGLat( ggi,    16,        108,       1024,   84.46801 )
        call SetGLat( ggi,    17,        120,       1024,   84.11681 )
        call SetGLat( ggi,    18,        120,       1024,   83.76560 )
        call SetGLat( ggi,    19,        125,       1024,   83.41440 )
        call SetGLat( ggi,    20,        135,       1024,   83.06319 )
        call SetGLat( ggi,    21,        144,       1024,   82.71198 )
        call SetGLat( ggi,    22,        150,       1024,   82.36077 )
        call SetGLat( ggi,    23,        160,       1024,   82.00956 )
        call SetGLat( ggi,    24,        160,       1024,   81.65835 )
        call SetGLat( ggi,    25,        180,       1024,   81.30714 )
        call SetGLat( ggi,    26,        180,       1024,   80.95593 )
        call SetGLat( ggi,    27,        180,       1024,   80.60471 )
        call SetGLat( ggi,    28,        192,       1024,   80.25350 )
        call SetGLat( ggi,    29,        192,       1024,   79.90229 )
        call SetGLat( ggi,    30,        200,       1024,   79.55107 )
        call SetGLat( ggi,    31,        216,       1024,   79.19986 )
        call SetGLat( ggi,    32,        216,       1024,   78.84864 )
        call SetGLat( ggi,    33,        216,       1024,   78.49743 )
        call SetGLat( ggi,    34,        225,       1024,   78.14621 )
        call SetGLat( ggi,    35,        240,       1024,   77.79500 )
        call SetGLat( ggi,    36,        240,       1024,   77.44378 )
        call SetGLat( ggi,    37,        243,       1024,   77.09256 )
        call SetGLat( ggi,    38,        250,       1024,   76.74135 )
        call SetGLat( ggi,    39,        256,       1024,   76.39013 )
        call SetGLat( ggi,    40,        270,       1024,   76.03891 )
        call SetGLat( ggi,    41,        270,       1024,   75.68770 )
        call SetGLat( ggi,    42,        288,       1024,   75.33648 )
        call SetGLat( ggi,    43,        288,       1024,   74.98526 )
        call SetGLat( ggi,    44,        288,       1024,   74.63405 )
        call SetGLat( ggi,    45,        300,       1024,   74.28283 )
        call SetGLat( ggi,    46,        300,       1024,   73.93161 )
        call SetGLat( ggi,    47,        320,       1024,   73.58040 )
        call SetGLat( ggi,    48,        320,       1024,   73.22918 )
        call SetGLat( ggi,    49,        320,       1024,   72.87796 )
        call SetGLat( ggi,    50,        324,       1024,   72.52674 )
        call SetGLat( ggi,    51,        360,       1024,   72.17552 )
        call SetGLat( ggi,    52,        360,       1024,   71.82431 )
        call SetGLat( ggi,    53,        360,       1024,   71.47309 )
        call SetGLat( ggi,    54,        360,       1024,   71.12187 )
        call SetGLat( ggi,    55,        360,       1024,   70.77065 )
        call SetGLat( ggi,    56,        360,       1024,   70.41944 )
        call SetGLat( ggi,    57,        375,       1024,   70.06822 )
        call SetGLat( ggi,    58,        375,       1024,   69.71700 )
        call SetGLat( ggi,    59,        384,       1024,   69.36578 )
        call SetGLat( ggi,    60,        384,       1024,   69.01456 )
        call SetGLat( ggi,    61,        400,       1024,   68.66334 )
        call SetGLat( ggi,    62,        400,       1024,   68.31213 )
        call SetGLat( ggi,    63,        400,       1024,   67.96091 )
        call SetGLat( ggi,    64,        432,       1024,   67.60969 )
        call SetGLat( ggi,    65,        432,       1024,   67.25847 )
        call SetGLat( ggi,    66,        432,       1024,   66.90725 )
        call SetGLat( ggi,    67,        432,       1024,   66.55603 )
        call SetGLat( ggi,    68,        432,       1024,   66.20482 )
        call SetGLat( ggi,    69,        450,       1024,   65.85360 )
        call SetGLat( ggi,    70,        450,       1024,   65.50238 )
        call SetGLat( ggi,    71,        450,       1024,   65.15116 )
        call SetGLat( ggi,    72,        480,       1024,   64.79994 )
        call SetGLat( ggi,    73,        480,       1024,   64.44872 )
        call SetGLat( ggi,    74,        480,       1024,   64.09750 )
        call SetGLat( ggi,    75,        480,       1024,   63.74629 )
        call SetGLat( ggi,    76,        480,       1024,   63.39507 )
        call SetGLat( ggi,    77,        486,       1024,   63.04385 )
        call SetGLat( ggi,    78,        500,       1024,   62.69263 )
        call SetGLat( ggi,    79,        500,       1024,   62.34141 )
        call SetGLat( ggi,    80,        500,       1024,   61.99019 )
        call SetGLat( ggi,    81,        512,       1024,   61.63897 )
        call SetGLat( ggi,    82,        512,       1024,   61.28776 )
        call SetGLat( ggi,    83,        540,       1024,   60.93654 )
        call SetGLat( ggi,    84,        540,       1024,   60.58532 )
        call SetGLat( ggi,    85,        540,       1024,   60.23410 )
        call SetGLat( ggi,    86,        540,       1024,   59.88288 )
        call SetGLat( ggi,    87,        540,       1024,   59.53166 )
        call SetGLat( ggi,    88,        576,       1024,   59.18044 )
        call SetGLat( ggi,    89,        576,       1024,   58.82922 )
        call SetGLat( ggi,    90,        576,       1024,   58.47800 )
        call SetGLat( ggi,    91,        576,       1024,   58.12679 )
        call SetGLat( ggi,    92,        576,       1024,   57.77557 )
        call SetGLat( ggi,    93,        576,       1024,   57.42435 )
        call SetGLat( ggi,    94,        600,       1024,   57.07313 )
        call SetGLat( ggi,    95,        600,       1024,   56.72191 )
        call SetGLat( ggi,    96,        600,       1024,   56.37069 )
        call SetGLat( ggi,    97,        600,       1024,   56.01947 )
        call SetGLat( ggi,    98,        600,       1024,   55.66825 )
        call SetGLat( ggi,    99,        640,       1024,   55.31703 )
        call SetGLat( ggi,   100,        640,       1024,   54.96581 )
        call SetGLat( ggi,   101,        640,       1024,   54.61460 )
        call SetGLat( ggi,   102,        640,       1024,   54.26338 )
        call SetGLat( ggi,   103,        640,       1024,   53.91216 )
        call SetGLat( ggi,   104,        640,       1024,   53.56094 )
        call SetGLat( ggi,   105,        640,       1024,   53.20972 )
        call SetGLat( ggi,   106,        640,       1024,   52.85850 )
        call SetGLat( ggi,   107,        648,       1024,   52.50728 )
        call SetGLat( ggi,   108,        675,       1024,   52.15606 )
        call SetGLat( ggi,   109,        675,       1024,   51.80484 )
        call SetGLat( ggi,   110,        675,       1024,   51.45362 )
        call SetGLat( ggi,   111,        675,       1024,   51.10241 )
        call SetGLat( ggi,   112,        675,       1024,   50.75119 )
        call SetGLat( ggi,   113,        675,       1024,   50.39997 )
        call SetGLat( ggi,   114,        720,       1024,   50.04875 )
        call SetGLat( ggi,   115,        720,       1024,   49.69753 )
        call SetGLat( ggi,   116,        720,       1024,   49.34631 )
        call SetGLat( ggi,   117,        720,       1024,   48.99509 )
        call SetGLat( ggi,   118,        720,       1024,   48.64387 )
        call SetGLat( ggi,   119,        720,       1024,   48.29265 )
        call SetGLat( ggi,   120,        720,       1024,   47.94143 )
        call SetGLat( ggi,   121,        720,       1024,   47.59021 )
        call SetGLat( ggi,   122,        720,       1024,   47.23899 )
        call SetGLat( ggi,   123,        729,       1024,   46.88778 )
        call SetGLat( ggi,   124,        729,       1024,   46.53656 )
        call SetGLat( ggi,   125,        750,       1024,   46.18534 )
        call SetGLat( ggi,   126,        750,       1024,   45.83412 )
        call SetGLat( ggi,   127,        750,       1024,   45.48290 )
        call SetGLat( ggi,   128,        750,       1024,   45.13168 )
        call SetGLat( ggi,   129,        750,       1024,   44.78046 )
        call SetGLat( ggi,   130,        768,       1024,   44.42924 )
        call SetGLat( ggi,   131,        768,       1024,   44.07802 )
        call SetGLat( ggi,   132,        768,       1024,   43.72680 )
        call SetGLat( ggi,   133,        768,       1024,   43.37558 )
        call SetGLat( ggi,   134,        800,       1024,   43.02436 )
        call SetGLat( ggi,   135,        800,       1024,   42.67315 )
        call SetGLat( ggi,   136,        800,       1024,   42.32193 )
        call SetGLat( ggi,   137,        800,       1024,   41.97071 )
        call SetGLat( ggi,   138,        800,       1024,   41.61949 )
        call SetGLat( ggi,   139,        800,       1024,   41.26827 )
        call SetGLat( ggi,   140,        800,       1024,   40.91705 )
        call SetGLat( ggi,   141,        800,       1024,   40.56583 )
        call SetGLat( ggi,   142,        810,       1024,   40.21461 )
        call SetGLat( ggi,   143,        810,       1024,   39.86339 )
        call SetGLat( ggi,   144,        864,       1024,   39.51217 )
        call SetGLat( ggi,   145,        864,       1024,   39.16095 )
        call SetGLat( ggi,   146,        864,       1024,   38.80973 )
        call SetGLat( ggi,   147,        864,       1024,   38.45851 )
        call SetGLat( ggi,   148,        864,       1024,   38.10730 )
        call SetGLat( ggi,   149,        864,       1024,   37.75608 )
        call SetGLat( ggi,   150,        864,       1024,   37.40486 )
        call SetGLat( ggi,   151,        864,       1024,   37.05364 )
        call SetGLat( ggi,   152,        864,       1024,   36.70242 )
        call SetGLat( ggi,   153,        864,       1024,   36.35120 )
        call SetGLat( ggi,   154,        864,       1024,   35.99998 )
        call SetGLat( ggi,   155,        864,       1024,   35.64876 )
        call SetGLat( ggi,   156,        864,       1024,   35.29754 )
        call SetGLat( ggi,   157,        864,       1024,   34.94632 )
        call SetGLat( ggi,   158,        900,       1024,   34.59510 )
        call SetGLat( ggi,   159,        900,       1024,   34.24388 )
        call SetGLat( ggi,   160,        900,       1024,   33.89266 )
        call SetGLat( ggi,   161,        900,       1024,   33.54145 )
        call SetGLat( ggi,   162,        900,       1024,   33.19023 )
        call SetGLat( ggi,   163,        900,       1024,   32.83901 )
        call SetGLat( ggi,   164,        900,       1024,   32.48779 )
        call SetGLat( ggi,   165,        900,       1024,   32.13657 )
        call SetGLat( ggi,   166,        900,       1024,   31.78535 )
        call SetGLat( ggi,   167,        900,       1024,   31.43413 )
        call SetGLat( ggi,   168,        900,       1024,   31.08291 )
        call SetGLat( ggi,   169,        960,       1024,   30.73169 )
        call SetGLat( ggi,   170,        960,       1024,   30.38047 )
        call SetGLat( ggi,   171,        960,       1024,   30.02925 )
        call SetGLat( ggi,   172,        960,       1024,   29.67803 )
        call SetGLat( ggi,   173,        960,       1024,   29.32681 )
        call SetGLat( ggi,   174,        960,       1024,   28.97559 )
        call SetGLat( ggi,   175,        960,       1024,   28.62438 )
        call SetGLat( ggi,   176,        960,       1024,   28.27316 )
        call SetGLat( ggi,   177,        960,       1024,   27.92194 )
        call SetGLat( ggi,   178,        960,       1024,   27.57072 )
        call SetGLat( ggi,   179,        960,       1024,   27.21950 )
        call SetGLat( ggi,   180,        960,       1024,   26.86828 )
        call SetGLat( ggi,   181,        960,       1024,   26.51706 )
        call SetGLat( ggi,   182,        960,       1024,   26.16584 )
        call SetGLat( ggi,   183,        960,       1024,   25.81462 )
        call SetGLat( ggi,   184,        960,       1024,   25.46340 )
        call SetGLat( ggi,   185,        960,       1024,   25.11218 )
        call SetGLat( ggi,   186,        960,       1024,   24.76096 )
        call SetGLat( ggi,   187,        960,       1024,   24.40974 )
        call SetGLat( ggi,   188,        960,       1024,   24.05852 )
        call SetGLat( ggi,   189,        960,       1024,   23.70731 )
        call SetGLat( ggi,   190,        960,       1024,   23.35609 )
        call SetGLat( ggi,   191,        972,       1024,   23.00487 )
        call SetGLat( ggi,   192,        972,       1024,   22.65365 )
        call SetGLat( ggi,   193,        972,       1024,   22.30243 )
        call SetGLat( ggi,   194,        972,       1024,   21.95121 )
        call SetGLat( ggi,   195,        972,       1024,   21.59999 )
        call SetGLat( ggi,   196,       1000,       1024,   21.24877 )
        call SetGLat( ggi,   197,       1000,       1024,   20.89755 )
        call SetGLat( ggi,   198,       1000,       1024,   20.54633 )
        call SetGLat( ggi,   199,       1000,       1024,   20.19511 )
        call SetGLat( ggi,   200,       1000,       1024,   19.84389 )
        call SetGLat( ggi,   201,       1000,       1024,   19.49267 )
        call SetGLat( ggi,   202,       1000,       1024,   19.14145 )
        call SetGLat( ggi,   203,       1000,       1024,   18.79023 )
        call SetGLat( ggi,   204,       1000,       1024,   18.43902 )
        call SetGLat( ggi,   205,       1000,       1024,   18.08780 )
        call SetGLat( ggi,   206,       1000,       1024,   17.73658 )
        call SetGLat( ggi,   207,       1000,       1024,   17.38536 )
        call SetGLat( ggi,   208,       1000,       1024,   17.03414 )
        call SetGLat( ggi,   209,       1000,       1024,   16.68292 )
        call SetGLat( ggi,   210,       1000,       1024,   16.33170 )
        call SetGLat( ggi,   211,       1000,       1024,   15.98048 )
        call SetGLat( ggi,   212,       1024,       1024,   15.62926 )
        call SetGLat( ggi,   213,       1024,       1024,   15.27804 )
        call SetGLat( ggi,   214,       1024,       1024,   14.92682 )
        call SetGLat( ggi,   215,       1024,       1024,   14.57560 )
        call SetGLat( ggi,   216,       1024,       1024,   14.22438 )
        call SetGLat( ggi,   217,       1024,       1024,   13.87316 )
        call SetGLat( ggi,   218,       1024,       1024,   13.52194 )
        call SetGLat( ggi,   219,       1024,       1024,   13.17073 )
        call SetGLat( ggi,   220,       1024,       1024,   12.81951 )
        call SetGLat( ggi,   221,       1024,       1024,   12.46829 )
        call SetGLat( ggi,   222,       1024,       1024,   12.11707 )
        call SetGLat( ggi,   223,       1024,       1024,   11.76585 )
        call SetGLat( ggi,   224,       1024,       1024,   11.41463 )
        call SetGLat( ggi,   225,       1024,       1024,   11.06341 )
        call SetGLat( ggi,   226,       1024,       1024,   10.71219 )
        call SetGLat( ggi,   227,       1024,       1024,   10.36097 )
        call SetGLat( ggi,   228,       1024,       1024,   10.00975 )
        call SetGLat( ggi,   229,       1024,       1024,    9.65853 )
        call SetGLat( ggi,   230,       1024,       1024,    9.30731 )
        call SetGLat( ggi,   231,       1024,       1024,    8.95609 )
        call SetGLat( ggi,   232,       1024,       1024,    8.60487 )
        call SetGLat( ggi,   233,       1024,       1024,    8.25365 )
        call SetGLat( ggi,   234,       1024,       1024,    7.90244 )
        call SetGLat( ggi,   235,       1024,       1024,    7.55122 )
        call SetGLat( ggi,   236,       1024,       1024,    7.20000 )
        call SetGLat( ggi,   237,       1024,       1024,    6.84878 )
        call SetGLat( ggi,   238,       1024,       1024,    6.49756 )
        call SetGLat( ggi,   239,       1024,       1024,    6.14634 )
        call SetGLat( ggi,   240,       1024,       1024,    5.79512 )
        call SetGLat( ggi,   241,       1024,       1024,    5.44390 )
        call SetGLat( ggi,   242,       1024,       1024,    5.09268 )
        call SetGLat( ggi,   243,       1024,       1024,    4.74146 )
        call SetGLat( ggi,   244,       1024,       1024,    4.39024 )
        call SetGLat( ggi,   245,       1024,       1024,    4.03902 )
        call SetGLat( ggi,   246,       1024,       1024,    3.68780 )
        call SetGLat( ggi,   247,       1024,       1024,    3.33658 )
        call SetGLat( ggi,   248,       1024,       1024,    2.98536 )
        call SetGLat( ggi,   249,       1024,       1024,    2.63415 )
        call SetGLat( ggi,   250,       1024,       1024,    2.28293 )
        call SetGLat( ggi,   251,       1024,       1024,    1.93171 )
        call SetGLat( ggi,   252,       1024,       1024,    1.58049 )
        call SetGLat( ggi,   253,       1024,       1024,    1.22927 )
        call SetGLat( ggi,   254,       1024,       1024,    0.87805 )
        call SetGLat( ggi,   255,       1024,       1024,    0.52683 )
        call SetGLat( ggi,   256,       1024,       1024,    0.17561 )

      case ( 320 )
        !                  latitude   reduced     regular  latitude
        !                   number     points      points
        !                  -------    -------     -------  --------
        call SetGLat( ggi,     1,         18,       1280,  89.78488 )
        call SetGLat( ggi,     2,         25,       1280,  89.50620 )
        call SetGLat( ggi,     3,         32,       1280,  89.22588 )
        call SetGLat( ggi,     4,         40,       1280,  88.94519 )
        call SetGLat( ggi,     5,         45,       1280,  88.66436 )
        call SetGLat( ggi,     6,         50,       1280,  88.38346 )
        call SetGLat( ggi,     7,         60,       1280,  88.10252 )
        call SetGLat( ggi,     8,         64,       1280,  87.82156 )
        call SetGLat( ggi,     9,         72,       1280,  87.54058 )
        call SetGLat( ggi,    10,         72,       1280,  87.25959 )
        call SetGLat( ggi,    11,         75,       1280,  86.97859 )
        call SetGLat( ggi,    12,         81,       1280,  86.69759 )
        call SetGLat( ggi,    13,         90,       1280,  86.41658 )
        call SetGLat( ggi,    14,         96,       1280,  86.13557 )
        call SetGLat( ggi,    15,        100,       1280,  85.85456 )
        call SetGLat( ggi,    16,        108,       1280,  85.57355 )
        call SetGLat( ggi,    17,        120,       1280,  85.29253 )
        call SetGLat( ggi,    18,        120,       1280,  85.01151 )
        call SetGLat( ggi,    19,        125,       1280,  84.73049 )
        call SetGLat( ggi,    20,        135,       1280,  84.44947 )
        call SetGLat( ggi,    21,        144,       1280,  84.16845 )
        call SetGLat( ggi,    22,        144,       1280,  83.88742 )
        call SetGLat( ggi,    23,        150,       1280,  83.60640 )
        call SetGLat( ggi,    24,        160,       1280,  83.32538 )
        call SetGLat( ggi,    25,        180,       1280,  83.04435 )
        call SetGLat( ggi,    26,        180,       1280,  82.76333 )
        call SetGLat( ggi,    27,        180,       1280,  82.48230 )
        call SetGLat( ggi,    28,        192,       1280,  82.20128 )
        call SetGLat( ggi,    29,        192,       1280,  81.92025 )
        call SetGLat( ggi,    30,        200,       1280,  81.63923 )
        call SetGLat( ggi,    31,        216,       1280,  81.35820 )
        call SetGLat( ggi,    32,        216,       1280,  81.07717 )
        call SetGLat( ggi,    33,        216,       1280,  80.79615 )
        call SetGLat( ggi,    34,        225,       1280,  80.51512 )
        call SetGLat( ggi,    35,        240,       1280,  80.23409 )
        call SetGLat( ggi,    36,        240,       1280,  79.95306 )
        call SetGLat( ggi,    37,        240,       1280,  79.67204 )
        call SetGLat( ggi,    38,        250,       1280,  79.39101 )
        call SetGLat( ggi,    39,        256,       1280,  79.10998 )
        call SetGLat( ggi,    40,        270,       1280,  78.82895 )
        call SetGLat( ggi,    41,        270,       1280,  78.54792 )
        call SetGLat( ggi,    42,        288,       1280,  78.26689 )
        call SetGLat( ggi,    43,        288,       1280,  77.98587 )
        call SetGLat( ggi,    44,        288,       1280,  77.70484 )
        call SetGLat( ggi,    45,        300,       1280,  77.42381 )
        call SetGLat( ggi,    46,        300,       1280,  77.14278 )
        call SetGLat( ggi,    47,        320,       1280,  76.86175 )
        call SetGLat( ggi,    48,        320,       1280,  76.58072 )
        call SetGLat( ggi,    49,        320,       1280,  76.29969 )
        call SetGLat( ggi,    50,        324,       1280,  76.01867 )
        call SetGLat( ggi,    51,        360,       1280,  75.73764 )
        call SetGLat( ggi,    52,        360,       1280,  75.45661 )
        call SetGLat( ggi,    53,        360,       1280,  75.17558 )
        call SetGLat( ggi,    54,        360,       1280,  74.89455 )
        call SetGLat( ggi,    55,        360,       1280,  74.61352 )
        call SetGLat( ggi,    56,        360,       1280,  74.33249 )
        call SetGLat( ggi,    57,        375,       1280,  74.05146 )
        call SetGLat( ggi,    58,        375,       1280,  73.77043 )
        call SetGLat( ggi,    59,        384,       1280,  73.48940 )
        call SetGLat( ggi,    60,        384,       1280,  73.20837 )
        call SetGLat( ggi,    61,        400,       1280,  72.92734 )
        call SetGLat( ggi,    62,        400,       1280,  72.64631 )
        call SetGLat( ggi,    63,        405,       1280,  72.36528 )
        call SetGLat( ggi,    64,        432,       1280,  72.08426 )
        call SetGLat( ggi,    65,        432,       1280,  71.80323 )
        call SetGLat( ggi,    66,        432,       1280,  71.52220 )
        call SetGLat( ggi,    67,        432,       1280,  71.24117 )
        call SetGLat( ggi,    68,        450,       1280,  70.96014 )
        call SetGLat( ggi,    69,        450,       1280,  70.67911 )
        call SetGLat( ggi,    70,        450,       1280,  70.39808 )
        call SetGLat( ggi,    71,        480,       1280,  70.11705 )
        call SetGLat( ggi,    72,        480,       1280,  69.83602 )
        call SetGLat( ggi,    73,        480,       1280,  69.55499 )
        call SetGLat( ggi,    74,        480,       1280,  69.27396 )
        call SetGLat( ggi,    75,        480,       1280,  68.99293 )
        call SetGLat( ggi,    76,        486,       1280,  68.71190 )
        call SetGLat( ggi,    77,        500,       1280,  68.43087 )
        call SetGLat( ggi,    78,        500,       1280,  68.14984 )
        call SetGLat( ggi,    79,        500,       1280,  67.86881 )
        call SetGLat( ggi,    80,        512,       1280,  67.58778 )
        call SetGLat( ggi,    81,        512,       1280,  67.30675 )
        call SetGLat( ggi,    82,        540,       1280,  67.02572 )
        call SetGLat( ggi,    83,        540,       1280,  66.74469 )
        call SetGLat( ggi,    84,        540,       1280,  66.46366 )
        call SetGLat( ggi,    85,        540,       1280,  66.18263 )
        call SetGLat( ggi,    86,        540,       1280,  65.90160 )
        call SetGLat( ggi,    87,        576,       1280,  65.62057 )
        call SetGLat( ggi,    88,        576,       1280,  65.33954 )
        call SetGLat( ggi,    89,        576,       1280,  65.05851 )
        call SetGLat( ggi,    90,        576,       1280,  64.77748 )
        call SetGLat( ggi,    91,        576,       1280,  64.49645 )
        call SetGLat( ggi,    92,        576,       1280,  64.21542 )
        call SetGLat( ggi,    93,        600,       1280,  63.93439 )
        call SetGLat( ggi,    94,        600,       1280,  63.65336 )
        call SetGLat( ggi,    95,        600,       1280,  63.37233 )
        call SetGLat( ggi,    96,        600,       1280,  63.09130 )
        call SetGLat( ggi,    97,        625,       1280,  62.81027 )
        call SetGLat( ggi,    98,        625,       1280,  62.52924 )
        call SetGLat( ggi,    99,        625,       1280,  62.24821 )
        call SetGLat( ggi,   100,        625,       1280,  61.96718 )
        call SetGLat( ggi,   101,        625,       1280,  61.68615 )
        call SetGLat( ggi,   102,        640,       1280,  61.40512 )
        call SetGLat( ggi,   103,        640,       1280,  61.12409 )
        call SetGLat( ggi,   104,        648,       1280,  60.84306 )
        call SetGLat( ggi,   105,        648,       1280,  60.56203 )
        call SetGLat( ggi,   106,        675,       1280,  60.28100 )
        call SetGLat( ggi,   107,        675,       1280,  59.99997 )
        call SetGLat( ggi,   108,        675,       1280,  59.71894 )
        call SetGLat( ggi,   109,        675,       1280,  59.43791 )
        call SetGLat( ggi,   110,        720,       1280,  59.15688 )
        call SetGLat( ggi,   111,        720,       1280,  58.87585 )
        call SetGLat( ggi,   112,        720,       1280,  58.59482 )
        call SetGLat( ggi,   113,        720,       1280,  58.31379 )
        call SetGLat( ggi,   114,        720,       1280,  58.03276 )
        call SetGLat( ggi,   115,        720,       1280,  57.75173 )
        call SetGLat( ggi,   116,        720,       1280,  57.47070 )
        call SetGLat( ggi,   117,        720,       1280,  57.18967 )
        call SetGLat( ggi,   118,        720,       1280,  56.90864 )
        call SetGLat( ggi,   119,        729,       1280,  56.62761 )
        call SetGLat( ggi,   120,        750,       1280,  56.34658 )
        call SetGLat( ggi,   121,        750,       1280,  56.06555 )
        call SetGLat( ggi,   122,        750,       1280,  55.78452 )
        call SetGLat( ggi,   123,        750,       1280,  55.50349 )
        call SetGLat( ggi,   124,        768,       1280,  55.22246 )
        call SetGLat( ggi,   125,        768,       1280,  54.94143 )
        call SetGLat( ggi,   126,        768,       1280,  54.66040 )
        call SetGLat( ggi,   127,        768,       1280,  54.37937 )
        call SetGLat( ggi,   128,        800,       1280,  54.09834 )
        call SetGLat( ggi,   129,        800,       1280,  53.81731 )
        call SetGLat( ggi,   130,        800,       1280,  53.53628 )
        call SetGLat( ggi,   131,        800,       1280,  53.25525 )
        call SetGLat( ggi,   132,        800,       1280,  52.97422 )
        call SetGLat( ggi,   133,        800,       1280,  52.69319 )
        call SetGLat( ggi,   134,        810,       1280,  52.41216 )
        call SetGLat( ggi,   135,        810,       1280,  52.13113 )
        call SetGLat( ggi,   136,        864,       1280,  51.85009 )
        call SetGLat( ggi,   137,        864,       1280,  51.56906 )
        call SetGLat( ggi,   138,        864,       1280,  51.28803 )
        call SetGLat( ggi,   139,        864,       1280,  51.00700 )
        call SetGLat( ggi,   140,        864,       1280,  50.72597 )
        call SetGLat( ggi,   141,        864,       1280,  50.44494 )
        call SetGLat( ggi,   142,        864,       1280,  50.16391 )
        call SetGLat( ggi,   143,        864,       1280,  49.88288 )
        call SetGLat( ggi,   144,        864,       1280,  49.60185 )
        call SetGLat( ggi,   145,        864,       1280,  49.32082 )
        call SetGLat( ggi,   146,        864,       1280,  49.03979 )
        call SetGLat( ggi,   147,        900,       1280,  48.75876 )
        call SetGLat( ggi,   148,        900,       1280,  48.47773 )
        call SetGLat( ggi,   149,        900,       1280,  48.19670 )
        call SetGLat( ggi,   150,        900,       1280,  47.91567 )
        call SetGLat( ggi,   151,        900,       1280,  47.63464 )
        call SetGLat( ggi,   152,        900,       1280,  47.35361 )
        call SetGLat( ggi,   153,        900,       1280,  47.07258 )
        call SetGLat( ggi,   154,        900,       1280,  46.79155 )
        call SetGLat( ggi,   155,        960,       1280,  46.51052 )
        call SetGLat( ggi,   156,        960,       1280,  46.22949 )
        call SetGLat( ggi,   157,        960,       1280,  45.94846 )
        call SetGLat( ggi,   158,        960,       1280,  45.66743 )
        call SetGLat( ggi,   159,        960,       1280,  45.38640 )
        call SetGLat( ggi,   160,        960,       1280,  45.10537 )
        call SetGLat( ggi,   161,        960,       1280,  44.82434 )
        call SetGLat( ggi,   162,        960,       1280,  44.54331 )
        call SetGLat( ggi,   163,        960,       1280,  44.26228 )
        call SetGLat( ggi,   164,        960,       1280,  43.98125 )
        call SetGLat( ggi,   165,        960,       1280,  43.70022 )
        call SetGLat( ggi,   166,        960,       1280,  43.41919 )
        call SetGLat( ggi,   167,        960,       1280,  43.13816 )
        call SetGLat( ggi,   168,        960,       1280,  42.85713 )
        call SetGLat( ggi,   169,        972,       1280,  42.57610 )
        call SetGLat( ggi,   170,        972,       1280,  42.29507 )
        call SetGLat( ggi,   171,       1000,       1280,  42.01404 )
        call SetGLat( ggi,   172,       1000,       1280,  41.73301 )
        call SetGLat( ggi,   173,       1000,       1280,  41.45198 )
        call SetGLat( ggi,   174,       1000,       1280,  41.17094 )
        call SetGLat( ggi,   175,       1000,       1280,  40.88991 )
        call SetGLat( ggi,   176,       1000,       1280,  40.60888 )
        call SetGLat( ggi,   177,       1000,       1280,  40.32785 )
        call SetGLat( ggi,   178,       1000,       1280,  40.04682 )
        call SetGLat( ggi,   179,       1024,       1280,  39.76579 )
        call SetGLat( ggi,   180,       1024,       1280,  39.48476 )
        call SetGLat( ggi,   181,       1024,       1280,  39.20373 )
        call SetGLat( ggi,   182,       1024,       1280,  38.92270 )
        call SetGLat( ggi,   183,       1024,       1280,  38.64167 )
        call SetGLat( ggi,   184,       1024,       1280,  38.36064 )
        call SetGLat( ggi,   185,       1080,       1280,  38.07961 )
        call SetGLat( ggi,   186,       1080,       1280,  37.79858 )
        call SetGLat( ggi,   187,       1080,       1280,  37.51755 )
        call SetGLat( ggi,   188,       1080,       1280,  37.23652 )
        call SetGLat( ggi,   189,       1080,       1280,  36.95549 )
        call SetGLat( ggi,   190,       1080,       1280,  36.67446 )
        call SetGLat( ggi,   191,       1080,       1280,  36.39343 )
        call SetGLat( ggi,   192,       1080,       1280,  36.11240 )
        call SetGLat( ggi,   193,       1080,       1280,  35.83137 )
        call SetGLat( ggi,   194,       1080,       1280,  35.55034 )
        call SetGLat( ggi,   195,       1080,       1280,  35.26931 )
        call SetGLat( ggi,   196,       1080,       1280,  34.98828 )
        call SetGLat( ggi,   197,       1080,       1280,  34.70725 )
        call SetGLat( ggi,   198,       1080,       1280,  34.42622 )
        call SetGLat( ggi,   199,       1125,       1280,  34.14519 )
        call SetGLat( ggi,   200,       1125,       1280,  33.86416 )
        call SetGLat( ggi,   201,       1125,       1280,  33.58313 )
        call SetGLat( ggi,   202,       1125,       1280,  33.30210 )
        call SetGLat( ggi,   203,       1125,       1280,  33.02107 )
        call SetGLat( ggi,   204,       1125,       1280,  32.74004 )
        call SetGLat( ggi,   205,       1125,       1280,  32.45901 )
        call SetGLat( ggi,   206,       1125,       1280,  32.17797 )
        call SetGLat( ggi,   207,       1125,       1280,  31.89694 )
        call SetGLat( ggi,   208,       1125,       1280,  31.61591 )
        call SetGLat( ggi,   209,       1125,       1280,  31.33488 )
        call SetGLat( ggi,   210,       1125,       1280,  31.05385 )
        call SetGLat( ggi,   211,       1125,       1280,  30.77282 )
        call SetGLat( ggi,   212,       1125,       1280,  30.49179 )
        call SetGLat( ggi,   213,       1152,       1280,  30.21076 )
        call SetGLat( ggi,   214,       1152,       1280,  29.92973 )
        call SetGLat( ggi,   215,       1152,       1280,  29.64870 )
        call SetGLat( ggi,   216,       1152,       1280,  29.36767 )
        call SetGLat( ggi,   217,       1152,       1280,  29.08664 )
        call SetGLat( ggi,   218,       1152,       1280,  28.80561 )
        call SetGLat( ggi,   219,       1152,       1280,  28.52458 )
        call SetGLat( ggi,   220,       1152,       1280,  28.24355 )
        call SetGLat( ggi,   221,       1152,       1280,  27.96252 )
        call SetGLat( ggi,   222,       1200,       1280,  27.68149 )
        call SetGLat( ggi,   223,       1200,       1280,  27.40046 )
        call SetGLat( ggi,   224,       1200,       1280,  27.11943 )
        call SetGLat( ggi,   225,       1200,       1280,  26.83840 )
        call SetGLat( ggi,   226,       1200,       1280,  26.55737 )
        call SetGLat( ggi,   227,       1200,       1280,  26.27634 )
        call SetGLat( ggi,   228,       1200,       1280,  25.99531 )
        call SetGLat( ggi,   229,       1200,       1280,  25.71428 )
        call SetGLat( ggi,   230,       1200,       1280,  25.43325 )
        call SetGLat( ggi,   231,       1200,       1280,  25.15222 )
        call SetGLat( ggi,   232,       1200,       1280,  24.87119 )
        call SetGLat( ggi,   233,       1200,       1280,  24.59016 )
        call SetGLat( ggi,   234,       1200,       1280,  24.30913 )
        call SetGLat( ggi,   235,       1200,       1280,  24.02810 )
        call SetGLat( ggi,   236,       1200,       1280,  23.74706 )
        call SetGLat( ggi,   237,       1200,       1280,  23.46603 )
        call SetGLat( ggi,   238,       1200,       1280,  23.18500 )
        call SetGLat( ggi,   239,       1200,       1280,  22.90397 )
        call SetGLat( ggi,   240,       1215,       1280,  22.62294 )
        call SetGLat( ggi,   241,       1215,       1280,  22.34191 )
        call SetGLat( ggi,   242,       1215,       1280,  22.06088 )
        call SetGLat( ggi,   243,       1215,       1280,  21.77985 )
        call SetGLat( ggi,   244,       1215,       1280,  21.49882 )
        call SetGLat( ggi,   245,       1215,       1280,  21.21779 )
        call SetGLat( ggi,   246,       1215,       1280,  20.93676 )
        call SetGLat( ggi,   247,       1250,       1280,  20.65573 )
        call SetGLat( ggi,   248,       1250,       1280,  20.37470 )
        call SetGLat( ggi,   249,       1250,       1280,  20.09367 )
        call SetGLat( ggi,   250,       1250,       1280,  19.81264 )
        call SetGLat( ggi,   251,       1250,       1280,  19.53161 )
        call SetGLat( ggi,   252,       1250,       1280,  19.25058 )
        call SetGLat( ggi,   253,       1250,       1280,  18.96955 )
        call SetGLat( ggi,   254,       1250,       1280,  18.68852 )
        call SetGLat( ggi,   255,       1250,       1280,  18.40749 )
        call SetGLat( ggi,   256,       1250,       1280,  18.12646 )
        call SetGLat( ggi,   257,       1250,       1280,  17.84543 )
        call SetGLat( ggi,   258,       1250,       1280,  17.56440 )
        call SetGLat( ggi,   259,       1250,       1280,  17.28337 )
        call SetGLat( ggi,   260,       1250,       1280,  17.00234 )
        call SetGLat( ggi,   261,       1250,       1280,  16.72131 )
        call SetGLat( ggi,   262,       1250,       1280,  16.44028 )
        call SetGLat( ggi,   263,       1250,       1280,  16.15925 )
        call SetGLat( ggi,   264,       1250,       1280,  15.87822 )
        call SetGLat( ggi,   265,       1250,       1280,  15.59718 )
        call SetGLat( ggi,   266,       1250,       1280,  15.31615 )
        call SetGLat( ggi,   267,       1280,       1280,  15.03512 )
        call SetGLat( ggi,   268,       1280,       1280,  14.75409 )
        call SetGLat( ggi,   269,       1280,       1280,  14.47306 )
        call SetGLat( ggi,   270,       1280,       1280,  14.19203 )
        call SetGLat( ggi,   271,       1280,       1280,  13.91100 )
        call SetGLat( ggi,   272,       1280,       1280,  13.62997 )
        call SetGLat( ggi,   273,       1280,       1280,  13.34894 )
        call SetGLat( ggi,   274,       1280,       1280,  13.06791 )
        call SetGLat( ggi,   275,       1280,       1280,  12.78688 )
        call SetGLat( ggi,   276,       1280,       1280,  12.50585 )
        call SetGLat( ggi,   277,       1280,       1280,  12.22482 )
        call SetGLat( ggi,   278,       1280,       1280,  11.94379 )
        call SetGLat( ggi,   279,       1280,       1280,  11.66276 )
        call SetGLat( ggi,   280,       1280,       1280,  11.38173 )
        call SetGLat( ggi,   281,       1280,       1280,  11.10070 )
        call SetGLat( ggi,   282,       1280,       1280,  10.81967 )
        call SetGLat( ggi,   283,       1280,       1280,  10.53864 )
        call SetGLat( ggi,   284,       1280,       1280,  10.25761 )
        call SetGLat( ggi,   285,       1280,       1280,   9.97658 )
        call SetGLat( ggi,   286,       1280,       1280,   9.69555 )
        call SetGLat( ggi,   287,       1280,       1280,   9.41452 )
        call SetGLat( ggi,   288,       1280,       1280,   9.13349 )
        call SetGLat( ggi,   289,       1280,       1280,   8.85246 )
        call SetGLat( ggi,   290,       1280,       1280,   8.57143 )
        call SetGLat( ggi,   291,       1280,       1280,   8.29040 )
        call SetGLat( ggi,   292,       1280,       1280,   8.00937 )
        call SetGLat( ggi,   293,       1280,       1280,   7.72833 )
        call SetGLat( ggi,   294,       1280,       1280,   7.44730 )
        call SetGLat( ggi,   295,       1280,       1280,   7.16627 )
        call SetGLat( ggi,   296,       1280,       1280,   6.88524 )
        call SetGLat( ggi,   297,       1280,       1280,   6.60421 )
        call SetGLat( ggi,   298,       1280,       1280,   6.32318 )
        call SetGLat( ggi,   299,       1280,       1280,   6.04215 )
        call SetGLat( ggi,   300,       1280,       1280,   5.76112 )
        call SetGLat( ggi,   301,       1280,       1280,   5.48009 )
        call SetGLat( ggi,   302,       1280,       1280,   5.19906 )
        call SetGLat( ggi,   303,       1280,       1280,   4.91803 )
        call SetGLat( ggi,   304,       1280,       1280,   4.63700 )
        call SetGLat( ggi,   305,       1280,       1280,   4.35597 )
        call SetGLat( ggi,   306,       1280,       1280,   4.07494 )
        call SetGLat( ggi,   307,       1280,       1280,   3.79391 )
        call SetGLat( ggi,   308,       1280,       1280,   3.51288 )
        call SetGLat( ggi,   309,       1280,       1280,   3.23185 )
        call SetGLat( ggi,   310,       1280,       1280,   2.95082 )
        call SetGLat( ggi,   311,       1280,       1280,   2.66979 )
        call SetGLat( ggi,   312,       1280,       1280,   2.38876 )
        call SetGLat( ggi,   313,       1280,       1280,   2.10773 )
        call SetGLat( ggi,   314,       1280,       1280,   1.82670 )
        call SetGLat( ggi,   315,       1280,       1280,   1.54567 )
        call SetGLat( ggi,   316,       1280,       1280,   1.26464 )
        call SetGLat( ggi,   317,       1280,       1280,   0.98361 )
        call SetGLat( ggi,   318,       1280,       1280,   0.70258 )
        call SetGLat( ggi,   319,       1280,       1280,   0.42155 )
        call SetGLat( ggi,   320,       1280,       1280,   0.14052 )

      case default
        write (*,'("ERROR - unsupported gg resolution : ",i4)') ggi%N
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    
!    ! * lat of cell center for a row
!    allocate( ggi%clat(ggi%nrow) )
!    do j = 1, ggi%nrow
!      if ( j <= ggi%nlat ) then
!        ggi%clat(j) = ggi%glat(j)
!      else
!        ggi%clat(j) = ggi%glat(ggi%nrow+1-j)
!      end if
!    end do
    
    ! * lat of cell bound for a row
    allocate( ggi%blat_deg(0:ggi%nlat) )
    ggi%blat_deg(0) = 90.0
    do j = 1, ggi%nlat-1
      ggi%blat_deg(j) = ( ggi%lat_deg(j) + ggi%lat_deg(j+1) )/2
    end do
    ggi%blat_deg(ggi%nlat) = -90.0
    ! rad
    allocate( ggi%blat(0:ggi%nlat) )
    ggi%blat = ggi%blat_deg * deg2rad
    
    ! * dlat of a cell
    allocate( ggi%dlat_deg(1:ggi%nlat) )
    ggi%dlat_deg = ggi%blat_deg(0:ggi%nlat-1) - ggi%blat_deg(1:ggi%nlat)
    ! rad
    allocate( ggi%dlat(1:ggi%nlat) )
    ggi%dlat = ggi%dlat_deg * deg2rad
    
    ! * total number of cells
    ggi%np = sum( ggi%nlon )
    
    ! * start and end indices for each row
    allocate( ggi%i1(ggi%nlat) )
    allocate( ggi%im(ggi%nlat) )
    do j = 1, ggi%nlat
      ggi%i1(j) = sum(ggi%nlon(1:j-1)) + 1
      ggi%im(j) = sum(ggi%nlon(1:j))
    end do
    
!    ! * row number for each cell
!    allocate( ggi%j(ggi%nrow) )
!    do j = 1, ggi%nrow
!      ggi%j(ggi%i1(j):ggi%im(j)) = j
!    end do
!    
!    ! * lon of cell and center for each cell
!    allocate( ggi%clon(ggi%np) )
!    allocate( ggi%blon(0:ggi%np) )
!    do j = 1, ggi%nrow
!      dlon = 360.0/ggi%nlon(j)
!      ggi%clon(ggi%i1(j):ggi%im(j)) = (/ ((i-1)*dlon, i=1,ggi%nlon(j)) /)
!      ggi%clon(ggi%i1(j):ggi%im(j)) = (/ ((i-0.5)*dlon, i=0,ggi%nlon(j)) /)
!    end do

    ! * area of cell in a row
    ! rad^2 :
    allocate( ggi%area(ggi%nlat) )
    do j = 1, ggi%nlat
      ggi%area(j) = ll_area( 0.0, ggi%dlon(j), ggi%blat(j-1), ggi%blat(j) )
    end do
    ! m^2 :
    allocate( ggi%area_m2(ggi%nlat) )
    ggi%area_m2 = ggi%area * ae**2
    
    ! * mask to use only a few lines on the grid;
    !   by default, use all:
    ggi%latflag = .true.
    
    ! ok
    status = 0

  end subroutine ggi_Init
  
  
  ! =
  
  ! Set definition of Gaussian grid for two oposite latitudes.
  ! This sets on row of the tables in
  !   http://wms.ecmwf.int/documents/manuals/libraries
  !                          /interpolation/gaussianGridsFIS.html
  !
  !    latitude   reduced     regular  latitude
  !     number     points      points
  !    -------    -------     -------  --------
  !        1          18         320   89.14152
  !        2          25         320   88.02943  
  !        3          36         320   86.91077
  !        :
  !
  
  subroutine SetGLat( ggi, j, nlon, nlon_reg, lat_deg )
  
    use Binas, only : deg2rad
  
    ! --- in/out --------------------------------
    
    type(TggGridInfo), intent(inout)     ::  ggi
    integer, intent(in)                  ::  j
    integer, intent(in)                  ::  nlon
    integer, intent(in)                  ::  nlon_reg
    real, intent(in)                     ::  lat_deg
    
    ! --- begin ---------------------------------
    
    ggi%nlon_reg   = nlon_reg

    if ( ggi%reduced ) then
      ggi%nlon(j)    = nlon
    else
      ggi%nlon(j)    = nlon_reg
    end if
    ggi%dlon_deg(j) = 360.0/ggi%nlon(j)
    ggi%dlon(j)     = ggi%dlon_deg(j) * deg2rad
    ggi%lat_deg(j)  = lat_deg
    ggi%lat(j)      = lat_deg * deg2rad
    
    ggi%nlon(ggi%nlat+1-j)     = ggi%nlon(j)
    ggi%dlon_deg(ggi%nlat+1-j) = ggi%dlon_deg(j)
    ggi%dlon(ggi%nlat+1-j)     = ggi%dlon(j)
    ggi%lat_deg(ggi%nlat+1-j)  = -lat_deg
    ggi%lat(ggi%nlat+1-j)      = -lat_deg * deg2rad
    
  end subroutine SetGLat
  
  
  ! ===
  
  
  subroutine ggi_Done( info, status )
  
    ! --- in/out ---------------------------------
    
    type(TggGridInfo), intent(inout)     ::  info
    integer, intent(out)               ::  status
    
    ! --- const ----------------------------------
    
    !character(len=*), parameter  ::  rname = mname//', ggi_Init'
    
    
    ! --- begin ---------------------------------
    
    ! free memory
    deallocate( info%lat )
    deallocate( info%lat_deg )
    deallocate( info%latflag )
    deallocate( info%blat )
    deallocate( info%blat_deg )
    deallocate( info%dlat )
    deallocate( info%dlat_deg )
    deallocate( info%nlon )
    deallocate( info%dlon )
    deallocate( info%dlon_deg )
    deallocate( info%i1 )
    deallocate( info%im )
    deallocate( info%area )
    deallocate( info%area_m2 )

    ! ok
    status = 0

  end subroutine ggi_Done
  
  
  ! ===
  
  
  subroutine GetLons( ggi, j, lons )
  
    use Binas, only : deg2rad
  
    ! --- in/out ------------------------------
    
    type(TggGridInfo), intent(in)        ::  ggi
    integer, intent(in)                  ::  j
    real, intent(out)                    ::  lons(:)
    
    ! --- local ----------------------------
    
    integer        ::  i
    integer        ::  nlon
    
    ! --- begin ------------------------------
    
    if ( j<1 .or. j>ggi%nlat ) then
      print *, 'GetLons : j=',j,'not in range 1,',ggi%nlat
      stop
    end if
    
    nlon = ggi%nlon(j)
    
    if ( size(lons) /= nlon ) then
      print *, 'GetLons : size(lons)=',size(lons), 'while nlon(',j,')=',nlon
      stop
    end if
    
    do i = 1, nlon
      lons(i) = 360.0*(i-1)/nlon  * deg2rad
    end do

  end subroutine GetLons    
  

  ! ===
  
  
  subroutine gg_Longitudes( ggi, gg )
  
    use Binas, only : deg2rad
  
    ! --- in/out ------------------------------
    
    type(TggGridInfo), intent(in)        ::  ggi
    real, intent(out)                    ::  gg(:)  ! deg
    
    ! --- local ----------------------------
    
    integer        ::  j
    
    ! --- begin ------------------------------
    
    call Check( ggi, gg )
    
    ! loop over rows
    do j = 1, ggi%nlat
      ! fill lons for this row:
      call GetLons( ggi, j, gg(ggi%i1(j):ggi%im(j)) )  ! rad
    end do
    
    ! convert to degrees:
    gg = gg / deg2rad      ! deg
    
  end subroutine gg_Longitudes    
  

  ! ===
  
  
  subroutine gg_Latitudes( ggi, gg )
  
    ! --- in/out ------------------------------
    
    type(TggGridInfo), intent(in)        ::  ggi
    real, intent(out)                    ::  gg(:)  ! deg
    
    ! --- local ----------------------------
    
    integer        ::  j
    
    ! --- begin ------------------------------
    
    call Check( ggi, gg )
    
    ! loop over rows
    do j = 1, ggi%nlat
      ! fill latitude for this row:
      gg(ggi%i1(j):ggi%im(j)) = ggi%lat_deg(j)   ! deg
    end do
    
  end subroutine gg_Latitudes    
  

  ! =============================================================
  
  
  subroutine gg_Check( ggi, gg )
  
    ! --- in/out ----------------------------------
    
    type(TggGridInfo), intent(in)    ::  ggi
    real, intent(in)                 ::  gg(:)
    
    ! --- begin ----------------------------------
    
    ! check size of data:
    if ( size(gg) /= ggi%np ) then
      print *, 'ggrid_Check : data size',size(gg),'while expected',ggi%np
      stop
    end if

  end subroutine gg_Check


  ! =====================================================
  
  
  subroutine gg_AreaOper( ggi, gg, oper, unit )
  
    ! --- in/out ----------------------------------
    
    type(TggGridInfo), intent(in)           ::  ggi
    real, intent(inout)                     ::  gg(:)
    character(len=*), intent(in)            ::  unit, oper
    
    ! --- local --------------------------------
    
    integer   ::  j
    integer   ::  i1, im
    real      ::  cell_area
    
    ! --- begin ----------------------------------
    
    call Check( ggi, gg )
    
    do j = 1, ggi%nlat
      
      ! select correct area for cells in this row:
      select case ( unit )
        case ( 'rad2' )
          cell_area = ggi%area(j)
        case ( 'm2' )
          cell_area = ggi%area_m2(j)
        case default
          print *, 'gg_AreaAverage : unknown unit "'//trim(unit)//'"'
          stop
      end select
      
      ! range of cells in this row:
      i1 = ggi%i1(j)
      im = ggi%im(j)
      
      ! assign or modify with cell area:
      select case ( oper )
        case ( '=' )
          gg(i1:im) = cell_area
        case ( '/' )
          gg(i1:im) = gg(i1:im) / cell_area
        case ( '*' )
          gg(i1:im) = gg(i1:im) * cell_area
        case default
          print *, 'gg_AreaAverage : unknown operation "'//trim(oper)//'"'
          stop
      end select

    end do
    
  end subroutine gg_AreaOper


  ! ===================================================
  
  !
  ! NOTE:
  !   o order : lat, lon !
  !   o unit  : degrees !
  !
  
  real function gg_Eval_lat_lon( ggi, gg, lat, lon )
  
    use Binas, only : deg2rad, pi
    
    use Num, only : Interp_Lin, CircInterp_Lin

    ! --- in/out -------------------------
    
    type(TggGridInfo), intent(in)   ::  ggi
    real, intent(in)                ::  gg(:)
    real, intent(in)                ::  lon, lat      ! rad
    
    ! --- local --------------------------
    
    real                  ::  lonX
    integer               ::  nlon, nlon_max
    integer               ::  nlat
    integer               ::  i1, im
    
    real, allocatable     ::  lons(:)
    
    real, allocatable     ::  gl(:)
    real, allocatable     ::  gl_lat(:)
    
    integer               ::  i, j
    integer               ::  j1, jm
    integer               ::  ilast
    integer               ::  status
    
    real                  ::  res
    
    ! --- begin --------------------------------
    
    call Check( ggi, gg )
    
    ! lon in [0,2pi)
    lonX = lon
    if ( lonX < 0 ) lonX = lonX + 2*pi

    nlat = ggi%nlat
    nlon_max = maxval(ggi%nlon)
    
    ! ll in lon, gg in lat
    allocate( gl(0:nlat+1) ); gl = 0.0

    ! latitudes from northpole to southpole:
    allocate( gl_lat(0:nlat+1) )
    gl_lat(0)      =  90.0 * deg2rad       ! north pole (rad)
    gl_lat(1:nlat) = ggi%lat               ! (rad)
    gl_lat(nlat+1) = -90.0 * deg2rad       ! north pole (rad)

    ! row in gg grid; doubled from -360.0 to 360.0
    allocate( lons(nlon_max) ); lons = 0.0
    
    ! select first and last Gaussian lat:
    j1 = nlat
    do
      if ( (j1 == 1) .or. (ggi%lat(j1) > lat) ) exit
      j1 = j1 - 1
    end do
    jm = 1
    do
      if ( (jm == nlat) .or. (ggi%lat(jm) < lat) ) exit
      jm = jm + 1
    end do

    ! loop over Gaussian latitudes, from north to south!
    do j = j1, jm

      ! number of lon points at this latitude:
      nlon = ggi%nlon(j)

      ! start and end indices in grid point array      
      i1 = ggi%i1(j)
      im = ggi%im(j)

      ! lons in [0,2pi)
      do i = 1, nlon
        lons(i) = (i-1)*360.0/nlon
      end do
      lons(1:nlon) = lons(1:nlon) * deg2rad

      ! set north pole  (j=1 in gg, j=nlat+1 in gl)
      if ( j == 1 ) then
        gl(0) = sum(gg(i1:im))/nlon
      end if
      
      ! set south pole
      if ( j == nlat ) then
        gl(nlat+1) = sum(gg(i1:im))/nlon
      end if

      ! Interpolate over lon (circular arrays):
      ilast = 1
      call CircInterp_Lin( lons(1:nlon), 360.0*deg2rad, gg(i1:im), &
                           lonX, gl(j), ilast, status )
      if (status/=0) stop 'ERROR in gg_Eval_lat_lon'
    end do

    ! Linear interpolation over lat; 
    ! negate lat axis to have increasing ax:
    ilast = jm
    call Interp_Lin( -1.0*gl_lat, gl, -1.0*lat, res, ilast, status )
    if (status/=0) stop 'ERROR in gg_Eval_lat_lon'

    gg_Eval_lat_lon = res

    ! free memory
    deallocate( gl )
    deallocate( gl_lat )
    deallocate( lons )

  end function gg_Eval_lat_lon
  
    
  ! ===================================================
  
  
  ! Compute nabla.gg

  !
  !  d/dx for irregular discrete x:
  !
  !
  ! ----+-------+----+-------------------
  !     x-h1     x   x+h2
  !     
  !  
  !                             2   2
  !              h2           h2 -h1         h1
  !            - -- f(x-h1) + ------- f(x) + -- f(x+h2)
  !              h1            h1 h2         h2
  !   f'(x) ~  ----------------------------------------
  !                         h1 + h2
  !
  !                2            2   2          2
  !            - h2 f(x-h1) +(h2 -h1 )f(x) + h1 f(x+h2)
  !   f'(x) ~  ----------------------------------------
  !                         h1 h2 (h1+h2)
  !
  !              2                   2
  !            h2 [f(x)-f(x-h1)] + h1 [f(x+h2)-f(x)]
  !   f'(x) ~  -------------------------------------
  !                         h1 h2 (h1+h2)
  !
  !            h2                 h1  
  !            --[f(x)-f(x-h1)] + -- [f(x+h2)-f(x)]
  !            h1                 h2
  !   f'(x) ~  -------------------------------------
  !                          h1 + h2
  !
  !              h2  f(x)-f(x-h1)    h1   f(x+h2)-f(x)
  !   f'(x) ~  ----- ------------ + ----- ------------
  !            h1+h2      h1        h1+h2       h2
  !


!  subroutine Nabla_gg( ggi, gg, nab )
!  
!    ! --- in/out -----------------------------------
!    
!    type(TggGridInfo), intent(in)    ::  ggi
!    real, intent(in)                 ::  gg(:)
!    real, intent(out)                ::  nab(:,:)
!    
!    ! --- local -----------------------------------
!    
!    integer     ::  i, j
!    integer     ::  i1, im
!    real        ::  h1, h2
!    integer     ::  i1n, imn, nin
!    
!    ! --- begin -----------------------------------
!    
!    if ( size(nab,2) /= 2 ) then
!      print *, 'size nab should be (np,2)'
!      stop 'error in Nabla_gg'
!    end if
!    call Check( ggi, gg )
!    call Check( ggi, nab(:,1) )
!    call Check( ggi, nab(:,2) )
!    
!    ! loop from north to south (!)
!    do j = 1, ggi%nlat
!    
!      i1 = ggi%i1(j)
!      im = ggi%im(j)
!    
!      ! dlon
!      nab(i1,1) = ( gg(i1+1) - gg(im) ) / (2*dlon)
!      do i = i1+1, im-1
!        nab(i,1) = ( gg(i+1) - gg(i-1) ) / (2*dlon)
!      end do
!      nab(im,1) = ( gg(i1) - gg(im-1) ) / (2*dlon)
!      
!      ! dlat
!      stop 'poles'
!
!      ! previous row
!      if ( j == 1 ) then
!        ...
!      else
!        h1 = ggi%glat(j) - ggi%glat(j-1)
!        i1p = ggi%i1(j-1)
!        imp = ggi%im(j-1)
!        nip = imp - i1p + 2
!        xx1(0:nin) = (/ -ggi%lon(i1n+1), ggi%lon(i1n:imn) /)
!        yy1(0:nin) = (/       gg(imn)  ,      gg(i1n:imn) /)
!      end if
!
!      ! next row      
!      if ( j == ggi%nlat ) then
!        ...
!      then
!        h2 = ggi%glat(j+1) - ggi%glat(j)
!        i1n = ggi%i1(j+1)
!        imn = ggi%im(j+1)
!        nin = imn - i1n + 2
!        xx2(0:nin) = (/ -ggi%lon(i1n+1), ggi%lon(i1n:imn) /)
!        yy2(0:nin) = (/ gg(imn), gg(i1n:imn) /)
!      end if
!      
!      do i = i1, im
!
!        y = gg(i)
!
!        ! interpolate previous and next row:
!        call linterp( xx1, yy1, lon(i), y1, ilast_prev )
!        call linterp( xx2, yy2, lon(i), y2, ilast_next )
!
!        !            h2                 h1  
!        !            --[f(x)-f(x-h1)] + -- [f(x+h2)-f(x)]
!        !            h1                 h2
!        !   f'(x) ~  -------------------------------------
!        !                          h1 + h2
!        !
!    
!        nab(i,2) = ( (y-y1)*h2/h1 + (y2-y)*h1/h2 ) / (h1+h2)
!        
!      end do
!      
!    end do
!    
!  end subroutine Nabla_gg
!  
  
  !
  !  (d/dx,d/dy) . (u,v)  ->  du/dx + dv/dy
  !
  
  
  subroutine Divergence_gg( ggi, u, v, div )
  
    use Binas, only : deg2rad
    use Binas, only : ae
    use Num, only : Interp_Lin

    ! --- in/out -------------------------------
    
    type(TggGridInfo), intent(in)  ::  ggi
    real, intent(in)               ::  u(:), v(:)
    real, intent(out)              ::  div(:)

    ! --- local --------------------------------
    
    integer              ::  status
    integer              ::  j
    integer              ::  i1, im
    integer              ::  i, im1, ip1
    integer              ::  nlon, nlon_max, n
    integer              ::  nlon_curr
    real                 ::  lat
    real, allocatable    ::  uu(:)
    real                 ::  dlon
    real, allocatable    ::  lons(:), lons_temp(:)
    real, allocatable    ::  vv_prev(:), vv_curr(:), vv_next(:), vv_temp(:)
    real                 ::  dlat_prev, dlat_next
    real                 ::  du_dx, dv_dy
    
    ! --- begin ---------------------------------
    
    call Check( ggi, u )
    call Check( ggi, v )
    call Check( ggi, div )
    
    nlon_max = maxval( ggi%nlon )
    allocate( lons(nlon_max) )
    allocate( lons_temp(nlon_max) )
    allocate( uu(0:nlon_max+1) )
    allocate( vv_prev(nlon_max) )
    allocate( vv_curr(nlon_max) )
    allocate( vv_next(nlon_max) )
    allocate( vv_temp(nlon_max) )

    do j = 1, ggi%nlat
    
      nlon = ggi%nlon(j)
      i1 = ggi%i1(j)
      im = ggi%im(j)
      
      lat = ggi%lat(j)

      call GetLons( ggi, j, lons(1:nlon) )

      ! current row of u
      !uu(0:nlon+1) = (/ u(im), u(i1:im), u(i1) /)
      !>>> adhoc fix after xlf90 in debug mode complained:
      uu(0)      = u(im)
      uu(1:nlon) = u(i1:im)
      uu(nlon+1) = u(i1)
      !<<<
      
      ! current row of v
      vv_curr(1:nlon) = v(i1:im)
      
      ! previous row of v
      if ( j == 1 ) then
        n = 1
        vv_prev(1:nlon) = sum(vv_curr(1:nlon))/nlon
        dlat_prev = 90.0*deg2rad - ggi%lat(j)
      else
        n = ggi%nlon(j-1)
        vv_temp(1:n) = v(ggi%i1(j-1):ggi%im(j-1))
        if ( n == nlon ) then
          vv_prev(1:nlon) = vv_temp(1:n)
        else
          call GetLons( ggi, j-1, lons_temp(1:n) )
          call Interp_Lin( (/ lons_temp(1:n), lons_temp(1) /), &
                           (/ vv_temp(1:n), vv_temp(1) /), &
                           lons(1:nlon), vv_prev(1:nlon), status )
          if (status/=0) stop 'ERROR in Divergence_gg'
        end if
        dlat_prev = ggi%lat(j-1) - ggi%lat(j)
      end if
      
      ! next row of v
      if ( j == ggi%nlat ) then
        n = 1
        vv_next(1:nlon) = sum(vv_curr(1:nlon))/nlon
        dlat_next = ggi%lat(j) + 90.0*deg2rad
      else
        n = ggi%nlon(j+1)
        vv_temp(1:n) = v(ggi%i1(j+1):ggi%im(j+1))
        if ( n == nlon ) then
          vv_next(1:nlon) = vv_temp(1:n)
        else
          call GetLons( ggi, j+1, lons_temp(1:n) )
          call Interp_Lin( (/ lons_temp(1:n), lons_temp(1) /), &
                           (/ vv_temp(1:n), vv_temp(1) /), &
                           lons(1:nlon), vv_next(1:nlon), status )
          if (status/=0) stop 'ERROR in Divergence_gg'
        end if
        dlat_next = ggi%lat(j) - ggi%lat(j+1)
      end if
      
      ! loop over all cells in this row
      do i = 1, nlon
      
        dlon = ggi%dlon(j)
        du_dx = ( uu(i+1) - uu(i-1) ) / (2*dlon*cos(lat)*ae)
        
        !            h2                 h1  
        !            --[f(x)-f(x-h1)] + -- [f(x+h2)-f(x)]
        !            h1                 h2
        !   f'(x) ~  -------------------------------------
        !                          h1 + h2

        dv_dy = ( (vv_curr(i)-vv_prev(i))*dlat_next/dlat_prev + &
                  (vv_next(i)-vv_curr(i))*dlat_prev/dlat_next ) &
                 / ( (dlat_prev+dlat_next)*ae )

        div(i1-1+i) = du_dx + dv_dy

      end do
      
    end do

    deallocate( lons )
    deallocate( lons_temp )
    deallocate( uu )
    deallocate( vv_prev )
    deallocate( vv_curr )
    deallocate( vv_next )
    deallocate( vv_temp )

  end subroutine Divergence_gg
  
    

  !
  !  (d/dx,d/dy) . g   ->  (dg/dx,dg/dy)
  !
  
  
  subroutine Nabla_gg( ggi, gg, dgdx, dgdy )
  
    use Binas, only : deg2rad
    use Binas, only : ae
    use Num, only : Interp_Lin

    ! --- in/out -------------------------------
    
    type(TggGridInfo), intent(in)  ::  ggi
    real, intent(in)               ::  gg(:)
    real, intent(out)              ::  dgdx(:), dgdy(:)

    ! --- local --------------------------------
    
    integer              ::  status
    integer              ::  j
    integer              ::  i1, im
    integer              ::  i, im1, ip1
    integer              ::  nlon, nlon_max, n
    integer              ::  nlon_curr
    real                 ::  lat
    real, allocatable    ::  uu(:)
    real                 ::  dlon
    real, allocatable    ::  lons(:), lons_temp(:)
    real, allocatable    ::  vv_prev(:), vv_curr(:), vv_next(:), vv_temp(:)
    real                 ::  dlat_prev, dlat_next
    real                 ::  du_dx, dv_dy
    
    ! --- begin ---------------------------------
    
    call Check( ggi, gg )
    call Check( ggi, dgdx )
    call Check( ggi, dgdy )
    
    nlon_max = maxval( ggi%nlon )
    allocate( lons(nlon_max) )
    allocate( lons_temp(nlon_max) )
    allocate( uu(0:nlon_max+1) )
    allocate( vv_prev(nlon_max) )
    allocate( vv_curr(nlon_max) )
    allocate( vv_next(nlon_max) )
    allocate( vv_temp(nlon_max) )

    do j = 1, ggi%nlat
    
      nlon = ggi%nlon(j)
      i1 = ggi%i1(j)
      im = ggi%im(j)
      
      lat = ggi%lat(j)

      call GetLons( ggi, j, lons(1:nlon) )

      ! current row of u
      uu(0:nlon+1) = (/ gg(im), gg(i1:im), gg(i1) /)
      
      ! current row of v
      vv_curr(1:nlon) = gg(i1:im)
      
      ! previous row of v
      if ( j == 1 ) then
        n = 1
        vv_prev(1:nlon+1) = sum(vv_curr(1:nlon))/nlon
        dlat_prev = 90.0*deg2rad - ggi%lat(j)
      else
        n = ggi%nlon(j-1)
        vv_temp(1:n) = gg(ggi%i1(j-1):ggi%im(j-1))
        if ( n == nlon ) then
          vv_prev(1:nlon) = vv_temp(1:n)
        else
          call GetLons( ggi, j-1, lons_temp(1:n) )
          call Interp_Lin( (/ lons_temp(1:n), lons_temp(1) /), &
                           (/ vv_temp(1:n), vv_temp(1) /), &
                           lons(1:nlon), vv_prev(1:nlon), status )
          if (status/=0) stop 'ERROR in Nabla_gg'
        end if
        dlat_prev = ggi%lat(j-1) - ggi%lat(j)
      end if
      
      ! next row of v
      if ( j == ggi%nlat ) then
        n = 1
        vv_next(1:nlon) = sum(vv_curr(1:nlon))/nlon
        dlat_next = ggi%lat(j) + 90.0*deg2rad
      else
        n = ggi%nlon(j+1)
        vv_temp(1:n) = gg(ggi%i1(j+1):ggi%im(j+1))
        if ( n == nlon ) then
          vv_next(1:nlon) = vv_temp(1:n)
        else
          call GetLons( ggi, j+1, lons_temp(1:n) )
          call Interp_Lin( (/ lons_temp(1:n), lons_temp(1) /), &
                           (/ vv_temp(1:n), vv_temp(1) /), &
                           lons(1:nlon), vv_next(1:nlon), status )
          if (status/=0) stop 'ERROR in Nabla_gg'
        end if
        dlat_next = ggi%lat(j) - ggi%lat(j+1)
      end if
      
      ! loop over all cells in this row
      do i = 1, nlon
      
        dlon = ggi%dlon(j)
        du_dx = ( uu(i+1) - uu(i-1) ) / (2*dlon*cos(lat)*ae)
        
        !            h2                 h1  
        !            --[f(x)-f(x-h1)] + -- [f(x+h2)-f(x)]
        !            h1                 h2
        !   f'(x) ~  -------------------------------------
        !                          h1 + h2

        dv_dy = ( (vv_curr(i)-vv_prev(i))*dlat_next/dlat_prev + &
                  (vv_next(i)-vv_curr(i))*dlat_prev/dlat_next ) &
                 / ( (dlat_prev+dlat_next)*ae )

        dgdx(i1-1+i) = du_dx
        dgdy(i1-1+i) = dv_dy

      end do
      
    end do

    deallocate( lons )
    deallocate( lons_temp )
    deallocate( uu )
    deallocate( vv_prev )
    deallocate( vv_curr )
    deallocate( vv_next )
    deallocate( vv_temp )

  end subroutine Nabla_gg
  
    

  
end module grid_type_gg
    

