!
!  Module and data types:
!
!    use grid_type_hyb
!
!    type(TLevelInfo)    ::  levi
!    type(TLevelInfo)    ::  levi_parent
!
!    integer              ::  status
!
!  Initialize in one of the following ways:
!
!   o by half level coeff  a (Pa)  and   b (0-1)
!
!      call Init( levi, 60, (/a0,..,a60/), (/b0,..,b60/), status )
!
!   o by a character key; yet supported:
!       'ec19'   : echam 19 levels
!       'ec31'   : ecmwf 31 levels
!       'ec40'   : ecmwf 40 levels
!       'ec60'   : ecmwf 60 levels
!       'ec62'   : ecmwf 62 levels
!       'ec91'   : ecmwf 91 levels
!       'tm31'   : tm 31 levels (reversed ec31)
!       'tm40'   : tm 40 levels (reversed ec40)
!       'tm60'   : tm 60 levels (reversed ec60)
!       'tm62'   : tm 62 levels (reversed ec62)
!       'tm91'   : tm 91 levels (reversed ec91)
!      call Init( levi, 'ec60', status )
!
!   o define a selection of half levels:
!
!      call Init( levi_parent, 'ec60', status )
!      call Init( levi, levi_parent, (/0,2,..,58,60/), status )
!
!  Copy 3D array llX into ll,
!  How levels are combined is specified by a key:
!   'bottom'     :  use bottom value (most close to the ground)
!   'top'        :  use top value (most close to the model top)
!   'sum'        :  sum values
!   'aver'       :  average of all levels
!   'mass-aver'  : mass weighted average
!  Eventually reverse.
!
!      call FillLevels( levi, ll, leviX, llX, combine_key, status )
!
!
!

module grid_type_hyb

  implicit none

  ! --- in/out --------------------------

  private

  public  ::  TLevelInfo
  public  ::  Init, Done
  public  ::  Check
  public  ::  Compare
  public  ::  HPressure, FPressure
  public  ::  FillLevels

  ! --- const ---------------------------------------

  character(len=*), parameter  ::  mname = 'grid_type_hyb'


  ! --- types ---------------------------

  type TLevelInfo
    ! name used for messages:
    character(len=32)                   ::  name
    ! number of levels
    integer                             ::  nlev
    ! hybride half level coeff; indices 0, 1, ..., nlev
    real, pointer                       ::  a(:), b(:)
    ! hybride full level coeff; indices 1, ..., nlev
    real, pointer                       ::  fa(:), fb(:)
    ! half level coeff per full level (1:2,1:nlev)
    real, pointer                       ::  fa_bounds(:,:), fb_bounds(:,:)
    ! standard pressures and mass (surface pressure 1e5)
    real, pointer                       ::  p0(:), fp0(:), m0(:)
    ! upwards od downwards ?
    character(len=1)                    ::  updo
    ! hybride increments; indices 1, ..., nlev
    real, pointer                       ::  da(:), db(:)
    ! -- parent levels
    logical                             ::  selection
    integer, pointer                    ::  hlevs(:)
    integer, pointer                    ::  flevs(:,:)
    integer                             ::  nlev_parent
    real, pointer                       ::  a_parent(:), b_parent(:)
  end type TLevelInfo


  ! --- const ------------------------------------------------------------------

  ! *** ec19 **************************************************************

  ! Hybrid coordinate coefficients at half level interfaces,
  ! specifying 19 vertical ECHAM levels.
  ! Coefficient a is given in Pa, b in [0,1] .

  real, parameter :: a_ec19(0:19) = (/ &
             0.000   ,  2000.000   ,  4000.000   ,  6046.111   , &
          8267.928   , 10609.510   , 12851.100   , 14698.500   , &
         15861.130   , 16116.240   , 15356.920   , 13621.460   , &
         11101.560   ,  8127.144   ,  5125.142   ,  2549.969   , &
           783.195   ,     0.000   ,     0.000   ,     0.000        /)

  real, parameter :: b_ec19(0:19) = (/ &
             0.000000000, 0.00000000, 0.00000000, 0.0003389933, &
             0.003357187, 0.01307004, 0.03407715, 0.0706498300, &
             0.125916700, 0.20119540, 0.29551960, 0.4054092000, &
             0.524932200, 0.64610790, 0.75969840, 0.8564376000, &
             0.928746900, 0.97298520, 0.99228150, 1.0000000000    /)


  ! *** ec31 **************************************************************

  ! Hybrid coordinate coefficients at half level interfaces,
  ! specifying 31  vertical ECMWF levels.
  ! Coefficient a is given in Pa, b in [0,1] .

  real, parameter :: a_ec31(0:31) = (/ &
        0.0         ,  2000.0         ,  4000.0         ,  6000.0         , &
     8000.0         ,  9976.135361    , 11820.539617    , 13431.393926    , &
    14736.356909    , 15689.207458    , 16266.610500    , 16465.005734    , &
    16297.619332    , 15791.598604    , 14985.269630    , 13925.517858    , &
    12665.291662    , 11261.228878    ,  9771.406290    ,  8253.212096    , &
     6761.341326    ,  5345.914240    ,  4050.717678    ,  2911.569385    , &
     1954.805296    ,  1195.889791    ,   638.148911    ,   271.626545    , &
       72.063577    ,     0.0         ,     0.0         ,     0.0             /)

  real, parameter :: b_ec31(0:31) = (/ &
        0.0         , 0.0             , 0.0             , 0.0             , &
        0.0         , 0.0003908582    , 0.0029197006    , 0.0091941320    , &
        0.0203191555, 0.0369748598    , 0.0594876397    , 0.0878949492    , &
        0.1220035886, 0.1614415235    , 0.2057032385    , 0.2541886223    , &
        0.3062353873, 0.3611450218    , 0.4182022749    , 0.4766881754    , &
        0.5358865832, 0.5950842740    , 0.6535645569    , 0.7105944258    , &
        0.7654052430, 0.8171669567    , 0.8649558510    , 0.9077158297    , &
        0.9442132326, 0.9729851852    , 0.9922814815    , 1.0                 /)

  ! *** ec40 **************************************************************

  ! Hybrid coordinate coefficients at half level interfaces,
  ! specifying 40 vertical ECMWF levels.
  ! Coefficient a is given in [Pa] .

  real, parameter  :: a_ec40(0:40) = (/   &
             0.000000,  2000.000000,  4000.000000,  6000.000000, &
          8000.000000,  9988.882813, 11914.523438, 13722.941406, &
         15369.730469, 16819.476563, 18045.183594, 19027.695313, &
         19755.109375, 20222.205078, 20429.863281, 20384.480469, &
         20097.402344, 19584.330078, 18864.750000, 17961.357422, &
         16899.468750, 15706.447266, 14411.124023, 13043.218750, &
         11632.758789, 10209.500977,  8802.356445,  7438.803223, &
          6144.314941,  4941.778320,  3850.913330,  2887.696533, &
          2063.779785,  1385.912598,   855.361755,   467.333588, &
           210.393890,    65.889244,     7.367743,     0.000000, &
             0.000000/)

  real, parameter  :: b_ec40(0:40) = (/   &
             0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,  &
             0.0000000000, 0.0001971156, 0.0015112918, 0.0048841573,  &
             0.0110761747, 0.0206778906, 0.0341211632, 0.0516904071,  &
             0.0735338330, 0.0996747017, 0.1300225258, 0.1643843055,  &
             0.2024759650, 0.2439331412, 0.2883229256, 0.3351548910,  &
             0.3838921785, 0.4339629412, 0.4847716093, 0.5357099175,  &
             0.5861684084, 0.6355474591, 0.6832686067, 0.7287858129,  &
             0.7715966105, 0.8112534285, 0.8473749161, 0.8796569109,  &
             0.9078838825, 0.9319403172, 0.9518215060, 0.9676452279,  &
             0.9796627164, 0.9882701039, 0.9940194488, 0.9976301193,  &
             1.0000000000 /)


  ! *** ec60 **************************************************************

  ! Hybrid coordinate coefficients at half level interfaces,
  ! specifying 60  vertical ECMWF levels.
  ! Coefficient a is given in Pa, b in [0,1] .

  real, parameter  :: a_ec60(0:60) = (/   &
             0.000000,    20.000000,    38.425343,    63.647804, &
            95.636963,   134.483307,   180.584351,   234.779053, &
           298.495789,   373.971924,   464.618134,   575.651001, &
           713.218079,   883.660522,  1094.834717,  1356.474609, &
          1680.640259,  2082.273926,  2579.888672,  3196.421631, &
          3960.291504,  4906.708496,  6018.019531,  7306.631348, &
          8765.053711, 10376.126953, 12077.446289, 13775.325195, &
         15379.805664, 16819.474609, 18045.183594, 19027.695313, &
         19755.109375, 20222.205078, 20429.863281, 20384.480469, &
         20097.402344, 19584.330078, 18864.750000, 17961.357422, &
         16899.468750, 15706.447266, 14411.124023, 13043.218750, &
         11632.758789, 10209.500977,  8802.356445,  7438.803223, &
          6144.314941,  4941.778320,  3850.913330,  2887.696533, &
          2063.779785,  1385.912598,   855.361755,   467.333588, &
           210.393890,    65.889244,     7.367743,     0.000000, &
             0.000000/)

  real, parameter  :: b_ec60(0:60) = (/   &
             0.00000000, 0.00000000, 0.00000000, 0.00000000,  &
             0.00000000, 0.00000000, 0.00000000, 0.00000000,  &
             0.00000000, 0.00000000, 0.00000000, 0.00000000,  &
             0.00000000, 0.00000000, 0.00000000, 0.00000000,  &
             0.00000000, 0.00000000, 0.00000000, 0.00000000,  &
             0.00000000, 0.00000000, 0.00000000, 0.00000000,  &
             0.00007582, 0.00046139, 0.00181516, 0.00508112,  &
             0.01114291, 0.02067788, 0.03412116, 0.05169041,  &
             0.07353383, 0.09967469, 0.13002251, 0.16438432,  &
             0.20247590, 0.24393314, 0.28832296, 0.33515489,  &
             0.38389215, 0.43396294, 0.48477158, 0.53570992,  &
             0.58616841, 0.63554746, 0.68326861, 0.72878581,  &
             0.77159661, 0.81125343, 0.84737492, 0.87965691,  &
             0.90788388, 0.93194032, 0.95182151, 0.96764523,  &
             0.97966272, 0.98827010, 0.99401945, 0.99763012,  &
             1.00000000 /)

  ! *** ec62 **************************************************************

  real, parameter  :: a_ec62(0:62) = (/   &
                                                                       0.000000, &
         988.835876,   1977.676270,   2966.516602,   3955.356934,   4944.197266, &
        5933.037598,   6921.870117,   7909.441406,   8890.707031,   9860.528320, &
       10807.783203,  11722.749023,  12595.006836,  13419.463867,  14192.009766, &
       14922.685547,  15638.053711,  16329.560547,  16990.623047,  17613.281250, &
       18191.029297,  18716.968750,  19184.544922,  19587.513672,  19919.796875, &
       20175.394531,  20348.916016,  20434.158203,  20426.218750,  20319.011719, &
       20107.031250,  19785.357422,  19348.775391,  18798.822266,  18141.296875, &
       17385.595703,  16544.585938,  15633.566406,  14665.645508,  13653.219727, &
       12608.383789,  11543.166992,  10471.310547,   9405.222656,   8356.252930, &
        7335.164551,   6353.920898,   5422.802734,   4550.215820,   3743.464355, &
        3010.146973,   2356.202637,   1784.854614,   1297.656128,    895.193542, &
         576.314148,    336.772369,    162.043427,     54.208336,      6.575628, &
           0.003160,      0.000000 /)

  real, parameter  :: b_ec62(0:62) = (/   &
                                                                       0.000000, &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000000,      0.000000,      0.000013,      0.000087,      0.000275, &
           0.000685,      0.001415,      0.002565,      0.004187,      0.006322, &
           0.009035,      0.012508,      0.016860,      0.022189,      0.028610, &
           0.036227,      0.045146,      0.055474,      0.067316,      0.080777, &
           0.095964,      0.112979,      0.131935,      0.152934,      0.176091, &
           0.201520,      0.229315,      0.259554,      0.291993,      0.326329, &
           0.362203,      0.399205,      0.436906,      0.475016,      0.513280, &
           0.551458,      0.589317,      0.626559,      0.662934,      0.698224, &
           0.732224,      0.764679,      0.795385,      0.824185,      0.850950, &
           0.875518,      0.897767,      0.917651,      0.935157,      0.950274, &
           0.963007,      0.973466,      0.982238,      0.989153,      0.994204, &
           0.997630,      1.000000 /)


  ! *** ec91 **************************************************************

  ! Hybrid coordinate coefficients at half level interfaces,
  ! specifying 91 vertical ECMWF levels.
  ! Coefficient a is given in Pa, b in [0,1] .

  real, parameter  :: a_ec91(0:91) = (/   &
           0.000000,      2.000040,      3.980832,      7.387186,     12.908319, &
          21.413612,     33.952858,     51.746601,     76.167656,    108.715561, &
         150.986023,    204.637451,    271.356506,    352.824493,    450.685791, &
         566.519226,    701.813354,    857.945801,   1036.166504,   1237.585449, &
        1463.163940,   1713.709595,   1989.874390,   2292.155518,   2620.898438, &
        2976.302246,   3358.425781,   3767.196045,   4202.416504,   4663.776367, &
        5150.859863,   5663.156250,   6199.839355,   6759.727051,   7341.469727, &
        7942.926270,   8564.624023,   9208.305664,   9873.560547,  10558.881836, &
       11262.484375,  11982.662109,  12713.897461,  13453.225586,  14192.009766, &
       14922.685547,  15638.053711,  16329.560547,  16990.623047,  17613.281250, &
       18191.029297,  18716.968750,  19184.544922,  19587.513672,  19919.796875, &
       20175.394531,  20348.916016,  20434.158203,  20426.218750,  20319.011719, &
       20107.031250,  19785.357422,  19348.775391,  18798.822266,  18141.296875, &
       17385.595703,  16544.585938,  15633.566406,  14665.645508,  13653.219727, &
       12608.383789,  11543.166992,  10471.310547,   9405.222656,   8356.252930, &
        7335.164551,   6353.920898,   5422.802734,   4550.215820,   3743.464355, &
        3010.146973,   2356.202637,   1784.854614,   1297.656128,    895.193542, &
         576.314148,    336.772369,    162.043427,     54.208336,      6.575628, &
           0.003160,      0.000000 /)

  real, parameter  :: b_ec91(0:91) = (/   &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000000,      0.000000,      0.000000,      0.000000,      0.000000, &
           0.000014,      0.000055,      0.000131,      0.000279,      0.000548, &
           0.001000,      0.001701,      0.002765,      0.004267,      0.006322, &
           0.009035,      0.012508,      0.016860,      0.022189,      0.028610, &
           0.036227,      0.045146,      0.055474,      0.067316,      0.080777, &
           0.095964,      0.112979,      0.131935,      0.152934,      0.176091, &
           0.201520,      0.229315,      0.259554,      0.291993,      0.326329, &
           0.362203,      0.399205,      0.436906,      0.475016,      0.513280, &
           0.551458,      0.589317,      0.626559,      0.662934,      0.698224, &
           0.732224,      0.764679,      0.795385,      0.824185,      0.850950, &
           0.875518,      0.897767,      0.917651,      0.935157,      0.950274, &
           0.963007,      0.973466,      0.982238,      0.989153,      0.994204, &
           0.997630,      1.000000 /)


  ! *** nc28 **************************************************************

  ! Hybrid coordinate coefficients at half level interfaces,
  ! specifying the 28 sigma levels of the NCEP reanalysis (surf -> top).
  ! Coefficient a is given in Pa, b in [0,1] .

  real, parameter  :: a_nc28(0:28) = 0.0

  real, parameter  :: b_nc28(0:28) = (/   &
              1.0000  ,  0.9900  ,   0.9742  ,   0.9546 ,  &
              0.9304  ,  0.9014  ,   0.8662  ,   0.8254 ,  &
              0.7774  ,  0.7242  ,   0.6644  ,   0.6014 ,  &
              0.5348  ,  0.4686  ,   0.4028  ,   0.3412 ,  &
              0.2838  ,  0.2326  ,   0.1876  ,   0.1488 ,  &
              0.1164  ,  0.0892  ,   0.0672  ,   0.0488 ,  &
              0.0348  ,  0.0228  ,   0.0138  ,   0.0064 ,  &
              0.0000 /)

  ! *** cmam **************************************************************

  real, parameter  :: a_msc71(0:71) = (/   &
  0.0000000E+00, 0.9697381E-01, 0.1359287E+00, 0.1903444E+00, 0.2661948E+00, 0.3729337E+00, &
  0.5219892E+00, 0.7302384E+00, 0.1023965E+01, 0.1429865E+01, 0.1996367E+01, 0.2795496E+01, &
  0.3910756E+01, 0.5469573E+01, 0.7652655E+01, 0.1064398E+02, 0.1476761E+02, 0.2035368E+02, &
  0.2782417E+02, 0.3784165E+02, 0.5115196E+02, 0.6867558E+02, 0.9154660E+02, 0.1215257E+03, &
  0.1601752E+03, 0.2096810E+03, 0.2730070E+03, 0.3534054E+03, 0.4539333E+03, 0.5791536E+03, &
  0.7337209E+03, 0.9020061E+03, 0.1151447E+04, 0.1438797E+04, 0.1756787E+04, 0.2144013E+04, &
  0.2592025E+04, 0.3131841E+04, 0.3718082E+04, 0.4406751E+04, 0.5184419E+04, 0.6002543E+04, &
  0.6901672E+04, 0.7882596E+04, 0.8892869E+04, 0.9923001E+04, 0.1096051E+05, 0.1196177E+05, &
  0.1289793E+05, 0.1371837E+05, 0.1436708E+05, 0.1479985E+05, 0.1496916E+05, 0.1481930E+05, &
  0.1430370E+05, 0.1346313E+05, 0.1242185E+05, 0.1128402E+05, 0.1013239E+05, 0.9018101E+04, &
  0.7976487E+04, 0.7027693E+04, 0.6182180E+04, 0.5440920E+04, 0.4799065E+04, 0.4202355E+04, &
  0.3590442E+04, 0.2963440E+04, 0.2321461E+04, 0.1664614E+04, 0.1115094E+04, 0.0000000E+00  /)

  real, parameter  :: b_msc71(0:71) = (/    &
  0.0000000E+00, 0.2585571E-09, 0.7040398E-09, 0.1535705E-08, 0.3013138E-08, 0.5588928E-08, &
  0.9977033E-08, 0.1738750E-07, 0.2995330E-07, 0.5067781E-07, 0.8520696E-07, 0.1431463E-06, &
  0.2392853E-06, 0.3989986E-06, 0.6646806E-06, 0.1095722E-05, 0.1800122E-05, 0.2924593E-05, &
  0.4696273E-05, 0.7484670E-05, 0.1182433E-04, 0.1850007E-04, 0.2865571E-04, 0.4416054E-04, &
  0.6735843E-04, 0.1018461E-03, 0.1529120E-03, 0.2279373E-03, 0.3362286E-03, 0.4919699E-03, &
  0.7133741E-03, 0.9869111E-03, 0.1466181E-02, 0.2084517E-02, 0.2893929E-02, 0.4006976E-02, &
  0.5507053E-02, 0.7581515E-02, 0.1018475E-01, 0.1375098E-01, 0.1841276E-01, 0.2415572E-01, &
  0.3156661E-01, 0.4113112E-01, 0.5287339E-01, 0.6738057E-01, 0.8526936E-01, 0.1069703E+00, &
  0.1332617E+00, 0.1646431E+00, 0.2011737E+00, 0.2437835E+00, 0.2934350E+00, 0.3496911E+00, &
  0.4130113E+00, 0.4790453E+00, 0.5416321E+00, 0.5992498E+00, 0.6510819E+00, 0.6971171E+00, &
  0.7375002E+00, 0.7725258E+00, 0.8025842E+00, 0.8281592E+00, 0.8497921E+00, 0.8694991E+00, &
  0.8893561E+00, 0.9093621E+00, 0.9295158E+00, 0.9498163E+00, 0.9665902E+00, 0.1000000E+01  /)



  ! --- interface -----------------------

  interface Init
    module procedure levi_Init
    module procedure levi_Init_levi
    module procedure levi_Init_key
    module procedure levi_Init_select
  end interface

  interface Done
    module procedure levi_Done
  end interface

  interface Check
    module procedure levi_Check
    module procedure levi_Check_3d
  end interface

  interface Compare
    module procedure levi_Compare
  end interface

  interface HPressure
    module procedure levi_HPressure_1d
    module procedure levi_HPressure_2d
    module procedure levi_HPressure_3d
  end interface

  interface FPressure
    module procedure levi_FPressure_1d
    module procedure levi_FPressure_2d
    module procedure levi_FPressure_3d
  end interface

  interface FillLevels
    module procedure levi_FillLevels_1d
    module procedure levi_FillLevels_2d
    module procedure levi_FillLevels_3d
  end interface


contains


  ! =========================================================


  subroutine levi_Init( levi, nlev, a, b, status, name, revert )

    use Binas, only : p0, grav

    ! --- in/out ---------------------------------------

    type(TLevelInfo), intent(out)       ::  levi
    integer, intent(in)                 ::  nlev
    real, intent(in)                    ::  a(:)
    real, intent(in)                    ::  b(:)
    integer, intent(out)                ::  status

    character(len=*), intent(in), optional  ::  name
    logical, intent(in), optional           ::  revert

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Init'

    ! --- local --------------------------------------

    integer             ::  l
    logical             ::  do_revert

    ! --- begin ---------------------------------------

    ! check param
    if ( (size(a) /= nlev+1) .or. (size(b) /= nlev+1) ) then
      write (*,'("ERROR - invalid length of half level coeff:")')
      write (*,'("ERROR -   nlev     : ",i4)') nlev
      write (*,'("ERROR -   size(a)  : ",i4," (",i4,")")') size(a), nlev+1
      write (*,'("ERROR -   size(b)  : ",i4," (",i4,")")') size(b), nlev+1
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! revert ?
    do_revert = .false.
    if ( present(revert) ) do_revert = revert
    ! store name ?
    if ( present(name) ) then
      levi%name = name
    else
      levi%name = 'levs'
    end if

    ! number of levels
    levi%nlev = nlev

    ! storage half level hybride coeff
    allocate( levi%a(0:nlev) )
    allocate( levi%b(0:nlev) )
    ! store:
    if ( do_revert ) then
      do l = 0, nlev
        levi%a(l) = a(nlev-l)
        levi%b(l) = b(nlev-l)
      end do
    else
      levi%a = a
      levi%b = b
    end if

    ! hybride coeff increments:
    allocate( levi%da(nlev) )
    allocate( levi%db(nlev) )
    levi%da = levi%a(1:nlev) - levi%a(0:nlev-1)
    levi%db = levi%b(1:nlev) - levi%b(0:nlev-1)

    ! full level hybride coeff:
    allocate( levi%fa(nlev) )
    allocate( levi%fb(nlev) )
    levi%fa = ( levi%a(0:nlev-1) + levi%a(1:nlev) )/2.0
    levi%fb = ( levi%b(0:nlev-1) + levi%b(1:nlev) )/2.0

    ! half level coeff per full level:
    allocate( levi%fa_bounds(1:2,1:nlev) )
    allocate( levi%fb_bounds(1:2,1:nlev) )
    levi%fa_bounds(1,:) = levi%a(0:nlev-1)
    levi%fa_bounds(2,:) = levi%a(1:nlev  )
    levi%fb_bounds(1,:) = levi%b(0:nlev-1)
    levi%fb_bounds(2,:) = levi%b(1:nlev  )

    ! fill standard pressures:
    allocate( levi%p0(0:nlev) )
    levi%p0(0:nlev) = levi%a(0:nlev) + levi%b(0:nlev) * p0     ! Pa
    allocate( levi%fp0(1:nlev) )
    levi%fp0 = levi%fa + levi%fb * p0     ! Pa
    allocate( levi%m0(1:nlev) )
    levi%m0 = abs( levi%da + levi%db * p0 )/grav     ! kg air / m2

    ! upwards (decreasing pressure) or downwards (increasing pressure)
    if ( levi%p0(0) > levi%p0(nlev) ) then
      levi%updo = 'u'
    else
      levi%updo = 'd'
    end if

    ! no selection info
    levi%selection = .false.
    nullify( levi%hlevs )
    nullify( levi%flevs )
    nullify( levi%a_parent )
    nullify( levi%b_parent )

    ! ok
    status = 0

  end subroutine levi_Init


  ! ***


  subroutine levi_Init_levi( levi, levi2, status )

    ! --- in/out ---------------------------------------

    type(TLevelInfo), intent(out)       ::  levi
    type(TLevelInfo), intent(in)        ::  levi2
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Init_levi'

    ! --- local --------------------------------------

    integer             ::  l

    ! --- begin ---------------------------------------

    ! copy fields:
    call Init( levi, levi2%nlev, levi2%a, levi2%b, status, name=trim(levi2%name) )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0

  end subroutine levi_Init_levi


  ! ***


  subroutine levi_Init_key( levi, key, status )

    ! --- in/out ---------------------------------------

    type(TLevelInfo), intent(out)       ::  levi
    character(len=*), intent(in)        ::  key
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Init_key'

    ! --- local ---------------------------------------

    real, allocatable   ::  a(:), b(:)
    integer             ::  l

    ! --- begin ---------------------------------------

    select case ( key )

      case ( 'ec19' )
        call Init( levi, 19, a_ec19, b_ec19, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case ( 'ec31' )
        call Init( levi, 31, a_ec31, b_ec31, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case ( 'ec40' )
        call Init( levi, 40, a_ec40, b_ec40, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case ( 'ec60' )
        call Init( levi, 60, a_ec60, b_ec60, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case ( 'ec62' )
        call Init( levi, 62, a_ec62, b_ec62, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case ( 'ec91' )
        call Init( levi, 91, a_ec91, b_ec91, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case ( 'tm31' )
        allocate( a(0:31) )
        allocate( b(0:31) )
        do l = 0, 31
          a(l) = a_ec31(31-l)
          b(l) = b_ec31(31-l)
        end do
        call Init( levi, 31, a, b, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        deallocate( a )
        deallocate( b )

      case ( 'tm40' )
        allocate( a(0:40) )
        allocate( b(0:40) )
        do l = 0, 40
          a(l) = a_ec40(40-l)
          b(l) = b_ec40(40-l)
        end do
        call Init( levi, 40, a, b, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        deallocate( a )
        deallocate( b )

      case ( 'tm60' )
        allocate( a(0:60) )
        allocate( b(0:60) )
        do l = 0, 60
          a(l) = a_ec60(60-l)
          b(l) = b_ec60(60-l)
        end do
        call Init( levi, 60, a, b, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        deallocate( a )
        deallocate( b )

      case ( 'tm62' )
        allocate( a(0:62) )
        allocate( b(0:62) )
        do l = 0, 62
          a(l) = a_ec62(62-l)
          b(l) = b_ec62(62-l)
        end do
        call Init( levi, 62, a, b, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        deallocate( a )
        deallocate( b )

      case ( 'tm91' )
        allocate( a(0:91) )
        allocate( b(0:91) )
        do l = 0, 91
          a(l) = a_ec91(91-l)
          b(l) = b_ec91(91-l)
        end do
        call Init( levi, 91, a, b, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        deallocate( a )
        deallocate( b )

      case ( 'nc28' )
        call Init( levi, 28, a_nc28, b_nc28, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case ( 'msc71' )
        call Init( levi, 71, a_msc71, b_msc71, status, name=trim(key) )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      case default
        write (*,'("ERROR - unknown key `",a,"`")') trim(key)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select

    ! ok
    status = 0

  end subroutine levi_Init_key


  ! ***


  subroutine levi_Init_select( levi, levi_parent, hlev_select, status, name )

    ! --- in/out ---------------------------------------

    type(TLevelInfo), intent(out)       ::  levi
    type(TLevelInfo), intent(in)        ::  levi_parent
    integer, intent(in)                 ::  hlev_select(:)
    integer, intent(out)                ::  status

    character(len=*), intent(in), optional  ::  name

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Init_select'

    ! --- local --------------------------------------

    integer             ::  nlev, l

    ! --- begin ---------------------------------------

    ! check length of array:
    if ( (size(hlev_select) < 2) .or. (size(hlev_select) > levi_parent%nlev+1) ) then
      write (*,'("ERROR - strange length of array with selected levels:")')
      write (*,'("ERROR -   expected : ",i4,", ..,",i4)') 2, levi_parent%nlev+1
      write (*,'("ERROR -   found    : ",i4)') size(hlev_select)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! check range of values:
    if ( (minval(hlev_select) /= 0) .or. (maxval(hlev_select) /= levi_parent%nlev) ) then
      write (*,'("ERROR - invalid range of selected levels:")')
      write (*,'("ERROR -   expected : ",i4,", ..,",i4)') 0, levi_parent%nlev
      write (*,'("ERROR -   found    : ",i4,", ..,",i4)') minval(hlev_select), maxval(hlev_select)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! set number of full levels
    nlev = size(hlev_select) - 1

    ! copy coeff
    call Init( levi, nlev, levi_parent%a(hlev_select), levi_parent%b(hlev_select), status, name=name )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! store selection info:
    levi%selection = .true.
    ! structure will contain selection info:
    !  o levi_parent half levels corresponding to levi half levels
    allocate( levi%hlevs(0:nlev) )
    levi%hlevs = hlev_select
    !  o range of levi_parent full levels covered by a levi level
    allocate( levi%flevs(nlev,2) )
    do l = 1, nlev
      levi%flevs(l,1) = hlev_select(l+1) + 1
      levi%flevs(l,2) = hlev_select(l)
    end do
    !  o original hybride coeff:
    levi%nlev_parent = levi_parent%nlev
    allocate( levi%a_parent(0:levi%nlev_parent) )
    allocate( levi%b_parent(0:levi%nlev_parent) )
    levi%a_parent = levi_parent%a
    levi%b_parent = levi_parent%b

    ! ok
    status = 0

  end subroutine levi_Init_select


  ! ***


  subroutine levi_Done( levi, status )

    ! --- in/out ---------------------------------------

    type(TLevelInfo), intent(inout)     ::  levi
    integer, intent(out)                ::  status

    ! --- begin ---------------------------------------

    ! deallocate storage for hybride coeff
    deallocate( levi%a )
    deallocate( levi%b )

    ! deallocate storage for full level hybride coeff
    deallocate( levi%fa )
    deallocate( levi%fb )

    ! deallocate storage for half level coeff per full level:
    deallocate( levi%fa_bounds )
    deallocate( levi%fb_bounds )

    ! deallocate storage for standard pressures and mass:
    deallocate( levi%p0  )
    deallocate( levi%fp0 )
    deallocate( levi%m0  )

    ! deallocate storage for increment coeff
    deallocate( levi%da )
    deallocate( levi%db )

    ! selection info ?
    if ( levi%selection ) then
      deallocate( levi%hlevs )
      deallocate( levi%flevs )
      deallocate( levi%a_parent )
      deallocate( levi%b_parent )
    end if

    ! ok
    status = 0

  end subroutine levi_Done


  ! ***


  subroutine levi_Check( levi, status )

    ! --- in/out ------------------------------

    type(TLevelInfo), intent(in)     ::  levi
    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Check'

    ! --- begin -------------------------------

    if ( levi%nlev < 1 ) then
      write (*,'("ERROR - level info contains strange number of levels:")')
      write (*,'("ERROR -   levi%nlev  : ",i4)') levi%nlev
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    if ( (.not. associated(levi%a)) .or. (.not. associated(levi%b)) ) then
      write (*,'("ERROR - hybride coeffs in level info not allocated properly.")')
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    if ( (size(levi%a)/=levi%nlev+1) .or. (size(levi%b)/=levi%nlev+1) ) then
      write (*,'("ERROR - hybride coeffs in level info wrong size:")')
      write (*,'("ERROR -   nlev   : ",i4)') levi%nlev
      write (*,'("ERROR -   size a : ",i4)') size(levi%a)
      write (*,'("ERROR -   size b : ",i4)') size(levi%b)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    if ( (.not. associated(levi%fa)) .or. (.not. associated(levi%fb)) ) then
      write (*,'("ERROR - f hybride coeffs in level info not allocated properly.")')
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    if ( (size(levi%fa)/=levi%nlev) .or. (size(levi%fb)/=levi%nlev) ) then
      write (*,'("ERROR - f hybride coeffs in level info wrong size:")')
      write (*,'("ERROR -   nlev    : ",i4)') levi%nlev
      write (*,'("ERROR -   size fa : ",i4)') size(levi%fa)
      write (*,'("ERROR -   size fb : ",i4)') size(levi%fb)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine levi_Check


  ! ***


  subroutine levi_Check_3d( levi, ll, status )

    ! --- in/out ------------------------------

    type(TLevelInfo), intent(in)     ::  levi
    real, intent(in)                 ::  ll(:,:,:)
    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Check_3d'

    ! --- begin -------------------------------

    call Check( levi, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    if ( size(ll,3) /= levi%nlev ) then
      write (*,'("ERROR - size of 3D grid does not match with level info")')
      write (*,'("ERROR -   size(ll,3) : ",i4)') size(ll,3)
      write (*,'("ERROR -   levi%nlev  : ",i4)') levi%nlev
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine levi_Check_3d



  ! ##############################################################
  ! ###
  ! ### combining levels
  ! ###
  ! ##############################################################

  !
  ! How could levi be formed from leviX ?
  !
  !  o levi could be combined out of leviX levels (thus leviX is a superset)
  !  o levels (superset?) in leviX might be a reversed version of levi levels
  !


  subroutine levi_Compare( levi, leviX, verbose, fillkey, reverse, status )

    ! --- in/out ----------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    type(TLevelInfo), intent(in)         ::  leviX
    integer, intent(in)                  ::  verbose
    character(len=10), intent(out)       ::  fillkey
    logical, intent(out)                 ::  reverse
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Compare'

    ! --- local -----------------------------------

    integer                ::  l
    real, allocatable      ::  Xa_rev(:), Xb_rev(:)

    ! --- begin ------------------------------------

    ! reverse hybride coeff
    allocate( Xa_rev(0:leviX%nlev), Xb_rev(0:leviX%nlev) )
    do l = 0, leviX%nlev
      Xa_rev(l) = leviX%a(leviX%nlev-l)
      Xb_rev(l) = leviX%b(leviX%nlev-l)
    end do

    ! by default ...
    fillkey = 'interpol'
    reverse  = levi%updo /= leviX%updo

    ! same number of levels ?
    if ( leviX%nlev == levi%nlev ) then

      ! exact match or reverse ...
      if      ( all(abs(leviX%a-levi%a)<0.1) .and. all(abs(leviX%b-levi%b)<0.1) ) then
        fillkey = 'copy'
        reverse = .false.
      else if ( all(abs(Xa_rev -levi%a)<0.1) .and. all(abs(Xb_rev -levi%b)<0.1) ) then
        fillkey = 'copy'
        reverse = .true.
      end if

    ! selection ?
    else if ( levi%selection .and. (leviX%nlev == levi%nlev_parent) ) then

      ! exact match or reverse ...
      if      ( all(abs(leviX%a-levi%a_parent)<0.1) .and. all(abs(leviX%b-levi%b_parent)<0.1) ) then
        fillkey = 'combine'
        reverse = .false.
      else if ( all(abs( Xa_rev-levi%a_parent)<0.1) .and. all(abs( Xb_rev-levi%b_parent)<0.1) ) then
        fillkey = 'combine'
        reverse  = .true.
      end if

    end if

    ! info on matching ?
    !if ( verbose > 0 ) then
    !  write (*,'("DEBUG - try to match levels:")')
    !  write (*,'("DEBUG -   requested :")')
    !  write (*,'("DEBUG -        l        a              b        ")')
    !  write (*,'("DEBUG -     ----  --------------  --------------")')
    !  do l = 0, levi%nlev
    !    write (*,'("DEBUG -     ",i4,2f16.8)') l, levi%a(l), levi%b(l)
    !  end do
    !  if ( levi%selection ) then
    !    write (*,'("DEBUG -   parent :")')
    !    write (*,'("DEBUG -        l        a              b        ")')
    !    write (*,'("DEBUG -     ----  --------------  --------------")')
    !    do l = 0, levi%nlev_parent
    !      write (*,'("DEBUG -     ",i4,2f16.8)') l, levi%a_parent(l), levi%b_parent(l)
    !    end do
    !  end if
    !  write (*,'("DEBUG -   found :")')
    !  write (*,'("DEBUG -        l        a              b        ")')
    !  write (*,'("DEBUG -     ----  --------------  --------------")')
    !  do l = 0, leviX%nlev
    !    write (*,'("DEBUG -     ",i4,2f16.8)') l, leviX%a(l), leviX%b(l)
    !  end do
    !  write (*,'("DEBUG -   fillkey : ",a)') trim(fillkey)
    !  write (*,'("DEBUG -   reverse : ",l2)') reverse
    !  write (*,'("DEBUG in ",a)') rname; status = 1 ; return
    !end if

    ! done
    deallocate( Xa_rev, Xb_rev )

    ! ok
    status = 0

  end subroutine levi_Compare


  ! =====================================================================
  ! ===
  ! === pressure levels
  ! ===
  ! =====================================================================

  subroutine levi_HPressure_1d( levi, sp, ph, status )

    ! --- in/out -------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp        ! Pa
    real, intent(out)                    ::  ph(:)     ! Pa
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_HPressure_1d'

    ! --- begin ---------------------------------

    ! check ...
    if ( size(ph) /= levi%nlev+1 ) then
      write (*,'("ERROR - shape of output grid does not match definition:")')
      write (*,'("ERROR -   half levels  : ",i4)') levi%nlev+1
      write (*,'("ERROR -   ph shape     : ",i4)') shape(ph)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! half level pressure
    ph = levi%a + levi%b * sp

    ! ok
    status = 0

  end subroutine levi_HPressure_1d


  ! *


  subroutine levi_HPressure_2d( levi, sp, ph, status )

    ! --- in/out -------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp(:)       ! Pa
    real, intent(out)                    ::  ph(:,:)     ! Pa
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_HPressure_2d'

    ! --- local ---------------------------------

    integer     ::  i

    ! --- begin ---------------------------------

    ! check ...
    if ( any( shape(ph) /= (/shape(sp),levi%nlev+1/) ) ) then
      write (*,'("ERROR - shape of output grid does not match definition:")')
      write (*,'("ERROR -   sp shape     : ",i6)') shape(sp)
      write (*,'("ERROR -   half levels  : ",i4)') levi%nlev+1
      write (*,'("ERROR -   ph shape     : ",i6,i4)') shape(ph)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! half level pressure
    do i = 1, size(sp)
      ph(i,:) = levi%a + levi%b * sp(i)
    end do

    ! ok
    status = 0

  end subroutine levi_HPressure_2d


  ! *


  subroutine levi_HPressure_3d( levi, sp, ph, status )

    ! --- in/out -------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp(:,:)       ! Pa
    real, intent(out)                    ::  ph(:,:,:)     ! Pa
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_HPressure_3d'

    ! --- local ---------------------------------

    integer     ::  i, j

    ! --- begin ---------------------------------

    ! check ...
    if ( any( shape(ph) /= (/shape(sp),levi%nlev+1/) ) ) then
      write (*,'("ERROR - shape of output grid does not match definition:")')
      write (*,'("ERROR -   sp shape     : ",2i6)') shape(sp)
      write (*,'("ERROR -   half levels  : ",i4)') levi%nlev+1
      write (*,'("ERROR -   ph shape     : ",2i6,i4)') shape(ph)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! half level pressure
    do i = 1, size(sp,1)
      do j = 1, size(sp,2)
        ph(i,j,:) = levi%a + levi%b * sp(i,j)
      end do
    end do

    ! ok
    status = 0

  end subroutine levi_HPressure_3d




  ! ***


  subroutine levi_FPressure_1d( levi, sp, pf, status )

    ! --- in/out -------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp        ! Pa
    real, intent(out)                    ::  pf(:)     ! Pa
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_FPressure_1d'

    ! --- begin ---------------------------------

    ! check ...
    if ( size(pf) /= levi%nlev ) then
      write (*,'("ERROR - shape of output grid does not match definition:")')
      write (*,'("ERROR -   full levels  : ",i4)') levi%nlev
      write (*,'("ERROR -   pf shape     : ",i4)') shape(pf)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! full level pressure
    pf = levi%fa + levi%fb * sp

    ! ok
    status = 0

  end subroutine levi_FPressure_1d


  ! *


  subroutine levi_FPressure_2d( levi, sp, pf, status )

    ! --- in/out -------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp(:)       ! Pa
    real, intent(out)                    ::  pf(:,:)     ! Pa
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_FPressure_2d'

    ! --- local ---------------------------------

    integer     ::  i

    ! --- begin ---------------------------------

    ! check ...
    if ( any( shape(pf) /= (/shape(sp),levi%nlev/) ) ) then
      write (*,'("ERROR - shape of output grid does not match definition:")')
      write (*,'("ERROR -   sp shape     : ",i6)') shape(sp)
      write (*,'("ERROR -   full levels  : ",i4)') levi%nlev
      write (*,'("ERROR -   pf shape     : ",i6,i4)') shape(pf)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! full level pressure
    do i = 1, size(sp)
      pf(i,:) = levi%fa + levi%fb * sp(i)
    end do

    ! ok
    status = 0

  end subroutine levi_FPressure_2d


  ! *


  subroutine levi_FPressure_3d( levi, sp, pf, status )

    ! --- in/out -------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp(:,:)       ! Pa
    real, intent(out)                    ::  pf(:,:,:)     ! Pa
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_FPressure_3d'

    ! --- local ---------------------------------

    integer     ::  i, j

    ! --- begin ---------------------------------

    ! check ...
    if ( any( shape(pf) /= (/shape(sp),levi%nlev/) ) ) then
      write (*,'("ERROR - shape of output grid does not match definition:")')
      write (*,'("ERROR -   sp shape     : ",2i6)') shape(sp)
      write (*,'("ERROR -   full levels  : ",i4)') levi%nlev
      write (*,'("ERROR -   pf shape     : ",2i6,i4)') shape(pf)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! full level pressure
    do i = 1, size(sp,1)
      do j = 1, size(sp,2)
        pf(i,j,:) = levi%fa + levi%fb * sp(i,j)
      end do
    end do

    ! ok
    status = 0

  end subroutine levi_FPressure_3d



  ! =====================================================================
  ! ===
  ! === 3D conversion
  ! ===
  ! =====================================================================


  !
  ! Copy llX (defined by leviX) into ll (defined by levi),
  ! eventually by combining levels.
  ! How levels are combined is specified by a key:
  !  'bottom'     :  use bottom value (most close to the ground)
  !  'top'        :  use top value (most close to the model top)
  !  'sum'        :  sum values
  !  'aver'       :  average of all levels
  !  'mass-aver'  :  mass weighted average
  ! Eventually reverse.
  ! Surface pressure field is used to compute massed in case of 'mass-aver' combination.
  !

  subroutine levi_FillLevels_1d( levi, ps, ll, leviX, llX, combine_key, status )

    ! --- in/out ----------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  ps        ! Pa
    real, intent(out)                    ::  ll(:)
    type(TLevelInfo), intent(in)         ::  leviX
    real, intent(in)                     ::  llX(:)
    character(len=*), intent(in)         ::  combine_key
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_FillLevels_1d'

    ! --- local ------------------------------

    integer            ::  verbose
    character(len=10)  ::  fillkey
    logical            ::  reverse

    integer            ::  k, l
    integer            ::  le, nle, le1, le2
    real               ::  Xfp0, Xfp0_save

    ! --- begin -------------------------------

    ! for safety ...
    status = 1

    ! correct number of levels ?
    if ( size(ll) /= levi%nlev ) then
      write (*,'("ERROR - number of levels in output grid does not match definition:")')
      write (*,'("ERROR -   ll  levels  : ",i3)') size(ll)
      write (*,'("ERROR -   levi nlev   : ",i3)') levi%nlev
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! necessary to combine or reverse ?
    verbose = 1
    call Compare( levi, leviX, verbose, fillkey, reverse, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! fill given key:
    select case ( fillkey )

      !
      ! ** copy levels
      !

      case ( 'copy' )

        ! loop over all TM levels
        do l = 1, levi%nlev
          ! index of corresponding level in field 'llX' read from file:
          k = l
          if ( reverse ) k = levi%nlev + 1 - l
          ! copy level:
          ll(l) = llX(k)
        end do

      !
      ! ** combine levels
      !

      case ( 'combine' )

        ! loop over all TM levels
        ll = 0.0
        do l = 1, levi%nlev
          ! loop over range of parent levels:
          le1 = levi%flevs(l,1)
          le2 = levi%flevs(l,2)
          nle = le2 - le1 + 1
          Xfp0_save = -1.0
          do le = le1, le2
            ! index of corresponding level in field 'llX' read from file:
            k = le
            if ( reverse ) k = levi%nlev_parent + 1 - le
            ! standard pressure for level in llX :
            Xfp0 = leviX%fp0(k)
            ! based on combine key, add contribution of this level:
            select case ( combine_key )
              case ( 'sum' )
                ll(l) = ll(l) + llX(k)
              case ( 'aver' )
                ll(l) = ll(l) + llX(k)/(nle*1.0)
              case ( 'mass-aver' )
                !
                !  m = dp * A / g = |da+db*ps| * A/g
                !  ( X1*mX1 + X2*mX2 + .. ) / m
                !  ~  ( X1*|daX1+dbX1*ps| + X2*|daX2+dbX2*ps| + .. ) / |da+db*ps|
                !
                ll(l) = ll(l) + llX(k)*abs(leviX%da(k)+leviX%db(k)*ps) / &
                                       abs( levi%da(l)+ levi%db(l)*ps)
              case ( 'bottom' )
                if ( (Xfp0_save < 0.0) .or. (Xfp0 >= Xfp0_save) ) then
                  ll(l) = llX(k)
                  Xfp0_save = Xfp0
                end if
              case ( 'top' )
                if ( (Xfp0_save < 0.0) .or. (Xfp0 <= Xfp0_save) ) then
                  ll(l) = llX(k)
                  Xfp0_save = Xfp0
                end if
              case default
                write (*,'("ERROR - combine key `",a,"` not supported.")') trim(combine_key)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select
          end do
        end do

    end select

    ! ok
    status = 0

  end subroutine levi_FillLevels_1d


  ! ***


!  subroutine levi_FillLevels_2d( levi, ps, ll, leviX, llX, combine_key, status )
!
!    ! --- in/out ----------------------------------
!
!    type(TLevelInfo), intent(in)         ::  levi
!    real, intent(in)                     ::  ps(:)        ! Pa
!    real, intent(out)                    ::  ll(:,:)
!    type(TLevelInfo), intent(in)         ::  leviX
!    real, intent(in)                     ::  llX(:,:)
!    character(len=*), intent(in)         ::  combine_key
!    integer, intent(out)                 ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/levi_FillLevels_2d'
!
!    ! --- local ------------------------------
!
!    integer            ::  verbose
!    character(len=10)  ::  fillkey
!    logical            ::  reverse
!
!    integer            ::  k, l
!    integer            ::  le, nle, le1, le2
!    real               ::  Xfp0, Xfp0_save
!
!    ! --- begin -------------------------------
!
!    ! for safety ...
!    status = 1
!
!    ! correct number of levels ?
!    if ( size(ll,2) /= levi%nlev ) then
!      write (*,'("ERROR - number of levels in output grid does not match definition:")')
!      write (*,'("ERROR -   ll  levels  : ",i3)') size(ll,2)
!      write (*,'("ERROR -   levi nlev   : ",i3)') levi%nlev
!      write (*,'("ERROR in ",a)') rname; status=1; return
!    end if
!
!    ! same grid sizes of ll and llX ?
!    if ( (size(ll,1) /= size(llX,1)) ) then
!      write (*,'("ERROR - horizontal sizes do not match:")')
!      write (*,'("ERROR -   ll     : ",i3)') size(ll ,1)
!      write (*,'("ERROR -   llX    : ",i3)') size(llX,1)
!      write (*,'("ERROR in ",a)') rname; status=1; return
!    end if
!
!    ! correct size of surface pressure field ?
!    if ( combine_key == 'mass-aver' ) then
!      if ( (size(ps,1) /= size(ll,1)) ) then
!        write (*,'("ERROR - horizontal sizes do not match:")')
!        write (*,'("ERROR -   ps     : ",i3)') size(ps ,1)
!        write (*,'("ERROR -   ll     : ",i3)') size(ll ,1)
!        write (*,'("ERROR in ",a)') rname; status=1; return
!      end if
!    end if
!
!    ! necessary to combine or reverse ?
!    verbose = 1
!    call Compare( levi, leviX, verbose, fillkey, reverse, status )
!    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
!
!    ! fill given key:
!    select case ( fillkey )
!
!      !
!      ! ** copy levels
!      !
!
!      case ( 'copy' )
!
!        ! loop over all TM levels
!        do l = 1, levi%nlev
!          ! index of corresponding level in field 'llX' read from file:
!          k = l
!          if ( reverse ) k = levi%nlev + 1 - l
!          ! copy level:
!          ll(:,l) = llX(:,k)
!        end do
!
!      !
!      ! ** combine levels
!      !
!
!      case ( 'combine' )
!
!        ! loop over all TM levels
!        ll = 0.0
!        do l = 1, levi%nlev
!          ! loop over range of parent levels:
!          le1 = levi%flevs(l,1)
!          le2 = levi%flevs(l,2)
!          nle = le2 - le1 + 1
!          Xfp0_save = -1.0
!          do le = le1, le2
!            ! index of corresponding level in field 'llX' read from file:
!            k = le
!            if ( reverse ) k = levi%nlev_parent + 1 - le
!            ! standard pressure for level in llX :
!            Xfp0 = leviX%fp0(k)
!            ! based on combine key, add contribution of this level:
!            select case ( combine_key )
!              case ( 'sum' )
!                ll(:,l) = ll(:,l) + llX(:,k)
!              case ( 'aver' )
!                ll(:,l) = ll(:,l) + llX(:,k)/(nle*1.0)
!              case ( 'mass-aver' )
!                !
!                !  m = dp * A / g = |da+db*ps| * A/g
!                !
!                !  ( X1*mX1 + X2*mX2 + .. ) / m
!                !    =  ( X1*|daX1+dbX1*ps| + X2*|daX2+dbX2*ps| + .. ) / |da+db*ps|
!                !
!                ll(:,l) = ll(:,l) + llX(:,k)*abs(leviX%da(k)+leviX%db(k)*ps) / &
!                                             abs( levi%da(l)+ levi%db(l)*ps)
!              case ( 'bottom' )
!                if ( (Xfp0_save < 0.0) .or. (Xfp0 >= Xfp0_save) ) then
!                  ll(:,l) = llX(:,k)
!                  Xfp0_save = Xfp0
!                end if
!              case ( 'top' )
!                if ( (Xfp0_save < 0.0) .or. (Xfp0 <= Xfp0_save) ) then
!                  ll(:,l) = llX(:,k)
!                  Xfp0_save = Xfp0
!                end if
!              case default
!                write (*,'("ERROR - combine key `",a,"` not supported.")') trim(combine_key)
!                write (*,'("ERROR in ",a)') rname; status=1; return
!            end select
!          end do
!        end do
!
!      !
!      ! ** error ..
!      !
!
!      case default
!        write (*,'("ERROR - fill key `",a,"` not supported.")') trim(fillkey)
!        write (*,'("ERROR in ",a)') rname; status=1; return
!    end select
!
!    ! ok
!    status = 0
!
!  end subroutine levi_FillLevels_2d


  ! ***


  subroutine levi_FillLevels_2d( levi, nw, ps, ll, &
                                 leviX, llX, combine_key, status )

    use Num, only : IntervalSum, Interp_Lin

    ! --- in/out ----------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    character(len=*), intent(in)         ::  nw
    real, intent(in)                     ::  ps(:)        ! Pa
    real, intent(out)                    ::  ll(:,:)
    type(TLevelInfo), intent(in)         ::  leviX
    real, intent(in)                     ::  llX(:,:)
    character(len=*), intent(in)         ::  combine_key
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_FillLevels_2d'

    ! --- local ------------------------------

    integer            ::  verbose
    character(len=10)  ::  fillkey
    logical            ::  reverse

    integer            ::  j, k, l
    integer            ::  le, nle, le1, le2
    real               ::  Xfp0, Xfp0_save

    real, allocatable  ::  phalf(:)
    real, allocatable  ::  phalfX(:)
    real, allocatable  ::  dpX(:)
    real               ::  pfac
    integer            ::  ilast

    ! --- begin -------------------------------

    ! for safety ...
    status = 1

    ! correct number of levels ?
    if ( ((nw == 'n') .and. (size(ll,2) /= levi%nlev  )) .or. &
         ((nw == 'w') .and. (size(ll,2) /= levi%nlev+1)) ) then
      write (*,'("ERROR - number of levels in output grid does not match definition:")')
      write (*,'("ERROR -   levi nlev   : ",i3)') levi%nlev
      write (*,'("ERROR -   nw          : ",a )') nw
      write (*,'("ERROR -   ll  levels  : ",i3)') size(ll,2)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! correct number of levels ?
    if ( ((nw == 'n') .and. (size(llX,2) /= leviX%nlev  )) .or. &
         ((nw == 'w') .and. (size(llX,2) /= leviX%nlev+1)) ) then
      write (*,'("ERROR - number of levels in input grid does not match definition:")')
      write (*,'("ERROR -   leviX nlev   : ",i3)') leviX%nlev
      write (*,'("ERROR -   nw           : ",a )') nw
      write (*,'("ERROR -   llX  levels  : ",i3)') size(llX,2)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! same horizontal grid sizes of ll and llX ?
    if ( (size(ll,1) /= size(llX,1)) ) then
      write (*,'("ERROR - horizontal size does not match:")')
      write (*,'("ERROR -   ll     : ",i3)') size(ll,1)
      write (*,'("ERROR -   llX    : ",i3)') size(llX,1)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! necessary to combine or reverse ?
    verbose = 0
    call Compare( levi, leviX, verbose, fillkey, reverse, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! correct size of surface pressure field ?
    if ( (combine_key == 'mass-aver') .or. (fillkey == 'interpol') ) then
      if ( (size(ps) /= size(ll,1)) ) then
        write (*,'("ERROR - horizontal sizes do not match:")')
        write (*,'("ERROR -   ps          : ",i3)') size(ps)
        write (*,'("ERROR -   ll          : ",i3)') size(ll,1)
        write (*,'("ERROR -   combine_key : ",a)') combine_key
        write (*,'("ERROR -   fillkey     : ",a)') fillkey
        write (*,'("ERROR -   reverse     : ",l1)') reverse
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
    end if

    ! levels or half-levels ?
    select case ( nw )

      !
      ! === levels ===================================
      !

      case ( 'n' )

        ! fill given key:
        select case ( fillkey )

          !
          ! ** copy levels
          !

          case ( 'copy' )

            ! loop over all TM levels
            if ( reverse ) then
              do l = 1, levi%nlev
                ! copy corresponding level in field 'llX' read from file:
                ll(:,l) = llX(:,levi%nlev+1-l)
              end do
            else
              ll = llX
            end if

          !
          ! ** combine levels
          !

          case ( 'combine' )

            ! loop over all TM levels
            ll = 0.0
            do l = 1, levi%nlev
              ! loop over range of parent levels:
              le1 = levi%flevs(l,1)
              le2 = levi%flevs(l,2)
              nle = le2 - le1 + 1
              Xfp0_save = -1.0
              do le = le1, le2
                ! index of corresponding level in field 'llX' read from file:
                k = le
                if ( reverse ) k = levi%nlev_parent + 1 - le
                ! standard pressure for level in llX :
                Xfp0 = leviX%fp0(k)
                ! based on combine key, add contribution of this level:
                select case ( combine_key )
                  case ( 'sum' )
                    ll(:,l) = ll(:,l) + llX(:,k)
                  case ( 'aver' )
                    ll(:,l) = ll(:,l) + llX(:,k)/(nle*1.0)
                  case ( 'mass-aver' )
                    !
                    !  m = dp * A / g = |da+db*ps| * A/g
                    !  ( X1*mX1 + X2*mX2 + .. ) / m
                    !  ~  ( X1*|daX1+dbX1*ps| + X2*|daX2+dbX2*ps| + .. ) / |da+db*ps|
                    !
                    ll(:,l) = ll(:,l) + llX(:,k)*abs(leviX%da(k)+leviX%db(k)*ps) / &
                                                 abs( levi%da(l)+ levi%db(l)*ps)
                  case ( 'bottom' )
                    if ( (Xfp0_save < 0.0) .or. (Xfp0 >= Xfp0_save) ) then
                      ll(:,l) = llX(:,k)
                      Xfp0_save = Xfp0
                    end if
                  case ( 'top' )
                    if ( (Xfp0_save < 0.0) .or. (Xfp0 <= Xfp0_save) ) then
                      ll(:,l) = llX(:,k)
                      Xfp0_save = Xfp0
                    end if
                  case default
                    write (*,'("ERROR - combine key not supported:")')
                    write (*,'("ERROR -   combine key : ",a)') trim(combine_key)
                    write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
                    write (*,'("ERROR -   nw key      : ",a)') trim(nw)
                    write (*,'("ERROR in ",a)') rname; status=1; return
                end select
              end do
            end do

          !
          ! ** interpol
          !

          case ( 'interpol' )

            allocate( phalf (0: levi%nlev) )
            allocate( phalfX(0:leviX%nlev) )
            allocate(    dpX(  leviX%nlev) )

            ! phalfX should be increasing
            pfac = 1.0
            if ( leviX%updo == 'u' ) pfac = -1.0

            select case ( combine_key )

              case ( 'sum' )

                do j = 1, size(ll,1)

                  phalf  =  levi%a  +  levi%b  * ps(j)    ! Pa
                  phalfX = leviX%a  + leviX%b  * ps(j)    ! Pa

                  ilast = 1
                  do l = 1, levi%nlev
                    ! take partly sums of f(i)*m(i) within [phalf(l-1),phalf(l)]
                    ! NOTE: if phalf(l) < phalf(l-1), the result of IntervalSum is negative
                    call IntervalSum( phalfX*pfac, llX(j,:), &
                                      phalf(l-1)*pfac, phalf(l)*pfac, &
                                      ll(j,l), ilast, status )
                    if (status/=0) then;
                      write (*,'("ERROR - from IntervalSum (combine key `sum`)")')
                      write (*,'("ERROR -    j    : ",2i4   )')  j
                      write (*,'("ERROR -   ps      : ", f16.4)') ps(j)
                      write (*,'("ERROR -   leviX%a : ",es11.4)') leviX%a
                      write (*,'("ERROR -   leviX%b : ",es11.4)') leviX%b
                      write (*,'("ERROR in ",a)') rname; status=1; return
                    end if
                    ll(j,l) = ll(j,l) * sign(1.0,(phalf(l)-phalf(l-1))*pfac)
                  end do

                end do

              case ( 'mass-aver' )

                do j = 1, size(ll,1)

                  phalf  =  levi%a  +  levi%b  * ps(j)    ! Pa
                  phalfX = leviX%a  + leviX%b  * ps(j)    ! Pa
                  dpX = abs(leviX%da + leviX%db * ps(j))    ! Pa

                  ilast = 1
                  do l = 1, levi%nlev
                    ! take partly sums of f(i)*m(i) within [phalf(l-1),phalf(l)]
                    ! NOTE: if phalf(l) < phalf(l-1), the result of IntervalSum is negative
                    call IntervalSum( phalfX*pfac, llX(j,:)*dpX, &
                                      phalf(l-1)*pfac, phalf(l)*pfac, &
                                      ll(j,l), ilast, status )
                    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                    ll(j,l) = ll(j,l) / (phalf(l)-phalf(l-1))*pfac
                  end do

                end do

              case ( 'top' )

                do j = 1, size(ll,1)

                  phalf  =  levi%a  +  levi%b  * ps(j)    ! Pa, 0..
                  phalfX = leviX%a  + leviX%b  * ps(j)    ! Pa, 0..

                  ! loop over all levels:
                  do l = 1, levi%nlev
                    ! interpolate to top half level pressure:
                    call Interp_Lin( phalfX(1:leviX%nlev)*pfac, llX(j,:), &
                                     phalf (l)           *pfac,  ll(j,l), ilast, status )
                    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                  end do

                end do

              case ( 'bottom' )

                do j = 1, size(ll,1)

                  phalf  =  levi%a  +  levi%b  * ps(j)    ! Pa, 0..
                  phalfX = leviX%a  + leviX%b  * ps(j)    ! Pa, 0..

                  ! loop over all levels:
                  do l = 1, levi%nlev
                    ! interpolate to bottom half level pressure:
                    call Interp_Lin( phalfX(0:leviX%nlev-1)*pfac, llX(j,:), &
                                     phalf (l-1)           *pfac,  ll(j,l), ilast, status )
                    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                  end do

                end do

              case default
                write (*,'("ERROR - combine key not supported:")')
                write (*,'("ERROR -   combine key : ",a)') trim(combine_key)
                write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
                write (*,'("ERROR -   nw key      : ",a)') trim(nw)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

            deallocate( phalf  )
            deallocate( phalfX )
            deallocate(    dpX )

          !
          ! ** error ..
          !

          case default
            write (*,'("ERROR - fill key not supported:")')
            write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
            write (*,'("ERROR -   nw key      : ",a)') trim(nw)
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      !
      ! === half levels ========================================
      !

      case ( 'w' )

        ! fill given key:
        select case ( fillkey )

          !
          ! ** copy levels
          !

          case ( 'copy' )

            ! loop over all TM levels
            if ( reverse ) then
              do l = 1, levi%nlev+1
                ! copy corresponding level in field 'llX' read from file:
                ll(:,l) = llX(:,levi%nlev+2-l)
              end do
            else
              ll = llX
            end if

          !
          ! ** combine levels
          !

          case ( 'combine' )

            ll = 0.0

            ! loop over target half levels
            do l = 0, levi%nlev
              ! index of corresponding level in field 'llX' read from file:
              ! note: k=0,..,nlev
              k = levi%hlevs(l)
              if ( reverse ) k = levi%nlev_parent - k
              ! copy:
              ll(:,l+1) = llX(:,k+1)
            end do

          !
          ! ** interpol
          !

          case ( 'interpol' )

            allocate( phalf (0: levi%nlev) )
            allocate( phalfX(0:leviX%nlev) )
            allocate(    dpX(  leviX%nlev) )

            ! phalfX should be increasing
            pfac = 1.0
            if ( leviX%updo == 'u' ) pfac = -1.0

            select case ( combine_key )

              case ( 'top', 'bottom' )

                do j = 1, size(ll,1)

                  phalf  =  levi%a  +  levi%b  * ps(j)    ! Pa
                  phalfX = leviX%a  + leviX%b  * ps(j)    ! Pa

                  ! loop over all half levels:
                  do l = 0, levi%nlev
                    ! interpolate to half level pressure:
                    call Interp_Lin( phalfX   *pfac, llX(j,:  ), &
                                     phalf (l)*pfac,  ll(j,l+1), ilast, status )
                    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                  end do

                end do

              case default
                write (*,'("ERROR - combine key not supported:")')
                write (*,'("ERROR -   combine key : ",a)') trim(combine_key)
                write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
                write (*,'("ERROR -   nw key      : ",a)') trim(nw)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

            deallocate( phalf  )
            deallocate( phalfX )
            deallocate(    dpX )

          !
          ! ** error ..
          !

          case default
            write (*,'("ERROR - fill key not supported:")')
            write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
            write (*,'("ERROR -   nw key      : ",a)') trim(nw)
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      !
      ! === error ========================================
      !

      case default

        write (*,'("ERROR - nw key `",a,"` not supported.")') trim(nw)
        write (*,'("ERROR in ",a)') rname; status=1; return

    end select

    ! ok
    status = 0

  end subroutine levi_FillLevels_2d


  ! ***


  subroutine levi_FillLevels_3d( levi, nw, ps, ll, &
                                 leviX, llX, combine_key, status )

    use Num, only : IntervalSum, Interp_Lin

    ! --- in/out ----------------------------------

    type(TLevelInfo), intent(in)         ::  levi
    character(len=*), intent(in)         ::  nw
    real, intent(in)                     ::  ps(:,:)        ! Pa
    real, intent(out)                    ::  ll(:,:,:)
    type(TLevelInfo), intent(in)         ::  leviX
    real, intent(in)                     ::  llX(:,:,:)
    character(len=*), intent(in)         ::  combine_key
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_FillLevels_3d'

    ! --- local ------------------------------

    integer            ::  verbose
    character(len=10)  ::  fillkey
    logical            ::  reverse

    integer            ::  i, j, k, l
    integer            ::  le, nle, le1, le2
    real               ::  Xfp0, Xfp0_save

    real, allocatable  ::  phalf(:)
    real, allocatable  ::  phalfX(:)
    real, allocatable  ::  dpX(:)
    real               ::  pfac
    integer            ::  ilast

    ! --- begin -------------------------------

    ! for safety ...
    status = 1

    ! correct number of levels ?
    if ( ((nw == 'n') .and. (size(ll,3) /= levi%nlev  )) .or. &
         ((nw == 'w') .and. (size(ll,3) /= levi%nlev+1)) ) then
      write (*,'("ERROR - number of levels in output grid does not match definition:")')
      write (*,'("ERROR -   levi nlev   : ",i3)') levi%nlev
      write (*,'("ERROR -   nw          : ",a )') nw
      write (*,'("ERROR -   ll  levels  : ",i3)') size(ll,3)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! correct number of levels ?
    if ( ((nw == 'n') .and. (size(llX,3) /= leviX%nlev  )) .or. &
         ((nw == 'w') .and. (size(llX,3) /= leviX%nlev+1)) ) then
      write (*,'("ERROR - number of levels in input grid does not match definition:")')
      write (*,'("ERROR -   leviX nlev   : ",i3)') leviX%nlev
      write (*,'("ERROR -   nw           : ",a )') nw
      write (*,'("ERROR -   llX  levels  : ",i3)') size(llX,3)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! same horizontal grid sizes of ll and llX ?
    if ( (size(ll,1) /= size(llX,1)) .or. (size(ll,2) /= size(llX,2)) ) then
      write (*,'("ERROR - horizontal sizes do not match:")')
      write (*,'("ERROR -   ll     : ",i3," x ",i3)') size(ll,1), size(ll,2)
      write (*,'("ERROR -   llX    : ",i3," x ",i3)') size(llX,1), size(llX,2)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! necessary to combine or reverse ?
    verbose = 0
    call Compare( levi, leviX, verbose, fillkey, reverse, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! correct size of surface pressure field ?
    if ( (combine_key == 'mass-aver') .or. (fillkey == 'interpol') ) then
      if ( (size(ps,1) /= size(ll,1)) .or. (size(ps,2) /= size(ll,2)) ) then
        write (*,'("ERROR - horizontal sizes do not match:")')
        write (*,'("ERROR -   ps          : ",i3," x ",i3)') size(ps ,1), size(ps ,2)
        write (*,'("ERROR -   ll          : ",i3," x ",i3)') size(ll ,1), size(ll ,2)
        write (*,'("ERROR -   combine_key : ",a)') combine_key
        write (*,'("ERROR -   fillkey     : ",a)') fillkey
        write (*,'("ERROR -   reverse     : ",l1)') reverse
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
    end if

    ! levels or half-levels ?
    select case ( nw )

      !
      ! === levels ===================================
      !

      case ( 'n' )

        ! fill given key:
        select case ( fillkey )

          !
          ! ** copy levels
          !

          case ( 'copy' )

            ! loop over all TM levels
            if ( reverse ) then
              do l = 1, levi%nlev
                ! copy corresponding level in field 'llX' read from file:
                ll(:,:,l) = llX(:,:,levi%nlev+1-l)
              end do
            else
              ll = llX
            end if

          !
          ! ** combine levels
          !

          case ( 'combine' )

            ! loop over all TM levels
            ll = 0.0
            do l = 1, levi%nlev
              ! loop over range of parent levels:
              le1 = levi%flevs(l,1)
              le2 = levi%flevs(l,2)
              nle = le2 - le1 + 1
              Xfp0_save = -1.0
              do le = le1, le2
                ! index of corresponding level in field 'llX' read from file:
                k = le
                if ( reverse ) k = levi%nlev_parent + 1 - le
                ! standard pressure for level in llX :
                Xfp0 = leviX%fp0(k)
                ! based on combine key, add contribution of this level:
                select case ( combine_key )
                  case ( 'sum' )
                    ll(:,:,l) = ll(:,:,l) + llX(:,:,k)
                  case ( 'aver' )
                    ll(:,:,l) = ll(:,:,l) + llX(:,:,k)/(nle*1.0)
                  case ( 'mass-aver' )
                    !
                    !  m = dp * A / g = |da+db*ps| * A/g
                    !  ( X1*mX1 + X2*mX2 + .. ) / m
                    !  ~  ( X1*|daX1+dbX1*ps| + X2*|daX2+dbX2*ps| + .. ) / |da+db*ps|
                    !
                    ll(:,:,l) = ll(:,:,l) + llX(:,:,k)*abs(leviX%da(k)+leviX%db(k)*ps) / &
                                                       abs( levi%da(l)+ levi%db(l)*ps)
                  case ( 'bottom' )
                    if ( (Xfp0_save < 0.0) .or. (Xfp0 >= Xfp0_save) ) then
                      ll(:,:,l) = llX(:,:,k)
                      Xfp0_save = Xfp0
                    end if
                  case ( 'top' )
                    if ( (Xfp0_save < 0.0) .or. (Xfp0 <= Xfp0_save) ) then
                      ll(:,:,l) = llX(:,:,k)
                      Xfp0_save = Xfp0
                    end if
                  case default
                    write (*,'("ERROR - combine key not supported:")')
                    write (*,'("ERROR -   combine key : ",a)') trim(combine_key)
                    write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
                    write (*,'("ERROR -   nw key      : ",a)') trim(nw)
                    write (*,'("ERROR in ",a)') rname; status=1; return
                end select
              end do
            end do

          !
          ! ** interpol
          !

          case ( 'interpol' )

            allocate( phalf (0: levi%nlev) )
            allocate( phalfX(0:leviX%nlev) )
            allocate(    dpX(  leviX%nlev) )

            ! phalfX should be increasing
            pfac = 1.0
            if ( leviX%updo == 'u' ) pfac = -1.0

            select case ( combine_key )

              case ( 'sum' )

                do j = 1, size(ll,2)
                  do i = 1, size(ll,1)

                    phalf  =  levi%a  +  levi%b  * ps(i,j)    ! Pa
                    phalfX = leviX%a  + leviX%b  * ps(i,j)    ! Pa

                    ilast = 1
                    do l = 1, levi%nlev
                      ! take partly sums of f(i)*m(i) within [phalf(l-1),phalf(l)]
                      ! NOTE: if phalf(l) < phalf(l-1), the result of IntervalSum is negative
                      call IntervalSum( phalfX*pfac, llX(i,j,:), &
                                        phalf(l-1)*pfac, phalf(l)*pfac, &
                                        ll(i,j,l), ilast, status )
                      if (status/=0) then;
                        write (*,'("ERROR - from IntervalSum (combine key `sum`)")')
                        write (*,'("ERROR -   i, j    : ",2i4   )') i, j
                        write (*,'("ERROR -   ps      : ", f16.4)') ps(i,j)
                        write (*,'("ERROR -   leviX%a : ",es11.4)') leviX%a
                        write (*,'("ERROR -   leviX%b : ",es11.4)') leviX%b
                        write (*,'("ERROR in ",a)') rname; status=1; return
                      end if
                      ll(i,j,l) = ll(i,j,l) * sign(1.0,(phalf(l)-phalf(l-1))*pfac)
                    end do

                  end do
                end do

              case ( 'mass-aver' )

                do j = 1, size(ll,2)
                  do i = 1, size(ll,1)

                    phalf  =  levi%a  +  levi%b  * ps(i,j)    ! Pa
                    phalfX = leviX%a  + leviX%b  * ps(i,j)    ! Pa
                    dpX = abs(leviX%da + leviX%db * ps(i,j))    ! Pa

                    ilast = 1
                    do l = 1, levi%nlev
                      ! take partly sums of f(i)*m(i) within [phalf(l-1),phalf(l)]
                      ! NOTE: if phalf(l) < phalf(l-1), the result of IntervalSum is negative
                      call IntervalSum( phalfX*pfac, llX(i,j,:)*dpX, &
                                        phalf(l-1)*pfac, phalf(l)*pfac, &
                                        ll(i,j,l), ilast, status )
                      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                      ll(i,j,l) = ll(i,j,l) / (phalf(l)-phalf(l-1))*pfac
                    end do

                  end do
                end do

              case ( 'top' )

                do j = 1, size(ll,2)
                  do i = 1, size(ll,1)

                    phalf  =  levi%a  +  levi%b  * ps(i,j)    ! Pa, 0..
                    phalfX = leviX%a  + leviX%b  * ps(i,j)    ! Pa, 0..

                    ! loop over all levels:
                    do l = 1, levi%nlev
                      ! interpolate to top half level pressure:
                      call Interp_Lin( phalfX(1:leviX%nlev)*pfac, llX(i,j,:), &
                                       phalf (l)           *pfac,  ll(i,j,l), ilast, status )
                      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                    end do

                  end do
                end do

              case ( 'bottom' )

                do j = 1, size(ll,2)
                  do i = 1, size(ll,1)

                    phalf  =  levi%a  +  levi%b  * ps(i,j)    ! Pa, 0..
                    phalfX = leviX%a  + leviX%b  * ps(i,j)    ! Pa, 0..

                    ! loop over all levels:
                    do l = 1, levi%nlev
                      ! interpolate to bottom half level pressure:
                      call Interp_Lin( phalfX(0:leviX%nlev-1)*pfac, llX(i,j,:), &
                                       phalf (l-1)           *pfac,  ll(i,j,l), ilast, status )
                      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                    end do

                  end do
                end do

              case default
                write (*,'("ERROR - combine key not supported:")')
                write (*,'("ERROR -   combine key : ",a)') trim(combine_key)
                write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
                write (*,'("ERROR -   nw key      : ",a)') trim(nw)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

            deallocate( phalf  )
            deallocate( phalfX )
            deallocate(    dpX )

          !
          ! ** error ..
          !

          case default
            write (*,'("ERROR - fill key not supported:")')
            write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
            write (*,'("ERROR -   nw key      : ",a)') trim(nw)
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      !
      ! === half levels ========================================
      !

      case ( 'w' )

        ! fill given key:
        select case ( fillkey )

          !
          ! ** copy levels
          !

          case ( 'copy' )

            ! loop over all TM levels
            if ( reverse ) then
              do l = 1, levi%nlev+1
                ! copy corresponding level in field 'llX' read from file:
                ll(:,:,l) = llX(:,:,levi%nlev+2-l)
              end do
            else
              ll = llX
            end if

          !
          ! ** combine levels
          !

          case ( 'combine' )

            ll = 0.0

            ! loop over target half levels
            do l = 0, levi%nlev
              ! index of corresponding level in field 'llX' read from file:
              ! note: k=0,..,nlev
              k = levi%hlevs(l)
              if ( reverse ) k = levi%nlev_parent - k
              ! copy:
              ll(:,:,l+1) = llX(:,:,k+1)
            end do

          !
          ! ** interpol
          !

          case ( 'interpol' )

            allocate( phalf (0: levi%nlev) )
            allocate( phalfX(0:leviX%nlev) )
            allocate(    dpX(  leviX%nlev) )

            ! phalfX should be increasing
            pfac = 1.0
            if ( leviX%updo == 'u' ) pfac = -1.0

            select case ( combine_key )

              case ( 'top', 'bottom' )

                do j = 1, size(ll,2)
                  do i = 1, size(ll,1)

                    phalf  =  levi%a  +  levi%b  * ps(i,j)    ! Pa
                    phalfX = leviX%a  + leviX%b  * ps(i,j)    ! Pa

                    ! loop over all half levels:
                    do l = 0, levi%nlev
                      ! interpolate to half level pressure:
                      call Interp_Lin( phalfX   *pfac, llX(i,j,:  ), &
                                       phalf (l)*pfac,  ll(i,j,l+1), ilast, status )
                      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                    end do

                  end do
                end do

              case default
                write (*,'("ERROR - combine key not supported:")')
                write (*,'("ERROR -   combine key : ",a)') trim(combine_key)
                write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
                write (*,'("ERROR -   nw key      : ",a)') trim(nw)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

            deallocate( phalf  )
            deallocate( phalfX )
            deallocate(    dpX )

          !
          ! ** error ..
          !

          case default
            write (*,'("ERROR - fill key not supported:")')
            write (*,'("ERROR -   fill key    : ",a)') trim(fillkey)
            write (*,'("ERROR -   nw key      : ",a)') trim(nw)
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      !
      ! === error ========================================
      !

      case default

        write (*,'("ERROR - nw key `",a,"` not supported.")') trim(nw)
        write (*,'("ERROR in ",a)') rname; status=1; return

    end select

    ! ok
    status = 0

  end subroutine levi_FillLevels_3d



end module grid_type_hyb
