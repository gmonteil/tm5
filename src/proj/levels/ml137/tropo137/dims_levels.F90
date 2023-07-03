#include "tm5.inc"
!
! Levels.
!

module dims_levels

  use const_ec_v, only : lme, a_ec, b_ec
  use dims_grid, only : nregions_all
  
  implicit none
  
  ! --- in/out ------------------------------
  
  public
  
  
  ! --- const -------------------------------
  
  
  integer, parameter ::  lm(0:nregions_all) = 137

  ! Only lowest layers are used for convec calculations (sub files)
  ! to avoid memory problems.
#ifdef without_lmax_conv
  ! don't bother about computation time and storage ...
  integer,parameter :: lmax_conv = 137
#else
  ! Quote Michiel van Weele:
  !     "..., convectie kan wel tot ~70 hPa (19km) gaan in de tropen."
  integer,parameter :: lmax_conv = 87
#endif

  ! select ECMWF half levels;
  ! TM levels are number bottom-up:
  integer, parameter   ::  echlevs(0:137) = (/ &
                 137, 136, 135, 134, 133, 132, 131, 130, &
       129, 128, 127, 126, 125, 124, 123, 122, 121, 120, &
       119, 118, 117, 116, 115, 114, 113, 112, 111, 110, &
       109, 108, 107, 106, 105, 104, 103, 102, 101, 100, &
        99,  98,  97,  96,  95,  94,  93,  92,  91,  90, &
        89,  88,  87,  86,  85,  84,  83,  82,  81,  80, &
        79,  78,  77,  76,  75,  74,  73,  72,  71,  70, &
        69,  68,  67,  66,  65,  64,  63,  62,  61,  60, &
        59,  58,  57,  56,  55,  54,  53,  52,  51,  50, &
        49,  48,  47,  46,  45,  44,  43,  42,  41,  40, &
        39,  38,  37,  36,  35,  34,  33,  32,  31,  30, &
        29,  28,  27,  26,  25,  24,  23,  22,  21,  20, &
        19,  18,  17,  16,  15,  14,  13,  12,  11,  10, &
         9,   8,   7,   6,   5,   4,   3,   2,   1,   0 /)

  ! values filled in module geometry ...
  real  :: bt(1:137+1)
  real  :: at(1:137+1)


end module dims_levels
