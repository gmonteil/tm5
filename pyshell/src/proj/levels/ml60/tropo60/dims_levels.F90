!###############################################################################
!
! Levels.
!
!### macro's ###################################################################
!
#include "tm5.inc"
!
!###############################################################################

module dims_levels

  use const_ec_v, only : lme, a_ec, b_ec

  use dims_grid, only : nregions

  implicit none

  ! --- in/out ------------------------------

  public


  ! --- const -------------------------------

  integer, parameter ::  lm(nregions) = 60


  !integer,parameter :: lmax_conv = 40   ! Sourish Basu -- more or less arbitrarily set
  integer,parameter :: lmax_conv = lm(1)   ! Sourish Basu -- to test ECMWF diffusion coefficients in the stratosphere

  integer,parameter :: echlevs(0:lm(1)) = (/ &
                    60, 59, 58, 57, &
                    56, 55, 54, 53, &
                    52, 51, 50, 49, &
                    48, 47, 46, 45, &
                    44, 43, 42, 41, &
                    40, 39, 38, 37, &
                    36, 35, 34, 33, &
                    32, 31, 30, 29, &
                    28, 27, 26, 25, &
                    24, 23, 22, 21, &
                    20, 19, 18, 17, &
                    16, 15, 14, 13, &
                    12, 11, 10,  9, &
                     8,  7,  6,  5, &
                     4,  3,  2,  1, &
                     0 /)

  real, parameter  ::   bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  ::   at(1:lm(1)+1) = a_ec(echlevs)


end module dims_levels
