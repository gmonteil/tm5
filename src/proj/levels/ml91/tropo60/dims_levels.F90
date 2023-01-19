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

  integer,parameter :: lmax_conv = 40   ! can be lowered to tropopause...

  integer, parameter  ::  echlevs(0:lm(1)) = (/ 91, 90, 89, 88, 87, 86, 85, &
                                                  84, 83, 82, 80, 78, 76, 74, &
                                                  72, 70, 68, 66, 64, 62, 60, &
                                                  58, 57, 56, 55, 54, 53, 52, &
                                                  51, 50, 49, 48, 47, 46, 45, &
                                                  44, 43, 42, 40, 39, 38, 36, &
                                                  34, 32, 30, 28, 26, 24, 22, &
                                                  20, 18, 16, 14, 12, 10,  8, &
                                                   6,  4,  2,  1,  0 /)

  real, parameter  :: bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  :: at(1:lm(1)+1) = a_ec(echlevs)


end module dims_levels
