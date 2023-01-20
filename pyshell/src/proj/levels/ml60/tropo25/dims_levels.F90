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

  integer, parameter ::  lm(nregions) = 25


  integer,parameter :: lmax_conv = 19   ! can be lowered to tropopause...

  integer,parameter :: echlevs(0:lm(1)) = (/ &
                     60,     58,       &
                     56,     54,       &
                     52,     50,       &
                     48,     46,       &
                     44,     42,       &
                     40,     38,       &
                     36,     34,       &
                     32,     30,       &
                     28,     26,       &
                     24,     22,       &
                     20,               &
                     16,               &
                     12,               &
                      8,               &
                      4,               &
                      0 /)

  real, parameter  ::   bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  ::   at(1:lm(1)+1) = a_ec(echlevs)


end module dims_levels
