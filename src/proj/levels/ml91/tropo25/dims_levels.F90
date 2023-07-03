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


  integer,parameter ::  echlevs(0:lm(1)) = (/91, 89, 87,&
                         84, 81, 79,&
                         76, 74, 71,&
                         68, 66, 63,&
                         60, 57, 53,&
                         50, 46, 41,&
                         36, 32, 27,&
                         21, 16 ,12,&
                          9,  0/)


  real, parameter  ::   bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  ::   at(1:lm(1)+1) = a_ec(echlevs)

end module dims_levels
