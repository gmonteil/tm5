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

  integer, parameter ::  lm(nregions) = 34

  integer,parameter :: lmax_conv = 23   ! can be lowered to tropopause...

  integer, parameter   ::  echlevs(0:lm(1)) = (/91, 89, 86,  &
                            83, 80, 77,&
                            74, 71, 68,&
                            65, 62, 59,&
                            56, 54, 52,&
                            50, 48, 46,&
                            44, 42, 40,&
                            38, 36, 34,&
                            32, 29 ,26,&
                            23, 20 ,17,&
                            14, 11 ,8 ,&
                             5,  0/)

  real, parameter  :: bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  :: at(1:lm(1)+1) = a_ec(echlevs)


end module dims_levels
