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

  integer, parameter ::  lm(nregions) = 36


  integer,parameter :: lmax_conv = 24   ! can be lowered to tropopause...

  integer,parameter :: echlevs(0:lm(1)) = (/ &
                     60,     59,    58,       &
                     56,     55,    54,       &
                     52,     51,    50,       &
                     48,     47,    46,       &
                     44,     43,    42,       &
                     40,     39,    38,       &
                     36,     35,    34,       &
                     32,     31,    30,       &
                     28,     27,    26,       &
                     24,     22,    20,       &
                     18,     16,    12,       &
                     10,      8,     4,       &
                      0 /)

  real, parameter  ::   bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  ::   at(1:lm(1)+1) = a_ec(echlevs)


end module dims_levels
