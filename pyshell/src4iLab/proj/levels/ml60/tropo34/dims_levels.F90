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

  ! lmax_conv determines to what level the clouds and convection
  ! are stored in the meteo files...

#ifdef without_lmax_conv
  ! don't bother about computation time and storage ...
  integer,parameter :: lmax_conv = 34
#else
  ! Quote Michiel van Weele:
  !     "..., convectie kan wel tot ~70 hPa (19km) gaan in de tropen."
  integer,parameter :: lmax_conv = 28
#endif

  integer,parameter :: echlevs(0:lm(1)) = (/ &
       60,     58,       &
       56,     54,       &
       52,     50,       &
       48,     46,       &
       44,     42,       &
       40,     39,       &
       38,     37,       &
       36,     35,       &
       34,     33,       &
       32,     31,       &
       30,     29,       &
       28,     27,       &
       26,     25,       &
       24,     23,       &
       22,     20,       &
       16,               &
       12,               &
       8,                &
       4,                &
       0 /)

  ! values filled in module geometry ...

  real, parameter  ::   bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  ::   at(1:lm(1)+1) = a_ec(echlevs)

end module dims_levels
