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

  ! lmax_conv determines to what level the clouds and convection 
  ! are stored in the meteo files...
#ifdef without_lmax_conv
    ! don't bother about computation time and storage ...
  integer,parameter :: lmax_conv = 25
#else
  ! Quote Michiel van Weele:
  !     "..., convectie kan wel tot ~70 hPa (19km) gaan in de tropen."
  integer,parameter :: lmax_conv = 19
#endif

  ! selection of 25 layers as subset of ml137/tropo34a,
  ! similar as ml60/tropo25 was selected out of ml60/tropo34 :
  integer,parameter :: echlevs(0:lm(1)) = (/ &
      137, 134, 131, 126, 122, 117, 113, 109, 105, 102, &
       98,       94,       90,  88,            81,      &
       75,       70,       64,       57,       51,  44, &
       34,  26,  20,  14,   0 /)

  real, parameter  ::   bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  ::   at(1:lm(1)+1) = a_ec(echlevs)


end module dims_levels
