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

  ! Only lowest layers are used for convec calculations (sub files)
  ! to avoid memory problems.
#ifdef without_lmax_conv
  ! don't bother about computation time and storage ...
  integer,parameter :: lmax_conv = 34
#else
  ! Quote Michiel van Weele:
  !     "..., convectie kan wel tot ~70 hPa (19km) gaan in de tropen."
  integer,parameter :: lmax_conv = 23
#endif

  ! select ECMWF half levels;
  ! TM levels are number bottom-up:
  !
  ! Selection of 34 half levels that are most similar to:
  !   levels/ml60/tropo34/dims_levels__ml60_tropo34.F90
  !
  integer, parameter   ::  echlevs(0:lm(1)) = (/ &
       137, 134, 131, 126, 122, 117, 113, 109, 105, 102, &
        98,  96,  94,  92,  90,  88,  85,  83,  81,  78, &
        75,  73,  70,  67,  64,  61,  57,  54,  51,  44, &
        34,  26,  20,  14,   0 /)

  ! values filled in module geometry ...
  real, parameter  ::   bt(1:lm(1)+1) = b_ec(echlevs)
  real, parameter  ::   at(1:lm(1)+1) = a_ec(echlevs)


end module dims_levels

