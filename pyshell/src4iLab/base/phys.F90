!
!ProTeX: 1.14-AJS
!
!BOI
!
! !TITLE:        Phys  -  Physics library.
! !AUTHORS:      Arjo Segers
! !AFFILIATION:  KNMI
! !DATE:         12/08/2004
!
!
! !INTRODUCTION: Introduction
!
! The 'phys' library contains subroutines for model physica.
!
!
! !INTRODUCTION: Geo-potential-height
!
! \begin{itemize}
!
! \item Potential height of an air parcel given top and bottom pressure,
!   temperature, and humidity:
!   \bv
!     subroutine PotentialHeight( ptop, pdown, T, Q, dh )
!       real, intent(in)        ::  ptop, pdown   ! pressure bounds (Pa)
!       real, intent(in)        ::  T             ! temperature (K)
!       real, intent(in)        ::  Q             ! specific humidity (kg h2o/kg air)
!       real, intent(out)       ::  dh            ! geo.pot.height bounds (m)
!     end subroutine PotentialHeight
!   \ev
!
! \item Geo-potential-height at half level boundaries:
!   \bv
!     subroutine GeoPotentialHeightB( lm, pb, T, Q, h0, gphb )
!       integer, intent(in)            ::  lm         ! number of levels
!       real, intent(in)               ::  pb(lm+1)   ! pressure bounds (Pa)
!       real, intent(in)               ::  T(lm)      ! temperature (K)
!       real, intent(in)               ::  Q(lm)      ! specific humidity (kg h2o/kg air)
!       real, intent(in)               ::  h0         ! initial height  (m)
!       real, intent(out)              ::  gphb(lm+1) ! geo.pot.height bounds (m)
!   \ev
!
! \item Geo-potential-height at full levels.
!   A full level has a pressure in between the top and bottom pressure of a layer.
!   Temperature and humidity are assumed to be constant in a layer.
!   \bv  
!     subroutine GeoPotentialHeight( updo, lm, pb, T, Q, h0, gph )
!       integer, intent(in)            ::  lm
!       real, intent(in)               ::  pb(lm+1)   ! pressure bounds (Pa)
!       real, intent(in)               ::  T(lm)      ! temperature (K)
!       real, intent(in)               ::  Q(lm)      ! specific humidity (kg h2o/kg air)
!       real, intent(in)               ::  h0         ! initial height  (m)
!       real, intent(out)              ::  gph(lm)    ! geo.pot.height (m)
!   \ev
!
! \end{itemize}
!
!
! !INTRODUCTION: Humidity functions
!
! \begin{itemize}
!
! \item Saturation specific humidity:
!   \bv
!     real function QSat( T, p )
!       real, intent(in)   ::  T     ! temperature (K)
!       real, intent(in)   ::  p     ! pressure (Pa)
!     end function QSat
!   \ev
!
! \item Derivative of saturation specific humidity with
!            respect to temperature:
!   \bv
!     real function dQSat_dT( T, p )
!       real, intent(in)   ::  T     ! temperature (K)
!       real, intent(in)   ::  p     ! pressure (Pa)
!     end function dQSat_dT
!   \ev
!
! \item Specific humidity given relative humidity, temperature, and pressure:
!   \bv
!     real function QH2O( R, T, p )
!       real, intent(in)  ::  R     ! rel. humidity (%)   
!       real, intent(in)  ::  T     ! temperature (in K)  
!       real, intent(in)  ::  p     ! pressure (Pa)       
!     end function QH2O
!   \ev
!
!
! !INTRODUCTION: Convective clouds
!
!   \bv
!     subroutine ConvCloudDim( updo, lm, detu, entd, &
!                              iclbas, ictop, icllfs, &
!                              status )
!
!       character(len=1), intent(in)   ::  updo
!       integer, intent(in)            ::  lm
!
!       real, intent(in)              ::  detu(lm)
!       real, intent(in)              ::  entd(lm)
!
!       character(len=1), intent(in)  ::  updo
!
!       ! cloud base, top, level of free sinking
!       integer, intent(out)          ::  iclbas
!       integer, intent(out)          ::  ictop
!       integer, intent(out)          ::  icllfs
!
!       integer, intent(out)          ::  status
!
!   \ev
!
!
! !INTRODUCTION: Cloud cover
!
! \bv
!   subroutine cf_overhead( nlev, yclfr, wccro, scheme, eps )
!
!     integer, intent(in)                     ::  nlev
!     real, intent(in)                        ::  yclfr(nlev)
!     real, intent(out)                       ::  wccro(nlev)
!     character(len=*), intent(in), optional  ::  scheme
!     real, intent(in), optional              ::  eps
!
!   input: 
!     nlev  : number of vertical levels
!     yclfr : cloud fraction (cc) per cell (0-1)
!  
!   output:
!     wccro: overhead cloud fraction
!  
!   optional arguments:
!     scheme='ecmwf'  : 'ecmwf' -> iovln=1
!                       'other' -> iovln=0
!     eps=1.0e-4      : cltres
!  
!   parameters:
!     iovln : switch
!       1 = ecmwf (maximum random overlap assumption) scheme
!       0 = another scheme 
!     cltres : threshold (minimum) cloud fraction used
!              for numerical stability (division by zero
!              and to eliminate small unrealistic cloud fractions
!  
!   Notes: 
!   - Index=1 of arrays (yclfr) corresponds to model top
!   - The clouds are supposed to be distributed homogeneously
!     in the vertical in each layer.
! \ev
!
! 
! !INTRODUCTION: Virtual temperature
!
! \begin{itemize}
!
! \item Convert from real temperature and specific humidity to virtual temperature:
!   \bv
!     elemental function VirtualTemperature( T, Q )
!       real                  ::  VirtualTemperature   ! K
!       real, intent(in)      ::  T                    ! real temperature (K)
!       real, intent(in)      ::  Q                    ! specific humidity (kg H2O / kg air)
!     end function VirtualTemperature 
!   \ev
!
! \item Convert from virtual temperature and specific humidity to temperature:
!   \bv
!     elemental function RealTemperature( Tv, Q )
!       real                  ::  RealTemperature    ! K
!       real, intent(in)      ::  Tv                 ! virtual temper (K)
!       real, intent(in)      ::  Q                  ! specific humidity (kg H2O / kg air)
!     end function RealTemperature
!   \ev
!
! \end{itemize}
!
!EOI
!

!
! Physics library.
!

module Phys

  use Phys_GPH
  use Phys_Humidity
  use Phys_Convec
  use Phys_CloudCover
  use Phys_Tv
  use Phys_sza
  
  public

end module Phys
