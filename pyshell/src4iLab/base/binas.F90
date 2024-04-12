module Binas

  implicit none

  public

  !
  !ProTeX: 1.14-AJS
  !
  !BOI
  !
  ! !TITLE:        Binas - geometrical and physical constants
  ! !AUTHORS:      Arjo Segers
  ! !AFFILIATION:  KNMI
  ! !DATE:         \today
  !
  ! !INTRODUCTION: Introduction
  !
  !   'BINAS' is the name an in The Netherlands common used table-book
  !   with scientific constants and formulae.
  !
  !
  ! !INTRODUCTION: Constants
  !
  !BOC

  ! ---------------------------------------------------------------
  ! gonio
  ! ---------------------------------------------------------------

  ! defintions for pi :
  !  o old definition:
  !real, parameter         ::  pi    = 3.1415927
  !  o EMOS definition (emos/interpolation/strlat.F, parameter PPI)
  real, parameter         ::  pi = 3.14159265358979

  ! two pi :
  real, parameter         ::  twopi = 2*pi

  ! factor to convert to radians from degrees:
  real, parameter         ::  deg2rad = pi/180.0     ! rad/deg


  ! ---------------------------------------------------------------
  ! earth
  ! ---------------------------------------------------------------

  ! Radius of earth as used in EMOS library (ECMWF model),
  ! see for example "jvod2uv.F"
  ! NOTE: the value 6.375e6 was used in TM !
  real, parameter         ::  ae = 6.371e6     ! m

  ! acceleration of gravity:
  !real, parameter         ::  grav = 9.81       ! m/s2
  real, parameter         ::  grav = 9.80665    ! m/s2

  ! Earth's angular speed of rotation
  !  Omega =  2 * pi * (365.25/364.25) / (86400.0)
  real, parameter         ::  Omega = 7.292e-5    ! rad/s


  ! ---------------------------------------------------------------
  ! molecules, mols, etc
  ! ---------------------------------------------------------------

  ! Avogadro number
  real, parameter        ::  Avog = 6.02205e23      ! mlc/mol

  ! Dobson units:
  real,parameter         ::  Dobs = 2.68668e16    ! (mlc/cm2) / DU


  !
  ! molar weights of components
  !

  ! naming convention:
  !  o old names 'xm***' are in g/mol
  !  o new names 'xm_***' are in kg/mol
  !

  ! atomic weights:
  real, parameter        ::  xm_h     =    1.00790e-3     ! kg/mol
  real, parameter        ::  xm_n     =   14.00670e-3     ! kg/mol
  real, parameter        ::  xm_c     =   12.01115e-3     ! kg/mol
  real, parameter        ::  xm_s     =   32.06400e-3     ! kg/mol
  real, parameter        ::  xm_o     =   15.99940e-3     ! kg/mol
  real, parameter        ::  xm_cl    =   35.45300e-3     ! kg/mol
  real, parameter        ::  xm_rn222 =  222.0e-3         ! kg/mol
  real, parameter        ::  xm_pb210 =  210.0e-3         ! kg/mol

  ! molecule weights:
  real, parameter        ::  xm_h2o   =  xm_h*2 + xm_o    ! kg/mol
  real, parameter        ::  xm_o3    =  xm_o * 3         ! kg/mol

  ! mass of air
  real, parameter        ::  xm_air   =  28.964e-3        ! kg/mol
  real, parameter        ::  xmair    =  28.94            ! g/mol; old name!

  ! atomic mass unit
  real, parameter        :: amu = 1.66053892D-27          ! kilograms

  !                  mlc/mol
  ! [cdob] = ------------------------ =   DU / (kg/m2)
  !          kg/mol cm2/m2 mlc/cm2/DU
  !

  real, parameter :: cdob_o3 = Avog / ( xm_o3 * 1.0e4 * Dobs )  ! DU/(kg/m2)

  ! ---------------------------------------------------------------
  ! gas
  ! ---------------------------------------------------------------

  ! gas constant
  real, parameter        ::  Rgas = 8.3144     ! J/mol/K

  ! gas constant for dry air
  !real, parameter        ::  rgas_x  = 287.05
  ! NEW:
  !   Rgas_air = Rgas / xmair = 287.0598
  real, parameter        ::  Rgas_air = Rgas / xm_air    ! J/kg/K

  ! water vapour
  !real,parameter         ::  rgasv = 461.51
  real, parameter        ::  Rgas_h2o = Rgas / xm_h2o   ! J/kg/K

  ! standard pressure
  real, parameter        ::  p0 = 1.0e5    ! Pa
  !real, parameter        ::  p0 = 1.01325e5    ! Pa  <-- suggestion Bram Bregman

  ! global mean pressure:
  real,parameter         ::  p_global = 98500.0   ! Pa

  ! specific heat of dry air at constant pressure
  !real, parameter        ::  cp0 = 1004.0           ! J/kg/K
  real, parameter        ::  cp_air = 1004.0           ! J/kg/K

  ! Latent heat of evaporation
  real, parameter        ::  lvap = 2.5e6     !  [J kg-1]

  ! Latent heat of condensation at 0 deg Celcius
  ! (heat (J) necesarry to evaporate 1 kg of water)
  real, parameter       ::  Lc = 22.6e5           ! J/kg

  ! kappa = R/cp = 2/7
  real, parameter        ::  kappa = 2.0/7.0
  ! 'kapa' is probably 'kappa' ....
  !real, parameter        ::  kapa = 0.286

  ! Von Karman constant (dry_dep)
  real, parameter        ::  vkarman = 0.4

  ! Boltzmann constant:
  real, parameter           ::  kbolz = 1.38066e-23    ! J/K

  ! ---------------------------------------------------------------
  ! virtual temperature :  Tv = T * ( 1 + eps1*q )
  ! ---------------------------------------------------------------

  real, parameter        ::  eps  = Rgas_air / Rgas_h2o
  real, parameter        ::  eps1 = ( 1.0 - eps )/eps


  ! ---------------------------------------------------------------
  ! other
  ! ---------------------------------------------------------------

  ! melting point
  real, parameter        ::  T0 = 273.16    ! K

  ! Rv/Rd
  real, parameter        ::  gamma = 6.5e-3

  ! density of pure water at 0 deg C
  real, parameter                  ::  rol   = 1000.         ! kg/m^3

  ! density of ice water at 0 deg C
  real, parameter                  ::  roi   = 917.          ! kg/m^3

  ! Planck times velocity of light
  real, parameter                  ::  hc = 6.626176e-34 * 2.997924580e8  ! Jm


  ! ---------------------------------------------------------------
  ! end
  ! ---------------------------------------------------------------

  !EOC

end module Binas
