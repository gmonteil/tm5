!
! Humidity functions.
!
! Copied from TMPP module 'msub_subg' (tmpp_sub_subg.f90)
!
! Note by Bas Henzing: 
!   Do not change the constants Rgas etc !
!   Thus, do not use them from Binas for example.
!   Reason: the constants c1 etc in the parameterisation are fits
!   given the coded values of Rgas etc .
!

module phys_humidity

  implicit none
  
  ! --- in/out ------------------------------------
  
  private
  
  public  ::  QSat
  public  ::  dQSat_dT
  public  ::  QH2O
  
  
contains


  !
  ! calculate saturation specific humidity
  !

  real function QSat( T, p )

    ! --- in/out ----------------------------

    real, intent(in)   ::  T     ! temperature (K)
    real, intent(in)   ::  p     ! pressure (Pa)

    ! --- const ----------------------------

    real, parameter    ::  rgasd = 287.05
    real, parameter    ::  rgasv = 461.51
    real, parameter    ::  eps = rgasd/rgasv

    real, parameter    ::  T0 = 273.16

    real, parameter    ::  c1  = 610.78
    real, parameter    ::  c3a =  17.269
    real, parameter    ::  c3b =  21.875
    real, parameter    ::  c4a =  35.86
    real, parameter    ::  c4b =   7.66

    ! --- local -------------------------------

    real              ::  es

    ! --- begin -----------------------------

    if ( p <= 0.0 ) then
      QSat = 0.0
      return
    end if

    ! set function 'qsat' equal 0 for temperatures t < 9K
    ! to ensure numerical stability 
    ! (the argument of the following exponential 
    ! function would otherwise exceeds maximum)

    if ( T < 9.0 ) then
      QSat = 0.
      return
    end if

    if ( T >= t0 ) then
      es = c1*exp(c3a*(T-T0)/(T-c4a))
    else
      es = c1*exp(c3b*(T-T0)/(T-c4b))
    end if

    QSat = eps / ( (p/es) - (1.0-eps) )

  end function QSat



  !
  ! calculate derivative of saturation specific humidity with
  ! respect to temperature
  !

  real function dQSat_dT( T, p )

    ! --- in/out -------------------------------

    real, intent(in)   ::  T     ! temperature (K)
    real, intent(in)   ::  p     ! pressure (Pa)

    ! --- const --------------------------------

    real, parameter    :: rgasd = 287.05
    real, parameter    :: rgasv = 461.51
    real, parameter    :: eps = rgasd/rgasv

    real, parameter    :: T0  = 273.16
    real, parameter    :: c1  = 610.78
    real, parameter    :: c3a =  17.269
    real, parameter    :: c3b =  21.875
    real, parameter    :: c4a =  35.86
    real, parameter    :: c4b =   7.66

    ! --- local ------------------------------

    real         ::  es
    real         ::  qsat

    ! --- begin -------------------------------

    if ( p < 0.0 ) then
      dQSat_dT = 0.0
      return
    endif

    ! set function 'dqsatdt' equal 0 for temperatures t less than 9K to ensure
    ! numerical stability (the argument of the following exponential 
    ! function would otherwise exceeds maximum)

    if ( t < 9.0 ) then
       dQSat_dT = 0.0
       return
    end if

    if ( T >= T0 ) then
      es = c1*exp(c3a*(T-T0)/(T-c4a))
      qsat = eps/((p/es)-(1-eps))
      dQSat_dT = c3a*(T0-c4a)*qsat/((T-c4a)*(T-c4a)*(1.0-(1.0-eps)*es/p))
    else
      es = c1*exp(c3b*(T-T0)/(T-c4b))
      qsat = eps/((p/es)-(1-eps))
      dQSat_dT = c3a*(T0-c4a)*qsat/((T-c4a)*(T-c4a)*(1.0-(1.0-eps)*es/p))
    end if

  end function dQSat_dT



  !
  ! calculate specific humidity
  !

  real function QH2O( r, T, p )

    ! --- in/out -------------------------------

    real, intent(in)    ::  R     ! rel. humidity (%)
    real, intent(in)    ::  T     ! temperature (in K)
    real, intent(in)    ::  p     ! pressure (Pa)

    ! --- const --------------------------------

    real, parameter    :: rgasd = 287.05
    real, parameter    :: rgasv = 461.51
    real, parameter    :: eps = rgasd/rgasv

    real, parameter    :: T0  = 273.16
    real, parameter    :: c1  = 610.78
    real, parameter    :: c3a =  17.269
    real, parameter    :: c3b =  21.875
    real, parameter    :: c4a =  35.86
    real, parameter    :: c4b =   7.66

    ! --- local ------------------------------

    real     ::  es

    ! --- begin -------------------------------

    if ( p <= 0.0 ) then
      QH2O = 0.0
      return
    endif

    ! set function 'qh2o' equal 0 for temperatures t less than 9K to ensure
    ! numerical stability (the argument of the following exponential 
    ! function would otherwise exceeds maximum)

    if ( T < 9.0 ) then
       qh2o = 0.0
       return
    endif

    if ( T >= T0 ) then
      es = c1*exp(c3a*(T-T0)/(T-c4a))
    else
      es = c1*exp(c3b*(T-T0)/(T-c4b))
    endif

    QH2O = (R/100.0)*eps / ((p/es)-(1.0-eps))

  end function QH2O


end module phys_humidity

      

