!
! virtual temperature calcualtions
!
!     xm_air  = 28.964 e-3  kg air/mol
!     xm_h2o  = 18.0   e-3  kg h2o/mol
!
!     eps1 = (kg air/mol) / (kg h2o/mol) - 1 = 0.609
!
!    virtual temperature :  Tv = T * ( 1 + eps1*q )
!

module phys_tv

  implicit none
  
  ! --- in/out -----------------------------
  
  private
  
  public  ::  VirtualTemperature
  public  ::  RealTemperature
  
  
  
contains


  ! convert from real temperature to virtual temperature
  
  elemental function VirtualTemperature( T, Q )

    use Binas, only : xm_air   ! 28.964 e-3  kg air/mol
    use Binas, only : xm_h2o   ! 18.0   e-3  kg h2o/mol
    
    ! --- in/out --------------------------------
    
    real                  ::  VirtualTemperature   ! K
    real, intent(in)      ::  T                    ! real temperature (K)
    real, intent(in)      ::  Q                    ! specific humidity (kg H2O / kg air)
    
    ! --- const ---------------------------------
    
    !   eps1 = (kg air/mol) / (kg h2o/mol) - 1 = 0.609
    real, parameter  ::  eps1 = xm_air/xm_h2o - 1.0
    
    ! --- begin ---------------------------------
    
    ! Tv =  T  * ( 1 + eps1*q )
    ! T  =  Tv / ( 1 + eps1*q )
    
    VirtualTemperature = T * ( 1.0 + eps1 * Q )   ! K

  end function VirtualTemperature


  ! convert from virtual temperature to temperature
  
  elemental function RealTemperature( Tv, Q )

    use Binas, only : xm_air   ! 28.964 e-3  kg air/mol
    use Binas, only : xm_h2o   ! 18.0   e-3  kg h2o/mol
    
    ! --- in/out --------------------------------
    
    real                  ::  RealTemperature    ! K
    real, intent(in)      ::  Tv                 ! virtual temper (K)
    real, intent(in)      ::  Q                  ! specific humidity (kg H2O / kg air)
    
    ! --- const ---------------------------------
    
    !   eps1 = (kg air/mol) / (kg h2o/mol) - 1 = 0.609
    real, parameter  ::  eps1 = xm_air/xm_h2o - 1.0
    
    ! --- begin ---------------------------------
    
    ! Tv =  T  * ( 1 + eps1*q )
    ! T  =  Tv / ( 1 + eps1*q )
    
    RealTemperature = Tv / ( 1.0 + eps1 * Q )   ! K

  end function RealTemperature   


end module phys_tv
