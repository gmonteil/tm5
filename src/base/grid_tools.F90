!
! Grid tools.
!
! may 2002, Arjo Segers
!

module grid_tools


  use Binas, only : deg2rad, pi, ae
  
  implicit none
  
  
  ! --- in/out --------------------------------------
  
  private
  
  public  ::  deg2rad, pi, ae
  
  public  ::  ll_area
  public  ::  ll_area_frac
  
  

contains


  ! ===================================================================
    
  
  ! given rectangle [west,east]x[south,north] in radians,
  ! compute area in rad^2
  
  real function ll_area( west, east, south, north )
  
    ! --- in/out -------------------------------------
    
    real, intent(in)  ::  west, east, south, north    ! rad
    
    ! --- begin ------------------------------------

    ll_area = ( sin(max(north,south)) - sin(min(north,south)) ) * abs(east-west)
    
  end function ll_area
  
  
  ! ===
  
  
  ! Compute fraction of rectangle 1 that covers rectangle 2,
  ! both defined in radians.
  ! Fraction is equal to ratio of shaded area and area 1.
  !
  !   +-----------------+         |
  !   |           +-----+---------+--
  !   |    1      |/////|         |
  !   +-----------------+   2     |
  !               |               |
  !             --+---------------+--
  !               |               |
  !
  
  real function ll_area_frac( west1, east1, south1, north1, &
                               west2, east2, south2, north2 )
    
    ! --- in/out ----------------------------
    
    real, intent(in)      ::  west1, east1, south1, north1   ! rad
    real, intent(in)      ::  west2, east2, south2, north2   ! rad

    ! --- local ---------------------------

    real        ::  xwest2, xeast2
    
    ! --- begin -----------------------------

    ! check ..    
    if ( east1 < west1 .or. east2 < west2 .or. &
         north1 < south1 .or. north2 < south2 ) then
      print *, 'found strange area:'
      print *, '  1:', west1, east1, south1, north1
      print *, '  2:', west2, east2, south2, north2
      stop 'FATAL ERROR IN ll_area_frac'
    end if
    
    ! shift rect 2 over 360 deg if it does not cover rect 1 at all;
    ! if they still not cover, the fraction will be set to zero later on:
    if ( west2 > east1 ) then
      ! cell 2 is east of cell 1; try to shift 360.0 westwards
      xwest2 = west2 - 2*pi
      xeast2 = east2 - 2*pi
    else if ( east2 < west1 ) then
      ! cell 2 is west of cell 1; try to shift 360.0 eastwards
      xwest2 = west2 + 2*pi
      xeast2 = east2 + 2*pi
    else
      ! just copy ...
      xwest2 = west2
      xeast2 = east2
    end if

    ! compute fraction; zero if rectangles do not cover:
    if ( (xeast2 <= west1 ) .or. (xwest2 >= east1 ) .or. &
         (north2 <= south1) .or. (south2 >= north1) ) then
      ll_area_frac = 0.0
    else
      ll_area_frac = &
         ll_area( max(west1 ,xwest2), min(east1 ,xeast2), &
                  max(south1,south2), min(north1,north2) ) &
          / ll_area( west1, east1, south1, north1 )
    end if

  end function ll_area_frac
  

end module grid_tools
  
