!
! Solar zenith stuff.
!

module phys_sza

  implicit none
  
  ! --- in/out ------------------------------
  
  private
  
  public  ::  cos_sza
  
  
contains


  ! =======================================================
  
  
  !
  ! Return cos( solar_zenith_angle ) given :
  !  o daynr : januari 1 = 1, ..., december 31 = 365|366
  !  o hour, minutes, seconds
  !  o lon, lat in degrees
  !

  real function cos_sza( daynr, hour, minu, sec, lon, lat )

    use binas, only : pi, deg2rad
    
    ! --- in/out --------------------------------

    integer, intent(in)   ::  daynr
    integer, intent(in)   ::  hour, minu, sec
    real, intent(in)      ::  lon        ! deg  [-180,180]
    real, intent(in)      ::  lat        ! deg  [ -90, 90]
    
    ! --- const ---------------------------------
    
    real, parameter   ::  piby  = pi / 180.0
    real, parameter   ::  obliq = 23.45 * piby

    ! --- local ---------------------------------
    
    real    ::  deday, delta
    real    ::  time
    real    ::  lonnoon
    
    ! --- begin -------------------------------------
    
    deday = 4.88 + 2.0*pi*real(daynr)/365.0 
    delta = asin( sin(obliq) * sin(deday) )
    
    ! seconds after 00:00
    time = hour*3600.0 + minu*60.0 + sec

    ! longitude at which it is noon at current time:
    !    00:00    pi      180 E
    !    06:00    pi/2     90 E
    !    12:00     0        0
    !    18:00   -pi/2     90 W
    !    24:00   -pi      180 W
    lonnoon = pi - 2.0*pi * real(time)/86400.0

    ! result:
    cos_sza = sin(delta)*sin(lat*deg2rad) + cos(delta)*cos(lat*deg2rad)*cos(lon*deg2rad-lonnoon)

  end function cos_sza
  

end module phys_sza
