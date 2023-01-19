module Phys_GPH

  implicit none
  
  ! --- in/out ---------------------------------
  
  private
  
  public  ::  PotentialHeight
  public  ::  GeoPotentialHeightB
  public  ::  GeoPotentialHeight
  
  
contains


  ! ===============================================
  

  !
  ! Compute height of air parcel between pressures p0 and p1:
  !
  !  h+dh  +----+  p 
  !        |    |                                                   
  !        |    |        temperature constant within parcel 
  !        |    |        specific humidity Q  (kg h2o / kg air)
  !        |    |                                                   
  !  h     +----+  p+dp  
  !         A = 1 m2   
  !
  ! Constants:
  !
  !   g   : gravity
  !
  !  eps1  =  xm_air/xm_h2o - 1  
  !        =  28.964/18.0 - 1.0  =  0.609
  !
  ! Mass of parcel:
  !
  !   dm  =  A dp/g    [kg air]
  !
  ! Number of mol:
  !
  !   #mol =  #mol_dry_air + #mol_h2o
  !
  !        =  (kg dry air) / (kg dry air/mol dry air) + (kg h2o) / (kg h2o/mol h2o)
  !
  !        =    dm (1-Q)   /        xm_air            +   dm Q   /      xm_h2o
  !
  !        =    dm   ( (1-Q)/xm_air + Q/xm_h2o   )
  !
  !        =  A dp/g ( 1 - Q + Q xm_air/xm_h2o   ) / xm_air
  !
  !        =  A dp/g ( 1 + (xm_air/xm_h2o - 1) Q ) / xm_air
  !
  !        =  A dp/g ( 1 +         eps1        Q ) / xm_air
  !
  ! General gas law:
  !
  !      p   V    =      R  T           #mol
  ! 
  !      p  A dh  =      R  T  A dp/g (1 + eps1 Q)/xm_air
  !
  !           dh  =      R/xm_air T(1+eQ)/p dp/g
  !
  ! Integrated height:
  !
  !   h = int dh  =  int R/xm_air T(1+eQ)/p dp/g 
  !
  !                                      p_down
  !               =  R/xm_air T(1+eQ)  (  int   1/p dp ) / g
  !                                      p=p_up
  !
  !                                             p_down
  !               =  R/xm_air T(1+eQ)  [ ln(p) ]
  !                                             p=p_up
  !
  !               =  R/M_air T(1+eQ)  ln(p_down/p_up) / g
  !
  !
  ! Usefull initial value:
  !
  !   h0  = orography/g       with orography in m*[g]=m2/s2
  !
  !
  ! Result:               GeoPotentialHeightB     GeoPotentialHeight
  !
  !         pb(n+1)    +-----+  gph(n+1)
  !                    |     |
  !  T(n), Q(n)        |     |                      gph(n)
  !                    |     |
  !         pb(n)      +-----+  gph(n)
  !                    |     |
  !  T(n-1), Q(n-1)    |     |                      gph(n-1)
  !                    |     |
  !         pb(n-1)    +-----+  gph(n-1)
  !                    |     |
  !                    :     :
  !                    |     |
  !         pb(2)      +-----+  gph(2)
  !                    |     |
  !                    |     |                      gph(1)
  !                    |     |
  !         pb(1)      +-----+  gph(1) = h0
  !   
  !
  
  
  subroutine PotentialHeight( ptop, pdown, T, Q, dh )
  
    use Binas, only : grav     !  9.80665    m/s2
    use Binas, only : xm_air   ! 28.964 e-3  kg air/mol
    use Binas, only : xm_h2o   ! 18.0   e-3  kg h2o/mol
    use Binas, only : Rgas     !  8.3144     J/mol/K
    
    ! --- in/out ---------------------------------------
    
    real, intent(in)               ::  ptop, pdown   ! pressure bounds (Pa)
    real, intent(in)               ::  T             ! temperature (K)
    real, intent(in)               ::  Q             ! specific humidity (kg h2o/kg air)
    real, intent(out)              ::  dh            ! geo.pot.height bounds (m)
    
    ! --- const ----------------------------------------
    
    !   eps1 = (kg air/mol) / (kg h2o/mol) - 1 = 0.609
    real, parameter  ::  eps1 = xm_air/xm_h2o - 1.0
    
    ! gas constant for dry air
    !   Rgas_air = Rgas / xm_air = 287.0598
    real, parameter  ::  Rgas_air = Rgas / xm_air 
  
    ! --- begin -----------------------------------------
    
    ! geo potential height:
    dh = T * ( 1.0 + eps1 * Q ) * Rgas_air * log(pdown/ptop) / grav
    
  end subroutine PotentialHeight
  
  
  ! ***

  
  ! GPH at pressure half levels
  
  subroutine GeoPotentialHeightB( lm, pb, T, Q, h0, gphb )
  
    ! --- in/out ---------------------------------------
    
    integer, intent(in)            ::  lm
    real, intent(in)               ::  pb(lm+1)   ! pressure bounds (Pa)
    real, intent(in)               ::  T(lm)      ! temperature (K)
    real, intent(in)               ::  Q(lm)      ! specific humidity (kg h2o/kg air)
    real, intent(in)               ::  h0         ! initial height  (m)
    real, intent(out)              ::  gphb(lm+1) ! geo.pot.height bounds (m)
    
    ! --- local ----------------------------------------

    integer         ::  topb_gaslaw
    integer         ::  l
    logical         ::  inftop
    real            ::  dh

    ! --- begin -----------------------------------------
    
    if ( pb(1) >= pb(lm+1) ) then
    
      !
      ! surface -> top
      !

      ! zero top pressure ? then gph is infinite there ...
      inftop = pb(lm+1) < tiny(1.0)

      ! maximum boundary for which gph can be computed from general gas law:
      if ( inftop ) then
        topb_gaslaw = lm
      else
        topb_gaslaw = lm+1
      end if

      ! initial height at bottom:
      gphb(1) = h0  ! m

      ! loop over all boundaries for which gph can be computed
      ! from general gas law:
      do l = 2, topb_gaslaw

        ! potential height of this layer:
        call PotentialHeight( pb(l), pb(l-1), T(l-1), Q(l-1), dh )

        ! geo potential height at this boundary:
        gphb(l) = gphb(l-1) + dh

      end do

      ! infinte top ? set height to double of mid dh:
      if ( inftop ) then
        call PotentialHeight( 0.5*pb(lm), pb(lm), T(lm), Q(lm), dh )
        gphb(lm+1) = gphb(lm) + 2*dh
      end if

    else
    
      !
      ! top -> surface
      !

      ! zero top pressure ? then gph is infinite there ...
      inftop = pb(1) < tiny(1.0)

      ! maximum boundary for which gph can be computed from general gas law:
      if ( inftop ) then
        topb_gaslaw = 2
      else
        topb_gaslaw = 1
      end if

      ! initial height at bottom:
      gphb(lm+1) = h0  ! m

      ! loop over all boundaries for which gph can be computed
      ! from general gas law:
      do l = lm, topb_gaslaw, -1

        ! potential height of this layer:
        call PotentialHeight( pb(l), pb(l+1), T(l), Q(l), dh )

        ! geo potential height at this boundary:
        gphb(l) = gphb(l+1) + dh

      end do

      ! infinte top ? set height to double of mid dh:
      if ( inftop ) then
        call PotentialHeight( 0.5*pb(2), pb(2), T(1), Q(1), dh )
        gphb(1) = gphb(2) + 2*dh
      end if

    end if
  end subroutine GeoPotentialHeightB
  
  
  ! ***

  
  ! GPH at pressure full levels
  
  subroutine GeoPotentialHeight( lm, pb, T, Q, h0, gph )
  
    ! --- in/out ---------------------------------------
    
    integer, intent(in)            ::  lm
    real, intent(in)               ::  pb(lm+1)   ! pressure bounds (Pa)
    real, intent(in)               ::  T(lm)      ! temperature (K)
    real, intent(in)               ::  Q(lm)      ! specific humidity (kg h2o/kg air)
    real, intent(in)               ::  h0         ! initial height  (m)
    real, intent(out)              ::  gph(lm)    ! geo.pot.height (m)
    
    ! --- local ----------------------------------------

    integer         ::  l
    real            ::  pf
    real            ::  dh

    ! --- begin -----------------------------------------
    
    if ( pb(1) >= pb(lm+1) ) then
    
      !
      ! surface -> top
      !

      ! initial height at bottom:
      pf = ( pb(1) + pb(2) )/2.0
      call PotentialHeight( pf, pb(1), T(1), Q(1), dh ) ! m
      gph(1) = h0 + dh
    
      ! loop over remaining full levels:   2,lm,1  or  lm-1,1,-1
      do l = 2, lm

        ! potential height of second half previous layer:
        call PotentialHeight( pb(l), pf, T(l-1), Q(l-1), dh )
        gph(l) = gph(l-1) + dh

        ! new full pressure:
        pf = ( pb(l) + pb(l+1) )/2.0

        ! potential height of first half of current layer:
        call PotentialHeight( pf, pb(l), T(l), Q(l), dh )
        gph(l) = gph(l) + dh

      end do
    
    else
    
      !
      ! top -> surface
      !

      ! initial height at bottom:
      pf = ( pb(lm+1) + pb(lm) )/2.0
      call PotentialHeight( pf, pb(lm+1), T(lm), Q(lm), dh ) ! m
      gph(lm) = h0 + dh
    
      ! loop over remaining full levels:   2,lm,1  or  lm-1,1,-1
      do l = lm-1, 1, -1

        ! potential height of second half previous layer:
        call PotentialHeight( pb(l+1), pf, T(l+1), Q(l+1), dh )
        gph(l) = gph(l+1) + dh

        ! new full pressure:
        pf = ( pb(l+1) + pb(l) )/2.0

        ! potential height of first half of current layer:
        call PotentialHeight( pf, pb(l+1), T(l), Q(l), dh )
        gph(l) = gph(l) + dh

      end do
    
    end if
    
  end subroutine GeoPotentialHeight
  
    

end module Phys_GPH


! ################################################################
!
!
!program test
!
!  use phys_gph
!
!  integer, parameter  ::  lm = 4
!  
!  real       ::  pb(lm+1)
!  real       ::  T(lm), Q(lm)
!  real       ::  gph(lm), gphb(lm+1)
!  
!  pb = (/ 1000.0, 500.0, 100.0, 50.0, 0.0 /)
!  T = (/ 293.0, 273.0, 260.0, 280.0 /)
!  Q = 0.0
!  
!  print *, ' '
!  print *, 'gph u'
!  call GeoPotentialHeight( lm, pb, T, Q, h0, gph )
!  print *, ' '
!  print *, 'gph u'
!  call GeoPotentialHeightB( lm, pb, T, Q, h0, gphb )
!
!  print *, ' '
!  do l = 1, lm
!    print *, gphb(l)
!    print *, '              ', gph(l)
!  end do
!  print *, gphb(lm+1)
!
!
!  pb = (/ 0.0, 50.0, 100.0, 500.0, 1000.0 /)
!  T = (/ 280.0, 260.0, 273.0, 293.0 /)
!  Q = 0.0
!  
!  print *, ' '
!  print *, 'gph d'
!  call GeoPotentialHeight( lm, pb, T, Q, h0, gph )
!  print *, ' '
!  print *, 'gphb d'
!  call GeoPotentialHeightB( lm, pb, T, Q, h0, gphb )
!
!  print *, ' '
!  do l = 1, lm
!    print *, gphb(l)
!    print *, '              ', gph(l)
!  end do
!  print *, gphb(lm+1)
!
!
!end program test
!
