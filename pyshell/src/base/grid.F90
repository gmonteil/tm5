!
!ProTeX: 1.15
!
!BOI
!
! !TITLE:        Grid  -  grid definitions and transformations
! !AUTHORS:      Arjo Segers
! !AFFILIATION:  KNMI
! !DATE:         21/04/2004
!
! !INTRODUCTION: Introduction
!
! \bv
!
! NAME
!   grid  -  grid definitions, interpolations, etc
!
! LAT/LON GRIDS
!
!   ! conect grid module:
!   use grid
!
!   ! declare grid definition and grid:
!   type(TllGridInfo)	 ::  lli
!   real, allocatable	 ::  ll(:,:)
!
!   ! define global grid of 2.5x2.5 degrees:
!   ! (mid point west, dlon, nlon, mid point south, dlat, nlat)
!   call Init( lli, -180.0, 2.5, 144, -90.0, 2.5, 72 )
!
!   ! storage of grid:
!   allocate( ll(lli%im,lli%jm) )
!
!   ! check if grid size is ok:
!   call Check( lli, ll )
!
!   ! devide or multiply all cells with their area:
!   call AreaOper( lli, ll, '/', 'rad2' )     ! [ll]/rad^2
!   call AreaOper( lli, ll, '/', 'm2'   )     ! [ll]/m^2
!   call AreaOper( lli, ll, '*', 'rad2' )     ! [ll] rad^2
!   call AreaOper( lli, ll, '*', 'm2'   )     ! [ll] m^2
!
!   !
!   ! Fill ll (defined by lli,nuv) with values from llX (defined by lliX,nuvX)
!   !
!   call FillGrid( lli, nuv, ll, lliX, nuvX, llX, combkey, status )
!   !
!   !   Key to identify data positions:
!   !      'n'  :  value valid for cell (center)            ll(1:nlon  ,1:nlat  )
!   !      'u'  :  value valid for east/west boundaries     ll(1:nlon+1,1:nlat  )
!   !      'v'  :  value valid for north/south boundaries   ll(1:nlon  ,1:nlat+1)
!   !
!   !   Coverage of lli by lliX :
!   !    o lliX is larger than or equal to lli   ->  all cells in ll changed
!   !    o lliX is smaller than lli              ->  only part of ll is changed
!   !
!   !   Create new ll from llX:
!   !    o llX is superset  ->  copy values from llX into ll
!   !    o llX is fine      ->  fill ll by combining cells in llX 
!   !                           (average/sum/etc given the combine key)
!   !
!   !   Combine keys: 'aver', 'sum'
!   !
!
!   ! clear memory:
!   deallocate( ll )
!   call Done( lli )
!
! HIARCHY
!
!   binas                    # general constants (pi etc)
!
!   singleton                # FFT routine
! 
!   grid_tools               # area of a rectangle etc
!
!     grid_type_ll           # define lat/lon grid
!     grid_type_gg           # define Gaussian grid
!     grid_type_sh           # define spherical harmonic grid
!
!        grid_interpol       # interpolate between grid types
!
!     grid_type_hyb          # hybride level definition
!
!        grid_3d             # transform 3D fields
!
!          grid              # main module, collects all grid stuff
!
!
! CHANGES
!   Arjo Segers
!
! \ev
!
!EOI
!

module grid

  use grid_tools

  use grid_type_ll
  use grid_type_gg
  use grid_type_sh
  use grid_interpol

  use grid_type_hyb
  
  use grid_3d

  implicit none
  
  public
  
end module grid

