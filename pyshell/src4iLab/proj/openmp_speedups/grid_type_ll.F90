!
! may 2002, Arjo Segers
!

! STAND ALONE VERSION FOR LAT/LON GRIDS

!
! Apply factor over region.
!
!   use grid, only : TRegion
!
!   type(TllRegion)   ::  llreg
!
!   ! define region:
!   call Init( llreg, west_deg, east_deg, south_deg, north_deg, status )
!
!   ! apply factor to data set within box, or to the complement;
!   ! if a grid cell only partly covers the region, the factor
!   ! is applied according to the area ratio;
!   ! x is at least 2d, if it is 3d the factor is applied for all levels:
!   call Region_Apply_Factor( lli, x, llreg, fac, status [,complement=.false.] )
!
!   ! done
!   call Done( llreg, status )
!


module grid_type_ll


  implicit none

  ! --- in/out --------------------------------

  private

  public  ::  TllGridInfo

  public  ::  Init, Done
  public  ::  Check

  public  ::  AreaOper
  public  ::  InterpolFractions
  public  ::  Interpol

  public  ::  BalanceMassFluxes
  public  ::  CheckMassBalance
  !public  ::  Poisson_init, Poisson_Done

  public  ::  GetRefinement
  public  ::  Match

  public  ::  FillGrid

  public  ::  EquivLat

  ! ~~ region

  public  ::  TllRegion
  public  ::  Region_Apply_Factor


  ! --- const ---------------------------------

  character(len=*), parameter  ::  mname = 'grid_type_ll'

  ! --- types ---------------------------------

  ! *** location, size, etc

  type TllGridInfo
    character(len=32)  ::  name
    ! * radius
    real               ::  R      ! m
    ! * spacing
    real               ::  dlon_deg, dlat_deg         ! degrees
    real               ::  dlon, dlat                 ! radians
    ! * size
    integer            ::  im, nlon
    integer            ::  jm, nlat
    ! * coordinates of gridpoint (cell center)
    !   indices 1, 2, ...
    real, pointer      ::  lon_deg(:), lat_deg(:)        ! degrees
    real, pointer      ::  lon(:), lat(:)                ! rad
    ! * coordinates of boundaries for cell around grid point;
    !   indices 0, 1, 2, ...
    real, pointer      ::  lon_bounds_deg(:,:), lat_bounds_deg(:,:)      ! degrees
    real, pointer      ::  lon_bounds(:,:), lat_bounds(:,:)              ! rad
    ! * boundaries in a rank-2 array:
    real, pointer      ::  blon_deg(:), blat_deg(:)      ! degrees
    real, pointer      ::  blon(:), blat(:)              ! rad
    ! * area for cell in certain row
    real, pointer      ::  area(:)                       ! rad^2
    real, pointer      ::  area_m2(:)                    ! m^2
    ! * cell length in m
    real, pointer      ::  dx(:), bdx(:)                 ! m
    real               ::  dy, bdy                       ! m
  end type TllGridInfo


  ! ~~ region


  type TllRegion
    ! region boundaries in degrees:
    real    ::  west_deg, east_deg, south_deg, north_deg
    ! idem in radians:
    real    ::  west, east, south, north
  end type TllRegion

  ! ~~ something called 'fac' needed for SolvePoissonEq_Zoom
  !type Tpoisson_fac
    !real, allocatable   :: fac(:,:)
  !end type Tpoisson_fac

  !type(Tpoisson_fac), allocatable   :: poisson_fac(:)

  ! --- interfaces ----------------------------

  interface Init
    module procedure llgridinfo_Init
    module procedure llreg_Init
    module procedure llgrid_PCTM_Init
  end interface

  interface Done
    module procedure llgridinfo_Done
    module procedure llreg_Done
  end interface

  interface Check
    module procedure llgrid_Check_i
    module procedure llgrid_Check_r
  end interface

  interface AreaOper
    module procedure llgrid_AreaOper_2d
    module procedure llgrid_AreaOper_3d
  end interface

  interface InterpolFractions
    module procedure llgrid_InterpolFractions
  end interface

  interface Interpol
    module procedure llgrid_Eval_2d
    module procedure llgrid_Eval_3d
  end interface

  interface Match
    module procedure llgrid_Match
  end interface

  interface BalanceMassFluxes
    module procedure BalanceMassFluxes_sm
    module procedure BalanceMassFluxes_m
  end interface

  interface EquivLat
    module procedure llgrid_EquivLat
    module procedure llgrid_EquivLat_sort
  end interface

  ! ~~ region

  interface Region_Apply_Factor
    module procedure llreg_Region_Apply_Factor_2d
    module procedure llreg_Region_Apply_Factor_3d
  end interface



contains


  ! ========================================================
  ! ===
  ! === Init, Done
  ! ===
  ! ========================================================


  !
  !  blat(j) +-------+
  !          |   |   |
  !   lat(j) |---o---|
  !          |   |   |
  !          +-------+
  !            lon(i)
  !               blon(i)

  !subroutine Poisson_init

    !use dims,   only : nregions, im, jm
    !use binas,  only : pi

    !implicit none

    !integer         :: region, i, j

    !allocate(poisson_fac(nregions))
    !do region = 1, nregions
        !allocate(poisson_fac(region)%fac(im(region), jm(region)))
        !do j = 1, jm(region)
            !do i = 1, im(region)
                !poisson_fac(region)%fac(i,j) = 2.0*(cos(2*pi*(i-1)/im(region))+cos(2*pi*(j-1)/jm(region))-2.)
            !end do
        !end do
        !poisson_fac(region)%fac(1,1) = 1.0   !to avoid division by zero
    !end do

  !end subroutine Poisson_init

  !subroutine Poisson_Done

    !implicit none

    !deallocate(poisson_fac)

  !end subroutine Poisson_Done

  subroutine llgridinfo_Init( lli, west_deg,  dlon_deg, im, &
                                   south_deg, dlat_deg, jm, status, name  )

    use Grid_Tools, only : ll_area
    use grid_tools, only : deg2rad, ae

    ! --- in/out ---------------------------------

    type(TllGridInfo), intent(out)     ::  lli
    real, intent(in)                   ::  west_deg, dlon_deg
    integer, intent(in)                ::  im
    real, intent(in)                   ::  south_deg, dlat_deg
    integer, intent(in)                ::  jm
    integer, intent(out)               ::  status

    character(len=*), intent(in), optional  ::  name

    ! --- local ---------------------------------

    integer          ::  i, j

    ! --- begin ---------------------------------

    ! store name ?
    if ( present(name) ) then
      lli%name = name
    else
      lli%name = 'anonymous'
    end if

    ! *** radius

    lli%R = ae

    ! *** grid spacing

    lli%dlon_deg = dlon_deg
    lli%dlon = dlon_deg*deg2rad

    lli%dlat_deg = dlat_deg
    lli%dlat = dlat_deg*deg2rad

    ! *** grid range

    lli%im = im
    lli%nlon = im

    lli%jm = jm
    lli%nlat = jm

    ! *** grid points

    ! * coor of grid points

    ! east-west
    allocate( lli%lon_deg(im) )
    do i = 1, im
      lli%lon_deg(i) = west_deg + (i-1)*dlon_deg
    end do
    allocate( lli%lon(im) )
    lli%lon = lli%lon_deg * deg2rad

    ! south-north
    allocate( lli%lat_deg(jm) )
    do j = 1, jm
      lli%lat_deg(j) = south_deg + (j-1)*dlat_deg
    end do
    allocate( lli%lat(jm) )
    lli%lat = lli%lat_deg * deg2rad

    ! *** cells with grid point in center;
    !     grid point at pole is at top of triangle !

    ! * bounds

    ! west-east
    allocate( lli%blon_deg(0:im) )
    do i = 0, im
      lli%blon_deg(i) = west_deg + (i-0.5)*dlon_deg
    end do
    allocate( lli%blon(0:im) )
    lli%blon = lli%blon_deg * deg2rad

    ! south-north
    allocate( lli%blat_deg(0:jm) )
    do j = 0, jm
      lli%blat_deg(j) = south_deg + (j-0.5)*dlat_deg
    end do
    if ( lli%blat_deg(0)  < -90.0 ) lli%blat_deg(0)  = -90.0
    if ( lli%blat_deg(0)  >  90.0 ) lli%blat_deg(0)  =  90.0
    if ( lli%blat_deg(jm) < -90.0 ) lli%blat_deg(jm) = -90.0
    if ( lli%blat_deg(jm) >  90.0 ) lli%blat_deg(jm) =  90.0
    allocate( lli%blat(0:jm) )
    lli%blat = lli%blat_deg * deg2rad

    ! * bounds in a rank-2 array:

    ! west-east
    allocate( lli%lon_bounds_deg(2,im) )
    lli%lon_bounds_deg(1,:) = lli%blon_deg(0:im-1)
    lli%lon_bounds_deg(2,:) = lli%blon_deg(1:im)
    allocate( lli%lon_bounds(2,im) )
    lli%lon_bounds = lli%lon_bounds_deg * deg2rad

    ! south-north
    allocate( lli%lat_bounds_deg(2,jm) )
    lli%lat_bounds_deg(1,:) = lli%blat_deg(0:jm-1)
    lli%lat_bounds_deg(2,:) = lli%blat_deg(1:jm)
    allocate( lli%lat_bounds(2,jm) )
    lli%lat_bounds = lli%lat_bounds_deg * deg2rad

    ! * area of cell in lat band
    ! rad^2 :
    allocate( lli%area(jm) )
    do j = 1, jm
      lli%area(j) = ll_area( 0.0, lli%dlon, lli%blat(j-1), lli%blat(j) )
    end do
    ! m^2 :
    allocate( lli%area_m2(jm) )
    lli%area_m2 = lli%area * lli%R**2

    ! * length in m
    ! east-west, mid latitude of cell
    allocate( lli%dx(jm) )
    lli%dx = lli%dlon * lli%R * cos(lli%lat)
    ! east-west, boundaries
    allocate( lli%bdx(0:jm) )
    lli%bdx = lli%dlon * lli%R * cos(lli%blat)
    ! north-south, the same for each longitude
    lli%dy  = lli%dlat * lli%R
    lli%bdy = lli%dlat * lli%R

    ! ok
    status = 0

  end subroutine llgridinfo_Init


  ! subroutine for defining a MERRA grid, useful if we want to work with PCTM/SiB

  subroutine llgrid_PCTM_Init(lli, im, jm, status, name)

    use Grid_Tools, only : ll_area
    use grid_tools, only : deg2rad, ae

    ! --- in/out ---------------------------------

    type(TllGridInfo), intent(out)     ::  lli
    integer, intent(in)                ::  im, jm
    integer, intent(out)               ::  status
    character(len=*), intent(in), optional  ::  name

    ! --- constant ------------------------------

    character(len=*), parameter :: rname = mname//'/llgrid_PCTM_Init'

    ! --- local ---------------------------------

    integer          ::  i, j

    ! --- begin ---------------------------------

    ! store name ?
    if ( present(name) ) then
      lli%name = name
    else
      lli%name = 'anonymous'
    end if

    ! *** radius

    lli%R = ae

    ! *** grid spacing, should be set to nonsense in the latitudinal direction because of the half degrees at the poles

    lli%dlon_deg = 360.0/im
    lli%dlon = lli%dlon_deg*deg2rad

    lli%dlat_deg = -999999.0
    lli%dlat = -999999.0

    ! *** grid range

    lli%im = im
    lli%nlon = im

    lli%jm = jm
    lli%nlat = jm

    ! We only do 288x181 grids for now, raise an error if this is not so
    if (im /= 288 .or. jm /= 181) then
        write(*,'(a, ": You requested a ", i3, " x ", i3, " PCTM grid")') rname, im, jm
        write(*,'(a, ": I only know how to form a 288 x 181 PCTM grid")') rname
        status = 1
        return
    end if

    ! *** grid points

    ! * bounds

    ! west-east
    allocate( lli%blon_deg(0:im) )
    ! the longitude boundaries are -180.625, -179.375, -178.125, ... 179.375
    do i = 0, im
      lli%blon_deg(i) = -180.625 + i*1.25
    end do
    allocate( lli%blon(0:im) )
    lli%blon = lli%blon_deg * deg2rad

    ! south-north
    allocate( lli%blat_deg(0:jm) )
    ! the latitude boundaries are -90.0, -89.5, -88.5, -87.5, ... 88.5, 89.5, 90.0
    lli%blat_deg(0) = -90.0
    lli%blat_deg(jm) = 90.0
    do j = 1, jm-1
      lli%blat_deg(j) = -89.5 + (j-1)*1.0
    end do
    allocate( lli%blat(0:jm) )
    lli%blat = lli%blat_deg * deg2rad

    ! * coor of grid points

    ! east-west
    allocate( lli%lon_deg(im) )
    do i = 1, im
      lli%lon_deg(i) = 0.5*(lli%blon_deg(i-1) + lli%blon_deg(i))
    end do
    allocate( lli%lon(im) )
    lli%lon = lli%lon_deg * deg2rad

    ! south-north
    allocate( lli%lat_deg(jm) )
    do j = 1, jm
      lli%lat_deg(j) = 0.5*(lli%blat_deg(j-1) + lli%blat_deg(j))
    end do
    allocate( lli%lat(jm) )
    lli%lat = lli%lat_deg * deg2rad

    ! *** cells with grid point in center;

    ! * bounds in a rank-2 array:

    ! west-east
    allocate( lli%lon_bounds_deg(2,im) )
    lli%lon_bounds_deg(1,:) = lli%blon_deg(0:im-1)
    lli%lon_bounds_deg(2,:) = lli%blon_deg(1:im)
    allocate( lli%lon_bounds(2,im) )
    lli%lon_bounds = lli%lon_bounds_deg * deg2rad

    ! south-north
    allocate( lli%lat_bounds_deg(2,jm) )
    lli%lat_bounds_deg(1,:) = lli%blat_deg(0:jm-1)
    lli%lat_bounds_deg(2,:) = lli%blat_deg(1:jm)
    allocate( lli%lat_bounds(2,jm) )
    lli%lat_bounds = lli%lat_bounds_deg * deg2rad

    ! * area of cell in lat band
    ! rad^2 :
    allocate( lli%area(jm) )
    do j = 1, jm
      ! ll_area( west, east, south, north ), in radians
      lli%area(j) = ll_area( 0.0, lli%dlon, lli%blat(j-1), lli%blat(j) )
    end do
    ! m^2 :
    allocate( lli%area_m2(jm) )
    lli%area_m2 = lli%area * lli%R**2

    ! * length in m
    ! east-west, mid latitude of cell
    allocate( lli%dx(jm) )
    lli%dx = lli%dlon * lli%R * cos(lli%lat)
    ! east-west, boundaries
    allocate( lli%bdx(0:jm) )
    lli%bdx = lli%dlon * lli%R * cos(lli%blat)
    ! north-south, the same for each longitude
    ! does not hold for this grid, because lli%dlat is -999999
    lli%dy  = lli%dlat * lli%R
    lli%bdy = lli%dlat * lli%R

    ! ok
    status = 0

  end subroutine llgrid_PCTM_Init

  ! ===


  subroutine llgridinfo_Done( lli, status )

    ! --- in/out ---------------------------------

    type(TllGridInfo), intent(inout)   ::  lli
    integer, intent(out)               ::  status

    ! --- begin ---------------------------------

    ! free memory
    deallocate( lli%lon_deg )
    deallocate( lli%lat_deg )
    deallocate( lli%lon )
    deallocate( lli%lat )

    deallocate( lli%blon_deg )
    deallocate( lli%blat_deg )
    deallocate( lli%blon )
    deallocate( lli%blat )

    deallocate( lli%lon_bounds_deg )
    deallocate( lli%lat_bounds_deg )
    deallocate( lli%lon_bounds )
    deallocate( lli%lat_bounds )

    deallocate( lli%area )
    deallocate( lli%area_m2 )

    deallocate( lli%dx  )
    deallocate( lli%bdx )

    ! ok
    status = 0

  end subroutine llgridinfo_Done


  ! =============================================================


  subroutine llgrid_Check_i( lli, nuv, ll, status )

    ! --- in/out ----------------------------------

    type(TllGridInfo), intent(in)  ::  lli
    character(len=*), intent(in)   ::  nuv
    integer, intent(in)            ::  ll(:,:)
    integer, intent(out)           ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_Check_i'

    ! --- begin ----------------------------------

    ! check shape of target grid:
    if ( ((nuv == 'n') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(ll,1) /= lli%nlon+1) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat+1))) ) then
      write (*,'("ERROR - target array does not match with grid definition:")')
      write (*,'("ERROR -   lli    : ",i3," x ",i3)') lli%nlon, lli%nlat
      write (*,'("ERROR -   nuv    : ",a          )') nuv
      write (*,'("ERROR -   ll     : ",i3," x ",i3)') shape(ll)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine llgrid_Check_i


  ! ***


  subroutine llgrid_Check_r( lli, nuv, ll, status )

    ! --- in/out ----------------------------------

    type(TllGridInfo), intent(in)  ::  lli
    character(len=*), intent(in)   ::  nuv
    real, intent(in)               ::  ll(:,:)
    integer, intent(out)           ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_Check_r'

    ! --- begin ----------------------------------

    ! check shape of target grid:
    if ( ((nuv == 'n') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(ll,1) /= lli%nlon+1) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat+1))) ) then
      write (*,'("ERROR - target array does not match with grid definition:")')
      write (*,'("ERROR -   lli    : ",i3," x ",i3)') lli%nlon, lli%nlat
      write (*,'("ERROR -   nuv    : ",a          )') nuv
      write (*,'("ERROR -   ll     : ",i3," x ",i3)') shape(ll)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine llgrid_Check_r


  ! =====================================================

  !  call AreaOper( lli, ll, '/' | '*' | '=', 'rad2' | 'm2' )

  subroutine llgrid_AreaOper_2d( lli, ll, oper, unit, status )

    ! --- in/out ----------------------------------

    type(TllGridInfo), intent(in)           ::  lli
    real, intent(inout)                     ::  ll(:,:)
    character(len=*), intent(in)            ::  unit, oper
    integer, intent(out)                    ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_AreaOper_2d'

    ! --- local --------------------------------

    integer   ::  j
    real      ::  cell_area

    ! --- begin ----------------------------------

    ! check ...
    if ( size(ll,2) /= lli%nlat ) then
      write (*,'("ERROR - unexpected size of ll grid:")')
      write (*,'("ERROR -   shape(ll) : ",i4," x ",i4)') shape(ll)
      write (*,'("ERROR -   lli%nlat  : ",i4)') lli%nlat
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! loop over latitudes:
    do j = 1, lli%nlat

      ! select correct area for cells in this row:
      select case ( unit )
        case ( 'rad2' )
          cell_area = lli%area(j)
        case ( 'm2' )
          cell_area = lli%area_m2(j)
        case default
          write (*,'("ERROR - unknown unit : ",a)') trim(unit)
          write (*,'("ERROR in ",a)') rname; status=1; return
      end select

      ! assign/mult/div by cell area:
      select case ( oper )
        case ( '=' )
          ll(:,j) = cell_area
        case ( '/' )
          ll(:,j) = ll(:,j) / cell_area
        case ( '*' )
          ll(:,j) = ll(:,j) * cell_area
        case default
          write (*,'("ERROR - unknown operation : ",a)') trim(oper)
          write (*,'("ERROR in ",a)') rname; status=1; return
      end select

    end do

    ! ok
    status = 0

  end subroutine llgrid_AreaOper_2d


  ! *


  subroutine llgrid_AreaOper_3d( lli, ll, oper, unit, status )

    ! --- in/out ----------------------------------

    type(TllGridInfo), intent(in)           ::  lli
    real, intent(inout)                     ::  ll(:,:,:)
    character(len=*), intent(in)            ::  oper
    character(len=*), intent(in)            ::  unit
    integer, intent(out)                    ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_AreaOper_3d'

    ! --- local --------------------------------

    integer   ::  l

    ! --- begin ----------------------------------

    ! loop over layers
    do l = 1, size(ll,3)
      ! apply 2d operator:
      call AreaOper( lli, ll(:,:,l), oper, unit, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end do

    ! ok
    status = 0

  end subroutine llgrid_AreaOper_3d


  ! =====================================

  !
  ! Interpolation to (lon,lat) in deg.
  !

  subroutine llgrid_InterpolFractions( lli, lon, lat, ii, jj, ff )

    ! --- in/out ---------------------------

    type(TllGridInfo), intent(in)    ::  lli
    real, intent(in)                 ::  lon, lat    ! deg
    integer, intent(out)             ::  ii(4)
    integer, intent(out)             ::  jj(4)
    real, intent(out)                ::  ff(4)

    ! --- local -----------------------------

    real          ::  lonX, latX
    real          ::  ir
    integer       ::  i1, i2
    real          ::  i1f, i2f
    real          ::  jr
    integer       ::  j1, j2
    real          ::  j1f, j2f

    ! --- begin -----------------------------

    ! bring  lon  in  [-180,180.0]
    lonX = modulo(lon,360.0)
    if ( lonX > 180.0 ) lonX = lonX - 360.0

    ! check lat
    latX = lat
    if ( latX < -90.0 .or. latX > 90.0 ) then
      write (*,'("ERROR - invalid lat (deg) :",f12.4)') latX
      write (*,'("ERROR in ",a)') 'llgrid_InterpolFractions'; stop
    end if

    !
    !    1  2
    !    3  4
    !

    ! i fractions ; circular
    !
    !               *
    !        [----+----][----+----] .. [----+----]
    !  ir        1.0        2.0           120.0
    !
    ir = ( lonX - lli%lon_deg(1) ) / lli%dlon_deg + 1.0
    i1 = floor(ir)
    i2 = i1 + 1
    !
    i2f = ( ir - i1 ) / ( i2 - i1 )
    i1f = 1.0 - i2f
    !
    if ( i1 < 1        ) i1 = lli%nlon
    if ( i2 > lli%nlon ) i2 = 1

    ! j fractions ; constant in half cells at poles
    jr = ( latX - lli%lat_deg(1) ) / lli%dlat_deg + 1.0
    j1 = floor(jr)
    j2 = j1 + 1
    !
    if ( j1 < 1 ) then
      j2f = ( jr - 0.5 ) / ( j2 - 0.5 )
      j1f = 1.0 - j2f
    else if ( j2 > lli%nlat ) then
      j2f = ( jr - j1 ) / ( lli%nlat+0.5 - j1 )
      j1f = 1.0 - j2f
    else
      j2f = ( jr - j1 ) / ( j2 - j1 )
      j1f = 1.0 - j2f
    end if

    ! fill output
    !
    !    j2   3    4
    !    j1   1    2
    !
    !        i1   i2
    !
    ii = (/ i1,      i2,      i1,      i2      /)
    jj = (/     j1,      j1,      j2,      j2  /)
    ff = (/ i1f*j1f, i2f*j1f, i1f*j2f, i2f*j2f /)

!    write (*,'("   ",2i5)') ii(3), ii(4)
!    write (*,'(i4," ",f4.2," ",f4.2," ",i4)') jj(3),ff(3),ff(4),jj(4)
!    write (*,'(i4," ",f4.2," ",f4.2," ",i4)') jj(1),ff(1),ff(2),jj(2)
!    write (*,'("   ",2i5)') ii(1), ii(2)

  end subroutine llgrid_InterpolFractions


  ! ***


  subroutine llgrid_Eval_2d( lli, ll, lon, lat, res )

    ! --- in/out ---------------------------

    type(TllGridInfo), intent(in)    ::  lli
    real, intent(in)                 ::  ll(:,:)
    real, intent(in)                 ::  lon, lat    ! deg
    real, intent(out)                ::  res

    integer ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_Eval_2d'

    ! --- local -----------------------------

    integer        ::  ii(4)
    integer        ::  jj(4)
    real           ::  ff(4)

    integer        ::  k
    real           ::  value

    ! --- begin -----------------------------

    ! check ...
    call Check( lli, 'n', ll, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if

    ! indices and fractions
    call InterpolFractions( lli, lon, lat, ii, jj, ff )

    ! init zero
    res = 0.0

    ! add contributions:
    do k = 1, 4

      ! handle poles
      if ( jj(k) < 1 ) then
        value = sum(ll(:,1))/lli%nlon
      else if ( jj(k) > lli%nlat ) then
        value = sum(ll(:,lli%nlat))/lli%nlon
      else
        value = ll(ii(k),jj(k))
      end if

      ! add fraction
      res = res + value * ff(k)

    end do

  end subroutine llgrid_Eval_2d


  ! ***


  subroutine llgrid_Eval_3d( lli, ll, lon, lat, res )

    ! --- in/out ---------------------------

    type(TllGridInfo), intent(in)    ::  lli
    real, intent(in)                 ::  ll(:,:,:)
    real, intent(in)                 ::  lon, lat    ! deg
    real, intent(out)                ::  res(size(ll,3))

    integer        ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_Eval_3d'

    ! --- local -----------------------------

    integer        ::  ii(4)
    integer        ::  jj(4)
    real           ::  ff(4)

    integer        ::  k
    real           ::  value(size(ll,3))

    ! --- begin -----------------------------

    ! check ...
    call Check( lli, 'n', ll(:,:,1), status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if

    ! indices and fractions
    call InterpolFractions( lli, lon, lat, ii, jj, ff )

    ! init zero
    res = 0.0

    ! add contributions:
    do k = 1, 4

      ! handle poles
      if ( jj(k) < 1 ) then
        value = sum(ll(:,1,:),1) / lli%nlon
      else if ( jj(k) > lli%nlat ) then
        value = sum(ll(:,lli%nlat,:),1) / lli%nlon
      else
        value = ll(ii(k),jj(k),:)
      end if

      ! add fraction
      res = res + value * ff(k)

    end do

  end subroutine llgrid_Eval_3d


  ! =================================================================
  ! ===
  ! ===  match fine grid with coarse grid
  ! ===
  ! =================================================================


  subroutine GetRefinement( cgi, fgi, &
                              refine_i, refine_j, &
                              cg_i1, cg_i2, cg_j1, cg_j2, status )

    ! --- in/out ------------------------------

    type(TllGridInfo), intent(in)     ::  cgi
    type(TllGridInfo), intent(in)     ::  fgi
    integer, intent(out)              ::  refine_i, refine_j
    integer, intent(out)              ::  cg_i1, cg_i2, cg_j1, cg_j2
    integer, intent(out)              ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/GetRefinement'

    ! --- local -----------------------------

    integer         ::  i, j

    ! --- begin -----------------------------

    ! *** determine refinement

    refine_i = nint( cgi%dlon_deg / fgi%dlon_deg )
    refine_j = nint( cgi%dlat_deg / fgi%dlat_deg )


    ! *** position of fine grid within coarse grid:

    ! search column in coarse grid with same west bound as fine grid:
    cg_i1 = 0
    do i = 1, cgi%im
      if ( cgi%blon_deg(i-1) == fgi%blon_deg(0) ) then
        cg_i1 = i
        exit
      end if
    end do
    if ( cg_i1 < 1 ) then
      write (*,'("ERROR - could not match west bound of fine grid:")')
      write (*,'("ERROR -   cgi%blon : ",f12.4)') cgi%blon_deg
      write (*,'("ERROR -   target   : ",f12.4)') fgi%blon_deg(0)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! search column in coarse grid with same east bound as fine grid:
    cg_i2 = 0
    do i = 1, cgi%im
      if ( cgi%blon_deg(i) == fgi%blon_deg(fgi%im) ) then
        cg_i2 = i
        exit
      end if
    end do
    if ( cg_i2 < 1 ) then
      write (*,'("ERROR - could not match east bound of fine grid")')
      write (*,'("ERROR -   cgi%blon : ",f12.4)') cgi%blon_deg
      write (*,'("ERROR -   target   : ",f12.4)') fgi%blon_deg(fgi%im)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! check ...
    if ( (cg_i2-cg_i1+1)*refine_i /= fgi%im ) then
      write (*,'("ERROR - i refinement not ok:")')
      write (*,'("ERROR -   coarse cells : ",f12.4)') cg_i1, cg_i2
      write (*,'("ERROR -   refinement   : ",f12.4)') refine_i
      write (*,'("ERROR -   fine cells   : ",f12.4)') fgi%im
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! search row in coarse grid with same south bound as fine grid:
    cg_j1 = 0
    do j = 1, cgi%jm
      if ( cgi%blat_deg(j-1) == fgi%blat_deg(0) ) then
        cg_j1 = j
        exit
      end if
    end do
    if ( cg_j1 < 1 ) then
      write (*,'("ERROR - could not match south bound of fine grid")')
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! search row in coarse grid with same south bound as fine grid:
    cg_j2 = 0
    do j = 1, cgi%jm
      if ( cgi%blat_deg(j) == fgi%blat_deg(fgi%jm) ) then
        cg_j2 = j
        exit
      end if
    end do
    if ( cg_j2 < 1 ) then
      write (*,'("ERROR - could not match north bound of fine grid")')
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! check ...
    if ( (cg_j2-cg_j1+1)*refine_j /= fgi%jm ) then
      write (*,'("ERROR - j refinement not ok:")')
      write (*,'("ERROR -   coarse cells : ",2i5)') cg_j1, cg_j2
      write (*,'("ERROR -   refinement   : ",i5)') refine_j
      write (*,'("ERROR -   fine cells   : ",i5)') fgi%jm
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine GetRefinement


  ! ***

  ! Relate fine grid and coarse grid:
  !
  !   call Match_cell( action, key, cgi, cg, fgi, fg )
  !
  !  1. The coarse grid covers the same area as the fine grid:
  !
  !     action:
  !       'combine'   : fill coarse grid by combining cells of fine grid
  !     key:
  !       'sum'       : sum values in fine grid
  !       'aver'      : aver values in fine grid
  !       'area-aver' : idem, weighted with area
  !
  !  2. The fine grid is a subset of the coarse grid:
  !
  !     action:
  !       'subset'    :  fill fine grid as a subset of coarse grid
  !     key:
  !       not used
  !
  !  3. The fine grid is a zooming area of the coarse grid:
  !
  !     action:
  !       'match'     : adjust values in fine grid to match coarse grid
  !     key:
  !       'sum'       : sum values in fine grid
  !       'aver'      : aver values in fine grid
  !       'area-aver' : idem, weighted with area
  !

  subroutine llgrid_Match( key, nuv, pgi, pg, tgi, tg, status )

    ! --- in/out ------------------------------

!    character(len=*), intent(in)      ::  action
    character(len=*), intent(in)      ::  key
    character(len=1), intent(in)      ::  nuv
    type(TllGridInfo), intent(in)     ::  pgi
    real, intent(in)                  ::  pg(pgi%im,pgi%jm)
    type(TllGridInfo), intent(in)     ::  tgi
    real, intent(inout)               ::  tg(tgi%im,tgi%jm)
    integer, intent(out)              ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_Match'

    ! --- begin -----------------------------

    ! call nuv specific routine
    select case ( nuv )
      case ( 'n' )
        call Match_cell( key, pgi, pg, tgi, tg, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      case ( 'u' )
        call Match_u( key, pgi, pg, tgi, tg, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      case ( 'v' )
        call Match_v( key, pgi, pg, tgi, tg, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      case default
        write (*,'("ERROR - unsupported nuv `",a,"`")') nuv
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select

    ! ok
    status = 0

  end subroutine llgrid_Match


  ! ***


  subroutine Match_cell( key, pgi, pg, tgi, tg, status )

    ! --- in/out ------------------------------

!    character(len=*), intent(in)      ::  action
    character(len=*), intent(in)      ::  key
    type(TllGridInfo), intent(in)     ::  pgi
    real, intent(in)                  ::  pg(pgi%im,pgi%jm)
    type(TllGridInfo), intent(in)     ::  tgi
    real, intent(inout)               ::  tg(tgi%im,tgi%jm)
    integer, intent(out)              ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/Match_cell'

    ! --- local -----------------------------

    integer      ::  refine_i, refine_j
    integer      ::  cg_i1, cg_i2, cg_j1, cg_j2
    integer      ::  ci, cj
    integer      ::  fg_i1, fg_i2, fg_j1, fg_j2
    integer      ::  fi, fj
    real         ::  fsum

    ! --- begin -----------------------------

!    select case ( action )
!
!      !
!      ! match fine grid with coarse cells:
!      !
!      case ( 'match' )

        ! determine refinement
        call GetRefinement( pgi, tgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! loop over cells in coarse grid covering fine grid:
        do cj = cg_j1, cg_j2
          do ci = cg_i1, cg_i2

            fg_i1 = (ci-cg_i1)*refine_i + 1  ;  fg_i2 = fg_i1-1 + refine_i
            fg_j1 = (cj-cg_j1)*refine_j + 1  ;  fg_j2 = fg_j1-1 + refine_j

            ! sum over cells in fine grid:
            select case ( key )
              case ( 'sum' )
                fsum = sum( tg(fg_i1:fg_i2,fg_j1:fg_j2) )
                ! distribute difference equally over all cells in fine grid:
                !                                     (/6,4/)             + (  14     - 10 )/  2      = (/8,6/)
                tg(fg_i1:fg_i2,fg_j1:fg_j2) = tg(fg_i1:fg_i2,fg_j1:fg_j2) + (pg(ci,cj)-fsum)/(refine_j*refine_i)
                ! cmk  corrected: divide by (refine_i * refine_j)
              case ( 'aver' )
                fsum = sum( tg(fg_i1:fg_i2,fg_j1:fg_j2) )/(refine_i*refine_j)
                ! add difference in averages to all cells in fine grid:
                !                                     (/6,4/)             + (   7     -  5 )          = (/8,6/)
                tg(fg_i1:fg_i2,fg_j1:fg_j2) = tg(fg_i1:fg_i2,fg_j1:fg_j2) + (pg(ci,cj)-fsum)
              case ( 'area-aver' )
                fsum = 0.0
                do fj = fg_j1, fg_j2
                  fsum = fsum + sum(tg(fg_i1:fg_i2,fj))*tgi%area(fj)
                end do
                fsum = fsum / pgi%area(cj)
                ! add difference in averages to all cells in fine grid:
                !                                     (/6,4/)             + (   7     -  5 )          = (/8,6/)
                tg(fg_i1:fg_i2,fg_j1:fg_j2) = tg(fg_i1:fg_i2,fg_j1:fg_j2) + (pg(ci,cj)-fsum)
              case default
                write (*,'("ERROR - unsupported key for match action:")')
!                write (*,'("ERROR -   action : ",a)') trim(action)
                write (*,'("ERROR -   key    : ",a)') trim(key)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

!            ! match fine grid with coarse cell:
!            if ( fsum == 0.0 ) then
!              if ( abs(pg(ci,cj)) > 1.0e-5 ) then
!                write (*,'("ERROR - zero sum over coarse cell (",i3,",",i3,")")') ci, cj
!                write (*,'("ERROR -   coarse cell : ",es12.4)') pg(ci,cj)
!                write (*,'("ERROR -   refine      : ",2i6)') refine_i, refine_j
!                write (*,'("ERROR -   fine grid   : ",es12.4)') tg(fg_i1:fg_i2,fg_j1:fg_j2)
!                write (*,'("ERROR in ",a)') rname; status=1; return
!              end if
!            else
!              ! scale towards value of coarse grid:
!              tg(fg_i1:fg_i2,fg_j1:fg_j2) = tg(fg_i1:fg_i2,fg_j1:fg_j2) * pg(ci,cj)/fsum
!            end if

          end do
        end do

!      !
!      ! fine grid is subset of coarse grid:
!      !
!      case ( 'subset' )
!
!        ! determine refinement
!        call GetRefinement( pgi, tgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2 )
!
!        if ( (refine_i /= 1) .or. (refine_j /= 1) ) then
!          write (*,'("ERROR - for this action the fine grid should a subset of a coarse grid:")')
!          write (*,'("ERROR -   action     : ",a)') trim(action)
!          write (*,'("ERROR -   refinement : ",2i6)') refine_i, refine_j
!          write (*,'("ERROR in ",a)') rname; status=1; return
!        end if
!
!        ! loop over cells in coarse grid covering fine grid:
!        do cj = cg_j1, cg_j2
!          do ci = cg_i1, cg_i2
!
!            fg_i1 = (ci-cg_i1) + 1
!            fg_j1 = (cj-cg_j1) + 1
!
!            ! fine grid is subset of coarse grid
!            tg(fg_i1,fg_j1) = pg(ci,cj)
!
!          end do
!        end do
!
!      !
!      ! collect cells in fine grid to coarse grid
!      !
!      case ( 'combine' )
!
!        ! determine refinement
!        ! NOTE: parent grid is fine now, target is coarse !
!        call GetRefinement( tgi, pgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2 )
!
!        if ( (cg_i1 /= 1) .or. (cg_i2 /= tgi%im) .or. &
!             (cg_j1 /= 1) .or. (cg_j2 /= tgi%jm) ) then
!          write (*,'("ERROR - for this action the fine grid should cover the complete coarse grid:")')
!          write (*,'("ERROR -   action  : ",a)') trim(action)
!          write (*,'("ERROR -   covered : [",i4,",",i4,"] x [",i4,",","]")') cg_i1,cg_i2,cg_j1,cg_j2
!          write (*,'("ERROR -   of      : ",2i5)') tgi%im, tgi%jm
!          write (*,'("ERROR in ",a)') rname; status=1; return
!        end if
!
!        ! loop over cells in coarse grid covering fine grid:
!        do cj = cg_j1, cg_j2
!          do ci = cg_i1, cg_i2
!
!            fg_i1 = (ci-cg_i1)*refine_i + 1  ;  fg_i2 = fg_i1-1 + refine_i
!            fg_j1 = (cj-cg_j1)*refine_j + 1  ;  fg_j2 = fg_j1-1 + refine_j
!
!            ! sum over cells in fine grid:
!            select case ( key )
!              case ( 'sum' )
!                fsum = sum( pg(fg_i1:fg_i2,fg_j1:fg_j2) )
!              case ( 'aver' )
!                fsum = sum( pg(fg_i1:fg_i2,fg_j1:fg_j2) )/(refine_i*refine_j)
!              case ( 'area-aver' )
!                fsum = 0.0
!                do fj = fg_j1, fg_j2
!                  fsum = fsum + sum(pg(fg_i1:fg_i2,fj))*pgi%area(fj)
!                end do
!                fsum = fsum / tgi%area(cj)
!              case default
!                write (*,'("ERROR - unsupported key for match action:")')
!                write (*,'("ERROR -   action : ",a)') trim(action)
!                write (*,'("ERROR -   key    : ",a)') trim(key)
!                write (*,'("ERROR in ",a)') rname; status=1; return
!            end select
!
!            ! collect cells in fine parent grid to target coarse grid
!            tg(ci,cj) = fsum
!
!          end do
!        end do
!
!      case default
!        write (*,'("ERROR - unknown action `",a,"`")') trim(action)
!        write (*,'("ERROR in ",a)') rname; status=1; return
!    end select

    ! ok
    status = 0

  end subroutine Match_cell


  ! ***

  ! flux through east/west boundaries

  subroutine Match_u( key, pgi, pg, tgi, tg, status )

    ! --- in/out ------------------------------

!    character(len=*), intent(in)      ::  action
    character(len=*), intent(in)      ::  key
    type(TllGridInfo), intent(in)     ::  pgi
    real, intent(in)                  ::  pg(0:pgi%im,pgi%jm)
    type(TllGridInfo), intent(in)     ::  tgi
    real, intent(inout)               ::  tg(0:tgi%im,tgi%jm)
    integer, intent(out)              ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/Match_u'

    ! --- local -----------------------------

    integer      ::  refine_i, refine_j
    integer      ::  cg_i1, cg_i2, cg_j1, cg_j2
    integer      ::  ci, cj
    integer      ::  fg_i, fg_j1, fg_j2
    integer      ::  fi, fj
    real         ::  fsum

    ! --- begin -----------------------------

!    select case ( action )
!
!      !
!      ! match fine grid with coarse cells:
!      !
!      case ( 'match' )

        ! determine refinement
        call GetRefinement( pgi, tgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! loop over cells in coarse grid covering fine grid:
        do cj = cg_j1, cg_j2
          do ci = cg_i1-1, cg_i2

            fg_i  = (ci-(cg_i1-1))*refine_i
            fg_j1 = (cj-cg_j1)*refine_j + 1  ;  fg_j2 = fg_j1-1 + refine_j

            ! sum over cells in fine grid:
            select case ( key )
              case ( 'sum' )
                fsum = sum( tg(fg_i,fg_j1:fg_j2) )
                ! distribute difference equally over all cells in fine grid:
                !                             (/6,4/)       + (  14     - 10 )/  2      = (/8,6/)
                tg(fg_i,fg_j1:fg_j2) = tg(fg_i,fg_j1:fg_j2) + (pg(ci,cj)-fsum)/refine_j
              case ( 'aver' )
                fsum = sum( tg(fg_i,fg_j1:fg_j2) )/(refine_j)
                ! add difference in averages to all cells in fine grid:
                !                             (/6,4/)       + (   7     -  5 )          = (/8,6/)
                tg(fg_i,fg_j1:fg_j2) = tg(fg_i,fg_j1:fg_j2) + (pg(ci,cj)-fsum)
              case default
                write (*,'("ERROR - unsupported key for match action:")')
!                write (*,'("ERROR -   action : ",a)') trim(action)
                write (*,'("ERROR -   key    : ",a)') trim(key)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

!            ! match fine grid with coarse cell:
!            if ( fsum == 0.0 ) then
!              if ( pg(ci,cj) /= 0.0 ) then
!                write (*,'("ERROR - zero sum over coarse cell (",i3,",",i3,")")') ci, cj
!                write (*,'("ERROR -   coarse cell : ",es12.4)') pg(ci,cj)
!                write (*,'("ERROR -   refine      : ",2i6)') refine_i, refine_j
!                write (*,'("ERROR -   target grid : ",es12.4)') tg(fg_i,fg_j1:fg_j2)
!                write (*,'("ERROR in ",a)') rname; status=1; return
!              end if
!            else
!              ! scale towards value of coarse grid:
!              tg(fg_i,fg_j1:fg_j2) = tg(fg_i,fg_j1:fg_j2)/fsum * pg(ci,cj)
!            end if

          end do
        end do

!      !
!      ! fine grid is subset of coarse grid:
!      !
!      case ( 'subset' )
!
!        ! determine refinement
!        call GetRefinement( pgi, tgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2 )
!
!        if ( (refine_i /= 1) .or. (refine_j /= 1) ) then
!          write (*,'("ERROR - for this action the fine grid should a subset of a coarse grid:")')
!          write (*,'("ERROR -   action     : ",a)') trim(action)
!          write (*,'("ERROR -   refinement : ",2i6)') refine_i, refine_j
!          write (*,'("ERROR in ",a)') rname; status=1; return
!        end if
!
!        ! loop over cells in coarse grid covering fine grid:
!        do cj = cg_j1, cg_j2
!          do ci = cg_i1-1, cg_i2
!
!            fg_i  = (ci-(cg_i1-1))
!            fg_j1 = (cj-cg_j1) + 1
!
!            ! fine grid is subset of coarse grid
!            tg(fg_i,fg_j1) = pg(ci,cj)
!
!          end do
!        end do
!
!      !
!      ! collect cells in fine grid to coarse grid
!      !
!      case ( 'combine' )
!
!        ! determine refinement
!        ! NOTE: parent grid is fine, targer is coarse !
!        call GetRefinement( tgi, pgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2 )
!
!        if ( (cg_i1 /= 1) .or. (cg_i2 /= tgi%im) .or. &
!             (cg_j1 /= 1) .or. (cg_j2 /= tgi%jm) ) then
!          write (*,'("ERROR - for this action the fine grid should cover the complete coarse grid:")')
!          write (*,'("ERROR -   action  : ",a)') trim(action)
!          write (*,'("ERROR -   covered : [",i4,",",i4,"] x [",i4,",","]")') cg_i1,cg_i2,cg_j1,cg_j2
!          write (*,'("ERROR -   of      : ",2i5)') tgi%im, tgi%jm
!          write (*,'("ERROR in ",a)') rname; status=1; return
!        end if
!
!        ! loop over cells in coarse grid covering fine grid:
!        do cj = cg_j1, cg_j2
!          do ci = cg_i1-1, cg_i2
!
!            fg_i  = (ci-(cg_i1-1))*refine_i
!            fg_j1 = (cj-cg_j1)*refine_j + 1  ;  fg_j2 = fg_j1-1 + refine_j
!
!            ! sum over cells in fine grid:
!            select case ( key )
!              case ( 'sum' )
!                fsum = sum( pg(fg_i,fg_j1:fg_j2) )
!              case ( 'aver' )
!                fsum = sum( pg(fg_i,fg_j1:fg_j2) )/(refine_j)
!              case default
!                write (*,'("ERROR - unsupported key for match action:")')
!                write (*,'("ERROR -   action : ",a)') trim(action)
!                write (*,'("ERROR -   key    : ",a)') trim(key)
!                write (*,'("ERROR in ",a)') rname; status=1; return
!            end select
!
!            ! collect cells in fine parent grid to target coarse grid
!            tg(ci,cj) = fsum
!
!          end do
!        end do
!
!      case default
!        write (*,'("ERROR - unknown action `",a,"`")') trim(action)
!        write (*,'("ERROR in ",a)') rname; status=1; return
!    end select

    ! ok
    status = 0

  end subroutine Match_u


  ! ***

  ! flux through north/south boundaries

  subroutine Match_v( key, pgi, pg, tgi, tg, status )

    ! --- in/out ------------------------------

!    character(len=*), intent(in)      ::  action
    character(len=*), intent(in)      ::  key
    type(TllGridInfo), intent(in)     ::  pgi
    real, intent(in)                  ::  pg(pgi%im,0:pgi%jm)
    type(TllGridInfo), intent(in)     ::  tgi
    real, intent(inout)               ::  tg(tgi%im,0:tgi%jm)
    integer, intent(out)              ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/Match_v'

    ! --- local -----------------------------

    integer      ::  refine_i, refine_j
    integer      ::  cg_i1, cg_i2, cg_j1, cg_j2
    integer      ::  ci, cj
    integer      ::  fg_i1, fg_i2, fg_j
    integer      ::  fi, fj
    real         ::  fsum

    ! --- begin -----------------------------

!    select case ( action )
!
!      !
!      ! match fine grid with coarse cells:
!      !
!      case ( 'match' )

        ! determine refinement
        call GetRefinement( pgi, tgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! loop over cells in coarse grid covering fine grid:
        do cj = cg_j1-1, cg_j2
          do ci = cg_i1, cg_i2

            fg_i1 = (ci-cg_i1)*refine_i + 1  ;  fg_i2 = fg_i1-1 + refine_i
            fg_j  = (cj-(cg_j1-1))*refine_j

            ! sum over cells in fine grid:
            select case ( key )
              case ( 'sum' )
                fsum = sum( tg(fg_i1:fg_i2,fg_j) )
                ! distribute difference equally over all cells in fine grid:
                !                             (/6,4/)       + (  14     - 10 )/  2      = (/8,6/)
                tg(fg_i1:fg_i2,fg_j) = tg(fg_i1:fg_i2,fg_j) + (pg(ci,cj)-fsum)/refine_i
              case ( 'aver' )
                fsum = sum( tg(fg_i1:fg_i2,fg_j) )/(refine_i)
                ! add difference in averages to all cells in fine grid:
                !                             (/6,4/)       + (   7     -  5 )          = (/8,6/)
                tg(fg_i1:fg_i2,fg_j) = tg(fg_i1:fg_i2,fg_j) + (pg(ci,cj)-fsum)
              case default
                write (*,'("ERROR - unsupported key for match action:")')
!                write (*,'("ERROR -   action : ",a)') trim(action)
                write (*,'("ERROR -   key    : ",a)') trim(key)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

!            ! match fine grid with coarse cell:
!            if ( fsum == 0.0 ) then
!              if ( pg(ci,cj) /= 0.0 ) then
!                write (*,'("ERROR - zero sum over coarse cell (",i3,",",i3,")")') ci, cj
!                write (*,'("ERROR -   coarse cell : ",es12.4)') pg(ci,cj)
!                write (*,'("ERROR -   refine      : ",2i6)') refine_i, refine_j
!                write (*,'("ERROR -   fine grid   : ",es12.4)') tg(fg_i1:fg_i2,fg_j)
!                write (*,'("ERROR in ",a)') rname; status=1; return
!              end if
!            else
!              ! scale towards value of coarse grid:
!              tg(fg_i1:fg_i2,fg_j) = tg(fg_i1:fg_i2,fg_j)/fsum * pg(ci,cj)
!            end if

          end do
        end do

!      !
!      ! fine grid is subset of coarse grid:
!      !
!      case ( 'subset' )
!
!        ! determine refinement
!        call GetRefinement( pgi, tgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2 )
!
!        if ( (refine_i /= 1) .or. (refine_j /= 1) ) then
!          write (*,'("ERROR - for this action the fine grid should a subset of a coarse grid:")')
!          write (*,'("ERROR -   action     : ",a)') trim(action)
!          write (*,'("ERROR -   refinement : ",2i6)') refine_i, refine_j
!          write (*,'("ERROR in ",a)') rname; status=1; return
!        end if
!
!        ! loop over cells in coarse grid covering fine grid:
!        do cj = cg_j1-1, cg_j2
!          do ci = cg_i1, cg_i2
!
!            fg_i1 = (ci-cg_i1) + 1
!            fg_j  = (cj-(cg_j1-1))
!
!            ! fine grid is subset of coarse grid
!            tg(fg_i1,fg_j) = pg(ci,cj)
!
!          end do
!        end do
!
!      !
!      ! collect cells in fine grid to coarse grid
!      !
!      case ( 'combine' )
!
!        ! determine refinement
!        ! NOTE: parent grid is fine, target is coarse !
!        call GetRefinement( tgi, pgi, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2 )
!
!        if ( (cg_i1 /= 1) .or. (cg_i2 /= tgi%im) .or. &
!             (cg_j1 /= 1) .or. (cg_j2 /= tgi%jm) ) then
!          write (*,'("ERROR - for this action the fine grid should cover the complete coarse grid:")')
!          write (*,'("ERROR -   action  : ",a)') trim(action)
!          write (*,'("ERROR -   covered : [",i4,",",i4,"] x [",i4,",","]")') cg_i1,cg_i2,cg_j1,cg_j2
!          write (*,'("ERROR -   of      : ",2i5)') tgi%im, tgi%jm
!          write (*,'("ERROR in ",a)') rname; status=1; return
!        end if
!
!        ! loop over cells in coarse grid covering fine grid:
!        do cj = cg_j1-1, cg_j2
!          do ci = cg_i1, cg_i2
!
!            fg_i1 = (ci-cg_i1)*refine_i + 1  ;  fg_i2 = fg_i1-1 + refine_i
!            fg_j  = (cj-(cg_j1-1))*refine_j
!
!            ! sum over cells in fine grid:
!            select case ( key )
!              case ( 'sum' )
!                fsum = sum( pg(fg_i1:fg_i2,fg_j) )
!              case ( 'aver' )
!                fsum = sum( pg(fg_i1:fg_i2,fg_j) )/(refine_i)
!              case default
!                write (*,'("ERROR - unsupported key for match action:")')
!                write (*,'("ERROR -   action : ",a)') trim(action)
!                write (*,'("ERROR -   key    : ",a)') trim(key)
!                write (*,'("ERROR in ",a)') rname; status=1; return
!            end select
!
!            ! collect cells in fine parent grid to target coarse grid
!            tg(ci,cj) = fsum
!
!          end do
!        end do
!
!      case default
!        write (*,'("ERROR - unknown action `",a,"`")') trim(action)
!        write (*,'("ERROR in ",a)') rname; status=1; return
!    end select

    ! ok
    status = 0

  end subroutine Match_v



  ! ========================================================
  ! ===
  ! === fill grid from other grid
  ! ===
  ! ========================================================


  !
  ! NUV
  !
  !   Key to identify data positions:
  !    'n'  :  value valid for cell (center)            ll(1:nlon  ,1:nlat  )
  !    'u'  :  value valid for east/west boundaries     ll(1:nlon+1,1:nlat  )
  !    'v'  :  value valid for north/south boundaries   ll(1:nlon  ,1:nlat+1)
  !
  ! ROUTINES
  !
  !   call FillGrid( lli, nuv, ll, lliX, nuvX, llX, combkey, status [,llX_w] )
  !
  !     Fill ll (defined by lli,nuv) with values from llX (defined by lliX,nuvX)
  !
  !     Coverage of lli by lliX :
  !      o lliX is larger than or equal to lli   ->  all cells in ll changed
  !      o lliX is smaller than lli              ->  only part of ll is changed
  !
  !     Create new ll from llX:
  !      o llX is superset  ->  copy values from llX into ll
  !      o llX is fine      ->  fill ll by combining cells in llX
  !                             (average/sum/etc given the combine key)
  !
  !        Combine keys:
  !
  !          'sum'         :   sum_i llX_i
  !
  !          'aver'        :   sum_i llX_i          /  sum_i i
  !
  !          'area-aver'   :   sum_i llX_i A_i      /  sum_i A_i
  !
  !          'weight'      :   sum_i llX_i llX_w_i  /  sum_i llX_w_i
  !            (only for nuv='n')
  !
  !      Return status:
  !         0 : ok
  !        -1 : only part of ll is filled
  !

  !
  !   AdjustGrid  NOT IMPLEMENTED YET
  !
  !     Adjust ll given llX:
  !      o llX is coarse   ->  adjust ll such that average/sum/.. of lli matches
  !

  subroutine FillGrid( lli, nuv, ll, lliX, nuvX, llX, combkey, status, llX_w )

    ! --- in/out --------------------------------

    type(TllGridInfo), intent(in)        ::  lli
    character(len=*), intent(in)         ::  nuv
    real, intent(out)                    ::  ll(:,:)
    type(TllGridInfo), intent(in)        ::  lliX
    character(len=*), intent(in)         ::  nuvX
    real, intent(in)                     ::  llX(:,:)
    character(len=*), intent(in)         ::  combkey
    integer, intent(out)                 ::  status
    real, intent(in), optional           ::  llX_w(:,:)

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/FillGrid'

    ! --- local ---------------------------------

    character(len=10)   ::  action
    integer             ::  di, dj
    integer             ::  i1, i2, j1, j2
    integer             ::  i, j
    integer             ::  i1X, i2X, j1X, j2X
    integer             ::  iX, jX
    integer             ::  diX, djX, nX
    real                ::  res, resw
    integer             ::  ia, ib, ja, jb
    integer             ::  iaX, ibX, jaX, jbX
    real                ::  llX_ab

    real, allocatable   ::  wwX(:,:)
    logical             ::  wwdiv

    logical             ::  only_part_of_ll

    ! --- begin ---------------------------------

    ! check input ...
    if ( nuv /= nuvX ) then
      write (*,'("ERROR - nuv keys should be equal:")') combkey
      write (*,'("ERROR -   nuv    : `",a,"`")') nuv
      write (*,'("ERROR -   nuvX   : `",a,"`")') nuvX
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! determine how lliX is related to lli:
    !   i1, i2, j1, j2       :  cell ranges in lli covered by cells of lliX
    !   i1X, i2X, j1X, j2X   :  cell ranges in lliX covering cells of lliX
    !   action               :  'copy', 'combine'
    call Relate( lli , i1 , i2 , j1 , j2 , &
                 lliX, i1X, i2X, j1X, j2X, &
                 action, status )
    select case ( status )
      case ( -1 )
        only_part_of_ll = .true.
      case ( 0 )
        only_part_of_ll = .false.
      case default
        write (*,'("ERROR in ",a)') rname; return
    end select

    ! what to do with cells in llX?
    select case ( action )

      !
      ! * copy
      !

      ! 'copy' : lli and lliX define same resolution;
      ! the cells from llX area [i1X,i2X] x [j1X,j2X]
      ! should be copied into ll area [i1,i2] x [j1,j2]
      case ( 'copy'  )

        select case ( nuv )

          case ( 'n' )

            ! loop over (selection of) cells of target grid lli:
            !OpenMP does not make sense here, it is just a copy ...
            !xOMP PARALLEL &
            !xOMP   default (none) &
            !xOMP   shared  (i1, i2, j1, j2, i1x, j1x, llx, ll) &
            !xOMP   private (i, j, ix, jx)
            !xOMP   DO
            do j = j1, j2
              do i = i1, i2
                ! source cell in llX:
                iX = i1X + i-i1
                jX = j1X + j-j1
                ! copy cell:
                ll(i,j) = llX(iX,jX)
              end do
            end do
            !xOMP   END DO
            !xOMP END PARALLEL

          case ( 'u' )

            ! loop over (selection of) cells of target grid lli:
            do j = j1, j2
              do i = i1, i2+1
                ! source cell in llX:
                iX = i1X + i-i1
                jX = j1X + j-j1
                ! copy cell:
                ll(i,j) = llX(iX,jX)
              end do
            end do

          case ( 'v' )

            ! loop over (selection of) cells of target grid lli:
            do j = j1, j2+1
              do i = i1, i2
                ! source cell in llX:
                iX = i1X + i-i1
                jX = j1X + j-j1
                ! copy cell:
                ll(i,j) = llX(iX,jX)
              end do
            end do

          case default
            write (*,'("ERROR - unsupported nuv `",a,"`")') nuv
            write (*,'("ERROR -   action : `",a,"`")') action
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      !
      ! * combine
      !

      ! 'combine' : lliX defines a fine resolution;
      ! the cells from llX area [i1X,i2X] x [j1X,j2X]
      ! should be combined and copied into ll area [i1,i2] x [j1,j2]
      case ( 'combine' )

        ! resolution of fine cells (lliX) in cell of lli,
        ! for example 3 x 2 :
        diX = (i2X-i1X+1) / (i2-i1+1)
        djX = (j2X-j1X+1) / (j2-j1+1)
        nX = diX * djX

        select case ( nuv )

          case ( 'n' )

            !
            ! set weight:
            !
            !   'weight'    :  (sum_i  f_i w_i) / (sum_i w_i)
            !
            !   'sum'       :  (sum_i f_i)                    :  w = 1.0
            !   'aver'      :  (sum_i f_i) / (sum_i)          :  w = 1.0,  wwdiv
            !   'area-aver' :  (sum_i f_i A_i) / (sum_i A_i)  :  w = A_i,  wwdiv

            ! same size as input grid:
            allocate( wwX(size(llX,1),size(llX,2)) )

            ! weight provided as argument ?
            if ( combkey == 'weight' ) then
              if ( .not. present(llX_w) ) then
                write (*,'("ERROR - combkey `weight` but llX_w not present ...")')
                write (*,'("ERROR in ",a)') rname; status=1; return
              end if
              if ( any( shape(llX_w) /= shape(llX) ) ) then
                write (*,'("ERROR - weight should have shape of input grid:")')
                write (*,'("ERROR -   shape(llX_w) : `",2i6)') shape(llX_w)
                write (*,'("ERROR -   shape(llX)   : `",2i6)') shape(llX)
                write (*,'("ERROR in ",a)') rname; status=1; return
              end if
              wwX = llX_w
              wwdiv = .true.
            else
              if ( present(llX_w) ) then
                write (*,'("ERROR - llX_w pressent but no combkey `weight` ...")')
                write (*,'("ERROR in ",a)') rname; status=1; return
              end if
              ! fill weight given combkey:
              select case ( combkey )
                case ( 'sum' )
                  wwX = 1.0
                  wwdiv = .false.
                case ( 'aver' )
                  wwX = 1.0
                  wwdiv = .true.
                case ( 'area-aver' )
                  call AreaOper( lliX, wwX, '=', 'm2', status )
                  if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                  wwdiv = .true.
                case default
                  write (*,'("ERROR - unsupported combkey `",a,"`")') combkey
                  write (*,'("ERROR -   nuv    : `",a,"`")') nuv
                  write (*,'("ERROR -   action : `",a,"`")') action
                  write (*,'("ERROR in ",a)') rname; status=1; return
              end select
            end if

            ! loop over (selection of) cells of target grid lli:
            !$OMP PARALLEL &
            !$OMP   default (none) &
            !$OMP   shared  (i1, i2, j1, j2, i1x, j1x ) &
            !$OMP   shared  (djx, dix, llx, wwx, wwdiv, ll ) &
            !$OMP   private (i, j, res, resw, ix, jx)
            !$OMP   DO
            do j = j1, j2
              do i = i1, i2

                ! start with zero result:
                res = 0.0
                resw = 0.0
                ! loop over source cells in llX:
                do jX = j1X + (j-j1)*djX, j1X + (j-j1+1)*djX-1
                  do iX = i1X + (i-i1)*diX, i1X + (i-i1+1)*diX-1
                    res  = res  + llX(iX,jX) * wwX(iX,jX)
                    resw = resw + wwX(iX,jX)
                  end do
                end do
                ! fill result:
                if ( wwdiv ) then
                  ll(i,j) = res / resw
                else
                  ll(i,j) = res
                end if

              end do
            end do
            !$OMP   END DO
            !$OMP END PARALLEL

            ! clear
            deallocate( wwX )

          case ( 'u' )

            ! loop over (selection of) cells of target grid lli:
            do j = j1, j2
              do i = i1, i2+1

                ! start with zero result:
                res = 0.0
                ! loop over points on east/west boundary in llX:
                iX = i1X + (i-i1)*diX
                do jX = j1X + (j-j1)*djX, j1X + (j-j1+1)*djX-1
                  ! add contribution of this llX cell:
                  select case ( combkey )
                    case ( 'sum' )
                      res = res + llX(iX,jX)
                    case ( 'aver', 'area-aver' )
                      res = res + llX(iX,jX)/real(djX)
                    case default
                      write (*,'("ERROR - unsupported combkey `",a,"`")') combkey
                      write (*,'("ERROR -   nuv    : `",a,"`")') nuv
                      write (*,'("ERROR -   action : `",a,"`")') action
                      write (*,'("ERROR in ",a)') rname; status=1; return
                  end select
                end do

                ! fill result:
                ll(i,j) = res

              end do
            end do

          case ( 'v' )

            ! loop over (selection of) cells of target grid lli:
            do j = j1, j2+1
              do i = i1, i2

                ! start with zero result:
                res = 0.0
                ! loop over points on north/south boundary in llX:
                jX = j1X + (j-j1)*djX
                do iX = i1X + (i-i1)*diX, i1X + (i-i1+1)*diX-1
                  ! add contribution of this llX cell:
                  select case ( combkey )
                    case ( 'sum' )
                      res = res + llX(iX,jX)
                    case ( 'aver', 'area-aver' )
                      res = res + llX(iX,jX)/real(diX)
                    case default
                      write (*,'("ERROR - unsupported combkey `",a,"`")') combkey
                      write (*,'("ERROR -   nuv    : `",a,"`")') nuv
                      write (*,'("ERROR -   action : `",a,"`")') action
                      write (*,'("ERROR in ",a)') rname; status=1; return
                  end select
                end do

                ! fill result:
                ll(i,j) = res

              end do
            end do

          case default
            write (*,'("ERROR - unsupported nuv `",a,"`")') nuv
            write (*,'("ERROR -   action : `",a,"`")') action
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      !
      ! * distribute
      !

      ! lliX defines a coarse resolution;
      ! the cells from llX area [i1X,i2X] x [j1X,j2X]
      ! should be sampled onto ll area [i1,i2] x [j1,j2]
      !
      ! Note: it is not possible to produce weighted distributions
      ! to have different values within a coarse grid,
      ! for example based on area, since what to do with a cell with
      ! zero area ?

      case ( 'distribute' )

        ! resolution of fine cells (lli) in cell of lliX,
        ! for example 3 x 2 :
        di = (i2-i1+1) / (i2X-i1X+1)
        dj = (j2-j1+1) / (j2X-j1X+1)

        select case ( nuv )

          case ( 'n' )

            ! loop over (selection of) coarse cells in coarse source grid lliX:
            do jX = j1X, j2X
              do iX = i1X, i2X

                ! loop over cells in fine target grid covered by coarse cell:
                do j = j1 + (jX-j1X)*dj, j1 + (jX-j1X)*dj + dj-1
                  do i = i1 + (iX-i1X)*di, i1 + (iX-i1X)*di + di-1

                    ! copy or take fraction of coarse value:
                    select case ( combkey )
                      case ( 'sum' )
                        ! coarse cell is sum of finer; take fraction:
                        !ll(i,j) = llX(iX,jX) / real(di*dj)
                        ! coarse cell is sum of finer; take area fractions:
                        ll(i,j) = llX(iX,jX) * lli%area(j) / lliX%area(jX)
                      case ( 'aver', 'area-aver', 'weight' )
                        ! coarse cell is aver of finer; take all the same:
                        ll(i,j) = llX(iX,jX)
                      case default
                        write (*,'("ERROR - unsupported combkey `",a,"`")') combkey
                        write (*,'("ERROR -   nuv    : `",a,"`")') nuv
                        write (*,'("ERROR -   action : `",a,"`")') action
                        write (*,'("ERROR in ",a)') rname; status=1; return
                    end select

                  end do
                end do

              end do
            end do

          case ( 'u' )

            !
            ! coarse      iaX         ibX
            !              |           |           |
            !          |   |   |   |   |   |   |   |   |
            ! fine        ia          ib
            !                  i
            !

            ! loop over (selection of) u boundaries of target grid lli:
            do j = j1, j2
              do i = i1, i2+1

                ! coarse cell in j direction:
                jX = j1X + floor(real(j-j1)/real(dj))

                ! left and right bound in fine grid that are covered by coarse bound:
                ia = i1 +   floor(real(i-i1)/real(di))*di
                ib = i1 + ceiling(real(i-i1)/real(di))*di

                ! corresponding coarse boundaries:
                iaX = i1X + (ia-i1)/di
                ibX = i1X + (ib-i1)/di

                ! interpolation in i direction of surounding lliX boundaries:
                if ( iaX == ibX ) then
                  llX_ab = llX(iaX,jX)
                else
                  llX_ab = llX(iaX,jX) * real(ib-i)/real(di) + &
                           llX(ibX,jX) * real(i-ia)/real(di)
                end if

                ! copy or take fraction of coarse value:
                select case ( combkey )
                  case ( 'sum' )
                    ! coarse cell is sum of finer; take fraction:
                    ll(i,j) = llX_ab / real(dj)
                  case ( 'aver', 'area-aver' )
                    ! coarse cell is aver of finer; take all the same:
                    ll(i,j) = llX_ab
                  case default
                    write (*,'("ERROR - unsupported combkey `",a,"`")') combkey
                    write (*,'("ERROR -   nuv    : `",a,"`")') nuv
                    write (*,'("ERROR -   action : `",a,"`")') action
                    write (*,'("ERROR in ",a)') rname; status=1; return
                end select

              end do
            end do

          case ( 'v' )

            !
            ! coarse       fine
            !
            !           --
            !    jaX -- -- ja
            !           -- j
            !           --
            !    jbX -- -- jb
            !           --
            !

            ! loop over (selection of) u boundaries of target grid lli:
            do i = i1, i2
              do j = j1, j2+1

                ! coarse cell in i direction:
                iX = i1X + floor(real(i-i1)/real(di))

                ! left and right bound in fine grid that are covered by coarse bound:
                ja = j1 +   floor(real(j-j1)/real(dj))*dj
                jb = j1 + ceiling(real(j-j1)/real(dj))*dj

                ! corresponding coarse boundaries:
                jaX = j1X + (ja-j1)/dj
                jbX = j1X + (jb-j1)/dj

                ! interpolation in j direction of surounding lliX boundaries:
                if ( jaX == jbX ) then
                  llX_ab = llX(iX,jaX)
                else
                  llX_ab = llX(iX,jaX) * real(jb-j)/real(dj) + &
                           llX(iX,jbX) * real(j-ja)/real(dj)
                end if

                ! copy or take fraction of coarse value:
                select case ( combkey )
                  case ( 'sum' )
                    ! coarse cell is sum of finer; take fraction:
                    ll(i,j) = llX_ab / real(di)
                  case ( 'aver', 'area-aver' )
                    ! coarse cell is aver of finer; take all the same:
                    ll(i,j) = llX_ab
                  case default
                    write (*,'("ERROR - unsupported combkey `",a,"`")') combkey
                    write (*,'("ERROR -   nuv    : `",a,"`")') nuv
                    write (*,'("ERROR -   action : `",a,"`")') action
                    write (*,'("ERROR in ",a)') rname; status=1; return
                end select

              end do
            end do

          case default
            write (*,'("ERROR - unsupported nuv `",a,"`")') nuv
            write (*,'("ERROR -   action : `",a,"`")') action
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      !
      ! * error ....
      !

      case default
        write (*,'("ERROR - unsupported action `",a,"`")') action
        write (*,'("ERROR -   lliX lon : ",f7.2,f6.2,i4)') lliX%lon_deg(1), lliX%dlon_deg, lliX%nlon
        write (*,'("ERROR -        lat : ",f7.2,f6.2,i4)') lliX%lat_deg(1), lliX%dlat_deg, lliX%nlat
        write (*,'("ERROR -    lli lon : ",f7.2,f6.2,i4)')  lli%lon_deg(1),  lli%dlon_deg,  lli%nlon
        write (*,'("ERROR -        lat : ",f7.2,f6.2,i4)')  lli%lat_deg(1),  lli%dlat_deg,  lli%nlat
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select

    ! ok
    if ( only_part_of_ll ) then
      status = -1
    else
      status = 0
    end if

  end subroutine FillGrid



  !
  ! determine how lliX is related to lli:
  !   i1, i2, j1, j2       :  cell ranges in lli covered by cells of lliX
  !   i1X, i2X, j1X, j2X   :  cell ranges in lliX covering cells of lliX
  !   action               :  'copy', 'combine'
  !

  subroutine Relate( lli , i1 , i2 , j1 , j2 , &
                     lliX, i1X, i2X, j1X, j2X, &
                     action, status )

    ! --- in/out --------------------------------

    type(TllGridInfo), intent(in)        ::  lli
    integer, intent(out)                 ::  i1, i2, j1, j2
    type(TllGridInfo), intent(in)        ::  lliX
    integer, intent(out)                 ::  i1X, i2X, j1X, j2X
    character(len=10), intent(out)       ::  action
    integer, intent(out)                 ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  name = mname//'/Relate'

    ! --- local ---------------------------------

    real           ::  dlon, dlonX
    integer        ::  nlon, nlonX
    real           ::  dlat, dlatX
    integer        ::  nlat, nlatX
    real           ::  west, westX
    real           ::  east, eastX
    real           ::  south, southX
    real           ::  north, northX
    real           ::  r, rX

    logical        ::  only_part_of_ll

    ! --- begin ---------------------------------

    ! shorthands ...
    dlon  = lli%dlon_deg         ;  dlonX  = lliX%dlon_deg
    nlon  = lli%nlon             ;  nlonX  = lliX%nlon
    dlat  = lli%dlat_deg         ;  dlatX  = lliX%dlat_deg
    nlat  = lli%nlat             ;  nlatX  = lliX%nlat
    west  = lli%blon_deg(0)      ;  westX  = lliX%blon_deg(0)
    east  = lli%blon_deg(nlon)   ;  eastX  = lliX%blon_deg(nlonX)
    south = lli%blat_deg(0)      ;  southX = lliX%blat_deg(0)
    north = lli%blat_deg(nlat)   ;  northX = lliX%blat_deg(nlatX)

    ! same spacing ?
    if ( (dlonX == dlon) .and. (dlatX == dlat) ) then
      ! cells from lliX should be copied to lli
      action = 'copy'
    else if ( (dlonX <= dlon) .and. (dlatX <= dlat) ) then
      ! cells from lliX should be combined:
      action = 'combine'
    else if ( (dlonX >= dlon) .and. (dlatX >= dlat) ) then
      ! distribute coarse cells of lliX over fine cells of lli:
      action = 'distribute'
    else
      write (*,'("ERROR - do not know how to relate grids:")')
      write (*,'("ERROR -   lli  dlon,dlat :",2f10.4)') dlon , dlat
      write (*,'("ERROR -   lliX dlon,dlat :",2f10.4)') dlonX, dlatX
      write (*,'("ERROR in ",a)') name; status=1; return
    end if

    ! assume by default that all is filled
    only_part_of_ll = .false.

    ! west boundary
    r  = abs(west - westX) / dlon
    rX = abs(west - westX) / dlonX
    ! relate lliX to lli:
    if ( (westX <= west) .and. (rX == nint(rX)) .and. (nint(rX)+1 <= nlonX) ) then
      ! at this side, all of lli is covered by lliX:
      i1  = 1
      i1X = nint(rX) + 1
    else if ( (westX > west) .and. (r == nint(r)) .and. (nint(r)+1 <= nlon) ) then
      ! at this side, only a part of lli is covered by lliX:
      i1  = nint(r) + 1
      i1X = 1
      only_part_of_ll = .true.
    else
      write (*,'("ERROR - do not know how to relate west bounds:")')
      write (*,'("ERROR -   lli  bound, spacing, number :",2f10.4,i6)') west , dlon , nlon
      write (*,'("ERROR -   lliX bound, spacing, number :",2f10.4,i6)') westX, dlonX, nlonX
      write (*,'("ERROR in ",a)') name; status=1; return
    end if

    ! east boundary
    r  = abs(east - eastX) / dlon
    rX = abs(east - eastX) / dlonX
    ! relate lliX to lli:
    if ( (eastX >= east) .and. (rX == nint(rX)) .and. (nlonX-nint(rX) >= 1) ) then
      ! at this side, all of lli is covered by lliX:
      i2  = nlon
      i2X = nlonX - nint(rX)
    else if ( (eastX < east) .and. (r == nint(r)) .and. (nlon-nint(r) >= 1) ) then
      ! at this side, only a part of lli is covered by lliX:
      i2  = nlon  - nint(r)
      i2X = nlonX
      only_part_of_ll = .true.
    else
      write (*,'("ERROR - do not know how to relate east bounds:")')
      write (*,'("ERROR -   lli  bound, spacing, number :",2f10.4,i6)') east , dlon , nlon
      write (*,'("ERROR -   lliX bound, spacing, number :",2f10.4,i6)') eastX, dlonX, nlonX
      write (*,'("ERROR in ",a)') name; status=1; return
    end if

    ! south boundary
    r  = abs(south - southX) / dlat
    rX = abs(south - southX) / dlatX
    ! relate lliX to lli:
    if ( (southX <= south) .and. (rX == nint(rX)) .and. (nint(rX)+1 <= nlatX) ) then
      ! at this side, all of lli is covered by lliX:
      j1  = 1
      j1X = nint(rX) + 1
    else if ( (southX > south) .and. (r == nint(r)) .and. (nint(r)+1 <= nlat) ) then
      ! at this side, only a part of lli is covered by lliX:
      j1  = nint(r) + 1
      j1X = 1
      only_part_of_ll = .true.
    else
      write (*,'("ERROR - do not know how to relate south bounds:")')
      write (*,'("ERROR -   lli  bound, spacing, number :",2f10.4,i6)') south , dlat , nlat
      write (*,'("ERROR -   lliX bound, spacing, number :",2f10.4,i6)') southX, dlatX, nlatX
      write (*,'("ERROR in ",a)') name; status=1; return
    end if

    ! north boundary
    r  = abs(north - northX) / dlat
    rX = abs(north - northX) / dlatX
    ! relate lliX to lli:
    if ( (northX >= north) .and. (rX == nint(rX)) .and. (nlatX-nint(rX) >= 1) ) then
      ! at this side, all of lli is covered by lliX:
      j2  = nlat
      j2X = nlatX - nint(rX)
    else if ( (northX < north) .and. (r == nint(r)) .and. (nlat-nint(r) >= 1) ) then
      ! at this side, only a part of lli is covered by lliX:
      j2  = nlat  - nint(r)
      j2X = nlatX
      only_part_of_ll = .true.
    else
      write (*,'("ERROR - do not know how to relate north bounds:")')
      write (*,'("ERROR -   lli  bound, spacing, number :",2f10.4,i6)') north , dlat , nlat
      write (*,'("ERROR -   lliX bound, spacing, number :",2f10.4,i6)') northX, dlatX, nlatX
      write (*,'("ERROR in ",a)') name; status=1; return
    end if

    ! ok
    if ( only_part_of_ll ) then
      status = -1
    else
      status = 0
    end if

  end subroutine Relate


  ! ========================================================
  ! ===
  ! === grid data operations
  ! ===
  ! ========================================================


!  subroutine Divergence_ll( lli, u, v, div )
!
!    ! --- in/out -------------------------------
!
!    type(TllGridInfo), intent(in)  ::  lli
!    type(TllGrid), intent(in)      ::  u, v
!    type(TllGrid), intent(out)     ::  div
!
!    ! --- local --------------------------------
!
!    integer  ::  i, j
!    integer  ::  im1, ip1
!
!    real     ::  wuim1(0:lme)
!    real     ::  wuip1(0:lme)
!    real     ::  wvjm1(0:lme)
!    real     ::  wvjp1(0:lme)
!
!    ! --- begin ---------------------------------
!
!    if ( lli%i1 == lli%im ) then
!      stop 'at least grid of 2 columns'
!    end if
!
!    if ( lli%j1 == lli%jm ) then
!      stop 'at least grid of 2 rows'
!    end if
!
!   stop 'circular ?'
!
!    do j = lli%j1, lli%jm
!
!      do i = lli%i1, lli%im
!
!        if ( j==lli%j1 .or. j==lli%jm ) then
!
!          div%d(i,j) = 0.0
!
!        else
!
!          if ( i==lli%i1 ) then
!            im1 = lli%im
!            ip1 = i + 1
!          else if ( i==lli%im ) then
!            im1 = lli%im
!            ip1 = i + 1
!          else
!            im1 = i - 1
!            ip1 = i + 1
!          end if
!
!          wuim1 = wu(im1,j,:)
!          wuip1 = wu(ip1,j,:)
!
!          wvjm1 = wv(i,j-1,:)
!          wvjp1 = wv(i,j+1,:)
!
!          dphi2 = eclon(2) - eclon(0)
!          dlat2 = eclat(j+1) - eclat(j-1)
!
!          qam(i,j,:) = ( wuip1 - wuim1 ) / (ceclat(j)*dphi2*ae) + &
!                       ( wvjp1 - wvjm1 ) / (dlat2*ae)
!
!        end if
!
!      end do
!    end do
!
!  end subroutine divergence_uv


  ! ==================================================================
  !
  ! mass flux balance
  !
  ! ==================================================================



  !
  ! (copied from TMPP;
  !  code should be cleaned ..)
  !


  !
  ! lli  : fine grid definition
  ! cgi  : coarse grid covering fine grid.
  !
  ! On input, fine grid fields pu/pv/pw should match with
  ! their coarse grid equivalents.
  !

  subroutine BalanceMassFluxes_sm( lli, lme, pu, pv, ww, dmdt, cgi, b_ec, status )

    ! --- in/out -------------------------------------

    type(TllGridInfo), intent(in)  ::  lli
    integer, intent(in)            ::  lme
    real,intent(inout)             ::  pu(0:lli%im,lli%jm,lme)
    real,intent(inout)             ::  pv(lli%im,0:lli%jm,lme)
    real, intent(in)               ::  ww(lli%im,lli%jm,lme)
    real, intent(in)               ::  dmdt(lli%im,lli%jm)
    type(TllGridInfo), intent(in)  ::  cgi
    real, intent(in)               ::  b_ec(0:lme)
    integer, intent(out)           ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/BalanceMassFluxes_sm'

    ! --- local ----------------------------------

    integer               ::  refine_i, refine_j
    integer               ::  cg_i1, cg_i2, cg_j1, cg_j2
    integer               ::  ci, cj
    integer               ::  fg_i1, fg_i2, fg_j1, fg_j2
    integer               ::  fi, fj

    integer               ::  i,j,l

    integer               ::  ifail
    real                  ::  massc(lli%im,lli%jm)
    real                  ::  dm(lli%im,lli%jm)
    real, allocatable     ::  cqu(:,:), cqv(:,:)

    ! --- begin -----------------------------------

    !print *, '<pascha_zoom>'

    !print *, '  correct layers 1 to',lme,' with (eta dp/deta) and dps/dt '

    ! determine refinement
    call GetRefinement( cgi, lli, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! correction fluxes:
    allocate( cqu(0:refine_i,refine_j) )
    allocate( cqv(refine_i,0:refine_j) )

    do l = 1, lme

      ! total horizontal massconvergence
      massc = 0.0
      massc = massc + sum(pu(1:lli%im  ,:,1:l),3)   ! leaving through east
      massc = massc - sum(pu(0:lli%im-1,:,1:l),3)   ! enter through west
      massc = massc + sum(pv(:,1:lli%jm  ,1:l),3)   ! leaving through north
      massc = massc - sum(pv(:,0:lli%jm-1,1:l),3)   ! enter through south

      ! set difference
      dm = (-massc) - ww(:,:,l) -  b_ec(l)*dmdt

      ! loop over cells in coarse grid covering fine grid:
      do cj = cg_j1, cg_j2
        do ci = cg_i1, cg_i2

          fg_i1 = (ci-cg_i1)*refine_i + 1  ;  fg_i2 = fg_i1-1 + refine_i
          fg_j1 = (cj-cg_j1)*refine_j + 1  ;  fg_j2 = fg_j1-1 + refine_j

          ! do FFT
          call SolvePoissonEq_zoom( refine_i, refine_j, &
                                    dm(fg_i1:fg_i2,fg_j1:fg_j2), &
                                    cqu, cqv, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

          ! add corrections to original flux:
          pu(fg_i1-1:fg_i2,fg_j1  :fg_j2,l) = pu(fg_i1-1:fg_i2,fg_j1  :fg_j2,l) + cqu
          pv(fg_i1  :fg_i2,fg_j1-1:fg_j2,l) = pv(fg_i1  :fg_i2,fg_j1-1:fg_j2,l) + cqv

        end do
      end do

    end do

    ! done
    deallocate( cqu )
    deallocate( cqv )

    ! ok
    status = 0

  end subroutine BalanceMassFluxes_sm


  ! ***

  !
  !  pw       :  vertical flux in direction of increasing level
  !

  subroutine BalanceMassFluxes_m( lli, pu, pv, pw, dm_dt, cgi, dt_sec, status )

    ! --- in/out -------------------------------------

    type(TllGridInfo), intent(in)  ::  lli
    real,intent(inout)             ::  pu(:,:,:)
    real,intent(inout)             ::  pv(:,:,:)
    real, intent(in)               ::  pw(:,:,:)
    real, intent(in)               ::  dm_dt(:,:,:)
    type(TllGridInfo), intent(in)  ::  cgi
    real, intent(in)               ::  dt_sec
    integer, intent(out)           ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/BalanceMassFluxes_m'

    ! --- local ----------------------------------

    integer             ::  refine_i, refine_j
    integer             ::  cg_i1, cg_i2, cg_j1, cg_j2
    integer             ::  ci, cj
    integer             ::  fg_i1, fg_i2, fg_j1, fg_j2
    integer             ::  fi, fj

    integer             ::  i,j,l
    integer             ::  nlev
    real                ::  massc(lli%im,lli%jm)
    real                ::  dm(lli%im,lli%jm)
    real, allocatable   ::  cqu(:,:), cqv(:,:)

    ! --- begin -----------------------------------

    ! number of levels
    nlev = size(pu,3)

    ! check ...
    if ( (size(pu   ,1) /= lli%nlon+1) .or. (size(pu   ,2) /= lli%nlat  ) .or. &
         (size(pv   ,1) /= lli%nlon  ) .or. (size(pv   ,2) /= lli%nlat+1) .or. &
         (size(pw   ,1) /= lli%nlon  ) .or. (size(pw   ,2) /= lli%nlat  ) .or. &
         (size(dm_dt,1) /= lli%nlon  ) .or. (size(dm_dt,2) /= lli%nlat  ) .or. &
         (size(pv   ,3) /= nlev) .or. (size(pw,3) /= nlev+1) .or. &
         (size(dm_dt,3) /= nlev) ) then
      write (*,'("ERROR - input arrays do not match with grid definition or levels:")')
      write (*,'("ERROR -   lli   : ",i3," x ",i3)') lli%nlon, lli%nlat
      write (*,'("ERROR -   pu    : ",i3," x ",i3," x ",i3)') shape(pu)
      write (*,'("ERROR -   pv    : ",i3," x ",i3," x ",i3)') shape(pv)
      write (*,'("ERROR -   pw    : ",i3," x ",i3," x ",i3)') shape(pw)
      write (*,'("ERROR -   dm_dt : ",i3," x ",i3," x ",i3)') shape(dm_dt)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! determine refinement:
    !   refine_i, refine_j
    !   cg_i1, cg_i2, cg_j1, cg_j2 : cells in coarse grid cgi covering fine grid lli:
    call GetRefinement( cgi, lli, refine_i, refine_j, cg_i1, cg_i2, cg_j1, cg_j2, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        !
        ! determine how lliX (finer?) is related to lli:
        !   i1, i2, j1, j2       :  cell ranges in lli covered by cells of lliX
        !   i1X, i2X, j1X, j2X   :  cell ranges in lliX covering cells of lliX
        !   action               :  'copy', 'combine'   to map lliX to lli
        !
        !
        !subroutine Relate( lli , i1 , i2 , j1 , j2 , &
        !                   lliX, i1X, i2X, j1X, j2X, &
        !                   action, status )

    ! correction fluxes:
    allocate( cqu(0:refine_i,refine_j) )
    allocate( cqv(refine_i,0:refine_j) )

    ! loop over levels
    do l = 1, nlev

      ! total horizontal mass convergence
      massc = 0.0
      massc = massc + pu(1:lli%im  ,:,l)       ! enter   through west
      massc = massc - pu(2:lli%im+1,:,l)       ! leaving through east
      massc = massc + pv(:,1:lli%jm  ,l)       ! enter   through south
      massc = massc - pv(:,2:lli%jm+1,l)       ! leaving through north
      massc = massc + pw(:,:,l  )              ! enter   through bottom
      massc = massc - pw(:,:,l+1)              ! leaving through top

      ! set difference
      dm = massc - dm_dt(:,:,l)

      ! loop over cells in coarse grid covering fine grid:
      do cj = cg_j1, cg_j2
        do ci = cg_i1, cg_i2

          fg_i1 = (ci-cg_i1)*refine_i + 1  ;  fg_i2 = fg_i1-1 + refine_i
          fg_j1 = (cj-cg_j1)*refine_j + 1  ;  fg_j2 = fg_j1-1 + refine_j

          ! do FFT
          call SolvePoissonEq_zoom( refine_i, refine_j, &
                                    dm(fg_i1:fg_i2,fg_j1:fg_j2), &
                                    cqu, cqv, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

          ! add corrections to original flux:
          pu(fg_i1:fg_i2+1,fg_j1:fg_j2  ,l) = pu(fg_i1:fg_i2+1,fg_j1:fg_j2  ,l) + cqu
          pv(fg_i1:fg_i2  ,fg_j1:fg_j2+1,l) = pv(fg_i1:fg_i2  ,fg_j1:fg_j2+1,l) + cqv

        end do  ! loop over coarse cells i
      end do  ! loop over coarse cells j

    end do  ! loop over levels

    ! done
    deallocate( cqu )
    deallocate( cqv )

    ! ok
    status = 0

  end subroutine BalanceMassFluxes_m


  ! ***


  subroutine CheckMassBalance( lli, pu, pv, sp1, sp2, dt, &
                                 rms_max, mav_max, status )

    use Binas, only : grav

    ! --- in/out -------------------------------------

    type(TllGridInfo), intent(in)  ::  lli
    real,intent(inout)             ::  pu(:,:,:)
    real,intent(inout)             ::  pv(:,:,:)
    real, intent(in)               ::  sp1(:,:)
    real, intent(in)               ::  sp2(:,:)
    real, intent(in)               ::  dt
    real, intent(in)               ::  rms_max
    real, intent(in)               ::  mav_max
    integer, intent(out)           ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/CheckMassBalance'

    ! --- local ----------------------------------

    real                ::  dmuv_dt(lli%im,lli%jm)
    real                ::  dmsp_dt(lli%im,lli%jm)
    real                ::  dm_dt_diff(lli%im,lli%jm)
    real                ::  rms_diff
    real                ::  mav_diff

    ! --- begin -----------------------------------

    ! check ...
    if ( (size(pu ,1) /= lli%nlon+1) .or. (size(pu ,2) /= lli%nlat  ) .or. &
         (size(pv ,1) /= lli%nlon  ) .or. (size(pv ,2) /= lli%nlat+1) .or. &
         (size(sp1,1) /= lli%nlon  ) .or. (size(sp1,2) /= lli%nlat  ) .or. &
         (size(sp2,1) /= lli%nlon  ) .or. (size(sp2,2) /= lli%nlat  ) .or. &
         (size(pu,3) /= size(pv,3)) ) then
      write (*,'("ERROR - input arrays do not match with grid definition or levels:")')
      write (*,'("ERROR -   lli   : ",i3," x ",i3)') lli%nlon, lli%nlat
      write (*,'("ERROR -   pu    : ",i3," x ",i3," x ",i3)') shape(pu)
      write (*,'("ERROR -   pv    : ",i3," x ",i3," x ",i3)') shape(pv)
      write (*,'("ERROR -   sp1   : ",i3," x ",i3)') shape(sp1)
      write (*,'("ERROR -   sp2   : ",i3," x ",i3)') shape(sp2)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! total horizontal mass convergence (kg/s)
    dmuv_dt = 0.0
    dmuv_dt = dmuv_dt + sum(pu(1:lli%nlon  ,:,:),3)       ! enter   through west
    dmuv_dt = dmuv_dt - sum(pu(2:lli%nlon+1,:,:),3)       ! leaving through east
    dmuv_dt = dmuv_dt + sum(pv(:,1:lli%nlat  ,:),3)       ! enter   through south
    dmuv_dt = dmuv_dt - sum(pv(:,2:lli%nlat+1,:),3)       ! leaving through north

    ! convert from kg/s to Pa/s
    dmuv_dt = dmuv_dt * grav                     ! Pa/s/m2
    call AreaOper( lli, dmuv_dt, '/', 'm2', status )     ! Pa/s
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! same according to surface pressures:
    dmsp_dt = ( sp2 - sp1 ) / dt

    ! imbalance relative to average surface pressure:
    dm_dt_diff = ( dmuv_dt - dmsp_dt ) / ( (sp1+sp2)/2.0 )

    ! statistics
    rms_diff = sqrt(  sum(dm_dt_diff**2) / size(dm_dt_diff) )
    mav_diff = maxval( abs( dm_dt_diff ) )

    ! check ...
    if ( (rms_diff > rms_max) .or. (mav_diff > mav_max) ) then
      write (*,'("ERROR - found relative surface pressure imbalance;")')
      write (*,'("ERROR - difference  [dmuv/dt]/[(sp1+sp2)/2]  and  [dsp/dt]/[(sp1+sp2)/2]  :")')
      write (*,'("ERROR -   rms      : ",es12.4," (",es12.4,")")') rms_diff, rms_max
      write (*,'("ERROR -   max abs. : ",es12.4," (",es12.4,")")') mav_diff, mav_max
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine CheckMassBalance


  ! ***


  ! Return flux field (u,v) such that netto flux in a cell matches m:
  !
  !    du   dv
  !    -- + --  =  M
  !    dx   dy
  !
  ! Fluxes at grid boundaries set to zero.
  !
  ! IN/OUT
  !   input:
  !     m(1:im,1:jm)   : target netto flux
  !   output:
  !     u(0:im,1:jm)   : u flux through cell boundaries
  !     v(1:im,0:jm)   : v vlux
  !
  !
  !                   +.......+.......+.......+.......+
  !                   :       :       :       :       :
  !                   :   o   :   o   :   o   :   o   :
  !                   :       :       :       :       :
  !       jm  +.......+---v---+---v---+---v---+---v---+.......+
  !           :       |       |       |       |       |       :
  !           :   o   u   *   u   *   u   *   u   *   u   o   :
  !           :       |       |       |       |       |       :
  !        :  +.......+---v---+---v---+---v---+---v---+.......+
  !           :       |       |       |       |       |       :
  !           :   o   u   *   u   *   u   *   u   *   u   o   :
  !           :       |       |       |       |       |       :
  !        1  +.......+---v---+---v---+---v---+---v---+.......+
  !           :       |       |       |       |       |       :
  !           :   o   u   *   u   *   u   *   u   *   u   o   :
  !           :       |       |       |       |       |       :
  !    j = 0  +.......+---v---+---v---+---v---+---v---+.......+
  !                   :       :       :       :       :
  !                   :   o   :   o   :   o   :   o   :
  !                   :       :       :       :       :
  !                   +.......+.......+.......+.......+
  !
  !              i=   0       1       2       ..     im
  !
  !      * =  m,  Psi
  !      o =      Psi (periodic)
  !
  !
  ! ALGORITHM
  !
  !  1. Solve the Poisson equation:
  !
  !         d^2    d^2
  !       ( ---- + ---- ) Psi = M
  !         dx^2   dy^2
  !
  !     and return  F = (u,v) = ( dPsi/dx, dPsi/dy )
  !     The differential equations is actually a difference equation:
  !     sum of fluxes is prescribed, we assume that the fluxes are differences
  !     of a discrete potential Psi.
  !
  !  2. The sollution is periodic in i and j.
  !     To obtain the zero boundary conditions, substract the boundaries
  !     from each row and column.
  !


  subroutine SolvePoissonEq_zoom( im, jm, m, u, v, status )

    use Binas, only : pi
    use grid_singleton, only : fftkind, fft

    ! --- in/out -------------------------------------

    integer, intent(in)            ::  im, jm
    real, intent(in)               ::  m(im,jm)
    real, intent(out)              ::  u(0:im,jm)
    real, intent(out)              ::  v(im,0:jm)
    integer, intent(out)           ::  status

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/SolvePoissonEq_zoom'

    ! --- local ----------------------------------

    integer            ::  i,j

    real               ::  fac(im,jm)
    complex(fftkind)   ::  A(im,jm)
    complex(fftkind)   ::  X(im,jm)

    real               ::  psi(im,jm)
    real               ::  row(im), col(jm)

    ! --- begin -----------------------------------

    !write(*, '("im = ", i3, ", jm = ", i3)') im, jm

    ! precalculate cos/sin array fac
    !!$OMP PARALLEL &
    !!$OMP   default (none) &
    !!$OMP   shared  (im, jm, fac) &
    !!$OMP   private (i, j)
    !!$OMP   DO
    do j = 1, jm
      do i=1,im
        fac(i,j) = 2.0*(cos(2*pi*(i-1)/im)+cos(2*pi*(j-1)/jm)-2.)
      end do
    end do
    !!$OMP   END DO
    !!$OMP END PARALLEL
    fac(1,1) = 1.0   !to avoid division by zero

    ! do FFT

    A = cmplx( m, 0.0 )
    X = fft( A, stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - from first fft; stat = ",i6)') status
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! compute fourier coefficients of correction potential
    A = cmplx( real(X)/fac, -aimag(X)/fac )
    A(1,1) = (0.0,0.0)

    ! get correction potential
    X = fft( A, stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - from second fft; stat = ",i6)') status
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    psi = real( X )

    ! correction flux in lon direction:
    do j = 1, jm
      !u(:,j) = (/ psi(1,j)-psi(im,j), psi(2:im,j)-psi(1:im-1,j), psi(1,j)-psi(im,j) /)
      u(0     ,j) = psi(1   ,j)-psi(im    ,j)
      u(1:im-1,j) = psi(2:im,j)-psi(1:im-1,j)
      u(im    ,j) = psi(1   ,j)-psi(im    ,j)
    end do
    ! substract flux over east/west boundary:
    col = u(0,:)
    do i = 0, im
      u(i,:) = u(i,:) - col
    end do

    ! correction flux in lat direction;
    do i = 1, im
      !v(i,:) = (/ psi(i,1)-psi(i,jm), psi(i,2:jm)-psi(i,1:jm-1), psi(i,1)-psi(i,jm) /)
      v(i,0     ) = psi(i,1   ) - psi(i,jm    )
      v(i,1:jm-1) = psi(i,2:jm) - psi(i,1:jm-1)
      v(i,jm    ) = psi(i,1   ) - psi(i,jm    )
    end do
    ! substract flux over north/west boundary:
    row = v(:,0)
    do j = 0, jm
      v(:,j) = v(:,j) - row
    end do

    ! Correction for flux over poles.
    ! The equation is solved with a 2D Fourrier transform, which
    ! assumes that the grid is part of an into infinity periodical
    ! extended grid; the solution is thus periodic from east into west
    ! but also from southpole into north pole (donut configuration).
    ! To obtain a zero flux over the poles, first the orginal flux over
    ! the poles is computed:
    !   vpole = psi(:,1) - psi(:,jm)
    ! This flux is substracted from each of the lat fluxes:
    !   v(:,j) := v(:,j) - vpole
    ! In this way, the netto mass convergence is maintained,
    ! but (u,v) is not a potential flow anymore!
    ! We will not bother about this small 'error', since we are dealing
    ! with a correction to a much larger error.

  end subroutine SolvePoissonEq_zoom


  ! =================================================================
  ! ===
  ! ===  equivalent latitudes
  ! ===
  ! =================================================================

  !
  ! eqvlatb(1:jm+1)   : bounds of eqv.lat. bins  (/-90.0,..,90.0/)
  !
  ! inds(im,jm)       : cell (i,j) is in   eqvlatb( inds(i,j), inds(i,j)+1 )
  !

  subroutine llgrid_EquivLat( lli, ff, eqvlatb, inds )

    use Binas, only : deg2rad, pi
    use Num, only : Interval, Interp_Lin
    use grid_tools, only : ll_area

    ! --- in/out -----------------------------

    type(TllGridInfo), intent(in)    ::  lli
    real, intent(in)                 ::  ff(:,:)
    real, intent(out)                ::  eqvlatb(:)
    integer, intent(out)             ::  inds(:,:)

    integer  ::  status

    ! --- const ------------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_EquivLat'

    ! --- local ------------------------------

    real                  ::  fmin, fmax, f

    integer               ::  nbins
    real, allocatable     ::  bin_inds(:)
    real, allocatable     ::  bin_fval(:)
    real, allocatable     ::  bin_area(:)
    real, allocatable     ::  bin_areacum(:)

    real                  ::  indf

    real                  ::  area_tot
    real                  ::  area_glob, area_south

    integer               ::  ilast
    integer               ::  i, j

    ! --- begin ------------------------------

    !if ( (maxval(lli%blon_deg)-minval(lli%blon_deg) < 360.0) .or. &
    !     (maxval(lli%blat_deg)-minval(lli%blat_deg) < 180.0)      ) then
    !  write (*,'("ERROR - equivalent lats for global grids only yet")')
    !  write (*,'("ERROR -   lons : ",f8.2,"-",f8.2)') minval(lli%blon_deg), maxval(lli%blon_deg)
    !  write (*,'("ERROR -   lats : ",f8.2,"-",f8.2)') minval(lli%blat_deg), maxval(lli%blat_deg)
    !  write (*,'("ERROR in ",a)') name; stop
    !end if

    call Check( lli, 'n', ff, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if
    call Check( lli, 'n', inds, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if

    ! for the moment ...
    nbins = lli%jm
    if ( size(eqvlatb) /= nbins+1 ) then
      write (*,'("ERROR - routine returns jm+1 bounds of jm eqv.lat bins:")')
      write (*,'("ERROR -   size(eqvlatb)   : ",i6)') size(eqvlatb)
      write (*,'("ERROR -   nbins+1 (=jm+1) : ",i6)') nbins+1
      write (*,'("ERROR in ",a)') rname; stop
    end if

    ! number of bins, field values, total area
    allocate( bin_inds(nbins) )
    allocate( bin_fval(nbins) )
    allocate( bin_area(nbins) )
    allocate( bin_areacum(nbins) )

    ! linear ax of values in field:
    fmin = minval(ff)
    fmax = maxval(ff)
    do j = 1, nbins
      bin_inds(j) = j*1.0
      bin_fval(j) = fmin*(nbins-j)/real(nbins-1) + fmax*(j-1)/real(nbins-1)
    end do

    ! add weights for each cell to corresponding bin:
    ilast = 0
    bin_area = 0.0
    do i = 1, lli%im
      do j = 1, lli%jm
        ! search index in bins:
        call Interp_Lin( bin_fval, bin_inds, ff(i,j), indf, ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if
        ! store index:
        inds(i,j) = nint(indf)
        ! add area to corresponding bin:
        bin_area(inds(i,j)) = bin_area(inds(i,j)) + lli%area(j)   ! rad^2
      end do
    end do

    ! total area:
    area_tot = lli%im * sum(lli%area)

    ! check ...
    if ( (sum(bin_area)-area_tot)/area_tot > 0.01 ) then
      write (*,'("ERROR - total area and areas in bins do not match:")')
      write (*,'("ERROR -   total area     : ",es12.2)') area_tot
      write (*,'("ERROR -   collected area : ",es12.2)') sum(bin_area)
      write (*,'("ERROR in ",a)') rname; stop
    end if

    ! cumulative bin area
    bin_areacum(1) = bin_area(1)
    do j = 2, nbins
      bin_areacum(j) = bin_areacum(j-1) + bin_area(j)
    end do

    !
    ! blat(nlat)  -------------------------------  90
    !                 |          |          \
    !                 |----------|           |
    !                 |//////////|           |
    !                 |----------|           | glob
    !                 |          | \         |
    !                 |          |  | south  |
    !                 |          | /        /
    ! blat(nlat)  -------------------------------  -90
    !              blon(0)   blon(nlon)

    ! global area and southern part (all in rad!)
    area_glob  = ll_area( 0.0, lli%blon(lli%nlon)-lli%blon(0), -0.5*pi, 0.5*pi      )  ! rad^2
    area_south = ll_area( 0.0, lli%blon(lli%nlon)-lli%blon(0), -0.5*pi, lli%blat(0) )  ! rad^2

    ! convert cumulative ax to (part of) [-1,1]
    bin_areacum = (area_south + bin_areacum ) / area_glob  ! [ 0,1]
    bin_areacum = -1.0 + 2.0*bin_areacum                   ! [-1,1]

    where ( bin_areacum < -1.0 ) bin_areacum = -1.0
    where ( bin_areacum >  1.0 ) bin_areacum =  1.0

    ! convert to latitude:
    eqvlatb(1) = lli%blat_deg(0)  ! deg
    do j = 1, nbins
      eqvlatb(j+1) = asin( bin_areacum(j) ) / deg2rad
    end do

       !print *, 0, area_south/area_glob, eqvlatb(1)
       !do j=1,nbins
       !  print *, j, bin_areacum(j), eqvlatb(j+1)
       !end do
       !stop

    ! done
    deallocate( bin_inds )
    deallocate( bin_fval )
    deallocate( bin_area )
    deallocate( bin_areacum )

  end subroutine llgrid_EquivLat



  !
  ! eqvlat1, eqvlat2 : lower and upper bound of eqv.lat.
  !

  subroutine llgrid_EquivLat_sort( lli, ff, eqvlatb1, eqvlatb2, status )

    use Binas, only : deg2rad, pi
    use Num, only : Sort
    use grid_tools, only : ll_area

    ! --- in/out -----------------------------

    type(TllGridInfo), intent(in)    ::  lli
    real, intent(in)                 ::  ff(:,:)
    real, intent(out)                ::  eqvlatb1(:,:)
    real, intent(out)                ::  eqvlatb2(:,:)
    integer, intent(out)             ::  status

    ! --- const ------------------------------

    character(len=*), parameter  ::  rname = mname//'/llgrid_EquivLat_sort'

    ! --- local ------------------------------

    real            ::  ffsort(lli%nlon,lli%nlat)
    integer         ::  ijsort(lli%nlon,lli%nlat,2)

    real            ::  areacum
    real            ::  bin_areacum(lli%nlon,lli%nlat,2)

    real            ::  area_tot
    real            ::  area_glob, area_south

    integer         ::  i, j
    integer         ::  is, js

    ! --- begin ------------------------------

    ! check input
    call Check( lli, 'n', ff, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! sort input array
    call Sort( ff, ffsort, ijsort, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! loop over sorted values to compute cumulative areas:
    areacum = 0.0
    do js = 1, lli%nlat
      do is = 1, lli%nlon

        ! original indices
        i = ijsort(is,js,1)
        j = ijsort(is,js,2)

        ! store lower cumulative area:
        bin_areacum(i,j,1) = areacum      ! rad^2

        ! add area of this cell:
        areacum = areacum + lli%area(j)    ! rad^2

        ! store upper cumulative area:
        bin_areacum(i,j,2) = areacum      ! rad^2

      end do
    end do

    ! total area:
    area_tot = lli%nlon * sum(lli%area)

    ! check ...
    if ( (bin_areacum(lli%nlon,lli%nlat,2)-area_tot)/area_tot > 0.01 ) then
      write (*,'("ERROR - total area and areas in bins do not match:")')
      write (*,'("ERROR -   total area     : ",es12.2)') area_tot
      write (*,'("ERROR -   collected area : ",es12.2)') bin_areacum(lli%nlon,lli%nlat,2)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    !
    ! blat(nlat)  -------------------------------  90
    !                 |          |          \
    !                 |----------|           |
    !                 |//////////|           |
    !                 |----------|           | glob
    !                 |          | \         |
    !                 |          |  | south  |
    !                 |          | /        /
    ! blat(nlat)  -------------------------------  -90
    !              blon(0)   blon(nlon)

    ! global area and southern part (all in rad!)
    area_glob  = ll_area( 0.0, lli%blon(lli%nlon)-lli%blon(0), -0.5*pi, 0.5*pi      )  ! rad^2
    area_south = ll_area( 0.0, lli%blon(lli%nlon)-lli%blon(0), -0.5*pi, lli%blat(0) )  ! rad^2

    ! convert cumulative ax to (part of) [-1,1]
    bin_areacum = (area_south + bin_areacum ) / area_glob  ! [ 0,1]
    bin_areacum = -1.0 + 2.0*bin_areacum                   ! [-1,1]

    ! adhoc truncations ...
    where ( bin_areacum < -1.0 ) bin_areacum = -1.0
    where ( bin_areacum >  1.0 ) bin_areacum =  1.0

    ! convert to latitude:
    eqvlatb1 = asin( bin_areacum(:,:,1) ) / deg2rad
    eqvlatb2 = asin( bin_areacum(:,:,2) ) / deg2rad

    ! ok
    status = 0

  end subroutine llgrid_EquivLat_sort


  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ~~~
  ! ~~~ region box
  ! ~~~
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  subroutine llreg_Init( llreg, west_deg, east_deg, south_deg, north_deg, status )

    use Binas, only : deg2rad

    ! --- in/out --------------------------------

    type(TllRegion), intent(inout)    ::  llreg
    real, intent(in)                  ::  west_deg, east_deg, south_deg, north_deg
    integer, intent(out)              ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = 'llreg_Init'

    ! --- begin --------------------------------

    ! store region boundaries:
    llreg%west_deg  =  west_deg
    llreg%east_deg  =  east_deg
    llreg%south_deg = south_deg
    llreg%north_deg = north_deg

    ! idem in radians:
    llreg%west  =  west_deg * deg2rad
    llreg%east  =  east_deg * deg2rad
    llreg%south = south_deg * deg2rad
    llreg%north = north_deg * deg2rad

    ! ok
    status = 0

  end subroutine llreg_Init


  ! ***


  subroutine llreg_Done( llreg, status )

    ! --- in/out --------------------------------

    type(TllRegion), intent(inout)    ::  llreg
    integer, intent(out)              ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = 'llreg_Done'

    ! --- begin --------------------------------

    ! nothing to be done

    ! ok
    status = 0

  end subroutine llreg_Done


  ! ***


  subroutine llreg_Region_Apply_Factor_2d( lli, x, llreg, fac, status, complement )

    use grid_tools, only : ll_area_frac

    ! --- in/out ---------------------------

    type(TllGridInfo), intent(in)   ::  lli
    real, intent(inout)             ::  x(:,:)
    type(TllRegion), intent(in)     ::  llreg
    real, intent(in)                ::  fac
    integer, intent(out)            ::  status
    logical, intent(in), optional   ::  complement

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = 'llreg_Region_Apply_Factor_2d'

    ! --- local --------------------------------

    logical      ::  do_complement
    real         ::  frac
    integer      ::  i, j

    ! --- begin --------------------------------

    ! apply to interior of region or to complement ?
    do_complement = .false.
    if ( present(complement) ) do_complement = complement

    ! check ...
    call Check( lli, 'n', x, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! loop over all cells in grid:
    do j = 1, lli%nlat
      do i = 1, lli%nlon

        ! fraction of grid cell covered by region;
        ! input to routine in radians:
        frac = ll_area_frac( lli%blon(i-1), lli%blon(i), lli%blat(j-1), lli%blat(j), &
                             llreg%west   , llreg%east , llreg%south  , llreg%north    )

        ! apply to interior or to complement ?
        if ( do_complement ) then
          ! apply factor if cell is (partly) outside region box:
          if ( frac < 1.0 ) then
            x(i,j) = fac * ( 1.0 - frac ) * x(i,j)
          end if
        else
          ! apply factor if cell is (partly) inside region box:
          if ( frac > 0.0 ) then
            x(i,j) = fac * frac * x(i,j)
          end if
        end if

      end do  ! i
    end do  ! j

    ! ok
    status = 0

  end subroutine llreg_Region_Apply_Factor_2d



  ! ***


  subroutine llreg_Region_Apply_Factor_3d( lli, x, llreg, fac, status, complement )

    use grid_tools, only : ll_area_frac

    ! --- in/out ---------------------------

    type(TllGridInfo), intent(in)   ::  lli
    real, intent(inout)             ::  x(:,:,:)
    type(TllRegion), intent(in)     ::  llreg
    real, intent(in)                ::  fac
    integer, intent(out)            ::  status
    logical, intent(in), optional   ::  complement

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = 'llreg_Region_Apply_Factor_3d'

    ! --- local --------------------------------

    logical      ::  do_complement
    real         ::  frac
    integer      ::  i, j

    ! --- begin --------------------------------

    ! apply to interior of region or to complement ?
    do_complement = .false.
    if ( present(complement) ) do_complement = complement

    ! check ...
    call Check( lli, 'n', x(:,:,1), status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! loop over all cells in grid:
    do j = 1, lli%nlat
      do i = 1, lli%nlon

        ! fraction of grid cell covered by region;
        ! input to routine in radians:
        frac = ll_area_frac( lli%blon(i-1), lli%blon(i), lli%blat(j-1), lli%blat(j), &
                             llreg%west   , llreg%east , llreg%south  , llreg%north    )

        ! apply to interior or to complement ?
        if ( do_complement ) then
          ! apply factor if cell is (partly) outside region box:
          if ( frac < 1.0 ) then
            x(i,j,:) = fac * ( 1.0 - frac ) * x(i,j,:)
          end if
        else
          ! apply factor if cell is (partly) inside region box:
          if ( frac > 0.0 ) then
            x(i,j,:) = fac * frac * x(i,j,:)
          end if
        end if

      end do  ! i
    end do  ! j

    ! ok
    status = 0

  end subroutine llreg_Region_Apply_Factor_3d



end module grid_type_ll



