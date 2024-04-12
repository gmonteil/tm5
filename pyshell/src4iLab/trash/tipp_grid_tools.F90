!----------------------------------------------------------
! Module tipp_grid_tools contains routines for converting data
! from one grid to the other.
!
! There are two main possibilities:
!
! 1. (re)distribution: this is typically applied for emission fields
!       These routines are strictly mass conserving.
!
!    o dist_vertical: distribute lat-lon field between two heights
!                     in a 3D field
!
!    o dist_latlon:   redistribute lat-lon field to different lat-lon grid
!                     (both coarsening and refining is supported)
!
! 2. interpolation: this is typically applied for concentration fields
!
!    o intp_3d:       do 3D linear interpolation from one lat-lon+hybrid
!                     grid to the other
!----------------------------------------------------------
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!----------------------------------------------------------

module tipp_grid_tools

  ! --- modules ----------------------------

  use GO, only : gol, goPr, goErr

  use grid, only                  : TllGridInfo, TLevelInfo, Init, Done

  implicit none

  ! --- in/out -----------------------------

  private

  public :: dist_vertical
  public :: dist_latlon
!  public :: intp_3d
  public :: oro_ecmwf

  ! --- const ------------------------------

  character(len=*), parameter  :: mname = 'module tipp_grid_tools'

contains

  subroutine dist_vertical( fld2d, fld3d, lli, at, bt, hlow, hhigh, hunit, hdef, status )
    !----------------------------------------------------------
    ! Distribute a field between two heights
    ! Some simple assumptions are made (thus no meteo needed):
    !  1. MSL pressure is assumed to be 1013 hPa
    !  2. Surface pressure is calculated from ECMWF orography
    !  2a. If orography is not available, flat surface at h=0 is assumed
    !  3. Fixed (virtual) temperature profile assumed
    !
    !  o fld2d(nx,ny,n)    :  2D field to be distributed [in]
    !                          'n' is extra dimension to process multiple fields
    !                          (e.g. months) at once
    !  o fld3d(nx,ny,nz,n) :  3D field to which fld2d has to be added [inout]
    !  o lli               :  latlon grid structure for fld2d and fld3d [in]
    !  o at, bt(nz+1)      :  hybrid pressure coefficients for levels in fld3d [in]
    !  o hlow, hhigh       :  heights between which fld2d has to be distributed [in]
    !                           hlow is closest to the surface
    !  o hunit              : unit of heights: 'm' (meters) or 'p':pascal [in]
    !  o hdef              :  height definition (only used if hunit='m'):
    !                          'sfc' = relative to surface
    !                          'msl' = relative to mean sea level
    !                           The 'msl' option gives possibility that part of
    !                           the field is below surface. This part is added to
    !                           the lowest level in fld3d.
    !  o status            :  0=success; 1=failure
    !----------------------------------------------------------

    ! --- modules --------------------------------

    use binas, only                 : T0, p0, xm_air, grav, Rgas

    ! --- in/out ----------------------------------------------

    real, intent(in)               :: fld2d(:,:,:)
    real, intent(inout)            :: fld3d(:,:,:,:)
    type(TllGridInfo), intent(in)  :: lli
    real, intent(in)               :: at(:), bt(:)
    real, intent(in)               :: hlow, hhigh
    character(len=1), intent(in)   :: hunit
    character(len=3), intent(in)   :: hdef
    integer, intent(out)           :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: rname = mname//'/dist_vertical'
    real, parameter                :: Tav = T0   ! average temperature
    real, parameter                :: c = xm_air*grav/(Rgas*Tav) ! a constant

    ! --- local -----------------------------------------------

    integer                        :: nx, ny, nz, n, i, j, l, hfac
    real                           :: f, p1b, p1t, p2b, p2t
    real                           :: tot2d, tot3db, tot3da, totbs
    real, dimension(:,:), allocatable :: oro, psurf
    real, dimension(:,:), allocatable :: logplow, logphigh
    real, dimension(:,:), allocatable :: logpdown, logpup

    ! --- begin -----------------------------------------------

    ! Dimensions
    nx = size( fld3d, 1 )  ! longitude dimension
    ny = size( fld3d, 2 )  ! latitude dimension
    nz = size( fld3d, 3 )  ! vertical dimension
    n = size( fld3d, 4 )   ! extra dimension

    ! Check dimensions
    if ( size( fld2d, 1 ) /= nx .or. size( fld2d, 2 ) /= ny .or. &
         size( fld2d, 3 ) /= n .or. &
         lli%im /= nx .or. lli%jm /= ny .or. &
         size( at ) /= nz+1 .or. size( bt ) /= nz+1 ) then
       write (*,'("ERROR - inconsistent dimensions in input")')
       write (*,'("      - fld2d (nx,ny,n):    ", 3i5)') &
            size(fld2d,1), size(fld2d,2), size(fld2d,3)
       write (*,'("      - fld3d (nx,ny,nz,n): ", 4i5)') &
            nx, ny, nz, n
       write (*,'("      - lli (nx,ny):        ", 2i5)') &
            lli%im, lli%jm
       write (*,'("      - at, bt (nz+1):      ", 2i5)') &
            size( at ), size( bt )
       TRACEBACK; status=1; return
    end if

    ! Allocations
    allocate( oro(nx,ny) )
    allocate( psurf(nx,ny) )
    allocate( logplow(nx,ny) )
    allocate( logphigh(nx,ny) )
    allocate( logpdown(nx,ny) )
    allocate( logpup(nx,ny) )

    ! Read surface altitude field
    call oro_ecmwf( lli, oro, status )
    if ( status /= 0 ) then
       write (*,'("WARNING - orography could not be read; assuming flat surface")')
       write (*,'("WARNING in ",a)') rname
       status = 0
       oro = 0.
    end if

    ! Surface pressure
    psurf = p0*exp( -c*oro )

    ! Height in meters or Pa
    select case( hunit )

    case( 'm' )

       ! Check on input
       if ( hlow < 0.0 .or. hhigh <= 0.0 .or. hhigh <= hlow ) then
          write (*,'("ERROR - input heights out of range")')
          write (*,'("ERROR - hlow, hhigh: ",2f10.3)') hlow, hhigh
          TRACEBACK; status=1; return
       end if

       ! defined wrt surface or mean sea level?
       select case( hdef )
       case( 'sfc' )
          hfac = 1
       case( 'msl' )
          hfac = 0
       case default
          write (*,'("ERROR - unknown height definition `",a,"`")') trim(hdef)
          TRACEBACK; status=1; return
       end select

       ! pressure at lower and higher level (logs for log interpolations)
       logplow  = log( p0 ) - c*(hfac*oro+hlow )
       logphigh = log( p0 ) - c*(hfac*oro+hhigh)

    case( 'p' )

       ! Check on input
       if ( hlow < 0.0 .or. hhigh <= 0.0 .or. hhigh >= hlow ) then
          write (*,'("ERROR - input pressures out of range")')
          write (*,'("ERROR - hlow, hhigh: ",2f10.3)') hlow, hhigh
          TRACEBACK; status=1; return
       end if

       ! log pressure at lower and higher level
       logplow  = alog( hlow )
       logphigh = alog( hhigh )

    case default

       write (*,'("ERROR - unknown height unit `",a,"`")') trim(hunit)
       TRACEBACK; status=1; return

    end select

    ! Checks of totals before and after
    tot2d = sum( fld2d )
    tot3db = sum( fld3d )
    ! Keep track of total below surface (possible in case hdef='msl')
    totbs = 0.

    ! First handle below surface
    logpup = alog( at(1) + bt(1)*psurf )
    do j = 1, ny
       do i = 1, nx
          ! Log of pressure at lower and higher level
          p2t = logplow(i,j)
          p2b = logphigh(i,j)
          ! Log of pressure at bottom and top of model layer
          p1t = alog( 2.e5 )   ! save below surface
          p1b = logpup(i,j)
          ! Determine fraction of field to be added in this layer
          f = max( ( min( p1t, p2t ) - max( p1b, p2b ) ) / ( p2t - p2b ), 0. )
          ! Add to output field in lowest layer
          fld3d(i,j,1,:) = fld3d(i,j,1,:) + f*fld2d(i,j,:)
          totbs = totbs + f*sum( fld2d(i,j,:) )
       end do
    end do

    ! loop over model layers (l outer loop is most efficient)
    do l = 2, nz
       logpdown = logpup
       logpup   = log( max( 1.0e-10, at(l+1) + bt(l+1)*psurf ) )
       do j = 1, ny
          do i = 1, nx
             ! Log of pressure at lower and higher level
             p2t = logplow(i,j)
             p2b = logphigh(i,j)
             ! Log of pressure at bottom and top of model layer
             p1t = logpdown(i,j)
             p1b = logpup(i,j)
             ! Determine fraction of field to be added in this layer
             f = max( ( min( p1t, p2t ) - max( p1b, p2b ) ) / ( p2t - p2b ), 0. )
             ! Add to output field
             fld3d(i,j,l,:) = fld3d(i,j,l,:) + f*fld2d(i,j,:)
          end do
       end do
    end do

    tot3da = sum( fld3d )

    if ( abs( (tot3da-tot3db)-tot2d ) > 1e-8*tot2d ) then
       write (*,'("ERROR - distribution of field not mass-conservative")')
       TRACEBACK; status=1; return
    end if

    if ( tot2d /= 0. .and. totbs /= 0. ) &
         write (*,'(a,": fraction of field below surface is: ",f9.6)') rname, totbs/tot2d

    deallocate( oro )
    deallocate( psurf )
    deallocate( logplow )
    deallocate( logphigh )
    deallocate( logpdown )
    deallocate( logpup )

    status = 0

  end subroutine dist_vertical



  subroutine dist_latlon( lli1, fld1, lli2, fld2, hcomb, status )
    !----------------------------------------------------------
    ! Redistribute field to different lat-lon grid
    ! Based on TM3 f77 routines 'fraction' and 'fill' but more general
    ! because target grid may be both coarser or finer than source grid
    !
    ! Input /output:
    !  o lli1, lli2       :  input/output lon-lat grid (see TMkit)
    !  o fld1, fld2       :  input/output fields
    !                         three dimensions: first two are lon-lat;
    !                         third is extra dimension, e.g. height
    !  o hcomb            :  horizontal combination key (see TMkit)
    !                         allowed values: 'sum' and 'aver'
    !      Note: this is somewhat complicated:
    !        - 'sum' means simply summing from fine to coarse grid, but
    !             it means 'area-weighted distribution' from coarse to fine grid
    !        - 'aver' means area-weighted averaging from fine to coarse grid, and
    !             it means 'uniform distribution' from coarse to fine grid
    !  o status           :  0=success; 1=failure
    !----------------------------------------------------------

    ! --- modules --------------------------------

    ! --- in/out ----------------------------------------------

    type(TllGridInfo), intent(in)  :: lli1
    real, intent(in)               :: fld1(:,:,:)
    type(TllGridInfo), intent(in)  :: lli2
    real, intent(out)              :: fld2(:,:,:)
    character(len=*), intent(in)   :: hcomb
    integer, intent(out)           :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: rname = mname//'/dist_latlon'

    ! --- local -----------------------------------------------

    integer                        :: nx1, ny1, nx2, ny2, n
    integer                        :: i1, i2, j1, j2
    real                           :: x1b, x1t, y1b, y1t, x2b, x2t, y2b, y2t
    real                           :: area, frac, ftot, tot1, tot2
    type(TllGridInfo)              :: slli1, slli2
    real, allocatable              :: fx(:,:)
    real, allocatable              :: fy(:,:)
    real, allocatable              :: tot(:)
    logical                        :: xswap, yswap

    ! --- begin -----------------------------------------------

    ! Dimensions
    nx1 = size( fld1, 1 )
    ny1 = size( fld1, 2 )
    n = size( fld1, 3 )
    nx2 = size( fld2, 1 )
    ny2 = size( fld2, 2 )

    ! storage:
    allocate( fx(nx1,nx2) )
    allocate( fy(ny1,ny2) )
    allocate( tot(n) )

    ! Check dimensions
    if ( lli1%im /= nx1 .or. lli1%jm /= ny1 ) then
       write (*,'("ERROR - input field and ll grid have different dimensions")')
       TRACEBACK; status=1; return
    end if
    if ( lli2%im /= nx2 .or. lli2%jm /= ny2 ) then
       write (*,'("ERROR - output field and ll grid have different dimensions")')
       TRACEBACK; status=1; return
    end if
    if ( size( fld2, 3 ) /= n ) then
       write (*,'("ERROR - third dimension of input and output fields differs")')
       TRACEBACK; status=1; return
    end if

    ! Define grid boundaries. Swap grids such that grid 2 is coarser than grid 1
    if ( nx2 <= nx1 ) then
       xswap = .false.
       if ( ny2 <= ny1 ) then
          yswap = .false.
          ! Note: simply slli1=lli1 etc is NOT allowed here!
          call Init( slli1, lli1%lon_deg(1), lli1%dlon_deg, lli1%nlon, &
                            lli1%lat_deg(1), lli1%dlat_deg, lli1%nlat, status )
          IF_NOTOK_RETURN(status=1)
          call Init( slli2, lli2%lon_deg(1), lli2%dlon_deg, lli2%nlon, &
                            lli2%lat_deg(1), lli2%dlat_deg, lli2%nlat, status )
          IF_NOTOK_RETURN(status=1)
       else
          yswap = .true.
          call Init( slli1, lli1%lon_deg(1), lli1%dlon_deg, lli1%nlon, &
                            lli2%lat_deg(1), lli2%dlat_deg, lli2%nlat, status )
          IF_NOTOK_RETURN(status=1)
          call Init( slli2, lli2%lon_deg(1), lli2%dlon_deg, lli2%nlon, &
                            lli1%lat_deg(1), lli1%dlat_deg, lli1%nlat, status )
          IF_NOTOK_RETURN(status=1)
       end if
    else
       xswap = .true.
       if ( ny2 <= ny1 ) then
          yswap = .false.
          call Init( slli1, lli2%lon_deg(1), lli2%dlon_deg, lli2%nlon, &
                            lli1%lat_deg(1), lli1%dlat_deg, lli1%nlat, status )
          IF_NOTOK_RETURN(status=1)
          call Init( slli2, lli1%lon_deg(1), lli1%dlon_deg, lli1%nlon, &
                            lli2%lat_deg(1), lli2%dlat_deg, lli2%nlat, status )
          IF_NOTOK_RETURN(status=1)
       else
          yswap = .true.
          call Init( slli1, lli2%lon_deg(1), lli2%dlon_deg, lli2%nlon, &
                            lli2%lat_deg(1), lli2%dlat_deg, lli2%nlat, status )
          IF_NOTOK_RETURN(status=1)
          call Init( slli2, lli1%lon_deg(1), lli1%dlon_deg, lli1%nlon, &
                            lli1%lat_deg(1), lli1%dlat_deg, lli1%nlat, status )
          IF_NOTOK_RETURN(status=1)
       end if
    end if
    if ( status /= 0. ) then
       write (*,'("ERROR - failed to swap lat-lon grids")')
       TRACEBACK; status=1; return
    end if

    ! Calculate fraction of longitude cell ix1 contained in cell ix2
    do i2 = 1, slli2%im
       x2b = slli2%blon_deg(i2-1)
       x2t = slli2%blon_deg(i2)
       do i1 = 1, slli1%im
          x1b = slli1%blon_deg(i1-1)
          x1t = slli1%blon_deg(i1)
          ! bring longitudes within +-180 deg from each other
          if ( x1b - x2b < -180. ) then
             x1b = x1b + 360.
             x1t = x1t + 360.
          end if
          if ( x1b - x2b > 180. ) then
             x1b = x1b - 360.
             x1t = x1t - 360.
          end if
          frac = max( ( min( x1t, x2t ) - max( x1b, x2b ) ) / ( x2t - x2b ), 0. )
          if ( xswap ) then
             fx(i2,i1) = frac
          else
             fx(i1,i2) = frac
          end if
       end do
    end do

    ! Check sum of fractions (only makes sense when grids span same domain)
    do i2 = 1, slli2%im
       if ( xswap ) then
          ftot = sum( fx(i2,:) )
       else
          ftot = sum( fx(:,i2) )
       end if
       if ( abs( ftot - 1. ) > 1.e-5 ) then
          !          write (*,'("ERROR - sum of longitude fractions /= 1")')
          !          write(*,*) i2, ftot
          !          TRACEBACK; status=1; return
       end if
    end do

    ! Calculate fraction of latitude cell iy1 contained in cell iy2
    do j2 = 1, slli2%jm
       y2b = slli2%blat_deg(j2-1)
       y2t = slli2%blat_deg(j2)
       do j1 = 1, slli1%jm
          y1b = slli1%blat_deg(j1-1)
          y1t = slli1%blat_deg(j1)
          frac = max( ( min( y1t, y2t ) - max( y1b, y2b ) ) / ( y2t - y2b ), 0. )
          if ( yswap ) then
             fy(j2,j1) = frac
          else
             fy(j1,j2) = frac
          end if
       end do
    end do

    ! Check sum of fractions (only makes sense when grids span same domain)
    do j2 = 1, slli2%jm
       if ( xswap ) then
          ftot = sum( fy(j2,:) )
       else
          ftot = sum( fy(:,j2) )
       end if
       if ( abs( ftot - 1. ) > 1.e-5 ) then
          !          write (*,'("ERROR - sum of latitude fractions /= 1")')
          !          write(*,*) j2, ftot
          !          TRACEBACK; status=1; return
       end if
    end do

    ! Redistribute fld1 on grid 1 to fld2 on grid 2
    ! Take latitude-dependent grid box area into account
    fld2 = 0.

    select case ( hcomb )
    case( 'sum' )
       do j2 = 1, ny2
          do i2 = 1, nx2
             tot = 0.
             area = 0.
             do j1 = 1, ny1
                if ( fy(j1,j2) == 0. ) cycle
                do i1 = 1, nx1
                   if ( fx(i1,i2) == 0. ) cycle
                   area = area + fx(i1,i2)*fy(j1,j2)*lli1%area_m2(j1)
                   tot = tot + fx(i1,i2)*fy(j1,j2)*fld1(i1,j1,:)
                end do
             end do
             fld2(i2,j2,:) = tot/area*lli2%area_m2(j2)
          end do
       end do
    case( 'aver' )
       do j2 = 1, ny2
          do i2 = 1, nx2
             tot = 0.
             area = 0.
             do j1 = 1, ny1
                if ( fy(j1,j2) == 0. ) cycle
                do i1 = 1, nx1
                   if ( fx(i1,i2) == 0. ) cycle
                   area = area + fx(i1,i2)*fy(j1,j2)*lli1%area_m2(j1)
                   tot = tot + fx(i1,i2)*fy(j1,j2)*fld1(i1,j1,:)*lli1%area_m2(j1)
                end do
             end do
             fld2(i2,j2,:) = tot/area
          end do
       end do
    case default
       write (*,'("ERROR - unsupported horizonal combination key:",a)') hcomb
       TRACEBACK; status=1; return
    end select

    !>>> AJS: segmentation faults, but check was not applied in the end
    !! Check totals/averages (only makes sense if both grids span same domain)
    !! Small differences can occur because of different surface areas grids 1 and 2
    !select case( hcomb )
    !  case( 'sum' )
    !   tot1 = sum( fld1 )
    !   tot2 = sum( fld2 )
    !   !       write (*,'(a,": Sum of input field  = ",e15.8)') name, tot1
    !   !       write (*,'(a,": Sum of output field = ",e15.8)') name, tot2
    !  case( 'aver' )
    !   tot1 = sum(sum(sum(fld1,dim=3),dim=1)*lli1%area_m2)/(nx1*sum(lli1%area_m2)*n)
    !   tot2 = sum(sum(sum(fld2,dim=3),dim=1)*lli2%area_m2)/(nx2*sum(lli2%area_m2)*n)
    !   !       write (*,'(a,": Mean of input field  = ",e15.8)') name, tot1
    !   !       write (*,'(a,": Mean of output field = ",e15.8)') name, tot2
    !end select
    !if ( abs( tot2-tot1 ) > 1.e-3*tot1 ) then
    !   !       write (*,'("WARNING - Difference in-out larger than 0.1%")')
    !   !       write (*,'("WARNING in ",a)') name
    !end if
    !<<<

    ! Deallocate swapped grids
    call Done( slli1, status )
    IF_NOTOK_RETURN(status=1)
    call Done( slli2, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( fx )
    deallocate( fy )
    deallocate( tot )

    status = 0

  end subroutine dist_latlon


!  subroutine intp_3d( fldX, lliX, leviX, fld, lli, levi, vtype, status )
!    !----------------------------------------------------------
!    ! Do a 3D interpolation from one lat-lon / hybrid grid
!    !  to the other.
!    ! In horizontal direction use bi-linear interpolation.
!    ! In vertical direction use either linear or hermite cubic interpolation.
!    ! Works also for the case nx(X)=1 or ny(X)=1.
!    !
!    ! To calculate pressures on hybrid grids, get surface pressure
!    !  field from ECMWF orography, similarly as in dist_vertical
!    !
!    !  o fldX(nxX,nyX,nzX,n)  : source field [in]
!    !                             'n' is extra dimension to process multiple fields
!    !                             (e.g. months) at once
!    !  o lliX, leviX          : latlon, hybrid source grids [in]
!    !  o fld(nx,ny,nz,n)      : target field [out]
!    !  o lli, levi            : latlon, hybrid target grids [in]
!    !  o vtype                : type of vertical interpolation [in]
!    !                             'lin' for linear, 'muherm' for hermite cubic
!    !  o status               : 0=success; 1=failure
!    !----------------------------------------------------------
!
!    ! --- modules --------------------------------
!
!    use binas,                 only : T0, p0, xm_air, grav, Rgas
!    use num,                   only : Interp_lin, CircInterp_lin, Interp_muherm
!
!    ! --- in/out ----------------------------------------------
!
!    real, intent(in)               :: fldX(:,:,:,:)
!    type(TllGridInfo), intent(in)  :: lliX
!    type(TLevelInfo), intent(in)   :: leviX
!    real, intent(out)              :: fld(:,:,:,:)
!    type(TllGridInfo), intent(in)  :: lli
!    type(TLevelInfo), intent(in)   :: levi
!    character(len=*), intent(in)   :: vtype
!    integer, intent(out)           :: status
!
!    ! --- const -----------------------------------------------
!
!    character(len=*), parameter    :: rname = mname//'/intp_3d'
!    real, parameter                :: Tav = T0   ! average temperature
!    real, parameter                :: c = xm_air*grav/(Rgas*Tav) ! a constant
!
!    ! --- local -----------------------------------------------
!
!    integer                        :: nxX, nyX, nzX, nx, ny, nz, n
!    integer                        :: i, ix, iy, iz
!    integer                        :: lastx, lasty, lastz
!    logical                        :: zswapX, zswap
!    real, dimension(:,:), allocatable     :: oroX, psurfX, oro, psurf
!    real, dimension(:), allocatable       :: faX, fbX, fa, fb
!    real, dimension(:), allocatable       :: pres3X
!    real, dimension(:,:), allocatable     :: pres2X, fld3X
!    real, dimension(:,:,:), allocatable   :: pres1X, fld2X, pres1
!    real, dimension(:,:,:,:), allocatable :: fld1X, fld1
!
!    ! --- begin -----------------------------------------------
!
!    ! Dimensions
!    nxX = size( fldX, 1 )
!    nyX = size( fldX, 2 )
!    nzX = size( fldX, 3 )
!    n = size( fldX, 4 )
!    nx = size( fld, 1 )
!    ny = size( fld, 2 )
!    nz = size( fld, 3 )
!
!    ! Allocations
!    allocate( oroX(nxX,nyX) )
!    allocate( psurfX(nxX,nyX) )
!    allocate( oro(nx,ny) )
!    allocate( psurf(nx,ny) )
!    allocate( faX(nzX) )
!    allocate( fbX(nzX) )
!    allocate( fa(nz) )
!    allocate( fb(nz) )
!    allocate( pres3X(nzX) )
!    allocate( fld3X(n,nzX) )
!    allocate( pres2X(nzX,nyX) )
!    allocate( fld2X(n,nzX,nyX) )
!    allocate( pres1X(nzX,nyX,nxX) )
!    allocate( fld1X(n,nzX,nyX,nxX) )
!    allocate( pres1(nz,ny,nx) )
!    allocate( fld1(n,nz,ny,nx) )
!
!    ! Read surface altitude fields
!    call oro_ecmwf( lliX, oroX, status )
!    if ( status /= 0 ) then
!       write (*,'("WARNING - orography source grid could not be read; assuming flat surface")')
!       write (*,'("WARNING in ",a)') rname
!       status = 0
!       oroX = 0.
!    end if
!    call oro_ecmwf( lli, oro, status )
!    if ( status /= 0 ) then
!       write (*,'("WARNING - orography source grid could not be read; assuming flat surface")')
!       write (*,'("WARNING in ",a)') rname
!       status = 0
!       oro = 0.
!    end if
!
!    ! Surface pressure fields
!    psurfX = p0*exp( -c*oroX )
!    psurf = p0*exp( -c*oro )
!
!    ! Determine whether vertical grids must be swapped
!    ! (for muherm pressure should be in increasing order)
!    zswapX = .false.
!    faX = leviX%fa
!    fbX = leviX%fb
!    if ( faX(1)+p0*fbX(1) > faX(nzX)+p0*fbX(nzX) ) then
!       zswapX = .true.
!       faX(1:nzX) = faX(nzX:1:-1)
!       fbX(1:nzX) = fbX(nzX:1:-1)
!    end if
!
!    zswap = .false.
!    fa = levi%fa
!    fb = levi%fb
!    if ( fa(1)+p0*fb(1) > fa(nz)+p0*fb(nz) ) then
!       zswap = .true.
!       fa(1:nz) = fa(nz:1:-1)
!       fb(1:nz) = fb(nz:1:-1)
!    end if
!
!    write(*,*) 'zswapX', zswapX
!    write(*,*) 'faX', faX
!    write(*,*) 'fbX', fbX
!    write(*,*) 'zswap ', zswap
!    write(*,*) 'fa', fa
!    write(*,*) 'fb', fb
!
!    ! Pressure on full grids (revert order for use in Interp_lin)
!    do ix = 1, nxX
!       do iy = 1, nyX
!          pres1X(:,iy,ix) = faX + psurfX(ix,iy)*fbX
!       end do
!    end do
!    do ix = 1, nx
!       do iy = 1, ny
!          pres1(:,iy,ix) = fa + psurf(ix,iy)*fb
!       end do
!    end do
!
!    ! Do the interpolation
!
!    ! Revert order of fldX for use in linear interpolation routines
!    ! and possibly swap vertical direction !
!    do ix = 1, nxX
!       do iy = 1, nyX
!          do iz = 1, nzX
!             do i = 1, n
!                if ( zswapX ) then
!                   fld1X(i,iz,iy,ix) = fldX(ix,iy,nzX+1-iz,i)
!                else
!                   fld1X(i,iz,iy,ix) = fldX(ix,iy,iz,i)
!                end if
!             end do
!          end do
!       end do
!    end do
!
!    lastx = 1
!
!    do ix = 1, nx
!
!       if ( nx /= 1 ) then
!
!          ! Interpolate to this longitude
!          call CircInterp_Lin( lliX%lon_deg, 360., fld1X, lli%lon_deg(ix), fld2X, lastx, status )
!          if ( status /= 0. ) then
!             TRACEBACK; status=1; return
!          end if
!
!          call CircInterp_Lin( lliX%lon_deg, 360., pres1X, lli%lon_deg(ix), pres2X, lastx, status )
!          if ( status /= 0. ) then
!             TRACEBACK; status=1; return
!          end if
!
!       else  ! 'interpolation' does not work in case of just one longitude
!
!          fld2X = fld1X(:,:,:,1)
!          pres2X = pres1X(:,:,1)
!
!       end if
!
!       lasty = 1
!
!       do iy = 1, ny
!
!          if ( ny /= 1 ) then
!
!             ! Interpolate to this latitude
!             call Interp_lin( lliX%lat_deg, fld2X, lli%lat_deg(iy), fld3X, lasty, status )
!             if ( status /= 0. ) then
!                TRACEBACK; status=1; return
!             end if
!
!             call Interp_lin( lliX%lat_deg, pres2X, lli%lat_deg(iy), pres3X, lasty, status )
!             if ( status /= 0. ) then
!                TRACEBACK; status=1; return
!             end if
!
!          else  ! 'interpolation' does not work in case of just one latitude
!
!             fld3X = fld2X(:,:,1)
!             pres3X = pres2X(:,1)
!
!          end if
!
!          ! Interpolate to target pressure grid
!
!          select case ( vtype )
!
!          case ( 'lin' )
!
!             lastz = 1
!             do iz = 1, nz
!                call Interp_lin( pres3X, fld3X, pres1(iz,iy,ix), fld1(:,iz,iy,ix), &
!                     lastz, status )
!                if ( status /= 0. ) then
!                   TRACEBACK; status=1; return
!                end if
!             end do
!
!          case ( 'muherm' )
!
!             ! (no 'extended' version of muherm implemented yet)
!             do i = 1, n
!                call Interp_muherm( pres3X, fld3X(i,:), pres1(:,iy,ix), &
!                     fld1(i,:,iy,ix), status )
!                if ( status /= 0. ) then
!                   TRACEBACK; status=1; return
!                end if
!             end do ! i
!
!          case default
!
!             if ( status /= 0 ) then
!                write (*,'("ERROR - unknown vertical interpolation type: ",a)') trim(vtype)
!                TRACEBACK; status=1; return
!             end if
!
!          end select
!
!       end do ! iy
!
!    end do ! ix
!
!    ! Copy to output field and possibly swap vertical direction of target field
!    do ix = 1, nx
!       do iy = 1, ny
!          do iz = 1, nz
!             do i = 1, n
!                if ( zswap ) then
!                   fld(ix,iy,iz,i) = fld1(i,nz+1-iz,iy,ix)
!                else
!                   fld(ix,iy,iz,i) = fld1(i,iz,iy,ix)
!                end if
!             end do
!          end do
!       end do
!    end do
!
!    ! Deallocate
!    deallocate( oroX )
!    deallocate( psurfX )
!    deallocate( oro )
!    deallocate( psurf )
!    deallocate( faX )
!    deallocate( fbX )
!    deallocate( fa )
!    deallocate( fb )
!    deallocate( pres3X )
!    deallocate( fld3X )
!    deallocate( pres2X )
!    deallocate( fld2X )
!    deallocate( pres1X )
!    deallocate( fld1X )
!    deallocate( pres1 )
!    deallocate( fld1 )
!
!    status = 0
!
!  end subroutine intp_3d


  subroutine oro_ecmwf( lli, oro, status )

    !----------------------------------------------------------
    ! Read ECMWF model orography on 1x1 deg global resolution,
    ! and convert to requested resolution (lli)
    !
    ! Input:
    !   o lli    - latlon grid
    !
    ! Output:
    !   o oro    - surface altitude on lli grid (in meters)
    !   o status
    !----------------------------------------------------------

    ! --- modules ---------------------------------------------

    use file_hdf
    use TIPP_Settings, only : tipp_ecmwf_dir
    use os_specs,      only : MAX_FILENAME_LEN

    ! --- in/out ----------------------------------------------

    type(TllGridInfo), intent(in)  :: lli
    real, intent(out)              :: oro(:,:)
    integer, intent(out)           :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: rname = mname//'/oro_ecmwf'

    character(len=*), parameter    :: sdsname = 'oro'
    integer, parameter             :: nx_f = 360, ny_f = 180

    ! --- local -----------------------------------------------

    integer                        :: nx, ny
    real                           :: field0(nx_f,ny_f)
    real                           :: field1(nx_f,ny_f,1)
    real, allocatable              :: field2(:,:,:)
    type(THdfFile)                 :: hdf
    type(TSds)                     :: sds
    type(TllGridInfo)              :: lliX
    character(len=MAX_FILENAME_LEN):: hdfname

    ! --- begin -----------------------------------------------

    ! dummy ...
    oro = 0.0

    ! storage:
    allocate( field2(size(oro,1),size(oro,2),1) )

    ! Shortcuts
    nx = size( oro, 1 )
    ny = size( oro, 2 )

    if ( lli%im /= nx .or. lli%jm /= ny ) then
       write (*,'("ERROR - inconsistent dimensions in input")')
       write (*,'("      - oro (nx,ny):     ", 2i5)') nx, ny
       write (*,'("      - lli (nx,ny):     ", 2i5)') lli%im, lli%jm
       TRACEBACK; status=1; return
    endif

    ! full path:
    write (hdfname,'(a,"/oro-glb1x1.hdf")') trim(tipp_ecmwf_dir)

    ! Read orography from HDF file
    call Init( hdf, trim(hdfname), 'read', status )
    IF_NOTOK_RETURN(status=1)

    call Init( sds, hdf, sdsname, status )
    IF_NOTOK_RETURN(status=1)

    call ReadData( sds, field0, status )
    IF_NOTOK_RETURN(status=1)
    field1(:,:,1) = field0

    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    call Done( hdf, status )
    IF_NOTOK_RETURN(status=1)

    ! Latlon grid in orography file
    call Init( lliX, -180.0+180.0/nx_f, 360.0/nx_f, nx_f, &
                      -90.0+ 90.0/ny_f, 180.0/ny_f, ny_f, status )

    ! Convert to target latlon resolution
    call dist_latlon( lliX, field1, lli , field2, 'aver', status )
    IF_NOTOK_RETURN(status=1)

    ! copy result:
    oro = field2(:,:,1)

    ! Finish grid
    call Done( lliX, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( field2 )

    ! ok
    status = 0

  end subroutine oro_ecmwf

end module tipp_grid_tools
