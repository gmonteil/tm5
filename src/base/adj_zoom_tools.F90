!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_zoom_tools

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- interface ---

  private

  public :: adj_mix_edges, adj_update_parent

  ! --- const ---

  character(len=*), parameter  ::  mname = 'adj_zoom_tools'


contains

!===========================================================================================================
!===========================================================================================================

  !=============================================================
  ! - adj_mix_edges is self-adjoint
  ! - however: reconstruction of required varialbles necessary
  !=============================================================

  subroutine adj_mix_edges( region, status )

    use zoom_tools, only:  mix_edges

    !__IO____________________________________________________________________

    integer,intent(in)            ::  region
    integer,intent(out)           ::  status

    !__CONST_________________________________________________________________

    character(len=*), parameter  ::  rname = mname//'adj_mix_edges'

    !__LOCAL_VARIABLES_______________________________________________________


    !__START_SUBROUTINE______________________________________________________

    call mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    call reconstruct_required_var_1( region, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_mix_edges

!===========================================================================================================
!===========================================================================================================

  !
  ! "A": exact reconstruction of the forward air masses in the edge (and halo) cells
  ! of a zoom region, for the advection step that is about to be adjointed.
  !
  ! Called from adj_mix_edges at the start of an adjoint x/y/z advection step of
  ! 'region'. On exit, the edge cells (region_dat%edge == -1) and the halo cells
  ! (0, im+1, 0, jm+1) of m_dat(region)%data hold the forward air masses as they
  ! were just after the forward advection loop of that step, before mix_edges
  ! (for x and y: after uncompress). The core cells are not touched: the adjoint
  ! already has them exactly.
  !
  ! Method. Within one parent half step ("block" of n_operators*tref_ child steps),
  ! the forward air mass of the child edge cells is fully determined by:
  !   o the parent put_xedges/put_yedges, done before the child block starts,
  !     from the parent air masses at that moment. Those are recovered from the
  !     current parent air mass by undoing the parent x/y/z mass fluxes:
  !       parent order xyz : put_x before x, put_y after x  (then z)
  !       parent order zyx : put_y after z, put_x after z,y
  !     (during the adjoint of the child block, the parent air mass is the forward
  !     one at the end of its advection triple);
  !   o the child advection steps of the block up to the current one, applied to
  !     the edge cells only, with the same rules as the forward code: scopes,
  !     compress/uncompress (x,y), zeroing of the outer mass fluxes, and
  !     mix_edges after each step.
  ! Mass changes only depend on the mass fluxes, never on neighbouring masses, so
  ! the core cells are not needed (they are kept at zero in the work array).
  !
  ! Assumptions, checked at run time (the adjoint stops with an error otherwise;
  ! these configurations could not be tested):
  !   o the zoom region is not cyclic / 360 degrees wide, does not touch a pole,
  !     does not cross the date line, and is not refined vertically (zref = 1);
  !   o no reduced grid in the zoom region, nor in the parent rows around it
  !     (see below);
  !   o the parent cells around the zoom region (one cell beyond its edges) are
  !     in the core of the parent: not in the parent's own zoom edges, not in
  !     the pole rows or across the date line of the global region.
  !
  ! Reduced grid (region.<name>.redgrid.* in the rc file):
  !   o SUPPORTED: a reduced grid in the global (parent) region, as long as it
  !     only covers latitudes away from the zoom region, i.e. the parent rows
  !     jbeg-1 ... jend+1 around the zoom region are NOT reduced. This is the
  !     usual case: the reduced grid is only applied close to the poles.
  !   o NOT SUPPORTED (run stops with an error):
  !       - a reduced grid in a zoom region itself (nred(region) /= 0): its
  !         x-advection then works on merged cells, not on the cells assumed here;
  !       - a reduced parent row next to (or under) the zoom region, i.e. the zoom
  !         region extends to latitudes where the parent grid is reduced: the
  !         parent air-mass change in such a row is not the plain flux difference
  !         used below to recover the parent masses at put_xedges/put_yedges.
  !     In both cases, either remove the reduced grid in the rows concerned, or
  !     move the zoom region away from those latitudes.
  !
  ! Date line: zoom regions crossing the date line (ibeg > iend) are not supported.
  ! The forward model does not support them consistently either: the 'zoomed' mask
  ! built in initexit.F90 (blank_zoom_region) marks no parent cell as covered when
  ! ibeg > iend, so the parent and the child would both apply the physics there.
  !

  subroutine reconstruct_required_var_1( region, status )

    use dims,        only : im, jm, lm
    use dims,        only : parent, xref, yref, zref, tref
    use dims,        only : ibeg, iend, jbeg, jend
    use dims,        only : xcyc, touch_sp, touch_np, zoom2D
    use dims,        only : nsplitsteps, n_operators, splitorderzoom
    use dims,        only : rstatus => status
    use global_data, only : wind_dat, region_dat
    use MeteoData  , only : m_dat
    use redgridZoom, only : grid_reduced, nred, imred

    implicit none

    !__IO____________________________________________________________________

    integer, intent(in)              ::  region
    integer, intent(out)             ::  status

    !__CONST_________________________________________________________________

    character(len=*), parameter      ::  rname = mname//'/reconstruct_required_var_1'

    !__LOCAL_VARIABLES_______________________________________________________

    real, dimension(:,:,:), pointer  ::  mc, mp, amp, bmp, cmp
    integer, dimension(:,:), pointer ::  edge
    ! work array: air mass of the child in its edge and halo cells
    real, allocatable                ::  mass_edge(:,:,:)
    ! parent air mass at the moment the parent copies it into the child's
    ! west/east edge cells (put_xedges) and south/north edge cells (put_yedges):
    real, allocatable                ::  parent_mass_x_edge(:,:,:)
    real, allocatable                ::  parent_mass_y_edge(:,:,:)
    ! change of the parent air mass in its own x, y and z advection steps:
    real, allocatable                ::  parent_mass_change_x(:,:,:)
    real, allocatable                ::  parent_mass_change_y(:,:,:)
    real, allocatable                ::  parent_mass_change_z(:,:,:)
    integer                          ::  my_parent, xref_, yref_, zref_, tref_
    integer                          ::  imc, jmc, lmc, nblk, p, p0, k
    integer                          ::  i, j, l, ip, jp
    integer                          ::  ib, ie_, jb, je_
    character(len=1)                 ::  op
    character(len=n_operators)       ::  stencil
    logical                          ::  xyz
    integer                          ::  grandparent, xrefp, yrefp
    integer                          ::  ipmin, ipmax, jpmin, jpmax

    !__START_SUBROUTINE______________________________________________________

    status = 0

    my_parent = parent(region)
    if ( region == 1 .or. my_parent == 0 ) return

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    zref_ = zref(region)/zref(my_parent)
    tref_ = tref(region)/tref(my_parent)

    ib = ibeg(region) ; ie_ = iend(region)
    jb = jbeg(region) ; je_ = jend(region)

    ! --- check the assumptions (configurations that could not be tested) ---

    if ( (xcyc(region) == 1) .or. (im(region)/xref(region) == im(1)) ) then
      write (gol,'("zoom region ",i0," is cyclic / 360 degrees wide: not supported")') region; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (touch_sp(region) == 1) .or. (touch_np(region) == 1) ) then
      write (gol,'("zoom region ",i0," touches a pole: not supported")') region; call goErr
      TRACEBACK; status=1; return
    end if
    if ( ib > ie_ ) then
      write (gol,'("zoom region ",i0," crosses the date line: not supported")') region; call goErr
      TRACEBACK; status=1; return
    end if
    if ( zref_ /= 1 ) then
      write (gol,'("zoom region ",i0," is refined vertically (zref=",i0,"): not supported")') region, zref_; call goErr
      TRACEBACK; status=1; return
    end if
    ! reduced grid in the zoom region itself (see "Reduced grid" in the header):
    if ( grid_reduced .and. (nred(region) /= 0) ) then
      write (gol,'("zoom region ",i0," uses a reduced grid: not supported")') region; call goErr
      write (gol,'("  remove the reduced grid for this region (region.<name>.redgrid.* in the rc file)")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! parent cells used: one cell beyond the zoom region on each side;
    ! they must be in the core of the parent, advected over the full scope:
    grandparent = parent(my_parent)
    if ( grandparent == 0 ) then
      ! global parent: no wrap around the date line, no pole rows
      ipmin = 1 ; ipmax = im(my_parent)
      jpmin = 2 ; jpmax = jm(my_parent)-1
    else
      ! zoom parent: outside its own edge cells
      xrefp = xref(my_parent)/xref(grandparent)
      yrefp = yref(my_parent)/yref(grandparent)
      ipmin = xrefp+1 ; ipmax = im(my_parent)-xrefp
      jpmin = yrefp+1 ; jpmax = jm(my_parent)-yrefp
    end if
    if ( (ib-1 < ipmin) .or. (ie_+1 > ipmax) .or. (jb-1 < jpmin) .or. (je_+1 > jpmax) ) then
      write (gol,'("zoom region ",i0," is too close to the edge of its parent ",i0,": not supported")') &
             region, my_parent; call goErr
      write (gol,'("  parent cells used : i ",i0,":",i0,", j ",i0,":",i0)') ib-1, ie_+1, jb-1, je_+1; call goErr
      write (gol,'("  parent core       : i ",i0,":",i0,", j ",i0,":",i0)') ipmin, ipmax, jpmin, jpmax; call goErr
      TRACEBACK; status=1; return
    end if
    ! reduced parent rows around the zoom region (see "Reduced grid" in the header);
    ! a reduced grid elsewhere in the parent (e.g. near the poles) is fine:
    if ( grid_reduced .and. (nred(my_parent) /= 0) ) then
      do jp = jb-1, je_+1
        if ( jp > size(imred,1) ) cycle
        if ( imred(jp,my_parent) /= im(my_parent) ) then
          write (gol,'("parent ",i0," of zoom region ",i0," has a reduced grid in row ",i0,": not supported")') &
                 my_parent, region, jp; call goErr
          write (gol,'("  the parent rows ",i0,":",i0," around the zoom region must not be reduced;")') &
                 jb-1, je_+1; call goErr
          write (gol,'("  limit the reduced grid of the parent to other latitudes, or move the zoom region")'); call goErr
          TRACEBACK; status=1; return
        end if
      end do
    end if

    imc = im(region) ; jmc = jm(region) ; lmc = lm(region)

    mc   => m_dat(region)%data
    mp   => m_dat(my_parent)%data
    amp  => wind_dat(my_parent)%am_t
    bmp  => wind_dat(my_parent)%bm_t
    cmp  => wind_dat(my_parent)%cm_t
    edge => region_dat(region)%edge

    ! forward position of the step being adjointed (adjoint counter is one ahead),
    ! and first position of the child block (= one parent half step):
    p    = rstatus(region) - 1
    nblk = n_operators * tref_
    p0   = (p/nblk) * nblk

    ! order of the parent advection triple = order of the first sub-block of the child block:
    do i = 1, n_operators
      stencil(i:i) = splitorderzoom(region, p0+i)
    end do
    xyz = ( index(stencil,'x') < index(stencil,'z') )

    ! --- parent air masses at the moments of put_xedges / put_yedges ---

    allocate( parent_mass_change_x(ib-1:ie_+1, jb-1:je_+1, lmc) )
    allocate( parent_mass_change_y(ib-1:ie_+1, jb-1:je_+1, lmc) )
    allocate( parent_mass_change_z(ib-1:ie_+1, jb-1:je_+1, lmc) )
    allocate( parent_mass_x_edge(ib-1:ie_+1, jb-1:je_+1, lmc) )
    allocate( parent_mass_y_edge(ib-1:ie_+1, jb-1:je_+1, lmc) )

    do l = 1, lmc
      do jp = jb-1, je_+1
        do ip = ib-1, ie_+1
          ! mass change of the parent cell in its x, y and z advection steps:
          parent_mass_change_x(ip,jp,l) = amp(ip-1,jp,l) - amp(ip,jp,l)
          parent_mass_change_y(ip,jp,l) = bmp(ip,jp,l)   - bmp(ip,jp+1,l)
          parent_mass_change_z(ip,jp,l) = cmp(ip,jp,l-1) - cmp(ip,jp,l)
        end do
      end do
    end do

    if ( xyz ) then
      ! parent: put_x | x | put_y | y | z  ; now = after z
      parent_mass_x_edge = mp(ib-1:ie_+1, jb-1:je_+1, 1:lmc) - parent_mass_change_z - parent_mass_change_y - parent_mass_change_x
      parent_mass_y_edge = mp(ib-1:ie_+1, jb-1:je_+1, 1:lmc) - parent_mass_change_z - parent_mass_change_y
    else
      ! parent: z | put_y | y | put_x | x  ; now = after x
      parent_mass_y_edge = mp(ib-1:ie_+1, jb-1:je_+1, 1:lmc) - parent_mass_change_x - parent_mass_change_y
      parent_mass_x_edge = mp(ib-1:ie_+1, jb-1:je_+1, 1:lmc) - parent_mass_change_x
    end if

    ! --- child edge air masses just after the puts (same formulas as put_x/yedges) ---

    allocate( mass_edge(0:imc+1, 0:jmc+1, 1:lmc) )
    mass_edge = 0.0

    if ( xyz ) then
      call set_xedges_from_parent( region, parent_mass_x_edge, mass_edge )
      call set_yedges_from_parent( region, parent_mass_y_edge, mass_edge )
    else
      call set_yedges_from_parent( region, parent_mass_y_edge, mass_edge )
      call set_xedges_from_parent( region, parent_mass_x_edge, mass_edge )
    end if

    ! --- child advection steps of the block, up to and including the current one ---

    do k = p0, p
      op = splitorderzoom(region, k+1)
      select case ( op )
        case ( 'x' )
          call advect_edges_x( region, k, mass_edge )
        case ( 'y' )
          call advect_edges_y( region, k, mass_edge )
        case ( 'z' )
          call advect_edges_z( region, k, mass_edge )
        case default
          cycle   ! no air-mass change
      end select
      ! mix_edges follows every advection step; for the current step the
      ! adjoint needs the state before it:
      if ( k < p ) call mix_edges_mass( region, mass_edge )
    end do

    ! --- store edge and halo cells ---

    do l = 1, lmc
      do j = 1, jmc
        do i = 1, imc
          if ( edge(i,j) == -1 ) mc(i,j,l) = mass_edge(i,j,l)
        end do
      end do
      mc(0      ,1:jmc,l) = mass_edge(0      ,1:jmc,l)
      mc(imc+1  ,1:jmc,l) = mass_edge(imc+1  ,1:jmc,l)
      mc(1:imc  ,0    ,l) = mass_edge(1:imc  ,0    ,l)
      mc(1:imc  ,jmc+1,l) = mass_edge(1:imc  ,jmc+1,l)
    end do

    deallocate( mass_edge, parent_mass_x_edge, parent_mass_y_edge )
    deallocate( parent_mass_change_x, parent_mass_change_y, parent_mass_change_z )
    nullify( mc, mp, amp, bmp, cmp, edge )

  end subroutine reconstruct_required_var_1


  !
  ! Helper routines of reconstruct_required_var_1. They act on the work array
  ! mass_edge(0:im+1,0:jm+1,1:lm) holding the air mass of the edge and halo cells of zoom
  ! region 'region' (core cells are not used), and reproduce, for the air mass
  ! only, what the forward routines named in their comments do.
  !

  ! put_xedges (advectx.F90): halo columns and x-wall cells from the parent air
  ! mass at the moment of the put (parent_mass_x_edge, parent indices)

  subroutine set_xedges_from_parent( region, parent_mass_x_edge, mass_edge )

    use dims, only : im, jm, lm, parent, xref, yref, zref
    use dims, only : ibeg, iend, jbeg

    integer, intent(in)     ::  region
    real, intent(in)        ::  parent_mass_x_edge(ibeg(region)-1:,jbeg(region)-1:,:)
    real, intent(inout)     ::  mass_edge(0:,0:,:)

    integer   ::  imc, jmc, lmc, xref_, yref_, zref_, ib, ie_, jb
    integer   ::  j, l, jp

    imc = im(region) ; jmc = jm(region) ; lmc = lm(region)
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    ib = ibeg(region) ; ie_ = iend(region) ; jb = jbeg(region)

    do l = 1, lmc
      do j = 1, jmc
        jp = jb + (j-1)/yref_
        mass_edge(0    ,j,l) = parent_mass_x_edge(ib-1 ,jp,l) / real(yref_*zref_)
        mass_edge(imc+1,j,l) = parent_mass_x_edge(ie_+1,jp,l) / real(yref_*zref_)
        mass_edge(1:xref_        ,j,l) = parent_mass_x_edge(ib ,jp,l) / real(xref_*yref_*zref_)
        mass_edge(imc-xref_+1:imc,j,l) = parent_mass_x_edge(ie_,jp,l) / real(xref_*yref_*zref_)
      end do
    end do

  end subroutine set_xedges_from_parent


  ! put_yedges (advecty.F90): halo rows and y-wall cells from the parent air
  ! mass at the moment of the put (parent_mass_y_edge, parent indices)

  subroutine set_yedges_from_parent( region, parent_mass_y_edge, mass_edge )

    use dims, only : im, jm, lm, parent, xref, yref, zref
    use dims, only : ibeg, jbeg, jend

    integer, intent(in)     ::  region
    real, intent(in)        ::  parent_mass_y_edge(ibeg(region)-1:,jbeg(region)-1:,:)
    real, intent(inout)     ::  mass_edge(0:,0:,:)

    integer   ::  imc, jmc, lmc, xref_, yref_, zref_, ib, jb, je_
    integer   ::  i, j, l, ip

    imc = im(region) ; jmc = jm(region) ; lmc = lm(region)
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    ib = ibeg(region) ; jb = jbeg(region) ; je_ = jend(region)

    do l = 1, lmc
      do i = 1, imc
        ip = ib + (i-1)/xref_
        mass_edge(i,0    ,l) = parent_mass_y_edge(ip,jb-1 ,l) / real(xref_*zref_)
        mass_edge(i,jmc+1,l) = parent_mass_y_edge(ip,je_+1,l) / real(xref_*zref_)
        do j = 1, yref_
          mass_edge(i,j      ,l) = parent_mass_y_edge(ip,jb ,l) / real(xref_*yref_*zref_)
          mass_edge(i,jmc+1-j,l) = parent_mass_y_edge(ip,je_,l) / real(xref_*yref_*zref_)
        end do
      end do
    end do

  end subroutine set_yedges_from_parent


  ! x advection at forward position k (advectxzoom + dynamu, air mass only)

  subroutine advect_edges_x( region, k, mass_edge )

    use dims,        only : im, jm, lm, parent, xref, yref, zref, tref
    use dims,        only : ibeg, iend, jbeg
    use dims,        only : zoom2D, nsplitsteps, n_operators, splitorderzoom
    use global_data, only : wind_dat

    integer, intent(in)     ::  region
    integer, intent(in)     ::  k
    real, intent(inout)     ::  mass_edge(0:,0:,:)

    real, dimension(:,:,:), pointer  ::  amc, amp
    integer   ::  imc, jmc, lmc, xref_, yref_, zref_, tref_, ib, ie_, jb
    integer   ::  n, q, i, j, l, js, je, ls, le, is, ie, jp
    logical   ::  enc, zero_outer
    character ::  dir
    real      ::  fw, fe

    imc = im(region) ; jmc = jm(region) ; lmc = lm(region)
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    tref_ = tref(region)/tref(parent(region))
    ib = ibeg(region) ; ie_ = iend(region) ; jb = jbeg(region)
    amc => wind_dat(region)%am_t
    amp => wind_dat(parent(region))%am_t

    ! scope (advectxzoom):
    js = 1 ; je = jmc ; ls = 1 ; le = lmc
    q = k/((nsplitsteps/2)*tref_)
    enc = .false.
    do n = 1, n_operators
      dir = splitorderzoom(region, q*(nsplitsteps/2)*tref_+n)
      select case ( dir )
        case ( 'x' )
          enc = .true.
        case ( 'y' )
          if ( .not. enc ) then
            js = 1 ; je = jmc
          else
            js = yref_+1 ; je = jmc-yref_
          end if
        case ( 'z' )
          if ( .not. enc ) then
            ls = 1 ; le = lmc
          else
            ls = zref_+1 ; le = lmc-zref_
          end if
          if ( zoom2D ) then
            ls = 1 ; le = lmc
          end if
      end select
    end do

    ! slopes: no coarse fluxes through the outer x-walls after the first sub-step
    ! (advectxzoom, only with slopes and without secmom):
    zero_outer = .false.
#ifdef slopes
#ifndef secmom
    zero_outer = ( mod(k,(nsplitsteps/2)*tref_) >= n_operators )
#endif
#endif

    is = xref_ ; ie = imc-xref_+1
    do l = ls, le
      do j = js, je
        ! compress_xedges:
        if ( xref_ > 1 ) then
          mass_edge(xref_      ,j,l) = sum(mass_edge(1:xref_,j,l))
          mass_edge(xref_-1    ,j,l) = mass_edge(0,j,l)
          mass_edge(imc-xref_+1,j,l) = sum(mass_edge(imc-xref_+1:imc,j,l))
          mass_edge(imc-xref_+2,j,l) = mass_edge(imc+1,j,l)
        end if
        ! mass update (dynamu):
        ! The fluxes through the outer x-walls of the child (amc(xref_-1) and
        ! amc(imc-xref_+1)) are taken from the parent, as set by coarsen_ambm
        ! (advect_tools.F90): the adjoint may have zeroed them temporarily in
        ! amc for the current step (adj_advectxzoom), while earlier steps of
        ! the block did use them.
        jp = jb + (j-1)/yref_
        do i = is, ie
          fw = amc(i-1,j,l)
          fe = amc(i  ,j,l)
          if ( i-1 == xref_-1     ) fw = amp(ib-1,jp,l) / real(yref_*zref_)
          if ( i   == imc-xref_+1 ) fe = amp(ie_ ,jp,l) / real(yref_*zref_)
          if ( zero_outer ) then
            if ( i-1 == xref_-1     ) fw = 0.0
            if ( i   == imc-xref_+1 ) fe = 0.0
          end if
          mass_edge(i,j,l) = mass_edge(i,j,l) + fw - fe
        end do
        ! uncompress_xedges:
        if ( xref_ > 1 ) then
          mass_edge(1:xref_        ,j,l) = mass_edge(xref_      ,j,l) / real(xref_)
          mass_edge(imc-xref_+1:imc,j,l) = mass_edge(imc-xref_+1,j,l) / real(xref_)
        end if
      end do
    end do

    nullify( amc, amp )

  end subroutine advect_edges_x


  ! y advection at forward position k (advectyzoom + dynamv, air mass only)

  subroutine advect_edges_y( region, k, mass_edge )

    use dims,        only : im, jm, lm, parent, xref, yref, zref, tref
    use dims,        only : ibeg, jbeg, jend
    use dims,        only : zoom2D, nsplitsteps, n_operators, splitorderzoom
    use global_data, only : wind_dat

    integer, intent(in)     ::  region
    integer, intent(in)     ::  k
    real, intent(inout)     ::  mass_edge(0:,0:,:)

    real, dimension(:,:,:), pointer  ::  bmc, bmp
    integer   ::  imc, jmc, lmc, xref_, yref_, zref_, tref_, ib, jb, je_
    integer   ::  n, q, i, j, l, is, ie, ls, le, js, je, ip
    logical   ::  enc, zero_outer
    character ::  dir
    real      ::  fs, fn

    imc = im(region) ; jmc = jm(region) ; lmc = lm(region)
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    tref_ = tref(region)/tref(parent(region))
    ib = ibeg(region) ; jb = jbeg(region) ; je_ = jend(region)
    bmc => wind_dat(region)%bm_t
    bmp => wind_dat(parent(region))%bm_t

    ! scope (advectyzoom):
    is = 1 ; ie = imc ; ls = 1 ; le = lmc
    q = k/((nsplitsteps/2)*tref_)
    enc = .false.
    do n = 1, n_operators
      dir = splitorderzoom(region, q*(nsplitsteps/2)*tref_+n)
      select case ( dir )
        case ( 'x' )
          if ( .not. enc ) then
            is = 1 ; ie = imc
          else
            is = xref_+1 ; ie = imc-xref_
          end if
        case ( 'y' )
          enc = .true.
        case ( 'z' )
          if ( .not. enc ) then
            ls = 1 ; le = lmc
          else
            ls = zref_+1 ; le = lmc-zref_
          end if
          if ( zoom2D ) then
            ls = 1 ; le = lmc
          end if
      end select
    end do

    ! same for the outer y-walls (advectyzoom, only with slopes and without secmom):
    zero_outer = .false.
#ifdef slopes
#ifndef secmom
    zero_outer = ( mod(k,(nsplitsteps/2)*tref_) >= n_operators )
#endif
#endif

    js = yref_ ; je = jmc-yref_+1
    do l = ls, le
      do i = is, ie
        ! compress_yedges:
        if ( yref_ > 1 ) then
          mass_edge(i,yref_      ,l) = sum(mass_edge(i,1:yref_,l))
          mass_edge(i,yref_-1    ,l) = mass_edge(i,0,l)
          mass_edge(i,jmc-yref_+1,l) = sum(mass_edge(i,jmc-yref_+1:jmc,l))
          mass_edge(i,jmc-yref_+2,l) = mass_edge(i,jmc+1,l)
        end if
        ! mass update (dynamv):
        ! fluxes through the outer y-walls (bmc(:,yref_), bmc(:,jmc-yref_+2))
        ! from the parent, as set by coarsen_ambm (see advect_edges_x):
        ip = ib + (i-1)/xref_
        do j = js, je
          fs = bmc(i,j  ,l)
          fn = bmc(i,j+1,l)
          if ( j   == yref_       ) fs = bmp(ip,jb   ,l) / real(xref_*zref_)
          if ( j+1 == jmc-yref_+2 ) fn = bmp(ip,je_+1,l) / real(xref_*zref_)
          if ( zero_outer ) then
            if ( j   == yref_       ) fs = 0.0
            if ( j+1 == jmc-yref_+2 ) fn = 0.0
          end if
          mass_edge(i,j,l) = mass_edge(i,j,l) + fs - fn
        end do
        ! uncompress_yedges:
        if ( yref_ > 1 ) then
          mass_edge(i,1:yref_        ,l) = mass_edge(i,yref_      ,l) / real(yref_)
          mass_edge(i,jmc-yref_+1:jmc,l) = mass_edge(i,jmc-yref_+1,l) / real(yref_)
        end if
      end do
    end do

    nullify( bmc, bmp )

  end subroutine advect_edges_y


  ! z advection at forward position k (advectzzoom + dynamw, air mass only)

  subroutine advect_edges_z( region, k, mass_edge )

    use dims,        only : im, jm, lm, parent, xref, yref, tref
    use dims,        only : nsplitsteps, n_operators, splitorderzoom
    use global_data, only : wind_dat

    integer, intent(in)     ::  region
    integer, intent(in)     ::  k
    real, intent(inout)     ::  mass_edge(0:,0:,:)

    real, dimension(:,:,:), pointer  ::  cmc
    integer   ::  imc, jmc, lmc, xref_, yref_, tref_
    integer   ::  n, q, i, j, l, is, ie, js, je
    logical   ::  enc
    character ::  dir

    imc = im(region) ; jmc = jm(region) ; lmc = lm(region)
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    tref_ = tref(region)/tref(parent(region))
    cmc => wind_dat(region)%cm_t

    ! scope (advectzzoom):
    is = 1 ; ie = imc ; js = 1 ; je = jmc
    q = k/((nsplitsteps/2)*tref_)
    enc = .false.
    do n = 1, n_operators
      dir = splitorderzoom(region, q*(nsplitsteps/2)*tref_+n)
      select case ( dir )
        case ( 'x' )
          if ( .not. enc ) then
            is = 1 ; ie = imc
          else
            is = xref_+1 ; ie = imc-xref_
          end if
        case ( 'y' )
          if ( .not. enc ) then
            js = 1 ; je = jmc
          else
            js = yref_+1 ; je = jmc-yref_
          end if
        case ( 'z' )
          enc = .true.
      end select
    end do

    do l = 1, lmc
      do j = js, je
        do i = is, ie
          mass_edge(i,j,l) = mass_edge(i,j,l) + cmc(i,j,l-1) - cmc(i,j,l)
        end do
      end do
    end do

    nullify( cmc )

  end subroutine advect_edges_z


  ! mix_edges (zoom_tools.F90), air mass only

  subroutine mix_edges_mass( region, mass_edge )

    use dims,        only : im, jm, lm, parent, xref, yref
    use global_data, only : region_dat

    integer, intent(in)     ::  region
    real, intent(inout)     ::  mass_edge(0:,0:,:)

    integer, dimension(:,:), pointer  ::  edge
    real, allocatable                 ::  tp(:,:)
    integer   ::  imc, jmc, lmc, xref_, yref_, imp, jmp
    integer   ::  i, j, l, ip, jp
    real      ::  xyref

    imc = im(region) ; jmc = jm(region) ; lmc = lm(region)
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    edge => region_dat(region)%edge

    xyref = real(xref_*yref_)
    imp = imc/xref_ ; jmp = jmc/yref_
    allocate( tp(imp,jmp) )

    do l = 1, lmc
      tp = 0.0
      do j = 1, jmc
        jp = 1 + (j-1)/yref_
        do i = 1, imc
          if ( edge(i,j) /= -1 ) cycle
          ip = 1 + (i-1)/xref_
          tp(ip,jp) = tp(ip,jp) + mass_edge(i,j,l)/xyref
        end do
      end do
      if ( xref_ > 1 ) then          ! x-walls (skipped by mix_edges when xref_ == 1)
        do j = 1, jmc
          jp = 1 + (j-1)/yref_
          mass_edge(1:xref_        ,j,l) = tp(1  ,jp)
          mass_edge(imc-xref_+1:imc,j,l) = tp(imp,jp)
        end do
      end if
      do j = 1, yref_                ! south y-wall
        do i = 1, imc
          mass_edge(i,j,l) = tp(1+(i-1)/xref_,1)
        end do
      end do
      do j = jmc-yref_+1, jmc        ! north y-wall
        do i = 1, imc
          mass_edge(i,j,l) = tp(1+(i-1)/xref_,jmp)
        end do
      end do
    end do

    deallocate( tp )
    nullify( edge )

  end subroutine mix_edges_mass


!===========================================================================================================
!===========================================================================================================


  subroutine adj_update_parent(region)
    !
    ! update rm of its parent in correspondence with itself
    ! written by mike botchev, march-june 1999
    !
    use dims
    use global_data  ,only : mass_dat
    use MeteoData   , only : m_dat
    use chem_param,   only : ntracet
    use ParTools  ,   only : ntracetloc
    use toolbox,      only : escape_tm

    implicit none

    !__IO____________________________________________________________________

    integer,intent(in)   :: region

    !__LOCAL_VARIABLES_______________________________________________________


    real,dimension(:,:,:,:),pointer          :: rm,rxm,rym,rzm
    real,dimension(:,:,:,:),pointer          :: rmp,rxmp,rymp,rzmp
    real,dimension(:,:,:),pointer            :: m,mp
    real,dimension(:,:,:),allocatable        :: to_parent

    integer :: my_parent
    integer :: xref_, yref_, zref_
    integer :: i, j, l
    integer :: ip, jp, lp
    integer :: n
    integer :: imp, jmp, lmp
    integer :: imr, jmr, lmr
    integer :: nt


    !__START_SUBROUTINE______________________________________________________


    if (region == 1) return

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    my_parent = parent(region)

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    zref_ = zref(region)/zref(my_parent)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    mp => m_dat(my_parent)%data
    rmp => mass_dat(my_parent)%rm_t
    rxmp => mass_dat(my_parent)%rxm_t
    rymp => mass_dat(my_parent)%rym_t
    rzmp => mass_dat(my_parent)%rzm_t

    nt=ntracetloc


    if ( okdebug ) then
       print *,'adj_update_parent: my_parent=',my_parent, &
            ' x-,y-,zref_: ',xref_,yref_,zref_
    end if

    imp = im(region)/xref_
    jmp = jm(region)/yref_
    lmp = lmr     ! CMK changed: when paralel over levels lmp can differ from lmr

    ! check imp
    if ( ibeg(region) < iend(region) .and. imp /= iend(region)-ibeg(region)+1 )  &
         call escape_tm('update_parent: stop')

    ! check jmp
    if ( jmp /= jend(region)-jbeg(region)+1 ) &
         call escape_tm('update_parent: stop')


    allocate(to_parent(imp,jmp,lmp))

    ! reconstruct masses
    call reconstruct_required_var_2(region)


    do n=1,nt

      do l=1,lmp
        do j=1,jmp
          do i=1,imp

            to_parent(i,j,l)  =  rmp(ibeg(region)-1+i,jbeg(region)-1+j,l,n)
            rmp(ibeg(region)-1+i,jbeg(region)-1+j,l,n) = 0.0
          enddo
        enddo
      enddo

      do l=1,lmr
        do j=1,jmr
          jp = 1 + (j-1)/yref_
          do i=1,imr
            ip = 1 + (i-1)/xref_
            rm(i,j,l,n) = rm(i,j,l,n) + to_parent(ip,jp,l)
          enddo
        enddo
      enddo  !l


    end do ! n

    ! from Maarten,s adjoint version ?????
    ! call mix_edges(region)

    deallocate(to_parent)

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

    nullify(mp)
    nullify(rmp)
    nullify(rxmp)
    nullify(rymp)
    nullify(rzmp)


  end subroutine adj_update_parent

!===========================================================================================================
!===========================================================================================================


  subroutine reconstruct_required_var_2(region)

    ! from Maarten's adjoint version:
    ! provide the correct air masses to the child...
    ! e.g
    ! region 1      |x|yz         |vsc
    ! region 2           xyzvsccsv
    !          in backward processing, the adj_update should provide the correct masses in
    !          the interface cells of the chlid (thus: at the end of the child operations).
    !          the 'internal' cell masses are conforming the parent masses, but
    !          the calls on the edges (0, im+1, jm+1) are sometines 'partly' advected.
    ! The situation is different for the x-interface and the y-interface (no z-zooming).
    ! In the above example, the x cells (0, imr+1) have 'seen' in the forward:
    !           o put_xedges puts the mass of the parent in the x-interface
    !           o edges remain the same during the full advection
    ! The y-cells (0, jmr+1), however.....
    !           o put_ydeges puts FULLY X-advected masses in 0, jmr+1
    ! In summary:
    !           the adj_update_parent now has to put the following masses in the cells (case parent = XYZ...)
    !           (i = 0, imr+1 (j=1,jm)) : no advection
    !           (j = 0, jmr+1 (i=1,im)) : X
    !
    ! To complicate the things even further, the situation is different when the order is different:
    ! region 1     csvz|y|x      |
    ! region 2             csvzyx
    !           the adj_update_parent now has to put the following masses in the cells (case parent = ZYX...)
    !           (i = 0, imr+1 (j=1,jm)) : ZY
    !           (j = 0, jmr+1 (i=1,im)) : Z
    !
    ! Furthermore, the total masses of the child should match the mass of the underlying parent.....
    !

    use dims
    use chem_param,   only : ntracet
    use ParTools  ,   only : ntracetloc
    use global_data,  only : mass_dat, wind_dat
    use MeteoData   , only : m_dat
    use toolbox,      only : escape_tm

    implicit none

    !__IO____________________________________________________________________

    integer,intent(in)            ::  region

    !__LOCAL_VARIABLES_______________________________________________________


    real,dimension(:,:,:,:),pointer          :: rm,rmp
    real,dimension(:,:,:),pointer            :: m,mp,amp,bmp,cmp

    real ::dummy
    integer i,iend_loop,iip,ip,j,jp,l,lp,my_parent,n,xref_,yref_,zref_,lb,le,ib,ie,jb,je
    integer :: imp,jmp,lmp,imr,jmr,lmr,q
    real m_value,mxval,sum_mass
    integer,dimension(3) :: mxloc
    character(len=n_operators) :: stencil
    logical                    :: xyz


    !__START_SUBROUTINE______________________________________________________

    if (region==1) return

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    my_parent = parent(region)

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    zref_ = zref(region)/zref(my_parent)


    if (okdebug) then
       print *,' adj_updating parent=',my_parent,' x-,y-,zref_: ',xref_, yref_, zref_
    endif

    imp = im(region)/xref_
    jmp = jm(region)/yref_
    lmp = lm(region)/zref_

    if(ibeg(region) < iend(region) .and. imp /= iend(region)-ibeg(region)+1)  &
         call escape_tm('stopped in update_parent')
    if(jmp.ne.jend(region)-jbeg(region)+1) call escape_tm('stopped in update_parent')
    if(lmp.ne.lend(region)-lbeg(region)+1) call escape_tm('stopped in update_parent')

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t

    amp => wind_dat(my_parent)%am_t
    bmp => wind_dat(my_parent)%bm_t
    cmp => wind_dat(my_parent)%cm_t
    mp => m_dat(my_parent)%data
    rmp => mass_dat(my_parent)%rm_t

    ! put mass parent in the child cells, 0, imr+1

    q=(status(my_parent)-1)/((nsplitsteps/2))      ! find q - the place in the splitorderzoom of the parent
    do i=1,n_operators
       stencil(i:i)=splitorderzoom(my_parent,q*(nsplitsteps/2)+i)
    enddo

    ! find out if parent has order xyz, or zyx...
    if( scan(stencil,'x') >  scan(stencil,'z')) then
       xyz = .false.
    else
       xyz = .true.
    endif

    if(xcyc(region)==0 ) then
       do l=1,lmr
          do j=1,jmr
             jp = jbeg(region) + (j-1)/yref_
             if(xyz) then ! x-edges: -  (YZX backwards)
                m(0,j,l) = (mp(ibeg(region)-1,jp,l)   &
                     - bmp(ibeg(region)-1,jp,l  ) + bmp(ibeg(region)-1,jp+1,l) &
                     - cmp(ibeg(region)-1,jp,l-1) + cmp(ibeg(region)-1,jp  ,l) &
                     - amp(ibeg(region)-2,jp,l  ) + amp(ibeg(region)-1,jp  ,l)  )/yref_
                m(imr+1,j,l) = (mp(iend(region)+1,jp,l)  &
                     - bmp(iend(region)+1,jp,l  ) + bmp(iend(region)+1,jp+1,l) &
                     - cmp(iend(region)+1,jp,l-1) + cmp(iend(region)+1,jp  ,l) &
                     - amp(iend(region)  ,jp,l  ) + amp(iend(region)+1,jp  ,l)  )/yref_

             else ! x-edges YZ  (X backwards)
                m(0,j,l) = (mp(ibeg(region)-1,jp,l)       + amp(ibeg(region)-1,jp,l) - amp(ibeg(region)-2,jp,l)  )/yref_
                m(imr+1,j,l) = (mp(iend(region)+1,jp,l)   + amp(iend(region)+1,jp,l) - amp(iend(region)  ,jp,l)  )/yref_
             endif
          enddo
       enddo
    endif !if(xcyc(region)==0)

    if(xcyc(region)==0) then
       do l=1,lmr
          do i=1,imr
             ip = ibeg(region) + (i-1)/xref_

             if(xyz) then    ! y-edges: X  (ZY backwards)

                if(touch_sp(region)==0) then
                   m(i,0,l) = (mp(ip, jbeg(region)-1,l) &
                        - cmp(ip, jbeg(region)-1,l-1) + cmp(ip, jbeg(region)-1,l) &
                        + bmp(ip, jbeg(region)  ,l  ) - bmp(ip, jbeg(region)-1,l)  )/xref_
                endif
                if(touch_np(region)==0) then
                   m(i,jmr+1,l) = (mp(ip, jend(region)+1,l)    &
                        - cmp(ip, jend(region)+1,l-1) + cmp(ip, jend(region)+1,l) &
                        - bmp(ip, jend(region)+1,l)   + bmp(ip, jend(region)+2,l)  )/xref_
                endif

             else  ! y-edges Z  (XY backwards)

                if(touch_sp(region)==0) then
                   m(i,0,l) =  (mp(ip, jbeg(region)-1,l) &
                        - amp(ip-1, jbeg(region)-1,l) + amp(ip, jbeg(region)-1,l) &
                        + bmp(ip,   jbeg(region),  l) - bmp(ip, jbeg(region)-1,l)  )/xref_
                endif
                if(touch_np(region)==0) then
                   m(i,jmr+1,l) = (mp(ip, jend(region)+1,l)  &
                        - amp(ip-1, jend(region)+1,l) + amp(ip, jend(region)+1,l) &
                        - bmp(ip,   jend(region)+1,l) + bmp(ip, jend(region)+2,l)  )/xref_
                endif
             endif
          enddo
       enddo
    endif !if(xcyc(region)==0)


  end subroutine reconstruct_required_var_2

!===========================================================================================================
!===========================================================================================================

end module adj_zoom_tools

!===========================================================================================================
!===========================================================================================================
