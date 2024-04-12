!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module advectx

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public :: AdvectX_Init, AdvectX_Done
  public :: advectxzoom, put_xedges


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'AdvectX'


  ! --- local ------------------------------------

  integer    ::  itim_zoom
  integer    ::  itim_work
  integer    ::  itim_dynam
  integer    ::  itim_loop
  integer    ::  itim_compress
  integer    ::  itim_uncompress
  integer    ::  itim_put_xedges


contains


  ! ====================================================================


  subroutine AdvectX_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/AdvectX_Init'

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_zoom, 'advectx zoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_work, 'advectx work', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynam, 'advectx dynam', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_loop, 'advectx loop', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_compress, 'advectx compress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_uncompress, 'advectx uncompress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_put_xedges, 'advectx put xedges', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine AdvectX_Init


  ! ***


  subroutine AdvectX_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/AdvectX_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine AdvectX_Done


  ! ***


  ! set parameters for advectx
  ! written by patrick berkvens and mike botchev, march-june 1999
  ! updated and modified by MK, dec 2002

  subroutine advectxzoom( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,          only : xref, yref, zref, tref, im, jm, lm
    use dims,          only : rstatus => status, nsplitsteps, n_operators, zoom2D
    use dims,          only : touch_sp, touch_np, okdebug
    use dims,          only : parent, splitorderzoom, zero
#ifdef with_budgets
    use budget_global, only : sum_advection
#endif
    use global_data  , only : mass_dat, wind_dat
    use toolbox,       only : escape_tm
#ifdef MPI
    use mpi_const,     only : my_real, com_lev, ierr
    use mpi_const,     only : root_t, mpi_sum, pe_first_tracer
    use mpi_comm,      only : barrier_t
#endif
    use partools     , only : myid, lmloc, ntracetloc
    use chem_param   , only : ntracet

    ! input/output
    integer,intent(in)                  :: region
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/advectxzoom'

    ! local
    real,dimension(:,:,:),  pointer     :: am

    integer                             :: js,je,ls,le,n,q
    integer                             :: imr,jmr,lmr,tref_,xref_,yref_,zref_
    real,dimension(:,:),allocatable     :: am0,am1  ! for slopes only
    logical                             :: x_encountered
    character(len=1)                    :: dir
#ifdef with_budgets
    real                                :: sum_old(ntracet), sum_new(ntracet)
    real                                :: sum_old_all(ntracet), sum_new_all(ntracet)
    integer                             :: itr
#endif

    ! start

    if ( ntracetloc == 0 ) return !WP! only PE's with nonzero lmloc proceed

    call GO_Timer_Start( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    ! write BC to cildren
    call put_xedges( region, status )
    IF_NOTOK_RETURN(status=1)

    tref_ = tref(region)/tref(parent(region))
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    allocate(am0(0:jm(region)+1, 0:lmr+1))
    allocate(am1(0:jm(region)+1, 0:lmr+1))

    am => wind_dat(region)%am_t

#ifdef with_budgets
    do itr = 1, ntracet
        sum_old(itr) = sum(mass_dat(region)%rm_t(1:im(region),1:jm(region),1:lmr,itr))
    end do
#endif

    ! determine the scope for advectx:

    if ( region == 1 ) then
       yref_ = 0
       zref_ = 0  ! to have js/je and ls/le properly computed
    end if
    q=rstatus(region)/((nsplitsteps/2)*tref_)
    ! find q - the place in the splitorderzoom
    ! corresponding to the begining of the last
    ! triple x-y-z of the parent
    x_encountered=.false.

    do n=1,n_operators
       ! now track four following steps from q !wp! changed to n_op.. steps

       dir=splitorderzoom(region,q*(nsplitsteps/2)*tref_+n)
       select case(dir)
       case('x')
          x_encountered=.true.
       case('y')
          if ( .not. x_encountered) then
             js=1                             ! y-substep is before x =>
             je=jm(region)                    ! full j-scope
          else
             if( touch_sp(region) == 1 ) then
                js=1                             !CMK special case...touching SP
             else
                js=yref_+1                       ! x-substep is before y =>
             end if
             if( touch_np(region) == 1 ) then
                je=jm(region)                    ! CMK special case: touching NP
             else
                je=jm(region)-yref_              ! restricted j-scope
             end if
          end if
       case('z')
          if  ( .not. x_encountered ) then
             ls=1                             ! z-substep is before x =>
             le=lmr                    ! full l-scope
          else
             ls=zref_+1                       ! x-substep is before z =>
             le=lmr-zref_              ! restricted l-scope
          end if
          if (zoom2D) then
             ls=1; le=lmr    !WP! this is always the case
          end if
       case ('c')
       case ('v')
       case ('d')
       case ('s')
       case default
          print *,'advectxzoom: strange value in splitorderzoom(',region,',',   &
               q*(nsplitsteps/2)*tref_+n,'): ',q
          call escape_tm( 'advectxzoom: error in x-advection ')
       end select

    end do  ! n=1,n_operators

#ifdef slopes
#ifndef secmom
    if ( ( mod(rstatus(region),(nsplitsteps/2)*tref_) >= n_operators ) .and. &
         ( im(region)/xref(region) < im(1) ) ) then
       ! IF (1) more than fourth substep, (2) slope scheme, and
       ! (3) the region is not [0;360] degrees wide THEN
       ! at this substep no coarse fluxes will be applied

       ! for slopes only: zero out velocities via the x-edges of the region
       ! to assure that no fluxes will be applied

       if ( okdebug ) print *,'advectxzoom: zeroing out outer mass fluxes'
       am0 = am(    xref_-1,:,:)
       am(    xref_-1,:,:) = zero
       am1 = am(imr-xref_+1,:,:)
       am(imr-xref_+1,:,:) = zero

    end if
#endif
#endif

    call advectx_work( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)

#ifdef slopes
#ifndef secmom
    if ( ( mod(rstatus(region),(nsplitsteps/2)*tref_) >= n_operators ).and. &
         ( im(region)/xref(region) < im(1) ) ) then

       ! for slopes only: recreate velocities via the x-edges
       am(    xref_-1,:,:) = am0
       am(imr-xref_+1,:,:) = am1

    end if
#endif
#endif

#ifdef with_budgets
    do itr = 1, ntracet
        sum_new(itr) = sum(mass_dat(region)%rm_t(1:im(region),1:jm(region),1:lmr,itr))
    end do
#ifdef MPI
    if( myid == pe_first_tracer ) &
         sum_advection(region,:) = sum_advection(region,:) + sum_new - sum_old
    if( okdebug .and. myid == root_t ) &
         print *, 'advectxzoom: region, sum_advection ',&
         region,sum_advection(region,:),sum_new,sum_old
#else
    sum_advection(region,:) = sum_advection(region,:) + sum_new - sum_old
    if ( okdebug ) print *, 'advectxzoom: region, sum_advection ', &
         region,sum_advection(region,:),sum_new,sum_old
#endif
#endif

    deallocate(am0)
    deallocate(am1)
    nullify(am)

    call GO_Timer_End( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectxzoom


  !
  ! makes reduced grid pre-/postprocessing and switches between
  ! dynamu and dynamu1
  ! written by mike botchev, march-june 1999
  !

  subroutine advectx_work( region, js,je,ls,le, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : im, jm, lm
    use redgridZoom, only : grid_reduced, nred, uni2red, uni2red_mf, red2uni
    use global_data, only : wind_dat
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use partools   , only : lmloc
    use chem_param,  only : ntracet

    ! input
    integer,intent(in) :: region
    integer,intent(in) :: js
    integer,intent(in) :: je
    integer,intent(in) :: ls
    integer,intent(in) :: le
    integer,intent(out) :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/advecty_work'

    ! local
    real,dimension(:,:,:),pointer     :: m
    real,dimension(:,:,:),pointer     :: am
    real,dimension(:,:,:),allocatable :: m_uni   ! for reduced grid...
    real,dimension(:,:,:),allocatable :: am_uni  ! for reduced grid...
    real                 :: mxval
    integer,dimension(3) :: mxloc
    integer              :: imr,jmr,lmr

    ! start

#ifndef slopes
    write (gol,'("slopes advection not defined")'); call goErr
    TRACEBACK; status=1; return
#endif

    call GO_Timer_Start( itim_work, status )
    IF_NOTOK_RETURN(status=1)

    imr = im(region) ; jmr = jm(region) ; lmr = lm(region)

    allocate(m_uni(-1:im(region)+2,-1:jm(region)+2, lmr))
    allocate(am_uni(0:im(region)+1, 0:jm(region)+1, 0:lmr+1))

    am => wind_dat(region)%am_t
    m  => m_dat(region)%data

    ! transform to the reduced grid:
    ! check for reduced grid in region

    if ( grid_reduced .and. (nred(region) /= 0) ) then
       ! save non-reduced m and am in m_uni and am_uni:
       m_uni = m
       am_uni = am
       ! reduce m,rm,rxm,rym,rzm:
       call uni2red(region)
       ! reduce am:
       call uni2red_mf(region)
    end if

    call dynamu( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    ! transform from the reduced grid:

    if ( grid_reduced .and. nred(region) /= 0 ) then
       ! advection on uniform grid:
       m_uni(1:imr,1:jmr,1:lmr) = m_uni(1:imr,1:jmr,1:lmr) + &
            am_uni(0:imr-1,1:jmr,1:lmr) - am_uni(1:imr,1:jmr,1:lmr)

       ! redistribute rm,rxm,rym,rzm by using m_uni and m
       call red2uni(region,m_uni)
       ! recreate rm/rxm/rym/rzm by using equal masses (not m_uni)
       !call red2uni_em(region)
       ! recreate am:
       am = am_uni
    end if

    nullify(am)
    nullify(m)

    deallocate(am_uni)
    deallocate(m_uni)

    call GO_Timer_End( itim_work, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectx_work


  !-----------------------------------------------------------------------
  !
  !****   dynamu  - east-west tracer transport  - v 9.1
  !
  ! programmed by       mh  mpi HH      23-feb-1995
  !
  ! purpose
  ! -------
  ! calculate amount of tracer moved in an east-west advection
  ! substep
  !
  ! interface
  ! ---------
  ! call dynamu
  !
  ! method
  ! ------
  ! slopes scheme
  !
  ! externals
  ! ---------
  ! none
  !
  ! reference
  ! ---------
  ! Russel and Lerner, 1979
  !-----------------------------------------------------------------------
  ! fixed bug in calculating maximum courant number
  !                    mh, Thu, Feb 24, 1994 12:49:27
  !     included code for limits of slopes to prevent negative tracer
  !     masses                   mh, 20-jun-1994
  !     zoom version written by mike botchev, march-june 1999
  !-----------------------------------------------------------------------

  subroutine dynamu( region, js,je,ls,le, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : okdebug, nregions, parent
    use dims,        only : im, jm, lm
    use dims,        only : xref, yref, zref, tref
    use dims,        only : zero, one
    use dims,        only : nloop_max, limits, xcyc
    use redgridZoom, only : grid_reduced, nred, imred
    use global_data, only : wind_dat
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use AdvectM_CFL, only : advectx_get_nloop
    use zoom_tools,  only : mix_edges
    use toolbox,     only : escape_tm
#ifdef with_budgets
    use budget_global, only : budget_flux, iflux1, iflux2
    use budget_global, only : apply_budget_global_flux => apply_budget_global
    use partools     , only : offsetn
    use chem_param   , only : ra
#endif
#ifdef MPI
    use mpi_const,   only : com_trac
    use mpi_const,   only : my_real,mpi_max,mpi_min
    use mpi_comm,    only : barrier_t
#endif
    use partools   , only : myid, root_t
    use partools   , only : ntracetloc

    ! input/output
    integer,intent(in) :: region
    integer,intent(in) :: js
    integer,intent(in) :: je
    integer,intent(in) :: ls
    integer,intent(in) :: le
    integer,intent(out) :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/dynamu'

    ! local
    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer    :: m,am
    real,dimension(:),allocatable              :: f,pf,fy,fz, mnew
    integer   ::  i,ie,is,j,l,n, nglob
    integer   ::  imr,jmr,lmr,tref_,xref_,yref_,zref_,ic
    real      ::  min_one,max_one,min_all,max_all
    real      ::  alpha
    real      ::  x, rmold
    ! cfl loop:
    integer            ::  iloop, nloop
    ! Openmp parameters
    integer            ::  iie, ns, ne
    integer            ::  nfail
    logical            ::  special_grid

    ! start

    call GO_Timer_Start( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_limits
    if ( .not. limits ) stop 'macro "with_limits" defined but variable "limits" not'
#else
    if ( limits ) stop 'macro "with_limits" not defined but variable "limits" is'
#endif

    if ( okdebug ) print *,'dynamu: region=',region, &
         ' js,je,ls,le=',js,je,ls,le
    if ( (region < 0) .or. (region > nregions) ) &
         call escape_tm( 'dynamu: STOP, illegal number of region !!!')

    ! prepare x-edges to allow for uniform work:

    call compress_xedges( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    ! compute refinement factors with respect to the parent
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    tref_ = tref(region)/tref(parent(region))
    imr=im(region);jmr=jm(region);lmr=lm(region)

    am => wind_dat(region)%am_t
    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    ! check js,je,ls,le:
    if((js/=yref_+1)  .and.(js/=1))   stop 'dynamu: Wrong value for JS '
    if((je/=jmr-yref_).and.(je/=jmr)) stop 'dynamu: Wrong value for JE '
    if((ls/=zref_+1)  .and.(ls/=1))   stop 'dynamu: Wrong value for LS '
    if((le/=lmr-zref_).and.(le/=lmr)) stop 'dynamu: Wrong value for LE '

    ! compute is/ie -- cells is:ie,js:je:ls:le will be updated
    is = xref_; ie = imr-xref_+1
    if ( xcyc(region) == 1 ) then
       ! periodic boundary condition
       is = 1; ie = imr
    end if

    ! reduced grid ?
    special_grid = grid_reduced .and. nred(region) /= 0

    if ( okdebug ) then

       max_one=maxval(m(1:imr,1:jmr,1:lmr))
       min_one=minval(m(1:imr,1:jmr,1:lmr))
#ifdef MPI
       call barrier_t
       call mpi_reduce(max_one,max_all,1,my_real,mpi_max,root_t,com_trac,ierr)
       call mpi_reduce(min_one,min_all,1,my_real,mpi_min,root_t,com_trac,ierr)
       if ( myid == root_t ) print *,'dynamu: m_min,m_max: ',min_all,max_all
       max_one=maxval(rm(1:imr,1:jmr,1:lmr,:))
       min_one=minval(rm(1:imr,1:jmr,1:lmr,:))
       call mpi_reduce(max_one,max_all,1,my_real,mpi_max,root_t,com_trac,ierr)
       call mpi_reduce(min_one,min_all,1,my_real,mpi_min,root_t,com_trac,ierr)
       if ( myid == root_t ) print *,'dynamu: rm_min,rm_max: ',min_all,max_all
#else
       print *,'dynamu: m_min,m_max: ',min_one, max_one
       max_one=maxval(rm(1:imr,1:jmr,1:lmr,:))
       min_one=minval(rm(1:imr,1:jmr,1:lmr,:))
       print *,'dynamu: rm_min,rm_max: ',min_one,max_one
#endif

    end if

    ! loop over vertical layers and latitudes

    call GO_Timer_Start( itim_loop, status )
    IF_NOTOK_RETURN(status=1)

    ! no failues yet:
    nfail = 0

    !afds$dsakfaOMP  shared  ( xref, xcyc ) &
    !$OMP PARALLEL &
    !$OMP  default (none) &
#ifdef with_budgets
    !$OMP  shared  (apply_budget_global_flux, iflux1, iflux2, budget_flux,ra) &
    !$OMP  shared  (offsetn) &
    !$OMP  private (nglob) &
#endif
    !$OMP  shared  ( region ) &
    !$OMP  shared  ( imr ) &
    !$OMP  shared  ( is, ie, js, je, ls, le ) &
    !$OMP  shared  ( special_grid, nred, imred ) &
    !$OMP  shared  ( m, am, rm, rxm, rym, rzm ) &
    !$OMP  reduction ( + : nfail ) &
    !$OMP  private ( status ) &
    !$OMP  shared  ( ntracetloc ) &
    !$OMP  private ( iloop, nloop ) &
    !$OMP  private ( i, j, l, iie ) &
    !$OMP  private ( alpha, f, pf, fy, fz, mnew )

    ! cells:
    allocate( mnew(imr) )
    ! interfaces:
    allocate(  f(0:imr) )
    allocate( pf(0:imr) )
    allocate( fy(0:imr) )
    allocate( fz(0:imr) )

    ! loop over levels:
    !$OMP   DO
    do l = ls, le

      ! loop over rows:
      do j = js, je

        ! reduced grid: less points will be handled
        iie = ie
        if ( special_grid ) iie = imred(j,region)

        ! get number of loops,
        ! divide am by number of loops:
        call advectx_get_nloop( is,iie, xcyc(region)==1, &
                                 m(is-1:iie+1,j,l), &
                                 am(is-1:iie,j,l), &
                                 nloop, status )
        if (status/=0) then
          nfail = nfail + 1
          cycle
        end if

        ! CFL loop
        do iloop = 1, nloop

          ! calculate new air mass distribution
          mnew(is:iie)=m(is:iie,j,l) + am(is-1:iie-1,j,l)-am(is:iie,j,l)

          if ( xcyc(region) == 1 ) then
            m  (0    ,j,l)   = m  (iie,j,l) !periodic BCs...
            m  (iie+1,j,l)    = m  (1,  j,l)
          end if

          ! loop over tracers
          do n=1,ntracetloc

            ! if requested limit zonal slopes such that no negative
            ! tracer masses should occur
#ifdef with_limits
            do i = is - 1, iie + 1
              rxm(i,j,l,n) = max( min( rxm(i,j,l,n), rm(i,j,l,n)), -rm(i,j,l,n) )
            end do
#endif

            ! periodic boundary conditions
            if ( xcyc(region) == 1 ) then
              rm (0,    j,l,n) = rm (iie,j,l,n)
              rxm(0,    j,l,n) = rxm(iie,j,l,n)
              rym(0,    j,l,n) = rym(iie,j,l,n)
              rzm(0,    j,l,n) = rzm(iie,j,l,n)
              rm (iie+1,j,l,n) = rm (1,  j,l,n)
              rxm(iie+1,j,l,n) = rxm(1,  j,l,n)
              rym(iie+1,j,l,n) = rym(1,  j,l,n)
              rzm(iie+1,j,l,n) = rzm(1,  j,l,n)
            end if

            ! calculate fluxes for rm,rxm,rym,rzm
            do i=is-1,iie
              if (am(i,j,l)>=zero)  then
                alpha=am(i,j,l)/m(i,j,l)
                f(i)=alpha*(rm(i,j,l,n)+(one-alpha)*rxm(i,j,l,n))
                pf(i)=am(i,j,l)*(alpha*alpha*rxm(i,j,l,n) - 3.*f(i))
                fy(i)=alpha*rym(i,j,l,n)
                fz(i)=alpha*rzm(i,j,l,n)
              else
                alpha=am(i,j,l)/m(i+1,j,l)
                f(i)=alpha*(rm(i+1,j,l,n)-(one+alpha)*rxm(i+1,j,l,n))
                pf(i)=am(i,j,l)*(alpha*alpha*rxm(i+1,j,l,n) - 3.*f(i))
                fy(i)=alpha*rym(i+1,j,l,n)
                fz(i)=alpha*rzm(i+1,j,l,n)
              end if
            end do

            ! redefine left fluxes using periodicity:
            if ( xcyc(region) == 1 ) then
              f(0)  = f (iie)
              pf(0) = pf(iie)
              fy(0) = fy(iie)
              fz(0) = fz(iie)
            end if

#ifdef with_budgets
            !
            ! add up flux budget!
            ! note that this fails if grid is reduced (iie check)
            !
            if ( apply_budget_global_flux ) then
              nglob = offsetn + n
              if ( (iflux1(region)-1 >= is) .and. (iflux1(region)-1 <= iie) ) then
                 budget_flux(region)%flux_x1(j,l,nglob) =  budget_flux(region)%flux_x1(j,l,nglob) &
                    + f(iflux1(region)-1)*1e3/ra(nglob)   ! moles
              endif
              if ( (iflux2(region)-1 >= is) .and. (iflux2(region)-1 <= iie) ) then
                 budget_flux(region)%flux_x2(j,l,nglob) =  budget_flux(region)%flux_x2(j,l,nglob) &
                    + f(iflux2(region)-1)*1e3/ra(nglob)
              endif
            end if
#endif
            !
            ! calculate new tracer mass, and tracer mass slopes
            !
            do i=is,iie
              rm(i,j,l,n)=rm(i,j,l,n)+(f(i-1)-f(i))
              rxm(i,j,l,n)=rxm(i,j,l,n)+(pf(i-1)-pf(i)       &
                   -(am(i-1,j,l)-am(i,j,l))*rxm(i,j,l,n)     &
                   +3.*((am(i-1,j,l)+am(i,j,l))*rm(i,j,l,n)  &
                   -(f(i-1)+f(i))*m(i,j,l)))/mnew(i)
#ifdef with_limits
              !CMK: apply limits again: might be in nloop!
              rxm(i,j,l,n) = max(min(rxm(i,j,l,n),rm(i,j,l,n)),-rm(i,j,l,n))
#endif
              rym(i,j,l,n)=rym(i,j,l,n)+(fy(i-1)-fy(i))
              rzm(i,j,l,n)=rzm(i,j,l,n)+(fz(i-1)-fz(i))
            end do !forall

            ! periodic boundary conditions
            if ( xcyc(region) == 1 ) then
              rm (    0,j,l,n) = rm (iie,j,l,n)
              rxm(    0,j,l,n) = rxm(iie,j,l,n)
              rym(    0,j,l,n) = rym(iie,j,l,n)
              rzm(    0,j,l,n) = rzm(iie,j,l,n)
              rm (iie+1,j,l,n) = rm (  1,j,l,n)
              rxm(iie+1,j,l,n) = rxm(  1,j,l,n)
              rym(iie+1,j,l,n) = rym(  1,j,l,n)
              rzm(iie+1,j,l,n) = rzm(  1,j,l,n)
            end if

          end do  ! end of n-loop

          ! store new air mass in m array
          m(is:iie,j,l)=mnew(is:iie)

        end do  ! cfl loop

        ! restore 'old' am
        am(is-1:ie,j,l) = am(is-1:ie,j,l)*nloop

      end do   ! l-loop
    end do    ! j-loop
    !$OMP   END DO

    ! clear:
    deallocate(f)
    deallocate(pf)
    deallocate(fy)
    deallocate(fz)
    deallocate(mnew)

    !$OMP END PARALLEL

    ! check ...
    if ( nfail > 0 ) then
      write (gol,'("failures from x advection : ",i6)') nfail; call goErr
      TRACEBACK; status=1; return
    end if

    call GO_Timer_End( itim_loop, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    nullify(am)
    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

    ! *

    call uncompress_xedges( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    ! *

    call mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    ! *

    call GO_Timer_End( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine dynamu


  ! ***


  !
  ! this is for slope only: condense data at the x-edges of the zoom region
  ! to allow for uniform work in dynamu/dynamu1
  ! written by mike botchev, march-june 1999
  !

  subroutine compress_xedges( region, js,je,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : xref, parent, im, okdebug
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use partools   , only : ntracetloc
    use chem_param,  only : ntracet

    ! --- in/out ---------------------------------

    integer, intent(in)   ::  region
    integer, intent(in)   ::  js
    integer, intent(in)   ::  je
    integer, intent(in)   ::  ls
    integer, intent(in)   ::  le
    integer, intent(out)  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/compress_xedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer    :: m
    integer :: imr,j,l,n,xref_

    ! --- start ----------------------------------

    call GO_Timer_Start( itim_compress, status )
    IF_NOTOK_RETURN(status=1)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    xref_ = xref(region)/xref(parent(region))
    imr = im(region)
    if ( (xref_ == 1 ) .or. ( im(region)/xref(region) == im(1) ) ) then
       if ( okdebug ) print *,"compress_xedges:", &
            " no refinement or periodic bc's, nothing to do"
       call GO_Timer_End( itim_compress, status )
       IF_NOTOK_RETURN(status=1)
       return
    end if
    if ( okdebug ) then
       print *,'compress_xedges: region=',region,' imr=',imr
       print *,'compress_xedges:        cells to be updated: ', &
            xref_-1,xref_,imr-xref_+1,imr-xref_+2
    end if

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( imr ) &
    !$OMP   shared ( xref_ ) &
    !$OMP   shared ( js, je, ls, le ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m, rm, rxm, rym, rzm ) &
    !$OMP   private ( j, l, n )
    !$OMP   DO
    do j=js,je
       do l=ls,le
          !sum_m = sum(m(1:xref_,j,l),1);  m (xref_,  j,l) = sum_m
          m (xref_,  j,l) = sum(m(1:xref_,j,l))
          m (xref_-1,j,l) = m (0,j,l)
          m (imr-xref_+1,j,l) = sum(m(imr-xref_+1:imr,j,l))
          m (imr-xref_+2,j,l) = m (imr+1,j,l)

          do n=1,ntracetloc
             rm (xref_,  j,l,n) = sum(rm(1:xref_,j,l,n))
             rm (xref_-1,j,l,n) = rm (0,j,l,n)
             rxm(xref_-1,j,l,n) = rxm(0,j,l,n)
             rym(xref_-1,j,l,n) = rym(0,j,l,n)
             rzm(xref_-1,j,l,n) = rzm(0,j,l,n)

             rm (imr-xref_+1,j,l,n) = sum(rm(imr-xref_+1:imr,j,l,n))
             rm (imr-xref_+2,j,l,n) = rm (imr+1,j,l,n)
             rxm(imr-xref_+2,j,l,n) = rxm(imr+1,j,l,n)
             rym(imr-xref_+2,j,l,n) = rym(imr+1,j,l,n)
             rzm(imr-xref_+2,j,l,n) = rzm(imr+1,j,l,n)
          end do
       end do
    end do
    !$OMP   END DO
    !$OMP END PARALLEL

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

    call GO_Timer_End( itim_compress, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine compress_xedges


  ! ***


  !
  ! distribute data at the x-edges of the zoom region
  ! over the whole interface cell
  ! written by mike botchev, march-june 1999
  !

  subroutine uncompress_xedges( region, js,je,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : xref, parent, im, okdebug
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use partools   , only : ntracetloc
    use chem_param,  only : ntracet

    ! --- in/out ---------------------------------

    integer, intent(in)   ::  region
    integer, intent(in)   ::  js
    integer, intent(in)   ::  je
    integer, intent(in)   ::  ls
    integer, intent(in)   ::  le
    integer, intent(out)  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/uncompress_xedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer    :: m
    integer :: imr,j,l,n,xref_
    real    :: m_edge,rm_edge

    ! --- start ----------------------------------

    call GO_Timer_Start( itim_uncompress, status )
    IF_NOTOK_RETURN(status=1)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    xref_ = xref(region)/xref(parent(region))
    imr = im(region)

    if ( (xref_==1) .or. (im(region)/xref(region)==im(1)) ) then
       if ( okdebug ) print *,"uncompress_xedges:", &
            " no refinement or periodic bc's, nothing to do"
       call GO_Timer_End( itim_uncompress, status )
       IF_NOTOK_RETURN(status=1)
       return
    end if
    if ( okdebug ) then
       print *,'uncompress_xedges: region=',region,' imr=',imr
       print *,'uncompress_xedges: cells to be updated: ',1,':', &
            xref_,' ',imr-xref_+1,':',imr
    end if

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( imr ) &
    !$OMP   shared ( xref_ ) &
    !$OMP   shared ( js, je, ls, le ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m, rm, rxm, rym, rzm ) &
    !$OMP   private ( j, l, n ) &
    !$OMP   private ( m_edge, rm_edge )
    !$OMP   DO
    do j=js,je
       do l=ls,le
          m_edge = m (xref_,j,l)
          m (1:xref_,j,l) = m_edge/xref_

          m_edge = m (imr-xref_+1,j,l)
          m (imr-xref_+1:imr,j,l) = m_edge/xref_
          !
          do n=1,ntracetloc
             rm_edge = rm (xref_,j,l,n)
             rm (1:xref_,j,l,n) = rm_edge/xref_
             !
             !*      rxm(1:xref_-1,j,l,n) = rxm(xref_,j,l,n)
             !*      rym(1:xref_-1,j,l,n) = rym(xref_,j,l,n)
             !*      rzm(1:xref_-1,j,l,n) = rzm(xref_,j,l,n)
             rxm(1:xref_,j,l,n) = rxm(xref_,j,l,n)/xref_
             rym(1:xref_,j,l,n) = rym(xref_,j,l,n)/xref_
             rzm(1:xref_,j,l,n) = rzm(xref_,j,l,n)/xref_
             !
             rm_edge = rm (imr-xref_+1,j,l,n)
             rm (imr-xref_+1:imr,j,l,n) = rm_edge/xref_
             !
             !*      rxm(imr-xref_+2:imr,j,l,n) = rxm(imr-xref_+1,j,l,n)
             !*      rym(imr-xref_+2:imr,j,l,n) = rym(imr-xref_+1,j,l,n)
             !*      rzm(imr-xref_+2:imr,j,l,n) = rzm(imr-xref_+1,j,l,n)
             rxm(imr-xref_+1:imr,j,l,n) = rxm(imr-xref_+1,j,l,n)/xref_
             rym(imr-xref_+1:imr,j,l,n) = rym(imr-xref_+1,j,l,n)/xref_
             rzm(imr-xref_+1:imr,j,l,n) = rzm(imr-xref_+1,j,l,n)/xref_
          end do
       end do
    end do
    !$OMP   END DO
    !$OMP END PARALLEL

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

    call GO_Timer_End( itim_uncompress, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine uncompress_xedges


  !
  ! passes values along the x-boundaries to all the children
  ! written by mike botchev, march-june 1999
  ! structure implemented by MK, dec 2002
  !

  subroutine put_xedges( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,          only : im, jm, lm, xref, yref, zref
    use dims,          only : ibeg, jbeg, lbeg, iend, nregions
    use dims,          only : children, okdebug
    use global_data,   only : mass_dat
    use MeteoData    , only : m_dat
#ifdef with_budgets
    use budget_global, only : sum_advection
#endif
    use chem_param,    only : ra, ntracet
#ifdef MPI
    use mpi_const
#endif
    use partools     , only : myid
    use partools     , only : lmloc, offsetl
    use partools     , only : ntracetloc, offsetn
    use partools     , only : which_par, previous_par
    use ParTools     , only : Par_Check_Domain

    ! in/out
    integer,intent(in)   :: region
    integer, intent(out)             ::  status

    ! const
    character(len=*), parameter  ::  rname = mname//'/put_xedges'

    ! local
    real,dimension(:,:,:,:),pointer     :: rm,rxm,rym,rzm,rmc,rxmc,rymc,rzmc
    real,dimension(:,:,:),  pointer     :: m,mc
    integer            :: child,ichild,j,jp,l,lp,jmc,lmc,n,i
    integer            :: xref_,yref_,zref_
    real               :: yzref,xyzref
#ifdef with_budgets
    real               :: sum_old(ntracet),sum_new(ntracet),sum_old_all(ntracet),sum_new_all(ntracet)
    integer            :: nzne, nzne_v,nznem
#endif
    real,dimension(:,:),allocatable :: toc0,toc1,tocm,tocm1
    integer :: communicator,root_id,nglob,lmr,nt

    ! start

    call GO_Timer_Start( itim_put_xedges, status )
    IF_NOTOK_RETURN(status=1)

    which_par=previous_par(region)

    !WP! make sure parents are on same domain
    call Par_Check_Domain( region, 'c', which_par )

    if ( which_par == 'tracer' .and. ntracetloc == 0 ) return
    if ( which_par == 'levels' .and. lmloc == 0 ) return

    ichild = 0
    do while(ichild<children(region,0))
       ichild = ichild + 1
       child = children(region,ichild)
       xref_ = xref(child)/xref(region)
       yref_ = yref(child)/yref(region)
       zref_ = zref(child)/zref(region)

       if (im(child)/xref(child)==im(1)) then
          ! this child is [0,360] degrees wide - skip it since
          ! periodic boundary conditions will be used
          if (okdebug)  &
               print *,'put_xedges: child ',child,' skipped'
          cycle
       endif

       if ( which_par == 'tracer' ) then

          m => m_dat(region)%data
          rm => mass_dat(region)%rm_t
          rxm => mass_dat(region)%rxm_t
          rym => mass_dat(region)%rym_t
          rzm => mass_dat(region)%rzm_t

          mc => m_dat(child)%data
          rmc => mass_dat(child)%rm_t
          rxmc => mass_dat(child)%rxm_t
          rymc => mass_dat(child)%rym_t
          rzmc => mass_dat(child)%rzm_t

          lmr = lm(region)
          nt=ntracetloc
#ifdef with_budgets
#ifdef MPI
          communicator=com_trac  !WP! assign com_trac as communicator
          root_id=root_t
#endif
#endif

       else if ( which_par == 'levels' ) then

          stop 'advecx not with new meteo yet'
          !m => mass_dat(region)%m_k
          rm => mass_dat(region)%rm_k
          rxm => mass_dat(region)%rxm_k
          rym => mass_dat(region)%rym_k
          rzm => mass_dat(region)%rzm_k

          !mc => mass_dat(child)%m_k
          rmc => mass_dat(child)%rm_k
          rxmc => mass_dat(child)%rxm_k
          rymc => mass_dat(child)%rym_k
          rzmc => mass_dat(child)%rzm_k

          lmr = lmloc
          nt=ntracet
#ifdef with_budgets
#ifdef MPI
          communicator=com_lev  !WP! assign com_lev as communicator
          root_id=root_k
#endif
#endif
       end if

       jmc = jm(child)
       lmc = lm(child)

       allocate(toc0(jmc,lmc))
       allocate(toc1(jmc,lmc))
       allocate(tocm(jmc,lmc))
       allocate(tocm1(jmc,lmc))

       yzref = 1./(yref_*zref_)
       xyzref = 1./(xref_*yref_*zref_)

       ! ~~ mass (all levels on each pe)

       ! loop through the cells of x-walls of the child
       do l=1,lmc
          lp = lbeg(child) + (l-1)/zref_
          do j=1,jmc
             jp = jbeg(child) + (j-1)/yref_
             toc0(j,l) = m(ibeg(child)-1,jp,lp)*yzref
             toc1(j,l) = m(ibeg(child)  ,jp,lp)*xyzref
             tocm(j,l) = m(iend(child)  ,jp,lp)*xyzref
             tocm1(j,l) = m(iend(child)+1,jp,lp)*yzref
          end do
       end do
       !     copy info to the child ....
       mc(          0,1:jmc,1:lmc) = toc0
       mc(im(child)+1,1:jmc,1:lmc) = tocm1
       do i=1,xref_
          mc(i            , 1:jmc , 1:lmc ) = toc1
          mc(im(child)+1-i, 1:jmc , 1:lmc ) = tocm
       end do
       nullify(mc)

       ! ~~ tracers

#ifdef with_budgets
       sum_old = 0.0
       sum_old_all = 0.0
       sum_new = 0.0
       sum_new_all = 0.0
#endif

       ! loop over tracers:
       do n = 1, nt

         ! global tracer index:
         nglob = offsetn + n

         do l=1,lmc
             lp = lbeg(child) + (l-1)/zref_
             do j=1,jmc
                jp = jbeg(child) + (j-1)/yref_
                toc0(j,l) = rm(ibeg(child)-1,jp,lp,n)*yzref
                toc1(j,l) = rm(ibeg(child)  ,jp,lp,n)*xyzref
                tocm(j,l) = rm(iend(child)  ,jp,lp,n)*xyzref
                tocm1(j,l) =rm(iend(child)+1,jp,lp,n)*yzref
             end do
          end do

          rmc(          0,1:jmc,1:lmc,n) = toc0
          rmc(im(child)+1,1:jmc,1:lmc,n) = tocm1

          ! calculate the mass that is associated with the 'implicit'
          !advection that is done here....

          do l=1,lmc
             do j=1,jmc
#ifdef with_budgets
                sum_new(nglob) = sum_new(nglob) +  xref_*(toc1(j,l) + tocm(j,l))
#endif
                do i=1,xref_
#ifdef with_budgets
                   sum_old(nglob) = sum_old(nglob) + rmc(i,j,l,n) + rmc(im(child)-i+1,j,l,n)
#endif
                   rmc(i             ,j,l,n) = toc1(j,l)
                   rmc(im(child)+1-i ,j,l,n) = tocm(j,l)
                end do
             end do
          end do

       end do

#ifdef with_budgets

#ifdef MPI
       !WP! only on one when in tracer
       if ( which_par == 'tracer' .and. myid /= pe_first_tracer ) sum_new=0.0
       !WP! only on one when in tracer
       if ( which_par == 'tracer' .and. myid /= pe_first_tracer ) sum_old=0.0
#endif

       !WP! When in levels, each sum_new and sum_old carries the full sum
       !WP! of each level for tracer one.
       !WP! When in tracer,sum_new and sum_old is zero on all PE's
       !WP! except PE_first_tracer

#ifdef MPI
       !add all
       call mpi_allreduce(sum_new,sum_new_all,1,my_real,mpi_sum,communicator,ierr)
       call mpi_allreduce(sum_old,sum_old_all,1,my_real,mpi_sum,communicator,ierr)
#else
       sum_new_all = sum_new
       sum_old_all = sum_old
#endif
       sum_advection(child,:) = sum_advection(child,:) + sum_new_all - sum_old_all

#ifdef MPI
       if ( okdebug .and. myid == root_id ) &
            print *, 'put_xedges: region, sum_advection ', &
            child,sum_advection(child,:)
#else
       if ( okdebug ) print *, 'put_xedges: region, sum_advection ', &
            child,sum_advection(child,:)
#endif

#endif

       nullify(rmc)

       do n=1,nt
          do l=1,lmc
             lp = lbeg(child) + (l-1)/zref_
             do j=1,jmc
                jp = jbeg(child) + (j-1)/yref_
                toc0(j,l) = rxm(ibeg(child)-1,jp,lp,n)*yzref
                toc1(j,l) = rxm(ibeg(child)  ,jp,lp,n)*yzref   !note the yzref...
                tocm(j,l) = rxm(iend(child)  ,jp,lp,n)*yzref
                tocm1(j,l) =rxm(iend(child)+1,jp,lp,n)*yzref
             end do
          end do

          rxmc(          0,1:jmc,1:lmc,n) = toc0
          rxmc(im(child)+1,1:jmc,1:lmc,n) = tocm1
          do i=1,xref_
             rxmc(i            , 1:jmc , 1:lmc ,n) = toc1
             rxmc(im(child)+1-i, 1:jmc , 1:lmc ,n) = tocm
          end do

       end do

       nullify(rxmc)

       do n=1,nt
          do l=1,lmc
             lp = lbeg(child) + (l-1)/zref_
             do j=1,jmc
                jp = jbeg(child) + (j-1)/yref_
                toc0(j,l) = rym(ibeg(child)-1,jp,lp,n)*yzref
                toc1(j,l) = rym(ibeg(child)  ,jp,lp,n)*yzref   !note the yzref...
                tocm(j,l) = rym(iend(child)  ,jp,lp,n)*yzref
                tocm1(j,l) =rym(iend(child)+1,jp,lp,n)*yzref
             end do
          end do
          rymc(          0,1:jmc,1:lmc,n) = toc0
          rymc(im(child)+1,1:jmc,1:lmc,n) = tocm1
          do i=1,xref_
             rymc(i            , 1:jmc , 1:lmc ,n) = toc1
             rymc(im(child)+1-i, 1:jmc , 1:lmc ,n) = tocm
          end do
       end do
       nullify(rymc)

       do n=1,nt
          do l=1,lmc
             lp = lbeg(child) + (l-1)/zref_
             do j=1,jmc
                jp = jbeg(child) + (j-1)/yref_
                toc0(j,l) = rzm(ibeg(child)-1,jp,lp,n)*yzref
                toc1(j,l) = rzm(ibeg(child)  ,jp,lp,n)*yzref   !note the yzref...
                tocm(j,l) = rzm(iend(child)  ,jp,lp,n)*yzref
                tocm1(j,l) =rzm(iend(child)+1,jp,lp,n)*yzref
             end do
          end do

          rzmc(          0,1:jmc,1:lmc,n) = toc0
          rzmc(im(child)+1,1:jmc,1:lmc,n) = tocm1
          do i=1,xref_
             rzmc(i            , 1:jmc , 1:lmc ,n) = toc1
             rzmc(im(child)+1-i, 1:jmc , 1:lmc ,n) = tocm
          end do

       end do
       nullify(rzmc)
       deallocate(toc0)
       deallocate(toc1)
       deallocate(tocm)
       deallocate(tocm1)

       nullify(rm)
       nullify(rxm)
       nullify(rym)
       nullify(rzm)
       nullify(m)
    end do

    call GO_Timer_End( itim_put_xedges, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine put_xedges



end module advectx
