!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module advecty

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------------

  private

  public :: AdvectY_Init, AdvectY_Done
  public :: advectyzoom, put_yedges


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'AdvectY'


  ! --- local ------------------------------------

  integer    ::  itim_zoom
  integer    ::  itim_dynam
  integer    ::  itim_compress
  integer    ::  itim_uncompress
  integer    ::  itim_put_yedges


contains


  ! ====================================================================


  subroutine AdvectY_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/AdvectY_Init'

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_zoom, 'advecty zoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynam, 'advecty dynam', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_compress, 'advecty compress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_uncompress, 'advecty uncompress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_put_yedges, 'advecty put yedges', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine AdvectY_Init


  ! ***


  subroutine AdvectY_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/AdvectY_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine AdvectY_Done


  ! ***


  !
  ! set parameters for advecty
  ! written by patrick berkvens and mike botchev, march-june 1999
  ! updated and modified by MK, dec 2002
  !

  subroutine advectyzoom( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,          only : xref, yref, zref, tref, im, jm, lm
    use dims,          only : rstatus => status, nsplitsteps, n_operators, zoom2D
    use dims,          only : touch_sp, touch_np, okdebug
    use dims,          only : parent, splitorderzoom, zero, xcyc
#ifdef with_budgets
    use budget_global, only : sum_advection
    use chem_param,    only : ntracet
#endif
    use global_data,   only : wind_dat, mass_dat
    use toolbox,       only : escape_tm
#ifdef MPI
    use mpi_const,     only : my_real,com_lev,ierr, &
                              root_t,mpi_sum,pe_first_tracer
#endif
    use partools     , only : myid, lmloc, ntracetloc

    ! --- in/out ---------------------------------

    integer,intent(in)        ::  region
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/advectyzoom'

    ! --- local ----------------------------------

    real,dimension(:,:,:),pointer               :: bm
    real,dimension(:,:,:,:),pointer             :: rm
    real,dimension(:,:),allocatable             :: bm0,bm1

    integer            :: is,ie,ls,le,n,q
    integer            :: imr,jmr,lmr,tref_,xref_,yref_,zref_
    logical            :: y_encountered
    character(len=1)   :: dir
#ifdef with_budgets
    real               :: sum_old(ntracet), sum_new(ntracet), sum_old_all(ntracet), sum_new_all(ntracet)
    integer            :: itr
#endif

    ! --- start ----------------------------------

    !WP! only PE's with nonzero lmloc proceed
    if ( ntracetloc==0) return

    call GO_Timer_Start( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    call put_yedges( region, status )
    IF_NOTOK_RETURN(status=1)

    tref_ = tref(region)/tref(parent(region))
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region);  jmr = jm(region);  lmr = lm(region)


    allocate(bm0(0:im(region)+1,0:lmr+1))
    allocate(bm1(0:im(region)+1,0:lmr+1))

    bm => wind_dat(region)%bm_t
    rm =>  mass_dat(region)%rm_t

#ifdef with_budgets
    do itr = 1, ntracet
        sum_old(itr) = sum(rm(1:imr,1:jmr,1:lmr,itr))
    end do
#endif

    ! determine the scope for advecty:

    if ( region == 1 ) then
       xref_ = 0; zref_ = 0  ! to have is/ie and ls/le properly computed
    end if
    ! find q - the place in the splitorderzoom
    q=rstatus(region)/((nsplitsteps/2)*tref_)
    ! corresponding to the begining of the last
    ! triple x-y-z of the parent
    y_encountered=.false.

    ! now track four following steps from q     !wp! changed to four steps

    do n=1,n_operators

       dir=splitorderzoom(region,q*(nsplitsteps/2)*tref_+n)

       select case(dir)
       case('x')
          if  ((.not.y_encountered).or.(xcyc(region) == 1)) then
             is=1                             ! x-substep is before y =>
             ie=im(region)                    ! full i-scope
          else
             is=xref_+1                       ! y-substep is before x =>
             ie=im(region)-xref_              ! restricted i-scope
          end if
       case('y')
          y_encountered=.true.
       case('z')
          if (.not.y_encountered) then
             ls=1                             ! z-substep is before y =>
             le=lmr                    ! full l-scope
          else
             ls=zref_+1                       ! y-substep is before z =>
             le=lmr-zref_              ! restricted l-scope
          end if
          if (zoom2D) then
             ls=1; le=lmr   !WP! this is always the case
          end if
       case('c')
       case('v')
       case('d')
       case('s')
       case default
          print *,'advectyzoom: strange value in splitorderzoom(',region,',',  &
               q*(nsplitsteps/2)*tref_+n,'): ',q
          call escape_tm('advectyzoom: Error in do_next3')
       end select

    end do

#ifdef slopes
#ifndef secmom
    if ( mod(rstatus(region),(nsplitsteps/2)*tref_) >= n_operators )  then
       ! if (1) more than fourth substep and (2) slope scheme then
       ! at this substep no coarse fluxes will be applied

       ! for slopes only: zero out velocities via the y-edges of the region
       ! to assure that no fluxes will be applied

       if ( okdebug ) print *,'advectyzoom: zeroing out outer ', &
            'mass fluxes (yref)',yref_
       if ( touch_sp(region) /= 1) then
          bm0 = bm(:,      yref_,:);   bm(:,      yref_,:) = zero
       else
          if ( okdebug) print *,'advectyzoom: skipping zeroing outer ', &
               'mass flux at SP'
       end if
       if ( touch_np(region) /= 1) then
          bm1 = bm(:,jmr-yref_+2,:);   bm(:,jmr-yref_+2,:) = zero
       else
          if ( okdebug) print *,'advectyzoom: skipping zeroing outer ', &
               'mass flux at NP'
       end if

    end if
#endif
#endif

    call dynamv( region, is,ie,ls,le, status )
    IF_NOTOK_RETURN(status=1)

#ifdef slopes
#ifndef secmom
    if ( mod(rstatus(region),(nsplitsteps/2)*tref_) >= n_operators ) then

       ! for slopes only: recreate velocities via the x-edges
       if ( touch_sp(region) /= 1 ) bm(:,      yref_,:) = bm0
       if ( touch_np(region) /= 1 ) bm(:,jmr-yref_+2,:) = bm1

    end if
#endif
#endif

#ifdef with_budgets
    do itr = 1, ntracet
        sum_new(itr) = sum(rm(1:imr,1:jmr,1:lmr,itr))
    end do
#ifdef MPI
    if ( myid==pe_first_tracer ) &
         sum_advection(region,:) = sum_advection(region.:) + sum_new - sum_old
    if ( okdebug .and. myid == root_t ) &
         print *, 'advectyzoom: advect_y, region, sum_advection ', &
         region,sum_advection(region,:),sum_new,sum_old
#else
    sum_advection(region,:) = sum_advection(region,:) + sum_new - sum_old
    if ( okdebug) print *, 'advectyzoom: advect_y, region, sum_advection ', &
         region,sum_advection(region,:),sum_new,sum_old
#endif
#endif

    deallocate(bm0)
    deallocate(bm1)
    nullify(bm)

    call GO_Timer_End( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectyzoom


  !-----------------------------------------------------------------------
  !
  !****   dynamv          - south-north tracer transport  v 9.1
  !
  !       programmed by           mh      mpi HH          23-feb-1994
  !
  !       purpose
  !       -------
  !       calculate amount of tracer moved in a south-north advection
  !       substep
  !
  !       interface
  !       ---------
  !       call dynamv
  !
  !       method
  !       ------
  !       slopes scheme
  !
  !       externals
  !       ---------
  !       none
  !
  !       reference
  !       ---------
  !       Russel and Lerner, 1979
  !-----------------------------------------------------------------------
  !       fixed bug in calculating maximum courant number
  !                               mh, Thu, Feb 24, 1994 12:49:27
  !     included code for limits of slopes to prevent negative tracer
  !     masses                   mh, 20-jun-1994
  !
  !     zoom version written by mike botchev, march-june 1999
  !-----------------------------------------------------------------------

  subroutine dynamv( region, is,ie,ls,le, status )

    use omp_lib
    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : xref, yref, zref
    use dims,        only : im, jm, lm
    use dims,        only : okdebug, parent, nregions
!    use dims,        only : xi, nxi
    use dims,        only : touch_sp, touch_np, limits, nloop_max, zero, one
    use redgridZoom, only : nred, jred, imredj, clustsize
    use global_data, only : wind_dat
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use zoom_tools,  only : mix_edges
    use toolbox,     only : escape_tm
#ifdef with_budgets
    use budget_global, only : budget_flux, jflux1, jflux2
    use budget_global, only : apply_budget_global_flux => apply_budget_global
    use partools     , only : tracer_loc
#endif
    use partools   , only : ntracetloc, myid
    use chem_param,  only : ntracet, ra

    ! --- in/out ---------------------------------

    integer,intent(in)  ::  region
    integer,intent(in)  ::  is
    integer,intent(in)  ::  ie
    integer,intent(in)  ::  ls
    integer,intent(in)  ::  le
    integer,intent(out) ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dynamv'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer           :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer           :: m,bm
    real,dimension(:,:),allocatable           :: mnew
    real,dimension(:,:),allocatable           :: f,pf,fx,fz
    integer                                   :: i,j,je,js,l,n,iee,iss
    integer                                   :: imr,imr2,jmr,lmr
    integer                                   :: xref_,yref_,zref_
    real                                      :: sfs,sfzs,sfn,sfzn
    real                                      :: beta,mxval
    integer,dimension(3)                      :: mxloc
    integer                                   :: iloop, nloop, nglob, offsetn
    integer,parameter                         :: max_nloop = 6
    real,dimension(:,:), allocatable :: mx
    logical                          :: cfl_ok
    integer        :: lrg, redfact, ixe, ixs
    real           :: summ

    ! --- start ----------------------------------

    call GO_Timer_Start( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_limits
    if ( .not. limits ) stop 'macro "with_limits" defined but variable "limits" not'
#else
    if ( limits ) stop 'macro "with_limits" not defined but variable "limits" is'
#endif

    if ( okdebug ) print *,'dynamv: region=',region,' is,ie,ls,le=',is,ie,ls,le
    if ( ( region < 0 ) .or. ( region > nregions ) )  then
       call escape_tm('dynamv: STOP, illegal number of region !!!')
    end if

    ! compute refinement factors with respect to the parent

    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))

    imr=im(region);jmr=jm(region);lmr=lm(region)

    call compress_yedges( region, is,ie,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    bm => wind_dat(region)%bm_t

    ! check is,ie,ls,le:
    if ( ( is /= xref_+1 ) .and. ( is /= 1 ) ) &
         call escape_tm( 'dynamv: Wrong value for IS in dynamv')
    if ( ( ie /= imr-xref_ ) .and. ( ie /= imr ) ) &
         call escape_tm( 'dynamv: Wrong value for IE in dynamv')
    if ( ( ls /= zref_+1 ) .and. ( ls /= 1 ) )  &
         call escape_tm( 'dynamv: Wrong value for LS in dynamv')
    if ( ( le /= lmr-zref_ ) .and. ( le /= lmr ) ) &
         call escape_tm( 'dynamv: Wrong value for LE in dynamv')

    ! compute js/je -- cells is:ie,js:je:ls:le will be updated
    if (region==1) then
       js = 2
       je = jmr-1
    else
       js = yref_
       je = jmr-yref_+1
       if ( touch_sp(region) == 1) js = 2    !cmk....new option.. reduced grid
       if ( touch_np(region) == 1) je = jmr-1  !exclude poles
    end if

    call GO_Timer_Start( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over tracers and vertical layers

    !$XOMP   shared ( touch_sp, touch_np ) &
    !$OMP PARALLEL &
    !$OMP   default (none) &
#ifdef with_budgets
    !$OMP   shared  (apply_budget_global_flux, budget_flux, jflux1, jflux2,ra) &
    !$OMP   shared  (tracer_loc) &
    !$OMP   private (nglob) &
#endif
    !$OMP   shared  (region ) &
    !$OMP   shared  (imr, jmr) &
    !$OMP   shared  (is, ie, js, je, ls, le) &
    !$OMP   shared  (okdebug) &
    !$OMP   shared  (ntracetloc) &
    !$OMP   shared  (myid) &
    !$OMP   shared  (rm, rym, rxm, rzm ) &
    !$OMP   shared  (m, bm ) &
    !$OMP   private (i, j, l, f, pf, fx, fz) &
    !$OMP   private (mnew, beta)

    allocate (f(im(region),0:jm(region))) ; f = 0.0
    allocate (pf(im(region),0:jm(region))) ; pf = 0.0
    allocate (fx(im(region),0:jm(region))) ; fx = 0.0
    allocate (fz(im(region),0:jm(region))) ; fz = 0.0
    allocate (mnew(im(region),jm(region))) ; mnew = 0.0

    ! loop over levels:
    !$OMP   DO
    do l=ls,le
      ! calculate new air mass distribution

      if (region==1) then
         mnew(1:imr,1:jmr) = m(1:imr,1:jmr,l) + &
              bm(1:imr,1:jmr,l) - bm(1:imr,2:jmr+1,l)
      else if ( touch_sp(region) == 1) then
         mnew(1:imr,1:je) = m(1:imr,1:je,l) + &
              bm(1:imr,1:je,l) - bm(1:imr,2:je+1,l)
      else if ( touch_np(region) == 1) then
         mnew(1:imr,js:jmr) = m(1:imr,js:jmr,l) + &
              bm(1:imr,js:jmr,l) - bm(1:imr,js+1:jmr+1,l)
      else
         mnew(is:ie,js:je) = m(is:ie,js:je,l) + &
              bm(is:ie,js  :je  ,l) - bm(is:ie,js+1:je+1,l)
      end if

      do n=1,ntracetloc

        ! if requested limit meridional slopes such that no negative
#ifdef with_limits
        rym(:,:,l,n) = max( min( rym(:,:,l,n), rm(:,:,l,n) ), -rm(:,:,l,n) )
#endif

        do j=js+1,je
            do i=is,ie
               if (bm(i,j,l)>=zero) then
                  beta=bm(i,j,l)/m(i,j-1,l)
                  f(i,j)=beta*(rm(i,j-1,l,n)+(one-beta)*rym(i,j-1,l,n))
                  pf(i,j)=bm(i,j,l)*(beta*beta*rym(i,j-1,l,n)-3.*f(i,j))
                  fx(i,j)=beta*rxm(i,j-1,l,n)
                  fz(i,j)=beta*rzm(i,j-1,l,n)
               else
                  beta=bm(i,j,l)/m(i,j,l)
                  f(i,j)=beta*(rm(i,j,l,n)-(one+beta)*rym(i,j,l,n))
                  pf(i,j)=bm(i,j,l)*(beta*beta*rym(i,j,l,n)-3.*f(i,j))
                  fx(i,j)=beta*rxm(i,j,l,n)
                  fz(i,j)=beta*rzm(i,j,l,n)
               end if
            end do
         end do

         !compute boundary fluxes
         if ( region == 1.or.touch_sp(region) == 1) then
            do i=is,ie
               fz(i,1)=zero
               fx(i,1)=zero
               pf(i,1)=zero
               f(i,1)=zero
               if (bm(i,2,l)>=zero) then
                  beta=bm(i,2,l)/m(i,1,l)
                  f(i,2)=beta*rm(i,1,l,n)
                  pf(i,2)=-3.*bm(i,2,l)*f(i,2)
                  fx(i,2)=zero
                  fz(i,2)=beta*rzm(i,1,l,n)
               else
                  beta=bm(i,2,l)/m(i,2,l)
                  f(i,2)=beta*(rm(i,2,l,n)-(one+beta)*rym(i,2,l,n))
                  pf(i,2)=bm(i,2,l)*(beta*beta*rym(i,2,l,n)-3.*f(i,2))
                  fx(i,2)=beta*rxm(i,2,l,n)
                  fz(i,2)=beta*rzm(i,2,l,n)
               end if
            end do
         else   !zoom region not touching south pole
            j = js
            do i=is,ie   !no reduced grid allowed
               if (bm(i,j,l)>=zero) then
                  beta=bm(i,j,l)/m(i,j-1,l)
                  f(i,j)=beta*(rm(i,j-1,l,n)+(one-beta)*rym(i,j-1,l,n))
                  pf(i,j)=bm(i,j,l)*(beta*beta*rym(i,j-1,l,n)-3.*f(i,j))
                  fx(i,j)=beta*rxm(i,j-1,l,n)
                  fz(i,j)=beta*rzm(i,j-1,l,n)
               else
                  beta=bm(i,j,l)/m(i,j,l)
                  f(i,j)=beta*(rm(i,j,l,n)-(one+beta)*rym(i,j,l,n))
                  pf(i,j)=bm(i,j,l)*(beta*beta*rym(i,j,l,n)-3.*f(i,j))
                  fx(i,j)=beta*rxm(i,j,l,n)
                  fz(i,j)=beta*rzm(i,j,l,n)
               end if
            end do
         end if ! compute boundary fluxes south pole...
         if ( region == 1.or.touch_np(region) == 1) then       ! north pole
            do i=is,ie
               if (bm(i,jmr,l)>=zero) then
                  beta=bm(i,jmr,l)/m(i,jmr-1,l)
                  f(i,jmr)=beta*(rm(i,jmr-1,l,n)+(one-beta)*rym(i,jmr-1,l,n))
                  pf(i,jmr)=bm(i,jmr,l)*(beta*beta*rym(i,jmr-1,l,n) - &
                       3.*f(i,jmr))
                  fx(i,jmr)=beta*rxm(i,jmr-1,l,n)
                  fz(i,jmr)=beta*rzm(i,jmr-1,l,n)
               else
                  beta=bm(i,jmr,l)/m(i,jmr,l)
                  ! for solid-body test and polar mixing
                  f(i,jmr)=beta*rm(i,jmr,l,n)
                  pf(i,jmr)=-3.*bm(i,jmr,l)*f(i,jmr)
                  fx(i,jmr)=zero
                  ! for solid-body test and polar mixing
                  fz(i,jmr)=beta*rzm(i,jmr,l,n)
               end if
            end do
         else   !zoom region not touching north pole
            j = je+1
            do i=is,ie   !no reduced grid allowed
               if (bm(i,j,l)>=zero) then
                  beta=bm(i,j,l)/m(i,j-1,l)
                  f(i,j)=beta*(rm(i,j-1,l,n)+(one-beta)*rym(i,j-1,l,n))
                  pf(i,j)=bm(i,j,l)*(beta*beta*rym(i,j-1,l,n)-3.*f(i,j))
                  fx(i,j)=beta*rxm(i,j-1,l,n)
                  fz(i,j)=beta*rzm(i,j-1,l,n)
               else
                  beta=bm(i,j,l)/m(i,j,l)
                  f(i,j)=beta*(rm(i,j,l,n)-(one+beta)*rym(i,j,l,n))
                  pf(i,j)=bm(i,j,l)*(beta*beta*rym(i,j,l,n)-3.*f(i,j))
                  fx(i,j)=beta*rxm(i,j,l,n)
                  fz(i,j)=beta*rzm(i,j,l,n)
               end if
            end do
         end if ! compute boundary fluxes north pole...

#ifdef with_budgets
         !cmk finished computing fluxes...now apply...
         ! first compute INCOMING fluxes from the south...
         if (apply_budget_global_flux) then
            nglob = tracer_loc(n)
            do i=is,ie
               ! AJS: something wrong here: values from f are used that are not set ...
               ! AJS: adhoc solution: set f to zero after allocation
               budget_flux(region)%flux_y1(i,l,nglob) = budget_flux(region)%flux_y1(i,l,nglob) + &
                 f(i,jflux1(region))*1e3/ra(nglob)
               budget_flux(region)%flux_y2(i,l,nglob) = budget_flux(region)%flux_y2(i,l,nglob) + &
                 f(i,jflux2(region))*1e3/ra(nglob)
            enddo
         endif
#endif
         do j=js,je
            do i=is,ie
               rm(i,j,l,n) =rm(i,j,l,n)+(f(i,j)-f(i,j+1))
               rym(i,j,l,n)=rym(i,j,l,n)+(pf(i,j)-pf(i,j+1)            &
                    -       (bm(i,j,l)-bm(i,j+1,l))*rym(i,j,l,n)       &
                    +       3.*((bm(i,j,l)+bm(i,j+1,l))*rm(i,j,l,n)    &
                    -       (f(i,j)+f(i,j+1))*m(i,j,l)))/mnew(i,j)
#ifdef with_limits
               rym(i,j,l,n) = max( min(rym(i,j,l,n), rm(i,j,l,n)), -rm(i,j,l,n) )
#endif
               rxm(i,j,l,n)=rxm(i,j,l,n)+(fx(i,j)-fx(i,j+1))
               rzm(i,j,l,n)=rzm(i,j,l,n)+(fz(i,j)-fz(i,j+1))
            end do
         end do !forall

         ! poles were excluded.....
         if ( region==1 .or. touch_sp(region) == 1 ) then
            rm(1:imr,  1,l,n) = rm(1:imr,  1,l,n) - f(1:imr,2)
            rzm(1:imr,1  ,l,n) = rzm(1:imr,  1,l,n) - fz(1:imr,2)
         end if
         if ( region==1 .or. touch_np(region) == 1 ) then
            rm(1:imr,jmr,l,n) = rm(1:imr,jmr,l,n) + f(1:imr,jmr)
            rzm(1:imr,jmr,l,n) = rzm(1:imr,jmr,l,n) + fz(1:imr,jmr)
         end if

      end do !      & n-loop

      ! store new air mass in m array

      if ( region == 1) then
         m(1:imr,1:jmr,l) = mnew(1:imr,1:jmr)
      else if ( touch_sp(region) == 1) then
         m(1:imr,1:je,l) = mnew(1:imr,1:je)
      else if ( touch_np(region) == 1) then
         m(1:imr,js:jmr,l) = mnew(1:imr,js:jmr)
      else
         m(is:ie,js:je,l) = mnew(is:ie,js:je)
      end if
      !

    end do! end of l-loop over vertical layers....
    !$OMP   END DO

    deallocate(mnew)
    deallocate( f)
    deallocate(pf)
    deallocate(fx)
    deallocate(fz)

    !$OMP END PARALLEL

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(bm)

    call GO_Timer_End( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    ! *

    call uncompress_yedges( region, is,ie,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    ! *

    call mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    ! *

    call GO_Timer_End( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine dynamv


  !
  ! this is for slope only: condense data at the y-edges of the zoom region
  ! to allow for uniform work in dynamv
  ! written by mike botchev, march-june 1999
  ! modified by MK, dec 2002
  !

  subroutine compress_yedges( region, is,ie,ls,le, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : yref, parent, jm, okdebug, touch_sp, touch_np
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use partools   , only : ntracetloc
    use chem_param,  only : ntracet

    ! input/output
    integer,intent(in)    ::  region
    integer,intent(in)    ::  is
    integer,intent(in)    ::  ie
    integer,intent(in)    ::  ls
    integer,intent(in)    ::  le
    integer, intent(out)  ::  status

    ! conts
    character(len=*), parameter  ::  rname = mname//'/compress_yedges'

    ! local
    real,dimension(:,:,:,:),pointer           :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer           :: m
    integer                                   :: jmr,i,l,n,yref_

    ! start

    call GO_Timer_Start( itim_compress, status )
    IF_NOTOK_RETURN(status=1)

    yref_ = yref(region)/yref(parent(region))
    jmr = jm(region)
    if ((yref_==1).or.(region==1)) then
       if ( okdebug) print *,'compress_yedges: no refinement, nothing to do'
       call GO_Timer_End( itim_compress, status )
       IF_NOTOK_RETURN(status=1)
       return
    end if
    if ( okdebug) then
       print *,'compress_yedges: region=',region,' jmr=',jmr
       print *,'compress_yedges: cells to be updated: ', &
            yref_-1,yref_,jmr-yref_+1,jmr-yref_+2
       if ( touch_sp(region) == 1) &
            print *,'compress_yedges: SP will be skipped here'
       if ( touch_np(region) == 1) &
            print *,'compress_yedges: NP will be skipped here'
    end if

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    !$XsafdOMP   shared ( touch_sp, touch_np ) &
    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( jmr ) &
    !$OMP   shared ( yref_ ) &
    !$OMP   shared ( is, ie, ls, le ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m, rm, rxm, rym, rzm ) &
    !$OMP   private ( i, l, n )
    !$OMP   DO
    do l=ls,le
       do i=is,ie

          if ( touch_sp(region) /= 1) then
             m (i,yref_,  l) = sum(m(i,1:yref_,l))
             m (i,yref_-1,l) = m (i,0,l)
          end if
          if ( touch_np(region) /= 1) then
             m (i,jmr-yref_+1,l) = sum(m(i,jmr-yref_+1:jmr,l))
             m (i,jmr-yref_+2,l) = m (i,jmr+1,l)
          end if
          do n=1,ntracetloc
             if ( touch_sp(region) /= 1) then
                rm (i,yref_,  l,n) = sum(rm(i,1:yref_,l,n))
                rm (i,yref_-1,l,n) = rm (i,0,l,n)
                rxm(i,yref_-1,l,n) = rxm(i,0,l,n)
                rym(i,yref_-1,l,n) = rym(i,0,l,n)
                rzm(i,yref_-1,l,n) = rzm(i,0,l,n)
             end if
             if ( touch_np(region) /= 1) then
                rm (i,jmr-yref_+1,l,n) = sum(rm(i,jmr-yref_+1:jmr,l,n))
                rm (i,jmr-yref_+2,l,n) = rm (i,jmr+1,l,n)
                rxm(i,jmr-yref_+2,l,n) = rxm(i,jmr+1,l,n)
                rym(i,jmr-yref_+2,l,n) = rym(i,jmr+1,l,n)
                rzm(i,jmr-yref_+2,l,n) = rzm(i,jmr+1,l,n)
             end if
          end do !forall

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

  end subroutine compress_yedges


  ! ***


  !
  ! distribute data at the y-edges of the zoom region
  ! over the whole interface cell
  ! written by mike botchev, march-june 1999
  !

  subroutine uncompress_yedges( region, is,ie,ls,le, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : yref, parent, jm, okdebug, touch_sp, touch_np
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use partools   , only : ntracetloc
    use chem_param,  only : ntracet

    ! input/output
    integer, intent(in)   ::  region
    integer, intent(in)   ::  is
    integer, intent(in)   ::  ie
    integer, intent(in)   ::  ls
    integer, intent(in)   ::  le
    integer, intent(out)  ::  status

    ! conts
    character(len=*), parameter  ::  rname = mname//'/uncompress_yedges'

    ! local
    real,dimension(:,:,:,:),pointer           :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer           :: m
    integer                                   :: jmr,i,l,n,yref_
    real                                      :: m_edge,rm_edge

    ! start

    call GO_Timer_Start( itim_uncompress, status )
    IF_NOTOK_RETURN(status=1)

    yref_ = yref(region)/yref(parent(region))
    jmr = jm(region)
    if ((yref_==1).or.(region==1)) then
       if ( okdebug) print *,'uncompress_yedges: no refinement, nothnig to do'
       call GO_Timer_End( itim_uncompress, status )
       IF_NOTOK_RETURN(status=1)
       return
    end if
    if ( okdebug) then
       print *,'uncompress_yedges: region=',region,' jmr=',jmr
       print *,'uncompress_yedges: cells to be updated: ',1,':', &
            yref_,' ',jmr-yref_+1,':',jmr
       if ( touch_sp(region) == 1) &
            print *,'uncompress_yedges: SP will be skipped here'
       if ( touch_np(region) == 1) &
            print *,'uncompress_yedges: NP will be skipped here'
    end if

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    !$XOMP   shared ( touch_sp, touch_np ) &
    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( jmr ) &
    !$OMP   shared ( yref_ ) &
    !$OMP   shared ( is, ie, ls, le ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m, rm, rxm, rym, rzm ) &
    !$OMP   private ( i, l, n ) &
    !$OMP   private ( m_edge, rm_edge )
    !$OMP   DO
    do l=ls,le
       do i=is,ie
          if ( touch_sp(region) /= 1) then
             m_edge = m (i,yref_,l)
             m (i,1:yref_,l) = m_edge/yref_
          end if
          if ( touch_np(region) /= 1) then
             m_edge = m (i,jmr-yref_+1,l)
             m (i,jmr-yref_+1:jmr,l) = m_edge/yref_
          end if

          do n=1,ntracetloc

            if ( touch_sp(region) /= 1 ) then
                rm_edge = rm (i,yref_,l,n)
                rm (i,1:yref_,l,n) = rm_edge/yref_
                !*      rxm(i,1:yref_-1,l,n) = rxm(i,yref_,l,n)
                !*      rym(i,1:yref_-1,l,n) = rym(i,yref_,l,n)
                !*      rzm(i,1:yref_-1,l,n) = rzm(i,yref_,l,n)
                rxm(i,1:yref_,l,n) = rxm(i,yref_,l,n)/yref_
                rym(i,1:yref_,l,n) = rym(i,yref_,l,n)/yref_
                rzm(i,1:yref_,l,n) = rzm(i,yref_,l,n)/yref_

             end if
             if ( touch_np(region) /= 1 ) then
                rm_edge = rm (i,jmr-yref_+1,l,n)
                rm (i,jmr-yref_+1:jmr,l,n) = rm_edge/yref_

                !*      rxm(i,jmr-yref_+2:jmr,l,n) = rxm(i,jmr-yref_+1,l,n)
                !*      rym(i,jmr-yref_+2:jmr,l,n) = rym(i,jmr-yref_+1,l,n)
                !*      rzm(i,jmr-yref_+2:jmr,l,n) = rzm(i,jmr-yref_+1,l,n)
                rxm(i,jmr-yref_+1:jmr,l,n) = rxm(i,jmr-yref_+1,l,n)/yref_
                rym(i,jmr-yref_+1:jmr,l,n) = rym(i,jmr-yref_+1,l,n)/yref_
                rzm(i,jmr-yref_+1:jmr,l,n) = rzm(i,jmr-yref_+1,l,n)/yref_
             end if
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

  end subroutine uncompress_yedges


  !
  ! passes values along the y-boundaries to all the children
  ! written by mike botchev, march-june 1999
  ! adapted for cray-run Maarten Krol, march, 2000
  ! touch_np, touch_sp implemented...MK, nov, 2000
  !

  subroutine put_yedges( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,          only : xref, yref, zref, im, jm, lm
    use dims,          only : touch_sp, touch_np, okdebug, children
    use dims,          only : ibeg, jbeg, lbeg, jend, nregions
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
    use partools     , only : Par_Check_Domain

    ! input/output
    integer,intent(in)                  :: region
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/put_yedges'

    ! local
    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm,rmc,rxmc,rymc,rzmc
    real,dimension(:,:,:),  pointer    :: m, mc

    real,dimension(:,:),allocatable    :: toc0,toc1,tocm,tocm1

    integer :: child, ichild, j, ip, l, lp, imc, lmc, n, i
    integer :: xref_, yref_, zref_
    real    :: xzref, xyzref
#ifdef with_budgets
    real    :: sum_old(ntracet), sum_new(ntracet), sum_old_all(ntracet), sum_new_all(ntracet)
    integer :: nglob
    integer :: communicator, root_id
#endif
    integer :: lmr, nt

    ! start

    call GO_Timer_Start( itim_put_yedges, status )
    IF_NOTOK_RETURN(status=1)

    which_par=previous_par(region)
    !WP! make sure parents are on same domain
    call Par_Check_Domain( region, 'c', which_par )
    if ( which_par == 'tracer' .and. ntracetloc == 0 ) then
      ! GM, 13 may 2024: if returning without doing anything is allowed here, then
      ! the timer should be ended
      call GO_Timer_End( itim_put_yedges, status )
      return
    endif
    if ( which_par == 'levels' .and. lmloc == 0 ) return  !WP!

    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    m => m_dat(region)%data

    ichild = 0
    do while(ichild<children(region,0))
       ichild = ichild + 1
       child = children(region,ichild)
       xref_ = xref(child)/xref(region)
       yref_ = yref(child)/yref(region)
       zref_ = zref(child)/zref(region)

       if ( touch_sp(child) == 1) then
          if (okdebug)  &
               print *,'put_yedges: child ',child,' sp  skipped'
       end if
       if ( touch_np(child) == 1) then
          if (okdebug)  &
               print *,'put_yedges: child ',child,' np  skipped'
       end if

       if ( which_par=='tracer') then

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
       else if ( which_par=='levels') then

          stop 'no advecty with new meteo yet'
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

       imc = im(child)
       lmc = lm(child)

       allocate(toc0(imc,lmc))
       allocate(toc1(imc,lmc))
       allocate(tocm(imc,lmc))
       allocate(tocm1(imc,lmc))

       ! ~~ mass (all levels on each pe)

       xzref = 1./(xref_*zref_)
       xyzref = 1./(xref_*yref_*zref_)

       ! loop through the cells of y-walls of the child
       do l=1,lmc
          lp = lbeg(child) + (l-1)/zref_
          do i=1,imc
             ip = mod(ibeg(child)-1 + (i-1)/xref_,im(region)) + 1
             toc0(i,l) =  m(ip,jbeg(child)-1,lp)*xzref
             toc1(i,l) =  m(ip,jbeg(child)  ,lp)*xyzref
             tocm(i,l) =  m(ip,jend(child)  ,lp)*xyzref
             tocm1(i,l) = m(ip,jend(child)+1,lp)*xzref
          end do
       end do

       if ( touch_sp(child) == 0)  mc(1:imc,0          ,1:lmc) = toc0
       if ( touch_np(child) == 0)  mc(1:imc,jm(child)+1,1:lmc) = tocm1
       do j=1,yref_
          if ( touch_sp(child) == 0)  mc( 1:imc , j,            1:lmc ) = toc1
          if ( touch_np(child) == 0)  mc( 1:imc , jm(child)+1-j,1:lmc ) = tocm
       end do

       nullify(mc)

       ! ~~ tracers

#ifdef with_budgets
       sum_old_all = 0.0
       sum_new = 0.0
       sum_new_all = 0.0
       sum_old = 0.0
#endif

       ! loop over tracers:
       do n = 1, nt

#ifdef with_budgets
         ! global tracer index:
         nglob = offsetn + n
#endif

          do l=1,lmc
             lp = lbeg(child) + (l-1)/zref_
             do i=1,imc
                ip = mod(ibeg(child)-1 + (i-1)/xref_,im(region)) + 1
                toc0(i,l) = rm(ip,jbeg(child)-1,lp,n)*xzref
                toc1(i,l) = rm(ip,jbeg(child)  ,lp,n)*xyzref
                tocm(i,l) = rm(ip,jend(child)  ,lp,n)*xyzref
                tocm1(i,l) =rm(ip,jend(child)+1,lp,n)*xzref
             end do
          end do

          if ( touch_sp(child) == 0)  rmc(1:imc,0          ,1:lmc,n) = toc0
          if ( touch_np(child) == 0)  rmc(1:imc,jm(child)+1,1:lmc,n) = tocm1

          do l=1,lmc
             do i=1,imc

#ifdef with_budgets
                if ( touch_sp(child) == 0) sum_new(nglob) = sum_new(nglob) + yref_*toc1(i,l)
                if ( touch_np(child) == 0) sum_new(nglob) = sum_new(nglob) + yref_*tocm(i,l)
#endif

                do j=1,yref_
                   if ( touch_sp(child) == 0)  then
#ifdef with_budgets
                      sum_old(nglob) = sum_old(nglob) + rmc( i, j, l , n)
#endif
                      rmc(i,j,l,n) = toc1(i,l)
                   end if
                   if ( touch_np(child) == 0 )  then
#ifdef with_budgets
                      sum_old(nglob) = sum_old(nglob) + rmc(i,jm(child)+1-j,l,n)
#endif
                      rmc(i,jm(child)+1-j,l,n) = tocm(i,l)
                   end if
                end do !j
             end do !i
          end do !l
       end do  !n

#ifdef with_budgets

#ifdef MPI
       !WP! only on one when in tracer
       if ( which_par == 'tracer' .and. myid /= pe_first_tracer ) sum_new=0.0
       if ( which_par == 'tracer' .and. myid /= pe_first_tracer ) sum_old=0.0
#endif

       !WP! When in levels, each sum_new and sum_old carries the full
       !WP! sum of each level for tracer one.
       !WP! When in tracer,sum_new and sum_old is zero on all PE's
       !WP! except PE_first_tracer

#ifdef MPI
       call mpi_allreduce(sum_new,sum_new_all,1,my_real, &
            mpi_sum,communicator,ierr)
       call mpi_allreduce(sum_old,sum_old_all,1,my_real, &
            mpi_sum,communicator,ierr)
#else
       sum_new_all = sum_new
       sum_old_all = sum_old
#endif

       sum_advection(child,:) = sum_advection(child,:) + sum_new_all - sum_old_all

#ifdef MPI
       if ( okdebug .and. myid==root_id ) &
            print *, 'put_yedges: region, sum_advection ', &
            child,sum_advection(child,:)
#else
       if ( okdebug) print *, 'put_yedges: region, sum_advection ', &
            child,sum_advection(child,:)
#endif

#endif

       nullify(rmc)

       do n=1,nt
         do l=1,lmc
           lp = lbeg(child) + (l-1)/zref_
           do i=1,imc
              ip = mod(ibeg(child)-1 + (i-1)/xref_,im(region)) + 1
              toc0(i,l) = rxm(ip,jbeg(child)-1,lp,n)*xzref
              toc1(i,l) = rxm(ip,jbeg(child)  ,lp,n)*xzref   !note the xzref...
              tocm(i,l) = rxm(ip,jend(child)  ,lp,n)*xzref
              tocm1(i,l) =rxm(ip,jend(child)+1,lp,n)*xzref
           end do
         end do

         if ( touch_sp(child) == 0 )  rxmc(1:imc,0          ,1:lmc,n) = toc0
         if ( touch_np(child) == 0 )  rxmc(1:imc,jm(child)+1,1:lmc,n) = tocm1
         do j=1,yref_
            if ( touch_sp(child) == 0 )  &
                 rxmc( 1:imc , j            ,1:lmc ,n) = toc1
            if ( touch_np(child) == 0 )  &
                 rxmc( 1:imc , jm(child)+1-j,1:lmc ,n) = tocm
         end do

       end do  !n
       nullify(rxmc)


       do n=1,nt
          do l=1,lmc
             lp = lbeg(child) + (l-1)/zref_
             do i=1,imc
                ip = mod(ibeg(child)-1 + (i-1)/xref_,im(region)) + 1
                toc0(i,l) = rym(ip,jbeg(child)-1,lp,n)*xzref
                toc1(i,l) = rym(ip,jbeg(child)  ,lp,n)*xzref   !note the xzref...
                tocm(i,l) = rym(ip,jend(child)  ,lp,n)*xzref
                tocm1(i,l) =rym(ip,jend(child)+1,lp,n)*xzref
             end do
          end do

          if ( touch_sp(child) == 0)  rymc(1:imc,0          ,1:lmc,n) = toc0
          if ( touch_np(child) == 0)  rymc(1:imc,jm(child)+1,1:lmc,n) = tocm1
          do j=1,yref_
             if ( touch_sp(child) == 0)  &
                  rymc( 1:imc , j            ,1:lmc ,n) = toc1
             if ( touch_np(child) == 0)  &
                  rymc( 1:imc , jm(child)+1-j,1:lmc ,n) = tocm
          end do

       end do  !n
       nullify(rymc)


       do n=1,nt
          do l=1,lmc
             lp = lbeg(child) + (l-1)/zref_
             do i=1,imc
                ip = mod(ibeg(child)-1 + (i-1)/xref_,im(region)) + 1
                toc0(i,l) = rzm(ip,jbeg(child)-1,lp,n)*xzref
                toc1(i,l) = rzm(ip,jbeg(child)  ,lp,n)*xzref   !note the xzref...
                tocm(i,l) = rzm(ip,jend(child)  ,lp,n)*xzref
                tocm1(i,l) =rzm(ip,jend(child)+1,lp,n)*xzref
             end do
          end do

          if ( touch_sp(child) == 0)  rzmc(1:imc,0          ,1:lmc,n) = toc0
          if ( touch_np(child) == 0)  rzmc(1:imc,jm(child)+1,1:lmc,n) = tocm1
          do j=1,yref_
             if ( touch_sp(child) == 0)  &
                  rzmc( 1:imc , j            ,1:lmc ,n) = toc1
             if ( touch_np(child) == 0)  &
                  rzmc( 1:imc , jm(child)+1-j,1:lmc ,n) = tocm
          end do

       end do  !n

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

     end do   ! region tree

    call GO_Timer_End( itim_put_yedges, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine put_yedges


end module advecty
