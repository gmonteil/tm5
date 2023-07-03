!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module advectz

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------------

  private

  public :: AdvectZ_Init, AdvectZ_Done
  public :: advectzzoom


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'AdvectZ'


  ! --- local ------------------------------------

  integer    ::  itim_dynamw


contains


  ! ====================================================================


  subroutine AdvectZ_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/AdvectZ_Init'

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_dynamw, 'advectz dynam', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine AdvectZ_Init


  ! ***


  subroutine AdvectZ_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/AdvectZ_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine AdvectZ_Done


  ! ***


  subroutine advectzzoom( region, status )

    use dims,        only : im, jm, lm, xref, yref, zref, tref
    use dims,        only : parent, zoom2D, touch_sp, touch_np, xcyc
    use dims,        only : splitorderzoom, n_operators, rstatus => status, nsplitsteps
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_const,   only : my_real,com_trac,ierr,mpi_sum
    use mpi_comm,    only : barrier_t
#endif
    use partools  , only : myid, root_t, ntracetloc

    ! --- in/out ---------------------------------

    integer,intent(in)        ::  region
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/advectzzoom'

    ! --- local ----------------------------------

    integer            :: is,ie,js,je,n,q
    integer            :: imr,jmr,lmr,tref_,xref_,yref_,zref_
    logical            :: z_encountered
    character(len=1)   :: dir
    real               :: sum_old, sum_new,sum_old_all,sum_new_all

    ! --- begin ----------------------------------

    !WP! only PE's with nonzero lmloc proceed
    if ( ntracetloc == 0 ) return

    ! write cells along the z-walls of each child
    ! to rm(0,...),rm(imr+1,...) and to interface cells of the child
    !
    ! this option should be there: zooming in z disabled....
    if ( .not. zoom2D ) then
       call escape_tm('advectzzoom: ERROR zoom2D should be true:'// &
            ' zooming along z not allowed')
       !CMKCALL put_zedges(region,rm,rxm,rym,rzm,m,  &
       !CMK                        rmgl,rxmgl,rymgl,rzmgl,mgl)
    end if

    tref_ = tref(region)/tref(parent(region))
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region);  jmr = jm(region);  lmr = lm(region)

    ! determine the scope for advectz:

    if ( region == 1 ) then
       xref_ = 0; yref_ = 0  ! to have is/ie and js/je properly computed
    end if
    ! find q - the place in the splitorderzoom
    q=rstatus(region)/((nsplitsteps/2)*tref_)
    ! corresponding to the begining of the last
    ! triple x-y-z of the parent
    z_encountered=.false.

    do n=1,n_operators          ! now track following steps from q

       dir=splitorderzoom(region,q*(nsplitsteps/2)*tref_+n)

       select case(dir)
       case('x')
          if ( (.not.z_encountered) .or. (xcyc(region) == 1) ) then
             is=1                          ! x-substep is before z =>
             ie=im(region)                 ! full i-scope
          else
             is=xref_+1                    ! z-substep is before x =>
             ie=im(region)-xref_           ! restricted i-scope
          end if
       case('y')
          if (.not.z_encountered) then
             js=1                          ! y-substep is before z =>
             je=jm(region)                 ! full j-scope
          else
             if(touch_sp(region) == 1 ) then
                js=1                       ! sp is southern edge--->full scope.
             else
                js=yref_+1                 ! x-substep is before y =>
             end if
             if ( touch_np(region) == 1 ) then
                je=jm(region)              ! np is northern edge--->full scope
             else
                je=jm(region)-yref_        ! restricted j-scope
             end if
          end if
       case('z')
          z_encountered=.true.
       case ('c')
       case ('v')
       case ('d')
       case ('s')
       case default
          print *,'advectzzoom: strange value in splitorderzoom(',region,',',  &
               q*(nsplitsteps/2)*tref_+n,'): ',q
          call escape_tm('advectzzoom: Error in advectz')
       end select

    end do

    call dynamw( region, is,ie,js,je, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectzzoom


  !-----------------------------------------------------------------------
  !
  !****   dynamw          - vertical tracer transport     v 9.1
  !
  ! programmed by         mh      mpi HH          23-feb-1995
  !
  !       purpose
  !       -------
  !       calculate amount of tracer moved in a vertical advection
  !       substep
  !
  !       interface
  !       ---------
  !       call dynamw
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
  !       changed order of loops for increased performance
  !                               mh, 11-jun-1994
  !     included code for limits of slopes to prevent negative tracer
  !     masses                   mh, 20-jun-1994
  !
  !     zoom version written by mike botchev, march-june 1999
  !
  !-----------------------------------------------------------------------!

  subroutine dynamw( region, is,ie,js,je, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : im, jm, lm
    use dims,        only : xref, yref
    use dims,        only : okdebug
    use dims,        only : nregions
    use dims,        only : parent
    use dims,        only : limits
    use global_data, only : wind_dat
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use partools   , only : myid, root_t
    use partools   , only : ntracetloc
#ifdef with_budgets
    use budget_global, only : budget_flux, lflux1, lflux2
    use budget_global, only : apply_budget_global_flux => apply_budget_global
    use partools     , only : offsetn
    use chem_param   , only : ra
#endif
    use zoom_tools,  only : mix_edges
    use toolbox,     only : escape_tm

    ! --- in/out ---------------------------------

    integer,intent(in) :: region
    integer,intent(in) :: is
    integer,intent(in) :: ie
    integer,intent(in) :: js
    integer,intent(in) :: je
    integer,intent(out) :: status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dynamw'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer   :: rm,rxm,rym,rzm
    real,dimension(:,:,:),pointer     :: m,cm
    real, allocatable                 :: f(:,:)
    integer                           :: i,j,l
    integer                           :: n
    integer                           :: imr,jmr,lmr
    integer                           :: xref_,yref_
    real                              :: max_one
    integer                           :: nglob

    ! --- start ----------------------------------

#ifdef with_limits
    if ( .not. limits ) stop 'macro "with_limits" defined but variable "limits" not'
#else
    if ( limits ) stop 'macro "with_limits" not defined but variable "limits" is'
#endif

    imr=im(region) ; jmr=jm(region) ; lmr=lm(region)

    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    m => m_dat(region)%data
    cm => wind_dat(region)%cm_t

    if ( okdebug ) print *,'dynamw: region=',region,' is,ie,js,je=',is,ie,js,je
    if ( ( region < 0 ) .or. ( region > nregions ) ) &
         call escape_tm( 'dynamw: STOP, illegal number of region !!!')

    ! compute refinement factors with respect to the parent
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))

    ! check is,ie,js,je:
    if ( ( is /= xref_+1 ) .and. ( is /= 1 ) ) &
         call escape_tm('dynamw: Wrong value for IS in dynamw')
    if ( ( ie /= imr-xref_ ) .and. ( ie /= imr ) ) &
         call escape_tm('dynamw: Wrong value for IE in dynamw')
    if ( ( js /= yref_+1 ) .and. ( js /= 1 ) ) &
         call escape_tm('dynamw: Wrong value for JS in dynamw')
    if ( ( je /= jmr-yref_ ) .and. ( je /= jmr ) ) &
         call escape_tm('dynamw: Wrong value for JE in dynamw')

    ! check ...
    if ( any(cm(:,:,0) /= 0.0) .or. any(cm(:,:,lmr) /= 0.0) ) then
      write (gol,'("found non-zero halo interfaces in cm")'); call goErr
      TRACEBACK; stop
    end if

    call GO_Timer_Start( itim_dynamw, status )
    IF_NOTOK_RETURN(status=1)

    ! calculate new air mass distribution

    !$OMP PARALLEL &
    !$OMP  default (none) &
#if defined (with_budgets)
    !$OMP  shared  (apply_budget_global_flux, budget_flux, lflux1, lflux2,ra) &
    !$OMP  shared  (offsetn) &
    !$OMP  shared  (region) &
    !$OMP  private (n,nglob) &
#endif
    !$OMP  shared  (is, ie, js, je, lmr ) &
    !$OMP  shared  (ntracetloc) &
    !$OMP  shared  (rm, rxm, rym, rzm) &
    !$OMP  shared  (m, cm) &
    !$OMP  private (i, j) &
    !$OMP  private (f)

    ! local storage:
    allocate(  f(0:lmr,1:ntracetloc) ) ;  f = 0.0

    !$OMP DO
    do j = js, je
      do i = is, ie

        ! column advection;
        ! update m and rm/rxm/rym/rzrm ;
        ! flux f could  be used for budgets:
        call dynamw_1d( lmr, cm(i,j,0:lmr), m(i,j,1:lmr), &
                         ntracetloc, f, &
                         rm(i,j,1:lmr,:), &
                         rxm(i,j,1:lmr,:), rym(i,j,1:lmr,:), rzm(i,j,1:lmr,:) )

#ifdef with_budgets
        ! calculate flux flowing in from below:
        if ( apply_budget_global_flux ) then
          ! loop over tracers
          do n = 1, ntracetloc
             ! global tracer index:
             nglob = offsetn + n
             ! add flux to budget:
             budget_flux(region)%flux_z1(i,j,nglob) = &
                 budget_flux(region)%flux_z1(i,j,nglob) + f(lflux1(region)-1,n)*1e3/ra(nglob)  ! moles
             budget_flux(region)%flux_z2(i,j,nglob) = &
                 budget_flux(region)%flux_z2(i,j,nglob) + f(lflux2(region)-1,n)*1e3/ra(nglob)  ! moles
          end do ! loop over tracers
        end if
#endif

      end do  ! i
    end do  ! j
    !$OMP END DO

    ! clear:
    deallocate(f)

    !$OMP END PARALLEL

    if ( okdebug .and. myid == root_t ) then
       max_one = maxval(rm(1:imr,1:jmr,1:lmr,1)/m(1:imr,1:jmr,1:lmr))
       print *, 'dynamw: z Maximum value mixing ratio',max_one
    end if

    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(m)
    nullify(cm)

    call GO_Timer_End( itim_dynamw, status )
    IF_NOTOK_RETURN(status=1)

    call mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

  end subroutine dynamw


  ! ***


  pure subroutine dynamw_1d( lmr, cm, m, &
                        ntr, f, rm, rxm, rym, rzm )

    use dims         , only : zero, one

    ! --- in/out ---------------------------------

    integer, intent(in)   ::  lmr
    real, intent(in)      ::  cm(0:lmr)
    real, intent(inout)   ::  m(1:lmr)
    integer, intent(in)   ::  ntr
    real, intent(out)     ::  f(0:lmr,1:ntr)   ! flux through boundaries;
                                               ! used for budgets
    real, intent(inout)   ::  rm(1:lmr,1:ntr)
    real, intent(inout)   ::  rxm(1:lmr,1:ntr)
    real, intent(inout)   ::  rym(1:lmr,1:ntr)
    real, intent(inout)   ::  rzm(1:lmr,1:ntr)

    ! --- local ----------------------------------

    integer              ::  l, n
    real                 ::  mnew(1:lmr)
    real                 ::  pf(0:lmr)
    real                 ::  fx(0:lmr)
    real                 ::  fy(0:lmr)
    real                 ::  gamma

    ! --- begin ----------------------------------

    ! loop over layers; cm at bottom (0) and top (lmr) are zero:
    do l = 1, lmr
      ! mid: two-sided update of the mass
      mnew(l) = m(l) + cm(l-1) - cm(l)
    end do

    ! if requested limit vertical slopes such that no non-negative
    ! tracer masses should occur
#ifdef with_limits
     rzm = max( min( rzm, rm ), -rm )
#endif

    ! loop over tracers
    do n = 1, ntr

      ! loop over interfaces:
      do l = 0, lmr
        ! compute f, pf, fx and fy
        if ( cm(l) == 0.0 ) then
          f(l,n) = 0.0
          pf(l)  = 0.0
          fx(l)  = 0.0
          fy(l)  = 0.0
        else if ( cm(l) > zero ) then
          gamma  = cm(l)/m(l)
          f(l,n) = gamma*(rm(l,n)+(one-gamma)*rzm(l,n))
          pf(l)  = cm(l)*(gamma*gamma*rzm(l,n)-3.*f(l,n))
          fx(l)  = gamma*rxm(l,n)
          fy(l)  = gamma*rym(l,n)
        else
          gamma  = cm(l)/m(l+1)
          f(l,n) = gamma*(rm(l+1,n)-(one+gamma)*rzm(l+1,n))
          pf(l)  = cm(l)*(gamma*gamma*rzm(l+1,n)-3.*f(l,n))
          fx(l)  = gamma*rxm(l+1,n)
          fy(l)  = gamma*rym(l+1,n)
        end if
      end do  ! l

      ! calculate new tracer mass, and tracer mass slopes
      ! update rm, rzm, rxm and rym in interior layers of the column

      do l = 1, lmr
        rm(l,n) = rm(l,n) + f(l-1,n)-f(l,n)
        rzm(l,n) = rzm(l,n) + &
            ( pf(l-1)-pf(l) &
              - (cm(l-1)-cm(l))* rzm(l,n)               &
              + 3.0*((cm(l-1)+cm(l))* rm(l,n)            &
              - (f(l-1,n)+f(l,n))* m(l) ) ) / mnew(l)
#ifdef with_limits
        rzm(l,n) = max(min(rzm(l,n),rm(l,n)),-rm(l,n))
#endif
        rxm(l,n) = rxm(l,n) + (fx(l-1)-fx(l))
        rym(l,n) = rym(l,n) + (fy(l-1)-fy(l))
      end do

    end do ! loop over tracers

    ! store new air mass in m array
    m = mnew

  end subroutine dynamw_1d


end module advectz
