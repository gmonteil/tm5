!###############################################################################
!
! module to check for CFL
!
! Exactly the same operations are performed on m as
! in the 'real' advection routines.
! However, rm, rxm, rym, rzm are not updated
! The 'old' m is saved using store_masses and
! restored with 'restore_masses'
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module advectm_cfl

  use GO,    only: gol, goErr, goPr
  use dims,  only: nregions, maxref, zoom_mode
  use dims,  only: revert
  use dims,  only: okdebug

  ! ------------ interface ---------------------
  private

  ! public routines:
  public :: check_cfl, init_cfl, done_cfl
  public :: Setup_MassFlow
  public :: advectx_get_nloop

  ! --- const --------------------------------

  character(len=*), parameter  ::  mname = 'advectm_cfl'

  integer, parameter             :: max_global_iteration=16
  integer                        :: global_iteration
  integer, save                  :: ndyn_save
  integer, parameter             :: n_operators = 3  ! xyz
  integer, parameter             :: nsplitsteps  = 6  ! xyzzyx
  !logical                        :: okdebug = .true.
  integer                        :: cfl_outputstep
  integer,dimension(nregions)    :: regionm_status
  character,dimension(nregions,maxref*nsplitsteps) :: splitorderzoom
  character,dimension(nsplitsteps),parameter       :: splitorder = &
       (/'x','y','z','z','y','x'/)

  ! local CFL counters:
  ! o Courant number: for each interface:
  !      xi = (mass flux [kg/time]) / (mass in source cell [kg])
  !   thus fraction [1/time] of mass leaving the adjacent cell
  !   throught the interface.
  !   Array 'xim' contains the maximum per region and x/y/z direction:
  real                             ::  xim(nregions,3) = 0.0
  ! o number of cells with |xi|>1 ?
  integer                          ::  nxim(nregions,3) = 0
  ! o maximum number of cfl loops per region and direction:
  integer                          ::  mloop_max(nregions,3) = 0

  type oldmass_data
    real,dimension(:,:,:),  pointer   :: m
    real,dimension(:,:,:),  pointer   :: am
    real,dimension(:,:,:),  pointer   :: bm
    real,dimension(:,:,:),  pointer   :: cm
  end type oldmass_data
  type(oldmass_data), dimension(nregions)   :: oldmass_dat

  ! timing:
  integer                             ::  itim_check_cfl
  integer                             ::  itim_advectmxzoom
  integer                             ::  itim_advectmyzoom
  integer                             ::  itim_advectmzzoom
  integer                             ::  itim_advectmx
  integer                             ::  itim_dynamum
  integer                             ::  itim_dynamvm
  integer                             ::  itim_dynamwm
  integer                             ::  itim_mixm_edges
  integer                             ::  itim_updatem_parent


contains


  subroutine init_cfl( status )

    use GO         , only : GO_Timer_Def
    use GO         , only : ReadRc
    use global_data, only : rcF
    use dims       , only : nregions, im, jm, lm

    ! --- in/out ---------------------------------

    integer, intent(out)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/init_cfl'

    ! --- local ----------------------------------

    integer           ::  region, imr, jmr, lmr
    logical           ::  found

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_check_cfl, 'check_cfl', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_advectmxzoom, 'advectmxzoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_advectmyzoom, 'advectmyzoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_advectmzzoom, 'advectmzzoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_advectmx, 'advectmx', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynamum, 'dynamum', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynamvm, 'dynamvm', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynamwm, 'dynamwm', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_mixm_edges, 'mixm_edges', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_updatem_parent, 'updatem_parent', status )
    IF_NOTOK_RETURN(status=1)

    ! settings:
    call ReadRc( rcF, 'cfl.outputstep', cfl_outputstep, status)
    IF_NOTOK_RETURN(status=1)

    ! temporary storage:
    do region = 1, nregions
       imr = im(region)
       jmr = jm(region)
       lmr = lm(region)
       allocate (oldmass_dat(region)%m(-1:imr+2, -1:jmr+2, lmr))
       allocate (oldmass_dat(region)%am(0:imr+1, 0:jmr+1, 0:lmr+1))
       allocate (oldmass_dat(region)%bm(0:imr+1, 0:jmr+1, 0:lmr+1))
       allocate (oldmass_dat(region)%cm(0:imr+1, 0:jmr+1, 0:lmr+1))
    end do

    ! something to do with operator splittting:
    call define_splitorderzoom

    ! ok
    status = 0

  end subroutine init_cfl


  ! ***


  subroutine Check_CFL( t1, t2, n, status )

    use GO         , only : gol, goErr, goPr
    use GO         , only : TDate, IncrDate, operator(+), operator(-), rTotal, wrtgol
    use global_data, only : wind_dat
    use toolbox,     only : escape_tm
    use dims,        only : ndyn_max
    use dims,        only : revert
    use datetime,    only : new_valid_timestep
    use partools   , only : myid, root

    ! --- in/out -------------------------------

    type(TDate), intent(in)  ::  t1, t2
    integer, intent(out)     ::  n
    integer, intent(out)     ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Check_CFL'

    ! --- local -------------------------------

    type(TDate)              ::  tr(2)
    integer                  ::  i
    integer                  ::  ndyn, ndyn_old
    integer                  ::  region
    real                     ::  fraction
    logical                  ::  cfl_ok

    ! --- begin -------------------------------

    ! init with allowed maximum:
    ndyn = ndyn_max

    ! n is the number of dynamic intervals within the
    ! time interval for which the meteo has been setup;
    ! first guess from allowed maximum:
    n = ceiling( abs(rTotal(t2-t1,'sec')) / real(ndyn) )

    ! info ...
!    call wrtgol( rname//': time range ', t1, ' - ', t2 ); call goPr
    write (gol,'(a,": initial timestep of  ",i5," sec")') rname, ndyn; call goPr
!    write (gol,'(a,": initial nr. of steps ",i0)') rname, n; call goPr

    ! check ..
    if ( n < 1 ) then
      write (gol,'("unrealistic number of initial steps: ",i0)') n; call goErr
      TRACEBACK; status=1; return
    end if

    ! increase number of time steps until cfl ok or maximum exceeded
    global_iteration = 0
    do

      ! increase number of time steps:
      global_iteration = global_iteration + 1

      ! info ...
!      write (gol,'(a,": attempt ",i0," ...")') rname, global_iteration; call goPr

      ! check ...
      if ( global_iteration > max_global_iteration ) then
        write (gol,'("exceeded maximum number of time steps")'); call goErr
        TRACEBACK; status=1; return
      end if

      ! loop over number of ndyn timesteps;
      ! apply advection over complete (large) time interval
      do i = 1, n

        ! small time intervals in positive direction;
        ! for reverse run, the 'store_masses' will set the initial mass
        ! to the absolute begin of the interval:
        if ( revert > 0 ) then
          tr(1) = t1 + IncrDate(sec=(i-1)*ndyn)
          tr(2) = t1 + IncrDate(sec= i   *ndyn)
        else
          tr(1) = t2 + IncrDate(sec=(i-1)*ndyn)
          tr(2) = t2 + IncrDate(sec= i   *ndyn)
        end if

        ! info ...
!        write (gol,'(a,":   timestep ",i0," / ",i0)') rname, i, n; call goPr
!        call wrtgol( rname//':     tr ', tr(1), ' - ', tr(2) ); call goPr
        !write (gol,'(a,":     ndyn ",i0," sec")') rname, ndyn; call goPr

        ! fill am/bm/cm with balanced mass flows, eventually time interpolated
        call Setup_MassFlow( tr, ndyn, status )
        IF_ERROR_RETURN(status=1)

        ! initial step?
        if ( i == 1 ) then
          ! save current air masses in 'oldmass_dat' ;
          ! in case of revert run, also advect mass backward in time,
          ! which needs the just computed am/bm/cm :
          call store_masses( n, status )
        IF_NOTOK_RETURN(status=1)
        end if

        ! first half step of advection:
        call determine_cfl_iter( 1, cfl_ok, status )     ! xyz ...
        IF_ERROR_RETURN(status=1)
        ! second half step if first is ok:
        if ( cfl_ok ) then
          call determine_cfl_iter( 1, cfl_ok, status )  ! ... zyx
          IF_ERROR_RETURN(status=1)
        end if

        ! info ...
!        write (gol,'(a,"     cfl_ok : ",l1)') rname, cfl_ok; call goPr

        ! cfl not ok ? then leave loop and decrease time step:
        if ( .not. cfl_ok ) exit

        ! flag ...
        regionm_status(:) = 0

      end do

      ! restore massa at begin of time interval
      call restore_masses( status )
      IF_NOTOK_RETURN(status=1)

      if ( .not. cfl_ok ) then
        ndyn_old = ndyn
        ! new ndyn should be a denominator of e.g. 3 hour
        call new_valid_timestep( ndyn, 3*3600, cfl_outputstep )
        ! info ...
        write (gol,'(a,":   reducing timestep to ",i5," sec")') rname, ndyn; call goPr
        ! update number of dynamic time steps:
        n = nint( real(n*ndyn_old)/real(ndyn) ) ! should work!
        ! update time integrated mass fluxes:
        fraction = real(ndyn)/real(ndyn_old)
        do region = 1, nregions
          wind_dat(region)%am_t = wind_dat(region)%am_t*fraction
          wind_dat(region)%bm_t = wind_dat(region)%bm_t*fraction
          wind_dat(region)%cm_t = wind_dat(region)%cm_t*fraction
        enddo
      else
        exit
      end if

    end do   ! global_iteration

    ! ok
    status = 0

  end subroutine Check_CFL


  ! ***


  ! Fill am/bm/cm :
  !  o interpolate balanced mass fluxes pu/pv in time (if necessary)
  !  o multiply with dt; store in am/bm/cm
  !  o handle edges of zoom regions

  subroutine Setup_MassFlow( tr, ndyn, status )

    use GO          , only : TDate, operator(<)
    !use GO          , only : PrintDate2
    use MeteoData   , only : pu_dat, pv_dat
    use MeteoData   , only : TimeInterpolation
    use zoom_tools  , only : coarsen_region
    use advect_tools, only : coarsen_ambm, dynam0

    ! --- in/out -----------------------------------

    type(TDate), intent(in)    ::  tr(2)
    integer, intent(in)        ::  ndyn
    integer, intent(out)       ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Setup_MassFlow'

    ! --- local ------------------------------------

    integer                 ::  n

    ! --- begin ----------------------------------

    !call PrintDate2( trim(rname)//': ', tr(1), ' - ', tr(2) )

    !
    ! time interpolation
    !

    ! loop over regions:
    do n = 1, nregions
      !
      call TimeInterpolation( pu_dat(n), tr, status )
      IF_ERROR_RETURN(status=1)
      !
      call TimeInterpolation( pv_dat(n), tr, status )
      IF_ERROR_RETURN(status=1)
      !
    end do

    !
    ! fill am/bm/cm
    !

    ! edges etc
    do n = 1, nregions
      ! pu/pv to am/bm/cm
      call dynam0( n, ndyn, status )
    IF_NOTOK_RETURN(status=1)
      ! update coarser parent grids from fine child grids (am,bm,cm,m)
      call coarsen_region( n )
    end do

    ! copy coarse mass-fluxes from the coarse to the fine region
    do n = 1, nregions
      call coarsen_ambm( n )
    end do

    !
    ! done
    !

    ! ok
    status = 0

  end subroutine Setup_MassFlow


  ! ****************************************************************
  ! ***
  ! *** mass data
  ! ***
  ! ****************************************************************


  subroutine store_masses( n, status )
    use dims,        only: nregions, tref, im, jm, lm
    use global_data, only: mass_dat, wind_dat
    use MeteoData  , only : m_dat
    use parTools,    only   : myid, root

    ! I/O
    integer, intent(in)  :: n    ! number of ndym timesteps...
    integer, intent(out) :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/store_masses'
    ! local
    real,dimension(:,:,:),pointer      :: m,am,bm,cm
    integer   :: region
    integer   :: imr
    integer   :: jmr
    integer   :: lmr
    integer   :: ntimes
    ! begin
    do region = 1, nregions
      ! store:
      oldmass_dat(region)%m = m_dat(region)%data
       oldmass_dat(region)%am = wind_dat(region)%am_t
       oldmass_dat(region)%bm = wind_dat(region)%bm_t
       oldmass_dat(region)%cm = wind_dat(region)%cm_t
      ! set flag:
       regionm_status(region) = 0
      ! reverse run?
       if (revert == -1) then
        ! set pointers:
        m => m_dat(region)%data
          am => wind_dat(region)%am_t
          bm => wind_dat(region)%bm_t
          cm => wind_dat(region)%cm_t
        ! local dims:
          imr = im(region)
          jmr = jm(region)
          lmr = lm(region)

        ! backward advection (two times since xyz..zyx) ;
        ! number of time steps; fluxes valid for ndyn/2/tref(region) sec
        ntimes = n * 2*tref(region)

        ! apply:
          m(1:imr,1:jmr,1:lmr) =  m(1:imr,1:jmr,1:lmr) &
                        - ntimes*am(0:imr-1,1:jmr  ,1:lmr  ) &   !east
                        + ntimes*am(1:imr  ,1:jmr  ,1:lmr  ) &   !west
                        - ntimes*bm(1:imr  ,1:jmr  ,1:lmr  ) &   !south
                        + ntimes*bm(1:imr  ,2:jmr+1,1:lmr  ) &   !north
                        - ntimes*cm(1:imr  ,1:jmr  ,0:lmr-1) &   !lower
                        + ntimes*cm(1:imr  ,1:jmr  ,1:lmr  )     !upper

        ! clear:
          nullify(am)
          nullify(bm)
          nullify(cm)
          nullify(m)

      end if  ! reverse run

    end do ! retions

    ! ok
    status = 0

  end subroutine store_masses

  ! *

  subroutine restore_masses( status )

    use dims,        only: nregions
    use global_data, only: mass_dat, wind_dat
    use MeteoData  , only : m_dat

    ! I/O
    integer, intent(out) :: status

    ! local:
    integer :: region

    ! begin

    ! loop over regions:
    do region = 1, nregions
      ! copy:
       m_dat(region)%data   = oldmass_dat(region)%m
       wind_dat(region)%am_t = oldmass_dat(region)%am
       wind_dat(region)%bm_t = oldmass_dat(region)%bm
       wind_dat(region)%cm_t = oldmass_dat(region)%cm
    enddo

    ! ok
    status = 0

  end subroutine restore_masses


  ! subroutine to determine a 'global' iteration
  ! that is needed to prevent CFL violations.
  ! Normally the operator splitted sequence is:
  ! xyz        vvzyx
  !    xyzvvzyx     vzyxxyzv
  ! If anywhere will occur a CFL violation, the global
  ! timestep is reduced and the whole sequence is called
  ! two (or more) times.
  ! The advantage is that 'masses' are largely restored by
  ! advection from neighbouring cells

  recursive subroutine determine_cfl_iter( region, cfl_ok, status )

    use dims,   only : parent, tref, revert
    use dims,   only : children
    use global_data, only: mass_dat

    implicit none

    ! input/output
    integer, intent(in)                 ::  region
    logical, intent(out)                ::  cfl_ok
    integer, intent(out)                ::  status

    ! const
    character(len=*), parameter  ::  rname = mname//'/determine_cfl_iter'

    ! local
    integer                            :: child, i, ichild, tref_, n_children
    integer, dimension(3)              :: ll,uu

    ! determine refinement factor with respect to the parent
    tref_ = tref(region)/tref(parent(region))

    do i=1,tref_

       call do_steps( region, cfl_ok, status )
       IF_NOTOK_RETURN(status=1)
       if ( .not. cfl_ok ) then; status=0; return; end if

       ! CALL advect_region for all the children (IF there are any)
       ichild = 0
       do while ( ichild < children(region,0) )
          ichild = ichild + 1
          child = children(region,ichild)
          call determine_cfl_iter( child, cfl_ok, status )
          IF_NOTOK_RETURN(status=1)
          if ( .not. cfl_ok ) then; status=0; return; end if
       end do
       !do the remaining steps if necessary...!
       if (mod(regionm_status(region),n_operators) /= 0 ) then
          call do_steps( region, cfl_ok, status )
          IF_NOTOK_RETURN(status=1)
          if ( .not. cfl_ok ) then; status=0; return; end if
       end if
    end do

    if ( region /= 1 ) then
      call updatem_parent( region, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine determine_cfl_iter


  ! ***


  subroutine do_steps( region, cfl_ok, status )

    use toolbox,        only : escape_tm

    ! input/output
    integer,intent(in)                  :: region
    logical, intent(out)                :: cfl_ok
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/advectxzoom'

    ! local
    integer          :: child, i123, ichild, tref_child, reg, rgi, j
    character        :: tobedone
    character(len=1) :: next_step, prev_step

    if ( okdebug ) print *,'do_steps: region ',region

    do i123=1,n_operators

       next_step = splitorderzoom(region,regionm_status(region)+1)
       if ( regionm_status(region) /= 0 ) then
          prev_step = splitorderzoom(region,regionm_status(region))
       else
          prev_step = ' '
       end if

       if (okdebug) then
             ! it's only to make work of the code visible -- step to be done
             ! will be printed with capital letter (X,Y or Z)
             do reg=1,region
                tobedone = ' '
                if ( reg == region ) tobedone = upper(next_step)
                print *, 'do_steps: ',reg,': ', &
                     splitorderzoom(reg,1:regionm_status(reg)),tobedone
             enddo
       endif
       tobedone = upper(next_step)
       select case(next_step)
       case('x')
          call advectmxzoom( region, cfl_ok, status )
          IF_NOTOK_RETURN(status=1)
       case('y')
          call advectmyzoom( region, cfl_ok, status )
          IF_NOTOK_RETURN(status=1)
       case('z')
          call advectmzzoom( region, cfl_ok, status )
          IF_NOTOK_RETURN(status=1)
       case default
          print *,'do_steps:  strange value in splitorderzoom: ',  &
               splitorderzoom(region,regionm_status(region))
          print *,'do_steps:  (must be x, y or z)'
          call escape_tm('do_steps: Error')
       end select

       if ( .not. cfl_ok ) then; status=0; return; end if

       regionm_status(region) = regionm_status(region)+1
           ! these statements are needed when
           ! more processes are involved that change rm
           ! for instance emissions.
           ! for m-advection only just leave when ready
       if ( mod(regionm_status(region),n_operators) == 0 ) then
          exit ! e.g.after zyx or xyz
       end if
    end do

    ! ok
    status = 0

  end subroutine do_steps

  character function upper(xyz)

    implicit none

    character(1),intent(in) :: xyz

    if (xyz=='x') then
       upper = 'X'
    else if (xyz=='y') then
       upper = 'Y'
    else if (xyz=='z') then
       upper = 'Z'
    else
       upper = '_'
    end if

  end function upper

 subroutine define_splitorderzoom
    !
    ! splitorderzoom is splitorder specified for each region
    ! the subroutine is normally called once per run of TM5
    ! ATTN: subroutine was written based on assumption that a region can not
    ! have number less than number of its parent: region > parent(region)
    ! written by patrick berkvens and mike botchev, march-june 1999
    !
    use dims,        only : zoom_mode, tref, parent

    implicit none
    ! local
    integer :: i, j, j0, region, tref_
    !character,dimension(3):: reverse

    if ( okdebug ) print *,'define_splitorderzoom: '
    splitorderzoom = ' '
    splitorderzoom(1,1:nsplitsteps) = splitorder
    if ( zoom_mode == 1 ) then
       ! zoom_mode==1 means:
       !  x y z       | z y x
       !  x y z z y x | z y x x y z
       !  ...
       do region=2,nregions
          tref_ = tref(region)/tref(parent(region))
          if ( (tref_ > 4 ) .or. ( tref_ == 3 ) ) then
             print *,'define_splitorderzoom: ERROR region = ',region
             print *,'define_splitorderzoom: ', &
                  'wrong value for tref(region): ',tref(region)
             print *,'define_splitorderzoom: ', &
                  'should be 1,2 or 4 times bigger than tref of its parent'
             print *,'define_splitorderzoom: ', &
                  'use another value for parameter zoom_mode'
             stop
          end if
          j0 = 1  ! step counter for parent
          j = 1   ! step counter for region
          do while ( splitorderzoom(parent(region),j0) /= ' ' )
             ! use step triple of parent tref_ times:
             do i=1,tref_
                if (mod(i,2)==1) then
                   ! copy step triple of parent:
                   !WP! changed 2 into n_oper
                   splitorderzoom(region,j:j+(n_operators-1)) = &
                        splitorderzoom(parent(region),j0:j0+(n_operators-1))
                else
                   ! step triple of parent in reverse order:
                   !    WP! changed 2 into n_op
                   splitorderzoom(region,j:j+(n_operators-1)) = &
                        reverse(splitorderzoom(parent(region), &
                        j0:j0+(n_operators-1)))
                end if
                !WP! changed 3 into n_ope...
                j = j + (n_operators)
             end do
             ! next step triple of parent  !WP! changed 3 into n_oper..
             j0 = j0 + (n_operators)
          end do
       end do
    else if ( zoom_mode == 2 ) then

       ! zoom_mode==2 means:
       !  x y z       | z y x
       !  x y z x y z | z y x z y x
       !  ...
       do region=2,nregions
          tref_ = tref(region)/tref(parent(region))
          j0 = 1  ! step counter for parent
          j = 1   ! step counter for region
          do while (splitorderzoom(parent(region),j0)/=' ')
             ! replicate step triple of parent tref_ times:
             do i=1,tref_
                splitorderzoom(region,j:j+(n_operators-1)) = &
                     splitorderzoom(parent(region),j0:j0+(n_operators-1))
                j = j + (n_operators)
             end do
             ! next step  of parent
             j0 = j0 + (n_operators)
          end do
       end do
    else
       print *,'define_splitorderzoom: wrong value for zoom_mode ',zoom_mode
       stop
    end if

  contains

    function reverse(str) result(inv_str)
      !WP! changed to n_operators....
      character,dimension(n_operators):: str,inv_str
      integer ::i

      do i=1,(n_operators)
         inv_str(i)=str(n_operators+1-i)
      end do

    end function reverse

  end subroutine define_splitorderzoom


  ! ================================================================================
  ! ===
  ! === advect m x
  ! ===
  ! ================================================================================


  ! set parameters for advectx
  ! written by patrick berkvens and mike botchev, march-june 1999
  ! updated and modified by MK, dec 2002

  subroutine advectmxzoom( region, cfl_ok, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,          only : xref, yref, zref, tref, im, jm, lm
    use dims,          only : zoom2D
    use dims,          only : touch_sp, touch_np, adv_scheme
    use dims,          only : parent, zero
    use global_data  , only : mass_dat, wind_dat
    use toolbox,       only : escape_tm

    ! input/output
    integer,intent(in)                  :: region
    logical, intent(out)                :: cfl_ok
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

    ! start

    call GO_Timer_Start( itim_advectmxzoom, status )
    IF_NOTOK_RETURN(status=1)

    call putm_xedges(region)   ! write BC to cildren
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
    ! determine the scope for advectx:
    if ( region == 1 ) then
       yref_ = 0
       zref_ = 0  ! to have js/je and ls/le properly computed
    end if
    q=regionm_status(region)/((nsplitsteps/2)*tref_)
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
       case default
          print *,'advectmxzoom: strange value in splitorderzoom(',region,',',   &
               q*(nsplitsteps/2)*tref_+n,'): ',q
          call escape_tm( 'advectmxzoom: error in x-advection ')
       end select

    end do  ! n=1,n_operators

    if ( ( mod(regionm_status(region),(nsplitsteps/2)*tref_) >= n_operators ) .and. &
         ( adv_scheme == 'slope' ) .and. &
         ( im(region)/xref(region) < im(1) ) ) then
       ! IF (1) more than fourth substep, (2) slope scheme, and
       ! (3) the region is not [0;360] degrees wide THEN
       ! at this substep no coarse fluxes will be applied

       ! for slopes only: zero out velocities via the x-edges of the region
       ! to assure that no fluxes will be applied

       if (okdebug) print *,'advectmxzoom: zeroing out outer mass fluxes'
       am0 = am(    xref_-1,:,:)
       am(    xref_-1,:,:) = zero
       am1 = am(imr-xref_+1,:,:)
       am(imr-xref_+1,:,:) = zero

    end if

    call advectmx( region, js,je,ls,le, cfl_ok, status )
    IF_NOTOK_RETURN(status=1)

    if ( ( mod(regionm_status(region),(nsplitsteps/2)*tref_) >= n_operators ).and. &
         ( adv_scheme == 'slope' ) .and. &
         ( im(region)/xref(region) < im(1) ) ) then
       am(    xref_-1,:,:) = am0
       am(imr-xref_+1,:,:) = am1
    end if

    deallocate(am0)
    deallocate(am1)
    nullify(am)

    call GO_Timer_End( itim_advectmxzoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectmxzoom


  ! ***


  !
  ! makes reduced grid pre-/postprocessing and switches between
  ! dynamu and dynamu1
  ! written by mike botchev, march-june 1999
  !

  subroutine advectmx( region, js,je,ls,le, cfl_ok, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : adv_scheme, im, jm, lm
    use redgridZoom, only : grid_reduced, nred, uni2red_mf
    use global_data, only : mass_dat, wind_dat
    use MeteoData  , only : m_dat
    use toolbox,     only : escape_tm

    ! input
    integer,intent(in) :: region
    integer,intent(in) :: js
    integer,intent(in) :: je
    integer,intent(in) :: ls
    integer,intent(in) :: le
    logical, intent(out)                ::  cfl_ok
    integer,intent(out)                 ::  status

    ! const
    character(len=*), parameter  ::  rname = mname//'/advectmx'

    ! local
    real,dimension(:,:,:),pointer     :: m
    real,dimension(:,:,:),pointer     :: am
    real,dimension(:,:,:),allocatable :: m_uni   ! for reduced grid...
    real,dimension(:,:,:),allocatable :: am_uni  ! for reduced grid...
    integer              :: imr,jmr,lmr

    ! start

    call GO_Timer_Start( itim_advectmx, status )
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
       ! reduce m (special local routine)
       call uni2red(region)
       ! reduce am:
       call uni2red_mf(region)
    end if

    call dynamum( region, js,je,ls,le, cfl_ok, status )
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

    call GO_Timer_End( itim_advectmx, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectmx


  ! ***


  subroutine uni2red(region)
    !
    ! transforms data from uniform grid to reduced grid
    ! written by mike botchev, march-june 1999
    ! modified by Maarten Krol, dec 2002
    !
    use dims, only       : im,jm,lm
    use redgridZoom, only: nred, clustsize, jred, imredj
    use global_data, only: mass_dat
    use MeteoData  , only : m_dat
    implicit none
    ! input
    integer,intent(in)  :: region
    ! local
    real,dimension(:,:,:),  pointer         :: m
    integer i,ie,is,j,l,lrg,redfact,lmr
    real summ
    ! start
    m => m_dat(region)%data
    lmr=lm(region)
    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  1,imredj(lrg,region)
             ! the is:ie  array section will be reduced to i
             is = (i-1)*redfact + 1
             ie = i*redfact
             summ = sum(m(is:ie,j,l))
             m(i,j,l) = summ
          end do  !i
       end do  !l
    end do   !redgrid...
    nullify(m)
  end subroutine uni2red


  ! ***


  subroutine red2uni(region,m_uni)
    !
    ! transforms data from reduced grid back to uniform grid
    ! written by mike botchev, march-june 1999
    !
    use dims
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use redgridZoom, only : nred, clustsize, jred, imredj
    implicit none
    ! input
    integer,intent(in) :: region
    real,dimension(-1:im(region)+2,-1:jm(region)+2,lm(region)), &
         intent(in) :: m_uni
    ! m_uni: same declaration as m!
    ! local
    real,dimension(:,:,:)  ,pointer    :: m
    integer :: i, ie, ii, is, j, l, lrg, n, redfact, lmr
    real    :: hi, mass, mass_coord, rmm, slope, m_old
    ! start
    lmr=lm(region)
    m => m_dat(region)%data
    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  imredj(lrg,region),1,-1
             ! the i cell will be distributed within the is:ie array section
             is = (i-1)*redfact + 1
             ie = i*redfact
             !m_uni is the mass-distribution in the non-reduced grid/divided by
             !the reduced_grid mass. This is used as distribution function!....
             mass=m(i,j,l); m(is:ie,j,l)= m_uni(is:ie,j,l)
          end do
          ! update cell(0,...) according to the periodic bc's:
          m(0,j,l) = m(im(region),j,l)
       end do
    end do
    nullify(m)

  end subroutine red2uni


  ! ***


  subroutine red2uni_em(region)
    !
    ! transforms data from reduced grid back to uniform grid
    ! written by mike botchev, march-june 1999
    !
    use dims, only       : im,jm,lm
    use redgridZoom, only: nred, clustsize, jred, imredj
    use global_data, only: mass_dat
    use MeteoData  , only : m_dat
    implicit none
    ! input
    integer,intent(in) :: region
    ! local
    real,dimension(:,:,:)  ,pointer    :: m
    integer i,ie,is,j,l,lrg,redfact,lmr
    real mass
    ! start
    lmr=lm(region)
    m => m_dat(region)%data
    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  imredj(lrg,region),1,-1
             is = (i-1)*redfact + 1
             ie = i*redfact
             mass=m(i,j,l); m(is:ie,j,l)= mass/(ie-is+1)
          enddo
          m(0,j,l) = m(im(region),j,l)
       end do
    end do
    nullify(m)
  end subroutine red2uni_em


  ! ***


  subroutine dynamum( region, js,je,ls,le, cfl_ok, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : nregions, parent
    use dims,        only : im, jm, lm, xref, yref, zref, tref
    use dims,        only : zero, one, limits, xcyc
    use redgridZoom, only : grid_reduced, nred, imred
    use global_data, only : mass_dat, wind_dat
    use MeteoData  , only : m_dat
    use toolbox,     only : escape_tm

    ! input/output
    integer,intent(in) :: region
    integer,intent(in) :: js
    integer,intent(in) :: je
    integer,intent(in) :: ls
    integer,intent(in) :: le
    logical, intent(out)                :: cfl_ok
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/dynamum'

    ! local
    real,dimension(:,:,:),  pointer    :: m,am
    real,dimension(:,:,:),allocatable  :: mnew
    integer i,ie,is,j,l,n
    integer imr,jmr,lmr,tref_,xref_,yref_,zref_,ic,nloop, iloop
    real :: max_alpha
    real alpha,x,rmold
    real,dimension(:), allocatable :: mx
    logical                        :: cfl_okl
    integer, parameter             :: max_nloop = 10
    ! openmp varaiables:
    integer   ::  l_mloop_max
    real      ::  l_xim

    call GO_Timer_Start( itim_dynamum, status )
    IF_NOTOK_RETURN(status=1)

    ! prepare x-edges to allow for uniform work:
    call compressm_xedges(region,js,je,ls,le)

    ! compute refinement factors with respect to the parent
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    tref_ = tref(region)/tref(parent(region))

    ! region dimensions:
    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    allocate( mnew(im(region),jm(region),lmr) )

    am => wind_dat(region)%am_t
    m => m_dat(region)%data

    ! compute is/ie -- cells is:ie,js:je:ls:le will be updated
    is = xref_; ie = imr-xref_+1

    if ((region)/xref(region)==im(1)) then
       ! periodic boundary condition
       is = 1; ie = imr
    end if

    ! init recuced variables:
    l_mloop_max = 0
    l_xim       = 0.0
    cfl_okl     = .true.

    ! loop over vertical layers:
    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( js, je, ls, le ) &
    !$OMP   shared ( xcyc ) &
    !$OMP   shared ( imr, nred, imred ) &
    !$OMP   reduction(.and.:cfl_okl) &
    !$OMP   reduction(max:l_mloop_max) &
    !$OMP   reduction(max:l_xim) &
    !$OMP   shared ( m, mnew ) &
    !$OMP   shared ( am ) &
    !$OMP   private ( is, ie ) &
    !$OMP   private ( i, j, l ) &
    !$OMP   private ( nloop ) &
    !$OMP   private ( alpha, max_alpha ) &
    !$OMP   private ( mx )
    !$OMP   DO
    do l = ls, le
      ! loop over latitudes:
      do j = js, je
        ! default range in lon direction:
        is = 1; ie = imr
        ! reduced grid: less points will be handled
        if (grid_reduced.and.nred(region).ne.0) ie = imred(j,region)

         ! Determine number of loops required to avoid CFLs...
         nloop = 1   ! default run the loop one time!
         alpha = 2.0
         allocate(mx(is-1:ie+1))
         !CMK added oct2003: nloop limit..
         max_alpha = 0.0
         do
           ! no cfl violations ? then leave:
           if ( abs(alpha) < one ) exit
           ! reached maximum ? then leave:
           if ( nloop >= max_nloop ) exit

           ! copy mass to temp array
           mx(is-1:ie+1) = m(is-1:ie+1,j,l)
           if(xcyc(region)==1) then
             mx(0) = mx(ie)
             mx(ie+1) = mx(1)
           end if

           xloop: do iloop = 1, nloop
             ! According to xlf error messages:
             !   "The reduction variable cfl_ok must be present
             !   on the right hand side of the reduction statement."
             cfl_okl = cfl_okl .and. .true.
             ! loop over longitude cells:
             do i=is-1,ie
               if (am(i,j,l)>=zero)  then
                 alpha=am(i,j,l)/mx(i)
               else
                 alpha=am(i,j,l)/mx(i+1)
               end if
               if((abs(alpha)>=one)) then                   !PB
                 cfl_okl = cfl_okl .and. .false.
                 exit xloop
               end if
               max_alpha = max(max_alpha, abs(alpha))
             end do  ! i
             ! update mass
             ! note: only is:ie updated, is-1, ie+1 are BC
             mx(is:ie) = mx(is:ie) + am(is-1:ie-1,j,l) - am(is:ie,j,l)
             if ( xcyc(region) == 1 ) then   !except for periodic BC...
               mx(0) = mx(ie)
               mx(ie+1) = mx(1)
             end if
           end do  xloop

           if(.not.cfl_okl) then
             ! reduce mass flux
             am(is-1:ie,j,l) = am(is-1:ie,j,l)*nloop/(nloop+1)
             !if ( okdebug ) then
             !  print *, 'dynamu: ',j,l, 'nloop increased to', nloop + 1, alpha
             !end if
             nloop = nloop + 1
             max_alpha = 0.0
             if (nloop == max_nloop) then
               ! not good solution found: set flag and return
               ! note that am is restored in the calling routines
               ! and m will be restored also!
               ! According to xlf error messages:
               !   "The reduction variable cfl_ok must be present
               !   on the right hand side of the reduction statement."
               cfl_okl = cfl_okl .and. .false.
               exit
             endif
          end if

        end do !while alpha>1

        ! clear:
        deallocate(mx)
        ! ok ?
        if ( cfl_okl ) then
          ! update loop variables:
          l_mloop_max = max( l_mloop_max, nloop )
          l_xim = max( l_xim, max_alpha )
          do iloop = 1,nloop   ! CFL loop
             ! calculate new air mass distribution
             mnew(is:ie,j,l)=m(is:ie,j,l) + am(is-1:ie-1,j,l)-am(is:ie,j,l)
             if ( xcyc(region) == 1 ) then
                m  (0   ,j,l)   = m  (ie,j,l) !periodic BCs...
                m  (ie+1,j,l)   = m  ( 1,j,l)
             end if
             m(is:ie,j,l)=mnew(is:ie,j,l)
          end do  ! iloop
          ! restore 'old' am
          am(is-1:ie,j,l) = am(is-1:ie,j,l)*nloop
        end if

      end do  ! l-loop
    end do   ! j-loop
    !$OMP   END DO
    !$OMP END PARALLEL

    ! update reduced variables:
    mloop_max(region,1) = max( mloop_max(region,1), l_mloop_max )
    xim      (region,1) = max( xim      (region,1), l_xim       )
    cfl_ok              = cfl_okl

    ! restore some arrays:
    if ( cfl_ok ) then
      if ( xcyc(region) == 1 ) then
        ! periodic boundary condition
        do j=js,je
          ! reduced grid
          if (grid_reduced.and.nred(region).ne.0) imr = imred(j,region)
          m  (0    ,j,:)   = m  (imr,j,:)
          m  (imr+1,j,:)   = m  (1,  j,:)
        end do
      end if
    end if

    ! clear:
    deallocate(mnew)
    nullify(am)
    nullify(m)

    ! interface cells:
    call uncompressm_xedges(region,js,je,ls,le)

    ! interface cells:
    call mixm_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    call GO_Timer_End( itim_dynamum, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine dynamum


  ! ***


  ! return status:
  !   1 : max_nloop exceeded

  pure subroutine advectx_get_nloop( is,ie, is_cyclic, m, am, nloop, status )

    use dims, only : zero, one

    ! --- in/out ---------------------------------

    integer, intent(in)               ::  is, ie
    logical, intent(in)               ::  is_cyclic
    real, intent(in)                  ::  m(is-1:ie+1)
    real, intent(inout)               ::  am(is-1:ie)
    integer, intent(out)              ::  nloop
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    integer, parameter             ::  max_nloop = 50

    ! --- local ----------------------------------

    real              ::  alpha
    real              ::  mx(is-1:ie+1)
    integer           ::  i
    integer           ::  iloop
    logical           ::  cfl_ok

    ! --- begin ----------------------------------

    ! Determine number of loops required to avoid CFLs...
    nloop = 1   ! default run the loop one time!
    alpha = 2.0

    !CMK added oct2003: nloop limit..
    do while ( abs(alpha) >= one .and. nloop < max_nloop )

       ! copy mass to temp array
       mx(is-1:ie+1) = m(is-1:ie+1)
       if ( is_cyclic ) then
         mx(0) = mx(ie)
         mx(ie+1) = mx(1)
       end if

       xloop: do iloop = 1, nloop
          cfl_ok = .true.
          do i=is-1,ie
             if (am(i)>=zero)  then
                alpha=am(i)/mx(i)
             else
                alpha=am(i)/mx(i+1)
             end if
             if((abs(alpha)>=one)) then                   !PB
                cfl_ok = .false.
                exit xloop
             end if
          end do
          ! update mass
          ! note: only is:ie updated, is-1, ie+1 are BC
          mx(is:ie) = mx(is:ie) + am(is-1:ie-1) - am(is:ie)
          if ( is_cyclic ) then   !except for periodic BC...
             mx(0) = mx(ie)
             mx(ie+1) = mx(1)
          end if
       end do xloop

       ! not ok yet ?
       if ( .not. cfl_ok ) then
          ! reduce mass flux
          am(is-1:ie) = am(is-1:ie)*nloop/(nloop+1)
          nloop = nloop + 1
          ! check ...
          if ( nloop == max_nloop ) then
            ! max_nloop exceeded
            status=1; return
          end if
       end if

     end do ! while alpha>1

    ! ok
    status = 0

  end subroutine advectx_get_nloop


  ! ***


  subroutine compressm_xedges(region,js,je,ls,le)
    !
    ! this is for slope only: condense data at the x-edges of the zoom region
    ! to allow for uniform work in dynamu/dynamu1
    ! written by mike botchev, march-june 1999
    !
    use dims,        only : xref, parent, im
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat

    implicit none

    ! input
    integer,intent(in) :: region
    integer,intent(in) :: js
    integer,intent(in) :: je
    integer,intent(in) :: ls
    integer,intent(in) :: le

    ! local
    real,dimension(:,:,:),  pointer    :: m
    integer :: imr,j,l,xref_
    ! start
    m => m_dat(region)%data
    xref_ = xref(region)/xref(parent(region))
    imr = im(region)
    if ( (xref_ == 1 ) .or. ( im(region)/xref(region) == im(1) ) ) then
       if ( okdebug ) print *,"compressm_xedges:", &
            " no refinement or periodic bc's, nothing to do"
       return
    end if
    if ( okdebug ) then
       print *,'compressm_xedges: region=',region,' imr=',imr
       print *,'compressm_xedges:        cells to be updated: ', &
            xref_-1,xref_,imr-xref_+1,imr-xref_+2
    end if

    do j=js,je
       do l=ls,le
          m (xref_,  j,l) = sum(m(1:xref_,j,l))
          m (xref_-1,j,l) = m (0,j,l)
          m (imr-xref_+1,j,l) = sum(m(imr-xref_+1:imr,j,l))
          m (imr-xref_+2,j,l) = m (imr+1,j,l)
       end do
    end do !forall

    nullify(m)

  end subroutine compressm_xedges

  subroutine uncompressm_xedges(region,js,je,ls,le)
    !
    ! distribute data at the x-edges of the zoom region
    ! over the whole interface cell
    ! written by mike botchev, march-june 1999
    !
    use dims,        only : xref, parent, im
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat

    implicit none

    ! input/output
    integer,intent(in) :: region
    integer,intent(in) :: js
    integer,intent(in) :: je
    integer,intent(in) :: ls
    integer,intent(in) :: le

    ! local
    real,dimension(:,:,:),  pointer    :: m
    integer :: imr,j,l,n,xref_
    real    :: m_edge

    ! start

    m => m_dat(region)%data

    xref_ = xref(region)/xref(parent(region))
    imr = im(region)

    if ( (xref_==1) .or. (im(region)/xref(region)==im(1)) ) then
       if ( okdebug ) print *,"uncompressm_xedges:", &
            " no refinement or periodic bc's, nothing to do"
       return
    end if
    if ( okdebug ) then
       print *,'uncompressm_xedges: region=',region,' imr=',imr
       print *,'uncompressm_xedges: cells to be updated: ',1,':', &
            xref_,' ',imr-xref_+1,':',imr
    end if

    ! forall(j=js:je,l=ls:le)
    do j=js,je
       do l=ls,le
          m_edge = m (xref_,j,l)
          m (1:xref_,j,l) = m_edge/xref_
          m_edge = m (imr-xref_+1,j,l)
          m (imr-xref_+1:imr,j,l) = m_edge/xref_
       end do !forall
    end do !forall

    nullify(m)
  end subroutine uncompressm_xedges



  subroutine putm_xedges(region)
    !
    ! passes values along the x-boundaries to all the children
    ! written by mike botchev, march-june 1999
    ! structure implemented by MK, dec 2002
    !
    use dims,          only : im, jm, lm, xref, yref, zref
    use dims,          only : ibeg, jbeg, lbeg, iend, nregions
    use dims,          only : children
    use global_data,   only : mass_dat
    use MeteoData    , only : m_dat

    implicit none

    ! in/out
    integer,intent(in)   :: region

    ! local
    real,dimension(:,:,:),  pointer     :: m,mc
    integer            :: child,ichild,j,jp,l,lp,jmc,lmc,n,i
    integer            :: xref_,yref_,zref_
    real               :: yzref,xyzref
    real,dimension(:,:),allocatable :: toc0,toc1,tocm,tocm1
    integer :: communicator,root_id,lglob,nglob,lmr,nt

    ! start
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
               print *,'putm_xedges: child ',child,' skipped'
          cycle
       endif
       m  => m_dat(region)%data
       mc => m_dat(child )%data
       lmr = lm(region)
       ! loop through the cells of x-walls of the child
       jmc = jm(child)
       lmc = lmr
       yzref = 1./(yref_*zref_)
       xyzref = 1./(xref_*yref_*zref_)
       allocate(toc0(jmc,lmc))
       allocate(toc1(jmc,lmc))
       allocate(tocm(jmc,lmc))
       allocate(tocm1(jmc,lmc))
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
       deallocate(toc0)
       deallocate(toc1)
       deallocate(tocm)
       deallocate(tocm1)
       nullify(m)
    end do
  end subroutine putm_xedges


  ! ================================================================================
  ! ===
  ! === advect m tools
  ! ===
  ! ================================================================================


  !
  ! distribute data at the x-edges of the zoom region
  ! over the whole interface cell
  ! written by mike botchev, march-june 1999
  ! adapted for the CRAY by Maarten Krol, 3-2000
  !

  subroutine mixm_edges( region, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : im,jm,lm, xref, yref, parent
    use dims,        only : touch_np, touch_sp, zref
    use global_data, only : mass_dat, region_dat
    use MeteoData  , only : m_dat
    use toolbox,     only : escape_tm

    ! input
    integer,intent(in)            ::  region
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/mixm_edges'

    ! local
    real,dimension(:,:,:),  pointer    :: m
    integer,dimension(:,:), pointer    :: edge

    integer    :: imr, i, j, l, n, jmr, lmr, ip, jp, imp, jmp, nt
    integer    :: xref_, yref_, zref_, xyzref_, xyref_
    real       :: m_edge, rm_edge, rxm_edge, rym_edge, rzm_edge
    logical    :: xedges , ynp,ysp

    ! relaxation factor for slopes; relax=0.0 means upwinding
    real,parameter:: relax=1.0
    real,dimension(:,:),allocatable   :: to_parent
    real,dimension(:,:),allocatable   :: to_parentx
    real,dimension(:,:),allocatable   :: to_parenty
    real,dimension(:,:),allocatable   :: to_parentz

    ! start

    call GO_Timer_Start( itim_mixm_edges, status )
    IF_NOTOK_RETURN(status=1)

    xedges = .true.
    ysp = .true.  !mix also NP/SP?
    ynp = .true.

    imr = im(region)
    jmr = jm(region)

    m => m_dat(region)%data
    lmr = lm(region)

    edge => region_dat(region)%edge

    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))

    xyref_ = xref_*yref_

    if ( okdebug ) then
       print *,'mix_edges: region=',region
    end if
    if ( (xref_ == 1) .or. (im(region)/xref(region) == im(1)) ) then
       if ( okdebug ) then
         print *, 'mix_edges: no refinement or periodic bc, skipping x-walls'
       end if
       xedges = .false.
       if ( region == 1 ) then
          if ( okdebug ) then
            print *, 'mix_edges: region = 1, returning'
          end if
          call GO_Timer_End( itim_mixm_edges, status )
          IF_NOTOK_RETURN(status=1)
          return
       end if
    end if

    if ( touch_sp(region) == 1) then
       ysp = .false.
       if ( okdebug ) then
         print *, 'mix_edges: SP touching...skip SP y-walls'
       end if
    end if
    if(touch_np(region) == 1) then
       ynp = .false.
       if ( okdebug ) then
         print *, 'mix_edges: NP touching...skip NP y-walls'
       end if
    end if
    imp = imr/xref_   !in 'coarse' resolution
    jmp = jmr/yref_

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( imr, jmr, lmr ) &
    !$OMP   shared ( imp, jmp ) &
    !$OMP   shared ( xref_, yref_, xyref_ ) &
    !$OMP   shared ( edge ) &
    !$OMP   shared ( xedges, ysp, ynp ) &
    !$OMP   shared ( m ) &
    !$OMP   private ( i, j, l ) &
    !$OMP   private ( ip, jp ) &
    !$OMP   private ( to_parent )
    !$OMP   DO
    do l=1,lmr !WP! depends on tracers or levels
       allocate(to_parent(imp,jmp))
       to_parent = 0.0
       do j=1,jmr
          jp = 1 + (j-1)/yref_
          do i=1,imr
             if ( edge(i,j) /= -1 ) cycle   !do only for edge of child...
             ip = 1 + (i-1)/xref_
             to_parent(ip,jp) = to_parent(ip,jp) + m(i,j,l)/xyref_
          end do
       end do
       if ( xedges ) then     !mix the x-wall of the zoom region.
          do j = 1,jmr
             jp = 1 + (j-1)/yref_
             do i = 1, xref_
                m(i,j,l) = to_parent(1,jp)
             end do
             do i = imr-xref_+1,imr
                m(i,j,l) = to_parent(imp,jp)
             end do
          end do
       end if    !xedges
       if ( ysp ) then
          do j = 1,yref_
             do i = 1, imr
                ip = 1 + (i-1)/xref_
                m(i,j,l) = to_parent(ip,1)
             end do
          end do
       end if
       if ( ynp ) then
          do j = jmr-yref_+1,jmr
             do i = 1, imr
                ip = 1 + (i-1)/xref_
                m(i,j,l) = to_parent(ip,jmp)
             end do
          end do
       end if
       deallocate(to_parent)
    end do  !l
    !$OMP   END DO
    !$OMP END PARALLEL

    nullify(m)
    nullify(edge)

    call GO_Timer_End( itim_mixm_edges, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine mixm_edges

  ! ***

  !
  ! update rm of its parent in correspondence with itself
  ! written by mike botchev, march-june 1999
  !

  subroutine updatem_parent( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,         only : im,jm,lm, xref, yref, zref, parent
    use dims,         only : jbeg, jend, ibeg, iend
    use global_data  ,only : mass_dat
    use MeteoData   , only : m_dat
    use toolbox,      only : escape_tm

    ! input
    integer,intent(in)   :: region
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/updatem_parent'

    ! local
    real,dimension(:,:,:),pointer            :: m,mp
    real,dimension(:,:,:),allocatable        :: to_parent

    real    :: dummy
    real    :: sum_parent, sum_to_parent, sum_parent_all, sum_to_parent_all
    integer :: i, iend_loop, iip, ip, j, jp, l, lp, n
    integer :: lb, le, ib, ie, jb, je
    integer :: my_parent, xref_, yref_, zref_
    integer :: imp, jmp, lmp, imr, jmr, lmr, nt
    real    :: m_value, mxval, sum_mass
    integer,dimension(3) :: mxloc

    ! start

    call GO_Timer_Start( itim_updatem_parent, status )
    IF_NOTOK_RETURN(status=1)

    if (region == 1) return
    imr = im(region)
    jmr = jm(region)
    my_parent = parent(region)
    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    zref_ = zref(region)/zref(my_parent)

    m  => m_dat(region   )%data
    mp => m_dat(my_parent)%data

    lmr = lm(region)
    if ( okdebug ) then
       print *,'updatem_parent: my_parent=',my_parent, &
            ' x-,y-,zref_: ',xref_,yref_,zref_
    end if

    imp = im(region)/xref_
    jmp = jm(region)/yref_
    !lmp = lm(region)/zref_
    lmp = lmr     ! CMK changed: when paralel over levels lmp can differ from lmr

    allocate(to_parent(imp,jmp,lmp))
    to_parent = 0.0
    do i=1,imr
       ip = 1 + (i-1)/xref_
       do j=1,jmr
          jp = 1 + (j-1)/yref_
          do l=1,lmr
             lp = 1 + (l-1)/zref_
             to_parent(ip,jp,lp) = to_parent(ip,jp,lp) + m(i,j,l)
          end do
       end do
    end do
    if ( ibeg(region) < iend(region) ) then   !no dateline crossing!
       mp(ibeg(region):iend(region),jbeg(region):jend(region),1:lmr) = &
            to_parent
    else
       mp(ibeg(region):im(region),jbeg(region):jend(region),1:lmr)  = &
            to_parent(1:im(region)-ibeg(region)+1,:,:)
       mp(1:iend(region),jbeg(region):jend(region),1:lmr)  = &
            to_parent(im(region)-ibeg(region)+2:,:,:)
    end if

    deallocate(to_parent)

    nullify(m)
    nullify(mp)

    call GO_Timer_End( itim_updatem_parent, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine updatem_parent


  ! ================================================================================
  ! ===
  ! === advect m y
  ! ===
  ! ================================================================================


  !
  ! set parameters for advecty
  ! written by patrick berkvens and mike botchev, march-june 1999
  ! updated and modified by MK, dec 2002
  !

  subroutine advectmyzoom( region, cfl_ok, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,          only : xref, yref, zref, tref, im, jm, lm
    use dims,          only : zoom2D
    use dims,          only : touch_sp, touch_np, adv_scheme, okdebug
    use dims,          only : parent, zero, xcyc
    use global_data,   only : wind_dat, mass_dat
    use toolbox,       only : escape_tm

    ! input/output
    integer,intent(in) :: region
    logical, intent(out)                :: cfl_ok
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/advectyzoom'

    ! local
    real,dimension(:,:,:),pointer               :: bm
    real,dimension(:,:),allocatable             :: bm0,bm1

    integer            :: is,ie,ls,le,n,q
    integer            :: imr,jmr,lmr,tref_,xref_,yref_,zref_
    logical            :: y_encountered
    character(len=1)   :: dir

    ! start

    call GO_Timer_Start( itim_advectmyzoom, status )
    IF_NOTOK_RETURN(status=1)

    call putm_yedges(region)

    tref_ = tref(region)/tref(parent(region))
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region);  jmr = jm(region);  lmr = lm(region)

    allocate(bm0(0:im(region)+1,0:lmr+1))
    allocate(bm1(0:im(region)+1,0:lmr+1))

    bm => wind_dat(region)%bm_t

    ! determine the scope for advecty:

    if ( region == 1 ) then
       xref_ = 0; zref_ = 0  ! to have is/ie and ls/le properly computed
    end if
    ! find q - the place in the splitorderzoom
    q=regionm_status(region)/((nsplitsteps/2)*tref_)
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
       case default
          print *,'advectmyzoom: strange value in splitorderzoom(',region,',',  &
               q*(nsplitsteps/2)*tref_+n,'): ',q
          call escape_tm('advectmyzoom: Error in do_next3')
       end select

    end do

    if ( ( mod(regionm_status(region),(nsplitsteps/2)*tref_) >= n_operators ) .and. &
         ( adv_scheme == 'slope' ) ) then
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

    call dynamvm( region, is,ie,ls,le, cfl_ok, status )
    IF_NOTOK_RETURN(status=1)

    if ( (mod(regionm_status(region),(nsplitsteps/2)*tref_) >= n_operators ) .and. &
         (adv_scheme=='slope') ) then

       ! for slopes only: recreate velocities via the x-edges
       if ( touch_sp(region) /= 1 ) bm(:,      yref_,:) = bm0
       if ( touch_np(region) /= 1 ) bm(:,jmr-yref_+2,:) = bm1

    end if
    deallocate(bm0)
    deallocate(bm1)
    nullify(bm)

    call GO_Timer_End( itim_advectmyzoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectmyzoom


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

  subroutine dynamvm( region, is,ie,ls,le, cfl_ok, status )

    use omp_lib
    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : xref, yref, zref, tref, im, jm, lm
    use dims,        only : parent, nregions
    use dims,        only : touch_sp, touch_np, limits, zero, one
    use redgridZoom, only : nred, jred, imredj, clustsize
    use global_data, only : mass_dat, wind_dat
    use MeteoData  , only : m_dat
    use zoom_tools,  only : mix_edges
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_comm,    only : barrier_k
    use mpi_const,   only : ntracetloc
#endif
    use chem_param,  only : ntracet

    ! input/output
    integer,intent(in) :: region
    integer,intent(in) :: is
    integer,intent(in) :: ie
    integer,intent(in) :: ls
    integer,intent(in) :: le
    logical, intent(out)                :: cfl_ok
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/dynamvm'

    ! local
    real,dimension(:,:,:),  pointer           :: m,bm
    real,dimension(:,:,:),allocatable         :: mnew
    integer                                   :: i,j,je,js,l,n,iee,iss
    integer                                   :: imr,imr2,jmr,lmr
    integer                                   :: tref_,xref_,yref_,zref_
    real                                      :: sfs,sfzs,sfn,sfzn
    real                                      :: beta, max_beta
    real,dimension(:,:), allocatable :: mx
    logical                          :: cfl_okl
    integer        :: lrg, redfact, ixe, ixs
    real           :: summ
    ! openmp:
    integer        ::  l_mloop_max
    real           ::  l_xim

    ! start

    call GO_Timer_Start( itim_dynamvm, status )
    IF_NOTOK_RETURN(status=1)

    ! compute refinement factors with respect to the parent

    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    tref_ = tref(region)/tref(parent(region))
    imr=im(region);jmr=jm(region);lmr=lm(region)

    call compressm_yedges(region,is,ie,ls,le)


    allocate(mnew(im(region),jm(region),lmr))

    m => m_dat(region)%data
    bm => wind_dat(region)%bm_t


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

    ! reduction variables:
    l_mloop_max = 0
    l_xim = 0.0
    cfl_okl = .true.

    ! CMK splitted the loop and simplified code because of recurring OMP problems at huygens
    ! CMK loop was executed by all threads.
    ! loop over tracers and vertical layers
    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( imr, jmr ) &
    !$OMP   shared ( is, ie, js, je, ls, le ) &
    !$OMP   shared ( touch_sp, touch_np ) &
    !$OMP   shared ( m, mnew, bm ) &
    !$OMP   reduction ( .and. : cfl_okl ) &
    !$OMP   reduction ( max : l_mloop_max ) &
    !$OMP   reduction ( max : l_xim ) &
    !$OMP   private ( i, j, l ) &
    !$OMP   private ( beta, max_beta ) &
    !$OMP   private ( mx )
    allocate(mx(is:ie,js-1:je+1))
    !$OMP   DO
    xloop: do l=ls,le
       !print *, l, omp_get_thread_num()+1, 'out of ', omp_get_num_threads()
       ! compute all inner fluxes
       ! f(*,j,*) is flux entering j-th cell from below?cmk
       beta = 2.0
       max_beta = 0.0
       mx(is:ie,js-1:je+1) = m(is:ie,js-1:je+1,l)
       do j=js+1,je
          do i=is,ie
             if (bm(i,j,l)>=zero) then
                beta=bm(i,j,l)/mx(i,j-1)
             else
                beta=bm(i,j,l)/mx(i,j)
             end if
             if (abs(beta)>=one) then
                cfl_okl = cfl_okl .and. .false.
             end if
             max_beta = max(max_beta, abs(beta))
          end do  !i
       end do !j
       if ( region == 1.or.touch_sp(region) == 1) then
          do i=is,ie
             if (bm(i,2,l)>=zero) then
                beta=bm(i,2,l)/mx(i,1)
             else
                beta=bm(i,2,l)/mx(i,2)
             end if
             if (abs(beta)>=one) then
                cfl_okl = cfl_okl .and. .false.
             end if
             max_beta = max(max_beta, abs(beta))
          end do
       else   !zoom region not touching south pole
          j = js
          do i=is,ie   !no reduced grid allowed
             if (bm(i,j,l)>=zero) then
                beta=bm(i,j,l)/mx(i,j-1)
             else
                beta=bm(i,j,l)/mx(i,j)
             end if
             if (abs(beta)>=one) then
                cfl_okl = cfl_okl .and. .false.
             end if
             max_beta = max(max_beta, abs(beta))
          end do
       end if ! compute boundary fluxes south pole...
       if ( region == 1.or.touch_np(region) == 1) then       ! north pole
          do i=is,ie
             if (bm(i,jmr,l)>=zero) then
                beta=bm(i,jmr,l)/mx(i,jmr-1)
             else
                beta=bm(i,jmr,l)/mx(i,jmr)
             end if
             if (abs(beta)>=one) then
                cfl_okl = cfl_okl .and. .false.
             end if
             max_beta = max(max_beta, abs(beta))
           end do
       else   !zoom region not touching north pole
          j = je+1
          do i=is,ie   !no reduced grid allowed
             if (bm(i,j,l)>=zero) then
                beta=bm(i,j,l)/mx(i,j-1)
             else
                beta=bm(i,j,l)/mx(i,j)
             end if
             if (abs(beta)>=one) then
                cfl_okl = cfl_okl .and. .false.
             end if
             max_beta = max(max_beta, abs(beta))
          end do
       end if ! compute boundary fluxes north pole...
       l_xim = max( l_xim, max_beta )

    end do xloop ! end of l-loop over vertical layers....
    !$OMP   END DO
    deallocate(mx)
    !$OMP END PARALLEL
    if ( cfl_okl ) then
       !$OMP PARALLEL &
       !$OMP   default ( none ) &
       !$OMP   shared ( region ) &
       !$OMP   shared ( imr, jmr ) &
       !$OMP   shared ( is, ie, js, je, ls, le ) &
       !$OMP   shared ( touch_sp, touch_np ) &
       !$OMP   shared ( m, mnew, bm ) &
       !$OMP   private ( l )
       !$OMP   DO
       do l = ls, le



            ! calculate new air mass distribution

            if (region==1) then
               mnew(1:imr,1:jmr,l) = m(1:imr,1:jmr,l) + &
                    bm(1:imr,1:jmr,l) - bm(1:imr,2:jmr+1,l)
            else if ( touch_sp(region) == 1) then
               mnew(1:imr,1:je,l) = m(1:imr,1:je,l) + &
                    bm(1:imr,1:je,l) - bm(1:imr,2:je+1,l)
            else if ( touch_np(region) == 1) then
               mnew(1:imr,js:jmr,l) = m(1:imr,js:jmr,l) + &
                    bm(1:imr,js:jmr,l) - bm(1:imr,js+1:jmr+1,l)
            else
               mnew(is:ie,js:je,l) = m(is:ie,js:je,l) + &
                    bm(is:ie,js  :je  ,l) - bm(is:ie,js+1:je+1,l)
            end if

            ! store new air mass in m array

            if ( region == 1) then
               m(1:imr,1:jmr,l) = mnew(1:imr,1:jmr,l)
            else if ( touch_sp(region) == 1) then
               m(1:imr,1:je,l) = mnew(1:imr,1:je,l)
            else if ( touch_np(region) == 1) then
               m(1:imr,js:jmr,l) = mnew(1:imr,js:jmr,l)
            else
               m(is:ie,js:je,l) = mnew(is:ie,js:je,l)
            end if
            !


       enddo
       !$OMP   END DO
       !$OMP END PARALLEL
    end if  ! cfl_okl


    mloop_max(region,2) = 1
    xim      (region,2) = l_xim
    cfl_ok              = cfl_okl

    ! clear:
         deallocate(mnew)
    nullify(m)
    nullify(bm)

    if ( cfl_ok ) then
      call uncompressm_yedges( region, is,ie,ls,le )

      call mixm_edges( region, status )
      IF_NOTOK_RETURN(status=1)

    end if

    call GO_Timer_End( itim_dynamvm, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine dynamvm


  ! ***


  subroutine putm_yedges(region)
    !
    ! passes values along the y-boundaries to all the children
    ! written by mike botchev, march-june 1999
    ! adapted for cray-run Maarten Krol, march, 2000
    ! touch_np, touch_sp implemented...MK, nov, 2000
    !

    use dims,          only : xref, yref, zref, im, jm, lm
    use dims,          only : touch_sp, touch_np, children
    use dims,          only : ibeg, jbeg, lbeg, jend, nregions
    use global_data,   only : mass_dat
    use MeteoData    , only : m_dat

    implicit none

    ! input/output
    integer,intent(in)   :: region
    ! local
    real,dimension(:,:,:),  pointer    :: m, mc
    real,dimension(:,:),allocatable    :: toc0,toc1,tocm,tocm1

    integer :: child, ichild, j, ip, l, lp, imc, lmc, n, i
    integer :: xref_, yref_, zref_
    real    :: xzref, xyzref
    integer :: lmr

    ! start

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
               print *,'putm_yedges: child ',child,' sp  skipped'
       end if
       if ( touch_np(child) == 1) then
          if (okdebug)  &
               print *,'putm_yedges: child ',child,' np  skipped'
       end if

       m  => m_dat(region)%data
       mc => m_dat(child )%data
       lmr = lm(region)
       ! loop through the cells of y-walls of the child
       imc = im(child)
       lmc = lmr
       xzref = 1./(xref_*zref_)
       xyzref = 1./(xref_*yref_*zref_)
       allocate(toc0(imc,lmc))
       allocate(toc1(imc,lmc))
       allocate(tocm(imc,lmc))
       allocate(tocm1(imc,lmc))
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
       deallocate(toc0)
       deallocate(toc1)
       deallocate(tocm)
       deallocate(tocm1)
       nullify(m)
    end do

  end subroutine putm_yedges

  subroutine compressm_yedges(region,is,ie,ls,le)
    !
    ! this is for slope only: condense data at the y-edges of the zoom region
    ! to allow for uniform work in dynamv
    ! written by mike botchev, march-june 1999
    ! modified by MK, dec 2002
    !
    use dims,        only : yref, parent, jm, touch_sp, touch_np
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat

    implicit none

    ! input/output
    integer,intent(in) :: region
    integer,intent(in) :: is
    integer,intent(in) :: ie
    integer,intent(in) :: ls
    integer,intent(in) :: le

    ! local
    real,dimension(:,:,:),  pointer           :: m
    integer                                   :: jmr,i,l,n,yref_

    ! start

    yref_ = yref(region)/yref(parent(region))
    jmr = jm(region)
    if ((yref_==1).or.(region==1)) then
       if ( okdebug) print *,'compressm_yedges: no refinement, nothing to do'
       return
    end if
    if ( okdebug) then
       print *,'compressm_yedges: region=',region,' jmr=',jmr
       print *,'compressm_yedges: cells to be updated: ', &
            yref_-1,yref_,jmr-yref_+1,jmr-yref_+2
       if ( touch_sp(region) == 1) &
            print *,'compress_yedges: SP will be skipped here'
       if ( touch_np(region) == 1) &
            print *,'compress_yedges: NP will be skipped here'
    end if

    m => m_dat(region)%data

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
       end do
    end do

    nullify(m)

  end subroutine compressm_yedges



  subroutine uncompressm_yedges(region,is,ie,ls,le)
    !
    ! distribute data at the y-edges of the zoom region
    ! over the whole interface cell
    ! written by mike botchev, march-june 1999
    !
    use dims,        only : yref, parent, jm, touch_sp, touch_np
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat

    implicit none

    ! input/output
    integer,intent(in) :: region
    integer,intent(in) :: is
    integer,intent(in) :: ie
    integer,intent(in) :: ls
    integer,intent(in) :: le

    ! local
    real,dimension(:,:,:),  pointer           :: m
    integer                                   :: jmr,i,l,n,yref_
    real                                      :: m_edge,rm_edge

    ! start

    yref_ = yref(region)/yref(parent(region))
    jmr = jm(region)
    if ((yref_==1).or.(region==1)) then
       if ( okdebug) print *,'uncompressm_yedges: no refinement, nothnig to do'
       return
    end if
    if ( okdebug) then
       print *,'uncompressm_yedges: region=',region,' jmr=',jmr
       print *,'uncompressm_yedges: cells to be updated: ',1,':', &
            yref_,' ',jmr-yref_+1,':',jmr
       if ( touch_sp(region) == 1) &
            print *,'uncompressm_yedges: SP will be skipped here'
       if ( touch_np(region) == 1) &
            print *,'uncompressm_yedges: NP will be skipped here'
    end if

    m => m_dat(region)%data

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
       end do
    end do

    nullify(m)

  end subroutine uncompressm_yedges


  ! ================================================================================
  ! ===
  ! === advect m z
  ! ===
  ! ================================================================================


  !
  ! set parameters for advectz
  ! written by patrick berkvens and mike botchev, march-june 1999
  !

  subroutine advectmzzoom( region, cfl_ok, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : im, jm, lm, xref, yref, zref, tref
    use dims,        only : parent, zoom2D, touch_sp, touch_np, xcyc
    use toolbox,     only : escape_tm

    !input/output
    integer,intent(in) :: region
    logical, intent(out)                :: cfl_ok
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/advectzzoom'

    ! local
    integer            :: is,ie,js,je,n,q
    integer            :: imr,jmr,lmr,tref_,xref_,yref_,zref_
    logical            :: z_encountered
    character(len=1)   :: dir

    ! start

    call GO_Timer_Start( itim_advectmzzoom, status )
    IF_NOTOK_RETURN(status=1)

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
    q=regionm_status(region)/((nsplitsteps/2)*tref_)
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
       case default
          print *,'advectmzzoom: strange value in splitorderzoom(',region,',',  &
               q*(nsplitsteps/2)*tref_+n,'): ',q
          call escape_tm('advectmzzoom: Error in advectz')
       end select

    end do

    call dynamwm( region, is,ie,js,je, cfl_ok, status )
    IF_NOTOK_RETURN(status=1)

    call GO_Timer_End( itim_advectmzzoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advectmzzoom



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
  ! AJS: removed everything that does not take into account that:
  !  o levels are always 1:lmr
  !  o advz_n = 1 (thus, just check if for the current values cfl is ok)
  !
  !-----------------------------------------------------------------------!

  subroutine dynamwm( region, is,ie,js,je, cfl_ok, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims,        only : im, jm, lm
    use dims,        only : zero, one
    use dims,        only : nregions
    use global_data, only : wind_dat, mass_dat
    use MeteoData  , only : m_dat
    use toolbox,     only : escape_tm

    ! input,output
    integer,intent(in) :: region
    integer,intent(in) :: is
    integer,intent(in) :: ie
    integer,intent(in) :: js
    integer,intent(in) :: je
    logical, intent(out)                :: cfl_ok
    integer,intent(out)                 :: status

    ! const
    character(len=*), parameter  ::  rname = mname//'/dynamwm'

    ! local
    real,dimension(:,:,:),pointer     :: m,cm
    real                              ::  gamma
    real                              ::  l_gamma
    integer                           :: i,j,l
    integer                           :: imr,jmr,lmr

    ! start

    call GO_Timer_Start( itim_dynamwm, status )
    IF_NOTOK_RETURN(status=1)

    if ( okdebug ) print *,'dynamwm: region=',region,' is,ie,js,je=',is,ie,js,je
    if ( ( region < 0 ) .or. ( region > nregions ) ) &
         call escape_tm( 'dynamw: STOP, illegal number of region !!!')

    ! local grid size:
    imr=im(region) ; jmr=jm(region) ; lmr=lm(region)

    ! point to region arrays:
    m  => m_dat(region)%data
    cm => wind_dat(region)%cm_t

    ! check bottom:
    if ( any(cm(:,:,0) /= 0.0) ) then
      write (gol,'("cm bottom flux should be zero")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! check top:
    if ( any(cm(:,:,lmr) /= 0.0) ) then
      write (gol,'("cm  top flux should be zero")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! reset gamma
    l_gamma = 0.0

    ! compute f, pf, fx and fy
    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( imr, jmr, lmr ) &
    !$OMP   shared ( m, cm ) &
    !$OMP   reduction (max:l_gamma) &
    !$OMP   private ( gamma ) &
    !$OMP   private ( i, j, l )
    !$OMP   DO
    do l = 0, lmr      ! 1,lmm1
       do j = 1,jmr
          do i = 1,imr
             if ( cm(i,j,l) == zero ) then
                gamma = 0.0
             else if ( cm(i,j,l) > zero ) then
                gamma = cm(i,j,l)/m(i,j,l)
             else
                gamma = cm(i,j,l)/m(i,j,l+1)
             end if
             ! update maximum:
             l_gamma = max( l_gamma, abs(gamma) )
             ! leave after first violation:
             if ( l_gamma >= one ) exit
          end do  ! i
          ! leave after first violation:
          if ( l_gamma >= one ) exit
       end do  ! j
    end do  ! l
    !$OMP   END DO
    !$OMP END PARALLEL

    ! update maximum Courrant number:
    xim(region,3) = max( xim(region,3), l_gamma )
    ! update counter (only one extra cfl loop step):
    mloop_max(region,3) = max( mloop_max(region,3), 1 )

    ! do not allow iterations in Z
    if ( l_gamma >= one ) then

      ! reset global flag:
      cfl_ok = .false.

    else

      ! new air mass:
      ! Sourish Basu: write a loop that is vectorizable
      m(is:ie, js:je, 1:lmr) = m(is:ie, js:je, 1:lmr) + cm(is:ie, js:je, 0:lmr-1) - cm(is:ie, js:je, 1:lmr)
!      !$OMP PARALLEL &
!      !$OMP   default ( none ) &
!      !$OMP   shared ( is, ie, js, je, lmr ) &
!      !$OMP   shared ( m, cm ) &
!      !$OMP   private ( i, j, l )
!      !$OMP   DO
!      do l = 1, lmr
!        do j = js, je
!          do i = is, ie
!            m(i,j,l) = m(i,j,l) + cm(i,j,l-1) - cm(i,j,l)
!          end do
!        end do
!      end do
!      !$OMP   END DO
!      !$OMP END PARALLEL

      call mixm_edges( region, status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! clear:
    nullify(m)
    nullify(cm)

    call GO_Timer_End( itim_dynamwm, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine dynamwm


  ! ***


  subroutine done_cfl
    use dims, only : nregions
    use ParTools, only : myid, root
    implicit none
    integer :: region
    if (myid == root) then
        print*,' '
        print*,'CFL info from program advectm:'
    endif
    do region=1, nregions
      deallocate(oldmass_dat(region)%m) !WP! deallocation on all PE's
      deallocate(oldmass_dat(region)%am)
      deallocate(oldmass_dat(region)%bm)
      deallocate(oldmass_dat(region)%cm)
      if (myid == root) then
          print '(a,3i4,f10.4)', ' exitus: x: region, nxim, mloop_max, xim',  &
             region, nxim(region,1), &
             mloop_max(region,1), xim(region,1)
          print '(a,3i4,f10.4)', ' exitus: y: region, nxim, mloop_max, xim',  &
             region, nxim(region,2), &
             mloop_max(region,2), xim(region,2)
          print '(a,3i4,f10.4)', ' exitus: z: region, nxim, mloop_max, xim',  &
             region, nxim(region,3), &
             mloop_max(region,3), xim(region,3)
       endif
    end do

    end subroutine done_cfl

end module advectm_cfl
