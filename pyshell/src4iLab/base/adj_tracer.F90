!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_Tracer

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------

  private

  public :: Adj_Tracer_Init, Adj_Tracer_Done
  public :: Adj_Tracer_Model


  ! --- const ----------------------------------

  ! module name
  character(len=*), parameter  ::  mname = 'adj_Tracer'


  ! --- local ------------------------------------

  ! timers:
  integer     ::  itim_tracer
  integer     ::  itim_timestep_init, itim_timestep_run, itim_timestep_user_output
  integer     ::  itim_model_init


contains


  ! ============================================================================


  subroutine Adj_Tracer_Init( status )

    use GO, only : GO_Timer_Def
    use Meteo, only : Meteo_Def_Timers

    ! --- in/out -------------------------------

    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Adj_Tracer_Init'

    ! --- begin --------------------------------

    ! define timers:
    call GO_Timer_Def( itim_tracer               , 'adj tracer model'        , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_model_init           , 'adj model init'          , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_timestep_init        , 'adj timestep init'       , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_timestep_run         , 'adj timestep run'        , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_timestep_user_output , 'adj timestep user output', status )
    IF_NOTOK_RETURN(status=1)

    ! initialize meteo module:
    call Meteo_Def_Timers( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Adj_Tracer_Init


  ! ***


  subroutine Adj_Tracer_Done( status )

    ! --- in/out -------------------------------

    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Adj_Tracer_Done'

    ! --- begin --------------------------------

    ! ok
    status = 0

  end subroutine Adj_Tracer_Done


  ! ***


  subroutine adj_Tracer_Model( status )

    use GO                  , only : goLabel
    use GO                  , only : GO_Timer_Start, GO_Timer_End
    use GO                  , only : TDate, NewDate, TIncrDate, IncrDate
    use GO                  , only : operator(-), operator(+), operator(*), operator(==), max
    use GO                  , only : rTotal
    use GO                  , only : wrtgol
    use GO                  , only : SystemDate
    use Dims,                 only : nregions
    use Dims,                 only : region_status => status
    use Dims,                 only : revert
    use Dims,                 only : nsec_read
    use Dims,                 only : newsrun, newmonth, newday
    use Dims,                 only : itau, itaui, itaue, itaur
    use Dims                , only : idate, idatei, idatee
    use Dims                , only : tref
    use Dims,                 only : ndyn
    use Dims,                 only : kmain, ndyn_max
    use Meteo               , only : Meteo_Init, Meteo_Done
    use Meteo               , only : Meteo_Setup_Mass
    use AdvectM_CFL         , only : Setup_MassFlow
    use Meteo               , only : Meteo_Setup_Other
    use adj_sources_sinks,    only : adj_trace0
    use datetime,             only : inctime, tstamp
    use Zoom_Tools          , only : Zoom_Tools_Init, Zoom_Tools_Done
    use zoom_tools,           only : coarsen_region
    use advect,               only : dynam0
    use advect_tools,         only : coarsen_ambm
    use var4d,                only : steps_region
    use advectm_cfl,          only : check_cfl
    use toolbox,              only : escape_tm
    use os_specs,             only : WRITE_STR_LEN
#ifndef without_chemistry
    use chemistry,            only : read_chemistry_fields
#endif

    use InitExit,             only : Start_TM5, Exit_TM5
    use initexit,             only : exit_TM5
    use adj_modelIntegration, only : adj_proces_region
    use adj_modelIntegration, only : adj_Proces_Init, adj_Proces_Done
    use adj_modelIntegration, only : adj_Proces_Setup
    use adj_modelIntegration, only : output_after_step

    use ParTools            , only : myid, root

    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Adj_Tracer'

    ! --- local -------------------------------

    integer             ::  region
    integer             ::  n
    type(TDate)         ::  t
    type(TDate)         ::  t0, t1
    type(TIncrDate)     ::  dt
    type(TDate)         ::  tdyn, tend, tr(2)
    integer             ::  nhalf
    type(TDate)         ::  tread1, tread2
    logical             ::  check_pressure
    character(len=WRITE_STR_LEN)  ::  line

    integer   ::  istep

    ! --- begin ---------------------------------------

    ! start label:
    call goLabel( timer=itim_tracer )

    !---------------!
    ! model startup !
    !---------------!

    ! start label:
    call GO_Timer_Start( itim_model_init, status )
    IF_NOTOK_RETURN(status=1)

    ! tracer for backward run:
    revert=-1

    ! initialize zoom tools:
    call Zoom_Tools_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! setup meteo input
    call Meteo_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init processes
    call adj_Proces_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! set-up and read user input;
    ! return time interval for which meteo was read:
    call Start_TM5( tread1, tread2, status )
    IF_NOTOK_RETURN(status=1)

    ! end label:
    call GO_Timer_End( itim_model_init, status )
    IF_NOTOK_RETURN(status=1)

    !-----------!
    ! main loop !
    !-----------!
    ! current time (begin of dynamics step)
    tdyn = NewDate( time6=idate  )
    tend = NewDate( time6=idatee )

    ! init half-step counter:
    nhalf = 0

    !! counter ...
    !istep = 0

    t0 = SystemDate()
    call wrtgol( 'Wall time before     adjoint run =  ', t0 ); call goPr

    ! time loop
    do
       if ( (revert*itau) >= (revert*itaue) ) exit

       ! * timestep init

       ! start label:
       call GO_Timer_Start( itim_timestep_init, status )
       IF_NOTOK_RETURN(status=1)

      ! next half step
      nhalf = modulo(nhalf,2) + 1

      !! next step:
      !istep = istep + 1

      ! info ...
!      if ( nhalf == 1 ) then
!        write (gol,'("")'); call goPr
!        call wrtgol( '>>> timestep from ', tdyn ); call goPr
!        write (gol,'("")'); call goPr
!        call wrtgol( 'current data valid for ', tread1, ' to ', tread2 ); call goPr
!      end if

      ! first set-up of chemistry related work...repeated if newmonth
      if ( newmonth .and. (.not. newsrun) ) then
        call adj_trace0
      end if

#ifndef without_chemistry
      if (newday .or. newsrun) then
        call read_chemistry_fields(-1, status)
        IF_NOTOK_RETURN(status=1)
      end if
#endif

      ! reached begin of time interval for which meteo is valid ? then setup new meteo:
      if ( tdyn == tread2 ) then

        ! setup meteo data for next interval;
        ! nread is the length (in seconds) of the interval in which
        ! surface pressure is interpolated (and mass fluxes are constant)
        tread1 = tdyn
        tread2 = max( tend, tdyn - IncrDate(sec=nsec_read) )

        ! setup mass and mass fluxes:
        !  o skip first time; already called in 'initexit/start'
        !  o check pressure implied by advection if advection is applied
#ifdef without_advection
        check_pressure = .false.
#else
        check_pressure = .true.
#endif
        call Meteo_Setup_Mass( tread1, tread2, status, check_pressure=check_pressure )
        IF_NOTOK_RETURN(status=1)

#ifndef without_advection
        ! determine dynamic timestep ndyn for this interval [tread1,tread2] ;
        ! the initial number of time steps n is increased until no cfl
        ! violations occure
        call Check_CFL( tread1, tread2, n, status )
        IF_NOTOK_RETURN(status=1)
#endif
        ! dynamic timestep in seconds:
        ndyn = nint( abs(rTotal(tread2-tread1,'sec')) / n )  ! sec

      end if  ! begin of new inteval

      ! setup meteo for dynamic step tdyn+[0,ndyn]
      if ( nhalf == 1 ) then
        ! time range of full (twice half) dynamic step:
        tr(1) = tdyn
        tr(2) = tdyn + revert * IncrDate( sec=ndyn )
        ! info ...
        call wrtgol( '- fullstep from ', tr(1), ' to ', tr(2) ); call goPr

#ifndef without_advection
        ! convert pu/pv to am/bm/cm, eventually time interpolated
        call Setup_MassFlow( tr, ndyn, status )
        IF_NOTOK_RETURN(status=1)
#endif
        ! setup (interpolate?) other meteo:
        call Meteo_Setup_Other( tr(1), tr(2), status )
       IF_NOTOK_RETURN(status=1)

        ! recalculate proces dependend fields if necessary
        call adj_Proces_Setup( tr(1), tr(2), status )
       IF_NOTOK_RETURN(status=1)

      end if

       ! end label:
      call GO_Timer_End( itim_timestep_init, status )
       IF_NOTOK_RETURN(status=1)

       ! * timestep run

       ! start label:
       call GO_Timer_Start( itim_timestep_run, status )
       IF_NOTOK_RETURN(status=1)


      ! full timestep operations:

       if ( ndyn > 0 .and. (mod(itau-itaui,ndyn) == 0) ) then
          ! reset the process status counters...
          region_status(1:nregions) = steps_region(1:nregions)
       end if

      ! operations on a half time step basis:

       if ( ndyn > 0 .and. (mod(itau-itaui,ndyn/2) == 0) ) then

        !! info ...
        !if ( myid == root ) then
        !  call tstamp(kmain,itau,'adj_tracer: Start processing main region')
        !end if

          itaur(:) = itau !synchronize time-count regions....

        ! time range of dynamic half=step:
        tr(1) = tdyn
        tr(2) = tdyn + revert * IncrDate( sec=nint(ndyn/2.0) )

        ! info ..
        !write (line,'("-- halfstep ",i0,"/2")') nhalf
        !call wrtgol( trim(line)//' from ', tr(1), ' to ', tr(2) ); call goPr

        call adj_Proces_Region( 1, tr, status )
          IF_NOTOK_RETURN(status=1)

          ! check times ...
          do region=2,nregions
             if ( itaur(region) /= itaur(1) ) then
                call escape_tm('tracer: exit of routine proces_region '//&
                     'with non-synchronized clocks')
             end if
          end do

       end if


      ! advance the model time with ndyn/2 seconds:
      call inctime( ndyn, status )
      IF_NOTOK_RETURN(status=1)
      ! update dynamic timestep:
      tdyn = tdyn + IncrDate( sec=revert*nint(ndyn/2.0) )

       ! end label:
       call GO_Timer_End( itim_timestep_run, status )
       IF_NOTOK_RETURN(status=1)

    end do   ! main loop

    t1 = SystemDate()
     call wrtgol( 'Wall time after     adjoint run =  ', t1 ); call goPr

    dt = t1 - t0
    call wrtgol( 'Wall time needed by adjoint run = ', dt ); call goPr

    ! *** done

    ! done processes
    call adj_Proces_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with zoom tools:
    call Zoom_Tools_Done( status )
    IF_NOTOK_RETURN(status=1)

    call exit_TM5(status)
    IF_NOTOK_RETURN(status=1)

    ! close meteo files etc;
    ! call this after 'exit_TM5' since mass array is saved:
    call Meteo_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ***

    ! end label:
    call goLabel()

    ! ok
    status = 0

  end subroutine adj_Tracer_Model


end module adj_Tracer
