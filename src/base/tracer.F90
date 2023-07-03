!#################################################################
!
!****   Global Atmospheric Tracer Model TM5, Version 1.0 , oct 2001
!
!       programmed by mh,  mpi HH,  1-oct-1991 (TM2 version)
!       programmed by Mike Botchev, Frank Dentener, Ad Jeuken,
!                   Patrick Berkvens and Maarten Krol
!
!       purpose
!       -------
!       solves the tracer continuitiy equation on an eulerian grid for
!       an arbitrary no. of tracers.
!       Allow zoom regions with higher resolution.
!       Perform advection, vertical transport, emissions,
!               chemistry, deposition, ....
!
!       interface
!       ---------
!       main program
!
!       method
!       ------
!       initialize the model and THEN cycle continuously through
!       the main model time loop.
!
!       externals
!       ---------
!
!       reference
!       ---------
!       see manual which does not exist.
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

module Tracer

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------

  private

  public :: Tracer_Init, Tracer_Done
  public :: Tracer_Model


  ! --- const ------------------------------------

  ! module name
  character(len=*), parameter  ::  mname = 'Tracer'


  ! --- types ------------------------------------

  ! timers:
  integer     ::  itim_tracer
  integer     ::  itim_timestep_init, itim_timestep_run, itim_timestep_done
  integer     ::  itim_model_init
  integer     ::  itim_read_meteo
  integer     ::  itim_dynam0, itim_coarsen_region, itim_coarsen_ambm, itim_check_cfl

  ! -- locals ------------------------------------

  logical     ::  restart_save_monthly, restart_save_weekly
  integer     ::  restart_save_weekday

contains


  ! ============================================================================


  subroutine Tracer_Init( status )

    use GO              , only : GO_Timer_Def, ReadRc
    use Meteo           , only : Meteo_Def_Timers
    use InitExit        , only : Init_TM5
    use global_data     , only : rcF

    ! --- in/out -------------------------------

    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Tracer_Init'

    ! --- begin --------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! define timers:
    call GO_Timer_Def( itim_tracer        , 'tracer model'   , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_model_init    , 'model init'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_timestep_init , 'timestep init'  , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_timestep_run  , 'timestep run'   , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_timestep_done , 'timestep done'  , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_read_meteo    , 'read meteo'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynam0        , 'dynam0'         , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_coarsen_region, 'coarsen_region' , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_coarsen_ambm  , 'coarsen_ambm'   , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_check_cfl     , 'check_cfl'      , status )
    IF_NOTOK_RETURN(status=1)

    ! initialize meteo module:
    call Meteo_Def_Timers( status )
    IF_NOTOK_RETURN(status=1)

    ! init model:
    !   paths, regions, operator splitting
    !write (gol,'("  init TM5 ...")'); call goPr
    call Init_TM5(status)
    IF_NOTOK_RETURN(status=1)

    ! Should we save restart files in the middle of a run?
    call ReadRc(rcF, 'restart.save.monthly', restart_save_monthly, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    call ReadRc(rcF, 'restart.save.weekly', restart_save_weekly, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    if (restart_save_weekly) then
        call ReadRc(rcF, 'restart.save.weekday', restart_save_weekday, status)
        IF_NOTOK_RETURN(status=1)
    end if

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Tracer_Init


  ! ***


  subroutine Tracer_Done( status )

    use InitExit  , only : Free_TM5

    ! --- in/out -------------------------------

    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Tracer_Done'

    ! --- begin --------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! clear model arrays etc:
    call Free_TM5(status)
    IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Tracer_Done


  ! ***


  subroutine Tracer_Model( status )

    use GO              , only : goLabel
    use GO              , only : GO_Timer_Start, GO_Timer_End
    use GO              , only : TDate, NewDate, TIncrDate, IncrDate
    use GO              , only : operator(-), operator(+), operator(*), operator(==), min
    use GO              , only : rTotal
    use GO              , only : wrtgol
    use GO              , only : SystemDate
    use Dims            , only : nregions
    use Dims            , only : region_status => status
    use Dims            , only : revert
    use Dims            , only : newsrun, newmonth, newday
    use Dims            , only : itau, itaui, itaue, itaur
    use Dims            , only : idate, idatei, idatee
    use Dims            , only : tref
    use Dims            , only : nsec_read
    use Dims            , only : ndyn
    use Dims            , only : kmain, ndyn_max
    use Meteo           , only : Meteo_Init, Meteo_Done
    use Meteo           , only : Meteo_Setup_Mass
    use AdvectM_CFL     , only : Setup_MassFlow
    use Meteo           , only : Meteo_Setup_Other
    use sources_sinks,    only : trace0
    use datetime,         only : inctime, tstamp
    use Zoom_Tools      , only : Zoom_Tools_Init, Zoom_Tools_Done
    use zoom_tools,       only : coarsen_region
    use advect,           only : dynam0
    use advect_tools,     only : coarsen_ambm
    use modelIntegration, only : Proces_Setup
    use modelIntegration, only : proces_region
    use advectm_cfl,      only : check_cfl
    use toolbox,          only : escape_tm
#ifndef without_chemistry
    use chemistry,        only : read_chemistry_fields
#endif

    use InitExit        , only : Start_TM5, Exit_TM5
    use Dims            , only : revert
    use ModelIntegration, only : Proces_Init, Proces_Done
    use ParTools,         only : Par_Init, Par_Done
    use ParTools,         only : Par_Barrier, Par_Check_Domain
    use ParTools,         only : myid, root
    use os_specs,         only : WRITE_STR_LEN
    use Go,               only : day_of_week
    use Restart,          only : Restart_Save

    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Tracer_Model'

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
    logical             ::  restart_save_now

    ! --- begin ---------------------------------------

    ! start label:
    call goLabel( timer=itim_tracer )

    !---------------!
    ! model startup !
    !---------------!

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! start label:
    call GO_Timer_Start( itim_model_init, status )
    IF_NOTOK_RETURN(status=1)

    ! initialize zoom tools:
    call Zoom_Tools_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! setup meteo input
    call Meteo_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init processes
    call Proces_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! tracer for forward run:
    revert = 1

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

    t0 = SystemDate()
    call wrtgol( 'Wall time before    forward run = ', t0 ); call goPr

    ! time loop
    do
       if ( (revert*itau) >= (revert*itaue) ) exit

       ! * timestep init

       ! save model status alternatingly on units 1 and 2
       ! note the staggered (nread/2) shift..
       ! IF ((nwrite>0).and.(mod(itau,nwrite)==0).and.myid==root) call  bisave

       ! start label:
       call GO_Timer_Start( itim_timestep_init, status )
       IF_NOTOK_RETURN(status=1)

      ! next half step
      nhalf = modulo(nhalf,2) + 1

      ! info at first half:
!      if ( nhalf == 1 ) then
!        write (gol,'("")'); call goPr
!        call wrtgol( '>>> timestep from ', tdyn ); call goPr
!        write (gol,'("")'); call goPr
!        call wrtgol( 'current data valid for ', tread1, ' to ', tread2 ); call goPr
!      end if

       ! first set-up of chemistry related work...repeated if newmonth
       if ( newmonth .and. (.not. newsrun) ) then
          call trace0( status )
          IF_NOTOK_RETURN(status=1)
       end if

#ifndef without_chemistry
       ! read chemistry-related fields, once per day
       if (newday .or. newsrun) then
          call read_chemistry_fields(1, status)
          IF_NOTOK_RETURN(status=1)
       end if
#endif

      ! reached end of time interval for which meteo is valid ? then setup new meteo:
      if ( tdyn == tread2 ) then

        ! setup meteo data for next interval;
        ! nread is the length (in seconds) of the interval in which
        ! surface pressure is interpolated (and mass fluxes are constant)
        tread1 = tdyn
        tread2 = min( tdyn + IncrDate(sec=nsec_read), tend )

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

      end if

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
        call Proces_Setup( tr(1), tr(2), status )
          IF_NOTOK_RETURN(status=1)

        ! should we save restart file?
        restart_save_now = .false.
        if (.not. newsrun) then
            if (restart_save_monthly .and. newmonth) restart_save_now = .true.
            if (restart_save_weekly .and. newday .and. day_of_week(tr(1)) == restart_save_weekday) restart_save_now = .true.
        end if
        if (restart_save_now) then
            call Restart_Save(status)
            IF_NOTOK_RETURN(status=1)
        end if

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
          region_status(1:nregions) = 0
       end if

       ! operations on a half time step basis:

       if ( ndyn > 0 .and. (mod(itau-itaui,ndyn/2) == 0) ) then

        !! info ...
        !if ( myid == root ) then
        !   call tstamp(kmain,itau,'tracer: Start processing main region')
        !end if

          itaur(:) = itau !synchronize time-count regions....
          call Par_Barrier

        ! time range of dynamic half=step:
        tr(1) = tdyn
        tr(2) = tdyn + IncrDate( sec=revert*nint(ndyn/2.0) )

        ! info ..
        !write (line,'("-- halfstep ",i0,"/2")') nhalf
        !call wrtgol( trim(line)//' from ', tr(1), ' to ', tr(2) ); call goPr

          !
          ! start recursive process for the main region = 1
          ! Within the advection, all the other processes are incorporated...
          ! Also the regions 'below' are called recursively....
          !
          !write (*,'("WARNING - skip processes ...")')
        call Proces_Region( 1, tr, status )
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

       ! * timestep done

       ! start label:
       call GO_Timer_Start( itim_timestep_done, status )
       IF_NOTOK_RETURN(status=1)

!       if (mod(itau,ndyn_max) == 0) then
!          call user_output_mean( status )
!          IF_NOTOK_RETURN(status=1)
!       end if

       ! end label:
       call GO_Timer_End( itim_timestep_done, status )
       IF_NOTOK_RETURN(status=1)

    end do   ! main loop

    t1 = SystemDate()
    call wrtgol( 'Wall time after     forward run = ', t1 ); call goPr

    dt = t1 - t0
    call wrtgol( 'Wall time needed by forward run = ', dt ); call goPr

    ! *** done

    ! done processes
    call Proces_Done( status )
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

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! end label:
    call goLabel()

    ! ok
    status = 0

  end subroutine Tracer_Model


end module Tracer
