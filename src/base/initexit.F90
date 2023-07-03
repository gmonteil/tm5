!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module initexit

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out ------------------------------

  private

  public :: TM5_Arguments
  public :: init_TM5
  public :: start_TM5
  public :: free_TM5
  public :: exit_TM5
  public :: TM5_OkFile

  character(len=*), parameter         :: mname = 'InitExit'

contains


  ! ===================================================================
  ! ===
  ! === arguments
  ! ===
  ! ===================================================================


  subroutine TM5_Arguments( status )

    use GO         , only : goArgCount, goGetArg
    use global_data, only : rcfile
    use partools   , only : myid, root, Par_Broadcast_Status, Par_Broadcast
    use os_specs   , only : LONG_STR_LEN

    ! --- in/out ----------------------------------

    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/TM5_Arguments'

    ! --- local -----------------------------------

    integer               ::  narg
    integer               ::  iarg
    character(len=LONG_STR_LEN)   ::  line

    ! --- begin -----------------------------------

    ! on root only, since some mpirun version do not parse
    ! all arguments to each executable:

    ! number of arguments:
    if (myid==root) call goArgCount( narg, status )
    call Par_Broadcast_Status(status,root)
    IF_NOTOK_RETURN(status=1)
    call Par_Broadcast( narg, root, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( narg == 0 ) then
      write (gol,'("no arguments found ...")') narg; call goErr
      TRACEBACK; status=1; return
    end if

    ! defaults:
    rcfile = 'None'

    ! loop over arguments:
    iarg = 0
    do
      ! next:
      iarg = iarg + 1
      ! get argument:
      if (myid==root) call goGetArg( iarg, line, status )
      call Par_Broadcast_Status(status,root)
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( line, root, status )
      IF_NOTOK_RETURN(status=1)
      ! specials ...
      select case ( trim(line) )
        ! arguments added by MPICH/mpirun :
        case ( '-p4pg', '-p4wd' )
          ! skip next argument:
          iarg = iarg + 1
        ! other ...
        case default
          ! not filled yet ?
          if ( trim(rcfile) == 'None' ) then
            rcfile = trim(line)
          else
            write (gol,'("unsupported argument : ",a)') trim(line); call goErr
            TRACEBACK; status=1; return
          end if
      end select
      ! last one is processed now ?
      if ( iarg == narg ) exit
    end do

    ! ok
    status = 0

  end subroutine TM5_Arguments


  ! ***


  subroutine TM5_Print_Usage( status )

    ! --- in/out ---------------------------------

    integer, intent(out)        ::  status

    ! --- begin ----------------------------------

    ! display usage line:
    write (*,'("Usage: tm5.x <rcfile>")')

    ! ok
    status = 0

  end subroutine TM5_Print_Usage



  ! ===================================================================
  ! ===
  ! === init/free of (multiple) TM5 runs
  ! ===
  ! ===================================================================


  !
  ! *** initialization of (multiple) TM5 runs
  !

  subroutine init_TM5(status)


    use GO,            only : ReadRc
    use dims,          only : okdebug, okdebug_tmm
    use dims,          only : nregions
    use dims,          only : istart
    use dims,          only : nsec_read
    use global_data,   only : rcF
    use global_data,   only : inputdir, outdir
    use global_data,   only : declare_fields
    use global_data,   only : region_dat
    use user_output,   only : user_output_init
    use TM5_Geometry , only : TM5_Geometry_Init
    use advectm_cfl,   only : init_cfl
    use ParTools     , only : Par_Init
    !use Sources_Sinks, only : Sources_Sinks_Init
    !use Emission     , only : Emission_Init
    use dims,          only : czeta, czetak
    use restart,       only : Restart_Init

    !__IO___________________________________________________________________

    integer, INTENT(OUT)             :: status

    !__CONST________________________________________________________________

    character(len=*), parameter      :: rname = mname//', init_TM5'

    !__LOCAL_VARIABLES______________________________________________________

    integer             :: region

    !__START_SUBROUTINE______________________________________________________

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    !===============================
    ! init parallelisation
    !===============================
    call Par_Init( status )
    IF_NOTOK_RETURN(status=1)

    !===============================
    ! set default control parameters
    !===============================

    call default( rcF, status )


    !=============================
    ! read parameters from rcfile
    !=============================

!    call ReadRc(rcF, 'czetak', czetak, status)
!    call ReadRc(rcF, 'czeta_blh', czeta_blh, status)

    call ReadRc( rcF, 'input.dir', inputdir, status )
    IF_NOTOK_RETURN(status=1)

    call ReadRc( rcF, 'output.dir', outdir, status )
    IF_NOTOK_RETURN(status=1)

    call ReadRc(rcF, 'istart', istart, status)
    IF_NOTOK_RETURN(status=1)

    call ReadRc( rcF, 'okdebug', okdebug, status, default=.false.)
    IF_ERROR_RETURN(status=1)

    call ReadRc( rcF, 'okdebug.tmm', okdebug_tmm, status, default=okdebug)
    IF_ERROR_RETURN(status=1)

    call ReadRc( rcF, 'timestep.read', nsec_read, status)
    IF_NOTOK_RETURN(status=1)

    !=================
    ! allocate memory
    !=================

    ! allocate memory
    call declare_fields( status )
    IF_NOTOK_RETURN(status=1)

!PB    call init_cfl          !initialise CFL preventing routine


    !=================
    ! geometry
    !=================

    write (gol,'(a,":   init model geometry ...")') rname; call goPr

    write (gol,'(a,":     determine children etc ...")') rname; call goPr
    call determine_children_etc

    write (gol,'(a,":     blank zoom regions ...")') rname; call goPr
    do region = 1, nregions
       call blank_zoom_region(region,region_dat(region))
    enddo
    !call print_zoom(region_dat)

    ! init geometry:
    call TM5_Geometry_Init( status )
    IF_NOTOK_RETURN(status=1)


    !=================
    ! dynamics
    !=================

    ! Scaling factor for convection
    call ReadRc(rcF, 'scale.convec', czeta, status, default=1.0)
    IF_ERROR_RETURN(status=1)

    ! Scaling factor for vertical diffusion
    call ReadRc(rcF, 'scale.vdiff', czetak, status, default=1.0)
    IF_ERROR_RETURN(status=1)

    ! Setting czeta and czetak to <1 would slow down vertical transport

    ! init operator splitting:
    write (gol,'(a,":   split order ...")') rname; call goPr
    call define_splitorderzoom

    ! init time window:
    write (gol,'(a,":   init time window ...")') rname; call goPr
    call Init_Time_Window( status )
    IF_NOTOK_RETURN(status=1)

    ! init time steps:
    write (gol,'(a,":   init time steps ...")') rname; call goPr
    call Init_Time_Steps( status )
    IF_NOTOK_RETURN(status=1)

    call Restart_Init(status)
    IF_NOTOK_RETURN(status=1)

    !! init emission arrays:
    !write (gol,'(a,":   init emissions ...")') rname; call goPr
    !call Emission_Init( status )
    !IF_NOTOK_RETURN(status=1)

    !! init sources/sinks:
    !write (gol,'(a,":   init sources/sinks ...")') rname; call goPr
    !call Sources_Sinks_Init( status )
    !IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok:
    status = 0

  end subroutine init_TM5


  !
  ! *** set default values of control variables
  !

  subroutine default( rcF, status )

    use GO         , only : TrcFile, ReadRc
    use dims       , only : icalendo
    use dims       , only : ndyn
    use dims       , only : nstep
    use dims       , only : istart
    use dims       , only : itaui, idatei
    use dims       , only : itaue, idatee
    use dims       , only : itaut, idatet
    use dims       , only : newsrun
    use dims       , only : czeta, czetak
    use dims       , only : cdebug, okdebug
    use dims       , only : limits
    use dims       , only : revert
    use dims       , only : nwrite, ndiag, ntrans, ninst, ndiagp1, ndiagp2
    use dims       , only : ncheck
    use datetime,    only : tau2date
    use toolbox,     only : escape_tm
    use var4d,       only : nasim

    !__IO___________________________________________________________________

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    !__CONST________________________________________________________________

    character(len=*), parameter      :: rname = mname//'/default'

    !__LOCAL_VARIABLES______________________________________________________

    integer :: i,k,n


    !__START_SUBROUTINE______________________________________________________


    icalendo=2  ! (real calendar, iyear0=1980)

    ndyn=5400

    ! window for observation representation:
    !nasim=3*3600

    istart=10   ! invalid value; specify in rcfile

    itaui=0
    newsrun=.true.
    call tau2date(itaui,idatei)
    if ( mod(idatei(4),3) /= 0 ) &
         call escape_tm('default: GMT start time should be multiple of 3')
    itaue=itaui
    call tau2date(itaue,idatee)
    itaut=0
    call tau2date(itaut,idatet)

    czeta=1.0    ! full convection

!    czetak = 1.0     ! no scaling of vertical diffusion
!    czeta_blh = 1.0  ! no scaling of TM5 blh

    !    limits=.true.
    limits=.false.

    nstep=0    ! initialize counter of model steps

    cdebug=.false.
    okdebug=.false.
    revert = 1

    !===================
    ! currently not used
    !===================
    nwrite=-1
    ndiag=0
    ntrans=0
    ninst=0
    ndiagp1=-2
    ndiagp2=-2
    czetak=1.     ! scaling factor for vertical diffusion
    ncheck=6      ! checking interval

    ! ok
    status = 0

  end subroutine default


  ! ***


  subroutine Init_Time_Window( status )

    use GO                 , only : ReadRc
    use GO                 , only : Get
    use GO                 , only : goTranslate
    use GO                 , only : Time_Window_Init
    use datetime           , only : time_window
    use global_data        , only : rcF

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'Init_Time_Window'

    ! --- local ----------------------------------

    character(len=32)                   ::  stime
    integer                             ::  time6_1(6)
    integer                             ::  time6_2(6)

    ! --- begin ----------------------------------
    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! * start time
    ! read from rcfile setting, yyyy-mm-dd hh:mn:sc format:
    call ReadRc( rcF, 'jobstep.timerange.start' , stime, status )
    IF_NOTOK_RETURN(status=1)
    ! replace '-', ':', etc by white space:
    call goTranslate( stime, '/-:', ' ', status )
    IF_NOTOK_RETURN(status=1)
    ! read 6 numbers:
    read (stime,*,iostat=status) time6_1
    if ( status /= 0 ) then
      write (gol,'("could not read start time from : ",a)') trim(stime); call goErr
      TRACEBACK; status=1; return
    end if

    ! * end time:
    ! read from rcfile setting, yyyy-mm-dd hh:mn:sc format:
    call ReadRc( rcF, 'jobstep.timerange.end' , stime, status )
    IF_NOTOK_RETURN(status=1)
    ! replace '-', ':', etc by white space:
    call goTranslate( stime, '/-:', ' ', status )
    IF_NOTOK_RETURN(status=1)
    ! read 6 numbers:
    read (stime,*,iostat=status) time6_2
    if ( status /= 0 ) then
      write (gol,'("could not read end time from : ",a)') trim(stime); call goErr
      TRACEBACK; status=1; return
    end if

    ! * fill time window:
    call Time_Window_Init( time_window, time6_1, time6_2, status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Init_Time_Window


  ! ***


  subroutine Init_Time_Steps( status )

    use GO,            only : ReadRc
    use dims,          only : ndyn, ndyn_max
    use global_data,   only : rcF

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'Init_Time_Steps'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! read dynamic timestep:
    call ReadRc(rcF, 'ndyn', ndyn, status)
    IF_NOTOK_RETURN(status=1)

    ! copy:
    ndyn_max = ndyn

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Init_Time_Steps


  ! ***


  ! final termination of (multiple) TM5 runs

  subroutine free_TM5( status )

    use ParTools,         only : Par_Done
    use Sources_Sinks   , only : Sources_Sinks_Done
    use Emission        , only : Emission_Done
    use TM5_Geometry    , only : TM5_Geometry_Done

    !__IO____________________________________________________________________

    integer, INTENT(OUT)             :: status

    !__CONST_________________________________________________________________

    character(len=*), parameter      :: rname = mname//'/free_TM5'

    !__LOCAL_VARIABLES_______________________________________________________

    !__START_SUBROUTINE______________________________________________________

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! done with emission arrays:
    call Emission_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with sources/sinks:
    call Sources_Sinks_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with geometry:
    call TM5_Geometry_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done parallelisation
    call Par_Done( status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr


  end subroutine free_TM5


  ! ***


  ! Write a dummy text file to indicate a normal end of the run.
  ! Necessary since some routines simply break with a 'stop' statement,
  ! which does not lead to model exit with error status.

  subroutine TM5_OkFile( status )

    use GO, only : goGetFU, ReadRc
    use global_data,only : rcF
    use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

    ! --- in/out ---------------------------------

    integer, intent(out)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5_OkFile'

    ! --- local ----------------------------------

    integer   ::  fu
    character(len=MAX_FILENAME_LEN) :: ok_file_name
    character(len=MAX_RCKEY_LEN)    :: output_dir

    ! --- begin ----------------------------------

    ! get free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! get output dir
    call ReadRc(rcF, 'output.dir', output_dir, status)
    IF_NOTOK_RETURN(status=1)

    ! write ok file:
    write(ok_file_name, '(a, "/tm5.ok")') trim(output_dir)
    open( unit=fu, form='formatted', status='unknown', file=trim(ok_file_name) )
    write (fu,'("Program terminated normally")')
    close( fu )

    ! ok
    status = 0

  end subroutine TM5_OkFile


  ! ===================================================================
  ! ===
  ! === start/exit of single TM5 runs
  ! ===
  ! ===================================================================


  !
  ! *** initialization of particular TM5 run
  !


  subroutine start_TM5( tread1, tread2, status)

    use GO                 , only : ReadRc
    use GO                 , only : Get
    use GO                 , only : TDate, NewDate, IncrDate
    use GO                 , only : operator(+), operator(-), operator(*), rTotal

    use dims,                only : revert
    use dims,                only : idate, idatei, idatee
    use dims,                only : itau, itaui, itaue
    use dims,                only : iyear0, julday0
    use dims,                only : ndyn, ndyn_max
    use dims,                only : nsec_read
    use dims,                only : newyr, newday, newmonth, newsrun
    use dims,                only : nstep, nstep0
    use dims,                only : nregions, region_name
    use dims,                only : tref
    use dims,                only : adv_scheme
    use dims,                only : cpu0
    use dims,                only : kmain
    use dims,                only : im, jm

    use global_data,         only : rcF
    use global_data,         only : region_dat
    use global_data,         only : declare_fields
    use global_data,         only : inputdir, outdir

    use global_data        , only : mass_dat

    use sources_sinks,       only : trace0
    use adj_sources_sinks,   only : adj_trace0

!#ifndef without_chemistry
    !use Chemistry          , only : Chemistry_Init
!#endif

    use Meteo              , only : Meteo_Setup_Mass
    use AdvectM_CFL        , only : Setup_MassFlow
    use Meteo              , only : Meteo_Setup_Other

    use geometry,            only : geomtryh, calc_dxy
    use user_output,         only : user_output_init
    use adj_user_output,     only : adj_user_output_init
    use datetime,            only : date2tau, tau2date, julday, tstamp
    use zoom_tools,          only : coarsen_region
    use advect,              only : dynam0
    use advect_tools,        only : advect_m, m2phlb, m2phlb1, coarsen_ambm
    use advectm_cfl,         only : init_cfl, check_cfl
    use parTools,            only : myid

    use datetime                   , only : time_window
    use var4d,                       only : runid, run_mode, iter

    use Var4D_IO_Mass              , only : restore_masses !, restore_pressure

!    use Emission                   , only : emissions
!    use Emission_Data              , only : ref_emissions, ref_emissions_apri
!    use Emission_Data              , only : n_cat

    !__IO___________________________________________________________________

    type(TDate), intent(out)    ::  tread1, tread2
    integer, intent(out)        ::  status

    !__CONST_______________________________________________________________

    character(len=*), parameter      :: rname = mname//'/start_TM5'

    !__LOCAL_VARIABLES______________________________________________________

    integer             :: region
    integer             :: n
    integer             :: i_cat
    integer             :: i_period
    integer             :: i,j
    type(TDate)         ::  tdyn, tr(2)
    integer             ::  nsec

    !__START_SUBROUTINE______________________________________________________


    write (gol,'(a," : entering")') trim(rname) ; call goPr

    !============
    ! set runid
    !============

!    call set_runid

    !============
    ! start time
    !============

    ! * 6-valued start and end arrays:
    select case ( revert )
      case ( 1 )
        call Get( time_window%t1, time6=idatei )
        call Get( time_window%t2, time6=idatee )
      case ( -1 )
        call Get( time_window%t2, time6=idatei )
        call Get( time_window%t1, time6=idatee )
      case default
        write (gol,'("unsupported revert : ",i6)') revert; call goErr
        TRACEBACK; status=1; return
    end select

    !iyear0 = idatei(1)-1     ! itau runs from beginning of PREVIOUS year
    iyear0 = time_window%t1%year-1  ! ensure that itau remains positive also in adjoint run

    julday0=julday(1,1,iyear0)

    call date2tau(idatei,itaui)
    call date2tau(idatee,itaue)

    itau=itaui
    call tau2date(itau,idate)

    ! set time flags
    newyr=.true.
    newday=.true.
    newmonth=.true.
    newsrun=.true.

    ! step counter
    nstep0=0
    nstep=0

    ! *

    ! initialise CFL loop:

    call init_cfl( status )
    IF_NOTOK_RETURN(status=1)

    !.........................................................
    ! new
    !.........................................................

    write (gol,'(a,":   initial meteo ...")') rname; call goPr

    ! safety ...
    if ( any( idate(4:6) /= 0 ) ) then
      write (gol,'("current implementation expects start at 00:00")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! setup from start time to end of interval  [k*nsec_read,(k+1)*nsec_read] ;
    ! step towards end of first interval from or around 00:00:
    nsec = (floor(idate(4)*3600.0/nsec_read)+1)*nsec_read
    ! set times:
    tread1 = NewDate(time6=idate)
    tread2 = NewDate(time6=(/idate(1:3),0,0,0/)) + revert * IncrDate(sec=nsec)

    !! setup for interval [-staggerd,nreadX]:
    !tread1 = NewDate(time6=idate) - IncrDate(sec=staggered)
    !tread2 = tread1 + IncrDate(sec=nreadX)

    ! setup pressure and mass fields
    !  o do not check pressure implied by advection
    !  o for adjoint, restore airmass and compute phlb and sp from this
    call Meteo_Setup_Mass( tread1, tread2, status, &
                             isfirst=.true., &
                             restore_airmass=(revert==-1) )
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

    ! * setup meteo for dynamic step tdyn+[0,ndyn]

    ! current time (begin of dynamics step)
    tdyn = NewDate( time6=idate  )

    ! time range of full (two half) dynamic step:
    tr(1) = tdyn
    tr(2) = tdyn + revert * IncrDate( sec=ndyn )

    ! convert pu/pv to am/bm/cm, eventually time interpolated
    call Setup_MassFlow( tr, ndyn, status )
    IF_NOTOK_RETURN(status=1)

    ! setup (interpolate?) other meteo:
    call Meteo_Setup_Other( tr(1), tr(2), status, isfirst=.true. )
      IF_NOTOK_RETURN(status=1)

    !=========================
    ! Concentrations
    !=========================

    call Start_TM5_Concentrations( status )
    IF_NOTOK_RETURN(status=1)


    !=========================
    ! initialize chemistry (OH), surface fields, datetime stuff etc.
    ! [Must be before reading emissions to have sec_year available.]
    !=========================

!#ifndef without_chemistry
    !! start chemistry:
    !call Chemistry_Init( status )
    !IF_NOTOK_RETURN(status=1)
!#endif

    select case ( revert )
      case ( 1 )
        call trace0( status )
        IF_NOTOK_RETURN(status=1)
      case ( -1 )
        call adj_trace0
      case default
        write (gol,'("unsupported revert : ",i6)') revert; call goErr
        TRACEBACK; status=1; return
    end select

    !=========================
    ! Emissions
    !=========================

    call Start_TM5_Emissions( status )
    IF_NOTOK_RETURN(status=1)

    !===================================
    ! initialise user-specified output
    !===================================

    if ( revert == 1 ) then
       call user_output_init(status)
       IF_NOTOK_RETURN(status=1)
    else
       call adj_user_output_init(status)
       IF_NOTOK_RETURN(status=1)
    end if


    !===================================
    ! original timing
    !===================================

    ! cputim
    call cputim( cpu0 )

    !===================================
    ! closing message
    !===================================

    write (gol,'("=======================")'); call goPr
    write (gol,'("+++++++++++++++++++++++")'); call goPr
    write (gol,'("model startup complete ")'); call goPr
    call tstamp(kmain,itau,'idate')
    write (gol,'("+++++++++++++++++++++++")'); call goPr
    write (gol,'("=======================")'); call goPr

    ! ok
    status = 0

    write (gol,'(a," : done")') trim(rname) ; call goPr

  end subroutine start_TM5


  ! ***


  subroutine Start_TM5_Emissions( status )

    use GO                , only : ReadRc
    use dims,               only : revert
    use global_data,        only : rcF
    use Emission,           only : Emission_Fwd_Setup
    use Emission,           only : Emission_Adj_Setup

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'Start_TM5_Emissions'

    ! --- local ----------------------------------


    ! --- begin ----------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    select case ( revert )

      ! forward run
      case ( 1 )

          write (gol,'("start: getting emissions for whole assimilation period")'); call goPr
          call Emission_Fwd_Setup( status )
          IF_NOTOK_RETURN(status=1)


      ! adjoint run
      case ( -1 )

        ! set all adj_emissions = 0.0
        call Emission_Adj_Setup( status )
        IF_NOTOK_RETURN(status=1)

      ! error ...
      case default

        write (gol,'("unsupported revert : ",i6)') revert; call goErr
        TRACEBACK; status=1; return

    end select

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Start_TM5_Emissions


  ! ***


  subroutine Start_TM5_Concentrations( status )

    use GO                 , only : ReadRc
    use dims,                only : revert
    use dims,                only : istart
    use dims,                only : nregions, region_name
    use dims,                only : adv_scheme
    use global_data,         only : rcF
    use global_data,         only : mass_dat
    use chem_param         , only : fscale
    use zoom_tools,          only : update_parent
    use sources_sinks,       only : trace1
    use adj_sources_sinks,   only : adj_trace1
    use io_save,             only : readnetcdf, save_filename
    use restart,             only : Restart_Read
    use os_specs,            only : MAX_FILENAME_LEN

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Start_TM5_Concentrations'

    ! --- local ----------------------------------

    integer                         :: region
    character(len=MAX_FILENAME_LEN) :: fname

    ! --- begin ----------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    if ( revert == 1 ) then

      ! JFM: need to read initial fields always, because of slopes.
      !   Tracer masses themselves may be reset in 4d-var / grad test.

      write(*,*) 'start: getting initial concentration field'

      !==============================
      ! get initial tracer fields...
      !  only in first iteration or in case of gradient test / 4dvar
      !  when initial concentrations are not in the state vector !!
      !==============================

      select case(istart)

        !=======================================
        ! initial tracer fields are set to zero
        !=======================================
        case(1)

          do region = 1,nregions

             mass_dat(region)%rm_t = 0.0

             if ( adv_scheme == 'slope' ) then

                mass_dat(region)%rxm_t=0.
                mass_dat(region)%rym_t=0.
                mass_dat(region)%rzm_t=0.
             end if
          end do


        !======================================
        ! initial tracer fields are calculated
        ! in trace1
        !======================================
        case(2)

          call trace1( status )
          IF_NOTOK_RETURN(status=1)

        !====================================
        ! initial tracer fields are obtained
        ! from a save file of a previous model run
        ! m is also obtained!
        !====================================
        case(3)

          call ReadRc(rcF, 'start.3.filename', fname, status)
             IF_NOTOK_RETURN(status=1)
          call readnetcdf(trim(fname), status)
          IF_NOTOK_RETURN(status=1)

        case(31)
          ! read in the save file from a previous run
          call save_filename(fname, status)
          IF_NOTOK_RETURN(status=1)
          call readnetcdf(trim(fname), status)
          IF_NOTOK_RETURN(status=1)

        case(33)

          !call Restart_Read( status, tracer_mass=.true., tendencies=.true. , air_mass=.true.) ! For cy3/chemistry code
          call Restart_Read( status, surface_pressure=.true., pressure=.true., tracer_mass=.true., air_mass=.true.)
          IF_NOTOK_RETURN(status=1)

        !====================================
        case default
        !====================================

          write (gol,'("start: illegal value of parameter istart : ",i6)') istart; call goErr
          TRACEBACK; status=1; return

       end select


       do region = nregions, 2, -1
          if ( revert == 1 ) call update_parent( region )
       end do

       !       end if

    else

       call adj_trace1  ! set all rm = 0.0

    end if

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Start_TM5_Concentrations


  !
  ! *** final termination of single TM5 run
  !

  subroutine Exit_TM5(status)

    use GO,           only: ReadRc
    use GO,           only: pathsep

    use dims,           only: nregions
    use dims,           only: nxi
    use dims,           only: xi
    use dims,           only: nloop_max
    use dims,           only: kmain
    use dims,           only: idate, itau
    use dims,           only: cpu0
    use dims,           only: nstep
    use dims,           only: kdebug
    use dims,           only: revert

    use datetime,       only : tstamp

    !use global_data,    only : free_fields

    use io_save,        only : savenetcdf, save_filename

    use sources_sinks,  only : trace_end

    use user_output,    only : user_output_done
    use adj_user_output, only : adj_user_output_done

    use advectm_cfl,    only : done_cfl

    use toolbox,        only : escape_tm

    use global_data,  only: rcF
    use global_data,  only: outdir
   use global_data,    only: mass_dat

    use Var4D_IO_Mass  , only : save_masses ! save_pressure
    use misctools,       only : check_dir
    use restart,         only : Restart_Save
    use os_specs,        only : MAX_FILENAME_LEN

!#ifndef without_chemistry
    !use Chemistry       , only : Chemistry_Done
!#endif

    !__IO___________________________________________________________________

    integer, INTENT(OUT)              :: status


    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//'/Exit_TM5'

    integer :: region
    real    :: cpu3, xx, gg
    integer :: istat

    character(len=MAX_FILENAME_LEN) :: save_subdir, fname

    !__START_SUBROUTINE______________________________________________________

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    call save_filename(fname, status)
    IF_NOTOK_RETURN(status=1)
    call check_dir(trim(fname))

    if(revert == 1) then

       ! savefile is written into sub-directory of output:
       !call ReadRc(rcF, 'save.output.subdir', save_subdir, status)
       !IF_NOTOK_RETURN(status=1)
       ! base name of save file; region name and extension are added by 'savehdf':
!       write( fname, '(4a,"save_",i4,4i2.2)' ) trim(outdir), pathsep, &
!                                       trim(save_subdir), pathsep, idate(1:5)
!       ! put out:
!       call savehdf( 'successful completion of run', trim(fname) )
       ! Also write a netcdf file
       !write( fname, '(4a,"save_",i4,4i2.2,".nc4")' ) trim(outdir), pathsep, &
                                       !trim(save_subdir), pathsep, idate(1:5)

       call savenetcdf(trim(fname), status)
       IF_NOTOK_RETURN(status=1)

       call Restart_Save(status, extra=.false., isfirst=.false.)
       IF_NOTOK_RETURN(status=1)
       ! JFM and MK        Feb 2007
       ! In zoom mode the original air masses are not reproducable from
       ! the surface pressure (problem related to routine coarsen_region).
       ! Therefore air masses are now saved instead of surface pressure.
       !
       call save_masses( status )
       IF_NOTOK_RETURN(status=1)

       ! complete user-specified output
       call user_output_done( 'successful completion of run', status )

       ! finish up chemistry/print budget, output shored-lived compounds
       call trace_end('successful completion of run','no-save-files-written-from-here',status)
       IF_NOTOK_RETURN(status=1)

    else

       ! savefile is written into sub-directory of output:
       !call ReadRc(rcF, 'save.output.subdir', save_subdir, status)
       !IF_NOTOK_RETURN(status=1)
       ! base name of save file; region name and extension are added by 'savehdf':
       !write( fname, '(4a,"adj_save_",i4,4i2.2,".nc4")' ) trim(outdir), pathsep, &
                                       !trim(save_subdir), pathsep, idate(1:5)
       call savenetcdf(trim(fname), status)
       IF_NOTOK_RETURN(status=1)

       call adj_user_output_done(status)
       IF_NOTOK_RETURN(status=1)

    endif

    ! CMKadj: leave memory allocated???
    !    call free_fields     !free memory of this run

    write (gol,'("")'); call goPr
    write (gol,'("exitus: CFL info:")'); call goPr
    do region=1, nregions
       write (gol,'(a,3i4,f10.4)') 'exitus: x: region, nxi, nloop_max, xi',  &
            region, nxi(region,1), &
            nloop_max(region,1), xi(region,1); call goPr
       write (gol,'(a,3i4,f10.4)') 'exitus: y: region, nxi, nloop_max, xi',  &
            region, nxi(region,2), &
            nloop_max(region,2), xi(region,2); call goPr
       write (gol,'(a,3i4,f10.4)') 'exitus: z: region, nxi, nloop_max, xi',  &
            region, nxi(region,3), &
            nloop_max(region,3), xi(region,3); call goPr
    end do

    call done_cfl

!#ifndef without_chemistry
    !! finish chemistry:
    !call Chemistry_Done( status )
    !IF_NOTOK_RETURN(status=1)
!#endif

    gol = ''; call goPr
    gol = repeat('-',80); call goPr
    gol = repeat('-',80); call goPr
    gol = ''; call goPr

    write (gol,'("exitus: program has terminated normally.")'); call goPr

    call tstamp(kmain,itau,'exitus: final time')
    call cputim(cpu3)
    cpu3 = cpu3-cpu0
    write (gol,'(" exitus: no of timesteps : ",i16)') nstep; call goPr
    write (gol,'(" exitus: cpu (s)         : ",f16.2)') cpu3; call goPr
    write (gol,'(" exitus: cpu/step        : ",f16.8)') cpu3/nstep; call goPr

    gol = ''; call goPr
    gol = repeat('-',80); call goPr
    gol = repeat('-',80); call goPr
    gol = ''; call goPr

    close(kdebug)

    write (gol,'(a," : done")') trim(rname) ; call goPr

  end subroutine exit_TM5


     !===========================================================================================================
     !===========================================================================================================

     subroutine cputim(time )
       !
       ! ***********************************************************************
       ! cpu time
       ! ***********************************************************************
       !
       implicit none

       real,intent(out) :: time

       integer :: clock,clockrate

       call system_clock(clock,clockrate)
       time=clock*(1.0/clockrate)
       !time=second()

     end subroutine cputim

     !===========================================================================================================
     !===========================================================================================================


     ! the following routines are called at start-up to calculate some
     ! useful information for the model, like array adresses and
     ! the lay-out of the defined regions.


     subroutine determine_children_etc
       !
       ! determines children and ibeg,iend,jbeg,jend,lbeg,lend
       ! the subroutine is normally called once per run of TM5
       ! written by patrick berkvens and mike botchev, march-june 1999
       !
       use dims
       implicit none
       ! local
       integer :: icells, jcells, lcells, jj, region
       integer :: xref_, yref_, zref_

       if ( okdebug ) print *,'determine_children_etc: '

       children = 0
       do region=2,nregions
          ! increase children's counter of the parent:
          children(parent(region),0) = children(parent(region),0) + 1
          jj = children(parent(region),0)
          ! this child is region:
          children(parent(region),jj) = region
       end do
       ! for example:
       !children(0,0) = 1  ! means: there is one child in the region 0
       !children(0,1) = 1  ! means: this child is region 1
       !children(1,0) = 1  ! means: there is one child in region 1
       !children(1,1) = 2  ! means: this child is region 1
       !children(2,0) = 0  ! no children in region 2

       do region=2,nregions
          ibeg(region) = nint((xbeg(region)-xbeg(parent(region)))    &
               /         (dx/xref(parent(region)))) + 1
          jbeg(region) = nint((ybeg(region)-ybeg(parent(region)))    &
               /         (dy/yref(parent(region)))) + 1
          lbeg(region) = nint((zbeg(region)-zbeg(parent(region)))    &
               /         (dz/zref(parent(region)))) + 1

          iend(region) = nint((xend(region)-xbeg(parent(region)))  &
               /         (dx/xref(parent(region))))
          jend(region) = nint((yend(region)-ybeg(parent(region)))  &
               /         (dy/yref(parent(region))))
          lend(region) = nint((zend(region)-zbeg(parent(region)))  &
               /         (dz/zref(parent(region))))

          ! test
          xref_ = xref(region)/xref(parent(region))
          yref_ = yref(region)/yref(parent(region))
          zref_ = zref(region)/zref(parent(region))

          ! compute icells - how many cells of parent in x-direction a region occupies
          if (iend(region)>ibeg(region)) then
             icells = iend(region)-ibeg(region)+1
          else
             icells = im(parent(region)) - ibeg(region) + 1 + iend(region)
          end if
          if (icells*xref_/=im(region)) then
             print *, 'determine_children_etc: ', &
                  'check xbeg/xend, im and ibeg/iend for region ',region
             print *, 'determine_children_etc: xbeg,xend,im,ibeg,iend:', &
                  xbeg(region),xend(region),im(region),ibeg(region),iend(region)
             stop
          end if
          jcells = jend(region)-jbeg(region)+1
          if (jcells*yref_/=jm(region)) then
             print *, 'determine_children_etc: ', &
                  'check ybeg/yend, jm and jbeg/jend for region ',region
             print *, 'determine_children_etc: ybeg,yend,jm,jbeg,jend:', &
                  ybeg(region),yend(region),jm(region),jbeg(region),jend(region)
             stop
          end if
          lcells = lend(region)-lbeg(region)+1
          if (lcells*zref_/=lm(region)) then
             print *, 'determine_children_etc: ', &
                  'check zbeg/zend, lm and lbeg/lend for region ',region
             print *, 'determine_children_etc: zbeg,zend,lm,lbeg,lend:', &
                  zbeg(region),zend(region),lm(region),lbeg(region),lend(region)

          end if

       end do
       if ( okdebug ) then
          print *,'determine_children_etc: ibeg(2:nregions): ',ibeg
          print *,'determine_children_etc: iend(2:nregions): ',iend
          print *,'determine_children_etc: jbeg(2:nregions): ',jbeg
          print *,'determine_children_etc: jend(2:nregions): ',jend
          print *,'determine_children_etc: lbeg(2:nregions): ',lbeg
          print *,'determine_children_etc: lend(2:nregions): ',lend
          print *,'determine_children_etc: parents: ',parent
          do region=1,nregions
             print *,'determine_children_etc: children(',region, &
                  ',0:nregions)=',children(region,:)
          end do
       end if
       ! determine isr,ier,jsr,jer: the scope of the cells to be
       ! processed per region.
       ! processes like vertical mixing, chemistry, sources/sinks
       ! are processed on this scope:
       ! do i=isr(region),ier(region)
       ! do j=jsr(region),jer(region)

       isr(1) = 1
       ier(1) = im(1)
       jsr(1) = 1
       jer(1) = jm(1)
       do region = 2, nregions
          xref_ = xref(region)/xref(parent(region))
          yref_ = yref(region)/yref(parent(region))
          isr(region) = xref_ + 1
          ier(region) = im(region)-xref_
          jsr(region) = yref_+1
          jer(region) = jm(region)-yref_
          if(xcyc(region)/=0) then
             isr(region) = 1
             ier(region) = im(region)
          end if
          if(touch_sp(region)==1) jsr(region) = 1
          if(touch_np(region)==1) jer(region) = jm(region)
       end do

     end subroutine determine_children_etc

     !===========================================================================================================
     !===========================================================================================================

     subroutine define_splitorderzoom
       !
       ! splitorderzoom is splitorder specified for each region
       ! the subroutine is normally called once per run of TM5
       ! ATTN: subroutine was written based on assumption that a region can not
       ! have number less than number of its parent: region > parent(region)
       ! written by patrick berkvens and mike botchev, march-june 1999
       !
       use dims
       use var4d,                only : steps_region

       implicit none
       ! local
       integer :: i, j, j0, region, tref_
       !character,dimension(3):: reverse

       if ( okdebug ) print *,'define_splitorderzoom: '
       splitorderzoom = ' '
       splitorderzoom(1,1:nsplitsteps) = splitorder
       !PBi
       steps_region(1) = nsplitsteps !ADJ
       !PBf
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
             !PBi
             steps_region(region) = j-1  ! store the total number of steps for this region   !ADJ
             !PBf
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

     !===========================================================================================================
     !===========================================================================================================

     subroutine print_zoom(region_dat)
       !
       !
       use dims
       implicit none

       ! input
       type(region_data),dimension(nregions),intent(in)   :: region_dat

       ! const
       character(len=*), parameter      :: rname = mname//', print_zoom'

       ! local
       integer            ::  region, j
       character(len=30)  ::  fmt
       integer            ::  iostat

       ! start

       print *,'print_zoom'
       print *,'print_zoom: nregions: ',nregions
       print *,'print_zoom: xbeg: ',xbeg
       print *,'print_zoom: xend: ',xend
       print *,'print_zoom: ybeg: ',ybeg
       print *,'print_zoom: yend: ',yend
       print *,'print_zoom: zbeg: ',zbeg
       print *,'print_zoom: zend: ',zend
       print *,'print_zoom: xref: ',xref
       print *,'print_zoom: yref: ',yref
       print *,'print_zoom: zref: ',zref
       print *,'print_zoom: im: ',im
       print *,'print_zoom: jm: ',jm
       print *,'print_zoom: lm: ',lm
       print *,'print_zoom: isr: ',isr
       print *,'print_zoom: ier: ',ier
       print *,'print_zoom: jsr: ',jsr
       print *,'print_zoom: jer: ',jer

       print *,'print_zoom: ibeg(2:nregions): ',ibeg
       print *,'print_zoom: iend(2:nregions): ',iend
       print *,'print_zoom: jbeg(2:nregions): ',jbeg
       print *,'print_zoom: jend(2:nregions): ',jend
       print *,'print_zoom: lbeg(2:nregions): ',lbeg
       print *,'print_zoom: lend(2:nregions): ',lend

       print *,'print_zoom: parents: ',parent
       do region=1,nregions
          print *,'print_zoom: children(',region,',0:nregions)=', &
               children(region,:)
       end do
       do region=1,nregions
          print *,'print_zoom: splitorderzoom(',region,',:)=', &
               splitorderzoom(region,:),' '
       end do

       do region = 1,nregions
          print *, 'print_zoom: Zoomed ... region', region
          write( fmt, '("(a13,i3,a3,",i3.3,"i1)")' ) im(region)
          do j=jm(region),1,-1
             print trim(fmt),' print_zoom: ', j,'---',&
                  region_dat(region)%zoomed(:,j)
          end do
          print *, 'print_zoom: Edge   ... region', region
          do j=jm(region),1,-1
             write (*,fmt=trim(fmt),iostat=iostat) ' print_zoom: ', j,'---', max(0,region_dat(region)%edge(:,j))
             if (iostat/=0) then
               write (gol,'("writing edge:")'); call goErr
               write (gol,'("writing edge to format : ",a)') trim(fmt); call goErr
               TRACEBACK; stop
              end if
          end do
       end do

     end subroutine print_zoom

     !===========================================================================================================
     !===========================================================================================================

     subroutine blank_zoom_region(region,region_dat)
       !
       ! determine which cells in a region are calculated in an overlying
       ! zoom region. Only the core of a zoom region is flagged. The
       ! interface cells are not flagged, since they are treated coarse
       ! ref: Krol et al, IAMAS paper, 2001
       !
       ! Example    region1:
       ! Zoomed  I2 represent the core-zoom region
       ! 1111111111111111111
       ! 1111111111111111111
       ! 1111222111111111111
       ! 1111222111111111111
       ! 1111222111111111111
       ! 1111222111111111111
       ! 1111111111111111111
       ! 1111111111111111111
       ! Edge     !indicates the interface cells
       ! 0000000000000000000
       ! 0003222300000000000
       ! 0001000100000000000
       ! 0001000100000000000
       ! 0001000100000000000
       ! 0001000100000000000
       ! 0003222300000000000
       ! 0000000000000000000
       !     3 means x+y, 2 means y, and 1 means only x interface cell.
       !
       use dims
       implicit none

       ! input/output
       integer, intent(in)               :: region
       type(region_data),intent(inout)   :: region_dat

       ! local
       integer                           :: i, j, n, n1, xref_, yref_
       integer,dimension(:,:),pointer    :: zoomed, edge

       zoomed => region_dat%zoomed
       edge   => region_dat%edge

       zoomed = region
       edge = 0
       do n=1,children(region,0)  !loop over the children
          n1 = children(region,n)    !child region

          do i=1,im(region)
             if( (i <= ibeg(n1) .or. i >= iend(n1)) .and. (xcyc(n1) /= 1) ) cycle
             do j=1,jm(region)
                if( (j <= jbeg(n1)) .and. (touch_sp(n1)/=1) ) cycle
                if( (j >= jend(n1)) .and. (touch_np(n1)/=1) ) cycle
                zoomed(i,j) = n1
             end do
          end do

          if ( xcyc(n1) /= 1 ) then
             edge(ibeg(n1),jbeg(n1):jend(n1)) = 1
             edge(iend(n1),jbeg(n1):jend(n1)) = 1
          end if
          if ( touch_sp(n1) /= 1 ) then
             edge(ibeg(n1):iend(n1),jbeg(n1)) = &
                  edge(ibeg(n1):iend(n1),jbeg(n1)) + 2
          end if
          if ( touch_np(n1) /= 1 ) then
             edge(ibeg(n1):iend(n1),jend(n1)) = &
                  edge(ibeg(n1):iend(n1),jend(n1)) + 2
          end if

       end do

       ! put the edge of the region itself to -1.....
       if ( parent(region) /= 0 ) then
          xref_ = xref(region)/xref(parent(region))
          yref_ = yref(region)/yref(parent(region))
          if(xcyc(region)/=1) then
             edge(1:xref_,:) = -1
             edge(im(region)-xref_+1:im(region),:) = -1
          end if
          if(touch_sp(region)/=1) edge(:,1:yref_) = -1
          if(touch_np(region)/=1) edge(:,jm(region)-yref_+1:jm(region)) = -1
       end if

       nullify(zoomed)
       nullify(edge)

     end subroutine blank_zoom_region

     !===========================================================================================================
     !===========================================================================================================

   end module initexit
