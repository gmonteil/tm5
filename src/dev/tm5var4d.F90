!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module TM5var4D

  use GO, only : gol, goPr, goErr
  use os_specs, only : MAX_FILENAME_LEN

  implicit none


  ! --- in/out -----------------------------------

  private

  public    ::  T_TM5var4D
  public    ::  TM5var4D_Init, TM5var4D_Done
  public    ::  TM5var4D_Run


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'TM5var4D'


  ! --- types ------------------------------------

  type T_TM5var4D
    ! timers:
    integer               ::  itim_runmode_forward
    integer               ::  itim_runmode_backward
    ! output directory
    character(len=MAX_FILENAME_LEN)   ::  outdir
  end type T_TM5var4D


contains


  ! ====================================================================
  ! ===
  ! === main
  ! ===
  ! ====================================================================


  subroutine TM5var4D_Init( TM5var4D, status )

    use GO        , only : GO_Init
    use Tracer    , only : Tracer_Init
    use Adj_Tracer, only : Adj_Tracer_Init
    use Var4D     , only : Init_Var4D
    use datetime,   only : system_clock_value_start, wall_clock_time

    ! --- in/out ---------------------------------

    type(T_TM5var4D), intent(out)      ::  TM5var4D
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    system_clock_value_start = wall_clock_time()

    ! init GO modules:
    write (*,'("  init GO modules ...")')
    call GO_Init( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from GO_Init : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! init settings:
    write (*,'("  init TM5/4D-var settings ...")')
    call TM5var4D_Settings_Init( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from TM5var4D_Settings_Init : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! init messages;
    write (*,'("  init messages ...")')
    call TM5var4D_Messages_Init( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from TM5var4D_Messages_Init : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! * for messages, use gol/goPr from here *

    ! init timing:
    write (gol,'("  init timing ...")'); call goPr
    call TM5var4D_Timing_Init( TM5var4D, status )
    IF_NOTOK_RETURN(status=1)

    ! init tracer model (module tracer):
    !   timers
    !   call Init_TM5 (from module initexit)
    !     paths, regions, operator splitting ...
    write (gol,'("  init tracer model ...")'); call goPr
    call Tracer_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init part of adjoint tracer model (module adj_tracer)
    !   timers
    write (gol,'("  init adjoint tracer model ...")'); call goPr
    call Adj_Tracer_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init 4D-var stuff (module var4d):
    call Init_Var4D( status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine TM5var4D_Init


  ! ***


  subroutine TM5var4D_Done( TM5var4D, status )

    use GO        , only : GO_Done
    use InitExit  , only : TM5_OkFile
    use Tracer    , only : Tracer_Done
    use Adj_Tracer, only : Adj_Tracer_Done
    use Var4D     , only : Free_Var4D
    use datetime,   only : system_clock_value_start, wall_clock_time

    ! --- in/out ---------------------------------

    type(T_TM5var4D), intent(inout)     ::  TM5var4D
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Done'

    ! --- local ----------------------------------

    double precision    :: system_clock_value_end

    ! --- begin ----------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! done with 4D-var:
    call Free_Var4D( status )
    IF_NOTOK_RETURN(status=1)

    ! done with tracer model:
    call Tracer_Done( status )
    IF_NOTOK_RETURN(status=1)

    call Adj_Tracer_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with timing:
    call TM5var4D_Timing_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! write ok file to indicate proper model end;
    ! (this routines uses GO modules, therefore not the final action):
    call TM5_OkFile( status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr

    system_clock_value_end = wall_clock_time()
    write(gol, '("Wall clock time taken for model run = ", f11.2, " seconds")') system_clock_value_end - system_clock_value_start ; call goPr

    ! * do not use gol/goPr from here *

    ! done with messages:
    call TM5var4D_Messages_Done( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from TM5var4D_Messages_Done : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! done with settings:
    call TM5var4D_Settings_Done( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from TM5var4D_Settings_Done : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! done with GO modules:
    call GO_Done( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from GO_Done : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! ok
    status = 0

  end subroutine TM5var4D_Done


  ! ***


  subroutine TM5var4D_Run( tm5var4d, status )

    use dims,                  only : revert, kdebug, icalendo, julday0, iyear0, ndyn
    use dims,                  only : nregions
    use dims,                  only : itau
    use GO,                    only : gol, goErr, goPr
    use GO,                    only : goLabel
    use GO,                    only : Init, Done, ReadRc
    use GO,                    only : pathsep
    use GO                   , only : goGetFU
    use global_data,           only : rcfile, rcF
    use Tracer               , only : Tracer_Model
    use adj_Tracer           , only : Adj_Tracer_Model
    use var4d,                 only : RUN_FORWARD, RUN_BACKWARD
    use var4d,                 only : init_var4d, free_var4d

    use initexit,              only : init_TM5, free_TM5, start_TM5, exit_TM5
    use datetime,              only : tau2date, date2tau, julday
    use toolbox,               only : escape_tm

    use Meteo,                 only : Meteo_Init, Meteo_Done

    use var4d,                 only : run_mode
    use dims,                  only : im, jm

    ! --- in/out ---------------------------------

    type(T_TM5var4D), intent(inout)     ::  tm5var4d
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter    :: rname = mname//'/TM5var4D_Run'

    ! --- local ----------------------------------

    integer             :: region
    character(len=MAX_FILENAME_LEN)  :: subdir    ! output subdirectory

    real                :: cost
    integer             :: n_iter    ! total number of iterations performed
    character(len = 20) :: fmt       ! format string

    integer             :: istat

    integer             :: iter_inner_max

    character(len=MAX_FILENAME_LEN)  :: FFileName
    integer             :: kstatus

    integer             :: i_cat
    integer             :: i_period
    integer             :: i,j, l

    ! --- begin ----------------------------------

    ! first label:
    call goLabel( label='main' )

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! output directory:
    call ReadRc(rcF, 'output.dir', tm5var4d%outdir, status)
    IF_NOTOK_RETURN(status=1)

    !================
    ! select run mode
    !================

    write (gol,'("")'); call goPr
    write (gol,'("run mode : ",i3)') run_mode; call goPr
    write (gol,'("")'); call goPr

    select case ( run_mode )


      ! ====================================================================
      case ( RUN_FORWARD )  ! FORWARD MODE
      ! ====================================================================

        ! start label:
        call goLabel( timer=tm5var4d%itim_runmode_forward )

        call Tracer_Model( status )
        IF_NOTOK_RETURN(status=1)

        ! end label:
        call goLabel()

      ! ====================================================================
      case ( RUN_BACKWARD )     ! BACKWARD MODE
      ! ====================================================================

        ! start label:
        call goLabel( timer=tm5var4d%itim_runmode_backward )

        call Adj_Tracer_Model( status )
        IF_NOTOK_RETURN(status=1)

        ! end label:
        call goLabel()


      ! ====================================================================
      case default
      ! ====================================================================

         write (gol,'("Invalid run mode : ",i6)') run_mode; call goErr
         TRACEBACK; status=1; return

    end select

    !==========================
    ! end
    !==========================

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! final label:
    call goLabel()

  end subroutine TM5var4D_Run


  ! ====================================================================
  ! ===
  ! === settings
  ! ===
  ! ====================================================================


  subroutine TM5var4D_Settings_Init( status )

    use GO         , only : Init
    use Global_Data, only : rcfile, rcF
    use InitExit   , only : TM5_Arguments

    ! --- in/out ---------------------------------

    integer, intent(inout)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Settings_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! read name of rcfile from arguments:
    call TM5_Arguments( status )
    IF_NOTOK_RETURN(status=1)

    ! read settings:
    call Init( rcF, trim(rcfile), status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TM5var4D_Settings_Init


  ! ***


  subroutine TM5var4D_Settings_Done( status )

    use GO         , only : Done
    use Global_Data, only : rcF

    ! --- in/out ---------------------------------

    integer, intent(inout)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Settings_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! done with settings:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TM5var4D_Settings_Done


  ! ====================================================================
  ! ===
  ! === mesages
  ! ===
  ! ====================================================================


  subroutine TM5var4D_Messages_Init( status )

    use GO, only : GO_Print_Init
    use GO, only : GO_Label_Init

    ! --- in/out ---------------------------------

    integer, intent(inout)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Messages_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! enable messages:
    call GO_Print_Init( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from GO_Print_Init : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! enable labels:
    call GO_Label_Init( status, trace=.true. )
    IF_NOTOK_RETURN(status=0)

    ! intro messsage:
    write (gol,'("")'); call goPr
    write (gol,'("================================")'); call goPr
    write (gol,'(" start main program             ")'); call goPr
    write (gol,'("================================")'); call goPr
    write (gol,'("")'); call goPr

    ! ok
    status = 0

  end subroutine TM5var4D_Messages_Init


  ! ***


  subroutine TM5var4D_Messages_Done( status )

    use GO, only : GO_Print_Done
    use GO, only : GO_Label_Done

    ! --- in/out ---------------------------------

    integer, intent(inout)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Messages_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! (almost?) final message:
    write (gol,'("===============================")'); call goPr
    write (gol,'("TM5 4DVAR terminated normally  ")'); call goPr
    write (gol,'("===============================")'); call goPr

    ! done with labels:
    call GO_Label_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with messages:
    call GO_Print_Done( status )
    if (status/=0) then
      write (*,'("ERROR - non-zero return status from  GO_Print_Done : ", i6)') status
      write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
      status=1; return
    end if

    ! ok
    status = 0

  end subroutine TM5var4D_Messages_Done


  ! ====================================================================
  ! ===
  ! === Timing
  ! ===
  ! ====================================================================


  subroutine TM5var4D_Timing_Init( tm5var4d, status )

    use GO, only : GO_Timer_Init, GO_Timer_Def

    ! --- in/out ---------------------------------

    type(T_TM5var4D), intent(inout)     ::  tm5var4d
    integer, intent(inout)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Timing_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! enable timing:
    call GO_Timer_Init( status )
    IF_NOTOK_RETURN(status=0)

    ! define timers:
    call GO_Timer_Def( tm5var4d%itim_runmode_forward , 'runmode forward'   , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( tm5var4d%itim_runmode_backward, 'runmode backward'  , status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TM5var4D_Timing_Init


  ! ***


  subroutine TM5var4D_Timing_Done( status )

    use GO, only : pathsep
    use GO, only : ReadRc
    use GO, only : GO_Timer_Done
    use misctools, only : check_dir

    use Global_Data, only : rcfile, rcF

    ! --- in/out ---------------------------------

    integer, intent(inout)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/TM5var4D_Timing_Done'

    ! --- local ----------------------------------

    integer               ::  l, slash_pos
    character(len=MAX_FILENAME_LEN)   ::  outdir
    character(len=MAX_FILENAME_LEN)    ::  subdir
    character(len=MAX_FILENAME_LEN)   ::  timing_file

    ! --- begin ----------------------------------

    ! output directory:
    call ReadRc( rcF, 'output.dir', outdir, status )
    IF_NOTOK_RETURN(status=1)
    ! timing subdirectory:
    call ReadRc( rcF, 'timing.output.subdir', subdir, status)
    IF_NOTOK_RETURN(status=1)

    ! output file for time profile:
    ! Sourish Basu, 20/10/2011
    ! The rcfile has the full path, so we need to extract the basename
    slash_pos = index(rcfile, pathsep, .true.) ! slash_pos refers to the index of the last '/'
    l = len_trim(rcfile)
    write (timing_file,'(6a)') trim(outdir), pathsep, trim(subdir), pathsep, rcfile(slash_pos+1:l-3), '.prf'

    ! done with timers; write profile to standard output and file:
    call check_dir(trim(timing_file))
    call GO_Timer_Done( status, file=trim(timing_file) )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TM5var4D_Timing_Done


end module TM5var4D

