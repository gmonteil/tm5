!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module user_output
  !--------------------------------------------------------------------------
  ! contains calls to user-specific output routines, e.g.
  ! instantaneous mix files, station output, output of flight tracks etc.
  !--------------------------------------------------------------------------

  use go,          only: gol, goErr


  implicit none

  private

  public :: user_output_init
  public :: user_output_step
  public :: user_output_mean
  public :: user_output_done

!  public :: station_data
!  public :: satellite_data
!  public :: flask_data

  character(len=*), parameter         :: mname = 'user_output'

  ! ------------ Settings ------------

  logical           :: model_data = .false.   ! signal for model output
!  logical           :: satellite_data      ! flag for station output
  logical           :: run_apri = .true.   ! set for first run...
!  logical           :: flask_data = .false.! signal for flask data, for Wouter's subroutine

  ! ---------- End settings ----------


contains

!===========================================================================================================
!===========================================================================================================

  subroutine user_output_init(status)
    !-------------------------------------------------------------
    ! user_output_init:
    !      Initialise user-specified model output (all regions)
    !-------------------------------------------------------------

    use GO,                    only : ReadRc
    use global_data,           only : rcF
    use dims,                  only : revert
    use User_Output_Model    , only : User_Output_Model_Init
!    use user_output_satellite, only : init_satellitedata
!    use user_output_flask,     only : user_output_flask_init



    implicit none

    !__IO____________________________________________________________________

    integer, intent(out) ::  status


    !__LOCAL_VARIABLES_______________________________________________________

    character(len=*), parameter      :: rname = mname//', user_output_init'

    character(len=20)   :: sflag

    !__START_SUBROUTINE______________________________________________________


    status = -1

    call ReadRc( rcF, 'model.output', model_data, status)
    IF_NOTOK_RETURN(status=1)
!    call ReadRc( rcF, 'output.satellite', satellite_data, status)
!    IF_NOTOK_RETURN(status=1)
!    call ReadRc( rcF, 'output.point', flask_data, status)
!    IF_NOTOK_RETURN(status=1)


    ! model data:
    if ( model_data ) then
      call User_Output_Model_Init( status )
      IF_NOTOK_RETURN(status=1)
    end if

!    if (satellite_data) call init_satellitedata

    ! flask data
!    if (flask_data) call user_output_flask_init(rcF)

    ! CMK: first run ....write, second run check:
!    if (run_apri) then
!       call obs_operator_control_write_init(status)
!       IF_NOTOK_RETURN(status=1)
!    else
!       call obs_operator_control_read(status)
!       IF_NOTOK_RETURN(status=1)
!    endif



    status = 0

  end subroutine user_output_init

!===========================================================================================================
!===========================================================================================================

  subroutine user_output_step ( region , status )
    !-------------------------------------------------------------
    ! user_output_step:
    !      Define user-specified model output for the region given
    !      Called every time step
    !-------------------------------------------------------------

    use dims,                 only : itaur
    use datetime,             only : tau2date
!    use user_output_satellite,only : output_satellitedata
!    use user_output_flask,    only : user_output_flask_sample

    implicit none

    !__IO___________________________________________________________________

    integer,intent(inout)           :: status
    integer,intent(in)              :: region

    !__LOCAL_VARIABLES_______________________________________________________

    character(len=*), parameter     :: rname = mname//', user_output_step'

    integer,dimension(6)            :: idate_f
    integer                         :: i


    !__START_SUBROUTINE______________________________________________________


    call tau2date(itaur(region),idate_f)

!    if(satellite_data) then
!       call output_satellitedata(region, status)
!       IF_NOTOK_RETURN(status=1)
!    endif
!    if(flask_data)     call user_output_flask_sample(region,itaur(region))

  end subroutine user_output_step


!===========================================================================================================
!===========================================================================================================

  subroutine user_output_mean( status )

    !__IO___________________________________________________________________

    integer, intent(out)          :: status

    !__CONST________________________________________________________________

    character(len=*), parameter  ::  rname = mname//'/user_output_mean'

    !__LOCAL_VARIABLES______________________________________________________

    !__START_SUBROUTINE______________________________________________________

    ! ok
    status = 0
    
  end subroutine user_output_mean

!===========================================================================================================
!===========================================================================================================

  subroutine user_output_done ( msg, status )
    !-------------------------------------------------------------
    ! user_output_done:
    !      Finalise user-specified model output for the region given
    !-------------------------------------------------------------

    use dims,                   only : nregions
    use dims,                   only : revert
    use var4d,                  only : iter

    use User_Output_Model     , only : User_Output_Model_Done
!    use user_output_satellite,  only : free_satellitedata
!    use user_output_flask,      only : user_output_flask_done


    implicit none

    !__IO___________________________________________________________________

    character(len=*),intent(in)  :: msg
    integer, INTENT(OUT)         :: status


    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//', user_output_done'


    !__START_SUBROUTINE______________________________________________________

    status = -1

!          ! a priori (forward) run
!    if (run_apri) then
!        call obs_operator_control_write_done(status)
!        IF_NOTOK_RETURN(status=1)
!    else

!        call obs_operator_control_read_done(status)
!        IF_NOTOK_RETURN(status=1)
!
!    endif




    ! model data:
    if ( model_data ) then
      call User_Output_Model_Done( status )
      IF_NOTOK_RETURN(status=1)
    end if

!    if (satellite_data) call free_satellitedata

!    if (flask_data) call user_output_flask_done

    status = 0

  end subroutine user_output_done

!===========================================================================================================
!===========================================================================================================

end module user_output
