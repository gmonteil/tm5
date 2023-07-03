!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_user_output


!  use user_output, only : satellite_data, flask_data


  use go,          only: gol, goErr

  implicit none

  private

  public :: adj_user_output_init
  public :: adj_user_output_step
  public :: adj_user_output_mean
  public :: adj_user_output_done

  character(len=*), parameter         :: mname = 'adj_user_output'


contains

!===========================================================================================================
!===========================================================================================================

  subroutine adj_user_output_init(status)

    use GO,                        only : ReadRc
    use global_data,               only : rcF

!    use adj_user_output_satellite, only : adj_init_satellitedata
!    use adj_user_output_flask,     only : adj_user_output_flask_init

!    use observation_operator,      only : obs_operator_control_read

    implicit none

    !__IO___________________________________________________________________

    integer, intent(OUT)            :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//', adj_user_output_init'

    !__START_SUBROUTINE______________________________________________________

!    call ReadRc( rcf, 'output.satellite', satellite_data , status)
!    IF_NOTOK_RETURN(status=1)
!    call ReadRc( rcF, 'output.point', flask_data, status)
!    IF_NOTOK_RETURN(status=1)

!    if (satellite_data) then
!       call adj_init_satellitedata(status)
!       IF_NOTOK_RETURN(status=1)
!    endif

!    if (flask_data) call adj_user_output_flask_init(rcF)

    status = 0

  end subroutine adj_user_output_init

!===========================================================================================================
!===========================================================================================================

  subroutine adj_user_output_step(region, status)

!    use adj_user_output_satellite, only : adj_output_satellitedata
!    use adj_user_output_flask,     only : adj_user_output_flask_step

    implicit none

    !__IO___________________________________________________________________

    integer,intent(in)              :: region
    integer,intent(inout)           :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//', adj_user_output_step'

    !__START_SUBROUTINE_____________________________________________________


!    if (flask_data) call adj_user_output_flask_step(region)
!    if(satellite_data) call adj_output_satellitedata(region, status)
!    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_user_output_step


!===========================================================================================================
!===========================================================================================================

  subroutine adj_user_output_mean

    use dims,                     only: itau
    use dims,                     only: revert

    use var4d,                    only: nasim

    use toolbox,                  only: escape_tm


!    use adj_observation_operator, only: adj_obs_operator


    implicit none

    !__IO___________________________________________________________________

    !__LOCAL_VARIABLES______________________________________________________


    integer ::  status

    !__START_SUBROUTINE_____________________________________________________


    !=====================
    ! observation operator
    !=====================
!    if(revert == -1) then
!       if(mod(itau, nasim) == 0) then
!          call adj_obs_operator(station_data, column_data, status)
!
!          if(status /= 0) then
!             call escape_tm('adj_user_output_mean: error adj_obs_operator')
!          endif
!
!       endif
!    endif

  end subroutine adj_user_output_mean

!===========================================================================================================
!===========================================================================================================

  subroutine adj_user_output_done(status)

!    use adj_user_output_satellite, only : adj_free_satellitedata
!    use adj_user_output_flask,     only : adj_user_output_flask_done

    implicit none

    !__IO___________________________________________________________________

    integer, intent(OUT)            :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//', adj_user_output_done'


    !__START_SUBROUTINE______________________________________________________

!    if(satellite_data) then
!      call adj_free_satellitedata( status )
!      IF_NOTOK_RETURN(status=1)
!    endif

!    if (flask_data)    call adj_user_output_flask_done

!    call obs_operator_control_read_done(status)
!    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine adj_user_output_done

!===========================================================================================================
!===========================================================================================================

end module adj_user_output


