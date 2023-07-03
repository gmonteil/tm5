!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module adj_user_output

  use go,          only: gol, goErr

  implicit none

  private

  public :: adj_user_output_init
  public :: adj_user_output_step
  public :: adj_user_output_done

  character(len=*), parameter         :: mname = 'adj_user_output'

  logical           :: flask_data = .false.     ! signal for flask data, for Wouter's subroutine
  logical           :: satellite_data = .false. ! signal for satellite data

contains

!===========================================================================================================
!===========================================================================================================

  subroutine adj_user_output_init(status)

    use GO,                        only : ReadRc
    use global_data,               only : rcF

    use adj_user_output_satellite, only : adj_user_output_satellite_init
    use adj_user_output_flask,     only : adj_user_output_flask_init

    implicit none

    !__IO___________________________________________________________________

    integer, intent(OUT)            :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//', adj_user_output_init'

    !__START_SUBROUTINE______________________________________________________

    call ReadRc( rcf, 'adjoint.input.satellite', satellite_data , status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'adjoint.input.point', flask_data, status)
    IF_NOTOK_RETURN(status=1)

    if (flask_data) call adj_user_output_flask_init(status)
    IF_NOTOK_RETURN(status=1)

    if (satellite_data) call adj_user_output_satellite_init(status)
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine adj_user_output_init

!===========================================================================================================
!===========================================================================================================

  subroutine adj_user_output_step(region, tr, status)

    !use dims,                 only : itaur
    !use datetime,             only : tau2date
    use Go,                         only : TDate
    use adj_user_output_flask,      only : adj_user_output_flask_step
    use adj_user_output_satellite,  only : adj_user_output_satellite_step

    implicit none

    !__IO___________________________________________________________________

    integer,intent(in)              :: region
    type(TDate), intent(in)         :: tr(2)
    integer,intent(inout)           :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//'/adj_user_output_step'
    !integer,dimension(6)             :: idate_f

    !__START_SUBROUTINE_____________________________________________________

    if (flask_data) call adj_user_output_flask_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)

    if (satellite_data) call adj_user_output_satellite_step(region, tr, status)
    !if (satellite_data) call adj_user_output_satellite_step(region, status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_user_output_step

!===========================================================================================================
!===========================================================================================================

  subroutine adj_user_output_done(status)

    use adj_user_output_flask,     only : adj_user_output_flask_done
    use adj_user_output_satellite, only : adj_user_output_satellite_done

    implicit none

    !__IO___________________________________________________________________

    integer, intent(OUT)            :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//', adj_user_output_done'


    !__START_SUBROUTINE______________________________________________________


    if (flask_data)    call adj_user_output_flask_done(status)
    IF_NOTOK_RETURN(status=1)
    if (satellite_data) call adj_user_output_satellite_done(status)
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine adj_user_output_done

!===========================================================================================================
!===========================================================================================================

end module adj_user_output
