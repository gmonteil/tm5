!###############################################################################
!
! the "test" project does not perform chemistry
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

module chemistry

  use GO, only : gol, goErr, goPr

  implicit none

  private

  public :: chemie, chemistry_init, chemistry_done, chemistry_apply

  logical :: chemistry_apply

  character(len=*), parameter :: mname = 'chemistry'

contains


  subroutine chemie( region, status )

    implicit none

    ! input/output
    integer,intent(in)                :: region
    integer,intent(out)               :: status

    ! ok
    status = 0

  end subroutine chemie

  subroutine chemistry_init( status )

    use GO, only : ReadRc
    use global_data, only : rcF

    integer,intent(out)               :: status
    character(len=*), parameter       :: rname = mname//'/chemistry_init'

    call ReadRc(rcF, 'proces.chemistry', chemistry_apply, status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine chemistry_init

  subroutine chemistry_done

    implicit none

  end subroutine chemistry_done


end module chemistry
