!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_chemistry

  use chemistry,                 only : chemistry_apply
  use go,                        only : gol, goErr

  implicit none

  ! --- in/out -----------------------------

  private

  public :: adj_Chemistry_Init
  public :: adj_Chemie
  public :: adj_Chemistry_Done

  ! --- interfaces--------------------------

  interface adj_Chemie
     module procedure adj_Chemistry_Apply
  end interface

  ! --- const ------------------------------

  character(len=*), parameter        :: mname = 'module adj_chemistry'

  ! timers:
  integer             ::  itim_chemistry


contains

  subroutine adj_Chemistry_Init( status )

    !---------------------------------------------------------
    ! Initialize CH4 chemistry:
    !  o calculate reaction rate
    !  o read OH field
    !---------------------------------------------------------

    ! --- modules ------------------------------
    use GO, only : ReadRc
    use global_data, only : rcF

    ! --- in/out ----------------------------------------------

    integer, intent(out)             :: status

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//', Init'

    !--- local ------------------------------------------------

    call ReadRc(rcF, 'proces.chemistry', chemistry_apply, status)
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine Adj_Chemistry_Init

  subroutine adj_Chemistry_Apply( region, status )

    ! --- in/out ----------------------------------------------

    integer, intent(in)              :: region
    integer, intent(out)             :: status

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//', Apply'

    ! ok
    status = 0

  end subroutine Adj_Chemistry_Apply



  subroutine adj_Chemistry_Done

    ! --- modules ------------------------------

    ! --- in/out ----------------------------------------------

    ! --- const ------------------------------

    character(len=*), parameter      :: rname = mname//', Done'

  end subroutine Adj_Chemistry_Done

end module adj_chemistry
