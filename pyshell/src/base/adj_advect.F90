!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_Advect

  use GO, only : gol, goPr, goErr

  implicit none

  ! ------------ interface ---------------------

  private

  public :: adj_Advect_Init, adj_Advect_Done

  
  ! --- const ------------------------------

  character(len=*), parameter ::  mname = 'adj_Advect'


contains

  
  ! ====================================================================


  subroutine adj_Advect_Init( status )

    use adj_AdvectX, only : adj_AdvectX_Init
    use adj_AdvectY, only : adj_AdvectY_Init
    use adj_AdvectZ, only : adj_AdvectZ_Init

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/adj_Advect_Init'

    ! --- begin ----------------------------------

    call adj_AdvectX_Init( status )
    IF_NOTOK_RETURN(status=1)

    call adj_AdvectY_Init( status )
    IF_NOTOK_RETURN(status=1)

    call adj_AdvectZ_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_Advect_Init
  

  ! ***


  subroutine adj_Advect_Done( status )

    use adj_AdvectX, only : adj_AdvectX_Done
    use adj_AdvectY, only : adj_AdvectY_Done
    use adj_AdvectZ, only : adj_AdvectZ_Done

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/adj_Advect_Done'

    ! --- begin ----------------------------------

    call adj_AdvectX_Done( status )
    IF_NOTOK_RETURN(status=1)

    call adj_AdvectY_Done( status )
    IF_NOTOK_RETURN(status=1)

    call adj_AdvectZ_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_Advect_Done
  

end module adj_Advect
