!###############################################################################
!
! Advection process.
!
! Module hierchy:
!   advect
!     advect[xyz]
!       advectm_cfl
!         advect_tools
!
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

module advect

  use GO, only : gol, goPr, goErr

  use Advect_Tools, only : dynam0

  implicit none

  ! ------------ interface ---------------------

  private

  public :: Advect_Init, Advect_Done
  public :: dynam0

  
  ! --- const ------------------------------

  character(len=*), parameter ::  mname = 'Advect'


contains

  
  ! ====================================================================


  subroutine Advect_Init( status )

    use AdvectX, only : AdvectX_Init
    use AdvectY, only : AdvectY_Init
    use AdvectZ, only : AdvectZ_Init

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Advect_Init'

    ! --- begin ----------------------------------

    call AdvectX_Init( status )
    IF_NOTOK_RETURN(status=1)

    call AdvectY_Init( status )
    IF_NOTOK_RETURN(status=1)

    call AdvectZ_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Advect_Init
  

  ! ***


  subroutine Advect_Done( status )

    use AdvectX, only : AdvectX_Done
    use AdvectY, only : AdvectY_Done
    use AdvectZ, only : AdvectZ_Done

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Advect_Done'

    ! --- begin ----------------------------------

    call AdvectX_Done( status )
    IF_NOTOK_RETURN(status=1)

    call AdvectY_Done( status )
    IF_NOTOK_RETURN(status=1)

    call AdvectZ_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Advect_Done

end module advect
