!###############################################################################
!
!       purpose
!       -------
!       perform emissions needed for TM5 4DVAR
!
!       interface
!       ---------
!       call emission_Init
!       call emission_Apply
!       call emission_Done
!
!       method
!       ------
!       subroutine emission_Init    is called from trace0
!       subroutine emission_Apply   is called from source1
!       subroutine emission_Done    is called from trace_end
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

module Emission

  use GO, only : gol, goErr, goPr

!  use Emission_Fwd , only : emissions
  use Emission_Fwd , only : Emission_Fwd_Setup
  use Emission_Fwd , only : Emission_Fwd_Apply

!  use Emission_Adj , only : adj_emissions
  use Emission_Adj , only : Emission_Adj_Setup
  use Emission_Adj , only : Emission_Adj_Apply

  use Emission_Data, only : tracers_em_info
  use Emission_Data, only : optim_emis_type
  use Emission_Data, only : ref_emissions_apri, ref_emissions

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  Emission_Init, Emission_Done

!  public  ::  emissions
  public  ::  Emission_Fwd_Setup
  public  ::  Emission_Fwd_Apply

!  public  ::  adj_emissions
  public  ::  Emission_Adj_Setup
  public  ::  Emission_Adj_Apply

  public  ::  optim_emis_type
  public  ::  tracers_em_info
  public  ::  ref_emissions, ref_emissions_apri


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'Emission'



contains


  ! ====================================================================


  subroutine Emission_Init( status )

    ! --- modules --------------------------------

    use Emission_Data, only : Emission_Data_Init
    use Emission_Fwd , only : Emission_Fwd_Init
    use Emission_Adj , only : Emission_Adj_Init

    ! --- in/out ---------------------------------

    integer, intent(out)             :: status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Emission_Init'

    !--- local -----------------------------------

    !--- begin -----------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! init shared emission data:
    call Emission_Data_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init emission data for forward run:
    call Emission_Fwd_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init emission data for adjoint run:
    call Emission_Adj_Init( status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Emission_Init


  ! ***


  subroutine emission_Done( status )

    ! --- modules --------------------------------

!    use Emission_Data, only : Emission_Data_Done
    use Emission_Fwd , only : Emission_Fwd_Done
    use Emission_Adj , only : Emission_Adj_Done

    ! --- in/out ---------------------------------

    integer, intent(out)             :: status

    ! --- const ----------------------------------

    character(len=*), parameter      :: rname = mname//'/emission_Done'

    !--- local -----------------------------------

    !--- begin -----------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! done with emission data for forward run:

    call Emission_Fwd_Done( status )
    IF_NOTOK_RETURN(status=1)


    ! done with emission data for adjoint run:

    call Emission_Adj_Done( status )
    IF_NOTOK_RETURN(status=1)


    ! done shared emission data: ==> deallocation is automatic at the end of the run ...
    ! call Emission_Data_Done( status )
    ! IF_NOTOK_RETURN(status=1)

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Emission_Done


end module emission
