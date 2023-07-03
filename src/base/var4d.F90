!###############################################################################
!
! Module containing data and data types that are 'global' in 4DVAR routines
! and general 4DVAR subroutines
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

module var4d

  use GO, only : gol, goErr, goPr, goLabel

  ! use all objects from the data module:
  use Var4D_Data

  implicit none


  ! --- in/out -----------------------------------

  ! all public:
  public


  ! --- const ------------------------------------

  character(len=*), parameter                 :: mname = 'var4d'



contains


  !==============================================================================================
  !==============================================================================================


  subroutine init_var4d( status )

    use go,                        only : ReadRc
    use global_data,               only : rcF
    use Var4D_IO_Mass            , only : Var4D_IO_Mass_Init
    use Var4D_State_Hori         , only : Var4D_State_Hori_Init

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)               :: status

    !__CONST________________________________________________________________

    character(len=*), parameter        :: rname = mname//'/init_var4d'

    !__LOCAL_VARIABLES______________________________________________________

    !__START_SUBROUTINE______________________________________________________

    call goLabel(rname)

    write (gol,'("==========================")'); call goPr
    write (gol,'("  init_var4d              ")'); call goPr
    write (gol,'("==========================")'); call goPr

    ! runmode:
    call ReadRc( rcF, 'my.runmode', run_mode, status)
    IF_NOTOK_RETURN(status=1)

    ! setup mass output:
    call Var4D_IO_Mass_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! init state vector:
    call Var4D_State_Hori_Init( status )
    IF_NOTOK_RETURN(status=1)

    call goLabel()

    ! ok
    status = 0

  end subroutine init_var4d


  !==============================================================================================
  !==============================================================================================

  subroutine free_var4d( status )
    use Var4D_IO_Mass              , only : Var4D_IO_Mass_Done
    use Var4D_State_Hori           , only : Var4D_State_Hori_Done

    !__IO___________________________________________________________________

    integer, intent(out)               :: status

    !__CONST________________________________________________________________

    character(len=*), parameter        :: rname = mname//'/free_var4d'

    !__LOCAL_VARIABLES______________________________________________________

    !__START_SUBROUTINE______________________________________________________

    ! done with mass io:
    call Var4D_IO_Mass_Done( status )
    IF_NOTOK_RETURN(status=1)

    !
    call Var4D_State_Hori_Done( status )
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine free_var4d


  !=============================================================================================
  !=============================================================================================


  subroutine set_runid

    use dims,                      only : revert

    !__IO___________________________________________________________________

    !__LOCAL_VARIABLES______________________________________________________

    !__START_SUBROUTINE______________________________________________________


       if ( revert == 1 ) then   ! forward

          select case( iter )
          case ( 0 )
             runid = 'APRIF'   ! forward
          case ( -1 )
             runid = 'APOSF'   ! forward
          case ( -2 )
             runid = 'DAPOS'   ! uncertainty step after iterations
          case default
             write(runid,'(i4.4,a1)') iter, 'F'
          end select

       else if ( revert == -1 ) then   ! adjoint

          select case ( iter )
          case ( 0 )
             runid = 'APRIA'   ! adjoint
          case ( -1 )
             runid = 'APOSA'   ! adjoint
          case default
             write(runid,'(i4.4,a1)') iter, 'A'
          end select

       endif



  end subroutine set_runid


end module var4d
