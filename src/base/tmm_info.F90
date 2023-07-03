!###############################################################################
!
! Info on meteo field:
!   o param name
!   o unit
!   o history
!   o ...
!
! use tmm_info
!
! type(TMeteoInfo)   ::  tmi
!
! ! Initialise with name and unit;
! ! history is empty or optionally filled from existing info's.
! call Init( tmi, 'name', 'unit', status, (/tmi1,tmi2,.../) )
!
! ! .. or init as copy of existing info:
! call Init( tmi, tmi2, status )
!
! ! Add text to history:
! call AddHistory( tmi, 'type==od', status )
!
! ! Add existing history to history:
! call AddHistory_tmi( tmi, tmi2, status )
!
! ! extract fields:
! call Get( tmi, status, name=name, unit=unit, history=history )
!
! ! ok
! call Done( tmi, status )
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module tmm_info

  use os_specs, only : HISTORY_LEN

  implicit none


  ! --- in/out --------------------------------------

  private

  public   ::  TMeteoInfo
  public   ::  Init, Done, Get
  public   ::  SetHistory
  public   ::  AddHistory


  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'tmm_info'


  ! --- types --------------------------------------

  type TMeteoInfo
    character(len=16)    ::  name=''
    character(len=16)    ::  unit=''
    character(len=HISTORY_LEN)  ::  history=''
  end type TMeteoInfo


  ! --- interfaces --------------------------------

  interface Init
    module procedure tmi_Init
    module procedure tmi_Init_copy
  end interface

  interface Done
    module procedure tmi_Done
  end interface

  interface Get
    module procedure tmi_Get
  end interface

  interface AddHistory
    module procedure tmi_AddHistory
    module procedure tmi_AddHistory_tmi
  end interface

  interface SetHistory
    module procedure tmi_SetHistory_tmi
  end interface


contains


  ! =============================================================


  subroutine tmi_Init( tmi, name, unit, status, tmis )

    ! --- in/out -----------------------------------

    type(TMeteoInfo), intent(out)            ::  tmi
    character(len=*), intent(in)             ::  name
    character(len=*), intent(in)             ::  unit
    integer, intent(out)                     ::  status

    type(TMeteoInfo), intent(in), optional   ::  tmis(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tmi_Init'

    ! --- local -------------------------------------

    integer         ::  i

    ! --- begin ------------------------------------

    ! store name and unit:
    tmi%name = trim(name)
    tmi%unit = trim(unit)

    ! start with empty history:
    tmi%history = ''

    ! include input histories ?
    if ( present(tmis) ) then
      ! loop over input info's :
      do i = 1, size(tmis)
        call AddHistory( tmi, tmis(i), status )
      end do
    end if

    ! add name and unit to history:
    call AddHistory( tmi, 'name=='//trim(tmi%name), status )
    call AddHistory( tmi, 'unit=='//trim(tmi%unit), status )

    ! ok
    status = 0

  end subroutine tmi_Init


  ! ***


  subroutine tmi_Init_copy( tmi, tmi2, status )

    ! --- in/out -----------------------------------

    type(TMeteoInfo), intent(out)            ::  tmi
    type(TMeteoInfo), intent(in)             ::  tmi2
    integer, intent(out)                     ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tmi_Init_copy'

    ! --- begin ------------------------------------

    ! copy:
    tmi = tmi2

    ! ok
    status = 0

  end subroutine tmi_Init_copy


  ! ***


  subroutine tmi_Done( tmi, status )

    ! --- in/out -----------------------------------

    type(TMeteoInfo), intent(inout)      ::  tmi
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tmi_Done'

    ! --- begin ------------------------------------

    ! ok
    status = 0

  end subroutine tmi_Done


  ! ***


  subroutine tmi_AddHistory( tmi, history, status )

    ! --- in/out -----------------------------------

    type(TMeteoInfo), intent(inout)           ::  tmi
    character(len=*), intent(in)              ::  history
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tmi_AddHistory'

    ! --- begin ------------------------------------

    ! add item to history, close with ';;'
    if ( len_trim(tmi%history)+len_trim(history)+2 < len(tmi%history) ) then
      tmi%history = trim(tmi%history)//trim(history)//';;'
    else
      !write (*,'("WARNING - history buffer too small; increase size in ",a)') mname
    end if

    ! ok
    status = 0

  end subroutine tmi_AddHistory


  ! ***


  subroutine tmi_AddHistory_tmi( tmi, tmi2, status )

    ! --- in/out -----------------------------------

    type(TMeteoInfo), intent(inout)           ::  tmi
    type(TMeteoInfo), intent(in)              ::  tmi2
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tmi_AddHistory_tmi'

    ! --- begin ------------------------------------

    ! add extra ';;', represents a new line:
    call AddHistory( tmi, trim(tmi2%history)//';;', status )

    ! ok
    status = 0

  end subroutine tmi_AddHistory_tmi


  ! ***


  subroutine tmi_SetHistory_tmi( tmi, tmi2, status )

    ! --- in/out -----------------------------------

    type(TMeteoInfo), intent(inout)           ::  tmi
    type(TMeteoInfo), intent(in)              ::  tmi2
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tmi_SetHistory_tmi'

    ! --- begin ------------------------------------

    ! replace existing history:
    tmi%history = tmi2%history

    ! ok
    status = 0

  end subroutine tmi_SetHistory_tmi


  ! ***


  subroutine tmi_Get( tmi, status, name, unit, history )

    ! --- in/out -----------------------------------

    type(TMeteoInfo), intent(in)              ::  tmi
    integer, intent(out)                      ::  status

    character(len=*), intent(out), optional   ::  name
    character(len=*), intent(out), optional   ::  unit
    character(len=*), intent(out), optional   ::  history

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tmi_Set'

    ! --- begin ------------------------------------

    ! fill name ?
    if ( present(name) ) name = trim(tmi%name)

    ! fill unit ?
    if ( present(unit) ) unit = trim(tmi%unit)

    ! fill history ?
    if ( present(history) ) history = trim(tmi%history)

    ! ok
    status = 0

  end subroutine tmi_Get



end module tmm_info
