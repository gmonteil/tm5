!
! Print routine labels
!
! Usage:
!
!   use GO_Label
!
!   ! Enable printing of routine labels:
!   call GO_Label_Init( status, trace=.true. )
!   if (status/=0) stop
!
!   ! start new labeled block of messages:
!   call goLabel( 'work in progress' )
!
!   ! end block:
!   call goLabel()
!
!   ! done with labels:
!   call GO_Label_Done( status )
!   if (status/=0) stop
!   
!
! Example:
!
!   The code:
!
!     call goLabel( 'main' )
!     write (gol,'("start of the main routine")'); call goPr
!
!       call goLabel('sub task')
!       write (gol,'("something to be done")'); call goPr
!       call goLabel()
!
!     call goLabel()
!
!   will produce the following output:
!
!     <main>
!       start of the main routine
!       <sub task>
!         something to be done
!       </sub task>
!     </main>
!
!
!#######################################################################
!
#define TRACEBACK write (gol,'("ERROR in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
!#######################################################################

module GO_Label

  use GO_Print, only : gol, goPr, goErr, goBug

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  GO_Label_Init, GO_Label_Done
  public  ::  goLabel


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'GO_Label'


  ! --- var --------------------------------------
  
  ! flags etc
  logical              ::  pr_trace
  
  ! stack with labels:
  integer, parameter   ::  mstack = 400
  character(len=64)    ::  labels(0:mstack)
  integer              ::  timers(0:mstack)
  integer              ::  istack = 0
  

contains


  ! ***


  subroutine GO_Label_Init( status, trace )

    ! --- in/out -------------------------

    integer, intent(out)              ::  status
    logical, intent(in), optional     ::  trace

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Label_Init'

    ! --- local --------------------------

    ! --- begin --------------------------

    ! trace labels ?
    pr_trace = .false.
    if ( present(trace) ) pr_trace = trace
    
    ! init stacks:
    labels(0) = '<no-label>'
    timers(0) = -1
    istack = 0
    
    !! used to be in GO_Print_Init, why was this necessary ?
    !if ( .not. pr_trace ) call GO_Print_DeIndent()
    
    ! ok
    status = 0

  end subroutine GO_Label_Init


  ! ***


  subroutine GO_Label_Done( status )

    ! --- in/out -------------------------

    integer, intent(out)                      ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Label_Done'

    ! --- local --------------------------

    ! --- begin --------------------------

    ! ok
    status = 0

  end subroutine GO_Label_Done


  ! ***************************************************************************
  ! ***
  ! *** routine labels
  ! ***
  ! ***************************************************************************


  subroutine goLabel( label, timer )
  
    use GO_Print, only : GO_Print_Indent, GO_Print_DeIndent
    use GO_Timer, only : GO_Timer_Get, GO_Timer_Start, GO_Timer_End
  
    ! --- in/out -------------------------------
    
    character(len=*), intent(in), optional  ::  label
    integer, intent(in), optional           ::  timer

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goLabel'
    
    ! --- local --------------------------------

    integer  :: status
    
    ! --- begin --------------------------------
    
    ! add new label to stack ?
    if ( present(label) .or. present(timer) ) then

      ! ** new labeled block **
      
      ! increase stack:
      istack = istack + 1
      ! check ...
      if ( istack > mstack ) then
        write (gol,'("BUG - stack too small; please increase mstack in go_label")'); call goBug
        TRACEBACK; stop
      else
        ! add label stack:
        if ( present(label) ) then
          ! copy the provided label:
          labels(istack) = trim(label)
        else if ( present(timer) ) then
          ! copy name:
          call GO_Timer_Get( timer, status, name=labels(istack) )
          IF_NOTOK_RETURN(stop)
        else
          ! dummy ...
          labels(istack) = 'dummy label'
        end if
        ! timer provided ?
        if ( present(timer) ) then
          ! store timer index:
          timers(istack) = timer
        else
          ! dummy ...
          timers(istack) = -1
        end if
      end if
      
      ! print labels ?
      if (pr_trace) then
        write (gol,'("<",a,">")') trim(labels(istack)); call goPr
      end if
      
      ! increase indent for new messages:
      call GO_Print_Indent()
      
      ! start timing ?
      if ( timers(istack) >= 0  ) then
        call GO_Timer_Start( timers(istack), status )
        IF_NOTOK_RETURN(stop)
      end if

    else

      ! ** end of labeled block **
      
      ! stop timing ?
      if ( timers(istack) >= 0 ) then
        call GO_Timer_End( timers(istack), status )
        IF_NOTOK_RETURN(stop)
      end if
      
      ! decrease indent:
      call GO_Print_DeIndent()
      
      ! print end label ?
      if (pr_trace) then
        write (gol,'("</",a,">")') trim(labels(istack)); call goPr
      end if
      
      ! remove from stack:
      istack = max( 0, istack - 1 )

    end if
    
  end subroutine goLabel
    

end module GO_Label



! #############################################################################
! ###
! ### test program
! ###
! #############################################################################
!
!
!module testmod
!
!  implicit none
!  
!  public
!  
!contains
!
!  subroutine subr( i, status )
!  
!    use go_print, only : goLabel, gol, goPr, goErr
!    
!    ! --- in/out ----------------------------------------
!    
!    integer, intent(in)           ::  i
!    integer, intent(out)          ::  status
!    
!    ! --- begin -----------------------------------------
!    
!    call goLabel( 'subr' )
!    
!    write (gol,'("welcome to subr !")'); call goPr
!    
!    select case ( i )
!    
!      case ( 0 )
!        write (gol,'("testing i : ",i2)') i; call goPr
!      
!      case ( 1 )
!        call subr2( 0, status )
!        if (status/=0) then; call goErr; status=1; return; end if
!      
!      case ( 2 )
!        call subr2( 1, status )
!        if (status/=0) then; call goErr; status=1; return; end if
!      
!      case default
!        write (gol,'("unsupported i : ",i2)') i; call goErr
!        call goErr; status=1; return
!        
!    end select
!        
!    call goLabel(); status=0
!    
!  end subroutine subr
!    
!
!  ! ***
!
!    
!  subroutine subr2( i, status )
!  
!    use go_print, only : goLabel, gol, goPr, goErr
!    
!    ! --- in/out ----------------------------------------
!    
!    integer, intent(in)           ::  i
!    integer, intent(out)          ::  status
!    
!    ! --- begin -----------------------------------------
!    
!    call goLabel('subr2')
!
!    write (gol,'("testing subr2")'); call goPr
!    
!    select case ( i )
!      case ( 0 )
!      case default
!        write (gol,'("wrong i : ",i2)') i; call goErr
!        call goErr; status=1; return
!    end select
!    
!    call goLabel; status=0
!    
!  end subroutine subr2
!    
!
!
!end module testmod
!
!
! ################################################################
!
!
!program test
!
!  use go_label
!  use testmod
!  
!  ! --- local -----------------------------------------
!  
!  integer           ::  status
!  
!  ! --- begin ------------------------------------------
!
!  call GO_Print_Init( status, trace=.false. )
!  call goLabel('test prog')
!  
!  write (gol,'("begin of program")'); call goPr
!  
!  call Subr( 2, status )
!  if (status/=0) then; call goErr; call exit(1); end if
!  
!  write (gol,'("end of program")'); call goPr
!
!  call goLabel()
!  call GO_Print_Done( status )
!
!end program test
!

