!###############################################################################
!
! NAME
!   go_SystemMod  -  machine and/or compiler specific stuff
!
! DESCRIPTION
!   The module go_SystemMod provides some basic constants for the
!   current compiler. In addition, some interfaces are defined
!   to routines for system calls, setting of exit statuses etc, 
!   which are non-standard Fortran, but often provided by the
!   vendor of the compiler. 
!   Since both constants and system routines differ from compiler
!   to compiler, this GO module is available in a number of copies,
!   each valid for a single compiler. If for some compiler a 
!   certain constant or system routine could not be filled,
!   a dummy value is used or a warning is issued.
!
!   The following system routines are defined:
!
!    o call goSystem( command, status )
!        Perform a system command, return exit status.
!
!    o call goExit( status )
!        Stop execution, set the exit status.
!
!    o call goArgCount( narg, status )
!        Count number of command line arguments.
!
!    o call goGetArg( nr, value, status )
!        Returns command line argument 'nr' in character string 'value'.
!
!    o call goSleep( nsec, status )
!        Wait for some seconds.
!
! ONLINE MANUALS
!
!   Latest known locations:
!
!   GNU gfortan
!       gcc.gnu.org/onlinedocs/
!
!   Intel Fortran Compiler
!       www.intel.com
!         Sitemap, Software, Find Product : Intel Compilers
!           Select 'Product Documentation' from the side menu
!       http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/index.htm
!    
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
!###############################################################################

module GO_System

  use GO_Print, only : gol, goPr, goErr
  
  implicit none

  ! --- in/out ------------------------------
  
  private
  
  public   ::  goSystem
  public   ::  goExit
  public   ::  goArgCount, goGetArg
  public   ::  goSleep
  
  public   ::  pathsep
  public   ::  goPathJoin
  
  
  ! --- const ---------------------------------
  
  ! module name
  character(len=*), parameter  ::  mname = 'GO_System'
  
  ! path seperation:
  character(len=1), parameter  ::  pathsep = '/'
    

contains



  ! ############################################################################
  ! ###
  ! ### goSystem
  ! ###
  ! ############################################################################


  ! Execute a system command, return exit status.
  
  subroutine goSystem( command, status )
  
#ifdef __INTEL_COMPILER
    use IFPort, only : System
    use IFPort, only : iErrNo, E2BIG, ENOENT, ENOEXEC, ENOMEM
#endif
  
    ! --- in/out -----------------------------------------------
    
    character(len=*), intent(in)       ::  command
    integer, intent(inout)             ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goSystem'
    
    ! --- local --------------------------------------------------
    
#ifdef __INTEL_COMPILER
    integer(4)     ::  stat
    integer(4)     ::  errno
#endif

    ! --- begin --------------------------------------------------
    
#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    stat = System( command )
    
    ! trap errors in ifort system command
    if ( stat == -1 ) then
      write (gol,'("error in call to IFort Portability command `system`:")'); call goErr
      errno = iErrNo()
      select case ( errno )
        case ( E2BIG )
          write (gol,'("  ",a)') 'The argument list is too long.'; call goErr
        case ( ENOENT )
          write (gol,'("  ",a)') 'The command interpreter cannot be found.'; call goErr
        case ( ENOEXEC )
          write (gol,'("  ",a)') 'The command interpreter file has an invalid format and is not executable.'; call goErr
        case ( ENOMEM )
          write (gol,'("  ",a)') 'Not enough system resources are available to execute the command.'; call goErr
        case default
          write (gol,'("  unknown iErrNo ",i)') errno; call goErr
      end select
      TRACEBACK; status=stat; return
    end if
    
    ! if the shell command exit status is 'n',
    ! then the number returned by 'system' is 256 * n
    status = stat / 256

#else
#ifdef __GFORTRAN__ 

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GFortran Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call System( command, status )

#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !write (gol,'("could not evaluate system command : ",a)') trim(command); call goErr
    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    !TRACEBACK; status=1; return

    ! just try ...
    call System( command, status )

#endif
#endif
  
  end subroutine goSystem
  
  
  ! ############################################################################
  ! ###
  ! ### goExit
  ! ###
  ! ############################################################################


  ! Stop execution, set exit status.
  
  subroutine goExit( status )

#ifdef __IBMC__
    use XLFUtility, only : Exit_
#endif
  
#ifdef __INTEL_COMPILER
    use IFPort, only : Exit
#endif

    ! --- in/out --------------------------------------------
    
    integer, intent(in)    ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goExit'
    
    ! --- begin --------------------------------------------
    
#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call Exit( int(status,kind=4) )

#else
#ifdef __GFORTRAN__
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GFortran compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call Exit( int(status,kind=4) )
    
#else
!#ifdef __xlc__
!    
!    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    ! IBM XLF compiler
!    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    call Exit_( int(status,kind=4) )
!    
!#else
!
!    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    ! error ...
!    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    
!    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
!    
!    !! this is an emergency, so for one time, the Fortran stop is allowed ...
!    !stop 'Fortran STOP in GO_System/goExit'
!
!    ! try this one ...
!    call Exit( int(status,kind=4) )
!
!#endif
    ! just stop, do not mind the error status ...
    stop
#endif
#endif
  
  end subroutine goExit


  ! ############################################################################
  ! ###
  ! ### goArgC
  ! ###
  ! ############################################################################


  ! Return number of command line arguments
  
  subroutine goArgCount( narg, status )
  
    ! --- in/out --------------------------------------------
    
    integer, intent(out)    ::  narg
    integer, intent(out)    ::  status

    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goArgCount'
    
    ! --- external ----------------------------------------------

#ifdef __INTEL_COMPILER
!    ! intrinsic function
!    interface
!      integer(4) function nArgs()
!      end function nArgs
!    end interface
#else
    integer, external  ::  iArgC
#endif

    ! --- begin -------------------------------------------------
    
#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !! isn't this the default ?
    !narg = iArgC()
    
    !! from the IFPort user manual:
    !narg = nArgs()
    
#if __INTEL_COMPILER == 1110
    ! from the online manual:
    narg = Command_Argument_Count()
#else    
    ! try this, often works:
    narg = iArgC()
#endif
    
#else


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !! always assign something ...
    !narg = -1

    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    !TRACEBACK; status=1; return
    
    ! try this, often works:
    narg = iArgC()

#endif

    ! ok
    status = 0
  
  end subroutine goArgCount
  
  
  ! ############################################################################
  ! ###
  ! ### goGetArg
  ! ###
  ! ############################################################################



  ! Return a command line argument
  
  subroutine goGetArg( pos, value, status )

    ! --- in/out --------------------------------------------------
    
    integer, intent(in)             ::  pos
    character(len=*), intent(out)   ::  value
    integer, intent(inout)          ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goGetArg'
    
    ! --- external ------------------------------------------------

!#ifdef __INTEL_COMPILER
!    ! intrinsic function
!    interface
!      subroutine GetArg( n, buffer, status )
!        integer, intent(in)               ::  n
!        character(len=*), intent(out)     ::  buffer
!        integer, intent(out), optional    ::  status
!      end subroutine GetArg
!    end interface
!#endif
    
    ! --- local -----------------------------------------------------
    
    integer     ::  n
    
    ! --- begin -----------------------------------------------------

#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#if __INTEL_COMPILER == 1110
    ! from the online manual:
    call Get_Command_Argument( pos, value=value, status=status )
#else    
    ! call system routine:    
    call GetArg( pos, value, status )
    ! Number of characters set in buffer is returned in status:
    if ( status <= 0 ) then
      status = 1
    else
      status = 0
    end if
#endif
    
#else
#ifdef __GFORTRAN__

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! call system routine:    
    call GetArg( pos, value )
    ! no status returned ...
    status = 0
    
#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !! use arguments to avoid compilation warnings:
    !status = pos
    !! always assign something ...
    !value = '?'
    
    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    !TRACEBACK; status=1; return

    ! try this ...
    call GetArg( pos, value )
    status = 0

#endif
#endif
  
  end subroutine goGetArg
  
  
  ! ############################################################################
  ! ###
  ! ### goSleep
  ! ###
  ! ############################################################################


  ! wait some seconds ...
  
  subroutine goSleep( nsec, status )

#ifdef __INTEL_COMPILER
    use IFPort, only : Sleep
#endif

    ! --- in/out --------------------------------------------------
    
    integer, intent(in)             ::  nsec
    integer, intent(out)            ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goSleep'
    
    ! --- begin -----------------------------------------------------

#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! call system routine:    
    call Sleep( nsec )

#else
#ifdef __GFORTRAN__

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! gfortran Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! call system routine:    
    call Sleep( nsec )

#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! use arguments to avoid compilation warnings:
    status = nsec
    
    write (gol,'("subroutine not implemented for this compiler")'); call goErr
    TRACEBACK; status=1; return

#endif
#endif

    ! ok
    status = 0
  
  end subroutine goSleep
  
  
  ! ***
  
  ! join two paths :  p = p1/p2
  
  subroutine goPathJoin( p1, p2, p, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)    ::  p1, p2
    character(len=*), intent(out)   ::  p
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goPathJoin'
    
    ! --- begin ----------------------------------
    
    ! combine paths:
    write (p,'(a,a,a)') trim(p1), pathsep, trim(p2)
    
    ! ok
    status = 0
    
  end subroutine goPathJoin
  
  


end module GO_System



