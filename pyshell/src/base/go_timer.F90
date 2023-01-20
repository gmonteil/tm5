!#######################################################################
!
! NAME
!   GO_Timer  -  General Objects : Timing routines
!
! USAGE
!
!  use GO_Timer
!
!  ! timer id's:
!  integer     ::  itim1, itim2, itim2a, itim2b
!
!  ! start timing:
!  call GO_Timer_Init( status )
!
!  ! define timer names, return timer id's:
!  call GO_Timer_Def( itim1 , 'part1' , status )
!  call GO_Timer_Def( itim2 , 'part2' , status )
!  call GO_Timer_Def( itim2a, 'part2a', status )
!  call GO_Timer_Def( itim2b, 'part2b', status )
!
!  ! first task:
!  call GO_Timer_Start(itim1,status)
!  ! ...
!  call GO_Timer_End(itim1,status)
!
!  ! second task:
!  call GO_Timer_Start(itim2,status)
!  ! ...
!    ! child tasks:
!    call GO_Timer_Start(itim2a,status)
!    ! ...
!    call GO_Timer_End(itim2a,status)
!    call GO_Timer_Start(itim2b,status)
!    ! ...
!    call GO_Timer_End(itim2b,status)
!  call GO_Timer_End(itim2,status)
!
!  ! stop timing, print profile to standard output;
!  ! if an output file name is provided, the timing data is written
!  ! to this file with the profile in the header:
!  call GO_Timer_Done( status  [,'profile.dat'] )
!
!
! HISTORY
!
!   2008 apr, Arjo Segers, TNO
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
! code compiled together with other GO modules ...
#define with_go
!
!#######################################################################

module GO_Timer

#ifdef with_go
  use GO_Print, only : gol, goPr, goErr
#endif

  implicit none


  ! --- in/out ------------------------

  private

  public    ::  GO_Timer_Init, GO_Timer_Done
  public    ::  GO_Timer_Def, GO_Timer_Start, GO_Timer_End
  public    ::  GO_Timer_Get


  ! --- const --------------------------

  character(len=*), parameter   ::  mname = 'GO_Timer'

  ! maximum number of times:
  integer, parameter    ::  maxtimer = 500

  ! real kind returned by cpu_time etc
  integer, parameter    ::  rknd = 8

  ! integer kind returned by system_clock etc
  integer, parameter    ::  iknd = 4

  ! --- types --------------------------

  type T_Timer
    ! label:
    character(len=64)   ::  name
    ! total time:
   !real(rknd)          ::  total_cpu
    real(rknd)          ::  total_sys
  end type T_Timer

  type T_Stopwatch
    !! timing using 'cpu_time' routine:
    !real(rknd)          ::  start_cpu
    !real(rknd)          ::  end_cpu
    !real(rknd)          ::  total_cpu
    ! timing using 'system_clock' routine:
    integer(iknd)       ::  start_sys
    integer(iknd)       ::  end_sys
    real(rknd)          ::  total_sys
  end type T_Stopwatch


  ! --- var ----------------------------

  ! list of timers:
  type(T_Timer)       ::  Timers(0:maxtimer)

  ! currently in use:
  integer             ::  ntimer

  ! root timer:
  integer             ::  itim_root

  ! parent-child relations:
  !logical             ::  child(0:maxtimer,maxtimer)

  ! StopWatch for each parent/child pair:
  type(T_Stopwatch)   ::  StopWatch(0:maxtimer,maxtimer)

  ! stack of current timers:
  integer             ::  stack(0:maxtimer)
  integer             ::  top

  ! parameters of system_clock :
  integer(iknd)       ::  sysclock_count_rate       ! clock ticks per second
  integer(iknd)       ::  sysclock_count_max        ! maximum number of ticks
  real(rknd)          ::  sysclock_tick2sec

#ifndef with_go
  ! message line:
  use os_specs, only : LONG_STR_LEN
  character(len=LONG_STR_LEN) ::  gol
#endif


contains


#ifndef with_go

  ! ********************************************************************
  ! ***
  ! *** GO surrogate
  ! ***
  ! ********************************************************************

  ! substitutes for message routines from GO modules

  ! display message:
  subroutine goPr
    write (*,'(a)') trim(gol)
  end subroutine goPr

  ! display error message:
  subroutine goErr
    write (*,'("ERROR - ",a)') trim(gol)
  end subroutine goErr

  ! free file unit:
  subroutine goGetFU( fu, status )
    integer, intent(out)    ::  fu
    integer, intent(out)    ::  status
    logical                 ::  opened
    fu = 456
    do
      inquire( unit=fu, opened=opened )
      if ( .not. opened ) exit
      fu = fu + 1
    end do
    status = 0
  end subroutine goGetFU

#endif


  ! ********************************************************************
  ! ***
  ! *** GO Timer Routines
  ! ***
  ! ********************************************************************


  subroutine GO_Timer_Init( status )

    ! --- in/out -------------------------

    integer, intent(out)      ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Timer_Init'

    ! --- local --------------------------

    integer(iknd)   ::  sysclock_count
    integer         ::  itimer, ichild

    ! --- begin --------------------------

    ! init system clock parameters:
    call system_clock( sysclock_count, sysclock_count_rate, sysclock_count_max )
    ! conversion from clock ticks to seconds:
    sysclock_tick2sec = 1.0/real(sysclock_count_rate,8)

    ! no timers defined yet:
    ntimer = 0

    ! dummy name for base, which might be used as parent:
    Timers(0)%name = '0'

    ! no children yet:
    !child = .false.

    ! no StopWatch yet:
    do itimer = 0, ntimer
      do ichild = 1, ntimer
        ! set accumulated time to zero:
       !StopWatch(itimer,itimer)%total_cpu = 0.0
        StopWatch(itimer,itimer)%total_sys = 0.0
      end do
    end do

    ! empty stack:
    stack = 0
    top = 0

    ! define root timer:
    call GO_Timer_Def( itim_root, 'root', status )
    IF_NOTOK_RETURN(status=1)

    ! start root:
    call GO_Timer_Start( itim_root, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine GO_Timer_Init


  ! ***


  subroutine GO_Timer_Done( status, file, verbose )

#ifdef with_go
    use GO_File, only : goGetFU
#endif

    ! --- in/out -------------------------

    integer, intent(out)                      ::  status
    character(len=*), intent(in), optional    ::  file
    logical, intent(in), optional             ::  verbose

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Timer_Done'

    ! --- external ---------------------------------

#ifdef _OPENMP
    integer, external  ::  OMP_GET_THREAD_NUM
    integer, external  ::  OMP_GET_NUM_THREADS
#endif

    ! --- local --------------------------

    logical             ::  do_verbose
    integer             ::  itimer, ichild
    character(len=40)   ::  label, child_label
    real(rknd)          ::  total, child_total
    real(rknd)          ::  children_total
    real                ::  frac
    integer             ::  fu
    real(rknd)          ::  child_totals(maxtimer)

    ! --- begin --------------------------

    ! arguments:
    do_verbose = .false.
    if ( present(verbose) ) do_verbose = verbose

    ! stop root:
    call GO_Timer_End( itim_root, status )
    IF_NOTOK_RETURN(status=1)

    ! also to file ?
    if ( present(file) ) then
      ! free file unit:
      call goGetFU( fu, status )
      IF_NOTOK_RETURN(status=1)
      ! open file:
      open( fu, file=trim(file), form='formatted', iostat=status )
      if (status/=0) then
        write (gol,'("opening timer output file : ",a)') trim(file); call goPr
        TRACEBACK; status=1; return
      end if
    end if

    ! print table
    if ( do_verbose ) then
      write (gol,'(" ")'); call goPr
      write (gol,'("------------------------------------------  ------------ ---------")'); call goPr
      write (gol,'("timer                                       system_clock       (%)")'); call goPr
      write (gol,'("------------------------------------------  ------------ ---------")'); call goPr
    end if
    ! also to file ?
    if ( present(file) ) then
      write (fu,'("#")')
      write (fu,'("# ------------------------------------------  ------------ ---------")')
      write (fu,'("# timer                                       system_clock       (%)")')
      write (fu,'("# ------------------------------------------  ------------ ---------")')
    end if

    ! loop over all timers:
    do itimer = 1, ntimer

      ! current values:
      label = trim(timers(itimer)%name)
     !total = timers(itimer)%total_cpu
      total = timers(itimer)%total_sys

      ! display ?
      if ( do_verbose ) then
        write (gol,'(" ")'); call goPr
        write (gol,'(a40,"  ",1("  ",f12.2,"          "))') label, total; call goPr
      end if
      ! also to file ?
      if ( present(file) ) then
        write (fu,'("#")')
        write (fu,'("# ",a40,"  ",1("  ",f12.2,"          "))') label, total
      end if

      ! loop over children:
     !children_total = 0.0
      children_total = 0.0
      do ichild = 1, ntimer

        ! child values:
        child_label = trim(timers(ichild)%name)
       !child_total = StopWatch(itimer,ichild)%total_cpu
        child_total = StopWatch(itimer,ichild)%total_sys

        ! no time spend here ? then skip:
        if ( child_total <= 0.0 ) cycle

        ! set fraction:
        if ( total > 0.0 ) then
          frac = child_total / total
        else
          frac = 1.0
        endif

        ! display?
        if ( do_verbose ) then
          write (gol,'("  ",a40,1("  ",f12.2," (",f5.1," %)"))') child_label, child_total, frac*100.0; call goPr
        end if
        ! also to file ?
        if ( present(file) ) then
          write (fu,'("#   ",a40,1("  ",f12.2," (",f5.1," %)"))') child_label, child_total, frac*100.0
        end if

        ! update sum:
        children_total = children_total + child_total

      end do ! child

      ! other ?
      if ( children_total > 0.0 ) then

        ! 'child' values:
        child_label = 'other'
        child_total = total - children_total

        ! check ...
        if ( child_total < 0.0 ) then
          ! tell the user to check the code ...
          write (gol,'("WARNING - total of children exceeds time spent by parent, probably a wrong start/end pair somewhere!")')
          ! next timer:
          cycle
        end if

        ! set fraction:
        if ( total > 0.0 ) then
          frac = child_total / total
        else
          frac = 1.0
        endif

        ! display:
        if ( do_verbose ) then
          write (gol,'("  ",a40,1("  ",f12.2," (",f5.1," %)"))') child_label, child_total, frac*100.0; call goPr
        end if
        ! also to file ?
        if ( present(file) ) then
          write (fu,'("#   ",a40,1("  ",f12.2," (",f5.1," %)"))') child_label, child_total, frac*100.0
        end if

      end if

    end do  ! timers

    ! close table:
    if ( do_verbose ) then
      write (gol,'(" ")'); call goPr
      write (gol,'("------------------------------------------  ------------ ---------")'); call goPr
      write (gol,'(" ")'); call goPr
    end if
    ! also to file ?
    if ( present(file) ) then
      write (fu,'("#")')
      write (fu,'("# ------------------------------------------  ------------ ---------")')
      write (fu,'("#")')
    end if

    ! write all data to the file:
    if ( present(file) ) then
      ! format number:
      write (fu,'("# format number:")')
      write (fu,'("[format]")')
      write (fu,'("2.0")')
      ! parallel:
      write (fu,'("# number of OpenMP thread; zero if not enabled:")')
      write (fu,'("[number-of-threads]")')
#ifdef _OPENMP
      !$OMP PARALLEL
      if ( OMP_GET_THREAD_NUM() == 0 ) then
        write (fu,*) OMP_GET_NUM_THREADS()
      end if
      !$OMP END PARALLEL
#else
      write (fu,*) 0
#endif
      ! numbers:
      write (fu,'("# number of timers:")')
      write (fu,'("[number-of-timers]")')
      write (fu,*) ntimer
      ! all data:
      write (fu,'("# index, total time, name")')
      write (fu,'("[timer-table]")')
      do itimer = 1, ntimer
        write (fu,'(i4,f12.4," ",a)') itimer, timers(itimer)%total_sys, trim(timers(itimer)%name)
      end do
      ! parent-child table:
      write (fu,'("# for each timer, total times spent on child processes")')
      write (fu,'("[timer-matrix]")')
      do itimer = 1, ntimer
        ! collect child totals:
        child_totals = 0.0
        do ichild = 1, ntimer
          child_totals(ichild) = StopWatch(itimer,ichild)%total_sys
        end do
        write (fu,'(1000f12.4)') child_totals(1:ntimer)
      end do
    end if

    ! close file if necessary:
    if ( present(file) ) then
      ! close:
      close( fu, iostat=status )
      if (status/=0) then
        write (gol,'("closing timer output file : ",a)') trim(file); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! ok
    status = 0

  end subroutine GO_Timer_Done


  ! ***


  subroutine GO_Timer_Def( itimer, name, status )

    ! --- in/out -------------------------

    integer, intent(out)            ::  itimer
    character(len=*), intent(in)    ::  name
    integer, intent(out)            ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Timer_Def'

    ! --- local --------------------------

    integer     ::  k

    ! --- begin --------------------------

    ! check if already defined;
    ! loop over currently defined timers:
    do k = 1, ntimer
      ! check ...
      if ( trim(name) == trim(timers(k)%name) ) then
        ! info ...
        write (gol,'("WARNING - timer `",a,"` already defined ...")') trim(name); call goPr
        ! set output handle:
        itimer = k
        ! ok
    status = 0
        return
      end if
    end do

    ! new number:
    ntimer = ntimer + 1

    ! check ...
    if ( ntimer > maxtimer ) then
      write (gol,'("could not define timer for `",a,"` ;")') trim(name); call goPr
      write (gol,'("reached maximum number of timers:")'); call goPr
      do k = 1, maxtimer
        write (gol,'("  ",i6," ",a)') k, trim(timers(k)%name); call goPr
      end do
      write (gol,'("increase value of parameter `maxtimer` in module `",a,"`")') trim(mname); call goPr
      TRACEBACK; status=1; return
    end if

    !! debug ...
    !print *, 'TTT def timer : ', ntimer, ' ', trim(name)

    ! current number:
    itimer = ntimer

    ! store:
    timers(itimer)%name = trim(name)

    ! init totals:
   !timers(itimer)%total_cpu = 0.0
    timers(itimer)%total_sys = 0.0

    ! ok:
    status = 0

  end subroutine GO_Timer_Def


  ! ***


  subroutine GO_Timer_Get( itimer, status, name )

    ! --- in/out -------------------------

    integer, intent(in)             ::  itimer
    integer, intent(out)            ::  status
    character(len=*), optional      ::  name

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Timer_Name'

    ! --- local --------------------------

    ! --- begin --------------------------

    ! extract values
    if ( present(name) ) name = trim(timers(itimer)%name)

    ! ok:
    status = 0

  end subroutine GO_Timer_Get


  ! ***


  subroutine GO_Timer_Start( itimer, status )

    ! --- in/out -------------------------

    integer, intent(in)             ::  itimer
    integer, intent(out)            ::  status

    ! --- local --------------------------

    integer       ::  i
    integer       ::  iparent

    ! --- begin --------------------------

    ! check ...
    if ( itimer < 1 ) then
      write (gol,'("timer id < 1 ; not defined ?")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( top == size(stack) ) then
      write (gol,'("timer stack out of bounds:")'); call goErr
      do i = 1, top
        write (gol,'(i6," : ",i6," `",a,"`")') i, stack(i), trim(Timers(i)%name); call goErr
      end do
      write (gol,'("probably bug in start/end calls, please check ...")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( top < 0 ) then
      write (gol,'("stack could not be lower than zero, but top is now : ",i6)') top; call goErr
      TRACEBACK; status=1; return
    end if

    ! add to stack:
    top = top + 1
    stack(top) = itimer

    ! current timer is on top of stack;
    ! parent code has timer stack(top-1):
    iparent = stack(top-1)

    ! set flag that parent calls this part of the code:
    !child(iparent,itimer) = .true.

    !! store time:
    !call cpu_time( StopWatch(iparent,itimer)%start_cpu )

    ! store ticks:
    call system_clock( StopWatch(iparent,itimer)%start_sys )

    ! ok:
    status = 0

  end subroutine GO_Timer_Start


  ! ***


  subroutine GO_Timer_End( itimer, status )

    ! --- in/out -------------------------

    integer, intent(in)             ::  itimer
    integer, intent(out)            ::  status

    ! --- local --------------------------

    integer          ::  iparent
   !real(rknd)       ::  dt_cpu
    real(rknd)       ::  dt_sys

    ! --- begin --------------------------

    ! check ..
    if ( stack(top) /= itimer ) then
      write (gol,'("end timer id not the same as start timer id:")'); call goErr
      write (gol,'("  start (top of stack) : ",i6," `",a,"`")') stack(top), trim(Timers(stack(top))%name); call goErr
      write (gol,'("  end                  : ",i6," `",a,"`")') itimer, trim(Timers(itimer)%name); call goErr
      write (gol,'("check if each timer start is followed by a correct timer end")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( top < 1 ) then
      write (gol,'("timer end but stack empty ...")'); call goErr
      write (gol,'("check if each call to timer_end has a corresponding call to timer_start")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! current timer is on top of stack;
    ! parent code has timer stack(top-1):
    iparent = stack(top-1)

    !! store time:
    !call cpu_time( stopwatch%end_cpu )
    !! add time increment:
    !dt_cpu = stopwatch%end_cpu - stopwatch%start_cpu
    !! add time increments:
    !StopWatch(iparent,itimer)%total_cpu = StopWatch(iparent,itimer)%total_cpu + dt_cpu
    !Timers   (        itimer)%total_cpu = Timers   (        itimer)%total_cpu + dt_cpu

    ! store time:
    call system_clock( StopWatch(iparent,itimer)%end_sys )
    ! trap reset:
    if ( StopWatch(iparent,itimer)%end_sys < StopWatch(iparent,itimer)%start_sys ) then
      ! set time increment:
      dt_sys = ( StopWatch(iparent,itimer)%end_sys + ( sysclock_count_max - StopWatch(iparent,itimer)%start_sys ) ) * sysclock_tick2sec
    else
      ! set time increment:
      dt_sys = ( StopWatch(iparent,itimer)%end_sys - StopWatch(iparent,itimer)%start_sys ) * sysclock_tick2sec
    end if
    ! add time increments:
    StopWatch(iparent,itimer)%total_sys = StopWatch(iparent,itimer)%total_sys + dt_sys
    Timers   (        itimer)%total_sys = Timers   (        itimer)%total_sys + dt_sys

    ! debugging ...
    !write (*,'("xxx added ",f6.2," to timer `",a,"`; called from `",a,"`")') dt_sys, trim(Timers(itimer)%name), trim(Timers(iparent)%name)

    ! pop from stack:
    top = top - 1

    ! ok:
    status = 0

  end subroutine GO_Timer_End


  ! ***


end module GO_Timer


!! ##########################################################
!! ###
!! ### test
!! ###
!! ##########################################################
!
!program test
!
!  use GO_Timer
!
!  implicit none
!
!  ! timer id's:
!  integer     ::  itim1, itim2, itim2a, itim2b
!
!  ! local:
!  integer     ::  status
!
!  ! start timing:
!  call GO_Timer_Init( status )
!  if (status/=0) stop 'ERROR from GO_Timer_Init'
!
!  ! define timer names, return timer id's:
!  call GO_Timer_Def( itim1 , 'part1' , status )
!  call GO_Timer_Def( itim2 , 'part2' , status )
!  call GO_Timer_Def( itim2a, 'part2a', status )
!  call GO_Timer_Def( itim2b, 'part2b', status )
!
!  ! first task:
!  call GO_Timer_Start(itim1,status)
!  ! ...
!  call Sleep( 1 )
!  ! ...
!  call GO_Timer_End(itim1,status)
!
!  ! second task:
!  call GO_Timer_Start(itim2,status)
!  ! ...
!  call Sleep( 2 )
!  ! ...
!    ! child tasks:
!    call GO_Timer_Start(itim2a,status)
!    ! ...
!    call Sleep( 2 )
!    ! ...
!    call GO_Timer_End(itim2a,status)
!    call GO_Timer_Start(itim2b,status)
!    ! ...
!    call Sleep( 3 )
!    ! ...
!    call GO_Timer_End(itim2b,status)
!  ! ...
!  call Sleep( 1 )
!  ! ...
!  call GO_Timer_End(itim2,status)
!
!  ! stop timing, print profile
!  call GO_Timer_Done( status )
!  if (status/=0) stop 'ERROR from GO_Timer_Done'
!
!end program test
!
