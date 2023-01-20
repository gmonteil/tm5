!###############################################################################
!
! go_date  -  module to manipulate date structures
!
! TYPES
!
!###############################################################################
!
! go_date  -  module to manipulate date structures
!
! TYPES
!
!  A structure is provided to store a date:
!
!   ! declare date types:
!   type(TDate)     ::  t0, t1, t, dt
!   type(TIncrDate) ::  dt
!
!  with fields:
!
!    character(len=4)   ::  calender     ! see 'CALENDERS'
!
!    integer        ::  year, month, day, hour, min, sec, mili
!
!    integer        ::  zone       ! minutes; add this to obtain GMT
!
!
! CALENDERS
!
!   A number of different calender types is supported:
!
!     'wall'   : wall clock time, including time zone
!
!     'greg'   : Gregorian calender, some years have a Februari 29
!     '366'    : every year has a Februari 29
!     '365'    : a year has never a  Februari 29
!     '360'    : every month has 30 days
!
!     'incr'   : incremental time step: year=0, month=0, day >= 0
!
!   The 'incr' type is a special calender which has no year
!    or month but might have any number of days.
!   Note that day==1 has the interpretation of 24 hours for an 'incr',
!   but means 'first' or 0 hours for one of the regular calenders.
!
!   Use the calender '360' if only operations on years and months are required.
!
!
! CREATING DATE STRUCTURES
!
!   To initialize a new date structure, a few routines are available.
!
!   Use routine 'NewDate' to initialize some fields and to fill
!   the rest with zero's. If no calender is specified,
!   the default value 'greg' is used (see also DEFAULTS).
!
!     t = NewDate( calender='greg', year=2000, month=1, ... )
!
!   Use routine 'IncrDate' to create a new increment:
!
!     dt = IncrDate( year=2000, month=1 )
!
!   Fill the time from the system clock in a date structure:
!
!     t = go_SystemDate()
!
!
! FIELD MANIPULATION
!
!   Use 'Set' to fill some specific fields of a date structure.
!   Special arrays:
!     time4 = (/year,month,day,hour/)
!     time5 = (/year,month,day,hour,min/)
!     time6 = (/year,month,day,hour,min,sec/)
!   Example:
!
!     call Set( t [,year=2000] [,month=1] [,day=2] ... &
!                 [,time4=time4] [,time5=time5] [,time6=time6])
!
!   Use 'Get' to obtain some specific fields of a date structure.
!
!     call Get( t [,year=year] [,month=month] ... &
!                 [,time4=time4] [,time5=time5] [,time6=time6] )
!
!   Check contents of a date structure:
!
!     call Check( t )
!
!   Normalize hours to {0,..,23}, minutes to {0,..,59}, etc:
!
!     call Normalize( t )
!
!   Return minimum or maximum of two dates:
!
!     tmin = min( t1, t2 )
!     tmax = max( t1, t2 )
!
!   Retrun begin|end of year|month|day specified by t :
!
!     t1 = Get_Begin_Of( t, 'year'  )    ! 00:00 at 01 Jan of this year
!     t1 = Get_Begin_Of( t, 'month' )    ! 00:00 at first day of this month
!     t1 = Get_Begin_Of( t, 'day'   )    ! 00:00 at this day
!
!     t2 = Get_End_Of( t1, 'year'  )    ! 00:00 at 01 Jan of next year
!     t2 = Get_End_Of( t1, 'month' )    ! 00:00 at first day of next month
!     t2 = Get_End_Of( t1, 'day'   )    ! 00:00 at next day
!
!
! INQUIRY FUNCTIONS
!
!   A few inquiry functions are provided.
!
!   The logical function 'LeapYear' tells you if the year
!   has a Februari 29 :
!
!     l = LeapYear( t )
!
!   Two integer functions are provided to count the total number
!   of days in a month or a year:
!
!     i = Days_in_Month( t )
!     i = Days_in_Year( t )
!
!   An integer function is provided to return the day number,
!   counting from 1 (Januari 1) to 360, 365, or 366 (last of December):
!
!     i = DayNumber( t )
!
!   The logical function 'Midnight' is true for times 00:00 :
!
!     l = Midnight( t )
!
!   The logical function 'Precisely' is true for times
!   that are precisely 1.5 hour etc:
!
!     l = Precisely( t,  1.5, 'hour' )    ! every 1.5 hour
!     l = Precisely( t, 24.0, 'hour' )    ! midnight
!
!
! OPERATORS
!
!   Operators '+' and '-' are redefined to perform operations
!   between two date structures.
!   Both should be of same calender type, or one should be
!   an increment:
!
!     t = t1 + t2
!     t = t1 - t2
!
!   Operators '*' and '/' are redefined for multiplication with
!   or division by a real or an integer:
!
!     t = t1 + dt * 2
!     t = t1 + dt * 3.1415
!     t = t1 + dt / 3.1415
!
!
! LOGICAL OPERATORS
!
!   Operators '==', '/=', '<', '<=', '>', '>=' are defined to
!   compare two dates.
!
!
! SUMMATION ROUTINES
!
!   The total number in a certain unit is returned by 'rTotal'
!   (real value) or 'iTotal' (integer, error if not possible).
!   Currently supported units are 'year', 'month', 'day',
!   'hour', 'min', 'sec', and 'mili'. If the total number is
!   not wel defined for a certain date, for example the
!   total number of years of today, an error message is produced.
!
!     r = rTotal( t, 'year'|'month'|... )
!     i = iTotal( t, 'year'|'month'|... )
!
!
! INTERPOLATION
!
!   For t in [t1,t2], return real coefficients alfa1 and alf2 such that:
!       t  =  alfa1 * t1  +  alfa2 * t2
!   Usefull for linear interpolation:
!       f(t)  ~  alfa1 * f(t1)  +  alfa2 * f(t2)
!
!     call InterpolFractions( t, t1, t2, alfa1, alfa2, status )
!
!   Return interval [tt(1),tt(2)] around t that matches with time resolution;
!   resolution specified by a step and unit:
!      3, 'hour'   # 00:00, 03:00, 06:00, ...
!
!     call Get_Surrounding_Interval( t, 3, 'hour', tt, status )
!
!
! OUTPUT
!
!   To obtain a pretty formatted print of the value of a date,
!   the 'Pretty' routine is provided. Output differs based on
!   the calender type.
!
!     print *, 't = '//trim(Pretty(t))
!
!   Special routines are available using the GO_Print facility:
!
!     call wrtgol( 'time             = ', t ); call goPr
!     call wrtgol( 'time range       = ', t1, ' - ', t2 ); call goPr
!     call wrtgol( 'time range (ref) = ', t1, ' - ', t2, ';', t0 ); call goPr
!
!
! DEFAULTS
!
!   For setting some default values, the subroutine 'go_DateDefaults'
!   is available. All arguments are optional:
!
!     call go_DateDefaults( [calender='greg'] )
!
!
!### macro's ##################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define PRINT_LEN 200
!
!###############################################################################


module GO_Date

  use GO_Print, only : gol, goErr, goPr

  implicit none

  ! --- in/out ---------------------------

  private

  public      ::  TDate, TIncrDate

  public      ::  goDateDefaults

  public      ::  NewDate, IncrDate, AnyDate, SystemDate
  public      ::  Get_Begin_Of, Get_End_Of

  public      ::  Set, Get
  public      ::  Check
  public      ::  Normalize

  public      ::  LeapYear, LeapDay
  public      ::  days_in_month
  public      ::  day_of_week
  public      ::  days_in_year
  public      ::  DayNumber
  public      ::  SecondNumber
  public      ::  Midnight
  public      ::  Precisely

  public      ::  operator(+)
  public      ::  operator(-)
  public      ::  operator(*)
  public      ::  operator(/)

  public      ::  IsAnyDate
  public      ::  isSameDay, isSameMonth
  public      ::  operator(==)
  public      ::  operator(/=)
  public      ::  operator(>)
  public      ::  operator(<)
  public      ::  operator(>=)
  public      ::  operator(<=)

  public      ::  assignment(=)

  public      ::  min, max
  public      ::  date_min, date_max

  public      ::  rTotal, iTotal

  public      ::  InterpolFractions
!  public      ::  IntervalFractions

  public      ::  Get_Surrounding_Interval

  public      ::  wrtgol
  public      ::  Pretty

  ! --- const -----------------------------------

  character(len=*), parameter  ::  mname = 'GO_Date'

  ! --- types -------------------------------------

  ! Strucure with fields to store year, month, day,
  ! hour and minute.
  ! Operators for assignment (=), adding (+),
  ! and comparission (==,<,>,>= and <=)
  ! have been defined for operations between
  ! instances of this type.

  type TDate
    ! type of calender: 'greg', '365', '360'
    character(len=4)    ::  calender
    ! year, month etc:
    integer             ::  year, month, day, hour, min, sec, mili
    ! difference with Coordinated Universal Time (UTC)
    integer             ::  zone   ! minutes
    ! error status
    integer             ::  status ! = 1
  end type TDate


  type TIncrDate
    ! days, hours, etc:
    integer             ::  day, hour, min, sec, mili
    ! error status
    integer             ::  status ! = 1
  end type TIncrDate

  ! testing ...
  type DateTime_Interval
    ! boundaries:
    type(TDate)         ::  left, right
    ! open ? default is closed:
    logical             ::  left_open, right_open
  end type DateTime_Interval


  ! --- var --------------------------------

  ! default calender type
  character(len=4)      ::  default_calender = 'greg'


  ! --- interface ---------------------------

  interface Check
    module procedure date_Check
    module procedure incrdate_Check
  end interface

  interface Pretty
    module procedure date_Pretty
    module procedure incrdate_Pretty
  end interface

  ! *

  interface LeapYear
    module procedure date_LeapYear
  end interface

  interface LeapDay
    module procedure date_LeapDay
  end interface

  interface days_in_month
    module procedure date_days_in_month
  end interface

  interface days_in_year
    module procedure date_days_in_year
    module procedure calc_days_in_year_greg
  end interface

  interface SecondNumber
    module procedure calc_SecondNumber_i
    module procedure calc_SecondNumber_r
  end interface

  interface DayNumber
    module procedure date_DayNumber
  end interface

  interface Midnight
    module procedure date_Midnight
  end interface

  interface Precisely
    module procedure date_Precisely
  end interface

  ! *

  interface Set
    module procedure date_Set
    module procedure incrdate_Set
  end interface

  interface Get
    module procedure date_Get
    module procedure incrdate_Get
  end interface

  ! *

  interface NewDate
    module procedure date_NewDate
  end interface

  interface AnyDate
    module procedure date_AnyDate
  end interface

  interface IncrDate
    module procedure incrdate_IncrDate
  end interface

  interface SystemDate
    module procedure date_SystemDate
  end interface

  interface Get_Begin_Of
    module procedure date_Get_Begin_Of
  end interface

  interface Get_End_Of
    module procedure date_Get_End_Of
  end interface

  ! * operators

  interface Normalize
    module procedure date_Normalize
    module procedure incrdate_Normalize
  end interface

  interface operator(+)
    module procedure t_plus_t
    module procedure t_plus_dt
    module procedure dt_plus_dt
  end interface

  interface operator(-)
    module procedure t_min_t
    module procedure t_min_dt
    module procedure dt_min_dt
  end interface

  interface operator(*)
    module procedure dt_times_r
    module procedure  r_times_dt
    module procedure dt_times_i
    module procedure  i_times_dt
  end interface

  interface operator(/)
    module procedure dt_div_r
    module procedure dt_div_i
  end interface

  ! assignments

  interface assignment(=)
    module procedure date_set_date
    module procedure incrdate_set_incrdate
  end interface

  ! * logical operators

  interface IsAnyDate
    module procedure date_IsAnyDate
  end interface

  interface operator(==)
    module procedure date_eq_date
  end interface

  interface operator(/=)
    module procedure date_ne_date
  end interface

  interface operator(>)
    module procedure date_gt_date
  end interface

  interface operator(<)
    module procedure date_lt_date
  end interface

  interface operator(>=)
    module procedure date_ge_date
  end interface

  interface operator(<=)
    module procedure date_le_date
  end interface

  ! *

  interface min
    module procedure date_min
  end interface

  interface max
    module procedure date_max
  end interface

  ! *

  interface rTotal
    module procedure date_rTotal
    module procedure incr_rTotal
  end interface

  interface iTotal
    module procedure date_iTotal
    module procedure incrdate_iTotal
  end interface

  ! *

  interface InterpolFractions
    module procedure date_InterpolFractions
  end interface

!  interface date_IntervalFractions
!    module procedure date_IntervalFractions
!  end interface

  ! *

  interface wrtgol
    module procedure wrtgol_t
    module procedure wrtgol_dt
    module procedure wrtgol_t1_t2
    module procedure wrtgol_t1_t2_t3
  end interface



contains


  ! ****************************************************
  ! ***
  ! *** set defaults
  ! ***
  ! ****************************************************


  subroutine goDateDefaults( calender )

    ! --- in/out --------------------------------

    character(len=*), intent(in), optional    ::  calender

    ! --- begin ----------------------------------

    if ( present(calender) ) default_calender = calender

  end subroutine goDateDefaults


  ! ****************************************************
  ! ***
  ! *** check
  ! ***
  ! ****************************************************

  !
  ! Check fields of a date:
  !   range etc
  !

  subroutine date_Check( t, status )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------------

    type(TDate), intent(in)        ::  t
    integer, intent(out)           ::  status

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Check'

    ! --- begin -----------------------------------

    ! already error status ? then leave immediatelly:
    if ( t%status /= 0 ) then
      write (gol,'("found error status in date")'); call goErr
      write (gol,'("  year,month,day : ",3i6)') t%year, t%month, t%day; call goErr
      write (gol,'("  hour,minu,sec,mili : ",4i6)') t%hour, t%min, t%sec, t%mili; call goErr
      TRACEBACK; status=1; return
    end if

    ! calender specific
    select case ( t%calender )
      case ( 'any' )
        ! always ok ...
        status = 0
        return
      case ( 'wall' )
        ! no special tests
      case ( 'greg', '366', '365', '360' )
        ! check month
        if ( t%month<1 .or. t%month>12 ) then
          call wrtgol( 'strange month in ', t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! check day
        if ( t%day<1 .or. t%day>days_in_month(t) ) then
          call wrtgol( 'strange day in ', t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! zone should be zero:
        if ( t%zone /= 0 ) then
          call wrtgol( 'expecting zero zone in date ', t ); call goErr
          TRACEBACK; status=1; return
        end if
      case default
        write (gol,'("unknown calender type: `",a,"`")') t%calender; call goErr
        write (gol,'("  year etc : ",6i5)') t%year, t%month, t%day, t%hour, t%min, t%sec; call goErr
        TRACEBACK; status=1; return
    end select

    ! check minutes
    if ( t%min<0 .or. t%min>59 ) then
      call wrtgol( 'found strange minutes in ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! check seconds
    if ( t%sec<0 .or. t%sec>59 ) then
      call wrtgol( 'found strange seconds in ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! check mili
    if ( t%mili<0 .or. t%mili>999 ) then
      call wrtgol( 'found strange mili seconds in ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine date_Check


  ! ***


  subroutine incrdate_Check( dt, status )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------------

    type(TIncrDate), intent(in)    ::  dt
    integer, intent(out)           ::  status

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_Check'

    ! --- begin -----------------------------------

    ! already error status ? then leave immediatelly:
    if ( dt%status /= 0 ) then
      write (gol,'("found error status in incrdate")'); call goErr
      write (gol,'("  day, hour,minu,sec,mili : ",5i6)') dt%day, dt%hour, dt%min, dt%sec, dt%mili; call goErr
      TRACEBACK; status=1; return
    end if

    ! every value is allowed for increments ...

    ! ok
    status = 0

  end subroutine incrdate_Check



  ! ****************************************************
  ! ***
  ! *** computation
  ! ***
  ! ****************************************************

  integer function day_of_week(t)
    ! Monday is 1, Tuesday is 2, etc.
    use go_print,   only : gol, goErr

    ! --- in/out -------------------------------
    type(TDate), intent(in)     :: t

    ! --- const -----------------------------
    character(len=*), parameter :: rname = mname//'/day_of_week'

    ! --- local -----------------------------
    type(TDate)                 :: ref_date
    type(TIncrDate)             :: dt

    ! --- begin ----------------

    ref_date = NewDate(time6=(/1900, 1, 1, 0, 0, 0/), calender='greg') ! 1900-01-01 was a Monday, or day number 1

    if (t < ref_date) then
        write (gol,'("Date must be on or after ", a)') trim(Pretty(ref_date)) ; call goErr
        TRACEBACK; stop
    end if

    dt = t - ref_date
    day_of_week = mod(dt%day, 7) + 1

  end function day_of_week

  ! Does this year have a 29 feb ?

  logical function calc_LeapYear( year )

    ! --- in/out -------------------------------

    integer, intent(in)           ::  year

    ! --- begin --------------------------------

    calc_LeapYear = ( (mod(year,4)==0) .and. .not.(mod(year,100)==0) ) &
                                        .or. (mod(year,400)==0)

  end function calc_LeapYear


  ! ***


  ! days per month

  integer function calc_days_in_month( calender, year, month )

    use GO_Print, only : gol, goErr

    ! --- in/out ---------------------------

    character(len=*), intent(in)   ::  calender
    integer, intent(in)            ::  year, month

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_days_in_month'

    ! --- const -----------------------------

    ! days in a month                      1  2  3  4  5  6  7  8  9 10 11 12
    integer, parameter :: days365(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! normal
    integer, parameter :: days366(12) = (/31,29,31,30,31,30,31,31,30,31,30,31/) ! leap year
    integer, parameter :: days360(12) = (/30,30,30,30,30,30,30,30,30,30,30,30/) ! fixed month

    ! --- begin ----------------

    select case ( calender )
      case ( 'wall', 'greg' )
        if ( calc_LeapYear(year) ) then
          calc_days_in_month = days366(month)
        else
          calc_days_in_month = days365(month)
        end if
      case ( '366' )
        calc_days_in_month = days366(month)
      case ( '365' )
        calc_days_in_month = days365(month)
      case ( '360' )
        calc_days_in_month = days360(month)
      case ( 'any' )
        calc_days_in_month = 0
      case default
        calc_days_in_month = -1
        write (gol,'("unknown calender type: ",a)') calender; call goErr
        TRACEBACK; stop
    end select

  end function calc_days_in_month


  ! ***

  subroutine calc_SecondNumber_i( t, seconds )
    ! number of seconds since the start of the day

    ! --- in/out ----------------
    type(TDate), intent(in)     :: t
    integer, intent(out)        :: seconds

    ! --- const -----------------------------
    character(len=*), parameter ::  rname = mname//'/calc_SecondNumber_i'

    ! --- begin ----------------

    seconds = t%hour * 3600 + t%min * 60 + t%sec

  end subroutine calc_SecondNumber_i

  subroutine calc_SecondNumber_r( t, seconds )
    ! number of seconds since the start of the day

    ! --- in/out ----------------
    type(TDate), intent(in)     :: t
    real, intent(out)           :: seconds

    ! --- const -----------------------------
    character(len=*), parameter ::  rname = mname//'/calc_SecondNumber_r'

    ! --- begin ----------------

    seconds = t%hour * 3600.0 + t%min * 60.0 + t%sec + t%mili/1000.0

  end subroutine calc_SecondNumber_r

  ! days per year

  integer function calc_days_in_year( calender, year )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------

    character(len=*), intent(in)   ::  calender
    integer, intent(in)            ::  year

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_days_in_year'

    ! --- begin ----------------

    select case ( calender )
      case ( 'wall', 'greg' )
        if ( calc_LeapYear(year) ) then
          calc_days_in_year = 366
        else
          calc_days_in_year = 365
        end if
      case ( '366' )
        calc_days_in_year = 366
      case ( '365' )
        calc_days_in_year = 365
      case ( '360' )
        calc_days_in_year = 360
      case ( 'any' )
        calc_days_in_year = 0
      case default
        write (gol,'("unknown calender type: ",a)') calender; call goErr
        TRACEBACK; stop
    end select

  end function calc_days_in_year


  ! by default the gregorian calender:

  integer function calc_days_in_year_greg( year )

    ! --- in/out ----------------

    integer, intent(in)            ::  year

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_days_in_year_greg'

    ! --- begin ----------------

    ! extract for Gregorian calender:
    calc_days_in_year_greg = calc_days_in_year( 'greg', year )

  end function calc_days_in_year_greg


  ! ***


  ! Returns the number of the day spedified by the date iy/im/id.
  ! The existence of a february 29 is checked.
  !
  ! ndays( 1995,  1,  1 ) =   1
  ! ndays( 1995, 12, 31 ) = 365
  ! ndays( 1996, 12, 31 ) = 366        29 feb every 4 year ...
  ! ndays( 1900, 12, 31 ) = 365         except every 100 year ...
  ! ndays( 2000, 12, 31 ) = 366          except every 400 year ...

  integer function calc_DayNumber( calender, year, month, day )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    character(len=*), intent(in)   ::  calender
    integer, intent(in)            ::  year, month, day

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_DayNumber'

    ! --- local -----------------------------

    integer            ::  imonth

    ! --- begin ----------------------------

    select case ( calender )
      case ( 'wall', 'greg', '366', '365', '360' )
        calc_DayNumber = day
        do imonth = 1, month-1
          calc_DayNumber = calc_DayNumber + calc_days_in_month(calender,year,imonth)
        end do
      case ( 'any' )
        calc_DayNumber = 0
      case default
        write (gol,'("unknown calender type: ",a)') calender; call goErr
        TRACEBACK; stop
    end select

  end function calc_DayNumber


  ! **********************************************


  logical function date_LeapYear( t )

    ! --- in/out -------------------------------

    type(TDate), intent(in)         ::  t

    ! --- begin --------------------------------

    ! calender specific
    select case ( t%calender )
      case ( 'wall', 'greg' )
        ! see above ...
        date_LeapYear = calc_LeapYear( t%year )
      case default
        ! no leap years ...
        date_LeapYear = .false.
    end select

  end function date_LeapYear


  ! ***


  ! is 29 feb ?

  logical function date_LeapDay( t )

    ! --- in/out -------------------------------

    type(TDate), intent(in)         ::  t

    ! --- begin --------------------------------

    ! calender specific
    select case ( t%calender )
      case ( 'wall', 'greg' )
        ! see above ...
        date_LeapDay = (t%month == 2) .and. (t%day == 29)
      case default
        ! no leap years ...
        date_LeapDay = .false.
    end select

  end function date_LeapDay


  ! ***


  integer function date_days_in_month( t )

    ! --- in/out ----------------

    type(TDate), intent(in)         ::  t

    ! --- begin ----------------

    date_days_in_month = calc_days_in_month( t%calender, t%year, t%month )

  end function date_days_in_month


  ! ***


  integer function date_days_in_year( t )

    ! --- in/out ----------------

    type(TDate), intent(in)         ::  t

    ! --- begin ----------------

    date_days_in_year = calc_days_in_year( t%calender, t%year )

  end function date_days_in_year


  ! ***


  integer function date_DayNumber( t )

    ! --- in/out ----------------

    type(TDate), intent(in)        ::  t

    ! --- begin ----------------

    date_DayNumber = calc_DayNumber( t%calender, t%year, t%month, t%day )

  end function date_DayNumber


  ! ***


  logical function date_Midnight( t )

    ! --- in/out -------------------------------

    type(TDate), intent(in)         ::  t

    ! --- begin --------------------------------

    date_Midnight = all( (/t%hour,t%min,t%sec,t%mili/) == 0 )

  end function date_Midnight


  ! ***


  logical function date_Precisely( t, val, unit )

    ! --- in/out -------------------------------

    type(TDate), intent(in)         ::  t
    real, intent(in)                ::  val
    character(len=*), intent(in)    ::  unit

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Precisely'

    ! --- local --------------------------------

    real      ::  rtot

    ! --- begin --------------------------------

    ! only some units supported yet:
    select case ( unit )

      ! every .. minute
      case ( 'min', 'minu' )
        ! check ...
        if ( (val <= 0.0) .or. (val > 60.0) ) then
          write (gol,'("strange value for minu : ",f8.2)') val; call goErr
          TRACEBACK; stop
        end if
        ! total minute fraction in this day:
        rtot = rTotal( IncrDate(0,t%hour,t%min,t%sec,t%mili), 'min' )
        ! modulo given value ?
        date_Precisely = modulo(rtot,val) == 0.0

      ! every .. hour
      case ( 'hour' )
        ! check ...
        if ( (val <= 0.0) .or. (val > 24.0) ) then
          write (gol,'("strange value for hour : ",f8.2)') val; call goErr
          TRACEBACK; stop
        end if
        ! total hour fraction in this day:
        rtot = rTotal( IncrDate(0,t%hour,t%min,t%sec,t%mili), 'hour' )
        ! modulo given value ?
        date_Precisely = modulo(rtot,val) == 0.0

      ! error ...
      case default
        write (gol,'("unsupported unit : ",a)') unit; call goErr
        TRACEBACK; stop

    end select

  end function date_Precisely



  ! ****************************************************
  ! ***
  ! *** Set/Get fields in a structure TDate
  ! ***
  ! ****************************************************

  !
  ! Fill fields in a 'TDate' structure.
  !

  subroutine date_Set( date, year, month, day, hour, min, sec, mili, &
                             zone, calender, time4, time5, time6 )

    ! --- in/out ------------------------------------

    type(TDate), intent(inout)               ::  date
    integer, intent(in), optional            ::  year
    integer, intent(in), optional            ::  month
    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour
    integer, intent(in), optional            ::  min
    integer, intent(in), optional            ::  sec
    integer, intent(in), optional            ::  mili
    integer, intent(in), optional            ::  zone
    character(len=*), intent(in), optional   ::  calender
    integer, intent(in), optional            ::  time4(4)
    integer, intent(in), optional            ::  time5(5)
    integer, intent(in), optional            ::  time6(6)

    ! --- local ----------------------------------

    if ( present(calender) ) date%calender  = calender

    if ( present(time4) ) then
      date%year  = time4(1)
      date%month = time4(2)
      date%day   = time4(3)
      date%hour  = time4(4)
    end if

    if ( present(time5) ) then
      date%year  = time5(1)
      date%month = time5(2)
      date%day   = time5(3)
      date%hour  = time5(4)
      date%min   = time5(5)
    end if

    if ( present(time6) ) then
      date%year  = time6(1)
      date%month = time6(2)
      date%day   = time6(3)
      date%hour  = time6(4)
      date%min   = time6(5)
      date%sec   = time6(6)
    end if

    if ( present(year ) ) date%year  = year
    if ( present(month) ) date%month = month
    if ( present(day  ) ) date%day   = day
    if ( present(zone ) ) date%zone  = zone
    if ( present(hour ) ) date%hour  = hour
    if ( present(min  ) ) date%min   = min
    if ( present(sec  ) ) date%sec   = sec
    if ( present(mili ) ) date%mili  = mili

  end subroutine date_Set


  ! *

  subroutine incrdate_Set( date, day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    type(TIncrDate), intent(inout)             ::  date
    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour
    integer, intent(in), optional            ::  min
    integer, intent(in), optional            ::  sec
    integer, intent(in), optional            ::  mili

    ! --- local ----------------------------------

    if ( present(day  ) ) date%day   = day
    if ( present(hour ) ) date%hour  = hour
    if ( present(min  ) ) date%min   = min
    if ( present(sec  ) ) date%sec   = sec
    if ( present(mili ) ) date%mili  = mili

  end subroutine incrdate_Set


  !
  ! Obtain fields from a 'TDate' structure.
  !

  subroutine date_Get( date, &
                       year, month, day, hour, min, sec, mili, &
                       zone, calender, time4, time5, time6 )

    ! --- in/out ------------------------------------

    type(TDate), intent(in)                  ::  date
    integer, intent(out), optional           ::  year, month, day
    integer, intent(out), optional           ::  hour, min, sec, mili
    integer, intent(out), optional           ::  zone
    integer, intent(out), optional           ::  time4(4)
    integer, intent(out), optional           ::  time5(5)
    integer, intent(out), optional           ::  time6(6)
    character(len=*), intent(out), optional  ::  calender

    ! --- local ----------------------------------

    if ( present(calender) ) calender = date%calender

    if ( present(year)  ) year  = date%year
    if ( present(month) ) month = date%month
    if ( present(day  ) ) day   = date%day
    if ( present(zone ) ) zone  = date%zone
    if ( present(hour ) ) hour  = date%hour
    if ( present(min  ) ) min   = date%min
    if ( present(sec  ) ) sec   = date%sec
    if ( present(mili ) ) mili  = date%mili

    if ( present(time4) ) time4 = (/ date%year, date%month, date%day, date%hour /)

    if ( present(time5) ) time5 = (/ date%year, date%month, date%day, &
                                     date%hour, date%min   /)

    if ( present(time6) ) time6 = (/ date%year, date%month, date%day, &
                                     date%hour, date%min  , date%sec /)

  end subroutine date_Get

  ! *

  subroutine incrdate_Get( date, day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    type(TIncrDate), intent(in)              ::  date
    integer, intent(out), optional           ::  day
    integer, intent(out), optional           ::  hour, min, sec, mili

    ! --- local ----------------------------------

    if ( present(day  ) ) day   = date%day
    if ( present(hour ) ) hour  = date%hour
    if ( present(min  ) ) min   = date%min
    if ( present(sec  ) ) sec   = date%sec
    if ( present(mili ) ) mili  = date%mili

  end subroutine incrdate_Get


  ! ****************************************************
  ! ***
  ! *** Return a date structure
  ! ***
  ! ****************************************************

  !
  ! Set fields to zero or fill some of them.
  !

  type(TDate) function date_NewDate( year, month, day, hour, min, sec, mili, &
                                     zone, calender, time4, time5, time6 )

    ! --- in/out ------------------------------------

    integer, intent(in), optional            ::  year, month, day
    integer, intent(in), optional            ::  hour, min, sec, mili
    integer, intent(in), optional            ::  zone
    character(len=*), intent(in), optional   ::  calender
    integer, intent(in), optional            ::  time4(4)
    integer, intent(in), optional            ::  time5(5)
    integer, intent(in), optional            ::  time6(6)

    ! --- local ----------------------------------

    ! set default calender type:
    date_NewDate%calender = default_calender

    ! Fields are zero by default:
    date_NewDate%year  = 0
    date_NewDate%month = 0
    date_NewDate%day   = 0
    date_NewDate%zone  = 0
    date_NewDate%hour  = 0
    date_NewDate%min   = 0
    date_NewDate%sec   = 0
    date_NewDate%mili  = 0

    ! Optionally, change some of them:
    if ( present(year    ) ) call Set( date_NewDate, year=year )
    if ( present(month   ) ) call Set( date_NewDate, month=month )
    if ( present(day     ) ) call Set( date_NewDate, day=day )
    if ( present(hour    ) ) call Set( date_NewDate, hour=hour )
    if ( present(min     ) ) call Set( date_NewDate, min=min )
    if ( present(sec     ) ) call Set( date_NewDate, sec=sec )
    if ( present(mili    ) ) call Set( date_NewDate, mili=mili )
    if ( present(zone    ) ) call Set( date_NewDate, zone=zone )
    if ( present(calender) ) call Set( date_NewDate, calender=calender )
    if ( present(time4   ) ) call Set( date_NewDate, time4=time4 )
    if ( present(time5   ) ) call Set( date_NewDate, time5=time5 )
    if ( present(time6   ) ) call Set( date_NewDate, time6=time6 )

    ! normalize too small/too large values:
    if ( date_NewDate%year /= 0000 ) call Normalize( date_NewDate )

    ! data filled, thus probably no error ...
    date_NewDate%status = 0

  end function date_NewDate


  ! ***


  type(TIncrDate) function incrdate_IncrDate( day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour, min, sec, mili

    ! --- local ----------------------------------

    ! Fields are zero by default:
    incrdate_IncrDate%day   = 0
    incrdate_IncrDate%hour  = 0
    incrdate_IncrDate%min   = 0
    incrdate_IncrDate%sec   = 0
    incrdate_IncrDate%mili  = 0

    ! Optionally, change some of them:
    if ( present(day     ) ) call Set( incrdate_IncrDate, day=day )
    if ( present(hour    ) ) call Set( incrdate_IncrDate, hour=hour )
    if ( present(min     ) ) call Set( incrdate_IncrDate, min=min )
    if ( present(sec     ) ) call Set( incrdate_IncrDate, sec=sec )
    if ( present(mili    ) ) call Set( incrdate_IncrDate, mili=mili )

    ! normalize too small/too large values:
    call Normalize( incrdate_IncrDate )

    ! data filled, thus probably no error ...
    incrdate_IncrDate%status = 0

  end function incrdate_IncrDate


  ! ***


  !
  ! Set fields to zero, special calender
  !

  type(TDate) function date_AnyDate()

    ! --- local ----------------------------------

    ! Set some fields, other are automatically zero:
    date_AnyDate = NewDate( calender='any' )

  end function date_AnyDate


  ! ***


  ! Fill with system time

  type(TDate) function date_SystemDate()

    ! --- in/out ------------------------------

    ! none ...

    ! --- local ------------------------------

    integer           ::  values(8)

    ! --- begin ------------------------------

    !
    ! Optional character output of Date_and_Time:
    !
    !   date   '20020812'
    !   time   '211757.314'
    !   zone   '+0200'
    !

    ! obtain system date and time:
    call Date_and_Time( values=values )

    ! fill fields in structure:
    call Set( date_SystemDate, calender='wall', &
                year=values(1), month=values(2), day=values(3), &
                zone=values(4), hour=values(5), &
                min=values(6), sec=values(7), mili=values(8) )

    ! no error ...
    date_SystemDate%status = 0

  end function date_SystemDate


  ! ***


  ! set to begin of : 'year', 'month', 'day'

  type(TDate) function date_Get_Begin_Of( t, what )

    ! --- in/out ------------------------------

    type(TDate), intent(in)           ::  t
    character(len=*), intent(in)      ::  what

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Get_Begin_Of'

    ! --- local ------------------------------

    ! --- begin ------------------------------

    ! end of what ?
    select case ( what )

      case ( 'year' )

        ! begin of year:
        date_Get_Begin_Of = NewDate( t%year, 01, 01, 00, 00, 00 )

      case ( 'month' )

        ! begin of month:
        date_Get_Begin_Of = NewDate( t%year, t%month, 01, 00, 00, 00 )

      case ( 'day' )

        ! begin of day:
        date_Get_Begin_Of = NewDate( t%year, t%month, t%day, 00, 00, 00 )

      case default

        write (gol,'("end of `",a,"` not supported yet")') what; call goErr
        TRACEBACK; stop

    end select

    ! no error ...
    date_Get_Begin_Of%status = 0

  end function date_Get_Begin_Of


  ! ***


  ! set to end of : 'month', ...

  type(TDate) function date_Get_End_Of( t, what )

    ! --- in/out ------------------------------

    type(TDate), intent(in)           ::  t
    character(len=*), intent(in)      ::  what

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Get_End_Of'

    ! --- local ------------------------------

    ! --- begin ------------------------------

    ! first set to begin:
    date_Get_End_Of = Get_Begin_Of( t, what )

    ! end of what ?
    select case ( what )
      case ( 'year' )
        ! add number of days:
        date_Get_End_Of = date_Get_End_Of + IncrDate(day=Days_in_Year(t))
      case ( 'month' )
        ! add number of days:
        date_Get_End_Of = date_Get_End_Of + IncrDate(day=Days_in_Month(t))
      case ( 'day' )
        ! add single day:
        date_Get_End_Of = date_Get_End_Of + IncrDate(day=1)
      case default
        write (gol,'("end of `",a,"` not supported yet")') what; call goErr
        TRACEBACK; stop
    end select

    ! no error ...
    date_Get_End_Of%status = 0

  end function date_Get_End_Of


  ! ************************************************
  ! ***
  ! *** operators
  ! ***
  ! ************************************************


  subroutine date_Normalize( t )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(inout)       ::  t

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Normalize'

    ! --- begin ---------------------------------

    ! mili seconds
    do
      if ( t%mili >= 0 ) exit
      t%sec = t%sec - 1
      t%mili = t%mili + 1000
    end do
    do
      if ( t%mili <= 999 ) exit
      t%mili = t%mili - 1000
      t%sec = t%sec + 1
    end do

    ! seconds
    do
      if ( t%sec >= 0 ) exit
      t%min = t%min - 1
      t%sec = t%sec + 60
    end do
    do
      if ( t%sec <= 59 ) exit
      t%sec = t%sec - 60
      t%min = t%min + 1
    end do

    ! minutes
    do
      if ( t%min >= 0 ) exit
      t%hour = t%hour - 1
      t%min = t%min + 60
    end do
    do
      if ( t%min <= 59 ) exit
      t%min = t%min - 60
      t%hour = t%hour + 1
    end do

    ! hours
    do
      if ( t%hour >= 0 ) exit
      t%day = t%day - 1
      t%hour = t%hour + 24
    end do
    do
      if ( t%hour <= 23 ) exit
      t%hour = t%hour - 24
      t%day = t%day + 1
    end do

    ! days, months, year
    select case ( t%calender )
      case ( 'wall', 'greg', '366', '365', '360' )
        do
          if ( t%day >= 1 ) exit
          t%month = t%month - 1
          do
            if ( t%month >= 1 ) exit
            t%year = t%year - 1
            t%month = t%month + 12
          end do
          t%day = t%day + days_in_month(t)
        end do
        do
          if ( t%day <= days_in_month(t) ) exit
          t%day = t%day - days_in_month(t)
          t%month = t%month + 1
          do
            if ( t%month <= 12 ) exit
            t%month = t%month - 12
            t%year = t%year + 1
          end do
        end do
      case default
        write (gol,'("unsupported calender type: ",a)') t%calender; call goErr
        TRACEBACK; stop
    end select

  end subroutine date_Normalize


  subroutine incrdate_Normalize( dt )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TIncrDate), intent(inout)    :: dt

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_Normalize'

    ! --- begin ---------------------------------

    ! mili seconds
    do
      if ( dt%mili >= 0 ) exit
      dt%sec = dt%sec - 1
      dt%mili = dt%mili + 1000
    end do
    do
      if ( dt%mili <= 999 ) exit
      dt%mili = dt%mili - 1000
      dt%sec = dt%sec + 1
    end do

    ! seconds
    do
      if ( dt%sec >= 0 ) exit
      dt%min = dt%min - 1
      dt%sec = dt%sec + 60
    end do
    do
      if ( dt%sec <= 59 ) exit
      dt%sec = dt%sec - 60
      dt%min = dt%min + 1
    end do

    ! minutes
    do
      if ( dt%min >= 0 ) exit
      dt%hour = dt%hour - 1
      dt%min = dt%min + 60
    end do
    do
      if ( dt%min <= 59 ) exit
      dt%min = dt%min - 60
      dt%hour = dt%hour + 1
    end do

    ! hours
    do
      if ( dt%hour >= 0 ) exit
      dt%day = dt%day - 1
      dt%hour = dt%hour + 24
    end do
    do
      if ( dt%hour <= 23 ) exit
      dt%hour = dt%hour - 24
      dt%day = dt%day + 1
    end do

  end subroutine incrdate_Normalize


  ! *** date = t1 + t2 ************************


  !
  !        t1             +   t2             ->   t1+t2
  !
  !        greg               incr                greg
  !        366                incr                366
  !        365                incr                365
  !
  !        360                360                 360
  !        360                incr                360
  !
  !        incr               greg                greg
  !        incr               366                 366
  !        incr               365                 365
  !        incr               360                 360
  !        incr               incr                incr
  !

  type(TDate) function t_plus_t( t1, t2 )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_plus_t'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t1, status )
    IF_NOTOK_STOP
    call Check( t2, status )
    IF_NOTOK_STOP

    ! any date ? return any date ..
    if ( (t1%calender == 'any') .or. (t2%calender == 'any') ) then
      t_plus_t = AnyDate()
      return
    end if

    ! calenders should be the same:
    if ( t1%calender /= t2%calender ) then
      write (gol,'("calenders should be the same : ")'); call goPr
      write (gol,'("  t1 : ",a)') trim(t1%calender); call goPr
      write (gol,'("  t2 : ",a)') trim(t2%calender); call goPr
      TRACEBACK; stop
    end if

    ! add all fields;
    t_plus_t = NewDate( calender=t1%calender, &
                        year  = t1%year  + t2%year  , &
                        month = t1%month + t2%month , &
                        day   = t1%day   + t2%day   , &
                        hour  = t1%hour  + t2%hour  , &
                        zone  = t1%zone  + t2%zone  , &
                        min   = t1%min   + t2%min   , &
                        sec   = t1%sec   + t2%sec   , &
                        mili  = t1%mili  + t2%mili     )

  end function t_plus_t


  ! *


  type(TDate) function t_plus_dt( t, dt )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t
    type(TIncrDate), intent(in)   ::  dt

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_plus_dt'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t, status )
    IF_NOTOK_STOP
    call Check( dt, status )
    IF_NOTOK_STOP

    ! any date ? return any date ..
    if ( t%calender == 'any' ) then
      t_plus_dt = AnyDate()
      return
    end if

    ! add fields; normalization is applied in routine:
    t_plus_dt = NewDate( calender = t%calender, &
                         year     = t%year             , &
                         month    = t%month            , &
                         day      = t%day   + dt%day   , &
                         hour     = t%hour  + dt%hour  , &
                         zone     = t%zone             , &
                         min      = t%min   + dt%min   , &
                         sec      = t%sec   + dt%sec   , &
                         mili     = t%mili  + dt%mili     )

  end function t_plus_dt


  ! *


  type(TIncrDate) function dt_plus_dt( dt1, dt2 )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt1
    type(TIncrDate), intent(in)   ::  dt2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_plus_dt'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( dt1, status )
    IF_NOTOK_STOP
    call Check( dt2, status )
    IF_NOTOK_STOP

    ! add fields:
    dt_plus_dt = IncrDate( day  = dt1%day  + dt2%day   , &
                           hour = dt1%hour + dt2%hour  , &
                           min  = dt1%min  + dt2%min   , &
                           sec  = dt1%sec  + dt2%sec   , &
                           mili = dt1%mili + dt2%mili     )

  end function dt_plus_dt


  ! *** date = t1 - t2

  !
  !        t1        ->  t2        ->   t1-t2      action
  !
  !        greg          greg           incr       difference
  !        greg          incr           greg       minus
  !
  !        366           366            incr       difference
  !        366           incr           366        minus
  !
  !        365           365            incr       difference
  !        365           incr           365        minus
  !
  !        360           360            360        difference
  !        360           incr           360        minus
  !
  !        incr          incr           incr       difference
  !


  type(TIncrDate) function t_min_t( t1, t2 )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_t'

    ! --- local ---------------------------------

    integer             ::  status
    integer             ::  ndays
    type(TDate)         ::  t

    ! --- begin ---------------------------------

    ! difference between two dates;
    ! algorithm implemented for positive differences only:

    ! difference should be positive:
    if ( t1 > t2 ) then

      ! call work routine:
      t_min_t = t_min_t_work( t1, t2 )

    else

      ! negative
      t_min_t = t_min_t_work( t2, t1 ) * (-1.0)

    end if

    call Normalize(t_min_t)

  end function t_min_t


  ! *

  type(TIncrDate) function t_min_t_work( t1, t2 )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_t_work'

    ! --- local ---------------------------------

    integer             ::  status
    integer             ::  ndays
    type(TDate)         ::  t

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t1, status )
    IF_NOTOK_STOP
    call Check( t2, status )
    IF_NOTOK_STOP

    ! any dates ? something wrong ...
    if ( (t1%calender == 'any') .or. (t2%calender == 'any') ) then
      write (gol,'("do not know how to compute difference between `any` dates ...")')
      TRACEBACK; stop
    end if

    ! calenders should be the same:
    if ( t1%calender /= t2%calender ) then
      write (gol,'("calenders should be the same : ")'); call goPr
      write (gol,'("  t1 : ",a)') trim(t1%calender); call goPr
      write (gol,'("  t2 : ",a)') trim(t2%calender); call goPr
      TRACEBACK; stop
    end if

    ! difference between two dates;
    ! algorithm implemented for positive differences only:

    ! difference should be positive:
    if ( t1 < t2 ) then
      write (gol,'("expect t1 to exceed t2 :")'); call goErr
      call wrtgol( '  t1 : ', t1 ); call goErr
      call wrtgol( '  t2 : ', t2 ); call goErr
      TRACEBACK; stop
    end if

    ! determine number of days between t1 and t2:
    t = t1
    ndays = daynumber(t) - 1
    do
      if ( t%year==t2%year ) exit
      t%year = t%year - 1
      ndays = ndays + days_in_year(t)
    end do
    ndays = ndays - (daynumber(t2)-1)

    ! store result:
    t_min_t_work = IncrDate( day  = ndays, &
                        hour = t1%hour - t2%hour, &
                        min  = t1%min  - t2%min , &
                        sec  = t1%sec  - t2%sec , &
                        mili = t1%mili - t2%mili           )

  end function t_min_t_work


  ! *


  type(TDate) function t_min_dt( t, dt )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t
    type(TIncrDate), intent(in)   ::  dt

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_dt'

    ! --- local ---------------------------------

    integer             ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t, status )
    IF_NOTOK_STOP
    call Check( dt, status )
    IF_NOTOK_STOP

    ! any date ? return any date ..
    if ( t%calender == 'any' ) then
      t_min_dt = AnyDate()
      return
    end if

    ! result is of same type as t
    ! ~ first guess:
    t_min_dt = NewDate( calender = t%calender      , &
                        year     = t%year          , &
                        month    = t%month         , &
                        day      = t%day  -dt%day  , &
                        hour     = t%hour -dt%hour , &
                        zone     = t%zone          , &
                        min      = t%min  -dt%min  , &
                        sec      = t%sec  -dt%sec  , &
                        mili     = t%mili -dt%mili         )
    ! ~ normalize from negative or too large values into regular values:
    call Normalize( t_min_dt )

  end function t_min_dt


  ! *


  type(TIncrDate) function dt_min_dt( dt1, dt2 )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)       ::  dt1
    type(TIncrDate), intent(in)       ::  dt2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_min_dt'

    ! --- local ---------------------------------

    integer             ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( dt1, status )
    IF_NOTOK_STOP
    call Check( dt2, status )
    IF_NOTOK_STOP

    ! fill result:
    dt_min_dt = IncrDate( day  = dt1%day  - dt2%day , &
                          hour = dt1%hour - dt2%hour, &
                          min  = dt1%min  - dt2%min , &
                          sec  = dt1%sec  - dt2%sec , &
                          mili = dt1%mili - dt2%mili    )

  end function dt_min_dt


  ! *** date = t * r ************************************************


  ! multiply time with a real factor;
  ! use round for fractions

  type(TIncrDate) function dt_times_r( dt, r )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
    real, intent(in)              ::  r

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_times_r'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! multiply each of the parts with r, round
    dt_times_r = IncrDate( day  = nint( dt%day  * r ), &
                           hour = nint( dt%hour * r ), &
                           min  = nint( dt%min  * r ), &
                           sec  = nint( dt%sec  * r ), &
                           mili = nint( dt%mili * r )      )

  end function dt_times_r


  ! *


  type(TIncrDate) function r_times_dt( r, dt )

    ! --- in/out --------------------------------

    real, intent(in)              ::  r
    type(TIncrDate), intent(in)   ::  dt

    ! --- begin ---------------------------------

    r_times_dt = dt * r

  end function r_times_dt


  ! *


  type(TIncrDate) function dt_times_i( dt, i )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
    integer, intent(in)           ::  i

    ! --- begin ---------------------------------

    dt_times_i = dt * (i*1.0)

  end function dt_times_i


  ! *


  type(TIncrDate) function i_times_dt( i, dt )

    ! --- in/out --------------------------------

    integer, intent(in)           ::  i
    type(TIncrDate), intent(in)   ::  dt

    ! --- begin ---------------------------------

    i_times_dt = dt * i

  end function i_times_dt


  ! *** dt = dt / r ************************************************


  type(TIncrDate) function dt_div_r( dt, r )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
    real, intent(in)              ::  r

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_div_r'

    ! --- local ---------------------------------

    integer    ::  status
    real       ::  rat
    integer    ::  intg
    real       ::  frac

    ! --- begin ---------------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! days:
    rat = dt%day / r
    intg = floor( rat )
    frac = rat - intg
    dt_div_r = IncrDate( day=intg )

    ! hours:
    rat = dt%hour / r + frac*24
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, hour=intg )

    ! mins:
    rat = dt%min / r + frac*60
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, min=intg )

    ! seconds:
    rat = dt%sec / r + frac*60
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, sec=intg )

    ! miliseconds:
    rat = dt%mili / r + frac*1000
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, mili=intg )

  end function dt_div_r

  ! *

  type(TIncrDate) function dt_div_i( dt, i )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
    integer, intent(in)           ::  i

    ! --- begin ---------------------------------

    dt_div_i = dt / (i*1.0)

  end function dt_div_i


  ! ************************************************
  ! ***
  ! *** logical operators
  ! ***
  ! ************************************************


  logical function date_IsAnyDate( t )

    ! --- in/out -------------------------------

    type(TDate), intent(in)        ::  t

    ! --- begin --------------------------------

    date_IsAnyDate = t%calender == 'any'

  end function date_IsAnyDate

  subroutine date_set_date(t_out, t_in)
    ! t1 = t2
    ! --- in/out --------------------------------

    type(TDate), intent(in)     :: t_in
    type(TDate), intent(out)    :: t_out

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_set_date'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( t_in, status )
    IF_NOTOK_STOP

    t_out%year  = t_in%year
    t_out%month = t_in%month
    t_out%day   = t_in%day
    t_out%hour  = t_in%hour
    t_out%min   = t_in%min
    t_out%sec   = t_in%sec
    t_out%mili  = t_in%mili
    t_out%zone  = t_in%zone

    t_out%status = 0
    t_out%calender = t_in%calender

  end subroutine date_set_date

  subroutine incrdate_set_incrdate(t_out, t_in)
    ! dt1 = dt2
    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)     :: t_in
    type(TIncrDate), intent(out)    :: t_out

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_set_incrdate'

    ! --- begin ---------------------------------

    t_out%day   = t_in%day
    t_out%hour  = t_in%hour
    t_out%min   = t_in%min
    t_out%sec   = t_in%sec
    t_out%mili  = t_in%mili

    call Normalize(t_out)

    t_out%status = 0

  end subroutine incrdate_set_incrdate

  logical function isSameDay(t1, t2)

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1, t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/isSameDay'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    ! compare values
    isSameDay = &
        ( t1%year  == t2%year  ) .and. &
        ( t1%month == t2%month ) .and. &
        ( t1%day   == t2%day   )

  end function isSameDay

  logical function isSameMonth(t1, t2)

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1, t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/isSameMonth'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    ! compare values
    isSameMonth = &
        ( t1%year  == t2%year  ) .and. &
        ( t1%month == t2%month )

  end function isSameMonth

  ! ***  date1 == date2

  logical function date_eq_date( t1, t2 )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_eq_date'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( t1, status )
    if (status /=0) then
        write(gol, '(a, ": Error checking LHS date")') rname ; call goErr
        IF_NOTOK_STOP
    end if

    call Check( t2, status )
    if (status /=0) then
        write(gol, '(a, ": Error checking RHS date")') rname ; call goErr
        IF_NOTOK_STOP
    end if

    ! any date ? always equal
    if ( (t1%calender == 'any') .or. (t2%calender == 'any') ) then
      date_eq_date = .true.
      return
    end if

    ! compare values
    date_eq_date = &
        ( t1%year  == t2%year  ) .and. &
        ( t1%month == t2%month ) .and. &
        ( t1%day   == t2%day   ) .and. &
        ( t1%zone  == t2%zone  ) .and. &
        ( t1%hour  == t2%hour  ) .and. &
        ( t1%min   == t2%min   ) .and. &
        ( t1%sec   == t2%sec   ) .and. &
        ( t1%mili  == t2%mili  )

  end function date_eq_date


  ! ***  date1 /= date2

  logical function date_ne_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_ne_date'

    ! --- begin ---------------------------------

    date_ne_date = .not. ( t1 == t2 )

  end function date_ne_date


  ! ***  date1 > date2


  logical function date_gt_date( t1, t2 )

    use go_print, only : gol, goErr

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_gt_date'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( t1, status )
    IF_NOTOK_STOP
    call Check( t2, status )
    IF_NOTOK_STOP

    ! any date ? always true
    if ( (t1%calender == 'any') .or. (t2%calender == 'any') ) then
      date_gt_date = .true.
      return
    end if

    if ( t1%year > t2%year ) then
      date_gt_date = .true.
      return
    else if ( t1%year < t2%year ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%month > t2%month ) then
      date_gt_date = .true.
      return
    else if ( t1%month < t2%month ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%day > t2%day ) then
      date_gt_date = .true.
      return
    else if ( t1%day < t2%day ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%hour > t2%hour ) then
      date_gt_date = .true.
      return
    else if ( t1%hour < t2%hour ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%min > t2%min ) then
      date_gt_date = .true.
      return
    else if ( t1%min < t2%min ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%sec > t2%sec ) then
      date_gt_date = .true.
      return
    else if ( t1%sec < t2%sec ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%mili > t2%mili ) then
      date_gt_date = .true.
      return
    else if ( t1%mili < t2%mili ) then
      date_gt_date = .false.
      return
    end if

    ! all fields are equal ...
    date_gt_date = .false.

  end function date_gt_date


  ! ***  date1 < date2


  logical function date_lt_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_lt_date = (.not.( ( t1 == t2 ) .or. ( t1 > t2 ) )) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_lt_date


  ! ***  date1 >= date2 ************************


  logical function date_ge_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_ge_date = ( t1 == t2 ) .or. ( t1 > t2 ) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_ge_date


  ! ***  date1 <= date2 ************************


  logical function date_le_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_le_date = (.not. ( t1 > t2 )) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_le_date


  ! ***********************************************
  ! ***
  ! *** min/max
  ! ***
  ! ***********************************************


  function date_min( t1, t2 )

    ! --- in/out ---------------------------------

    type(TDate)                   ::  date_min
    type(TDate), intent(in)       ::  t1, t2

    ! --- begin ----------------------------------

    if ( t1 <= t2 ) then
      date_min = t1
    else
      date_min = t2
    end if

  end function date_min


  ! ***


  function date_max( t1, t2 )

    ! --- in/out ---------------------------------

    type(TDate)                   ::  date_max
    type(TDate), intent(in)       ::  t1, t2

    ! --- begin ----------------------------------

    if ( t1 >= t2 ) then
      date_max = t1
    else
      date_max = t2
    end if

  end function date_max


  ! ***********************************************
  ! ***
  ! *** totals
  ! ***
  ! ***********************************************


  real function date_rTotal( t, unit )

    use go_print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TDate), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_rTotal'

    ! --- local -----------------------------------

    integer            ::  status
    real               ::  nday
    integer            ::  iyear

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOTOK_STOP

    ! not all arguments are possible ...
    select case ( t%calender )
      case ( 'wall', 'greg', '366', '365' )
        select case ( unit )
          case ( 'year' )
            if ( any( (/t%month,t%day,t%hour,t%min,t%sec,t%mili/) /= 0 ) ) then
              write (gol,'("do not know how to count total:")'); call goErr
              write (gol,'("  unit : ",a)') unit; call goErr
              call wrtgol( '  t    : ', t ); call goErr
              TRACEBACK; stop
            end if
          case ( 'month' )
            if ( any( (/t%day,t%hour,t%min,t%sec,t%mili/) /= 0 ) ) then
              write (gol,'("do not know how to count total:")'); call goErr
              write (gol,'("  unit : ",a)') unit; call goErr
              call wrtgol( '  t    : ', t ); call goErr
              TRACEBACK; stop
            end if
        end select
      case ( 'incr' )
        select case ( unit )
          case ( 'year', 'month' )
            write (gol,'("do not know how to count total in incremental date:")') unit; call goErr
            write (gol,'("  unit : ",a)') unit; call goErr
            call wrtgol( '  t    : ', t ); call goErr
            TRACEBACK; stop
        end select
    end select

    ! precount total number of days for some of the units:
    select case ( unit )
      case ( 'day', 'hour', 'min', 'sec', 'mili' )
        nday = 0.0
        do iyear = 1, t%year-1
          nday  = nday  + calc_days_in_year(t%calender,iyear)
        end do
        nday  = nday + DayNumber( t ) - 1
    end select

    ! count time units:
    select case ( unit )
      case ( 'year' )
        ! set 'nday' to a reference length of the year;
        ! if this length is not constant during the years, the
        ! values of t%month etc have been checked to be zero:
        nday = days_in_year(t) * 1.0
        ! count fractional years:
        date_rTotal = t%year                                           + &
                      t%month / 12.0                                   + &
                      t%day   / nday                                   + &
                      t%hour  / nday / 24.0                            + &
                      t%min   / nday / 24.0 / 60.0                     + &
                      t%sec   / nday / 24.0 / 60.0 / 60.0              + &
                      t%mili  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'month' )
        ! set 'nday' to a reference length of the month;
        ! if this length is not constant during the years, the
        ! values of t%day etc been checked to be zero:
        nday = days_in_month(t) * 1.0
        ! count fractional months:
        date_rTotal = t%year  * 12.0                                   + &
                      t%month                                          + &
                      t%day   / nday                                   + &
                      t%hour  / nday / 24.0                            + &
                      t%min   / nday / 24.0 / 60.0                     + &
                      t%sec   / nday / 24.0 / 60.0 / 60.0              + &
                      t%mili  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        date_rTotal = nday                                             + &
                      t%hour         / 24.0                            + &
                      t%min          / 24.0 / 60.0                     + &
                      t%sec          / 24.0 / 60.0 / 60.0              + &
                      t%mili         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        date_rTotal = nday           * 24.0                            + &
                      t%hour                                           + &
                      t%min                 / 60.0                     + &
                      t%sec                 / 60.0 / 60.0              + &
                      t%mili                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        date_rTotal = nday           * 24.0 * 60.0                     + &
                      t%hour                * 60.0                     + &
                      t%min                                            + &
                      t%sec                        / 60.0              + &
                      t%mili                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        date_rTotal = nday           * 24.0 * 60.0 * 60.0              + &
                      t%hour                * 60.0 * 60.0              + &
                      t%min                        * 60.0              + &
                      t%sec                                            + &
                      t%mili                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        date_rTotal = nday           * 24.0 * 60.0 * 6.00 * 1000.0     + &
                      t%hour                * 60.0 * 60.0 * 1000.0     + &
                      t%min                        * 60.0 * 1000.0     + &
                      t%sec                               * 1000.0     + &
                      t%mili
      case default
        write (gol,'("do not know how to count time in unit : ",a)') trim(unit); call goErr
        TRACEBACK; stop
    end select

  end function date_rTotal


  ! ***


  real function incr_rTotal( dt, unit )

    use go_print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TIncrDate), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incr_rTotal'

    ! --- local -----------------------------------

    integer            ::  status

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! count time units:
    select case ( unit )
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        incr_rTotal = dt%day                                            + &
                      dt%hour         / 24.0                            + &
                      dt%min          / 24.0 / 60.0                     + &
                      dt%sec          / 24.0 / 60.0 / 60.0              + &
                      dt%mili         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        incr_rTotal = dt%day          * 24.0                            + &
                      dt%hour                                           + &
                      dt%min                 / 60.0                     + &
                      dt%sec                 / 60.0 / 60.0              + &
                      dt%mili                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        incr_rTotal = dt%day          * 24.0 * 60.0                     + &
                      dt%hour                * 60.0                     + &
                      dt%min                                            + &
                      dt%sec                        / 60.0              + &
                      dt%mili                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        incr_rTotal = dt%day          * 24.0 * 60.0 * 60.0              + &
                      dt%hour                * 60.0 * 60.0              + &
                      dt%min                        * 60.0              + &
                      dt%sec                                            + &
                      dt%mili                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        incr_rTotal = dt%day          * 24.0 * 60.0 * 6.00 * 1000.0     + &
                      dt%hour                * 60.0 * 60.0 * 1000.0     + &
                      dt%min                        * 60.0 * 1000.0     + &
                      dt%sec                               * 1000.0     + &
                      dt%mili
      case default
        write (gol,'("do not know how to count time in unit : ",a)') trim(unit); call goErr
        TRACEBACK; stop
    end select

  end function incr_rTotal


  ! ***


  integer function date_iTotal( t, unit )

    use go_print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TDate), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_iTotal'

    ! --- local -----------------------------

    integer            ::  status
    real               ::  rtot
    integer            ::  itot

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOTOK_STOP

    ! determine total some as a real value:
    rtot = rTotal( t, unit )

    ! round to integer value:
    itot = nint(rtot)

    ! result should be pure integer ....
    if ( itot*1.0 == rtot ) then
      date_iTotal = itot
    else
      write (gol,'("date does not contain integer total:")'); call goErr
      write (gol,'("  unit : ",a)') trim(unit); call goErr
      call wrtgol( '  t    : ', t ); call goErr
      TRACEBACK; stop
    end if

  end function date_iTotal


  ! ***


  integer function incrdate_iTotal( dt, unit )

    use go_print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TIncrDate), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_iTotal'

    ! --- local -----------------------------

    integer            ::  status
    real               ::  rtot
    integer            ::  itot

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! determine total some as a real value:
    rtot = rTotal( dt, unit )

    ! round to integer value:
    itot = nint(rtot)

    ! result should be pure integer ....
    if ( itot*1.0 == rtot ) then
      incrdate_iTotal = itot
    else
      write (gol,'("date does not contain integer total:")'); call goErr
      write (gol,'("  unit : ",a)') trim(unit); call goErr
      call wrtgol( '  dt   : ', dt ); call goErr
      TRACEBACK; stop
    end if

  end function incrdate_iTotal


  ! ***********************************************
  ! ***
  ! *** interpolation
  ! ***
  ! ***********************************************

  !
  ! Return coeff such that
  !   t = alfa1 * t1 + alfa2 * t2
  !

  subroutine date_InterpolFractions( t, t1, t2, alfa1, alfa2, status )

    use go_print, only : gol, goErr

    ! --- in/out -----------------------------

    type(TDate), intent(in)    ::  t
    type(TDate), intent(in)    ::  t1
    type(TDate), intent(in)    ::  t2
    real, intent(out)          ::  alfa1
    real, intent(out)          ::  alfa2
    integer, intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_InterpolFractions'

    ! --- local ------------------------------

    real    ::  ds, ds1

    ! --- begin ------------------------------

    ! check ...
    if ( ((t1 < t2) .and. ((t < t1) .or. (t2 < t))) .or. &
         ((t2 < t1) .and. ((t < t2) .or. (t1 < t)))      ) then
      write (gol,'("t not in [t1,t2] :")'); call goErr
      call wrtgol( '  t  = ', t  ); call goErr
      call wrtgol( '  t1 = ', t1 ); call goErr
      call wrtgol( '  t2 = ', t2 ); call goErr
      TRACEBACK; status=1; return
    end if

    ! compute differences in seconds:
    ds  = rTotal( t2 - t1, 'sec' )
    ds1 = rTotal( t  - t1, 'sec' )

    ! (almost) zero length interval?
    if ( abs(ds) < tiny(ds) ) then
      ! zero length, interpolate to middle:
      alfa2 = 0.5
    else
      ! fraction for second value:
      alfa2 = ds1 / ds
    end if
    ! fraction for first value:
    alfa1 = 1.0 - alfa2

    ! ok
    status = 0

  end subroutine date_InterpolFractions


  ! ***


!  subroutine date_IntervalFractions( t0, dt, t1, t2, n, ii, ff, status )
!
!    ! --- in/out ---------------------------------
!
!    type(TDate), intent(in)             ::  t0
!    type(TIncrDate), intent(in)         ::  dt
!    type(TDate), intent(in)             ::  t1, t2
!    integer, intent(out)                ::  n
!    integer, pointer, intent(inout)     ::  ii(:)
!    real, pointer, intent(inout)        ::  ff(:)
!    integer, intent(out)                ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/date_IntervalFractions'
!
!    ! --- local ------------------------------
!
!    real          ::  dt_sec
!    integer       ::  i1, i2
!    type(TDate)   ::  st1, st2
!    integer       ::  k
!
!    ! --- begin ----------------------------------
!
!    ! check ...
!    if ( t1 > t2 ) then
!      call wrtgol( 'no valid interval : ', t1, ' - ', t2 ); call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! check ...
!    if ( (t1 < t0) .or. (t2 < t0) ) then
!      call wrtgol( 'interval (partly) before start time : ' ); call goErr
!      call wrtgol( '  start time : ', t0 ); call goErr
!      call wrtgol( '  interval   : ', t1, ' - ', t2 ); call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! time step in seconds:
!    dt_sec = rTotal( dt, 'sec' )
!
!    ! intervals with t1 and t2:
!    i1 = ceiling( rTotal( t1-t0, 'sec' ) / dt_sec )
!    i2 = ceiling( rTotal( t2-t0, 'sec' ) / dt_sec )
!
!    ! total number:
!    n = i2 - i1 + 1
!
!    ! set arrays:
!    if ( associated(ii) ) deallocate( ii )
!    allocate( ii(n) )
!    if ( associated(ff) ) deallocate( ff )
!    allocate( ff(n) )
!
!    ! loop over intervals:
!    k = 0
!    do i = i1, i2
!      ! joint sub interval:
!      st1 = max( t1, t0+(i-1)*dt )
!      st2 = min( t2, t0+    i*dt )
!      ! store index and fraction:
!      k = k + 1
!      ii(k) = i
!      ff(k) = rTotal(st2-st1,'sec') / dt_sec
!    end do
!
!    ! ok
!    status = 0
!
!  end subroutine date_IntervalFractions


  ! ***********************************************
  ! ***
  ! *** interpolation
  ! ***
  ! ***********************************************

  !
  ! Return interval [tt(1),tt(2)] around t that matches with time resolution;
  ! resolution specified by a step and unit:
  !    3, 'hour'   # 00:00, 03:00, 06:00, ...
  ! Interval boundaries are the same if t is exactly on a step.
  !

  subroutine Get_Surrounding_Interval( t, step, units, tt, status )

    ! --- in/out -----------------------------

    type(TDate), intent(in)       ::  t
    integer, intent(in)           ::  step
    character(len=*), intent(in)  ::  units
    type(TDate), intent(out)      ::  tt(2)
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Get_Surrounding_Interval'

    ! --- local ------------------------------

    real          ::  r
    type(TDate)   ::  t0

    ! --- begin ------------------------------

    ! which units ?
    select case ( units )

      ! x hourly
      case ( 'hour' )
        ! begin of day:
        t0 = Get_Begin_Of( t, 'day' )
        ! hour fraction:
        r = rTotal( t - t0, 'hour' )
        ! start of interval:
        tt(1) = t0 + IncrDate( hour=int(floor(r/step))*step )
        ! end:
        if ( t == tt(1) ) then
          tt(2) = tt(1)
        else
          tt(2) = tt(1) + IncrDate(hour=step)
        end if

      ! unknown ...
      case default
        write (gol,'("unsupported units `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine Get_Surrounding_Interval


  ! ***********************************************
  ! ***
  ! *** print
  ! ***
  ! ***********************************************

  character(len=PRINT_LEN) function date_Pretty( t )

    ! --- in/out -------------------------

    type(TDate), intent(in)       ::  t

    ! --- const --------------------------

    character(len=3), parameter  ::  month_name(12) = &
                                      (/ 'jan','feb','mar','apr','may','jun', &
                                         'jul','aug','sep','oct','nov','dec' /)

    ! --- local --------------------------

    integer                   ::  zone_abs, zone_hour, zone_min
    character(len=1)          ::  zone_sign
    character(len=PRINT_LEN)  ::  s

    ! --- begin --------------------------

    select case ( t%calender )
      case ( 'wall' )
        if ( t%zone < 0 ) then
          zone_sign = '-'
        else
          zone_sign = '+'
        end if
        zone_abs  = abs(t%zone)
        zone_hour = floor(zone_abs/60.0)
        zone_min  = zone_abs - zone_hour*60
        write (s,'(i2,":",i2.2,":",i2.2,":",i3.3,       &
                  & " ",i2.2," ",a3," ",i4.4,           &
                  & " (GMT",a1,i2.2,":",i2.2,")")')     &
          t%hour, t%min, t%sec, t%mili, &
          t%day, month_name(t%month), t%year, &
          zone_sign, zone_hour, zone_min
      case ( 'greg', '366', '365', '360', 'any' )
        ! only print seconds and milliseconds if they are non-zero
        if (t%mili == 0) then
            if (t%sec == 0) then
                write (s,'(i4.4,"/",i2.2,"/",i2.2," ",i2,":",i2.2)') &
                    t%year, t%month, t%day, t%hour, t%min
            else
                write (s,'(i4.4,"/",i2.2,"/",i2.2," ",i2,":",i2.2,":",i2.2)') &
                    t%year, t%month, t%day, t%hour, t%min, t%sec
            end if
        else
            write (s,'(i4.4,"/",i2.2,"/",i2.2," ",i2,":",i2.2,":",i2.2,".",i3.3)') &
                t%year, t%month, t%day, t%hour, t%min, t%sec, t%mili
        end if
      case default
        s = 'no-calender'
    end select

    date_Pretty = s

  end function date_Pretty


  ! *


  character(len=PRINT_LEN) function incrdate_Pretty( dt )

    ! --- in/out -------------------------

    type(TIncrDate), intent(in)       ::  dt

    ! --- local --------------------------

    character(len=PRINT_LEN)       ::  s

    ! --- begin --------------------------

    if (dt%mili == 0) then
        if (dt%sec == 0) then
            write (s,'(i5," days ",i2,":",i2.2)') dt%day, dt%hour, dt%min
        else
            write (s,'(i5," days ",i2,":",i2.2,":",i2.2)') dt%day, dt%hour, dt%min, dt%sec
        end if
    else
        write (s,'(i5," days ",i2,":",i2.2,":",i2.2,".",i3.3)') dt%day, dt%hour, dt%min, dt%sec, dt%mili
    end if

    incrdate_Pretty = s

  end function incrdate_Pretty


  ! *


  subroutine wrtgol_t( msg, t )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t

    ! --- local ---------------------------------

    character(len=PRINT_LEN)     ::  s

    ! --- begin -----------------------------------

    s = date_Pretty( t )
    write (gol,'(a,a)') msg, trim(s)

  end subroutine wrtgol_t


  ! *


  subroutine wrtgol_dt( msg, dt )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TIncrDate), intent(in)     ::  dt

    ! --- local ---------------------------------

    character(len=PRINT_LEN)     ::  s

    ! --- begin -----------------------------------

    s = incrdate_Pretty( dt )
    write (gol,'(a,a)') msg, trim(s)

  end subroutine wrtgol_dt


  ! *


  subroutine wrtgol_t1_t2( msg, t, msg2, t2 )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(TDate), intent(in)         ::  t2

    ! --- local ---------------------------------

    character(len=PRINT_LEN)     ::  s
    character(len=PRINT_LEN)     ::  s2

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = date_Pretty( t2 )
    write (gol,'(a,a,a,a)') msg, trim(s), msg2, trim(s2)

  end subroutine wrtgol_t1_t2


  ! *


  subroutine wrtgol_t1_t2_t3( msg, t, msg2, t2, msg3, t3 )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(TDate), intent(in)         ::  t2
    character(len=*), intent(in)    ::  msg3
    type(TDate), intent(in)         ::  t3

    ! --- local ---------------------------------

    character(len=PRINT_LEN)     ::  s
    character(len=PRINT_LEN)     ::  s2
    character(len=PRINT_LEN)     ::  s3

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = date_Pretty( t2 )
    s3 = date_Pretty( t3 )
    write (gol,'(a,a,a,a,a,a)') msg, trim(s), msg2, trim(s2), msg3, trim(s3)

  end subroutine wrtgol_t1_t2_t3


end module GO_Date


