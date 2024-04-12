module datetime_for

implicit none

private

public :: mktime, format_datetime, time_tm, time_dt, utctime, isleap, nextDay, nextMonth
public :: assignment(=), operator(+), operator(-)
public :: kind_mktime
public :: SEC_PER_DAY, SEC_PER_HR, SEC_PER_MIN

type time_tm
    integer :: tm_year, tm_mon, tm_mday, tm_hour, tm_min, tm_sec, tm_wday, tm_yday
end type time_tm

type time_dt
    integer :: days, seconds
end type time_dt

integer, parameter      :: YEAR_ZERO = 1900, SEC_PER_DAY = 86400, SEC_PER_HR = 3600, SEC_PER_MIN = 60
integer, parameter      :: MIN_PER_HR = 60, HR_PER_DAY = 24, MON_PER_YR = 12
integer                 :: mon_lengths(2,12) = reshape((/ 31, 31, 28, 29, 31, 31, 30, 30, 31, 31, &
                           30, 30, 31, 31, 31, 31, 30, 30, 31, 31, 30, 30, 31, 31 /), (/ 2, 12 /))
integer                 :: year_lengths(2) = (/ 365, 366 /)
integer, parameter      :: kind_mktime = selected_int_kind(10)

interface assignment (=)
    module procedure array_to_tm
    module procedure tm_to_array
    module procedure array_to_tm_short
    module procedure tm_to_array_short
    module procedure array_to_dt
    module procedure dt_to_array
    module procedure tm_to_tm
end interface

interface operator (+)
    module procedure add_times_t_dt
    module procedure add_times_dt_t
end interface

interface operator(-)
    module procedure subtract_times
    module procedure subtract_dt
end interface

contains

function nextMonth(time)
    ! Returns the starting point of the next month. For example, nextMonth(12:30, May 23, 2004) gives (0:0, Jun 1, 2004),
    ! while nextMonth(0:0, Jun 1, 2004) gives (0:0, Jul 1, 2004).
    type(time_tm), intent(in)   :: time
    type(time_tm)               :: nextMonth

    nextMonth%tm_mon = time%tm_mon + 1

    if (nextMonth%tm_mon > 12) then
        nextMonth%tm_mon = 1
        nextMonth%tm_year = time%tm_year + 1
    else
        nextMonth%tm_year = time%tm_year
    end if

    nextMonth%tm_mday = 1
    nextMonth%tm_hour = 0
    nextMonth%tm_min = 0
    nextMonth%tm_sec = 0

end function nextMonth

function nextDay(time)
    ! Identical to nextMonth, but returns the start of the next day
    type(time_tm), intent(in)   :: time
    type(time_tm)               :: nextDay
    integer                     :: li

    li = isleap(time%tm_year)

    nextDay%tm_mday = time%tm_mday + 1

    if (nextDay%tm_mday > mon_lengths(li,time%tm_mon)) then
        nextDay%tm_mday = 1
        nextDay%tm_mon = time%tm_mon + 1
    else
        nextDay%tm_mon = time%tm_mon
    end if

    if (nextDay%tm_mon > 12) then
        nextDay%tm_mon = 1
        nextDay%tm_year = time%tm_year + 1
    else
        nextDay%tm_year = time%tm_year
    end if

    nextDay%tm_hour = 0
    nextDay%tm_min = 0
    nextDay%tm_sec = 0

end function nextDay

function subtract_times(time1, time2) result (dt)

    type(time_tm), intent(in) :: time1, time2
    type(time_dt)             :: dt
    integer(kind=8)           :: t1, t2
    type(time_tm)             :: s_time, e_time

    ! assume that time1 and time2 are sanitized
    s_time = time2
    e_time = time1
    t1 = mktime(e_time)
    t2 = mktime(s_time)
    dt%days = (t1-t2)/SEC_PER_DAY
    dt%seconds = mod(t1-t2, int(SEC_PER_DAY, kind(t1)))

end function subtract_times

function subtract_dt(time_e, dt) result (time_s)

    type(time_tm), intent(in) :: time_e
    type(time_dt), intent(in) :: dt
    type(time_tm)             :: time_s, scratch_time
    integer(kind=8)           :: start_time, time_delta

    ! assume that time_e is already sanitized
    scratch_time = time_e
    start_time = mktime(scratch_time)
    time_delta = SEC_PER_DAY * dt%days + dt%seconds
    if (time_delta > start_time) then
        time_s = (/ -1, -1, -1, -1, -1, -1, -1, -1 /)
        return
    end if
    start_time = start_time - time_delta
    time_s = utctime(start_time)

end function subtract_dt

function add_times_dt_t(timedelta, timep) result (fin_time)

    type(time_tm), intent(in) :: timep
    type(time_dt), intent(in) :: timedelta
    type(time_tm)             :: fin_time

    ! assume that timep is already sanitized
    fin_time = timep
    fin_time%tm_sec = fin_time%tm_sec + timedelta%seconds
    call sanitize(fin_time)
    fin_time%tm_mday = fin_time%tm_mday + timedelta%days
    call sanitize(fin_time)

end function add_times_dt_t

function add_times_t_dt(timep, timedelta) result (fin_time)

    type(time_tm), intent(in) :: timep
    type(time_dt), intent(in) :: timedelta
    type(time_tm)             :: fin_time

    ! assume that timep is already sanitized
    fin_time = timep
    fin_time%tm_sec = fin_time%tm_sec + timedelta%seconds
    call sanitize(fin_time)
    fin_time%tm_mday = fin_time%tm_mday + timedelta%days
    call sanitize(fin_time)

end function add_times_t_dt

subroutine tm_to_tm(ttarget,tsource)

    type(time_tm), intent(in)   :: tsource
    type(time_tm), intent(out)  :: ttarget

    ttarget%tm_year = tsource%tm_year
    ttarget%tm_mon = tsource%tm_mon
    ttarget%tm_mday = tsource%tm_mday
    ttarget%tm_hour = tsource%tm_hour
    ttarget%tm_min = tsource%tm_min
    ttarget%tm_sec = tsource%tm_sec
    ttarget%tm_wday = tsource%tm_wday
    ttarget%tm_yday = tsource%tm_yday

end subroutine tm_to_tm

subroutine array_to_dt(dt,cdate)

    integer, intent(in)        :: cdate(2)
    type(time_dt), intent(out) :: dt

    dt%days = cdate(1)
    dt%seconds = cdate(2)

end subroutine array_to_dt

subroutine dt_to_array(cdate,dt)

    type(time_dt), intent(in) :: dt
    integer, intent(out)      :: cdate(2)

    cdate(1) = dt%days
    cdate(2) = dt%seconds

end subroutine dt_to_array

subroutine array_to_tm(t,cdate)

    integer, intent(in)        :: cdate(:)
    type(time_tm), intent(out) :: t

    if (size(cdate) .lt. 3) then
        t%tm_year = -1
        t%tm_mon = -1
        t%tm_mday = -1
        t%tm_hour = -1
        t%tm_min = -1
        t%tm_sec = -1
        t%tm_wday = -1
        t%tm_yday = -1
        return
    end if
    t%tm_hour = 0
    t%tm_min = 0
    t%tm_sec = 0
    t%tm_year = cdate(1)
    t%tm_mon  = cdate(2)
    t%tm_mday = cdate(3)
    if (size(cdate) .ge. 4) t%tm_hour = cdate(4)
    if (size(cdate) .ge. 5) t%tm_min  = cdate(5)
    if (size(cdate) .ge. 6) t%tm_sec  = cdate(6)
    call sanitize(t)

end subroutine array_to_tm

subroutine array_to_tm_short(t,cdate)

    integer(2), intent(in)     :: cdate(:)
    type(time_tm), intent(out) :: t

    if (size(cdate) .lt. 3) then
        t%tm_year = -1
        t%tm_mon = -1
        t%tm_mday = -1
        t%tm_hour = -1
        t%tm_min = -1
        t%tm_sec = -1
        t%tm_wday = -1
        t%tm_yday = -1
        return
    end if
    t%tm_hour = 0
    t%tm_min = 0
    t%tm_sec = 0
    t%tm_year = cdate(1)
    t%tm_mon  = cdate(2)
    t%tm_mday = cdate(3)
    if (size(cdate) .ge. 4) t%tm_hour = cdate(4)
    if (size(cdate) .ge. 5) t%tm_min  = cdate(5)
    if (size(cdate) .ge. 6) t%tm_sec  = cdate(6)
    ! set some useful arrays
    mon_lengths(1,:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    mon_lengths(2,:) = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    call sanitize(t)

end subroutine array_to_tm_short

subroutine tm_to_array(cdate,t)

    type(time_tm), intent(in) :: t
    integer, intent(out)      :: cdate(6)

    cdate(1) = t%tm_year
    cdate(2) = t%tm_mon
    cdate(3) = t%tm_mday
    cdate(4) = t%tm_hour
    cdate(5) = t%tm_min
    cdate(6) = t%tm_sec

end subroutine tm_to_array

subroutine tm_to_array_short(cdate,t)

    type(time_tm), intent(in) :: t
    integer(2), intent(out)   :: cdate(6)

    cdate(1) = t%tm_year
    cdate(2) = t%tm_mon
    cdate(3) = t%tm_mday
    cdate(4) = t%tm_hour
    cdate(5) = t%tm_min
    cdate(6) = t%tm_sec

end subroutine tm_to_array_short

function isleap(y)
    integer, intent(in) :: y
    integer             :: isleap
    logical             :: leap

    isleap = 1
    leap = (mod(y,4) == 0 .and. mod(y,100) /= 0) .or. mod(y,400) == 0
    if (leap) isleap = 2
end function isleap

function format_datetime(timep) result (date_str)

    type(time_tm), intent(in) :: timep
    character(len=80)         :: date_str

    write(date_str, '(i0.2,a,i0.2,a,i0.4,a,i0.2,a,i0.2,a,i0.2)') &
        timep%tm_mday,'/',timep%tm_mon,'/',timep%tm_year,' ',&
        timep%tm_hour,':',timep%tm_min,':',timep%tm_sec
    date_str = trim(adjustl(date_str))

end function format_datetime

subroutine sanitize(timep)

    ! changes timep (if needed) so that 1 <= months <= 12, 0 <= hours <= 23, etc.
    type(time_tm), intent(inout) :: timep

    timep%tm_min = timep%tm_min + timep%tm_sec/SEC_PER_MIN
    timep%tm_sec = mod(timep%tm_sec, SEC_PER_MIN)
    if (timep%tm_sec < 0) then
        timep%tm_sec = timep%tm_sec + SEC_PER_MIN
        timep%tm_min = timep%tm_min - 1
    end if
    timep%tm_hour = timep%tm_hour + timep%tm_min/MIN_PER_HR
    timep%tm_min = mod(timep%tm_min, MIN_PER_HR)
    if (timep%tm_min < 0) then
        timep%tm_min = timep%tm_min + MIN_PER_HR
        timep%tm_hour = timep%tm_hour - 1
    end if
    timep%tm_mday = timep%tm_mday + timep%tm_hour/HR_PER_DAY
    timep%tm_hour = mod(timep%tm_hour, HR_PER_DAY)
    if (timep%tm_hour < 0) then
        timep%tm_hour = timep%tm_hour + HR_PER_DAY
        timep%tm_mday = timep%tm_mday - 1
    end if
    timep%tm_year = timep%tm_year + (timep%tm_mon-1)/MON_PER_YR
    timep%tm_mon = mod(timep%tm_mon, MON_PER_YR)
    if (timep%tm_mon == 0) timep%tm_mon = MON_PER_YR
    ! do not account for negative months now
    do while (timep%tm_mday > mon_lengths(isleap(timep%tm_year), timep%tm_mon))
        timep%tm_mday = timep%tm_mday - mon_lengths(isleap(timep%tm_year), timep%tm_mon)
        timep%tm_mon = timep%tm_mon + 1
        timep%tm_year = timep%tm_year + (timep%tm_mon-1)/MON_PER_YR
        timep%tm_mon = mod(timep%tm_mon, MON_PER_YR)
        if (timep%tm_mon == 0) timep%tm_mon = MON_PER_YR
    end do
    timep%tm_year = timep%tm_year + (timep%tm_mon-1)/MON_PER_YR

end subroutine sanitize

function utctime(t) result(timep)

    ! converts seconds since the start of YEAR_ZERO into timep
    type(time_tm)               :: timep
    integer(kind=8), intent(in) :: t
    integer(kind=8)             :: seconds
    integer                     :: year, day, month

    year = YEAR_ZERO
    seconds = t
    do while (seconds > SEC_PER_DAY*year_lengths(isleap(year)))
        seconds = seconds - SEC_PER_DAY*year_lengths(isleap(year))
        year = year + 1
    end do
    timep%tm_year = year
    month = 1
    do while (seconds > SEC_PER_DAY*mon_lengths(isleap(year),month))
        seconds = seconds - SEC_PER_DAY*mon_lengths(isleap(year),month)
        month = month + 1
    end do
    timep%tm_mon = month
    timep%tm_mday = (seconds/SEC_PER_DAY) + 1
    seconds = mod(seconds, int(SEC_PER_DAY, kind(seconds)))
    timep%tm_hour = seconds/SEC_PER_HR
    seconds = mod(seconds, int(SEC_PER_HR, kind(seconds)))
    timep%tm_min = seconds/SEC_PER_MIN
    timep%tm_sec = mod(seconds, int(SEC_PER_MIN, kind(seconds)))

end function utctime

function mktime(timep) result (seconds)

    ! converts timep into seconds since YEAR_ZERO
    type(time_tm), intent(inout) :: timep
    integer                      :: year, month
    integer(kind=8)              :: seconds, day

    call sanitize(timep)
    year = YEAR_ZERO
    day = 0
    do while (year < timep%tm_year)
        day = day + year_lengths(isleap(year))
        year = year + 1
    end do
    month = 1
    do while (month < timep%tm_mon)
        day = day + mon_lengths(isleap(year),month)
        month = month + 1
    end do
    day = day + timep%tm_mday - 1
    seconds = day * SEC_PER_DAY + timep%tm_hour * SEC_PER_HR + timep%tm_min * SEC_PER_MIN + timep%tm_sec

end function mktime

end module datetime_for
