!###############################################################################
!
! module contains all routines that deal with date/time calculation
!
! subroutine  chardate(idate6,cdate)
! subroutine  inctime
! subroutine  tau2date(itaux,idatex)
! subroutine  date2tau(idatex,itaux)
! subroutine  calc_sm( mlen, sec_day, sec_month, sec_year )
! subroutine  dayl(day,daylen,jdim,lat_start,dlat)
! subroutine  caldat(julian,mm,id,iyyy)
! integer function  julday(mm,id,iy)
! integer function get_day(mm,dd,mlen)
! subroutine  tstamp(kunit,itaux,msg)
! integer function get_num_days(tau_start, tau_end)
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

module datetime

  use GO, only : gol, goErr, goPr
  use GO, only : T_Time_Window

  implicit none


  ! --- in/out -------------------------------------

  private

  public :: inctime, tau2date, date2tau, dayl, idate2ddate
!  public :: calc_sm, get_day
  public :: julday, new_valid_timestep
  public :: tstamp, wrtgol_tstamp
  public :: get_num_days, isDayBoundary
  public :: solar_zenith_angle
  public :: system_clock_value_start
  public  ::  time_window, wall_clock_time

  interface date2tau
    module procedure date2tau_int
    module procedure date2tau_short
  end interface date2tau

  ! --- var --------------------------------------

  double precision            ::  system_clock_value_start
  type(T_Time_Window)         ::  time_window


contains

  function wall_clock_time()

    implicit none

    integer   :: clock_max, clock_reading
    real      :: clock_rate, wall_clock_time

    call system_clock(clock_reading, clock_rate, clock_max)
    wall_clock_time = real(clock_reading)/clock_rate

  end function wall_clock_time

  function solar_zenith_angle(int_date, lat, lon)
    !-----------------------------------------------------------------------
    ! Given the integer date as (year, month, day, hour, minute, second), and
    ! the latitude and longitude in degrees, calculate the solar zenith angle.
    ! https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
    ! Checked against https://www.esrl.noaa.gov/gmd/grad/antuv/SolarCalc.jsp
    !-----------------------------------------------------------------------
    use binas,  only        : pi

    integer, intent(in)     :: int_date(6)
    real, intent(in)        :: lat, lon
    real                    :: solar_zenith_angle

    real                    :: frac_year, eqtime, decl, time_offset, tst, ha, gamma, cos_phi

    ! calculate fractional year
    frac_year = idate2ddate(int_date)
    gamma = frac_year - floor(frac_year)
    gamma = 2.0 * pi * gamma ! convert to radians

    eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) - 0.014615 * cos(2*gamma) - 0.040849 * sin(2*gamma)) ! minutes
    decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - 0.006758 * cos(2*gamma) + 0.000907 * sin(2*gamma) - 0.002697 * cos(3*gamma) + 0.00148 * sin(3*gamma) ! radians
    time_offset = eqtime + 4.0*lon ! minute
    tst = int_date(4) * 60.0 + int_date(5) + int_date(6)/60.0 + time_offset ! minute
    ha = (pi/180.0) * (tst/4.0 - 180.0) ! radians
    cos_phi = sin(pi*lat/180.0) * sin(decl) + cos(pi*lat/180.0) * cos(decl) * cos(ha)
    solar_zenith_angle = 180.0 * acos(cos_phi)/pi

  end function solar_zenith_angle

  function get_num_days(tau_start, tau_end)
    !-----------------------------------------------------------------------
    !
    !   Given two integers, tau_start and tau_end, return the number of days
    !   that span that period, including partial days. For example, if
    !   tau_start -> 2009/3/12 18:30:00 and tau_end -> 2009/3/25 07:20:00,
    !   then the return value should be 14, i.e., the number of days from
    !   2009/3/12 to 2009/3/25, both inclusive.
    !
    !-----------------------------------------------------------------------

    integer, intent(in)  :: tau_start, tau_end
    integer :: get_num_days

    integer :: temp_date(6), another_temp_date(6)
    integer :: tau_start_copy, tau_end_copy
    logical :: adjoint

    get_num_days = 0
    adjoint = (tau_end < tau_start) ! switch to detect an adjoint run

    ! Is the starting time on a day boundary?
    call tau2date(tau_start, temp_date)
    if (all(temp_date(4:6) == (/0,0,0/))) then
        tau_start_copy = tau_start
    else
        get_num_days = get_num_days + 1
        another_temp_date = (/ temp_date(1), temp_date(2), temp_date(3), 0, 0, 0/)
        call date2tau(another_temp_date, tau_start_copy)
        if (.not. adjoint) tau_start_copy = tau_start_copy + 86400
    end if

    ! Is the ending time on a day boundary?
    call tau2date(tau_end, temp_date)
    if (all(temp_date(4:6) == (/0,0,0/))) then
        tau_end_copy = tau_end
    else
        get_num_days = get_num_days + 1
        another_temp_date = (/ temp_date(1), temp_date(2), temp_date(3), 0, 0, 0/)
        call date2tau(another_temp_date, tau_end_copy)
        if (adjoint) tau_end_copy = tau_end_copy + 86400
    end if

    ! Now count the whole number of days
    if (adjoint) then
        get_num_days = get_num_days + (tau_start_copy-tau_end_copy)/86400
    else
        get_num_days = get_num_days + (tau_end_copy-tau_start_copy)/86400
    end if

  end function get_num_days

  ! ==============================================================


  !-----------------------------------------------------------------------
  !
  !     put date/time in character string
  !
  !     on input:  idate6 contains date/time
  !     on output: cdata contains date/time string
  !
  !-----------------------------------------------------------------------

  subroutine chardate(idate6,cdate)

    integer,dimension(6),intent(in) :: idate6  ! "chardate" input date
    character(24)                   :: cdate   ! "chardate" output string
    !
    integer :: k
    integer :: iostat
    character(3),dimension(12),parameter :: mon= &
         (/'Jan','Feb','Mar','Apr','May','Jun', &
         'Jul','Aug','Sep','Oct','Nov','Dec'/)
    !
    ! put date/time in string
    !
    write(cdate,'(i4.4,"-",a3,"-",i2.2," ",i2.2,":",i2.2,":",i2.2,"   ")',iostat=iostat) &
         idate6(1),mon(idate6(2)),idate6(3:6)
    !if (iostat/=0) then
    !  write (*,*) 'ERROR - writing date values : ', idate6(1),mon(idate6(2)),idate6(3:6)
    !  stop 'ERROR in datetime/chardate'
    !end if
    !
  end subroutine chardate


  !-----------------------------------------------------------------------
  !**** inctime
  !
  !     purpose
  !     -------
  !     increment time and set newday/newmonth/newyr switches
  !
  !     influences
  !     ----------
  !     nstep,itau,idate,newyr,newmonth,newday,newsrun
  !
  !     externals
  !     ---------
  !     subroutines: tau2date
  !                  tstamp
  !-----------------------------------------------------------------------

  subroutine inctime( ndyn, status )

    !use dims, only : ndyn
    use dims, only : tref,idate,itau,nstep,revert,nregions
    use dims, only : cdebug,kdebug,newyr,newmonth,newday,newsrun, newhour

    ! in/out
    integer, intent(in)   ::  ndyn
    integer, intent(out)  ::  status

    ! local
    integer              :: ninc, k, itautmp
    integer,dimension(6) :: idtmp
    !
    ! add time step of ninc seconds

    ninc=ndyn/(2*tref(1))  !cmk   !region = 1 is master
    !
    nstep=nstep+1
    itau=itau+revert*ninc
!    do k = 1, 6
!       idtmp(k) = idate(k)
!    end do
    call tau2date(itau,idate)

    ! JFM: change to work properly in adjoint mode
    itautmp = itau - ninc
    call tau2date(itautmp,idtmp)
    !
    ! set switches
    !
    newyr = ( idate(1) /= idtmp(1) )
    newmonth = ( idate(2) /= idtmp(2) )
    newday = ( idate(3) /= idtmp(3) )
    newhour(:)  = ( idate(4) /= idtmp(4) )
    newsrun = .false.

    if ( cdebug ) then
       call tstamp(kdebug,itau,'time incremented')
    endif

    ! ok
    status = 0

  end subroutine inctime

  function isDayBoundary(itaux)
    !-----------------------------------------------------------------------
    !
    !     purpose
    !     -------
    !     given an itau, return whether it is at a day boundary, i.e., whether
    !     UTC time is 00:00:00
    !
    !-----------------------------------------------------------------------

    implicit none

    ! input/output
    integer,intent(in)              :: itaux
    logical                         :: isDayBoundary

    ! local
    integer, dimension(6)           :: idatex

    call tau2date(itaux, idatex)
    isDayBoundary = all(idatex(4:6) .eq. (/0, 0, 0/))

  end function isDayBoundary

  subroutine tau2date(itaux,idatex)
    !-----------------------------------------------------------------------
    !**** tau2date
    !
    !     purpose
    !     -------
    !     calculate date from given time in seconds
    !
    !     parameters
    !     ----------
    !     on input : itaux contains date/time in seconds
    !     on output: idatex contains date in year,month,day,hour,min,sec
    !
    !     dependencies
    !     ------------
    !     icalendo determines the type calendar used for the conversion
    !     iyear0 is the reference year for the calculation
    !     julday0 is the reference julian day for the calculation
    !
    !     externals
    !     ---------
    !     subroutines: caldat
    !     funtions:    julday
    !-----------------------------------------------------------------------
    use dims, only : icalendo, iyear0, julday0

    implicit none

    ! input/output
    integer,intent(in)               :: itaux
    integer,dimension(6),intent(out) :: idatex

    ! local
    integer :: julian
    integer :: idayy,iyeary,idummy
    !
    ! compute time (hour,min,sec) and number of days
    !
    idatex(6)=mod(itaux,60)
    idatex(5)=mod(itaux/60,60)
    idatex(4)=mod(itaux/3600,24)
    idayy=itaux/86400
    !
    ! permanent 360 year calendar with 30 days in each month
    !
    if ( icalendo == 1 ) then
       idatex(3)=mod(idayy,30)+1
       idatex(2)=mod(idayy/30,12)+1
       idatex(1)=iyear0+idayy/360
       !
       ! real calendar
       !
    else if ( icalendo == 2 ) then
       julian=julday0+idayy
       call caldat(julian,idatex(2),idatex(3),idatex(1))
       !
       ! permanent 365 day year calendar
       !
    else if ( icalendo == 3 ) then
       iyeary=idayy/365
       idatex(1)=iyear0+iyeary
       ! use Jan 1, 1981 as a year containing 365 days and add doy
       julian=julday(1,1,1981)+idayy-iyeary*365
       call caldat(julian,idatex(2),idatex(3),idummy)
       !
       ! permanent leap year calendar
       !
    else if ( icalendo == 4 ) then
       iyeary=idayy/366
       idatex(1)=iyear0+iyeary
       ! use Jan 1, 1980 as a year containing 366 days and add doy
       julian=julday(1,1,1980)+idayy-iyeary*366
       call caldat(julian,idatex(2),idatex(3),idummy)
       !
       ! illegal option icalendo
       !
    else
       write(*,*) ' tau2date: ERROR while computing date'
       write(*,*) ' tau2date: Illegal calendar type'
       write(*,*) '           icalendo = ',icalendo
       stop
    end if

  end subroutine tau2date

  subroutine date2tau_int(idatex,itaux)
    !-----------------------------------------------------------------------
    !**** date2tau
    !
    !     purpose
    !     -------
    !     calculate time in seconds from given date
    !
    !     parameters
    !     ----------
    !     on input : idatex contains date in year,month,day,hour,min,sec
    !     on output: itaux contains date/time in seconds
    !
    !     dependencies
    !     ------------
    !     icalendo determines the type calendar used for the conversion
    !     iyear0 is the reference year for the calculation
    !     julday0 is the reference julian day for the calculation
    !
    !     externals
    !     ---------
    !     funtions: julday
    !-----------------------------------------------------------------------
    use dims, only : icalendo, iyear0, julday0, kmain

    implicit none

    ! input/output
    integer,dimension(6),intent(in) :: idatex
    integer,intent(out)             :: itaux

    ! local
    integer :: idaysec
    !
    ! compute the seconds the day is old
    !
    idaysec=idatex(6)+idatex(5)*60+idatex(4)*3600
    !
    ! permanent 360 year calendar with 30 days in each month
    !
    if ( icalendo == 1 ) then
       itaux=idaysec+(idatex(3)-1)*86400+(idatex(2)-1)*2592000 &
            +(idatex(1)-iyear0)*31104000
       !
       ! real calendar
       !
    else if ( icalendo == 2 ) then
       itaux=86400*(julday(idatex(2),idatex(3),idatex(1))-julday0)+idaysec
       !
       ! permanent 365 day year calendar
       !
    else if ( icalendo == 3 ) then
       itaux=86400*(julday(idatex(2),idatex(3),1981)-julday(1,1,1981))  &
            +(idatex(1)-iyear0)*365*86400+idaysec
       !
       ! permanent leap year calendar
       !
    else if ( icalendo == 4 ) then
       itaux=86400*(julday(idatex(2),idatex(3),1980)-julday(1,1,1980))  &
            +(idatex(1)-iyear0)*366*86400+idaysec
       !
       ! illegal option icalendo
       !
    else
       write(kmain,*) ' date2tau: ERROR while computing time'
       write(kmain,*) ' date2tau: Illegal calendar type'
       write(kmain,*) '           icalendo = ',icalendo
       stop
    end if

  end subroutine date2tau_int

  subroutine date2tau_short(idatex_short,itaux)
    !-----------------------------------------------------------------------
    !**** date2tau
    !
    !     purpose
    !     -------
    !     calculate time in seconds from given date
    !
    !     parameters
    !     ----------
    !     on input : idatex contains date in year,month,day,hour,min,sec
    !     on output: itaux contains date/time in seconds
    !
    !     dependencies
    !     ------------
    !     icalendo determines the type calendar used for the conversion
    !     iyear0 is the reference year for the calculation
    !     julday0 is the reference julian day for the calculation
    !
    !     externals
    !     ---------
    !     funtions: julday
    !-----------------------------------------------------------------------
    use dims, only : icalendo, iyear0, julday0, kmain

    implicit none

    ! input/output
    integer(2),dimension(6),intent(in) :: idatex_short
    integer,intent(out)             :: itaux

    ! local
    integer :: idaysec
    integer,dimension(6) :: idatex

    idatex(:) = idatex_short(:)
    !
    ! compute the seconds the day is old
    !
    idaysec=idatex(6)+idatex(5)*60+idatex(4)*3600
    !
    ! permanent 360 year calendar with 30 days in each month
    !
    if ( icalendo == 1 ) then
       itaux=idaysec+(idatex(3)-1)*86400+(idatex(2)-1)*2592000 &
            +(idatex(1)-iyear0)*31104000
       !
       ! real calendar
       !
    else if ( icalendo == 2 ) then
       itaux=86400*(julday(idatex(2),idatex(3),idatex(1))-julday0)+idaysec
       !
       ! permanent 365 day year calendar
       !
    else if ( icalendo == 3 ) then
       itaux=86400*(julday(idatex(2),idatex(3),1981)-julday(1,1,1981))  &
            +(idatex(1)-iyear0)*365*86400+idaysec
       !
       ! permanent leap year calendar
       !
    else if ( icalendo == 4 ) then
       itaux=86400*(julday(idatex(2),idatex(3),1980)-julday(1,1,1980))  &
            +(idatex(1)-iyear0)*366*86400+idaysec
       !
       ! illegal option icalendo
       !
    else
       write(kmain,*) ' date2tau: ERROR while computing time'
       write(kmain,*) ' date2tau: Illegal calendar type'
       write(kmain,*) '           icalendo = ',icalendo
       stop
    end if

  end subroutine date2tau_short

!  integer function get_day(mm,dd,mlen)
!    !
!    ! returns day number (from 1 January)
!    !
!    implicit none
!    ! input/output
!    integer,intent(in) :: mm   ! month in year
!    integer,intent(in) :: dd   ! day in month, year
!    integer,intent(in),dimension(12)  :: mlen  !lengt of months (days)
!    ! local
!    integer    :: m

!    get_day = 0
!    do m=1,mm-1
!       get_day = get_day + mlen(m)
!    enddo
!    get_day = get_day + dd

!  end function get_day



  integer function julday(mm,id,iy)
    !-----------------------------------------------------------------------
    !**** julday
    !
    !    purpose
    !    -------
    !    calculate julian day from given date
    !
    !    parameters
    !    ----------
    !    on input : mm, id, iyyy contain month, day and year
    !    on output: julday contains the julian day
    !
    !    dependencies
    !    ------------
    !    julday0 is the reference julian day for the calculation
    !
    !    reference
    !    ---------
    !    J. Meeuws, "Astronomical formulea for calculators" 19xx
    !-----------------------------------------------------------------------
    implicit none

    ! input, output
    integer,intent(in) :: mm  ! month
    integer,intent(in) :: id  ! day
    integer,intent(in) :: iy  ! year

    ! local
    integer,parameter  :: igreg=15+31*(10+12*1582)
    integer            :: julday0, jy, jm, ja, iyyy

    ! handle dates before 0 AD
    !
    iyyy=iy
    if ( iy == 0 ) then
       stop 'julday:  ERROR invalid year 0 AD'
    end if
    if ( iy < 0 ) then
       iyyy=iy+1
    end if
    !
    !calculate julian day from date in gregorian calendar
    !
    if ( mm > 2 ) then
       jy=iyyy
       jm=mm+1
    else
       jy=iyyy-1
       jm=mm+13
    end if
    julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
    !
    !handle julian calender
    !
    if ( id+31*(mm+12*iyyy) >= igreg ) then
       ja=int(0.01*jy)
       julday=julday+2-ja+int(0.25*ja)
    end if

  end function julday



  subroutine caldat(julian,mm,id,iyyy)
    !-----------------------------------------------------------------------
    !**** caldat
    !
    !     purpose
    !     -------
    !     calculate date from given julian day
    !
    !     parameters
    !     ----------
    !     on input : julday contains the julian day
    !     on output: mm, id, iyyy contain month, day and year
    !
    !     dependencies
    !     ------------
    !     julday0 is the reference julian day for the calculation
    !
    !     reference
    !     ---------
    !     J. Meeuws, "Astronomical formulea for calculators" 19xx
    !-----------------------------------------------------------------------
    implicit none

    ! input/output
    integer,intent(in)  :: julian
    integer,intent(out) :: mm
    integer,intent(out) :: id
    integer,intent(out) :: iyyy

    ! local
    integer,parameter   :: igreg=2299161
    integer             :: jalpha, ja, jb, jc, jd, je
    !
    ! handle gregorian and julian date
    !
    if ( julian >= igreg )then
       jalpha=int(((julian-1867216)-0.25)/36524.25)
       ja=julian+1+jalpha-int(0.25*jalpha)
    else
       ja=julian
    end if
    jb=ja+1524
    jc=int(6680.+((jb-2439870)-122.1)/365.25)
    jd=365*jc+int(0.25*jc)
    je=int((jb-jd)/30.6001)
    id=jb-jd-int(30.6001*je)
    mm=je-1
    if ( mm > 12 ) mm=mm-12
    iyyy=jc-4715
    if ( mm > 2 ) iyyy=iyyy-1
    !
    ! handle dates before 0 AD
    !
    if ( iyyy <= 0 ) iyyy=iyyy-1

  end subroutine caldat



!  subroutine calc_sm( mlen, sec_day, sec_month, sec_year )
!    !
!    !
!    !
!    use dims, only : icalendo, idate

!    implicit none

!    ! input/output
!    real,intent(out) :: sec_day     ! # seconds in day
!    real,intent(out) :: sec_month   ! # seconds in current month
!    real,intent(out) :: sec_year    ! # seconds in current year
!    integer,intent(out),dimension(12):: mlen ! days per month (current year)

!    ! start

!    sec_day=86400.
!    mlen(1)=31
!    mlen(2)=28
!    mlen(3)=31
!    mlen(4)=30
!    mlen(5)=31
!    mlen(6)=30
!    mlen(7)=31
!    mlen(8)=31
!    mlen(9)=30
!    mlen(10)=31
!    mlen(11)=30
!    mlen(12)=31 ! only for regular year

!    !
!    ! calender option
!    !
!    sec_year=365.*sec_day
!    if ( icalendo == 1 ) then
!       mlen(:)=30
!       sec_year=360.*sec_day
!    end if
!    if ( icalendo == 4 ) then
!       mlen(2)=29
!       sec_year=366.*sec_day
!    end if
!    if ( icalendo == 2 .and. (mod(idate(1),4) == 0) .and.   &
!         (mod(idate(1),100) /= 0) .or. (mod(idate(1),400) == 0) ) then
!       mlen(2)=29
!       sec_year=366.*sec_day
!    end if
!    sec_month=sec_day*mlen(idate(2))

!    !write(*,*) 'calc_sm: sec_month',sec_month

!  end subroutine calc_sm

  function idate2ddate(idate)
    ! Calculates fractional years
    implicit none

    integer, intent(in) :: idate(6)
    real                :: idate2ddate
    real                :: seconds_in_year, seconds_elapsed
    integer             :: itau_ref, itau

    if (isleap(idate(1)) == 2) then
        seconds_in_year = 366.0 * 86400.0
    else
        seconds_in_year = 365.0 * 86400.0
    end if

    call date2tau((/idate(1), 1, 1, 0, 0, 0/), itau_ref)
    call date2tau(idate, itau)

    seconds_elapsed = dble(itau-itau_ref)

    idate2ddate = dble(idate(1)) + seconds_elapsed/seconds_in_year

  end function idate2ddate

  function isleap(y)
    integer, intent(in) :: y
    integer             :: isleap
    logical             :: leap

    isleap = 1
    leap = (mod(y,4) == 0 .and. mod(y,100) /= 0) .or. mod(y,400) == 0
    if (leap) isleap = 2
  end function isleap


  subroutine dayl(day,daylen,jdim,lat_start,dlat)
    !
    !***  calculates daylength (hours) depending on
    !***  latitude (phi) and day of year (day)
    !
    !     programmed by:
    !     implemented by: fd IMAU Tue Feb 27 17:51:57 MET 1997
    !     modified by MK for zoom version may 2001

    !     purpose
    !     -------
    !     calculates daylength
    !
    !     interface
    !     ---------
    !     call dayl(day,daylen,jdim,lat_start,dlat)
    !
    !       day          : the day of the year based on 365 days a year
    !       daylen(jdim) : length of day at jdim latitudes
    !       lat_start    : first latitude (degrees from -90 to +90)
    !                      at the southernmost edge
    !       dlat         : increment
    !
    !     method
    !     ------
    !     none
    !
    !     external
    !     ---------
    !     none
    !
    !     reference
    !     ---------
    !
    !-------------------------------
    implicit none

    ! input/output
    integer,intent(in)               :: day    ! the day of the year
    integer,intent(in)               :: jdim   ! dimension of daylen
    real,intent(out),dimension(jdim) :: daylen ! length of day at jdim lats
    real,intent(in)                  :: lat_start ! first latitude
    real,intent(in)                  :: dlat   ! latitude increment

    ! local
    real                             :: nj,dj,phi,td,a,phix,xh,pi
    integer                          :: j,idayy=365

    ! start
    pi = acos(-1.0)
    dj= dlat

    do j=1,jdim
       phix = (lat_start+(j-0.5)*dj)*pi/180.
       td = -float(mod(day+10,idayy))*2.*pi/idayy
       a = cos(td)*pi/180.*23.45
       xh = tan(a)*tan(phix)
       if ( abs(xh) <= 1. ) then
          ! CMK BUG: N-S reversal removed dec2004
          daylen(j) = 24*(1-acos(-xh)/pi)
       else
          if ( xh <= -1 ) daylen(j) = 24.0
          if ( xh >=  1 ) daylen(j) = 0.0
       end if
    end do ! j

  end subroutine dayl



  !----------------------------------------------------------------------
  ! write time stamp and msg on unit kunit
  !----------------------------------------------------------------------

  subroutine tstamp( kunit, itaux, msg )

    use GO, only : gol, goPr

    ! --- in/out -----------------------

    integer,intent(in)          :: kunit  ! unit to write "tstamp" to  <--- ignored

    integer,intent(in)          :: itaux  ! "tstamp" time
    character(len=*),intent(in) :: msg    ! "tstamp" message

    ! --- begin ---------------------

    call wrtgol_tstamp( itaux, msg ); call goPr

  end subroutine tstamp


  ! ***


  subroutine wrtgol_tstamp( itaux, msg )

    use GO, only : gol, goPr

    ! --- in/out -----------------------

    integer,intent(in)          :: itaux  ! "tstamp" time
    character(len=*),intent(in) :: msg    ! "tstamp" message

    ! --- local ----------------------

    integer,dimension(6) :: idatex
    character(len=24)    :: cdate

    ! --- begin ---------------------

    ! convert from seconds to year/month/etc:
    call tau2date(itaux,idatex)

    ! write in characters:
    call chardate(idatex,cdate)

    ! depricated ...
    !write (kunit,'(a1,a24,a1,a)') ' ',cdate,' ',msg

    ! display:
    write (gol,'(a24," ",a)') cdate, trim(msg)

  end subroutine wrtgol_tstamp


  ! ***


  subroutine new_valid_timestep( dtime, nread, cfl_outputstep)

    use dims, only        : nregions, tref

    ! --- in/out -------------------------

    integer, intent(inout)  :: dtime   ! current timestep
                                       ! to be replaced with a valid new timestep
    integer, intent(in)     :: nread   ! nread (e.g. 3hr) should be a
                                       ! multiple of the new 'valid' timestep

    integer, intent(in)     :: cfl_outputstep   ! choose times

    ! --- begin -----------------------------

    ! loop until time is largest multiple of nread, clf_outputstep, and 2*tref
    do
      dtime = dtime-1
      if (mod(nread,dtime) == 0 .and. mod(cfl_outputstep,dtime) == 0 .and.&
          mod(dtime,maxval(2*tref(1:nregions))) == 0) exit
      if (dtime < maxval(2*tref(1:nregions))) then
        write (gol,'("no valid timestep found:")'); call goPr
        write (gol,'("  dtime (s)          : ",i8)') dtime; call goErr
        write (gol,'("  nread (s)          : ",i8)') nread; call goErr
        write (gol,'("  cfl_outputstep (s) : ",i8)') cfl_outputstep; call goErr
        write (gol,'("STOP in ",a)') 'datetime/new_valid_timestep'; call goErr; stop
      end if
    enddo

  end subroutine new_valid_timestep


end module datetime
