!###############################################################################
!
! time tools:
!  o time window
!  o time profile : partitioning of a window into intervals
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!###############################################################################

module GO_Time

  use GO_Print, only : gol, goPr, goErr
  use GO_Date , only : TDate

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  T_Time_Window
  public  ::  Time_Window_Init, Time_Window_Done

  public  ::  wrtgol

  public  ::  T_Time_Profile, T_Time_Period
  public  ::  Time_Profile_Init, Time_Profile_Done
  public  ::  Time_Profile_Index
  public  ::  Time_Profile_Covers
  public  ::  Time_Profile_Length
  public  ::  Time_Profile_Overlap_Matrix
  public  ::  operator(.in.), operator(.overlaps.)


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'GO_Time'


  ! --- types ------------------------------------

  ! time window:
  type T_Time_Window
    ! start and end time:
    type(TDate)                         ::  t1, t2
  end type T_Time_Window

  ! time description of a single period:
  type T_Time_Period
    ! start, end:
    type(TDate)                        ::  t1, t2
    ! mid:
    type(TDate)                        ::  tmid
  end type T_Time_Period

  ! total profile:
  type T_Time_Profile
    ! time window:
    type(TDate)                         ::  t1, t2
    ! resolution key:
    character(len=32)                   ::  reskey
    ! number of time periods within window:
    integer                             ::  n_period
    ! information per period:
    type(T_Time_Period), pointer        ::  period(:)
  end type T_Time_Profile


  ! --- interfaces -------------------------------

  interface wrtgol
    module procedure Time_Period_wrtgol
  end interface

  interface operator(.in.)
    module procedure Time_Profile_point_in
    module procedure Time_Profile_window_in
  end interface

  interface operator(.overlaps.)
    module procedure Time_Window_overlap
  end interface

  interface Time_Window_Init
    module procedure Time_Window_Init_int
    module procedure Time_Window_Init_date
  end interface

  interface Time_Profile_Init
    module procedure Time_Profile_Init_edges
    module procedure Time_Profile_Init_res
  end interface

contains


  ! ********************************************************************
  ! ***
  ! *** time window
  ! ***
  ! ********************************************************************


  subroutine Time_Window_Init_int( tw, time6_1, time6_2, status )

    use GO_Date, only : NewDate

    ! --- in/out ---------------------------------

    type(T_Time_Window), intent(out)          ::  tw
    integer, intent(in)                       ::  time6_1(6), time6_2(6)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Window_Init_int'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! store:
    tw%t1 = NewDate( time6=time6_1 )
    tw%t2 = NewDate( time6=time6_2 )

    ! ok:
    status = 0

  end subroutine Time_Window_Init_int

  subroutine Time_Window_Init_date( tw, time1, time2, status )

    use GO_Date, only : NewDate

    ! --- in/out ---------------------------------

    type(T_Time_Window), intent(out)          ::  tw
    type(TDate), intent(in)                   ::  time1, time2
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Window_Init_date'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! store:
    tw%t1 = time1
    tw%t2 = time2

    ! ok:
    status = 0

  end subroutine Time_Window_Init_date


  ! ***


  subroutine Time_Window_Done( tw, status )

    ! --- in/out ---------------------------------

    type(T_Time_Window), intent(inout)          ::  tw
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Window_Done'

    ! --- begin ----------------------------------

    ! ok:
    status = 0

  end subroutine Time_Window_Done


  ! ********************************************************************
  ! ***
  ! *** time period
  ! ***
  ! ********************************************************************


  subroutine Time_Period_wrtgol( msg, tp )

    use GO_Date, only : wrtgol

    ! --- in/out ---------------------------------

    character(len=*), intent(in)                ::  msg
    type(T_Time_Profile), intent(inout)         ::  tp

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Period_wrtgol'

    ! --- begin ----------------------------------

    ! fill gol with  :  msg [ t1, t2 ]
    call wrtgol( msg, tp%t1, ', ', tp%t2 )

  end subroutine Time_Period_wrtgol


  ! ********************************************************************
  ! ***
  ! *** time profile
  ! ***
  ! ********************************************************************

  subroutine Time_Profile_Init_edges(tp, tbeg, tend, status)

    use Go_Date,    only : TDate, operator(-), operator(+), operator(/), TIncrDate

    type(T_Time_Profile), intent(out)   :: tp
    type(TDate), intent(in)             :: tbeg(:), tend(:)
    integer, intent(out)                :: status

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_Init_edges'

    integer         :: np, i
    type(TIncrDate) :: dt

    np = size(tbeg, 1)

    tp%t1 = tbeg(1)
    tp%t2 = tend(np)

    tp%n_period = np
    allocate(tp%period(np))
    do i = 1, np
      tp%period(i)%t1 = tbeg(i)
      tp%period(i)%t2 = tend(i)
      dt = tend(i) - tbeg(i)
      dt = dt/2
      tp%period(i)%tmid = tbeg(i) + dt
    end do
    ! set the resolution to something that does not matter
    tp%reskey = 'arbitrary'

    status = 0

  end subroutine Time_Profile_Init_edges

  subroutine Time_Profile_Init_res( tp, tw, reskey, status )

    use GO_Date, only : NewDate, IncrDate, MidNight, Get_End_Of
    use GO_Date, only : operator(+), operator(-), operator(/), operator(>=)
    use GO_Date, only : rTotal
    use GO_String, only : goSplitLine

    use GO_Date, only : wrtgol

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(out)         ::  tp
    type(T_Time_Window), intent(in)           ::  tw
    character(len=*), intent(in)              ::  reskey   ! time resolution key: 'monthly' ...
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_Init_res'

    ! --- local ----------------------------------

    character(len=20)   ::  sres, sn
    integer             ::  nres
    integer             ::  i, j
    integer             ::  np

    ! --- begin ----------------------------------

    ! copy time window:
    tp%t1 = tw%t1
    tp%t2 = tw%t2

    ! store resolution key:
    tp%reskey = trim(reskey)

    ! format is:  <resolution>[+<n>]
    ! with <n>=1 by default; examples: "monthly", "daily", "daily+8"

    ! split if necessary:
    call goSplitLine( trim(reskey), sres, '+', sn, status )
    IF_NOTOK_RETURN(status=1)
    if ( len_trim(sn) > 0 ) then
      read (sn,'(i4)',iostat=status) nres
      if (status/=0) then
        write (gol,'("could not read integer number from time resolution key : ",a)') trim(reskey); call goErr
        TRACEBACK; status=1; return
      end if
    else
      nres = 1
    end if

    ! swith between different profile resolutions:
    select case ( trim(sres) )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    case ( 'monthly' )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! first guess of number of time periods; round to upper value:
      np = ceiling( real( (tw%t2%year-tw%t1%year) * 12 + (tw%t2%month-tw%t1%month+1) )/real(nres) )

      ! storage for time info for each period:
      allocate( tp%period(np) )

      ! loop over periods:
      i = 0
      do
        ! increase counter:
        i = i + 1
        ! fill start time:
        if ( i == 1 ) then
          tp%period(i)%t1 = tw%t1
        else
          tp%period(i)%t1 = tp%period(i-1)%t2
        end if
        ! fill end time:
        tp%period(i)%t2 = Get_End_Of( tp%period(i)%t1, 'month' )
        if ( nres > 1 ) then
          do j = 1, nres-1
            tp%period(i)%t2 = Get_End_Of( tp%period(i)%t2, 'month' )
          end do
        end if
        ! fill current number:
        tp%n_period = i
        ! finished ?
        if ( tp%period(i)%t2 >= tp%t2 ) then
          ! reset to end of window:
          tp%period(i)%t2 = tp%t2
          ! leave:
          exit
        end if
      end do

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    case ( 'daily' )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! first guess of number of time periods; round to upper value:
      np = ceiling( rTotal( tw%t2-tw%t1, 'day' )/real(nres) )

      ! storage for time info for each period:
      allocate( tp%period(np) )

      ! loop over periods:
      i = 0
      do
        ! increase counter:
        i = i + 1
        ! fill start time:
        if ( i == 1 ) then
          tp%period(i)%t1 = tw%t1
        else
          tp%period(i)%t1 = tp%period(i-1)%t2
        end if
        ! fill end time:
        if ( i == 1 ) then
          tp%period(i)%t2 = Get_End_Of( tp%period(i)%t1, 'day' )
          ! if ( nres > 1 ) tp%period(i)%t2 = tp%period(i)%t1 + IncrDate(day=nres-1)
          if ( nres > 1 ) tp%period(i)%t2 = tp%period(i)%t2 + IncrDate(day=nres-1)
        else
          tp%period(i)%t2 = tp%period(i)%t1 + IncrDate(day=nres)
        end if
        ! fill current number:
        tp%n_period = i
        ! finished ?
        if ( tp%period(i)%t2 >= tp%t2 ) then
          ! reset to end of window:
          tp%period(i)%t2 = tp%t2
          ! leave:
          exit
        end if
      end do

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      write (gol,'("unsupported time profile key : ",a)') trim(tp%reskey); call goErr
      TRACEBACK; status=1; return

    end select

    ! set mid times:
    do i = 1, tp%n_period
      ! mid of interval:
      tp%period(i)%tmid = tp%period(i)%t1 + ( tp%period(i)%t2 - tp%period(i)%t1 )/2
    end do

    !! debug ...
    !write (gol,'("xxx ",a)') trim(reskey); call goPr
    !do i = 1, tp%n_period
    !  call wrtgol( 'xxx : ', tp%period(i)%t1, '  -  ', tp%period(i)%t2, '   ;   ', tp%period(i)%tmid ); call goPr
    !end do

    ! ok
    status = 0

  end subroutine Time_Profile_Init_res


  ! ***


  subroutine Time_Profile_Done( tp, status )

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(inout)       ::  tp
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! clear:
    deallocate( tp%period )

    ! ok
    status = 0

  end subroutine Time_Profile_Done


  ! ***


  logical function Time_Profile_point_in( t, tp )

    use GO_Date, only : operator(<=)

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(in)          ::  tp
    type(TDate), intent(in)                   ::  t

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_point_in'

    ! --- begin ----------------------------------

    ! in window ?
    Time_Profile_point_in = (tp%t1 <= t) .and. (t <= tp%t2)

  end function Time_Profile_point_in

  logical function Time_Profile_window_in( tp, tw ) ! if time profile fits inside time window

    use GO_Date, only : operator(<=)

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(in)          ::  tp
    type(T_Time_Window), intent(in)           ::  tw

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_window_in'

    ! --- begin ----------------------------------

    ! in window ?
    Time_Profile_window_in = (tw%t1 <= tp%t1) .and. (tp%t2 <= tw%t2)

  end function Time_Profile_window_in

  logical function Time_Window_overlap(tw1, tw2)
    ! If two time windows overlap
    use GO_Date, only   : operator(<=)

    type(T_Time_Window), intent(in) :: tw1, tw2

    character(len=*), parameter   ::  rname = mname//'/Time_Window_overlap'

    if (tw1%t2 <= tw2%t1 .or. tw2%t2 <= tw1%t1) then
      Time_Window_overlap = .false.
    else
      Time_Window_overlap = .true.
    end if

  end function Time_Window_overlap

  ! ***


  ! search for first period in profile that includes t ;
  ! by default, boundaries match too, thus  t in [t1,t]_i ,
  ! unless 'at_left_side' flag is set, then t in [t,t2]_i+1

  subroutine Time_Profile_Index( tp, t, ind, status, at_left_side )

    use GO_Date, only : operator(>=), operator(<=), operator(==), wrtgol

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(in)          ::  tp
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  ind
    integer, intent(out)                      ::  status
    logical, intent(in), optional             ::  at_left_side

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_Index'

    ! --- local ----------------------------------

    integer                  ::  i

    ! --- begin ----------------------------------

    ! check ...
    if ( tp%n_period <= 0 ) then
      write (gol,'("no periods found; using undefined time profile ?")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! init result:
    ind = -1

    ! loop over periods:
    do i = 1, tp%n_period
      ! match ?
      if ( (t >= tp%period(i)%t1) .and. (t <= tp%period(i)%t2) ) then
        ! store:
        ind = i
        ! leave:
        exit
      end if
    end do

    ! check ...
    if ( ind < 0 ) then
      write (gol,'("could not find time period:")'); call goErr
      call wrtgol( '  target time : ', t ); call goErr
      write (gol,'("  periods : ",i6)') tp%n_period; call goErr
      do i = 1, tp%n_period
        call wrtgol( '    ', tp%period(i)%t1, ' - ', tp%period(i)%t2  ); call goErr
      end do
      TRACEBACK; status=1; return
    end if

    ! at right border ?
    if ( t == tp%period(ind)%t2 ) then
      ! argument provided ?
      if ( present(at_left_side) ) then
        ! should be at left border if possible ?
        if ( at_left_side ) then
          ! set to next interval if possible:
          if ( ind < tp%n_period ) ind = ind+1
        end if
      end if
    end if

    ! ok
    status = 0

  end subroutine Time_Profile_Index

  ! ***

  ! Given two time profiles, tp1 with n1 intervals and tp2 with n2 intervals,
  ! return a n1xn2 matrix whose element (i,j) is
  ! (a) If per_unit_time:
  !     1 if tp1(i) has any overlap with tp2(j), 0 otherwise
  ! (b) If not per_unit_time:
  !     Fraction of tp1(i) that is covered by tp2(j)

  subroutine Time_Profile_Overlap_Matrix(tp1, tp2, ov_mat, per_unit_time)

    use GO_Date, only : rTotal
    use GO_Date, only : operator(<=), operator(>=), operator(-), min, max

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(in)    :: tp1, tp2
    real, intent(out)                   :: ov_mat(tp1%n_period, tp2%n_period)
    logical, intent(in), optional       :: per_unit_time

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_Overlap_Matrix'

    ! --- local ----------------------------------

    integer                             :: i1, i2
    logical                             :: per_time
    type(TDate)                         :: t_lb, t_ub
    real                                :: dt_num, dt_den

    if (present(per_unit_time)) then
      per_time = per_unit_time
    else
      per_time = .true.
    end if

    ov_mat = 0.0

    do i1 = 1, tp1%n_period
      do i2 = 1, tp2%n_period
        ! calculate how much of tp1(i1) is covered by tp2(i2)
        if (tp2%period(i2)%t2 <= tp1%period(i1)%t1 .or. tp2%period(i2)%t1 >= tp1%period(i1)%t2) cycle
        t_lb = max(tp2%period(i2)%t1, tp1%period(i1)%t1)
        t_ub = min(tp2%period(i2)%t2, tp1%period(i1)%t2)
        dt_num = rTotal(t_ub - t_lb, 'sec')
        if (per_time) then
          dt_den = rTotal(tp2%period(i2)%t2 - tp2%period(i2)%t1, 'sec')
        else
          dt_den = rTotal(tp1%period(i1)%t2 - tp1%period(i1)%t1, 'sec')
        end if
        ov_mat(i1, i2) = dt_num/dt_den
      end do
    end do

  end subroutine Time_Profile_Overlap_Matrix

  ! ***


  !
  ! return info on how to distribute an interval [t1,t2] over the profile:
  !   n       :  number of periods (partly) covered by the interval
  !   ip(n)   :  period indices
  !   ff(n)   :  fraction of interval [t1,t2] that covers the period
  !   gg(n)   :  fraction of the period that is covered by [t1,t2]
  !

  subroutine Time_Profile_Covers( tp, t1, t2, n, ip, ff, gg, status )

    use GO_Date, only : operator(<=), operator(==), operator(-), min, max
    use GO_Date, only : rTotal

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(in)          ::  tp
    type(TDate), intent(in)                   ::  t1, t2
    integer, intent(out)                      ::  n
    integer, pointer, intent(out)             ::  ip(:)
    real, pointer, intent(out)                ::  ff(:)
    real, pointer, intent(out)                ::  gg(:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_Covers'

    ! --- local ----------------------------------

    integer                  ::  i1, i2, i_period, i
    real                     ::  dt12, dtp, dts

    ! --- begin ----------------------------------

    ! check ...
    if ( tp%n_period <= 0 ) then
      write (gol,'("no periods found; using undefined time profile ?")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! no coverage at al, or empty interval ?
    if ( (t2 <= tp%t1) .or. (tp%t2 <= t1) .or. (t1==t2) ) then
      ! no matches:
      n = 0
    else
      ! period including t1:
      if ( t1 <= tp%t1 ) then
        ! first period:
        i1 = 1
      else
        ! search interval including t1 or with t1 at left side:
        call Time_Profile_Index( tp, t1, i1, status, at_left_side=.true. )
        IF_NOTOK_RETURN(status=1)
      end if
      ! period including t2:
      if ( tp%t2 <= t2 ) then
        i2 = tp%n_period
      else
        ! search interval including t2 or with t2 at right side:
        call Time_Profile_Index( tp, t2, i2, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! total number:
      n = i2 - i1 + 1
      ! storage:
      if ( associated(ip) ) deallocate(ip)
      allocate( ip(n) )
      if ( associated(ff) ) deallocate(ff)
      allocate( ff(n) )
      if ( associated(gg) ) deallocate(gg)
      allocate( gg(n) )
      ! length of interval:
      dt12 = rTotal( t2-t1, 'sec' )
      ! loop:
      do i_period = i1, i2
        ! index in output arrays:
        i = i_period - i1 + 1
        ! store index:
        ip(i) = i_period
        ! length of perid:
        dtp = rTotal( tp%period(i_period)%t2 - tp%period(i_period)%t1, 'sec' )
        ! length of shared interval:
        dts = rTotal( min(tp%period(i_period)%t2,t2) - max(tp%period(i_period)%t1,t1), 'sec' )
        ! faction of [t1,t2] covering period:
        ff(i) = dts/dt12
        ! fraction of period covered by [t1,t2] :
        gg(i) = dts/dtp
      end do
    end if

    ! ok
    status = 0

  end subroutine Time_Profile_Covers


  ! ***


  !
  ! Return length of window or single period in specified unit ('day','sec',...)
  !

  subroutine Time_Profile_Length( tp, unit, r, status, period )

    use GO_Date, only : rTotal, operator(-)

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(in)    ::  tp
    character(len=*), intent(in)        ::  unit
    real, intent(out)                   ::  r
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  period

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Time_Profile_Length'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! single period or total window ?
    if ( present(period) ) then
      ! check ...
      if ( (period < 1) .or. (period > tp%n_period) ) then
        write (gol,'("requested period ",i6," out of range 1 .. ",i6)') period, tp%n_period; call goErr
        TRACEBACK; status=1; return
      end if
      ! return difference:
      r = rTotal( tp%period(period)%t2 - tp%period(period)%t1, unit )
    else
      ! return difference:
      r = rTotal( tp%t2 - tp%t1, unit )
    end if

    ! ok
    status = 0

  end subroutine Time_Profile_Length


end module GO_Time
