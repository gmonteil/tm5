!
!   num_quad  -  numerical quadrature
!

module num_quad

  implicit none
  
  ! --- in/out -------------------
  
  private
  
  public   ::  IntervalQuad_Lin
  public   ::  IntervalQuad_Cos_Lin
  public   ::  IntervalQuad_Const
  public   ::  IntervalSum
  
  ! --- const ---------------------------------------
  
  character(len=*), parameter  ::  mname = 'num_quad'
  
  
contains


  !-------------------------------------------------------------------
  !
  ! Syntax:
  !
  !   call IntervalQuad_Lin( x, y, a, b, c, ilast )
  !
  ! Compute an approximation to the integral:
  !
  !    b
  !   int  y(x) dx
  !   x=a
  !
  ! Input is formed by two vectors x and y: y contains the function 
  ! values to the abscissa values contained in x. 
  ! The result, stored in c is computed using the trapezoid formula 
  ! throughout the domain of integration. 
  ! If a<x(1) or b>x(n) then the function is extended to a, resp. b, 
  ! assuming a constant value of y(1), resp. y(n).
  ! The index ilast should be set to 1 in the first call. If subsequent
  ! integrals on the same vector are to be computed, it should be kept
  ! unaltered between calls. This will speed up the lookup of values
  ! if the various integrals are ascending in x.
  !
  !--------------------------------------------------------------------

  subroutine IntervalQuad_Lin( x,y, a,b, c, ilast, status )

    use GO       , only : gol, goErr
    use num_tools, only : interval

    ! --- in/out ---------------------

    real, intent(in)       ::  x(:), y(:)
    real, intent(in)       ::  a, b
    real, intent(out)      ::  c
    integer, intent(inout) ::  ilast
    integer, intent(out)   ::  status

    ! --- const ---------------------

    character(len=*), parameter  ::  rname = mname//'/IntervalQuad_Lin'

    ! --- local -------------------

    integer         ::  n
    real            ::  ya, yb
    integer         ::  i, i1, i2, iflag

    ! --- begin ---------------------

    n = size(x)
    if ( size(y) /= n ) then
      write (gol,'("x and y should have same lengths:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(y)  : ",i6)') size(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    i1=ilast

    call interval(x,a,i1,iflag)
    if (iflag == -1) then              ! a < x(1) case
      c=y(1)*(x(1)-a)
    else if (iflag == 1) then          ! a > x(n) case
      c=-y(n)*(a-x(n))
    else
      ya=y(i1)+(a-x(i1))*(y(i1+1)-y(i1))/(x(i1+1)-x(i1))
      c=-0.5*(a-x(i1))*(y(i1)+ya)
    end if

    i2=i1
    call interval(x,b,i2,iflag)
    if (iflag == -1) then              ! b < x(1) case
      c=c-y(1)*(x(1)-b)
    else if (iflag == 1) then          ! b > x(n) case
      c=c+y(n)*(b-x(n))
    else
      yb=y(i2)+(b-x(i2))*(y(i2+1)-y(i2))/(x(i2+1)-x(i2))
      c=c+0.5*(b-x(i2))*(y(i2)+yb)
    end if

    if (i1 < i2) then
      do i=i1,i2-1
        c=c+0.5*(x(i+1)-x(i))*(y(i+1)+y(i))
      end do
    else
      do i=i2,i1-1
        c=c-0.5*(x(i+1)-x(i))*(y(i+1)+y(i))
      end do
    end if

    ilast=i2

    ! ok
    status = 0

  end subroutine IntervalQuad_Lin


  ! ===

  !-------------------------------------------------------------------
  !
  ! Syntax:
  !
  !   call lquad_cos( x, y, a, b, c, ilast )
  !
  ! Compute an approximation to the integral:
  !
  !    b
  !   int  y(x) cos(x) dx
  !   x=a
  !
  ! Input is formed by two vectors x and y: y contains the function 
  ! values to the abscissa values contained in x. 
  ! The result, stored in c is computed using the trapezoid formula 
  ! throughout the domain of integration. 
  ! If a<x(1) or b>x(n) then the function is extended to a, resp. b, 
  ! assuming a constant value of y(1), resp. y(n).
  ! The index ilast should be set to 1 in the first call. If subsequent
  ! integrals on the same vector are to be computed, it should be kept
  ! unaltered between calls. This will speed up the lookup of values
  ! if the various integrals are ascending in x.
  !
  ! Assuming y(x) = c + dx in the gridbox :
  !
  !    b                                                         b
  !   int (c +dx) cos(x) dx = [ sin(x) + d( x sin(x) + cos(x) ) ]
  !   x=a                                                        x=a
  !
  !--------------------------------------------------------------------

  subroutine IntervalQuad_Cos_Lin( x, y, a, b, c, ilast, status )

    use GO       , only : gol, goErr
    use num_tools, only : interval

    ! --- in/out --------------------------------

    real, intent(in)         ::  x(:), y(:)
    real, intent(in)         ::  a, b
    real, intent(out)        ::  c
    integer, intent(inout)   ::  ilast
    integer, intent(out)     ::  status

    ! --- const ---------------------

    character(len=*), parameter  ::  rname = mname//'/IntervalQuad_Cos_Lin'

    ! --- local ----------------------------------

    integer   ::  n
    real      ::  s, co
    integer   ::  i, i1, i2
    integer   ::  iflag

    ! --- begin ----------------------------------

    n = size(x)
    if ( size(y) /= n ) then
      write (gol,'("x and y should have same lengths:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(y)  : ",i6)') size(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    i1=ilast

    call interval(x,a,i1,iflag)
    if (iflag == -1) then              ! a < x(1) case
      c = y(1)*(sin(x(1))-sin(a))
    else if (iflag == 1) then          ! a > x(n) case
      c = -y(n)*(sin(a)-sin(x(n)))
    else
      ! int (co + s*x) cos(x) dx = INt (co cos(x) dx + Int s*x cos(x) dx 
      ! co sin(x) |  + s x sin(x) | + s cosx |
      s = (y(i1+1)-y(i1))/(x(i1+1)-x(i1))
      co = y(i1) - x(i1)*s
      c = (co + s*a)    *sin(a) - &
          (co + s*x(i1))*sin(x(i1)) + &
           s*(cos(a)-cos(x(i1)))
      c = -c   ! negative because outside a-b
    end if

    i2=i1
    call interval(x,b,i2,iflag)
    if (iflag == -1) then              ! b < x(1) case
      c = c - y(1)*(sin(x(1))-sin(b))
    else if (iflag == 1) then          ! b > x(n) case
      c = c + y(n)*(sin(b)-sin(x(n)))
    else
      s = (y(i2+1)-y(i2))/(x(i2+1)-x(i2))
      co = y(i2) - x(i2)*s
      c = c + (co + s*b)*sin(b) - &
          (co + s*x(i2))    *sin(x(i2)) + &
           s*(cos(b)-cos(x(i2)))
    end if

    if (i1 < i2) then
      do i=i1,i2-1
        s = (y(i+1)-y(i))/(x(i+1)-x(i))
        co = y(i) - x(i)*s
        c = c + (co + s*x(i+1))*sin(x(i+1)) - &
                (co + s*x(i))  *sin(x(i))   + &
           s*(cos(x(i+1))-cos(x(i)))
      end do
    else
      do i=i2,i1-1
        s = (y(i+1)-y(i))/(x(i+1)-x(i))
        co = y(i) - x(i)*s
        c = c - (co + s*x(i+1))*sin(x(i+1)) + &
                (co + s*x(i))  *sin(x(i))   - &
           s*(cos(x(i+1))-cos(x(i)))
      end do
    end if

    ilast=i2

    ! ok
    status = 0

  end subroutine IntervalQuad_Cos_Lin

      
  !-------------------------------------------------------------------
  ! 
  ! NAME
  !   Integral_const
  !
  ! INTERFACE
  !   subroutine IntervalQuad_Const( x, y, a, b, c, ilast )
  !     real, intent(in)        ::  x(:), y(:)
  !     real, intent(in)        ::  a, b
  !     real, intent(out)       ::  c
  !     integer, intent(inout)  ::  ilast
  !
  ! DESCRIPTION
  !   Compute integral over [a,b] for function specified by:
  !     o interval boundaries x(1),..,x(n),x(n+1)
  !     o function values y(1),..,y(n) ; constant in interval
  !
  !                  y(1)       y(2)        y(3)
  !                        +-----o-----+
  !                        |/ / / / / /+-----o----+
  !              +----o----+ / / / / / / / / / /|
  !                     | / / / / / / / / / / / |
  !                     |/ / / / / / / / / / / /|
  !          ----+------|--+-----------+--------|-+---
  !             x(1)      x(2)        x(3)       x(4)
  !                     a                       b
  !
  !   The result is stored in c.
  !   If a<x(1) or b>x(n+1) then the function is extended to a, resp. b, 
  !   assuming a constant value of y(1), resp. y(n).
  !   The index ilast should be set to 1 in the first call. 
  !   If subsequent integrals on the same vector are to be computed, 
  !   it should be kept unaltered between calls. This will speed up the 
  !   lookup of values if the various integrals are ascending in x.
  !
  !--------------------------------------------------------------------
  

  subroutine IntervalQuad_Const( x, y, a,b, c, ilast, status )

    use GO       , only : gol, goErr
    use num_tools, only : interval

    ! --- in/out ---------------------

    real, intent(in)       ::  x(:), y(:)
    real, intent(in)       ::  a, b
    real, intent(out)      ::  c
    integer, intent(inout) ::  ilast
    integer, intent(out)   ::  status

    ! --- const ---------------------

    character(len=*), parameter  ::  rname = mname//'/IntervalQuad_Const'

    ! --- local -------------------

    integer         ::  n
    real            ::  ya, yb
    integer         ::  i, i_a, i_b, iflag

    ! --- begin ---------------------

    ! check array sizes:
    n = size(y)
    if ( size(x) /= n+1 ) then
      write (gol,'("x should have one element more than y:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(y)  : ",i6)') size(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! check interval:    
    if ( b < a ) then
      write (gol,'("found strange interval [a,b] :")'); call goErr
      write (gol,'("  a  : ",es12.4)') a; call goErr
      write (gol,'("  b  : ",es12.4)') b; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if
    
    ! fill c with contribution of interval including a:
    i_a = ilast
    call Interval( x, a, i_a, iflag )
    select case ( iflag )
      case ( -1 )
        ! a < x(1)  ;  extend y to the left:
        c = y(1) * ( x(1) - a )
      case ( 0 )
        ! x(i_a) < a < x(i_a+1)
        c = y(i_a) * ( x(i_a+1) - a )
      case ( 1 )
        ! a > x(n+1)  ;  extend y to the right
        ! negative contribution to integral:
        c = - y(n) * ( a - x(n+1) )
      case default
        write (gol,'("unsupported iflag from call to Interval : ",i6)') iflag
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! add contributions of interval including b:
    i_b = i_a
    call Interval( x, b, i_b, iflag )
    select case ( iflag )
      case ( -1 )
        ! b < x(1)  ;  negative contribution
        c = c - y(1) * ( x(1) - b )
      case ( 0 )
        ! x(i_b) < b < x(i_b+1)
        c = c + y(i_b) * ( b - x(i_b) )
      case ( 1 )
        ! b > x(n+1)
        c = c + y(n) * ( b - x(n+1) )
      case default
        write (gol,'("unsupported iflag from call to Interval : ",i6)') iflag
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! add contributions of intermediate intervals
    do i = i_a+1, i_b-1
      c = c + y(i) * (x(i+1)-x(i))
    end do

    ! set index of last interval (including b):
    ilast = i_b
    
    ! ok
    status = 0

  end subroutine IntervalQuad_Const
  
  
  
  !-------------------------------------------------------------------
  ! 
  ! NAME
  !   IntervalSum  -  add contributions of small intervals
  !
  ! INTERFACE
  !   subroutine IntervalSum( x, y, a, b, c, ilast, fac )
  !     real, intent(in)              ::  x(:), y(:)
  !     real, intent(in)              ::  a, b
  !     real, intent(out)             ::  c
  !     integer, intent(inout)        ::  ilast
  !     real, intent(out), optional   ::  fac(:)
  !
  ! DESCRIPTION
  !   Compute sum of values y over invterval [a,b],
  !   for function y specified by:
  !     o interval boundaries x(1),..,x(n),x(n+1)
  !     o function values y(1),..,y(n) ; constant in interval
  !
  !                  y(1)       y(2)        y(3)
  !
  !                            100%
  !                        +-----o-----+    60%
  !                  30%   |     ^     +-----o----+
  !              +----o----+     |              ^
  !                     ^        |              |
  !                     |        |              |
  !                     |        |              |
  !          ----+------|--+-----------+--------|-+---
  !             x(1)      x(2)        x(3)       x(4)
  !                     a                       b
  !
  !   Each y(i) contributes to the sum for the fraction of [x(i),x(i+1)]
  !   covered by [a,b].
  !   The result is stored in c.
  !   If a<x(1) or b>x(n+1) an error is issued.
  !
  !   The index ilast should be set to 1 in the first call. 
  !   If subsequent integrals on the same vector are to be computed, 
  !   it should be kept unaltered between calls. This will speed up the 
  !   lookup of values if the various integrals are ascending in x.
  !
  !   If b < a, the result is negative.
  !
  !   If present, the array 'fac' is filt with the factors applied to y(:).
  !
  !--------------------------------------------------------------------
  

  subroutine IntervalSum( x, y, a_in, b_in, c, ilast, status, fac )

    use GO       , only : gol, goErr
    use num_tools, only : interval

    ! --- in/out ---------------------

    real, intent(in)              ::  x(:), y(:)
    real, intent(in)              ::  a_in, b_in
    real, intent(out)             ::  c
    integer, intent(inout)        ::  ilast
    integer, intent(inout)        ::  status
    real, intent(out), optional   ::  fac(:)

    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/IntervalSum'
    
    ! --- local -------------------

    integer         ::  n
    real            ::  ya, yb
    real            ::  f
    integer         ::  i, i_a, i_b, iflag
    
    real            ::  a, b, plusmin

    ! --- begin ---------------------
    
    ! check array sizes:
    n = size(y)
    if ( size(x) /= n+1 ) then
      write (gol,'("x should have one element more than y:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(y)  : ",i6)') size(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! check size of factor array if present:
    if ( present(fac) ) then
      if ( size(fac) /= n ) then
        write (gol,'("ERROR in ",a)') rname; status=1; return
        write (gol,'("fac should have same size as y:")'); call goErr
        write (gol,'("  size(y)   : ",i6)') size(y); call goErr
        write (gol,'("  size(fac) : ",i6)') size(fac); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
    end if
    
    ! increasing or decreasing interval ?
    if ( a_in <= b_in ) then
      a = a_in
      b = b_in
      plusmin = 1.0
    else
      a = b_in
      b = a_in
      plusmin = -1.0
    end if

    ! init output factors:    
    if ( present(fac) ) fac = 0.0

    ! fill c with contribution of interval including a:
    i_a = ilast
    call interval(x,a,i_a,iflag)
    select case ( iflag )
      case ( -1 )
        ! a < x(1)
        write (gol,'("interval is partly less than x :")'); call goErr
        write (gol,'("  a       : ",es12.4)')    a; call goErr
        write (gol,'("  x(1)    : ",es12.4)') x(1); call goErr
        write (gol,'("ERROR in ",a)') rname; status=1; return
      case ( 0 )
        ! x(i_a) < a < x(i_a+1)
        f = ( x(i_a+1) - a ) / ( x(i_a+1) - x(i_a) )
        c = y(i_a) * f
        if ( present(fac) ) fac(i_a) = fac(i_a) + f
      case ( 1 )
        ! a > x(n+1)
        write (gol,'("interval partly exceeds x :")'); call goErr
        write (gol,'("  a       : ",es12.4)')      a; call goErr
        write (gol,'("  x(n+1)  : ",es12.4)') x(n+1); call goErr
        write (gol,'("ERROR in ",a)') rname; status=1; return
      case default
        write (gol,'("unsupported iflag from call to Interval : ",i6)') iflag
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! add contributions of interval including b:
    i_b = i_a
    call interval( x, b, i_b, iflag )
    select case ( iflag )
      case ( -1 )
        ! b < x(1)
        write (gol,'("interval is outside x :")'); call goErr
        write (gol,'("  b       : ",es12.4)')    a; call goErr
        write (gol,'("  x(1)    : ",es12.4)') x(1); call goErr
        write (gol,'("ERROR in ",a)') rname; status=1; return
      case ( 0 )
        ! x(i_b) < b < x(i_b+1)
        if ( i_b > i_a ) then
          ! b in other interval; add contrib [x(i_b),b]
          f = ( b - x(i_b) ) / ( x(i_b+1) - x(i_b) )
        else
          ! a and b in same interval; substract contrib [b,x(i_b+1)]
          f = - ( x(i_b+1) - b ) / ( x(i_b+1) - x(i_b) )
        end if
        c = c + y(i_b) * f
        if ( present(fac) ) fac(i_b) = fac(i_b) + f
      case ( 1 )
        ! b > x(n+1)
        write (gol,'("interval exceeds x :")'); call goErr
        write (gol,'("  b       : ",es12.4)')      b; call goErr
        write (gol,'("  x(n+1)  : ",es12.4)') x(n+1); call goErr
        write (gol,'("ERROR in ",a)') rname; status=1; return
      case default
        write (gol,'("unsupported iflag from call to Interval : ",i6)') iflag; call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! add contributions of intermediate intervals
    do i = i_a+1, i_b-1
      c = c + y(i)
      if ( present(fac) ) fac(i) = fac(i) + 1.0
    end do

    ! set index of last interval (including b):
    ilast = i_b
    
    ! apply factor for b<a
    c = plusmin * c
    
    ! ok
    status = 0

  end subroutine IntervalSum
  
  
end module num_quad
