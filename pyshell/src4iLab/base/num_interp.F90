module num_interp

  implicit none
  
  ! --- in/out -------------------
  
  private
  
  public   ::  Interp_Lin
  public   ::  Interp_Lin_Weight
  public   ::  CircInterp_Lin
  public   ::  Interp_MuHerm


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'num_interp'
  
  
  ! --- interfaces -----------------
  
  interface Interp_Lin
    module procedure linterp_1
    module procedure linterp_1_row
    module procedure linterp_2
    module procedure linterp_3
  end interface
  
  
  interface CircInterp_Lin
    module procedure linterp_1_circ
  end interface

  
contains


  ! =====================================================



  !----------------------------------------------------------------------
  ! NAME
  !   linterp - makes a linear interpolation.
  !
  ! SYNOPSIS
  !   call linterp( x,y, x0,y0, last )
  !
  ! ARGUMENTS
  !    x     (in)    the array of x-values
  !    y     (in)    the array of y-values
  !    x0    (in)    the location where interpolation is desired
  !    y0    (out)   The interpolated value
  !    last  (inout) A first guess of the index of the left hand
  !                  edge of the interval, in which x0 lies. On return
  !                  this variable has the exact index. 
  !
  ! DESCRIPTION
  !   This version uses linear interpolation only in the interior of
  !   the array x. If x0<x(1), or if x0>x(end) then y0 will be 
  !   x(1) or x(end) resp.
  !
  !   The modified Boor routine INTERVAL is used to locate x0 within 
  !   the array x.
  !
  !-----------------------------------------------------------------------

  subroutine linterp_1( x,y, x0, y0, last, status )

    use GO, only : gol, goErr
    use num_tools, only : interval

    ! --- in/out ------------------------

    real, intent(in)        ::  x(:), y(:)
    real, intent(in)        ::  x0
    real, intent(out)       ::  y0
    integer, intent(inout)  ::  last
    integer, intent(out)    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/linterp_1'

    ! --- local -------------------------

    integer    ::  mflag
    integer    ::  n

    ! --- begin -----------------------

    ! length of input arrays
    n = size(x)
    if ( size(y) /= n ) then
      write (gol,'("x and y should have same lengths:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(y)  : ",i6)') size(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! determine interval in x that contains x0 :
    call Interval( x, x0, last, mflag )

    if ( mflag == 0 ) then
      y0 = y(last) + (x0-x(last))*(y(last+1)-y(last))/(x(last+1)-x(last))
    else if ( mflag .lt. 0 ) then
      ! x0 is <x(1)
      y0 = y(1)
    else
      ! x0 is >x(n)
      y0 = y(n)
    endif
    
    ! ok
    status = 0

  end subroutine linterp_1


  ! ===

  ! similar, but now y is a rank 2 array;
  ! last dimension corresponds to 'x'

  subroutine linterp_2( x, y, x0, y0, last, status )

    use GO       , only : gol, goErr
    use num_tools, only : interval

    ! --- in/out ------------------------

    real, intent(in)        ::  x(:), y(:,:)
    real, intent(in)        ::  x0
    real, intent(out)       ::  y0(:)
    integer, intent(inout)  ::  last
    integer, intent(out)    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/linterp_2'

    ! --- local -------------------------

    integer    ::  mflag
    integer    ::  n

    ! --- begin -----------------------

    ! length of input arrays
    n = size(x)
    if ( size(y,2) /= n ) then
      write (gol,'("dim 2 of y should have same length as x :")'); call goErr
      write (gol,'("  size(x)   : ", i6)') size(x); call goErr
      write (gol,'("  shape(y)  : ",2i6)') shape(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
      stop
    end if
    if ( size(y0) /= size(y,1) ) then
      write (gol,'("output y0 should have same number of elements as first dim of y:")'); call goErr
      write (gol,'("  size(y0)  : ", i6)') size(y0); call goErr
      write (gol,'("  shape(y)  : ",2i6)') shape(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
      stop
    end if

    ! determine position of x0 in array x :
    call Interval( x, x0, last, mflag )

    ! inter- or extra-polation based on position of x0:
    if ( mflag == 0 ) then
      y0 = y(:,last) + (x0-x(last))*(y(:,last+1)-y(:,last))/(x(last+1)-x(last))
    else if ( mflag .lt. 0 ) then
      ! x0 is <x(1)
      y0 = y(:,1)
    else
      ! x0 is >x(n)
      y0 = y(:,n)
    endif
    
    ! ok
    status = 0

  end subroutine linterp_2


  ! ===

  ! similar, but now y is a rank 3 array;
  ! last dimension corresponds to 'x'

  subroutine linterp_3( x, y, x0, y0, last, status )

    use GO       , only : gol, goErr
    use num_tools, only : interval

    ! --- in/out ------------------------

    real, intent(in)        ::  x(:), y(:,:,:)
    real, intent(in)        ::  x0
    real, intent(out)       ::  y0(size(y,1),size(y,2))
    integer, intent(inout)  ::  last
    integer, intent(out)    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/linterp_3'

    ! --- local -------------------------

    integer    ::  n
    integer    ::  mflag

    ! --- begin -----------------------

    ! length of input arrays
    n = size(x)
    if ( size(y,3) /= n ) then
      write (gol,'("dim 3 of y should have same length as x :")'); call goErr
      write (gol,'("  size(x)   : ", i6)') size(x); call goErr
      write (gol,'("  shape(y)  : ",3i6)') shape(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
      stop
    end if
    if ( (size(y0,1) /= size(y,1)) .or. (size(y0,2) /= size(y,2)) ) then
      write (gol,'("output y0 should have same shape as first dims of y:")'); call goErr
      write (gol,'("  shape(y0) : ",2i6)') shape(y0); call goErr
      write (gol,'("  shape(y)  : ",3i6)') shape(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
      stop
    end if

    ! determine position of x0 in array x :
    call Interval( x, x0, last, mflag )

    ! inter- or extra-polation based on position of x0:
    if ( mflag == 0 ) then
      ! x0 in [x(last),x(last+1)]
      y0 = y(:,:,last) + (x0-x(last))*(y(:,:,last+1)-y(:,:,last))/(x(last+1)-x(last))
    else if ( mflag .lt. 0 ) then
      ! x0 is <x(1)
      y0 = y(:,:,1)
    else
      ! x0 is >x(n)
      y0 = y(:,:,n)
    endif

    ! ok
    status = 0

  end subroutine linterp_3


  ! ===
  
  !
  ! weights such that:
  !  
  !  y(:) defined on axis x(:)
  !  x0 on axis x(:) 
  !  y interpolated to x0 = sum( y(:) * w(:) )
  !

  subroutine Interp_Lin_Weight( x, x0, w, status )

    use GO       , only : gol, goErr
    use num_tools, only : interval

    ! --- in/out ------------------------

    real, intent(in)        ::  x(:)
    real, intent(in)        ::  x0
    real, intent(out)       ::  w(:)
    integer, intent(out)    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/Interp_Lin_Weight'

    ! --- local -------------------------

    integer    ::  last
    integer    ::  mflag
    integer    ::  n

    ! --- begin -----------------------

    ! length of input arrays
    n = size(x)
    if ( size(w) /= n ) then
      write (gol,'("arrays x and w should have same length:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(w)  : ",i6)') size(w); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if
    
    ! initial zero
    w = 0.0

    last = 1
    call Interval( x, x0, last, mflag )

    if ( mflag == 0 ) then
      ! y0 = y(last) + (x0-x(last))*(y(last+1)-y(last))/(x(last+1)-x(last))
      !    = [1-(x0-x(last))/(x(last+1)-x(last))]*y(last) + [(x0-x(last))/(x(last+1)-x(last))]*y(last+1)
      w(last+1) = (x0-x(last))/(x(last+1)-x(last))
      w(last  ) = 1.0 - w(last+1)
    else if ( mflag .lt. 0 ) then
      ! x0 is <x(1)
      !y0 = y(1)
      w(1) = 1.0
    else
      ! x0 is >x(n)
      !y0 = y(n)
      w(n) = 1.0
    endif
    
    ! ok
    status = 0

  end subroutine Interp_Lin_Weight


  ! ===

  ! idem with x periodic with periode p
  !
  ! 25-04-2002
  !   Debugged: division by zero if x(n)==x(1)
  !

  subroutine linterp_1_circ( x, p, y, x0, y0, last, status )

    use GO       , only : gol, goErr
    use num_tools, only : Interval

    ! --- in/out ------------------------

    real, intent(in)        ::  x(:), y(:)
    real, intent(in)        ::  p
    real, intent(in)        ::  x0
    real, intent(out)       ::  y0
    integer, intent(inout)  ::  last
    integer, intent(out)    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/linterp_1_circ'

    ! --- local -------------------------

    integer    ::  mflag
    integer    ::  n, i
    real       ::  x0p

    ! --- begin -----------------------

    ! length of input arrays
    n = size(x)
    if ( size(y) /= n ) then
      write (gol,'("x and y should have same lengths:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(y)  : ",i6)') size(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! check position of x0
    if ( x0 >= x(1) .and. x0 <= x(n) ) then
      ! x0 in [x(1),x(n)]
      x0p = x0
    else
      ! set x0p in [x(1),x(1)+p)
      x0p = x(1) + modulo(x0-x(1),p)
    end if

    ! deterimine interval with x0p :
    call Interval( x, x0p, last, mflag )

    if ( mflag == 0 ) then
      ! default interpolation
      y0 = y(last) + (x0p-x(last))*(y(last+1)-y(last))/(x(last+1)-x(last))
    else if ( mflag < 0 ) then
      ! impossible that x0 < x(1) ...
      write (gol,'("internal error: x0 is smaller than x(1) :")'); call goErr
      write (gol,'(" x0   : ",e16.4)') x0; call goErr
      write (gol,'(" p    : ",e16.4)') p; call goErr
      write (gol,'(" x0p  : ",e16.4)') x0p; call goErr
      write (gol,'(" x    : ",e16.4)') x; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    else
      ! mflag > 0, thus x0p is between x(n) and first x(i) > x(n)
      i = 1
      do
        if ( x(n) < p+x(i) ) exit
        i = i + 1
      end do
      y0 = y(n) + (x0p-x(n))*(y(i)-y(n))/(modulo(x(i),p)-modulo(x(n),p))
    end if
    
    ! ok
    status = 0

  end subroutine linterp_1_circ


  ! ===


  subroutine linterp_1_row( x,y, x0, y0, status )

    use GO, only : gol, goErr

    ! --- in/out ------------------------

    real, intent(in)        ::  x(:), y(:)
    real, intent(in)        ::  x0(:)
    real, intent(out)       ::  y0(:)
    integer, intent(out)    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/linterp_1_row'

    ! --- local -------------------------

    integer    ::  n
    integer    ::  i
    integer    ::  last

    ! --- begin -----------------------

    ! check input shapes:
    n = size(x0)
    if ( size(y0) /= n ) then
      write (gol,'("x and y should have same lengths:")'); call goErr
      write (gol,'("  size(x)  : ",i6)') size(x); call goErr
      write (gol,'("  size(y)  : ",i6)') size(y); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! loop over interpolation points:
    last = 0
    do i = 1, n
      call Interp_lin( x, y, x0(i), y0(i), last, status )
      if (status/=0) then; write (gol,'("ERROR in ",a)') rname; status=1; return; end if
    end do
    
    ! ok
    status = 0

  end subroutine linterp_1_row


  ! =======================================================================


  !----------------------------------------------------------------------
  ! Hermite cubic interpolation from one vector to annother
  ! Input:        x       the array of x-values
  !               y       the array of y-values
  !               n       the no. of x and y values
  !               x0      the locations where interpolation is desired
  ! Output:       y0      the interpolated values
  ! Input:  m       no of values to interpolate
  ! 
  !       both arrays x(n) and x0(m) are assumed to be strictly monotonically
  !       increasing
  !
  !       If x0.lt.x(1), or if x0.gt.x(n) then y0 will be extrapolated
  ! linearly
  !
  !-----------------------------------------------------------------------

  ! original name: muherm

  subroutine Interp_MuHerm( x, y, x0, y0, status )

    ! --- in/out ----------------------------------

    real, intent(in)        ::  x(:)         
    real, intent(in)        ::  y(size(x))   
    real, intent(in)        ::  x0(:)        
    real, intent(out)       ::  y0(size(x0)) 
    integer, intent(out)    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/Interp_MuHerm'

    ! --- local -----------------------------------

    integer     ::  n, m
    integer     ::  last, mflag
    real        ::  fx(size(x))     

    real        ::  a1, a2, a3, a4    
    real        ::  g1, g2, g3 
    real        ::  hp, hm
    real        ::  xi
    integer     ::  i, k
    integer     ::  il

    ! --- begin -------------------------------------

    n = size(x)
    m = size(x0)

    ! calculate slopes first
    ! use 3pt formula at borders
    g1=(2*x(1)-x(2)-x(3))/(x(1)-x(2))/(x(1)-x(3))
    g2=(x(1)-x(3))/(x(2)-x(1))/(x(2)-x(3))
    g3=(x(1)-x(2))/(x(3)-x(1))/(x(3)-x(2))
    fx(1)=g1*y(1)+g2*y(2)+g3*y(3)

    g1=(x(n)-x(n-1))/(x(n-2)-x(n-1))/(x(n-2)-x(n))
    g2=(x(n)-x(n-2))/(x(n-1)-x(n-2))/(x(n-1)-x(n))
    g3=(2*x(n)-x(n-1)-x(n-2))/(x(n)-x(n-2))/(x(n)-x(n-1))
    fx(n)=g1*y(n-2)+g2*y(n-1)+g3*y(n)

    do i = 2, n-1
      hp=x(i+1)-x(i)
      hm=x(i)-x(i-1)
      g1=-hp/(hm*(hp+hm))
      g2=(hp-hm)/(hp*hm)  
      g3=hm/(hp*(hp+hm))
      fx(i)=y(i-1)*g1+y(i)*g2+y(i+1)*g3
    end do

    ! do now the interpolation

    ! index of left border of interval where x0(k) is located
    il=1

    do k = 1, m
      if(x0(k).lt.x(1)) then
        y0(k)=y(1)+fx(1)*(x0(k)-x(1))
      else if(x0(k).gt.x(n)) then
        y0(k)=y(n)+fx(n)*(x0(k)-x(n))
      else
        do
          if ( x0(k) .gt. x(il+1) ) then
            il=il+1
            cycle
          else
            xi=(x0(k)-x(il))/(x(il+1)-x(il))
            a1=f1(1-xi)
            a2=f1(xi)
            a3=-f2(1-xi)
            a4=f2(xi)
            y0(k)=a1*y(il)+a2*y(il+1)+(x(il+1)-x(il))*(a3*fx(il)+a4*fx(il+1))
            exit
          end if
        end do
      endif
    end do

    ! ok
    status = 0
    return

  contains

    real function  f1(u)
      real, intent(in)  ::  u
      f1 = u*u*(-2*u+3)
    end function f1

    real function f2(u)
      real, intent(in)  :: u
      f2 = u*u*(u-1)
    end function f2

  end subroutine Interp_MuHerm
      
      
end module num_interp
