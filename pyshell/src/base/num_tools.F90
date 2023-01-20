module num_tools

  implicit none
  
  ! --- in/out -------------------
  
  private
  
  public   ::  Interval
  public   ::  Swap
  

  
contains


  !------------------------------------------------------------------
  !
  ! Computes ileft = max( i, 1 <= i <= n, .and. x(i) < a )
  !
  ! Input:
  !       x:      a real sequence of length n, assumed nondecreasing.
  !       a:      the point whose location with tespect to the sequence
  !               x is to be determined.
  !       ileft:  a 'first guess' estimate of the index to be found.
  !               this variable should be set to 1 for the first call to
  !               INTERVAL.
  !               Interpretate 'ileft' as the interval number {1,..,n-1}
  !
  ! Output:
  !       ileft   iflag     meaning
  !       ------  --------  -------------------------------------
  !        1       -1        a  <  x(1)
  !        i        0              x(i)  <=  a  <  x(i+1)
  !        n-1      0                        a  =  x( n )
  !        n        1                              x( n )  <  a
  !
  ! Method:
  !  Same as de Boor's INTERV sr, but not using a local index variable to
  !  store the last index found. Instead the 'first guess' index LEFT is
  !  used to initiate a search for the 'true' LEFT.
  !
  !------------------------------------------------------------------

  subroutine Interval( x, a, ileft, iflag )

    ! --- in/out ---------------------------------

    real, intent(in)        ::  x(:)
    real, intent(in)        ::  a
    integer, intent(inout)  ::  ileft
    integer, intent(out)    ::  iflag

    ! --- local ----------------------------------

    integer     ::  iright, istep, middle
    integer     ::  n

    ! --- begin -----------------------------------

    n = size(x)

    ! check correctness of first guess:
    if ( (ileft < 1) .or. (ileft>n) ) ileft = 1

    iright = ileft + 1
    if ( iright >= n ) then
      if ( a == x(n) ) then
        iflag = 0
        ileft = n-1
        return
      else if ( a > x(n) ) then
        iflag = 1
        ileft = n
        return
      end if
      if ( n <= 1 ) then
        iflag = -1
        ileft = 1
        return
      end if
      ileft = n-1
      iright = n
    end if

    if ( a >= x(iright) ) then

      ! now a >= x(iright). increase iright to capture a.
      istep = 1
      do
        ileft = iright
        iright = ileft + istep
        if ( iright >= n ) then
          if ( a == x(n) ) then
            iflag = 0
            ileft = n-1
            return
          else if ( a > x(n) ) then
            iflag = 1
            ileft = n
            return
          end if
          iright = n
          exit
        end if
        if ( a < x(iright) ) exit
        istep = istep*2
      end do

    else if ( a >= x(ileft) ) then

      iflag = 0
      return

    else

      ! now a  <  x(ileft).  decrease ileft to capture a

      istep = 1
      do
        iright = ileft
        ileft = iright - istep
        if ( ileft <= 1 ) then
          ileft = 1
          if ( a < x(1) ) then
            iflag = -1
            ileft = 1
            return
          end if
          exit
        end if
        if ( a >= x(ileft) ) exit
        istep = istep*2
      end do

    end if

    ! now x(ileft) <= x(iright). narrow the interval.

    do
      middle = ( ileft + iright ) / 2
      if ( middle == ileft ) then
        iflag = 0
        return
      end if
      if ( a < x(middle) ) then
        iright = middle
      else
        ileft = middle
      end if
    end do

  end subroutine Interval


  ! ==========================================================
  
  !
  ! call Swap( lx, xx )
  !
  !   Swaps the elements in array xx with length lx .
  !  

  subroutine Swap( lx, xx )

    ! --- in/out -------------------------------
    
    integer, intent(in)       ::  lx
    real, intent(inout)       ::  xx(lx)
    
    ! --- local --------------------------------
    
    real        ::  swp
    integer     ::  l
    
    ! --- begin --------------------------------

    do l = 1, lx/2
      swp = xx(l)
      xx(l) = xx(lx-l+1)
      xx(lx-l+1) = swp
    end do

  end subroutine Swap


end module num_tools
