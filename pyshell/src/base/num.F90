!
!ProTeX: 1.14-AJS
!
!BOI
!
! !TITLE:        Num  -  numerical interpolation, quadrature, etc
! !AUTHORS:      Arjo Segers
! !AFFILIATION:  KNMI
! !DATE:         \today
!
! !INTRODUCTION: Tools
!
! \begin{itemize}
!
! \item Search interval that contains a value:
!
!   \bv
!     subroutine Interval( x, a, ileft, iflag )
!       real, intent(in)        ::  x(:)
!       real, intent(in)        ::  a
!       integer, intent(inout)  ::  ileft
!       integer, intent(out)    ::  iflag
!     end subroutine Interval
!
!   Computes:
!     ileft = max( i, 1 <= i <= n, and x(i) < a )
!
!   Input:
!     x:      a real sequence of length n, assumed nondecreasing.
!     a:      the point whose location with tespect to the sequence
!             x is to be determined.
!     ileft:  a 'first guess' estimate of the index to be found.
!             this variable should be set to 1 for the first call to
!             INTERVAL.
!             Interpretate 'ileft' as the interval number {1,..,n-1}
!
!   Output:
!     ileft   iflag     meaning
!     ------  --------  -------------------------------------
!      1       -1        a  <  x(1)
!      i        0              x(i)  <=  a  <  x(i+1)
!      n-1      0                        a  =  x( n )
!      n        1                              x( n )  <  a
!
!   Method:
!     Same as de Boor's INTERV sr, but not using a local index variable to
!     store the last index found. Instead the 'first guess' index LEFT is
!     used to initiate a search for the 'true' LEFT.
!
!   \ev
!
! \item Swaps the elements in array 'xx' with length 'lx' :
!
!   \bv
!     subroutine Swap( lx, xx )
!       integer, intent(in)       ::  lx
!       real, intent(inout)       ::  xx(lx)
!     end subroutine Swap
!   \ev
!
! \end{itemize}
!
!
! !INTRODUCTION: Interpolation
!
! \begin{itemize}
!
! \item Linear interpolation, 1D arrays.
!   \bv
!     subroutine Interp_Lin( x,y, x0, y0, last, status )
!       real, intent(in)        ::  x(n), y(n)
!       real, intent(in)        ::  x0
!       real, intent(out)       ::  y0
!       integer, intent(inout)  ::  last
!       integer, intent(out)    ::  status
!     end subroutine Interp_Lin
!   \ev
!   Interpolate 'y' defined at positions 'x' linear to 'x0',
!   return result in 'y0'. Index 'last' is a help variable used
!   to speedup search algorithm in repated calls to this routine.
!
! \item Linear interpolation in 2D array 'y'; last dimension corresponds to 'x'.
!   \bv
!     subroutine Interp_Lin( x,y, x0, y0, last, status )
!       real, intent(in)        ::  x(n2), y(n1,n2)
!       real, intent(in)        ::  x0
!       real, intent(out)       ::  y0(n1)
!       integer, intent(inout)  ::  last
!       integer, intent(out)    ::  status
!     end subroutine Interp_Lin
!   \ev
!
! \item Linear interpolation in 3D array 'y'; last dimension corresponds to 'x'.
!   \bv
!     subroutine Interp_Lin( x,y, x0, y0, last, status )
!       real, intent(in)        ::  x(n3), y(n1,n2,n3)
!       real, intent(in)        ::  x0
!       real, intent(out)       ::  y0(n1,n2)
!       integer, intent(inout)  ::  last
!       integer, intent(out)    ::  status
!     end subroutine Interp_Lin
!   \ev
!
! \item Compute interpolation weights. 
!   \bv
!     subroutine Interp_Lin_Weight( x, x0, w, status )
!       real, intent(in)        ::  x(n)
!       real, intent(in)        ::  x0
!       real, intent(out)       ::  w(n)
!       integer, intent(out)    ::  status
!     end subroutine InterpolWeight
!   \ev
!   Let array 'y' be defined on axis points 'x',
!   and let 'x0' be a point on the axis. The weights 'w' are then
!   computed such that:
!   \bv
!     y interpolated to x0 = sum( y(:) * w(:) )
!   \ev
!
! \item Linear interpolation, 1D arrays, 'x' cyclic with period 'p' :
!   \bv
!     subroutine CircInterp_Lin( x, p, y, x0, y0, last, status )
!       real, intent(in)        ::  x(:), y(:)
!       real, intent(in)        ::  p
!       real, intent(in)        ::  x0
!       real, intent(out)       ::  y0
!       integer, intent(inout)  ::  last
!       integer, intent(out)    ::  status
!     end subroutine Interp_Lin
!   \ev
!
! \item Linear interpolation in 1D arrays to multiple points 'x0'.
!   \bv
!     subroutine Interp_Lin( x,y, x0, y0, status )
!       real, intent(in)        ::  x(n), y(n)
!       real, intent(in)        ::  x0(q)
!       real, intent(out)       ::  y0(q)
!       integer, intent(out)    ::  status
!     end subroutine Interp_Lin
!   \ev
!
! \item Hermite-cubic interpolation in 1D arrays to multiple points 'x0'.
!   \bv
!     subroutine Interp_MuHerm( x,y, x0, y0, status )
!       real, intent(in)        ::  x(n), y(n)
!       real, intent(in)        ::  x0(q)
!       real, intent(out)       ::  y0(q)
!       integer, intent(out)    ::  status
!     end subroutine Interp_MuHerm
!   \ev
!
! \end{itemize}
!
!
!
! !INTRODUCTION: Quadrature
!
! \begin{itemize}
!
!   \item Qudrature assuming linear interpolation.
!     \bv
!       subroutine IntervalQuad_Lin( x,y, a,b, c, ilast, status )
!         real, intent(in)       ::  x(n), y(n)
!         real, intent(in)       ::  a, b
!         real, intent(out)      ::  c
!         integer, intent(inout) ::  ilast
!         integer, intent(out)   ::  status
!       end subroutine IntervalQuad_Lin
!     \ev
!     Compute an approximation to the integral:
!     \bv
!        b
!       int  y(x) dx
!       x=a
!     \ev
!     Input is formed by two vectors 'x' and 'y' .
!     Vector 'y' contains the function 
!     values to the abscissa values contained in 'x'. 
!     The result, stored in 'c', is computed using the trapezoid formula 
!     throughout the domain of integration. 
!     If 'a < x(1)' or 'b > x(n)' then the function is extended to 'a', resp. 'b', 
!     assuming a constant value of 'y(1)', resp. 'y(n)'.
!
!     The index 'ilast' should be set to 1 in the first call. If subsequent
!     integrals on the same vector are to be computed, it should be kept
!     unaltered between calls. This will speed up the lookup of values
!     if the various integrals are ascending in 'x'.
!
!   \item Qudrature of cosined function assuming linear interpolation.
!     \bv
!       subroutine IntervalQuad_Cos_Lin( x, y, a, b, c, ilast, status )
!         real, intent(in)         ::  x(n), y(n)
!         real, intent(in)         ::  a, b
!         real, intent(out)        ::  c
!         integer, intent(inout)   ::  ilast
!         integer, intent(out)     ::  status
!     \ev
!     Compute an approximation to the integral:
!     \bv
!        b
!       int  y(x) cos(x) dx
!       x=a
!     \ev
!     See routine 'IntervalQuad_Lin' for method.
!
!   \item Qudrature of piecewise constant function.
!     \bv
!       subroutine IntervalQuad_Const( x, y, a,b, c, ilast, status )
!         real, intent(in)       ::  x(n+1), y(n)
!         real, intent(in)       ::  a, b
!         real, intent(out)      ::  c
!         integer, intent(inout) ::  ilast
!         integer, intent(out)   ::  status
!       end subroutine IntervalQuad_Const
!     \ev
!     Compute integral over '[a,b]' for function specified by:
!     \begin{itemize}
!       \item interval boundaries 'x(1),..,x(n),x(n+1)'
!       \item function values 'y(1),..,y(n)' ; constant in interval
!     \end{itemize}
!     \bv
!  
!                    y(1)       y(2)        y(3)
!                          +-----o-----+
!                          |/ / / / / /+-----o----+
!                +----o----+ / / / / / / / / / /|
!                       | / / / / / / / / / / / |
!                       |/ / / / / / / / / / / /|
!            ----+------|--+-----------+--------|-+---
!               x(1)      x(2)        x(3)       x(4)
!                       a                       b
!     \ev
!     See routine 'IntervalQuad_Lin' for method.
!
!   \item Sum values defined in interval.
!     \bv
!       subroutine IntervalSum( x, y, a_in, b_in, c, ilast, status, fac )
!         real, intent(in)              ::  x(n+1), y(n)
!         real, intent(in)              ::  a_in, b_in
!         real, intent(out)             ::  c
!         integer, intent(inout)        ::  ilast
!         integer, intent(inout)        ::  status
!         real, intent(out), optional   ::  fac(n)
!       end subroutine IntervalSum
!     \ev
!     Compute sum of values 'y' over interval '[a,b]',
!     \begin{itemize}
!       \item interval boundaries 'x(1),..,x(n),x(n+1)'
!       \item function values 'y(1),..,y(n)' ; constant in interval
!     \end{itemize}
!     \bv
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
!     \ev
!     Each 'y(i)' contributes to the sum for the fraction 
!     of '[x(i),x(i+1)]' covered by '[a,b]'.
!     If optional output array 'fac' is provided it is filled
!     with the contributions for each element of 'y'.
!     The result is stored in 'c'.
!     If 'a<x(1)' or 'b>x(n+1)' an error is issued.
!
! \end{itemize}
!
! !INTRODUCTION: Sorting
!
! \begin{itemize}
!
! \item Sort rank-1 array:
!
!   \bv
!     subroutine Sort( x, xsort, isort, status )
!       real, intent(in)      ::  x    (n)
!       real, intent(out)     ::  xsort(n)
!       integer, intent(out)  ::  isort(n)
!       integer, intent(out)  ::  status
!     end subroutine
!   \ev
!
!   Sort the values of an array 'x' in increasing order.
!   Sorted 'x' is stored in 'xsort', indices in 'isort' such that:
!   \bv
!     xsort(i) = x(isort(i))
!   \ev
!
! \item Sort rank-2 array:
!
!   \bv
!     subroutine Sort( x, xsort, isort, status )
!       real, intent(in)      ::  x    (n1,n2)
!       real, intent(out)     ::  xsort(n1,n2)
!       integer, intent(out)  ::  isort(n1,n2,2)
!       integer, intent(out)  ::  status
!     end subroutine Sort
!   \ev
!
!   Sort the values of an array 'x' in increasing order.
!   Sorted 'x' is stored in 'xsort', indices in 'isort' such that:
!   \bv
!     xsort(i,j) = x(isort(i,j,1),isort(i,j,2))
!   \ev
!
! \end{itemize}
!
! !INTRODUCTION: Random numbers
!
! \begin{itemize}
!
! \item Gaussian random deviate:
!
!   \bv
!     x = GasDev()
!   \ev
!
! \end{itemize}
!
!
!EOI
!

module num

  use num_tools
  use num_sort
  use num_minimum
  use num_interp
  use num_quad
  use num_matrix
  use num_random
  
  implicit none
  
  public

end module num
