!
! Sorting.
!
!
! call Sort( x, xsort, isort, status )
!
!   Sorts the values of x(:) in increasing order.
!   Sorted x is stored in xsort,
!   indices in isort such that x(isort) = xsort .
!


module num_sort

  implicit none
  
  
  ! --- in/out ---------------------------------
  
  private
  
  public  ::  Sort
  
  
  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'num_sort'
  
  
  ! --- interfaces ------------------------------
  
  interface Sort
    module procedure sort_r_1d
    module procedure sort_r_2d
  end interface
  
  
contains

  
  ! ===================================================
  
  !
  ! call Sort( x, xsort, isort, status )
  !
  !   Sorts the values of x(:) in increasing order.
  !   Sorted x is stored in xsort,
  !   indices in isort such that:
  !      xsort(i) = x(isort(i))
  !
  
  
  subroutine sort_r_1d( x, xsort, isort, status )
  
    use GO       , only : gol, goErr
    use Num_Tools, only : Interval
  
    ! --- in/out --------------------------------
    
    real, intent(in)      ::  x(:)
    real, intent(out)     ::  xsort(:)
    integer, intent(out)  ::  isort(:)
    integer, intent(out)  ::  status
    
    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/sort_r_1d'

    ! --- local ---------------------------------
    
    integer  ::  ns
    real     ::  a
    integer  ::  i
    integer  ::  ileft
    integer  ::  iflag
  
    ! --- begin ---------------------------------
    
    ! check input
    if ( any(shape(xsort) /= shape(x)) .or. any(shape(isort) /= shape(x)) ) then
      write (gol,'("shape of output arrays should be same as input:")'); call goErr
      write (gol,'("  shape(x)      : ",i6)') shape(x); call goErr
      write (gol,'("  shape(xsort)  : ",i6)') shape(xsort); call goErr
      write (gol,'("  shape(isort)  : ",i6)') shape(isort); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if
    
    ! init
    xsort(1) = x(1)    ! sorted input
    isort(1) = 1       ! sorted indices
    ns       = 1       ! number of yet sorted values

    ! no interval searched yet
    ileft = 1
    
    ! loop over other elements
    do i = 2, size(x)

      ! current value
      a = x(i)

      ! where in sorted array ?
      call Interval( xsort(1:ns), a, ileft, iflag )

      if ( iflag < 0 ) then
        isort(1:ns+1) = (/ i, isort(1:ns) /)
        xsort(1:ns+1) = (/ a, xsort(1:ns) /)
      else if ( iflag > 0 ) then
        isort(1:ns+1) = (/ isort(1:ns), i /)
        xsort(1:ns+1) = (/ xsort(1:ns), a /)
      else
        isort(1:ns+1) = (/ isort(1:ileft), i, isort(ileft+1:ns) /)
        xsort(1:ns+1) = (/ xsort(1:ileft), a, xsort(ileft+1:ns) /)
      end if

      ns = ns + 1

    end do

    ! ok
    status = 0
    
  end subroutine sort_r_1d


  !
  ! call Sort( x, xsort, isort, status )
  !
  !   Sorts the values of x(:) in increasing order.
  !   Sorted x is stored in xsort,
  !   indices in isort such that:
  !      xsort(i,j) = x(isort(i,j,1),isort(i,j,2))
  !
  
  
  subroutine sort_r_2d( x, xsort, isort, status )
  
    use GO       , only : gol, goErr
    use Num_Tools, only : Interval
  
    ! --- in/out --------------------------------
    
    real, intent(in)      ::  x(:,:)
    real, intent(out)     ::  xsort(:,:)
    integer, intent(out)  ::  isort(:,:,:)
    integer, intent(out)  ::  status
    
    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/sort_r_2d'

    ! --- local ---------------------------------
    
    real     ::  x_1d(size(x))
    real     ::  xsort_1d(size(x))
    integer  ::  isort_1d(size(x))
    
    integer  ::  ns
    real     ::  a
    integer  ::  i
    integer  ::  ileft
    integer  ::  iflag

    ! --- begin ---------------------------------
    
    ! check xsort
    if ( any(shape(xsort) /= shape(x)) ) then
      write (gol,'("shape of sorted output should be same as input:")'); call goErr
      write (gol,'("  shape(x)      : ",2i6)') shape(x); call goErr
      write (gol,'("  shape(xsort)  : ",2i6)') shape(xsort); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if
    
    ! check isort
    if ( any(shape(isort) /= (/size(x,1),size(x,2),2/)) ) then
      write (gol,'("shape of isort be (/x1,x2,2/) :")'); call goErr
      write (gol,'("  shape(x)      : ",2i6)') shape(x); call goErr
      write (gol,'("  shape(isort)  : ",3i6)') shape(isort); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if
    
    ! store in 1d arrays:
    x_1d     = reshape(x    ,(/size(x)/))
    xsort_1d = reshape(xsort,(/size(x)/))
    isort_1d = reshape(isort,(/size(x)/))
    
    ! sort 1d array:
    call Sort( x_1d, xsort_1d, isort_1d, status )
    if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
    
    ! fill sorted x :
    xsort = reshape(xsort_1d,shape(x))
    
    ! fill indices
    isort(:,:,1) =  modulo(reshape(isort_1d,shape(x))-1,size(x,1))+1
    isort(:,:,2) = ceiling(reshape(isort_1d,shape(x))/real(size(x,1)))

    ! ok
    status = 0
    
  end subroutine sort_r_2d


end module num_sort


! ######################################################################


!program test
!
!  use num_sort, only : Sort
!  
!  implicit none
!
!  integer, parameter  ::  n1=10, n2=5
!  real      ::  x(n1,n2)
!  real      ::  xsort(n1,n2)
!  integer   ::  isort(n1,n2,2)
!  integer   ::  status
!
!  integer   ::  i, j
!  
!  call Random_Seed()
!  call Random_Number( x )
!  
!  call Sort( x, xsort, isort, status )
!  if (status/=0) stop 'error from sort'
!    
!  do j=1, n2
!    do i = 1, n1
!      write (*,'(2i4,f8.4,2i4,2f8.4)') i, j, x(i,j), isort(i,j,:), xsort(i,j), x(isort(i,j,1),isort(i,j,2))
!    end do
!  end do
!
!end program test
