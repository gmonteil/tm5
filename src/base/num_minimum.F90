!
!  type(TMinimumInfo)   ::  mini
!  real                 ::  x, y
!
!  call Init( mini, x0, dx, maxiter=10, atol=0.01, verbose=.true. )
!  do
!    call InitIter( mini, x )
!
!    ...
!    y = ... func( x )
!    ...
!
!    call DoneInter( mini, y, status )
!    if ( status ==  0 ) exit    ! ok, converged
!    if ( status == -1 ) cycle   ! try next iteration
!    if ( status ==  1 ) stop 'exceeded maximum number of iterations'
!    if ( status ==  2 ) stop 'no clear minimum'
!    stop 'error in iteration'
!  end do
!  call Done( mini )
!   
!
module num_minimum

  implicit none
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  TMinimumInfo
  public  ::  Init, Done
  public  ::  InitIter, DoneIter
  
  
  ! --- types ----------------------------
  
  type TMinimumInfo
    integer          ::  iter
    real             ::  dx
    real             ::  x1, x2, x3
    real             ::  y1, y2, y3
    real, pointer    ::  x(:)
    real, pointer    ::  y(:)
    !
    integer          ::  maxiter
    real             ::  atol
    !
    logical          ::  verbose
  end type TMinimumInfo
  
  
  ! --- interfaces --------------------------
  
  interface Init
    module procedure mini_Init
  end interface
  
  interface Done
    module procedure mini_Done
  end interface
  
  interface InitIter
    module procedure mini_InitIter
  end interface
  
  interface DoneIter
    module procedure mini_DoneIter
  end interface
  
    

contains


  ! ===========================================================

  
  subroutine mini_Init( mini, x0, dx, maxiter, atol, verbose )
  
    ! --- in/out ----------------------------
    
    type(TMinimumInfo), intent(out)    ::  mini
    real, intent(in)                   ::  x0, dx
    integer, intent(in)                ::  maxiter
    real, intent(in)                   ::  atol
    logical, intent(in)                ::  verbose
    
    ! --- begin -------------------------------
    
    mini%maxiter = maxiter
    mini%atol = atol
    
    mini%iter = 0

    mini%dx  = dx
    
    allocate( mini%x(mini%maxiter) )
    mini%x    = 0.0
    mini%x(1) = x0
    mini%x(2) = x0 + dx
    
    allocate( mini%y(mini%maxiter) )
    mini%y = 0.0
    
    mini%verbose = verbose

  end subroutine mini_Init


  ! ***


  subroutine mini_Done( mini )
  
    ! --- in/out ----------------------------
    
    type(TMinimumInfo), intent(inout)    ::  mini
    
    ! --- begin -------------------------------
    
    deallocate( mini%x )
    deallocate( mini%y )

  end subroutine mini_Done
  
  
  ! ===================================================


  subroutine mini_InitIter( mini, x )
  
    ! --- in/out ----------------------------
    
    type(TMinimumInfo), intent(inout)  ::  mini
    real, intent(out)                  ::  x
    
    ! --- local -----------------------------
    
    integer         ::  i1, i2, i3
    
    ! --- begin -------------------------------
    
    ! next interation step:
    mini%iter = mini%iter + 1

    ! leave if this is initialisation:
    if ( mini%iter <= 2 ) then
      x = mini%x(mini%iter)
      return
    end if
    
    ! previous and prev-previous:
    i1 = mini%iter-2
    i2 = mini%iter-1
    i3 = mini%iter
    
    mini%x1 = mini%x(i1)
    mini%x2 = mini%x(i2)
    
    mini%y1 = mini%y(i1)
    mini%y2 = mini%y(i2)
    
    do

      if ( mini%y2 <= mini%y1 ) then

        !
        !         o
        !             o
        !                 ?
        !  -------+---+---+------- x
        !         1   2   3
        !

        ! take same step ...

      else if ( mini%y2 > mini%y1 ) then

        !
        !             o
        !     o        
        !         ?        
        !  ---+---+---+----------- x
        !     1   3   2    
        !

        ! reverse search direction, decrase step:
        mini%dx   = - mini%dx/2.0

      end if

      ! next point:
      mini%x3 = mini%x2 + mini%dx
      
      ! already processed ? take x and y and try again:
      if ( any( mini%x(1:i2) == mini%x3 ) ) then
        do i3 = 1, i2
          if ( mini%x(i3) == mini%x3 ) then
            mini%y3 = mini%y(i3)
            exit
          end if
        end do
        mini%x1 = mini%x2  ; mini%y1 = mini%y2
        mini%x2 = mini%x3  ; mini%y2 = mini%y3
      else
        exit
      end if
      
    end do

    ! store next point:
    mini%x(mini%iter) = mini%x3

    ! set output:
    x = mini%x3

  end subroutine mini_InitIter


  ! ***

  !
  ! status   description
  !  -1       try next iter
  !   0       converged
  !   1       reached max iter
  !   2       error; no minium in interval ?
  !  

  subroutine mini_DoneIter( mini, y, status )
  
    ! --- in/out ----------------------------
    
    type(TMinimumInfo), intent(inout)    ::  mini
    real, intent(in)                     ::  y
    integer, intent(out)                 ::  status
    
    ! --- local -----------------------------
    
    integer         ::  i1, i2, i3
    
    ! --- begin -------------------------------

    ! previous and prev-previous:
    i1 = mini%iter-2
    i2 = mini%iter-1
    i3 = mini%iter

    ! store evaluation:
    mini%y(i3) = y
    
    if ( mini%verbose ) then
      print *, '        iter:', mini%iter, mini%x(i3), mini%y(i3)
      !
      !print *, 'iter', mini%iter
      !print *, '  x=', mini%x
      !print *, '  y=', mini%y
      !print *, ' '
    end if
    
    if ( mini%iter <= 2 ) then
      ! next iter
      status = -1
      return
    end if

!print *, '          --->',    mini%x(i3), mini%x(i2), abs(mini%x(i3)-mini%x(i2)), mini%atol
    if ( abs(mini%x(i3)-mini%x(i2)) <= mini%atol ) then
      !if ( mini%verbose ) then
      !  print *, 'iter: minimal difference :', mdif
      !  print *, '                atol * y :', mini%atol * abs(mini%y(2))
      !end if
      status = 0
      return
    end if

    if ( mini%iter == mini%maxiter ) then
      if ( mini%verbose ) then
        print *, 'iter: reached maximum number of iterations'
      end if
      status = 1
      return
    end if
      
    ! by default next iter:
    status = -1
    
  end subroutine mini_DoneIter
  
  


end module num_minimum
