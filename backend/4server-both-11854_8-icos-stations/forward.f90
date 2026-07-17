!======================================================================
!
!   Name:               forward.f90
!
!======================================================================
program forward
  ! Deklaration der Variablen
  implicit none
  integer ( kind = 4 )    :: n     ! dim of control space
  integer ( kind = 4 )    :: m     ! dim of obs space
  integer ( kind = 4 )    :: k     ! dim of target space
  real ( kind = 8 ), allocatable :: x0(:) ! prior control vector and uncertainty
  real ( kind = 8 ), allocatable :: y(:) ! concentration simulated from initial
  real ( kind = 8 ), allocatable :: yinit(:) ! concentration, simulated from initial
  real ( kind = 8 ), allocatable :: oj(:,:) ! sensitivity of obs to control vector (obs jacobian)
  ! set dimensions of control, obs, and target spaces
  call setdim2 (n, m, k )
  print*, 'n, m, k = ', n, m, k
  print*, 'Jac size [GB] = ', n*m*8/1.e9
  print*, 'Jac size [kB] = ', n*m*8/1.e3
  ! allocate fields
  allocate (x0(n),y(m),yinit(m),oj(m,n))
  ! set values of fields
  call setvalf (n, m, k, x0, yinit, oj)
  ! Forward run
  y = yinit + matmul (oj,x0)
  call wconc(m,y)
  call wfmat(6,'dummy.dat',1,m,y,'y simulated') ! this is just for development
  call wfmat(10,'output/c.dat',1,m,y,'y simulated') ! this is just for development
  
end program forward
