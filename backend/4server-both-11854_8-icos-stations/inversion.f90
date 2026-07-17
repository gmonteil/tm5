!======================================================================
!
!   Name:               inversion.f90
!
!======================================================================
program inversion
  ! Deklaration der Variablen
  implicit none
  integer ( kind = 4 )    :: n     ! dim of control space
  integer ( kind = 4 )    :: m     ! dim of obs space
  integer ( kind = 4 )    :: k     ! dim of target space
  real ( kind = 8 ), allocatable :: x0(:) ! prior control vector and uncertainty
  real ( kind = 8 ), allocatable :: sx0(:) ! its uncertainty (1 sigma)
  real ( kind = 8 ), allocatable :: xp(:) ! perturbed prior, only needed for identical twin
  real ( kind = 8 ), allocatable :: x(:) ! posterior control vector
  real ( kind = 8 ), allocatable :: deltax(:) ! control vector increment (post minus prior)
  real ( kind = 8 ), allocatable :: y(:), sy(:) ! obs and their uncertainty
  real ( kind = 8 ), allocatable :: y0(:) ! obs equivalents simulated from prior
  real ( kind = 8 ), allocatable :: yinit(:) ! simulated from initial [ppb]
  real ( kind = 8 ), allocatable :: deltay(:) ! ypost-y0
  real ( kind = 8 ), allocatable :: oj(:,:) ! sensitivity of obs to control vector (obs jacobian)
  real ( kind = 8 ), allocatable :: t(:), t0(:) ! targets post and prior 
  real ( kind = 8 ), allocatable :: deltat(:) ! targets increment (post minus prior)
  real ( kind = 8 ), allocatable :: tj(:,:) ! Sensitivity of targets to control vector
  real ( kind = 8 ), allocatable :: ct0(:,:) ! Covariance of prior target uncertainty
  real ( kind = 8 ), allocatable :: ct(:,:) ! Covariance of posterior target uncertainty
  real ( kind = 8 ), allocatable :: w(:), u(:,:), v(:,:) ! sigular vectors and values
  real ( kind = 8 ), allocatable :: tju(:,:), tjum(:,:) ! aux fields
  logical :: matu = .true., matv = .true. 
  logical :: twin = .false.
  integer ( kind = 4) :: i, j, l, ierr
  real ( kind = 8 ) :: eps
  ! set dimensions of control, obs, and target spaces
  call setdim2 (n, m, k )
  print*, 'n, m, k = ', n, m, k
  print*, 'Jac size [GB] = ', n*m*8/1.e9
  print*, 'Jac size [kB] = ', n*m*8/1.e3
  ! allocate fields
  allocate (x0(n),sx0(n),x(n),xp(n),deltax(n))
  allocate (y(m),sy(m),y0(m),yinit(m),oj(m,n),deltay(m))
  allocate (t(k), t0(k), deltat(k), tj(k,n), ct0(k,k), ct(k,k))
  allocate (w(m), u(n,m), v(m,m), tju(k,m), tjum (k,m))
  ! set values of fields
  call setval2 (n, m, k, x0, sx0, xp, y, sy, yinit, oj, tj, twin)
  ! SVD  
  eps = sqrt(epsilon(eps))*10
  ierr = 0
  print*, 'starting SVD'
  call svd ( n, m, transpose(oj), w, matu, u, matv, v, ierr )
  tju = matmul(tj,u)
  do l = 1, k
     tjum(l,:)= tju(l,:) * w(:)/(1._8+w(:)**2)
  enddo
  y0 = matmul (oj,x0)
  deltat = matmul(tjum,matmul(transpose(v),y-y0))
  do l = 1, k
     tju(l,:)= tju(l,:) * sqrt(w(:)**2/(1._8+w(:)**2)) 
  enddo
  ct0 = matmul(tj,transpose(tj))
  ct = ct0 - matmul(tju,transpose(tju))
  t0 = matmul (tj, x0)
  ! output
  call wtarget(k,t0,deltat,ct0,ct)
  do i = 1, n
     u(i,:)= u(i,:) * w(:)/(1._8+w(:)**2)
  enddo
  deltax = matmul(u,matmul(transpose(v),y-y0))
  call wfmat(6,'dummy.dat',1,n,(deltax+x0)*sx0,'x post in MtCH4/cell') ! this is just for development
  call wemission(n,(deltax+x0)*sx0)
  print*, 'abs delta x: max, avg = ', maxval(abs(deltax)), sum(abs(deltax))/n

  call wts('output/cpr.dat', m, yinit+sy*y0)
  call wts('output/cobs.dat', m, yinit+sy*y)
  deltay=matmul(oj,deltax)
  call wts('output/cpost.dat', m, yinit+sy*(y0+deltay))
  print*, 'RMSE cprior = ', sqrt(sum((sy*(y0-y))**2)/m)
  print*, 'RMSE cpost = ', sqrt(sum((sy*(y0+deltay-y))**2)/m)
  write(20,'(a8,6a16)') 'conc #', 'RMSE to prior', 'RMSE to post' , 'reduction' 
  do j = 1,m 
     write(20,'(i8,6(2x,e14.8))') j, sy(j)*sqrt((y0(j)-y(j))**2), sy(j)*sqrt((y0(j)+deltay(j)-y(j))**2), &
          sy(j)*sqrt((y0(j)-y(j))**2)-sy(j)*sqrt((y0(j)+deltay(j)-y(j))**2)
  enddo
  call wconcpost(m,yinit+sy*y0,yinit+sy*(y0+deltay),yinit+sy*y)
  if (twin) then 
     print*, 'twin diagnostics '
     do j = 1, min(5,m)
        print'(a10,i3,6(2x,e14.8))', 'delta (y)', j, y(j), y0(j), y(j)-y0(j), deltay(j),(y(j)-y0(j)-deltay(j))/deltay(j)
     enddo
     print*, 'absdiff delta y: max, avg = ', maxval((abs(y-y0-deltay)/deltay)), sum(abs((y-y0-deltay)/deltay))/m
     print*
     do i = 1, min(5,n)
        print'(a10,i3,6(2x,e14.8))', 'delta (x)', i, deltax(i), xp(i)-x0(i), (xp(i)-x0(i)-deltax(i))/deltax(i)
     enddo
     print*, 'absdiff delta x: max, avg = ', maxval(abs((xp-x0-deltax)/(xp-x0))), sum(abs((xp-x0-deltax)/(xp-x0)))/n
     print* 
     print*, 'obs signal of diff to truth = ', maxval(abs(matmul(oj,(xp-x0-deltax)))), sum(abs(matmul(oj,(xp-x0-deltax))))
     print*, 'relativ to obs = ', maxval(abs(matmul(oj,(xp-x0-deltax))/y)), sum(abs(matmul(oj,(xp-x0-deltax))/y))
  endif
end program inversion
