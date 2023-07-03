module linalg_interface

public :: eigvals, inverse, matmul_fast, eigensystem

contains

subroutine matmul_fast(mat1,mat2,res)

    real, intent(in) :: mat1(:,:), mat2(:,:)
    real, intent(inout) :: res(:,:)
    integer :: dim1, dim2, dim3 ! mat1 = dim1xdim2, mat2 = dim2xdim3
    real, allocatable :: aux(:)

    dim1 = size(mat1,1)
    dim2 = size(mat1,2)
    dim3 = size(mat2,2)

    if ((size(res,1) /= dim1) .or. (size(res,2) /= dim3)) write(0,*) 'matmul_fast :: the sizes are all wrong!'

!    call dgemms(mat1, dim1, 'N', mat2, dim2, 'N', res, dim1, dim1, dim2, dim3, aux, 0)
    call dgemul(mat1, dim1, 'N', mat2, dim2, 'N', res, dim1, dim1, dim2, dim3)
!    call dgemm('N', 'N', dim1, dim3, dim2, 1.0D0, mat1, dim1, mat2, dim2, 0.0D0, res, dim1)

end subroutine matmul_fast

subroutine eigensystem (A, V)

    real, intent(inout) :: A(:,:), V(:,:)
    integer :: mat_size, i
    real, allocatable :: eigenvalues(:)

    mat_size = size(A,1)
    if ((size(A,2) /= mat_size) .or. (size(V,1) /= mat_size) .or. (size(V,2) /= mat_size)) write(0,*) 'eigensystem :: matrix sizes do not match'

    allocate(eigenvalues(mat_size))

    call eigvals(mat_size, A, eigenvalues, check_positive=.false.)
    ! A should be a diagonal matrix of eigenvalues
    ! V should have the eigenvectors
    V(:,:) = A(:,:)
    A = 0.0
    do i = 1, mat_size
        A(i,i) = eigenvalues(i)
    end do

    deallocate(eigenvalues)

end subroutine eigensystem

subroutine inverse (mat_size, A)

    integer, intent(in) :: mat_size
    real, dimension(mat_size,mat_size), intent(inout) :: A
    integer :: istat
    real, allocatable :: A_copy(:,:), work(:), ipiv(:)
    integer, allocatable :: ipvt(:)

    allocate(A_copy(mat_size, mat_size))
    A_copy = A

    allocate(ipvt(mat_size))
    allocate(work(mat_size))

    call DGETRF(mat_size, mat_size, A_copy, mat_size, ipvt, istat)
    if (istat /= 0) then
        write(0,*) 'DGETRF failed'
        stop
    end if
    call DGETRI(mat_size, A_copy, mat_size, ipvt, work, 0, istat)
    if (istat /= 0) then
        write(0,*) 'DGETRI failed'
        stop
    end if

    deallocate(ipvt, work)

    A = A_copy
    deallocate(A_copy)

end subroutine inverse

subroutine eigvals (mat_size, A, evals, check_positive)

    integer, intent(in) :: mat_size
    real, dimension(mat_size,mat_size), intent(inout) :: A
    real, dimension(mat_size), intent(out) :: evals
    logical, intent(in), optional :: check_positive

    real, allocatable :: A_copy(:,:), work(:), evecs(:,:)
    integer :: istat, len_work, i, conv_evals
    integer, allocatable :: iwork(:), failed_indices(:)
    logical :: check = .false.

    allocate(A_copy(mat_size, mat_size))
    allocate(evecs(mat_size,mat_size))
    A_copy = A

    if (present(check_positive)) check = check_positive

    if (check) then
        call DPOTRF('L', mat_size, A_copy, mat_size, istat)
        if ( istat /= 0 ) then
           write (0,'("WARNING - Cholesky decomposition failed with status ",i5)') istat
           if ( istat > 0 ) then
              write (0,'("ERROR - this leading minor is not positive definite")')
              do i = 1, istat
                 write (*,*) i, A(i,1:istat)
              end do
           end if
        end if
        A_copy = A
    end if

    allocate(iwork(5*mat_size))
    allocate(work(mat_size))
    allocate(failed_indices(mat_size))
    call DSYEVX('V', 'A', 'L', mat_size, A_copy, mat_size, 0.0, 0.0, 0, 0, 1.0D-15, conv_evals, evals, evecs, mat_size, work, 0, iwork, failed_indices, istat)
    if (istat /= 0) then
        write(0,'(i5, " eigenvalues failed to converge")') istat
    end if
    deallocate(iwork, work, failed_indices)

    A = evecs
    deallocate(A_copy, evecs)

end subroutine eigvals

subroutine test_inverse(mat_size, A, A_inv)

    integer, intent(in) :: mat_size
    real, dimension(mat_size,mat_size), intent(in) :: A, A_inv

    real, allocatable :: Identity(:,:)
    real :: diag_sum, off_diag_sum
    integer :: i,j

    allocate(Identity(mat_size,mat_size))

    Identity = matmul(A, A_inv)

    diag_sum = 0.0
    off_diag_sum = 0.0
    do i=1,mat_size
        diag_sum = diag_sum + abs(Identity(i,i))
        do j=1,i-1
            off_diag_sum = off_diag_sum + abs(Identity(i,j))
        end do
        do j=i+1,mat_size
            off_diag_sum = off_diag_sum + abs(Identity(i,j))
        end do
    end do
    write(*,'(a, f14.6, a, f14.8, a)') 'Sum of absolute diagonals = ', diag_sum, ', error of ', abs(mat_size-diag_sum), ' from theoretical.'
    write(*,'(a, f14.8)') 'Sum of absolute off-diagonals = ', off_diag_sum

    deallocate(Identity)

end subroutine test_inverse

subroutine test_eigvals(mat_size, A, eval, evec)

    integer, intent(in) :: mat_size
    real, dimension(mat_size), intent(in) :: eval
    real, dimension(mat_size,mat_size), intent(in) :: A, evec

    real, allocatable :: eig_residues(:,:)
    integer :: i

    allocate(eig_residues(mat_size,mat_size))
    eig_residues = matmul(A, evec)
    do i=1,mat_size
        ! construct residues for eigenvector/eigenvalue #i, which are evec_lapack(:,i) and eval_lapack(i)
        eig_residues(:,i) = eig_residues(:,i) - eval(i) * evec(:,i)
    end do
    write(*,'(a, f20.12)') 'Sum of absolute residuals = ', sum(abs(eig_residues))
    deallocate(eig_residues)

end subroutine test_eigvals

end module linalg_interface
