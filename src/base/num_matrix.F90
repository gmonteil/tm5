!
! NAME
!   num_matrix   -  matrix tools
!
! USAGE
!
!   use num_matrix
!
!   type(TFullMatrix)   ::  A, U, Vt
!   type(TDiagMatrix)   ::  S
!   type(TBlockMatrix)  ::  B, C
!
!
!   ! matrix-matrix multiplications:
!   !
!   !   A * B  =:  C
!   !
!   call Multiply( A, B, C )
!
!   ! singular value decomposition:
!   !
!   !   A  =  U S V'
!   !
!   call SVD( A, U, S, Vt )
!
!   ! eigen value decomposition:
!   !
!   !   A  =  Q D Q'
!   !
!   call EVD( A, D, Q [,vrange=] [,irange=] )
!

! ***********
!
!   matrix clasifications:
!     rectangular | square
!     symmetric | asymmetric
!     band
!       diagonal
!         identity
!     upper triangular | lower triangular
!


module num_matrix

  use num_matrix_full
  use num_matrix_block
  use num_matrix_diag

  implicit none


  ! --- in/out ------------------------------------

  private

  public    ::  Init, Done
  public    ::  SetStorage, ClearStorage
  public    ::  IsZero

  public    ::  TFullMatrix
  public    ::  SetFull

  public    ::  TBlockMatrix
  public    ::  SetBlock

  public    ::  TDiagMatrix
  public    ::  SetDiag

  public    ::  Multiply
!  public    ::  SVD
  public    ::  EVD



  ! --- const ----------------------------------

  character(len=*), parameter  ::  mname = 'module num_matrix'


  ! --- interfaces ------------------------------

!  interface SVD
!    module procedure mat_SVD
!  end interface

  interface EVD
    module procedure mat_EVD
  end interface

  interface Multiply
    module procedure mat_Multiply_full_block
  end interface

!  interface MultiplyX
!    module procedure mat_MultiplyX_full_block
!  end interface




contains


  ! ==============================================================

  !
  !   A * B  +  0.0 * C  =:  C
  !
  !               B    j1  j2
  !                 [         ]
  !             i1  [  xxxxx  ]
  !             i2  [  xxxxx  ]
  !                 [         ]
  !    A          C
  !      [xxxxxxxx] [  xxxxx  ]
  !      [xxxxxxxx] [  xxxxx  ]
  !
  !                 A(:,i1:i2) B(i1:i2,j1:j2) =  C(:,j1:j2)


  subroutine mat_Multiply_full_block( A, B, C )

    ! --- in/out -------------------------------

    type(TFullMatrix), intent(in)      ::  A
    type(TBlockMatrix), intent(in)     ::  B
    type(TBlockMatrix), intent(inout)  ::  C

    ! --- const ----------------------------------

    character(len=*), parameter  ::  name = mname//', mat_Multiply_full_block'

    ! --- begin --------------------------------

    ! set shape of output matrix and ranges of non-zero block:
    call SetBlock( C, A%m, B%n, 1, A%m, B%j1, B%j2 )

#ifdef with_lapack
    ! call BLAS level 3 routine:
    !
    !   1.0 * A * B  +  0.0 * C  =:  C
    !
    select case ( A%knd )
      case ( 4 )
        call sGeMM( 'N', 'N', A%m, B%bn, B%bm, &
                    1.0, A%a(1,B%i1), size(A%a,1), B%a, size(B%a,1), &
                    0.0, C%a, size(C%a,1) )
      case ( 8 )
        call dGeMM( 'N', 'N', A%m, B%bn, B%bm, &
                    1.0, A%a(1,B%i1), size(A%a,1), B%a, size(B%a,1), &
                    0.0, C%a, size(C%a,1) )
      case default
        write (*,'("ERROR - blas routine ?GeMM not implemented for kind ",i2)') A%knd
        write (*,'("ERROR in ",a)') name; stop
    end select
#else
    write (*,'("ERROR - program should be compiled with lapack")')
    write (*,'("ERROR in ",a)') name; stop
#endif

  end subroutine mat_Multiply_full_block


  !


!  subroutine mat_MultiplyX_full_block( A, B, C )
!
!    ! --- in/out -------------------------------
!
!    type(TFullMatrix), intent(in)      ::  A
!    type(TBlockMatrix), intent(in)     ::  B
!    type(TBlockMatrix), intent(inout)  ::  C
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter  ::  name = mname//', mat_Multiply_full_block'
!
!    ! --- begin --------------------------------
!
!    write (1357,*) '                Multiply:'
!    write (1357,*) '                  A : ',A%m, A%n
!    write (1357,*) '                      ',shape(A%a)
!    write (1357,*) '                  B : ',B%m, B%n
!    write (1357,*) '                      ',B%i1, B%i2
!    write (1357,*) '                      ',B%j1, B%j2
!    write (1357,*) '                      ',shape(B%a)
!
!    write (1357,*) '                set C ...'
!    ! set shape of output matrix and ranges of non-zero block:
!    call SetBlock( C, A%m, B%n, 1, A%m, B%j1, B%j2 )
!
!    ! call BLAS level 3 routine:
!    !
!    !   1.0 * A * B  +  0.0 * C  =:  C
!    !
!    select case ( A%knd )
!      case ( 4 )
!        write (1357,*) '                sGeMM ...'
!        call sGeMM( 'N', 'N', A%m, B%bn, B%bm, &
!                    1.0, A%a(1,B%i1), size(A%a,1), B%a, size(B%a,1), &
!                    0.0, C%a, size(C%a,1) )
!        write (1357,*) '                  ok'
!      case ( 8 )
!        write (1357,*) '                dGeMM ...'
!        call dGeMM( 'N', 'N', A%m, B%bn, B%bm, &
!                    1.0, A%a(1,B%i1), size(A%a,1), B%a, size(B%a,1), &
!                    0.0, C%a, size(C%a,1) )
!        write (1357,*) '                  ok'
!      case default
!        write (*,'("ERROR - blas routine ?GeMM not implemented for kind ",i)') A%knd
!        write (*,'("ERROR in ",a)') name; stop
!    end select
!    write (1357,*) '                done'
!
!  end subroutine mat_MultiplyX_full_block


  ! ================================================================

  !
  ! Singular value decomposition:
  !
  !     A    =   U      S       V'
  !
  !   m x n    m x s  s x s   s x n
  !

!  subroutine mat_SVD( AA, U, S, Vt )

!    ! --- in/out ---------------------------------

!    type(TFullMatrix), intent(in)       ::  AA
!    type(TFullMatrix), intent(out)      ::  U
!    type(TDiagMatrix), intent(out)      ::  S
!    type(TFullMatrix), intent(out)      ::  Vt

!    ! --- const ----------------------------------

!    character(len=*), parameter  ::  name = mname//', mat_SVD'

!    ! --- local ----------------------------------

!    integer               ::  m, n, ns

!    type(TFullMatrix)     ::  A
!    real, allocatable     ::  work(:)
!    integer, allocatable  ::  iwork(:)

!    integer               ::  info

!    ! --- begin ----------------------------------

!    ! copy of input matrix:
!    call Init( A, '(mat_SVD) A' )
!    call SetStorage( A, AA%m, AA%n )
!    call SetFull( A, AA%a, 'N' )

!    ! input size:
!    m = A%m
!    n = A%n

!    ! number of singular values:
!    ns = max(1,min(m,n))

!    ! setup output:
!    call SetStorage( U , m , ns )
!    call SetStorage( S , ns     )
!    call SetStorage( Vt, ns, n  )

!    ! work space
!    allocate( work(4*(ns)**2+max(m,n)+9*ns) )
!    allocate( iwork(8*ns) )

!    !
!    ! decompose
!    !
!#ifdef with_lapack
!    select case ( A%knd )
!      case ( 4 )
!        call sGeSdD( 'S', m, n, A%a, size(A%a,1), &
!                     S%a, U%a, size(U%a,1), Vt%a, size(Vt%a,1), &
!                     work, size(work), iwork, info )
!      case ( 8 )
!        call dGeSdD( 'S', m, n, A%a, size(A%a,1), &
!                      S%a, U%a, size(U%a,1), Vt%a, size(Vt%a,1), &
!                      work, size(work), iwork, info )
!      case default
!        write (*,'("ERROR - lapack routine ?GeSvdD not implemented for kind ",i2)') A%knd
!        write (*,'("ERROR in ",a)') name; stop
!    end select
!    if ( info /= 0 ) then
!      write (*,'("ERROR - from ?GeSvdD; info=",i6)') info
!      write (*,'("ERROR in ",a)') name; stop
!    end if
!#else
!    write (*,'("ERROR - program should be compiled with lapack")')
!    write (*,'("ERROR in ",a)') name; stop
!#endif

!    ! output filled now ...
!    S%zero  = .false.
!    U%zero  = .false.
!    Vt%zero = .false.

!    ! done
!    call Done( A )
!    deallocate(  work )
!    deallocate( iwork )

!  end subroutine mat_SVD


  ! ====================================================


  subroutine mat_EVD( AA, DD, QQ, vrange, irange )

    ! --- in/out ---------------------------------

    type(TFullMatrix), intent(in)      ::  AA
    type(TDiagMatrix), intent(out)     ::  DD
    type(TFullMatrix), intent(out)     ::  QQ
    real, intent(in), optional         ::  vrange(2)
    integer, intent(in), optional      ::  irange(2)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  name = mname//', mat_EVD'

    ! --- local ----------------------------------

    real, allocatable     ::  A(:,:)

    integer               ::  n

    character(len=1)      ::  range
    real                  ::  vl, vu
    integer               ::  il, iu

    real                  ::  abstol

    integer               ::  neigval

    real, allocatable     ::  work(:)
    integer, allocatable  ::  iwork(:)
    integer, allocatable  ::  ifail(:)

    integer               ::  info

    ! --- begin ----------------------------------

    ! check ...
    if ( (AA%m < 1) .or. (AA%n < 1) .or. (AA%m /= AA%n) ) then
      write (*,'("ERROR - matrix strange or not square:")')
      write (*,'("ERROR -   m,n : ",2i6)') AA%m, AA%n
      write (*,'("ERROR in ",a)') name; stop
    end if

    ! input size:
    n = AA%n

    ! copy input
    allocate( A(n,n) )
    A = AA%a

    ! eigenvalues and vectors:
    call SetStorage( DD, n )
    call SetStorage( QQ, n, n )

    ! work space
    allocate(  work(8*n) )  ;  work  = 0.0
    allocate( iwork(5*n) )  ;  iwork = 0
    allocate( ifail(  n) )  ;  ifail = 0

    if ( present(vrange) ) then
      range = 'V'
      vl = vrange(1)
      vu = vrange(2)
      il = 1
      iu = n
    else if ( present(irange) ) then
      range = 'I'
      vl = -huge(1.0)
      vu =  huge(1.0)
      il = irange(1)
      iu = irange(2)
    else
      range = 'A'
      vl = -huge(1.0)
      vu =  huge(1.0)
      il = 1
      iu = n
    end if

    ! default tolerance
    abstol = -1.0

    ! decompose:
#ifdef with_lapack
    select case ( kind(A) )
      case ( 4 )
        call sSyEvX ( 'V', range, 'U', n, A, size(A,1), vl, vu, il, iu, abstol, &
                      neigval, DD%a, QQ%a, size(QQ%a,1), &
                      work, size(work), iwork, ifail, info)
      case ( 8 )
        call dSyEvX ( 'V', range, 'U', n, A, size(A,1), vl, vu, il, iu, abstol, &
                      neigval, DD%a, QQ%a, size(QQ%a,1), &
                      work, size(work), iwork, ifail, info)
      case default
        write (*,'("ERROR - lapack routine ?SyEvX not implemented for kind ",i6)') kind(A)
        write (*,'("ERROR in ",a)') name; stop
    end select
    if ( info /= 0 ) then
      write (*,'("ERROR - from ?SyEvX; info=",i6)') info
      write (*,'("ERROR -   ifail=",i6)') ifail(1:neigval)
      write (*,'("ERROR in ",a)') name; stop
    end if
#else
    write (*,'("ERROR - program should be compiled with lapack")')
    write (*,'("ERROR in ",a)') name; stop
#endif

    ! 'truncate' diagonal matrix (acutal storage is unchanged):
    DD%m = neigval
    DD%n = neigval

    ! 'truncate' eigenvector matrix (acutal storage is unchanged):
    QQ%m = n
    QQ%n = neigval

    ! ok
    deallocate( A )
    deallocate( work  )
    deallocate( iwork )
    deallocate( ifail )

  end subroutine mat_EVD



end module num_matrix
