!
! NAME
!   num_matrix_diag   -  defines diag matrix
!
! USAGE
!
!   use num_matrix_diag
!
!   type(TDiagMatrix)   ::  A
!
!   ! initialize matrix:
!   !  o no memory allocated
!   !  o zero flag
!   call Init( A, 'matrix A' )
!
!   ! (re)allocate memory to store contents:
!   call SetStorage( A, am )
!   call ClearStorage( A )
!
!   ! zero flag ?
!   if ( IsZero(A) ) ...
!  
!   ! define diag matrix:
!   call SetDiag( A, m, n )
!   call SetDiag( A, arr(:) )    ! define as square matrix
!  
!   ! done
!   call Done( A )
!



module num_matrix_diag

  implicit none
  
  
  ! --- in/out ------------------------------------
  
  private
  
  public    ::  TDiagMatrix
  public    ::  Init, Done
  public    ::  SetStorage, ClearStorage
  public    ::  IsZero
  public    ::  SetDiag
  
  
  ! --- const ----------------------------------

  character(len=*), parameter  ::  mname = 'module num_matrix_diag'


  ! --- types ------------------------------------
  
  type TDiagMatrix
    ! key name for error mesages
    character(len=20)    ::  key
    ! logical matrix dimension:
    integer              ::  m, n
    ! contents flags
    logical              ::  zero
    ! physical storage:
    real, pointer        ::  a(:)
    integer              ::  am
    integer              ::  knd
  end type TDiagMatrix
  
  
  ! --- interfaces ------------------------------
  
  interface Init
    module procedure mat_Init
  end interface
  
  interface Done
    module procedure mat_Done
  end interface
  
  interface SetStorage
    module procedure mat_SetStorage
  end interface
  
  interface ClearStorage
    module procedure mat_ClearStorage
  end interface
  
  interface IsZero
    module procedure mat_IsZero
  end interface
  
  interface SetDiag
    module procedure mat_SetDiag_range
    module procedure mat_SetDiag_array
  end interface
  
  
  
  
  
contains


  ! ==============================================================
  
  
  subroutine mat_Init( mat, key )
  
    ! --- in/out ------------------------------
    
    type(TDiagMatrix), intent(out)    ::  mat
    character(len=*), intent(in)      ::  key
   
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_Init'
    
    ! --- begin -------------------------------
    
    ! store key
    mat%key = key
    
    ! dummy size:
    mat%m = 0
    mat%n = 0
    
    ! start with zero matrix:
    mat%zero     = .true.

    ! initialize physical storage
    nullify( mat%a )
    mat%knd = kind(1.0)
    
  end subroutine mat_Init
  
  
  ! ***
  
  
  subroutine mat_Done( mat )
  
    ! --- in/out ------------------------------
    
    type(TDiagMatrix), intent(inout)     ::  mat
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_Done'
    
    ! --- begin -------------------------------
    
    ! clear memory ...
    call ClearStorage( mat )
    
  end subroutine mat_Done



  ! =========================================================  


  subroutine mat_SetStorage( mat, am )
  
    ! --- in/out ------------------------------
    
    type(TDiagMatrix), intent(inout)   ::  mat
    integer, intent(in)                ::  am
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetStorage'
    
    ! --- local -------------------------------
    
    integer        ::  stat
    
    ! --- begin -------------------------------
    
    ! check ...
    if ( (am < 1)  ) then
      write (*,'("ERROR - strange storage definition:")')
      write (*,'("ERROR -   storage      : ",i6)') am
      write (*,'("ERROR -   matrix key   : ",a)') mat%key
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! set maximum shape:
    mat%am = am
    
    ! allocate memory:
    if ( associated(mat%a) ) then
      if ( size(mat%a) /= mat%am ) then
        deallocate( mat%a, stat=stat )
        if ( stat /= 0 ) then
          write (*,'("ERROR - error during deallocation of matrix;")')
          write (*,'("ERROR -   status       : ",i6)') stat
          write (*,'("ERROR -   matrix key   : ",a)') mat%key
          write (*,'("ERROR in ",a)') name; stop
        end if
      end if
    end if
    if ( .not. associated(mat%a) ) then
      allocate( mat%a(mat%am), stat=stat )
      if ( stat /= 0 ) then
        write (*,'("ERROR - error during allocation of matrix;")')
        write (*,'("ERROR -   status       : ",i6)') stat
        write (*,'("ERROR -   matrix key   : ",a)') mat%key
        write (*,'("ERROR in ",a)') name; stop
      end if
    end if
        
    ! not zero anymore
    mat%zero = .false.
    
  end subroutine mat_SetStorage
  
  
  ! ***
  
  
  subroutine mat_ClearStorage( mat )
  
    ! --- in/out ------------------------------
    
    type(TDiagMatrix), intent(inout)   ::  mat
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_ClearStorage'
    
    ! --- local -------------------------------
    
    integer        ::  stat
    
    ! --- begin -------------------------------
    
    ! remove physical storage ...
    if ( associated(mat%a) ) then
      deallocate( mat%a, stat=stat )
      if ( stat /= 0 ) then
        write (*,'("ERROR - error during deallocation of matrix;")')
        write (*,'("ERROR -   status       : ",i6)') stat
        write (*,'("ERROR -   matrix key   : ",a)') mat%key
        write (*,'("ERROR in ",a)') name; stop
      end if
      nullify( mat%a )
    end if

  end subroutine mat_ClearStorage
  
  
  ! =========================================================  
  
  
  logical function mat_IsZero( mat )
  
    ! --- in/out ------------------------------
    
    type(TDiagMatrix), intent(in)       ::  mat
    
    ! --- begin -------------------------------
    
    mat_IsZero = mat%zero
    
  end function mat_IsZero
  
  

  ! =========================================================  


  subroutine mat_SetDiag_range( mat, m, n )
  
    ! --- in/out ------------------------------
    
    type(TDiagMatrix), intent(inout)   ::  mat
    integer, intent(in)                ::  m, n
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetDiag_range'
    
    ! --- begin -------------------------------
    
    ! memory allocated ?
    if ( .not. associated(mat%a) ) then
      write (*,'("ERROR - no storage allocated")')
      write (*,'("ERROR -   matrix key   : ",a)') mat%key
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! check ...
    if ( (m < 1) .or. (n < 1) .or. &
         (min(m,n) > mat%am) ) then
      write (*,'("ERROR - strange diag matrix definition:")')
      write (*,'("ERROR -   matrix shape : ",2i6)') m, n
      write (*,'("ERROR -   storage      : ",2i6)') mat%am
      write (*,'("ERROR -   matrix key   : ",a)') mat%key
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! logical shape of the matrix:
    mat%m = m
    mat%n = n
    
  end subroutine mat_SetDiag_range
  
  
  ! ***
  
  
  subroutine mat_SetDiag_array( mat, arr )
  
    ! --- in/out ------------------------------
    
    type(TDiagMatrix), intent(inout)   ::  mat
    real, intent(in)                   ::  arr(:)
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetDiag_array'
    
    ! --- local --------------------------------
    
    integer       ::  m, n
    
    ! --- begin -------------------------------
    
    ! extract shape; square matrix with arr(:) as diagonal
    m = size(arr)
    n = size(arr)
    
    ! set shape etc, checks included:
    call SetDiag( mat, m, n )
    
    ! fill contents
    mat%a(1:m) = arr
    
  end subroutine mat_SetDiag_array
  
  
end module num_matrix_diag
