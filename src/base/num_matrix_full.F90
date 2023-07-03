!
! NAME
!   num_matrix_full   -  defines full matrix
!
! USAGE
!
!   use num_matrix_full
!
!   type(TFullMatrix)   ::  A
!
!   ! initialize matrix:
!   !  o no memory allocated
!   !  o zero flag
!   call Init( A )
!
!   ! (re)allocate memory to store contents:
!   call SetStorage( A, am, an )
!   call ClearStorage( A )
!
!   ! zero flag ?
!   if ( IsZero(A) ) ...
!  
!   ! define full matrix:
!   call SetFull( A, m, n )
!   call SetFull( A, arr(:,:) )
!  
!   ! done
!   call Done( A )
!



module num_matrix_full

  implicit none
  
  
  ! --- in/out ------------------------------------
  
  private
  
  public    ::  TFullMatrix
  public    ::  Init, Done
  public    ::  SetStorage, ClearStorage
  public    ::  IsZero
  public    ::  SetFull
  
  
  ! --- const ----------------------------------

  character(len=*), parameter  ::  mname = 'module num_matrix_full'


  ! --- types ------------------------------------
  
  type TFullMatrix
    ! key name for error mesages
    character(len=20)    ::  key
    ! logical matrix dimension:
    integer              ::  m, n
    ! contents flags
    logical              ::  zero
    ! physical storage:
    real, pointer        ::  a(:,:)
    integer              ::  am, an
    integer              ::  knd
  end type TFullMatrix
  
  
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
  
  interface SetFull
    module procedure mat_SetFull_range
    module procedure mat_SetFull_array
  end interface
  
  
  
  
  
contains


  ! ==============================================================
  
  
  subroutine mat_Init( mat, key )
  
    ! --- in/out ------------------------------
    
    type(TFullMatrix), intent(out)    ::  mat
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
    
    type(TFullMatrix), intent(inout)     ::  mat
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_Done'
    
    ! --- begin -------------------------------
    
    ! clear memory ...
    call ClearStorage( mat )
    
  end subroutine mat_Done



  ! =========================================================  


  subroutine mat_SetStorage( mat, am, an )
  
    ! --- in/out ------------------------------
    
    type(TFullMatrix), intent(inout)   ::  mat
    integer, intent(in)                ::  am, an
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetStorage'
    
    ! --- local -------------------------------
    
    integer        ::  stat
    
    ! --- begin -------------------------------
    
    ! check ...
    if ( (am < 1) .or. (an < 1)  ) then
      write (*,'("ERROR - strange storage definition:")')
      write (*,'("ERROR -   storage      : ",2i6)') am, an
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! set maximum shape:
    mat%am = am
    mat%an = an
    
    ! allocate memory:
    if ( associated(mat%a) ) then
      if ( any( shape(mat%a) /= (/mat%am,mat%an/) ) ) then
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
      allocate( mat%a(mat%am,mat%an), stat=stat )
      if ( stat /= 0 ) then
        write (*,'("ERROR - error during allocation of matrix;")')
        write (*,'("ERROR -   status       : ",i6)') stat
        write (*,'("ERROR -   matrix key   : ",a)') mat%key
        write (*,'("ERROR in ",a)') name; stop
      end if
    end if
        
  end subroutine mat_SetStorage
  
  
  ! ***
  
  
  subroutine mat_ClearStorage( mat )
  
    ! --- in/out ------------------------------
    
    type(TFullMatrix), intent(inout)   ::  mat
    
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
    
    type(TFullMatrix), intent(in)       ::  mat
    
    ! --- begin -------------------------------
    
    mat_IsZero = mat%zero
    
  end function mat_IsZero
  
  

  ! =========================================================  


  subroutine mat_SetFull_range( mat, m, n )
  
    ! --- in/out ------------------------------
    
    type(TFullMatrix), intent(inout)   ::  mat
    integer, intent(in)                ::  m, n
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetFull_range'
    
    ! --- begin -------------------------------
    
    ! memory allocated ?
    if ( .not. associated(mat%a) ) then
      write (*,'("ERROR - no storage allocated")')
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! check ...
    if ( (m < 1) .or. (n < 1) .or. &
         (m > mat%am) .or. (n > mat%an) ) then
      write (*,'("ERROR - strange full matrix definition:")')
      write (*,'("ERROR -   matrix shape : ",2i6)') m, n
      write (*,'("ERROR -   storage      : ",2i6)') mat%am, mat%an
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! logical shape of the matrix:
    mat%m = m
    mat%n = n
    
  end subroutine mat_SetFull_range
  
  
  ! ***
  
  
  subroutine mat_SetFull_array( mat, A, transA )
  
    ! --- in/out ------------------------------
    
    type(TFullMatrix), intent(inout)   ::  mat
    real, intent(in)                   ::  A(:,:)
    character(len=1), intent(in)       ::  transA
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetFull_array'
    
    ! --- local --------------------------------
    
    integer       ::  m, n
    
    ! --- begin -------------------------------
    
    ! extract shape
    m = size(A,1)
    n = size(A,2)
    
    if ( transA == 'N' ) then
    
      ! set shape etc, checks included:
      call SetFull( mat, m, n )
    
      ! fill contents
      mat%a(1:mat%m,1:mat%n) = A
    
    else if ( transA == 'T' ) then
    
      ! set shape etc, checks included:
      call SetFull( mat, n, m )
    
      ! fill contents
      mat%a(1:mat%m,1:mat%n) = transpose(A)
    
    else
    
      write (*,'("ERROR - unknown key for normal/transposed : ",a)') transA
      write (*,'("ERROR in ",a)') name; stop

    end if
    
    ! not zero anymore
    mat%zero = .false.
    
  end subroutine mat_SetFull_array
  
  
end module num_matrix_full
