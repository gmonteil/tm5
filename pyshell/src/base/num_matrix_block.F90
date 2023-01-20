!
! NAME
!   num_matrix_block   -  defines block matrix
!
! USAGE
!
!   use num_matrix_block
!
!   type(TBlockMatrix)   ::  A
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
!   ! define block:
!   call SetBlock( A, m, n, i1,i2, j1,j2 )
!   call SetBlock( A, m, n, i1, j1, arr(:,:) )
!  
!   ! done
!   call Done( A )
!



module num_matrix_block

  implicit none
  
  
  ! --- in/out ------------------------------------
  
  private
  
  public    ::  TBlockMatrix
  public    ::  Init, Done
  public    ::  SetStorage, ClearStorage
  public    ::  IsZero
  public    ::  SetBlock
  
  
  ! --- const ----------------------------------

  character(len=*), parameter  ::  mname = 'module num_matrix_block'


  ! --- types ------------------------------------
  
  type TBlockMatrix
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
    !
    ! ** block location and size:
    integer              ::  i0, j0
    integer              ::  bm, bn
    integer              ::  i1, i2, j1, j2
  end type TBlockMatrix
  
  
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
  
  interface SetBlock
    module procedure mat_SetBlock_range
    module procedure mat_SetBlock_array
  end interface
  
  
  
  
  
contains


  ! ==============================================================
  
  
  subroutine mat_Init( mat, key )
  
    ! --- in/out ------------------------------
    
    type(TBlockMatrix), intent(out)    ::  mat
    character(len=*), intent(in)       ::  key
    
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
    
    ! dummy block 
    mat%i0 = 0
    mat%j0 = 0
    mat%i1 = 0
    mat%i2 = 0
    mat%j1 = 0
    mat%j2 = 0
    mat%bm = 0
    mat%bn = 0
    
  end subroutine mat_Init
  
  
  ! ***
  
  
  subroutine mat_Done( mat )
  
    ! --- in/out ------------------------------
    
    type(TBlockMatrix), intent(inout)     ::  mat
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_Done'
    
    ! --- begin -------------------------------
    
    ! clear memory ...
    call ClearStorage( mat )
    
  end subroutine mat_Done



  ! =========================================================  


  subroutine mat_SetStorage( mat, am, an )
  
    ! --- in/out ------------------------------
    
    type(TBlockMatrix), intent(inout)  ::  mat
    integer, intent(in)                ::  am, an
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetStorage'
    
    ! --- local -------------------------------
    
    integer        ::  stat
    
    ! --- begin -------------------------------
    
    ! check ...
    if ( (am < 1) .or. (an < 1)  ) then
      write (*,'("ERROR - strange storage definition:")')
      write (*,'("ERROR -   storage      : ",2i7)') am, an
      write (*,'("ERROR -   matrix key   : ",a)') mat%key
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! allocate memory:
    if ( associated(mat%a) ) then
      if ( any( shape(mat%a) /= (/am,an/) ) ) then
        deallocate( mat%a, stat=stat )
        if ( stat /= 0 ) then
          write (*,'("ERROR - error during deallocation of matrix;")')
          write (*,'("ERROR -   matrix key   : ",a)') mat%key
          write (*,'("ERROR in ",a)') name; stop
        end if
      end if
    end if
    if ( .not. associated(mat%a) ) then
      allocate( mat%a(am,an), stat=stat )
      if ( stat /= 0 ) then
        write (*,'("ERROR - error during allocation of matrix;")')
        write (*,'("ERROR -   matrix key   : ",a)') mat%key
        write (*,'("ERROR in ",a)') name; stop
      end if
    end if
        
    ! set maximum shape:
    mat%am = am
    mat%an = an
    
  end subroutine mat_SetStorage
  
  
  ! ***
  
  
  subroutine mat_ClearStorage( mat )
  
    ! --- in/out ------------------------------
    
    type(TBlockMatrix), intent(inout)  ::  mat
    
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
        write (*,'("ERROR -   matrix key   : ",a)') mat%key
        write (*,'("ERROR in ",a)') name; stop
      end if
      nullify( mat%a )
    end if

  end subroutine mat_ClearStorage
  
  
  ! =========================================================  
  
  
  logical function mat_IsZero( mat )
  
    ! --- in/out ------------------------------
    
    type(TBlockMatrix), intent(in)       ::  mat
    
    ! --- begin -------------------------------
    
    mat_IsZero = mat%zero
    
  end function mat_IsZero
  
  

  ! =========================================================  


  subroutine mat_SetBlock_range( mat, m, n, i1, i2, j1, j2 )
  
    ! --- in/out ------------------------------
    
    type(TBlockMatrix), intent(inout)  ::  mat
    integer, intent(in)                ::  m, n
    integer, intent(in)                ::  i1, i2, j1, j2
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetBlock_range'
    
    ! --- begin -------------------------------
    
    ! memory allocated ?
    if ( .not. associated(mat%a) ) then
      write (*,'("ERROR - no storage allocated")')
      write (*,'("ERROR -   matrix key   : ",a)') mat%key
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! check ...
    if ( (m < 1) .or. (n < 1) .or. &
         (i1 < 1) .or. (i1 > m) .or. (i2 < i1) .or. (i2 > m) .or. &
         (j1 < 1) .or. (j1 > n) .or. (j2 < j1) .or. (j2 > n) .or. &
         (i2-i1+1 > mat%am) .or. (j2-j1+1 > mat%an) ) then
      write (*,'("ERROR - strange block definition:")')
      write (*,'("ERROR -   matrix shape : ",2i6)') m, n
      write (*,'("ERROR -   block i1, i2 : ",2i6)') i1, i2
      write (*,'("ERROR -         j1, j2 : ",2i6)') j1, j2
      write (*,'("ERROR -   storage      : ",2i6)') mat%am, mat%an
      write (*,'("ERROR -   matrix key   : ",a)') mat%key
      write (*,'("ERROR in ",a)') name; stop
    end if
    
    ! logical shape of the matrix:
    mat%m = m
    mat%n = n
    
    ! block range:
    mat%i1 = i1
    mat%i2 = i2
    mat%j1 = j1
    mat%j2 = j2
    
    ! block base position
    mat%i0 = mat%i1 - 1
    mat%j0 = mat%j1 - 1
    
    ! block size
    mat%bm = mat%i2 - mat%i1 + 1
    mat%bn = mat%j2 - mat%j1 + 1
    
  end subroutine mat_SetBlock_range
  
  
  ! ***
  
  
  subroutine mat_SetBlock_array( mat, m, n, i1, j1, arr )
  
    ! --- in/out ------------------------------
    
    type(TBlockMatrix), intent(inout)  ::  mat
    integer, intent(in)                ::  m, n
    integer, intent(in)                ::  i1, j1
    real, intent(in)                   ::  arr(:,:)
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  name = mname//', mat_SetBlock_array'
    
    ! --- local -------------------------------
    
    integer           ::  i2, j2
    
    ! --- begin -------------------------------
    
    ! end of ranges:
    i2 = i1 + size(arr,1) - 1
    j2 = j1 + size(arr,2) - 1
    
    ! set block position, checks included:
    call SetBlock( mat, m, n, i1, i2, j1, j2 )
    
    ! fill contents
    mat%a(1:mat%bm,1:mat%bn) = arr
    
    ! not zero anymore
    mat%zero = .false.
    
  end subroutine mat_SetBlock_array
  
  
end module num_matrix_block
