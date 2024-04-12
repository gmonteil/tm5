!
! PArray
!
!
! Template for PArray routines with integer arguments.
!
! To generate kind specific versions, use:
!
!    sed -e 's/2/1/g' parray_iwp.f90 > parray_i1.f90
!    sed -e 's/2/2/g' parray_iwp.f90 > parray_i2.f90
!    sed -e 's/2/4/g' parray_iwp.f90 > parray_i4.f90
!    sed -e 's/2/8/g' parray_iwp.f90 > parray_i8.f90
!

module PArray_i2

  implicit none
  
  ! --- in/out -------------------------
  
  private
  
  public    :: pa_Init, pa_Done, pa_SetShape, pa_SetCopy
  
  
  ! --- interfaces ---------------------------
  
  interface pa_Init
    module procedure pa_Init_i2_1d
    module procedure pa_Init_i2_2d
    module procedure pa_Init_i2_3d
    module procedure pa_Init_i2_4d
    module procedure pa_Init_i2_5d
    module procedure pa_Init_i2_6d
    module procedure pa_Init_i2_7d
  end interface
  
  interface pa_Done
    module procedure pa_Done_i2_1d
    module procedure pa_Done_i2_2d
    module procedure pa_Done_i2_3d
    module procedure pa_Done_i2_4d
    module procedure pa_Done_i2_5d
    module procedure pa_Done_i2_6d
    module procedure pa_Done_i2_7d
  end interface
  
  interface pa_SetShape
    module procedure pa_SetShape_i2_1d_shp
    module procedure pa_SetShape_i2_1d_n
    module procedure pa_SetShape_i2_2d_shp
    module procedure pa_SetShape_i2_2d_n
    module procedure pa_SetShape_i2_3d_shp
    module procedure pa_SetShape_i2_3d_n
    module procedure pa_SetShape_i2_4d_shp
    module procedure pa_SetShape_i2_4d_n
    module procedure pa_SetShape_i2_5d_shp
    module procedure pa_SetShape_i2_5d_n
    module procedure pa_SetShape_i2_6d_shp
    module procedure pa_SetShape_i2_6d_n
    module procedure pa_SetShape_i2_7d_shp
    module procedure pa_SetShape_i2_7d_n
  end interface
  
  interface pa_SetCopy
    module procedure pa_SetCopy_i2_1d
    module procedure pa_SetCopy_i2_2d
    module procedure pa_SetCopy_i2_3d
    module procedure pa_SetCopy_i2_4d
    module procedure pa_SetCopy_i2_5d
    module procedure pa_SetCopy_i2_6d
    module procedure pa_SetCopy_i2_7d
  end interface
  
  
contains


  ! =========================================================
  ! ===
  ! === integer(2)
  ! ===
  ! =========================================================
  
  
  ! *******************************************
  ! ***
  ! *** integer(2) 1D
  ! ***
  ! *******************************************
  
  
  subroutine pa_Init_i2_1d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:)
    
    ! --- begin ---------------------------
    
    nullify( x )
    
  end subroutine pa_Init_i2_1d


  ! ***

  
  subroutine pa_Done_i2_1d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) deallocate( x )
    
  end subroutine pa_Done_i2_1d


  ! ***
  
  
  subroutine pa_SetShape_i2_1d_shp( x, n )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:)
    integer, intent(in)   ::  n(1)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( any( shape(x) /= n ) ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n(1)) )
    
  end subroutine pa_SetShape_i2_1d_shp


  ! ***


  subroutine pa_SetShape_i2_1d_n( x, n1 )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:)
    integer, intent(in)   ::  n1
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( size(x) /= n1 ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n1) )
    
  end subroutine pa_SetShape_i2_1d_n


  ! ***
  
  
  subroutine pa_SetCopy_i2_1d( x, y )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:)
    integer(2), intent(in)  ::  y(:)
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, shape(y) )
    x = y
    
  end subroutine pa_SetCopy_i2_1d


  ! *******************************************
  ! ***
  ! *** integer(2) 2D
  ! ***
  ! *******************************************
  
  
  subroutine pa_Init_i2_2d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:)
    
    ! --- begin ---------------------------
    
    nullify( x )
    
  end subroutine pa_Init_i2_2d


  ! ***
  
  
  subroutine pa_Done_i2_2d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) deallocate( x )
    
  end subroutine pa_Done_i2_2d


  ! ***
  
  
  subroutine pa_SetShape_i2_2d_shp( x, n )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:)
    integer, intent(in)   ::  n(2)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( any( shape(x) /= n ) ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n(1),n(2)) )
    
  end subroutine pa_SetShape_i2_2d_shp


  ! ***
  
  
  subroutine pa_SetShape_i2_2d_n( x, n1, n2 )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:)
    integer, intent(in)   ::  n1, n2
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, (/n1,n2/) )
    
  end subroutine pa_SetShape_i2_2d_n


  ! ***
  
  
  subroutine pa_SetCopy_i2_2d( x, y )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:)
    integer(2), intent(in)  ::  y(:,:)
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, shape(y) )
    x = y
    
  end subroutine pa_SetCopy_i2_2d


  ! *******************************************
  ! ***
  ! *** integer(2) 3D
  ! ***
  ! *******************************************
  
  
  subroutine pa_Init_i2_3d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:)
    
    ! --- begin ---------------------------
    
    nullify( x )
    
  end subroutine pa_Init_i2_3d


  ! ***
  
  
  subroutine pa_Done_i2_3d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) deallocate( x )
    
  end subroutine pa_Done_i2_3d


  ! ***
  
  
  subroutine pa_SetShape_i2_3d_shp( x, n )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:)
    integer, intent(in)   ::  n(3)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( any( shape(x) /= n ) ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n(1),n(2),n(3)) )
    
  end subroutine pa_SetShape_i2_3d_shp


  ! ***
  
  
  subroutine pa_SetShape_i2_3d_n( x, n1, n2, n3 )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:)
    integer, intent(in)   ::  n1, n2, n3
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, (/n1,n2,n3/) )
    
  end subroutine pa_SetShape_i2_3d_n


  ! ***
  
  
  subroutine pa_SetCopy_i2_3d( x, y )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:)
    integer(2), intent(in)  ::  y(:,:,:)
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, shape(y) )
    x = y
    
  end subroutine pa_SetCopy_i2_3d


  ! *******************************************
  ! ***
  ! *** integer(2) 4D
  ! ***
  ! *******************************************
  
  
  subroutine pa_Init_i2_4d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:)
    
    ! --- begin ---------------------------
    
    nullify( x )
    
  end subroutine pa_Init_i2_4d


  ! ***
  
  
  subroutine pa_Done_i2_4d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) deallocate( x )
    
  end subroutine pa_Done_i2_4d


  ! ***
  
  
  subroutine pa_SetShape_i2_4d_shp( x, n )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:)
    integer, intent(in)   ::  n(4)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( any( shape(x) /= n ) ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n(1),n(2),n(3),n(4)) )
    
  end subroutine pa_SetShape_i2_4d_shp


  ! ***
  
  
  subroutine pa_SetShape_i2_4d_n( x, n1, n2, n3, n4 )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:)
    integer, intent(in)   ::  n1, n2, n3, n4
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, (/n1,n2,n3,n4/) )
    
  end subroutine pa_SetShape_i2_4d_n


  ! ***
  
  
  subroutine pa_SetCopy_i2_4d( x, y )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:)
    integer(2), intent(in)  ::  y(:,:,:,:)
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, shape(y) )
    x = y
    
  end subroutine pa_SetCopy_i2_4d


  ! *******************************************
  ! ***
  ! *** integer(2) 5D
  ! ***
  ! *******************************************
  
  
  subroutine pa_Init_i2_5d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    nullify( x )
    
  end subroutine pa_Init_i2_5d


  ! ***
  
  
  subroutine pa_Done_i2_5d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) deallocate( x )
    
  end subroutine pa_Done_i2_5d


  ! ***
  
  
  subroutine pa_SetShape_i2_5d_shp( x, n )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:,:)
    integer, intent(in)   ::  n(5)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( any( shape(x) /= n ) ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n(1),n(2),n(3),n(4),n(5)) )
    
  end subroutine pa_SetShape_i2_5d_shp


  ! ***
  
  
  subroutine pa_SetShape_i2_5d_n( x, n1, n2, n3, n4, n5 )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:,:)
    integer, intent(in)   ::  n1, n2, n3, n4, n5
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, (/n1,n2,n3,n4,n5/) )
    
  end subroutine pa_SetShape_i2_5d_n


  ! ***
  
  
  subroutine pa_SetCopy_i2_5d( x, y )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:)
    integer(2), intent(in)  ::  y(:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, shape(y) )
    x = y
    
  end subroutine pa_SetCopy_i2_5d


  ! *******************************************
  ! ***
  ! *** integer(2) 6D
  ! ***
  ! *******************************************
  
  
  subroutine pa_Init_i2_6d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    nullify( x )
    
  end subroutine pa_Init_i2_6d


  ! ***
  
  
  subroutine pa_Done_i2_6d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) deallocate( x )
    
  end subroutine pa_Done_i2_6d


  ! ***
  
  
  subroutine pa_SetShape_i2_6d_shp( x, n )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:,:,:)
    integer, intent(in)   ::  n(6)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( any( shape(x) /= n ) ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n(1),n(2),n(3),n(4),n(5),n(6)) )
    
  end subroutine pa_SetShape_i2_6d_shp


  ! ***
  
  
  subroutine pa_SetShape_i2_6d_n( x, n1, n2, n3, n4, n5, n6 )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:,:,:)
    integer, intent(in)   ::  n1, n2, n3, n4, n5, n6
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, (/n1,n2,n3,n4,n5,n6/) )
    
  end subroutine pa_SetShape_i2_6d_n


  ! ***
  
  
  subroutine pa_SetCopy_i2_6d( x, y )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:,:)
    integer(2), intent(in)  ::  y(:,:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, shape(y) )
    x = y
    
  end subroutine pa_SetCopy_i2_6d


  ! *******************************************
  ! ***
  ! *** integer(2) 7D
  ! ***
  ! *******************************************
  
  
  subroutine pa_Init_i2_7d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    nullify( x )
    
  end subroutine pa_Init_i2_7d


  ! ***
  
  
  subroutine pa_Done_i2_7d( x )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) deallocate( x )
    
  end subroutine pa_Done_i2_7d


  ! ***
  
  
  subroutine pa_SetShape_i2_7d_shp( x, n )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:,:,:,:)
    integer, intent(in)   ::  n(7)
    
    ! --- begin ---------------------------
    
    if ( associated(x) ) then
      if ( any( shape(x) /= n ) ) deallocate( x )
    end if
    if ( .not. associated(x) ) allocate( x(n(1),n(2),n(3),n(4),n(5),n(6),n(7)) )
    
  end subroutine pa_SetShape_i2_7d_shp


  ! ***
  
  
  subroutine pa_SetShape_i2_7d_n( x, n1, n2, n3, n4, n5, n6, n7 )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer  ::  x(:,:,:,:,:,:,:)
    integer, intent(in)   ::  n1, n2, n3, n4, n5, n6, n7
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, (/n1,n2,n3,n4,n5,n6,n7/) )
    
  end subroutine pa_SetShape_i2_7d_n


  ! ***
  
  
  subroutine pa_SetCopy_i2_7d( x, y )
  
    ! --- in/out ---------------------------
    
    integer(2), pointer     ::  x(:,:,:,:,:,:,:)
    integer(2), intent(in)  ::  y(:,:,:,:,:,:,:)
    
    ! --- begin ---------------------------
    
    call pa_SetShape( x, shape(y) )
    x = y
    
  end subroutine pa_SetCopy_i2_7d




end module PArray_i2
