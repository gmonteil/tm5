!
!ProTeX: 1.14-AJS
!
!BOI
!
! !TITLE:        PArray - Pointer Array tools
! !AUTHORS:      Arjo Segers, Gijs van Soest
! !AFFILIATION:  KNMI
! !DATE:         \today
!
! !INTRODUCTION: Description
!
!  Module to facilitate allocation and deallocation of pointer arrays.
!
!  Use the module by:
!  \bv
!    use PArray, only : ...
!  \ev
!  First define a pointer array of some type and rank:
!  \bv
!    <type>, pointer      ::  x(:,:,...) 
!  \ev
!  Yet supported:
!  \bv
!    types   : integer(1), integer(2), integer(4), integer(8)
!              real(4), real(8)
!              complex(4), complex(8)
!              logical
!              character(len=<n>)     ! <-- <n> should be a constant !
!              
!    ranks   : 1 - 7
!  \ev
!  Initialize with 'pa_Init', which in fact only nullifies the pointer:
!  \bv
!    call pa_Init( x )
!  \ev
!  Allocate memory for 'x' using one of the following calls:
!  \bv
!    ! define the shape with an integer array:
!    call pa_SetShape( x, (/2,3,4/) )
! 
!    ! define the shape with a series of integer arguments:
!    call pa_SetShape( x, 2, 3, 4 )
! 
!    ! make 'x' a copy of an array 'y' with the same shape and contents;
!    ! 'x' and 'y' should have the same rank:
!    call pa_SetCopy( x, y )
!  \ev
!  Finally, deallocate the memory using the 'pa_Done' routine:
!  \bv
!    call pa_Done( x )
!  \ev
!
!
! !INTRODUCTION: Source files
!
!   The parray module consists of a main and a number of sub modules
!   for different types and kinds.
!   Each type has its own sub module; for types with multiple kinds,
!   sub modules are generated from a template for each individual kind.
!   \bv
!     configure             : script to create source files from template
!                             and Makefile
!
!     parray.f90.in         : template for main module
!
!     parray_iwp.f90.in     : template for integer type
!     parray_rwp.f90.in     : template for real    type
!     parray_cwp.f90.in     : template for complex type
!     parray_l.f90          : sub module for logical type
!     parray_s.f90          : sub module for character string type
!   \ev
!
!
! !INTRODUCTION: History
!
!   See CVS log.
!
!EOI
!

module PArray

  use parray_i1
  use parray_i2
  use parray_i4
  use parray_i8
  use parray_r4
  use parray_r8
  use parray_c4
  use parray_c8
  use parray_l
  use parray_s

  implicit none
  
  ! --- in/out -------------------------
  
  private
  
  public    :: pa_Init, pa_Done, pa_SetShape, pa_SetCopy
  
  
end module PArray
