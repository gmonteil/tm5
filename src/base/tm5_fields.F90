!###############################################################################
!
! Field types and routines.
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module TM5_Fields

  implicit none

  ! --- in/out ------------------------------

  public


  ! --- const ------------------------------

  character(len=*), parameter        :: mname = 'TM5_Fields'


  ! --- types -------------------------------

  ! general 2D field
  type T_Field_2D
    real, allocatable                   :: field(:,:)
  end type T_Field_2D

  ! general 3D field
  type T_Field_3D
    real, allocatable                   :: field(:,:,:) ! (nx,ny,nt)
  end type T_Field_3D

  ! general 4D field:
  type T_Field_4D
     real, allocatable                  ::  field(:,:,:,:)   ! (nx,ny,nz,nt)
  end type T_Field_4D

  ! 4D field on all regions:
  type T_Fields_4D
    type(T_Field_4D), allocatable       ::  cat(:)      ! (ncat)
  end type T_Fields_4D

  type T_Fields_5D
      type(T_Fields_4D), allocatable    ::  tracer(:)
  end type T_Fields_5D

contains


  ! ====================================================================


  subroutine Fields_4D_Init( f, region, layers, ncat, nt, status )

    use Dims, only : im, jm, lm

    ! --- in/out ---------------------------------

    type(T_Fields_4D), intent(out)      ::  f
    integer, intent(in)                 ::  region
    logical, intent(in)                 ::  layers
    integer, intent(in)                 ::  ncat
    integer, intent(in), dimension(ncat)::  nt
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter      :: rname = mname//'/Fields_4D_Init'

    ! --- local ----------------------------------

    integer               ::  i

    ! --- begin ----------------------------------

    ! allocate storage for all regions:
    allocate( f%cat(ncat) )

    ! allocate each cat:
    do i = 1, ncat
      ! allocate layers ?
      if ( layers ) then
        ! allocate time series of 3D fields:
        allocate( f%cat(i)%field(im(region),jm(region),lm(region),nt(i)) )
      else
        ! allocate time series of surface fields:
        allocate( f%cat(i)%field(im(region),jm(region),1    ,nt(i)) )
      end if
    end do

    ! ok
    status = 0

  end subroutine Fields_4D_Init


  ! ***


  subroutine Fields_4D_Done( f, ncat, status )


    ! --- in/out ---------------------------------

    type(T_Fields_4D), intent(inout)    ::  f
    integer, intent(in)                 ::  ncat
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter      :: rname = mname//'/Fields_4D_Done'

    ! --- local ----------------------------------

    integer               ::  i

    ! --- begin ----------------------------------

    ! clear each region:
    do i = 1, ncat
      deallocate( f%cat(i)%field )
    end do

    ! clear region storage:
    deallocate( f%cat )

    ! ok
    status = 0

  end subroutine Fields_4D_Done



end module TM5_Fields
