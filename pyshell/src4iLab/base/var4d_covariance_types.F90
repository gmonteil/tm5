!###############################################################################
!
! Definition of 4D-var covariance data structures.
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

module Var4D_Covariance_Types

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------

  private

  public  ::  Tcov_data
  public  ::  init_covariance, done_covariance


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'Var4D_Covariance_Types'
  
  
  ! --- types ------------------------------------

  ! covariance matrix stored as its eigenvalue decomposition
  type Tcov_data
     real, dimension(:,:), pointer            :: P         ! eigenvectors
     real, dimension(:), pointer              :: sqrt_lam  ! square root of eigenvalues
     real                                     :: corlen    ! decorrelation length
     character(len=1)                         :: choice    ! type of correlation
     integer                                  :: n         ! dimension of the matrix
  end type Tcov_data



  ! --- var --------------------------------------


contains


  !==============================================================================================


  !-----------------------------------------------------------------------
  ! Initialize a covariance structure
  ! Optional input:
  !   o corlen  - correlation length
  !   o choice  - type of correlation
  !-----------------------------------------------------------------------

  subroutine init_covariance( cov, n, status, corlen, choice )

    ! --- modules ---------------------------------------------

    ! --- in/out ----------------------------------------------

    type(Tcov_data), intent(inout)     :: cov
    integer, intent(in)                :: n
    real, intent(in), optional         :: corlen
    character(len=1), intent(in), optional :: choice
    integer, intent(out)               :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter        :: rname = mname//'/init_covariance'

    ! --- local -----------------------------------------------

    ! --- begin -----------------------------------------------

    allocate( cov%P(n,n) )

    allocate( cov%sqrt_lam(n) )

    cov%n = n

    if ( present(corlen) ) then
       cov%corlen = corlen
    else
       cov%corlen = -1.
    end if

    if ( present(choice) ) then
       cov%choice = choice
    else
       cov%choice = '-'
    end if

    status = 0

  end subroutine init_covariance


  !=============================================================================================


  !-----------------------------------------------------------------------
  ! Deallocate a covariance structure
  !-----------------------------------------------------------------------

  subroutine done_covariance( cov, status )

    ! --- modules ---------------------------------------------

    ! --- in/out ----------------------------------------------

    type(Tcov_data), intent(inout)     :: cov
    integer, intent(out)               :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter        :: rname = mname//'/done_covariance'

    ! --- local -----------------------------------------------

    ! --- begin -----------------------------------------------

    deallocate( cov%P )

    deallocate( cov%sqrt_lam )

    status = 0

  end subroutine done_covariance


end module Var4D_Covariance_Types

