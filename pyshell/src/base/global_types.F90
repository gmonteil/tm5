!###############################################################################
!
! Type definitions.
!
! Extracted from 'dims.F90'.
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

module global_types

  implicit none

  ! --- in/out ------------------------------

  public


  ! --- types -------------------------------

  ! define structures to contain dimensionolized data:

  type region_data
     integer,dimension(:,:),pointer  :: zoomed
     integer,dimension(:,:),pointer  :: edge
     real,dimension(:),pointer       :: dxyp    ! area (m2)
     integer,dimension(:,:),allocatable  :: hor_idx ! list of (i,j) for a region
     integer                               :: n_hor
  end type region_data


  type wind_data

#ifdef MPI
     !
     ! am_k       westward mass flux through east (kg/s)   (halo = 1)
     ! bm_k       northward mass flux through south (kg/s) (halo = 1)
     ! cm_k       upward flux through bottom (kg/s) (halo = 1)
     ! pu_k,pv_k  to read in the pre-calculated fluxes (halo = 1)
     !
     real,dimension(:,:,:),pointer   :: am_k
     real,dimension(:,:,:),pointer   :: bm_k
     real,dimension(:,:,:),pointer   :: cm_k
#endif
     !
     ! am_t       westward mass flux through east (kg/s)   (halo = 1)
     ! bm_t       northward mass flux through south (kg/s) (halo = 1)
     ! cm_t       upward flux through bottom (kg/s) (halo = 1)
     ! pu_t,pv_t  to read in the pre-calculated fluxes (halo = 1)
     !
     real,dimension(:,:,:),pointer   :: am_t
     real,dimension(:,:,:),pointer   :: bm_t
     real,dimension(:,:,:),pointer   :: cm_t
  end type wind_data

  type mass_data

     ! m_k     air mass (kg)   (halo = 2)
     ! p       surface pressure (Pa)   (halo = 2)
     ! phlb_k  pressure at the grid boundaries (Pa) :
     !         l=1 means bottom lowest box (halo = 0)
     ! rm_k    tracer masses (kg)    (halo = 2)
     ! rxm_k   tracer slopes (only transported) in kg/(halfgridsize) (halo = 1)
     ! rym_k   tracer slopes (only transported) in kg/(halfgridsize) (halo = 1)
     ! rzm_k   tracer slopes (only transported) in kg/(halfgridsize) (halo = 1)
     ! m_t     air mass (kg)   (halo = 2)
     ! phlb_t  pressure at the grid boundaries (Pa) :
     !         l=1 means bottom lowest box (halo = 0)
     ! rm_t    tracer masses (kg)    (halo = 2)
     ! rxm_t   tracer slopes (only transported) in kg/(halfgridsize) (halo = 1)
     ! rym_t   tracer slopes (only transported) in kg/(halfgridsize) (halo = 1)
     ! rzm_t   tracer slopes (only transported) in kg/(halfgridsize) (halo = 1)
     !
     real,dimension(:,:,:,:),pointer   :: rm_k
     real,dimension(:,:,:,:),pointer   :: rxm_k
     real,dimension(:,:,:,:),pointer   :: rym_k
     real,dimension(:,:,:,:),pointer   :: rzm_k
#ifdef secmom
     real,dimension(:,:,:,:),pointer   :: rxxm_k
     real,dimension(:,:,:,:),pointer   :: rxym_k
     real,dimension(:,:,:,:),pointer   :: rxzm_k
     real,dimension(:,:,:,:),pointer   :: ryym_k
     real,dimension(:,:,:,:),pointer   :: ryzm_k
     real,dimension(:,:,:,:),pointer   :: rzzm_k
#endif
     real,dimension(:,:,:,:),pointer   :: rm_t
     real,dimension(:,:,:,:),pointer   :: rxm_t
     real,dimension(:,:,:,:),pointer   :: rym_t
     real,dimension(:,:,:,:),pointer   :: rzm_t
#ifdef secmom
     real,dimension(:,:,:,:),pointer   :: rxxm_t
     real,dimension(:,:,:,:),pointer   :: rxym_t
     real,dimension(:,:,:,:),pointer   :: rxzm_t
     real,dimension(:,:,:,:),pointer   :: ryym_t
     real,dimension(:,:,:,:),pointer   :: ryzm_t
     real,dimension(:,:,:,:),pointer   :: rzzm_t
#endif

  end type mass_data


  type conv_data
     !
     ! cloud_base  bottom updraft
     ! cloud_top   top updraft
     ! cloud_lfs   level of free sinking (downdraft)
     !
     integer, dimension(:,:), allocatable   :: cloud_base
     integer, dimension(:,:), allocatable   :: cloud_top
     integer, dimension(:,:), allocatable   :: cloud_lfs
     !
     ! dkg         diffusion coefs.
     ! blh         boudary layer height
     !
     real, dimension(:,:,:), allocatable    :: dkg
     real, dimension(:,:), allocatable      :: blh
     real, dimension(:,:,:), allocatable    :: kvh

  end type conv_data


end module global_types
