!#################################################################
!
! Grids.
!
!### macro's #####################################################
!
module dims_grid

  implicit none
  
  ! --- in/out ------------------------------
  
  public
  
  
  ! --- const -------------------------------
  
  
  ! Basic model definition: resolution etc. including some routines
  ! to fill the data structure.

  ! basic (coarsest) resolution in degrees for x and y (dz default 1.0)

  real, parameter     ::  dx = 6.0
  real, parameter     ::  dy = 4.0
  real, parameter     ::  dz = 1.0


  ! Maximum number of zoom regions, 
  ! including the basic (coarsest grid) region;
  ! arrays are allocated for each of these regions:
  integer, parameter  ::  nregions_max = 1
  
  ! Actual number of zoom regions,
  ! during testing this could be set to 1 to quickly run the model.
  integer, parameter :: nregions = 1

  ! region_name is used to recognise the METEO files
  ! region_name is also used in the HDF output file name
  ! region 1 should always be the global domain

  integer, parameter  ::  len_region_name = 10
  character(len=len_region_name), parameter  ::  region_name(1:nregions) = &
       (/ 'glb600x400'/)

  ! coordinates (in degrees) for each region:
  ! xcyc = 1 if the region has cyclic x-boundary conditions
  ! touch_np = 1 if region touches the north pole
  ! touch_sp = 1 if region touches the south pole
  ! xbeg : the westmost border of the region
  ! xend : the eastmost border of the region
  ! ybeg : the southmost border of the region
  ! yend : the northmost border of the region

  integer, parameter  ::  xcyc    (nregions) = (/    1 /)
  integer, parameter  ::  touch_np(nregions) = (/    1 /)
  integer, parameter  ::  touch_sp(nregions) = (/    1 /)
  integer, parameter  ::  xbeg    (nregions) = (/ -180 /)
  integer, parameter  ::  xend    (nregions) = (/  180 /)
  integer, parameter  ::  ybeg    (nregions) = (/  -90 /)
  integer, parameter  ::  yend    (nregions) = (/   90 /)
  integer, parameter  ::  im      (nregions) = (/   60 /)
  integer, parameter  ::  jm      (nregions) = (/   45 /)


  ! maximum refinement factor (can be arbitrary in principle):

  integer, parameter :: maxref = 11

  ! refinement factors for each region (<= maxref)
  ! tref may differ from xref/yref. In the current 
  ! implementation it should be 2,2,4,6,...

  integer, parameter  :: xref(0:nregions) = (/ 1, 1 /)
  integer, parameter  :: yref(0:nregions) = (/ 1, 1 /)
  integer, parameter  :: zref(0:nregions) = (/ 1, 1 /)
  integer, parameter  :: tref(0:nregions) = (/ 1, 1 /)

  ! Define the parent of each region. 
  ! Global region 2 should have parent 0 (globe single cell);
  ! global surface region should have parent 2 (global region).
  integer, parameter  ::  parent(nregions) = (/ 0 /)

end module dims_grid
