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
  integer, parameter  ::  nregions_max = 3
  
  ! Actual number of zoom regions,
  ! during testing this could be set to 1 to quickly run the model.
  integer, parameter :: nregions = 3

  ! region_name is used to recognise the METEO files
  ! region_name is also used in the HDF output file name
  ! region 1 should always be the global domain

  integer, parameter  ::  len_region_name = 10
  !-- FIT-IC: high resolution (1x1) for
  !           Germany, Netherlands, and Suisse,
  !           which fit into rectangular bounds 2E-16E x 44N-56N
  !           In order not to intersect with the 6x4 global grid-cells
  !           the following bounds are used: 0E-18E x 42N-58N
  character(len=len_region_name), parameter  ::  region_name(1:nregions) = &
       (/ 'glb600x400', 'eur300x200', 'gns100x100'/)

  ! coordinates (in degrees) for each region:
  ! xcyc = 1 if the region has cyclic x-boundary conditions
  ! touch_np = 1 if region touches the north pole
  ! touch_sp = 1 if region touches the south pole
  ! xbeg : the westmost border of the region
  ! xend : the eastmost border of the region
  ! ybeg : the southmost border of the region
  ! yend : the northmost border of the region

  integer, parameter  ::  xcyc    (nregions) = (/    1,   0,  0 /)
  integer, parameter  ::  touch_np(nregions) = (/    1,   0,  0 /)
  integer, parameter  ::  touch_sp(nregions) = (/    1,   0,  0 /)
  integer, parameter  ::  xbeg    (nregions) = (/ -180, -36,  0 /)
  integer, parameter  ::  xend    (nregions) = (/  180,  54, 18 /)
  integer, parameter  ::  ybeg    (nregions) = (/  -90,  22, 42 /)
  integer, parameter  ::  yend    (nregions) = (/   90,  74, 58 /)
  integer, parameter  ::  im      (nregions) = (/   60,  30, 18 /)
  integer, parameter  ::  jm      (nregions) = (/   45,  26, 16 /)


  ! maximum refinement factor (can be arbitrary in principle):

  integer, parameter :: maxref = 11

  ! refinement factors for each region (<= maxref)
  ! tref may differ from xref/yref. In the current 
  ! implementation it should be 2,2,4,6,...
  !MVO-NOTE::looking at the code (initexit.F90, l1314ff) the
  !          refinement factor is w.r.t. to the basic resolution (here: 6x4)
  integer, parameter  :: xref(0:nregions) = (/ 1, 1, 2, 6/)
  integer, parameter  :: yref(0:nregions) = (/ 1, 1, 2, 4/)
  integer, parameter  :: zref(0:nregions) = (/ 1, 1, 1, 1/)
  ! integer, parameter  :: tref(0:nregions) = (/ 1, 1, 2, 1/)!MVO:unclear whether this must be changed as well
  integer, parameter  :: tref(0:nregions) = (/ 1, 1, 1, 1/)!MVO:unclear whether this must be changed as well

  ! Define the parent of each region. 
  ! Global region 2 should have parent 0 (globe single cell);
  ! global surface region should have parent 2 (global region).
  integer, parameter  ::  parent(nregions) = (/ 0, 1, 2/)

end module dims_grid
