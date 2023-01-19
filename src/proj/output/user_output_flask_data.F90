#include "tm5.inc"

module user_output_flask_data

  use chem_param,   only : tracer_name_len

  public

  integer, parameter    :: VERT_COORD_ALT = 1
  integer, parameter    :: VERT_COORD_PRES = 2

  integer, parameter    :: POINT_INTERPOLATION_GRIDBOX = 1
  integer, parameter    :: POINT_INTERPOLATION_SLOPES = 2
  integer, parameter    :: POINT_INTERPOLATION_LINEAR = 3

  integer, parameter    :: POINT_ERROR_GRADIENT = 1
  integer, parameter    :: POINT_ERROR_NEIGHBOR = 2

  integer, parameter    :: lowest_sample_layer = 2 ! for stations in the first model layer, this is the layer in which we'll sample them

  ! "id" is meant to be a flask identifier guaranteed to be unique, for the
  ! purpose of matching input file records with output file records.  In
  ! CarbonTracker, it is generated when the observations netCDF file is generated,
  ! and is the record dimension. In PyShell 4DVAR, it is a record dimension,
  ! going from 1 to n_obs.
  !
  ! "station_id" is an integer descriptor meant to give humans an indication of
  ! where each observation comes from. For CarbonTracker, this can mean NOAA
  ! database eventnumbers or strings made up to represent 4-hour averages from
  ! continuous instruments.  It is only present to help identify observations.

  type flask_sample
     integer    :: id  ! unique indentifier placed in input file
     real       :: lat ! sample latitude
     real       :: lon ! sample longitude
     real       :: alt, pres ! sample altitude (meters above sea level, or pressure in Pascals)
     integer    :: vert_coord_type ! whether vertical coordinate is in m.a.s.l. or Pa
     integer    :: itau_start ! sampling start time
     integer    :: itau_center ! sampling center time
     integer    :: itau_end ! sampling end time
     integer    :: sampling_strategy ! 1 for 4-hour averages, 2 for instantaneous
     integer    :: station_id ! another sample identifier; could be NOAA's CCG database event number
     real       :: mix ! mixing ratio of interest
     real       :: var ! variance of the sampled gridbox and nearest neighbors
     integer    :: nsamples ! number of samples accumulated
     real       :: weight   ! unlike nsamples, this accounts for different time step sizes
     integer    :: region ! zoom region index from which this flask is sampled
     integer    :: ifr, jfr ! i,j region indices for flask's grid cell
     integer    :: ifn, jfn ! i,j region indices for flask's "next" grid cell
     real       :: rif, rjf, rlf ! fractions from center of ifr,jfr box
     real       :: surface_height ! surface height in meters
     integer    :: lfr ! vertical level number of the sample
     real       :: wcx, wcy ! x and y weighting factors for slopes interpolation
     logical    :: evaluated ! whether or not the average has been computed from mix/nsamples
     real       :: u,v,blh,q,pressure,temperature ! meteorological variables
     logical    :: below_surface_warning = .false.
  end type flask_sample

  type flask_output_data
    integer, allocatable, dimension(:)  :: flask_id, nsamples, station_id, sampling_strategy
    real, allocatable, dimension(:)     :: avg_time, weight
    real, allocatable, dimension(:)     :: u, v, blh, q, pressure, temperature
    real, allocatable, dimension(:)     :: mix_ratio, mix_ratio_sigma
    integer                             :: n_obs
  end type flask_output_data

  type flask_region_counter
    integer, allocatable :: iflask(:) ! for a region, stores the 'iflask's corresponding to that region
    integer :: n_obs
  end type flask_region_counter

  type tracer_inout
    character(len=tracer_name_len)       :: name
    integer                              :: nflasks
    logical                              :: flask_active=.true.
    real                                 :: flask_minerror
    type(flask_sample), dimension(:), allocatable      :: flasks
    type(flask_output_data), dimension(:), allocatable :: flask_output
    type(flask_region_counter), dimension(:), allocatable :: flask_counter
  end type tracer_inout

  type flask_tracer_forcing
    character(len=tracer_name_len)  :: name                     ! name of the tracer
    logical                         :: has_data                 ! whether there are any forcings for this tracer within a specific zoom region
    real, allocatable               :: rif(:), rjf(:)           ! interpolation coefficients calculated during the forward run
    integer, allocatable            :: ifr(:), jfr(:)           ! gridbox indices for each observation
    integer, allocatable            :: ifn(:), jfn(:)           ! gridbox indices for the "next" grid cell
    real, allocatable               :: wcx(:), wcy(:)           ! x and y weighting factors for slopes interpolation
    real, allocatable               :: lat(:), lon(:), alt(:)   ! latitude, longitude and altitude (meters) of sampling
    integer, allocatable            :: itau_start(:)            ! sampling start time
    integer, allocatable            :: itau_center(:)           ! sampling center time
    integer, allocatable            :: itau_end(:)              ! sampling end time
    integer, allocatable            :: time_window_length(:)    ! custom time window length for averaging
    integer, allocatable            :: times(:,:)               ! 6-component date and time of sampling
    real, allocatable               :: forcing(:)               ! mixing ratio forcing of interest
    integer, allocatable            :: nsamples(:)              ! number of samples accumulated
    real, allocatable               :: total_wt(:)              ! total weight accumulated during the forward run
    integer, allocatable            :: sampling_strategy(:)     ! sampling strategy (1 for 4-hour average, 2 for instantaneous samples, 3 for dT sampling, 4 for custom)
  end type flask_tracer_forcing

  type flask_region_forcing
    logical                                 :: has_data ! whether a particular region has at least one sample
    type(flask_tracer_forcing), allocatable :: io_tracer(:) ! dimension ntracet
  end type flask_region_forcing

end module user_output_flask_data
