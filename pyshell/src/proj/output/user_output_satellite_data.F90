#include "tm5.inc"

module user_output_satellite_data

  use chem_param,   only : tracer_name_len

  public

  integer, parameter    :: SAT_INTERPOLATION_GRIDBOX = 1
  integer, parameter    :: SAT_INTERPOLATION_SLOPES = 2
  integer, parameter    :: SAT_INTERPOLATION_LINEAR = 3

  type T_output_fields
      real, allocatable         :: sampled_profiles(:,:), std_sampled_profiles(:,:)
      !real, allocatable         :: sampled_psurf(:), modeled_mixing(:), std_modeled_mixing(:), observed_mixing(:), std_observed_mixing(:)
      real, allocatable         :: sampled_psurf(:), modeled_mixing(:), std_modeled_mixing(:)
      integer, allocatable      :: input_positions(:), instr_positions(:)
      integer(2), allocatable   :: sampling_strategy(:), nsamples(:)
      real, allocatable         :: weight(:)
      integer                   :: nobs ! # of observations
      real, allocatable         :: t(:,:), gph(:,:), q(:,:) ! temperature, geopotential height, specific humidity
  end type T_output_fields

  type T_satellite_sample
    integer     :: serial ! unique ID per sample
    integer     :: instr_serial ! sequence within instrument
    real        :: lat ! latitude
    real        :: lon ! longitude
    integer     :: sample_itau ! sampling time in TM5 time coordinate
    integer     :: itau_start, itau_end ! time boundaries for when this measurement should be assimilated
    integer(2)  :: sampling_strategy ! 2 for instantaneous, 3 for ndyn
    !real        :: obs_mix, obs_std, mod_mix, mod_var ! total column mixing ratio and standard deviation, as measured and modeled
    real        :: mod_mix, mod_var ! total column mixing ratio and standard deviation, as modeled
    integer(2)  :: nsamples ! # of times a particular sample is sampled
    real        :: weight
    integer     :: region ! zoom region
    integer     :: ifr, jfr ! i,j indices for sample's grid cell
    integer     :: ifn, jfn ! i,j indices for the neighboring grid cell
    real        :: rif, rjf ! fractional distance from the center of cell ifr,jfr
    real        :: wcx, wcy ! x and y weighting factors for slopes interpolation
    logical     :: evaluated ! whether or not the average has been computed from mix/nsamples
    real        :: pressure ! surface pressure
    real, allocatable :: mod_profile(:), var_mod_profile(:) ! modeled mixing ratio
    real, allocatable :: gph(:), t(:), q(:) ! geopotential height, temperature, specific humidity
  end type T_satellite_sample

  type T_satellite_region_counter
    integer, allocatable :: iobs(:) ! for a region, stores the 'iobs's corresponding to that region
    integer :: n_obs
  end type T_satellite_region_counter

  type T_tracer_inout
    character(len=tracer_name_len)                  :: name ! name of the tracer
    logical                                         :: has_data ! whether there are any data for this tracer
    logical                                         :: output_meteo ! whether to output meteo fields as well
    integer                                         :: n_obs
    type(T_output_fields), allocatable              :: out_field(:)
    type(T_satellite_sample), allocatable           :: sat_obs(:)
    type(T_satellite_region_counter), allocatable   :: sat_reg_counter(:)
  end type T_tracer_inout

  type sat_tracer_forcing
    character(len=tracer_name_len)  :: name  ! the name of the tracer
    logical                         :: has_data ! whether a particular region has at least one sample for this tracer
    real, allocatable               :: rif(:), rjf(:) ! interpolation coefficients calculated during the forward run
    integer, allocatable            :: ifr(:), jfr(:) ! gridbox indices for each observation
    integer, allocatable            :: ifn(:), jfn(:) ! gridbox indices for neighbors
    real, allocatable               :: wcx(:), wcy(:) ! x and y weighting factors for slopes interpolation
    real, allocatable               :: lat(:), lon(:) ! latitude and longitude of sampling
    integer, allocatable            :: sample_itau(:) ! sampling time
    integer, allocatable            :: itau_start(:), itau_end(:) ! time interval during which this sample should be counted
    real, allocatable               :: forcing(:,:) ! mixing ratio forcing of interest
    integer(2), allocatable         :: nsamples(:) ! number of samples accumulated
    real, allocatable               :: total_wt(:) ! time-weighted weight of each sample
    integer(2), allocatable         :: sampling_strategy(:) ! sampling strategy (1 for 4-hour average, 2 for instantaneous samples, etc.)
    integer                         :: nobs ! # of observations
    integer(2), allocatable         :: idates(:,:) ! only for temporary storage
  end type sat_tracer_forcing

  type sat_region_forcing
    logical                                 :: has_data ! whether a particular region has at least one sample
    type(sat_tracer_forcing), allocatable   :: io_tracer(:)
  end type sat_region_forcing

end module user_output_satellite_data
