!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module user_output_satellite

  ! module included to write columns along a satellite track
  ! called from user_output
  ! Sourish Basu, 2011

  ! The input data are expected to be in files called inputfile_YYYYMM.nc, which are
  ! in a folder called 'satellites' within the input folder. This can be changed in
  ! init_satellitedata. variable indir_satellite. Each inputfile_YYYYMM.nc needs the
  ! following fields:
  !
  ! Dimensions
  ! ==========
  !
  ! n_obs    :: The number of observations stored in the file
  ! n_levels :: Number of levels over which the prior profile and averaging kernel is specified
  ! n_time   :: Number of integers needed to specify time, currently 6
  !
  ! Variables
  ! =========
  !
  ! float32 latitude(n_obs)
  ! float32 longitude(n_obs)
  ! int16 cdate(n_obs, n_time)
  !     Lat, lon and time for the satellite samples
  ! float64 column_mixing(n_obs)
  ! float64 sigma_column_mixing(n_obs)
  !     Measured total column mixing ratio and its uncertainty
  ! int16 sampling_strategy(n_obs)
  !     2 :: Sample only within the dynamic timestep containing the measurement time. So if the
  !          dynamic timestep is reduced to 11.25 minutes due to CFL check, the measurement counts
  !          only over that time step.
  !     3 :: Sample throughout the window of ndyn_max seconds containing the measurement. So if
  !          the measurement is at 12:13 PM, then it is sampled at all dynamic timesteps between
  !          12:00 and 13:30 (if ndyn_max is 5400).
  !
  ! The input file is expected to contain averaging kernels and prior profiles as well. However, they
  ! are not read by the Fortran code. They are only used by the Python code that generates the
  ! adjoint forcings later. The output from this Fortran code are stored in files called
  ! sat-track_YYYYMM.nc4.
  !
  ! The non-obvious rc keys are:
  !
  ! output.satellite.sample.parent (T/F)   :: Whether to sample interface cells in the parent zoom region or not
  ! output.satellite.interpolate (T/F)     :: Whether to interpolate tracer mixing ratio to satellite coordinates or use the gridbox value
  ! output.satellite.errors                :: Set to 'neighbors' to estimate representation errors by evaluating neighboring grid cells
  ! output.satellite.meteo (T/F)           :: Whether to output meteo data such as temperature and geopotential height
  ! output.satellite.dryair.mixratio (T/F) :: Whether to correct airmass for moisture. If true, the airmass is assumed to be m*(1-q)
  ! output.satellite.split.period          :: Set to 'm' if there is one input file per month, 'd' if one per day (for large datasets such as CarbonSAT)

  use dims,                         only : nregions
  use toolbox,                      only : escape_tm
  use go,                           only : gol, goErr
  use orderpack,                    only : searchSorted, mrgrnk
  use user_output_satellite_data,   only : T_tracer_inout, T_output_fields, T_satellite_sample, T_satellite_region_counter
  use user_output_satellite_data,   only : SAT_INTERPOLATION_GRIDBOX, SAT_INTERPOLATION_SLOPES, SAT_INTERPOLATION_LINEAR
  use file_netcdf
  use datetime_for

  implicit none

  ! interface

  private

  public :: user_output_satellite_step
  public :: user_output_satellite_init
  public :: user_output_satellite_done

  logical                               :: period_has_data, output_satellite_verbose
  logical, allocatable, dimension(:)    :: period_opened ! keeps track of which month files have been already opened
  logical, allocatable, dimension(:)    :: period_exists ! keeps track of which months have satellite data
  logical, allocatable, dimension(:)    :: mask ! mask observations outside the model run time
  character(len=200)                    :: outdir_satellite, indir_rc, indir_satellite
  character(len=300)                    :: output_filename
  integer                               :: satellite_interpolation
  integer(2), allocatable               :: timlist(:,:)
  integer, allocatable                  :: region_rank(:)
  character(len=80)                     :: module_name = 'user_output_satellite'
  character(len=80)                     :: satellite_error_choice
  logical                               :: sample_in_parent = .false.
  integer                               :: parent_sample_choice
!  logical                               :: sample_dry_air = .false. ! this will replicate carbontracker sampling
  character(len=1)                      :: split_period ! 'd' if there is one file per day, 'm' if one per month
  logical, pointer                      :: new_period ! points to newday if there is one file per day, newmonth if there is one file per month

  type(T_tracer_inout), allocatable     :: IO_tracer(:)
  type(T_output_fields), allocatable    :: out_field(:)
  type(T_satellite_sample), allocatable :: sat_obs(:)
  type(T_satellite_region_counter), allocatable :: sat_reg_counter(:)

  character(len=*), parameter       ::  mname = 'user_output_satellite'

contains

subroutine user_output_satellite_init(status)

    use dims,        only: im, jm, lm, xbeg, xend, ybeg, yend, xref, yref
    use chem_param,  only: ntracet, names
    use dims,        only: nregions, newmonth, idatei, idatee, newday
    use datetime,    only: date2tau
    use omp_lib
    use GO,         only : ReadRc
    use global_data,only : rcF

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    integer              :: region, i, n_period, itr
    character(len=200)   :: fname
    logical              :: file_exist
    integer              :: tau_beg, tau_end, idate_temp(6)

    character(len=*), parameter         :: rname = mname//'/user_output_satellite_init'

    !__START_SUBROUTINE______________________________________________________

    call ReadRc(rcF, 'output.dir', outdir_satellite, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'input.dir', indir_rc, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.sample.parent', sample_in_parent, status, default=.true.)
    IF_ERROR_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.interpolation', satellite_interpolation, status, default=SAT_INTERPOLATION_LINEAR)
    IF_ERROR_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.errors', satellite_error_choice, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.verbose', output_satellite_verbose, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.split.period', split_period, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.output.directory', indir_satellite, status)
    IF_NOTOK_RETURN(status=1)

    parent_sample_choice = 0
    if (sample_in_parent) parent_sample_choice = 1

    ! rank the regions in inverse order of refinement (most refined first)
    allocate(region_rank(nregions))
    call mrgrnk(xref(1:nregions) * yref(1:nregions), region_rank)
    region_rank = region_rank(nregions:1:-1)
    ! allocate pos_counter and out_field
    ! during forward run, idatee > idatei
    select case (split_period)
        case ('m')
            n_period = (idatee(1)-idatei(1))*12 + (idatee(2)-idatei(2)+1)
            new_period => newmonth
        case ('d')
            idate_temp = idatei
            idate_temp(4) = 0
            idate_temp(5) = 0
            idate_temp(6) = 0
            call date2tau(idate_temp, tau_beg)
            idate_temp = idatee
            idate_temp(4) = 0
            idate_temp(5) = 0
            idate_temp(6) = 0
            call date2tau(idate_temp, tau_end)
            n_period = (tau_end-tau_beg)/SEC_PER_DAY + 1
            new_period => newday
        case default
            write(0,*) 'Wrong split period selected'
            IF_NOTOK_RETURN(status=1)
    end select

    allocate(period_opened(n_period))
    allocate(period_exists(n_period))
    do i=1,n_period
        period_opened(i) = .false.
        period_exists(i) = .true. ! assume data exist for all months, later set the empty months to false
    end do
    period_has_data = .false.

    allocate(IO_tracer(ntracet))
    do itr = 1, ntracet
        call ReadRc(rcF, 'output.satellite.meteo.'//trim(names(itr)), IO_tracer(itr)%output_meteo, status)
        IF_NOTOK_RETURN(status=1)
    end do

    status = 0

end subroutine user_output_satellite_init

subroutine user_output_satellite_done(status)

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter         :: rname = mname//'/user_output_satellite_done'

    !__START_SUBROUTINE______________________________________________________

    call close_satellitedatafile(status)
    IF_NOTOK_RETURN(status=1)

    ! deallocate some stuff
    if(allocated(IO_tracer)) deallocate(IO_tracer)


    status = 0

end subroutine user_output_satellite_done

subroutine user_output_satellite_step(region, tr, status)

    use dims,           only : im, jm, lm, idatei, newsrun
    use datetime,       only : tau2date, date2tau
    use datetime_for,   only : SEC_PER_DAY
    use Go,             only : TDate, Get, operator(+), operator(-), operator(/)

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    !__CONST________________________________________________________________

    character(len=*), parameter :: rname = mname//'/user_output_satellite_step'

    !__LOCAL_VARIABLES______________________________________________________

    character(len=200)          :: fname
    integer                     :: i_period, tau_beg, tau_mid
    integer, dimension(6)       :: midpt_date, idate_temp
    type(TDate)                 :: tmid

    !__START_SUBROUTINE______________________________________________________

    ! calculate the period index, to keep track of which periods have been read already
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call Get(tmid, time6=midpt_date)

    select case (split_period)
        case ('m')
            i_period = (midpt_date(1)-idatei(1))*12 + (midpt_date(2)-idatei(2)+1)
        case ('d')
            idate_temp = idatei
            idate_temp(4:6) = 0
            call date2tau(idate_temp, tau_beg)
            call date2tau(midpt_date, tau_mid)
            i_period = (tau_mid - tau_beg)/SEC_PER_DAY + 1
    end select

    ! If this is a new period but not a new run, we need to close the previous
    ! period's files and open the next period's files
    if ( .not. period_opened(i_period) ) then
        if ( .not. newsrun .and. period_has_data ) then
            call close_satellitedatafile(status)
            IF_NOTOK_RETURN(status=1)
        end if
        if (period_exists(i_period)) then
            call read_samples(midpt_date, i_period, status) ! read the input track for this month, sets 'period_has_data' according to whether a track file exists for this month or not
            IF_NOTOK_RETURN(status=1)
        end if
    end if

    call user_output_satellite_sample(region, tr, status)
    IF_NOTOK_RETURN(status=1)

    status = 0

end subroutine user_output_satellite_step

function in_region(lat,lon,region)

    use dims, only : dx, dy, xref, yref, xbeg, ybeg, xend, yend, parent

    implicit none

    real(4), intent(in) :: lat, lon
    integer, intent(in) :: region
    logical :: in_region

    real :: dxp, dyp

    in_region = .false.

    ! if region is the global region, then return true
    if (parent(region) == 0) then
        in_region = .true.
        return
    end if

    ! else, go through the exercise
    dxp = dx/xref(parent(region))
    dyp = dy/yref(parent(region))

    if ((lon .gt. xbeg(region)+dxp*parent_sample_choice) .and. (lon .lt. xend(region)-dxp*parent_sample_choice) &
        .and. (lat .gt. ybeg(region)+dyp*parent_sample_choice) .and. (lat .lt. yend(region)-dyp*parent_sample_choice)) in_region = .true.

end function in_region

subroutine read_samples(midpt_date, i_period, status)

    use dims,          only : nregions, lm, dx, dy, xref, yref, itaui, itaue
    use dims,          only : xcyc, xbeg, ybeg, xend, yend, im, jm, tref, ndyn, ndyn_max
    use dims,          only : isr, ier, jsr, jer
    use chem_param,    only : ntracet, names
    use datetime,      only : date2tau, tau2date
    ! For writing status messages
    use dims,          only : kmain, itau
    use datetime,      only : tstamp
    use Go,            only : goPr, gol

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in) :: midpt_date(6), i_period
    integer, intent(out) :: status

    !__LOCAL_VARIABLES______________________________________________________

    integer                 :: i, j, itrac, dims2(2), region, nc_id, grp_id, idate_temp(6), n_obs, n_instruments, i_instr
    character(len=256)      :: fname
    logical                 :: file_exist
    type(time_tm)           :: dummy_time
    integer, allocatable    :: region_index(:), region_counter(:), input_pos(:), sample_time(:), instr_serial(:), instr_counter(:)
    integer(2), allocatable :: sampling_strategy(:)
    integer(1), allocatable :: instrument(:), instrument_nums(:)
    real(4), allocatable    :: latitude(:), longitude(:)
    real, allocatable       :: observed_mixing(:), std_observed_mixing(:)
    real                    :: dxr, dyr
    type(time_tm)           :: sampleTime, startTime, endTime
    character(len=*), parameter :: rname = mname//'/read_samples'
    ! For writing status messages
    character(len=256)      :: write_string

    !__START_SUBROUTINE______________________________________________________

    write(gol,'(a, " :: entering")') rname ; call goPr

    ! Open inputfile
    select case (split_period)
        case ('m')
            write(fname,"(a,a,i4.4,i2.2,a)") trim(indir_satellite), '/inputfile_', midpt_date(1:2), '.nc4'
        case ('d')
            write(fname,"(a,a,i4.4,2i2.2,a)") trim(indir_satellite), '/inputfile_', midpt_date(1:3), '.nc4'
    end select

    inquire( file=fname, exist=file_exist )
    if ( file_exist ) then

        nc_id = nc_open(fname, 'r', status)
        IF_NOTOK_RETURN(status=1)

        do itrac = 1, ntracet
            IO_tracer(itrac)%name = names(itrac)
            IO_tracer(itrac)%has_data = nc_grp_exists(nc_id, names(itrac))
            ! Sometimes the tracer group may exist but have no observations
            if (IO_tracer(itrac)%has_data) then
                grp_id = nc_get_group(nc_id, names(itrac), status)
                IF_NOTOK_RETURN(status=1)
                n_obs  = nc_get_dim(grp_id, 'n_obs')
                IO_tracer(itrac)%has_data = IO_tracer(itrac)%has_data .and. (n_obs > 0)
            end if

            if (.not. IO_tracer(itrac)%has_data) then
                write(*, '("File ", a, " has no data for tracer ", a)') trim(fname), trim(names(itrac))
                cycle
            end if

            ! Read longitudes/latitudes and sample times into memory
            latitude = nc_read_var(grp_id, 'latitude', status)
            IF_NOTOK_RETURN(status=1)
            longitude = nc_read_var(grp_id, 'longitude', status)
            IF_NOTOK_RETURN(status=1)

            ! Read the instrument number
            instrument = nc_read_var(grp_id, 'instrument', status)
            IF_NOTOK_RETURN(status=1)

            allocate(input_pos(n_obs))
            input_pos(:) = (/ (i, i=1, n_obs) /)

            ! For each instrument, make a unique index that points to the corresponding entry in the instrument group. For this to
            ! work, the input file needs to have the instrument-specific variables (such as the averaging kernels) sorted in the
            ! same order as the common variables (such as lat, lon, time). That is, if we subselect all the common variables for
            ! a specific instrument, they should be in the same order as the instrument-specific variables for that instrument.
            ! So, e.g., if instrument = [1, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 2, 2, 3, 1],
            ! then       instr_serial = [1, 2, 1, 3, 1, 2, 2, 4, 5, 3, 3, 4, 5, 4, 6],
            !                                     |           |              |
            !                               first 3        fourth 1          fifth 2
            ! i.e., instr_serial for any given sounding is the serial number *for that specific instrument*
            allocate(instr_serial(n_obs))
            ! Instruments need not be serially numbered from 1. E.g., there could be three instruments, 11, 5 and 36.
            n_instruments = nc_get_dim(grp_id, 'instruments') ! number of instruments
            instrument_nums = nc_read_var(grp_id, 'instrument_nums') ! the integer ID of each instrument

            allocate(instr_counter(n_instruments))
            instr_counter(:) = 1
            do i = 1, n_obs
                ! Which instrument does this obs belong to?
                do i_instr = 1, n_instruments
                    if (instrument(i) == instrument_nums(i_instr)) exit
                end do
                if (instrument(i) /= instrument_nums(i_instr)) then
                    ! exited the loop because the loop variable was exhausted, not because a match was found
                    write(*, '(a, " :: ", i10, "-th observation is from instrument ", i3)') rname, i, instrument(i)
                    write(*, '(a, " :: this instrument does not exist in the variable instrument_nums in file ", a)') rname, trim(fname)
                    status = 1
                    IF_NOTOK_RETURN(status=1)
                end if
                ! i-th obs belongs to i_instr-th instrument

                instr_serial(i) = instr_counter(i_instr)
                instr_counter(i_instr) = instr_counter(i_instr) + 1
            end do
            deallocate(instr_counter, instrument, instrument_nums)

            sampling_strategy = nc_read_var(grp_id, 'sampling_strategy', status)
            IF_NOTOK_RETURN(status=1)
            timlist = nc_read_var(grp_id, 'cdate', status)
            IF_NOTOK_RETURN(status=1)
            allocate(sample_time(n_obs))
            ! convert the 6 integers to a single number, the internal time coordinate of TM5
            do i=1,n_obs
                call date2tau(timlist(:,i), sample_time(i))
            end do
            deallocate(timlist)

            ! find the regions
            allocate(region_index(n_obs))
            do i = 1, n_obs
                do j=1,nregions
                    region = region_rank(j)
                    if (in_region(latitude(i), longitude(i), region)) exit
                end do
                ! In case the lat/lon are invalid, set region = -1, which will be fitlered out later
                if ((abs(latitude(i)) > 90.0) .or. (abs(longitude(i)) > 180.0)) then
                    region = -1
                    write(*,'(a, " :: Satellite sample ", i6, " has an invalid latitude and/or longitude, skipped")') &
                        rname, i
                end if

                region_index(i) = region
            end do
            ! filter out observations outside simulation window
            allocate(mask(n_obs))
            mask = (sample_time .ge. itaui) .and. (sample_time .lt. itaue)  ! sample only itaui <= obs < itaue
            ! also filter out observations with invalid regions
            mask = mask .and. (region_index > 0)

            n_obs = count(mask)  ! note that this changes the value of n_obs
            if (output_satellite_verbose) write(*,'(a, " :: ", i8, " satellite observations to be assimilated")') rname, n_obs

            if (n_obs == 0) then
                if (allocated(latitude)) deallocate(latitude)
                if (allocated(longitude)) deallocate(longitude)
                !if (allocated(observed_mixing)) deallocate(observed_mixing)
                !if (allocated(std_observed_mixing)) deallocate(std_observed_mixing)
                if (allocated(sampling_strategy)) deallocate(sampling_strategy)
                if (allocated(sample_time)) deallocate(sample_time)
                if (allocated(region_index)) deallocate(region_index)
                if (allocated(mask)) deallocate(mask)
                if (allocated(input_pos)) deallocate(input_pos)
                if (allocated(instr_serial)) deallocate(instr_serial)
                IO_tracer(itrac)%has_data = .false.
            end if

            if (.not. IO_tracer(itrac)%has_data) cycle ! go to the next tracer

            if (allocated(IO_tracer(itrac)%sat_obs)) deallocate(IO_tracer(itrac)%sat_obs)
            if (allocated(IO_tracer(itrac)%sat_reg_counter)) deallocate(IO_tracer(itrac)%sat_reg_counter)

            IO_tracer(itrac)%n_obs = n_obs
            allocate(IO_tracer(itrac)%sat_obs(n_obs))
            do i=1,n_obs
                IO_tracer(itrac)%sat_obs(i)%nsamples    = 0
                IO_tracer(itrac)%sat_obs(i)%weight      = 0.0
                IO_tracer(itrac)%sat_obs(i)%mod_mix     = 0.0
                IO_tracer(itrac)%sat_obs(i)%mod_var     = 0.0
                !IO_tracer(itrac)%sat_obs(i)%obs_mix     = 0.0
                !IO_tracer(itrac)%sat_obs(i)%obs_std     = 0.0
                IO_tracer(itrac)%sat_obs(i)%evaluated   = .false.
                IO_tracer(itrac)%sat_obs(i)%pressure    = 0.0
            end do

            ! now apply the mask
            IO_tracer(itrac)%sat_obs%serial = pack(input_pos, mask)
            IO_tracer(itrac)%sat_obs%instr_serial = pack(instr_serial, mask)
            IO_tracer(itrac)%sat_obs%lat = pack(latitude, mask)
            IO_tracer(itrac)%sat_obs%lon = pack(longitude, mask)
            !IO_tracer(itrac)%sat_obs%obs_mix = pack(observed_mixing, mask)
            !IO_tracer(itrac)%sat_obs%obs_std = pack(std_observed_mixing, mask)
            IO_tracer(itrac)%sat_obs%sampling_strategy = pack(sampling_strategy, mask)
            IO_tracer(itrac)%sat_obs%sample_itau = pack(sample_time, mask)
            IO_tracer(itrac)%sat_obs%region = pack(region_index, mask)

            ! initialize with bogus values
            IO_tracer(itrac)%sat_obs%ifr = -1
            IO_tracer(itrac)%sat_obs%jfr = -1
            IO_tracer(itrac)%sat_obs%ifn = -1
            IO_tracer(itrac)%sat_obs%jfn = -1
            IO_tracer(itrac)%sat_obs%wcx = -1e12
            IO_tracer(itrac)%sat_obs%wcy = -1e12

            ! deallocate the already packed arrays
            !deallocate(input_pos, latitude, longitude, observed_mixing, std_observed_mixing, sampling_strategy, sample_time, region_index, mask)
            deallocate(input_pos, instr_serial, latitude, longitude, sampling_strategy, sample_time, region_index, mask)

            ! Now go through the samples and calculate region, ifr, etc.
            do i = 1, n_obs
                region = IO_tracer(itrac)%sat_obs(i)%region
                ! calculate itau_start and itau_end
                select case (IO_tracer(itrac)%sat_obs(i)%sampling_strategy)
                    case (2) ! instantaneous
                        IO_tracer(itrac)%sat_obs(i)%itau_start = IO_tracer(itrac)%sat_obs(i)%sample_itau
                        IO_tracer(itrac)%sat_obs(i)%itau_end   = IO_tracer(itrac)%sat_obs(i)%sample_itau

                    case (3) ! sample within ndyn
                        IO_tracer(itrac)%sat_obs(i)%itau_start = IO_tracer(itrac)%sat_obs(i)%sample_itau - &
                            mod(IO_tracer(itrac)%sat_obs(i)%sample_itau - itaui, ndyn_max/tref(region)) - 1 ! 1 second leeway on either side
                        IO_tracer(itrac)%sat_obs(i)%itau_end   = IO_tracer(itrac)%sat_obs(i)%itau_start  + &
                            ndyn_max/tref(region) + 1 ! 1 second leeway on either side

                        if (output_satellite_verbose) then
                            call tau2date(IO_tracer(itrac)%sat_obs(i)%itau_start, idate_temp)
                            startTime = idate_temp(1:6)
                            call tau2date(IO_tracer(itrac)%sat_obs(i)%itau_end, idate_temp)
                            endTime = idate_temp(1:6)
                            call tau2date(IO_tracer(itrac)%sat_obs(i)%sample_itau, idate_temp)
                            sampleTime = idate_temp(1:6)
                            write(*,'(a, " :: observation over ", f6.2, ",", f7.2, " at ", a, " will be sampled between ", a, " and ", a)') &
                                rname, IO_tracer(itrac)%sat_obs(i)%lat, IO_tracer(itrac)%sat_obs(i)%lon, &
                                trim(adjustl(format_datetime(sampleTime))), trim(adjustl(format_datetime(startTime))), trim(adjustl(format_datetime(endTime)))
                        end if

                    case default
                        if (output_satellite_verbose) then
                            call tau2date(IO_tracer(itrac)%sat_obs(i)%sample_itau, idate_temp)
                            sampleTime = idate_temp(1:6)
                            write(*,'(a, " :: unknown sampling strategy ", i2, " for observation at ", a, " over ", f6.2, ", ", f7.2)') &
                                rname, IO_tracer(itrac)%sat_obs(i)%sampling_strategy, trim(adjustl(format_datetime(sampleTime))), &
                                IO_tracer(itrac)%sat_obs(i)%lat, IO_tracer(itrac)%sat_obs(i)%lon
                        end if
                end select

                dyr = dy/yref(region)
                dxr = dx/xref(region)
                IO_tracer(itrac)%sat_obs(i)%rif = (IO_tracer(itrac)%sat_obs(i)%lon-float(xbeg(region)))/dxr + 0.99999
                IO_tracer(itrac)%sat_obs(i)%rjf = (IO_tracer(itrac)%sat_obs(i)%lat-float(ybeg(region)))/dyr + 0.99999
                IO_tracer(itrac)%sat_obs(i)%ifr  = int(IO_tracer(itrac)%sat_obs(i)%rif)   ! i-index of grid cell in which observation is located
                IO_tracer(itrac)%sat_obs(i)%jfr  = int(IO_tracer(itrac)%sat_obs(i)%rjf)   ! j-index of grid cell in which observation is located
                !fraction from the center of the is-box  (-0.5---+0.5)
                IO_tracer(itrac)%sat_obs(i)%rif = IO_tracer(itrac)%sat_obs(i)%rif-IO_tracer(itrac)%sat_obs(i)%ifr-0.5
                !idem js
                IO_tracer(itrac)%sat_obs(i)%rjf = IO_tracer(itrac)%sat_obs(i)%rjf-IO_tracer(itrac)%sat_obs(i)%jfr-0.5
                !the neighbors for pressure interpolation
                if(IO_tracer(itrac)%sat_obs(i)%rif .gt. 0) then
                    IO_tracer(itrac)%sat_obs(i)%ifn = IO_tracer(itrac)%sat_obs(i)%ifr+1
                else
                    IO_tracer(itrac)%sat_obs(i)%ifn = IO_tracer(itrac)%sat_obs(i)%ifr-1
                end if
                if(IO_tracer(itrac)%sat_obs(i)%rjf .gt. 0) then
                    IO_tracer(itrac)%sat_obs(i)%jfn = IO_tracer(itrac)%sat_obs(i)%jfr+1
                else
                    IO_tracer(itrac)%sat_obs(i)%jfn = IO_tracer(itrac)%sat_obs(i)%jfr-1
                end if
                ! x- / y-weighting of grid cell in which observation is located
                IO_tracer(itrac)%sat_obs(i)%wcx = (1.0-abs(IO_tracer(itrac)%sat_obs(i)%rif))    ! 1.0 ... 0.5
                IO_tracer(itrac)%sat_obs(i)%wcy = (1.0-abs(IO_tracer(itrac)%sat_obs(i)%rjf))    ! 1.0 ... 0.5
                !=================================================================
                ! if index of neighbour is exceeding range of region set
                ! neighbour = current cell (i.e. no interpolation)
                ! in case of cyclic x-boundaries take corresponding cyclic i index
                !=================================================================

                ! y-direction:
                ! ~ original, used in older 'T..' codes and pyshell:
                !if (IO_tracer(itrac)%sat_obs(i)%jfn < 1         ) IO_tracer(itrac)%sat_obs(i)%jfn = 1
                !if (IO_tracer(itrac)%sat_obs(i)%jfn > jm(region)) IO_tracer(itrac)%sat_obs(i)%jfn = jm(region)
                ! ~ changed to T38 equivalent:
                if (IO_tracer(itrac)%sat_obs(i)%jfn < jsr(region) ) IO_tracer(itrac)%sat_obs(i)%jfn = jsr(region)
                if (IO_tracer(itrac)%sat_obs(i)%jfn > jer(region) ) IO_tracer(itrac)%sat_obs(i)%jfn = jer(region)

                ! x-direction, check on cyclic:
                if ( xcyc(region) == 0 ) then
                    ! non-cyclic boundaries
                    ! ~ original, used in older 'T..' codes and pyshell:
!                    if (IO_tracer(itrac)%sat_obs(i)%ifn < 1          ) IO_tracer(itrac)%sat_obs(i)%ifn = 1
!                    if (IO_tracer(itrac)%sat_obs(i)%ifn > im(region) ) IO_tracer(itrac)%sat_obs(i)%ifn = im(region)
                    ! ~ changed to T38 equivalent:
                    if (IO_tracer(itrac)%sat_obs(i)%ifn < isr(region) ) IO_tracer(itrac)%sat_obs(i)%ifn = isr(region)
                    if (IO_tracer(itrac)%sat_obs(i)%ifn > ier(region) ) IO_tracer(itrac)%sat_obs(i)%ifn = ier(region)
                else
                    ! cyclic x-boundaries
                    if (IO_tracer(itrac)%sat_obs(i)%ifn < 1)          IO_tracer(itrac)%sat_obs(i)%ifn = im(region)
                    if (IO_tracer(itrac)%sat_obs(i)%ifn > im(region)) IO_tracer(itrac)%sat_obs(i)%ifn = 1
                end if

                ! allocate mod_profile and var_mod_profile
                allocate(IO_tracer(itrac)%sat_obs(i)%mod_profile(lm(region)))
                allocate(IO_tracer(itrac)%sat_obs(i)%var_mod_profile(lm(region)))
                IO_tracer(itrac)%sat_obs(i)%mod_profile = 0.0
                IO_tracer(itrac)%sat_obs(i)%var_mod_profile = 0.0

                if (IO_tracer(itrac)%output_meteo) then
                    allocate(IO_tracer(itrac)%sat_obs(i)%q(lm(region)))
                    allocate(IO_tracer(itrac)%sat_obs(i)%gph(lm(region)+1))
                    allocate(IO_tracer(itrac)%sat_obs(i)%t(lm(region)))
                    IO_tracer(itrac)%sat_obs(i)%q = 0.0
                    IO_tracer(itrac)%sat_obs(i)%gph = 0.0
                    IO_tracer(itrac)%sat_obs(i)%t = 0.0
                end if
            end do ! i = 1, n_obs

            allocate(IO_tracer(itrac)%sat_reg_counter(nregions))
            do region=1,nregions
                IO_tracer(itrac)%sat_reg_counter(region)%n_obs = count(IO_tracer(itrac)%sat_obs(:)%region == region)
                allocate(IO_tracer(itrac)%sat_reg_counter(region)%iobs(IO_tracer(itrac)%sat_reg_counter(region)%n_obs))
            end do ! region
            allocate(region_counter(nregions))
            region_counter = 1
            do i=1,n_obs
                region = IO_tracer(itrac)%sat_obs(i)%region
                j = region_counter(region)
                IO_tracer(itrac)%sat_reg_counter(region)%iobs(j) = i
                region_counter(region) = region_counter(region) + 1
            end do ! i
            deallocate(region_counter)

            if (sum(IO_tracer(itrac)%sat_reg_counter(:)%n_obs) /= n_obs) write(*,'(a,a,a)') 'WARNING :: ', rname, ' :: not all satellite samples accounted for'
            write(write_string, '("For tracer ", a, ", ", i6, " samples read in")') trim(names(itrac)), n_obs
            call tstamp(kmain, itau, trim(write_string))

        end do ! itrac

        call nc_close(nc_id, status)
        IF_NOTOK_RETURN(status=1)

        period_opened(i_period) = .true.
        period_exists(i_period) = .true.

    else
        if (output_satellite_verbose) write(*,'(a, " :: satellite input file ", a, " does not exist")') rname, trim(fname)
        period_exists(i_period) = .false.
        IO_tracer(:)%has_data = .false.
    end if ! file_exist

    ! This period has data only if at least one tracer has data for this period
    period_has_data = any(IO_tracer(:)%has_data)

    if (period_has_data) then
        select case (split_period)
            case ('m')
                write(output_filename,'(a,a,i4.4,i2.2,a4)') trim(outdir_satellite),'/satellite/sat-track_',midpt_date(1:2),'.nc4'
            case ('d')
                write(output_filename,'(a,a,i4.4,2i2.2,a4)') trim(outdir_satellite),'/satellite/sat-track_',midpt_date(1:3),'.nc4'
        end select
    end if ! period_has_data

    status = 0

    write(gol,'(a, " :: done")') rname ; call goPr

end subroutine read_samples

subroutine user_output_satellite_sample(region, tr, status)

    use chem_param,  only : ntracet, fscale
    use dims,        only : lm
    use global_data, only : mass_dat
    use MeteoData,   only : m_dat, humid_dat, gph_dat, temper_dat, phlb_dat
    use datetime,    only : date2tau, tau2date
    use Go,          only : TDate, Pretty, Get, NewDate, rTotal, operator(-)

    implicit none

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    character(len=*), parameter  ::  rname = mname//'/user_output_satellite_sample'

    real, dimension(:,:,:), pointer   :: m, gph, T, phlb, q
    real, dimension(:,:,:,:), pointer :: rm, rxm, rym
    ! SB: I've been told by MK that TM5 will never have region-dependent vertical resolution, so I'm going
    ! to declare some fixed-dimension arrays for working with profiles.
    integer, parameter :: lm_fixed = lm(1)
    real :: profile(-1:1,-1:1,lm_fixed), rmf(lm_fixed), rmf_p1(lm_fixed+1), mean_mix(lm_fixed), var_mix(lm_fixed)
    integer :: num_neighbors = 9
    real :: col_2d(2,2)
    ! End of SB hard-coding

    logical     :: in_window
    integer     :: i, j, k, i_obs, i_region, itrac
    integer     :: ifr,jfr,lfr,ifn,jfn
    !! Debug
    !integer     :: obs_in_window, t_great, t_less
    !type(TDate) :: sat_t
    !integer     :: sat_t_int(6)
    !! End debug
    real        :: wcx, wcy, rif, rjf, pres, mod_var, sum_rm, sum_m, weight
    integer     :: itau_tr(2), idate_temp(6)
    type(TDate) :: t_temp

    status = 0

    ! if this period does not have data, return
    if (.not. period_has_data) return

    weight = rTotal(tr(2)-tr(1), 'sec')

    ! Pointers to global arrays
    ! In the pointer assignments below, sometimes we need to explicitly specify the bounds on the right hand side. The issue
    ! is that for an assignment like
    !
    ! m => m_dat(region)%data
    !
    ! where m_dat(region)%data has bounds -1:im(region)+2, -1:jm(region)+2, m will have the same bounds. However, for an
    ! assignment like
    !
    ! p => phlb_dat(region)%data(:,:,1)
    !
    ! even though phlb_dat(region)%data has bounds -1:im(region)+2, -1:jm(region)+2, p, pointing to a slice of the pointee
    ! array, will not have that information and will have bounds 1:im(region)+4, 1:jm(region)+4. So for p we will do pointer
    ! bounds remapping, supported since Fortran 2003. Or, we could just use phlb, and point it to the entire array.

    m    => m_dat(region)%data
    rm   => mass_dat(region)%rm_t
    rxm  => mass_dat(region)%rxm_t
    rym  => mass_dat(region)%rym_t
    q    => humid_dat(region)%data
    phlb => phlb_dat(region)%data

    ! convert tr to integer
    call Get(tr(1), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(1))
    call Get(tr(2), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(2))

    do itrac = 1, ntracet

        !obs_in_window = -1

        if (.not. IO_tracer(itrac)%has_data) cycle

        if (IO_tracer(itrac)%output_meteo) then
            gph  => gph_dat(region)%data
            t    => temper_dat(region)%data
        end if

        !! Debug
        !obs_in_window = 0
        !t_great = 0
        !t_less = 0
        !call tau2date(maxval(IO_tracer(itrac)%sat_obs(:)%itau_end), sat_t_int)
        !sat_t = NewDate(time6=sat_t_int)
        !write(*,'(a, " :: max val of itau in IO_tracer = ", a, ", tr2 = ", a)') rname, trim(Pretty(sat_t)), trim(Pretty(tr(2)))
        !! End debug

        !$OMP PARALLEL PRIVATE (i, j, pres, profile, mean_mix, var_mix, mod_var, col_2d, sum_rm, sum_m) &
        !$OMP PRIVATE (rmf, rmf_p1, i_region, i_obs, in_window, idate_temp) &
        !$OMP PRIVATE (ifr, ifn, jfr, jfn, rif, rjf, wcx, wcy) &
        !$OMP REDUCTION (+:status)
        !$OMP DO schedule(static,10)
        do i_region = 1, IO_tracer(itrac)%sat_reg_counter(region)%n_obs
            i_obs = IO_tracer(itrac)%sat_reg_counter(region)%iobs(i_region)

            !  1. Is model time in the sampling window?
            in_window = .false.
            ! if the dynamic timestep is 12:30 -- 13:15, then itaur(region) points to 12:30
            select case (IO_tracer(itrac)%sat_obs(i_obs)%sampling_strategy)

                case (2) ! instantaneous
                    if (itau_tr(1) .le. IO_tracer(itrac)%sat_obs(i_obs)%sample_itau .and. &
                    itau_tr(2) .gt. IO_tracer(itrac)%sat_obs(i_obs)%sample_itau) in_window = .true.

                case(3) ! within ndyn/tref
                    if (itau_tr(1) .ge. IO_tracer(itrac)%sat_obs(i_obs)%itau_start .and. &
                    itau_tr(2) .le. IO_tracer(itrac)%sat_obs(i_obs)%itau_end) in_window = .true.

                    !if (itau_tr(1) .lt. IO_tracer(itrac)%sat_obs(i_obs)%itau_start) t_less = t_less + 1
                    !if (itau_tr(2) .gt. IO_tracer(itrac)%sat_obs(i_obs)%itau_end) t_great = t_great + 1

                case default
                    if (output_satellite_verbose) then
                        call tau2date(IO_tracer(itrac)%sat_obs(i_obs)%sample_itau, idate_temp)
                        t_temp = NewDate(time6=idate_temp)
                        write(gol,'(a, " :: unknown sampling strategy ", i2, " for observation at ", a, " over ", f6.2, ", ", f7.2, " :: skipping sample")') &
                            rname, IO_tracer(itrac)%sat_obs(i_obs)%sampling_strategy, trim(Pretty(t_temp)), &
                            IO_tracer(itrac)%sat_obs(i_obs)%lat, IO_tracer(itrac)%sat_obs(i_obs)%lon
                    else
                        write(gol,*) 'Unknown sampling strategy, skipping flask'
                    end if
                    call goErr

            end select

            if (in_window) then

                !obs_in_window = obs_in_window + 1

                !  2. Use slopes to determine tracer concentrations at the site.
                ifr = IO_tracer(itrac)%sat_obs(i_obs)%ifr
                ifn = IO_tracer(itrac)%sat_obs(i_obs)%ifn
                jfr = IO_tracer(itrac)%sat_obs(i_obs)%jfr
                jfn = IO_tracer(itrac)%sat_obs(i_obs)%jfn
                rif = IO_tracer(itrac)%sat_obs(i_obs)%rif
                rjf = IO_tracer(itrac)%sat_obs(i_obs)%rjf
                wcx = IO_tracer(itrac)%sat_obs(i_obs)%wcx
                wcy = IO_tracer(itrac)%sat_obs(i_obs)%wcy

                if (output_satellite_verbose) then
                    write(*,'("========== ", a, " ==========")') rname
                    write(*,'(a, " :: Sample over ",f6.2,",",f7.2," in region ",i1," at gridbox ",i3,",",i3," between ",a," and ",a)') &
                        rname, IO_tracer(itrac)%sat_obs(i_obs)%lat, IO_tracer(itrac)%sat_obs(i_obs)%lon, region, jfr, ifr, &
                        trim(Pretty(tr(1))), trim(Pretty(tr(2)))
                end if

                select case (satellite_interpolation)

                    case (SAT_INTERPOLATION_GRIDBOX)
                        rmf = rm(ifr,jfr,:,itrac) / m(ifr,jfr,:) *fscale(itrac)

                        sum_rm = sum(rm(ifr,jfr,:,itrac))
                        sum_m = sum(m(ifr,jfr,:))

                        ! sample surface pressure
                        pres = phlb(ifr,jfr,1)

                    case (SAT_INTERPOLATION_SLOPES)
                        rmf = ( rm(ifr,jfr,:,itrac) + 2.0*(rif*rxm(ifr,jfr,:,itrac) + rjf*rym(ifr,jfr,:,itrac)) ) &
                            / m(ifr,jfr,:) * fscale(itrac)

                        sum_rm = sum(rm(ifr,jfr,:,itrac) + 2.0*(rif*rxm(ifr,jfr,:,itrac) + rjf*rym(ifr,jfr,:,itrac)))
                        sum_m = sum(m(ifr,jfr,:))

                        ! sample surface pressure
                        pres = phlb(ifr,jfr,1)*wcx*wcy + phlb(ifn,jfr,1)*(1.0-wcx)*wcy + &
                            phlb(ifr,jfn,1)*wcx*(1.0-wcy) + phlb(ifn,jfn,1)*(1.0-wcx)*(1.0-wcy)

                    case (SAT_INTERPOLATION_LINEAR)
                        rmf = ( &
                             wcx  *      wcy  * rm(ifr,jfr,:,itrac) / m(ifr,jfr,:)  + &
                        (1.0-wcx) *      wcy  * rm(ifn,jfr,:,itrac) / m(ifn,jfr,:)  + &
                             wcx  * (1.0-wcy) * rm(ifr,jfn,:,itrac) / m(ifr,jfn,:)  + &
                        (1.0-wcx) * (1.0-wcy) * rm(ifn,jfn,:,itrac) / m(ifn,jfn,:)) * fscale(itrac)

                        sum_rm = sum(wcx*wcy*rm(ifr,jfr,:,itrac) + (1.0-wcx)*wcy*rm(ifn,jfr,:,itrac) + &
                            wcx*(1.0-wcy)*rm(ifr,jfn,:,itrac) + (1.0-wcx)*(1.0-wcy)*rm(ifn,jfn,:,itrac))
                        sum_m = sum(wcx*wcy*m(ifr,jfr,:) + (1.0-wcx)*wcy*m(ifn,jfr,:) + &
                            wcx*(1.0-wcy)*m(ifr,jfn,:) + (1.0-wcx)*(1.0-wcy)*m(ifn,jfn,:))

                        ! sample surface pressure
                        pres = phlb(ifr,jfr,1)*wcx*wcy + phlb(ifn,jfr,1)*(1.0-wcx)*wcy + &
                            phlb(ifr,jfn,1)*wcx*(1.0-wcy) + phlb(ifn,jfn,1)*(1.0-wcx)*(1.0-wcy)

                    case default
                        write (gol,'("unsupported satellite interpolation index ",i6)') satellite_interpolation; call goErr
                        status = status+1

                end select

                IO_tracer(itrac)%sat_obs(i_obs)%mod_profile = IO_tracer(itrac)%sat_obs(i_obs)%mod_profile + rmf * weight

                ! store the modeled total column, without any averaging kernel
                IO_tracer(itrac)%sat_obs(i_obs)%mod_mix = IO_tracer(itrac)%sat_obs(i_obs)%mod_mix + weight * sum_rm*fscale(itrac)/sum_m

                ! calculate the representation error
                select case (trim(adjustl(satellite_error_choice)))

                    case ('2d')
                        ! Calculate the modeled total columns at sampling point and three nearest neighbors, then
                        ! take the standard deviation across all four cells. This will not fill var_mix, because it
                        ! does not calculate the error in the profile, only the error in the total column.
                        col_2d(1,1) = sum(rm(ifr,jfr,:,itrac)) / sum(m(ifr,jfr,:)) * fscale(itrac)
                        col_2d(1,2) = sum(rm(ifn,jfr,:,itrac)) / sum(m(ifn,jfr,:)) * fscale(itrac)
                        col_2d(2,1) = sum(rm(ifr,jfn,:,itrac)) / sum(m(ifr,jfn,:)) * fscale(itrac)
                        col_2d(2,2) = sum(rm(ifn,jfn,:,itrac)) / sum(m(ifn,jfn,:)) * fscale(itrac)

                        ! The calculation for mod_var used to be
                        ! mod_var = sum(col_2d**2)/4.0 - (sum(col_2d)/4.0)**2
                        ! However, if the elements of col_2d are really close, the above can lead to a small negative number,
                        ! even though it is theoretically impossible. E.g., try the four numbers
                        ! col_2d = [395.300031257851, 395.300035334569, 395.300032213274, 395.300027276924]
                        ! With the above formula, this will yield mod_var = -2.9103830456733704e-11 in both Fortran and Python.
                        ! To get around this numerical issue, we will first subtract the mean, and then take the sum of squares.
                        ! We can reuse mod_var as the scalar which stores the mean AND later the sum of squares.
                        mod_var = sum(col_2d)/4.0
                        col_2d = col_2d - mod_var

                        mod_var = sum(col_2d*col_2d)/4.0
                        var_mix = 0.0

                    case ('neighbors')
                        ! Calculate the modeled profiles at the cell itself and three nearest neighbors using linear interpolation,
                        ! then take the variance across those cells. Use the 'profile' array for this, no need for a new array.
                        profile(0,0,:) = rm(ifr,jfr,:,itrac) / m(ifr,jfr,:) * fscale(itrac)
                        profile(0,1,:) = rm(ifn,jfr,:,itrac) / m(ifn,jfr,:) * fscale(itrac)
                        profile(1,0,:) = rm(ifr,jfn,:,itrac) / m(ifr,jfn,:) * fscale(itrac)
                        profile(1,1,:) = rm(ifn,jfn,:,itrac) / m(ifn,jfn,:) * fscale(itrac)

                        mean_mix = sum(sum(profile(0:1,0:1,:), 1), 1)/4.0
                        do i = 0, 1
                            do j = 0, 1
                                profile(i,j,:) = profile(i,j,:) - mean_mix
                            end do
                        end do

                        var_mix = sum(sum(profile(0:1,0:1,:)**2, 1), 1)/4.0
                        mod_var = 0.0

                    case ('gradient')
                        ! calculate the modeled profiles at the midpoints of all neighboring cells, then take the variance
                        ! in that profile across all nine neighbors
                        do i = -1, 1
                            do j = -1, 1
                                profile(i,j,:) = (rm(ifr,jfr,:,itrac) + &
                                    2.0*(i*rxm(ifr,jfr,:,itrac) + j*rym(ifr,jfr,:,itrac))) / &
                                    m(ifr,jfr,:) * fscale(itrac)
                            end do ! j
                        end do ! i
                        mean_mix = sum(sum(profile, 1), 1)/num_neighbors
                        do i = -1, 1
                            do j = -1, 1
                                profile(i,j,:) = profile(i,j,:) - mean_mix
                            end do ! j
                        end do ! i
                        var_mix = sum(sum(profile * profile, 1), 1)/num_neighbors
                        mod_var = 0.0

                    case default
                        write(0,'(a, " : Error algorithm ", a, " not known")') rname, trim(adjustl(satellite_error_choice))
                        status = 1

                end select

                IO_tracer(itrac)%sat_obs(i_obs)%var_mod_profile = IO_tracer(itrac)%sat_obs(i_obs)%var_mod_profile + var_mix * weight * weight

                IO_tracer(itrac)%sat_obs(i_obs)%mod_var = IO_tracer(itrac)%sat_obs(i_obs)%mod_var + mod_var * weight * weight

                IO_tracer(itrac)%sat_obs(i_obs)%pressure = IO_tracer(itrac)%sat_obs(i_obs)%pressure + pres * weight

                ! if output_meteo, sample gph, t and q
                if (IO_tracer(itrac)%output_meteo) then
                    rmf = wcx*wcy*t(ifr,jfr,:) + (1.0-wcx)*wcy*t(ifn,jfr,:) + wcx*(1.0-wcy)*t(ifr,jfn,:) + (1.0-wcx)*(1.0-wcy)*t(ifn,jfn,:)
                    IO_tracer(itrac)%sat_obs(i_obs)%t = IO_tracer(itrac)%sat_obs(i_obs)%t + rmf * weight

                    rmf_p1 = wcx*wcy*gph(ifr,jfr,:) + (1.0-wcx)*wcy*gph(ifn,jfr,:) + wcx*(1.0-wcy)*gph(ifr,jfn,:) + (1.0-wcx)*(1.0-wcy)*gph(ifn,jfn,:)
                    IO_tracer(itrac)%sat_obs(i_obs)%gph = IO_tracer(itrac)%sat_obs(i_obs)%gph + rmf_p1 * weight

                    rmf = wcx*wcy*q(ifr,jfr,:) + (1.0-wcx)*wcy*q(ifn,jfr,:) + wcx*(1.0-wcy)*q(ifr,jfn,:) + (1.0-wcx)*(1.0-wcy)*q(ifn,jfn,:)
                    IO_tracer(itrac)%sat_obs(i_obs)%q = IO_tracer(itrac)%sat_obs(i_obs)%q + rmf * weight
                end if

                IO_tracer(itrac)%sat_obs(i_obs)%nsamples = IO_tracer(itrac)%sat_obs(i_obs)%nsamples + 1
                IO_tracer(itrac)%sat_obs(i_obs)%weight   = IO_tracer(itrac)%sat_obs(i_obs)%weight   + weight

            end if ! in_window

        end do ! i_region
        !$OMP END DO
        !$OMP END PARALLEL

        IF_NOTOK_RETURN(status=1)

        if (IO_tracer(itrac)%output_meteo) nullify(t, gph)

        !write(*,'(a, " :: for tracer ", i1, ", obs between ", a, " and ", a, " = ", i6, ", t_great = ", i6, ", t_less = ", i6)') &
            !rname, itrac, trim(Pretty(tr(1))), trim(Pretty(tr(2))), obs_in_window, t_great, t_less
    end do ! itrac

    nullify(m, phlb, q)
    nullify(rm, rxm, rym)

end subroutine user_output_satellite_sample

subroutine close_satellitedatafile(status)

    use dims,     only : lm, region_name, at, bt
    use chem_param, only : ntracet, names
    use datetime, only : tau2date
    use misctools, only : check_dir

    implicit none

    integer, intent(out)            :: status

    character(len=*), parameter     :: rname = mname//'/close_satellitedatafile'

    !__LOCAL_VARIABLES______________________________________________________
    integer                 :: group_id, tgrp_id, region, i, itrac, met_group_id, output_fid, nobs
    integer, allocatable    :: i_obs(:)
    real                    :: weight

    if (period_has_data) then
        do itrac = 1, ntracet
            if (.not. IO_tracer(itrac)%has_data) cycle
            ! first calculate the mean observations
            do i=1,IO_tracer(itrac)%n_obs
                if((.not. IO_tracer(itrac)%sat_obs(i)%evaluated) .and. (IO_tracer(itrac)%sat_obs(i)%nsamples .ge. 1)) then
                    weight = IO_tracer(itrac)%sat_obs(i)%weight

                    IO_tracer(itrac)%sat_obs(i)%mod_profile = IO_tracer(itrac)%sat_obs(i)%mod_profile / weight
                    IO_tracer(itrac)%sat_obs(i)%var_mod_profile = IO_tracer(itrac)%sat_obs(i)%var_mod_profile / weight / weight
                    IO_tracer(itrac)%sat_obs(i)%mod_var = IO_tracer(itrac)%sat_obs(i)%mod_var / weight / weight
                    IO_tracer(itrac)%sat_obs(i)%mod_mix = IO_tracer(itrac)%sat_obs(i)%mod_mix / weight
                    IO_tracer(itrac)%sat_obs(i)%pressure = IO_tracer(itrac)%sat_obs(i)%pressure / weight
                    if (IO_tracer(itrac)%output_meteo) then
                        IO_tracer(itrac)%sat_obs(i)%t = IO_tracer(itrac)%sat_obs(i)%t / weight
                        IO_tracer(itrac)%sat_obs(i)%gph = IO_tracer(itrac)%sat_obs(i)%gph / weight
                        IO_tracer(itrac)%sat_obs(i)%q = IO_tracer(itrac)%sat_obs(i)%q / weight
                    end if
                    IO_tracer(itrac)%sat_obs(i)%evaluated = .true.
                end if ! evaluated, nsamples
            end do ! i


            ! segregate according to region
            allocate(IO_tracer(itrac)%out_field(nregions))
            allocate(i_obs(nregions))
            do region = 1, nregions
                !write(*,'(a, " :: tracer ", a, " in region ", a, " has region = ", i6, " and ", i6, " evaluated, with total nsamples ", i6)') &
                    !rname, region_name(region), names(itrac), count(IO_tracer(itrac)%sat_obs(:)%region == region), &
                    !count(IO_tracer(itrac)%sat_obs(:)%evaluated), sum(IO_tracer(itrac)%sat_obs(:)%nsamples)
                nobs = count(IO_tracer(itrac)%sat_obs(:)%region == region .and. IO_tracer(itrac)%sat_obs(:)%evaluated)
                IO_tracer(itrac)%out_field(region)%nobs = nobs
                allocate(IO_tracer(itrac)%out_field(region)%sampled_profiles(lm(region),nobs))
                allocate(IO_tracer(itrac)%out_field(region)%std_sampled_profiles(lm(region),nobs))
                allocate(IO_tracer(itrac)%out_field(region)%sampled_psurf(nobs))
                allocate(IO_tracer(itrac)%out_field(region)%modeled_mixing(nobs))
                allocate(IO_tracer(itrac)%out_field(region)%std_modeled_mixing(nobs))
                !allocate(IO_tracer(itrac)%out_field(region)%observed_mixing(nobs))
                !allocate(IO_tracer(itrac)%out_field(region)%std_observed_mixing(nobs))
                allocate(IO_tracer(itrac)%out_field(region)%input_positions(nobs))
                allocate(IO_tracer(itrac)%out_field(region)%instr_positions(nobs))
                allocate(IO_tracer(itrac)%out_field(region)%nsamples(nobs))
                allocate(IO_tracer(itrac)%out_field(region)%weight(nobs))
                allocate(IO_tracer(itrac)%out_field(region)%sampling_strategy(nobs))
                if (IO_tracer(itrac)%output_meteo) then
                    allocate(IO_tracer(itrac)%out_field(region)%t(lm(region),nobs))
                    allocate(IO_tracer(itrac)%out_field(region)%gph(lm(region)+1,nobs))
                    allocate(IO_tracer(itrac)%out_field(region)%q(lm(region),nobs))
                end if
            end do ! region

            i_obs(:) = 1
            do i = 1, IO_tracer(itrac)%n_obs
                if (IO_tracer(itrac)%sat_obs(i)%evaluated) then
                    region = IO_tracer(itrac)%sat_obs(i)%region

                    IO_tracer(itrac)%out_field(region)%sampled_profiles(:,i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%mod_profile
                    IO_tracer(itrac)%out_field(region)%std_sampled_profiles(:,i_obs(region)) = sqrt(IO_tracer(itrac)%sat_obs(i)%var_mod_profile)
                    IO_tracer(itrac)%out_field(region)%sampled_psurf(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%pressure
                    IO_tracer(itrac)%out_field(region)%modeled_mixing(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%mod_mix
                    IO_tracer(itrac)%out_field(region)%std_modeled_mixing(i_obs(region)) = sqrt(IO_tracer(itrac)%sat_obs(i)%mod_var)
                    !IO_tracer(itrac)%out_field(region)%observed_mixing(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%obs_mix
                    !IO_tracer(itrac)%out_field(region)%std_observed_mixing(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%obs_std
                    IO_tracer(itrac)%out_field(region)%input_positions(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%serial
                    IO_tracer(itrac)%out_field(region)%instr_positions(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%instr_serial
                    IO_tracer(itrac)%out_field(region)%nsamples(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%nsamples
                    IO_tracer(itrac)%out_field(region)%weight(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%weight
                    IO_tracer(itrac)%out_field(region)%sampling_strategy(i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%sampling_strategy
                    if (IO_tracer(itrac)%output_meteo) then
                        IO_tracer(itrac)%out_field(region)%t(:,i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%t
                        IO_tracer(itrac)%out_field(region)%gph(:,i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%gph
                        IO_tracer(itrac)%out_field(region)%q(:,i_obs(region)) = IO_tracer(itrac)%sat_obs(i)%q
                    end if
                    i_obs(region) = i_obs(region) + 1
                end if ! IO_tracer(itrac)%sat_obs(i)%evaluated
            end do ! i

            do region=1,nregions
                if (i_obs(region)-IO_tracer(itrac)%out_field(region)%nobs /= 1) then
                    write(*,'(a,i1)') 'There is a problem in counting satellite samples in region ', region
                    IO_tracer(itrac)%out_field(region)%nobs = i_obs(region) - 1
                end if
            end do
            deallocate(i_obs)
        end do ! itrac

        ! now write the data
        call check_dir(output_filename)

        output_fid = nc_open(output_filename, 'c', status)
        IF_NOTOK_RETURN(status=1)

        call nc_create_dim(output_fid, 'date_components', 6)

        do region=1,nregions
            group_id = nc_create_group(output_fid, region_name(region))
            call nc_create_dim(group_id, 'n_lev', lm(region))
            if (any(IO_tracer(:)%output_meteo)) call nc_create_dim(group_id, 'n_lev_p1', lm(region)+1)
            do itrac = 1, ntracet
                ! If this tracer has no data, either globally or for this region, cycle to the next tracer
                if (.not. IO_tracer(itrac)%has_data) then
                    write(*, '(a, " :: tracer ", a, " has no data in region ", a)') rname, trim(names(itrac)), trim(region_name(region))
                    cycle
                end if
                if (IO_tracer(itrac)%out_field(region)%nobs == 0) then
                    write(*, '(a, " :: tracer ", a, " has nobs = 0 in region ", a)') rname, trim(names(itrac)), trim(region_name(region))
                    cycle
                end if
                tgrp_id = nc_create_group(group_id, names(itrac))

                call nc_create_dim(tgrp_id, 'n_obs', IO_tracer(itrac)%out_field(region)%nobs)
                nobs = IO_tracer(itrac)%out_field(region)%nobs

                call nc_dump_var(tgrp_id, 'input_positions', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%input_positions(1:nobs), &
                    (/'long_name'/), (/'position of observation in input file'/))
                call nc_dump_var(tgrp_id, 'instrument_positions', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%instr_positions(1:nobs), &
                    (/'long_name'/), (/'position of observation for this instrument (identical to input_positions for a single instrument)'/))
                call nc_dump_var(tgrp_id, 'profiles',  (/'n_lev', 'n_obs'/), IO_tracer(itrac)%out_field(region)%sampled_profiles(1:lm(region),1:nobs), &
                    (/'long_name'/), (/'sampled model profile for each satellite observation'/))
                call nc_dump_var(tgrp_id, 'std_profiles',  (/'n_lev', 'n_obs'/), IO_tracer(itrac)%out_field(region)%std_sampled_profiles(1:lm(region),1:nobs), &
                    (/'long_name'/), (/'standard deviation of sampled model profile'/))
                call nc_dump_var(tgrp_id, 'psurf', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%sampled_psurf(1:nobs), &
                    (/'long_name'/), (/'sampled surface pressures'/))
                !call nc_dump_var(tgrp_id, 'column_mixing', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%observed_mixing(1:nobs), &
                    !(/'long_name'/), (/'observed total column mixing ratio'/))
                !call nc_dump_var(tgrp_id, 'sigma_column_mixing', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%std_observed_mixing(1:nobs), &
                    !(/'long_name'/), (/'error in observed total column mixing ratio'/))
                call nc_dump_var(tgrp_id, 'model_column', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%modeled_mixing(1:nobs), &
                    (/'long_name'/), (/'modeled total column mixing ratio (pressure weighted sum of profile, no averaging kernel)'/))
                call nc_dump_var(tgrp_id, 'sigma_model_column', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%std_modeled_mixing(1:nobs), &
                    (/'long_name'/), (/'error in modeled total column mixing ratio'/))
                call nc_dump_var(tgrp_id, 'nsamples', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%nsamples(1:nobs), &
                    (/'long_name'/), (/'samples averaged to get modeled profile'/))
                call nc_dump_var(tgrp_id, 'total_weight', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%weight(1:nobs), &
                    (/'long_name'/), (/'sum of weights of all samples, weights proportional to dynamic time step length'/))
                call nc_dump_var(tgrp_id, 'sampling_strategy', (/'n_obs'/), IO_tracer(itrac)%out_field(region)%sampling_strategy(1:nobs), &
                    (/'long_name'/), (/'sampling strategy for satellite observation'/))

                if (IO_tracer(itrac)%output_meteo) then
                    met_group_id = nc_create_group(tgrp_id, 'meteo')
                    call nc_dump_var(met_group_id, 'temperature', (/'n_lev','n_obs'/), IO_tracer(itrac)%out_field(region)%t(1:lm(region),1:nobs), &
                        (/'long_name'/), (/'temperature of grid-box'/))
                    call nc_dump_var(met_group_id, 'geopotential_height', (/'n_lev_p1','n_obs   '/), IO_tracer(itrac)%out_field(region)%gph(1:lm(region)+1,1:nobs), &
                        (/'long_name'/), (/'geopotential height at layer boundaries'/))
                    call nc_dump_var(met_group_id, 'specific_humidity', (/'n_lev   ','n_obs   '/), IO_tracer(itrac)%out_field(region)%q(1:lm(region),1:nobs), &
                        (/'long_name'/), (/'specific humidity (grams/gram)'/))
                    call nc_dump_var(met_group_id, 'at', (/'n_lev_p1'/), at, (/'comment'/), (/'AT coefficient in Pa (P = AT + BT*Psurf)'/))
                    call nc_dump_var(met_group_id, 'bt', (/'n_lev_p1'/), bt, (/'comment'/), (/'BT coefficient (P = AT + BT*Psurf)'/))
                end if
            end do ! itrac
        end do ! region

        call nc_close(output_fid)

        do itrac = 1, ntracet
            if (allocated(IO_tracer(itrac)%out_field))          deallocate(IO_tracer(itrac)%out_field)
            if (allocated(IO_tracer(itrac)%sat_obs))            deallocate(IO_tracer(itrac)%sat_obs)
            if (allocated(IO_tracer(itrac)%sat_reg_counter))    deallocate(IO_tracer(itrac)%sat_reg_counter)
        end do

    end if ! period_has_data

    status = 0

end subroutine close_satellitedatafile

end module user_output_satellite
