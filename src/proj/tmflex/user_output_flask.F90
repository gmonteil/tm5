!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module user_output_flask
! !PUBLIC TYPES: none
!
! !PUBLIC MEMBER FUNCTIONS:
!    user_output_flask_init
!    user_output_flask_sample
!    user_output_flask_done
!
! !DESCRIPTION:
!
!   Routine for simulating flask samples in TM5.
!
!   "Flask" is shorthand for any x,y,z,t sampling of the model.  This
!   routine is intended for producing model output of low-frequency,
!   intermittent, or sporadic sampling (like NOAA Cooperative Sampling
!   Network flasks), but is also useful for other observations like
!   aircraft profiles or daily averages from continuous station data.
!
!   I/O is via netCDF files, using the MDF interface.  Observation
!   list from input file is truncated to begin/end times of the current run.
!
!   Input file format:
!     Routine expects to read a netCDF (v3) file with:
!
!     "id":  record dimension, contains unique integers
!     "id":  integer ID for each observation, for tracing back
!     "lat": 1-d real array, length id, of latitudes (degrees N)
!     "lon": 1-d real array, length id, of longitudes (degrees E)
!     "alt": 1-d real array, length id, of altitudes (m ASL)
!
!     "date_components": 2-D variable with dims (id,6),
!        with time of observations in 6-integer 'idate' format
!        (year, month, day, hour, minute, second in UTC time.
!
!     "sampling_strategy" 1-D integer array designed to represent
!        distinct schemes for sampling TM5 to represent the
!        observations in this file.  Add to the following list if
!        you define a new scheme:
!
!        VALUE                  MEANING
!
!          1     Four-hour averages, centered on observation time.
!          2     Instantaneous values, no temporal averaging
!
!     "station_id" 1-D integer array representing which station an
!        observation came from. For mobile sampling, such as from an
!        aircraft, this number can be random, or one per aircraft.
!
!   Controlling rc-file keys:
!
!     output.point  (logical, default F)
!         Whether or not this routine is active.  Note that
!         unlike other rc-file keys, this one is read in
!         user_output.F90.
!
!     output.point.infile  (string, required, no default value)
!         Path to input netCDF file, e.g. "/path/to/input/flask_obs.nc"
!
!     output.point.verbose (logical, optional, default F)
!         If true, some extra information is printed to stdout.
!
!     output.point.meteo (logical, optional, default F)
!         If true, output meteo variables averaged at flask sites.
!         These variables are u, v, blh, q, press, and temp.
!
!         (1) obs are in-window for this simulation if itaui <= obs < itaue,
!             where itaui and itaue are the initial and ending simulation
!             time values, and obs is the observation center time.  The
!             default behavior of this routine is instead to sample all
!             obs meeting the criterion itaui-window/2 <= obs <= itaue+window/2.
!             This will result in obs being sampled in multiple simulations
!             when their sampling windows cross itaui or itaue.  This
!             must be dealt with in post-processing.
!
!
! !REVISION HISTORY:
!
!   Sourish Basu, September 2014 : Store interpolation coefficients for using pressure as the vertical coordinate
!
!   Sourish Basu, April 2014 : Break up input and output into monthly/daily chunks
!
!   Sourish Basu, Oct 2013 : Modifying for multiple tracers
!
!   Sourish Basu, Sep 2011 (I think) : Modified for PyShell 4DVAR
!
!   Andy Jacobson, Aug 2010
!
!   Adapted from code in user_output_forecast, user_output_station,
!   and user_output_noaa, mostly if not all originally written by
!   Wouter Peters.  Also used Arjo Segers' file_MDF example code
!   for netCDF I/O.
!
!
!EOP
!-------------------------------------------------------------------------

  use go,                       only : gol, goErr

  use user_output_flask_data,   only : VERT_COORD_ALT, VERT_COORD_PRES
  use user_output_flask_data,   only : POINT_INTERPOLATION_GRIDBOX, POINT_INTERPOLATION_SLOPES, POINT_INTERPOLATION_LINEAR
  use user_output_flask_data,   only : tracer_inout

  implicit none

  private

  public                :: user_output_flask_init
  public                :: user_output_flask_done
  public                :: user_output_flask_step

  character(len=512)    :: indir_point, outfile_point, outdir_point
  character(len=1)      :: split_period

  logical               :: flask_sample_meteo = .false.
  logical               :: flask_verbose = .false.
  integer               :: flask_interpolation
!  logical               :: flask_sample_interpolate = .true.
  logical               :: sample_in_parent = .false.
  logical               :: flask_active = .true. ! if there's no flask data, set this to false

  integer               :: parent_sample_choice
  integer               :: assim_window
  real,parameter        :: flask_missing_value=-1.0e34
  character(len=80)     :: flask_error_choice
  logical               :: period_has_data

  character(len=*), parameter       ::  mname = 'user_output_flask'

  type(tracer_inout), dimension(:), allocatable :: IO_tracer

  integer, allocatable                  :: region_rank(:)
  logical, allocatable, dimension(:)    :: period_opened ! keeps track of which period files have been already opened
  logical, allocatable, dimension(:)    :: period_exists ! keeps track of which periods have point data

contains

  subroutine user_output_flask_init(status)

    use GO,             only : ReadRc
    use global_data,    only : rcF
    use orderpack,      only : mrgrnk
    use datetime,       only : date2tau
    use dims,           only : nregions, idatei, idatee, xref, yref
    use datetime_for,   only : SEC_PER_DAY

    implicit none

    ! --- in/out ---------------------------------
    integer, intent(out)    :: status

    ! local
    integer                 :: n_period, tau_beg, tau_end, idate_temp(6), i

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/user_output_flask_init'

    write(*,'(a, " :: entering")') rname

    call ReadRc( rcF, 'output.point.input.dir', indir_point, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.point.verbose', flask_verbose, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'output.point.meteo', flask_sample_meteo, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'output.point.errors', flask_error_choice, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.point.interpolation', flask_interpolation, status, default=POINT_INTERPOLATION_SLOPES)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'output.point.timewindow', assim_window, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.point.sample.parent', sample_in_parent, status, default=.true.)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'output.point.split.period', split_period, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.dir', outdir_point, status)
    IF_NOTOK_RETURN(status=1)

    parent_sample_choice = 0
    if (sample_in_parent) parent_sample_choice = 1

    if(flask_verbose) then
        write (*,'(a, " :: verbose output requested.")') rname
    end if

    ! rank the regions in inverse order of refinement (most refined first)
    allocate(region_rank(nregions))
    call mrgrnk(xref(1:nregions) * yref(1:nregions), region_rank)
    region_rank = region_rank(nregions:1:-1)
    ! allocate pos_counter and out_field
    ! during forward run, idatee > idatei
    select case (split_period)
        case ('a') ! the entire period in one go
            n_period = 1
        case ('m')
            n_period = (idatee(1)-idatei(1))*12 + (idatee(2)-idatei(2)+1)
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

    status = 0

    write(*,'(a, " :: done")') rname

  end subroutine user_output_flask_init

  subroutine user_output_flask_done(status)

    implicit none

    !__IO___________________________________________________________________
    integer, intent(out)                :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter         :: rname = mname//'/user_output_flask_done'

    !__START_SUBROUTINE______________________________________________________

    call user_output_flask_write

    ! deallocate some stuff
    if(allocated(IO_tracer)) deallocate(IO_tracer)

    status = 0

  end subroutine user_output_flask_done

  function in_region(lat,lon,region)

    use dims, only : dx, dy, xref, yref, xbeg, ybeg, xend, yend, parent

    implicit none

    real, intent(in) :: lat, lon
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
        .and. (lat .gt. ybeg(region)+dyp*parent_sample_choice) .and. (lat .lt. yend(region)-dyp*parent_sample_choice)) &
            in_region = .true.

  end function in_region

  subroutine user_output_flask_step(region, tr, status)

    use dims,           only : idatei, newsrun
    use datetime,       only : tau2date, date2tau
    use datetime_for,   only : SEC_PER_DAY
    use Go,             only : TDate, Get, operator(+), operator(-), operator(/)

    implicit none

    !__IO___________________________________________________________________
    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________
    character(len=*), parameter :: rname = mname//'/user_output_flask_step'

    integer                 :: midpt_date(6), i_period, idate_temp(6), tau_beg, tau_mid
    type(TDate)             :: tmid

    !__START_SUBROUTINE______________________________________________________

    ! calculate the period index, to keep track of which periods have been read already
    ! Which period number are we in?
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call Get(tmid, time6=midpt_date)

    select case (split_period)
        case ('a')
            i_period = 1
        case ('m')
            i_period = (midpt_date(1)-idatei(1))*12 + (midpt_date(2)-idatei(2)+1)
        case ('d')
            idate_temp = idatei
            idate_temp(4:6) = 0
            call date2tau(midpt_date, tau_mid)
            call date2tau(idate_temp, tau_beg)
            i_period = (tau_mid - tau_beg)/SEC_PER_DAY + 1
    end select

    ! If this is a new period but not a new run, we need to close the previous
    ! period's files and open the next period's files
    if ( .not. period_opened(i_period) ) then
        if ( .not. newsrun .and. flask_active ) call user_output_flask_write
        if (period_exists(i_period)) then
            call read_samples(midpt_date, i_period, status) ! read the input track for this month, sets 'period_has_data' according to whether a track file exists for this month or not
            IF_NOTOK_RETURN(status=1)
        end if
    end if

    call user_output_flask_sample(region, tr, status)
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine user_output_flask_step

  subroutine read_samples(midpt_date, i_period, status)

    use chem_param,     only : ntracet, names
    use file_netcdf
    use Go,             only : ReadRc
    use global_data,    only : rcF
    use go_date,        only : TDate, Pretty, NewDate
    use dims,           only : nregions, itaui, itaue, ndyn_max, tref, region_name, xcyc
    use dims,           only : idatei, idatee, xref, yref, dx, dy, xbeg, ybeg, im, jm
    use dims,           only : isr, ier, jsr, jer
    use datetime,       only : tau2date, date2tau

    implicit none

    !__IO___________________________________________________________________
    integer, intent(in)     :: midpt_date(6), i_period
    integer, intent(out)    :: status

    character(len=*), parameter :: rname = mname//'/read_samples'

    !__LOCAL_VARIABLES______________________________________________________
    integer, allocatable    :: itau_center(:), itau_start(:), itau_end(:), region_list(:)
    integer, allocatable    :: station_id(:), id(:), region_counter(:), time_window_length(:)
    logical, allocatable    :: mask(:)
    real, allocatable       :: lat(:), lon(:), alt(:), pres(:)
    integer(2), allocatable :: idate_f(:,:), sampling_strategy(:), vert_coord_type(:)

    character(len=256)      :: fname
    character(len=16)       :: idstr
    logical                 :: flask_file_exists
    integer                 :: hid, itrac, nflasks, iflask, region, i_region, idate_temp(6)
    integer                 :: grp_id, idatef(6)
    type(TDate)             :: windowBeg, windowEnd, windowMid
    real                    :: dxr, dyr

    !__START_SUBROUTINE______________________________________________________

    write(*,'(a, " :: entering")') rname

    flask_active = .true.

    ! Open inputfile
    select case (split_period)
        case ('a')
            write(fname,"(a,a)") trim(indir_point), '/point_input.nc4'
            write(outfile_point,"(a,a)") trim(outdir_point), '/point/point_output.nc4'
        case ('m')
            write(fname,"(a,a,i4.4,i2.2,a)") trim(indir_point), '/point_input_', midpt_date(1:2), '.nc4'
            write(outfile_point,"(a,a,i4.4,i2.2,a)") trim(outdir_point), '/point/point_output_', midpt_date(1:2), '.nc4'
        case ('d')
            write(fname,"(a,a,i4.4,2i2.2,a)") trim(indir_point), '/point_input_', midpt_date(1:3), '.nc4'
            write(outfile_point,"(a,a,i4.4,2i2.2,a)") trim(outdir_point), '/point/point_output_', midpt_date(1:3), '.nc4'
    end select

    inquire(file=trim(fname), exist=flask_file_exists)

    if (flask_file_exists) then

        if (flask_verbose) then
            write (*,'(a, " :: input from ",a,".")') rname, trim(fname)
            write (*,'(a, " :: output to ",a,".")') rname, trim(outfile_point)
        end if

        ! read input file

        ! get from nc input file:
        !    lat, lon, alt, date_components, sampling_strategy, station_id
        !   restrict to flasks between itaui and itaue
        !   for each flask, determine region, ifr, jfr, ...things that don't change.

        allocate(IO_tracer(ntracet))

        hid = nc_open(trim(fname), 'r', status)
        IF_NOTOK_RETURN(status=1)

        do itrac = 1, ntracet

            IO_tracer(itrac)%name = names(itrac)

            call ReadRc( rcF, 'output.point.'//trim(names(itrac))//'.minerror', IO_tracer(itrac)%flask_minerror, status)
            IF_NOTOK_RETURN(status=1)

            IO_tracer(itrac)%flask_active = nc_grp_exists(hid, names(itrac))
            if (.not. IO_tracer(itrac)%flask_active) cycle

            grp_id = nc_get_group(hid, names(itrac), status)
            IF_NOTOK_RETURN(status=1)
            nflasks = nc_get_dim(grp_id, 'id', status)
            IF_NOTOK_RETURN(status=1)

            if (flask_verbose) write (*,'(a, " :: ",i8," obs in input file for tracer ", a)') rname, nflasks, names(itrac)

            allocate(mask(nflasks))
            allocate(itau_start(nflasks))
            allocate(itau_center(nflasks))
            allocate(itau_end(nflasks))
            allocate(region_list(nflasks))

            id = nc_read_var(grp_id, 'id', status)
            IF_NOTOK_RETURN(status=1)
            idate_f = nc_read_var(grp_id, 'date_components', status) ! shape of idate_f is now 6 x id
            IF_NOTOK_RETURN(status=1)
            lat = nc_read_var(grp_id, 'lat', status)
            IF_NOTOK_RETURN(status=1)
            lon = nc_read_var(grp_id, 'lon', status)
            IF_NOTOK_RETURN(status=1)

            ! For the vertical coordinate, if only 'alt' or 'pres' exists as a variable, just read that
            ! and set vert_coord_type correspondingly. Otherwise, read both and also read vert_coord_type.
            if (nc_var_exists(grp_id, 'alt')) alt = nc_read_var(grp_id, 'alt', status)
            if (nc_var_exists(grp_id, 'pres')) pres = nc_read_var(grp_id, 'pres', status)

            if (allocated(alt) .and. allocated(pres)) then
                vert_coord_type = nc_read_var(grp_id, 'vert_coord_type', status)
            else
                allocate(vert_coord_type(size(lat)))
                if (allocated(alt)) vert_coord_type = VERT_COORD_ALT
                if (allocated(pres)) vert_coord_type = VERT_COORD_PRES
            end if

            sampling_strategy = nc_read_var(grp_id, 'sampling_strategy', status)
            IF_NOTOK_RETURN(status=1)
            station_id = nc_read_var(grp_id, 'station_id', status)
            IF_NOTOK_RETURN(status=1)

            ! We look for another variable called time_window_length, which is the number of seconds on each side of the sample
            ! time over which to average. In older input files this variable may not be there, in which case we fill it up with
            ! dummy values. In any case, this variable is only relevant for those samples for which sampling_strategy is 4.
            if (nc_var_exists(grp_id, 'time_window_length')) then
                time_window_length = nc_read_var(grp_id, 'time_window_length')
            else
                allocate(time_window_length(size(lat)))
                time_window_length = 0
            end if

            !$omp parallel private (iflask, i_region, region) &
            !$omp private (idate_temp, windowMid, windowBeg, windowEnd)
            !$omp do schedule(static)
            do iflask = 1, nflasks

                call date2tau(idate_f(:,iflask),itau_center(iflask))

                ! assign region
                do i_region = 1, nregions
                     region = region_rank(i_region)
                    if (in_region(lat(iflask), lon(iflask), region)) exit
                end do
                ! In case the lat/lon are invalid, set region = -1, which will be fitlered out later
                if ((abs(lat(iflask)) > 90.0) .or. (abs(lon(iflask)) > 180.0) .or. alt(iflask) < -100) then
                    region = -1
                    write(*,'(a, " :: Flask sample ", i6, " has an invalid latitude and/or longitude, skipped")') &
                        rname, iflask
                end if

                region_list(iflask) = region

                ! For invalid samples, no need to calculate time window
                if (region_list(iflask) <= 0) cycle

                select case (sampling_strategy(iflask))

                case (1) ! N-hour window

                    itau_start(iflask) = itau_center(iflask)-assim_window*3600
                    itau_end(iflask) = itau_center(iflask)+assim_window*3600

                    if (flask_verbose) then
                        if ((itau_center(iflask) .ge. itaui) .and. (itau_center(iflask) .lt. itaue)) then
                            call tau2date(itau_center(iflask), idate_temp)
                            windowMid = NewDate(time6=idate_temp)
                            call tau2date(itau_start(iflask), idate_temp)
                            windowBeg = NewDate(time6=idate_temp)
                            call tau2date(itau_end(iflask), idate_temp)
                            windowEnd = NewDate(time6=idate_temp)
                            write(*,'(a, " :: Flask sample at ", a, " in region ", a, " will be sampled from ", a, " to ", a)') rname, &
                                trim(Pretty(windowMid)), trim(region_name(region_list(iflask))), trim(Pretty(windowBeg)), trim(Pretty(windowEnd))
                        end if
                    end if

                case (2) ! instantaneous

                    itau_start(iflask) = itau_center(iflask)
                    itau_end(iflask) = itau_center(iflask)

                case (3) ! dT sampling

                    itau_start(iflask) = itau_center(iflask) - mod(itau_center(iflask) - itaui, ndyn_max/tref(region_list(iflask))) - 1 ! 1 second leeway on either side
                    itau_end(iflask)   = itau_start(iflask)  + ndyn_max/tref(region_list(iflask)) + 1 ! 1 second leeway on either side

                    if (flask_verbose) then
                        if ((itau_center(iflask) .ge. itaui) .and. (itau_center(iflask) .lt. itaue)) then
                            call tau2date(itau_center(iflask), idate_temp)
                            windowMid = NewDate(time6=idate_temp)
                            call tau2date(itau_start(iflask), idate_temp)
                            windowBeg = NewDate(time6=idate_temp)
                            call tau2date(itau_end(iflask), idate_temp)
                            windowEnd = NewDate(time6=idate_temp)
                            write(*,'(a, " :: Flask sample in region ", a, " at ", a, " will be sampled from ", a, " to ", a)') rname, &
                                trim(region_name(region_list(iflask))), trim(Pretty(windowMid)), trim(Pretty(windowBeg)), trim(Pretty(windowEnd))
                        end if
                    end if

                case (4) ! sampling time window provided per sample

                    itau_start(iflask) = itau_center(iflask) - time_window_length(iflask)
                    itau_end(iflask)   = itau_center(iflask) + time_window_length(iflask)

                    if (time_window_length(iflask) .le. 0) write(*,'(a, " :: WARNING: flask sample ", i6, " for tracer ", a, &
                        " has a non-positive window length")') rname, iflask, names(itrac)

                case default

                    write (0,'("[", a, "] Flask with id ",i10," and event number ",i10,":")') rname, id(iflask), station_id(iflask)
                    write (0, '("  Unknown sampling strategy = ",i10,".")') sampling_strategy(iflask)

                end select

            end do ! iflask
            !$omp end do
            !$omp end parallel

            mask = (itau_center .ge. itaui) .and. (itau_center .lt. itaue)  ! sample only itaui <= obs < itaue
            !mask=((itau_end .gt. itaui) .and. (itau_start .lt. itaue))  ! sample all obs whose sampling windows fall inside (itaui,itaue)

            ! Some basic sanity checks, such as latitudes between -90 and 90, longitudes between -180 and 180. If a
            ! latitude/longitude is invalid, in_region returns -1 as the region. Filter those out. In fact, region should always
            ! be a natural number, so zero is out as well.
            mask = mask .and. (region_list > 0)

            nflasks=count(mask)  ! note that this changes the value of nflasks
            if (flask_verbose) write (*,'(a, " :: ",i4," obs in time range ",i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2," to ",i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,".")') &
                rname, nflasks, idatei(1:5), idatee(1:5)

            if(nflasks .eq. 0) then
                if(allocated(idate_f)) deallocate(idate_f)
                if(allocated(mask)) deallocate(mask)
                if(allocated(id)) deallocate(id)
                if(allocated(itau_center)) deallocate(itau_center)
                if(allocated(itau_start)) deallocate(itau_start)
                if(allocated(itau_end)) deallocate(itau_end)
                if(allocated(time_window_length)) deallocate(time_window_length)
                if(allocated(lat)) deallocate(lat)
                if(allocated(lon)) deallocate(lon)
                if(allocated(alt)) deallocate(alt)
                if(allocated(pres)) deallocate(pres)
                if(allocated(vert_coord_type)) deallocate(vert_coord_type)
                if(allocated(region_list)) deallocate(region_list)
                if(allocated(sampling_strategy)) deallocate(sampling_strategy)
                if(allocated(station_id)) deallocate(station_id)
                IO_tracer(itrac)%flask_active = .false.
            end if

            if(allocated(IO_tracer(itrac)%flasks)) deallocate(IO_tracer(itrac)%flasks)  ! avoid double allocation
            IO_tracer(itrac)%nflasks = nflasks

            if (.not. IO_tracer(itrac)%flask_active) cycle

            allocate(IO_tracer(itrac)%flasks(nflasks))

            do iflask=1,nflasks
                IO_tracer(itrac)%flasks(iflask)%mix = 0.0
                IO_tracer(itrac)%flasks(iflask)%var = 0.0
                IO_tracer(itrac)%flasks(iflask)%nsamples = 0
                IO_tracer(itrac)%flasks(iflask)%weight = 0.0
                IO_tracer(itrac)%flasks(iflask)%evaluated = .false.
                if(flask_sample_meteo) then
                    IO_tracer(itrac)%flasks(iflask)%u = 0.0
                    IO_tracer(itrac)%flasks(iflask)%v = 0.0
                    IO_tracer(itrac)%flasks(iflask)%blh = 0.0
                    IO_tracer(itrac)%flasks(iflask)%q = 0.0
                    IO_tracer(itrac)%flasks(iflask)%pressure = 0.0
                    IO_tracer(itrac)%flasks(iflask)%temperature = 0.0
                end if
            end do

            IO_tracer(itrac)%flasks(:)%id=pack(id,mask)
            IO_tracer(itrac)%flasks(:)%lat=pack(lat,mask)
            IO_tracer(itrac)%flasks(:)%lon=pack(lon,mask)
            IO_tracer(itrac)%flasks(:)%vert_coord_type = pack(vert_coord_type,mask)

            if (allocated(alt)) then
                IO_tracer(itrac)%flasks(:)%alt = pack(alt,mask)
            else
                IO_tracer(itrac)%flasks(:)%alt = -9999.9
            end if

            if (allocated(pres)) then
                IO_tracer(itrac)%flasks(:)%pres = pack(pres,mask)
            else
                IO_tracer(itrac)%flasks(:)%pres = -9999.0
            end if

            IO_tracer(itrac)%flasks(:)%region = pack(region_list, mask)
            IO_tracer(itrac)%flasks(:)%station_id=pack(station_id,mask)
            IO_tracer(itrac)%flasks(:)%itau_start=pack(itau_start,mask)
            IO_tracer(itrac)%flasks(:)%itau_center=pack(itau_center,mask)
            IO_tracer(itrac)%flasks(:)%itau_end=pack(itau_end,mask)
            IO_tracer(itrac)%flasks(:)%sampling_strategy=pack(sampling_strategy,mask)

            ! initialize structure with default and undefined values
            IO_tracer(itrac)%flasks%ifr = -1
            IO_tracer(itrac)%flasks%jfr = -1
            IO_tracer(itrac)%flasks%ifn = -1
            IO_tracer(itrac)%flasks%jfn = -1
            IO_tracer(itrac)%flasks%wcx = -1e12
            IO_tracer(itrac)%flasks%wcy = -1e12

            !$omp parallel private (iflask,idstr,idatef,region,dxr,dyr)
            !$omp do schedule(static)
            do iflask=1,nflasks

                write(idstr,*) IO_tracer(itrac)%flasks(iflask)%id
                call tau2date(IO_tracer(itrac)%flasks(iflask)%itau_center, idatef)

                if(IO_tracer(itrac)%flasks(iflask)%itau_start .lt. itaui) then
                    write (*,'("[", a, "]  attention: flask ",a,": ",i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2," close to start; sampling continued from previous run.")') &
                        rname, trim(idstr), idatef(1:5)
                end if

                if(IO_tracer(itrac)%flasks(iflask)%itau_end .gt. itaue) then
                    write (*,'("[", a, "]  attention: flask ",a,": ",i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2," close to end; sampling will continue in subsequent run.")') &
                        rname, trim(adjustl(idstr)), idatef(1:5)
                end if

                ! compute indices and slopes weighting factors

                region = IO_tracer(itrac)%flasks(iflask)%region
                dyr = dy/yref(region)
                dxr = dx/xref(region)
                ! compute indices and weighting factors

                IO_tracer(itrac)%flasks(iflask)%rif = (IO_tracer(itrac)%flasks(iflask)%lon-float(xbeg(region)))/dxr + 0.99999
                IO_tracer(itrac)%flasks(iflask)%rjf = (IO_tracer(itrac)%flasks(iflask)%lat-float(ybeg(region)))/dyr + 0.99999

                IO_tracer(itrac)%flasks(iflask)%ifr  = int(IO_tracer(itrac)%flasks(iflask)%rif)   ! i-index of grid cell in which observation is located
                IO_tracer(itrac)%flasks(iflask)%jfr  = int(IO_tracer(itrac)%flasks(iflask)%rjf)   ! j-index of grid cell in which observation is located

                !fraction from the center of the is-box  (-0.5---+0.5)
                IO_tracer(itrac)%flasks(iflask)%rif = IO_tracer(itrac)%flasks(iflask)%rif-IO_tracer(itrac)%flasks(iflask)%ifr-0.5
                !idem js
                IO_tracer(itrac)%flasks(iflask)%rjf = IO_tracer(itrac)%flasks(iflask)%rjf-IO_tracer(itrac)%flasks(iflask)%jfr-0.5

                !the neighbour for pressure interpolation
                if(IO_tracer(itrac)%flasks(iflask)%rif .gt. 0) then
                    IO_tracer(itrac)%flasks(iflask)%ifn = IO_tracer(itrac)%flasks(iflask)%ifr+1
                else
                    IO_tracer(itrac)%flasks(iflask)%ifn = IO_tracer(itrac)%flasks(iflask)%ifr-1
                end if

                !the neighbour for y interpolation
                if(IO_tracer(itrac)%flasks(iflask)%rjf .gt. 0) then
                    IO_tracer(itrac)%flasks(iflask)%jfn = IO_tracer(itrac)%flasks(iflask)%jfr+1
                else
                    IO_tracer(itrac)%flasks(iflask)%jfn = IO_tracer(itrac)%flasks(iflask)%jfr-1
                end if

                ! x- / y-weighting of grid cell in which observation is located
                IO_tracer(itrac)%flasks(iflask)%wcx = (1.0-abs(IO_tracer(itrac)%flasks(iflask)%rif))    ! 1.0 ... 0.5
                IO_tracer(itrac)%flasks(iflask)%wcy = (1.0-abs(IO_tracer(itrac)%flasks(iflask)%rjf))    ! 1.0 ... 0.5

                !=================================================================
                ! if index of neighbour is exceeding range of region set
                ! neighbour = current cell (i.e. no interpolation)
                ! in case of cyclic x-boundaries take corresponding cyclic i index
                !=================================================================
                !if ( IO_tracer(itrac)%flasks(iflask)%jfn < 1) IO_tracer(itrac)%flasks(iflask)%jfn=1
                !if ( IO_tracer(itrac)%flasks(iflask)%jfn > jm(region) ) IO_tracer(itrac)%flasks(iflask)%jfn=jm(region)
                if ( IO_tracer(itrac)%flasks(iflask)%jfn < jsr(region)) IO_tracer(itrac)%flasks(iflask)%jfn = jsr(region)
                if ( IO_tracer(itrac)%flasks(iflask)%jfn > jer(region)) IO_tracer(itrac)%flasks(iflask)%jfn = jer(region)
                if ( xcyc(region) == 0 ) then
                    ! non-cyclic boundaries
                    !if ( IO_tracer(itrac)%flasks(iflask)%ifn < 1) IO_tracer(itrac)%flasks(iflask)%ifn=1
                    !if ( IO_tracer(itrac)%flasks(iflask)%ifn > im(region) ) IO_tracer(itrac)%flasks(iflask)%ifn=im(region)
                    if ( IO_tracer(itrac)%flasks(iflask)%ifn < isr(region)) IO_tracer(itrac)%flasks(iflask)%ifn = isr(region)
                    if ( IO_tracer(itrac)%flasks(iflask)%ifn > ier(region)) IO_tracer(itrac)%flasks(iflask)%ifn = ier(region)
                else
                    ! cyclic x-boundaries
                    if ( IO_tracer(itrac)%flasks(iflask)%ifn < 1 ) IO_tracer(itrac)%flasks(iflask)%ifn=im(region)
                    if ( IO_tracer(itrac)%flasks(iflask)%ifn > im(region) ) IO_tracer(itrac)%flasks(iflask)%ifn=1
                end if

            end do ! iflask
            !$omp end do
            !$omp end parallel

            if (flask_verbose) then
                write(*,'(a, " :: list of observations to be sampled during this simulation:")') rname
                write(*,'(a, " ::    flask    ID     region  longitude   (i)  latitude   (j)  altitude        date         time eventnumber")') rname
                do iflask = 1,nflasks
                    call tau2date(IO_tracer(itrac)%flasks(iflask)%itau_center,idatef)
                    write(*,'(a, " ::    ",i5," ",i5," ",a,"  ",f9.2," (",i3,")  ",f8.2," (",i3,")  ",f8.1,"  ",i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2, " UTC ",i11)') &
                    rname, iflask, IO_tracer(itrac)%flasks(iflask)%id, trim(region_name(IO_tracer(itrac)%flasks(iflask)%region)), &
                    IO_tracer(itrac)%flasks(iflask)%lon, IO_tracer(itrac)%flasks(iflask)%ifr, &
                    IO_tracer(itrac)%flasks(iflask)%lat, IO_tracer(itrac)%flasks(iflask)%jfr, &
                    IO_tracer(itrac)%flasks(iflask)%alt, idatef, IO_tracer(itrac)%flasks(iflask)%station_id
                end do
            end if

            ! divide by region
            allocate(IO_tracer(itrac)%flask_counter(nregions))
            do region=1,nregions
                IO_tracer(itrac)%flask_counter(region)%n_obs = count(IO_tracer(itrac)%flasks(:)%region == region)
                allocate(IO_tracer(itrac)%flask_counter(region)%iflask(IO_tracer(itrac)%flask_counter(region)%n_obs))
            end do ! region
            allocate(region_counter(nregions))
            region_counter = 1
            do iflask=1,nflasks
                region = IO_tracer(itrac)%flasks(iflask)%region
                i_region = region_counter(region)
                IO_tracer(itrac)%flask_counter(region)%iflask(i_region) = iflask
                region_counter(region) = region_counter(region) + 1
            end do ! iflask
            deallocate(region_counter)

            if (sum(IO_tracer(itrac)%flask_counter(:)%n_obs) /= nflasks) write(*,'(a,a,a)') 'WARNING :: ', rname, ' :: not all flasks accounted for'

            if(allocated(idate_f)) deallocate(idate_f)
            if(allocated(mask)) deallocate(mask)
            if(allocated(id)) deallocate(id)
            if(allocated(itau_center)) deallocate(itau_center)
            if(allocated(itau_start)) deallocate(itau_start)
            if(allocated(itau_end)) deallocate(itau_end)
            if(allocated(time_window_length)) deallocate(time_window_length)
            if(allocated(lat)) deallocate(lat)
            if(allocated(lon)) deallocate(lon)
            if(allocated(alt)) deallocate(alt)
            if(allocated(pres)) deallocate(pres)
            if(allocated(vert_coord_type)) deallocate(vert_coord_type)
            if(allocated(sampling_strategy)) deallocate(sampling_strategy)
            if(allocated(station_id)) deallocate(station_id)
            if(allocated(region_list)) deallocate(region_list)

        end do ! itrac

        call nc_close(hid)

        period_opened(i_period) = .true.
        period_exists(i_period) = .true.

        ! This period has data only if at least one tracer has data for this period
        flask_active = any(IO_tracer(:)%flask_active)

        ! If no tracer has data, deallocate IO_tracer
        if (.not. flask_active) deallocate(IO_tracer)

    else

        ! input file does not exist

        nflasks = 0
        flask_active = .false. ! just for this period
        period_exists(i_period) = .false.

    end if ! flask_file_exists

    status = 0

    write(*,'(a, " :: done")') rname

  end subroutine read_samples

  subroutine user_output_flask_evaluate
    ! average flask samples over nsamples/weights

    use datetime,       only : date2tau, tau2date
    use chem_param,     only : ntracet

    implicit none

    character(len=*), parameter ::  rname = mname//'/user_output_flask_evaluate'


    integer                 :: itrac
    integer                 :: iflask, nsamples
    character(len=15)       :: idstr
    integer, dimension(6)   :: idatef
    real                    :: weight

    if (.not. flask_active) return

    do itrac = 1, ntracet

        if (.not. IO_tracer(itrac)%flask_active) cycle

        do iflask = 1, IO_tracer(itrac)%nflasks
            if (flask_verbose) write(*,'("FLASK VERBOSE : evaluated = ", l1, " :: nsamples = ", i3)') IO_tracer(itrac)%flasks(iflask)%evaluated, IO_tracer(itrac)%flasks(iflask)%nsamples
            nsamples = IO_tracer(itrac)%flasks(iflask)%nsamples
            weight   = IO_tracer(itrac)%flasks(iflask)%weight

            if((.not. IO_tracer(itrac)%flasks(iflask)%evaluated) .and. (nsamples .ge. 1)) then
                IO_tracer(itrac)%flasks(iflask)%mix = IO_tracer(itrac)%flasks(iflask)%mix / weight
                IO_tracer(itrac)%flasks(iflask)%var = IO_tracer(itrac)%flasks(iflask)%var / weight / weight
                if (flask_sample_meteo) then
                    IO_tracer(itrac)%flasks(iflask)%u = IO_tracer(itrac)%flasks(iflask)%u / weight
                    IO_tracer(itrac)%flasks(iflask)%v = IO_tracer(itrac)%flasks(iflask)%v / weight
                    IO_tracer(itrac)%flasks(iflask)%q = IO_tracer(itrac)%flasks(iflask)%q / weight
                    IO_tracer(itrac)%flasks(iflask)%blh = IO_tracer(itrac)%flasks(iflask)%blh / weight
                    IO_tracer(itrac)%flasks(iflask)%pressure = IO_tracer(itrac)%flasks(iflask)%pressure / weight
                    IO_tracer(itrac)%flasks(iflask)%temperature = IO_tracer(itrac)%flasks(iflask)%temperature / weight
                end if
                IO_tracer(itrac)%flasks(iflask)%evaluated = .true.
            else
                write(idstr,*) IO_tracer(itrac)%flasks(iflask)%id
                call tau2date(IO_tracer(itrac)%flasks(iflask)%itau_center, idatef)
                write (*,'("[", a, "]  attention: flask ",a," (station ID ",i8,") at ",i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2," not evaluated; nsamples is ",i4,".")') &
                    rname, trim(adjustl(idstr)), IO_tracer(itrac)%flasks(iflask)%station_id, idatef(1:5), nsamples
                IO_tracer(itrac)%flasks(iflask)%mix = flask_missing_value
                IO_tracer(itrac)%flasks(iflask)%var = flask_missing_value
            end if

        end do ! iflask
    end do ! itrac

  end subroutine user_output_flask_evaluate

  subroutine user_output_flask_write

    use chem_param,     only : ntracet, fscale, names
    use dims,           only : itaui, itaue, nregions, region_name
    use datetime,       only : tau2date
    use misctools,      only : check_dir
    use file_netcdf

    implicit none

    ! --- const ------------------------------
    character(len=*), parameter         :: rname = mname//'/user_output_flask_write'

    ! local
    real, allocatable       :: mix_all(:), var_all(:), avetime(:), weight(:)
    integer, allocatable    :: nsamples(:), i_obs(:)

    integer                 :: hid, grp_id, tgrp_id, n_obs, met_grp_id
    integer                 :: dim_id, dim_ntracet, dim_tracer_name_len
    integer                 :: var_id, var_ntracet, var_tracer_name_len
    integer                 :: var_tracer_names, var_flask, itrac, nobs, nflasks
    integer                 :: var_nsamples, var_avetime, var_surface_height
    integer                 :: var_u,var_v,var_blh,var_q,var_pressure,var_temperature
    integer                 :: iflask, region, i, idatei(6), idatee(6)
    character(len=1024)     :: attstring
    integer :: status

    if (.not. flask_active) return

    call user_output_flask_evaluate

    ! write results to output structure

    do itrac = 1, ntracet
        if (.not. IO_tracer(itrac)%flask_active) cycle
        nflasks=IO_tracer(itrac)%nflasks

        allocate(avetime(nflasks))
        avetime = IO_tracer(itrac)%flasks(:)%itau_end-IO_tracer(itrac)%flasks(:)%itau_start

        allocate(nsamples(nflasks))
        nsamples = IO_tracer(itrac)%flasks(:)%nsamples

        allocate(weight(nflasks))
        weight = IO_tracer(itrac)%flasks(:)%weight

        allocate(mix_all(nflasks))
        allocate(var_all(nflasks))
        ! ...and fill them
        do iflask=1,nflasks
            mix_all(iflask) = IO_tracer(itrac)%flasks(iflask)%mix
            var_all(iflask) = IO_tracer(itrac)%flasks(iflask)%var
        end do

        ! segregate according to region
        ! we want to write id, tracer names, flask (sampled mixing ratio), nsamples, averaging time, surface height, station ID
        ! optionally, u,v, blh, q, pressure, temperature
        allocate(IO_tracer(itrac)%flask_output(nregions))
        allocate(i_obs(nregions))
        ! type flask_output_data
        !    integer, allocatable, dimension(:)  :: flask_id, nsamples, station_id
        !    real, allocatable, dimension(:)     :: avg_time, weight
        !    real, allocatable, dimension(:)     :: u, v, blh, q, pressure, temperature
        !    real, allocatable, dimension(:,:)   :: mix_ratio, mix_ratio_sigma
        !    integer                            :: n_obs
        ! end type flask_output_data
        do region = 1, nregions

            nobs = count(IO_tracer(itrac)%flasks(:)%region == region .and. IO_tracer(itrac)%flasks(:)%evaluated)
            IO_tracer(itrac)%flask_output(region)%n_obs = nobs
            allocate(IO_tracer(itrac)%flask_output(region)%flask_id(nobs))
            allocate(IO_tracer(itrac)%flask_output(region)%nsamples(nobs))
            allocate(IO_tracer(itrac)%flask_output(region)%weight(nobs))
            allocate(IO_tracer(itrac)%flask_output(region)%sampling_strategy(nobs))
            allocate(IO_tracer(itrac)%flask_output(region)%station_id(nobs))
            allocate(IO_tracer(itrac)%flask_output(region)%avg_time(nobs))
            allocate(IO_tracer(itrac)%flask_output(region)%mix_ratio(nobs))
            allocate(IO_tracer(itrac)%flask_output(region)%mix_ratio_sigma(nobs))
            if (flask_sample_meteo) then
                allocate(IO_tracer(itrac)%flask_output(region)%u(nobs))
                allocate(IO_tracer(itrac)%flask_output(region)%v(nobs))
                allocate(IO_tracer(itrac)%flask_output(region)%blh(nobs))
                allocate(IO_tracer(itrac)%flask_output(region)%q(nobs))
                allocate(IO_tracer(itrac)%flask_output(region)%pressure(nobs))
                allocate(IO_tracer(itrac)%flask_output(region)%temperature(nobs))
            end if
        end do

        i_obs(:) = 1
        do i = 1, nflasks
            if (IO_tracer(itrac)%flasks(i)%evaluated) then
                region = IO_tracer(itrac)%flasks(i)%region
                IO_tracer(itrac)%flask_output(region)%flask_id(i_obs(region))           = IO_tracer(itrac)%flasks(i)%id
                IO_tracer(itrac)%flask_output(region)%nsamples(i_obs(region))           = nsamples(i)
                IO_tracer(itrac)%flask_output(region)%weight(i_obs(region))             = weight(i)
                IO_tracer(itrac)%flask_output(region)%sampling_strategy(i_obs(region))  = IO_tracer(itrac)%flasks(i)%sampling_strategy
                IO_tracer(itrac)%flask_output(region)%station_id(i_obs(region))         = IO_tracer(itrac)%flasks(i)%station_id
                IO_tracer(itrac)%flask_output(region)%avg_time(i_obs(region))           = avetime(i)
                IO_tracer(itrac)%flask_output(region)%mix_ratio(i_obs(region))          = mix_all(i)
                IO_tracer(itrac)%flask_output(region)%mix_ratio_sigma(i_obs(region))    = sqrt(var_all(i))
                if (flask_sample_meteo) then
                    IO_tracer(itrac)%flask_output(region)%u(i_obs(region))              = IO_tracer(itrac)%flasks(i)%u
                    IO_tracer(itrac)%flask_output(region)%v(i_obs(region))              = IO_tracer(itrac)%flasks(i)%v
                    IO_tracer(itrac)%flask_output(region)%blh(i_obs(region))            = IO_tracer(itrac)%flasks(i)%blh
                    IO_tracer(itrac)%flask_output(region)%q(i_obs(region))              = IO_tracer(itrac)%flasks(i)%q
                    IO_tracer(itrac)%flask_output(region)%pressure(i_obs(region))       = IO_tracer(itrac)%flasks(i)%pressure
                    IO_tracer(itrac)%flask_output(region)%temperature(i_obs(region))    = IO_tracer(itrac)%flasks(i)%temperature
                end if
                i_obs(region) = i_obs(region) + 1
            end if
        end do
        do region=1,nregions
            if (i_obs(region)-IO_tracer(itrac)%flask_output(region)%n_obs /= 1) write(*,'(a,i1)') 'There is a problem in counting point samples in region ', region
        end do

        deallocate(mix_all)
        deallocate(var_all)
        deallocate(avetime)
        deallocate(nsamples)
        deallocate(weight)
        deallocate(i_obs)

    end do !itrac

    ! new file:
    call check_dir(outfile_point)
    hid = nc_open(trim(outfile_point), 'c', status)
    IF_NOTOK_RETURN(status=1)

    ! define dimensions:
    call nc_create_dim(hid, 'tracer', ntracet)

    call tau2date(itaui,idatei)
    call tau2date(itaue,idatee)

    write(attstring,'(i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2, " UTC")') idatei
    call nc_set_attrs(hid, "model_start_date", trim(attstring))
    write(attstring,'(i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2, " UTC")') idatee
    call nc_set_attrs(hid, "model_end_date", trim(attstring))

    ! now create the groups and write the variables
    do region = 1, nregions
        grp_id = nc_create_group(hid, region_name(region))
        do itrac = 1, ntracet
            if (.not. IO_tracer(itrac)%flask_active) cycle
            n_obs = IO_tracer(itrac)%flask_output(region)%n_obs
            if (n_obs > 0) then
                tgrp_id = nc_create_group(grp_id, names(itrac))
                call nc_create_dim(tgrp_id, 'samples', n_obs)
                call nc_dump_var(tgrp_id, 'id', (/'samples'/), IO_tracer(itrac)%flask_output(region)%flask_id(1:n_obs))
                call nc_dump_var(tgrp_id, 'mixing_ratio', (/'samples'/), IO_tracer(itrac)%flask_output(region)%mix_ratio(1:n_obs))
                call nc_dump_var(tgrp_id, 'mixing_ratio_sigma', (/'samples'/), IO_tracer(itrac)%flask_output(region)%mix_ratio_sigma(1:n_obs))
                call nc_dump_var(tgrp_id, 'nsamples', (/'samples'/), IO_tracer(itrac)%flask_output(region)%nsamples(1:n_obs))
                call nc_dump_var(tgrp_id, 'total_weight', (/'samples'/), IO_tracer(itrac)%flask_output(region)%weight(1:n_obs))
                call nc_dump_var(tgrp_id, 'sampling_strategy', (/'samples'/), IO_tracer(itrac)%flask_output(region)%sampling_strategy(1:n_obs))
                call nc_dump_var(tgrp_id, 'averaging_time', (/'samples'/), IO_tracer(itrac)%flask_output(region)%avg_time(1:n_obs))
                call nc_dump_var(tgrp_id, 'station_id', (/'samples'/), IO_tracer(itrac)%flask_output(region)%station_id(1:n_obs))
                call nc_dump_var(tgrp_id, 'ifn', (/'samples'/), IO_tracer(itrac)%flasks(:)%ifn)
                call nc_dump_var(tgrp_id, 'ifr', (/'samples'/), IO_tracer(itrac)%flasks(:)%ifr)
                call nc_dump_var(tgrp_id, 'jfn', (/'samples'/), IO_tracer(itrac)%flasks(:)%jfn)
                call nc_dump_var(tgrp_id, 'jfr', (/'samples'/), IO_tracer(itrac)%flasks(:)%jfr)
                call nc_dump_var(tgrp_id, 'rif', (/'samples'/), IO_tracer(itrac)%flasks(:)%rif)
                call nc_dump_var(tgrp_id, 'rjf', (/'samples'/), IO_tracer(itrac)%flasks(:)%rjf)
                call nc_dump_var(tgrp_id, 'lfr', (/'samples'/), IO_tracer(itrac)%flasks(:)%lfr)
                call nc_dump_var(tgrp_id, 'rlf', (/'samples'/), IO_tracer(itrac)%flasks(:)%rlf)
                if (flask_sample_meteo) then
                    met_grp_id = nc_create_group(tgrp_id, 'meteo')
                    call nc_dump_var(met_grp_id, 'boundary_layer_height', (/'samples'/), IO_tracer(itrac)%flask_output(region)%blh(1:n_obs))
                    call nc_dump_var(met_grp_id, 'specific_humidity', (/'samples'/), IO_tracer(itrac)%flask_output(region)%q(1:n_obs))
                    call nc_dump_var(met_grp_id, 'pressure', (/'samples'/), IO_tracer(itrac)%flask_output(region)%pressure(1:n_obs))
                    call nc_dump_var(met_grp_id, 'temperature', (/'samples'/), IO_tracer(itrac)%flask_output(region)%temperature(1:n_obs))
                    call nc_dump_var(met_grp_id, 'zonal_wind', (/'samples'/), IO_tracer(itrac)%flask_output(region)%u(1:n_obs))
                    call nc_dump_var(met_grp_id, 'meridional_wind', (/'samples'/), IO_tracer(itrac)%flask_output(region)%v(1:n_obs))
                end if ! flask_sample_meteo
            end if ! n_obs > 0
        end do ! itrac
    end do ! region

    ! close file:
    call nc_close(hid)

    ! free up the IO_tracer data structure
    deallocate(IO_tracer)

  end subroutine user_output_flask_write

  subroutine user_output_flask_sample(region, tr, status)

    use global_data,    only : mass_dat, region_dat, conv_dat
    use MeteoData,      only : gph_dat, phlb_dat, humid_dat, m_dat, temper_dat
    use dims,           only : lm, itaur, ndyn, tref
    use chem_param,     only : fscale, names, ntracet
    use datetime,       only : tau2date, date2tau
    use Go,             only : TDate, Pretty, Get, rTotal, operator(-)

    implicit none

    ! input/output
    integer, intent(in)             :: region
    type(TDate), intent(in)         :: tr(2)
    integer, intent(out)            :: status

    real,dimension(:,:,:), pointer          :: m,gph
    real,dimension(:,:,:,:), pointer        :: rm, rxm, rym, rzm
    real,dimension(:,:), pointer            :: blh ! boundary layer height [m]
    real,dimension(:,:,:), pointer          :: T ! temperature
    real,dimension(:,:,:), pointer          :: phlb ! pressure grid boundaries
    real,dimension(:,:,:), pointer          :: pu ! mass flux x-direction [kg/s]
    real,dimension(:,:,:), pointer          :: pv ! mass flux y-direction [kg/s]
    real,dimension(:,:,:), pointer          :: q ! specific humidity [kg/kg]

    real, dimension(0:lm(1))    :: height, pressure
    integer                     :: itrac, iflask, lfrt, lfrtn, lmr, lfrn
    integer                     :: l, itr, i_region
    integer                     :: n, i, j, k
    integer                     :: ifr,jfr,lfr,ifn,jfn
    integer(2)                  :: vchoice
    real                        :: alt,rlf, wcz, pres
    real                        :: rmf, mean_mix, var_mix, weight
    real                        :: wcx,wcy,rif,rjf
    real                        :: mix_n(-1:1,-1:1,-1:1)   ! static for OMP!
    logical                     :: in_window
    integer, dimension(6)       :: idate_temp
    integer                     :: itau_tr(2)

    ! --- const ------------------------------
    character(len=*), parameter             ::  rname = mname//'/user_output_flask_sample'

    status = 0

    if (.not. flask_active) return

    weight = rTotal(tr(2)-tr(1), 'sec')

    ! pointers to global arrays

    m    => m_dat(region)%data
    rm   => mass_dat(region)%rm_t
    rxm  => mass_dat(region)%rxm_t
    rym  => mass_dat(region)%rym_t
    rzm  => mass_dat(region)%rzm_t
    gph  => gph_dat(region)%data
    blh  => conv_dat(region)%blh
    t    => temper_dat(region)%data
    phlb => phlb_dat(region)%data
    q    => humid_dat(region)%data

    lmr = lm(region)

    ! convert tr to integer
    call Get(tr(1), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(1))
    call Get(tr(2), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(2))

    do itrac = 1, ntracet
        if (.not. IO_tracer(itrac)%flask_active) cycle

        !$omp parallel do schedule(static) &
        !$omp private (iflask, in_window, i_region, lfrn, wcz, height, pressure, pres) &
        !$omp private (alt, ifr, ifn, jfr, jfn, rif, rjf, wcx, wcy, lfr, rlf, vchoice) &
        !$omp private (itr, rmf, mix_n, i, j, k, l, mean_mix, var_mix) &
        !$omp reduction (+:status)
        do i_region = 1, IO_tracer(itrac)%flask_counter(region)%n_obs
           iflask = IO_tracer(itrac)%flask_counter(region)%iflask(i_region)
           ! Is model time in the sampling window?
           in_window = .false.
           ! itaur(region) points to the beginning of the dynamic time step

           select case (IO_tracer(itrac)%flasks(iflask)%sampling_strategy)

              case (1,3,4) ! 4-hour average, or dT sampling, or sampling with custom time window

                ! the dynamic timestep tr(1) --> tr(2) must be completely inside the interval (itau_start, itau_end)
                if ((itau_tr(1) .ge. IO_tracer(itrac)%flasks(iflask)%itau_start) .and. &
                    (itau_tr(2) .le. IO_tracer(itrac)%flasks(iflask)%itau_end)) in_window = .true.

              case (2) ! instantaneous

                if ((itau_tr(1) .le. IO_tracer(itrac)%flasks(iflask)%itau_center) .and. &
                    (itau_tr(2) .gt. IO_tracer(itrac)%flasks(iflask)%itau_center)) in_window = .true.

              case default

                 write(gol,*) 'Unknown sampling strategy, skipping flask'; call goErr

           end select

           if(in_window)  then
              ! Use slopes to determine tracer concentrations at the site.
              vchoice = IO_tracer(itrac)%flasks(iflask)%vert_coord_type
              alt = IO_tracer(itrac)%flasks(iflask)%alt
              pres = IO_tracer(itrac)%flasks(iflask)%pres
              ifr = IO_tracer(itrac)%flasks(iflask)%ifr
              ifn = IO_tracer(itrac)%flasks(iflask)%ifn
              jfr = IO_tracer(itrac)%flasks(iflask)%jfr
              jfn = IO_tracer(itrac)%flasks(iflask)%jfn
              rif = IO_tracer(itrac)%flasks(iflask)%rif
              rjf = IO_tracer(itrac)%flasks(iflask)%rjf
              wcx = IO_tracer(itrac)%flasks(iflask)%wcx
              wcy = IO_tracer(itrac)%flasks(iflask)%wcy


              ! interpolate the altitude to site position...
              lfr = 1 !layer

              select case (vchoice)

                case (VERT_COORD_ALT)
                  do l=0,lm(region)
                     height(l) = wcx * wcy  * gph(ifr,jfr,l+1) + (1.0-wcx) * wcy * gph(ifn,jfr,l+1) + &
                        wcx * (1.0-wcy) * gph(ifr,jfn,l+1) + (1.0-wcx) *(1.0-wcy) * gph(ifn,jfn,l+1)
                     if(l==0) IO_tracer(itrac)%flasks(iflask)%surface_height = height(0)
                  end do

                  do l=0,lm(region) ! selects layer , note that we start from second layer from surface
                     if(height(l) .gt. alt) exit
                  end do

                  select case(l)
                  case(0)
                        if (.not. IO_tracer(itrac)%flasks(iflask)%below_surface_warning) then
                            if (flask_verbose) then
                                write (*,'(a, " :: WARNING:  For flask ",i8,":")') rname, IO_tracer(itrac)%flasks(iflask)%id
                                write (*,'(a, " :: Sample altitude of ",f8.2," m is below surface height of ",f8.2," m.")') rname, alt, height(0)
                                write (*,'(a, " :: Will sample at surface.")') rname
                            end if
                            IO_tracer(itrac)%flasks(iflask)%below_surface_warning = .True.
                        end if
                        lfr = 1
                        rlf = -0.5  !surface...
                  case default
                        lfr = l  !the site layer
                        ! the offset from the center of the layer (-0.5--->+0.5)
                        ! (interpolation is in (m))
                        rlf = (alt-height(l-1))/(height(l)-height(l-1)) - 0.5
                  end select

                case (VERT_COORD_PRES)
                  do l=0,lm(region)
                     pressure(l) = wcx * wcy  * phlb(ifr,jfr,l+1) + (1.0-wcx) * wcy * phlb(ifn,jfr,l+1) + &
                        wcx * (1.0-wcy) * phlb(ifr,jfn,l+1) + (1.0-wcx) *(1.0-wcy) * phlb(ifn,jfn,l+1)
                  end do

                  do l=0,lm(region) ! selects layer , note that we start from second layer from surface
                     if(pressure(l) .gt. pres) exit
                  end do

                  select case(l)
                  case(0)
                        if (.not. IO_tracer(itrac)%flasks(iflask)%below_surface_warning) then
                            if (flask_verbose) then
                                write (*,'(a, " :: WARNING:  For flask ",i8,":")') rname, IO_tracer(itrac)%flasks(iflask)%id
                                write (*,'(a, " :: Sample pressure of ",f8.2," Pa is below surface pressure of ",f8.2," Pa.")') rname, pres, pressure(0)
                                write (*,'(a, " :: Will sample at surface.")') rname
                            end if
                            IO_tracer(itrac)%flasks(iflask)%below_surface_warning = .True.
                        end if
                        lfr = 1
                        rlf = -0.5  !surface...
                  case default
                        lfr = l  !the site layer
                        ! the offset from the center of the layer (-0.5--->+0.5)
                        ! (interpolation is in (m))
                        rlf = (pres-pressure(l-1))/(pressure(l)-pressure(l-1)) - 0.5
                  end select

              end select ! vchoice

              IO_tracer(itrac)%flasks(iflask)%lfr = lfr

              !=================================
              !the neighbour for z interpolation
              !=================================
              if ( rlf .gt. 0 ) then
                 lfrn = lfr+1
              else
                 lfrn = lfr-1
              end if
              ! z-weighting of grid cell in which observation is located
              wcz = (1.0-abs(rlf))  !.0 ... 0.5

              !=========================================================
              ! if vertical neighbor is 0 (which does not exist)
              ! take vertical layer with l=2 for EXTRApolation to ground
              !=========================================================

              IF(lfrn == 0) THEN
                 lfrn=2
                 wcz=1.0-rlf  ! 1.0 ... 1.5
              END IF
              IF(lfrn == lmr+1) THEN
                 !=========================================================
                 ! if vertical neighbor is lmr+1 (which does not exist)
                 ! -> no interpolation
                 !=========================================================
                 lfrn=lmr ! no interpolation
                 wcz=1.0
              END IF

              IO_tracer(itrac)%flasks(iflask)%rlf = rlf

              ! rm-value is obtained from rm + slopes.
              ! slope = rxm = (rm*dX/dx *deltaX/2)

              select case (flask_interpolation)
                case (POINT_INTERPOLATION_GRIDBOX)
                    rmf = rm(ifr,jfr,lfr,itrac) / m(ifr,jfr,lfr) * fscale(itrac)

                case (POINT_INTERPOLATION_SLOPES)
                    rmf = ( rm(ifr,jfr,lfr,itrac) + 2.0*(rif*rxm(ifr,jfr,lfr,itrac) + rjf*rym(ifr,jfr,lfr,itrac) + &
                        rlf*rzm(ifr,jfr,lfr,itrac)) ) / m(ifr,jfr,lfr) * fscale(itrac)

                case (POINT_INTERPOLATION_LINEAR)
                    rmf = ( &
                             wcx  *      wcy  *      wcz  * rm(ifr,jfr,lfr ,itrac) / m(ifr,jfr,lfr )  + &
                        (1.0-wcx) *      wcy  *      wcz  * rm(ifn,jfr,lfr ,itrac) / m(ifn,jfr,lfr )  + &
                             wcx  * (1.0-wcy) *      wcz  * rm(ifr,jfn,lfr ,itrac) / m(ifr,jfn,lfr )  + &
                        (1.0-wcx) * (1.0-wcy) *      wcz  * rm(ifn,jfn,lfr ,itrac) / m(ifn,jfn,lfr )  + &
                             wcx  *      wcy  * (1.0-wcz) * rm(ifr,jfr,lfrn,itrac) / m(ifr,jfr,lfrn)  + &
                        (1.0-wcx) *      wcy  * (1.0-wcz) * rm(ifn,jfr,lfrn,itrac) / m(ifn,jfr,lfrn)  + &
                             wcx  * (1.0-wcy) * (1.0-wcz) * rm(ifr,jfn,lfrn,itrac) / m(ifr,jfn,lfrn)  + &
                        (1.0-wcx) * (1.0-wcy) * (1.0-wcz) * rm(ifn,jfn,lfrn,itrac) / m(ifn,jfn,lfrn)) * fscale(itrac)

                case default
                    write (gol,'("unsupported point interpolation index ",i6)') flask_interpolation; call goErr
                    status = status+1

              end select

              if (flask_verbose) then
                write(*, '(a, " :: Sampled ", f10.4, " from tracer ", a4, " in cell ", i3, ",", i3, ",", i3, " between ", a, " and ", a)') &
                    rname, rmf, trim(names(itrac)), ifr, jfr, lfr, trim(Pretty(tr(1))), trim(Pretty(tr(2)))
                write(*, '(a, " :: rif = ", f9.6, ", rjf = ", f9.6, ", rlf = ", f9.6, ", rxm = ", f15.5, ", rym = ", f15.5, ", rzm = ", f15.5)') &
                    rname, rif, rjf, rlf, rxm(ifr,jfr,lfr,itrac), rym(ifr,jfr,lfr,itrac), rzm(ifr,jfr,lfr,itrac)
                write(*, '(a, " :: Tracer mass = ", es20.12, ", air mass = ", es20.12, ", fscale = ", es18.10)') &
                    rname, rm(ifr,jfr,lfr,itrac), m(ifr,jfr,lfr), fscale(itrac)
              end if

              IO_tracer(itrac)%flasks(iflask)%mix = IO_tracer(itrac)%flasks(iflask)%mix + rmf * weight

              ! calculate the representation error
              select case (trim(adjustl(flask_error_choice)))

                case ('neighbors')

                    do i = -1, 1
                        do j = -1, 1
                            do k = -1, 1
                                mix_n(i,j,k) = (rm(ifr,jfr,lfr,itrac) + &
                                    2.0*(i*rxm(ifr,jfr,lfr,itrac) + j*rym(ifr,jfr,lfr,itrac) + k*rzm(ifr,jfr,lfr,itrac))) / &
                                    m(ifr,jfr,lfr) * fscale(itrac)
                            end do
                        end do
                    end do
                    mean_mix = sum(mix_n)/size(mix_n)
                    mix_n = mix_n - mean_mix
                    var_mix = sum(mix_n * mix_n)/size(mix_n)

                case ('gradient')

                    var_mix = ((2.0*rxm(ifr,jfr,lfr,itrac))**2 + (2.0*rym(ifr,jfr,lfr,itrac))**2 + (2.0*rzm(ifr,jfr,lfr,itrac))**2) / &
                        m(ifr,jfr,lfr)**2 * fscale(itrac)**2

              end select

              IO_tracer(itrac)%flasks(iflask)%var = IO_tracer(itrac)%flasks(iflask)%var + &
                weight*weight * max(var_mix, IO_tracer(itrac)%flask_minerror**2)

              ! sample meteo
              if(flask_sample_meteo) then

                rmf =   wcx   *      wcy  * blh(ifr,jfr) + &
                    (1.0-wcx) *      wcy  * blh(ifn,jfr) + &
                        wcx   * (1.0-wcy) * blh(ifr,jfn) + &
                    (1.0-wcx) * (1.0-wcy) * blh(ifn,jfn)

                IO_tracer(itrac)%flasks(iflask)%blh = IO_tracer(itrac)%flasks(iflask)%blh + rmf * weight

                rmf =      wcx  *      wcy  *      wcz  * T(ifr,jfr,lfr)  + &
                      (1.0-wcx) *      wcy  *      wcz  * T(ifn,jfr,lfr)  + &
                           wcx  * (1.0-wcy) *      wcz  * T(ifr,jfn,lfr)  + &
                      (1.0-wcx) * (1.0-wcy) *      wcz  * T(ifn,jfn,lfr)  + &
                           wcx  *      wcy  * (1.0-wcz) * T(ifr,jfr,lfrn) + &
                      (1.0-wcx) *      wcy  * (1.0-wcz) * T(ifn,jfr,lfrn) + &
                           wcx  * (1.0-wcy) * (1.0-wcz) * T(ifr,jfn,lfrn) + &
                      (1.0-wcx) * (1.0-wcy) * (1.0-wcz) * T(ifn,jfn,lfrn)

                IO_tracer(itrac)%flasks(iflask)%temperature = IO_tracer(itrac)%flasks(iflask)%temperature + rmf * weight

                rmf =      wcx  *      wcy  *      wcz  * q(ifr,jfr,lfr)  + &
                      (1.0-wcx) *      wcy  *      wcz  * q(ifn,jfr,lfr)  + &
                           wcx  * (1.0-wcy) *      wcz  * q(ifr,jfn,lfr)  + &
                      (1.0-wcx) * (1.0-wcy) *      wcz  * q(ifn,jfn,lfr)  + &
                           wcx  *      wcy  * (1.0-wcz) * q(ifr,jfr,lfrn) + &
                      (1.0-wcx) *      wcy  * (1.0-wcz) * q(ifn,jfr,lfrn) + &
                           wcx  * (1.0-wcy) * (1.0-wcz) * q(ifr,jfn,lfrn) + &
                      (1.0-wcx) * (1.0-wcy) * (1.0-wcz) * q(ifn,jfn,lfrn)

                IO_tracer(itrac)%flasks(iflask)%q = IO_tracer(itrac)%flasks(iflask)%q + rmf * weight

                rmf = (((0.5-rlf) * phlb(ifr,jfr,lfr) + (0.5+rlf) * phlb(ifr,jfr,lfrn)) *      wcx  *      wcy  + &
                       ((0.5-rlf) * phlb(ifn,jfr,lfr) + (0.5+rlf) * phlb(ifn,jfr,lfrn)) * (1.0-wcx) *      wcy  + &
                       ((0.5-rlf) * phlb(ifr,jfn,lfr) + (0.5+rlf) * phlb(ifr,jfn,lfrn)) *      wcx  * (1.0-wcy) + &
                       ((0.5-rlf) * phlb(ifn,jfn,lfr) + (0.5+rlf) * phlb(ifn,jfn,lfrn)) * (1.0-wcx) * (1.0-wcy))

                IO_tracer(itrac)%flasks(iflask)%pressure = IO_tracer(itrac)%flasks(iflask)%pressure + rmf * weight


              end if ! flask_sample_meteo

              IO_tracer(itrac)%flasks(iflask)%nsamples  = IO_tracer(itrac)%flasks(iflask)%nsamples  + 1
              IO_tracer(itrac)%flasks(iflask)%weight    = IO_tracer(itrac)%flasks(iflask)%weight    + weight

           end if ! in_window
        end do ! i_region
        !$omp end parallel do

        IF_NOTOK_RETURN(status=1)

    end do ! itrac

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(gph)
    nullify(T)
    nullify(blh)
    nullify(phlb)
    nullify(q)

  end subroutine user_output_flask_sample

end module user_output_flask
