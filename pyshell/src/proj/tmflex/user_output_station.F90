!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#define CHECK_NCSTAT(nc_ret_code) if (nc_ret_code /= NF90_NOERR) then; status=1; IF_NOTOK_RETURN(status=1); end if
!
!
!###############################################################################
#include "tm5.inc"

module user_output_station

use GO,         only : gol, goErr, goPr
use go_date,    only : TDate
use os_specs,   only : DUMMY_STR_LEN

use user_output_flask_data, only : POINT_INTERPOLATION_GRIDBOX, POINT_INTERPOLATION_SLOPES, POINT_INTERPOLATION_LINEAR
use user_output_flask_data, only : POINT_ERROR_NEIGHBOR, POINT_ERROR_GRADIENT
use user_output_flask_data, only : lowest_sample_layer

implicit none

private

public :: user_output_station_init, user_output_station_step, user_output_station_done

character(len=*), parameter :: mname = 'user_output_station'
character(len=DUMMY_STR_LEN):: outdir_station, outfile_station
logical                     :: station_verbose, station_sample_meteo
logical                     :: sample_in_parent
integer                     :: station_interpolation
integer                     :: parent_sample_choice = 0
integer, allocatable        :: point_error_choice(:)
integer                     :: num_stations ! number of stations
integer                     :: io_nrecords ! number of dynamic timesteps during the run
integer                     :: netcdf_ID ! ID of the output file
integer(4), allocatable     :: midpoint_dates(:,:)
character(len=1)            :: split_period
logical, allocatable        :: period_done(:)
type(TDate), allocatable    :: period_bounds(:,:) ! n_period x 2

type T_stations
    integer                 :: site_num         ! numeric station ID
    character(len=15)       :: station_id       ! station identifier
    character(len=DUMMY_STR_LEN):: station_name ! station name
    real                    :: lat              ! latitude of station
    real                    :: lon              ! longitude of station
    real                    :: alt              ! height of station, MASL
    integer                 :: region           ! region of station
    integer(4), allocatable :: nsamples(:)      ! number of samples accumulated
    real, allocatable       :: weight(:)        ! time-interval-based weights
    real, allocatable       :: mix(:,:), var(:,:) ! mixing ratio and variance of mixing ratio
    integer                 :: ifr, jfr         ! i,j region indices for flask's grid cell
    integer                 :: lfr              ! vertical level number of the sample
    integer                 :: ifn, jfn         ! i,j region indices for flask's "next" grid cell
    real                    :: rif, rjf, rlf    ! fractions from center of ifr,jfr box
    character(len=3)        :: station_type     ! sampling type (FM = flask, CM = continuous, TM = tower, could be any three letter code)
    integer                 :: type_id          ! when stations are classified, an integer ID denoting the station type (FM, TM, ...)
    real                    :: surface_height   ! surface height in meters
    real                    :: wcx, wcy         ! x and y weighting factors for slopes interpolation
    real, allocatable       :: blh(:),q(:),pressure(:),temperature(:) ! meteorological variables
    logical                 :: below_surface_warning = .false.
end type T_stations

type T_station_region_counter
    integer, allocatable :: istat(:) ! for a region, stores the station indices corresponding to that region
    integer :: n_obs
end type T_station_region_counter

type(T_stations), allocatable :: stations(:)
type(T_station_region_counter), allocatable :: stat_counter(:)

contains

subroutine user_output_station_init(status)

    use GO,             only : gol, goErr, goPr, goSplitLine, goReadFromLine
    use GO,             only : ReadRc
    use global_data,    only : rcF
    use dims,           only : im, jm, lm, dx, dy, xref, yref, xbeg, ybeg, xend, yend, tref
    use dims,           only : nregions, region_name, idatei, idatee, xcyc
    use datetime,       only : date2tau
    use orderpack,      only : mrgrnk
    use datetime_for,   only : SEC_PER_DAY
    use go,             only : NewDate, Get_End_Of
    use chem_param,     only : names, ntracet
    use os_specs,       only : MAX_RCKEY_LEN

    implicit none

    ! --- in/out ---------------------------------
    integer, intent(out)    :: status

    ! --- local ----------------------------------
    integer                 :: sunit = 111, tau_end, tau_beg, itr, global_error_choice
    integer                 :: i, region, i_region, idate_temp(6), n_period
    integer, allocatable    :: region_rank(:), region_counter(:)
    real                    :: dxr, dyr
    type(TDate)             :: cur_date, next_date
    character(len=DUMMY_STR_LEN) :: dummy_str, stationlist_filename
    character(len=MAX_RCKEY_LEN) :: rc_key

    ! --- const ------------------------------
    character(len=*), parameter         :: rname = mname//'/user_output_station_init'

    write(*,'(a, " :: entering")') rname

    call ReadRc( rcF, 'output.station.timeseries.filename', stationlist_filename, status) ! must be the full path
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.station.verbose', station_verbose, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'output.station.meteo', station_sample_meteo, status, default=.false.)
    IF_ERROR_RETURN(status=1)

    call ReadRc( rcF, 'output.point.errors', global_error_choice, status)
    IF_NOTOK_RETURN(status=1)
    allocate(point_error_choice(ntracet))
    do itr = 1, ntracet
        write(rc_key, '("output.point.errors.", a)') trim(names(itr))
        call ReadRc( rcF, trim(rc_key), point_error_choice(itr), status, default=global_error_choice )
        IF_ERROR_RETURN(status=1)
    end do

    call ReadRc( rcF, 'output.point.interpolation', station_interpolation, status, default=POINT_INTERPOLATION_SLOPES)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'output.point.sample.parent', sample_in_parent, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.point.split.period', split_period, status, default='a')
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.dir', outdir_station, status)
    IF_NOTOK_RETURN(status=1)

    if (sample_in_parent) parent_sample_choice = 1

    ! now read the station list
    open(sunit, FORM='FORMATTED', FILE=stationlist_filename, STATUS='OLD')
    read(sunit, '(a)') dummy_str ! the first line is the header
    ! count the number of stations
    num_stations = 0
    do
        read(sunit, '(a)', END=100) dummy_str
        num_stations = num_stations + 1
    end do
100 close(sunit)
    allocate(stations(num_stations))
    ! read the station coordinates
    open(sunit, FORM='FORMATTED', FILE=stationlist_filename, STATUS='OLD')
    read(sunit, '(a)') dummy_str ! the first line is the header
    do i=1,num_stations
        read(sunit, '(a)') dummy_str
        dummy_str = adjustl(dummy_str) ! remove any leading whitespaces

        call  goReadFromLine(dummy_str, stations(i)%site_num, status, sep=' ')
        IF_NOTOK_RETURN(status=1)

        call goReadFromLine(dummy_str, stations(i)%station_id, status, sep=' ')
        IF_NOTOK_RETURN(status=1)

        call goReadFromLine(dummy_str, stations(i)%lat, status, sep=' ')
        IF_NOTOK_RETURN(status=1)

        call goReadFromLine(dummy_str, stations(i)%lon, status, sep=' ')
        IF_NOTOK_RETURN(status=1)

        call goReadFromLine(dummy_str, stations(i)%alt, status, sep=' ')
        IF_NOTOK_RETURN(status=1)

        call goReadFromLine(dummy_str, stations(i)%station_type, status, sep=' ')
        IF_NOTOK_RETURN(status=1)

        stations(i)%station_name = trim(dummy_str)
!        read(sunit, '(a6,3f8.2,1x,a2,1x,a60)') stations(i)%station_id, stations(i)%lat, stations(i)%lon, stations(i)%alt, stations(i)%station_type, stations(i)%station_name
    end do ! i
    close(sunit)

    ! initialize some variables
    stations%ifr = -1
    stations%jfr = -1
    stations%ifn = -1
    stations%jfn = -1
    stations%wcx = -1e12
    stations%wcy = -1e12

    ! calculate rif, rjf, wcx, wcy, etc.
    allocate(region_rank(nregions))
    call mrgrnk(xref(1:nregions) * yref(1:nregions), region_rank)
    region_rank = region_rank(nregions:1:-1) ! finest regions occur first
    do i=1,num_stations

        ! assign region
        do i_region=1,nregions
            region = region_rank(i_region)
            if (in_region(stations(i)%lat, stations(i)%lon, region)) exit
        end do ! region
        stations(i)%region = region

        dxr = dx/xref(region) ! dxr is the longitude grid size, so 1.0 for 1x1
        dyr = dy/yref(region) ! dyr is the latitude grid size, so 2.0 for 3x2

        stations(i)%rif = ( stations(i)%lon - float(xbeg(region)) )/dxr + 0.9999999
        stations(i)%rjf = ( stations(i)%lat - float(ybeg(region)) )/dyr + 0.9999999

        stations(i)%ifr  = int(stations(i)%rif)   ! i-index of grid cell in which observation is located
        stations(i)%jfr  = int(stations(i)%rjf)   ! j-index of grid cell in which observation is located

        !fraction from the center of the is-box and js-box  (-0.5---+0.5)
        stations(i)%rif = stations(i)%rif-stations(i)%ifr-0.5
        stations(i)%rjf = stations(i)%rjf-stations(i)%jfr-0.5

        !the neighbour for x interpolation
        if(stations(i)%rif .gt. 0) then
          stations(i)%ifn = stations(i)%ifr+1
        else
          stations(i)%ifn = stations(i)%ifr-1
        endif

        !the neighbour for y interpolation
        if(stations(i)%rjf .gt. 0) then
          stations(i)%jfn = stations(i)%jfr+1
        else
          stations(i)%jfn = stations(i)%jfr-1
        endif

        ! x- / y-weighting of grid cell in which observation is located
        stations(i)%wcx = (1.0-abs(stations(i)%rif))    ! 1.0 ... 0.5
        stations(i)%wcy = (1.0-abs(stations(i)%rjf))    ! 1.0 ... 0.5

        !=================================================================
        ! if index of neighbour is exceeding range of region set
        ! neighbour = current cell (i.e. no interpolation)
        ! in case of cyclic x-boundaries take corresponding cyclic i index
        !=================================================================
        if ( stations(i)%jfn < 1) stations(i)%jfn=1
        if ( stations(i)%jfn > jm(region) ) stations(i)%jfn=jm(region)
        if ( xcyc(region) == 0 ) then
          ! non-cyclic boundaries
          if ( stations(i)%ifn < 1) stations(i)%ifn=1
          if ( stations(i)%ifn > im(region) ) stations(i)%ifn=im(region)
        else
          ! cyclic x-boundaries
          if ( stations(i)%ifn < 1 ) stations(i)%ifn=im(region)
          if ( stations(i)%ifn > im(region) ) stations(i)%ifn=1
        endif

    end do ! i

    allocate(stat_counter(nregions))
    do region=1,nregions
        stat_counter(region)%n_obs = count(stations(:)%region == region)
        allocate(stat_counter(region)%istat(stat_counter(region)%n_obs))
    end do ! region
    allocate(region_counter(nregions))
    region_counter = 1
    do i=1,num_stations
        region = stations(i)%region
        i_region = region_counter(region)
        stat_counter(region)%istat(i_region) = i
        region_counter(region) = region_counter(region) + 1
    end do ! i
    deallocate(region_counter)
    if (sum(stat_counter(:)%n_obs) /= num_stations) write(*,'(a,a,a)') 'WARNING :: ', rname, ' :: not all stations accounted for'

    select case (split_period)
        case ('a') ! the entire period in one go
            n_period = 1
        case ('y')
            n_period = idatee(1)-idatei(1)+1
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

    allocate(period_done(n_period))
    period_done = .false.

    ! create arrays with starting and ending times per period
    allocate(period_bounds(n_period,2))
    select case (split_period)
    case ('a')
        period_bounds(1,1) = NewDate(time6=idatei)
        period_bounds(1,2) = NewDate(time6=idatee)
    case ('y')
        cur_date = NewDate(time6=idatei)
        do i = 1, n_period
            next_date = Get_End_Of(cur_date, 'year')
            period_bounds(i,1) = cur_date
            period_bounds(i,2) = next_date
            cur_date = next_date
        end do
        period_bounds(n_period,2) = NewDate(time6=idatee)
    case ('m')
        cur_date = NewDate(time6=idatei)
        do i = 1, n_period
            next_date = Get_End_Of(cur_date, 'month')
            period_bounds(i,1) = cur_date
            period_bounds(i,2) = next_date
            cur_date = next_date
        end do
        period_bounds(n_period,2) = NewDate(time6=idatee)
    case ('d')
        cur_date = NewDate(time6=idatei)
        do i = 1, n_period
            next_date = Get_End_Of(cur_date, 'day')
            period_bounds(i,1) = cur_date
            period_bounds(i,2) = next_date
            cur_date = next_date
        end do
        period_bounds(n_period,2) = NewDate(time6=idatee)
    end select

    status = 0

    write(*,'(a, " :: done")') rname

end subroutine user_output_station_init

subroutine init_station_struct(start_date, end_date, status)

    use go_date,    only : TIncrDate, iTotal, operator(-), operator(+), IncrDate, Get
    use dims,       only : ndyn_max
    use chem_param, only : ntracet

    implicit none

    type(TDate), intent(in) :: start_date, end_date
    integer, intent(out)    :: status

    character(len=*), parameter :: rname = mname//'/init_station_struct'

    type(TDate)             :: dummy_t
    type(TIncrDate)         :: dt
    integer                 :: total_seconds_in_interval, idate_temp(6), i

    ! total number of records?
    total_seconds_in_interval = iTotal(end_date - start_date, 'sec')
    io_nrecords = total_seconds_in_interval/ndyn_max
    allocate(midpoint_dates(6,io_nrecords))
    ! fill the midpoint dates
    do i=1,io_nrecords
        dt = IncrDate(sec = i*ndyn_max - ndyn_max/2)
        dummy_t = start_date + dt
        call Get(dummy_t, time6=idate_temp)
        midpoint_dates(:,i) = idate_temp(:)
    end do

    ! construct the output filename
    call Get(start_date, time6=idate_temp)
    select case (split_period)
    case ('a')
        write(outfile_station, '(a,"/stations/stations.nc4")') trim(outdir_station)
    case ('y')
        write(outfile_station, '(a,"/stations/stations_",i4.4,".nc4")') trim(outdir_station), idate_temp(1)
    case ('m')
        write(outfile_station, '(a,"/stations/stations_",i4.4,i2.2,".nc4")') trim(outdir_station), idate_temp(1:2)
    case ('d')
        write(outfile_station, '(a,"/stations/stations_",i4.4,2i2.2,".nc4")') trim(outdir_station), idate_temp(1:3)
    case default
        write(0,*) 'Wrong split period selected'
        IF_NOTOK_RETURN(status=1)
    end select

    do i = 1, num_stations

        allocate(stations(i)%mix(io_nrecords,ntracet))
        allocate(stations(i)%var(io_nrecords,ntracet))
        allocate(stations(i)%nsamples(io_nrecords))
        allocate(stations(i)%weight(io_nrecords))
        stations(i)%mix = 0.0
        stations(i)%var = 0.0
        stations(i)%nsamples = 0
        stations(i)%weight = 0.0
        stations(i)%type_id = -99 ! junk value to catch problems later, since they should all become assigned valid values later
        if (station_sample_meteo) then
            allocate(stations(i)%blh(io_nrecords))
            allocate(stations(i)%q(io_nrecords))
            allocate(stations(i)%pressure(io_nrecords))
            allocate(stations(i)%temperature(io_nrecords))
            stations(i)%blh = 0.0
            stations(i)%q = 0.0
            stations(i)%pressure = 0.0
            stations(i)%temperature = 0.0
        end if

    end do

    status = 0

end subroutine init_station_struct

subroutine user_output_station_done(status)

    implicit none

    integer, intent(out)        :: status

    character(len=*), parameter :: rname = mname//'/user_output_station_done'

    call user_output_station_write(status)

    deallocate(stations)
    deallocate(period_done, period_bounds)
    deallocate(point_error_choice)

    status = 0

end subroutine user_output_station_done

subroutine user_output_station_write(status)

    use dims,       only : region_name, idatei, idatee
    use chem_param, only : ntracet, names, mixrat_unit_name
    use go_date,    only : Pretty, NewDate
    use misctools,  only : check_dir, find_unique_strings
    use file_netcdf
    use netcdf

    implicit none

    integer, intent(out)            :: status

    character(len=*), parameter     :: rname = mname//'/user_output_station_write'

!    integer, allocatable            :: num_stations_by_type(:), type_gid(:), unique_type_idx(:)
!    character(len=2), allocatable   :: station_types(:)
    integer                         :: io_status, grp_id, meteo_grp_id, i, int_date(6), itr, total_types, i_type
    type(Tdate)                     :: form_date
    character(len=256)              :: tracer_names, unit_names
    character(len=64)               :: st_group_name

    ! first, evaluate the station data
    do i=1,num_stations
        do itr = 1, ntracet
            stations(i)%mix(:,itr) = stations(i)%mix(:,itr)/stations(i)%weight
            stations(i)%var(:,itr) = stations(i)%var(:,itr)/stations(i)%weight/stations(i)%weight
        end do ! itr
        if (station_sample_meteo) then
            stations(i)%blh = stations(i)%blh/stations(i)%weight
            stations(i)%q = stations(i)%q/stations(i)%weight
            stations(i)%pressure = stations(i)%pressure/stations(i)%weight
            stations(i)%temperature = stations(i)%temperature/stations(i)%weight
        end if
    end do ! i

!    ! count the different types of measurements
!    call find_unique_strings(stations(:)%station_type, total_types, unique_type_idx)
!    allocate(num_stations_by_type(total_types), type_gid(total_types), station_types(total_types))
!    do i_type = 1, total_types
!        station_types(i_type) = stations(unique_type_idx(i_type))%station_type
!        num_stations_by_type(i_type) = count(stations(:)%station_type == station_types(i_type))
!        do i = 1, num_stations
!            if (stations(i)%station_type == station_types(i_type)) stations(i)%type_id = i_type
!        end do
!    end do

!    do i = 1, num_stations
!        if (stations(i)%type_id < 0) then
!            write(gol,'(a," :: station ", a, "/", a, " has not been classified yet")') rname, trim(stations(i)%station_type), trim(stations(i)%station_id)
!            call goErr
!            status = 1
!            IF_NOTOK_RETURN(status=1)
!        end if
!    end do

    ! now open a netcdf file for writing, the filename is already correct at this point
    call check_dir(outfile_station)
    netcdf_ID = nc_open(outfile_station, 'c', status)
    IF_NOTOK_RETURN(status=1)

    ! write the starting and ending dates
    form_date = NewDate(time6 = idatei)
    io_status = nf90_put_att(netcdf_ID, NF90_GLOBAL, 'starting time', trim(Pretty(form_date))) ; CHECK_NCSTAT(io_status)
    form_date = NewDate(time6 = idatee)
    io_status = nf90_put_att(netcdf_ID, NF90_GLOBAL, 'ending time', trim(Pretty(form_date))) ; CHECK_NCSTAT(io_status)

    ! write the tracer names
    do itr = 1, ntracet
        write(tracer_names, '("tracer_", i3.3)') itr
        io_status = nf90_put_att(netcdf_ID, NF90_GLOBAL, trim(tracer_names), trim(names(itr))) ; CHECK_NCSTAT(io_status)
    end do
    !tracer_names = trim(names(1))

    ! make a string with unit names
    unit_names = mixrat_unit_name(1)
    do itr = 2, ntracet
        unit_names = trim(unit_names)//', '//mixrat_unit_name(itr)
    end do
    unit_names = 'mole fraction ('//trim(unit_names)//')'

    ! create groups and global dimensions
    call nc_create_dim(netcdf_ID, 'samples', NF90_UNLIMITED)
    call nc_create_dim(netcdf_ID, 'tracers', ntracet)
    call nc_create_dim(netcdf_ID, 'date_components', 6)

    call nc_dump_var(netcdf_ID, 'date_midpoints', (/'date_components', 'samples        '/), midpoint_dates(1:6,1:io_nrecords))

!    do i_type = 1, total_types
!        type_gid(i_type) = nc_create_group(netcdf_ID, station_types(i_type))
!    end do

    do i=1,num_stations

        ! create the groups
!        grp_id = nc_create_group(type_gid(stations(i)%type_id), stations(i)%station_id)
        write(st_group_name, '("station_", i5.5)') stations(i)%site_num
        grp_id = nc_create_group(netcdf_ID, trim(st_group_name))
        if (station_sample_meteo) meteo_grp_id = nc_create_group(grp_id, 'meteo')

        ! set group attributes
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'name', trim(stations(i)%station_name)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'latitude', stations(i)%lat) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'longitude', stations(i)%lon) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'altitude', stations(i)%alt) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'region', region_name(stations(i)%region)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'abbr', trim(stations(i)%station_type)//'/'//trim(stations(i)%station_id)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'code', trim(stations(i)%station_id)) ; CHECK_NCSTAT(io_status)

        ! write the variables, corresponding to mix, var and nsamples
        call nc_dump_var(grp_id, 'mixing_ratio', (/'samples','tracers'/), stations(i)%mix(1:io_nrecords,1:ntracet), (/'unit'/), (/trim(unit_names)/))
        call nc_dump_var(grp_id, 'mixing_ratio_err', (/'samples','tracers'/), sqrt(stations(i)%var(1:io_nrecords,1:ntracet)), (/'unit','type'/), (/ trim(unit_names), 'model representation error'/))
        call nc_dump_var(grp_id, 'nsamples', (/'samples'/), stations(i)%nsamples(1:io_nrecords), (/'unit'/), (/'number of samples within the dynamic timestep'/))
        if (station_sample_meteo) then
            ! write the variables corresponding to blh, q, pres, temp
            call nc_dump_var(meteo_grp_id, 'blh', (/'samples'/), stations(i)%blh(1:io_nrecords), (/'unit'/), (/'boundary layer height in meters'/))
            call nc_dump_var(meteo_grp_id, 'q', (/'samples'/), stations(i)%q(1:io_nrecords), (/'unit'/), (/'specific humidity, mass of water/mass of air'/))
            call nc_dump_var(meteo_grp_id, 'pres', (/'samples'/), stations(i)%pressure(1:io_nrecords), (/'unit'/), (/'pressure in Pascals'/))
            call nc_dump_var(meteo_grp_id, 'temp', (/'samples'/), stations(i)%temperature(1:io_nrecords), (/'unit'/), (/'temperature in Kelvins'/))
        end if

    end do ! i

    call nc_close(netcdf_ID)

    deallocate(midpoint_dates)

    do i = 1, num_stations

        deallocate(stations(i)%mix)
        deallocate(stations(i)%var)
        deallocate(stations(i)%nsamples)
        deallocate(stations(i)%weight)
        if (station_sample_meteo) then
            deallocate(stations(i)%blh)
            deallocate(stations(i)%q)
            deallocate(stations(i)%pressure)
            deallocate(stations(i)%temperature)
        end if

    end do

!    if (allocated(unique_type_idx)) deallocate(unique_type_idx)
!    if (allocated(num_stations_by_type)) deallocate(num_stations_by_type)
!    if (allocated(type_gid)) deallocate(type_gid)
!    if (allocated(station_types)) deallocate(station_types)

    status = 0

end subroutine user_output_station_write

subroutine user_output_station_sample(region, tr, i_tstep)

    use global_data,    only : mass_dat, region_dat, conv_dat
    use MeteoData,      only : gph_dat, m_dat, temper_dat, phlb_dat, humid_dat
    use chem_param,     only : fscale, ntracet
    use dims,           only : lm
    use Go,             only : TDate, rTotal, operator(-)

    implicit none

    ! input/output
    integer,intent(in)                  :: region, i_tstep
    type(TDate), intent(in)             :: tr(2)

    character(len=*), parameter         :: rname = mname//'/user_output_station_done'

    integer                             :: ifr,jfr,lfr,ifn,jfn,l,lfrn,itr,i,j,k,i_region,i_stat,lmr
    real                                :: alt, rlf, wcz, rmf, mean_mix, var_mix
    real                                :: wcx, wcy, rif, rjf, weight
    real,dimension(:,:,:), pointer      :: m,gph
    real,dimension(:,:,:,:), pointer    :: rm, rxm, rym, rzm
    real,dimension(:,:), pointer        :: blh ! boundary layer height [m]
    real,dimension(:,:,:), pointer      :: T ! temperature
    real,dimension(:,:,:), pointer      :: phlb ! pressure grid boundaries
    real,dimension(:,:,:), pointer      :: q ! specific humidity [kg/kg]
    real,dimension(0:lm(region))        :: height
    real                                :: mix_n(-1:1,-1:1,-1:1) ! static for OpenMP

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

    weight = rTotal(tr(2) - tr(1), 'sec') ! from emission_fwd.F90

    !$OMP PARALLEL private (mix_n, i, j, k, l, mean_mix, var_mix) &
    !$OMP private (i_stat, i_region, lfrn, wcz, height, rmf) &
    !$OMP private (alt, ifr, ifn, jfr, jfn, rif, rjf, wcx, wcy, lfr, rlf)
    !$OMP DO schedule(dynamic)
    do i_region = 1, stat_counter(region)%n_obs
        i_stat = stat_counter(region)%istat(i_region)

        alt = stations(i_stat)%alt
        ifr = stations(i_stat)%ifr
        ifn = stations(i_stat)%ifn
        jfr = stations(i_stat)%jfr
        jfn = stations(i_stat)%jfn
        rif = stations(i_stat)%rif
        rjf = stations(i_stat)%rjf
        wcx = stations(i_stat)%wcx
        wcy = stations(i_stat)%wcy

        ! interpolate the altitude to site position...
        lfr = 1 !layer

        do l=0,lm(region)
            height(l) = wcx       * wcy       * gph(ifr,jfr,l+1) + &
                        (1.0-wcx) * wcy       * gph(ifn,jfr,l+1) + &
                        wcx       * (1.0-wcy) * gph(ifr,jfn,l+1) + &
                        (1.0-wcx) * (1.0-wcy) * gph(ifn,jfn,l+1)
            if(l==0) stations(i_stat)%surface_height = height(0)
        end do

        do l = lowest_sample_layer, lm(region) ! selects layer, note that we skip the first layer and start from second layer from surface
            if (height(l) .gt. alt) exit
        end do

        select case(l)
            case(0) ! will never happen, because we started from l=2 above
                if (.not. stations(i_stat)%below_surface_warning) then
                    if (station_verbose) then
                        write (gol,'("WARNING:  For station ", a)') stations(i_stat)%station_name
                        call goPr
                        write (gol,'("  Sample altitude of ",f8.2,"m is below surface height of ",f8.2,"m.")') alt, height(0)
                        call goPr
                        write (gol,'("  Will sample at surface.")')
                        call goPr
                    end if
                    stations(i_stat)%below_surface_warning = .True.
                end if
                lfr = 1
                rlf = -0.5  !surface...
            case default
                lfr = l  !the site layer
                ! the offset from the center of the layer (-0.5--->+0.5)
                ! (interpolation is in (m))
                ! Here we need to distinguish two cases, (1) when alt is between height(l-1) and height(l), and
                ! (2) when alt is below height(l-1), but we are sampling in layer l because the site is too close to the ground.
                if (alt >= height(l-1)) then
                    rlf = (alt-height(l-1))/(height(l)-height(l-1)) - 0.5
                else
                    ! the altitude of the station is below height(l-1), but we're sampling in layer l=2 because we don't want to sample too close to the ground
                    rlf = 0.0
                end if
        end select
        stations(i_stat)%lfr = lfr

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

        if (lfrn < lowest_sample_layer) then
            lfrn = lowest_sample_layer
            wcz = 1.0 ! no interpolation in the z direction
        end if
        if (lfrn == lmr+1) then
            !=========================================================
            ! if vertical neighbor is lmr+1 (which does not exist)
            ! -> no interpolation
            !=========================================================
            lfrn = lmr ! no interpolation
            wcz = 1.0
        end if
        stations(i_stat)%rlf = rlf

        do itr = 1, ntracet

            select case (station_interpolation)
              case (POINT_INTERPOLATION_GRIDBOX)
                rmf = rm(ifr,jfr,lfr,itr) / m(ifr,jfr,lfr) * fscale(itr)

              case (POINT_INTERPOLATION_SLOPES)
                rmf = ( rm(ifr,jfr,lfr,itr) + 2.0*(rif*rxm(ifr,jfr,lfr,itr) + rjf*rym(ifr,jfr,lfr,itr) + rlf*rzm(ifr,jfr,lfr,itr)) ) &
                    / m(ifr,jfr,lfr) * fscale(itr)

              case (POINT_INTERPOLATION_LINEAR)

                rmf = ( &
                         wcx  *      wcy  *      wcz  * rm(ifr,jfr,lfr ,itr) / m(ifr,jfr,lfr )  + &
                    (1.0-wcx) *      wcy  *      wcz  * rm(ifn,jfr,lfr ,itr) / m(ifn,jfr,lfr )  + &
                         wcx  * (1.0-wcy) *      wcz  * rm(ifr,jfn,lfr ,itr) / m(ifr,jfn,lfr )  + &
                    (1.0-wcx) * (1.0-wcy) *      wcz  * rm(ifn,jfn,lfr ,itr) / m(ifn,jfn,lfr )  + &
                         wcx  *      wcy  * (1.0-wcz) * rm(ifr,jfr,lfrn,itr) / m(ifr,jfr,lfrn)  + &
                    (1.0-wcx) *      wcy  * (1.0-wcz) * rm(ifn,jfr,lfrn,itr) / m(ifn,jfr,lfrn)  + &
                         wcx  * (1.0-wcy) * (1.0-wcz) * rm(ifr,jfn,lfrn,itr) / m(ifr,jfn,lfrn)  + &
                    (1.0-wcx) * (1.0-wcy) * (1.0-wcz) * rm(ifn,jfn,lfrn,itr) / m(ifn,jfn,lfrn)) * fscale(itr)

              case default
                write (gol,'("unsupported station interpolation index ",i6)') station_interpolation; call goErr

            end select

            stations(i_stat)%mix(i_tstep,itr) = stations(i_stat)%mix(i_tstep,itr) + weight * rmf

            ! calculate the representation error
            select case (point_error_choice(itr))

                case (POINT_ERROR_NEIGHBOR)
                    do i = -1, 1
                        do j = -1, 1
                            do k = -1, 1
                                mix_n(i,j,k) = (rm(ifr,jfr,lfr,itr) + &
                                2.0*(i*rxm(ifr,jfr,lfr,itr) + j*rym(ifr,jfr,lfr,itr) + k*rzm(ifr,jfr,lfr,itr))) / &
                                m(ifr,jfr,lfr) * fscale(itr)
                            end do
                        end do
                    end do
                    mean_mix = sum(mix_n)/size(mix_n)
                    mix_n = mix_n - mean_mix
                    var_mix = sum(mix_n * mix_n)/size(mix_n)

                case (POINT_ERROR_GRADIENT)
                    var_mix = ((2.0*rxm(ifr,jfr,lfr,itr))**2 + (2.0*rym(ifr,jfr,lfr,itr))**2 + (2.0*rzm(ifr,jfr,lfr,itr))**2) / m(ifr,jfr,lfr)**2 * fscale(itr)**2

            end select
            stations(i_stat)%var(i_tstep,itr) = stations(i_stat)%var(i_tstep,itr) + weight * weight * var_mix
        end do ! itr

        ! sample meteo
        if(station_sample_meteo) then
            rmf = wcx       * wcy       * blh(ifr,jfr)  + &
                  (1.0-wcx) * wcy       * blh(ifn,jfr) + &
                  wcx       * (1.0-wcy) * blh(ifr,jfn) + &
                  (1.0-wcx) * (1.0-wcy) * blh(ifn,jfn)
            stations(i_stat)%blh(i_tstep) = stations(i_stat)%blh(i_tstep) + weight * rmf

            rmf =      wcx  *      wcy  *      wcz  * T(ifr,jfr,lfr)  + &
                  (1.0-wcx) *      wcy  *      wcz  * T(ifn,jfr,lfr)  + &
                       wcx  * (1.0-wcy) *      wcz  * T(ifr,jfn,lfr)  + &
                  (1.0-wcx) * (1.0-wcy) *      wcz  * T(ifn,jfn,lfr)  + &
                       wcx  *      wcy  * (1.0-wcz) * T(ifr,jfr,lfrn) + &
                  (1.0-wcx) *      wcy  * (1.0-wcz) * T(ifn,jfr,lfrn) + &
                       wcx  * (1.0-wcy) * (1.0-wcz) * T(ifr,jfn,lfrn) + &
                  (1.0-wcx) * (1.0-wcy) * (1.0-wcz) * T(ifn,jfn,lfrn)
            stations(i_stat)%temperature(i_tstep) = stations(i_stat)%temperature(i_tstep) + weight * rmf

            rmf =      wcx  *      wcy  *      wcz  * q(ifr,jfr,lfr)  + &
                  (1.0-wcx) *      wcy  *      wcz  * q(ifn,jfr,lfr)  + &
                       wcx  * (1.0-wcy) *      wcz  * q(ifr,jfn,lfr)  + &
                  (1.0-wcx) * (1.0-wcy) *      wcz  * q(ifn,jfn,lfr)  + &
                       wcx  *      wcy  * (1.0-wcz) * q(ifr,jfr,lfrn) + &
                  (1.0-wcx) *      wcy  * (1.0-wcz) * q(ifn,jfr,lfrn) + &
                       wcx  * (1.0-wcy) * (1.0-wcz) * q(ifr,jfn,lfrn) + &
                  (1.0-wcx) * (1.0-wcy) * (1.0-wcz) * q(ifn,jfn,lfrn)
            stations(i_stat)%q(i_tstep) = stations(i_stat)%q(i_tstep) + weight * rmf

            rmf = (((0.5-rlf) * phlb(ifr,jfr,lfr) + (0.5+rlf) * phlb(ifr,jfr,lfrn)) *      wcx  *      wcy  + &
                   ((0.5-rlf) * phlb(ifn,jfr,lfr) + (0.5+rlf) * phlb(ifn,jfr,lfrn)) * (1.0-wcx) *      wcy  + &
                   ((0.5-rlf) * phlb(ifr,jfn,lfr) + (0.5+rlf) * phlb(ifr,jfn,lfrn)) *      wcx  * (1.0-wcy) + &
                   ((0.5-rlf) * phlb(ifn,jfn,lfr) + (0.5+rlf) * phlb(ifn,jfn,lfrn)) * (1.0-wcx) * (1.0-wcy))
            stations(i_stat)%pressure(i_tstep) = stations(i_stat)%pressure(i_tstep) + weight * rmf
        end if ! station_sample_meteo

        stations(i_stat)%nsamples(i_tstep) = stations(i_stat)%nsamples(i_tstep) + 1
        stations(i_stat)%weight(i_tstep)   = stations(i_stat)%weight(i_tstep)   + weight

    end do ! i_region
    !$OMP END DO
    !$OMP END PARALLEL

    nullify(m, rm, rxm, rym, rzm)
    nullify(gph, T, phlb, blh, q)

end subroutine user_output_station_sample

subroutine user_output_station_step(region, tr, status)

    use datetime,       only : date2tau, tau2date
    use dims,           only : itaur, ndyn, tref, idatei, idatee, ndyn_max, newsrun
    use datetime_for,   only : SEC_PER_DAY
    use go,             only : Get, TDate, rTotal, operator(+), operator(-), operator(/), Pretty

    implicit none

    integer, intent(in)         :: region
    type(TDate), intent(in)     :: tr(2)
    integer, intent(out)        :: status

    character(len=*), parameter :: rname = mname//'/user_output_station_step'

    integer                     :: i_tstep, midpt_time, i_period, tau_beg, tau_mid, secs_in_window
    integer, dimension(6)       :: midpt_date, idate_temp
    type(TDate)                 :: tmid

    ! Which period number are we in?
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call Get(tmid, time6=midpt_date)

    select case (split_period)
        case ('a')
            i_period = 1
        case ('y')
            i_period = midpt_date(1)-idatei(1)+1
        case ('m')
            i_period = (midpt_date(1)-idatei(1))*12 + (midpt_date(2)-idatei(2)+1)
        case ('d')
            idate_temp = idatei
            idate_temp(4:6) = 0
            call date2tau(midpt_date, tau_mid)
            call date2tau(idate_temp, tau_beg)
            i_period = (tau_mid - tau_beg)/SEC_PER_DAY + 1
    end select

    ! If this is a new period but not a new run, we need to write out the
    ! previous period's file and re-initialize the stations data structure
    if (.not. period_done(i_period)) then
        ! period_done is a misnomer, because it doesn't mean that we are done with this period, but instead it simply
        ! means that we are done with writing the data from the *previous* period
        if (.not. newsrun) call user_output_station_write(status) ! write the previous period's data
        call init_station_struct(period_bounds(i_period,1), period_bounds(i_period,2), status)
        period_done(i_period) = .true.
    end if

    ! Within this period, which timestep should we fill?
    secs_in_window = nint(rTotal(tmid - period_bounds(i_period,1), 'sec'))
    i_tstep = 1 + secs_in_window/ndyn_max

    call user_output_station_sample(region, tr, i_tstep)

    status = 0

end subroutine user_output_station_step

function in_region(lat,lon,region)

    use dims, only : dx, dy, xref, yref, xbeg, ybeg, xend, yend, parent

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
        .and. (lat .gt. ybeg(region)+dyp*parent_sample_choice) .and. (lat .lt. yend(region)-dyp*parent_sample_choice)) in_region = .true.

end function in_region

end module user_output_station
