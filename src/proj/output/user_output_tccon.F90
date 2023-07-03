!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#define CHECK_NCSTAT(nc_ret_code) if (nc_ret_code /= NF90_NOERR) then; status=1; IF_NOTOK_RETURN(status=1); end if
!
!###############################################################################
#include "tm5.inc"

module user_output_tccon

use go_date,    only : TDate
use GO,         only : gol, goErr, goPr

use user_output_satellite_data,   only : SAT_INTERPOLATION_GRIDBOX, SAT_INTERPOLATION_SLOPES, SAT_INTERPOLATION_LINEAR

implicit none

private

public :: user_output_tccon_init, user_output_tccon_step, user_output_tccon_done

character(len=*), parameter :: mname = 'user_output_tccon'
character(len=300)          :: stationlist_filename, outdir_tccon, outfile_tccon
logical                     :: station_verbose, station_sample_meteo
integer                     :: tccon_interpolation
logical                     :: sample_in_parent, calc_profile_var
integer                     :: parent_sample_choice = 0
character(len=80)           :: station_error_choice
integer                     :: num_stations ! number of stations
integer                     :: io_nrecords ! number of dynamic timesteps during the run
integer                     :: netcdf_ID ! ID of the output file
integer(4), allocatable     :: midpoint_dates(:,:)
character(len=1)            :: split_period
logical, allocatable        :: period_done(:)
type(TDate), allocatable    :: period_bounds(:,:)

type T_tccon
    character(len=12)       :: station_id     ! station identifier
    character(len=60)       :: station_name      ! station name
    real                    :: lat       ! latitude of station
    real                    :: lon       ! longitude of station
    integer                 :: region ! region of station
    integer(4), allocatable :: nsamples(:) ! number of samples accumulated
    real, allocatable       :: mix(:,:,:), var(:,:,:) ! mixing ratio and variance of mixing ratio
    integer                 :: ifr, jfr ! i,j region indices for flask's grid cell
    integer                 :: ifn, jfn ! i,j region indices for flask's "next" grid cell
    integer                 :: rif, rjf ! fractions from center of ifr,jfr box
    real                    :: wcx, wcy ! x and y weighting factors for slopes interpolation
    real, allocatable       :: q(:,:), pressure(:,:), temperature(:,:), gph(:,:) ! meteorological variables
    real, allocatable       :: weight(:)    ! weight of each time step, depends on length of step
end type T_tccon

type T_tccon_region_counter
    integer, allocatable :: istat(:) ! for a region, stores the station indices corresponding to that region
    integer :: n_obs
end type T_tccon_region_counter

type(T_tccon), allocatable :: stations(:)
type(T_tccon_region_counter), allocatable :: stat_counter(:)

contains

subroutine init_tccon_struct(start_date, end_date, status)

    use go_date,    only : TIncrDate, iTotal, operator(-), operator(+), IncrDate, Get
    use dims,       only : ndyn_max, lm
    use chem_param, only : ntracet

    implicit none

    type(TDate), intent(in) :: start_date, end_date
    integer, intent(out)    :: status

    character(len=*), parameter :: rname = mname//'/init_tccon_struct'

    type(TDate)             :: dummy_t
    type(TIncrDate)         :: dt
    integer                 :: total_seconds_in_interval, idate_temp(6), i, region

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
        write(outfile_tccon, '(a,"/tccon/tccon.nc4")') trim(outdir_tccon)
    case ('m')
        write(outfile_tccon, '(a,"/tccon/tccon_",i4.4,i2.2,".nc4")') trim(outdir_tccon), idate_temp(1:2)
    case ('d')
        write(outfile_tccon, '(a,"/tccon/tccon_",i4.4,2i2.2,".nc4")') trim(outdir_tccon), idate_temp(1:3)
    case default
        write(0,*) 'Wrong split period selected'
        IF_NOTOK_RETURN(status=1)
    end select

    do i = 1, num_stations

        region = stations(i)%region

        allocate(stations(i)%mix(lm(region),io_nrecords,ntracet))

        allocate(stations(i)%nsamples(io_nrecords))
        allocate(stations(i)%pressure(lm(region)+1,io_nrecords))
        allocate(stations(i)%weight(io_nrecords))
        stations(i)%mix         = 0.0
        stations(i)%nsamples    = 0
        stations(i)%pressure    = 0.0
        stations(i)%weight      = 0.0
        if (calc_profile_var) then
            allocate(stations(i)%var(lm(region),io_nrecords,ntracet))
            stations(i)%var         = 0.0
        end if

        if (station_sample_meteo) then
            allocate(stations(i)%q(lm(region),io_nrecords))
            allocate(stations(i)%temperature(lm(region),io_nrecords))
            allocate(stations(i)%gph(lm(region)+1,io_nrecords))
            stations(i)%q           = 0.0
            stations(i)%temperature = 0.0
            stations(i)%gph         = 0.0
        end if

    end do

    status = 0

end subroutine init_tccon_struct

subroutine user_output_tccon_init(status)

    use GO,             only : ReadRc
    use go_date,        only : NewDate, Get_End_Of
    use global_data,    only : rcF, outdir
    use dims,           only : im, jm, lm, dx, dy, xref, yref, xbeg, ybeg, xend, yend, tref
    use dims,           only : nregions, region_name, xcyc, ndyn, ndyn_max, idatee, idatei
    use chem_param,     only : ntrace, ntracet
    use datetime,       only : date2tau, tau2date
    use orderpack,      only : mrgrnk
    use datetime_for,   only : SEC_PER_DAY

    implicit none

    ! --- in/out ---------------------------------
    integer, intent(out)    :: status

    ! --- local ----------------------------------
    integer                 :: sunit = 111, n_period, tau_beg, tau_end
    integer                 :: i, region, i_region, idate_temp(6)
    character(len=300)      :: dummy_str
    integer, allocatable    :: region_rank(:), region_counter(:)
    real                    :: dxr, dyr
    type(TDate)             :: cur_date, next_date

    ! --- const ------------------------------
    character(len=*), parameter         :: rname = mname//'/user_output_tccon_init'

    call ReadRc( rcF, 'input.tccon.filename', stationlist_filename, status) ! must be the full path
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.tccon.verbose', station_verbose, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.tccon.meteo', station_sample_meteo, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.satellite.errors', station_error_choice, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.satellite.interpolation', tccon_interpolation, status, default=SAT_INTERPOLATION_LINEAR)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'output.satellite.sample.parent', sample_in_parent, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.dir', outdir_tccon, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.satellite.split.period', split_period, status, default='a')
    IF_NOTOK_RETURN(status=1)

    if (sample_in_parent) parent_sample_choice = 1

    ! if station_error_choice is '2d', do not try to write or fill the profile variance
    select case (trim(adjustl(station_error_choice)))
        case ('neighbors', 'gradient')
            calc_profile_var = .true.
        case default
            calc_profile_var = .false.
    end select

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
        read(sunit, '(a3,2f8.2,1x,a60)') stations(i)%station_id, stations(i)%lat, stations(i)%lon, stations(i)%station_name
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

        dxr = dx/xref(region)
        dyr = dy/yref(region)

        stations(i)%rif = (stations(i)%lon-float(xbeg(region)))/dxr + 0.99999
        stations(i)%rjf = (stations(i)%lat-float(ybeg(region)))/dyr + 0.99999

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
        end if

        !the neighbour for y interpolation
        if(stations(i)%rjf .gt. 0) then
          stations(i)%jfn = stations(i)%jfr+1
        else
          stations(i)%jfn = stations(i)%jfr-1
        end if

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
        end if

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
    if (sum(stat_counter(:)%n_obs) /= num_stations) write(*,'(a,a,a)') 'WARNING :: ', rname, ' :: not all tccon stations accounted for'

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

    allocate(period_done(n_period))
    period_done = .false.

    ! create arrays with starting and ending times per period
    allocate(period_bounds(n_period,2))
    select case (split_period)
    case ('a')
        period_bounds(1,1) = NewDate(time6=idatei)
        period_bounds(1,2) = NewDate(time6=idatee)
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

end subroutine user_output_tccon_init

subroutine user_output_tccon_write(status)

    use dims,       only : region_name, idatei, idatee, lm
    use chem_param, only : ntracet, names, mixrat_unit_name, tracer_name_len
    use go_date,    only : Pretty, NewDate
    use misctools,  only : check_dir
    use file_netcdf
    use netcdf

    implicit none

    integer, intent(out)    :: status

    character(len=*), parameter :: rname = mname//'/user_output_tccon_write'

    integer                     :: io_status, grp_id, meteo_grp_id, i, int_date(6), levs, l, itr
    type(TDate)                 :: form_date
    character(len=256)          :: unit_names
    character(len=256)          :: tracer_names

    ! the number of levels should be same for all stations
    levs = lm(1)

    ! first, evaluate the tccon data
    !$omp parallel do schedule(dynamic) private(i,l,itr)
    do i=1,num_stations
        do l=1,levs
            do itr = 1, ntracet
                stations(i)%mix(l,:,itr) = stations(i)%mix(l,:,itr)/stations(i)%weight
                if (calc_profile_var) stations(i)%var(l,:,itr) = stations(i)%var(l,:,itr)/stations(i)%weight/stations(i)%weight ! variance needs to be divided by weight**2
            end do ! itr
        end do ! l
        do l=1,levs+1
            stations(i)%pressure(l,:) = stations(i)%pressure(l,:)/stations(i)%weight
        end do ! l
        if (station_sample_meteo) then
            do l=1,levs
                stations(i)%q(l,:) = stations(i)%q(l,:)/stations(i)%weight
                stations(i)%temperature(l,:) = stations(i)%temperature(l,:)/stations(i)%weight
            end do ! l
            do l=1,levs+1
                stations(i)%gph(l,:) = stations(i)%gph(l,:)/stations(i)%weight
            end do ! l
        end if
    end do ! i
    !$omp end parallel do

    ! now open a netcdf file for writing
    call check_dir(outfile_tccon)
    netcdf_ID = nc_open(outfile_tccon, 'c', status)
    IF_NOTOK_RETURN(status=1)

    ! write the tracer names
    tracer_names = trim(names(1))
    do itr = 2, ntracet
        tracer_names = tracer_names//','//trim(names(itr))
    end do
    call nc_set_attrs(netcdf_ID, 'tracer_names', trim(tracer_names))

    ! make a string with unit names
    unit_names = mixrat_unit_name(1)
    do itr = 2, ntracet
        unit_names = trim(unit_names)//', '//mixrat_unit_name(itr)
    end do
    unit_names = 'mole fraction ('//trim(unit_names)//')'

    ! write the starting and ending dates
    form_date = NewDate(time6 = idatei)
    io_status = nf90_put_att(netcdf_ID, NF90_GLOBAL, 'starting time', trim(Pretty(form_date))) ; CHECK_NCSTAT(io_status)
    form_date = NewDate(time6 = idatee)
    io_status = nf90_put_att(netcdf_ID, NF90_GLOBAL, 'ending time', trim(Pretty(form_date))) ; CHECK_NCSTAT(io_status)

    ! create groups and global dimensions
    call nc_create_dim(netcdf_ID, 'samples', NF90_UNLIMITED)
    call nc_create_dim(netcdf_ID, 'tracers', ntracet)
    call nc_create_dim(netcdf_ID, 'date_components', 6)
    call nc_create_dim(netcdf_ID, 'levels', levs)
    call nc_create_dim(netcdf_ID, 'nlev_p1', levs+1)

    call nc_dump_var(netcdf_ID, 'date_midpoints', (/'date_components', 'samples        '/), midpoint_dates)

    do i=1,num_stations

        ! create the station group
        grp_id = nc_create_group(netcdf_ID, stations(i)%station_id)
        if (station_sample_meteo) meteo_grp_id = nc_create_group(grp_id, 'meteo')

        ! set group attributes
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'name', trim(stations(i)%station_name)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'latitude', stations(i)%lat) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'longitude', stations(i)%lon) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(grp_id, NF90_GLOBAL, 'region', region_name(stations(i)%region)) ; CHECK_NCSTAT(io_status)

        ! write the variables, corresponding to mix, var, pressure boundaries and nsamples
        call nc_dump_var(grp_id, 'mixing_ratio',     (/'levels ','samples','tracers'/), stations(i)%mix,        (/'unit'/), (/ trim(unit_names) /))
        if (calc_profile_var) call nc_dump_var(grp_id, 'mixing_ratio_err', (/'levels ','samples','tracers'/), sqrt(stations(i)%var),  (/'unit','type'/), (/ trim(unit_names), 'model representation error'/))
        call nc_dump_var(grp_id, 'pres',             (/'nlev_p1','samples'/),           stations(i)%pressure,   (/'unit'/), (/'pressure at layer boundaries in Pascals'/))
        call nc_dump_var(grp_id, 'nsamples',         (/'samples'/),                     stations(i)%nsamples,   (/'unit'/), (/'number of samples within the maximum dynamic timestep'/))
        if (station_sample_meteo) then
            ! write the variables corresponding to q, gph, temp
            call nc_dump_var(meteo_grp_id, 'q',     (/'levels ','samples'/), stations(i)%q,             (/'unit'/), (/'specific humidity, mass of water/mass of air'/))
            call nc_dump_var(meteo_grp_id, 'gph',   (/'nlev_p1','samples'/), stations(i)%gph,           (/'unit'/), (/'geopotential height at layer boundaries in meters'/))
            call nc_dump_var(meteo_grp_id, 'temp',  (/'levels ','samples'/), stations(i)%temperature,   (/'unit'/), (/'temperature in Kelvins'/))
        end if

    end do ! i

    call nc_close(netcdf_ID)

    deallocate(midpoint_dates)

    do i = 1, num_stations

        deallocate(stations(i)%mix)
        if (calc_profile_var) deallocate(stations(i)%var)
        deallocate(stations(i)%nsamples)
        deallocate(stations(i)%pressure)
        deallocate(stations(i)%weight)

        if (station_sample_meteo) then
            deallocate(stations(i)%q)
            deallocate(stations(i)%temperature)
            deallocate(stations(i)%gph)
        end if

    end do
    status = 0

end subroutine user_output_tccon_write

subroutine user_output_tccon_done(status)

    implicit none

    character(len=*), parameter :: rname = mname//'/user_output_tccon_done'

    integer, intent(out)    :: status

    call user_output_tccon_write(status)
    IF_NOTOK_RETURN(status=1)

    deallocate(stations)
    deallocate(period_done, period_bounds)

    status = 0

end subroutine user_output_tccon_done

subroutine user_output_tccon_sample(region, tr, i_tstep)

    use dims,           only : lm
    use global_data,    only : mass_dat
    use MeteoData,      only : m_dat, phlb_dat, gph_dat, humid_dat, temper_dat
    use chem_param,     only : fscale, ntracet
    use Go,             only : TDate, rTotal, operator(-)

    implicit none

    ! input/output
    integer,intent(in)      :: region, i_tstep
    type(TDate), intent(in) :: tr(2)

    character(len=*), parameter :: rname = mname//'/user_output_tccon_sample'

    ! local
    real, pointer       :: m(:,:,:), T(:,:,:), phlb(:,:,:), q(:,:,:), gph(:,:,:)
    real, pointer       :: rm(:,:,:,:), rxm(:,:,:,:), rym(:,:,:,:)
    real, pointer       :: p(:,:) ! surface pressure
    integer, parameter  :: lm_fixed = lm(1)
    real                :: profile(-1:1,-1:1,lm_fixed)
    real                :: rmf(lm_fixed), rmf_p1(lm_fixed), mean_mix(lm_fixed), var_mix(lm_fixed)
    integer, parameter  :: num_neighbors = 9
    logical             :: in_window
    integer             :: i, itr, j, k, i_region, ifr, jfr, lfr, ifn, jfn
    integer             :: midpt_time, i_stat
    real                :: wcx, wcy, rif, rjf, pres, weight

    ! pointers to global arrays
    m       => m_dat(region)%data
    rm      => mass_dat(region)%rm_t
    rxm     => mass_dat(region)%rxm_t
    rym     => mass_dat(region)%rym_t
    t       => temper_dat(region)%data
    phlb    => phlb_dat(region)%data
    q       => humid_dat(region)%data
    gph     => gph_dat(region)%data

    weight = rTotal(tr(2) - tr(1), 'sec') ! from emission_fwd.F90

    do itr = 1, ntracet
        !$omp parallel do private (i_region, i_stat, ifr, ifn, jfr, jfn, rif, rjf, wcx, wcy) &
        !$omp private (rmf, rmf_p1, i, j, profile, mean_mix, var_mix) schedule(dynamic)
        do i_region = 1, stat_counter(region)%n_obs
            i_stat = stat_counter(region)%istat(i_region)

            ifr = stations(i_stat)%ifr
            ifn = stations(i_stat)%ifn
            jfr = stations(i_stat)%jfr
            jfn = stations(i_stat)%jfn
            rif = stations(i_stat)%rif
            rjf = stations(i_stat)%rjf
            wcx = stations(i_stat)%wcx
            wcy = stations(i_stat)%wcy

            select case (tccon_interpolation)
                case (SAT_INTERPOLATION_GRIDBOX)
                    rmf = rm(ifr,jfr,:,itr) / m(ifr,jfr,:) * fscale(itr)

                case (SAT_INTERPOLATION_SLOPES)
                    rmf = ( rm(ifr,jfr,:,itr) + 2.0*(rif*rxm(ifr,jfr,:,itr) + rjf*rym(ifr,jfr,:,itr)) ) &
                        / m(ifr,jfr,:) * fscale(itr)

                case (SAT_INTERPOLATION_LINEAR)
                    rmf = ( &
                             wcx  *      wcy  * rm(ifr,jfr,:,itr) / m(ifr,jfr,:)  + &
                        (1.0-wcx) *      wcy  * rm(ifn,jfr,:,itr) / m(ifn,jfr,:)  + &
                             wcx  * (1.0-wcy) * rm(ifr,jfn,:,itr) / m(ifr,jfn,:)  + &
                        (1.0-wcx) * (1.0-wcy) * rm(ifn,jfn,:,itr) / m(ifn,jfn,:)) * fscale(itr)

                case default
                    write (gol,'("unsupported TCCON interpolation index ",i6)') tccon_interpolation; call goErr

            end select

            stations(i_stat)%mix(:,i_tstep,itr) = stations(i_stat)%mix(:,i_tstep,itr) + weight * rmf

            ! calculate the representation error
            select case (trim(adjustl(station_error_choice)))

                case ('neighbors')
                        do i = -1, 1
                            do j = -1, 1
                                profile(i,j,:) = (rm(ifr,jfr,:,itr) + &
                                    2.0*(i*rxm(ifr,jfr,:,itr) + j*rym(ifr,jfr,:,itr))) / &
                                    m(ifr,jfr,:) * fscale(itr)
                            end do ! j
                        end do ! i
                        mean_mix = sum(sum(profile, 1), 1)/num_neighbors
                        do i = -1, 1
                            do j = -1, 1
                                profile(i,j,:) = profile(i,j,:) - mean_mix
                            end do ! j
                        end do ! i
                        var_mix = sum(sum(profile * profile, 1), 1)/num_neighbors

                    case ('gradient')
                        var_mix = ((2.0*rxm(ifr,jfr,:,itr))**2 + (2.0*rym(ifr,jfr,:,itr))**2) / m(ifr,jfr,:)**2 * fscale(itr)**2

            end select

            if (calc_profile_var) stations(i_stat)%var(:,i_tstep,itr) = stations(i_stat)%var(:,i_tstep,itr) + weight * weight * var_mix

            rmf_p1 = wcx*wcy*phlb(ifr,jfr,:) + (1.0-wcx)*wcy*phlb(ifn,jfr,:) + wcx*(1.0-wcy)*phlb(ifr,jfn,:) + (1.0-wcx)*(1.0-wcy)*phlb(ifn,jfn,:)
            stations(i_stat)%pressure(:,i_tstep) = stations(i_stat)%pressure(:,i_tstep) + weight * rmf_p1

            ! sample meteo
            if(station_sample_meteo) then

                rmf =      wcx  *      wcy  * T(ifr,jfr,:)  + &
                      (1.0-wcx) *      wcy  * T(ifn,jfr,:)  + &
                           wcx  * (1.0-wcy) * T(ifr,jfn,:)  + &
                      (1.0-wcx) * (1.0-wcy) * T(ifn,jfn,:)
                stations(i_stat)%temperature(:,i_tstep) = stations(i_stat)%temperature(:,i_tstep) + weight * rmf

                rmf =      wcx  *      wcy  * q(ifr,jfr,:)  + &
                      (1.0-wcx) *      wcy  * q(ifn,jfr,:)  + &
                           wcx  * (1.0-wcy) * q(ifr,jfn,:)  + &
                      (1.0-wcx) * (1.0-wcy) * q(ifn,jfn,:)
                stations(i_stat)%q(:,i_tstep) = stations(i_stat)%q(:,i_tstep) + weight * rmf

                rmf_p1 = wcx*wcy*gph(ifr,jfr,:) + (1.0-wcx)*wcy*gph(ifn,jfr,:) + wcx*(1.0-wcy)*gph(ifr,jfn,:) + (1.0-wcx)*(1.0-wcy)*gph(ifn,jfn,:)
                stations(i_stat)%gph(:,i_tstep) = stations(i_stat)%gph(:,i_tstep) + weight * rmf_p1

            end if ! sample_meteo

            stations(i_stat)%nsamples(i_tstep)  = stations(i_stat)%nsamples(i_tstep) + 1
            stations(i_stat)%weight(i_tstep)    = stations(i_stat)%weight(i_tstep) + weight

        end do ! i_region
        !$omp end parallel do
    end do ! itr

    nullify(m, rm, rxm, rym)
    nullify(T, phlb, q, gph)

end subroutine user_output_tccon_sample

subroutine user_output_tccon_step(region, tr, status)

    use dims,           only : idatei, newsrun, ndyn_max
    use datetime_for,   only : SEC_PER_DAY
    use datetime,       only : date2tau, time_window
    use go,             only : Get, TDate, rTotal, operator(+), operator(-), operator(/)

    implicit none

    character(len=*), parameter         :: rname = mname//'/user_output_tccon_step'

    ! input/output
    integer, intent(in)     :: region
    type(Tdate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    integer                 :: i_tstep, i_period, tau_beg, midpt_date(6), idate_temp(6), tau_mid, secs_in_window
    type(Tdate)             :: tmid

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

    ! If this is a new period but not a new run, we need to write out the
    ! previous period's file and re-initialize the stations data structure
    if (.not. period_done(i_period)) then
        if (.not. newsrun) then 
            call user_output_tccon_write(status) ! write the previous period's data
            IF_NOTOK_RETURN(status=1)
        end if
        call init_tccon_struct(period_bounds(i_period,1), period_bounds(i_period,2), status)
        period_done(i_period) = .true.
    end if

    ! Within this period, which timestep should we fill? The time resolution is ndyn_max
    secs_in_window = nint(rTotal(tmid - period_bounds(i_period,1), 'sec'))
    i_tstep = 1 + secs_in_window/ndyn_max

    call user_output_tccon_sample(region, tr, i_tstep)

    status = 0

end subroutine user_output_tccon_step

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

end module user_output_tccon
