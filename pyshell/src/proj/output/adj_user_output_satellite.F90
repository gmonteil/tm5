!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module adj_user_output_satellite

  ! module included to read model-satellite (=departure) mismatch output fields
  ! called from adj_user_output

  use dims,                         only : nregions
  use toolbox,                      only : escape_tm
  use go,                           only : gol, goErr
  use user_output_satellite_data,   only : sat_region_forcing
  use user_output_satellite_data,   only : SAT_INTERPOLATION_GRIDBOX, SAT_INTERPOLATION_SLOPES, SAT_INTERPOLATION_LINEAR
  use file_netcdf

  implicit none

  ! interface

  private

  public :: adj_user_output_satellite_init
  public :: adj_user_output_satellite_step
  public :: adj_user_output_satellite_done

  type(sat_region_forcing), allocatable :: sat_dep(:)
  logical, allocatable, dimension(:) :: period_exists ! keeps track of which months have departure files
  logical, allocatable, dimension(:,:) :: sat_data_read ! keeps track of which period/region departures have already been read, replaces isOpen of Pim
  logical :: output_satellite_verbose
!  logical :: sample_dry_air = .false. ! this will replicate carbontracker sampling
  character(len=1) :: split_period ! 'd' if there is one file per day, 'm' if one per month
  logical, pointer :: new_period ! points to newday if there is one file per day, newmonth if there is one file per month
  character(len=200) :: outdir_satellite
  integer :: satellite_interpolation

  character(len=*), parameter       ::  mname = 'adj_user_output_satellite'

contains

subroutine adj_user_output_satellite_init(status)

    use dims,           only : nregions, idatei, idatee, newmonth, newday
    use GO,             only : ReadRc
    use datetime,       only : date2tau
    use global_data,    only : rcF
    use datetime_for,   only : SEC_PER_DAY
    use chem_param,     only : ntracet

    implicit none

    integer, intent(out)    :: status

    integer :: region, n_period, idate_temp(6), tau_beg, tau_end
    character(len=*), parameter         :: rname = mname//'/adj_user_output_satellite_init'

    call ReadRc(rcF, 'output.satellite.interpolation', satellite_interpolation, status, default=SAT_INTERPOLATION_LINEAR)
    IF_ERROR_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.verbose', output_satellite_verbose, status, default=.false.)
    IF_ERROR_RETURN(status=1)
!    call ReadRc(rcF, 'output.satellite.dryair.mixratio', sample_dry_air, status)
!    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.dir', outdir_satellite, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.satellite.split.period', split_period, status, default='m')
    IF_NOTOK_RETURN(status=1)

    allocate(sat_dep(nregions))
    do region = 1, nregions
        allocate(sat_dep(region)%io_tracer(ntracet))
    end do

    ! during adjoint run, idatei > idatee
    select case (split_period)
        case ('m')
            n_period = (idatei(1)-idatee(1))*12 + (idatei(2)-idatee(2)+1)
            new_period => newmonth
        case ('d')
            idate_temp = idatei
            idate_temp(4:6) = (/ 0,0,0 /)
            call date2tau(idate_temp, tau_beg)
            idate_temp = idatee
            idate_temp(4:6) = (/ 0,0,0 /)
            call date2tau(idate_temp, tau_end)
            ! tau_beg > tau_end
            n_period = (tau_beg-tau_end)/SEC_PER_DAY + 1
            new_period => newday
        case default
            write(0,*) 'Wrong split period selected'
            IF_NOTOK_RETURN(status=1)
    end select

    allocate(period_exists(n_period))
    ! assume departures exist for all periods
    period_exists = .true.

    allocate(sat_data_read(n_period,nregions))
    ! no file has yet been read it
    sat_data_read = .false.

    status = 0

end subroutine adj_user_output_satellite_init

subroutine adj_user_output_satellite_step(region, tr, status)
    ! called every timestep

    use dims,           only : idatee
    use datetime,       only : tau2date, date2tau
    use datetime_for,   only : SEC_PER_DAY
    use Go,             only : TDate, Get, operator(+), operator(-), operator(/)

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    integer                     :: i_period, midpt_date(6), tau_beg, idate_temp(6), tau_mid
    character(len=*), parameter :: rname = mname//'/adj_output_satellitedata'
    type(TDate)                 :: tmid

    !__START_SUBROUTINE______________________________________________________

    ! during the adjoint run, tr(1) > tr(2)
    ! if the run period is Jan 1 to Jan 3, then idatei => Jan 3, idatee => Jan 1
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call Get(tmid, time6=midpt_date)

    select case (split_period)
        case ('m')
            i_period = (midpt_date(1)-idatee(1))*12 + (midpt_date(2)-idatee(2)+1)
        case ('d')
            idate_temp = idatee
            idate_temp(4:6) = 0
            call date2tau(idate_temp, tau_beg)
            call date2tau(midpt_date, tau_mid)
            i_period = (tau_mid - tau_beg)/SEC_PER_DAY + 1
    end select

    if (.not. sat_data_read(i_period,region)) call read_samples(region, midpt_date, i_period, status)
    IF_NOTOK_RETURN(status=1)

    if (period_exists(i_period)) call adj_user_output_satellite_addforcing(region, tr, status)
    IF_NOTOK_RETURN(status=1)

    status = 0

end subroutine adj_user_output_satellite_step

subroutine read_samples(region, midpt_date, i_period, status)

    use dims,        only : itaur, ndyn, tref, region_name, dx, dy, xref, yref, nregions
    use dims,        only : xbeg, xend, ybeg, yend, itaue, newsrun, ndyn_max
    use dims,        only : xcyc, im, jm
    use dims,        only : isr, ier, jsr, jer
    use datetime,    only : tau2date, date2tau
    use toolbox,     only : escape_tm
    use chem_param,  only : names, ntracet
    use Go,          only : NewDate, TDate, Pretty
    ! For writing status messages
    use dims,        only : kmain, itau
    use datetime,    only : tstamp

    implicit none

    !__IO___________________________________________________________________
    integer, intent(in)    :: region, midpt_date(6), i_period
    integer, intent(out)   :: status

    !__CONST________________________________________________________________
    character(len=*), parameter         :: rname = mname//'/read_samples'

    !__LOCAL_VARIABLES______________________________________________________
    integer             :: nc_id, group_id, tgrp_id, io_status
    integer             :: i, itrac, i_obs, idate_temp(6), nobs
    character(len=256)  :: fname
    logical             :: file_exist
    real                :: dxr, dyr
    type(TDate)         :: sampleTime, startTime, endTime

    ! For writing status messages
    character(len=256)      :: write_string

    !__START_SUBROUTINE______________________________________________________
    if (.not. newsrun .and. new_period) then
        ! deallocate arrays from last period
        do itrac = 1, ntracet
            if (allocated(sat_dep(region)%io_tracer(itrac)%rif))    deallocate(sat_dep(region)%io_tracer(itrac)%rif)
            if (allocated(sat_dep(region)%io_tracer(itrac)%rjf))    deallocate(sat_dep(region)%io_tracer(itrac)%rjf)
            if (allocated(sat_dep(region)%io_tracer(itrac)%ifr))    deallocate(sat_dep(region)%io_tracer(itrac)%ifr)
            if (allocated(sat_dep(region)%io_tracer(itrac)%jfr))    deallocate(sat_dep(region)%io_tracer(itrac)%jfr)
            if (allocated(sat_dep(region)%io_tracer(itrac)%ifn))    deallocate(sat_dep(region)%io_tracer(itrac)%ifn)
            if (allocated(sat_dep(region)%io_tracer(itrac)%jfn))    deallocate(sat_dep(region)%io_tracer(itrac)%jfn)
            if (allocated(sat_dep(region)%io_tracer(itrac)%wcx))    deallocate(sat_dep(region)%io_tracer(itrac)%wcx)
            if (allocated(sat_dep(region)%io_tracer(itrac)%wcy))    deallocate(sat_dep(region)%io_tracer(itrac)%wcy)
            if (allocated(sat_dep(region)%io_tracer(itrac)%lat))    deallocate(sat_dep(region)%io_tracer(itrac)%lat)
            if (allocated(sat_dep(region)%io_tracer(itrac)%lon))    deallocate(sat_dep(region)%io_tracer(itrac)%lon)
            if (allocated(sat_dep(region)%io_tracer(itrac)%itau_start))     deallocate(sat_dep(region)%io_tracer(itrac)%itau_start)
            if (allocated(sat_dep(region)%io_tracer(itrac)%itau_end))       deallocate(sat_dep(region)%io_tracer(itrac)%itau_end)
            if (allocated(sat_dep(region)%io_tracer(itrac)%sample_itau))    deallocate(sat_dep(region)%io_tracer(itrac)%sample_itau)
            if (allocated(sat_dep(region)%io_tracer(itrac)%idates))         deallocate(sat_dep(region)%io_tracer(itrac)%idates)
            if (allocated(sat_dep(region)%io_tracer(itrac)%forcing))        deallocate(sat_dep(region)%io_tracer(itrac)%forcing)
            if (allocated(sat_dep(region)%io_tracer(itrac)%nsamples))       deallocate(sat_dep(region)%io_tracer(itrac)%nsamples)
            if (allocated(sat_dep(region)%io_tracer(itrac)%total_wt))       deallocate(sat_dep(region)%io_tracer(itrac)%total_wt)
            if (allocated(sat_dep(region)%io_tracer(itrac)%sampling_strategy))  deallocate(sat_dep(region)%io_tracer(itrac)%sampling_strategy)
        end do ! itrac
    end if

    ! get the name of the correct departurefile
    select case (split_period)
        case ('m')
            write(fname,"(a,a,i4.4,i2.2,a4)") trim(outdir_satellite),'/satellite/sat-track_departures_', midpt_date(1:2),'.nc4'
        case ('d')
            write(fname,"(a,a,i4.4,2i2.2,a4)") trim(outdir_satellite),'/satellite/sat-track_departures_', midpt_date(1:3),'.nc4'
    end select

    ! By default, assume that no region has data for any tracer
    sat_dep(region)%has_data = .false.
    sat_dep(region)%io_tracer(:)%has_data = .false.

    inquire(file=fname, exist=file_exist)
    if (file_exist) then
        nc_id = nc_open(fname, 'r', status)
        IF_NOTOK_RETURN(status=1)

        if (nc_grp_exists(nc_id, region_name(region))) then ! departure file has, e.g., /nam300x200
            group_id = nc_get_group(nc_id, region_name(region), io_status)
            if(io_status /= 0) call escape_tm(' Error opening '//trim(adjustl(fname)))
            sat_dep(region)%has_data = .true.

            do itrac = 1, ntracet
                sat_dep(region)%io_tracer(itrac)%has_data = .false.
                if (.not. nc_grp_exists(group_id, names(itrac))) cycle
                sat_dep(region)%io_tracer(itrac)%name = names(itrac)
                tgrp_id = nc_get_group(group_id, names(itrac), io_status)

                sat_dep(region)%io_tracer(itrac)%nobs = nc_get_dim(tgrp_id, 'n_obs')
                if (sat_dep(region)%io_tracer(itrac)%nobs .ge. 1) then
                    sat_dep(region)%io_tracer(itrac)%has_data = .true.
                    sat_dep(region)%io_tracer(itrac)%idates = nc_read_var(tgrp_id, 'idate')
                    sat_dep(region)%io_tracer(itrac)%lat = nc_read_var(tgrp_id, 'lat')
                    sat_dep(region)%io_tracer(itrac)%lon = nc_read_var(tgrp_id, 'lon')
                    sat_dep(region)%io_tracer(itrac)%forcing = nc_read_var(tgrp_id, 'departures')
                    sat_dep(region)%io_tracer(itrac)%nsamples = nc_read_var(tgrp_id, 'nsamples')
                    sat_dep(region)%io_tracer(itrac)%total_wt = nc_read_var(tgrp_id, 'total_weight')
                    sat_dep(region)%io_tracer(itrac)%sampling_strategy = nc_read_var(tgrp_id, 'sampling_strategy')
                end if
            end do ! itrac
        end if ! nc_grp_exists(nc_id, region_name(region))
        call nc_close(nc_id)

        do itrac = 1, ntracet
            if (.not. sat_dep(region)%io_tracer(itrac)%has_data) cycle
            nobs = sat_dep(region)%io_tracer(itrac)%nobs
            ! convert the sample times to TM5 times
            allocate(sat_dep(region)%io_tracer(itrac)%sample_itau(nobs))
            do i = 1, nobs
                call date2tau(sat_dep(region)%io_tracer(itrac)%idates(:,i), sat_dep(region)%io_tracer(itrac)%sample_itau(i))
            end do

            allocate(sat_dep(region)%io_tracer(itrac)%itau_start(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%itau_end(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%rif(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%rjf(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%ifr(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%jfr(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%ifn(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%jfn(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%wcx(nobs))
            allocate(sat_dep(region)%io_tracer(itrac)%wcy(nobs))

            dxr = dx/xref(region)
            dyr = dy/yref(region)

            do i_obs = 1, nobs
                ! calculate itau_start and itau_end
                select case (sat_dep(region)%io_tracer(itrac)%sampling_strategy(i_obs))
                    case (2) ! instantaneous
                        sat_dep(region)%io_tracer(itrac)%itau_start(i_obs) = sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs)
                        sat_dep(region)%io_tracer(itrac)%itau_end(i_obs)   = sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs)
                    case (3) ! symmetric
                        sat_dep(region)%io_tracer(itrac)%itau_start(i_obs) = sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs) - &
                            mod(sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs) - itaue, ndyn_max/tref(region)) - 1 ! 1 second leeway on either side
                        sat_dep(region)%io_tracer(itrac)%itau_end(i_obs)   = sat_dep(region)%io_tracer(itrac)%itau_start(i_obs) + &
                            ndyn_max/tref(region) + 1 ! 1 second leeway on either side

                        if (output_satellite_verbose) then
                            call tau2date(sat_dep(region)%io_tracer(itrac)%itau_start(i_obs), idate_temp)
                            startTime = NewDate(time6=idate_temp)
                            call tau2date(sat_dep(region)%io_tracer(itrac)%itau_end(i_obs), idate_temp)
                            endTime = NewDate(time6=idate_temp)
                            call tau2date(sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs), idate_temp)
                            sampleTime = NewDate(time6=idate_temp)
                            write(*,'(a, " :: observation over ", f6.2, ",", f7.2, " at ", a, " will be sampled between ", a, " and ", a)') &
                                rname, sat_dep(region)%io_tracer(itrac)%lat(i_obs), sat_dep(region)%io_tracer(itrac)%lon(i_obs), &
                                trim(Pretty(sampleTime)), trim(Pretty(startTime)), trim(Pretty(endTime))
                        end if

                    case default
                        if (output_satellite_verbose) then
                            call tau2date(sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs), idate_temp)
                            sampleTime = NewDate(time6=idate_temp)
                            write(*,'(a, " :: unknown sampling strategy ", i2, " for observation at ", a, " over ", f6.2, ", ", f7.2)') &
                                rname, sat_dep(region)%io_tracer(itrac)%sampling_strategy(i_obs), trim(Pretty(sampleTime)), &
                                sat_dep(region)%io_tracer(itrac)%lat(i_obs), sat_dep(region)%io_tracer(itrac)%lon(i_obs)
                        end if
                end select

                ! now calculate ifr, jfr, etc.
                sat_dep(region)%io_tracer(itrac)%rif(i_obs) = (sat_dep(region)%io_tracer(itrac)%lon(i_obs)-float(xbeg(region)))/dxr + 0.99999
                sat_dep(region)%io_tracer(itrac)%rjf(i_obs) = (sat_dep(region)%io_tracer(itrac)%lat(i_obs)-float(ybeg(region)))/dyr + 0.99999
                sat_dep(region)%io_tracer(itrac)%ifr(i_obs) = int(sat_dep(region)%io_tracer(itrac)%rif(i_obs))   ! i-index of grid cell in which observation is located
                sat_dep(region)%io_tracer(itrac)%jfr(i_obs) = int(sat_dep(region)%io_tracer(itrac)%rjf(i_obs))   ! j-index of grid cell in which observation is located
                sat_dep(region)%io_tracer(itrac)%rif(i_obs) = sat_dep(region)%io_tracer(itrac)%rif(i_obs)-sat_dep(region)%io_tracer(itrac)%ifr(i_obs)-0.5
                sat_dep(region)%io_tracer(itrac)%rjf(i_obs) = sat_dep(region)%io_tracer(itrac)%rjf(i_obs)-sat_dep(region)%io_tracer(itrac)%jfr(i_obs)-0.5

                !the neighbors for pressure interpolation
                if(sat_dep(region)%io_tracer(itrac)%rif(i_obs) .gt. 0) then
                    sat_dep(region)%io_tracer(itrac)%ifn(i_obs) = sat_dep(region)%io_tracer(itrac)%ifr(i_obs) + 1
                else
                    sat_dep(region)%io_tracer(itrac)%ifn(i_obs) = sat_dep(region)%io_tracer(itrac)%ifr(i_obs) - 1
                end if

                if(sat_dep(region)%io_tracer(itrac)%rjf(i_obs) .gt. 0) then
                    sat_dep(region)%io_tracer(itrac)%jfn(i_obs) = sat_dep(region)%io_tracer(itrac)%jfr(i_obs) + 1
                else
                    sat_dep(region)%io_tracer(itrac)%jfn(i_obs) = sat_dep(region)%io_tracer(itrac)%jfr(i_obs) - 1
                end if

                !=================================================================
                ! if index of neighbour is exceeding range of region set
                ! neighbour = current cell (i.e. no interpolation)
                ! in case of cyclic x-boundaries take corresponding cyclic i index
                !=================================================================

                ! y-direction:
                ! ~ original, used in older 'T..' codes and pyshell:
                !if (sat_dep(region)%io_tracer(itrac)%jfn(i_obs) < 1         ) sat_dep(region)%io_tracer(itrac)%jfn(i_obs)=1
                !if (sat_dep(region)%io_tracer(itrac)%jfn(i_obs) > jm(region)) sat_dep(region)%io_tracer(itrac)%jfn(i_obs)=jm(region)
                ! ~ changed to T38 equivalent:
                if (sat_dep(region)%io_tracer(itrac)%jfn(i_obs) < jsr(region)) sat_dep(region)%io_tracer(itrac)%jfn(i_obs) = jsr(region)
                if (sat_dep(region)%io_tracer(itrac)%jfn(i_obs) > jer(region)) sat_dep(region)%io_tracer(itrac)%jfn(i_obs) = jer(region)

                ! x-direction, check on cyclic:
                if ( xcyc(region) == 0 ) then
                    ! non-cyclic boundaries
                    ! ~ original, used in older 'T..' codes and pyshell:
!                    if (sat_dep(region)%io_tracer(itrac)%ifn(i_obs) < 1         ) sat_dep(region)%io_tracer(itrac)%ifn(i_obs)=1
!                    if (sat_dep(region)%io_tracer(itrac)%ifn(i_obs) > im(region)) sat_dep(region)%io_tracer(itrac)%ifn(i_obs)=im(region)
                    ! ~ changed to T38 equivalent:
                    if (sat_dep(region)%io_tracer(itrac)%ifn(i_obs) < isr(region)) sat_dep(region)%io_tracer(itrac)%ifn(i_obs) = isr(region)
                    if (sat_dep(region)%io_tracer(itrac)%ifn(i_obs) > ier(region)) sat_dep(region)%io_tracer(itrac)%ifn(i_obs) = ier(region)
                else
                    ! cyclic x-boundaries
                    if (sat_dep(region)%io_tracer(itrac)%ifn(i_obs) < 1         ) sat_dep(region)%io_tracer(itrac)%ifn(i_obs)=im(region)
                    if (sat_dep(region)%io_tracer(itrac)%ifn(i_obs) > im(region)) sat_dep(region)%io_tracer(itrac)%ifn(i_obs)=1
                end if

                ! x- / y-weighting of grid cell in which observation is located
                sat_dep(region)%io_tracer(itrac)%wcx(i_obs) = (1.0-abs(sat_dep(region)%io_tracer(itrac)%rif(i_obs)))    ! 1.0 ... 0.5
                sat_dep(region)%io_tracer(itrac)%wcy(i_obs) = (1.0-abs(sat_dep(region)%io_tracer(itrac)%rjf(i_obs)))    ! 1.0 ... 0.5
            end do ! i_obs
        end do ! itrac
    else
        if (output_satellite_verbose) write(*,'(a, " :: departure file ",a," not found")') rname, trim(fname)
        period_exists(i_period) = .false.
    end if ! file_exist

    ! since I tried reading this period/region, set sat_data_read to true, even if the file is missing
    ! this way, we won't be trying to read the departure for a period/region repeatedly in case the file does not exist
    sat_data_read(i_period,region) = .true.

    status = 0

end subroutine read_samples

subroutine adj_user_output_satellite_addforcing(region, tr, status)

    !use dims,         only : nregions,im,isr,ier,jm,jsr,jer,lm,itaur,ndyn,tref
    use dims,         only : lm
    use global_data,  only : mass_dat
    use MeteoData,    only : m_dat, humid_dat
    use datetime,     only : tau2date, date2tau
    use toolbox,      only : escape_tm
    use chem_param,   only : ntracet, fscale
    use go,           only : TDate, Get, Pretty, NewDate, rTotal, operator(-)

    implicit none

    !__IO___________________________________________________________________
    character(len=*), parameter       :: rname = mname//'/adj_user_output_satellite_addforcing'

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    ! __CONST_______________________________________________________________

    !__LOCAL_VARIABLES______________________________________________________
    real, dimension(:,:,:,:), pointer :: rm, rxm, rym
    real, dimension(:,:,:), pointer   :: m, q
    logical                           :: in_window
    integer                           :: i_obs, itrac, ifr, jfr, ifn, jfn, idate_temp(6), itau_tr(2)
    real                              :: rif, rjf, wcx, wcy, weight, total_wt
    real, allocatable                 :: rmf(:)
    type(TDate)                       :: t_temp

    !__START_SUBROUTINE______________________________________________________

    status = 0

    ! if there's no data in this region, return
    if (.not. sat_dep(region)%has_data) return

    weight = rTotal(tr(1)-tr(2), 'sec')

    m   => m_dat(region)%data
    rm  => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    q   => humid_dat(region)%data

    !allocate(rmf(lm(region)))

    ! convert tr to integer
    call Get(tr(1), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(1))
    call Get(tr(2), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(2))

    do itrac = 1, ntracet
        if (.not. sat_dep(region)%io_tracer(itrac)%has_data) cycle

        ! Parallelizing the loop below does not work because the arrays rm, rxm, etc. are modified "in place", and different threads
        ! may overwrite the same elements. Consider the following two assignments:
        !
        ! rm(ifr,jfr,:,itrac) = rm(ifr,jfr,:,itrac) + rmf * wcx       * wcy       / m(ifr,jfr,:) * fscale(itrac)
        ! rm(ifn,jfr,:,itrac) = rm(ifn,jfr,:,itrac) + rmf * (1.0-wcx) * wcy       / m(ifn,jfr,:) * fscale(itrac)
        !
        ! For this to be the correct adjoint code, the 'rm' array on the RHS has to be identical for the two statements. However,
        ! in between the two statements, another thread could modify rm. I don't see a way around this. So for now we will not use
        ! OpenMP in the adjoint satellite code.

        !!$omp parallel private(i_obs, in_window, idate_temp, t_temp, ifr, jfr, rif, rjf, wcx, wcy, ifn, jfn, total_wt, rmf) reduction(+:status)
        allocate(rmf(lm(region)))
        !!$omp do schedule(guided) ordered
        do i_obs = 1, sat_dep(region)%io_tracer(itrac)%nobs

            in_window = .false.
            ! tr(1) points to the later time, tr(2) points to the earlier time

            select case(sat_dep(region)%io_tracer(itrac)%sampling_strategy(i_obs))

                case (2) ! instantaneous
                    if (itau_tr(1) .gt. sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs) .and. &
                        itau_tr(2) .le. sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs)) in_window = .true.

                case (3) ! within ndyn/tref
                    if (itau_tr(2) .ge. sat_dep(region)%io_tracer(itrac)%itau_start(i_obs) .and. &
                        itau_tr(1) .le. sat_dep(region)%io_tracer(itrac)%itau_end(i_obs)) in_window = .true.

                case default
                    if (output_satellite_verbose) then
                        call tau2date(sat_dep(region)%io_tracer(itrac)%sample_itau(i_obs), idate_temp)
                        t_temp = NewDate(time6=idate_temp)
                        write(gol,'(a, " :: unknown sampling strategy ", i2, " for observation at ", a, " over ", f6.2, ", ", f7.2, " :: skipping sample")') &
                            rname, sat_dep(region)%io_tracer(itrac)%sampling_strategy(i_obs), trim(Pretty(t_temp)), &
                            sat_dep(region)%io_tracer(itrac)%lat(i_obs), sat_dep(region)%io_tracer(itrac)%lon(i_obs)
                    else
                        write(gol,*) 'Unknown sampling strategy, skipping flask'
                    end if
                    call goErr

            end select

            if (.not. in_window) cycle ! go to the next observation

            ifr = sat_dep(region)%io_tracer(itrac)%ifr(i_obs)
            jfr = sat_dep(region)%io_tracer(itrac)%jfr(i_obs)
            rif = sat_dep(region)%io_tracer(itrac)%rif(i_obs)
            rjf = sat_dep(region)%io_tracer(itrac)%rjf(i_obs)
            wcx = sat_dep(region)%io_tracer(itrac)%wcx(i_obs)
            wcy = sat_dep(region)%io_tracer(itrac)%wcy(i_obs)
            ifn = sat_dep(region)%io_tracer(itrac)%ifn(i_obs)
            jfn = sat_dep(region)%io_tracer(itrac)%jfn(i_obs)

            if (output_satellite_verbose) then
                write(*,'("========== ", a, " ==========")') rname
                write(*,'(a, " :: forcing over ",f6.2,",",f7.2," in region ",i1," at gridbox ",i2,",",i2," between ",a," and ",a)') &
                    rname, sat_dep(region)%io_tracer(itrac)%lat(i_obs), sat_dep(region)%io_tracer(itrac)%lon(i_obs), region, jfr, ifr, &
                    trim(Pretty(tr(1))), trim(Pretty(tr(2)))
            end if

            total_wt = sat_dep(region)%io_tracer(itrac)%total_wt(i_obs)
            rmf = sat_dep(region)%io_tracer(itrac)%forcing(:,i_obs) * weight/total_wt

            select case (satellite_interpolation)

                case (SAT_INTERPOLATION_GRIDBOX)
                    rm(ifr,jfr,:,itrac)  = rm(ifr,jfr,:,itrac)  + rmf / m(ifr,jfr,:) * fscale(itrac)

                case (SAT_INTERPOLATION_SLOPES)
                    ! begin forward code
                    !
                    ! rmf = ( rm(ifr,jfr,:,itr) + 2.0*(rif*rxm(ifr,jfr,:,itr) + rjf*rym(ifr,jfr,:,itr)) ) / m(ifr,jfr,:)*fscale(itr)
                    ! sat_obs(i_obs)%mod_profile = sat_obs(i_obs)%mod_profile + rmf
                    ! ----------------------------------------------------------------------------------------------------------
                    ! sat_obs(i)%mod_profile = sat_obs(i)%mod_profile/sat_obs(i)%nsamples
                    !
                    ! end forward code
                    rm(ifr,jfr,:,itrac)  = rm(ifr,jfr,:,itrac)  + rmf / m(ifr,jfr,:) * fscale(itrac)
                    rxm(ifr,jfr,:,itrac) = rxm(ifr,jfr,:,itrac) + 2.0*rif*rmf / m(ifr,jfr,:) * fscale(itrac)
                    rym(ifr,jfr,:,itrac) = rym(ifr,jfr,:,itrac) + 2.0*rjf*rmf / m(ifr,jfr,:) * fscale(itrac)

                case (SAT_INTERPOLATION_LINEAR)
                    rm(ifr,jfr,:,itrac) = rm(ifr,jfr,:,itrac) + rmf * wcx       * wcy       / m(ifr,jfr,:) * fscale(itrac)
                    rm(ifn,jfr,:,itrac) = rm(ifn,jfr,:,itrac) + rmf * (1.0-wcx) * wcy       / m(ifn,jfr,:) * fscale(itrac)
                    rm(ifr,jfn,:,itrac) = rm(ifr,jfn,:,itrac) + rmf * wcx       * (1.0-wcy) / m(ifr,jfn,:) * fscale(itrac)
                    rm(ifn,jfn,:,itrac) = rm(ifn,jfn,:,itrac) + rmf * (1.0-wcx) * (1.0-wcy) / m(ifn,jfn,:) * fscale(itrac)

                case default
                    write (gol,'("unsupported satellite interpolation index ",i6)') satellite_interpolation; call goErr
                    status = status + 1

            end select

        end do ! i_obs
        !!$omp end do
        deallocate(rmf)
        !!$omp end parallel
        IF_NOTOK_RETURN(status=1)
    end do ! itrac

    !deallocate(rmf)

    nullify(m, rm, rxm, rym)

end subroutine adj_user_output_satellite_addforcing

subroutine adj_user_output_satellite_done(status)

    implicit none

    !__IO___________________________________________________________________
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________
    character(len=*), parameter         :: rname = mname//'/adj_free_satellitedata'

    !__START_SUBROUTINE______________________________________________________
    if (allocated(sat_dep)) deallocate(sat_dep)

    if (allocated(period_exists)) deallocate(period_exists)
    if (allocated(sat_data_read)) deallocate(sat_data_read)

    status = 0

end subroutine adj_user_output_satellite_done

end module adj_user_output_satellite
