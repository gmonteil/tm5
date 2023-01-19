!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

! !REVISION HISTORY:
!
!   Sourish Basu, April 2014 : Break up input and output into monthly/daily chunks
!
!   Sourish Basu, Oct 2013 : Modifying for multiple tracers
!
!   Sourish Basu, Sep 2011 (I think) : Modified for PyShell 4DVAR

module adj_user_output_flask

    use go,                     only : gol, goErr
    use user_output_flask_data, only : VERT_COORD_ALT, VERT_COORD_PRES
    use user_output_flask_data, only : POINT_INTERPOLATION_GRIDBOX, POINT_INTERPOLATION_SLOPES, POINT_INTERPOLATION_LINEAR
    use user_output_flask_data, only : flask_region_forcing

    implicit none

    private

    public              :: adj_user_output_flask_init
    public              :: adj_user_output_flask_done
    public              :: adj_user_output_flask_step

    real, parameter     :: flask_missing_value=-1.0e34
    integer             :: assim_window

    character(len=*), parameter :: mname = 'adj_user_output_flask'
    character(len=1)            :: split_period
    character(len=1024)         :: indir_point
    logical                     :: use_point_forcing ! if the departure file does not exist, set this to false and do not apply any forcing
    logical                     :: flask_verbose
    logical, allocatable        :: period_exists(:) ! keeps track of which periods have point data
    logical, allocatable        :: point_data_read(:) ! keeps track of which periods we've already tried reading the departures for
    integer                     :: flask_interpolation

    type(flask_region_forcing), allocatable     :: flasks(:)

contains

subroutine adj_user_output_flask_init(status)

    use Go,             only : ReadRc
    use global_data,    only : rcF
    use datetime_for,   only : SEC_PER_DAY
    use datetime,       only : date2tau
    use dims,           only : idatei, idatee, nregions

    implicit none

    ! --- in/out ---------------------------------
    integer, intent(out)        :: status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/adj_user_output_flask_init'

    ! --- local ------------------------------
    integer                     :: n_period, idate_temp(6), tau_beg, tau_end

    call ReadRc( rcF, 'output.point.verbose', flask_verbose, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.point.interpolation', flask_interpolation, status, default=POINT_INTERPOLATION_SLOPES)
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'output.point.timewindow', assim_window, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.dir', indir_point, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.point.split.period', split_period, status, default='a')
    IF_NOTOK_RETURN(status=1)

    ! during adjoint run, idatei > idatee
    select case (split_period)
        case ('a')
            n_period = 1
        case ('m')
            n_period = (idatei(1)-idatee(1))*12 + (idatei(2)-idatee(2)+1)
        case ('d')
            idate_temp = idatei
            idate_temp(4:6) = (/ 0,0,0 /)
            call date2tau(idate_temp, tau_beg)
            idate_temp = idatee
            idate_temp(4:6) = (/ 0,0,0 /)
            call date2tau(idate_temp, tau_end)
            ! tau_beg > tau_end
            n_period = (tau_beg-tau_end)/SEC_PER_DAY + 1
        case default
            write(0,*) 'Wrong split period selected'
            IF_NOTOK_RETURN(status=1)
    end select

    allocate(period_exists(n_period))
    ! assume departures exist for all periods
    period_exists = .true.

    allocate(point_data_read(n_period))
    ! no file has yet been read it
    point_data_read = .false.

    status = 0

end subroutine adj_user_output_flask_init

subroutine read_samples(midpt_date, i_period, status)

    use file_netcdf
    use go_date,        only : TDate, Set, Pretty
    use dims,           only : newsrun, nregions, region_name, dx, dy, xref, yref
    use dims,           only : itaue, ndyn_max, tref, xbeg, ybeg, im, jm, xcyc
    use dims,           only : isr, ier, jsr, jer
    use chem_param,     only : names, ntracet
    use datetime,       only : date2tau, tau2date

    implicit none

    ! --- i/o ----------------------
    integer, intent(in)     :: midpt_date(6), i_period
    integer, intent(out)    :: status

    ! --- const --------------------
    character(len=*), parameter :: rname = mname//'/read_samples'

    ! --- local --------------------
    logical                 :: flask_file_exists
    integer                 :: hid, region, grp_id, itrac, tgrp_id, nflasks, iflask, idate_temp(6)
    real                    :: dxr, dyr
    character(len=512)      :: fname
    type(TDate)             :: windowBeg, windowMid, windowEnd

    use_point_forcing = .false.

    ! Open departure file
    select case (split_period)
        case ('a')
            write(fname,"(a,a)") trim(indir_point), '/point/point_departures.nc4'
        case ('m')
            write(fname,"(a,a,i4.4,i2.2,a)") trim(indir_point), '/point/point_departures_', midpt_date(1:2), '.nc4'
        case ('d')
            write(fname,"(a,a,i4.4,2i2.2,a)") trim(indir_point), '/point/point_departures_', midpt_date(1:3), '.nc4'
    end select

    ! in case this is not a new run, deallocate the forcing array from the previous period
    if (.not. newsrun) then
        if (allocated(flasks)) deallocate(flasks)
    end if

    allocate(flasks(nregions))

    inquire(file=trim(fname), exist=flask_file_exists)

    if (flask_file_exists) then

        hid = nc_open(trim(fname), 'r', status)
        IF_NOTOK_RETURN(status=1)

        do region=1, nregions

            flasks(region)%has_data = nc_grp_exists(hid, region_name(region))

            if (.not. flasks(region)%has_data) cycle

            grp_id = nc_get_group(hid, region_name(region))

            allocate(flasks(region)%io_tracer(ntracet))

            do itrac = 1, ntracet

                flasks(region)%io_tracer(itrac)%has_data = nc_grp_exists(grp_id, names(itrac))

                if (.not. flasks(region)%io_tracer(itrac)%has_data) cycle

                flasks(region)%io_tracer(itrac)%name = names(itrac)
                tgrp_id = nc_get_group(grp_id, names(itrac)) ! tgrp_id is the id for the tracer group within a zoom region

                nflasks = nc_get_dim(tgrp_id, 'samples')

                flasks(region)%io_tracer(itrac)%nsamples            = nc_read_var(tgrp_id, 'nsamples')
                flasks(region)%io_tracer(itrac)%total_wt            = nc_read_var(tgrp_id, 'total_weight')
                flasks(region)%io_tracer(itrac)%sampling_strategy   = nc_read_var(tgrp_id, 'sampling_strategy')
                flasks(region)%io_tracer(itrac)%time_window_length  = nc_read_var(tgrp_id, 'time_window_length')
                flasks(region)%io_tracer(itrac)%forcing             = nc_read_var(tgrp_id, 'forcing')
                flasks(region)%io_tracer(itrac)%times               = nc_read_var(tgrp_id, 'date_components')
                flasks(region)%io_tracer(itrac)%alt                 = nc_read_var(tgrp_id, 'alt')
                flasks(region)%io_tracer(itrac)%lat                 = nc_read_var(tgrp_id, 'lat')
                flasks(region)%io_tracer(itrac)%lon                 = nc_read_var(tgrp_id, 'lon')

                allocate(flasks(region)%io_tracer(itrac)%ifr(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%jfr(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%rif(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%rjf(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%wcx(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%wcy(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%ifn(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%jfn(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%itau_start(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%itau_center(nflasks))
                allocate(flasks(region)%io_tracer(itrac)%itau_end(nflasks))

                if (flask_verbose) write(*,'(a, " :: ",i8," obs in input file for tracer ", a, " in region", a)') rname, nflasks, names(itrac), region_name(region)

                dxr = dx/xref(region)
                dyr = dy/yref(region)

                do iflask = 1,nflasks

                    call date2tau(flasks(region)%io_tracer(itrac)%times(:, iflask), flasks(region)%io_tracer(itrac)%itau_center(iflask))

                    select case (flasks(region)%io_tracer(itrac)%sampling_strategy(iflask))

                        case (1) ! 4-hour average

                            flasks(region)%io_tracer(itrac)%itau_start(iflask) = flasks(region)%io_tracer(itrac)%itau_center(iflask)-assim_window*3600
                            flasks(region)%io_tracer(itrac)%itau_end(iflask)   = flasks(region)%io_tracer(itrac)%itau_center(iflask)+assim_window*3600

                            if (flask_verbose) then
                                call tau2date(flasks(region)%io_tracer(itrac)%itau_center(iflask), idate_temp)
                                call Set(windowMid, time6=idate_temp)
                                call tau2date(flasks(region)%io_tracer(itrac)%itau_start(iflask), idate_temp)
                                call Set(windowBeg, time6=idate_temp)
                                call tau2date(flasks(region)%io_tracer(itrac)%itau_end(iflask), idate_temp)
                                call Set(windowEnd, time6=idate_temp)
                                write(*,'(a, " :: Flask sample at ", a, " will be sampled from ", a, " to ", a)') rname, &
                                    trim(Pretty(windowMid)), trim(Pretty(windowBeg)), trim(Pretty(windowEnd))
                            end if

                        case (2) ! instantaneous

                            flasks(region)%io_tracer(itrac)%itau_start(iflask) = flasks(region)%io_tracer(itrac)%itau_center(iflask)
                            flasks(region)%io_tracer(itrac)%itau_end(iflask)   = flasks(region)%io_tracer(itrac)%itau_center(iflask)

                        case (3) ! dT sampling

                            flasks(region)%io_tracer(itrac)%itau_start(iflask) = flasks(region)%io_tracer(itrac)%itau_center(iflask) &
                                - mod(flasks(region)%io_tracer(itrac)%itau_center(iflask) - itaue, ndyn_max/tref(region)) - 1 ! 1 second leeway on either side
                            flasks(region)%io_tracer(itrac)%itau_end(iflask)   = flasks(region)%io_tracer(itrac)%itau_start(iflask) &
                                + ndyn_max/tref(region) + 1 ! 1 second leeway on either side

                            if (flask_verbose) then
                                call tau2date(flasks(region)%io_tracer(itrac)%itau_center(iflask), idate_temp)
                                call Set(windowMid, time6=idate_temp)
                                call tau2date(flasks(region)%io_tracer(itrac)%itau_start(iflask), idate_temp)
                                call Set(windowBeg, time6=idate_temp)
                                call tau2date(flasks(region)%io_tracer(itrac)%itau_end(iflask), idate_temp)
                                call Set(windowEnd, time6=idate_temp)
                                write(*,'(a, " :: Flask sample in region ", a, " at ", a, " will be sampled from ", a, " to ", a)') &
                                    rname, trim(region_name(region)), trim(Pretty(windowMid)), trim(Pretty(windowBeg)), &
                                    trim(Pretty(windowEnd))
                            end if

                        case (4) ! custom time window length per sample

                            flasks(region)%io_tracer(itrac)%itau_start(iflask) = flasks(region)%io_tracer(itrac)%itau_center(iflask) &
                                - flasks(region)%io_tracer(itrac)%time_window_length(iflask)
                            flasks(region)%io_tracer(itrac)%itau_end(iflask)   = flasks(region)%io_tracer(itrac)%itau_center(iflask) &
                                + flasks(region)%io_tracer(itrac)%time_window_length(iflask)

                    end select

                    ! calculate ifr, jfr, rif, rjf, wcx, wcy
                    flasks(region)%io_tracer(itrac)%rif(iflask) = (flasks(region)%io_tracer(itrac)%lon(iflask)-float(xbeg(region)))/dxr + 0.99999
                    flasks(region)%io_tracer(itrac)%rjf(iflask) = (flasks(region)%io_tracer(itrac)%lat(iflask)-float(ybeg(region)))/dyr + 0.99999
                    flasks(region)%io_tracer(itrac)%ifr(iflask) = int(flasks(region)%io_tracer(itrac)%rif(iflask))   ! i-index of grid cell in which observation is located
                    flasks(region)%io_tracer(itrac)%jfr(iflask) = int(flasks(region)%io_tracer(itrac)%rjf(iflask))   ! j-index of grid cell in which observation is located
                    flasks(region)%io_tracer(itrac)%rif(iflask) = flasks(region)%io_tracer(itrac)%rif(iflask)-flasks(region)%io_tracer(itrac)%ifr(iflask)-0.5
                    flasks(region)%io_tracer(itrac)%rjf(iflask) = flasks(region)%io_tracer(itrac)%rjf(iflask)-flasks(region)%io_tracer(itrac)%jfr(iflask)-0.5

                    ! calculate ifn, jfn
                    if (flasks(region)%io_tracer(itrac)%rif(iflask) .gt. 0) then
                        flasks(region)%io_tracer(itrac)%ifn(iflask) = flasks(region)%io_tracer(itrac)%ifr(iflask)+1
                    else
                        flasks(region)%io_tracer(itrac)%ifn(iflask) = flasks(region)%io_tracer(itrac)%ifr(iflask)-1
                    end if
                    if (flasks(region)%io_tracer(itrac)%rjf(iflask) .gt. 0) then
                        flasks(region)%io_tracer(itrac)%jfn(iflask) = flasks(region)%io_tracer(itrac)%jfr(iflask)+1
                    else
                        flasks(region)%io_tracer(itrac)%jfn(iflask) = flasks(region)%io_tracer(itrac)%jfr(iflask)-1
                    end if

                    ! x- / y-weighting of grid cell in which observation is located
                    flasks(region)%io_tracer(itrac)%wcx(iflask) = (1.0-abs(flasks(region)%io_tracer(itrac)%rif(iflask)))    ! 1.0 ... 0.5
                    flasks(region)%io_tracer(itrac)%wcy(iflask) = (1.0-abs(flasks(region)%io_tracer(itrac)%rjf(iflask)))    ! 1.0 ... 0.5

                    !if (flasks(region)%io_tracer(itrac)%jfn(iflask) < 1) flasks(region)%io_tracer(itrac)%jfn(iflask)=1
                    !if (flasks(region)%io_tracer(itrac)%jfn(iflask) > jm(region)) flasks(region)%io_tracer(itrac)%jfn(iflask)=jm(region)
                    if (flasks(region)%io_tracer(itrac)%jfn(iflask) < jsr(region)) flasks(region)%io_tracer(itrac)%jfn(iflask) = jsr(region)
                    if (flasks(region)%io_tracer(itrac)%jfn(iflask) > jer(region)) flasks(region)%io_tracer(itrac)%jfn(iflask) = jer(region)
                    if ( xcyc(region) == 0 ) then
                        ! non-cyclic boundaries
                        !if (flasks(region)%io_tracer(itrac)%ifn(iflask) < 1) flasks(region)%io_tracer(itrac)%ifn(iflask)=1
                        !if (flasks(region)%io_tracer(itrac)%ifn(iflask) > im(region)) flasks(region)%io_tracer(itrac)%ifn(iflask)=im(region)
                        if (flasks(region)%io_tracer(itrac)%ifn(iflask) < isr(region)) flasks(region)%io_tracer(itrac)%ifn(iflask) = isr(region)
                        if (flasks(region)%io_tracer(itrac)%ifn(iflask) > ier(region)) flasks(region)%io_tracer(itrac)%ifn(iflask) = ier(region)
                    else
                        ! cyclic x-boundaries
                        if (flasks(region)%io_tracer(itrac)%ifn(iflask) < 1) flasks(region)%io_tracer(itrac)%ifn(iflask)=im(region)
                        if (flasks(region)%io_tracer(itrac)%ifn(iflask) > im(region)) flasks(region)%io_tracer(itrac)%ifn(iflask)=1
                    end if

                end do ! iflask = 1,nflasks

            end do ! itrac = 1, ntracet

            ! If no tracer has data for this region, set the regional has_data to false
            flasks(region)%has_data = any(flasks(region)%io_tracer(:)%has_data)

        end do ! region=1, nregions

        write(*,'(a," : read in departures from ",a)') rname, trim(fname)

        call nc_close(hid)
        period_exists(i_period) = .true.
        use_point_forcing = any(flasks(:)%has_data)

    else

        use_point_forcing = .false.
        period_exists(i_period) = .false.

    end if

    point_data_read(i_period) = .true. ! I have already tried reading this period, whether or not there are any data

    status=0

end subroutine read_samples

subroutine adj_user_output_flask_addforcing(region, tr, status)

    use global_data,        only : mass_dat
    use MeteoData,          only : m_dat, gph_dat
    use chem_param,         only : ntracet, fscale
    use dims,               only : lm
    use go,                 only : TDate, Get, Pretty, rTotal, operator(-)
    use datetime,           only : date2tau

    implicit none

    ! __IO__________________________________________________________________
    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    ! __CONST_______________________________________________________________
    character(len=*), parameter   ::  rname = mname//'/adj_user_output_flask_addforcing'

    ! __LOCAL_______________________________________________________________
    real, dimension(:,:,:), pointer     :: m, gph
    real, dimension(:,:,:,:), pointer   :: rm, rxm, rym, rzm
    real, dimension(0:lm(region))       :: height

    integer                             :: itrac, iflask, ifr, jfr, lfr, l, ifn, jfn, lfn, idate_temp(6), itau_tr(2)
    logical                             :: in_window
    real                                :: rmf, rif, rjf, rlf, wcx, wcy, wcz, alt, weight, total_wt

    ! __BEGIN_______________________________________________________________

    status = 0

    ! If we don't need to apply forcing, return
    if (.not. use_point_forcing) return

    ! unless there's flask data in this region, leave
    if (.not. flasks(region)%has_data) return

    weight = rTotal(tr(1)-tr(2), 'sec')

    m   => m_dat(region)%data
    rm  => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    gph => gph_dat(region)%data

    ! convert tr to integer
    call Get(tr(1), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(1))
    call Get(tr(2), time6=idate_temp)
    call date2tau(idate_temp, itau_tr(2))

    do itrac = 1, ntracet
        if (.not. flasks(region)%io_tracer(itrac)%has_data) cycle

        do iflask = 1, size(flasks(region)%io_tracer(itrac)%nsamples)

            in_window = .false.
            ! itaur(region) points to the end of the dynamic time step, i.e., the later time
            ! tr(1) points to the later time, tr(2) points to the earlier time

            select case (flasks(region)%io_tracer(itrac)%sampling_strategy(iflask))

               case (1,3,4) ! 4-hour average, or dT sampling, or sampling with custom time window

                    ! the dynamic timestep tr(2) --> tr(1) must be completely inside the interval (itau_start, itau_end)
                    if (itau_tr(2) .ge. flasks(region)%io_tracer(itrac)%itau_start(iflask) .and. &
                        itau_tr(1) .le. flasks(region)%io_tracer(itrac)%itau_end(iflask)) in_window = .true.

               case (2) ! instantaneous

                    if (itau_tr(1) .gt. flasks(region)%io_tracer(itrac)%itau_center(iflask) .and. &
                        itau_tr(2) .le. flasks(region)%io_tracer(itrac)%itau_center(iflask)) in_window = .true.

               case default

                  write(*,*) 'Unknown sampling strategy, skipping flask'
              end select

            if (in_window) then

                ifr = flasks(region)%io_tracer(itrac)%ifr(iflask)
                jfr = flasks(region)%io_tracer(itrac)%jfr(iflask)
                rif = flasks(region)%io_tracer(itrac)%rif(iflask)
                rjf = flasks(region)%io_tracer(itrac)%rjf(iflask)
                wcx = flasks(region)%io_tracer(itrac)%wcx(iflask)
                wcy = flasks(region)%io_tracer(itrac)%wcy(iflask)
                ifn = flasks(region)%io_tracer(itrac)%ifn(iflask)
                jfn = flasks(region)%io_tracer(itrac)%jfn(iflask)
                alt = flasks(region)%io_tracer(itrac)%alt(iflask)

                ! The XY grid is static, but the vertical layer for sampling can vary from timestep to timestep
                ! So we need to calculate lfr as in the forward code
                lfr = 1 !layer
                do l=0,lm(region)
                    height(l) = wcx*wcy*gph(ifr,jfr,l+1) + (1.0-wcx)*wcy*gph(ifn,jfr,l+1) + wcx*(1.0-wcy)*gph(ifr,jfn,l+1) + (1.0-wcx)*(1.0-wcy)*gph(ifn,jfn,l+1)
                enddo

                do l=0,lm(region) ! selects layer , note that we start from second layer from surface
                    if(height(l) .gt. alt) exit
                end do

                select case(l)
                    case(0)
                        lfr = 1
                        rlf = -0.5  !surface...
                    case default
                    lfr = l  !the site layer
                    ! the offset from the center of the layer (-0.5--->+0.5)
                    ! (interpolation is in (m))
                    rlf = (alt-height(l-1))/(height(l)-height(l-1)) - 0.5
                end select

                ! the neighbour for z interpolation
                if ( rlf .gt. 0 ) then
                   lfn = lfr+1
                else
                   lfn = lfr-1
                endif
                ! z-weighting of grid cell in which observation is located
                wcz = (1.0-abs(rlf))  !.0 ... 0.5

                ! edges ...
                if ( lfn == 0 ) then
                   ! if vertical neighbor is 0 (which does not exist)
                   ! take vertical layer with l=2 for EXTRApolation to ground
                   lfn=2
                   wcz=1.0-rlf  ! 1.0 ... 1.5
                else if( lfn == lm(region)+1 ) then
                   ! if vertical neighbor is lmr+1 (which does not exist)
                   ! -> no interpolation
                   lfn=lm(region) ! no interpolation
                   wcz=1.0
                end if

                total_wt = flasks(region)%io_tracer(itrac)%total_wt(iflask)
                rmf = flasks(region)%io_tracer(itrac)%forcing(iflask) * weight/total_wt

                select case (flask_interpolation)
                    case (POINT_INTERPOLATION_GRIDBOX)
                        rm(ifr,jfr,lfr,itrac)  = rm(ifr,jfr,lfr,itrac)  + rmf / m(ifr,jfr,lfr) * fscale(itrac)

                    case (POINT_INTERPOLATION_SLOPES)
                        ! begin forward code
                        ! rmf = ( rm(ifr,jfr,lfr,itr) + &
                        !       2.0 * (rif*rxm(ifr,jfr,lfr,itr) + rjf*rym(ifr,jfr,lfr,itr) + rlf*rzm(ifr,jfr,lfr,itr)) ) / &
                        !       m(ifr,jfr,lfr)*fscale(itr)
                        ! flasks(iflask)%mix(itr) = flasks(iflask)%mix(itr) + rmf
                        ! flasks(iflask)%nsamples = flasks(iflask)%nsamples + 1
                        ! ----------------------------------------------------------
                        ! flasks(iflask)%mix(itr)=flasks(iflask)%mix(itr)/flasks(iflask)%nsamples
                        ! end forward code
                        rm(ifr,jfr,lfr,itrac)  = rm(ifr,jfr,lfr,itrac)  + rmf / m(ifr,jfr,lfr) * fscale(itrac)
                        rxm(ifr,jfr,lfr,itrac) = rxm(ifr,jfr,lfr,itrac) + 2.0*rif*rmf / m(ifr,jfr,lfr) * fscale(itrac)
                        rym(ifr,jfr,lfr,itrac) = rym(ifr,jfr,lfr,itrac) + 2.0*rjf*rmf / m(ifr,jfr,lfr) * fscale(itrac)
                        rzm(ifr,jfr,lfr,itrac) = rzm(ifr,jfr,lfr,itrac) + 2.0*rlf*rmf / m(ifr,jfr,lfr) * fscale(itrac)

                    case (POINT_INTERPOLATION_LINEAR)
                        rm(ifr,jfr,lfr,itrac) = rm(ifr,jfr,lfr,itrac) + rmf * wcx       * wcy       * wcz       / m(ifr,jfr,lfr) * fscale(itrac)
                        rm(ifn,jfr,lfr,itrac) = rm(ifn,jfr,lfr,itrac) + rmf * (1.0-wcx) * wcy       * wcz       / m(ifn,jfr,lfr) * fscale(itrac)
                        rm(ifr,jfn,lfr,itrac) = rm(ifr,jfn,lfr,itrac) + rmf * wcx       * (1.0-wcy) * wcz       / m(ifr,jfn,lfr) * fscale(itrac)
                        rm(ifn,jfn,lfr,itrac) = rm(ifn,jfn,lfr,itrac) + rmf * (1.0-wcx) * (1.0-wcy) * wcz       / m(ifn,jfn,lfr) * fscale(itrac)
                        rm(ifr,jfr,lfn,itrac) = rm(ifr,jfr,lfn,itrac) + rmf * wcx       * wcy       * (1.0-wcz) / m(ifr,jfr,lfn) * fscale(itrac)
                        rm(ifn,jfr,lfn,itrac) = rm(ifn,jfr,lfn,itrac) + rmf * (1.0-wcx) * wcy       * (1.0-wcz) / m(ifn,jfr,lfn) * fscale(itrac)
                        rm(ifr,jfn,lfn,itrac) = rm(ifr,jfn,lfn,itrac) + rmf * wcx       * (1.0-wcy) * (1.0-wcz) / m(ifr,jfn,lfn) * fscale(itrac)
                        rm(ifn,jfn,lfn,itrac) = rm(ifn,jfn,lfn,itrac) + rmf * (1.0-wcx) * (1.0-wcy) * (1.0-wcz) / m(ifn,jfn,lfn) * fscale(itrac)

                    case default
                        write (gol,'("unsupported point interpolation index ",i6)') flask_interpolation; call goErr
                        IF_NOTOK_RETURN(status=1)

                end select

                if (flask_verbose) then
                    write(*, '(a, " :: Forced tracer ", i2, " by ", f10.4, " in cell ", i3, ",", i3, ",", i3, " between ", a, " and ", a)') &
                        rname, itrac, rmf, ifr, jfr, lfr, trim(Pretty(tr(1))), trim(Pretty(tr(2)))
                    write(*, '(a, " :: rif = ", f9.6, ", rjf = ", f9.6, ", rlf = ", f9.6, ", rxm = ", f15.5, ", rym = ", f15.5, ", rzm = ", f15.5)') &
                        rname, rif, rjf, rlf, rxm(ifr,jfr,lfr,itrac), rym(ifr,jfr,lfr,itrac), rzm(ifr,jfr,lfr,itrac)
                    write(*, '(a, " :: Adjoint tracer mass added = ", es20.12, ", air mass = ", es20.12, ", fscale = ", es18.10)') &
                        rname, rmf / m(ifr,jfr,lfr) * fscale(itrac), m(ifr,jfr,lfr), fscale(itrac)
                end if

            end if ! in_window

        end do ! iflask
    end do ! itrac

    nullify(m, gph)
    nullify(rm, rxm, rym, rzm)

end subroutine adj_user_output_flask_addforcing

subroutine adj_user_output_flask_step(region, tr, status)

    use datetime,       only : tau2date, date2tau
    use datetime_for,   only : SEC_PER_DAY
    use dims,           only : idatee
    use Go,             only : TDate, Get, operator(+), operator(-), operator(/)

    implicit none

    !__IO___________________________________________________________________
    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________
    integer                     :: i_period, midpt_date(6), tau_beg, idate_temp(6), tau_mid
    character(len=*), parameter :: rname = mname//'/adj_user_output_flask_step'
    type(TDate)                 :: tmid

    !__START_SUBROUTINE______________________________________________________

    ! during the adjoint run, tr(1) > tr(2)
    ! if the run period is Jan 1 to Jan 3, then idatei => Jan 3, idatee => Jan 1
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call Get(tmid, time6=midpt_date)

    select case (split_period)
        case ('a')
            i_period = 1
        case ('m')
            i_period = (midpt_date(1)-idatee(1))*12 + (midpt_date(2)-idatee(2)+1)
        case ('d')
            idate_temp = idatee
            idate_temp(4:6) = 0
            call date2tau(midpt_date, tau_mid)
            call date2tau(idate_temp, tau_beg)
            i_period = (tau_mid - tau_beg)/SEC_PER_DAY + 1
    end select

    if (.not. point_data_read(i_period)) then
        call read_samples(midpt_date, i_period, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (period_exists(i_period)) call adj_user_output_flask_addforcing(region, tr, status)
    IF_NOTOK_RETURN(status=1)

    status = 0

end subroutine adj_user_output_flask_step

subroutine adj_user_output_flask_done(status)

    implicit none
    integer, intent(out)    :: status

    ! local variables
    character(len=*), parameter :: rname = mname//'/adj_user_output_flask_done'

    if (allocated(flasks)) deallocate(flasks)
    if (allocated(period_exists)) deallocate(period_exists)
    if (allocated(point_data_read)) deallocate(point_data_read)

    status = 0

end subroutine adj_user_output_flask_done

end module adj_user_output_flask
