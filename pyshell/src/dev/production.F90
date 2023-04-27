module production

    use emission_data,  only : T_tracer_info
    use global_data,    only : rcf, region_dat
    use go,             only : readrc, tdate, gol, goerr
    use dims,           only : nregions, im, jm, lm, region_name
    use chem_param,     only : tracer_name_len, ntracet, tracer_name => names

    implicit none

    private

    public :: production_fwd

    type d2d
        real, dimension(:, :), allocatable :: data
    end type d2d

    type t_production_data
        character(len=tracer_name_len)  :: name                 ! tracer name, for convenience
        character(len=200)              :: filename             ! filename with production fields
        integer                         :: profile_type = 1     ! not used yet
        logical                         :: enabled = .false.    ! whether there is production enabled for this tracer
        type(tdate)                     :: data_start, data_end ! start and end of the data
        type(d2d), dimension(nregions)  :: field           ! production field, for each region
        integer                         :: current_index
        real, dimension(lm(1))          :: profile
    end type t_production_data

    type(t_production_data), dimension(:), allocatable, target    :: production_data

    character(len=*), parameter :: mname = 'production'
    character(len=30)           :: rname

contains

    subroutine production_init(status)

        ! This allocates and initialize the "production" data structure (for all tracers),
        ! and initializes it for the tracer(s) for which this is relevant.
        !
        ! The following rc-keys are read:
        ! - production.{tracer}.enabled (default: False)
        ! - production.{tracer}.filename ==> read if the previous one is "True". Assumed to be similar for all regions of a tracer

        integer, intent(out)                :: status
        integer                             :: itrac, ireg
        type(t_production_data), pointer    :: prod

        rname = mname // 'production_init'

        allocate(production_data(ntracet))

        do itrac = 1, ntracet
            prod => production_data(itrac)
            prod%name = tracer_name(itrac)
            call readrc(rcf, trim('production.' // trim(tracer_name(itrac)) // '.enabled'), prod%enabled, status, default=.false.)
            status = handle_err(status)
            if (prod%enabled) then
                call readrc(rcf, trim('production.' // trim(tracer_name(itrac)) // '.filename'), prod%filename, status)
                status = handle_err(status)
            end if
            do ireg = 1, nregions
                allocate(prod%field(ireg)%data(im(ireg), jm(ireg)))
            end do
            call load_profile(prod%profile)
        end do
    end subroutine production_init


    subroutine load_production(prod, start)

        ! Load production data from a netcdf file.
        ! The netCDF file should have:
        ! - one global "time" variable (integer, dimension(nt))
        ! - one group for each region, with a "production" variable in it (real, dimension(nt, nlat, nlon)).
        !
        ! The production should be in kg[tracer]/gridcell/s
        !
        ! The "time" values should be in days since a reference date, given as a (year, month, day) tuple by the "origin" attribute of the "time" variable.
        !
        ! There is currently no check on the end of the time interval, so in theory, a run in 2020 with a production
        ! file stopping in 2010 would work (it would just recycle the last step of that file ...).

        use netcdf
        use go_date,    only : tdate, newdate, rtotal, operator(-), operator(+), operator(>=), operator(<=), incrdate
        use dims,       only : idatee

        type(t_production_data), intent(inout)  :: prod  ! Structure to complete
        type(tdate), intent(in)                 :: start ! Start of the period for which we want the production field

        integer :: fid, gid, varid_prod, varid_time, dimid ! netcdf-related variables
        character(len=NF90_MAX_NAME)            :: dimname ! useless but required variable
        integer, dimension(:), allocatable      :: time    ! time coordinate of the data (days since refdate)
        integer, dimension(3)                   :: origin  ! ref date of the "time" coordinate (/year, month, day/)
        type(tdate) :: origin_date  ! reference date of the "time" coordinate (as a date type)
        integer     :: ntime        ! number of time step in the file
        integer     :: ireg         ! region iterator
        integer     :: status
        integer     :: istart       ! days since the reference time of the data
        integer     :: idx          ! index of the slice of data that we want to read

        rname = mname // '/load_production'

        ! Open file
        status = handle_err(nf90_open(prod%filename, NF90_NOWRITE, fid))

        ! Read the time coordinate -- common to all regions
        status = handle_err(nf90_inq_varid(fid, 'time', varid_time))            ! get id for "time" variable
        status = handle_err(nf90_inq_dimid(fid, 'time', dimid))                 ! get the "time" dimension id
        status = handle_err(nf90_inquire_dimension(fid, dimid, dimname, ntime)) ! get the number of timesteps
        allocate(time(ntime))
        status = handle_err(nf90_get_var(fid, varid_time, time))                ! read the data
        status = handle_err(nf90_get_att(fid, varid_time, 'origin', origin))    ! read the "origin" attribute
        origin_date = newdate(year=origin(1), month=origin(2), day=origin(3))   ! convert it to a date

        ! Convert TM5 date to the same time coordinates (i.e. "days since {origin_date}")
        istart = int(rtotal(start - origin_date, 'day'))

        ! Determine the time index that we want to read. It is the last value for which time_int <= istart is true:
        idx = findloc(time <= istart, .true., back=.true., dim=1)

        ! Read the production for each region
        do ireg = 1, nregions
            status = handle_err(nf90_inq_ncid(fid, region_name(ireg), gid))     ! get the group id
            status = handle_err(nf90_inq_varid(gid, 'production', varid_prod))  ! get the id for the "production" variable
            ! read the data:
            status = handle_err( &
                    nf90_get_var( &
                            gid, &
                            varid_prod, &
                            prod%field(ireg)%data, &
                            start = (/idx, 1, 1/), &
                            count = (/1, im(ireg), jm(ireg)/) &
                            ) &
                    )
        end do

        ! Store the start and end time of the data. This is common to all regions:
        prod%data_start = origin_date + incrdate(day=time(idx))
        if (idx < ntime) then
            prod%data_start = origin_date + incrdate(day=time(idx + 1))
        else
            prod%data_end = newdate(year=idatee(1), month=idatee(2), day=idatee(3), hour=idatee(4), min=idatee(5), sec=idatee(6))
        end if

        status = handle_err(nf90_close(fid))

    end subroutine load_production


    subroutine load_profile(profile)

        ! For now, just enforce stratosphere as layers 16 to 23

        real, dimension(:), intent(out) :: profile
        integer :: status

        rname = mname//'/load_profile'

        profile(:) = 0.
        select case (lm(1))
            case (25)
                profile(16: 23) = 1.
                profile = profile / sum(profile)
            case default
                write(0,'("I do not know where the stratosphere is for ", i2.2, " layers")') lm(1)
                write(0,'("Please code the proper case in ", a)') rname
                status = handle_err(1)
        end select

    end subroutine load_profile


    subroutine production_fwd(itrac, ireg, time_start, time_end)
        use go_date, only : operator(>=), operator(<)
        use dims,   only : isr, jsr, ier, jer
        use go,     only : rtotal, operator(-)
        use global_data, only : mass_dat

        integer, intent(in)                 :: itrac, ireg
        type(tdate), intent(in)             :: time_start, time_end
        type(t_production_data), pointer    :: prod
        integer :: i, j
        real    :: dtime

        prod => production_data(itrac)

        ! Gatekeep tracers without production
        if (.not. prod%enabled) return

        ! Load new data in memory if needed
        if (time_start >= prod%data_end) call load_production(prod, time_start)
        if (time_end < prod%data_start) call load_production(prod, time_start)

        ! timestep emissions
        dtime = abs(rtotal(time_end - time_start, 'sec' ))

        ! Apply the production
        do i = isr(ireg), ier(ireg)
            do j = jsr(ireg), jer(ireg)
                if (region_dat(ireg)%zoomed(i, j) /= ireg) cycle
                mass_dat(ireg)%rm_t(i, j, :, itrac) = mass_dat(ireg)%rm_t(i, j, :, itrac) + prod%profile * dtime * prod%field(ireg)%data(i, j)
            end do
        end do

    end subroutine production_fwd


    function handle_err(status_in) result(status_out)
        use go, only : gol, goerr

        integer, intent(in) :: status_in
        integer             :: status_out

        status_out = status_in
        if (status_in/=0) then
            write (gol,'("in ",a," (",a,", line",i5,")")') rname, __file__, __line__
            call goErr
            return
        end if
    end function handle_err

end module production
