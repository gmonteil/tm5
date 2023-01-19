!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module user_output_profile

    use go_date,     only : tdate, tincrdate, itotal, newdate, pretty
    use go_rc,       only : trcfile, readrc, init, done
    use global_data, only : rcf
    use go,          only : newdate

    implicit none

    private

    public :: user_output_profile_init, user_output_profile_step, user_output_profile_done

    type t_profile
        character(len=12)   :: id
        character(len=120)  :: name
        real                :: lat, lon
        integer             :: region
        real, allocatable   :: mix(:, :, :)
        integer             :: ifr, jfr ! i,j region indices for flask's grid cell
        integer             :: ifn, jfn ! i,j region indices for flask's "next" grid cell
        real                :: rif, rjf ! fractions from center of ifr,jfr box
        real                :: surface_height ! surface height in meters
        real                :: wcx, wcy ! x and y weighting factors for slopes interpolation
    end type t_profile

    type(t_profile), dimension(:), allocatable :: profiles

    character(len=*), parameter             :: mname = 'user_output_profiles'
    integer, dimension(:), allocatable      :: region_rank
    integer(4), dimension(:,:), allocatable :: midpoint_dates
    real, dimension(:), allocatable         :: levels
    real    :: zmin_profiles, zmax_profiles
    integer :: nlev_profiles
    integer :: ntstep_out, itstep
    integer :: nprof
    real    :: counter=0
    type(tdate) :: lastOutput

    contains

        ! Interfaces
        subroutine user_output_profile_init(status)
            integer, intent(out)    :: status

            call init_times
            call parse_config_file(status)
            call init_levels
            call init_profiles
            itstep = 1
            status = 0
        end subroutine user_output_profile_init


        subroutine user_output_profile_step(region, tr, status)
            use dims,       only : ndyn_max
            use go_date,    only : operator(+), operator(-), operator(/), rTotal
            use chem_param, only : ntracet

            integer, intent(out)        :: status
            integer, intent(in)         :: region
            type(TDate), intent(in)     :: tr(2)
            type(TDate)                 :: tmid
            integer :: iprof, ilev, nsec, dt
            real, dimension(ntracet,nlev_profiles) :: mix
            logical :: newt
            real    :: frac_t

            ! Determine if the time step should increase or not:
            tmid = tr(1) + (tr(2)-tr(1))/2.0
            nsec = nint(rTotal(tmid - lastOutput, 'sec'))
            newt = nsec >= ndyn_max

            ! If so, add the leftover concentrations and average :
            frac_t = 1.
            if (newt) then
                frac_t = 1-(nsec-ndyn_max)*1./(nint(rTotal(tr(2)-tr(1), 'sec')))
                do iprof = 1, nprof
                    if (profiles(iprof)%region /= region) cycle
                    ! Add the leftover concentration of the time step
                    call interp_mix(profiles(iprof), mix)
                    profiles(iprof)%mix(itstep, :, :) = profiles(iprof)%mix(itstep, :, :) + mix*frac_t
                    ! Average
                    profiles(iprof)%mix(itstep, :, :) = profiles(iprof)%mix(itstep, :, :)/(counter+frac_t)
                enddo
                ! reset the counter and lastOutput, and increase the time step :
                counter = 0.
                itstep = itstep+1
                lastOutput = tmid
                ! reverse frac_t
                frac_t = 1-frac_t
            endif

            ! Accumulate the concentrations
            do iprof = 1, nprof
                if (profiles(iprof)%region /= region) cycle
                call interp_mix(profiles(iprof), mix)
                profiles(iprof)%mix(itstep, :, :) = profiles(iprof)%mix(itstep, :, :) + mix*frac_t
            enddo
            counter = counter + frac_t

            status = 0
        end subroutine user_output_profile_step


        subroutine user_output_profile_done
            integer :: iprof
            ! Do the final averaging:
            do iprof = 1, nprof
                profiles(iprof)%mix(itstep, :, :) = profiles(iprof)%mix(itstep, :, :)/counter
            enddo
            ! Write to file
            call write_profiles
            ! Deallocate
            do iprof = 1, nprof
                deallocate(profiles(iprof)%mix)
            enddo
            deallocate(region_rank)
            deallocate(levels)
            deallocate(midpoint_dates)
            deallocate(profiles)
        end subroutine user_output_profile_done


        ! Internal subroutines
        subroutine parse_config_file(status)
            character(len=*), parameter :: rname = mname//'/parse_config_file'
            character(len=300)  :: profile_filename

            type(trcfile) :: config
            ! Read the rcfile containing the coordinate of the profiles
            integer, intent(out)    :: status
            character(len=11)       :: prefix
            integer                 :: iprof

            call ReadRc(rcF, 'output.profile.filename', profile_filename, status)
            IF_NOTOK_RETURN(status=1)

            call Init(config, trim(profile_filename), status)
            IF_NOTOK_RETURN(status=1)

            call readrc(config, 'profiles.n', nprof, status)
            allocate(profiles(nprof))
            do iprof = 1, nprof
                if (iprof < 10) then
                    write(prefix, '(a,i1,a)') 'profile.', iprof, '.'
                else if (iprof < 100) then
                    write(prefix, '(a,i2,a)') 'profile.', iprof, '.'
                else

                endif
                call readrc(config, trim(prefix)//'id', profiles(iprof)%id, status)
                call readrc(config, trim(prefix)//'lat', profiles(iprof)%lat, status)
                call readrc(config, trim(prefix)//'lon', profiles(iprof)%lon, status)
                call readrc(config, trim(prefix)//'name', profiles(iprof)%name, status)
            enddo
            call readrc(config, 'profiles.zmin', zmin_profiles, status)
            call readrc(config, 'profiles.zmax', zmax_profiles, status)
            call readrc(config, 'profiles.nlev', nlev_profiles, status)
            call done(config, status)
        end subroutine parse_config_file


        subroutine init_profiles
            use chem_param, only : ntracet
            ! Initialize the profiles structures for each site 
            integer :: iprof

            do iprof = 1, nprof
                call assign_region(profiles(iprof))
                call compute_coefficients(profiles(iprof))
                allocate(profiles(iprof)%mix(ntstep_out, ntracet, nlev_profiles))
                profiles(iprof)%mix = 0 
            enddo
        end subroutine init_profiles
        

        subroutine init_times

            use go_date, only : operator(-), IncrDate, operator(+), Get
            use dims,    only : ndyn_max, idatei, idatee

            integer               :: tt
            integer               :: total_seconds
            type(TIncrDate)       :: dt
            type(TDate)           :: dummy_t, start_date, end_date
            integer, dimension(6) :: idate_temp

            start_date    = NewDate(time6=idatei)
            end_date      = NewDate(time6=idatee)
            total_seconds = iTotal(end_date - start_date, 'sec')
            ntstep_out    = total_seconds/ndyn_max
            allocate(midpoint_dates(6,ntstep_out))
            do tt=1,ntstep_out
                dt = IncrDate(sec = tt*ndyn_max - ndyn_max/2)
                dummy_t = start_date + dt
                call Get(dummy_t, time6=idate_temp)
                midpoint_dates(:,tt) = idate_temp(:)
            enddo
            lastOutput = start_date
        end subroutine init_times

        subroutine assign_region(profile)
            use dims,      only : xref, yref, nregions
            use orderpack, only : mrgrnk
            type(t_profile), intent(inout) :: profile
            integer :: i_region, reg

            if (.not. allocated(region_rank)) then
                allocate(region_rank(nregions))
                call mrgrnk(xref(1:nregions) * yref(1:nregions), region_rank)
                region_rank = region_rank(nregions:1:-1) ! finest regions occur first
            endif

            do i_region=1,nregions
                reg = region_rank(i_region)
                if (in_region(profile%lat, profile%lon, reg)) exit
            enddo
            profile%region = reg
        end subroutine assign_region

        subroutine compute_coefficients(profile)
            use dims, only : xref, yref, dx, dy, im, jm, xcyc, xbeg, ybeg
            type(t_profile), intent(inout) :: profile
            profile%rif = (profile%lon-float(xbeg(profile%region)))*xref(profile%region)/dx + 0.99999
            profile%rjf = (profile%lat-float(ybeg(profile%region)))*yref(profile%region)/dy + 0.99999
            profile%ifr = int(profile%rif) ! i-index of grid cell in which observation is located
            profile%jfr = int(profile%rjf) ! j-index of grid cell in which observation is located

            !fraction from the center of the is-box and js-box  (-0.5---+0.5)
            profile%rif = profile%rif-profile%ifr-0.5
            profile%rjf = profile%rjf-profile%jfr-0.5

            !the neighbour for x interpolation
            if (profile%rif > 0) then
                profile%ifn = profile%ifr+1
            else
                profile%ifn = profile%ifr-1
            endif

            !the neighbour for y interpolation
            if (profile%rjf > 0) then
                profile%jfn = profile%jfr+1
            else
                profile%jfn = profile%jfr-1
            endif

            ! x- / y-weighting of grid cell in which observation is located
            profile%wcx = (1.-abs(profile%rif))
            profile%wcy = (1.-abs(profile%rjf))

            !=================================================================
            ! if index of neighbour is exceeding range of region set
            ! neighbour = current cell (i.e. no interpolation)
            ! in case of cyclic x-boundaries take corresponding cyclic i index
            !=================================================================
            profile%jfn = max(profile%jfn, 1)
            profile%jfn = min(profile%jfn, jm(profile%region))
            if ( xcyc(profile%region) == 0 ) then
                ! non-cyclic boundaries
                profile%ifn = max(profile%ifn, 1)
                profile%ifn = min(profile%ifn, im(profile%region))
            else 
                ! cyclic x-boundaries
                if (profile%ifn < 1) profile%ifn = im(profile%region)
                if (profile%ifn > im(profile%region)) profile%ifn = 1
            endif

        end subroutine compute_coefficients

        subroutine init_levels
            real    :: dlev
            integer :: ilev
            allocate(levels(nlev_profiles))
            dlev = (zmax_profiles-zmin_profiles)/float(nlev_profiles-1)
            do ilev = 1, nlev_profiles 
                levels(ilev) = zmin_profiles + (ilev-1)*dlev
            enddo
        end subroutine init_levels

        subroutine interp_mix(profile, mix)
            use dims,        only : lm
            use meteodata,   only : gph_dat, m_dat
            use chem_param,  only : fscale, ntracet
            use global_data, only : mass_dat
            type(t_profile), intent(inout)                       :: profile
            real, dimension(ntracet, nlev_profiles), intent(out) :: mix
            real, dimension(0:lm(profile%region))                :: height
            real, dimension(:,:,:), pointer                      :: gph, m
            real, dimension(:,:,:,:), pointer                    :: rm, rxm, rym, rzm
            real    :: wcx, wcy, wcz
            real    :: rmf
            integer :: ifr, ifn, jfr, jfn, lfr
            integer :: tracer
            integer :: ilevout, ilevin, ilevin_n
            real    :: rlf

            ! pointers to global arrays
            gph  => gph_dat(profile%region)%data
            m    => m_dat(profile%region)%data
            rm   => mass_dat(profile%region)%rm_t
            rxm  => mass_dat(profile%region)%rxm_t
            rym  => mass_dat(profile%region)%rym_t
            rzm  => mass_dat(profile%region)%rzm_t

            wcx = profile%wcx
            wcy = profile%wcy
            ifr = profile%ifr
            ifn = profile%ifn
            jfr = profile%jfr
            jfn = profile%jfn

            ! Determine the altitude profile at the site 
            do ilevin=0,lm(profile%region)
                height(ilevin) =       wcx  *      wcy  * gph(ifr,jfr,ilevin+1) & 
                                + (1.0-wcx) *      wcy  * gph(ifn,jfr,ilevin+1) &
                                +      wcx  * (1.0-wcy) * gph(ifr,jfn,ilevin+1) &
                                + (1.0-wcx) * (1.0-wcy) * gph(ifn,jfn,ilevin+1)
            enddo
            profile%surface_height = height(0)

            do ilevout = 1, nlev_profiles
                if (levels(ilevout) <= profile%surface_height) then
                    ilevin = 1
                    rlf = -.5
                else
                    do ilevin=0,lm(profile%region)
                        if (height(ilevin) > levels(ilevout)) exit
                    enddo
                    rlf = (levels(ilevout)-height(ilevin-1))/(height(ilevin)-height(ilevin-1)) - 0.5
                endif

                !=================================
                !the neighbour for z interpolation
                !=================================
                if (rlf > 0) then
                    ilevin_n = ilevin+1
                else 
                    ilevin_n = ilevin-1
                endif

                ! z-weighting of grid cell in which observation is located
                wcz = (1.0-abs(rlf))  !.0 ... 0.5

                !=========================================================
                ! if vertical neighbor is 0 (which does not exist)
                ! take vertical layer with l=2 for EXTRApolation to ground
                !=========================================================

                if (ilevin_n == 0) then
                    ilevin_n = 2
                    wcz = 1.-rlf
                endif

                !=========================================================
                ! if vertical neighbor is lmr+1 (which does not exist)
                ! -> no interpolation
                !=========================================================
                if (ilevin == lm(profile%region)+1) then
                    ilevin_n = lm(profile%region)
                    wcz = 1.
                endif
                
                do tracer = 1, ntracet
                    rmf = ( &
                         wcx  *      wcy  *      wcz  * rm(ifr,jfr,ilevin ,tracer) / m(ifr,jfr,ilevin )  + &
                    (1.0-wcx) *      wcy  *      wcz  * rm(ifn,jfr,ilevin ,tracer) / m(ifn,jfr,ilevin )  + &
                         wcx  * (1.0-wcy) *      wcz  * rm(ifr,jfn,ilevin ,tracer) / m(ifr,jfn,ilevin )  + &
                    (1.0-wcx) * (1.0-wcy) *      wcz  * rm(ifn,jfn,ilevin ,tracer) / m(ifn,jfn,ilevin )  + &
                         wcx  *      wcy  * (1.0-wcz) * rm(ifr,jfr,ilevin_n,tracer) / m(ifr,jfr,ilevin_n)  + &
                    (1.0-wcx) *      wcy  * (1.0-wcz) * rm(ifn,jfr,ilevin_n,tracer) / m(ifn,jfr,ilevin_n)  + &
                         wcx  * (1.0-wcy) * (1.0-wcz) * rm(ifr,jfn,ilevin_n,tracer) / m(ifr,jfn,ilevin_n)  + &
                    (1.0-wcx) * (1.0-wcy) * (1.0-wcz) * rm(ifn,jfn,ilevin_n,tracer) / m(ifn,jfn,ilevin_n)) * fscale(tracer)
                    !rmf = (rm(ifr, jfr, ilevin, tracer) + 2*( &
                    !          profile%rif*rxm(ifr, jfr, ilevin, tracer) & 
                    !        + profile%rjf*rym(ifr, jfr, ilevin, tracer) &
                    !        + rlf*rzm(ifr, jfr, ilevin, tracer)) &
                    !      ) / m(ifr, jfr, ilevin) * fscale(tracer)
                    mix(tracer, ilevout) = rmf
                enddo
            enddo
            nullify(gph, m, rm , rxm, rym, rzm)

        end subroutine interp_mix

        subroutine write_profiles
            
            use file_netcdf
            use netcdf
            use misctools,  only : check_dir
            use dims,       only : idatei, idatee
            use chem_param, only : ntracet, names
            
            integer            :: iprof, itr
            integer            :: ncf
            integer            :: status
            character(len=300) :: outpath
            type(tdate)        :: form_date
            character(len=8*ntracet)  :: tracer_names

            call readrc(rcf, 'output.dir', outpath, status)
            call check_dir(trim(outpath)//'/profiles')
            do iprof = 1, nprof
                ncf = nc_open(trim(outpath)//'/profiles/profile.'//trim(profiles(iprof)%id)//'.nc', 'c', status)
                call nc_create_dim(ncf, 'tracers', ntracet)
                call nc_create_dim(ncf, 'ncol_d', 6)
                call nc_create_dim(ncf, 'ntstep', ntstep_out)
                call nc_create_dim(ncf, 'nlev', nlev_profiles)
                call nc_dump_var(ncf, 'date_midpoints', (/'ncol_d', 'ntstep'/), midpoint_dates(1:6,:))
                call nc_dump_var(ncf, 'mix', (/'ntstep ', 'tracers', 'nlev   '/), profiles(iprof)%mix)
                call nc_dump_var(ncf, 'levels', (/'nlev'/), levels)
                ! write the starting and ending dates
                form_date = NewDate(time6 = idatei)
                call nc_set_attrs(ncf, 'starting time', trim(Pretty(form_date)))
                form_date = NewDate(time6 = idatee)
                ! write the tracer names
                tracer_names = trim(names(1))
                do itr = 2, ntracet
                    tracer_names = tracer_names//','//trim(names(itr))
                end do
                call nc_set_attrs(ncf, 'tracer_names', trim(tracer_names))
                call nc_set_attrs(ncf, 'ending time', trim(Pretty(form_date)))
                call nc_set_attrs(ncf, 'name', trim(profiles(iprof)%name))
                status = nf90_put_att(ncf, NF90_GLOBAL, 'latitude', profiles(iprof)%lat)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'longitude', profiles(iprof)%lon)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'region', profiles(iprof)%region)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'surface_height', profiles(iprof)%surface_height)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'ifn', profiles(iprof)%ifn)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'ifr', profiles(iprof)%ifr)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'jfn', profiles(iprof)%jfn)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'jfr', profiles(iprof)%jfr)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'rif', profiles(iprof)%rif)
                status = nf90_put_att(ncf, NF90_GLOBAL, 'rjf', profiles(iprof)%rjf)
                call nc_close(ncf)
            enddo

        end subroutine write_profiles

        function in_region(lat,lon,region)

            use dims, only : dx, dy, xref, yref, xbeg, ybeg, xend, yend, parent

            real, intent(in)    :: lat, lon
            integer, intent(in) :: region
            logical             :: in_region

            real    :: dxp, dyp

            in_region = .false.

            ! if region is the global region, then return true
            if (parent(region) == 0) then
                in_region = .true.
                return
            end if

            ! else, go through the exercise
            dxp = dx/xref(parent(region))
            dyp = dy/yref(parent(region))

            if ((lon .gt. xbeg(region)) .and. (lon .lt. xend(region)) .and. (lat .gt. ybeg(region)) .and. (lat .lt. yend(region))) in_region = .true.

        end function in_region

end module user_output_profile
