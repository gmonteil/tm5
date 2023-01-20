!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!###############################################################################
#include "tm5.inc"

module user_output_flux1x1

implicit none

private

public  :: user_output_flux1x1_init, user_output_flux1x1_step, user_output_flux1x1_done, calculate_flux1x1_indices
public  :: flux1x1_dhour, flux1x1_3d, time_index_flux1x1, write_flux1x1, grid_translate_1x1, calculate_grid_translate

character(len=512)                  :: outdir_flux1x1, outfile_name_prefix

type Tfluxdata
    real, allocatable       :: prod(:,:,:,:,:)    ! production, nx x ny x nz x nt x ntrace
    real, allocatable       :: p(:,:,:,:)         ! pressure at level boundaries, nx x ny x (nz+1) x nt
    real, allocatable       :: gph(:,:,:,:)       ! pressure at level boundaries, nx x ny x (nz+1) x nt
    real, allocatable       :: weight(:,:,:)      ! the weight per dynamic timestep, per 1x1 pixel
    integer, allocatable    :: times(:,:)         ! integer times, 6 x nt
    real, allocatable       :: area_m2(:,:)       ! area of each grid cell, in square meters
end type Tfluxdata

type Tfluxfile
    integer                 :: nc_id
    logical                 :: opened           ! flag to indicate whether file is open
end type Tfluxfile

! For each cell in each zoom region, we need to create a list of 1x1 cells it contributes to, and the fraction it
! contributes. For example, for the 3x2 grid cell between (2N, 4N) and (90E, 93E), it contributes to six 1x1 grid
! cells, with i-indices 271, 272 and 273, and j-indices 93 and 94. So we need a list such as [(271,93), (271, 94),
! (272,93), (272,94), (273,93), (273,94)]. For each of these six elements, we need the fractional contribution of
! the original cell to the 1x1 cell. In other words, how much of the cell ((2N, 4N), (90E, 93E)) falls within the
! the 1x1 cell (272,94)? For the particular example, that array of fractions will be [f1, f2, f1, f2, f1, f2],
! where 3(f1+f2) = 1.
type Tlinarray
    integer, allocatable    :: ilist(:)
    integer, allocatable    :: jlist(:)
    real, allocatable       :: frac(:)
    real, allocatable       :: inv_frac(:)
    integer                 :: N
end type Tlinarray

type Tgridmap
    type(Tlinarray), allocatable    :: cell(:,:) ! im(region) x jm(region)
end type Tgridmap

integer                     :: flux1x1_dhour        ! time interval in which to bin the flux
integer                     :: time_index_flux1x1   ! time index to fill during a time step
logical                     :: write_flux1x1        ! variable to tell emission_fwd to fill up flux1x1_3d
logical                     :: flux_is_3d           ! True if flux is 3D, i.e., with 3D production, False if only surface

type(Tfluxdata)             :: flux1x1_3d
type(Tfluxfile)             :: output_file
type(Tgridmap), allocatable :: grid_translate_1x1(:) ! dimension(nregions)

character(len=*), parameter :: mname = 'user_output_flux1x1'

contains

pure function calculate_rect_area(lat1, lat2, lon1, lon2)
    ! Given a rectangular region, calculate the surface area in square meters
    ! Formula :
    ! \Delta S = R_e^2 (lon_2 - lon_1) (sin(lat_2)-sin(lat_1))
    ! where lat_2 > lat_1 and lon_2 > lon_1. It is assumed that the latitudes and
    ! longitudes supplied are in radians.

    use binas,  only : ae

    implicit none

    real, intent(in) :: lat1, lon1, lat2, lon2
    double precision :: calculate_rect_area

    calculate_rect_area = ae * ae * (lon2-lon1) * (sin(lat2) - sin(lat1))

end function calculate_rect_area

subroutine calculate_grid_translate

    use dims,           only : im, jm, nregions, okdebug
    use tm5_geometry,   only : lli, lli_1x1

    implicit none

    integer             :: imr, jmr, region, i1, j1, ncell, i, j
    real                :: lat1, lat2, lon1, lon2, ovlap_area

    ! it is possible for this routine to be called multiple times, because C14 routines want to use this functionality without writing 1x1 fluxes
    if (allocated(grid_translate_1x1)) return

    allocate(grid_translate_1x1(nregions))

    do region = 1, nregions

        imr = im(region)
        jmr = jm(region)

        allocate(grid_translate_1x1(region)%cell(imr,jmr))

        do i = 1, imr
            do j = 1, jmr
                ! The cell (i,j) is between longitudes lli(region)%blon(i-1) and lli(region)%blon(i), and latitudes
                ! lli(region)%blat(j-1) and lli(region)%blat(j), all in radians. First, we count the number of 1x1 cells
                ! which overlap with this box.
                ncell = 0
                do i1 = 1, lli_1x1%nlon
                    ! If this 1x1 cell is completely outside cell (i,j), cycle
                    if (lli_1x1%blon(i1-1) .ge. lli(region)%blon(i) .or. lli_1x1%blon(i1) .le. lli(region)%blon(i-1)) cycle

                    do j1 = 1, lli_1x1%nlat
                        ! If this 1x1 cell is completely outside cell (i,j), cycle
                        if (lli_1x1%blat(j1-1) .ge. lli(region)%blat(j) .or. lli_1x1%blat(j1) .le. lli(region)%blat(j-1)) cycle
                        ncell = ncell + 1
                    end do ! j1
                end do ! i1

                allocate(grid_translate_1x1(region)%cell(i,j)%ilist(ncell))
                allocate(grid_translate_1x1(region)%cell(i,j)%jlist(ncell))
                allocate(grid_translate_1x1(region)%cell(i,j)%frac(ncell))
                allocate(grid_translate_1x1(region)%cell(i,j)%inv_frac(ncell))
                grid_translate_1x1(region)%cell(i,j)%N = ncell

                ! Now calculate the overlaps
                ncell = 0
                do i1 = 1, lli_1x1%nlon
                    if (lli_1x1%blon(i1-1) .ge. lli(region)%blon(i) .or. lli_1x1%blon(i1) .le. lli(region)%blon(i-1)) cycle
                    do j1 = 1, lli_1x1%nlat
                        if (lli_1x1%blat(j1-1) .ge. lli(region)%blat(j) .or. lli_1x1%blat(j1) .le. lli(region)%blat(j-1)) cycle
                        ! We have reached this far, means there is an overlap
                        ncell = ncell + 1

                        lon1 = max(lli_1x1%blon(i1-1), lli(region)%blon(i-1))
                        lon2 = min(lli_1x1%blon(i1), lli(region)%blon(i))
                        lat1 = max(lli_1x1%blat(j1-1), lli(region)%blat(j-1))
                        lat2 = min(lli_1x1%blat(j1), lli(region)%blat(j))
                        ovlap_area = calculate_rect_area(lat1, lat2, lon1, lon2)

                        grid_translate_1x1(region)%cell(i,j)%ilist(ncell) = i1
                        grid_translate_1x1(region)%cell(i,j)%jlist(ncell) = j1
                        ! fraction of grid cell i,j that lies within 1x1 grid cell i1,j1
                        grid_translate_1x1(region)%cell(i,j)%frac(ncell) = ovlap_area/lli(region)%area_m2(j)
                        ! fraction of grid cell i1,j1 that lies within grid cell i,j
                        if ((lli_1x1%blon(i1-1) .ge. lli(region)%blon(i-1)) .and. &
                            (lli_1x1%blon(i1)   .le. lli(region)%blon(i))   .and. &
                            (lli_1x1%blat(j1-1) .ge. lli(region)%blat(j-1)) .and. &
                            (lli_1x1%blat(j1)   .le. lli(region)%blat(j))) then
                            grid_translate_1x1(region)%cell(i,j)%inv_frac(ncell) = 1.0
                        else
                            grid_translate_1x1(region)%cell(i,j)%inv_frac(ncell) = ovlap_area/lli_1x1%area_m2(j1)
                        end if
                    end do ! j1
                end do ! i1

                ! Since every grid cell must be completely covered by a set of 1x1 cells, the fractions must add up
                ! to 1. However, due to rounding error this does not happen, and as a result the sum of flux1x1_3d%prod
                ! and budget_global, over one day, can be different by up to 1 part in 10^4, which is a far worse error
                ! than the 1 part in 10^15 expected from double precision numbers. Therefore, we explicitly make sure
                ! that the sum of fractions for each grid cell is 1.
                grid_translate_1x1(region)%cell(i,j)%frac = grid_translate_1x1(region)%cell(i,j)%frac/sum(grid_translate_1x1(region)%cell(i,j)%frac)

                ! Now print the overlaps for debugging
                if (okdebug) then
                    write(*,'("For region ", i1, ", cell ", 2i3, ", list of i = ")',advance="no") region, i, j
                    do ncell = 1, grid_translate_1x1(region)%cell(i,j)%N
                        write(*, '(i4)', advance="no") grid_translate_1x1(region)%cell(i,j)%ilist(ncell)
                    end do
                    write(*,*)
                    write(*,'("For region ", i1, ", cell ", 2i3, ", list of j = ")',advance="no") region, i, j
                    do ncell = 1, grid_translate_1x1(region)%cell(i,j)%N
                        write(*, '(i4)', advance="no") grid_translate_1x1(region)%cell(i,j)%jlist(ncell)
                    end do
                    write(*,*)
                    write(*,'("For region ", i1, ", cell ", 2i3, ", list of frac = ")',advance="no") region, i, j
                    do ncell = 1, grid_translate_1x1(region)%cell(i,j)%N
                        write(*, '(f5.2)', advance="no") grid_translate_1x1(region)%cell(i,j)%frac(ncell)
                    end do
                    write(*,'(" :: sum = 1 + ", es14.7)') sum(grid_translate_1x1(region)%cell(i,j)%frac) - 1.0
                end if

            end do ! j
        end do ! i

    end do ! region

end subroutine calculate_grid_translate

subroutine user_output_flux1x1_init(status)

    use dims,           only : lm
    use chem_param,     only : ntrace
    use GO,             only : ReadRc
    use global_data,    only : rcF, outdir
    use datetime_for,   only : SEC_PER_DAY
    use tm5_geometry,   only : lli_1x1

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter :: rname = mname//'/user_output_flux1x1_init'

    integer             :: lmr, j
    logical             :: dir_exist

    !__START_SUBROUTINE______________________________________________________

    ! read outdir from rcfile
    call ReadRc(rcF, 'output.dir', outdir_flux1x1, status)
    IF_NOTOK_RETURN(status=1)
    outdir_flux1x1 = trim(outdir_flux1x1)//'/flux'
    call ReadRc(rcF, 'output.flux1x1.tstep', flux1x1_dhour, status) ! actually, number of seconds per accumulation bin, so flux1x1_dhour is a misnomer
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.flux1x1.filename.prefix', outfile_name_prefix, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.flux1x1.3d', flux_is_3d, status, default=.false.)
    IF_ERROR_RETURN(status=1)

    write_flux1x1 = .true.

    ! check whether the ouput directory exists or not, and if not, create it
    inquire(file=trim(outdir_flux1x1)//'/.', exist=dir_exist)
    if (.not. dir_exist) call system('mkdir -p '//trim(outdir_flux1x1))

    lmr = lm(1) ! same vertical grid for all regions

    ! allocate array
    allocate(flux1x1_3d%prod(lli_1x1%nlon,lli_1x1%nlat,lmr,SEC_PER_DAY/flux1x1_dhour,ntrace))
    flux1x1_3d%prod = 0.0
    allocate(flux1x1_3d%times(6,SEC_PER_DAY/flux1x1_dhour))
    flux1x1_3d%times = 0
    allocate(flux1x1_3d%p(lli_1x1%nlon,lli_1x1%nlat,lmr+1,SEC_PER_DAY/flux1x1_dhour))
    flux1x1_3d%p = 0.0
    allocate(flux1x1_3d%gph(lli_1x1%nlon,lli_1x1%nlat,lmr+1,SEC_PER_DAY/flux1x1_dhour))
    flux1x1_3d%gph = 0.0
    allocate(flux1x1_3d%weight(lli_1x1%nlon,lli_1x1%nlat,SEC_PER_DAY/flux1x1_dhour))
    flux1x1_3d%weight = 0.0
    allocate(flux1x1_3d%area_m2(lli_1x1%nlon,lli_1x1%nlat))
    do j = 1, lli_1x1%nlat
        flux1x1_3d%area_m2(:,j) = lli_1x1%area_m2(j)
    end do

    output_file%nc_id = 0
    output_file%opened = .false.

    ! Now calculate grid_translate_1x1
    call calculate_grid_translate

    status = 0

end subroutine user_output_flux1x1_init

subroutine user_output_flux1x1_done(status)

    use file_netcdf

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter :: rname = mname//'/user_output_flux1x1_done'

    !__START_SUBROUTINE______________________________________________________

    ! write the last day's fluxes
    call output_fluxrecord

!    if (output_file%opened) then
!        call nc_close(output_file%nc_id)
!        output_file%opened = .false.
!    end if

    deallocate(flux1x1_3d%prod)
    deallocate(flux1x1_3d%p)
    deallocate(flux1x1_3d%gph)
    deallocate(flux1x1_3d%times)
    deallocate(flux1x1_3d%weight)
    deallocate(flux1x1_3d%area_m2)

    if (allocated(grid_translate_1x1)) deallocate(grid_translate_1x1)

    status = 0

end subroutine user_output_flux1x1_done

subroutine open_fluxdatafile(idate_f, status)

    use chem_param,    only : ntrace
    use datetime,      only : tau2date
    use datetime_for,  only : SEC_PER_DAY
    use dims,          only : lm
    use tm5_geometry,  only : lli_1x1
    use misctools,     only : check_dir
    use file_netcdf
    use netcdf,        only : NF90_UNLIMITED

    implicit none

    integer, intent(in)     :: idate_f(6)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len = 512)    :: FFilename
    integer                 :: n
    integer                 :: io_status
    
    character(len=*), parameter :: rname = mname//'/open_fluxdatafile'

    !__START_SUBROUTINE______________________________________________________

    if (.not. output_file%opened) then
        write(FFilename, '(a, "/", i4.4, "/", i2.2, "/", a, i4.4, 2(i2.2), ".nc4")') &
            trim(outdir_flux1x1), idate_f(1), idate_f(2), trim(outfile_name_prefix), idate_f(1:3)
        call check_dir(FFilename)

        output_file%nc_id = nc_open(FFilename, 'c', status)
        IF_NOTOK_RETURN(status=1)

        ! create the global dimensions
        call nc_create_dim(output_file%nc_id, 'times', NF90_UNLIMITED)
        !call nc_create_dim(output_file%nc_id, 'times', SEC_PER_DAY/flux1x1_dhour)
        call nc_create_dim(output_file%nc_id, 'tracers', ntrace)
        if (flux_is_3d) then
            call nc_create_dim(output_file%nc_id, 'levels', lm(1)) ! no zooming in the vertical direction in TM5
            call nc_create_dim(output_file%nc_id, 'boundaries', lm(1)+1)
        end if
        call nc_create_dim(output_file%nc_id, 'lat1x1', lli_1x1%nlat)
        call nc_create_dim(output_file%nc_id, 'lon1x1', lli_1x1%nlon)

        output_file%opened = .true.
    end if

    status = 0

end subroutine open_fluxdatafile

subroutine output_fluxrecord ! only called at the end of the day

    use chem_param,     only : ntrace, ra, names
    use tm5_geometry,   only : lli_1x1
    use netcdf,         only : nf90_put_att, NF90_GLOBAL
    use file_netcdf

    implicit none

    !__LOCAL_VARIABLES______________________________________________________

    integer                   :: n, nlev, itr, i, j, ntime
    integer                   :: io_status
    real, allocatable         :: sum_emis(:,:), daily_total(:)
    character(len=256)        :: tracer_names

    ! divide by nsamples to get average pressures and GPH's
    nlev = size(flux1x1_3d%prod, 3)
    do n = 1, nlev+1 ! over levels
        flux1x1_3d%p(:,:,n,:)   = flux1x1_3d%p(:,:,n,:)  /flux1x1_3d%weight
        flux1x1_3d%gph(:,:,n,:) = flux1x1_3d%gph(:,:,n,:)/flux1x1_3d%weight
    end do

    ! allocate sum_emis(ntime,ntrace) and make sums of emissions
    ntime = size(flux1x1_3d%prod, 4)
    allocate(sum_emis(ntime, ntrace))
    do itr = 1, ntrace
        do i = 1, ntime
            ! prod contains the total tracer mass added per time step, in Kg
            ! so sum_emis contains the total tracer mass added per time step over the entire atmosphere
            sum_emis(i, itr) = sum(flux1x1_3d%prod(:,:,:,i,itr))
        end do
    end do

    allocate(daily_total(ntrace))
    daily_total = sum(sum_emis, 1)

    ! Up to now prod contains the total tracer mass added per time step, in Kg
    ! Reduce it to micromoles tracer/m^2/sec
    do itr = 1, ntrace
        do i = 1, lli_1x1%nlon
            do j = 1, lli_1x1%nlat
                flux1x1_3d%prod(i,j,:,:,itr) = 1.0D9 * flux1x1_3d%prod(i,j,:,:,itr) / ra(itr) / flux1x1_3d%area_m2(i,j) / flux1x1_dhour
            end do
        end do
    end do

    ! write the tracer names
    do itr = 1, ntrace
        write(tracer_names, '("tracer_", i3.3)') itr
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, trim(tracer_names), trim(names(itr)))
    end do

    ! turn on netcdf compression
    nc_variables_deflate = .true.
    nc_deflate_level = 1

    if (flux_is_3d) then
        call nc_dump_var(output_file%nc_id, 'emission', (/'lon1x1 ', 'lat1x1 ', 'levels ', 'times  ', 'tracers'/), &
            flux1x1_3d%prod, (/'unit'/), (/'micromoles per square meter per second'/))
        call nc_dump_var(output_file%nc_id, 'pressure', (/'lon1x1    ', 'lat1x1    ', 'boundaries', 'times    '/), &
            flux1x1_3d%p, (/'unit'/), (/'layer boundary in Pascals'/))
        call nc_dump_var(output_file%nc_id, 'gph', (/'lon1x1    ', 'lat1x1    ', 'boundaries', 'times    '/), &
            flux1x1_3d%gph, (/'unit'/), (/'layer boundary in meters'/))
    else
        call nc_dump_var(output_file%nc_id, 'emission', (/'lon1x1 ', 'lat1x1 ', 'times  ', 'tracers'/), &
            sum(flux1x1_3d%prod, dim=3), (/'unit'/), (/'micromoles per square meter per second'/))
    end if

    call nc_dump_var(output_file%nc_id, 'area', (/'lon1x1', 'lat1x1'/), flux1x1_3d%area_m2, (/'unit'/), (/'square meters'/))

    ! now write the totals
    call nc_dump_var(output_file%nc_id, 'total_emis_per_timestep', (/'times  ', 'tracers'/), sum_emis, &
        (/'unit'/), (/'Kg tracer per time step'/))
    call nc_dump_var(output_file%nc_id, 'total_emis_per_day', (/'tracers'/), daily_total, &
        (/'unit'/), (/'Kg tracer per day, for this day'/))

    ! write out the molar masses as well
    call nc_dump_var(output_file%nc_id, 'molar_mass', (/'tracers'/), ra, (/'unit'/), (/'grams tracer/mole'/))

    ! turn off netcdf compression
    nc_variables_deflate = .false.

    flux1x1_3d%prod   = 0.0
    flux1x1_3d%p      = 0.0
    flux1x1_3d%gph    = 0.0
    flux1x1_3d%weight = 0.0

    deallocate(daily_total, sum_emis)

    call nc_close(output_file%nc_id)
    output_file%opened = .false.

end subroutine output_fluxrecord

subroutine calculate_flux1x1_indices(tr, time_index, weight)

    ! subroutine to be called from emission_fwd to calculate the time index and weight,
    ! prevents a whole bunch of extra 'use' statements in emission_fwd

    use Go,             only : TDate, operator(+), operator(-), operator(/), SecondNumber, rTotal

    implicit none

    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: time_index
    real, intent(out)       :: weight

    integer                 :: sec_in_day
    type(TDate)             :: tmid

    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call SecondNumber(tmid, sec_in_day)
    time_index = sec_in_day/flux1x1_dhour + 1

    weight = rTotal(tr(2) - tr(1), 'sec') ! from emission_fwd.F90

end subroutine calculate_flux1x1_indices

subroutine user_output_flux1x1_step(region, tr, status)

    use dims,           only : newsrun, nregions, revert
    use Go,             only : TDate, Get
    use go,             only : gol, goErr

    !__IO___________________________________________________________________

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter :: rname = mname//"/user_output_flux1x1_step"
    logical                     :: write_file
    integer, dimension(6)       :: t_temp

    !__START_SUBROUTINE______________________________________________________

    if (revert /= 1) then
        write(gol,'(a, " should only be called from a forward run")') rname ; call goErr
        status = 1 ; IF_NOTOK_RETURN(status=1)
    end if

    ! A new flux file is to be opened if tr(1) points to the beginning of a day, but only once for all regions
    call Get(tr(1), time6=t_temp)
    ! We can't use newday because that stays true during the entire nread seconds at the beginning of a day
    write_file = (region == nregions) .and. all(t_temp(4:6) == 0)
    if (write_file) then
        ! If we are not starting a new run, we need to close the previous day's file
        if (.not. newsrun) call output_fluxrecord
        call open_fluxdatafile(t_temp, status)
        IF_NOTOK_RETURN(status=1)
    end if

    status = 0

end subroutine user_output_flux1x1_step

end module user_output_flux1x1
