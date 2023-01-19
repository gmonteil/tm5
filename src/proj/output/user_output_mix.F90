!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#define CHECK_NCSTAT(nc_ret_code) if (nc_ret_code /= NF90_NOERR) then; status=1; IF_NOTOK_RETURN(status=1); end if
!
!###############################################################################
#include "tm5.inc"

module user_output_mix

use dims,  only : nregions

use datetime_for, only : SEC_PER_DAY

implicit none

private

public :: user_output_mix_init, user_output_mix_step, user_output_mix_done

character(len=512)                  :: outdir_mix, outfile_name_prefix

type mixdata
    real, allocatable               :: mix(:,:,:,:,:)     ! mixing ratio, nx x ny x nz x nt x ntrace
    real, allocatable               :: p(:,:,:)           ! surface pressure, nx x ny x nt
    real, allocatable               :: gph(:,:,:,:)       ! pressure at level boundaries, nx x ny x (nz+1) x nt
    real, allocatable               :: q(:,:,:,:)         ! specific humidity, nx x ny x nz x nt
    real, allocatable               :: t(:,:,:,:)         ! temperature, nx x ny x nz x nt
    integer, allocatable            :: nsamples(:)
    real, allocatable               :: tau(:)             ! should be an integer, but that causes problems while adding stuff
    real, allocatable               :: weight(:)          ! the weight per dynamic timestep
    integer, allocatable            :: times(:,:)         ! integer times, 6 x nt
end type mixdata

type mixfile
    integer                         :: nc_id
    logical                         :: opened  ! flag to indicate whether file is open
    integer, allocatable            :: grp_id(:)
    logical, allocatable            :: exists(:)
end type mixfile

integer                             :: mix_dhour           ! time interval between mixing ratio outputs
integer, allocatable                :: tod_index(:)        ! time-of-day index
logical                             :: mix_output_meteo    ! whether to write out t, q, etc.
integer                             :: deflate_lvl

type(mixdata), dimension(nregions)  :: mixf
type(mixfile)                       :: output_file

character(len=*), parameter         :: mname = 'user_output_mix'

contains

subroutine user_output_mix_init(status)

    use dims,        only: im,jm,lm
    use chem_param,  only: ntrace
    use dims,        only: nregions
    use GO,          only: gol, goErr, goPr
    use GO,          only: ReadRc
    use global_data, only: rcF, outdir

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter :: rname = mname//'/user_output_mix_init'

    integer             :: region
    integer             :: imr, jmr, lmr
    logical             :: dir_exist

    !__START_SUBROUTINE______________________________________________________

    ! read outdir from rcfile
    call ReadRc(rcF, 'output.dir', outdir_mix, status)
    IF_NOTOK_RETURN(status=1)
    outdir_mix = trim(outdir_mix)//'/mix'
    call ReadRc(rcF, 'output.mix.tstep', mix_dhour, status) ! actually, number of seconds per accumulation bin, so mix_dhour is a misnomer
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.mix.filename.prefix', outfile_name_prefix, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.mix.output.meteo', mix_output_meteo, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    call ReadRc(rcF, 'output.mix.deflate.level', deflate_lvl, status, default=1)
    IF_ERROR_RETURN(status=1)

    ! check whether the ouput directory exists or not, and if not, create it
    inquire(file=trim(outdir_mix)//'/.', exist=dir_exist)
    if (.not. dir_exist) call system('mkdir -p '//trim(outdir_mix))

    allocate(tod_index(nregions))

    do region=1,nregions

       imr = im(region)
       jmr = jm(region)
       lmr = lm(region)

       ! allocate array
       allocate(mixf(region)%mix(imr,jmr,lmr,SEC_PER_DAY/mix_dhour,ntrace))
       mixf(region)%mix = 0.0
       allocate(mixf(region)%nsamples(SEC_PER_DAY/mix_dhour))
       mixf(region)%nsamples = 0
       allocate(mixf(region)%weight(SEC_PER_DAY/mix_dhour))
       mixf(region)%weight = 0.0
       allocate(mixf(region)%tau(SEC_PER_DAY/mix_dhour))
       mixf(region)%tau = 0.0
       allocate(mixf(region)%times(6,SEC_PER_DAY/mix_dhour))
       mixf(region)%times = 0
       allocate(mixf(region)%p(imr,jmr,SEC_PER_DAY/mix_dhour))
       mixf(region)%p = 0.0
       allocate(mixf(region)%gph(imr,jmr,lmr+1,SEC_PER_DAY/mix_dhour))
       mixf(region)%gph = 0.0
       if (mix_output_meteo) then
          allocate(mixf(region)%q(imr,jmr,lmr,SEC_PER_DAY/mix_dhour))
          mixf(region)%q = 0.0
          allocate(mixf(region)%t(imr,jmr,lmr,SEC_PER_DAY/mix_dhour))
          mixf(region)%t = 0.0
       end if

    end do

    output_file%nc_id = 0
    output_file%opened = .false.
    allocate(output_file%exists(nregions))
    allocate(output_file%grp_id(nregions))
    do region = 1, nregions
        output_file%exists(region) = .false.
        output_file%grp_id(region) = 0
    end do

    status = 0

end subroutine user_output_mix_init

subroutine user_output_mix_done(status)

    use dims,           only : nregions
    use file_netcdf,    only : nc_close

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    integer   :: region

    !__START_SUBROUTINE______________________________________________________

    do region = 1, nregions
        call output_mixrecord(region)
        call close_mixdatafile(region)
    end do

    if (output_file%opened) then
        call nc_close(output_file%nc_id)
        output_file%opened = .false.
    end if

    do region=1,nregions
       deallocate(mixf(region)%mix)
       deallocate(mixf(region)%nsamples)
       deallocate(mixf(region)%weight)
       deallocate(mixf(region)%times)
       deallocate(mixf(region)%tau)
       deallocate(mixf(region)%p)
       deallocate(mixf(region)%gph)
       if (mix_output_meteo) then
          deallocate(mixf(region)%q)
          deallocate(mixf(region)%t)
       end if
       output_file%exists(region) = .false.
    end do

    deallocate(tod_index)

    status = 0

end subroutine user_output_mix_done

subroutine user_output_mix_step(region, tr, status)

    use dims,           only : itaur,itaui,ndyn,tref
    use dims,           only : newsrun
    use Go,             only : TDate, operator(+), operator(-), operator(/), SecondNumber, Get

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter :: rname = mname//'/user_output_mix_step'
    integer                     :: sec_in_day
    type(TDate)                 :: tmid
    integer, dimension(6)       :: t_temp

    !__START_SUBROUTINE______________________________________________________

    status = 0

    call Get(tr(1), time6=t_temp)
    ! If we are starting a new day, close the old file and open a new file
    if (all(t_temp(4:6) == 0)) then
        if (.not. newsrun) then
            call output_mixrecord(region)
            call close_mixdatafile(region)
        end if
        if (.not. output_file%opened) then
            call open_mixdatafile(t_temp, status)
            IF_NOTOK_RETURN(status=1)
        end if
    end if

    ! which time index to put the mixing ratios in?
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call SecondNumber(tmid, sec_in_day)
    tod_index(region) = sec_in_day/mix_dhour + 1

    call fill_mixrecord(region, tr)

end subroutine user_output_mix_step

subroutine close_mixdatafile(region)

  use file_netcdf,  only : nc_close

  implicit none

  integer, intent(in) :: region

  output_file%exists(region) = .false.
  ! if all the regions have been written, close the file
  if (.not. any(output_file%exists)) then
      call nc_close(output_file%nc_id)
      output_file%opened = .false.
  end if

end subroutine close_mixdatafile

subroutine open_mixdatafile(idate_f, status)

    use dims,           only : nregions, dx, dy, xref, yref, tref, im, jm, lm, xbeg, xend, ybeg, yend, region_name, at, bt
    use chem_param,     only : ntrace, names
    use misctools,      only : check_dir
    use tm5_geometry,   only : lli
    use netcdf,         only : nf90_put_att, NF90_GLOBAL, NF90_UNLIMITED, NF90_NOERR
    use file_netcdf

    implicit none

    !__IO___________________________________________________________________
    integer, intent(in)     :: idate_f(6)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len = 512)   :: FFilename, attr_name
    integer                :: n, region, io_status

    character(len=*), parameter         :: rname = mname//'/open_mixdatafile'

    !__START_SUBROUTINE______________________________________________________

    status = 0

    if (.not. output_file%opened) then
        write(FFilename, '(a, "/", i4.4, "/", i2.2, "/", a, i4.4, 2(i2.2), ".nc4")') &
            trim(outdir_mix), idate_f(1), idate_f(2), trim(outfile_name_prefix), idate_f(1:3)
        call check_dir(FFilename)

        output_file%nc_id = nc_open(FFilename, 'c', status)
        IF_NOTOK_RETURN(status=1)
        output_file%opened = .true.
        ! write the global attributes
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'nregions', nregions) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'dx', dx) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'dy', dy) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'xref', xref(0:nregions)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'yref', yref(0:nregions)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'tref', tref(0:nregions)) ; CHECK_NCSTAT(io_status)

        ! create the global dimensions
        call nc_create_dim(output_file%nc_id, 'times', NF90_UNLIMITED)
        !call nc_create_dim(output_file%nc_id, 'times', SEC_PER_DAY/mix_dhour)
        call nc_create_dim(output_file%nc_id, 'tracers', ntrace)
        call nc_create_dim(output_file%nc_id, 'idate', 6)
        call nc_create_dim(output_file%nc_id, 'levels', lm(1)) ! no zooming in the vertical direction in TM5
        call nc_create_dim(output_file%nc_id, 'boundaries', lm(1)+1)

        ! write the tracer names
        do n = 1, ntrace
            write(attr_name, '("tracer_name_", i3.3)') n
            io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, trim(attr_name), names(n)) ; CHECK_NCSTAT(io_status)
        end do

        ! dump the at and bt coefficients
        call nc_dump_var(output_file%nc_id, 'at', (/ 'boundaries' /), &
            at(:), (/'long_name'/), (/'hybrid sigma-pressure coefficients (P = at + bt * Psurf) in Pascals'/))
        call nc_dump_var(output_file%nc_id, 'bt', (/ 'boundaries' /), &
            bt(:), (/'long_name'/), (/'hybrid sigma-pressure coefficients (P = at + bt * Psurf)'/))

    end if

    do region = 1, nregions
        output_file%grp_id(region) = nc_create_group(output_file%nc_id, region_name(region))
        output_file%exists(region) = .true.
        ! write the region-specific attributes
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'im', im(region)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'jm', jm(region)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'xbeg', xbeg(region)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'ybeg', ybeg(region)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'xend', xend(region)) ; CHECK_NCSTAT(io_status)
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'yend', yend(region)) ; CHECK_NCSTAT(io_status)

        ! create the group dimensions
        call nc_create_dim(output_file%grp_id(region), 'latitude', jm(region))
        call nc_create_dim(output_file%grp_id(region), 'longitude', im(region))
        call nc_create_dim(output_file%grp_id(region), 'lat_edges', jm(region)+1)
        call nc_create_dim(output_file%grp_id(region), 'lon_edges', im(region)+1)

        ! write the latitudes and longitudes
        call nc_dump_var(output_file%grp_id(region), 'latitude', (/'latitude'/), lli(region)%lat_deg, &
            (/'description'/), (/'latitudes of grid cell centers (degrees north)'/))
        call nc_dump_var(output_file%grp_id(region), 'longitude', (/'longitude'/), lli(region)%lon_deg, &
            (/'description'/), (/'longitudes of grid cell centers (degrees east)'/))
        call nc_dump_var(output_file%grp_id(region), 'lat_edges', (/'lat_edges'/), lli(region)%blat_deg, &
            (/'description'/), (/'latitudes of grid cell edges (degrees north)'/))
        call nc_dump_var(output_file%grp_id(region), 'lon_edges', (/'lon_edges'/), lli(region)%blon_deg, &
            (/'description'/), (/'longitudes of grid cell edges (degrees east)'/))

    end do

end subroutine open_mixdatafile

subroutine fill_mixrecord(region, tr)

    use chem_param,   only : ntrace, fscale
    use dims,         only : im, jm, lm
    !use dims,         only : itaur, ndyn, tref, ndyn_max
    use global_data,  only : mass_dat
    use MeteoData   , only : m_dat, gph_dat, phlb_dat, humid_dat, temper_dat
    use datetime,     only : date2tau
    use Go,           only : TDate, rTotal, operator(-), operator(+), operator(/), Get

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)      :: region
    type(TDate), intent(in)  :: tr(2)

    !__LOCAL_VARIABLES______________________________________________________

    integer                   :: n
    integer                   :: imr,jmr,lmr, t_temp(6), local_tau
    real, dimension(:,:,:), pointer      :: m, phlb, gph, q, t
    real, dimension(:,:,:,:), pointer    :: rm
    real                      :: weight
    type(TDate)               :: tmid

    !__START_SUBROUTINE______________________________________________________

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    m       => m_dat(region)%data
    rm      => mass_dat(region)%rm_t
    gph     => gph_dat(region)%data
    phlb    => phlb_dat(region)%data
    q       => humid_dat(region)%data
    t       => temper_dat(region)%data

    weight = rTotal(tr(2) - tr(1), 'sec') ! from emission_fwd.F90

    do n=1,ntrace
       mixf(region)%mix(1:imr,1:jmr,1:lmr,tod_index(region),n) = mixf(region)%mix(1:imr,1:jmr,1:lmr,tod_index(region),n) + &
        weight * rm(1:imr,1:jmr,1:lmr,n) / m(1:imr,1:jmr,1:lmr) * fscale(n)
    end do

    mixf(region)%p(1:imr,1:jmr,tod_index(region)) = mixf(region)%p(1:imr,1:jmr,tod_index(region)) + &
        weight * phlb(1:imr,1:jmr,1)
    mixf(region)%gph(1:imr,1:jmr,1:lmr+1,tod_index(region)) = mixf(region)%gph(1:imr,1:jmr,1:lmr+1,tod_index(region)) + &
        weight * gph(1:imr,1:jmr,1:lmr+1)
    if (mix_output_meteo) then
        mixf(region)%q(1:imr,1:jmr,1:lmr,tod_index(region)) = mixf(region)%q(1:imr,1:jmr,1:lmr,tod_index(region)) + &
            weight * q(1:imr,1:jmr,1:lmr)
        mixf(region)%t(1:imr,1:jmr,1:lmr,tod_index(region)) = mixf(region)%t(1:imr,1:jmr,1:lmr,tod_index(region)) + &
            weight * t(1:imr,1:jmr,1:lmr)
    end if
    ! The 'tau' for this sample should be the midpoint of tr(1) and tr(2)
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call Get(tmid, time6=t_temp)
    call date2tau(t_temp, local_tau)

    mixf(region)%nsamples(tod_index(region)) = mixf(region)%nsamples(tod_index(region)) + 1
    mixf(region)%tau(tod_index(region)) = mixf(region)%tau(tod_index(region)) + weight*local_tau
    mixf(region)%weight(tod_index(region)) = mixf(region)%weight(tod_index(region)) + weight

    nullify(rm, m)   ! reset pointers
    nullify(gph, phlb, q, t)

end subroutine fill_mixrecord

subroutine output_mixrecord(region) ! only called at the end of the day

    use chem_param,   only : ntrace
    use dims,         only : im, jm, lm, at, bt
    use datetime,     only : tau2date
    use file_netcdf

    implicit none

    !__IO___________________________________________________________________

    integer,intent(in)    :: region

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter :: rname = mname//'/output_mixrecord'

    integer                   :: n
    integer                   :: imr,jmr,lmr
    ! integer                   :: io_status

    !__START_SUBROUTINE______________________________________________________

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    ! divide by nsamples to get average mixing ratio over mix_dhour seconds
    do n=1,size(mixf(region)%mix,4)
        mixf(region)%mix(:,:,:,n,:) = mixf(region)%mix(:,:,:,n,:)/mixf(region)%weight(n)
        mixf(region)%p(:,:,n) = mixf(region)%p(:,:,n)/mixf(region)%weight(n)
        mixf(region)%gph(:,:,:,n) = mixf(region)%gph(:,:,:,n)/mixf(region)%weight(n)
        if (mix_output_meteo) then
            mixf(region)%q(:,:,:,n) = mixf(region)%q(:,:,:,n)/mixf(region)%weight(n)
            mixf(region)%t(:,:,:,n) = mixf(region)%t(:,:,:,n)/mixf(region)%weight(n)
        end if
        mixf(region)%tau(n) = mixf(region)%tau(n)/mixf(region)%weight(n)
        call tau2date(int(mixf(region)%tau(n)), mixf(region)%times(:,n))
    end do

    if (deflate_lvl > 0) then
        ! turn on netcdf compression
        nc_variables_deflate = .true.
        nc_deflate_level = deflate_lvl
    end if

    call nc_dump_var(output_file%grp_id(region), 'mix', (/ 'longitude', 'latitude ', 'levels   ', 'times    ', 'tracers  ' /), &
        mixf(region)%mix(:,:,:,:,:), (/ 'long_name' /), (/ 'mixing ratio' /))
    call nc_dump_var(output_file%grp_id(region), 'pressure', (/ 'longitude', 'latitude ', 'times    ' /), &
        mixf(region)%p(:,:,:), (/'long_name'/), (/'surface pressure (in Pa)'/))
    call nc_dump_var(output_file%grp_id(region), 'gph', (/ 'longitude ', 'latitude  ', 'boundaries', 'times     ' /), &
        mixf(region)%gph(:,:,:,:), (/'long_name'/), (/'geopotential height (in meters) at the level boundaries'/))

    if (mix_output_meteo) then
        call nc_dump_var(output_file%grp_id(region), 'q', (/'longitude', 'latitude ', 'levels   ', 'times    '/), &
            mixf(region)%q(:,:,:,:), (/'long_name'/), (/'specific humidty (mass water vapor/mass air)'/))
        call nc_dump_var(output_file%grp_id(region), 't', (/'longitude', 'latitude ', 'levels   ', 'times    '/), &
            mixf(region)%t(:,:,:,:), (/'long_name'/), (/'temperature (K)'/))
    end if
    call nc_dump_var(output_file%grp_id(region), 'nsamples', (/'times'/), mixf(region)%nsamples(:))
    call nc_dump_var(output_file%grp_id(region), 'sample_times', (/'idate','times'/), mixf(region)%times(:,:))

    ! mixing ratios written, so the mixing ratio array can be re-initialized
    mixf(region)%mix = 0.0
    mixf(region)%p = 0.0
    mixf(region)%gph = 0.0
    mixf(region)%q = 0.0
    mixf(region)%t = 0.0
    mixf(region)%nsamples = 0
    mixf(region)%times = 0
    mixf(region)%tau = 0.0
    mixf(region)%weight = 0.0

    if (deflate_lvl > 0) then
        ! turn off netcdf compression that was previously turned on
        nc_variables_deflate = .false.
    end if

end subroutine output_mixrecord

end module user_output_mix
