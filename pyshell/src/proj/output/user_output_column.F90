!### macro's ###################################################################
!
#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module user_output_column

  ! module included to write instantaneous output fields
  ! called from user_output

  use dims,  only : nregions

  use toolbox,      only : escape_tm
  use datetime_for, only : SEC_PER_DAY

  implicit none

  ! interface

  private

  public :: user_output_column_step
  public :: user_output_column_init
  public :: user_output_column_done
!  public :: write_columnfiles
  public :: column_dhour

  logical                             :: write_columnfiles
  character(len=512)                  :: outdir_column, outfile_name_prefix

  type columndata
     real, allocatable    :: col(:,:,:,:)   ! mixing ratio, nlon x nlat x ntime x ntrace
     real, allocatable    :: pres(:,:,:)    ! surface pressure, nlon x nlat x ntime
     integer, allocatable :: nsamples(:)    ! the number of times a particular timestep is sampled
     real, allocatable    :: tau(:)         ! should be an integer, but that causes problems while adding stuff
     real, allocatable    :: weight(:)      ! the weight per dynamic timestep
     integer, allocatable :: times(:,:)     ! integer times, 6 x nt
  end type columndata

  type columnfile
      integer :: nc_id
      logical :: opened
      integer, allocatable :: grp_id(:) ! the group ID for each region
      logical, allocatable :: exists(:) ! whether column output file has a group for a certain region
  end type columnfile

  integer              :: column_dhour           ! time interval between mixing ratio outputs
  integer, allocatable :: tod_index(:)           ! time-of-day index

  type(columndata), dimension(nregions)  :: columnf
  type(columnfile)  :: output_file

  character(len=*), parameter         :: mname = 'user_output_column'

contains

  subroutine user_output_column_init(status)

    use dims,        only: im,jm,ndyn
    use chem_param,  only: ntrace
    use dims,        only: nregions
    use GO,          only: ReadRc
    use global_data, only: rcF

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    integer             :: region, dynamic_timestep
    integer             :: imr,jmr
    character(len=*), parameter      :: rname = mname//', user_output_column_init'
    logical             :: dir_exist


    !__START_SUBROUTINE______________________________________________________

    status = -1

    ! read outdir from rcfile
    call ReadRc(rcF, 'output.dir', outdir_column, status)
    IF_NOTOK_RETURN(status=1)
    outdir_column = trim(outdir_column)//'/col'
    call ReadRc(rcF, 'output.totalcol.tstep', column_dhour, status)
    IF_NOTOK_RETURN(status=1)
    call ReadRc(rcF, 'output.totalcol.filename.prefix', outfile_name_prefix, status)
    IF_NOTOK_RETURN(status=1)

    ! check whether the ouput directory exists or not, and if not, create it
    inquire(file=trim(outdir_column)//'/.', exist=dir_exist)
    if (.not. dir_exist) call system('mkdir -p '//trim(outdir_column))

    allocate(tod_index(nregions))

    do region=1,nregions

       imr = im(region)
       jmr = jm(region)

       ! allocate array
       allocate(columnf(region)%col(imr,jmr,SEC_PER_DAY/column_dhour,ntrace))
       allocate(columnf(region)%pres(imr,jmr,SEC_PER_DAY/column_dhour))
       allocate(columnf(region)%nsamples(SEC_PER_DAY/column_dhour))
       columnf(region)%nsamples = 0
       columnf(region)%col = 0.0
       columnf(region)%pres = 0.0
       allocate(columnf(region)%tau(SEC_PER_DAY/column_dhour))
       columnf(region)%tau = 0.0
       allocate(columnf(region)%times(6,SEC_PER_DAY/column_dhour))
       columnf(region)%times = 0
       allocate(columnf(region)%weight(SEC_PER_DAY/column_dhour))
       columnf(region)%weight = 0.0

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

  end subroutine user_output_column_init

  subroutine user_output_column_done(status)

    use dims,           only : nregions
    use file_netcdf,    only : nc_close

    implicit none

    !__IO___________________________________________________________________

    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    integer   :: region
    character(len=*), parameter :: rname = mname//'/user_output_column_done'

    !__START_SUBROUTINE______________________________________________________

    do region = 1, nregions
        call output_columnrecord(region)
        call close_columndatafile(region)
    end do

    if (output_file%opened) then
        call nc_close(output_file%nc_id)
        output_file%opened = .false.
    end if

    do region = 1, nregions
       deallocate(columnf(region)%col)
       deallocate(columnf(region)%pres)
       deallocate(columnf(region)%nsamples)
       deallocate(columnf(region)%tau)
       deallocate(columnf(region)%weight)
       deallocate(columnf(region)%times)
       output_file%exists(region) = .false.
    end do

    deallocate(tod_index)

    status = 0

  end subroutine user_output_column_done

  subroutine user_output_column_step(region, tr, status)

    use dims,           only : itaur, itaui, ndyn, tref
    use dims,           only : newsrun
    use Go,             only : TDate, operator(+), operator(-), operator(/), SecondNumber, Get

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter :: rname = mname//'/user_output_column_step'
    integer                     :: sec_in_day
    type(TDate)                 :: tmid
    integer, dimension(6)       :: t_temp

    !__START_SUBROUTINE______________________________________________________

    status = 0

    call Get(tr(1), time6=t_temp)
    ! If we are starting a new day, close the old file and open a new file
    if (all(t_temp(4:6) == 0)) then
        if (.not. newsrun) then
            call output_columnrecord(region)
            call close_columndatafile(region)
        end if
        if (.not. output_file%opened) then
            call open_columndatafile(t_temp, status)
            IF_NOTOK_RETURN(status=1)
        end if
    end if

    ! which time index to put the mixing ratios in?
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call SecondNumber(tmid, sec_in_day)
    tod_index(region) = sec_in_day/column_dhour + 1

    call fill_columnrecord(region, tr)

  end subroutine user_output_column_step

  subroutine close_columndatafile(region)

    use file_netcdf,  only : nc_close

    implicit none

    integer, intent(in) :: region
    character(len=*), parameter :: rname = mname//'/close_columndatafile'

    output_file%exists(region) = .false.
    ! if all the regions have been written, close the file
    if (.not. any(output_file%exists)) then
        call nc_close(output_file%nc_id)
        output_file%opened = .false.
    end if

  end subroutine close_columndatafile

  subroutine open_columndatafile(idate_f, status)

    use dims,          only : nregions, xref, yref, tref, dx, dy
    use dims,          only : isr, jsr, ier, jer, xbeg, xend, ybeg, yend
    use dims,          only : region_name, im, jm, ibeg, jbeg, iend, jend
    use chem_param,    only : ntrace, names
    use global_data,   only : rcF
    use GO,            only : ReadRc
    use misctools,     only : check_dir
    use tm5_geometry,  only : lli
    use netcdf,        only : nf90_put_att, NF90_GLOBAL, NF90_UNLIMITED
    use file_netcdf

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)     :: idate_f(6)
    integer, intent(out)    :: status

    !__LOCAL_VARIABLES______________________________________________________


    character(len = 200)      :: FFilename, attr_name
    integer                   :: n, region, io_status

    character(len=*), parameter         :: rname = mname//'/open_columndatafile'

    !__START_SUBROUTINE______________________________________________________

    status = 0

    if (.not. output_file%opened) then
        write(FFilename, '(a, "/", i4.4, "/", i2.2, "/", a, i4.4, 2(i2.2), ".nc4")') &
            trim(outdir_column), idate_f(1), idate_f(2), trim(outfile_name_prefix), idate_f(1:3)
        call check_dir(FFilename)

        output_file%nc_id = nc_open(FFilename, 'c', status)
        IF_NOTOK_RETURN(status=1)
        
        output_file%opened = .true.
        ! write the global attributes
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'nregions', nregions)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'dx', dx)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'dy', dy)
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'xref', xref(0:nregions))
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'yref', yref(0:nregions))
        io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, 'tref', tref(0:nregions))

        ! create the global dimensions
        call nc_create_dim(output_file%nc_id, 'times', NF90_UNLIMITED)
        !call nc_create_dim(output_file%nc_id, 'times', SEC_PER_DAY/column_dhour)
        call nc_create_dim(output_file%nc_id, 'tracers', ntrace)
        call nc_create_dim(output_file%nc_id, 'idate', 6)

        ! write the tracer names
        do n = 1, ntrace
            write(attr_name, '("tracer_name_", i3.3)') n
            io_status = nf90_put_att(output_file%nc_id, NF90_GLOBAL, trim(attr_name), names(n))
        end do

        if (io_status /= 0) then
            status = 1
            IF_NOTOK_RETURN(status=1)
        end if
    end if

    do region = 1, nregions
        output_file%grp_id(region) = nc_create_group(output_file%nc_id, region_name(region))
        output_file%exists(region) = .true.
        ! write the region-specific attributes
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'im', im(region))
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'jm', jm(region))
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'xbeg', xbeg(region))
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'ybeg', ybeg(region))
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'xend', xend(region))
        io_status = nf90_put_att(output_file%grp_id(region), NF90_GLOBAL, 'yend', yend(region))

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

        if (io_status /= 0) then
            status = 1
            IF_NOTOK_RETURN(status=1)
        end if

    end do

  end subroutine open_columndatafile

  subroutine fill_columnrecord(region, tr)

    use chem_param,   only : ntrace, fscale
    use dims,         only : im, jm, lm, itaur, ndyn, tref, ndyn_max
    use global_data,  only : mass_dat
    use MeteoData,    only : m_dat, phlb_dat, sp_dat
    use Go,           only : TDate, rTotal, operator(-), operator(+), operator(/), Get
    use datetime,     only : date2tau

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)     :: region
    type(TDate), intent(in) :: tr(2)

    !__LOCAL_VARIABLES______________________________________________________

    integer                   :: n
    integer                   :: imr,jmr,lmr, t_temp(6), local_tau
    real, dimension(:,:,:), pointer   :: m, phlb
    real, dimension(:,:,:,:), pointer :: rm
    real                      :: weight
    type(TDate)               :: tmid

    character(len=*), parameter :: rname = mname//'/fill_columnrecord'

    !__START_SUBROUTINE______________________________________________________

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    m    => m_dat(region)%data
    rm   => mass_dat(region)%rm_t
    phlb => phlb_dat(region)%data

    weight = rTotal(tr(2) - tr(1), 'sec') ! from emission_fwd.F90

    do n=1,ntrace
        columnf(region)%col(1:imr,1:jmr,tod_index(region),n) = columnf(region)%col(1:imr,1:jmr,tod_index(region),n) + &
            weight * sum(rm(1:imr,1:jmr,1:lmr,n), dim=3) * fscale(n) / sum(m(1:imr,1:jmr,1:lmr), dim=3)
    end do

    columnf(region)%pres(1:imr,1:jmr,tod_index(region)) = columnf(region)%pres(1:imr,1:jmr,tod_index(region)) + &
        weight * phlb(1:imr,1:jmr,1)
    columnf(region)%nsamples(tod_index(region)) = columnf(region)%nsamples(tod_index(region)) + 1

    ! The 'tau' for this sample should be the midpoint of tr(1) and tr(2)
    tmid = tr(1) + (tr(2)-tr(1))/2.0
    call Get(tmid, time6=t_temp)
    call date2tau(t_temp, local_tau)

    columnf(region)%tau(tod_index(region)) = columnf(region)%tau(tod_index(region)) + &
        weight * local_tau
    columnf(region)%weight(tod_index(region)) = columnf(region)%weight(tod_index(region)) + weight

    nullify(rm)   ! reset pointers
    nullify(m)    ! reset pointers
    nullify(phlb)

  end subroutine fill_columnrecord

  subroutine output_columnrecord(region) ! only called at the end of the day

    use chem_param,   only: ntrace
    use dims,         only: im, jm, lm
    use datetime,     only : tau2date
    use file_netcdf

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)      :: region

    !__LOCAL_VARIABLES______________________________________________________


    integer                   :: n, j
    integer                   :: imr,jmr
    integer, dimension(6)     :: idate_f
    integer                   :: io_status
    character(len=*), parameter :: rname = mname//'/output_columnrecord'

    !__START_SUBROUTINE______________________________________________________


    imr = im(region)
    jmr = jm(region)

    ! divide by nsamples to get average column mixing ratio over column_dhour seconds
    do n=1,size(columnf(region)%col,3)
        columnf(region)%col(:,:,n,:) = columnf(region)%col(:,:,n,:)/columnf(region)%weight(n)
        columnf(region)%pres(:,:,n) = columnf(region)%pres(:,:,n)/columnf(region)%weight(n)
        columnf(region)%tau(n) = columnf(region)%tau(n)/columnf(region)%weight(n)
        call tau2date(int(columnf(region)%tau(n)), columnf(region)%times(:,n))
    end do

    ! turn on netcdf compression
    nc_variables_deflate = .true.
    nc_deflate_level = 5

    call nc_dump_var(output_file%grp_id(region), 'column_mix', (/ 'longitude', 'latitude ', 'times    ', 'tracers  ' /), &
        columnf(region)%col(:,:,:,:), (/ 'long_name' /), (/ 'column-averaged mixing ratio' /))
    call nc_dump_var(output_file%grp_id(region), 'psurf', (/ 'longitude', 'latitude ', 'times    ' /), &
        columnf(region)%pres(:,:,:), (/ 'long_name' /), (/ 'surface pressure in Pa' /))
    call nc_dump_var(output_file%grp_id(region), 'nsamples', (/'times'/), &
        columnf(region)%nsamples(:), (/ 'long_name' /), (/ 'number of samples in each time step' /))
    call nc_dump_var(output_file%grp_id(region), 'sample_times', (/'idate','times'/), &
        columnf(region)%times(:,:))

    ! mixing ratios written, so the mixing ratio array can be re-initialized
    columnf(region)%col = 0.0
    columnf(region)%pres = 0.0
    columnf(region)%nsamples = 0
    columnf(region)%times = 0
    columnf(region)%tau = 0.0
    columnf(region)%weight = 0.0

    ! turn off netcdf compression
    nc_variables_deflate = .false.

  end subroutine output_columnrecord

end module user_output_column
