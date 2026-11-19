#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"


module chemistry

    use go,             only : TDate, readrc, T_Time_Window, operator(<), NewDate, TIncrDate, operator(+), IncrDate, date_min, operator(>=)
    use global_data,    only : rcf
    use chem_param,     only : ntracet, tracers, tracer_t, react_t, ntlow
    use dims,           only : im, jm, lm, isr, ier, jsr, jer, nregions, revert
    use meteo,          only : pclim_dat
    use file_netcdf
    use grid_type_ll,   only : init_grid => init, tllgridinfo
    use grid_type_hyb,  only : init_levels => init, tlevelinfo
    use grid_type_gg,   only : init_gg => init
    use grid_interpol,  only : aver
    use num,            only : intervalquad_const
    use grid_3d,        only : fill3d
    use go,             only : gol, goerr
    use tm5_geometry,   only : lli, levi
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

    implicit none

    public :: chemistry_init, chemistry_step, chemistry_done, read_chemistry_fields
    private

    character(len=*), parameter   ::  mname = 'chemistry'

    type oh_input
       integer                                  :: year = -1
       real, dimension(:, :, :, :), allocatable :: conc ![molec / cm**3]
       type(tllgridinfo)                        :: lli
       type(tlevelinfo)                         :: levi
       logical                                  :: initialized = .false.
    end type oh_input
    ! integer                                  :: cur_year=-1
    ! real, dimension(:, :, :, :), allocatable :: cur_ohfield_in
    ! type(tllgridinfo)                        :: cur_lli_in
    ! type(tlevelinfo)                         :: cur_levi_in
    type(oh_input) :: oh_cams
    type(oh_input) :: oh_spiva

    contains

        subroutine chemistry_init(status)
            integer, intent(out)    :: status
            status = 0
            ! TODO: this is duplicate of init_chem in chem_param.F90 ...
        end subroutine chemistry_init


        subroutine chemistry_step(region, period, status)

            use global_data,    only : mass_dat, region_dat
            use dims,           only : isr, ier, jsr, jer
            use go,             only : operator(-), rtotal

            integer, intent(in)                     :: region
            type(TDate), dimension(2), intent(in)   :: period
            integer, intent(out)                    :: status
            integer     :: itr
            real, dimension(:, :, :), pointer       :: rm, rxm, rym, rzm
            real, dimension(:, :, :), allocatable   :: loss_rate
            integer     :: is, ie, js, je
            real        :: dtime
            character(len=*), parameter             :: rname = &
                 mname//'/chemistry_step'

            dtime = abs(rtotal(period(2) - period(1), 'sec'))

            is = isr(region) ; ie = ier(region)
            js = jsr(region) ; je = jer(region)

            !-- MVO,2026-03-24:
            !   - allocate loss_rate once here, rather than for
            !     every tracer internally within get_total_loss_rate
            allocate(loss_rate(im(region), jm(region), lm(region)))
            
            do itr = 1, ntracet
                if (tracers(itr)%has_chem) then

                    call get_total_loss_rate(region, period, tracers(itr), loss_rate, status)
                    IF_NOTOK_RETURN(status=1)

                    rm => mass_dat(region)%rm_t(is : ie, js : je, :, itr)
                    
                    rxm => mass_dat(region)%rxm_t(is : ie, js : je, :, itr)
                    rym => mass_dat(region)%rym_t(is : ie, js : je, :, itr)
                    rzm => mass_dat(region)%rzm_t(is : ie, js : je, :, itr)

                    rm = rm * (1 - loss_rate(is : ie, js : je, :) * dtime)
#ifdef slopes
                    rxm = rxm * (1 - loss_rate(is: ie, js: je, :) * dtime)
                    rym = rym * (1 - loss_rate(is: ie, js: je, :) * dtime)
                    rzm = rzm * (1 - loss_rate(is: ie, js: je, :) * dtime)
#endif

                    nullify(rm, rxm, rym, rzm)
                end if
            end do
            !-- dispose memory 
            deallocate(loss_rate)
            status = 0
        end subroutine chemistry_step


        subroutine chemistry_done(status)
            integer, intent(out)    :: status
            status = 0
        end subroutine chemistry_done


        subroutine read_chemistry_fields(direction, status)
            integer, intent(in)     :: direction
            integer, intent(out)    :: status
            status = 0
        end subroutine read_chemistry_fields

        subroutine get_total_loss_rate(region, period, tracer, loss_rate, status)
            ! Get the total tracer loss (i.e. sum of reaction rate * mass * dtime) for a given tracer

            Type(TDate), dimension(2), intent(in)   :: period
            integer, intent(in)                     :: region
            type(tracer_t), intent(inout)           :: tracer
            ! real, dimension(:, :, :), allocatable, intent(out)   :: loss_rate
            real, dimension(:, :, :), intent(out)   :: loss_rate
            integer, intent(out)                    :: status

            integer                                 :: ireac
            real, dimension(:, :, :), allocatable   :: lrate
            character(len=*), parameter             :: rname = &
                 mname//'/get_total_loss_rate'

!            allocate(loss_rate(im(region), jm(region), lm(region)))

            loss_rate = 0

            do ireac = 1, tracer%nreact

                ! Loss rate (s-1) of this reaction in each grid cell
                call get_lossrate(tracer%reactions(ireac), period, region, lrate, status)
                IF_NOTOK_RETURN(status=1)

                ! Make sure the loss rate is set to 0 outside the region/vertical domain of application of the reaction
                call apply_l_domain(lrate, tracer%reactions(ireac), region, status)
                IF_NOTOK_RETURN(status=1)

                ! Accumulate
                loss_rate = loss_rate + lrate

            end do

        end subroutine get_total_loss_rate


        subroutine get_lossrate(reaction, period, region, lrate, status)
            ! Loss rate (s-1) of the tracer due to this reaction, in each grid cell of the region
            ! (not yet restricted to the domain of application of the reaction)
            type(react_t), intent(inout)                        :: reaction
            type(TDate), dimension(2), intent(in)               :: period
            integer, intent(in)                                 :: region
            real, dimension(:, :, :), allocatable, intent(out)  :: lrate
            integer, intent(out)                                :: status

            character(len=*), parameter                 :: rname = mname//'/get_lossrate'

            ! Read the input data (also tells us what they contain: loss rates or concentrations)
            call read_fields(reaction, period, status)
            IF_NOTOK_RETURN(status=1)

            if (reaction%lossrate) then
                ! The input data are already loss rates
                lrate = reaction%conc(region)%values * reaction%scalef
            else
                ! The input data are concentrations of the reactive species: loss rate = k(T) * concentration
                lrate = get_rrate_field(reaction, region) * reaction%conc(region)%values
            end if
        end subroutine get_lossrate


        function get_rrate_field(reaction, region) result(field3d)
            ! Return the actual reaction rate in each grid cell (which is a function of temperature)
            ! The reaction rate will be set to 0 outside the current region (so chemistry is computed only once)

            use meteodata,  only : temper_dat

            type(react_t), intent(in)               :: reaction
            integer, intent(in)                     :: region
            real, dimension(:, :, :), allocatable   :: field3d
            integer                                 :: ilon, ilat, ilev, itemp

            allocate(field3d(im(region), jm(region), lm(region)))
            field3d = 0.

            do ilev = 1, lm(region)
                do ilat = jsr(region), jer(region)
                    do ilon = isr(region), ier(region)
                        ! itemp = nint(temper_dat(region)%data(ilon, ilat, ilev) - real(ntlow))
                        field3d(ilon, ilat, ilev) = reaction%rate(int(temper_dat(region)%data(ilon, ilat, ilev))) * reaction%scalef
                    end do
                end do
            end do

        end function get_rrate_field


        subroutine apply_l_domain(field, reaction, region, status)
            ! domain :
            ! total  = total atmosphere
            ! tropo  = troposphere only
            ! strato = stratosphere only
            ! set loss rates outside specified domain to zero

            use tm5_geometry,   only : lli, levi
            use global_data,    only : region_dat
            use meteo,          only : pclim_dat    ! pressure climatology
            use go,             only : gol, goerr

            real, dimension(:, :, :), intent(inout)     :: field
            type(react_t), intent(in)                   :: reaction
            integer, intent(in)                         :: region
            integer, intent(out)                        :: status

            integer                                     :: ilon, ilat, ilev
            real                                        :: lat
            real                                        :: pres, pres_tropopause
            character(len=*), parameter                 :: rname = mname//'/Apply_L_Domain'

            do ilat = 1, jm(region)
                ! current latitude:
                lat = lli(region)%lat(ilat)  ! in radians
                ! parameterisation for tropopause:
                !~ lowest at 300 hPa at poles, highest at equator at 85 hPa :
                !pres_tropopause = 300.0e2 - 215.0e2*cos(lat)
                !~ MCK, 2017-09: fixed, should be cos**2 ?
                !  Lawrence et al., 2001, ACP, doi:10.5194/acp-1-37-2001
                pres_tropopause = 300.0e2 - 215.0e2 * cos(lat) ** 2
                ! loop over longitudes:

                do ilon = 1, im(region)

                    ! If we are outside the current region, set the field to 0 throughout the column
                    if (region_dat(region)%zoomed(ilon, ilat) /= region) then
                        field(ilon, ilat, :) = 0
                        cycle
                    end if

                    ! Otherwise, see if we are in the right vertical domain
                    if (trim(reaction%domain) == 'total') cycle

                    do ilev = 1, lm(region)
                        pres = levi%fa(ilev) + levi%fb(ilev) * pclim_dat(region)%data(ilon, ilat, 1)
                        select case (trim(reaction%domain))
                            case ('tropo')
                                if (pres < pres_tropopause) field(ilon, ilat, ilev) = 0
                            case ('strato')
                                if (pres >= pres_tropopause) field(ilon, ilat, ilev) = 0
                            case default
                                write (gol,'("unsuported domain :",a)') trim(reaction%domain); call goErr
                                TRACEBACK; status=1; return
                        end select
                    end do
                end do
            end do
        end subroutine apply_l_domain


        subroutine read_fields(reaction, period, status)
            ! Make sure that the input data of the reaction are in memory (reaction%conc) for the requested period.
            ! Depending on the format of the input file (reaction%version), these are either
            !   - concentrations of the species reacting with the tracer (reaction%lossrate = .false.), or
            !   - loss rates of the tracer, in s-1 (reaction%lossrate = .true.)
            ! The reader sets reaction%lossrate accordingly.
            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period
            integer, intent(out)                        :: status

            character(len=*), parameter                 :: rname = mname//'/read_fields'

            status = 0

            ! Check if we need to read new data (i.e. new period). Otherwise return what's already in memory
            if ((revert == 1) .and. (period(1) < reaction%data_period%t2)) return
            if ((revert == -1) .and. (period(2) >= reaction%data_period%t1)) return

            select case (reaction%version)
                case ('cams        ')
                    call get_conc_field_cams(reaction, period, status)
                case ('spivakovsky ')
                    call get_conc_field_spivakovsky(reaction, period, status)
                case ('bruehl_echam')
                    call get_lossrate_field_bruehl_echam(reaction, period, status)
                case default
                    write (gol,'("unsupported reaction version: ",a)') trim(reaction%version); call goErr
                    TRACEBACK; status=1; return
            end select
            IF_NOTOK_RETURN(status=1)
        end subroutine read_fields


        subroutine get_conc_field_cams(reaction, period, status)

            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period

            real, dimension(:, :, :, :), allocatable    :: ohfield_in
            real, dimension(:), allocatable             :: hyb_a, hyb_b
            character(len=*), parameter                 :: rname = mname//'/get_conc_field_cams'
            integer, intent(out)                        :: status

            integer             :: ncf
            character(len=256)  :: reacfile
            integer             :: nlay
            integer             :: region
            integer             :: i_month
            type(TDate)         :: period_start

            period_start = date_min(period(1), period(2))

            ! The rate constant is mandatory for these reactions (default in chem_param is NaN)
            if (ieee_is_nan(reaction%parameters(1))) then
                write (gol,'("reaction ",a," reads concentrations and requires a rate")') trim(reaction%name); call goErr
                TRACEBACK; status=1; return
            end if

            reaction%climatology = .false.
            reaction%data_timestep = '1m'
            reaction%data_period = get_new_period(period_start, reaction%data_timestep)

            !
            ! int64 month(month) ;
            ! double OH(month, pressure_level, lat, lon) ;
            !     OH:_FillValue = NaN ;
            !     OH:units = "molec / cm ** 3" ;
            !
            i_month = reaction%data_period%t1%month

            ! Load required data from the netCDF file
            ! assumption is that OH files are annual,
            ! only re-loading file in case the year has changed!
            if (oh_cams%year /= reaction%data_period%t1%year) then
                write(reacfile, '(a,i4,a)') trim(reaction%file)//'_', reaction%data_period%t1%year, '_monthly.nc'
                ncf = nc_open(reacfile, 'r', status)
                IF_NOTOK_RETURN(status=1)
                ohfield_in = nc_read_var(ncf, 'OH', status)     ! needs to be pre-converted in molec/cm3!
                IF_NOTOK_RETURN(status=1)
                hyb_a = nc_read_var(ncf, 'z_a_grid', status)
                IF_NOTOK_RETURN(status=1)
                hyb_b = nc_read_var(ncf, 'z_b_grid', status)
                IF_NOTOK_RETURN(status=1)
                call nc_close(ncf)

                ! First timestep, need to initialize the grid
                if (.not. oh_cams%initialized) then
                    ! Create the horizontal and vertical grids for the input file
                    nlay = size(hyb_a) - 1
                    call init_grid(oh_cams%lli, -179.5, 1.0, 360, -89.5, 1.0, 180, status)
                    IF_NOTOK_RETURN(status=1)
                    call init_levels(oh_cams%levi, nlay, hyb_a, hyb_b, status, name='OH')
                    IF_NOTOK_RETURN(status=1)
                    oh_cams%initialized = .true.
                endif

                !-- update (global) CAMS OH input
                oh_cams%year = reaction%data_period%t1%year
                !--  move CAMS OH to global variable
                call move_alloc(ohfield_in, oh_cams%conc)
            endif

            !-- MVO-TODO, 2026-03-24
            !   - next improvement would be to store the regionally
            !     regridded values to global fields,
            !     and then have only pointers to those at the
            !     tracers!
            !
            ! Do the regridding, for each region
            do region = 1, nregions
                if (.not. allocated(reaction%conc(region)%values)) allocate(reaction%conc(region)%values(im(region), jm(region), lm(region)))

                reaction%conc(region)%values = 0.
                ! print*, 'DEBUG@'//rname//' region=',region,&
                !      'i_period=',i_period, 'i_month=',i_month,&
                !      'mean(ohfield_in)=',&
                !      sum(ohfield_in(:,:,:,i_month))/size(ohfield_in(:,:,:,i_month))
                call fill3d( &
                    lli(region), levi, 'n', pclim_dat(region)%data(:, :, 1), &
                    reaction%conc(region)%values, oh_cams%lli, oh_cams%levi, &
                    oh_cams%conc(:, :, :, i_month), &
                    'mass-aver', status &
                )
                IF_NOTOK_RETURN(status=1)
            enddo

        end subroutine get_conc_field_cams


        subroutine get_conc_field_spivakovsky(reaction, period, status)

            ! Return the concentration field of a reactive species (in units of ...), during the requested time interval
            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period
            integer, intent(out)                        :: status

            real, dimension(:, :, :, :), allocatable    :: ohfield_in
            character(len=*), parameter                 :: rname = mname//'/get_conc_field_spivakovsky'
            integer             :: ncf
            integer             :: region
            type(TDate)         :: period_start
            
            period_start = date_min(period(1), period(2))
            
            ! Hard-coded stuff specific to the Spivakovsky OH field. Code below needs to be revised when more
            ! OH fields are implemented
            ! The rate constant is mandatory for these reactions (default in chem_param is NaN)
            if (ieee_is_nan(reaction%parameters(1))) then
                write (gol,'("reaction ",a," reads concentrations and requires a rate")') trim(reaction%name); call goErr
                TRACEBACK; status=1; return
            end if

            reaction%climatology = .true.
            reaction%data_timestep = '1m'
            reaction%data_period = get_new_period(period_start, '1m')

            ! GM, 14 oct 2024:
            ! Code adapted from https://sourceforge.net/p/tm5/cy3_4dvar/ci/d13c-ch4/tree/proj/tracer/d13CH4/src/import_chemistry_fields.F90#l232
            ! works well for the Spivakovski fields, but only for these, since the vertical coordinates need to follow that of TM5

            ! Spivakovski tropospheric OH in cm ** -3
            ! MACC OH in molec . cm ** -3
            ! Lrate(i,j,l) = rrates(kCH4OH, itemp(i,j,l) * OHconc(i,j,l))

            if (.not. oh_spiva%initialized) then
                ncf = nc_open(reaction%file, 'r', status)
                IF_NOTOK_RETURN(status=1)
                ohfield_in = nc_read_var(ncf, 'field', status)
                IF_NOTOK_RETURN(status=1)
                call nc_close(ncf)

                call  move_alloc(ohfield_in, oh_spiva%conc)

                ! Create coordinates for the OH field
                call init_grid(oh_spiva%lli, -179.5, 1.0, 360, -89.5, 1.0, 180, status)
                IF_NOTOK_RETURN(status=1)
                call init_levels(oh_spiva%levi, 'tm60', status)
                IF_NOTOK_RETURN(status=1)
                oh_spiva%initialized = .true.
            endif

            ! Apply regridding
            do region = 1, nregions
                if (.not. allocated(reaction%conc(region)%values)) allocate(reaction%conc(region)%values(im(region), jm(region), lm(region)))
                ! If we need to load new data:
                reaction%conc(region)%values = 0.

                call fill3d( &
                    lli(region), levi, 'n', pclim_dat(region)%data(:, :, 1), &
                    reaction%conc(region)%values, oh_spiva%lli, oh_spiva%levi, &
                    oh_spiva%conc(:, :, :, reaction%data_period%t1%month), &
                    'mass-aver', status)
                IF_NOTOK_RETURN(status=1)
            enddo

        end subroutine get_conc_field_spivakovsky


        subroutine get_lossrate_field_bruehl_echam(reaction, period, status)
            ! --- in/out ---------------------------------

            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period
            integer, intent(out)                        :: status

            ! --- const ----------------------------------

            character(len=*), parameter                 :: rname = mname//'/get_lossrate_field_bruehl_echam'

            !--- local -----------------------------------

            integer                          ::  region
            integer                          ::  month
            integer                          ::  ncf
            logical                          ::  exist
            character(len=32)                ::  units
            real(4), allocatable             ::  lats(:)
            real, allocatable                ::  lats_rev(:)
            real, allocatable                ::  field(:,:)
            real(4), allocatable             ::  field_3d(:,:,:)
            real, allocatable                ::  rates(:,:)
            real, allocatable                ::  fpres(:,:)
            real, allocatable                ::  hpres(:,:)
            real, allocatable                ::  rates_llz(:,:,:)
            real, allocatable                ::  hpres_llz(:,:,:)
            integer                          ::  imr, jmr, lmr
            integer                          ::  i, j, k
            integer                          ::  ilev
            real                             ::  p1, p2
            integer                          ::  ilast
            type(TDate)                      ::  period_start

            !--- begin -----------------------------------

            period_start = date_min(period(1), period(2))
            reaction%lossrate = .true.
            reaction%climatology = .true.
            reaction%data_timestep = '1m'
            reaction%data_period = get_new_period(period_start, reaction%data_timestep)

            ! extract current month:
            month = reaction%data_period%t1%month

            associate ( ds => reaction%lossrate_field_gg )

            ! check ...
            ds%filename = reaction%file
            inquire( file=trim(ds%filename), exist=exist )
            if ( .not. exist ) then
              write (gol,'("Bruehl output file not found : ",a)') trim(ds%filename); call goErr
              TRACEBACK; status=1; return
            end if

            ! open file:
            ncf = nc_open(trim(ds%filename), 'r', status)
            IF_NOTOK_RETURN(status=1)

            ! no longitudes read yet ?
            if ( .not. ds%grid_is_initialized ) then
              ! latitudes:
              ds%nlat = nc_get_dim(ncf, 'lat', status)
              IF_NOTOK_RETURN(status=1)
              ! read:
              lats = nc_read_var(ncf, 'lat', status)
              IF_NOTOK_RETURN(status=1)
              ! units:
              units = nc_get_attr(ncf, 'lat', 'units', status)
              IF_NOTOK_RETURN(status=1)
              ! convert if necessary:
              select case ( trim(units) )
                case ( 'degree', 'degree_north' )
                case default
                  write (gol,'("unsupported latitude units : ",a)') trim(units); call goErr
                  TRACEBACK; status=1; return
              end select
              ! gaussian grid, N=number of lats on hemisphere average,
              ! for N=64/2 it is to the values found in the ECHAM/MESSY files
              ! for guassian grids the order is always north->south, so need to switch:
              call init_gg( ds%ggi, ds%nlat/2, .false., status )
              IF_NOTOK_RETURN(status=1)
              ! reversed ax:
              allocate( lats_rev(ds%nlat) )
              do j = 1, ds%nlat
                lats_rev(j) = lats(ds%nlat+1-j)
              end do
              ! check ...
              if ( any( abs( ds%ggi%lat_deg - lats_rev ) > 1.0e-4 ) ) then
                write (gol,'("lats in assumed gaussian grid do not match with file ...")'); call goErr
                TRACEBACK; status=1; return
              end if
              ! clear:
              deallocate( lats )
              deallocate( lats_rev )

              ! levels:
              ds%nlev = nc_get_dim(ncf, 'level', status)
              IF_NOTOK_RETURN(status=1)

              ! reset flag:
              ds%grid_is_initialized = .true.
            end if

            ! storage:
            allocate( field(ds%nlat,ds%nlev) )

            ! read:
            field_3d = nc_read_var(ncf, 'Lrate', status)
            IF_NOTOK_RETURN(status=1)
            field = field_3d(:,:,month)
            ! units:
            units = nc_get_attr(ncf, 'Lrate', 'units', status)
            IF_NOTOK_RETURN(status=1)
            select case ( trim(units) )
              case ( 's-1' ) ;
                ! no conversion
              case default
                write (gol,'("unsupported rates units : ",a)') trim(units); call goErr
                TRACEBACK; status=1; return
            end select
            ! storage for rates:
            allocate( rates(ds%nlat,ds%nlev) )
            ! revert:
            do j = 1, ds%nlat
              rates(j,:) = field(ds%nlat+1-j,:)
            end do

            ! read:
            field_3d = nc_read_var(ncf, 'pressure', status)
            IF_NOTOK_RETURN(status=1)
            field = field_3d(:,:,month)
            ! units:
            units = nc_get_attr(ncf, 'pressure', 'units', status)
            IF_NOTOK_RETURN(status=1)
            select case ( trim(units) )
              case ( 'Pa' ) ;
                ! no conversion
              case ( 'hPa' )
                ! hPa to Pa:
                field = field * 1.0e2  ! Pa
              case default
                write (gol,'("unsupported pressure units : ",a)') trim(units); call goErr
                TRACEBACK; status=1; return
            end select
            ! storage for pressure, copy while reverting latitudes:
            allocate( fpres(ds%nlat,1:ds%nlev) )
            do j = 1, ds%nlat
              fpres(j,:) = field(ds%nlat+1-j,:)
            end do
            ! extend to half level pressures:
            allocate( hpres(ds%nlat,0:ds%nlev) )
            hpres(:,0) = fpres(:,1) - 0.5*(fpres(:,2)-fpres(:,1))
            do ilev = 1, ds%nlev-1
              hpres(:,ilev) = 0.5 * ( fpres(:,ilev) + fpres(:,ilev+1) )
            end do
            hpres(:,ds%nlev) = 0.0 

            ! close file:
            call nc_close(ncf)

            ! loop over regions:
            do region = 1, nregions
              ! dimensions:
              imr = im(region) ; jmr = jm(region) ; lmr = lm(region)
              if (.not. allocated(reaction%conc(region)%values)) allocate(reaction%conc(region)%values(imr, jmr, lmr))
              ! average to horizontal target grid, levels original:
              allocate( rates_llz(imr,jmr,ds%nlev) )
              do ilev = 1, ds%nlev
                call aver( ds%ggi, rates(:,ilev), lli(region), rates_llz(:,:,ilev) )
              end do
              ! idem for pressure:
              allocate( hpres_llz(imr,jmr,0:ds%nlev) )
              do ilev = 0, ds%nlev
                call aver( ds%ggi, hpres(:,ilev), lli(region), hpres_llz(:,:,ilev) )
              end do
              ! vertical interpolation:
              do j = 1, jmr
                do i = 1, imr
                  ilast = 1
                  ! loop over target layers;
                  ! use negative pressure to have increasing ax:
                  do k = 1, lmr
                    ! layer half level pressures:
                    !AJS: don't use this! phlb not initialized correctly for adjoint run ...
                    !     also problem that we use current pressue to interpolat monthly average,
                    !     thus differences between forward and adjoint run 
                    !     (might be zero if only stratosphere is used)
                    !p1 = -1.0 * mass_dat(region)%phlb_t(i,j,k  )
                    !p2 = -1.0 * mass_dat(region)%phlb_t(i,j,k+1)
                    ! climatology:
                    p1 = -1.0 * levi%p0(k-1)
                    p2 = -1.0 * levi%p0(k)
                    ! integral over interval [p1,p2], with rates constant withing layer:
                    call intervalquad_const( -1.0*hpres_llz(i,j,:), rates_llz(i,j,:), &
                                                p1, p2, reaction%conc(region)%values(i,j,k), &
                                                ilast, status )
                    IF_NOTOK_RETURN(status=1)
                    ! average:
                    reaction%conc(region)%values(i,j,k) = reaction%conc(region)%values(i,j,k) / (p2-p1)
                  end do  ! layer
                end do  ! i
              end do  ! j
              ! clear:
              deallocate( rates_llz )
              deallocate( hpres_llz )
            end do

            ! clear:
            deallocate( rates )
            deallocate( fpres )
            deallocate( hpres )
            deallocate( field )
            deallocate( field_3d )

            ! store current month:
            ds%month = month

            status = 0

            end associate

        end subroutine get_lossrate_field_bruehl_echam


        function get_new_period(date, tres) result(new_period)
            ! *Very* ad-hoc function to handle the dates in this module. Only the cases that are currently in 
            ! use are implemented ...

            type(TDate), intent(in)         :: date
            type(TIncrDate)                 :: dt
            character(len=2), intent(in)    :: tres
            type(T_Time_Window)             :: new_period

            select case (tres)
                case ('1m')
                    new_period%t1 = NewDate(year=date%year, month=date%month, day=1)
                    if (new_period%t1%month == 12) then
                        new_period%t2 = NewDate(year=date%year + 1, month=1, day=1)
                    else
                        new_period%t2 = NewDate(year=date%year, month=date%month + 1, day=1)
                    endif
                case ('3d')
                   ! dt = IncrDate(day=3)
                    new_period%t1 = NewDate(year=date%year, month=date%month, day=date%day)
                    new_period%t2 = new_period%t1 + IncrDate(day=3)
            end select

        end function get_new_period


end module chemistry
