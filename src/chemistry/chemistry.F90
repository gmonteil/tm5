#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"


module chemistry

    use go,             only : TDate, readrc, T_Time_Window, operator(<), NewDate, TIncrDate, operator(+), IncrDate
    use global_data,    only : rcf
    use chem_param,     only : ntracet, tracers, tracer_t, react_t, ntlow
    use dims,           only : im, jm, lm, isr, ier, jsr, jer, nregions
    use meteo,          only : pclim_dat
    use file_netcdf
    use grid_type_ll,   only : init_grid => init, tllgridinfo
    use grid_type_hyb,  only : init_levels => init, tlevelinfo
    use grid_3d,        only : fill3d
    use go,             only : gol, goerr
    use tm5_geometry,   only : lli, levi

    implicit none

    public :: chemistry_init, chemistry_step, chemistry_done, read_chemistry_fields
    private

    character(len=*), parameter   ::  mname = 'chemistry'

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
            
            dtime = abs(rtotal(period(2) - period(1), 'sec'))

            is = isr(region) ; ie = ier(region)
            js = jsr(region) ; je = jer(region)

            do itr = 1, ntracet
                if (tracers(itr)%has_chem) then

                    call get_total_loss_rate(region, period, tracers(itr), loss_rate)

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


        subroutine get_total_loss_rate(region, period, tracer, loss_rate)
            ! Get the total tracer loss (i.e. sum of reaction rate * mass * dtime) for a given tracer

            Type(TDate), dimension(2), intent(in)   :: period
            integer, intent(in)                     :: region
            type(tracer_t), intent(inout)           :: tracer
            real, dimension(:, :, :), allocatable, intent(out)   :: loss_rate
            integer                                 :: ireac
            real, dimension(:, :, :), allocatable   :: rrate

            allocate(loss_rate(im(region), jm(region), lm(region)))

            loss_rate = 0

            do ireac = 1, tracer%nreact

                ! Get the mass of the species the tracer reacts with
                call get_conc_field(tracer%reactions(ireac), period, region)

                ! Calculate the reaction rate
                rrate = get_rrate_field(tracer%reactions(ireac), region)

                ! Make sure the reaction rate is set to 0 outside the region/vertical domain of application of the reaction
                call apply_l_domain(rrate, tracer%reactions(ireac), region)

                ! Calculate the loss rate
                loss_rate = loss_rate + rrate * tracer%reactions(ireac)%conc(region)%values

            end do

        end subroutine get_total_loss_rate


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


        subroutine apply_l_domain(field, reaction, region)

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
            integer                                     :: status
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


        subroutine get_conc_field(reaction, period, region)

            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period
            integer, intent(in)                         :: region

            ! Check if we need to read new data (i.e. new period). Otherwise return what's already in memory
            if (period(1) < reaction%data_period%t2) return

            select case (reaction%version)
                case ('cams        ')
                    call get_conc_field_cams(reaction, period)
                case ('spivakovsky ')
                    call get_conc_field_spivakovsky(reaction, period)
            end select

        end subroutine get_conc_field


        subroutine get_conc_field_cams(reaction, period)

            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period

            real, dimension(:, :, :, :), allocatable    :: ohfield_in
            real, dimension(:), allocatable             :: hyb_a, hyb_b 
            character(len=*), parameter                 :: rname = mname//'/get_conc_field_cams'

            integer             :: ncf
            integer             :: status
            character(len=60)   :: reacfile
            type(tllgridinfo)   :: lli_in
            type(tlevelinfo)    :: levi_in
            integer             :: nlay
            integer             :: region

            print*, 'get cams'
            reaction%climatology = .false.
            reaction%data_timestep = '3d'

            print*, 'get periods'
            reaction%data_period = get_new_period(period(1), reaction%data_timestep)

            print*, 'read nc'
            write(reacfile, '(a,i4,a)') trim(reaction%file)//'_', period(1)%year, '.nc'
            print*, trim(reacfile)

            ncf = nc_open(reacfile, 'r', status)
            IF_NOTOK_RETURN(status=1)

            print*, 'read oh'
            ohfield_in = nc_read_var(ncf, 'oh', status)
            IF_NOTOK_RETURN(status=1)

            print*, 'read at'
            hyb_a = nc_read_var(ncf, 'z_a_grid', status)
            IF_NOTOK_RETURN(status=1)
            
            print*, 'read bt'
            hyb_b = nc_read_var(ncf, 'z_b_grid', status)
            IF_NOTOK_RETURN(status=1)
            
            call nc_close(ncf)

            print*, 'init grid'
            call init_grid(lli_in, -179.5, 1.0, 360, -89.5, 1.0, 180, status)
            IF_NOTOK_RETURN(status=1)

            print*, 'init levels'
            call init_levels(levi_in, nlay, hyb_a, hyb_b, status, name='OH', revert=.true.)
            IF_NOTOK_RETURN(status=1)

            print*, 'regrid'
            nlay = size(hyb_a) - 1

            do region = 1, nregions
                if (.not. allocated(reaction%conc(region)%values)) allocate(reaction%conc(region)%values(im(region), jm(region), lm(region)))

                print*, 'init data'
                reaction%conc(region)%values = 0.

          ! mole/(mole air) * mlc/mole * (mole air)/m3 * m3/cm3 = mlc/m3
          !                     Avog        p/(RT)     *  1e-6
          !  data_gp = data_gp * Avog * pres_gp/(Rgas*tmpr_gp) * 1e-6  ! mlc/cm3

            ! Interpolate vertically (and horizontally?)
            ! call init_grid(lli_in, 0, 0.75, 480, -90, 0.75, 241, status)

                print*, 'fill3d'
                call fill3d( &
                    lli(region), levi, 'n', pclim_dat(region)%data(:, :, 1), &
                    reaction%conc(region)%values, lli_in, levi_in, &
                    ohfield_in(:, :, :, reaction%data_period%t1%month), &
                    'mass-aver', status)
            enddo

            print*, 'deallocate'
            IF_NOTOK_RETURN(status=1)
            deallocate(ohfield_in)

            print*, 'read cams oh'

        end subroutine get_conc_field_cams


        subroutine get_conc_field_spivakovsky(reaction, period)! result(field3d)

            ! Return the concentration field of a reactive species (in units of ...), during the requested time interval
            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period
            real(4), dimension(:, :, :, :), allocatable :: tmp
            real, dimension(:, :, :, :), allocatable    :: field4d
            character(len=*), parameter                 :: rname = mname//'/get_conc_field_spivakovsky'
            type(tllgridinfo)   :: lli_in
            type(tlevelinfo)    :: levi_in
            integer             :: ncf
            integer             :: status
            integer             :: region
            
            ! Hard-coded stuff specific to the Spivakovsky OH field. Code below needs to be revised when more
            ! OH fields are implemented
            reaction%climatology = .true.
            reaction%data_timestep = '1m'
            reaction%data_period = get_new_period(period(1), '1m')

            ! GM, 14 oct 2024:
            ! Code adapted from https://sourceforge.net/p/tm5/cy3_4dvar/ci/d13c-ch4/tree/proj/tracer/d13CH4/src/import_chemistry_fields.F90#l232
            ! works well for the Spivakovski fields, but only for these, since the vertical coordinates need to follow that of TM5

            ! Spivakovski tropospheric OH in cm ** -3
            ! MACC OH in molec . cm ** -3
            ! Lrate(i,j,l) = rrates(kCH4OH, itemp(i,j,l) * OHconc(i,j,l))

            ncf = nc_open(reaction%file, 'r', status)
            IF_NOTOK_RETURN(status=1)
            field4d = nc_read_var(ncf, 'field', status)
            IF_NOTOK_RETURN(status=1)

            ! Create coordinates for the OH field
            call init_grid(lli_in, -179.5, 1.0, 360, -89.5, 1.0, 180, status)
            IF_NOTOK_RETURN(status=1)
            call init_levels(levi_in, 'tm60', status)
            IF_NOTOK_RETURN(status=1)

            ! Apply regridding
            do region = 1, nregions
                if (.not. allocated(reaction%conc(region)%values)) allocate(reaction%conc(region)%values(im(region), jm(region), lm(region)))
                ! If we need to load new data:
                reaction%conc(region)%values = 0.

                call fill3d( &
                    lli(region), levi, 'n', pclim_dat(region)%data(:, :, 1), &
                    reaction%conc(region)%values, lli_in, levi_in, &
                    field4d(:, :, :, reaction%data_period%t1%month), &
                    'mass-aver', status)
                IF_NOTOK_RETURN(status=1)
            enddo
            deallocate(field4d)

        end subroutine get_conc_field_spivakovsky

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
