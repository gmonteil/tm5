#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"


module chemistry

    use go,             only : TDate, readrc, T_Time_Window, operator(<), NewDate
    use global_data,    only : rcf
    use chem_param,     only : ntracet, tracers, tracer_t, react_t, ntlow
    use dims,           only : im, jm, lm, isr, ier, jsr, jer

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

            use global_data, only   : mass_dat
            use dims,                    only : isr, ier, jsr, jer

            integer, intent(in)                     :: region
            type(TDate), dimension(2), intent(in)   :: period
            integer, intent(out)                    :: status
            integer     :: itr
            real, dimension(:, :, :), pointer       :: rm, rxm, rym, rzm
            real, dimension(:, :, :), allocatable   :: loss
            integer     :: is, ie, js, je

            is = isr(region) ; ie = ier(region)
            js = jsr(region) ; je = jer(region)

            do itr = 1, ntracet
                if (tracers(itr)%has_chem) then

                    rm => mass_dat(region)%rm_t(is:ie, js:je, :, itr)
                    rxm => mass_dat(region)%rxm_t(is:ie, js:je, :, itr)
                    rym => mass_dat(region)%rym_t(is:ie, js:je, :, itr)
                    rzm => mass_dat(region)%rzm_t(is:ie, js:je, :, itr)

                    call get_loss(region, period, tracers(itr), loss)

                    rm = rm - loss
                    rxm = rxm * (1 - loss)
                    rym = rym * (1 - loss)
                    rzm = rzm * (1 - loss)

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


        subroutine get_loss(region, period, tracer, loss)
            use go, only : rtotal, operator(-)

            ! Get the total tracer loss (i.e. sum of reaction rate * mass * dtime) for a given tracer

            Type(TDate), dimension(2), intent(in)   :: period
            integer, intent(in)                     :: region
            type(tracer_t), intent(inout)           :: tracer
            real, dimension(:, :, :), allocatable, intent(out)   :: loss
            integer                                 :: ireac
            real                                    :: dtime
            real, dimension(:, :, :), allocatable   :: rrate

            allocate(loss(im(region), jm(region), lm(region)))

            ! time step for this region
            dtime = abs(rtotal(period(2) - period(1), 'sec'))

            loss = 0

            do ireac = 1, tracer%nreact

                ! Get the mass of the species the tracer reacts with
                call get_conc_field(tracer%reactions(ireac), period, region)

                ! Calculate the reaction rate
                rrate = get_rrate_field(tracer%reactions(ireac), region)

                ! Make sure the reaction rate is set to 0 outside the region/vertical domain of application of the reaction
                call apply_l_domain(rrate, tracer%reactions(ireac), region)

                ! Calculate the loss rate
!                print*, ireac, minval(rrate), maxval(rrate), dtime
                loss = loss + rrate * tracer%reactions(ireac)%data * dtime

            end do

        end subroutine get_loss


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
                        field3d(ilon, ilat, ilev) = reaction%rate(int(temper_dat(region)%data(ilon, ilat, ilev)))
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


        subroutine get_conc_field(reaction, period, region)! result(field3d)

            use meteo,          only : pclim_dat
            use file_netcdf
            use grid_type_ll,   only : init_grid => init, tllgridinfo
            use grid_type_hyb,  only : init_levels => init, tlevelinfo
            use grid_3d,        only : fill3d
            use go,             only : gol, goerr
            use tm5_geometry,   only : lli, levi

            ! Return the concentration field of a reactive species (in units of ...), during the requested time interval
            type(react_t), intent(inout)                :: reaction
            type(TDate), dimension(2), intent(in)       :: period
            integer, intent(in)                         :: region
            real(4), dimension(:, :, :, :), allocatable :: tmp
            real, dimension(:, :, :, :), allocatable    :: field4d
            character(len=*), parameter                 :: rname = mname//'/get_conc_field'
            type(tllgridinfo)   :: lli_in
            type(tlevelinfo)    :: levi_in
            integer             :: ncf
            integer             :: status
            
            ! Check if we need to read new data (i.e. new period). Otherwise return what's already in memory
            !if (period(1) < reaction%data_period%t1) return

            ! Hard-coded stuff specific to the Spivakovsky OH field. Code below needs to be revised when more
            ! OH fields are implemented
            reaction%climatology = .true.
            reaction%data_timestep = 'm'
            reaction%data_period = get_new_period(period(1), 'm')
            if (.not. allocated(reaction%data)) allocate(reaction%data(im(region), jm(region), lm(region)))

            ! If we need to load new data:
            reaction%data = 0.

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
            !field4d = tmp

            ! Create coordinates for the OH field
            call init_grid(lli_in, -179.5, 1.0, 360, -89.5, 1.0, 180, status)
            IF_NOTOK_RETURN(status=1)
            call init_levels(levi_in, 'tm60', status)
            IF_NOTOK_RETURN(status=1)
            call fill3d( &
                    lli(region), levi, 'n', pclim_dat(region)%data(:, :, 1), &
                    reaction%data, lli_in, levi_in, &
                    field4d(:, :, :, reaction%data_period%t1%month), &
                    'mass-aver', status)
            IF_NOTOK_RETURN(status=1)
            deallocate(field4d)

        end subroutine get_conc_field

        function get_new_period(date, tres) result(new_period)
            ! Very ad-hoc function to handle the dates in this module. Only the cases that are currently in 
            ! use are implemented ...

            type(TDate), intent(in)         :: date
            character(len=1), intent(in)    :: tres
            type(T_Time_Window)             :: new_period

            select case (tres)
                case ('m')
                    new_period%t1 = NewDate(year=date%year, month=date%month, day=1)
                    if (new_period%t1%month == 12) then
                        new_period%t2 = NewDate(year=date%year + 1, month=1, day=1)
                    else
                        new_period%t2 = NewDate(year=date%year, month=date%month + 1, day=1)
                    endif
            end select

        end function get_new_period


end module chemistry