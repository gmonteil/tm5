#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"


module chemistry

    use go,             only : TDate, readrc
    use global_data,    only : rcf
    use chem_param,     only : ntracet, tracers, tracer_t, react_t, ntlow
    use dims,           only : im, jm, lm

    implicit none

    public :: chemistry_init, chemistry_step, chemistry_done, read_chemistry_fields
    private

    contains

        subroutine chemistry_init(status)
            integer, intent(out)    :: status
            status = 0
            ! TODO: this is duplicate of init_chem in chem_param.F90 ...
        end subroutine chemistry_init


        subroutine chemistry_step(region, period, status)
            integer, intent(in)                     :: region
            type(TDate), dimension(2), intent(in)   :: period
            integer, intent(out)                    :: status
            integer     :: itr

            do itr = 1, ntracet
                print*, tracers(itr)%species
                print*, tracers(itr)%has_chem
                if (tracers(itr)%has_chem) then

                    rm => mass_dat(region)%rm_t(:, :, :, itr)
                    rxm => mass_dat(region)%rxm_t(:, :, :, itr)
                    rym => mass_dat(region)%rym_t(:, :, :, itr)
                    rzm => mass_dat(region)%rzm_t(:, :, :, itr)

                    loss_rate = get_loss_rate(region, period, itr)

                    rm = rm - loss_rate * dtime
                    rxm = rxm * (1 - loss_rate * dtime)
                    rym = rym * (1 - loss_rate * dtime)
                    rzm = rzm * (1 - loss_rate * dtime)

                    nullify(rm, rxm, rym, rzm)
                end if
            end do

            status = 10
        end subroutine chemistry_step


        subroutine chemistry_done(status)
            integer, intent(out)    :: status
            status = 10
        end subroutine chemistry_done


        subroutine read_chemistry_fields(direction, status)
            integer, intent(in)     :: direction
            integer, intent(out)    :: status
            status = 0
        end subroutine read_chemistry_fields


        function get_loss_rate(region, period, tracer) result(loss_rate)
            ! Get the total tracer loss rate (i.e. sum of reaction rate * mass * dtime) for a given tracer

            use go,     only : rtotal

            Type(TDate), dimension(2), intent(in)   :: period
            integer, intent(in)                     :: region
            type(tracer_t), intent(in)              :: tracer
            real, dimension(:, :, :), allocatable   :: loss_rate
            real, dimension(:, :, :), allocatable   :: conc     ! concentration of the species the tracer reacts with
            integer                                 :: ireac
            real                                    :: rtotal

            allocate(loss_rate(im(region), jm(region), lm(region)))

            ! time step for this region
            dtime = abs(rtotal(period(2) - period(1), 'sec'))

            do ireac = 1, tracer%nreac

                ! Get the mass of the species the tracer reacts with
                conc = get_conc_field(tracer%reactions(ireac), period)

                ! Calculate the reaction rate
                rrate = get_rrate_field(tracer%reactions(ireac), region)

                ! Make sure the reaction rate is set to 0 outside the region/vertical domain of application of the reaction
                call apply_l_domain(rrate, tracer%reactions(ireac), region)

                ! Calculate the loss rate
                loss_rate = loss_rate + rrate * conc * dtime
            end do

        end function get_loss_rate


        function get_rrate_field(reaction, region) result(field3d)
            ! Return the actual reaction rate in each grid cell (which is a function of temperature)
            ! The reaction rate will be set to 0 outside the current region (so chemistry is computed only once)

            use global_data,    only : meteo_dat

            type(react_t), intent(in)               :: reaction
            integer, intent(in)                     :: region
            real, dimension(:, :, :), allocatable   :: field3d
            integer                                 :: ilon, ilat, ilev, itemp

            allocate(field3d(im(region), jm(region), lm(region)))
            field3d = 0.

            do ilon = 1, im(region)
                do ilat = 1, jm(region)
                    do ilev = 1, lm(region)
                        itemp = nint(meteo_dat(region)%T(ilon, ilat, ilev) - real(ntlow))
                        field3d(ilon, ilat, ilev) = reaction%rate(itemp)
                    end do
                end do
            end do

        end function get_rrate_field


        subroutine apply_l_domain(field, reaction, region, status)

            use tm5_geometry,   only : lli
            use dims,           only : im, jm, lm
            use global_data,    only : region_dat

            real, dimension(:, :, :), intent(inout)     :: field
            type(react_t), intent(in)                   :: reaction
            integer, intent(in)                         :: region
            integer, intent(in)                         :: status
            integer                                     :: ilon, ilat, ilev
            real                                        :: lat
            real                                        :: pres_tropopause

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
                                write (gol,'("unsuported domain :",a)') trim(domain); call goErr
                                TRACEBACK; status=1; return
                        end select
                    end do
                end do
            end do

            status = 0

        end subroutine apply_l_domain


        function get_conc_field(filename, period) result(field3d)
            ! Return the concentration field of a reactive species (in units of ...), during the requested time interval
            character(len=*), intent(in)            :: filename
            type(TDate), intent(in)                 :: period
            real, dimension(:, :, :), allocatable   :: field3d

            allocate(field3d(im(region), jm(region), lm(region)))
            field3d = 0.

        end function get_conc_field


end module chemistry