#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"


module Emission_Adj

    use GO,             only : gol, goErr, goPr, tdate
    use emission_data,  only : source_apply, adj_emissions, tracers_em_info, t_tracer_info
    use emission_fwd,   only : emission_fwd_setup
    use dims,           only : nregions, isr, ier, jsr, jer, itau, region_name
    use chem_param,     only : ntracet
    use os_specs,       only : MAX_FILENAME_LEN
    use global_data,    only : region_dat, mass_dat
    use datetime,       only : tau2date

    implicit none

    private
    public  :: emission_adj_init
    public  :: emission_adj_done
    public  :: emission_adj_setup
    public  :: emission_adj_apply

    character(len=*), parameter        :: mname = 'Emission_Adj'

    ! Adjoint has the structure adj_emis(region)%tracer(itrac)%values(ilon, ilat)
    ! No need to bother with categories: the adjoint is the same for all categories (by definition)
    type t_surfemis
        ! dimensions: lon, lat, tracer
        real, dimension(:, :, :), allocatable   :: values
    end type t_surfemis
    type(t_surfemis), dimension(:), allocatable :: adj_emis

    character(len=MAX_FILENAME_LEN), dimension(:), allocatable     :: filename_adj

    contains

        subroutine emission_adj_init(status)
            integer, intent(out)    :: status
            integer                 :: ireg
            
            allocate(adj_emis(nregions))
            allocate(filename_adj(nregions))
            do ireg = 1, nregions
                allocate(adj_emis(ireg)%values(isr(ireg) : ier(ireg), jsr(ireg) : jer(ireg), ntracet))
            enddo
            status = 0
        end subroutine emission_adj_init


        subroutine emission_adj_done(status)
            integer, intent(out)    :: status
            status = 0
        end subroutine emission_adj_done


        subroutine emission_adj_setup(status)
            integer, intent(out)            :: status
            character(len=*), parameter     :: rname = mname//'/emission_adj_setup'
            call emission_fwd_setup(status)
        end subroutine emission_adj_setup


        subroutine emission_adj_apply(ireg, tr, status)
            
            use go, only : rtotal, operator(-)

            integer, intent(in)     :: ireg
            type(tdate), intent(in) :: tr(2)
            integer, intent(out)    :: status

            integer :: itrac
            character(len=*), parameter :: rname = mname//'/emission_adj_apply'
            real                        :: dtime

            ! timestep emissions
            dtime = abs(rtotal(tr(2) - tr(1), 'sec'))

            do itrac = 1, ntracet
                call emission_adj_apply_tracer(ireg, itrac, tracers_em_info(ireg)%tracer(itrac), dtime)
            enddo

            call write_adj_emis(ireg, dtime, status)
            IF_NOTOK_RETURN(status=1)

        end subroutine emission_adj_apply


        subroutine emission_adj_apply_tracer(ireg, itrac, tracer, dtime)

            integer, intent(in)                     :: ireg, itrac
            real, intent(in)                        :: dtime
            type(t_tracer_info), intent(in)         :: tracer
            real, dimension(:, :, :, :), pointer    :: adj_rm
            real, dimension(:, :, :, :), pointer    :: adj_rzm

            integer                         :: j, i

            adj_rm => mass_dat(ireg)%rm_t
            adj_rzm => mass_dat(ireg)%rzm_t

            do j = jsr(ireg), jer(ireg)
                do i = isr(ireg), ier(ireg)
                    if (region_dat(ireg)%zoomed(i, j) /= ireg) cycle
                    adj_emis(ireg)%values(i, j, itrac) = adj_emis(ireg)%values(i, j, itrac) + adj_rm(i, j, 1, itrac) * dtime
                    adj_emis(ireg)%values(i, j, itrac) = adj_emis(ireg)%values(i, j, itrac) - adj_rzm(i, j, 1, itrac) * dtime
                enddo
            enddo

            nullify(adj_rm, adj_rzm)

        end subroutine emission_adj_apply_tracer

        
        subroutine write_adj_emis(ireg, dtime, status)

            integer, intent(in)                 :: ireg
            integer, intent(out)                :: status
            real, intent(in)                    :: dtime

            ! integer                             :: iday
            integer, dimension(6)               :: idate_mid
            character(len=MAX_FILENAME_LEN)     :: filename

            status = 0

            ! TODO: this is probably not correct in adjoint mode!
            !iday = get_num_days(itaui, itaur(ireg) + ndyn / 4 / tref(ireg))
            !call tau2date(itaur(ireg) + ndyn / 4 / tref(ireg), idate_mid)
            call tau2date(itau - nint(dtime), idate_mid)

            ! One adjoint file should contain the emissions of all tracers for one given region and adjoint time step
            write(filename, '("adjemis.", a, ".", i4.4, 2i2.2, ".nc")'), trim(region_name(ireg)), idate_mid(1), idate_mid(2), idate_mid(3)

            if (filename == filename_adj(ireg)) then 
                ! If we are still in the same adjoint timestep, do nothing (i.e. keep accumulating the adjoint emissions)
                return
            endif

            ! Write the adjoint emissions to a netCDF file (TODO: write the code!)

            ! Reset the adjoint fields
            adj_emis(ireg)%values = 0

            ! Store the filename
            filename_adj(ireg) = filename

        end subroutine write_adj_emis

        
end module emission_adj