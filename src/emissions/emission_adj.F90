#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"


module Emission_Adj

    use GO,             only : gol, goErr, goPr, tdate, rtotal, operator(-)
    use emission_data,  only : source_apply, adj_emissions, tracers_em_info, t_tracer_info
    use emission_fwd,   only : emission_fwd_setup
    use dims,           only : nregions, isr, ier, jsr, jer, itau, region_name, itaui, itaur, ndyn, tref
    use chem_param,     only : ntracet
    use os_specs,       only : MAX_FILENAME_LEN
    use global_data,    only : region_dat, mass_dat
    use datetime,       only : tau2date, get_num_days
    use tm5_geometry,   only : lli

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

    logical :: adjoint_output_is_initialized = .false.

    contains

        subroutine emission_adj_init(status)
            integer, intent(out)    :: status
            integer                 :: ireg

            allocate(adj_emis(nregions))
            do ireg = 1, nregions
                allocate(adj_emis(ireg)%values(lli(ireg)%nlon, lli(ireg)%nlat, ntracet))
                !allocate(adj_emis(ireg)%values(isr(ireg) : ier(ireg), jsr(ireg) : jer(ireg), ntracet))
                adj_emis(ireg)%values = 0.
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

            integer, intent(in)     :: ireg
            type(tdate), intent(in) :: tr(2)
            integer, intent(out)    :: status

            integer :: itrac
            character(len=*), parameter :: rname = mname//'/emission_adj_apply'

            do itrac = 1, ntracet
                call emission_adj_apply_tracer(ireg, itrac, tracers_em_info(ireg)%tracer(itrac), tr)
            enddo

            call write_adj_emis(ireg, tr, status)
            IF_NOTOK_RETURN(status=1)
        end subroutine emission_adj_apply


        subroutine emission_adj_apply_tracer(ireg, itrac, tracer, tr)

            integer, intent(in)                     :: ireg, itrac
            type(tdate), dimension(2), intent(in)   :: tr
            type(t_tracer_info), intent(in)         :: tracer
            real, dimension(:, :, :, :), pointer    :: adj_rm
            real, dimension(:, :, :, :), pointer    :: adj_rzm

            integer :: j, i
            real    :: dtime

            ! timestep emissions
            dtime = abs(rtotal(tr(2) - tr(1), 'sec'))

            adj_rm => mass_dat(ireg)%rm_t
            adj_rzm => mass_dat(ireg)%rzm_t

            ! Loop over categories
            ! No loop over categories here, as all categories for one tracer share the same adjoint
            do j = jsr(ireg), jer(ireg)
                do i = isr(ireg), ier(ireg)
                    if (region_dat(ireg)%zoomed(i, j) /= ireg) cycle
                    adj_emis(ireg)%values(i, j, itrac) = adj_emis(ireg)%values(i, j, itrac) + adj_rm(i, j, 1, itrac) * dtime
                    adj_emis(ireg)%values(i, j, itrac) = adj_emis(ireg)%values(i, j, itrac) - adj_rzm(i, j, 1, itrac) * dtime
                enddo
            enddo
            print*, "em_adj_apply_tr", ireg, itrac, sum(adj_rm(:, :, :, itrac))

            nullify(adj_rm, adj_rzm)
        end subroutine emission_adj_apply_tracer


        subroutine write_adj_emis(ireg, tr, status)

          use file_netcdf,    only : nc_open, nc_close, nc_create_dim, nc_dump_var
          use global_data,    only : RcF, outdir
          use go,             only : pathsep, ReadRc
          use os_specs,       only : MAX_FILENAME_LEN, MAX_RCKEY_LEN
          use chem_param,     only : names
          use dims,           only : itaue, itaui, itau
          use datetime,       only : tau2date

          ! This is called by "emission_adj_apply" at the end of each timestep (?)
          !
          ! In the forward we read a new emis if the *current* filename is different from the *previous* filename
          ! Here we need to write a new adj_emis if the *current* filename is different than the *next* filename (otherwise, we accumulate)

          integer, intent(in)                     :: ireg
          integer, intent(out)                    :: status
          type(tdate), dimension(2), intent(in)   :: tr

          character(len=*), parameter             :: rname = mname//'/write_adj_emis'
          character(len=MAX_FILENAME_LEN)         :: filename
          character(len=MAX_RCKEY_LEN)            :: subdir_adjemis
          logical                                 :: dir_exist
          integer                                 :: fid
          real    :: nsec

          status = 0

          ! Unless the current step is exactly at midnight, just keep accumulating adjoint fields
          if ((tr(2)%hour + tr(2)%min + tr(2)%sec) > 0) return

          ! Calculate the number of sec since the start of the simulation
          nsec = itaue - itau - abs(rtotal(tr(2) - tr(1), 'sec'))

          call ReadRc(rcF, 'output.adjemis', &
               subdir_adjemis, status, default='adjemis')
          IF_ERROR_RETURN(status=1)
          !-- make sure directory exists
          inquire(file=trim(outdir)//pathsep//trim(subdir_adjemis), exist=dir_exist)
          if( .not. dir_exist ) &
               call system('mkdir -p '//trim(outdir)//pathsep//trim(subdir_adjemis))
          ! Write the adjoint emissions to a netCDF file (TODO: write the code!)
          write(filename, '(4a,"adjemis.", a, ".", i4.4, 2i2.2, ".nc")') &
               trim(outdir), pathsep, trim(subdir_adjemis), pathsep,&
               trim(region_name(ireg)), tr(2)%year, tr(2)%month, tr(2)%day
          
          fid = nc_open(trim(filename), 'c', status)
          IF_NOTOK_RETURN(status=1)

          ! Create dimensions
          call nc_create_dim(fid, "lat", lli(ireg)%nlat)
          !call nc_create_dim(fid, "lat", jer(ireg) - jsr(ireg) + 1)
          call nc_create_dim(fid, 'lon', lli(ireg)%nlon)
          !call nc_create_dim(fid, 'lon', ier(ireg) - isr(ireg) + 1)
          call nc_create_dim(fid, 'tracer', ntracet)
          call nc_create_dim(fid, 'time', 1)

          ! Create coordinates
          ! print*, ireg, lli(ireg)%lat_deg
          call nc_dump_var(fid, 'lat', ['lat'], lli(ireg)%lat_deg, ['units    '], ['degrees_north'])
          !call nc_dump_var(fid, 'lat', ['latitude'], lli(ireg)%lat_deg(jsr(ireg):jer(ireg)), ['units    '], ['degrees_north'])
          call nc_dump_var(fid, 'lon', ['lon'], lli(ireg)%lon_deg, ['units    '], ['degrees_east'])
          !call nc_dump_var(fid, 'lon', ['longitude'], lli(ireg)%lon_deg(isr(ireg):ier(ireg)), ['units    '], ['degrees_east'])
          call nc_dump_var(fid, 'time', ['time'], [nsec])

          ! Create and the adjoint emission field
          !call nc_dump_var(fid, 'adj_emis', ['tracer   ', 'latitude ', 'longitude'], adj_emis(ireg)%values)
          call nc_dump_var(fid, 'adj_emis', ['lon      ', 'lat      ', 'tracer   '], adj_emis(ireg)%values)
          ! print*, shape(adj_emis(ireg)%values)
          ! print*, lli(ireg)%nlon, lli(ireg)%nlat

          call nc_close(fid)

          !-- MVO-COMMENT:maybe better to move this into the calling context?
          !               (to not have side-effect in this I/O routine)
          ! Reset the adjoint fields
          adj_emis(ireg)%values = 0

        end subroutine write_adj_emis


end module emission_adj
