
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"

module emission_fwd

    use chem_param,     only : ntracet, tracer_name => names, tracer_name_len
    use dims,           only : nregions, region_name, itaui, itaue, itaur
    use os_specs,       only : DUMMY_STR_LEN, MAX_FILENAME_LEN
    use go,             only : readrc, tdate, goerr, gol, gopr
    use global_data,    only : rcf
    use emission_data,  only : tracers_em_info, ref_emissions_apri, source_apply, t_tracer_info
    use datetime,       only : time_window, get_num_days, tau2date
    use production,     only : production_fwd

    implicit none

    private

    public :: emission_fwd_init             ! ==> called by emission.F90
    public :: emission_fwd_setup            ! ==> called by modelIntegration.F90
    public :: emission_fwd_apply            ! ==> called by emission.F90
    public :: emission_fwd_done             ! ==> called by initexit.F90
    public :: emission_fwd_after_read       ! ==> called by sources_sinks.F90

    integer :: dailycycle_type

    character(len=*), parameter   :: mname = 'Emission_Fwd'

contains

    subroutine emission_fwd_init(status)
        ! Call chain:
        ! main -> tm5var4d/tm5var4d_Run -> tracer/tracer_model -> modelintegration/proces_init -> emission/emission_init -> emission_fwd_init
        ! ==> this is more or less identical to emission_fwd_setup, so putting everything there.
        integer, intent(out) :: status
        status = 0
    end subroutine emission_fwd_init


    subroutine emission_fwd_after_read(status)
        ! This routine is called from trace_after_read, which is called by both tracer and adj_tracer.
        ! Thus, there is no separate Emission_Adj_After_Read.
        integer, intent(out) :: status
        status = 0
    end subroutine emission_fwd_after_read


    subroutine emission_fwd_done(status)
        integer, intent(out) :: status
        status = 0
    end subroutine emission_fwd_done


    subroutine emission_fwd_setup(status)
        ! Setup the emissions:
        ! 1. Read the coarse-resolution emissions (and store them in memory)
        ! 2. Determine whether dailycycle are to be applied, for each tracer/category
        !
        ! Call chain: main -> tm5var4d/tm5var4d_Run -> tracer/tracer_model -> initexit/start_tm5 -> initexit/start_tm5_emissions -> emission_fwd_setup

        use emission_read_pyshell,  only : read_emissions_from_pyshell
        use dims,                   only : adv_scheme

        integer, intent(out)            :: status
        logical                         :: apply_dailycycle, dcycle_cat
        integer                         :: itrac, ireg, icat, ncat, nday, iday
        character(len=DUMMY_STR_LEN)    :: cat
        character(len=tracer_name_len)  :: tracer
        character(len=256)              :: emis_indir
        character(len=256)              :: dailycycle_pfx

        character(len=*), parameter     :: rname = mname//'/emission_fwd_setup'

        status = 0

        if (.not. source_apply) return

        ! Check that we are using the slope scheme
        if (trim(adv_scheme) /= 'slope') then
            write (gol,'("adv_scheme `",a,"` not supported")') trim(adv_scheme); call goErr
            status = 1
        end if
#ifndef slopes
        write (gol,'("adv_scheme "slopes" while macro slopes not defined")'); call goErr
        status = 1
#endif

        !call read_emissions_from_pyshell(tracers_em_info, time_window, ref_emissions_apri, status)
        !IF_NOTOK_RETURN(status=1)

    end subroutine emission_fwd_setup


    subroutine emission_fwd_apply(ireg, tr, status)
        ! Application of the emissions:
        ! 1. (optional) read the dailycycle files
        ! 2. apply the surface emissions
        ! 3. (optional) apply 3D emission fields (i.e. "production")
        !
        ! Call chain:
        ! main -> tm5var4d/tm5var4d_Run -> tracer/tracer_model -> modelIntegration/proces_region -> modelIntegration/do_steps -> source_sinks/source1 -> emission_fwd_apply
        integer, intent(in)     :: ireg
        type(tdate), intent(in) :: tr(2)
        integer, intent(out)    :: status
        integer :: itrac
        character(len=*), parameter :: rname = mname//'/emission_fwd_apply'

        do itrac = 1, ntracet
            call emission_fwd_apply_tracer(ireg, itrac, tracers_em_info(ireg)%tracer(itrac), tr, status)
            !call production_fwd(itrac, ireg, tr(1), tr(2))
        end do

    end subroutine emission_fwd_apply


    subroutine emission_fwd_apply_tracer(ireg, itrac, tracer, tr, status)

        use go,                     only : rtotal, newdate, operator(-), time_profile_index
        use dims,                   only : ndyn, tref, newsrun, im, jm, isr, jsr, ier, jer, itau
        use emission_read_pyshell,  only : close_dailycycle, read_dailycycle
        use global_data,            only : region_dat, mass_dat

        integer, intent(in)                 :: ireg, itrac
        type(T_tracer_info), intent(in)     :: tracer
        type(tdate), intent(in)             :: tr(2)
        integer, intent(out)                :: status

        logical :: apply_dailycycle
        integer :: iday, icat, iperiod_dcycle, iperiod_emis
        integer, dimension(6) :: idate_mid, idater

        real, dimension(:, :, :, :), pointer    :: em
        real, dimension(:, :, :), allocatable   :: emis

        real    :: dtime
        real    :: x
        integer :: i, j

        character(len=*), parameter :: rname = mname//'/emission_fwd_apply_tracer'

        ! find out the day number
        iday = get_num_days(itaui, itaur(ireg) + ndyn / 4 / tref(ireg))
        call tau2date(itaur(ireg) + ndyn / 4 / tref(ireg), idate_mid)

        emis = read_emis(idate_mid, ireg, tracer, status)

        ! timestep emissions
        dtime = abs(rtotal(tr(2) - tr(1), 'sec'))

        ! Loop over categories :
        do icat = 1, tracer%n_cat
            do j = jsr(ireg), jer(ireg)
                do i = isr(ireg), ier(ireg)
                    if (region_dat(ireg)%zoomed(i, j) /= ireg) cycle
                    x = emis(icat, i, j) * dtime
                    mass_dat(ireg)%rm_t(i, j, 1, itrac) = mass_dat(ireg)%rm_t(i, j, 1, itrac) + x
                    mass_dat(ireg)%rzm_t(i, j, 1, itrac) = mass_dat(ireg)%rzm_t(i, j, 1, itrac) - x
                end do
            end do
        end do

    end subroutine emission_fwd_apply_tracer


    function read_emis(idate_local, ireg, tracer, status) result(emis)
        use dims, only : im, jm
        use file_netcdf

        ! This is inefficient as the netCDF files get read multiple times, but this is just a proof of concept

        !__IO___________________________________________________________________
        integer, dimension(6), intent(in)       :: idate_local
        integer, intent(in)                     :: ireg
        type(t_tracer_info), intent(in)         :: tracer
        integer, intent(out)                    :: status
        real, dimension(:, :, :), allocatable   :: emis

        !__LOCAL_VARIABLES______________________________________________________
        character(len=MAX_FILENAME_LEN)     :: fname
        logical                             :: file_exists
        integer                             :: fid
        integer                             :: i_cat
        real, dimension(:, :), allocatable  :: tmp

        character(len=*), parameter     :: rname = mname//'/read_emis'

        ! Construct the file name
        ! should be something like {path.emissions}/tm5emis.{tracer}.{region}.{date}.nc

        fname = get_emis_filename(idate_local, ireg, trim(tracer%name))
        fid = nc_open(trim(fname), 'r', status)
        IF_NOTOK_RETURN(status=1)

        allocate(emis(tracer%n_cat, im(ireg), jm(ireg)))

        emis(:, :, :) = 0

        do i_cat = 1, tracer%n_cat
            tmp = nc_read_var(fid, tracer%cat(i_cat)%name, status)
            emis(i_cat, :, :) = tmp
            deallocate(tmp)
        end do

        call nc_close(fid)

    end function read_emis


    function get_emis_filename(idate, ireg, trname) result(filename)
        use dims, only : region_name

        integer, dimension(6), intent(in)   :: idate
        integer, intent(in)                 :: ireg
        character(len=*), intent(in)        :: trname
        character(len=MAX_FILENAME_LEN)     :: filename

        character(len=MAX_FILENAME_LEN)     :: prefix
        integer                             :: status

        call readrc(rcf, 'emissions.'// trname // '.prefix', prefix, status)
        write(filename, '(a, ".", a, ".", a, ".", i4.4, 2i2.2, ".nc")') trim(prefix), trname, trim(region_name(ireg)), idate(1), idate(2), idate(3)

    end function get_emis_filename

end module emission_fwd
