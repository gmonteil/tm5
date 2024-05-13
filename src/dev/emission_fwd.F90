
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"

module emission_fwd

    use chem_param,     only : ntracet, tracer_name => names, tracer_name_len
    use dims,           only : nregions, region_name, itaui, itaue, itaur
    use os_specs,       only : DUMMY_STR_LEN
    use go,             only : readrc, tdate, goerr, gol, gopr
    use global_data,    only : rcf
    use emission_data,  only : tracers_em_info, ref_emissions_apri, source_apply
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

        call read_emissions_from_pyshell(tracers_em_info, time_window, ref_emissions_apri, status)
        IF_NOTOK_RETURN(status=1)

        ! Check whether we want to apply daily cycle or not, per tracer and per region
        do itrac = 1, ntracet
            tracer = tracer_name(itrac)
            ! Figure out if this tracer needs sub-period chunking for any of its categories
            apply_dailycycle = .false.
            do ireg = 1, nregions
                ncat = tracers_em_info(ireg)%tracer(itrac)%n_cat
                do icat = 1, ncat
                    cat = tracers_em_info(ireg)%tracer(itrac)%cat(icat)%name
                    call readrc(rcf, trim(tracer) // '.' // trim(cat) // '.dailycycle', dcycle_cat, status, default=.false.)
                    IF_NOTOK_RETURN(status=1)
                    if (dcycle_cat) write(*,'("Tracer ",a," category ",a," has sub-period chunking")') trim(tracer), trim(cat)
                    apply_dailycycle = apply_dailycycle .or. dcycle_cat
                end do
            end do

            if (apply_dailycycle) then
                call readrc(rcf, trim(tracer) // '.dailycycle.prefix', dailycycle_pfx, status)
                IF_NOTOK_RETURN(status=1)
                call readrc(rcf, 'dailycycle.folder', emis_indir, status)
                IF_NOTOK_RETURN(status=1)
                call readrc(rcf, trim(tracer) // '.dailycycle.type', dailycycle_type, status) ! 0 for scaling, 1 for adding
                IF_NOTOK_RETURN(status=1)
                do ireg = 1, nregions
                    tracers_em_info(ireg)%tracer(itrac)%dailycycle%apply = .true.
                    ncat = tracers_em_info(ireg)%tracer(itrac)%n_cat
                    allocate(tracers_em_info(ireg)%tracer(itrac)%dailycycle%cycle_cat(ncat))
                    ! The apply component will be set by read_dailycycle
                    tracers_em_info(ireg)%tracer(itrac)%dailycycle%pfx = dailycycle_pfx
                    tracers_em_info(ireg)%tracer(itrac)%dailycycle%emis_indir = emis_indir
                    tracers_em_info(ireg)%tracer(itrac)%dailycycle%dtype = dailycycle_type
                    if (.not. allocated(tracers_em_info(ireg)%tracer(itrac)%dailycycle%day_opened)) then
                        ! How many days during the model run?
                        nday = get_num_days(itaui, itaue)
                        allocate(tracers_em_info(ireg)%tracer(itrac)%dailycycle%day_opened(nday))
                        tracers_em_info(ireg)%tracer(itrac)%dailycycle%day_opened(:) = .false.
                    end if ! allocated(day_opened)
                end do
            else
                ! emission_fwd_apply checks per category if there is a daily cycle. If the tracer does not have
                ! a daily cycle, this check is redundant, but triggers a segfault if dailycycle%cycle_cat is
                ! not allocated. We could catch that with an if statement, but I'd rather avoid extra if statements
                ! when possible. So here we allocate dailycycle%cycle_cat and set the apply components to false.
                do ireg = 1, nregions
                    ncat = tracers_em_info(ireg)%tracer(itrac)%n_cat
                    allocate(tracers_em_info(ireg)%tracer(itrac)%dailycycle%cycle_cat(ncat))
                    tracers_em_info(ireg)%tracer(itrac)%dailycycle%cycle_cat(:)%apply = .false.
                end do ! region`
            end if
        end do

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

        do itrac = 1, ntracet
            call emission_fwd_apply_tracer(ireg, itrac, tracers_em_info(ireg)%tracer(itrac), tr, status)
            call production_fwd(itrac, ireg, tr(1), tr(2))
        end do

    end subroutine emission_fwd_apply


    subroutine emission_fwd_apply_tracer(ireg, itrac, tracer, tr, status)

        use go,                     only : rtotal, newdate, operator(-), time_profile_index
        use dims,                   only : ndyn, tref, newsrun, im, jm, isr, jsr, ier, jer, itau
        use emission_read_pyshell,  only : close_dailycycle, read_dailycycle
        use emission_data,          only : T_tracer_info
        use global_data,            only : region_dat, mass_dat

        integer, intent(in)                 :: ireg, itrac
        type(T_tracer_info), intent(in)     :: tracer
        type(tdate), intent(in)             :: tr(2)
        integer, intent(out)                :: status

        logical :: apply_dailycycle
        integer :: iday, icat, iperiod_dcycle, iperiod_emis
        integer, dimension(6) :: idate_mid, idater

        real, dimension(:, :, :, :), pointer    :: em

        real    :: dtime
        real    :: x
        integer :: i, j

        character(len=*), parameter :: rname = mname//'/emission_fwd_apply_tracer'

        ! Check whether to read new dailycycle file
        if (tracer%dailycycle%apply) then
            ! find out the day number
            iday = get_num_days(itaui, itaur(ireg) + ndyn / 4 / tref(ireg))
            call tau2date(itaur(ireg) + ndyn / 4 / tref(ireg), idate_mid)

            if (.not. tracer%dailycycle%day_opened(iday)) then
                ! close_dailycycle only needs the date to write the emissions applied, in case verbose_debug_output is defined
                ! in that case, we should pass a date that's a little before itaur(region), otherwise the emission applied on day 1
                ! gets written to a file with a filename corresponding to day 2
                if (.not. newsrun) call close_dailycycle(ireg, itrac)
                call read_dailycycle(iday, idate_mid, ireg, itrac, status)
            end if
        end if

        ! timestep emissions
        dtime = abs(rtotal(tr(2) - tr(1), 'sec'))

        ! convert time to 6 valued array:
        call tau2date(itau, idater)

        ! Loop over categories :
        do icat = 1, tracer%n_cat
            em => ref_emissions_apri(ireg)%tracer(itrac)%cat(icat)%field !  kg/cl/s or factor around 0.0 or 1.0

            ! calculate the time step index for daily cycle, if it is to be applied
            if (tracer%dailycycle%cycle_cat(icat)%apply) then ! this is true only if the daily cycle file read in previously has a cycle for this category
                ! time index within one day:
                iperiod_dcycle = (idate_mid(4) * 3600 + idate_mid(5) * 60 + idate_mid(6)) / tracer%dailycycle%cycle_cat(icat)%dtime + 1
            end if

            ! get period index for emissions:
            call time_profile_index(tracer%cat(icat)%time_profile, newdate(time6=idater), iperiod_emis, status, at_left_side=.true.)
            IF_NOTOK_RETURN(status=1)

            ! em contains real emissions in kg/cl/s
            !   c(t+dt) = c(t) + em dt

            do j = jsr(ireg), jer(ireg)
                do i = isr(ireg), ier(ireg)
                    if (region_dat(ireg)%zoomed(i, j) /= ireg) cycle

                    if (tracer%dailycycle%cycle_cat(icat)%apply) then
                        if (dailycycle_type == 0) then
                            x = em(i, j, 1, iperiod_emis) * dtime * tracer%dailycycle%cycle_cat(icat)%scaling(i, j, iperiod_dcycle)
                        else if (dailycycle_type == 1) then
                            x = em(i, j, 1, iperiod_emis) + dtime * tracer%dailycycle%cycle_cat(icat)%anomaly(i, j, iperiod_dcycle)
                        end if
                    else
                        x = em(i, j, 1, iperiod_emis) * dtime
                    end if
                    mass_dat(ireg)%rm_t(i, j, 1, itrac) = mass_dat(ireg)%rm_t(i, j, 1, itrac) + x
                    mass_dat(ireg)%rzm_t(i, j, 1, itrac) = mass_dat(ireg)%rzm_t(i, j, 1, itrac) - x
                end do
            end do
        end do


    end subroutine emission_fwd_apply_tracer

end module emission_fwd