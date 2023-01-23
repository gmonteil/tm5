!###############################################################################
!
!  Emisisons for forward run.
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if

#include "tm5.inc"

module Emission_Fwd

    use GO, only : gol, goErr, goPr

    implicit none

    private

    public :: Emission_Fwd_Init
    public :: Emission_Fwd_Setup
    public :: Emission_Fwd_Apply
    public :: Emission_Fwd_Done
    public :: Emission_Fwd_After_Read

    character(len=*), parameter   :: mname = 'Emission_Fwd'

    ! timers:
    integer ::  itim_emission

contains

    subroutine Emission_Fwd_Init( status )

        use Emission_Fwd_CO2,   only : Init_CO2_emis => Emission_Fwd_Init
        use GO,                 only : GO_Timer_Def

        integer, intent(out)        :: status
        character(len=*), parameter :: rname = mname//'/Emission_Fwd_Init'

        ! define timers:
        call GO_Timer_Def( itim_emission, 'emission', status )
        IF_NOTOK_RETURN(status=1)

        call Init_CO2_emis(status)
        IF_NOTOK_RETURN(status=1)

    end subroutine Emission_Fwd_Init

    subroutine Emission_Fwd_After_Read(status)

        use Emission_Fwd_CO2,   only : fwd_after_read_CO2 => Emission_Fwd_After_Read

        character(len=*), parameter :: rname = mname//'/Emission_Fwd_After_Read'
        integer, intent(out)        :: status

        call fwd_after_read_CO2(status)
        IF_NOTOK_RETURN(status=1)

    end subroutine Emission_Fwd_After_Read

    subroutine Emission_Fwd_Done(status)

        use Emission_Fwd_CO2,   only : Done_emis_CO2 => Emission_Fwd_Done

        character(len=*), parameter :: rname = mname//'/Emission_Fwd_Done'
        integer, intent(out)        :: status

        call Done_emis_CO2(status)
        IF_NOTOK_RETURN(status=1)

    end subroutine Emission_Fwd_Done

    subroutine Emission_Fwd_Setup( status )

        ! --- modules ------------------------------

        use GO                     , only : Readrc, TDate, operator(-), rTotal
        use dims                   , only : nregions, region_name
        use dims                   , only : isr, ier, jsr, jer
        use dims                   , only : sec_normal_year, itaui, itaue
        use global_data            , only : rcF
        use global_data            , only : region_dat
        use datetime               , only : time_window, get_num_days
        use Emission_Data          , only : tracers_em_info, ntracet
        use Emission_Data          , only : ref_emissions_apri
        use Emission_Read_PyShell  , only : Read_Emissions_From_PyShell
        use chem_param             , only : emis_unit, emis_unit_name
        use os_specs               , only : DUMMY_STR_LEN

        ! --- in/out ----------------------------------------------

        integer, intent(out)             :: status

        ! --- const ------------------------------

        character(len=*), parameter      :: rname = mname//'/Emission_Fwd_Setup'

        !--- local ------------------------------------------------

        integer                         :: i_period, n_day, i_cat, n_cat, month, year
        integer                         :: region, i, j, itrac
        real                            :: tot
        integer,dimension(:,:), pointer :: zoomed
        character(len=8)                :: tracer_name
        character(len=256)              :: emis_indir ! where are the daily cycle files?
        logical                         :: apply_dailycycle, dcycle_cat ! whether or not to apply the diurnal cycle
        character(len=256)              :: dailycycle_pfx ! the prefix for the daily cycle file
        integer                         :: dailycycle_type ! 0 for scaling, 1 for adding daily cycles over smaller time steps
        real                            :: emis_scale, secs_per_period, total_seconds
        character(len=2)                :: emis_scale_name
        type(TDate)                     :: t1, t2
        character(len=DUMMY_STR_LEN)    :: cat_name

        !--- begin ------------------------------------------------

        !
        ! Read from file
        !

        call Read_Emissions_From_PyShell ( tracers_em_info, time_window, ref_emissions_apri, status )
        IF_NOTOK_RETURN(status=1)

        ! Check whether we want to apply daily cycle or not, per tracer and per region
        do itrac = 1, ntracet
            tracer_name = trim(tracers_em_info(1)%tracer(itrac)%name)
            ! Figure out if this tracer needs sub-period chunking for any of its categories
            apply_dailycycle = .false.
            do region = 1, nregions
                n_cat = tracers_em_info(region)%tracer(itrac)%n_cat
                do i_cat = 1, n_cat
                    cat_name = tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name
                    call ReadRc(rcf, trim(tracer_name)//'.'//trim(cat_name)//'.dailycycle', dcycle_cat, status, default=.false.)
                    IF_ERROR_RETURN(status=1)
                    if (dcycle_cat) write(*,'("Tracer ",a," category ",a," has sub-period chunking")') trim(tracer_name), trim(cat_name)
                    apply_dailycycle = apply_dailycycle .or. dcycle_cat
                end do ! i_cat
            end do ! region
            !call ReadRc(rcf, trim(tracer_name)//'.emission.dailycycle', apply_dailycycle, status)
            !IF_NOTOK_RETURN(status=1)
            if (apply_dailycycle) then
                call ReadRc(rcf, trim(tracer_name)//'.dailycycle.prefix', dailycycle_pfx, status)
                IF_NOTOK_RETURN(status=1)
                call ReadRc(rcF, 'dailycycle.folder', emis_indir, status)
                IF_NOTOK_RETURN(status=1)
                call ReadRc(rcf, trim(tracer_name)//'.dailycycle.type', dailycycle_type, status) ! 0 for scaling, 1 for adding
                IF_NOTOK_RETURN(status=1)
                do region = 1, nregions
                    tracers_em_info(region)%tracer(itrac)%dailycycle%apply = apply_dailycycle
                    n_cat = tracers_em_info(region)%tracer(itrac)%n_cat
                    allocate(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(n_cat))
                    ! The apply component will be set by read_dailycycle
                    tracers_em_info(region)%tracer(itrac)%dailycycle%pfx = dailycycle_pfx
                    tracers_em_info(region)%tracer(itrac)%dailycycle%emis_indir = emis_indir
                    tracers_em_info(region)%tracer(itrac)%dailycycle%dtype = dailycycle_type
                    if (.not. allocated(tracers_em_info(region)%tracer(itrac)%dailycycle%day_opened)) then
                        ! How many days during the model run?
                        n_day = get_num_days(itaui, itaue)
                        allocate(tracers_em_info(region)%tracer(itrac)%dailycycle%day_opened(n_day))
                        do i=1,n_day
                            tracers_em_info(region)%tracer(itrac)%dailycycle%day_opened(i) = .false.
                        end do ! n_day
                    end if ! allocated(day_opened)
                end do ! region
            else
                ! emission_fwd_apply checks per category if there is a daily cycle. If the tracer does not have
                ! a daily cycle, this check is redundant, but triggers a segfault if dailycycle%cycle_cat is
                ! not allocated. We could catch that with an if statement, but I'd rather avoid extra if statements
                ! when possible. So here we allocate dailycycle%cycle_cat and set the apply components to false.
                do region = 1, nregions
                    n_cat = tracers_em_info(region)%tracer(itrac)%n_cat
                    allocate(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(n_cat))
                    tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(:)%apply = .false.
                end do ! region
            end if ! apply_dailycycle
        end do ! itrac

        ! Check whether whole emission array has been filled
        do itrac = 1, ntracet
            do region = 1, nregions
                do i_cat = 1, tracers_em_info(region)%tracer(itrac)%n_cat
                    !print *, 'checking region:', region, 'i_cat', i_cat, 'nperiod:',&
                    !tracers_em_info(region)%tracer(itrac)%cat(i_cat)%time_profile%n_period
                    tot = 0.0
                    do i_period = 1, tracers_em_info(region)%tracer(itrac)%cat(i_cat)%time_profile%n_period
                        tot = tot + sum( ref_emissions_apri(region)%tracer(itrac)%cat(i_cat)%field(:,:,:,i_period) )
                    end do
                    !print *, "total over periods, for tracer ",tracers_em_info(region)%tracer(itrac)%name,":", tot
                    if ( tot == 0.0 ) then
                        write (gol,'("WARNING - found zero emissions in region for category ",i2,", region", i3)') i_cat, region; call goErr
                        TRACEBACK;
                    end if
                end do
            end do
        end do

        !
        ! Print total emissions for the entire simulation period
        !

        write(*,'(a8,3x,a10,3x,a20,3x,a12,3x,a20)') 'Tracer  ', 'Region    ', 'Category            ', 'Total/yr', 'Total (Kg)'
        write(*,'(a)') repeat('=', 82)

        do itrac = 1, ntracet
            emis_scale = emis_unit(itrac)
            emis_scale_name = emis_unit_name(itrac)
            do region = 1, nregions
                zoomed => region_dat(region)%zoomed
                do i_cat = 1, tracers_em_info(region)%tracer(itrac)%n_cat
                    tot = 0.0
                    total_seconds = 0.0
                    do i_period = 1, tracers_em_info(region)%tracer(itrac)%cat(i_cat)%time_profile%n_period
                        t1 = tracers_em_info(region)%tracer(itrac)%cat(i_cat)%time_profile%period(i_period)%t1
                        t2 = tracers_em_info(region)%tracer(itrac)%cat(i_cat)%time_profile%period(i_period)%t2
                        secs_per_period = rTotal(t2-t1, 'sec')
                        do j = jsr(region), jer(region)
                            do i = isr(region), ier(region)
                                if ( zoomed(i,j) /= region ) cycle
                                tot = tot + secs_per_period * ref_emissions_apri(region)%tracer(itrac)%cat(i_cat)%field(i,j,1,i_period)
                            end do ! i
                        end do ! j
                        total_seconds = total_seconds + secs_per_period
                    end do ! i_period
                    write(*,'(a8,3x,a10,3x,a20,3x,f9.2,1x,a2,3x,es20.13)') tracers_em_info(1)%tracer(itrac)%name, &
                            region_name(region), tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name, &
                            sec_normal_year * (tot/total_seconds) * emis_scale, emis_scale_name, tot
                end do ! i_cat
                nullify( zoomed )
            end do ! region
            write(*,'(a)') repeat('-', 82)
        end do ! itrac

        ! ok
        status = 0

    end subroutine Emission_Fwd_Setup

    subroutine Emission_Fwd_Apply(region, tr, status)

        use Emission_Fwd_CO2,       only : Apply_emis_CO2 => Emission_Fwd_Apply
        use GO,                     only : GO_Timer_Start, GO_Timer_End, TDate
        ! Need additional modules if we want to write 1x1 fluxes
        use user_output_flux1x1,    only : flux1x1_3d, write_flux1x1, calculate_flux1x1_indices, grid_translate_1x1
        use dims,                   only : nregions, isr, ier, jsr, jer
        use global_data,            only : region_dat
        use MeteoData,              only : phlb_dat, gph_dat

        character(len=*), parameter :: rname = mname//'/Emission_Fwd_Apply'

        integer, intent(out)        :: status
        type(TDate), intent(in)     :: tr(2)
        integer, intent(in)         :: region
        ! Need additional variables for writing out 1x1 fluxes
        real                        :: weight
        integer                     :: flux1x1_tidx, i_1x1, j_1x1, il, i, j
        integer, pointer            :: zoomed(:,:)
        real, pointer               :: gph(:,:,:), phlb(:,:,:)

        call GO_Timer_Start( itim_emission, status )
        IF_NOTOK_RETURN(status=1)

        if (write_flux1x1) then
            zoomed => region_dat(region)%zoomed
            gph    => gph_dat(region)%data
            phlb   => phlb_dat(region)%data

            call calculate_flux1x1_indices(tr, flux1x1_tidx, weight)
            ! The pressures and geopotential heights need not be added once per tracer/category, so do that addition here
            do j = jsr(region), jer(region)
                do i = isr(region), ier(region)
                    if ( zoomed(i,j) /= region ) cycle
                    do il = 1, grid_translate_1x1(region)%cell(i,j)%N
                        i_1x1 = grid_translate_1x1(region)%cell(i,j)%ilist(il)
                        j_1x1 = grid_translate_1x1(region)%cell(i,j)%jlist(il)

                        flux1x1_3d%p(i_1x1,j_1x1,:,flux1x1_tidx) = flux1x1_3d%p(i_1x1,j_1x1,:,flux1x1_tidx) + &
                                weight * phlb(i,j,:)
                        flux1x1_3d%gph(i_1x1,j_1x1,:,flux1x1_tidx) = flux1x1_3d%gph(i_1x1,j_1x1,:,flux1x1_tidx) + &
                                weight * gph(i,j,:)
                        flux1x1_3d%weight(i_1x1,j_1x1,flux1x1_tidx) = flux1x1_3d%weight(i_1x1,j_1x1,flux1x1_tidx) + weight
                    end do
                end do
            end do

            nullify(zoomed, gph, phlb)
        end if

        call Apply_emis_CO2(region, tr, status)
        IF_NOTOK_RETURN(status=1)

        call GO_Timer_End( itim_emission, status )
        IF_NOTOK_RETURN(status=1)

    end subroutine Emission_Fwd_Apply

end module Emission_Fwd
