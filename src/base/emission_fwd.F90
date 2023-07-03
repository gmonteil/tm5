!###############################################################################
!
!  Emisisons for forward run.
!
!   Note that this file contains a dummy copy, and does nothing tracer specific. In the multi-tracer code,
!   there are tracer-generic and tracer-specific things to be done by emission_adj. In this file, which is
!   a copy of the file from the project 'tracers/COCO2', all the tracer-specific parts have been commented
!   out. If you want to write a multi-tracer code, you need to uncomment and modify those parts.
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

        ! Include tracer-specific modules if needed, such as
        ! use emission_fwd_SF6, only : init_sf6_emis => emission_fwd_emis
        use GO,                 only : GO_Timer_Def
        use dims,               only : adv_scheme

        integer, intent(out)        :: status
        character(len=*), parameter :: rname = mname//'/Emission_Fwd_Init'

        ! I don't know why this check has to be in the emission code, but it's been here for historical reasons...
        select case ( trim(adv_scheme) )
        case ( 'slope' )
#ifndef slopes
            write (gol,'("adv_scheme `",a,"` while macro slopes not defined")') trim(adv_scheme); call goErr
            TRACEBACK; stop
#endif
        case default
            write (gol,'("adv_scheme `",a,"` not supported")') trim(adv_scheme); call goErr
            TRACEBACK; stop
        end select

        ! define timers:
        call GO_Timer_Def( itim_emission, 'emission', status )
        IF_NOTOK_RETURN(status=1)

        ! call tracer-specific routines, such as
        ! call init_sf6_emis(status)
        ! IF_NOTOK_RETURN(status=1)

    end subroutine Emission_Fwd_Init

    subroutine Emission_Fwd_After_Read(status)

        ! include tracer-specific modules, such as
        ! use emission_fwd_sf6, only : fwd_after_read_sf6 => emission_fwd_after_read

        character(len=*), parameter :: rname = mname//'/Emission_Fwd_After_Read'
        integer, intent(out)        :: status

        ! call tracer-specific routines, such as
        ! call fwd_after_read_sf6(status)
        ! IF_NOTOK_RETURN(status=1)

        status = 0

    end subroutine Emission_Fwd_After_Read

    subroutine Emission_Fwd_Done(status)

        ! import tracer-specific modules, such as
        ! use emission_fwd_sf6, only : done_emis_sf6 => emission_fwd_done

        character(len=*), parameter :: rname = mname//'/Emission_Fwd_Done'
        integer, intent(out)        :: status

        ! call tracer-specific routines, such as
        ! call done_emis_sf6(status)
        ! IF_NOTOK_RETURN(status=1)

        status = 0

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
        use chem_param             , only : emis_unit, emis_unit_name, tracer_name_len
        use os_specs               , only : MAX_FILENAME_LEN, SHORT_STR_LEN, DUMMY_STR_LEN

        ! --- in/out ----------------------------------------------

        integer, intent(out)             :: status

        ! --- const ------------------------------

        character(len=*), parameter      :: rname = mname//'/Emission_Fwd_Setup'

        !--- local ------------------------------------------------

        integer                         :: i_period, n_day, i_cat, n_cat, month, year
        integer                         :: region, i, j, itrac
        real                            :: tot
        integer,dimension(:,:), pointer :: zoomed
        character(len=tracer_name_len)  :: tracer_name
        character(len=MAX_FILENAME_LEN) :: emis_indir ! where are the daily cycle files?
        logical                         :: apply_dailycycle, dcycle_cat ! whether or not to apply the diurnal cycle
        character(len=SHORT_STR_LEN)    :: dailycycle_pfx ! the prefix for the daily cycle file
        integer                         :: dailycycle_type ! 0 for scaling, 1 for adding daily cycles over smaller time steps
        real                            :: emis_scale, secs_per_period, total_seconds
        character(len=2)                :: emis_scale_name
        type(TDate)                     :: t1, t2
        character(len=DUMMY_STR_LEN)    :: cat_name
        logical                         :: source_apply

        !--- begin ------------------------------------------------

        !
        ! Read from file
        !
        call ReadRc(rcf, 'proces.source', source_apply, status)
        if (.not. source_apply) return

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

            ! apply_dailycycle now says whether at least one of the categories of this tracer has a daily cycle,
            ! and if it does, that means that the tracer has a daily cycle
            do region = 1, nregions
                tracers_em_info(region)%tracer(itrac)%dailycycle%apply = apply_dailycycle
            end do

            if (apply_dailycycle) then
                call ReadRc(rcf, trim(tracer_name)//'.dailycycle.prefix', dailycycle_pfx, status)
                IF_NOTOK_RETURN(status=1)
                call ReadRc(rcF, 'dailycycle.folder', emis_indir, status)
                IF_NOTOK_RETURN(status=1)
                call ReadRc(rcf, trim(tracer_name)//'.dailycycle.type', dailycycle_type, status, default=1) ! 0 for scaling, 1 for adding
                IF_ERROR_RETURN(status=1)
                do region = 1, nregions
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
                        write(gol,'("WARNING - found zero emissions for tracer ", a, " in region ", i1, " for category ", i2)') trim(tracers_em_info(1)%tracer(itrac)%name), region, i_cat; call goPr
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

        ! import tracer-specific routines, such as
        ! use emission_fwd_sf6, only : apply_emis_sf6 => emission_fwd_apply
        use GO,                 only : GO_Timer_Start, GO_Timer_End, TDate
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

        ! At this point, we can call tracer-specific routines, such as
        ! call apply_emis_sf6(region, tr, status)
        ! IF_NOTOK_RETURN(status=1)
        !
        ! That is needed, e.g., when there is something tracer-specific about the emission code, such as calculating the high frequency ocean
        ! flux of CO2 from online winds, or injecting stuff high up such as for biomass burning CO. However, for a generic passive tracer
        ! with only surface fluxes, there is no need for a tracer-specific routine. So this module includes a generic, tracer-agnostic routine
        ! for injecting the emissions. If you call tracer-specific routines, either comment out the call to Apply_Emis_generic, or take care
        ! not to replicate that function in your tracer-specific routine.

        call Apply_Emis_generic(region, tr, status)
        IF_NOTOK_RETURN(status=1)

        call GO_Timer_End( itim_emission, status )
        IF_NOTOK_RETURN(status=1)

    end subroutine Emission_Fwd_Apply

    subroutine Apply_Emis_generic(region, tr, status)

        ! --- modules ------------------------------

        use GO,                     only : GO_Timer_Start, GO_Timer_End
        use GO,                     only : NewDate, Time_Profile_Index, TDate
        use GO,                     only : operator(+), operator(-), operator(/), rTotal
        use dims,                   only : tref, newsrun
        use dims,                   only : nregions, sec_day, ndyn
        use dims,                   only : im, jm, isr, ier, jsr, jer, lm
        use dims,                   only : itau, itaur, itaui
        use global_data,            only : mass_dat, region_dat
        use datetime,               only : tau2date, get_num_days
        use Var4D_Data,             only : calc_NL
        use Emission_Data,          only : tracers_em_info, ntracet
        use Emission_Data,          only : ref_emissions_apri
        use Emission_Read_PyShell,  only : read_dailycycle, close_dailycycle
        use chem_param,             only : ntracet
#ifdef with_budgets
    use budget_global,          only : budget_time_profile
    use budget_global,          only : budg_dat, nzon_vg, apply_budget_global, budemig, sum_emission
#endif
        ! variables for writing the 1x1 flux
        use user_output_flux1x1,    only : flux1x1_3d, write_flux1x1, calculate_flux1x1_indices, grid_translate_1x1

        ! --- in/out ----------------------------------------------

        integer, intent(in)              :: region
        type(TDate), intent(in)          :: tr(2)
        integer, intent(out)             :: status

        ! --- const -----------------------------------------------

        character(len=*), parameter      :: rname = mname//'/Apply_Emis_generic'

        !--- local ------------------------------------------------

        integer             ::  i, j, l, itr
        integer, pointer    ::  zoomed(:,:)
        real, pointer       ::  rm(:,:,:,:)
        real, pointer       ::  rzm(:,:,:,:)
        real, pointer       ::  em(:,:,:,:)
        real, allocatable   ::  ref_em_apri(:,:,:), emfac(:,:)
        real                ::  x, dtime
        integer             ::  idater(6), idate_mid(6)
        integer             ::  i_period, i_cat
        integer             ::  i_period_budget, i_day, i_period_dcycle
        integer             ::  nzone, nzone_v
        logical             ::  apply_dailycycle
        integer             ::  dailycycle_type
        ! some variables for writing the 1x1 flux
        real                :: weight, f_1x1
        integer             :: flux1x1_tidx, i_1x1, j_1x1, il

        !--- begin ------------------------------------------------

        ! short:
        zoomed => region_dat(region)%zoomed
        rm     => mass_dat(region)%rm_t
        rzm    => mass_dat(region)%rzm_t

        ! storage:
        allocate( emfac(im(region), jm(region)) )

        ! timestep emissions
        dtime = abs(rTotal( tr(2) - tr(1), 'sec' ))

        ! convert time to 6 valued array:
        call tau2date( itau, idater )

        ! get period index for budgets:
#ifdef with_budgets
    if ( apply_budget_global ) then
      call Time_Profile_Index( budget_time_profile, NewDate(time6=idater), i_period_budget, status , at_left_side = .true. )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

        ! Calculate some variables necessary for the 1x1 flux
        if (write_flux1x1) then
            call calculate_flux1x1_indices(tr, flux1x1_tidx, weight)
        end if

        do itr = 1, ntracet

            ! check whether to read new daily cycle file
            apply_dailycycle = tracers_em_info(region)%tracer(itr)%dailycycle%apply ! This T/F will come from an rc file, read in by emission_fwd/emission_fwd_setup

            if (apply_dailycycle) then
                dailycycle_type = tracers_em_info(region)%tracer(itr)%dailycycle%dtype
                ! find out the day number
                i_day = get_num_days(itaui, itaur(region) + ndyn/4/tref(region))
                call tau2date(itaur(region) + ndyn/4/tref(region), idate_mid)
                if ( .not. tracers_em_info(region)%tracer(itr)%dailycycle%day_opened(i_day) ) then
                    ! close_dailycycle only needs the date to write the emissions applied, in case verbose_debug_output is defined
                    ! in that case, we should pass a date that's a little before itaur(region), otherwise the emission applied on day 1
                    ! gets written to a file with a filename corresponding to day 2
                    if ( .not. newsrun ) call close_dailycycle(region, itr)
                    call read_dailycycle(i_day, idate_mid, region, itr, status)
                    IF_NOTOK_RETURN(status=1)
                end if
            end if ! apply_dailycycle

            ! ~ loop over categories:
            do i_cat = 1, tracers_em_info(region)%tracer(itr)%n_cat

                ! point to emission array:
                em => ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field   !  kg/cl/s or factor around 0.0 or 1.0

                ! calculate the time step index for daily cycle, if it is to be applied
                apply_dailycycle = tracers_em_info(region)%tracer(itr)%dailycycle%cycle_cat(i_cat)%apply ! this is true only if the daily cycle file read in previously has a cycle for this category
                if (apply_dailycycle) &
                        i_period_dcycle = (idate_mid(4) * 3600 + idate_mid(5) * 60 + idate_mid(6)) / tracers_em_info(region)%tracer(itr)%dailycycle%cycle_cat(i_cat)%dtime + 1
                ! this is the time index within one day

                ! get period index for emissions:
                call Time_Profile_Index( tracers_em_info(region)%tracer(itr)%cat(i_cat)%time_profile, NewDate(time6=idater), i_period, status, at_left_side=.true. )
                IF_NOTOK_RETURN(status=1)

                ! em contains real emissions in kg/cl/s
                !   c(t+dt) = c(t) + em dt
                emfac = em(:,:,1,i_period)   ! kg/cl/s

                !
                ! Emission operator:
                !
                !  c(t+dt) = c(t) + emfac(t) dt
                !          ~ c(t) + em ew(t) dt     ! linearized
                !
                ! where:
                !   c   : tracer mass [kg cl-1]
                !   em  : emission field [kg cl-1 s-1]
                !   ew  : emission weight for a specific time [1]
                !   dt  : time step (s)
                !
                ! In matrix form:
                !
                !  [ c  ]   [ I  ew*dt ] [ c  ]
                !  [    ] = [          ] [    ]
                !  [ em ]   [ O    I   ] [ em ]
                !
                ! Adjoint:
                !
                !  [ c*  ]   [   I   O ] [ c*  ]
                !  [     ] = [         ] [     ]
                !  [ em* ]   [ ew*dt I ] [ em* ]
                !
                ! where:
                !   c*   : adjoint concentration 1/[kg cl-1]
                !   em*  : adjoint emission field 1/[kg cl-1 s-1]
                !
                ! which gives the adjoint operator:
                !
                !   em* =  c* ew dt  +  em*
                !
                do j = jsr(region), jer(region)
                    do i = isr(region), ier(region)

                        if ( zoomed(i,j) /= region ) cycle

                        if (apply_dailycycle) then
                            select case (dailycycle_type)
                            case (0)
                                x = emfac(i,j) * dtime * &
                                        tracers_em_info(region)%tracer(itr)%dailycycle%cycle_cat(i_cat)%scaling(i,j,i_period_dcycle)
                            case (1)
                                x = (emfac(i,j) + &
                                        tracers_em_info(region)%tracer(itr)%dailycycle%cycle_cat(i_cat)%anomaly(i,j,i_period_dcycle)) * dtime
                            end select
                        else
                            x = emfac(i,j) * dtime   ! kg/cl/timestep
                        end if ! apply_dailycycle

                        rm(i,j,1,itr)  = rm(i,j,1,itr)  + x
#ifdef slopes
                    rzm(i,j,1,itr) = rzm(i,j,1,itr) - x
#endif
                        ! 1x1 flux
                        if (write_flux1x1) then
                            do il = 1, grid_translate_1x1(region)%cell(i,j)%N
                                i_1x1 = grid_translate_1x1(region)%cell(i,j)%ilist(il)
                                j_1x1 = grid_translate_1x1(region)%cell(i,j)%jlist(il)
                                f_1x1 = grid_translate_1x1(region)%cell(i,j)%frac(il)
                                flux1x1_3d%prod(i_1x1,j_1x1,1,flux1x1_tidx,itr) = flux1x1_3d%prod(i_1x1,j_1x1,1,flux1x1_tidx,itr) + &
                                        f_1x1 * x
                            end do
                        end if ! write_flux1x1

                        !=========
                        ! budget
                        !=========
#ifdef with_budgets
                    if ( apply_budget_global ) then

                        !global budget
                        nzone = budg_dat(region)%nzong(i,j)
                        nzone_v = nzon_vg(1)

                        budemig(nzone,nzone_v,itr,i_period_budget) = budemig(nzone,nzone_v,itr,i_period_budget) + x ! [kg]

                        sum_emission(region,itr) = sum_emission(region,itr) + x  ! [kg]
                    end if
#endif
                    end do ! i
                end do ! j
            end do ! i_cat
        end do ! itr

        ! unlink:
        nullify( zoomed )
        nullify( rm )
        nullify( rzm )
        nullify( em )

        ! clear:
        deallocate( emfac )

        ! ok:
        status = 0

    end subroutine Apply_Emis_generic

end module Emission_Fwd
