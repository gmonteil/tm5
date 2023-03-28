!###############################################################################
!
!  Emisisons for forward run.
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"

module Emission_Fwd_CO2

    use GO, only : gol, goErr, goPr

    use dims,          only : nregions
    use TM5_Fields,    only : T_Fields_4D

    implicit none

    private

    public :: Emission_Fwd_Init
    public :: Emission_Fwd_Apply
    public :: Emission_Fwd_Done
    public :: Emission_Fwd_After_Read

    ! --- const ------------------------------

    character(len=*), parameter        :: mname = 'Emission_Fwd_CO2'

    ! --- var ------------------------------

    ! some total number:
    real                                  ::  emissions_total

    contains

        ! allocate arrays for forward run
        subroutine Emission_Fwd_Init( status )
            integer, intent(out)             :: status
            character(len=*), parameter      :: rname = mname//'/Emission_Fwd_Init'
            status = 0
        end subroutine Emission_Fwd_Init

        subroutine Emission_Fwd_After_Read ( status )
            integer, intent(out)             :: status
            character(len=*), parameter      :: rname = mname//'/Emission_Fwd_After_Read'
            status = 0
        end subroutine Emission_Fwd_After_Read

        subroutine Emission_Fwd_Done( status )
            integer, intent(out)             :: status
            character(len=*), parameter      :: rname = mname//'/Emission_Fwd_Done'
            status = 0
        end subroutine Emission_Fwd_Done

        !---------------------------------------------------------
        ! o Apply emissions
        ! o Update budgets
        !---------------------------------------------------------
        subroutine Emission_Fwd_Apply( region, tr, itracer, status )
            use GO,                     only : GO_Timer_Start, GO_Timer_End
            use GO,                     only : NewDate, Time_Profile_Index, TDate
            use GO,                     only : operator(+), operator(-), operator(/), rTotal
            use dims,                   only : tref, newsrun
            use dims,                   only : nregions, sec_day, ndyn
            use dims,                   only : im, jm, isr, ier, jsr, jer, lm
            use dims,                   only : itau, itaur, itaui
            use dims,                   only : adv_scheme
            use global_data,            only : mass_dat, region_dat !, meteo_dat, conv_dat
            use datetime,               only : tau2date, get_num_days
            use Var4D_Data,             only : calc_NL
            use Emission_Data,          only : tracers_em_info
            use Emission_Data,          only : ref_emissions_apri
            use Emission_Read_PyShell,  only : read_dailycycle, close_dailycycle
        !    use chem_param,             only : ico2
#ifdef with_budgets
            use budget_global,          only : budget_time_profile
            use budget_global,          only : budg_dat, nzon_vg, apply_budget_global, budemig, sum_emission
#endif
            ! variables for writing the 1x1 flux
            use user_output_flux1x1,    only : flux1x1_3d, write_flux1x1, calculate_flux1x1_indices, grid_translate_1x1

            integer, intent(in)              :: region
            type(TDate), intent(in)          :: tr(2)
            integer, intent(in)              :: itracer
            integer, intent(out)             :: status
            character(len=*), parameter      :: rname = mname//'/Emission_Fwd_Apply'
            integer             ::  i, j, l, ct
            integer, pointer    ::  zoomed(:,:)
        !    real, pointer       ::  m(:,:,:)
            real, pointer       ::  rm(:,:,:,:)
            real, pointer       ::  rzm(:,:,:,:)
        !    real, pointer       ::  blh(:,:) ! boundary layer height [m]
            real, pointer       ::  em(:,:,:,:)
        !    real, pointer       ::  t(:,:,:)
            real, allocatable   ::  ref_em_apri(:,:,:)
            real, allocatable   ::  ref_em(:,:,:)
            real, allocatable   ::  frac(:,:,:)
            real, allocatable   ::  emfac(:,:,:)
            real                ::  x, kr
            integer             ::  idater(6), idate_mid(6)
            real                ::  dtime, temp
            integer             ::  i_period, i_cat, itemp
        !    integer             ::  year, month, day
            integer             ::  iw
            integer             ::  i_period_budget, i_day, i_period_dcycle
            integer             ::  nzone, nzone_v
            logical             ::  apply_dailycycle
            integer             ::  dailycycle_type
            ! some variables for writing the 1x1 flux
            real                :: weight, f_1x1
            integer             :: flux1x1_tidx, i_1x1, j_1x1, il

            ! check whether to read new daily cycle file
            apply_dailycycle = tracers_em_info(region)%tracer(itracer)%dailycycle%apply ! This T/F will come from an rc file, read in by emission_fwd/emission_fwd_setup

            if (apply_dailycycle) then
                dailycycle_type = tracers_em_info(region)%tracer(itracer)%dailycycle%dtype
                ! find out the day number
                i_day = get_num_days(itaui, itaur(region) + ndyn/4/tref(region))
                call tau2date(itaur(region) + ndyn/4/tref(region), idate_mid)
                if ( .not. tracers_em_info(region)%tracer(itracer)%dailycycle%day_opened(i_day) ) then
                    ! close_dailycycle only needs the date to write the emissions applied, in case verbose_debug_output is defined
                    ! in that case, we should pass a date that's a little before itaur(region), otherwise the emission applied on day 1
                    ! gets written to a file with a filename corresponding to day 2
                    if ( .not. newsrun ) call close_dailycycle(region, itracer)
                    call read_dailycycle(i_day, idate_mid, region, itracer, status)
                    IF_NOTOK_RETURN(status=1)
                end if
            end if

            ! check ...
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

            ! short:
            zoomed => region_dat(region)%zoomed
            rm     => mass_dat(region)%rm_t
            rzm    => mass_dat(region)%rzm_t

            !em     => emissions(region)%surf   !  kg/cl/s or factor around 0.0 or 1.0

            ! storage:
            allocate( emfac      (im(region),jm(region),1) )

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

            ! ~ loop over categories:
            do i_cat = 1, tracers_em_info(region)%tracer(itracer)%n_cat

                ! point to emission array:
                em     => ref_emissions_apri(region)%tracer(itracer)%cat(i_cat)%field   !  kg/cl/s or factor around 0.0 or 1.0

                ! calculate the time step index for daily cycle, if it is to be applied
                apply_dailycycle = tracers_em_info(region)%tracer(itracer)%dailycycle%cycle_cat(i_cat)%apply ! this is true only if the daily cycle file read in previously has a cycle for this category
                if (apply_dailycycle) then
                    i_period_dcycle = (idate_mid(4) * 3600 + idate_mid(5) * 60 + idate_mid(6)) / tracers_em_info(region)%tracer(itracer)%dailycycle%cycle_cat(i_cat)%dtime + 1 ! this is the time index within one day
                endif

                ! get period index for emissions:
                call Time_Profile_Index( tracers_em_info(region)%tracer(itracer)%cat(i_cat)%time_profile, NewDate(time6=idater), i_period, status, at_left_side=.true. )
                IF_NOTOK_RETURN(status=1)

                !
                ! The emission factor contains the release rate:
                !
                !  c(t+dt) = c(t) + emfac(t) dt
                !

                ! apply in-month profiles to a-priori emissions if necessary:
                ! ~ current time values:
                !      year  = idater(1)
                !      month = idater(2)
                !      day   = idater(3)
                !      ! weighted ?
                !      if ( associated(Em_nday_Info(year,i_cat)%p) ) then
                !        ! current weight index:
                !        iw = Em_nday_Info(year,i_cat)%p%iw(year,month,day)
                !        ! daily average is fraction ~1.0 of monthly average:
                !        frac = Em_nday_Info(year,i_cat)%p%w(iw)%distr(region)%field
                !        ! change from monthly average to daily average:
                !        ref_em_apri = ref_em_apri * frac  ! kg/cl/s
                !        ref_em      = ref_em      * frac  ! kg/cl/s
                !      end if

                ! em contains real emissions in kg/cl/s
                !   c(t+dt) = c(t) + em dt
                emfac = em(:,:,:,i_period)   ! kg/cl/s

                !write (gol,'("add emissions on ",i4,2("-",i2.2)," ",i2.2,2(":",i2.2)," over ",f8.2," sec")') idate, dtime; call goPr

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
                                    x = emfac(i,j,1) * dtime * &
                                        tracers_em_info(region)%tracer(itracer)%dailycycle%cycle_cat(i_cat)%scaling(i,j,i_period_dcycle)
                                case (1)
                                    x = (emfac(i,j,1) + &
                                        tracers_em_info(region)%tracer(itracer)%dailycycle%cycle_cat(i_cat)%anomaly(i,j,i_period_dcycle)) * dtime
                            end select
                        else
                            x = emfac(i,j,1) * dtime   ! kg/cl/timestep
                        end if ! apply_dailycycle
                        rm(i,j,1,itracer)  = rm(i,j,1,itracer) + x
#ifdef slopes
                        rzm(i,j,1,itracer) = rzm(i,j,1,itracer) - x
#endif
                        ! 1x1 flux
                        if (write_flux1x1) then
                            do il = 1, grid_translate_1x1(region)%cell(i,j)%N
                                i_1x1 = grid_translate_1x1(region)%cell(i,j)%ilist(il)
                                j_1x1 = grid_translate_1x1(region)%cell(i,j)%jlist(il)
                                f_1x1 = grid_translate_1x1(region)%cell(i,j)%frac(il)
                                flux1x1_3d%prod(i_1x1,j_1x1,1,flux1x1_tidx,itracer) = flux1x1_3d%prod(i_1x1,j_1x1,1,flux1x1_tidx,itracer) + &
                                    f_1x1 * x
                            end do
                        end if

                        !=========
                        ! budget
                        !=========
#ifdef with_budgets
                        if ( apply_budget_global ) then

                            !global budget
                            nzone = budg_dat(region)%nzong(i,j)
                            nzone_v = nzon_vg(1)

                            budemig(nzone,nzone_v,itracer,i_period_budget) = budemig(nzone,nzone_v,itracer,i_period_budget) + x ! [kg]
                            !budemi_dat(region)%emi(i,j,nzone_v,itracer) = &
                            !     budemi_dat(region)%emi(i,j,nzone_v,itracer) + x ! [kg]

                            sum_emission(region,itracer) = sum_emission(region,itracer) + x  ! [kg]
                        end if
#endif
                    end do ! i
                end do ! j
            end do ! i_cat
            ! unlink:
            nullify( zoomed )
        !    nullify( m )
            nullify( rm )
            nullify( rzm )
            nullify( em )
        !    nullify( t )
        !    nullify( blh )

            ! clear:
            deallocate( emfac       )

            ! ok:
            status = 0

        end subroutine Emission_Fwd_Apply

!  subroutine read_dailycycle(i_day, idate_local, region, itrac)

!    use dims, only : nregions, region_name, sec_day
!    use Emission_data, only : tracers_em_info
!    use file_netcdf

!    implicit none

!    !__IO___________________________________________________________________

!    integer, intent(in) :: i_day, idate_local(6), region, itrac

!    !__CONSTANTS_____________________________________________________________

!    character(len=*), parameter      :: rname = mname//'/read_dailycycle'

!    !__LOCAL_VARIABLES______________________________________________________

!    integer              :: nc_id, grp_id, i_cat, cat_grp_id
!    character(len=256)   :: fname
!    logical              :: file_exist
!    integer              :: n_cat, n_tstep, status
!    integer              :: dailycycle_type
!    character(len=256)   :: dailycycle_pfx, tracer_name, emis_indir
!    character(len=*), parameter :: rname = mname//'/read_dailycycle'

!    !__START_SUBROUTINE______________________________________________________

!    dailycycle_pfx = tracers_em_info(region)%tracer(itrac)%dailycycle%pfx
!    dailycycle_type = tracers_em_info(region)%tracer(itrac)%dailycycle%dtype
!    tracer_name = tracers_em_info(region)%tracer(itrac)%name
!    n_cat = tracers_em_info(region)%tracer(itrac)%n_cat
!    emis_indir = tracers_em_info(region)%tracer(itrac)%dailycycle%emis_indir

!    write(fname,"(a,a1,a,i4.4,2i2.2,a4)") trim(emis_indir), '/', trim(dailycycle_pfx), idate_local(1:3), '.nc4'

!    inquire( file=fname, exist=file_exist )
!    if ( file_exist ) then
!        nc_id = nc_open(trim(fname), 'r')
!        grp_id = nc_get_group(nc_id, region_name(region), status)
!        IF_NOTOK_RETURN(status=1)
!        do i_cat = 1,n_cat
!            if (nc_grp_exists(grp_id, tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name)) then
!                cat_grp_id = nc_get_group(grp_id, tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name, status)
!                IF_NOTOK_RETURN(status=1)
!                select case (dailycycle_type)
!                    case (0)
!                        tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%scaling = nc_read_var(cat_grp_id, 'emission_scaling_factor', status)
!                        IF_NOTOK_RETURN(status=1)
!                    case (1)
!                        tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%anomaly = nc_read_var(cat_grp_id, 'emission_anomaly', status)
!                        IF_NOTOK_RETURN(status=1)
!                end select
!                n_tstep = nc_get_dim(cat_grp_id, 'timesteps', status)
!                IF_NOTOK_RETURN(status=1)
!                tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%dtime = sec_day/n_tstep
!                tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%apply = .true.
!            else
!                tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%apply = .false.
!            end if ! category exists in dailycycle file
!        end do ! i_cat
!        call nc_close(nc_id)
!        tracers_em_info(region)%tracer(itrac)%dailycycle%day_opened(i_day) = .true.
!        !write(*,'(a,a,a)') ' Diurnal cycle file ', trim(fname), ' read in'
!    else
!        write(0,'(a,a,a)') ' ERROR :: Diurnal cycle file ', trim(fname), ' does not exist'
!        status=1
!        IF_NOTOK_RETURN(status=1)
!    end if

!  end subroutine read_dailycycle

!  subroutine close_dailycycle(region, itrac)

!    use Emission_Data, only : tracers_em_info

!    implicit none

!    integer, intent(in) :: region, itrac

!    character(len=*), parameter      :: rname = mname//'/close_dailycycle'
!    integer :: n_cat, i_cat

!    n_cat = tracers_em_info(region)%tracer(itrac)%n_cat

!    do i_cat = 1, n_cat
!        if (allocated(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%scaling)) &
!            deallocate(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%scaling)
!        if (allocated(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%anomaly)) &
!            deallocate(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%anomaly)
!        tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%apply = .false.
!    end do ! i_cat

!  end subroutine close_dailycycle

end module Emission_Fwd_CO2
