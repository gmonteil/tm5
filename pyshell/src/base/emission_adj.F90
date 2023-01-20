!###############################################################################
!
!       purpose
!       -------
!       perform adjoint emissions needed for TM5 4DVAR
!
!       interface
!       ---------
!       call Emission_Adj_Init
!       call Emission_Adj_Apply
!       call Emission_Adj_Done
!
!       method
!       ------
!       subroutine Emission_Adj_Init    is called from adj_trace0
!       subroutine Emission_Adj_Apply   is called from adj_source1
!       subroutine Emission_Adj_Done    is called from adj_trace_end
!
!       Note that this file contains a dummy copy, and does nothing tracer specific. In the multi-tracer code,
!       there are tracer-generic and tracer-specific things to be done by emission_adj. In this file, which is
!       a copy of the file from the project 'tracers/COCO2', all the tracer-specific parts have been commented
!       out. If you want to write a multi-tracer code, you need to uncomment and modify those parts.
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"

module Emission_Adj

  use GO, only : gol, goErr, goPr
  use Emission_data, only : adj_emissions

  implicit none

  private

  ! public routines
  public :: Emission_Adj_Init
  public :: Emission_Adj_Setup
  public :: Emission_Adj_Apply
  public :: Emission_Adj_Done

  ! --- const ------------------------------

  character(len=*), parameter        :: mname = 'Emission_Adj'


contains

  subroutine Emission_Adj_Init( status )

    ! import tracer-specific modules, such as
    ! use emission_adj_sf6, only : init_adjemis_sf6 => emission_adj_init
    use dims,               only : nregions
    use dims,               only : adv_scheme
    use Emission_Data,      only : tracers_em_info
    use chem_param,         only : ntracet
    use TM5_Fields,         only : Fields_4D_Init

    integer, intent(out)            :: status
    character(len=*), parameter     :: rname = mname//'/Emission_Adj_Init'
    integer                         :: region, itrac, i_cat
    integer, allocatable            :: nt(:)

    ! Why is this check here? Who knows...
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

    ! emission arrays:
    allocate( adj_emissions(nregions) )

    ! loop over var4d categories:
    do region = 1, nregions
       allocate( adj_emissions(region)%tracer(ntracet) )
       do itrac = 1, ntracet
          allocate(nt(tracers_em_info(region)%tracer(itrac)%n_cat))
          do i_cat = 1, tracers_em_info(region)%tracer(itrac)%n_cat
             nt(i_cat) = tracers_em_info(region)%tracer(itrac)%cat(i_cat)%time_profile%n_period
          end do
          call Fields_4D_Init( adj_emissions(region)%tracer(itrac), region, .false.,&
                               tracers_em_info(region)%tracer(itrac)%n_cat, nt , status )
          IF_NOTOK_RETURN(status=1)
          deallocate(nt)
       end do
    end do

    ! call tracer-specific routines, such as
    ! call init_adjemis_sf6(status)
    ! IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Emission_Adj_Init

  subroutine Emission_Adj_Done( status )

    ! import tracer-specific modules, such as
    ! use emission_adj_sf6, only : done_adjemis_sf6 => emission_adj_done
    use Emission_Data,  only : tracers_em_info, update_parent_emissions
    use dims,           only : revert, nregions, region_name, im, jm
    use global_data,    only : rcF
    use GO,             only : ReadRc
    use chem_param,     only : ntracet
    use TM5_Fields,     only : Fields_4D_Done
    use file_netcdf
    use os_specs,       only : MAX_FILENAME_LEN

    integer, intent(out)        :: status
    character(len=*), parameter :: rname = mname//'/Emission_Adj_Done'
    character(len=MAX_FILENAME_LEN) :: outdir_adj_emis, outFile
    integer                     :: nc_io, region, rgrp, itrac, tgrp, i_cat, ntime, cgrp
    character(len=40)           :: c_cat

    if (revert == -1) then

        ! First update parent adj_emissions
        ! For some reason Arjo has this call, but for me it causes the adjoint test to fail...
!        call update_parent_emissions( adj_emissions )

        call ReadRc(rcF, 'output.dir', outdir_adj_emis, status)
        IF_NOTOK_RETURN(status=1)

        write (outFile,'(a,"/adj_emissions.nc4")') trim(outdir_adj_emis)
        nc_io = nc_open(trim(outFile),'c',status)
        IF_NOTOK_RETURN(status=1)

        call nc_create_dim(nc_io,'nregions',nregions)
        call nc_create_dim(nc_io,'tracer',ntracet)
        do region = 1, nregions
           rgrp = nc_create_group(nc_io,region_name(region),status)
           do itrac = 1, ntracet
              tgrp = nc_create_group(rgrp,TRIM(tracers_em_info(region)%tracer(itrac)%name),status)
              call nc_create_dim(tgrp,'lon ',im(region),status)
              call nc_create_dim(tgrp,'lat ',jm(region),status)
              call nc_create_dim(tgrp,'ncat',tracers_em_info(region)%tracer(itrac)%n_cat,status)
              do i_cat = 1, tracers_em_info(region)%tracer(itrac)%n_cat
                 c_cat = tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name
                 ntime = tracers_em_info(region)%tracer(itrac)%cat(i_cat)%time_profile%n_period
                 cgrp = nc_create_group(tgrp,trim(c_cat),status)
                 call nc_create_dim(cgrp,'time',ntime,status)
                 call nc_dump_var(cgrp,'adj_emis',(/'lon ','lat ','time'/),  &
                      adj_emissions(region)%tracer(itrac)%cat(i_cat)%field(:,:,1,:))
              end do !i_cat
           end do ! itrac
        end do ! region
        call nc_close(nc_io, status)
        IF_NOTOK_RETURN(status=1)

    end if ! revert

    ! clear:
    do region = 1, nregions
       do itrac = 1, ntracet
          call Fields_4D_Done( adj_emissions(region)%tracer(itrac), tracers_em_info(region)%tracer(itrac)%n_cat,  status )
          IF_NOTOK_RETURN(status=1)
       end do
    end do
    deallocate( adj_emissions )

    ! call tracer-specific routines, such as
    ! call done_adjemis_sf6(status)
    ! IF_NOTOK_RETURN(status=1)

  end subroutine Emission_Adj_Done

  subroutine Emission_Adj_Setup( status )

    ! --- modules ------------------------------

    use Dims,                   only : nregions, itaui, itaue
    use datetime,               only : time_window, get_num_days
    use Emission_Data,          only : ref_emissions_apri, tracers_em_info, adj_emissions
    use Emission_Read_PyShell,  only : Read_Emissions_From_PyShell
    use chem_param,             only : ntracet, tracer_name_len
    use global_data,            only : rcF
    use GO,                     only : ReadRc
    use os_specs,               only : MAX_FILENAME_LEN, SHORT_STR_LEN, DUMMY_STR_LEN

    ! --- in/out ----------------------------------------------

    integer, intent(out)             :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/Emission_Adj_Setup'

    !--- local ------------------------------------------------

    integer                         :: i_cat, itrac, region, n_cat, n_day, i
    character(len=tracer_name_len)  :: tracer_name
    character(len=MAX_FILENAME_LEN) :: emis_indir ! where are the daily cycle files?
    logical                         :: apply_dailycycle, dcycle_cat ! whether or not to apply the diurnal cycle
    character(len=SHORT_STR_LEN)    :: dailycycle_pfx ! the prefix for the daily cycle file
    integer                         :: dailycycle_type ! 0 for scaling, 1 for adding daily cycles over smaller time steps
    character(len=DUMMY_STR_LEN)    :: cat_name

    !--- begin ------------------------------------------------

    ! loop over categories:
    ! loop over tracers:
    ! loop over regions:
    do region = 1, nregions
       do itrac = 1, ntracet
          do i_cat = 1, tracers_em_info(region)%tracer(itrac)%n_cat
             ! set adjoint emissions to zero:
             adj_emissions(region)%tracer(itrac)%cat(i_cat)%field = 0.0
          end do
       end do
    end do

    call Read_Emissions_From_PyShell ( tracers_em_info, time_window, &
                                        ref_emissions_apri, status )
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
        end if ! apply_dailycycle
    end do ! itrac

    ! ok
    status = 0

  end subroutine Emission_Adj_Setup

  subroutine Emission_Adj_Apply( region, tr, status )

    ! import tracer-specific routines, such as
    ! use emission_adj_sf6, only : apply_adjemis_sf6 => emission_adj_apply
    use Go,                 only : TDate

    integer, intent(in)              :: region
    type(TDate), intent(in)          :: tr(2)
    integer, intent(out)             :: status

    character(len=*), parameter      :: rname = mname//'/Emission_Adj_Apply'

    ! At this point, we can call tracer-specific routines, such as
    ! call apply_adjemis_sf6(region, tr, status)
    ! IF_NOTOK_RETURN(status=1)
    !
    ! That is needed, e.g., when there is something tracer-specific about the adjoint emission code, such as injecting adjoint tracer mass
    ! stuff high up such as for biomass burning CO. However, for a generic passive tracer with only surface fluxes, there is no need for a
    ! tracer-specific routine. So this module includes a generic, tracer-agnostic routine for injecting the adjoint emissions. If you call
    ! tracer-specific routines, either comment out the call to Apply_AdjEmis_generic,  or take care not to replicate that function in your
    ! tracer-specific routines.

    call Apply_AdjEmis_generic(region, tr, status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Emission_Adj_Apply

  subroutine Apply_AdjEmis_generic(region, tr, status)

    ! --- modules ------------------------------

    use GO                     , only : TDate, NewDate, Get
    use GO                     , only : Time_Profile_Index
    use GO                     , only : operator(+), operator(-), operator(/), rTotal
    use dims,                    only : tref, idatee, itaue, sec_day
    use dims,                    only : im, jm, isr, ier, jsr, jer, lm
    use dims,                    only : itau, itaur, ndyn, itaui, newsrun
    use global_data,             only : mass_dat, region_dat
    use datetime,                only : tau2date, get_num_days
    use Emission_Data          , only : tracers_em_info, ntracet
    use Emission_Read_PyShell  , only : read_dailycycle, close_dailycycle

    use Emission_Data          , only : ref_emissions, ref_emissions_apri

    ! --- in/out ----------------------------------------------

    integer, intent(in)              :: region
    type(TDate), intent(in)          :: tr(2)
    integer, intent(out)             :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/Apply_AdjEmis_generic'

    !--- local ------------------------------------------------

    integer             ::  i, j, l, itr
    integer, pointer    ::  zoomed(:,:)
    real, pointer       ::  adj_rm(:,:,:,:)
    real, pointer       ::  adj_rzm(:,:,:,:)
    real, pointer       ::  adj_em(:,:,:,:)
    real                ::  x
    type(TDate)         ::  tread
    integer             ::  idater(6)
    real                ::  dtime
    integer             ::  year, month, day
    integer             ::  iw
    integer             ::  i_cat, i_day, idate_mid(6)
    integer             ::  i_period, i_period_dcycle
    logical             ::  apply_dailycycle
    integer             ::  dailycycle_type

    !--- begin ------------------------------------------------

    ! short:
    zoomed  => region_dat(region)%zoomed
    adj_rm  => mass_dat(region)%rm_t
    adj_rzm => mass_dat(region)%rzm_t

    ! timestep emissions (should be tr(1)-tr(2), but the 'abs' takes care of that)
    dtime = abs(rTotal( tr(2) - tr(1), 'sec' ))

    ! mid time:
    call tau2date( itau - nint( dtime ), idater )
    tread = NewDate( time6=idater )

    ! current time values:
    call Get( tread, year=year, month=month, day=day )

    do itr = 1, ntracet
        ! read the daily cycle, if necessary
        apply_dailycycle = tracers_em_info(region)%tracer(itr)%dailycycle%apply
        dailycycle_type = tracers_em_info(region)%tracer(itr)%dailycycle%dtype

        if (apply_dailycycle) then
            ! find out the day number
            i_day = get_num_days(itaui, itaur(region) - ndyn/4/tref(region))
            call tau2date(itaur(region) - ndyn/4/tref(region), idate_mid)
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

            ! short:
            adj_em => adj_emissions(region)%tracer(itr)%cat(i_cat)%field   !  kg/cl/s or factor around 0.0 or 1.0

            ! calculate the time step index for daily cycle
            apply_dailycycle = tracers_em_info(region)%tracer(itr)%dailycycle%cycle_cat(i_cat)%apply
            if (apply_dailycycle) &
                i_period_dcycle = (idate_mid(4) * 3600 + idate_mid(5) * 60 + idate_mid(6)) / tracers_em_info(region)%tracer(itr)%dailycycle%cycle_cat(i_cat)%dtime + 1
                ! This points to the time index within one day

            ! get period index:
            call Time_Profile_Index( tracers_em_info(region)%tracer(itr)%cat(i_cat)%time_profile, tread, i_period, status, at_left_side=.true. )
            IF_NOTOK_RETURN(status=1)

            !
            ! The emission weight contains the factors that should be
            ! multplied with the emission field in the state to give
            ! the release rate:
            !
            !  c(t+dt) = c(t) + em emw(t) dt
            !
            ! em contains real emissions in kg/cl/s
            !   c(t+dt) = c(t) + em dt

            !
            ! Emission operator:
            !
            !  c(t+dt) = c(t) + emfac(t) dt
            !          ~ c(t) + em ew(t) dt     ! linearized
            !
            ! where :
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
            !  [ adj_c  ]   [   I   O ] [ adj_c  ]
            !  [        ] = [         ] [    ]
            !  [ adj_em ]   [ ew^dt I ] [ adj_em ]
            !
            ! where:
            !   adj_c   : adjoint concentration [cl kg-1]
            !   adj_em  : adjoint emission field [cl s kg-1]
            !
            ! which gives the adjoint operator:
            !
            !   adj_em =  adj_c ew dt  +  adj_em
            !

            do j = jsr(region), jer(region)
                do i = isr(region), ier(region)

                    if ( zoomed(i,j) /= region ) cycle

                    if (apply_dailycycle) then
                        select case (dailycycle_type)
                        case (0)
                            x = tracers_em_info(region)%tracer(itr)%dailycycle%cycle_cat(i_cat)%scaling(i,j,i_period_dcycle) * dtime
                        case (1)
                            x = dtime
                        end select
                    else
                        x = dtime
                    end if ! apply_dailycycle

                    adj_em(i,j,1,i_period) = adj_em(i,j,1,i_period) + adj_rm(i,j,1,itr) * x
#ifdef slopes
                    adj_em(i,j,1,i_period) = adj_em(i,j,1,i_period) - adj_rzm(i,j,1,itr) * x
#endif

                end do ! i
            end do ! j

            nullify( adj_em  )

        end do ! i_cat
    end do ! itr

    ! unlink:
    nullify( zoomed )
    nullify( adj_rm  )
    nullify( adj_rzm )

    ! ok:
    status = 0

  end subroutine Apply_AdjEmis_generic

end module Emission_Adj
