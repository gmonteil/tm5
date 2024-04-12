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

    use Emission_Adj_CO2,   only : Init_adjemis_CO2 => Emission_Adj_Init
    use dims,               only : nregions
    use Emission_Data,      only : tracers_em_info
    use chem_param,         only : ntracet
    use TM5_Fields,         only : Fields_4D_Init

    integer, intent(out)            :: status
    character(len=*), parameter     :: rname = mname//'/Emission_Adj_Init'
    integer                         :: region, itrac, i_cat
    integer, allocatable            :: nt(:)

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

    call Init_adjemis_CO2(status)
    IF_NOTOK_RETURN(status=1)

   ! ok
    status = 0

  end subroutine Emission_Adj_Init

  subroutine Emission_Adj_Done( status )

    use Emission_Adj_CO2,   only : Done_adjemis_CO2 => Emission_Adj_Done

    use Emission_Data,  only : tracers_em_info, update_parent_emissions
    use dims,           only : revert, nregions, region_name, im, jm
    use global_data,    only : rcF
    use GO,             only : ReadRc
    use chem_param,     only : ntracet
    use TM5_Fields,     only : Fields_4D_Done
    use file_netcdf

    integer, intent(out)        :: status
    character(len=*), parameter :: rname = mname//'/Emission_Adj_Done'
    character(len=512)          :: outdir_adj_emis, outFile
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
                 print*, i_cat, c_cat, sum(adj_emissions(region)%tracer(itrac)%cat(i_cat)%field(:,:,1,:))
              end do !i_cat
           end do ! itrac
        end do ! region
        call nc_close(nc_io, status)
        IF_NOTOK_RETURN(status=1)

    end if ! revert

    ! clear:
!    do region = 1, nregions
!       do itrac = 1, ntracet
!          call Fields_4D_Done( adj_emissions(region)%tracer(itrac), tracers_em_info(region)%tracer(itrac)%n_cat,  status )
!          IF_NOTOK_RETURN(status=1)
!       end do
!    end do
!    deallocate( adj_emissions )

    call Done_adjemis_CO2(status)
    IF_NOTOK_RETURN(status=1)

  end subroutine Emission_Adj_Done

  subroutine Emission_Adj_Setup( status )

    ! --- modules ------------------------------

    use Dims,                   only : nregions, itaui, itaue
    use datetime,               only : time_window, get_num_days
    use Emission_Data,          only : ref_emissions_apri, tracers_em_info, adj_emissions
    use Emission_Read_PyShell,  only : Read_Emissions_From_PyShell
    use chem_param,             only : ntracet
    use global_data,            only : rcF
    use GO,                     only : ReadRc
    use os_specs,               only : DUMMY_STR_LEN

    ! --- in/out ----------------------------------------------

    integer, intent(out)             :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/Emission_Adj_Setup'

    !--- local ------------------------------------------------

    integer                         :: i_cat, itrac, region, n_cat, n_day, i
    character(len=8)                :: tracer_name
    character(len=256)              :: emis_indir ! where are the daily cycle files?
    logical                         :: apply_dailycycle, dcycle_cat ! whether or not to apply the diurnal cycle
    character(len=256)              :: dailycycle_pfx ! the prefix for the daily cycle file
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
                if (dcycle_cat) write(gol,'(a, ": tracer ",a," category ",a," has sub-period chunking")') rname, trim(tracer_name), trim(cat_name) ; call goPr
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
                n_cat = tracers_em_info(region)%tracer(itrac)%n_cat
                allocate(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(n_cat))
                tracers_em_info(region)%tracer(itrac)%dailycycle%apply = apply_dailycycle
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

    use Emission_Adj_CO2,   only : Apply_adjemis_CO2 => Emission_Adj_Apply
    use Go,                 only : TDate

    integer, intent(in)              :: region
    type(TDate), intent(in)          :: tr(2)
    integer, intent(out)             :: status

    character(len=*), parameter      :: rname = mname//'/Emission_Adj_Apply'

    call Apply_adjemis_CO2(region, tr, status)
    IF_NOTOK_RETURN(status=1)

  end subroutine Emission_Adj_Apply

end module Emission_Adj
