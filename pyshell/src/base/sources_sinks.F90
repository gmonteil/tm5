!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module sources_sinks
    !-----------------------------------------------------------------------
    !       purpose
    !       -------
    !       - apply CH4 sources for CH4 FW run
    !       - read k_CH4 file (trace0)
    !       Peter Bergamaschi 12/2004
    !
    !       interface
    !       ---------
    !       call trace0
    !       call trace1
    !       call trace_after_read
    !       call source1
    !       call source2
    !       call trace_end
    !
    !       method
    !       ------
    !       subroutine trace0 is called from start, and is used
    !                to input tracer specific values.
    !       subroutine trace1 is called from start when istart = 2, and provides
    !                a way to do some initialization, that has to be done only
    !                at the beginning of a run.
    !       subroutine trace_after_read is called when new meteo comes available.
    !       subroutine source1 is called every nsrce seconds and provides
    !                the source/sink calculation.
    !       subroutine source2 is called every nsrce seconds and provides
    !                boundary conditions.
    !       subroutine trace_end is called at the end, and deallocates variables
    !                and provides budget output.
    !       in zoom mode, recursive calling of the children is performed
    !
    !-----------------------------------------------------------------------

    ! --- modules ------------------------------

    use dims,       only : nregions
    use chem_param, only : chem_data
    use GO,         only : gol, goErr, goPr
    use os_specs,   only : MAX_FILENAME_LEN, DUMMY_STR_LEN

    implicit none

    ! --- in/out -----------------------------

    private

    public :: Sources_Sinks_Init, Sources_Sinks_Done
    public :: trace0
    public :: trace1
    !  public :: trace_after_read
    public :: source1
    public :: source2
    public :: trace_end

    ! --- const ------------------------------

    character(len=*), parameter     ::  mname = 'module sources_sinks'

    ! --- var --------------------------------


    !  ! diffusion files:
    !  logical                         ::  with_diffusion_files
    !  character(len=MAX_FILENAME_LEN)              ::  diffusion_dir

    ! timers
    integer     ::  itim_trace0, itim_trace1
    integer     ::  itim_trace_after_read
    integer     ::  itim_calc_kzz
    integer     ::  itim_read_diffusion
    integer     ::  itim_write_diffusion
    integer     ::  itim_read_surface
    integer     ::  itim_dry_deposition


contains


    ! ============================================================================


    subroutine Sources_Sinks_Init( status )

        use GO, only : GO_Timer_Def
        use GO, only : ReadRc
        use global_data, only : rcF

        ! --- in/out -------------------------------

        integer, intent(out)              ::  status

        ! --- const ----------------------------

        character(len=*), parameter  ::  rname = mname//'/Sources_Sinks_Init'

        ! --- begin --------------------------------

        ! define timers:
        call GO_Timer_Def( itim_trace0          , 'trace0', status )
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_Def( itim_trace1          , 'trace1', status )
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_Def( itim_trace_after_read, 'trace_after_read', status )
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_Def( itim_calc_kzz        , 'calc Kzz', status )
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_Def( itim_read_diffusion  , 'read diffusion', status )
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_Def( itim_write_diffusion , 'write diffusion', status )
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_Def( itim_read_surface    , 'read surface file', status )
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_Def( itim_dry_deposition  , 'dry deposition', status )
        IF_NOTOK_RETURN(status=1)

        !    ! read previously stored diffusion fields ?
        !    call ReadRc( rcF, 'diffusion.files', with_diffusion_files, status )
        !    IF_NOTOK_RETURN(status=1)
        !    ! directory for temporary diffusion files:
        !    call ReadRc(rcF, 'diffusion.dir', diffusion_dir, status)
        !    IF_NOTOK_RETURN(status=1)

        ! ok
        status = 0

    end subroutine Sources_Sinks_Init


    ! ***


    subroutine Sources_Sinks_Done( status )

        ! --- in/out -------------------------------

        integer, intent(out)              ::  status

        ! --- const ----------------------------

        character(len=*), parameter  ::  rname = mname//'/Sources_Sinks_Done'

        ! --- begin --------------------------------

        ! ok
        status = 0

    end subroutine Sources_Sinks_Done


    !=====================================================================================================
    !=====================================================================================================


    subroutine trace0( status )
        !-----------------------------------------------------------------------
        !
        !       trace0
        !
        !
        !       purpose
        !       -------
        !       initialise data needed for chemistry,
        !       e.g. read emissions and boundary conditions
        !       routine is calculated at start and at beginning each month...
        !
        !       interface
        !       ---------
        !       call trace0
        !
        !       method
        !       ------
        !       subroutine trace0 is called from module initexit, and provides a way
        !                to input tracer specific values.
        !
        !       external
        !       --------
        !           calc_sm
        !           ini_zoneg
        !           diagbudg
        !
        !       reference
        !       ---------
        !       see above
        !---------------------------------------------------------------------

        use dims,          only : idate, dxy11, nlat180, newsrun
        !    use dims,          only : sec_day,sec_month,sec_year, mlen
        use dims,          only : newsrun
        use dims,          only : revert
        use geometry,      only : calc_dxy
#ifdef with_budgets
        use budget_global, only : diagbudg, ini_zoneg
        use budget_global, only : apply_budget_global
#endif
        use GO           , only : GO_Timer_Start, GO_Timer_End
        use GO,            only : ReadRc
        use global_data,   only : rcF
        use toolbox,       only : escape_tm
!#ifndef without_chemistry
        !use chemistry,     only : Chemistry_Init !, Chemistry_Setup
!#endif


#ifdef MPI
        use mpi_const,only : myid,lmloc
#endif

        implicit none

        !__IO___________________________________________________________________

        integer, intent(out)   :: status

        !__LOCAL_VARIABLES______________________________________________________

        integer                :: region

        !__CONST________________________________________________________________


        character(len=*), parameter      :: rname = mname//'/trace0'

        !__START_SUBROUTINE______________________________________________________

        call GO_Timer_Start( itim_trace0, status )
        IF_NOTOK_RETURN(status=1)

        !=======================
        ! general
        !=======================

        !    print * , 'trace0: newmonth ', newmonth
        !    print * , 'trace0: start month ', idate(2)

        ! calculate some conversion factors related to time...
        !    call calc_sm(mlen,sec_day,sec_month,sec_year)

        !=======================
        ! budget
        !=======================

#ifdef with_budgets
        if ( newsrun ) then

            call ReadRc(rcF, 'budget.global', apply_budget_global, status)
            IF_NOTOK_RETURN(status=1)


            if ( apply_budget_global ) then

                ! the global budget
                call ini_zoneg( status )
                IF_NOTOK_RETURN(status=1)

                do region = 1, nregions
                    call diagbudg(region,1,status)
                    IF_NOTOK_RETURN(status=1)      ! initialize budget files global
                end do

            end if

        end if
#endif

        !=======================
        ! surface fields
        !=======================

        !    PB now called from initexit
        !    ! fields for surface processes (1x1 --> coarsened)
        !    if ( newsrun ) call declare_surface_fields

        !=======================
        ! emissions
        ! JFM: this is now called from initexit: start_TM5
        !  because it depends on run_mode and iteration whether this has to be done.
        !=======================

        !    >>> now called from modelIntegration/Proces_Init
        !        and every timestep from modelIntegration/Process_Setup
!#ifndef without_chemistry
        !    !=======================
        !    ! chemistry (OH sink)
        !    !=======================
        !    ! first call ? overall init only:
        !    if (newsrun) then
        !        call chemistry_init( status )
        !        IF_NOTOK_RETURN(status=1)
        !    end if
        !    ! setup tasks to per forformed every time:
        !    call Chemistry_Setup( status )
        !    IF_NOTOK_RETURN(status=1)
!#endif
        !    <<<

        !    print * , 'trace0: end month ', idate(2)

        call GO_Timer_End( itim_trace0, status )
        IF_NOTOK_RETURN(status=1)

        ! ok
        status = 0

    end subroutine trace0


    !=====================================================================================================
    !=====================================================================================================

    subroutine trace1( status )
        !-----------------------------------------------------------------------
        !
        !       purpose
        !       -------
        !       this subroutine determines the initial values for tracer mass and
        !       and its slopes.  it is called from input when istart = 2.
        !
        !       interface
        !       ---------
        !       call trace1
        !
        !       method
        !       ------
        !
        !       reference
        !       ---------
        !       see above
        !
        !-----------------------------------------------------------------------

        use GO         , only : GO_Timer_Start, GO_Timer_End, ReadRc
        use dims,        only : im, jm, lm, adv_scheme, region_name
        use toolbox,     only : print_totalmass
        use global_data, only : mass_dat, rcf
        use MeteoData  , only : m_dat
        use zoom_tools,  only : update_parent
        use chem_param,  only : fscale, ntracet, ntrace, ntrace_chem
        use chem_param,  only : names
        use ParTools,    only : ntracetloc,ntracet_ar,myid
        use file_netcdf

        implicit none

        !__IO___________________________________________________________________

        integer, intent(out)     :: status

        !__LOCAL_VARIABLES______________________________________________________

        character(len=*), parameter       :: rname = mname//'/trace1'
        real,dimension(:,:,:,:),pointer   :: rm, rxm, rym, rzm
        !    real,dimension(:,:,:,:),pointer   :: rmc
        real,dimension(:,:,:),  pointer   :: m
        real, allocatable                 :: mix_ini_array(:,:,:)
        integer                           :: region
        integer                           :: n,nn
        real                              :: mix_ini
        logical                           :: iniconc_from_file
        character(len=MAX_FILENAME_LEN)   :: iniconc_file_name
        integer                           :: nc_id, grp_id, tgrp_id

        !__START_SUBROUTINE______________________________________________________

        call GO_Timer_Start( itim_trace1, status )
        IF_NOTOK_RETURN(status=1)

        call ReadRc(rcF, 'start.2.iniconc_from_file', iniconc_from_file, status)
        IF_NOTOK_RETURN(status=1)

        if (iniconc_from_file) then
            call ReadRc(rcF, 'start.2.iniconcfile', iniconc_file_name, status)
            IF_NOTOK_RETURN(status=1)
        end if

        do region=1,nregions

            m => m_dat(region)%data
            rm => mass_dat(region)%rm_t
            rxm => mass_dat(region)%rxm_t
            rym => mass_dat(region)%rym_t
            rzm => mass_dat(region)%rzm_t

            do n=1,ntracet
                if (iniconc_from_file) then
                    ! read the iniconc file
                    write(gol,'(a,a)') 'Reading initial mixing ratio from ', trim(adjustl(iniconc_file_name)) ; call goPr
                    nc_id = nc_open(trim(adjustl(iniconc_file_name)), 'r', status)
                    IF_NOTOK_RETURN(status=1)

                    grp_id = nc_get_group(nc_id, region_name(region))
                    tgrp_id = nc_get_group(grp_id, names(n))
                    mix_ini_array = nc_read_var(tgrp_id, 'mixing_ratio')
                    call nc_close(nc_id)
                    write(*,'(a,a)') 'Initial mixing ratio read in from ', trim(adjustl(iniconc_file_name))
                    rm(1:im(region),1:jm(region),1:lm(region),n) = mix_ini_array(1:im(region),1:jm(region),1:lm(region)) * m(1:im(region),1:jm(region),1:lm(region)) / fscale(n)
                else
                    call ReadRc(rcF, 'start.2.iniconc.'//trim(names(n)), mix_ini, status)
                    IF_NOTOK_RETURN(status=1)
                    rm(1:im(region),1:jm(region),1:lm(region),n) = mix_ini * m(1:im(region),1:jm(region),1:lm(region)) / fscale(n)
                end if
                if ( adv_scheme == 'slope' ) then
                    rxm(:,:,:,n) = 0.0
                    rym(:,:,:,n) = 0.0
                    rzm(:,:,:,n) = 0.0
                end if
                if (allocated(mix_ini_array)) deallocate(mix_ini_array)
            end do ! tracer

            nullify(m)
            nullify(rm)
            nullify(rxm)
            nullify(rym)
            nullify(rzm)

        end do ! region

        do region = nregions,2,-1
            call update_parent(region)
        end do

        do region = 1, nregions
            do n=1,ntracet
                call print_totalmass(region,n)
            end do
        end do

        call GO_Timer_End( itim_trace1, status )
        IF_NOTOK_RETURN(status=1)

        status = 0

    end subroutine trace1


    !---------------------------------------------------------------------
    !
    !       source1
    !
    !
    !       purpose
    !       -------
    !       this subroutine changes the tracer mass and its slopes
    !       by chemical sources
    !
    !       interface
    !       ---------
    !       call source1
    !
    !       method
    !       ------
    !
    !
    !
    !       reference
    !       --------
    !       see above
    !
    !---------------------------------------------------------------------

    subroutine source1( region, tr, status )

        use GO            , only : TDate
        use emission,       only : Emission_Fwd_Apply
#ifndef without_dry_deposition
        use dry_deposition, only : Dry_Deposition_Apply
#endif

        !__IO___________________________________________________________________

        integer, intent(in)        :: region
        type(TDate), intent(in)    :: tr(2)
        integer, intent(out)       :: status

        !__CONST________________________________________________________________

        character(len=*), parameter   ::  rname = mname//'/source1'

        !__LOCAL_VARIABLES______________________________________________________

        !__START_SUBROUTINE______________________________________________________


        call Emission_Fwd_Apply( region, tr, status )
        IF_NOTOK_RETURN(status=1)

#ifndef without_dry_deposition
        call Dry_Deposition_Apply( region, tr, status )
        IF_NOTOK_RETURN(status=1)
#endif

        ! ok:
        status = 0

    end subroutine source1

    !=====================================================================================================
    !=====================================================================================================

    subroutine source2(region)
        ! ---------------------------------------------------------------------
        !
        !       source2
        !
        !
        !       purpose
        !       -------
        !       apply boundary conditions to selected tracers
        !
        !       interface
        !       ---------
        !       call source2
        !
        !       method
        !       ------
        !
        !
        !       external
        !       --------
        !       none
        !
        !       reference
        !       --------
        !       see above
        !
        !
        !---------------------------------------------------------------------

        implicit none

        !__IO___________________________________________________________________

        integer, intent(in)        :: region

        !__LOCAL_VARIABLES______________________________________________________

        !__START_SUBROUTINE______________________________________________________



    end subroutine source2

    !=====================================================================================================
    !=====================================================================================================

    subroutine trace_end(msg, file_name, status)
        !
        ! deallocate the memory and do everything needed to finalise chemistry
        !
        use dims,          only : region_name, im, jm, lm
        use chem_param,    only : ntracet, ntrace_chem, names
        use global_data,   only : mass_dat
        use io_hdf,        only : io_write, DFACC_WRITE
#ifdef with_budgets
        use budget_global, only : diagbudg
        use budget_global, only : init_mass, apply_budget_global
        use budget_global, only : sum_chemistry, sum_update
        use budget_global, only : sum_advection, sum_deposition, sum_stratosphere
        use budget_global, only : sum_emission
#endif

        use emission,      only : emission_done

#ifdef MPI
        use mpi_const, only : my_real, mpi_sum, com_lev, root_t, myid
        use mpi_const, only : ntracetloc, lmloc, tracer_active, ierr, root_k
        use mpi_comm,  only : check_domain
        use mpi_comm,  only : gather_tracer_k, barrier, stopmpi
#endif

        implicit none

        ! input
        character(len=*),intent(in) :: file_name
        character(len=*),intent(in) :: msg
        integer, intent(out)        :: status

        character(len=*), parameter :: rname = mname//'/trace_end'

        ! local
        character(len=DUMMY_STR_LEN) :: frmt
        integer                  :: region, n, ind, io, imr, jmr, lmr
        integer                  :: sfstart, sfend, itr
        real                     :: mass(ntracet)

        ! start

#ifdef with_budgets
        if (apply_budget_global) then
            do region =1,nregions

                call diagbudg(region,3,status)    !add depdry
                IF_NOTOK_RETURN(status=1)
                call diagbudg(region,2,status)    !write output, free memory...
                IF_NOTOK_RETURN(status=1)

#ifdef MPI
                if ( myid == root_t ) then
#endif
                    do itr = 1, ntracet
                        mass(itr) = sum(mass_dat(region)%rm_t(1:im(region),1:jm(region),1:lm(region),itr))
                    end do

                    write(*,*) '--------------------------------------------------------'
                    write(*,'("Budget of tracers in region ", a, " (Kg)")') trim(region_name(region))
                    if (ntracet .le. 9) then
                        write(frmt, '(a,i1,"a20)")') "(a17,", ntracet
                    else
                        write(frmt, '(a,i2,"a20)")') "(a17,", ntracet
                    end if
                    write(*,frmt) "Tracer   : ", names(1:ntracet)
                    write(*,*) '--------------------------------------------------------'

                    if (ntracet .le. 9) then
                        write(frmt, '(a,i1,"es20.7)")') "(a17,", ntracet
                    else
                        write(frmt, '(a,i2,"es20.7)")') "(a17,", ntracet
                    end if

                    write(*,frmt) "Initial mass   : ", init_mass(region,:)
                    write(*,frmt) "emitted   : ", sum_emission(region,:)
                    write(*,frmt) "chemistry   : ", sum_chemistry(region,:)
                    write(*,frmt) "update_p   : ", sum_update(region,:)
                    write(*,frmt) "advection   : ", sum_advection(region,:)
                    write(*,frmt) "deposition   : ", sum_deposition(region,:)
                    write(*,frmt) "stratosphere   : ", sum_stratosphere(region,:)
                    write(*,frmt) "Final mass   : ", mass
                    write(*,frmt) "Budget error   : ", init_mass(region,:)  + sum_emission(region,:) + sum_chemistry(region,:) + sum_update(region,:) + sum_advection(region,:) &
                        + sum_deposition(region,:) + sum_stratosphere(region,:) - mass

                    write(*,*) '--------------------------------------------------------'
#ifdef MPI
                end if ! myid == root_t
#endif
            end do ! region
        endif   ! apply_budget_global
#endif

        status = 0

    end subroutine trace_end

    !  !=====================================================================================================
    !  !=====================================================================================================
    !
    !  subroutine read_diffusion( status )
    !
    !    !
    !    ! Read vertical diffusion from file and store in conv_dat%dkg
    !    !
    !    ! Output:
    !    !  o status      0 = ok
    !    !               -1 = file not available for at least one region
    !    !                1 = read error
    !    !
    !
    !    ! --- modules ------------------------------
    !
    !    use GO                     , only : GO_Timer_Start, GO_Timer_End
    !    use GO                     , only : pathsep
    !    use dims,                    only : itau, revert, region_name
    !    use global_data,             only : conv_dat
    !    use datetime,                only : tau2date
    !    use file_hdf,                only : Init, ReadData, Done, THdfFile, TSds
    !
    !    ! --- in/out ----------------------------------------------
    !
    !    integer, intent(out)             :: status
    !
    !    ! --- const -----------------------------------------------
    !
    !    character(len=*), parameter      :: rname = mname//', read_diffusion'
    !
    !    ! --- local -----------------------------------------------
    !
    !    character(len=MAX_FILENAME_LEN)  :: fname
    !    integer                          :: region
    !    integer, dimension(6)            :: idater
    !    logical                          :: exist
    !    type(THdfFile)                   :: hdf
    !    type(TSds)                       :: sds
    !
    !    ! --- begin -----------------------------------------------
    !
    !    ! start timing:
    !    call GO_Timer_Start( itim_read_diffusion, status )
    !    IF_NOTOK_RETURN(status=1)
    !
    !    if ( revert == 1 ) then
    !       call tau2date( itau, idater )
    !    else
    !       ! read field from 3 hours before
    !       call tau2date( itau - 3600*3, idater )
    !    end if
    !
    !    do region = 1, nregions
    !
    !       ! Create file name
    !       write(fname,'(2a,i4.4,a1,i2.2,a1,3a,i4,5i2.2,a)') &
    !            trim(diffusion_dir), pathsep, idater(1), pathsep, idater(2), pathsep, 'dkg_', &
    !            trim(region_name(region)), '_', idater, '.hdf'
    !
    !       inquire( file=trim(fname), exist=exist )
    !
    !       if ( .not. exist ) then
    !          write (gol,'(a,": dkg file not found for region",i2)') rname, region; call goPr
    !          write (gol,'(a,":   ",a)') rname, trim(fname); call goPr
    !          call GO_Timer_End( itim_read_diffusion, status )
    !          IF_NOTOK_RETURN(status=1)
    !          status = -1; return
    !       end if
    !
    !       ! info ...
    !       !write(*,'(a,": reading ", a)') rname, trim(fname)
    !
    !       ! Open file
    !       call Init( hdf, fname, 'read', status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       ! Read dkg
    !       call Init( sds, hdf, 'dkg', status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call ReadData( sds, conv_dat(region)%dkg, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call Done( sds, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       ! Read blh
    !       call Init( sds, hdf, 'blh', status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call ReadData( sds, conv_dat(region)%blh, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call Done( sds, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       ! Close file
    !       call Done( hdf, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !    end do
    !
    !    ! end timing:
    !    call GO_Timer_End( itim_read_diffusion, status )
    !    IF_NOTOK_RETURN(status=1)
    !
    !    ! ok
    !    status = 0
    !
    !  end subroutine read_diffusion
    !
    !  !=====================================================================================================
    !  !=====================================================================================================
    !
    !  subroutine write_diffusion( status )
    !
    !    !
    !    ! Write vertical diffusion to files (for all regions)
    !    !
    !    ! Output:
    !    !  o status      0 = ok, else not ok
    !    !
    !
    !    ! --- modules ------------------------------
    !
    !    use GO                     , only : GO_Timer_Start, GO_Timer_End
    !    use GO                     , only : pathsep
    !    use dims,                    only : itau, revert, region_name
    !    use global_data,             only : conv_dat
    !    use datetime,                only : tau2date
    !    use file_hdf,                only : Init, WriteData, Done, THdfFile, TSds
    !    use misctools,               only : check_dir
    !
    !    ! --- in/out ----------------------------------------------
    !
    !    integer, intent(out)             :: status
    !
    !    ! --- const -----------------------------------------------
    !
    !    character(len=*), parameter      :: rname = mname//', write_diffusion'
    !
    !    ! --- local -----------------------------------------------
    !
    !    character(len=MAX_FILENAME_LEN)  :: fname
    !    integer                          :: region
    !    integer, dimension(6)            :: idater
    !    type(THdfFile)                   :: hdf
    !    type(TSds)                       :: sds
    !
    !    ! --- begin -----------------------------------------------
    !
    !    ! start timing:
    !    call GO_Timer_Start( itim_write_diffusion, status )
    !    IF_NOTOK_RETURN(status=1)
    !
    !    if ( revert == 1 ) then
    !       call tau2date( itau, idater )
    !    else
    !       ! read field from 3 hours before
    !       call tau2date( itau - 3600*3, idater )
    !    end if
    !
    !    do region = 1, nregions
    !
    !       ! Create file name
    !       write(fname,'(2a,i4.4,a1,i2.2,a1,3a,i4,5i2.2,a)') &
    !            trim(diffusion_dir), pathsep, idater(1), pathsep, idater(2), pathsep, 'dkg_', &
    !            trim(region_name(region)), '_', idater, '.hdf'
    !
    !       ! info ...
    !       write (gol,'(a,": creating ", a)') rname, trim(fname); call goPr
    !       call check_dir(fname)
    !
    !       ! Create file
    !       call Init( hdf, fname, 'create', status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       ! Write dkg
    !       call Init( sds, hdf, 'dkg', shape(conv_dat(region)%dkg), 'real(8)', status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call WriteData( sds, conv_dat(region)%dkg, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call Done( sds, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !
    !       ! Write blh
    !       call Init( sds, hdf, 'blh', shape(conv_dat(region)%blh), 'real(8)', status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call WriteData( sds, conv_dat(region)%blh, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       call Done( sds, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !       ! Close file
    !       call Done( hdf, status )
    !       IF_NOTOK_RETURN(status=1)
    !
    !    end do
    !
    !    ! end timing:
    !    call GO_Timer_End( itim_write_diffusion, status )
    !    IF_NOTOK_RETURN(status=1)
    !
    !    status = 0
    !
    !  end subroutine write_diffusion
    !
    !  !=====================================================================================================
    !  !=====================================================================================================

end module sources_sinks
