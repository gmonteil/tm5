!###############################################################################
!
! Read routines for emission data.
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"

#define DIMS2D(my_array) lbound(my_array,1):ubound(my_array,1),lbound(my_array,2):ubound(my_array,2)

module Emission_Read_PyShell

    use GO, only : gol, goPr, goErr

    implicit none


    ! --- in/out -----------------------------------

    private

    public  ::  Read_Emissions_From_PyShell, read_dailycycle, close_dailycycle


    ! --- const ------------------------------------

    character(len=*), parameter        :: mname = 'Emission_Read_PyShell'


contains


    ! ================================================================================

    subroutine Read_Emissions_From_PyShell( em_cat_info, time_window, &
            ref_emissions_apri, status )

        ! --- modules ---------------------------------------------

        use GO                     , only : Readrc, Time_Window_Init, T_Time_Profile, Time_Profile_Done
        use GO                     , only : T_Time_Window, Time_Profile_Covers, Time_Profile_Init
        use GO                     , only : wrtgol, goPr, TDate, Pretty, NewDate
        use GO                     , only : operator(.in.), operator(.overlaps.)
        use GO_date,                 only : min, max
        use Dims                   , only : nregions, region_name
        use Dims                   , only : im, jm, lm
        use Global_Data            , only : rcF
        use Emission_Data          , only : T_Tracers_Info
        use TM5_Fields             , only : T_Fields_5D
        use chem_param,              only : ntracet
        use file_netcdf
        use os_specs               , only : MAX_FILENAME_LEN

        ! --- in/out ----------------------------------------------

        type(T_Tracers_Info), intent(inout) :: em_cat_info(nregions)
        type(T_Time_Window), intent(in)     :: time_window
        type(T_Fields_5D), intent(inout)    :: ref_emissions_apri(nregions)
        integer, intent(out)                :: status

        ! --- const -----------------------------------------------

        character(len=*), parameter        :: rname = mname//'/Read_Emissions_From_PyShell'

        ! --- local -----------------------------------------------

        character(len=MAX_FILENAME_LEN) :: PyShell_inputdir
        character(len=MAX_FILENAME_LEN) :: fname
        integer                         :: nc_id, group_id, region, i_cat, cgroup_id, tgroup_id
        real, allocatable               :: dummy_array(:,:,:)
        integer                         :: itr, n_cat

        integer(2), allocatable         :: dummy_int16(:,:)
        integer, allocatable            :: dummy_int(:,:)
        type(T_Time_Window)             :: emis_file_time_window, model_tw, file_tw
        type(T_Time_Profile)            :: emis_file_time_profile
        integer                         :: i1, i2, emis_file_ntime, itime, num_windows, etime
        type(TDate), allocatable        :: emis_file_tstart(:), emis_file_tend(:)
        integer, allocatable            :: overlap_index(:)
        logical                         :: zero_emis

        ! --- begin -----------------------------------------------

        call ReadRc( rcf, 'PyShell.em.filename', fname, status )
        IF_NOTOK_RETURN(status=1)

        write(gol,'(a, ": Reading ",a)') rname, trim(fname) ; call goPr

        nc_id = nc_open(trim(fname), 'r', status)
        IF_NOTOK_RETURN(status=1)

        do region = 1, nregions
            group_id = nc_get_group(nc_id, region_name(region))
            do itr = 1, ntracet
                tgroup_id = nc_get_group(group_id, em_cat_info(region)%tracer(itr)%name)
                n_cat = em_cat_info(region)%tracer(itr)%n_cat
                do i_cat = 1, n_cat
                    write(gol,'(a, ": reading emission for tracer ", a, ", category ", a)') rname, trim(em_cat_info(region)%tracer(itr)%name), &
                            em_cat_info(region)%tracer(itr)%cat(i_cat)%name ; call goPr
                    cgroup_id = nc_get_group(tgroup_id, em_cat_info(region)%tracer(itr)%cat(i_cat)%name)
                    dummy_array = nc_read_var(cgroup_id, 'emission',status)
                    IF_NOTOK_RETURN(status=1)

                    ! Check for NaN emission
                    if (any(dummy_array /= dummy_array)) then
                        write(0,'(a, " :: input emission array for tracer ", a, ", category ", a, " has NaN values")') &
                                rname, trim(em_cat_info(region)%tracer(itr)%name), trim(em_cat_info(region)%tracer(itr)%cat(i_cat)%name)
                        status=1
                        IF_NOTOK_RETURN(status=1)
                    end if

                    ! Check that the shape of the read-in emission is what is required
                    ! ref_emissions_apri(region)%cat(i_cat)%field is of dimension (nx,ny,nz,nt)
                    ! where nx = im(region), ny = jm(region), nz = 1, nt = em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile%n_period
                    if ((size(dummy_array,1) /= im(region)) .or. (size(dummy_array,2) /= jm(region))) then
                        write(0,'("      ", a, " :: expected horizontal resolution of emission input for region ", i1, " was ", i3, " x ", i3)') rname, region, im(region), jm(region)
                        write(0,'("      ", a, " :: instead, emission file gave ", i3, " x ", i3)') rname, size(dummy_array,1), size(dummy_array,2)
                        status = 1
                        IF_NOTOK_RETURN(status=1)
                    end if

                    if (em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile%n_period /= size(dummy_array,3)) then
                        ! This can happen, e.g., because the run job has been split up into multiple steps
                        ! So the emission file contains emissions for the entire run period, but this forward   .... (1)
                        ! run only knows about the split period.
                        emis_file_ntime = size(dummy_array, 3)
                        dummy_int16 = nc_read_var(cgroup_id, 'time_start', status)
                        IF_NOTOK_RETURN(status=1)
                        allocate(dummy_int(DIMS2D(dummy_int16)))
                        dummy_int(:,:) = dummy_int16(:,:)

                        allocate(emis_file_tstart(emis_file_ntime))
                        do itime = 1, emis_file_ntime
                            emis_file_tstart(itime) = NewDate(time6 = dummy_int(:,itime))
                        end do
                        deallocate(dummy_int, dummy_int16)

                        dummy_int16 = nc_read_var(cgroup_id, 'time_end', status)
                        IF_NOTOK_RETURN(status=1)
                        allocate(dummy_int(DIMS2D(dummy_int16)))
                        dummy_int(:,:) = dummy_int16(:,:)
                        allocate(emis_file_tend(emis_file_ntime))
                        do itime = 1, emis_file_ntime
                            emis_file_tend(itime) = NewDate(time6 = dummy_int(:,itime))
                        end do
                        deallocate(dummy_int, dummy_int16)

                        call Time_Window_Init(emis_file_time_window, emis_file_tstart(1), emis_file_tend(emis_file_ntime), status)
                        IF_NOTOK_RETURN(status=1)

                        ! For the situation in (1) to have happened, the emis_file_time_window must cover
                        ! em_cat_info(region)%cat(i_cat)%time_profile completely, so check for that
                        if (em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile .in. emis_file_time_window) then
                            ! Redefine em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile to match the time windows in the
                            ! emission file. Also redefine ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field to have the
                            ! correct shape.
                            call Time_Profile_Init(emis_file_time_profile, emis_file_tstart, emis_file_tend, status)
                            IF_NOTOK_RETURN(status=1)
                            ! Store the starting and ending times for the model
                            call Time_Window_Init(model_tw, em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile%t1, &
                                    em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile%t2, status)
                            IF_NOTOK_RETURN(status=1)
                            ! Destroy the old time profile
                            call Time_Profile_Done(em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile, status)
                            IF_NOTOK_RETURN(status=1)

                            ! Below, determine which of the time windows in the emission file have overlaps with the model time window.
                            ! We can reuse emis_file_tstart and emis_file_tend, because their purpose has been served. Suppose
                            !
                            ! model time window = .......... tm1 .......................................... tm2 ............
                            ! file time windows = tf1 .. tf2 ........ tf3 ..... tf4 ........... tf5 .. tf6 ...... tf7 .. tf8
                            !
                            ! In that case, num_windows = 5
                            ! emis_file_tstart(1:5) = [tm1, tf3, tf4, tf5, tf6]
                            ! emis_file_tend(1:5)   = [tf3, tf4, tf5, tf6, tm2]
                            ! overlap_index(1:5)    = [2,3,4,5,6]

                            allocate(overlap_index(emis_file_ntime))
                            overlap_index = -1
                            num_windows = 1
                            do itime = 1, emis_file_time_profile%n_period
                                call Time_Window_Init(file_tw, emis_file_time_profile%period(itime)%t1, &
                                        emis_file_time_profile%period(itime)%t2, status)
                                IF_NOTOK_RETURN(status=1)
                                if (file_tw .overlaps. model_tw) then
                                    emis_file_tstart(num_windows) = max(file_tw%t1, model_tw%t1)
                                    emis_file_tend(num_windows) = min(file_tw%t2, model_tw%t2)
                                    overlap_index(num_windows) = itime
                                    num_windows = num_windows + 1
                                end if
                            end do
                            num_windows = num_windows - 1
                            ! Now reinitialize em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile
                            call Time_Profile_Init(em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile, emis_file_tstart(1:num_windows), emis_file_tend(1:num_windows), status)
                            IF_NOTOK_RETURN(status=1)
                            ! Redefine ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field
                            deallocate(ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field)
                            allocate(ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field(im(region), jm(region), lm(region), num_windows))
                            ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field = 0.0

                            write(*,'("For tracer ", a, " in region ", a, ", category ", a, ", the time profile in the emission file does not match the one in the model")') &
                                    trim(em_cat_info(region)%tracer(itr)%name), trim(region_name(region)), trim(em_cat_info(region)%tracer(itr)%cat(i_cat)%name)
                            do itime = 1, num_windows
                                etime = overlap_index(itime)
                                ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field(:,:,1,itime) = dummy_array(:,:,etime)
                                ! Also print for debugging purposes
                                write(*,'("Using window ", i3, " (", a, " to ", a, ")", " to fill model time window ", &
                                        i3, " (", a, " to ", a, ")")'), etime, trim(Pretty(emis_file_time_profile%period(etime)%t1)), &
                                        trim(Pretty(emis_file_time_profile%period(etime)%t2)), itime, &
                                        trim(Pretty(em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile%period(itime)%t1)), &
                                        trim(Pretty(em_cat_info(region)%tracer(itr)%cat(i_cat)%time_profile%period(itime)%t2))
                            end do

                            deallocate(overlap_index)

                        else
                            write(0,'(a, " :: region = ", a, ", tracer = ", a, ", category = ", a)') rname, region, &
                                    trim(em_cat_info(region)%tracer(itr)%name), trim(em_cat_info(region)%tracer(itr)%cat(i_cat)%name)
                            call wrtgol(rname//" :: emission file goes from ", emis_file_time_window%t1, &
                                    " to ", emis_file_time_window%t2) ; call goPr
                            call wrtgol(rname//" :: simulation time window goes from ", time_window%t1, &
                                    " to ", time_window%t2) ; call goPr
                            status = 1
                            IF_NOTOK_RETURN(status=1)
                        end if

                        deallocate(emis_file_tstart, emis_file_tend)

                    else
                        ref_emissions_apri(region)%tracer(itr)%cat(i_cat)%field(:,:,1,:) = dummy_array(:,:,:)
                    end if

                    deallocate(dummy_array)

                    ! We should implement reading tf_diurnal from the emission file at a later stage
                    !em_cat_info(region)%tracer(itr)%cat(i_cat)%tf_diurnal = 1.0
                    !em_cat_info(region)%tracer(itr)%cat(i_cat)%use_tf_diurnal = .false.
                end do ! category
            end do ! tracer
        end do ! region

        status = 0

    end subroutine Read_Emissions_From_PyShell

    subroutine read_dailycycle(i_day, idate_local, region, itrac, status)

        use dims, only : nregions, region_name, sec_day
        use Emission_data, only : tracers_em_info
        use file_netcdf
        ! For writing status messages
        use dims,        only : kmain, itau
        use datetime,    only : tstamp
        use chem_param,  only : names, tracer_name_len
        use os_specs,    only : MAX_FILENAME_LEN, WRITE_STR_LEN, SHORT_STR_LEN

        implicit none

        !__IO___________________________________________________________________

        integer, intent(in)     :: i_day, idate_local(6), region, itrac
        integer, intent(out)    :: status

        !__CONSTANTS_____________________________________________________________

        character(len=*), parameter      :: rname = mname//'/read_dailycycle'

        !__LOCAL_VARIABLES______________________________________________________

        integer              :: nc_id, grp_id, i_cat, cat_grp_id
        character(len=MAX_FILENAME_LEN) :: fname, emis_indir
        logical              :: file_exist
        integer              :: n_cat, n_tstep
        integer              :: dailycycle_type
        character(len=tracer_name_len)  :: tracer_name
        character(len=SHORT_STR_LEN) :: dailycycle_pfx
        ! For writing status messages
        character(len=WRITE_STR_LEN) :: write_string

        !__START_SUBROUTINE______________________________________________________

        dailycycle_pfx = tracers_em_info(region)%tracer(itrac)%dailycycle%pfx
        dailycycle_type = tracers_em_info(region)%tracer(itrac)%dailycycle%dtype
        tracer_name = tracers_em_info(region)%tracer(itrac)%name
        n_cat = tracers_em_info(region)%tracer(itrac)%n_cat
        emis_indir = tracers_em_info(region)%tracer(itrac)%dailycycle%emis_indir

        ! The daily cycle files are called
        ! emis_indir/YYYY/MM/dailycycle_pfx_YYYYMMDD.nc4
        write(fname,'(a, "/", i4.4, "/", i2.2, "/", a, i4.4, 2i2.2, ".nc4")') &
                trim(emis_indir), idate_local(1), idate_local(2), trim(dailycycle_pfx), idate_local(1:3)

        inquire( file=fname, exist=file_exist )
        if ( file_exist ) then
            nc_id = nc_open(trim(fname), 'r', status)
            IF_NOTOK_RETURN(status=1)

            grp_id = nc_get_group(nc_id, region_name(region), status)
            IF_NOTOK_RETURN(status=1)
            do i_cat = 1,n_cat
                if (nc_grp_exists(grp_id, tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name)) then
                    cat_grp_id = nc_get_group(grp_id, tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name, status)
                    IF_NOTOK_RETURN(status=1)
                    select case (dailycycle_type)
                    case (0)
                        tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%scaling = nc_read_var(cat_grp_id, 'emission_scaling_factor', status)
                        IF_NOTOK_RETURN(status=1)
                    case (1)
                        tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%anomaly = nc_read_var(cat_grp_id, 'emission_anomaly', status)
                        IF_NOTOK_RETURN(status=1)
                    end select
                    n_tstep = nc_get_dim(cat_grp_id, 'timesteps', status)
                    IF_NOTOK_RETURN(status=1)
                    tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%dtime = sec_day/n_tstep
                    tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%apply = .true.
                    !                write(write_string, '(a, "/", a, "/", a, " : ", i2, " steps/day")') &
                    !                    trim(region_name(region)), trim(names(itrac)), trim(tracers_em_info(region)%tracer(itrac)%cat(i_cat)%name), n_tstep
                    !                call tstamp(kmain, itau, trim(write_string))
                else
                    tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%apply = .false.
                end if ! category exists in dailycycle file
            end do ! i_cat
            call nc_close(nc_id)
            tracers_em_info(region)%tracer(itrac)%dailycycle%day_opened(i_day) = .true.
        else
            write(0,'(a,a,a)') ' ERROR :: Diurnal cycle file ', trim(fname), ' does not exist'
            status=1
            IF_NOTOK_RETURN(status=1)
        end if

        status = 0

    end subroutine read_dailycycle

    subroutine close_dailycycle(region, itrac)

        use Emission_Data, only : tracers_em_info

        implicit none

        integer, intent(in) :: region, itrac

        character(len=*), parameter      :: rname = mname//'/close_dailycycle'
        integer :: n_cat, i_cat

        n_cat = tracers_em_info(region)%tracer(itrac)%n_cat

        do i_cat = 1, n_cat
            if (allocated(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%scaling)) &
                    deallocate(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%scaling)
            if (allocated(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%anomaly)) &
                    deallocate(tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%anomaly)
            tracers_em_info(region)%tracer(itrac)%dailycycle%cycle_cat(i_cat)%apply = .false.
        end do ! i_cat

    end subroutine close_dailycycle


end module Emission_Read_PyShell
