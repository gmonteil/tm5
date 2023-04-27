!###############################################################################
!
! Storage for data shared by 'Emission' and 'Adj_Emission' modules.
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module Emission_Data

    use GO, only : gol, goPr, goErr

    use GO        , only : T_Time_Window, T_Time_Profile
    use Dims      , only : nregions, region_name
    use chem_param, only : ntracet
    use TM5_Fields, only : T_Fields_5D
    use chem_param, only : tracer_name_len
    use os_specs,   only : MAX_FILENAME_LEN, SHORT_STR_LEN, MAX_RCKEY_LEN, MAX_RCVAL_LEN

    implicit none

    ! --- in/out -----------------------------------

    private

    public  ::  optim_emis_type
    public  ::  tracers_em_info
    public  ::  T_Tracers_info, ntracet, t_tracer_info
    public  ::  ref_emissions, ref_emissions_apri, adj_emissions
    public  ::  Emission_Data_Init, Emission_Data_Done, update_parent_emissions

    ! --- const ------------------------------------

    character(len=*), parameter        :: mname = 'Emission_Data'


    ! --- types ------------------------------------

    ! TIPP emission info:
    type T_TIPP_Info
        character(len=40)                       ::  category   ! TIPP subcategory name
        character(len=40)                       ::  source     ! TIPP subcategory source
        integer                                 ::  year       ! TIPP subcategory year
    end type T_TIPP_Info

    ! --- daily cycles------------------------------
    type dcycle_per_cat
        logical                             :: apply = .false. ! by default, do not apply a daily cycle
        ! within one tracer, one category may have a daily cycle while another may not, hence the
        ! above flag needs to be category specific, not tracer specific
        real, dimension(:,:,:), allocatable :: scaling ! n_tstep, n_lat, n_lon
        real, dimension(:,:,:), allocatable :: anomaly ! n_tstep, n_lat, n_lon
        integer                             :: dtime ! time resolution in seconds for this category
    end type dcycle_per_cat

    type diurnal_cycle
        logical                             :: apply = .false. ! by default, do not apply a diurnal cycle
        integer                             :: dtype
        character(len=SHORT_STR_LEN)        :: pfx
        character(len=MAX_FILENAME_LEN)     :: emis_indir
        logical, allocatable, dimension(:)  :: day_opened
        type(dcycle_per_cat), dimension(:), allocatable :: cycle_cat
    end type diurnal_cycle
    ! information on emission data (categories, months, unit)

    type T_Em_Cat_Info
        ! category name:
        character(len=40)                         ::  name         ! name of category
        ! temporal resolution:
        character(len=16)                         ::  treskey      ! 'monthly', ...
        type(T_Time_Profile)                      ::  time_profile ! time resolution, intervals, etc
        ! TIPP categories etc
        integer                                   ::  n_tipp       ! number of TIPP subcategories
        type(T_TIPP_Info), pointer                ::  tipp(:)
        !
        ! * 4D-var stuff
        !     Should be moved to 4D-var modules!
        !     For the moment kept here since all numbers are read from rcfile
        !
        ! variance:
!        real                                      ::  error       ! standard deviation (%)
        ! horizontal correlation:
!        real                                      ::  corlen      ! spatial correl length (km)
!        character(len=1)                          ::  corchoice   ! type of spatial corr
        ! temporal correlation:
!        real                                      ::  tau          ! temporal correl length (months)
!        character(len=1)                          ::  tauchoice    ! type of temporal corr
        !real,dimension(24)                        ::  tf_diurnal = 1.0 ! Read this in from read_emission_pyshell if needed
        !logical                                   ::  use_tf_diurnal = .false.
    end type T_Em_Cat_Info

    ! Information on tracer (name, emissions, emission categories)
    type T_tracer_info
        character(len=tracer_name_len)  :: name          ! tracer name
        integer                         :: n_cat         ! number of categories per region
        type(T_Em_Cat_Info),allocatable :: cat(:)
        type(diurnal_cycle)             :: dailycycle
    end type T_Tracer_info

    type T_tracers_info
        type(T_tracer_info),allocatable :: tracer(:)
    end type T_tracers_info

    ! --- var --------------------------------------

    type(T_tracers_info), allocatable      :: tracers_em_info(:)     ! contains information on all tracers
    integer                                :: optim_emis_type
    type(T_fields_5D), allocatable         :: ref_emissions(:)       ! () reference emission field ; used for linearization
    type(T_fieldS_5D), allocatable, target :: ref_emissions_apri(:)  ! (n_cat) reference emission field (a priori)

    ! adjoint model variables:
    type(T_Fields_5D), pointer            ::  adj_emissions(:)   ! (n_cat) dJ_obs/d(emissions) in adjoint run

    ! The hierarchy is region -> tracer -> category
    ! The first two levels are interchangeable, but since André already coded it this way, let's stick to it

contains


    ! ================================================================================


    subroutine Emission_Data_Init( status )

        use GO         , only : Get
        use GO         , only : TrcFile, ReadRc
        use TM5_Fields , only : Fields_4D_Init
        use global_data, only : rcF
        use datetime   , only : time_window
        use chem_param , only : names
        use go, only : time_profile_init

        ! --- in/out ---------------------------------

        integer, intent(out)                ::  status

        ! --- const ---------------------------------

        character(len=*), parameter      :: rname = mname//'/Emission_Data_Init'

        ! --- local ----------------------------------

        integer                             ::  i_cat, i_tra
        character(len=MAX_RCKEY_LEN)        ::  rckey
        character(len=MAX_RCVAL_LEN)        ::  line

        !    integer                             ::  year1, year2
        integer                             ::  year
        integer                             ::  np, region, ncat
        integer, dimension(:), allocatable  ::  nt

        character(20), dimension(:), allocatable :: catnames
        character(len=200) :: tmpstr


        ! --- begin ----------------------------------

        ! optimize emissions ? this should be known in order
        ! to determine how many super-categories have to be allocated:

        ! number of emissions categories:
        ! categories and prior information

        ! info ...
        write (gol,'("emission data init ...")'); call goPr

        allocate(tracers_em_info(nregions))

        ! emission category info: different for each tracer and region
        do region=1,nregions
            allocate(tracers_em_info(region)%tracer(ntracet))
            do i_tra=1,ntracet
                tracers_em_info(region)%tracer(i_tra)%name = names(i_tra)

                ! read number of categories from rc file.lines should have the structure:
                call ReadRc( rcF, 'emissions.'//trim(names(i_tra))//'.'//trim(region_name(region))//'.ncats', ncat, status)
                IF_NOTOK_RETURN(status=1)
                tracers_em_info(region)%tracer(i_tra)%n_cat = ncat

                ! emission category info: different for different regions and tracers!
                allocate(tracers_em_info(region)%tracer(i_tra)%cat(ncat))

                ! Read name of categories (comma-separated):
                allocate(catnames(ncat))
                call readrc(rcf, 'emissions.'//trim(names(i_tra))//'.'//trim(region_name(region))//'.categories', tmpstr, status)
                IF_NOTOK_RETURN(status=1)
                read(tmpstr, *) catnames

                ! loop over var4d categories: note MK: we allow for different categories in different regions
                do i_cat=1, ncat
!                    if (i_cat < 10) then
!                        write (rckey,'(".",a10,".category",i1)') region_name(region), i_cat
!                    else
!                        write (rckey,'(".",a10,".category",i2)') region_name(region), i_cat
!                    endif
!                    rckey = "emission."//trim(names(i_tra))//rckey
                    tracers_em_info(region)%tracer(i_tra)%cat(i_cat)%name = catnames(i_cat)

                    ! Read the key "emission.{tracer}.{region}.{category}.treskey.
                    ! Should be either "daily" or "monthly"
                    call readrc(&
                            rcf, &
                            'emissions.'//trim(names(i_tra))//'.'//region_name(region)//'.'//trim(catnames(i_cat)), &
                            tracers_em_info(region)%tracer(i_tra)%cat(i_cat)%treskey, &
                            status)
                    IF_NOTOK_RETURN(status=1)

                    call time_profile_init( &
                        tracers_em_info(region)%tracer(i_tra)%cat(i_cat)%time_profile, &
                        time_window, &
                        tracers_em_info(region)%tracer(i_tra)%cat(i_cat)%treskey, &
                        status)

!                    call ReadRc( rcF, trim(rckey), line, status)
!                    IF_NOTOK_RETURN(status=1)
!                    call Cat_Info_Init(tracers_em_info(region)%tracer(i_tra)%cat(i_cat), line, status)
!                    IF_NOTOK_RETURN(status=1)
                end do ! category
                deallocate(catnames)
            end do ! tracer
        end do ! region

        ! which optim type ?
        call ReadRc( rcF, 'var4d.optim_emis.type', optim_emis_type, status )
        IF_NOTOK_RETURN(status=1)

        ! Type of emission optimization
        select case ( optim_emis_type )
        case ( 1, 2, 3 )
            write (gol,'(a,": Type of emissions optimization is: ",i1)') rname, optim_emis_type; call goPr
        case default
            write(gol,'(a,": Wrong emission optimization type: ", i3)') rname, optim_emis_type; call goErr
            TRACEBACK; status=1; return
        end select

        ! emission arrays:
        allocate( ref_emissions(nregions) )  ! reference emission field
        allocate( ref_emissions_apri(nregions) )  ! reference emission field (a priori)

        ! loop over var4d categories:
        do region = 1, nregions
            allocate(ref_emissions(region)%tracer(ntracet))
            allocate(ref_emissions_apri(region)%tracer(ntracet))
            do i_tra = 1, ntracet
                ! allocate time series of surface fields:
                ncat = tracers_em_info(region)%tracer(i_tra)%n_cat
                allocate( nt(ncat) )

                do i_cat = 1, ncat
                    nt(i_cat) = tracers_em_info(region)%tracer(i_tra)%cat(i_cat)%time_profile%n_period
                end do
                write(gol, '(a, ": Allocating emissions for region ", i1)') rname, region ; call goPr
                call Fields_4D_Init( ref_emissions(region)%tracer(i_tra), region, .false., ncat, nt, status )
                IF_NOTOK_RETURN(status=1)
                call Fields_4D_Init( ref_emissions_apri(region)%tracer(i_tra), region, .false., ncat, nt, status )
                IF_NOTOK_RETURN(status=1)
                deallocate(nt)

                ! init to zero:
                do i_cat=1,ncat
                    ref_emissions(region)%tracer(i_tra)%cat(i_cat)%field = 0
                    ref_emissions_apri(region)%tracer(i_tra)%cat(i_cat)%field = 0
                end do ! category
            end do ! tracer
        end do ! region

        ! *

        ! extract year range from time window:
        !    call Get( time_window%t1, year=year1 )
        !    call Get( time_window%t2, year=year2 )

        status = 0

    end subroutine Emission_Data_Init


    ! ***


    subroutine Emission_Data_Done( status )

        use TM5_Fields, only : Fields_4D_Done

        ! --- in/out ---------------------------------

        integer, intent(out)                ::  status

        ! --- const ---------------------------------

        character(len=*), parameter      :: rname = mname//'/Emission_Data_Done'

        ! --- local ----------------------------------

        integer       ::  year
        integer       ::  i_cat, region, itr, ncat

        ! --- begin ----------------------------------

        ! emissions
        do region = 1, nregions
            do itr = 1,ntracet
                ncat = tracers_em_info(region)%tracer(itr)%n_cat
                call Fields_4D_Done( ref_emissions(region)%tracer(itr), ncat, status)
                IF_NOTOK_RETURN(status=1)
                call Fields_4D_Done( ref_emissions_apri(region)%tracer(itr), ncat, status)
                IF_NOTOK_RETURN(status=1)
            end do ! tracer
            deallocate(ref_emissions(region)%tracer)
            deallocate(ref_emissions_apri(region)%tracer)
        end do ! region

        ! clear:
        deallocate( ref_emissions      )
        deallocate( ref_emissions_apri )

        do region = 1, nregions
            do itr = 1, ntracet
                deallocate(tracers_em_info(region)%tracer(itr)%cat)
            end do
            deallocate(tracers_em_info(region)%tracer)
        end do

        deallocate(tracers_em_info)

        ! ok
        status = 0

    end subroutine Emission_Data_Done


    subroutine Cat_Info_Done( eci, status )

        use GO, only : Time_Profile_Done

        ! --- in/out ----------------------------------------------

        type(T_Em_Cat_Info), intent(inout)  ::  eci
        integer, intent(out)                ::  status

        ! --- const -----------------------------------------------

        character(len=*), parameter         ::  rname = mname//'/Cat_Info_Done'

        ! --- local -----------------------------------------------

        integer                             ::  j

        ! --- begin -----------------------------------------------

        ! done with time profile:
        call Time_Profile_Done( eci%time_profile, status )
        IF_NOTOK_RETURN(status=1)

        ! clear tipp info:
        deallocate( eci%tipp )

        ! ok:
        status = 0

    end subroutine Cat_Info_Done


    ! ============================================================================


    subroutine TIPP_Info_Init( tipp_info, tipp_string, status )

        use GO                    , only : goReadFromLine

        ! --- in/out ----------------------------------------------

        type(T_TIPP_Info), intent(out)      ::  tipp_info
        character(len=*), intent(inout)     ::  tipp_string
        integer, intent(out)                ::  status

        ! --- const -----------------------------------------------

        character(len=*), parameter         ::  rname = mname//'/TIPP_Info_Init'

        ! --- local -----------------------------------------------

        ! --- begin -----------------------------------------------

        ! tipp_string format  :  category[-[source][-[year]]]

        call goReadFromLine( tipp_string, tipp_info%category, status, sep='-' )
        IF_ERROR_RETURN(status=1)
        call goReadFromLine( tipp_string, tipp_info%source  , status, sep='-' )
        IF_ERROR_RETURN(status=1)
        call goReadFromLine( tipp_string, tipp_info%year    , status, sep='-', default=0000 )
        IF_ERROR_RETURN(status=1)

        ! ok:
        status = 0

    end subroutine TIPP_Info_Init


    ! ============================================================================

    ! copied from T38 var4d_state_emis.F90

    subroutine update_parent_emissions( emis )

        use dims,                    only : nregions
        use dims,                    only : isr, ier, jsr, jer, xref, yref
        use dims,                    only : ibeg, iend, jbeg, jend
        use dims,                    only : parent
        use global_data,             only : mass_dat
        use toolbox,                 only : escape_tm
        use TM5_Fields             , only : T_Fields_5D

        !__IO___________________________________________________________________

        type(T_Fields_5D), intent(inout), target   ::  emis(:)  ! (nregions)

        !__LOCAL_VARIABLES______________________________________________________

        integer                            :: region

        real, dimension(:,:,:,:), pointer  :: e
        real, dimension(:,:,:,:), pointer  :: ep
        real, dimension(:,:,:,:), allocatable  :: to_parent

        integer                            :: i_period, n_period
        integer                            :: i_cat, itr
        integer                            :: i, j
        integer                            :: ip, jp
        integer                            :: my_parent
        integer                            :: xref_, yref_
        integer                            :: imp, jmp

        !__START_SUBROUTINE______________________________________________________

        ! JFM: change 13/04/2007
        ! Emissions at region edges are taken from parent grid cell.
        ! Therefore, restrict updating of parent emissions here to
        ! the inner zoom region (i.e. without edges).
        ! In previous version, the edges were also updated.

        ! loop over regions:
        do region = nregions, 2, -1

            ! determine parent
            my_parent = parent(region)

            xref_ = xref(region)/xref(my_parent)
            yref_ = yref(region)/yref(my_parent)

            imp = (ier(region)-isr(region)+1)/xref_
            jmp = (jer(region)-jsr(region)+1)/yref_

            ! check calculated imp, jmp, lmp
            ! stop when zooming over dateline
            if(ibeg(region) >= iend(region)) call escape_tm('stopped in update_parent_emissions')
            if(imp /= iend(region)-ibeg(region)-1) call escape_tm('stopped in update_parent_emissions')
            if(jmp /= jend(region)-jbeg(region)-1) call escape_tm('stopped in update_parent_emissions')

            ! loop over tracers
            do itr = 1, ntracet

                ! loop over categories:
                do i_cat = 1, tracers_em_info(region)%tracer(itr)%n_cat

                    n_period = tracers_em_info(region)%tracer(itr)%cat(i_cat)%time_profile%n_period
                    allocate( to_parent(imp,jmp,1,n_period) )

                    e  => emis(region)%tracer(itr)%cat(i_cat)%field
                    ep => emis(my_parent)%tracer(itr)%cat(i_cat)%field

                    to_parent = 0.0
                    do i = isr(region), ier(region)
                        ip = 1 + (i-isr(region))/xref_
                        if ( ip < 1 .or. ip > imp ) then
                            write(*,*) 'update_parent_emissions: ip out of range'
                            write(*,*) 'ip, imp = ', ip, imp
                            call escape_tm('stopped in update_parent_emissions')
                        end if
                        do j = jsr(region), jer(region)
                            jp = 1 + (j-jsr(region))/yref_
                            if ( jp < 1 .or. jp > jmp ) then
                                write(*,*) 'update_parent_emissions: jp out of range'
                                write(*,*) 'jp, jmp = ', jp, jmp
                                call escape_tm('stopped in update_parent_emissions')
                            end if
                            ! add contribution:
                            to_parent(ip,jp,:,:) = to_parent(ip,jp,:,:) + e(i,j,:,:)
                        end do
                    end do

                    ! insert in parrent:
                    ep(ibeg(region)+1:iend(region)-1,jbeg(region)+1:jend(region)-1,:,:) = to_parent

                    deallocate( to_parent )

                    nullify( e )
                    nullify( ep )

                end do  ! categories
            end do ! tracers
        end do   ! regions

    end subroutine update_parent_emissions

end module Emission_Data

