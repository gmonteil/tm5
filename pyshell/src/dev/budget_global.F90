!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module budget_global
    !
    ! this budget routine defines the global budgets
    ! for instance, the budget in certain lattitude bands (SH, NH)
    ! and heights (BL, FT, STRAT)
    ! Apart from that, the budget in the finest zoom region can be defined
    ! There is a difficulty with the 'INTERFACE' regions. Here, the
    ! emissions, conv, chemistry are performed in the parent resolutions.
    ! A straightforward separation between advection and these
    ! processes becomes therefore problamatic. Therefore, the advection
    ! budget is not calculated for now.
    !
    ! The budget is defined according to the 'coarsest' resolution (e.g. 6x4)
    ! The zoom regions will be defined according to the defined regions.
    !
    ! 07/2006: removed dry deposition budget: is now in dry_deposition

    use dims,       only : nregions, im, jm, lm
    use dims,       only : dx, dy, xref, xbeg, xend, yref, ybeg, yend
    use dims,       only : itau, idate, idatei, idatee
    use dims,       only : region_name
    use chem_param, only : names, ra, nmark, njnum, jnam, nreac
    use chem_param, only : nreacw, rwnam, ntrace, ntracet
    use chem_param, only : marknam, ratnam
    use go,         only : gol, goErr
    use GO        , only : T_Time_Profile
    use os_specs,   only : MAX_FILENAME_LEN, MAX_RCVAL_LEN

    implicit none

    ! --- interface ---

    private

    ! public routines
    public :: ini_zoneg, diagbudg, budget_transportg, init_budget

    ! public types
    public :: buddep_data, budg_data, budflux_data

    ! public variables
    public :: init_mass, buddep_dat, budg_dat, nzon_vg, nbudg, nbud_vg
    public :: sum_deposition, sum_wet, sum_chemistry, sum_emission
    public :: sum_update, sum_advection, sum_stratosphere
    public :: iflux1, iflux2, jflux1, jflux2, lflux1, lflux2
    public :: budemig, budstratg, budrjg, budrrg, budrwg, budmarkg
    public :: budget_flux, apply_budget_global
    public :: budget_time_profile

    ! --- const ------------------------------

    character(len=*), parameter      :: mname = 'budget_global'


    ! --- module variables ---

    logical               :: apply_budget_global = .true.

    !  integer,parameter     :: nbudg=7, nbud_vg=4
    integer,parameter     :: nbudg=3, nbud_vg=4

    !  deposition data on gridbasis
    type buddep_data
        real,dimension(:,:,:,:),pointer    :: lsp     ! im,jm,nbud_vg,ntrace
        real,dimension(:,:,:,:),pointer    :: cp      ! im,jm,nbud_vg,ntrace
        real,dimension(:,:,:,:),pointer    :: cond    ! im,jm,nbud_vg,ntrace
        real,dimension(:,:,:),pointer      :: dry     ! im,jm,ntrace
    end type buddep_data

    type(buddep_data),dimension(nregions),target  :: buddep_dat, buddep_dat_all

    ! real,dimension(nbudg,nbud_vg,ntracet)        :: rmoldg,rmnewg
    real, dimension(:, :, :), allocatable :: rmoldg, rmnewg

    ! budget for tracer 1 for output checking.......
    real, dimension(:, :), allocatable :: init_mass, sum_deposition, sum_wet, sum_chemistry, sum_update, sum_advection, sum_stratosphere, sum_emission

    !    budrrg  : accumulator for chemical reaction fluxes, separated by tracer
    !    budrjg  : accumulator for photolysis reactions
    !    buddryg : accumulator for dry deposition losses
    !    budwetg : accumulator for wet deposition losses
    !    budemig : accumulator for emissions entering calculations
    !    budstratg : accumulator for stratospheric losses and gains
    !    budconcg  : average concentrations
    !    budchngg  : change of concentration
    !    budinig   : initial concentration at begin of run
    !    budconvg : convection in and out of region per tracer
    !
    !  real, dimension(nbudg,nbud_vg,nreac)  :: budrrg
    real, dimension(:,:,:,:,:), allocatable  :: budrrg
    real, dimension(:, :, :), allocatable :: budrrg_all, budrjg, budrjg_all
    real, dimension(:, :, :), allocatable :: budwetg_lsp, budwetg_cp, budwetg_lsp_all, budwetg_cp_all
    real, dimension(:,:,:,:), allocatable  :: budemig
    real, dimension(:, :, :), allocatable :: budstratg, budconcg, budconcg_all, budinig, budemig_all, budstratg_all
    real, dimension(:, :, :), allocatable :: budinig_all, budchngg, budchngg_all, budmarkg, budmarkg_all
    real, dimension(:, :, :), allocatable :: budrwg, budrwg_all
    real, dimension(:, :, :), allocatable :: budadvxg, budadvyg, budadvzg
    real, dimension(:, :, :), allocatable :: budadvxg_all, budadvyg_all, budadvzg_all
    real, dimension(:, :, :), allocatable :: budconvg, budconvg_all, budvdifg, budvdifg_all

    character(len=20), dimension(nbudg)   :: zone_nameg
    character(len=20), dimension(nbud_vg) :: zone_namvg

    type budg_data
        ! defines the horizontal region definition
        integer,dimension(:,:),pointer      :: nzong
    end type budg_data

    type(budg_data),dimension(nregions),target    :: budg_dat

    integer,dimension(lm(1))          :: nzon_vg
    integer                           :: lglob !WP! for offset on level domain

    ! added dec 2004: fluxes of all species through some predefined levels:
    type budflux_data
        real,dimension(:,:,:),pointer      :: flux_x1   ! jm,lm,ntracet
        real,dimension(:,:,:),pointer      :: flux_x2   ! jm,lm,ntracet
        real,dimension(:,:,:),pointer      :: flux_y1   ! im,lm,ntracet
        real,dimension(:,:,:),pointer      :: flux_y2   ! im,lm,ntracet
        real,dimension(:,:,:),pointer      :: flux_z1   ! im,jm,ntracet
        real,dimension(:,:,:),pointer      :: flux_z2   ! im,jm,ntracet
    end type budflux_data

    type(budflux_data),dimension(nregions),target  :: budget_flux, budget_flux_all

    integer, dimension(nregions)          :: iflux1, iflux2  ! calculate flux flowing in from west (2x)
    integer, dimension(nregions)          :: jflux1, jflux2  ! calculate flux flowing in from south (2x)
    integer, dimension(nregions)          :: lflux1, lflux2  ! calculate flux flowing in from below (2x)

    character(len=MAX_RCVAL_LEN) ::  budget_subdir

    type(T_Time_Profile)      ::  budget_time_profile


contains


    subroutine init_budget

        allocate(rmoldg(nbudg, nbud_vg, ntracet), rmnewg(nbudg, nbud_vg, ntracet))
        allocate(init_mass(nregions, ntracet), sum_deposition(nregions, ntracet), sum_wet(nregions, ntracet))
        allocate(sum_chemistry(nregions, ntracet), sum_update(nregions, ntracet), sum_advection(nregions, ntracet))
        allocate(sum_stratosphere(nregions, ntracet), sum_emission(nregions, ntracet))
        init_mass(:, :) = 0
        sum_deposition(:, :) = 0
        sum_wet(:, :) = 0
        sum_chemistry(:, :) = 0
        sum_update(:, :) = 0
        sum_advection(:, :) = 0
        sum_stratosphere(:, :) = 0
        sum_emission(:, :) = 0
        allocate(budwetg_lsp(nbudg, nbud_vg, ntrace), budwetg_lsp_all(nbudg, nbud_vg, ntrace))
        allocate(budwetg_cp(nbudg, nbud_vg, ntrace), budwetg_cp_all(nbudg, nbud_vg, ntrace))
        allocate(budstratg(nbudg, nbud_vg, ntrace), budstratg_all(nbudg, nbud_vg, ntrace))
        allocate(budconcg(nbudg, nbud_vg, ntrace), budconcg_all(nbudg, nbud_vg, ntrace))
        allocate(budinig(nbudg, nbud_vg, ntrace), budinig_all(nbudg, nbud_vg, ntrace))
        allocate(budemig_all(nbudg, nbud_vg, ntrace))
        allocate(budrrg_all(nbudg, nbud_vg, nreac))
        allocate(budrjg(nbudg, nbud_vg, njnum), budrjg_all(nbudg, nbud_vg, njnum))
        allocate(budchngg(nbudg, nbud_vg, ntrace), budchngg_all(nbudg, nbud_vg, ntrace))
        budchngg(:, :, :) = 0
        budchngg_all(:, :, :) = 0
        allocate(budmarkg(nbudg, nbud_vg, nmark), budmarkg_all(nbudg, nbud_vg, nmark))
        allocate(budrwg(nbudg, nbud_vg, nreacw), budrwg_all(nbudg, nbud_vg, nreacw))
        allocate(budadvxg(nbudg, nbud_vg, ntracet), budadvxg_all(nbudg, nbud_vg, ntracet))
        allocate(budadvyg(nbudg, nbud_vg, ntracet), budadvyg_all(nbudg, nbud_vg, ntracet))
        allocate(budadvzg(nbudg, nbud_vg, ntracet), budadvzg_all(nbudg, nbud_vg, ntracet))
        allocate(budconvg(nbudg, nbud_vg, ntracet), budconvg_all(nbudg, nbud_vg, ntracet))
        allocate(budvdifg(nbudg, nbud_vg, ntracet), budvdifg_all(nbudg, nbud_vg, ntracet))

    end subroutine init_budget

    subroutine ini_zoneg( status )
        !
        ! user defined structure of budget zones
        !
        use GO         , only : ReadRc
        use GO         , only : pathsep
        use GO         , only : Time_Profile_Init
        use global_data, only : outdir, rcF
        use datetime   , only : time_window
        !    use io_hdf     , only : io_write, DFACC_CREATE, DFNT_INT32
        use file_netcdf
        use toolbox    , only : lvlpress, escape_tm
        use ParTools   , only : root,myid
        use misctools  , only : check_dir

        implicit none

        ! --- in/out ----------------------------------------------

        integer, intent(out)   ::  status

        ! --- const -----------------------------------------------

        character(len=*), parameter        :: rname = mname//'/ini_zoneg'

        ! --- local -----------------------------------------------

        integer :: i,j,n, nix,niy,region
        integer :: io,sfstart,sfend,istat,sfsnatt
        integer :: ilvl1, ilvl2, lats, ilvl3
        real    :: psurf, plev, rlat, rlon
        integer :: dyr, dxr

        character(len=MAX_FILENAME_LEN) ::  fname
        character(len=MAX_RCVAL_LEN)    ::  budget_tres

        ! --- begin -----------------------------------------------

        call ReadRc(rcF, 'budget.global', apply_budget_global,status)
        IF_NOTOK_RETURN(status=1)
        if (.not. apply_budget_global) return

        call ReadRc(rcF, 'budget.output.subdir', budget_subdir, status)
        IF_NOTOK_RETURN(status=1)

        call ReadRc(rcF, 'budget.global.timechunk', budget_tres, status)
        IF_NOTOK_RETURN(status=1)

        ! init budget time profile:
        call Time_Profile_Init( budget_time_profile, time_window, trim(budget_tres), status )
        IF_NOTOK_RETURN(status=1)

        region = 1
        !
        !
        zone_namvg(1)='BL'
        zone_namvg(2)='FT'
        zone_namvg(3)='STRAT'
        !
        zone_nameg(1)='SH 90-30'
        zone_nameg(2)='30S-30N '
        zone_nameg(3)='NH 90-30'
        !    zone_nameg(4)='Europe'
        !    zone_nameg(5)='N-america'
        !    zone_nameg(6)='Asia'
        !    zone_nameg(7)='Africa'

        psurf = 98400.0
        plev  = 85000.0
        ilvl1 = lvlpress(region,plev, psurf)
        plev  = 50000.0
        ilvl2 = lvlpress(region,plev, psurf)
        plev  = 10000.0
        ilvl3 = lvlpress(region,plev, psurf)
        nzon_vg(1:ilvl1) = 1 !boundary layer
        nzon_vg(ilvl1+1:ilvl2)= 2 !850-->500 hPa
        nzon_vg(ilvl2+1:ilvl3)= 3 !500-->100 hPa
        nzon_vg(ilvl3+1:lm(1))= 4 !100 hPa-->TOA
        lflux1(:) = ilvl1+1
        lflux2(:) = ilvl2+1



        ! allocate arrays for horizontal decomposition of the regions:

        do n=1,nregions
            allocate( budg_dat(n)%nzong(im(n),jm(n)))
            allocate( buddep_dat(n)%lsp(im(n),jm(n),nbud_vg,ntrace))
            allocate( buddep_dat(n)% cp(im(n),jm(n),nbud_vg,ntrace))
            allocate( buddep_dat(n)% cond(im(n),jm(n),nbud_vg,ntrace))
            allocate( buddep_dat(n)%dry(im(n),jm(n),ntrace))
            allocate( buddep_dat_all(n)%lsp(im(n),jm(n),nbud_vg,ntrace))
            allocate( buddep_dat_all(n)% cp(im(n),jm(n),nbud_vg,ntrace))
            allocate( buddep_dat_all(n)% cond(im(n),jm(n),nbud_vg,ntrace))
            allocate( buddep_dat_all(n)%dry(im(n),jm(n),ntrace))
            ! allocate flux data:
            allocate( budget_flux(n)%flux_x1(jm(n),lm(n),ntracet))
            allocate( budget_flux(n)%flux_y1(im(n),lm(n),ntracet))
            allocate( budget_flux(n)%flux_z1(im(n),jm(n),ntracet))
            allocate( budget_flux(n)%flux_x2(jm(n),lm(n),ntracet))
            allocate( budget_flux(n)%flux_y2(im(n),lm(n),ntracet))
            allocate( budget_flux(n)%flux_z2(im(n),jm(n),ntracet))
            allocate( budget_flux_all(n)%flux_x1(jm(n),lm(n),ntracet))
            allocate( budget_flux_all(n)%flux_y1(im(n),lm(n),ntracet))
            allocate( budget_flux_all(n)%flux_z1(im(n),jm(n),ntracet))
            allocate( budget_flux_all(n)%flux_x2(jm(n),lm(n),ntracet))
            allocate( budget_flux_all(n)%flux_y2(im(n),lm(n),ntracet))
            allocate( budget_flux_all(n)%flux_z2(im(n),jm(n),ntracet))
            budget_flux(n)%flux_x1(:,:,:) = 0.0
            budget_flux(n)%flux_y1(:,:,:) = 0.0
            budget_flux(n)%flux_z1(:,:,:) = 0.0
            budget_flux(n)%flux_x2(:,:,:) = 0.0
            budget_flux(n)%flux_y2(:,:,:) = 0.0
            budget_flux(n)%flux_z2(:,:,:) = 0.0
        end do

        allocate( budemig(nbudg,nbud_vg,ntrace,budget_time_profile%n_period))
        allocate( budrrg(nbudg,nbud_vg,ntrace,nreac,budget_time_profile%n_period))

        ! copy zone information to all regions...

        do n = 1,nregions
            dyr = nint(dy/yref(n))
            dxr = nint(dx/xref(n))
            do j=1,jm(n)
                lats = ybeg(n) + (j-1)*dyr
                select case(lats)
                case(-90:-31)
                    budg_dat(n)%nzong(:,j) = 3  !SH 90-30
                case(-30:29)
                    budg_dat(n)%nzong(:,j) = 2  !-30->30
                case(30:90)
                    budg_dat(n)%nzong(:,j) = 1  !NH 30-90
                case default
                    call escape_tm('ini_zoneg: Wrong latitude in budget routine')
                end select
            end do

            !       do j=1,jm(n)
            !          rlat = ybeg(n) + (float(j)-0.5)*dyr
            !          if(  nint(ybeg(n) + float(j-1)*dyr) == 34) jflux1(n) = j
            !          if(  nint(ybeg(n) + float(j-1)*dyr) == 62) jflux2(n) = j
            !          !if(  nint(ybeg(n) + float(j)*dyr) == 62)   jflux2(n) = j !CMK corrected
            !          do i=1,im(n)
            !             rlon = xbeg(n) + (float(i)-0.5)*dxr
            !             if(  nint(xbeg(n) + float(i-1)*dxr) == -12) iflux1(n) = i
            !             if(  nint(xbeg(n) + float(i-1)*dxr) == 36)  iflux2(n) = i   ! corrected CMK
            !             !cmk if(  nint(xbeg(n) + float(i)*dxr) == 36)  iflux2(n) = i   ! corrected CMK
            !             if(rlon > -12. .and. rlon < 36. .and. rlat > 34. .and. rlat < 62.)  then
            !                budg_dat(n)%nzong(i,j) = 4   ! EUROPE
            !             elseif(rlon > -126. .and. rlon < -66. .and. rlat > 22. .and. rlat < 62.)  then
            !                budg_dat(n)%nzong(i,j) = 5   ! NAM
            !             elseif(rlon > 72. .and. rlon < 150. .and. rlat > -10. .and. rlat < 50.)  then
            !                budg_dat(n)%nzong(i,j) = 6   ! Asia
            !             elseif(rlon > -18. .and. rlon < 48. .and. rlat > -34. .and. rlat < 34.)  then
            !                budg_dat(n)%nzong(i,j) = 7   ! Africa
            !             endif
            !          enddo
            !       enddo

            if ( myid == root ) then
                write (fname,'(4a,"regionsg_",a,".nc")') trim(outdir), pathsep, trim(budget_subdir), pathsep, region_name(n)
                call check_dir(trim(fname))
                io = nc_open(trim(fname), 'c', status)
                IF_NOTOK_RETURN(status=1)

                call nc_create_dim(io, 'im', im(n))
                call nc_create_dim(io, 'jm', jm(n))
                call nc_create_dim(io, 'lm', lm(n))
                call nc_dump_var(io, 'region_global', (/'im', 'jm'/), budg_dat(n)%nzong)
                call nc_close(io)

                !          write (fname,'(4a,"regionsg_",a,".hdf")') trim(outdir), pathsep, trim(budget_subdir), pathsep, region_name(n)
                !          call check_dir(trim(fname))
                !          io = SFSTART( trim(fname), DFACC_CREATE )
                !          istat = sfsnatt(io,'im',    DFNT_INT32, 1, im(n))
                !          istat = sfsnatt(io,'jm',    DFNT_INT32, 1, jm(n))
                !          istat = sfsnatt(io,'lm',    DFNT_INT32, 1, lm(n))
                !          call io_write(io,im(n),'im',jm(n),'jm', budg_dat(n)%nzong,'region_global')
                !          print *,'ini_zoneg: closing output file',SFEND(io)
            end if

        end do

        ! ok
        status = 0

    end subroutine ini_zoneg



    subroutine diagbudg(region,k,status)
        !-----------------------------------------------------------------------
        !
        !
        !       interface
        !       ---------
        !       call call_diagbud(region,k)
        !       if k=1   reset budget variables budrr,budrj,buddry,
        !                budwet,budemi,budstrat,budconc,budchng,
        !                initialise budini, allocate memory
        !       if k=2   output of budget variables, free memory
        !       if k=3   accumulate deposition budget
        !
        !       method
        !       ------
        !       output of time accumulated budget files for postprocessing
        !
        !       externals
        !       ---------
        !       none
        !
        !       reference
        !       ---------
        !       none
        !-----------------------------------------------------------------------
        use global_data, only : outdir
        use global_data, only : region_dat, mass_dat
        use io_hdf,      only : io_write
        use io_hdf,      only : DFACC_CREATE, DFNT_INT32, DFNT_INT64, DFNT_CHAR
        use toolbox,     only : escape_tm
        use file_netcdf
        use string_functions
        use netcdf,      only : nf90_put_att, NF90_GLOBAL
        use misctools,   only : check_dir
        use chem_param,  only : ntracet
        use datetime,    only : tau2date
        use dims,        only : itaui, itaue

#ifdef MPI
    use mpi_const
#endif

        implicit none

        ! in/out
        integer,intent(in) :: region,k
        integer, intent(out) :: status

        ! local
        real,dimension(:,:,:,:),pointer    ::rm
        integer,dimension(:,:) ,pointer    ::zoomed,edge
        integer :: nzone,nzone_v
        integer :: i,j,l,n,nt,iy
        character(len=MAX_FILENAME_LEN) :: budget_global_fname

        integer,parameter       :: nproces=15     ! number of budeget processes

        character(len=10), dimension(nproces), parameter :: proces_name=(/&
                'ratnam    ','jnam      ','rwnam     ','marknam   ','names     ', &
                        'names     ','names     ','names     ','names     ','names     ', &
                        'names     ','names     ','names     ','names     ','names     '/)
        integer, dimension(nproces) :: proces_nr
        integer   :: sfstart, sfend, sfsnatt, sfscatt, io, istat, itr
        integer   :: communicator,offsetn,offsetl,lglob,nglob,root_id,lmr
        integer   :: idatei(6), idatee(6)

        character(len=*), parameter :: rname = mname//'/diagbudg'
        !save

        ! start
        if (.not. apply_budget_global) return

        proces_nr = (/ &
            nreac,  njnum,  nreacw,   nmark,  ntrace, &
            ntrace,  ntrace,  ntrace,  ntrace,  ntrace, &
            ntrace, ntracet, ntracet, ntracet, ntracet &
        /)

        !WP diagbud can be called during chemistry, which is parallel
        !   over levels. However, a version parallel over
        !WP tracers is also implemented for future changes.

#ifdef MPI
    which_par=previous_par(region)

    if ( which_par=='tracer'.and.ntracetloc==0) return
    if (which_par=='levels'.and.lmloc==0) return  !WP!
#endif

#ifdef MPI

    if (which_par=='tracer') then

       rm => mass_dat(region)%rm_t
       lmr = lm(region)
       nt=ntracetloc
       communicator=com_trac  !WP! assign com_trac as communicator
       root_id=root_t
       offsetn=sum(ntracet_ar(0:myid-1) )  !offset for global value of n
       offsetl=0  ! no offset for levels
    else if (which_par=='levels') then

       rm => mass_dat(region)%rm_k
       lmr = lmloc
       nt=ntracet
       communicator=com_lev  !WP! assign com_lev as communicator
       root_id=root_k
       offsetl=sum(lmar(0:myid-1) )  ! offset for global value of l
       offsetn=0    ! no offset for tracers

    end if

#else
        rm => mass_dat(region)%rm_t
        lmr = lm(region)
        nt=ntracet
        offsetn=0  !
        offsetl=0  ! no offset for levels
#endif

        select case(k)

        case(1) !               _____________________case 1: initialise

#ifdef MPI
       if ( which_par /= 'tracer' ) &
            call escape_tm( &
            'diagbudg: initializing budgets should be done on tracer domain')
#endif

            if ( region == 1 ) then     !only for region = 1....
                !write(kmain,*) 'diagbudg: reset budget variables'
                budrrg=0.
                budrrg_all=0.
                budrjg=0.
                budrjg_all=0.
                budemig=0.
                budemig_all=0.
                budstratg=0.
                budstratg_all=0.
                budconcg=0.
                budconcg_all=0.
                budinig=0.
                budinig_all=0.
                budmarkg=0.
                budmarkg_all=0.
                budrwg=0.
                budrwg_all=0.
                budconvg=0.0
                budconvg_all=0.0
                budvdifg=0.0
                budvdifg_all=0.0
                budadvxg=0.0
                budadvxg_all=0.0
                budadvyg=0.0
                budadvyg_all=0.0
                budadvzg=0.0
                budadvzg_all=0.0
                !
                do n=1,nregions
                    buddep_dat(n)%lsp = 0.0
                    buddep_dat(n)%cp  = 0.0
                    buddep_dat(n)%cond  = 0.0
                    buddep_dat(n)%dry = 0.0
                    buddep_dat_all(n)%lsp = 0.0
                    buddep_dat_all(n)%cp  = 0.0
                    buddep_dat_all(n)%cond  = 0.0
                    buddep_dat_all(n)%dry = 0.0
                end do
            end if   !only for region 1


            zoomed => region_dat(region)%zoomed
            edge   => region_dat(region)%edge

            ! tracer 1 budget...
#ifdef MPI
       if ( myid == pe_first_tracer ) &
#endif


            do itr = 1, ntracet
                init_mass(region,itr) = sum(rm(1:im(region),1:jm(region),1:lmr,itr))
            end do

            sum_emission(region,:)   = 0.0
            sum_deposition(region,:) = 0.0    !
            sum_wet(region,:)        = 0.0    !
            sum_chemistry(region,:)  = 0.0
            sum_update(region,:)     = 0.0
            sum_advection(region,:)  = 0.0
            !
            ! initialise budinig
            !
            do n=1,nt
                nglob=n+offsetn
                do l=1,lmr
                    nzone_v=nzon_vg(l)
                    do j=1,jm(region)
                        do i=1,im(region)
                            !                   !skip zoomed part of the region
                            !                   if ( zoomed(i,j) /= region .or. edge(i,j) >= 1 ) cycle
                            ! include the edges....parent takes care!
                            nzone=budg_dat(region)%nzong(i,j)
                            budinig(nzone,nzone_v,nglob) = &
                                    budinig(nzone,nzone_v,nglob) + &
                                            !                        rm(i,j,l,n)/ra(nglob)*1000.
                                            rm(i,j,l,n)   !kg tracer mass
                        end do !i
                    end do!j
                end do!l
            end do ! n

            nullify(rm)
            nullify(zoomed)
            nullify(edge)


        case(2) !____________________________write and close, free memory
            !


            zoomed => region_dat(region)%zoomed
            edge   => region_dat(region)%edge

            if ( region == 1 ) budconcg = 0.0
            do n=1,nt
                nglob=n+offsetn !WP! offset is zero on levels domain
                do l=1,lmr
                    lglob=l+offsetl !WP! offset is zero on tracer domain
                    nzone_v=nzon_vg(lglob)
                    do j=1,jm(region)
                        do i=1,im(region)
                            !                   if (zoomed(i,j)/=region .or. edge(i,j) >= 1) cycle
                            nzone=budg_dat(region)%nzong(i,j)
                            budconcg(nzone,nzone_v,nglob) = &
                                    budconcg(nzone,nzone_v,nglob) + &
                                            !                        rm(i,j,l,n)/ra(nglob)*1000.
                                            rm(i,j,l,n) !kg tracer mass
                        end do !i
                    end do!j
                end do!l
            end do !n
            !

            nullify(rm)
            nullify(zoomed)
            nullify(edge)


            ! calculate change and write....
            ! do this in the step for the last region,
            ! since all fields are de-allocated at the end
            if ( region == nregions ) then

                do  n=1,nt
                    nglob=n+offsetn  !WP! offset is zero on levels domain
                    do i=1,nbudg
                        do iy=1,nbud_vg
                            budchngg(i,iy,nglob) = budconcg(i,iy,nglob) - &
                                    budinig(i,iy,nglob)
                        end do ! iy
                    end do !i
                end do !n

#ifdef MPI
          call gather_budget
#endif

#ifdef MPI
          if (myid==root_id) then     !WP! only root writes
#endif

                !
                ! write to hdf file
                !
                if (apply_budget_global) then
                    call tau2date(itaui,idatei)
                    call tau2date(itaue,idatee)
                    write(budget_global_fname, '(a,"/",a,"/budget_global_",i4.4,3i2.2,"_",i4.4,3i2.2,".nc")') trim(outdir), trim(budget_subdir), idatei(1:4), idatee(1:4)
                    !budget_global_fname = trim(outdir)//'/'//trim(budget_subdir)//'/budget_global.nc'
                    call check_dir(budget_global_fname)
                    io = nc_open(trim(budget_global_fname), 'w', status)
                    IF_NOTOK_RETURN(status=1)

                    call nc_create_dim(io, 'n_period', budget_time_profile%n_period)
                    call nc_create_dim(io, 'nreac', nreac)
                    call nc_create_dim(io, 'nbud_v', nbud_vg)
                    call nc_create_dim(io, 'nbud', nbudg)
                    call nc_create_dim(io, 'ntrace', ntrace)
                    call nc_create_dim(io, 'ntracet', ntracet)

                    call nc_dump_var(io, 'budrr', (/'nbud    ','nbud_v  ','ntrace  ','nreac   ','n_period'/), budrrg(:,:,:,:,:))
                    call nc_dump_var(io, 'budemi', (/'nbud    ','nbud_v  ','ntrace  ','n_period'/), budemig(:,:,:,:))
                    call nc_dump_var(io, 'budconc', (/'nbud  ', 'nbud_v', 'ntrace'/), budconcg(:,:,:))
                    call nc_dump_var(io, 'budchng', (/'nbud  ', 'nbud_v', 'ntrace'/), budchngg(:,:,:))
                    call nc_dump_var(io, 'budini', (/'nbud  ', 'nbud_v', 'ntrace'/), budinig(:,:,:))
                    call nc_dump_var(io, 'budadvx', (/'nbud   ', 'nbud_v ', 'ntracet'/), budadvxg(:,:,:))
                    call nc_dump_var(io, 'budadvy', (/'nbud   ', 'nbud_v ', 'ntracet'/), budadvyg(:,:,:))
                    call nc_dump_var(io, 'budadvz', (/'nbud   ', 'nbud_v ', 'ntracet'/), budadvzg(:,:,:))
                    call nc_dump_var(io, 'budconv', (/'nbud   ', 'nbud_v ', 'ntracet'/), budconvg(:,:,:))
                    call nc_dump_var(io, 'budvdif', (/'nbud   ', 'nbud_v ', 'ntracet'/), budvdifg(:,:,:))

                    istat = nf90_put_att(io, NF90_GLOBAL, 'nregions', nregions)
                    istat = nf90_put_att(io, NF90_GLOBAL, 'itau', itau)
                    istat = nf90_put_att(io, NF90_GLOBAL, 'njnum', njnum)
                    istat = nf90_put_att(io, NF90_GLOBAL, 'nmark', nmark)
                    istat = nf90_put_att(io, NF90_GLOBAL, 'nreacw', nreacw)
                    istat = nf90_put_att(io, NF90_GLOBAL, 'idate', idate(1:6))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'idatei', idatei(1:6))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'idatee', idatee(1:6))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'proces_nr', proces_nr(1:nproces))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'proces_name', trim(string_concat(proces_name(1:nproces), ',')))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'zone_name', trim(string_concat(zone_nameg(1:nbudg), ',')))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'zone_namv', trim(string_concat(zone_namvg(1:nbud_vg), ',')))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'iflux1', iflux1(1:nregions))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'iflux2', iflux2(1:nregions))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'jflux1', jflux1(1:nregions))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'jflux2', jflux2(1:nregions))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'lflux1', lflux1(1:nregions))
                    istat = nf90_put_att(io, NF90_GLOBAL, 'lflux2', lflux2(1:nregions))

                    call nc_close(io)

                end if ! if (apply_budget_global) then

                do n=1,nregions   !write and deallocate fields
                    !leave allocated for later use! deallocate( budg_dat(n)%nzong)
                    deallocate( buddep_dat(n)%lsp)
                    deallocate( buddep_dat(n)% cp)
                    deallocate( buddep_dat(n)% cond)
                    deallocate( buddep_dat(n)%dry)
                    deallocate( buddep_dat_all(n)%lsp)
                    deallocate( buddep_dat_all(n)% cp)
                    deallocate( buddep_dat_all(n)% cond)
                    deallocate( buddep_dat_all(n)%dry)

                    deallocate( budget_flux(n)%flux_x1)
                    deallocate( budget_flux(n)%flux_y1)
                    deallocate( budget_flux(n)%flux_z1)
                    deallocate( budget_flux_all(n)%flux_x1)
                    deallocate( budget_flux_all(n)%flux_y1)
                    deallocate( budget_flux_all(n)%flux_z1)
                    deallocate( budget_flux(n)%flux_x2)
                    deallocate( budget_flux(n)%flux_y2)
                    deallocate( budget_flux(n)%flux_z2)
                    deallocate( budget_flux_all(n)%flux_x2)
                    deallocate( budget_flux_all(n)%flux_y2)
                    deallocate( budget_flux_all(n)%flux_z2)

                end do

                deallocate( budemig)
                deallocate( budrrg)

                !             print *, 'diagbudg: Close budget file region',region, sfend(io)

#ifdef MPI
          end if ! root_t
#endif

            end if !region == nregions...


        case(3)
            !
            ! add wet and dry deposition fields to budget
            !
            if (region == 1) then
                !          buddryg=0.
                budwetg_lsp=0.
                budwetg_cp=0.
                !          budwetg_nat=0.
            endif
            do n=1,nt
                nglob=n+offsetn !WP! offset is zero on levels domain
                do nzone_v=1,nbud_vg
                    do j=1,jm(region)
                        do i=1,im(region)
                            nzone=budg_dat(region)%nzong(i,j)
                            budwetg_lsp(nzone,nzone_v,nglob) = &
                                    budwetg_lsp(nzone,nzone_v,nglob) + &
                                            buddep_dat(region)%lsp(i,j,nzone_v,nglob)
                            budwetg_cp (nzone,nzone_v,nglob) = &
                                    budwetg_cp (nzone,nzone_v,nglob) + &
                                            buddep_dat(region)%cp (i,j,nzone_v,nglob)
                            !                    budwetg_nat (nzone,nzone_v,nglob) = &
                            !                        budwetg_nat (nzone,nzone_v,nglob) + &
                            !                        buddep_dat(region)%cond (i,j,nzone_v,nglob)
                        end do
                    end do!j
                end do!nzone_v
            end do !nt
            if (region==nregions) then
                do n=1,ntracet
                    ! compensate double counting wet deposition in convection
                    budconvg(:,:,n) = budconvg(:,:,n) + budwetg_cp(:,:,n)
                end do
            endif
            ! dry deposition: processed parallel over levels --> only
            !                 on PE processed the surface / other have zeros:
            !       do n=1,ntrace   !CMK changed....
            !         do j=1,jm(region)
            !            do i=1,im(region)
            !               nzone = budg_dat(region)%nzong(i,j)
            !               buddryg(nzone,1,n) = buddryg(nzone,1,n) + &
            !                    buddep_dat(region)%dry(i,j,n)
            !            end do
            !         end do!j
            !       end do !

        case default
            call escape_tm('diagbudg: Wrong k entered the budget routines ' )
        end select

#ifdef MPI
  contains


    subroutine gather_budget
      implicit none
      integer :: nsend, region

      nsend=nbudg*nbud_vg*nmark
      call mpi_allreduce(budmarkg,budmarkg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      nsend=nbudg*nbud_vg*nreacw
      call mpi_allreduce(budrwg,budrwg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      nsend=nbudg*nbud_vg*nreac
      call mpi_allreduce(budrrg,budrrg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      nsend=nbudg*nbud_vg*njnum
      call mpi_allreduce(budrjg,budrjg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      nsend=nbudg*nbud_vg*ntrace
      call mpi_allreduce(budinig,budinig_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
!      call mpi_allreduce(buddryg,buddryg_all,nsend,my_real, &
!           mpi_sum,communicator,ierr)
      call mpi_allreduce(budwetg_lsp,budwetg_lsp_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budwetg_cp,budwetg_cp_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
!      call mpi_allreduce(budwetg_nat,budwetg_nat_all,nsend,my_real, &
!          mpi_sum,communicator,ierr)
      call mpi_allreduce(budemig,budemig_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budstratg,budstratg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budconcg,budconcg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budchngg,budchngg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      nsend=nbudg*nbud_vg*ntracet
      call mpi_allreduce(budadvxg,budadvxg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budadvyg,budadvyg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budadvzg,budadvzg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budconvg,budconvg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      call mpi_allreduce(budvdifg,budvdifg_all,nsend,my_real, &
           mpi_sum,communicator,ierr)
      do region = 1, nregions
         nsend = im(region)*jm(region)*nbud_vg*ntrace
         call mpi_allreduce(buddep_dat(region)%lsp, &
              buddep_dat_all(region)%lsp, &
              nsend,my_real,mpi_sum,communicator,ierr)
         call mpi_allreduce(buddep_dat(region)%cp, &
              buddep_dat_all(region)%cp, &
              nsend,my_real,mpi_sum,communicator,ierr)
         nsend = im(region)*jm(region)*ntrace
         call mpi_allreduce(buddep_dat(region)%dry, &
              buddep_dat_all(region)%dry, &
              nsend,my_real,mpi_sum,communicator,ierr)
         nsend = jm(region)*lm(region)*ntracet
         call mpi_allreduce(budget_flux(region)%flux_x1, &
              budget_flux_all(region)%flux_x1, &
              nsend,my_real,mpi_sum,communicator,ierr)
         call mpi_allreduce(budget_flux(region)%flux_x2, &
              budget_flux_all(region)%flux_x2, &
              nsend,my_real,mpi_sum,communicator,ierr)
         nsend = im(region)*lm(region)*ntracet
         call mpi_allreduce(budget_flux(region)%flux_y1, &
              budget_flux_all(region)%flux_y1, &
              nsend,my_real,mpi_sum,communicator,ierr)
         call mpi_allreduce(budget_flux(region)%flux_y2, &
              budget_flux_all(region)%flux_y2, &
              nsend,my_real,mpi_sum,communicator,ierr)
         nsend = im(region)*jm(region)*ntracet
         call mpi_allreduce(budget_flux(region)%flux_z1, &
              budget_flux_all(region)%flux_z1, &
              nsend,my_real,mpi_sum,communicator,ierr)
         call mpi_allreduce(budget_flux(region)%flux_z2, &
              budget_flux_all(region)%flux_z2, &
              nsend,my_real,mpi_sum,communicator,ierr)
      end do

    end subroutine gather_budget
#endif

        status = 0

    end subroutine diagbudg



    subroutine budget_transportg(region,stat,proces,prev_step)
        !
        ! routine to evaluate the budget of advection
        ! routine is called before and after advection
        ! The difference is used to calculate the budget.
        ! Problems occur with zoome regions. Here advection is applied by
        ! writing x-adges and y-edges (in put_xedges, put_yedges).
        ! Knowledge of the previous step is required to
        ! calculate the budget....
        !
        use dims
        use global_data, only : region_dat,mass_dat
        use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_const
#endif

        implicit none

        ! in/out
        integer,intent(in)   :: stat   ! stat = 0 before, stat = 1 after process
        integer,intent(in)   :: region
        character(len=*)     :: proces       ! the process (xyzvcs)
        character(len=*)     :: prev_step    ! the previous process.

        ! const
        character(len=*), parameter  ::  rname = mname//'/budget_transportg'

        ! local
        real,dimension(:,:,:,:),pointer       :: rm
        integer,dimension(:,:) ,pointer       :: zoomed, edge
        integer :: i,j,l,n,nzone,nzone_v,lmr,nt,offsetn,nglob
        logical,dimension(:,:),allocatable    :: skip_flag

        ! start
        if (.not.apply_budget_global) return

#ifdef MPI
    if ( previous_par(region) /= 'tracer' )  &
         call escape_tm('budget_transportg: must be on tracer domain')
    if (ntracetloc==0) return
#endif

        rm => mass_dat(region)%rm_t
        rm => mass_dat(region)%rm_t
        lmr=lm(region)
#ifdef MPI
    nt=ntracetloc
    which_par='tracer'
    offsetn=sum(ntracet_ar(0:myid-1))
#else
        nt=ntracet
        offsetn=0
#endif

        zoomed => region_dat(region)%zoomed
        edge => region_dat(region)%edge

        if ( okdebug ) then
            print *, 'budget_transportg: region, stat, proces, prev_step = ',&
                    region, stat, proces,' ', prev_step
        end if

        allocate(skip_flag(im(region),jm(region)))

        skip_flag = .false.

        where ( zoomed /= region .or. edge >= 1 )
            skip_flag = .true.
        end where

        if ( proces(1:4) == 'advx' .and. prev_step /= 'y' ) then
            ! include y-edges if XYZ sequence
            where(edge >= 2)
                skip_flag = .false.    !only remaining problem: reduced grid!
            end where
        else if ( proces(1:4) == 'advy' .and. prev_step /= 'x' ) then
            where ( edge == 1 .or. edge == 3 )
                skip_flag = .false.
            end where
        else if ( proces(1:4) == 'advz' .and. prev_step /= 'y' ) then
            ! zyx sequence: include interface
            where ( edge >= 1 )
                skip_flag = .false.
            end where
        else if ( proces(1:4) == 'conv' ) then   ! always interface
            where ( edge >= 1 )
                skip_flag = .false.
            end where
        end if

        if ( stat == 0 ) then ! save current tracer concentrations
            rmoldg = 0.0
            do n=1,nt
                nglob=n+offsetn
                do l=1,lmr
                    nzone_v = nzon_vg(l)
                    do j=1,jm(region)
                        do i=1,im(region)
                            if (skip_flag(i,j)) cycle
                            nzone=budg_dat(region)%nzong(i,j)
                            rmoldg(nzone,nzone_v,nglob) = &
                                    rmoldg(nzone,nzone_v,nglob) + rm(i,j,l,n)
                        end do
                    end do
                end do
            end do
        else
            rmnewg = 0.0
            do n=1,nt  !WP! local value of ntracet
                nglob=n+offsetn
                do l=1,lmr
                    nzone_v = nzon_vg(l)
                    do j=1,jm(region)
                        do i=1,im(region)
                            if (skip_flag(i,j)) cycle
                            nzone=budg_dat(region)%nzong(i,j)
                            rmnewg(nzone,nzone_v,nglob) = &
                                    rmnewg(nzone,nzone_v,nglob) + rm(i,j,l,n)
                        end do
                    end do
                end do
            end do

            select CASE(proces(1:4))
            case('advx')
                do n=1,nt  !WP! local value of ntracet
                    nglob=n+offsetn
                    budadvxg(:,:,nglob) = budadvxg(:,:,nglob) + &
                            (rmnewg(:,:,nglob)-rmoldg(:,:,nglob))/ra(nglob)*1e3
                end do
            case('advy')
                do n=1,nt  !WP! local value of ntracet
                    nglob=n+offsetn
                    budadvyg(:,:,nglob) = budadvyg(:,:,nglob) + &
                            (rmnewg(:,:,nglob)-rmoldg(:,:,nglob))/ra(nglob)*1e3
                end do
            case('advz')
                do n=1,nt  !WP! local value of ntracet
                    nglob=n+offsetn
                    budadvzg(:,:,nglob) = budadvzg(:,:,nglob) + &
                            (rmnewg(:,:,nglob)-rmoldg(:,:,nglob))/ra(nglob)*1e3
                end do
            case('conv')
                do n=1,nt  !WP! local value of ntracet
                    nglob=n+offsetn
                    budconvg(:,:,nglob) = budconvg(:,:,nglob) + &
                            (rmnewg(:,:,nglob)-rmoldg(:,:,nglob))/ra(nglob)*1e3
                end do
            case('vdif')
                do n=1,nt  !WP! local value of ntracet
                    nglob=n+offsetn
                    budvdifg(:,:,nglob) = budvdifg(:,:,nglob) + &
                            (rmnewg(:,:,nglob)-rmoldg(:,:,nglob))/ra(nglob)*1e3
                end do
            case default
                print *,'budget_transportg: no valid value for transport budget,', &
                        ' budget not updated!'
                TRACEBACK; stop
            end select
        end if

        deallocate(skip_flag)

        nullify(rm)
        nullify(zoomed)
        nullify(edge)

    end subroutine budget_transportg


end module budget_global
