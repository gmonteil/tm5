#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if

module chem_param

    use binas,          only : xmair, Avog, grav, amu
    use global_types,   only : emis_data
    use global_data,    only : ntracet
    use os_specs,       only : MAX_FILENAME_LEN
    use go,             only : gol, goerr, T_Time_Window, NewDate

    implicit none

    type d3_data
        ! d3_data%d3   : 3D data, e.g. nox emissions
        real, dimension(:,:,:), allocatable :: d3
    end type d3_data

    type d4_data
        type(d3_data), allocatable :: field_t(:) ! time dimension added to d3_data
    end type d4_data

    type d23_data
        ! d23_data%d23  : lat/pres fields, e.g. o3 climatology in ppmv
        real, dimension(:,:), pointer     :: d23
    end type d23_data

    type d2_data
        ! d2_data%d2  : lat fields, e.g. hno3/o3 ratios at 10 hPa
        real, dimension(:), pointer       :: d2
    end type d2_data

    type isop_data
        ! isop_data%scalef_isop  : (jm,ntim) scalefactor isoprene emissions
        real, dimension(:,:), pointer   :: scalef_isop
    end type isop_data

    type chem_data
        ! chem_data%rm_k  : 'chemistry' tracers are parallel over layers...
        real, dimension(:,:,:,:), pointer :: rm_k
    end type chem_data

    ! definition of the chemistry: #reactions, order of species, etc.
    ! parameters needed for chemistry and rate constants

    ! nmark: number of 'marked' tracers
    integer, parameter :: nmark = 0
    character(len=8), dimension(nmark) :: marknam != [] !(/'CO2     '/)

    integer :: ntrace       ! number of tracers for chemistry
    integer :: ntrace_chem  ! number of non-transported tracers
    integer :: maxtrace     ! total number of tracers
    integer, parameter :: tracer_name_len = 20

    ! components numbers
!    integer, parameter :: ico2 = 1

    ! additional fields used in chemistry routine alone
    ! (more meteo-like files in units different from #/cm3)
    integer, parameter                  :: n_extra = 0
    integer, parameter                  :: nstd    = 1
!    integer, dimension(nstd), parameter :: istd    = (/ico2/)
    !
    ! species name
    character(len=tracer_name_len), dimension(:), allocatable :: names
    !
    ! some stuff for the budgets....
    integer, parameter           :: njnum = 1
    integer, parameter           :: nrdep = 0
    integer, parameter           :: nreac = 1
    integer, parameter           :: nreacw = 1
    character(len=tracer_name_len), dimension(njnum)  :: jnam   = (/'dummy   '/)
    character(len=tracer_name_len), dimension(nreacw) :: rwnam  = (/'dummy   '/)
    character(len=tracer_name_len), dimension(nreac)  :: ratnam = (/'dummy   '/)

    character(len=10), dimension(:), allocatable :: mixrat_unit_name
    character(len=10), dimension(:), allocatable :: emis_unit_name
    real, dimension(:), allocatable :: ra, mixrat_unit
    real, dimension(:), allocatable :: emis_unit, fscale
    real, parameter         :: factor_airdens = 1.0D-6/(xmair * amu * grav) ! without the 'D-6', we would get molec/m^3

    ! tracer indices on which dry depositions should be applied:
    ! deposition parameters:
    integer, parameter   :: ndep = 0

    integer, parameter  :: ntlow = 160
    integer, parameter  :: nthigh = 330
    integer, parameter  :: ntemp = nthigh - ntlow

    type react_t
        character(len=MAX_FILENAME_LEN) :: file         ! path to the file storing the reactive species concentration
        character(len=20)               :: domain       ! domain of application of the reaction: total, tropo or strato
        real, dimension(2)              :: parameters   ! reaction parameters
        character(len=50)               :: name         ! reaction name (as in the rc-file)
        real, dimension(ntlow:nthigh)   :: rate = 0     ! reaction rate as function of temperature
        logical                         :: climatology = .false.
        character(len=1)                        :: data_timestep
        type(T_Time_Window)                     :: data_period
        real, dimension(:, :, :), allocatable   :: data ! in-memory field

        !contains
        !    procedure init => init_reaction
    end type react_t

    type tracer_t
        character(len=tracer_name_len)  :: species      ! Name of the chemical species
        logical                         :: has_chem     ! Turn on or off chemistry for that tracer
        integer                         :: nreact
        type(react_t), dimension(:), allocatable    :: reactions
    end type tracer_t
    character(len=20), dimension(:), allocatable    :: reacnames_tmp    ! Temporary array for storing reaction names

    type(tracer_t), dimension(:), allocatable, target   :: tracers

    contains

    subroutine init_chem

        use global_data, only : rcf
        use go, only : readrc

        integer :: status
        integer :: itrac, ireac, nreac
        character(len=500) :: rcval
        character(len=*), parameter                 :: rname = 'chem_param/init_chem'

        call readrc(rcf, 'tracers.number', ntrace, status, default=1)
        if (ntrace == 0) return
        
        call readrc(rcf, 'tracers.names', rcval, status)

        allocate(names(ntrace))
        allocate(ra(ntrace))
        allocate(mixrat_unit(ntrace), mixrat_unit_name(ntrace))
        allocate(emis_unit(ntrace), emis_unit_name(ntrace))
        allocate(fscale(ntrace))

        allocate(tracers(ntrace))

        ntracet = ntrace
        ntrace_chem = ntrace - ntracet
        maxtrace = ntrace

        read(rcval, *), names

        do itrac = 1, ntrace
            call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.molar_mass', ra(itrac), status)
            IF_NOTOK_RETURN(status=1)

            ! mixrat_unit_value should be 1.e6 for ppm, 1.e9 for ppb, etc.
            call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.mixrat_unit_value', mixrat_unit(itrac), status)
            IF_NOTOK_RETURN(status=1)
            call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.mixrat_unit_name', mixrat_unit_name(itrac), status)
            IF_NOTOK_RETURN(status=1)
            ! emis_unit_value should enable conversion from kg(tracer) to unit(tracer), e.g. if we want PgC for a CO2
            ! tracer, then it should be ~(12/44)*1.e-12.
            ! The unit names can be arbitrary (no check on that ...)
            call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.emis_unit_value', emis_unit(itrac), status)
            IF_NOTOK_RETURN(status=1)
            call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.emis_unit_name', emis_unit_name(itrac), status)
            IF_NOTOK_RETURN(status=1)

            ! Initialize the "tracers" object
            call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.chemistry', tracers(itrac)%has_chem, status)
            IF_NOTOK_RETURN(status=1)
            call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.species', tracers(itrac)%species, status)
            IF_NOTOK_RETURN(status=1)

            if (tracers(itrac)%has_chem) then

                ! Number of reactions
                call readrc(rcf, 'tracers.' //trim(names(itrac)) // '.nreac', nreac, status, 0)
                IF_NOTOK_RETURN(status=1)
                tracers(itrac)%nreact = nreac
                allocate(reacnames_tmp(nreac))

                ! Reaction names
                call readrc(rcf, 'tracers.' //trim(names(itrac)) //'.reaction_names', rcval, status)
                IF_NOTOK_RETURN(status=1)
                read(rcval, *) reacnames_tmp

                ! Read reaction parameters
                allocate(tracers(itrac)%reactions(nreac))
                do ireac = 1, nreac
                    tracers(itrac)%reactions(ireac)%name = reacnames_tmp(ireac)
                    call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.' // trim(reacnames_tmp(ireac)) // '.file', tracers(itrac)%reactions(ireac)%file, status)
                    IF_NOTOK_RETURN(status=1)
                    call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.' // trim(reacnames_tmp(ireac)) // '.domain', tracers(itrac)%reactions(ireac)%domain, status)
                    IF_NOTOK_RETURN(status=1)
                    call readrc(rcf, 'tracers.' // trim(names(itrac)) // '.' // trim(reacnames_tmp(ireac)) // '.rate', tracers(itrac)%reactions(ireac)%parameters, status)
                    IF_NOTOK_RETURN(status=1)
                    call init_reaction(tracers(itrac)%reactions(ireac))
                end do

                deallocate(reacnames_tmp)
            end if
        end do

        !   fscale(ntrace): scaling factor for conversion of mixing ratios
        !   in kg tracer per kg air to practical mixing ratio units (e.g. ppm)
        !   In this version: ratio of molecular weight  of tracer to that of air
        fscale = xmair / ra(:) * mixrat_unit(:)

    end subroutine init_chem

    subroutine init_reaction(self)
        !class(react_t), intent(inout)   :: self
        type(react_t), intent(inout)   :: self
        integer                         :: itemp
        real                            :: ztrec
        do itemp = ntlow, nthigh
            ztrec = 1. / real(itemp)
            self%rate(itemp) = self%parameters(1) * exp(self%parameters(2) * ztrec)
        end do
        self%data_period%t1 = NewDate(year=1900, month=1, day=1)
        self%data_period%t2 = NewDate(year=1900, month=1, day=1)
    end subroutine init_reaction

end module chem_param
