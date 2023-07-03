!###############################################################################
!
! this module declares and deallocates the actual data
! dimensions are defined in module dims
! MK december 2002
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

module global_data

    use GO, only : gol, goPr, goErr
    use GO, only : TrcFile
    use os_specs, only : MAX_FILENAME_LEN

    use dims,       only : region_data
    use dims,       only : wind_data
    use dims,       only : mass_data
    use global_types, only : conv_data
    use dims,       only : nregions, im, jm, lm, nlon360, nlat180
    use dims,       only : okdebug, lmax_conv
  use chem_param, only : ntracet, emis_data
#ifdef MPI
  use mpi_const,  only : ntracetloc,lmloc,allocate_mass
#endif

    implicit none

    ! --- in/out -----------------------------------

    private

    ! flags
    public :: use_surf

    ! variables

    ! NOTE: ifort compiler requires that both type and variable are public ..

  public  ::  rcfile, rcF
  public  ::  inputdir
  public  ::  outdir

    public  :: region_data, region_dat
    public  :: wind_data, wind_dat
    public  :: mass_data, mass_dat
    public  :: conv_data, conv_dat
    public  :: emis_data
    public  :: albedo

    ! subroutines
    public  :: declare_fields, free_fields
    public  :: free_massdat, assign_massdat


    ! --- const ------------------------------------------

    ! module name
    character(len=*), parameter  ::  mname = 'global_data'


    ! --- var ---------------------------------------------

    ! name of rc file:
    character(len=MAX_FILENAME_LEN)  ::  rcfile
    ! rc file values:
    type(TrcFile)       :: rcF

    ! path to input files:
    character(len=MAX_FILENAME_LEN)  ::  inputdir = './'

    ! path to scripts:
    character(len=MAX_FILENAME_LEN)  ::  bindir   = './'

    ! path to output files:
    character(len=MAX_FILENAME_LEN)  ::  outdir        = './'

    logical,parameter :: use_surf   = .true. !flag to declare/read/use oro/albedo

    type(region_data),dimension(nregions),target   :: region_dat
    type(wind_data)  ,dimension(nregions),target   :: wind_dat
    type(mass_data)  ,dimension(nregions),target   :: mass_dat
    type(conv_data)  ,dimension(nregions),target   :: conv_dat

    type(emis_data)  ,dimension(nregions),target   :: albedo


contains


    subroutine declare_fields( status )
        !
        !    subroutine to allocate the memory for all the regions
        !    should be called at start of program
        !

        use GO  , only : NewDate

        implicit none

        ! in/out
        integer, intent(out)   ::  status

        ! const
        character(len=*), parameter  ::  rname = mname//'/declare_fields'

        ! local
        integer   :: region
        integer   :: imr,jmr,lmr

        ! start
        do region=1,nregions
            imr = im(region)
            jmr = jm(region)
            lmr = lm(region)
            allocate ( region_dat(region)%zoomed(imr,jmr))
            allocate ( region_dat(region)%edge(imr,jmr))
            allocate ( region_dat(region)%dxyp(jmr))

            allocate ( wind_dat(region)%am_t(0:imr+1,0:jmr+1,0:lmr+1) )
            allocate ( wind_dat(region)%bm_t(0:imr+1,0:jmr+1,0:lmr+1) )
            allocate ( wind_dat(region)%cm_t(0:imr+1,0:jmr+1,0:lmr+1) )

            ! diffusion coeff:
            allocate(conv_dat(region)%dkg(imr,jmr,lmax_conv))
            conv_dat(region)%dkg = 0.0
            allocate(conv_dat(region)%kvh(imr,jmr,lmr))
            conv_dat(region)%kvh = 0.0
            allocate(conv_dat(region)%blh(imr,jmr))

            allocate ( conv_dat(region)%cloud_base(imr,jmr))
            allocate ( conv_dat(region)%cloud_top (imr,jmr))
            allocate ( conv_dat(region)%cloud_lfs (imr,jmr))

            if ( use_surf ) then
                allocate ( albedo(region)%surf(imr,jmr))
                albedo(region)%surf = 0.05    ! default value
            end if

            call assign_massdat(region,1) ! WP! initialize mass_dat on tracer domain
#ifdef MPI
       ! WP! initialize mass_dat on level domain
       if ( .not. allocate_mass ) call assign_massdat(region,0)
#endif

            if ( okdebug ) print *, 'declare_fields: Fields for grid, mass, ', &
                    'and winds allocated, region:', region
        end do

        ! ok
        status = 0

    end subroutine declare_fields



    subroutine free_fields( status )
        !
        !    subroutine to free the memory for all the regions
        !    should be called at end of program
        !
#ifdef MPI
    use mpi_const,only:myid,allocate_mass
#endif

        implicit none

        ! in/out
        integer, intent(out)    ::  status

        ! const
        character(len=*), parameter  ::  rname = mname//'/free_fields'

        ! local
        integer   :: region

        ! start
        do region=1,nregions
#ifdef MPI
       if ( okdebug ) print *, 'free_fields: deallocate... ', &
            myid, region,im(region),jm(region)
#else
            if ( okdebug ) print *, 'free_fields: deallocate...', &
                    region,im(region),jm(region)
#endif

            deallocate ( region_dat(region)%zoomed )
            deallocate ( region_dat(region)%edge )
            deallocate ( region_dat(region)%dxyp )

            deallocate ( wind_dat(region)%am_t )
            deallocate ( wind_dat(region)%bm_t )
            deallocate ( wind_dat(region)%cm_t )

            deallocate ( conv_dat(region)%dkg )
            deallocate ( conv_dat(region)%kvh )
            deallocate ( conv_dat(region)%blh )

            deallocate ( conv_dat(region)%cloud_base )
            deallocate ( conv_dat(region)%cloud_top  )
            deallocate ( conv_dat(region)%cloud_lfs  )

            if (use_surf) then
                deallocate ( albedo(region)%surf)
            endif
            call free_massdat(region,0) ! WP! free mass_dat on tracer domain
#ifdef MPI
       ! WP! free mass_dat on level domain
       if ( .not. allocate_mass ) call free_massdat(region,1)
       if ( okdebug ) print *, 'free_fields: Fields for grid, mass, ', &
            'and winds deallocated, region:', myid, region
#else
            if ( okdebug ) print *, 'free_fields: Fields for grid, mass, ', &
                    'and winds deallocated, region:', region
#endif
        end do

        ! ok
        status = 0

    end subroutine free_fields



    subroutine free_massdat(reg,iaction)

        implicit none

        ! input/output
        integer,intent(in) :: reg
        integer,intent(in) :: iaction

        ! local
        integer            :: imr,jmr,lmr,nt

        ! start
        if ( iaction == 0 ) then ! swapping from tracer to k
            deallocate ( mass_dat(reg)%rm_t)
            deallocate ( mass_dat(reg)%rxm_t)
            deallocate ( mass_dat(reg)%rym_t)
            deallocate ( mass_dat(reg)%rzm_t)
#ifdef secmom
       deallocate ( mass_dat(reg)%rxxm_t)
       deallocate ( mass_dat(reg)%rxym_t)
       deallocate ( mass_dat(reg)%rxzm_t)
       deallocate ( mass_dat(reg)%ryym_t)
       deallocate ( mass_dat(reg)%ryzm_t)
       deallocate ( mass_dat(reg)%rzzm_t)
#endif
        else if( iaction == 1 ) then !WP! swapping from k to tracer
#ifdef MPI
       deallocate ( mass_dat(reg)%rm_k)
       deallocate ( mass_dat(reg)%rxm_k)
       deallocate ( mass_dat(reg)%rym_k)
       deallocate ( mass_dat(reg)%rzm_k)
#ifdef secmom
       deallocate ( mass_dat(reg)%rxxm_k)
       deallocate ( mass_dat(reg)%rxym_k)
       deallocate ( mass_dat(reg)%rxzm_k)
       deallocate ( mass_dat(reg)%ryym_k)
       deallocate ( mass_dat(reg)%ryzm_k)
       deallocate ( mass_dat(reg)%rzzm_k)
#endif
#endif
        end if

    end subroutine free_massdat



    subroutine assign_massdat(reg,iaction)

        implicit none

        ! input/output
        integer,intent(in) :: reg
        integer,intent(in) :: iaction

        ! local
        integer            :: imr,jmr,lmr,nt

        ! start
        if ( iaction == 0 ) then !WP! swapping from tracer to k
#ifdef MPI
       imr=im(reg) ; jmr=jm(reg) ; lmr= lmloc ; nt = ntracet
       allocate ( mass_dat(reg)%rm_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rxm_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rym_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rzm_k(-1:imr+2,-1:jmr+2,lmr,nt) )
#ifdef secmom
       allocate ( mass_dat(reg)%rxxm_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rxym_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rxzm_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%ryym_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%ryzm_k(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rzzm_k(-1:imr+2,-1:jmr+2,lmr,nt) )
#endif
       mass_dat(reg)%rm_k=0.0
       mass_dat(reg)%rxm_k=0.0
       mass_dat(reg)%rym_k=0.0
       mass_dat(reg)%rzm_k=0.0
#ifdef secmom
       mass_dat(reg)%rxxm_k=0.0
       mass_dat(reg)%rxym_k=0.0
       mass_dat(reg)%rxzm_k=0.0
       mass_dat(reg)%ryym_k=0.0
       mass_dat(reg)%ryzm_k=0.0
       mass_dat(reg)%rzzm_k=0.0
#endif
#endif
        else if( iaction == 1 ) then ! swapping from k to tracer
            imr=im(reg) ; jmr=jm(reg) ; lmr= lm(reg)
#ifdef MPI
       nt = ntracetloc
#else
            nt = ntracet
#endif
            allocate ( mass_dat(reg)%rm_t(-1:imr+2,-1:jmr+2,lmr,nt) )
            allocate ( mass_dat(reg)%rxm_t(-1:imr+2,-1:jmr+2,lmr,nt) )
            allocate ( mass_dat(reg)%rym_t(-1:imr+2,-1:jmr+2,lmr,nt) )
            allocate ( mass_dat(reg)%rzm_t(-1:imr+2,-1:jmr+2,lmr,nt) )
#ifdef secmom
       allocate ( mass_dat(reg)%rxxm_t(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rxym_t(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rxzm_t(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%ryym_t(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%ryzm_t(-1:imr+2,-1:jmr+2,lmr,nt) )
       allocate ( mass_dat(reg)%rzzm_t(-1:imr+2,-1:jmr+2,lmr,nt) )
#endif
            mass_dat(reg)%rm_t=0.0
            mass_dat(reg)%rxm_t=0.0
            mass_dat(reg)%rym_t=0.0
            mass_dat(reg)%rzm_t=0.0
#ifdef secmom
       mass_dat(reg)%rxxm_t=0.0
       mass_dat(reg)%rxym_t=0.0
       mass_dat(reg)%rxzm_t=0.0
       mass_dat(reg)%ryym_t=0.0
       mass_dat(reg)%ryzm_t=0.0
       mass_dat(reg)%rzzm_t=0.0
#endif
        end if

    end subroutine assign_massdat


end module global_data
