!#################################################################
!
! routine to perform main processing and
! handles time ordering of model processes
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

module modelIntegration

  use GO, only : gol, goPr, goErr
  use Grid, only : TllRegion
  use emission_data, only : source_apply

  implicit none


  ! --- in/out -----------------------------

  private

  public :: Proces_Init
  public :: Proces_Done
  public :: Proces_Setup

  public :: proces_region



  ! --- const ----------------------------------

  ! module name
  character(len=*), parameter  ::  mname = 'ModelIntegration'

  ! --- var ---------------------------------

  character(len=1)             ::  output_after_step
  logical                      ::  print_masses
  character(len=1)             ::  print_mass_step

  ! timer handles:
  integer ::  itim_advectx, itim_advecty, itim_advectz
  integer ::  itim_source_sink, itim_chemie
  integer ::  itim_put_edges
  integer ::  itim_output
  integer ::  itim_update_parent

  logical    ::  chemistry_apply
!  logical    ::  source_apply
!  logical    ::  convec_apply
!  logical    ::  vdiff_apply

  ! mask area:
  logical           ::  mask_apply
  real              ::  mask_region(4)  ! west, east, south, north (degees)
  type(TllRegion)   ::  mask_llr
  real              ::  mask_factor
  logical           ::  mask_complement


contains


  ! =====================================================


  subroutine Proces_Init( status )

    use GO                 , only : GO_Timer_Def
    use GO                 , only : ReadRc
    use Grid               , only : Init
    use global_data        , only : rcF

    use dims               , only : nregions, nregions_all, iglbsfc
    use MeteoData          , only : sp1_dat, sp2_dat, sp_dat, spm_dat
    use MeteoData          , only : tsp_dat
    use MeteoData          , only : phlb_dat, m_dat
    use MeteoData          , only : mfu_dat, mfv_dat, mfw_dat
    use MeteoData          , only : pu_dat, pv_dat, pw_dat
    use MeteoData          , only : omega_dat
    use MeteoData          , only : lwc_dat, iwc_dat, cc_dat, ccu_dat, cco_dat
    use MeteoData          , only : entu_dat, entd_dat, detu_dat, detd_dat
    use MeteoData          , only : temper_dat, humid_dat
    use MeteoData          , only : gph_dat
    use MeteoData          , only : oro_dat
    use MeteoData          , only : pclim_dat
    use MeteoData          , only : blh_dat
    use Meteo              , only : mdat_set
    use Meteo              , only : Meteo_Alloc
    use Advect             , only : Advect_Init
    use Sources_Sinks      , only : Sources_Sinks_Init
    use Emission           , only : Emission_Init
#ifndef without_dry_deposition
    use Dry_Deposition     , only : Dry_Deposition_Init
#endif
#ifndef without_chemistry
    use Chemistry          , only : Chemistry_Init
#endif
#ifndef without_convection
    use Convection         , only : Convection_Init
#ifndef without_diffusion
    use Diffusion          , only : Diffusion_Init
#endif
#endif

    implicit none

    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Proces_Init'

    ! --- local ---------------------------------

    integer           ::  n

    ! --- begin ---------------------------------
    write (gol,'(a," : entering")') trim(rname) ; call goPr

    !
    ! which processes ?
    !

!    call ReadRc( rcF, 'proces.convec'   , convec_apply , status   )
!    IF_NOTOK_RETURN(status=1)
!    call ReadRc( rcF, 'proces.vdiff'    , vdiff_apply  , status, default=convec_apply   )
!    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'proces.chemistry', chemistry_apply, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'proces.source'   , source_apply, status    )
    IF_NOTOK_RETURN(status=1)

    call ReadRc(rcF, 'do_steps.print.mass', print_masses, status, default=.false.)
    IF_ERROR_RETURN(status=1)

    if (print_masses) then
        call ReadRc(rcF, 'do_steps.print.mass.step', print_mass_step, status, default='a')
        IF_ERROR_RETURN(status=1)
    end if

    !
    ! area mask ?
    !

    ! apply mask ?
    call ReadRc( rcF, 'mask.apply', mask_apply, status )
    IF_NOTOK_RETURN(status=1)
    ! setup ?
    if ( mask_apply ) then
      ! read domain:
      call ReadRc( rcF, 'mask.region', mask_region, status )
      IF_NOTOK_RETURN(status=1)
      ! setup region structure:
      call Init( mask_llr, mask_region(1), mask_region(2), mask_region(3), mask_region(4), status )
      IF_NOTOK_RETURN(status=1)
      ! factor applied:
      call ReadRc( rcF, 'mask.factor', mask_factor, status )
      IF_NOTOK_RETURN(status=1)
      ! factor applied:
      call ReadRc( rcF, 'mask.complement', mask_complement, status )
      IF_NOTOK_RETURN(status=1)
    end if


    !
    ! init processes
    !

    ! NOTE: these inits should select the appropriate meteo

    call Emission_Init( status )
    IF_NOTOK_RETURN(status=1)

    call Sources_Sinks_Init( status )
    IF_NOTOK_RETURN(status=1)

    call Advect_Init( status )
    IF_NOTOK_RETURN(status=1)

#ifndef without_convection
    ! init convec:
    call Convection_Init( status )
    IF_NOTOK_RETURN(status=1)
#ifndef without_diffusion
    ! init diffusion:
    call Diffusion_Init( status )
    IF_NOTOK_RETURN(status=1)
#endif
#endif

#ifndef without_chemistry
    ! initialize chemistry:
    call Chemistry_Init( status )
    IF_NOTOK_RETURN(status=1)
#endif

#ifndef without_dry_deposition
    call Dry_Deposition_Init( status )
    IF_NOTOK_RETURN(status=1)
#endif


    !
    ! allocate meteo
    !

    !>>> should be done by processes in future!
    ! loop over zoom regions:
    do n = 1, nregions

      ! enable meteo:
      call mdat_set( sp1_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      call mdat_set( sp2_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      call mdat_set( sp_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call mdat_set( tsp_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call mdat_set( spm_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call mdat_set( phlb_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( m_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call mdat_set( mfu_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( mfv_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( mfw_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call mdat_set( pu_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( pv_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( pw_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call mdat_set( temper_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( humid_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call mdat_set( gph_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      !! enable meteo:
      !call mdat_set( lwc_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call mdat_set( iwc_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call mdat_set( cc_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call mdat_set( ccu_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call mdat_set( cco_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)

      ! convec fields enabled in Convection_Init ...

      ! enable meteo:
      call mdat_set( blh_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

    end do

    ! enable meteo:
    call mdat_set( oro_dat(iglbsfc), status, used=.true. )
    IF_NOTOK_RETURN(status=1)

    !<<<

    ! some meteo implies other ...
    do n = 1, nregions_all
      ! I want to use pclim for CO, and therefore I set pclim_dat(:)%used = .true. from within emission_fwd_CO/emission_fwd_init.
      ! Trouble is, that's called from tracer_init, after which tracer_model is called. Tracer_model calls meteo_init at each
      ! time step, which resets pclim_dat(:)%used to .false. Therefore, I need to set it back to .true. here explicitly. The same
      ! problem might occur for met fields needed for dry deposition, since they're set to .true. in dry_deposition_init. We'll
      ! cross that bridge when we get to it.
      call mdat_set( pclim_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! pclim computed from oro:
      if ( pclim_dat(n)%used ) then
        call mdat_set(    oro_dat(n), status, used=.true. )
      end if

      ! convec requires gph and omega:
      if ( entu_dat(n)%used ) then
        call mdat_set(    gph_dat(n), status, used=.true. )
        call mdat_set(  omega_dat(n), status, used=.true. )
      end if

      ! omega (Pa/s) is computed form mfw (kg/s):
      if ( omega_dat(n)%used ) then
        call mdat_set(    mfw_dat(n), status, used=.true. )
      end if
      ! gph is computed from oro/sp/temper/humid
      if ( gph_dat(n)%used ) then
        call mdat_set(    oro_dat(n), status, used=.true. )
        call mdat_set(     sp_dat(n), status, used=.true. )
        call mdat_set( temper_dat(n), status, used=.true. )
        call mdat_set(  humid_dat(n), status, used=.true. )
      end if
      ! sp is interpolated in time between sp1 and sp2:
      if ( sp_dat(n)%used ) then
        call mdat_set(    sp1_dat(n), status, used=.true. )
        call mdat_set(    sp2_dat(n), status, used=.true. )
      end if
      ! cloud covers and overhead/underfeet cloud covers
      if ( cc_dat(n)%used ) then
        call mdat_set( cco_dat(n), status, used=.true. )
        call mdat_set( ccu_dat(n), status, used=.true. )
      end if
    end do

    ! allocate used meteo
    call Meteo_Alloc( status )
    IF_NOTOK_RETURN(status=1)



    !
    ! timing
    !

    ! define ...
    call GO_Timer_Def( itim_advectx    , 'advectx'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_advecty    , 'advecty'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_advectz    , 'advectz'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_put_edges , 'put edges'     , status )
    IF_NOTOK_RETURN(status=1)
    if ( source_apply ) then
      call GO_Timer_Def( itim_source_sink, 'source_sink' , status )
      IF_NOTOK_RETURN(status=1)
    end if
    call GO_Timer_Def( itim_output , 'user_output'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_chemie , 'chemie'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_update_parent , 'update_parent'     , status )
    IF_NOTOK_RETURN(status=1)

    ! read how to sample the output:
    ! x,y,z,c,v,e   : always after this process
    ! a             : after all steps (testing, not recommended)
    ! o             : the 'old' way (default)

    call ReadRc( rcF, 'output.after.step', output_after_step, status)
    IF_ERROR_RETURN(status=1)
    write(gol, *) '****************************************' ; call goPr
    write(gol, *) 'Output will be collected after step :  ',output_after_step ; call goPr
    write(gol, *) '****************************************' ; call goPr

    write (gol,'(a," : done")') trim(rname) ; call goPr
    !
    ! done
    !

    ! ok
    status = 0

  end subroutine Proces_Init


  ! ***


  subroutine Proces_Done( status )

    use Grid  , only : Done
    use Advect, only : Advect_Done
#ifndef without_convection
    use Convection    , only : Convection_Done
#ifndef without_diffusion
    use Diffusion     , only : Diffusion_Done
#endif
#endif
#ifndef without_chemistry
    use Chemistry, only : Chemistry_Done
#endif
#ifndef without_dry_deposition
    use dry_deposition, only : Dry_Deposition_Done
#endif

    implicit none

    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Proces_Done'

    ! --- begin ---------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! done with advection:
    call Advect_Done( status )
    IF_NOTOK_RETURN(status=1)

#ifndef without_convection
    ! done with convec:
     call Convection_Done( status )
    IF_NOTOK_RETURN(status=1)
#ifndef without_diffusion
    ! done with diffusion:
    call Diffusion_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif
#endif

#ifndef without_chemistry
    ! done with chemistry:
    call Chemistry_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif

#ifndef without_dry_deposition
    call Dry_Deposition_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! mask defined ?
    if ( mask_apply ) then
      ! done with region structure:
      call Done( mask_llr, status )
      IF_NOTOK_RETURN(status=1)
    end if

    write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Proces_Done


  ! ***


  ! called every timestep, setup OH fields etc

  subroutine Proces_Setup( t1, t2, status )

    use GO            , only : TDate
#ifndef without_diffusion
    use Diffusion     , only : Diffusion_Setup
#endif
#ifndef without_dry_deposition
    use dry_deposition, only : Dry_Deposition_Calc
#endif
    use Emission_fwd  , only : Emission_fwd_After_Read
!#ifndef without_chemistry
    !use Chemistry     , only : Chemistry_Setup
!#endif
!#ifndef without_dry_deposition
    !use Dry_Deposition, only : Dry_Deposition_Setup
!#endif

    implicit none

    ! --- in/out --------------------------------

    type(TDate), intent(in)   ::  t1, t2
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Proces_Setup'

    ! --- begin ---------------------------------
    !write (gol,'(a," : entering")') trim(rname) ; call goPr
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! from 'sources_sinks/trace_after_read

#ifndef without_diffusion
    ! setup diffusion fields:
    call Diffusion_Setup( t1, t2, status )
    IF_NOTOK_RETURN(status=1)
#endif

#ifndef without_dry_deposition
    call Dry_Deposition_Calc( t1, t2, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! update emission fields (empty routine?)
    call Emission_Fwd_After_Read( status )
    IF_NOTOK_RETURN(status=1)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!#ifndef without_chemistry
    !! setup chemistry:
    !call Chemistry_Setup( t1, t2, status )
    !IF_NOTOK_RETURN(status=1)
!#endif

!#ifndef without_dry_deposition
    !call Dry_Deposition_Setup( status )
    !IF_NOTOK_RETURN(status=1)
!#endif
    !write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Proces_Setup


  ! ***


  !
  ! for a given region, performs tref(region)/tref(parent(region))
  ! advection steps for each direction x,y,z
  ! including all region's children, granchildren, etc.
  ! written by patrick berkvens and mike botchev, march-june 1999
  !WP! Also chemistry and sources_sinks are called during advection,
  !WP! wouter peters, july 2000
  !
  !example (cmk)
  !   XYZ                              ECCEZYX
  !      xyz          eccezyx                  cezyx          xyz          ec
  !         xyzeccezyx       ecxyzzyxec             cezyxxyzec   zyxeccexyz
  ! in this case the operations e and c should be executed in the
  ! interface region
  ! but not in the core of the zoom region (otherwise double counting)
  ! This results in the most simple algorithm, because you should leave
  ! the DO_STEPS routine AFTER
  ! (1) finishing all the steps...
  ! (2) XYZ
  ! MK 2000...
  !

  recursive subroutine proces_region( region, tr, status )

    use GO          , only : TDate, IncrDate, rTotal, operator(+), operator(-), wrtgol
    use GO          , only : GO_Timer_Start, GO_Timer_End
    use dims        , only : tref, revert, ndyn, itaur
    use dims        , only : parent, children
    use dims        , only : statusr => status
    use dims        , only : n_operators
    use user_output,  only : user_output_step
    use zoom_tools,   only : update_parent
    !use advect_tools, only : m2phlb, m2phlb1
    use advect_tools, only : mass_to_pressure

    implicit none

    ! --- in/out -------------------------------

    integer, intent(in)       ::  region
    type(TDate), intent(in)   ::  tr(2)
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Proces_Region'

    ! --- local ---------------------------------

    integer             ::  child, i, ichild, tref_, n_children, dtime
    integer             ::  nsec_
    type(TDate)         ::  tr_(2)
    character(len=64)   ::  line

    ! --- begin ----------------------------------
    !write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! determine refinement factor with respect to the parent
    tref_ = tref(region)/tref(parent(region))
    !main time step for all processes. note that all timesteps are default/2
    dtime = revert*ndyn/(2*tref(region))
    !e.g. xyzvec  (firts call) and cevzyx (second call)

    ! time step in region:
    nsec_ = nint( rTotal( tr(2) - tr(1), 'sec' ) / tref_ )
    ! check ...
    if ( nsec_ /= dtime ) then
      write (gol,'("unexpected time steps : nsec_ = ",i0,", while dtime = ",i0)') nsec_, dtime; call goErr
      !TRACEBACK; status=1; return
    end if
    ! loop over time steps in region:
    do i = 1, tref_

       ! timerange in region:
       tr_(1) = tr (1) + IncrDate( sec=nsec_*(i-1) )
       tr_(2) = tr_(1) + IncrDate( sec=nsec_ )

       call do_steps( region, tr_, status )
       IF_NOTOK_RETURN(status=1)

       ! CALL advect_region for all the children (IF there are any)

       ichild = 0
       do while ( ichild < children(region,0) )
          ichild = ichild + 1
          child = children(region,ichild)
          call proces_region( child, tr_, status )
          IF_NOTOK_RETURN(status=1)
       end do
       !do the remaining steps if necessary...!
       if ( mod(statusr(region),n_operators) /= 0 ) then
          call do_steps( region, tr_, status )
          IF_NOTOK_RETURN(status=1)
       end if
       itaur(region) = itaur(region) + dtime     !update region time.....

    end do

    ! update parent using child data ?
    ! compute pressure grid from (changed) air mass:

    call GO_Timer_Start( itim_update_parent, status )
    IF_NOTOK_RETURN(status=1)

    !if ( region /= 1 ) then
       !call update_parent( region )
       !call m2phlb( region )
    !else
       !call m2phlb1( 1 )
    !end if
    if (region > 1) call update_parent(region)
    call mass_to_pressure(region)

    call GO_Timer_End( itim_update_parent, status )
    IF_NOTOK_RETURN(status=1)

    !if(okdebug) call Par_Check_mass(region,'after proces_region')

    ! Accumulate output, write info at stations
    ! Note: the edges in the interface may not have been updated:
    !   this is the task of the parent.
    !write (gol,'(a," : done")') trim(rname) ; call goPr

  end subroutine proces_region


  ! ***


  subroutine do_steps( region, tr, status )

    !-------------------------------------------------------------
    ! performs next three (x-, y- and z-) advection steps
    ! written by patrick berkvens and mike botchev, march-june 1999
    !wp! added the call to chemistry. subroutine now does four
    !wp! steps (c-x-y-z) wp, july 2000
    !mk! changed the order for proper processing of the interfaces.
    !mk! feb-2001
    !-------------------------------------------------------------

    use GO            , only : GO_Timer_Start, GO_Timer_End
    use GO            , only : TDate
    use dims          , only : statusr => status
    use dims          , only : okdebug
    use dims          , only : splitorderzoom, n_operators
    !use redgridzoom
#ifndef without_chemistry
    use chemistry,     only : Chemistry_Step
#endif
    use sources_sinks, only : source1,source2
    use budget_global, only : budget_transportg,diagbudg,apply_budget_global
    use toolbox,       only : escape_tm
#ifdef MPI
    use mpi_comm,      only : stopmpi
    use mpi_const,     only : previous_par
#endif
    use ParTools     , only : myid, root
    use ParTools     , only : ntracetloc, ntracet_ar
    use ParTools     , only : Par_Check_Domain, Par_Check_Mass
    use global_data,   only : mass_dat
    use MeteoData,     only : m_dat
    use chem_param,    only : names, ntracet   !CMKTEMP
    use convection,    only : convec
    use advectx,    only : advectxzoom, put_xedges
    use advecty,    only : advectyzoom, put_yedges
    use advectz,    only : advectzzoom
    use user_output,only : user_output_step
    use global_data, only : mass_dat
    use dims,        only : isr, ier, jsr, jer, lm, im, jm
    use chem_param, only  : names, ntracet
    use os_specs,   only : WRITE_STR_LEN

    implicit none

    ! --- in/out -----------------------------------

    integer,intent(in)        ::  region
    type(TDate), intent(in)   ::  tr(2)
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/do_steps'

    ! --- local -----------------------------------

    integer          :: child, i123, ichild, tref_child, reg, rgi
    character        :: tobedone
    character(len=1) :: next_step, prev_step   !cmk
    real,dimension(:,:,:,:), pointer :: rm
    real,dimension(:,:,:), pointer   :: m
    real             :: tot_before(0:ntracet), tot_after(0:ntracet), tot_before_sr(0:ntracet), tot_after_sr(0:ntracet)
    integer          :: itr
    character(len=WRITE_STR_LEN) :: write_string

    ! --- begin -----------------------------------

    !write (gol,'(a," : entering")') trim(rname) ; call goPr
    if ( okdebug .and. myid==root ) print *,'do_steps: region ',region

    !example (cmk)
    !   XYZ          EC CEX
    !      xyzeccezyx      cexyzzyxec
    ! in this case the operations e and c should be executed in
    ! the interface region but not in the core of the zoom region
    ! (otherwise double counting)
    ! This results in the most simple algorithm, because you should leave
    ! the DO_STEPS routine AFTER a Z OR after finishing all steps...

    ! note that the parent is responsible for the interface with the children
    ! this means that the interface may not be 'update' in case of (e.g.)
    ! xyz        a
    !    xyzaazyx
    ! in this case, the edges of the child are not updated with the proces a.
    ! this has consequences for the mmix and save output....
    ! THE INTERFACE IS PART OF THE PARENT, NOT THE CHILD.....

    rm => mass_dat(region)%rm_t
    m  => m_dat(region)%data

    do i123=1,n_operators  !WP! changed 3 to n_operators

      next_step = splitorderzoom(region,statusr(region)+1)
      if ( statusr(region) /= 0 ) then
        prev_step = splitorderzoom(region,statusr(region))
      else
        prev_step = ' '
      end if

      if ( okdebug .and. myid == root ) then

        ! it's only to make work of the code visible -- step to be done
        ! will be printed with capital letter (X,Y or Z)
        do reg=1,region
          tobedone = ' '
          if ( reg == region ) tobedone = upper(next_step)
          print *, 'do_steps: ',reg,': ', &
               splitorderzoom(reg,1:statusr(reg)),tobedone
        end do

      end if ! okdebug

      if (print_masses .and. myid == root) then
        do itr = 1, ntracet
            tot_before_sr(itr) = sum(rm(isr(region):ier(region), jsr(region):jer(region), 1:lm(region), itr))
            tot_before(itr)    = sum(rm(1:im(region), 1:jm(region), 1:lm(region), itr))
        end do
        tot_before_sr(0) = sum(m(isr(region):ier(region), jsr(region):jer(region), 1:lm(region)))
        tot_before(0)    = sum(m(1:lm(region), 1:jm(region), 1:lm(region)))
      end if

      select case(next_step)

        case ( 'v' )

#ifndef without_convection
!          if ( convec_apply ) then

            !WP! data must be on tracers for convection
            call Par_Check_domain(region,'n','tracer')
            if ( okdebug ) call Par_Check_mass(region,'bef_convec')

#ifdef with_budgets
            if (apply_budget_global) then
               call budget_transportg(region,0,'conv',prev_step)
            endif
#endif
            if ( okdebug ) print*, 'do_steps: Call convection'

            call Convec( region, next_step, tr, status )
            IF_NOTOK_RETURN(status=1)

#ifdef with_budgets
            if (apply_budget_global) then
               call budget_transportg(region,1,'conv',prev_step)
            endif
#endif

            if ( okdebug ) call Par_Check_mass(region,'aft_convec')

!          end if
#endif

        case ( 'd' )

#ifndef without_diffusion
!          if ( vdiff_apply ) then

            !WP! data must be on tracers for convection
            call Par_Check_domain(region,'n','tracer')
            if ( okdebug ) call Par_Check_mass(region,'bef_vdiff')

#ifdef with_budgets
            if (apply_budget_global) then
               call budget_transportg(region,0,'vdif',prev_step)
            endif
#endif
            if ( okdebug ) print*, 'do_steps: Call vertical diffusion'

            call Convec( region, next_step, tr, status )
            IF_NOTOK_RETURN(status=1)

#ifdef with_budgets
            if (apply_budget_global) then
               call budget_transportg(region,1,'vdif',prev_step)
            endif
#endif
            if ( okdebug ) call Par_Check_mass(region,'aft_vdiff')

!          end if
#endif

        case ( 'c' )

#ifndef without_chemistry
          if ( chemistry_apply ) then

            call GO_Timer_Start( itim_chemie, status )
            IF_NOTOK_RETURN(status=1)
            !WP! data must be on levels for chemistry
            call Par_Check_domain(region,'n','levels')
            if ( okdebug ) call Par_Check_mass(region,'bef_chem')

            if ( okdebug ) print*, 'do_steps: Call Chemistry_Step'
            call Chemistry_Step( region, tr, status )
            IF_NOTOK_RETURN(status=1)

            if ( okdebug ) call Par_Check_mass(region,'aft_chem')

            call GO_Timer_End( itim_chemie, status )
            IF_NOTOK_RETURN(status=1)

          end if
#endif

        case( 's' )

          if ( source_apply ) then

            call GO_Timer_Start( itim_source_sink, status )
            IF_NOTOK_RETURN(status=1)

            !WP! data must be on tracers for convection
            call Par_Check_domain(region,'n','tracer')
            if ( okdebug ) call Par_Check_mass(region,'bef_source')

            if ( okdebug ) print*, 'do_steps: Call source'

            call source1( region, tr, status )
            IF_NOTOK_RETURN(status=1)

            call source2(region)

            if ( okdebug ) call Par_Check_mass(region,'aft_source')

            call GO_Timer_End( itim_source_sink, status )
            IF_NOTOK_RETURN(status=1)

          end if

        case( 'x' )

          call GO_Timer_Start( itim_advectx, status )
          IF_NOTOK_RETURN(status=1)

          !WP! data must be on tracers for X
          call Par_Check_domain(region,'c','tracer')
          if ( okdebug ) call Par_Check_mass(region,'befX_tracers')

#ifdef with_budgets
          if (apply_budget_global) then
             call budget_transportg(region,0,'advx',prev_step)
          endif
#endif

          call advectxzoom( region, status )
          IF_NOTOK_RETURN(status=1)

#ifdef with_budgets
          if (apply_budget_global) then
             call budget_transportg(region,1,'advx',prev_step)
          endif
#endif

          if ( okdebug ) call Par_Check_mass(region,'aftX_tracers')

          call GO_Timer_End( itim_advectx, status )
          IF_NOTOK_RETURN(status=1)

        case( 'y' )

          call GO_Timer_Start( itim_advecty, status )
          IF_NOTOK_RETURN(status=1)

          !WP! data must be on tracers for Y
          call Par_Check_domain(region,'c','tracer')
          if ( okdebug ) call Par_Check_mass(region,'befY_tracers')

#ifdef with_budgets
          if (apply_budget_global) then
             call budget_transportg(region,0,'advy',prev_step)
          endif
#endif

          call advectyzoom( region, status )
          IF_NOTOK_RETURN(status=1)

#ifdef with_budgets
          if (apply_budget_global) then
             call budget_transportg(region,1,'advy',prev_step)
          endif
#endif

          if ( okdebug ) call Par_Check_mass(region,'aftY_tracers')

          call GO_Timer_End( itim_advecty, status )
          IF_NOTOK_RETURN(status=1)

        case( 'z' )

          call GO_Timer_Start( itim_advectz, status )
          IF_NOTOK_RETURN(status=1)

          !WP! data must be on tracers for Z
          call Par_Check_domain(region,'n','tracer')
          if ( okdebug ) call Par_Check_mass(region,'befZ_tracers')

#ifdef with_budgets
          if (apply_budget_global) then
             call budget_transportg(region,0,'advz',prev_step)
          endif
#endif

          call advectzzoom( region, status )
          IF_NOTOK_RETURN(status=1)

#ifdef with_budgets
          if (apply_budget_global) then
             call budget_transportg(region,1,'advz',prev_step)
          endif
#endif

          if ( okdebug ) call Par_Check_mass(region,'aftZ_tracers')

          call GO_Timer_End( itim_advectz, status )
          IF_NOTOK_RETURN(status=1)

        case default

          print *,'do_steps:  strange value in splitorderzoom: ',  &
               splitorderzoom(region,statusr(region))
          print *,'do_steps:  (must be c, e, v, x, y or z)'
          print *,'do_steps:  (found ', next_step, ' instead)'
          call escape_tm('do_steps: Error')

      end select

      if (print_masses .and. myid == root) then
        do itr = 1, ntracet
            tot_after_sr(itr) = sum(rm(isr(region):ier(region), jsr(region):jer(region), 1:lm(region), itr))
            tot_after(itr)    = sum(rm(1:im(region), 1:jm(region), 1:lm(region), itr))
        end do

        tot_after_sr(0) = sum(m(isr(region):ier(region), jsr(region):jer(region), 1:lm(region)))
        tot_after(0)    = sum(m(1:im(region), 1:jm(region), 1:lm(region)))

        ! now print
        if (next_step == print_mass_step .or. print_mass_step == 'a') then
            write(gol,'(a, " :: air mass (im:jm) in region ", i1, " changed in process ", a1, " by ", f25.8)') &
                rname, region, next_step, tot_after(0)-tot_before(0) ; call goPr
            write(gol,'(a, " :: air mass (isr:jsr) in region ", i1, " changed in process ", a1, " by ", f25.8)') &
                rname, region, next_step, tot_after_sr(0)-tot_before_sr(0) ; call goPr

            do itr = 1, ntracet
                write(gol, '(a, " :: tracer ", a, " mass (im:jm) in region ", i1, " changed in process ", a1, " by ", f25.8)') &
                    rname, trim(names(itr)), region, next_step, tot_after(itr)-tot_before(itr) ; call goPr
                write(gol, '(a, " :: tracer ", a, " mass (isr:jsr) in region ", i1, " changed in process ", a1, " by ", f25.8)') &
                    rname, trim(names(itr)), region, next_step, tot_after_sr(itr)-tot_before_sr(itr) ; call goPr
            end do
        end if

      end if

      ! apply mask ?
      if ( mask_apply ) then
        ! mask concentration in a region:
        call MaskOut( region, status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! update region processing status:
      statusr(region) = statusr(region)+1

      ! flexible sample scheme ('a' means after all steps), o = old 'default',
      ! output_after_step = v,c,e,x,y,z,a
      if ((output_after_step == next_step) .or. (output_after_step == 'a')) then
        call GO_Timer_Start( itim_output, status )
        IF_NOTOK_RETURN(status=1)
        call user_output_step( region, tr, status)
        IF_NOTOK_RETURN(status=1)
        call GO_Timer_End( itim_output, status )
        IF_NOTOK_RETURN(status=1)
      endif

      ! e.g. xyzvec...........cevzyx
      ! exit after xyz or vec of cevzyx....
      ! if the exit is after vec...the interface of child should be updated...
      !
      if ( mod(statusr(region),n_operators) == 0 ) then
        ! after xyz....vec...update child interface CMK removed feb 2003...
        ! BUT replaced in march 2003 by MK: budgets were corrupted when
        ! the edges are NOT updated!!!  still to find out why???
        if ( next_step /= 'x' ) then
          call GO_Timer_Start( itim_put_edges, status )
          IF_NOTOK_RETURN(status=1)
          call put_xedges( region, status )
          IF_NOTOK_RETURN(status=1)
          call put_yedges( region, status )
          IF_NOTOK_RETURN(status=1)
          call GO_Timer_End( itim_put_edges, status )
          IF_NOTOK_RETURN(status=1)
        end if    ! if not .... zyxcev
        exit ! e.g. vceecvzyx
      end if


      !exit after 'yz' sequence to proces xyz children
      if ( next_step == 'z' ) then
         prev_step = splitorderzoom(region,statusr(region)-1)
         if ( prev_step == 'y' ) exit
      end if

    end do

    nullify(rm, m)

    !write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine do_steps


  ! ***


  character function upper(xyz)

    implicit none

    character(1),intent(in) :: xyz

    if (xyz=='x') then
       upper = 'X'
    else if (xyz=='y') then
       upper = 'Y'
    else if (xyz=='z') then
       upper = 'Z'
    else if (xyz=='c') then
       upper = 'C'
    else if (xyz=='s') then
       upper = 'S'
    else if (xyz=='v') then
       upper = 'V'
    else
       upper = '_'
    end if

  end function upper


  subroutine CheckMass( region )

    ! Routine to check whether mass converted to surface pressure
    ! converted back to mass returns indeed the same mass
    ! JFM  Feb 2007

    use binas, only       : grav
    use dims, only        : at, bt, im, jm, lm, idate
    use global_data, only : mass_dat, region_dat
    use MeteoData  , only : m_dat

    implicit none

    integer, intent(in)   :: region


    ! local
    real,dimension(:,:,:),pointer    :: m
    real,dimension(:,:,:),allocatable    :: phlb, mnew, mold
    real,dimension(:,:),allocatable      :: ps
    real,dimension(:),pointer        :: dxyp
    integer                          :: i,j,l
    character(len=80)                :: fname

    m    => m_dat(region)%data
    dxyp => region_dat(region)%dxyp

    allocate( ps(im(region),jm(region)) )
    allocate( phlb(im(region),jm(region),lm(region)+1) )
    allocate( mnew(im(region),jm(region),lm(region)) )
    allocate( mold(im(region),jm(region),lm(region)) )

    ps = 0.0

    do l=1,lm(region)
       do j=1,jm(region)
          do i=1,im(region)
             ps(i,j) = ps(i,j) + m(i,j,l)*grav/dxyp(j)
          end do
       end do
    end do

    do l=1,lm(region)+1
       do j=1,jm(region)
          do  i=1,im(region)
             phlb(i,j,l) = at(l)+bt(l)*ps(i,j)
          end do
       end do
    end do

    !
    !cmk ----
    !     compute m (kg), the mass of air in each box.  (at the poles, m
    !     is the air mass of a full cylindrical grid box. This same mass is
    !     placed in every cell for j=1 or j=jm)
    !----

    do l=1,lm(region)
       do j=1,jm(region)
          do  i=1,im(region)
             mnew(i,j,l)=(phlb(i,j,l)-phlb(i,j,l+1))*dxyp(j)/grav
          end do
       end do
    end do

    write(*,*) 'CheckMass: abs. difference ', region, &
         sum( abs( mnew - m(1:im(region),1:jm(region),1:lm(region)) ) )

    write( fname, '("m",i4,5i2.2,"_",i2.2)' ) idate, region
    call WriteAirMass( fname )

    mold = m(1:im(region),1:jm(region),1:lm(region))
    m(1:im(region),1:jm(region),1:lm(region)) = mnew

    write( fname, '("mnew",i4,5i2.2,"_",i2.2)' ) idate, region
    call WriteAirMass( fname )

    m(1:im(region),1:jm(region),1:lm(region)) = mold

    deallocate( ps )
    deallocate( phlb )
    deallocate( mnew )
    deallocate( mold )

    nullify( m )
    nullify( dxyp )

  end subroutine CheckMass


  ! ***


  subroutine WriteAirMass( fname_base )

    ! Writes air mass field to ascii file (one file per region).
    ! Useful for testing whether air mass in forward and adjoint
    ! mode is equivalent.

    use dims,                       only : nregions, im, jm, lm
    use global_data,                only : mass_dat
    use meteoData                 , only : m_dat

    implicit none

    character(len=*), intent(in)    :: fname_base

    integer, parameter              :: kout = 777
    integer                         :: region, i, j, l, imr, jmr, lmr
    integer, dimension(3)           :: lowbound, upbound
    real, dimension(:,:,:), pointer :: m
    character(len=80)               :: fname

    do region = 1, nregions
       imr = im(region)
       jmr = jm(region)
       lmr = lm(region)
       write( fname, '(a,"_",i2.2,".dat")' ) trim(fname_base), region
       open( unit=kout, form='formatted', file=fname )
       m => m_dat(region)%data
       lowbound = lbound(m)
       upbound = ubound(m)
       do i = 1, imr
          do j = 1, jmr
             do l = 1, lmr
                write( kout, '(3i3.3,e20.10)' ) i, j, l, m(i,j,l)
             end do
          end do
       end do
       nullify(m)
       close( kout )
    enddo

  end subroutine WriteAirMass


  ! ***


  subroutine MaskOut( region, status )

    use Grid, only : Region_Apply_Factor
    use dims       , only : im, jm, lm
    use chem_param , only : ntrace, ntracet
    use ParTools   , only : which_par, lmloc
    use global_data, only : mass_dat
    use TM5_Geometry, only : lli

    implicit none

    ! --- in/out ---------------------------------

    integer, intent(in)           ::  region
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MaskOut'

    ! --- local ----------------------------------

    real, dimension(:,:,:,:), pointer   :: rm                   ! tracer mass
#ifdef slopes
    real, dimension(:,:,:,:), pointer   :: rxm, rym, rzm        ! slopes of tracer mass
#ifdef secmom
    real, dimension(:,:,:,:), pointer   :: rxxm, rxym, rxzm, ryym, yzm, rzzm
#endif
#endif
    integer         ::  itr
    integer         ::  imr, jmr, lmr

    ! --- begin ----------------------------------

    ! region horizontal dimension:
    imr = im(region) ; jmr = jm(region)

    ! point to tracer arrays:
    if ( which_par == 'tracer' ) then
      ! all layers:
      lmr = lm(region)
      ! set pointers:
      rm  => mass_dat(region)%rm_t
#ifdef slopes
      rxm => mass_dat(region)%rxm_t
      rym => mass_dat(region)%rym_t
      rzm => mass_dat(region)%rzm_t
#ifdef secmom
      rxxm => mass_dat(region)%rxxm_t
      rxym => mass_dat(region)%rxym_t
      rxzm => mass_dat(region)%rxzm_t
      ryym => mass_dat(region)%ryym_t
      ryzm => mass_dat(region)%ryzm_t
      rzzm => mass_dat(region)%rzzm_t
#endif
#endif
   else
      ! local layers only:
      lmr = lmloc
      ! set pointers:
      rm  => mass_dat(region)%rm_k
#ifdef slopes
      rxm => mass_dat(region)%rxm_k
      rym => mass_dat(region)%rym_k
      rzm => mass_dat(region)%rzm_k
#ifdef secmom
      rxxm => mass_dat(region)%rxxm_k
      rxym => mass_dat(region)%rxym_k
      rxzm => mass_dat(region)%rxzm_k
      ryym => mass_dat(region)%ryym_k
      ryzm => mass_dat(region)%ryzm_k
      rzzm => mass_dat(region)%rzzm_k
#endif
#endif
    end if

    ! loop over transported species:
    do itr = 1, ntracet
      call Region_Apply_Factor( lli(region), rm (1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
#ifdef slopes
      call Region_Apply_Factor( lli(region), rxm(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
      call Region_Apply_Factor( lli(region), rym(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
      call Region_Apply_Factor( lli(region), rzm(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
#ifdef secmom
      call Region_Apply_Factor( lli(region), rxxm(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
      call Region_Apply_Factor( lli(region), rxym(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
      call Region_Apply_Factor( lli(region), rxzm(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
      call Region_Apply_Factor( lli(region), ryym(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
      call Region_Apply_Factor( lli(region), ryzm(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
      call Region_Apply_Factor( lli(region), rzzm(1:imr,1:jmr,1:lmr,itr), mask_llr, mask_factor, status, complement=mask_complement )
      IF_NOTOK_RETURN(status=1)
#endif
#endif
    end do

    ! any non transported tracers ?
    if ( ntrace > ntracet ) then
      write (gol,'("non transported tracers, but chem_dat not defined:")'); call goErr
      write (gol,'("  ntracet     : ",i6)') ntracet; call goErr
      write (gol,'("  ntrace      : ",i6)') ntrace; call goErr
      TRACEBACK; status=1; return
    end if

    ! done:
    nullify( rm )
#ifdef slopes
    nullify( rxm, rym, rzm )
#ifdef secmom
    nullify( rxxm, rxym, rxzm, ryym, ryzm, rzzm )
#endif
#endif

    ! ok
    status = 0

  end subroutine MaskOut


end module modelIntegration
