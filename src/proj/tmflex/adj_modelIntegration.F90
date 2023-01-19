!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_modelIntegration

  !*************************************************
  !* routine to perform main processing and        *
  !* handles time ordering of model processes      *
  !*************************************************
  !

  use go,          only: gol, goErr , goPr

  implicit none

  ! --- in/out -----------------------------

  private

  public :: adj_proces_region
  public :: adj_Proces_Init
  public :: adj_Proces_Done
  public :: adj_Proces_Setup

  public :: output_after_step


  ! --- const -------------------------------

  character(len=*), parameter         :: mname = 'module adj_modelIntegration'


  ! --- var ---------------------------------

  character(len=1)                    :: output_after_step

  logical    ::  advection_apply
  logical    ::  chemistry_apply
  logical    ::  source_apply
!  logical    ::  convec_apply
!  logical    ::  vdiff_apply

contains

  ! =====================================================


  subroutine adj_Proces_Init( status )

    use GO, only : ReadRc
    use global_data, only : rcF

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
    use Meteo              , only : Set
    use Meteo              , only : Meteo_Alloc

    use adj_Advect         , only : adj_Advect_Init
#ifndef without_convection
    use Convection         , only : Convection_Init
#ifndef without_diffusion
    use Diffusion          , only : Diffusion_Init
#endif
#endif
#ifndef without_chemistry
    use adj_Chemistry,   only : adj_Chemistry_Init
#endif
#ifndef without_dry_deposition
    use Dry_Deposition,  only : Dry_Deposition_Init
#endif

    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- local ---------------------------------

    character(len=*), parameter      :: rname = mname//', adj_Proces_Init'

    ! --- local ---------------------------------

    integer           ::  n

    ! --- begin ---------------------------------

    !
    ! which processes ?
    !

    call ReadRc( rcF, 'proces.advection', advection_apply, status )
    IF_NOTOK_RETURN(status=1)
!    call ReadRc( rcF, 'proces.convec'   , convec_apply, status    )
!    IF_NOTOK_RETURN(status=1)
!    call ReadRc( rcF, 'proces.vdiff'    , vdiff_apply  , status, default=convec_apply   )
!    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'proces.chemistry', chemistry_apply, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'proces.source'   , source_apply, status    )
    IF_NOTOK_RETURN(status=1)

    !
    ! init processes
    !

    ! NOTE: these inits should select the appropriate meteo

    call adj_Advect_Init( status )
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
!#ifndef without_diffusion
#ifndef without_chemistry
    ! init chemistry:
    call adj_Chemistry_Init( status )
    IF_NOTOK_RETURN(status=1)
#endif
!    ! init diffusion:
#ifndef without_dry_deposition
    ! dinit dry deposition:
    call Dry_Deposition_Init( status )
    IF_NOTOK_RETURN(status=1)
#endif


!#ifdef NEWMET
    !
    ! allocate meteo
    !

    !>>> should be done by processes in future!
    ! loop over zoom regions:
    do n = 1, nregions

!      ! enable meteo:
!      call Set( oro_dat(n), status, used=.true. )
!      IF_NOTOK_RETURN(status=1)
!
!      ! no lsm, only needed glbsfc for diffusion

      ! enable meteo:
      call Set( sp1_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      call Set( sp2_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      call Set( sp_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call Set( tsp_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call Set( spm_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call Set( phlb_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call Set( m_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call Set( mfu_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call Set( mfv_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call Set( mfw_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call Set( pu_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call Set( pv_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call Set( pw_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call Set( temper_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call Set( humid_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! enable meteo:
      call Set( gph_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! not used yet ...
      !! enable meteo:
      !call Set( lwc_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call Set( iwc_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call Set( cc_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call Set( ccu_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)
      !! enable meteo:
      !call Set( cco_dat(n), status, used=.true. )
      !IF_NOTOK_RETURN(status=1)

      ! convec fields enabled in Convection_Init ...

      ! enable meteo:
      call Set( blh_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

    end do

    ! enable meteo:
    call Set( oro_dat(iglbsfc), status, used=.true. )
    IF_NOTOK_RETURN(status=1)

    !<<<

    ! some meteo implies other ...
    do n = 1, nregions_all
      ! I want to use pclim for CO, and therefore I set pclim_dat(:)%used = .true. from within emission_fwd_CO/emission_fwd_init.
      ! Trouble is, that's called from tracer_init, after which tracer_model is called. Tracer_model calls meteo_init at each
      ! time step, which resets pclim_dat(:)%used to .false. Therefore, I need to set it back to .true. here explicitly. The same
      ! problem might occur for met fields needed for dry deposition, since they're set to .true. in dry_deposition_init. We'll
      ! cross that bridge when we get to it.
      call Set( pclim_dat(n), status, used=.true. )
      IF_NOTOK_RETURN(status=1)

      ! pclim computed from oro:
      if ( pclim_dat(n)%used ) then
        call Set(    oro_dat(n), status, used=.true. )
      end if
      ! convec requires gph and omega:
      if ( entu_dat(n)%used ) then
        call Set(    gph_dat(n), status, used=.true. )
        call Set(  omega_dat(n), status, used=.true. )
      end if
      ! omega (Pa/s) is computed from mfw (kg/s):
      if ( omega_dat(n)%used ) then
        call Set(    mfw_dat(n), status, used=.true. )
      end if
      ! gph is computed from oro/sp/temper/humid
      if ( gph_dat(n)%used ) then
        call Set(    oro_dat(n), status, used=.true. )
        call Set(     sp_dat(n), status, used=.true. )
        call Set( temper_dat(n), status, used=.true. )
        call Set(  humid_dat(n), status, used=.true. )
      end if
      ! sp is interpolated in time between sp1 and sp2:
      if ( sp_dat(n)%used ) then
        call Set(    sp1_dat(n), status, used=.true. )
        call Set(    sp2_dat(n), status, used=.true. )
      end if
      ! cloud covers and overhead/underfeet cloud covers
      if ( cc_dat(n)%used ) then
        call Set( cco_dat(n), status, used=.true. )
        call Set( ccu_dat(n), status, used=.true. )
      end if
    end do

    ! allocate used meteo
    call Meteo_Alloc( status )
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


    !
    ! done
    !

    ! ok
    status = 0

  end subroutine adj_Proces_Init


  ! ***


  subroutine adj_Proces_Done( status )

    use adj_Advect, only : adj_Advect_Done
#ifndef without_convection
    use Convection    , only : Convection_Done
#ifndef without_diffusion
    use Diffusion     , only : Diffusion_Done
#endif
#endif
!#ifndef without_diffusion
#ifndef without_chemistry
    use adj_Chemistry,   only : adj_Chemistry_Done
#endif
#ifndef without_dry_deposition
    use Dry_Deposition,  only : Dry_Deposition_Done
#endif

    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_Proces_Done'

    ! --- begin ---------------------------------

    ! done with advection:
    call adj_Advect_Done( status )
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
!#ifndef without_diffusion
    !! done with emissions:
    !call adj_Emission_Done( status )
    !IF_NOTOK_RETURN(status=1)
!    ! done with diffusion:
#ifndef without_chemistry
    ! done with chemistry:
    call adj_Chemistry_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif
!    call Diffusion_Done( status )
#ifndef without_dry_deposition
    ! done with dry deposition:
    call Dry_Deposition_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif
!    IF_NOTOK_RETURN(status=1)
    ! ok
    status = 0
!#endif
  end subroutine adj_Proces_Done


  ! ***


  ! called every timestep, setup OH fields etc

  subroutine adj_Proces_Setup( t1, t2, status )

    use GO            , only : TDate
#ifndef without_diffusion
    use Diffusion     , only : Diffusion_Setup
#endif
#ifndef without_dry_deposition
    use dry_deposition, only : Dry_Deposition_Calc
#endif
    use Emission_fwd  , only : Emission_fwd_After_Read
!#ifndef without_chemistry
    !use adj_Chemistry , only : adj_Chemistry_Setup
!#endif
!#ifndef without_dry_deposition
    !use Dry_Deposition, only : Dry_Deposition_Setup
!#endif

    ! --- in/out --------------------------------

    type(TDate), intent(in)   ::  t1, t2
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_Proces_Setup'

    ! --- begin ---------------------------------

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
    !call adj_Chemistry_Setup( t1, t2, status )
    !IF_NOTOK_RETURN(status=1)
!#endif

!#ifndef without_dry_deposition
    !call Dry_Deposition_Setup( status )
    !IF_NOTOK_RETURN(status=1)
!#endif

    ! ok
    status = 0

  end subroutine adj_Proces_Setup

!===========================================================================================================
!===========================================================================================================



  recursive subroutine adj_proces_region( region, tr, status )

    use GO              , only : TDate, IncrDate, rTotal, operator(+), operator(-), wrtgol
    use dims            , only : itaur, tref, ndyn, parent, children, revert
    use adj_user_output , only : adj_user_output_step
    use zoom_tools      , only : update_parent
    use advect_tools    , only : m2phlb, m2phlb1
    use adj_zoom_tools  , only : adj_update_parent

    implicit none

    !__IO____________________________________________________________________

    integer, intent(in)       ::  region
    type(TDate), intent(in)   ::  tr(2)
    integer, intent(out)      ::  status

    !__CONST________________________________________________________________

    character(len=*), parameter  ::  rname = mname//'/adj_proces_region'

    !__LOCAL_VARIABLES_______________________________________________________

    integer             ::  child, i, ichild, tref_, n_children, dtime
    integer             ::  nsec_
    type(TDate)         ::  tr_(2)
    character(len=64)   ::  line

    !__START_SUBROUTINE______________________________________________________

    !! info ...
    !write (line,'("        proces region ",i0)') region
    !call wrtgol( trim(line)//' from ', tr(1), ' to ', tr(2) ); call goPr

    ! determine refinement factor with respect to the parent
    tref_ = tref(region)/tref(parent(region))

    !main time step for all processes. note that all timesteps are default/2
    dtime = revert*ndyn/(2*tref(region))

    if ( region /= 1 ) then
      call adj_update_parent(region)
    end if

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

      call adj_do_steps( region, tr_, 1, status )
      IF_NOTOK_RETURN(status=1)

      ! CALL advect_region for all the children (IF there are any)

      ichild = children(region,0)
      do while (ichild > 0)

        child = children(region,ichild)

        call adj_proces_region( child, tr_, status )
        IF_NOTOK_RETURN(status=1)

        ichild = ichild-1

      end do

      !do the remaining steps if necessary...!
      !if ( mod(status(region),n_operators) /= 0 ) then
      call adj_do_steps( region, tr_, 2, status )
      IF_NOTOK_RETURN(status=1)
      !end if

      itaur(region) = itaur(region) + dtime     !update region time (dtime < 0)
    end do

    if ( region /= 1 ) then
       call m2phlb(region)
    else
       call m2phlb1(1)
    end if

    ! ok:
    status = 0

  end subroutine adj_proces_region

!===========================================================================================================
!===========================================================================================================


  subroutine adj_do_steps( region, tr, pass, status )

    use GO                  , only : TDate
    use dims,                 only: okdebug
    use dims,                 only: n_operators
    use dims,                 only: splitorderzoom
    use dims,                 only: statusr => status
    !mkadj
    use dims, only: itaur     , im, jm, itau

    use var4d,                only : steps_region
#ifndef without_chemistry
    use adj_chemistry,      only : adj_chemistry_step
#endif
    use adj_sources_sinks,  only : adj_source1, adj_source2
    use toolbox,        only : escape_tm
    use convection,     only : convec
    use adj_advectx, only : adj_advectxzoom
    use adj_advecty, only : adj_advectyzoom
    use adj_advectz, only : adj_advectzzoom
    !mkadj
    use global_data, only : mass_dat, wind_dat
    use datetime, only : tau2date
    use adj_user_output,  only : adj_user_output_step

    implicit none

    !__IO____________________________________________________________________

    integer, intent(in)       :: region
    type(TDate), intent(in)   :: tr(2)
    integer, intent(in)       :: pass    ! indicates the first or second pass
    integer, intent(out)      :: status

    !__CONST________________________________________________________________

    character(len=*), parameter  ::  rname = mname//'/adj_do_steps'

    !__LOCAL_VARIABLES_______________________________________________________
    integer          :: child, i123, ichild, tref_child, reg, rgi, j
    character        :: tobedone
    character(len=1) :: next_step, prev_step   !cmk

    !mkadj
    integer, dimension(6) :: idate_f


    !__START_SUBROUTINE______________________________________________________


    if (okdebug) print *,' adj_do_steps region/pass ',region, pass

    ! cmk.....for the 1D version
    ! the sequence in forward direction can be somethin like:
    !1 x                  cs  | scx
    !2  x      csscx               scx      x      cs
    !3   xcsscx     scxxcs            scxxcs xcsscx
    !
    ! with 2 regions:
    !1 x      cs | scx
    !2  xsccsx        scxxsc

    ! in the adjoint version we start from the right.
    ! Thus in this example, the first call:
    !    1. quits when is start with an x
    !    2. processes sc and then quits (thus also quits when x)
    ! The second call, consequently
    !    1. processes the whole xcs sequence (exits when 's')
    !    2. processes only 'x' (exit after 'x', if mod(status(region),n_operators) == 0)
    ! implementation : see below...

    noper: DO i123=1,n_operators  !WP! changed 3 to n_operators

      next_step = splitorderzoom( region, statusr(region) )
      if(statusr(region) /= steps_region(region) ) then
        prev_step = splitorderzoom( region, statusr(region)+1 )
      else
        prev_step = ' '
      endif

      if(pass == 1 .and. (next_step == 'z'.or.next_step == 'x')) exit noper

      IF (okdebug) THEN
        ! it's only to make work of the code visible -- step to be done
        ! will be printed with capital letter (X,Y or Z)
        DO reg=1,region
          tobedone = ' '
          IF (reg==region) tobedone = upper(next_step)
          print *,reg,': ',tobedone, splitorderzoom(reg,statusr(reg)+1:steps_region(reg))
          !print *,reg,': ',splitorderzoom(reg,1:status(reg)),tobedone
        ENDDO
      ENDIF

      tobedone = upper(next_step)

      if ((output_after_step == next_step) .or. (output_after_step == 'a')) then
         call adj_user_output_step( region, tr, status )
         IF_NOTOK_RETURN(status=1)
      end if

      select case(next_step)

        case('x')
           if ( advection_apply ) then
              call adj_advectxzoom( region, status )
              IF_NOTOK_RETURN(status=1)
           end if
        case('y')
           if ( advection_apply ) then
              call adj_advectyzoom( region, status )
              IF_NOTOK_RETURN(status=1)
           end if
        case('z')
           if ( advection_apply ) then
              call adj_advectzzoom( region, status )
              IF_NOTOK_RETURN(status=1)
           end if
        case('v')
#ifndef without_convection
           !if ( convec_apply ) then
              ! matrix will be transposed 'on the fly' depending on revert=-1
              call Convec( region, next_step, tr, status )
              IF_NOTOK_RETURN(status=1)
          !end if
#endif
        case('d')
#ifndef without_diffusion
           !if ( vdiff_apply ) then
              ! matrix will be transposed 'on the fly' depending on revert=-1
              call Convec( region, next_step, tr, status )
              IF_NOTOK_RETURN(status=1)
          !end if
#endif
        case('s')
          if ( source_apply ) then
            call adj_source2(region)
            call adj_source1( region, tr, status )
            IF_NOTOK_RETURN(status=1)
          end if

        case('c')
#ifndef without_chemistry
           if ( chemistry_apply ) then
              call adj_chemistry_step( region, tr, status )
              IF_NOTOK_RETURN(status=1)
           endif
#endif
        case default
          print *,'do_steps:  strange value in splitorderzoom: ',  splitorderzoom(region,statusr(region))
          print *,'do_steps:  (must be c, e, v, x, y or z)'
          call escape_tm('do_steps: Error')
       end select

 !cmkadj     write(2,'(i4,a4,f12.9)'), region,next_step, mass_dat(3)%m_t(10,10,1)*1e-12
      statusr(region) = statusr(region)-1
      if (pass == 2 .and. next_step == 'x' .and. mod(statusr(region),n_operators) == 0 ) exit noper

    ENDDO noper

    ! ok:
    status = 0

  end subroutine adj_do_steps

!===========================================================================================================
!===========================================================================================================

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

!===========================================================================================================
!===========================================================================================================

end module adj_modelIntegration

!===========================================================================================================
!===========================================================================================================

