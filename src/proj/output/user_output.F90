!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module user_output
  !--------------------------------------------------------------------------
  ! contains calls to user-specific output routines, e.g.
  ! instantaneous mix files, station output, output of flight tracks etc.
  !--------------------------------------------------------------------------

  use go,          only: gol, goErr


  implicit none

  private

  public :: user_output_init
  public :: user_output_step
  public :: user_output_done

  character(len=*), parameter         :: mname = 'user_output'

  ! ------------ Settings ------------

  logical           :: model_data     ! signal for model output
  logical           :: satellite_data ! flag for satellite output
  logical           :: flask_data     ! signal for flask data, for Wouter's subroutine
  logical           :: station_data   ! signal for writing station time series
  logical           :: tccon_data     ! flag for writing total column time series over certain points
  logical           :: mix_data       ! flag for writing out the 3D mixing ratio fields
  logical           :: column_data    ! flag for writing out the 2D column average mixing ratio fields
  logical           :: flux1x1_data   ! flag for writing 1x1 fluxes

  ! ---------- End settings ----------

contains

!===========================================================================================================
!===========================================================================================================

  subroutine user_output_init(status)
    !-------------------------------------------------------------
    ! user_output_init:
    !      Initialise user-specified model output (all regions)
    !-------------------------------------------------------------

    use GO,                    only : ReadRc, RcHasKey
    use global_data,           only : rcF
    use dims,                  only : revert
    use User_Output_Model    , only : User_Output_Model_Init
    use user_output_satellite, only : user_output_satellite_init
    use user_output_flask,     only : user_output_flask_init
    use user_output_station,   only : user_output_station_init
    use user_output_tccon,     only : user_output_tccon_init
    use user_output_mix,       only : user_output_mix_init
    use user_output_column,    only : user_output_column_init
    use user_output_flux1x1,   only : user_output_flux1x1_init

    implicit none

    !__IO____________________________________________________________________

    integer, intent(out) ::  status


    !__LOCAL_VARIABLES_______________________________________________________

    character(len=*), parameter :: rname = mname//', user_output_init'
    integer                     :: nfound

    !__START_SUBROUTINE______________________________________________________


    status = -1

    model_data      = .false.
    satellite_data  = .false.
    flask_data      = .false.
    station_data    = .false.
    tccon_data      = .false.
    mix_data        = .false.
    column_data     = .false.
    flux1x1_data    = .false.

    if (RcHasKey(rcF, 'model.output', nfound)) then
        call ReadRc( rcF, 'model.output', model_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (RcHasKey(rcF, 'output.point', nfound)) then
        call ReadRc( rcF, 'output.point', flask_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (RcHasKey(rcF, 'output.satellite', nfound)) then
        call ReadRc( rcF, 'output.satellite', satellite_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (RcHasKey(rcF, 'output.station.timeseries', nfound)) then
        call ReadRc( rcF, 'output.station.timeseries', station_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (RcHasKey(rcF, 'output.tccon', nfound)) then
        call ReadRc(rcF, 'output.tccon', tccon_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (RcHasKey(rcF, 'output.mix', nfound)) then
        call ReadRc(rcF, 'output.mix', mix_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (RcHasKey(rcF, 'output.totalcol', nfound)) then
        call ReadRc(rcF, 'output.totalcol', column_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    if (RcHasKey(rcF, 'output.flux1x1', nfound)) then
        call ReadRc(rcF, 'output.flux1x1', flux1x1_data, status)
        IF_NOTOK_RETURN(status=1)
    end if

    ! model data:
    if ( model_data ) then
      call User_Output_Model_Init( status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! satellite data
    if (satellite_data) call user_output_satellite_init(status)
    IF_NOTOK_RETURN(status=1)

    ! flask data
    if (flask_data)     call user_output_flask_init(status)
    IF_NOTOK_RETURN(status=1)

    ! station time series
    if (station_data)   call user_output_station_init(status)
    IF_NOTOK_RETURN(status=1)

    ! tccon stations
    if (tccon_data) call user_output_tccon_init(status)
    IF_NOTOK_RETURN(status=1)

    ! 4D mixing ratio fields
    if (mix_data) call user_output_mix_init(status)
    IF_NOTOK_RETURN(status=1)

    ! Column average mixing ratios
    if (column_data) call user_output_column_init(status)
    IF_NOTOK_RETURN(status=1)

    ! 1x1 fluxes
    if (flux1x1_data) call user_output_flux1x1_init(status)
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine user_output_init

!===========================================================================================================
!===========================================================================================================

  subroutine user_output_step ( region , tr, status )
    !-------------------------------------------------------------
    ! user_output_step:
    !      Define user-specified model output for the region given
    !      Called every time step
    !-------------------------------------------------------------

    use Go,                   only : TDate
    use dims,                 only : itaur
    use datetime,             only : tau2date
    use user_output_flask,    only : user_output_flask_step
    use user_output_station,  only : user_output_station_step
    use user_output_tccon,    only : user_output_tccon_step
    use user_output_mix,      only : user_output_mix_step
    use user_output_column,   only : user_output_column_step
    use user_output_satellite,only : user_output_satellite_step
    use user_output_flux1x1,  only : user_output_flux1x1_step

    implicit none

    !__IO___________________________________________________________________

    integer, intent(inout)          :: status
    integer, intent(in)             :: region
    type(TDate), intent(in)         :: tr(2)

    !__LOCAL_VARIABLES_______________________________________________________

    character(len=*), parameter     :: rname = mname//', user_output_step'

    integer,dimension(6)            :: idate_f
    integer                         :: i


    !__START_SUBROUTINE______________________________________________________

    if (satellite_data) call user_output_satellite_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)
    if (flask_data)     call user_output_flask_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)
    if (station_data)   call user_output_station_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)
    if (tccon_data)     call user_output_tccon_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)
    if (mix_data)       call user_output_mix_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)
    if (column_data)    call user_output_column_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)
    if (flux1x1_data)   call user_output_flux1x1_step(region, tr, status)
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine user_output_step

  subroutine user_output_done ( msg, status )
    !-------------------------------------------------------------
    ! user_output_done:
    !      Finalise user-specified model output for the region given
    !-------------------------------------------------------------

    use dims,                   only : nregions
    use dims,                   only : revert
    use var4d,                  only : iter

    use User_Output_Model     , only : User_Output_Model_Done
    use user_output_satellite,  only : user_output_satellite_done
    use user_output_flask,      only : user_output_flask_done
    use user_output_station,    only : user_output_station_done
    use user_output_tccon,      only : user_output_tccon_done
    use user_output_mix,        only : user_output_mix_done
    use user_output_column,     only : user_output_column_done
    use user_output_flux1x1,    only : user_output_flux1x1_done

    implicit none

    !__IO___________________________________________________________________

    character(len=*),intent(in)  :: msg
    integer, INTENT(OUT)         :: status


    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter      :: rname = mname//', user_output_done'


    !__START_SUBROUTINE______________________________________________________

    status = 0

    ! model data:
    if ( model_data ) then
      call User_Output_Model_Done( status )
      IF_NOTOK_RETURN(status=1)
    end if

    if (satellite_data) call user_output_satellite_done(status)
    IF_NOTOK_RETURN(status=1)
    if (flask_data)     call user_output_flask_done(status)
    IF_NOTOK_RETURN(status=1)
    if (station_data)   call user_output_station_done(status)
    IF_NOTOK_RETURN(status=1)
    if (tccon_data)     call user_output_tccon_done(status)
    IF_NOTOK_RETURN(status=1)
    if (mix_data)       call user_output_mix_done(status)
    IF_NOTOK_RETURN(status=1)
    if (column_data)    call user_output_column_done(status)
    IF_NOTOK_RETURN(status=1)
    if (flux1x1_data)   call user_output_flux1x1_done(status)
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine user_output_done

!===========================================================================================================
!===========================================================================================================

end module user_output
