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

module Emission_Adj_CO2

  use GO, only : gol, goErr, goPr

  use TM5_Fields,   only : T_Fields_4D
  use emission_data, only : adj_emissions

  implicit none

  ! --- in/out -----------------------------

  private

  ! public routines
  public :: Emission_Adj_Init
  public :: Emission_Adj_Apply
  public :: Emission_Adj_Done

  ! --- const ------------------------------

  character(len=*), parameter        :: mname = 'Emission_Adj'

contains

  !---------------------------------------------------------
  ! allocate adjoint emission fields
  !---------------------------------------------------------

  subroutine Emission_Adj_Init( status )

    ! --- modules ------------------------------

    ! --- in/out ----------------------------------------------

    integer, intent(out)             :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/Emission_Adj_Init'

    !--- local ------------------------------------------------

    !--- begin ------------------------------------------------

    ! ok
    status = 0

  end subroutine Emission_Adj_Init


  !---------------------------------------------------------
  ! Finish adjoint emission application
  !---------------------------------------------------------

  subroutine Emission_Adj_Done( status )

    ! --- modules ---------------------------------------------

    ! --- in/out ----------------------------------------------

    integer, intent(out)             :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/Emission_Adj_Done'

    !--- local ------------------------------------------------

    !--- begin ------------------------------------------------

    ! ok
    status = 0

  end subroutine Emission_Adj_Done

  !---------------------------------------------------------
  ! o Apply adjoint emissions
  !---------------------------------------------------------

  subroutine Emission_Adj_Apply( region, tr, status )

    ! --- modules ------------------------------

    use GO                     , only : TDate, NewDate, Get
    use GO                     , only : Time_Profile_Index
    use GO                     , only : operator(+), operator(-), operator(/), rTotal
    use dims,                    only : tref, idatee, itaue, sec_day
    use dims,                    only : im, jm, isr, ier, jsr, jer, lm
    use dims,                    only : itau, itaur, ndyn, itaui
    use dims,                    only : adv_scheme, newsrun
    use global_data,             only : mass_dat, region_dat
    use datetime,                only : tau2date, get_num_days
    use Emission_Data          , only : tracers_em_info
    use Emission_Read_PyShell  , only : read_dailycycle, close_dailycycle

    use Emission_Data          , only : ref_emissions, ref_emissions_apri
    use chem_param,              only : ico2

    ! --- in/out ----------------------------------------------

    integer, intent(in)              :: region
    type(TDate), intent(in)          :: tr(2)
    integer, intent(out)             :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/Emission_Adj_Apply'

    !--- local ------------------------------------------------

    integer             ::  i, j, l
    integer, pointer    ::  zoomed(:,:)
    real, pointer       ::  adj_rm(:,:,:,:)
    real, pointer       ::  adj_rzm(:,:,:,:)
    real, pointer       ::  adj_em(:,:,:,:)
!    real, allocatable   ::  ref_em_apri(:,:)
!    real, allocatable   ::  ref_em(:,:)
!    real, allocatable   ::  frac(:,:)
    real, allocatable   ::  emw(:,:)
    real                ::  x
    type(TDate)         ::  tread
    integer             ::  idater(6)
    real                ::  dtime
    integer             ::  year, month, day
    integer             ::  iw
    integer             ::  i_cat, i_day, idate_mid(6)
    integer             ::  i_period, i_period_dcycle
    logical             ::  apply_dailycycle
    integer             ::  dailycycle_type

    !--- begin ------------------------------------------------

    ! check ...
    select case ( trim(adv_scheme) )
      case ( 'slope' )
#ifndef slopes
        write (gol,'("adv_scheme `",a,"` while macro slopes not defined")') trim(adv_scheme); call goErr
        TRACEBACK; stop
#endif
      case default
        write (gol,'("adv_scheme `",a,"` not supported")') trim(adv_scheme); call goErr
        TRACEBACK; stop
    end select

    ! short:
    zoomed  => region_dat(region)%zoomed
    adj_rm  => mass_dat(region)%rm_t
    adj_rzm => mass_dat(region)%rzm_t

    ! storage:
!    allocate( ref_em_apri(im(region),jm(region)) )
!    allocate( ref_em     (im(region),jm(region)) )
    allocate( emw        (im(region),jm(region)) )
!    allocate( frac       (im(region),jm(region)) )

    ! timestep emissions (should be tr(1)-tr(2), but the 'abs' takes care of that)
    dtime = abs(rTotal( tr(2) - tr(1), 'sec' ))

    ! mid time:
    call tau2date( itau - nint( dtime ), idater )
    tread = NewDate( time6=idater )

    ! current time values:
    call Get( tread, year=year, month=month, day=day )

    ! read the daily cycle, if necessary
    apply_dailycycle = tracers_em_info(region)%tracer(ico2)%dailycycle%apply
    dailycycle_type = tracers_em_info(region)%tracer(ico2)%dailycycle%dtype

    if (apply_dailycycle) then
        ! find out the day number
        i_day = get_num_days(itaui, itaur(region) - ndyn/4/tref(region))
        call tau2date(itaur(region) - ndyn/4/tref(region), idate_mid)
        if ( .not. tracers_em_info(region)%tracer(ico2)%dailycycle%day_opened(i_day) ) then
            ! close_dailycycle only needs the date to write the emissions applied, in case verbose_debug_output is defined
            ! in that case, we should pass a date that's a little before itaur(region), otherwise the emission applied on day 1
            ! gets written to a file with a filename corresponding to day 2
            if ( .not. newsrun ) call close_dailycycle(region, ico2)
            call read_dailycycle(i_day, idate_mid, region, ico2, status)
            IF_NOTOK_RETURN(status=1)
        end if
    end if

    ! ~ loop over categories:
    do i_cat = 1, tracers_em_info(region)%tracer(ico2)%n_cat

      ! short:
      adj_em  => adj_emissions(region)%tracer(ico2)%cat(i_cat)%field   !  kg/cl/s or factor around 0.0 or 1.0

      ! calculate the time step index for daily cycle
      apply_dailycycle = tracers_em_info(region)%tracer(ico2)%dailycycle%cycle_cat(i_cat)%apply
      if (apply_dailycycle) i_period_dcycle = (idate_mid(4) * 3600 + idate_mid(5) * 60 + idate_mid(6)) / &
        tracers_em_info(region)%tracer(ico2)%dailycycle%cycle_cat(i_cat)%dtime + 1 ! This points to the time index within one day

      ! get period index:
      call Time_Profile_Index( tracers_em_info(region)%tracer(ico2)%cat(i_cat)%time_profile, tread, i_period, status, at_left_side=.true. )
!      print *, 'Emission adjoint:', region, i_cat, i_period
      IF_NOTOK_RETURN(status=1)

      !
      ! The emission weight contains the factors that should be
      ! multplied with the emission field in the state to give
      ! the release rate:
      !
      !  c(t+dt) = c(t) + em emw(t) dt
      !

      ! copy emisison fields:
!      ref_em_apri = ref_emissions_apri(region)%cat(i_cat)%field(:,:,1,i_period)   !  kg/cl/s
!      ref_em      = ref_emissions     (region)%cat(i_cat)%field(:,:,1,i_period)   !  kg/cl/s

!      ! weighted with in-month profiles ?
!      if ( associated(Em_nday_Info(year,i_cat)%p) ) then
!        ! current weight index:
!        iw = Em_nday_Info(year,i_cat)%p%iw(year,month,day)
!        ! daily average is fraction ~1.0 of monthly average:
!        frac = Em_nday_Info(year,i_cat)%p%w(iw)%distr(region)%field(:,:,1)
!        ! change from monthly average to daily average:
!        ref_em_apri(:,:) = ref_em_apri(:,:) * frac  ! kg/cl/s
!        ref_em     (:,:) = ref_em     (:,:) * frac  ! kg/cl/s
!      end if

      ! em contains real emissions in kg/cl/s
      !   c(t+dt) = c(t) + em dt
      emw = 1.0
      !
      ! Emission operator:
      !
      !  c(t+dt) = c(t) + emfac(t) dt
      !          ~ c(t) + em ew(t) dt     ! linearized
      !
      ! where :
      !   c   : tracer mass [kg cl-1]
      !   em  : emission field [kg cl-1 s-1]
      !   ew  : emission weight for a specific time [1]
      !   dt  : time step (s)
      !
      ! In matrix form:
      !
      !  [ c  ]   [ I  ew*dt ] [ c  ]
      !  [    ] = [          ] [    ]
      !  [ em ]   [ O    I   ] [ em ]
      !
      ! Adjoint:
      !
      !  [ adj_c  ]   [   I   O ] [ adj_c  ]
      !  [        ] = [         ] [    ]
      !  [ adj_em ]   [ ew^dt I ] [ adj_em ]
      !
      ! where:
      !   adj_c   : adjoint concentration [cl kg-1]
      !   adj_em  : adjoint emission field [cl s kg-1]
      !
      ! which gives the adjoint operator:
      !
      !   adj_em =  adj_c ew dt  +  adj_em
      !

         do j = jsr(region), jer(region)
           do i = isr(region), ier(region)

             if ( zoomed(i,j) /= region ) cycle

             if (apply_dailycycle) then
                select case (dailycycle_type)
                    case (0)
                        x = emw(i,j) * tracers_em_info(region)%tracer(ico2)%dailycycle%cycle_cat(i_cat)%scaling(i,j,i_period_dcycle) * dtime
                    case (1)
                        x = emw(i,j) * dtime
                end select
             else
                x = emw(i,j) * dtime
             end if

             adj_em(i,j,1,i_period) = adj_em(i,j,1,i_period) + adj_rm (i,j,1,ico2) * x
#ifdef slopes
             adj_em(i,j,1,i_period) = adj_em(i,j,1,i_period) - adj_rzm(i,j,1,ico2) * x
#endif

           end do ! i
         end do ! j

      nullify( adj_em  )

    end do ! i_cat

    ! unlink:
    nullify( zoomed )
    nullify( adj_rm  )
    nullify( adj_rzm )

    ! clear:
    deallocate( emw         )

    ! ok:
    status = 0

  end subroutine Emission_Adj_Apply

end module Emission_Adj_CO2
