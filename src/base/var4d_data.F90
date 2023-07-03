!###############################################################################
!
! Module containing data and data types that are 'global' in 4DVAR routines
! and general 4DVAR subroutines
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

module Var4D_Data

  use dims,                     only : nregions
  use GO,                       only : gol, goErr, goPr, goLabel

  implicit none

  ! all public:
  public
  

  ! -------- CONSTANTS --------------------

  integer, parameter               :: RUN_FORWARD = 1
  integer, parameter               :: RUN_BACKWARD = 2

  !integer, parameter             :: interpolation_method = 1    !gridbox_value
  !integer, parameter             :: interpolation_method = 2   !slopes
  integer, parameter               :: interpolation_method = 3   !linear


  ! -------- TYPES -------------------------------

  ! measurements of different types
  type data_4DVAR
     integer                    :: otype     ! 1 = station, c = column
     integer                    :: index     ! for station data only: station index
     real                       :: lat       ! latitude
     real                       :: lon       ! longitude 
     real                       :: height    ! height of station / surface elevation for SCIAMACHY pixels
     integer                    :: region    ! region of measurement
     integer                    :: is
     integer                    :: js
     integer                    :: ls        ! ls = -1 for column observation
     real                       :: y         ! value of observation
     real                       :: y_std     ! standard deviation of observation
     integer                    :: count     ! number of obs used for averaging
     real                       :: dy_obs    ! uncertainty of observation
     real                       :: dy_mod    ! uncertainty of corresponding model value
                                             ! (representativeness error)
     real                       :: dy_modh   ! 'horizontal representativeness error'
     real                       :: dy_modv   ! 'vertical representativeness error'
     real                       :: dy_tot    ! total uncertainty for data point (=obs+mod error)
     real                       :: y_M       ! model result corresponding to observation (Hx)
     real                       :: y_M_std   ! standard deviation of model result (within assimilation window)
     real                       :: res       ! residual (Hx - y)
     real                       :: res1      ! E-1 x (Hx - y) 
     integer                    :: ncount    ! number of averaging steps in this timeslot

  end type data_4DVAR


  ! -------- VARIABLES --------------------

  integer,dimension(nregions)      :: steps_region = 0   !ADJ

  character(len=5)                 :: runid  ! identifier for run
  integer                          :: iter
  integer                          :: calc_NL = 0   ! for emissions....non-linear..
  ! APRIF   a priori forward run (iter = 0)
  ! APOSF   a posteriori forward run (iter = -1)
  ! 0001F   forward run, iteration 1
  ! 0001A   adjoint run, iteration 1




  integer                :: run_mode                ! run mode

  integer                :: nasim                  ! length 4DVAR time slot

  real                   :: J_tot                  ! cost function - total
  real                   :: J_obs                  ! cost function - observations
  real                   :: J_para                 ! cost function - parameters


end module Var4D_Data
