!###############################################################################
!
! this module declares the values needed in MPI communications
! WP january 2003
! 
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################

module mpi_const

  use dims, only : nregions
  use chem_param, only: ntracet

  implicit none

  include 'mpif.h'

  integer :: my_real ! platform dependent reference to real values for MPI 
  integer :: myid    ! PE number in mpi_comm_world
  integer :: npes    ! number of PE's
  integer :: pe_first_tracer ! lowest myid involved in processes over tracers
  integer :: pe_first_l ! lowest myid involved in processes over levels
  integer :: ierr       ! return status of MPI routine calls
  integer :: com_trac   ! communicator with only PE's having nonzero ntracetloc
  integer :: com_lev    ! communicator with inly PE's having nonzero lmloc
  integer :: myid_t     ! PE number in com_trac 
                        ! (can differ from mpi_comm_world!)
  integer :: myid_k     ! PE number in com_lev
  integer :: root       ! myid of root in mpi_comm_world 
  integer :: root_k     ! myid of root in com_lev 
  integer :: root_t     ! myid root in com_trac
  character(len=6)        :: which_par  ! either 'levels' or 'tracer'
  ! previous_par  : previous paralel regime
  character(len=6),dimension(nregions) :: previous_par

  !integer,dimension(0:npes-1) :: lmar  ! number of levels assigned to each PE
  !integer,dimension(0:npes-1) :: ntracet_ar ! nr of transported tracers "  "
  integer, allocatable :: lmar(:)  ! number of levels assigned to each PE
  integer, allocatable :: ntracet_ar(:) ! nr of transported tracers "  "

  integer :: lmloc         ! nr of levels at this PE
  integer :: ntracetloc    ! nr of tracers and transported tracers at this PE
  ! tracer_active  : determines whether tracer is active on processer 
  logical,dimension(ntracet) :: tracer_active
  ! tracer_loc  : determines location in the local array 
  integer,dimension(ntracet) :: tracer_loc
  ! proc_tracer  : determines processor that handles itracer
  integer,dimension(ntracet) :: proc_tracer
  ! allocate_mass  : switch to allocate and deallocate mass after each swap
  logical,parameter          :: allocate_mass=.false.

end module mpi_const
