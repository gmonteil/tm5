!###############################################################################
!
! Parallelisation stuff.
! Dummy variables and routines for single processor case.
!
! npes                       1                  3
!
! myid                       0            0            1            2
!
! root                       0            0            0            0
! root_k                     0            0            0            0
! root_t                     0            0            0            0
!
! lm                        25           25           25           25
! lmloc                     25            8            8            9
! lmar(0:npes-1)          (/25/)      (/8,8,9/)    (/8,8,9/)    (/8,8,9/)
! offsetl                    0            0            8           16
!
! ntrace                    42           42           42           42
!
! ntracet                   26           26           26           26
! ntracetloc                26            8            9            9
! ntracet_ar(0:npes-1)    (/26/)      (/8,9,9/)    (/8,9,9/)    (/8,9,9/)
! offsetn                    0            0            8           17
!
! tracer_active(1:26)       T  1        T  1          F           F   
! tracer_loc(1:26)          :  :        :             :           :
!                           T  8        T  8          F           F 
!                           T  9        F             T  9        F 
!                           :  :        :             :  :        :
!                           T 17        F             T 17        F
!                           T 18        F             F           T 18
!                           :  :        :             :           :  :
!                           T 26        F             F           T 26
!
!
!
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

module ParTools

  use GO        , only : gol, goPr, goErr
#ifdef MPI
  use mpi_const, only : npes, myid
  use mpi_const, only : root, root_k, root_t
  use mpi_const, only : tracer_active, tracer_loc, ntracetloc, ntracet_ar
  use mpi_const, only : lmloc, lmar
  use mpi_const, only : which_par, previous_par
#else
  use dims      , only : nregions
#endif 
  use chem_param, only: ntracet
  
  implicit none

  ! --- in/out -----------------------------------
  
  private
  
  public  ::  npes, myid
  public  ::  procname
  public  ::  root, root_k, root_t
  public  ::  lmloc, lmar, offsetl
  public  ::  ntracetloc, ntracet_ar, offsetn
  public  ::  tracer_loc, tracer_active
  public  ::  which_par, previous_par

  public  ::  Par_Init, Par_Done
  public  ::  Par_Barrier
  public  ::  Par_Check_Domain
  public  ::  Par_Check_Mass
  public  ::  Par_StopMPI
  
  public  ::  Par_Broadcast_Status
  public  ::  Par_Broadcast

  
  ! --- const --------------------------------------

  character(len=*), parameter  ::  mname = 'ParTools'
  
  
  ! --- var --------------------------------------

#ifdef MPI
#else

  ! same as used from mpi const:
  integer           ::  npes
  integer           ::  myid     ! PE number in mpi_comm_world (always one if not MPI)
  integer           ::  root     ! myid of root in mpi_comm_world
  integer           ::  root_k   ! myid of root in com_lev
  integer           ::  root_t   ! myid root in com_trac

  ! nr of levels at this PE
  integer               ::  lmloc

  ! number of levels actually assigned to each PE
  integer, allocatable  ::  lmar(:)
  
  ! nr of tracers and transported tracers at this PE
  integer               ::  ntracetloc    

  ! nr of transported tracers to each PE
  integer, allocatable  ::  ntracet_ar(:)

  ! tracer_active  : determines whether tracer is active on processer 
  logical           ::  tracer_active(ntracet)

  ! tracer_loc  : determines location in the local array 
  integer           ::  tracer_loc(ntracet)

  ! either 'levels' or 'tracer'
  character(len=6)        ::  which_par
  character(len=6)        ::  previous_par(nregions)
 
#endif

  ! level offset: lglob = offsetl+lloc
  integer           ::  offsetl

  ! chemical offset: itracer_global = offsetn + itracer_local
  integer           ::  offsetn

  ! character keys for each processor
  character(len=6)  ::  procname

  
  ! --- interfaces -----------------------------------
  
  interface Par_Broadcast
    module procedure Par_Broadcast_i
    module procedure Par_Broadcast_s
    module procedure Par_Broadcast_r2
    module procedure Par_Broadcast_r3
  end interface

  
contains


  ! ===================================================
  
  !
  ! Inititialisation.
  ! Copied from main program.
  !
  
  subroutine Par_Init( status )
  
    use Dims, only : lm
#ifdef MPI
    use mpi_comm, only : startmpi
#endif

    ! --- in/out ---------------------------------
    
    integer, intent(out)          ::  status
    
    ! --- local ----------------------------------
    
    integer        ::  n
    
    ! --- begin ----------------------------------
    
#ifdef MPI

    call startmpi

#else

    ! single processor
    npes = 1

    ! dummy id for root processor
    root = 0

    ! only one processor, so this is always root ...
    myid = root
    root_k = root
    root_t = root 

    ! *** levels

    ! allocate arrays
    allocate( lmar(0:npes-1) )

    ! nr of levels at this PE
    ! (take levels from region 1)
    lmloc = lm(1)

    ! number of levels actually assigned to each PE
    ! (take levels from region 1)
    lmar(0) = lm(1)
  
    ! *** tracers

    ! allocate arrays
    allocate( ntracet_ar(0:npes-1) )

    ! tracer_active  : determines whether tracer is active on processer 
    tracer_active = .true.
 
    ! number of transported tracers
    ntracetloc = ntracet
 
    ! tracer location
    tracer_loc = -9
    do n = 1, ntracet
      tracer_loc(n) = n
    end do

    ! ***

    ! initially 'parallel' over tracers:
    which_par = 'tracer'
    previous_par(:) = 'tracer'

#endif

    ! level offset: lglob = offsetl+lloc
    offsetl = 0
    if ( myid > 0 ) offsetl = sum(lmar(0:myid-1))

    ! tracer offset: itracer_global = offsetn + itracer_local
    offsetn = 0
    if ( myid > 0 ) offsetn = sum(ntracet_ar(0:myid-1))

    ! set processor names: pe0000, pe0001, ...
    write (procname,'("pe",i4.4)') myid

    ! ok
    status = 0

  end subroutine Par_Init
  
  
  ! ***
  
  
  subroutine Par_Done( status )
  
#ifdef MPI
    use mpi_comm, only : stopmpi
#endif
    
    ! --- in/out ------------------------------
    
    integer, intent(out)          ::  status
    
    ! --- begin ------------------------------

#ifdef MPI
    
    !CMK ???? call stopmpi
    
#else
  
    ! deallocate arrays
    deallocate( lmar       )
    deallocate( ntracet_ar )

#endif
    
    ! ok
    status = 0

  end subroutine Par_Done
  

  ! ***
  
  
  subroutine Par_Barrier

#ifdef MPI
    use mpi_comm, only : barrier
#endif

    ! --- begin ------------------------------------
      
#ifdef MPI   
    call barrier
#endif

  end subroutine Par_Barrier
    


  ! ***
  
  
  subroutine Par_Check_Domain( region1, others, to_check )

#ifdef MPI
    use mpi_comm, only : check_domain    
#endif

    ! --- in/out -----------------------------------
    
    integer,intent(in)              ::  region1 
    character(len=*),intent(in)     ::  others
    character(len=*),intent(in)     ::  to_check

    ! --- begin ------------------------------------
      
#ifdef MPI   

    call check_domain( region1, others, to_check )

#endif

  end subroutine Par_Check_Domain


  ! ***
  

  subroutine Par_Check_Mass( region, text )

#ifdef MPI
    use mpi_comm, only : check_mass
#endif

    ! --- in/out -----------------------------------
    
    integer, intent(in)             ::  region
    character(len=*), intent(in)    ::  text

    ! --- start ------------------------------------

#ifdef MPI
    call check_mass( region, text )
#endif

  end subroutine Par_Check_Mass


  ! ***
  
  
  subroutine Par_StopMPI

#ifdef MPI
    use mpi_comm, only : stopmpi
#endif

    ! --- begin ------------------------------------
      
#ifdef MPI   

    call stopmpi

#endif

  end subroutine Par_StopMPI


  ! **************************************************************
  

  subroutine Par_Broadcast_Status( istat, id )
  
#ifdef MPI
    use mpi_const, only : MPI_INTEGER, localComm
#endif

    ! --- in/out -------------------------------------
    
    integer, intent(inout)            ::  istat
    integer, intent(in)               ::  id
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Par_Broadcast_Status'
    
    ! --- local ------------------------------------
    
    integer          ::  status
    
    ! --- begin ------------------------------------
    
    ! send the input status to all other processes:
    call Par_Broadcast( istat, id, status )
    if (status/=0) then
      write (gol,'("broadcasting status")'); call goErr
      TRACEBACK; istat=1; return
    end if
    
    ! each process has now return status 'istat' ...
    
  end subroutine Par_Broadcast_Status


  ! ***
  

  subroutine Par_Broadcast_i( i, id, status )
  
#ifdef MPI
    use mpi_const, only : MPI_INTEGER, localComm
#endif

    ! --- in/out -------------------------------------
    
    integer, intent(inout)            ::  i
    integer, intent(in)               ::  id
    integer, intent(out)              ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Par_Broadcast_i'
    
    ! --- begin ------------------------------------

    if ( npes > 1 ) then
#ifdef MPI
      call MPI_BCast( i, 1, MPI_INTEGER, id, localComm, status )
      IF_NOTOK_RETURN(status=1)
#else
      write (gol,'("please implement for non-mpi parallel library.")'); call goErr
      TRACEBACK; status=1; return
#endif
    end if
    
    ! ok
    status = 0
    
  end subroutine Par_Broadcast_i


  ! ***
  

  subroutine Par_Broadcast_s( s, id, status )
  
#ifdef MPI
    use mpi_const, only : MPI_CHARACTER, localComm
#endif

    ! --- in/out -------------------------------------
    
    character(len=*), intent(inout)   ::  s
    integer, intent(in)               ::  id
    integer, intent(out)              ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Par_Broadcast_s'
    
    ! --- begin ------------------------------------

    if ( npes > 1 ) then
#ifdef MPI
      call MPI_BCast( s, len(s), MPI_CHARACTER, id, localComm, status )
      if (status/=0) then; write (*,'("ERROR from MPI_BCAST ",a)') rname; status=1; return; end if
#else
      write (gol,'("please implement for non-mpi parallel library.")'); call goErr
      TRACEBACK; status=1; return
#endif
    end if
    
    ! ok
    status = 0
    
  end subroutine Par_Broadcast_s


  ! ***
  

  subroutine Par_Broadcast_r2( x, id, status )
  
#ifdef MPI
    use mpi_const, only : MPI_REAL8, localComm
#endif

    ! --- in/out -------------------------------------
    
    real(8), intent(inout)   ::  x(:,:)
    integer, intent(in)      ::  id
    integer, intent(out)     ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Par_Broadcast_r2'
    
    ! --- begin ------------------------------------

    if ( npes > 1 ) then
#ifdef MPI
      call mpi_bcast( x, size(x), MPI_REAL8, id, localComm, status )
      if (status/=0) then; TRACEBACK; status=1; return; end if
#else
      write (gol,'("please implement for non-mpi parallel library.")'); call goErr
      TRACEBACK; status=1; return
#endif
    end if
    
    ! synchronize ...
    call Par_Barrier
    
    ! ok
    status = 0
    
  end subroutine Par_Broadcast_r2


  ! ***
  
  
  subroutine Par_Broadcast_r3( x, id, status )
  
#ifdef MPI
    use mpi_const, only : MPI_REAL8, localComm
#endif

    ! --- in/out -------------------------------------
    
    real(8), intent(inout)   ::  x(:,:,:)
    integer, intent(in)      ::  id
    integer, intent(out)     ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Par_Broadcast_r3'
    
    ! --- begin ------------------------------------

    if ( npes > 1 ) then
#ifdef MPI
      call mpi_bcast( x, size(x), MPI_REAL8, id, localComm, status )
      if (status/=0) then; TRACEBACK; status=1; return; end if
#else
      write (gol,'("please implement for non-mpi parallel library.")'); call goErr
      TRACEBACK; status=1; return
#endif
    end if
    
    ! synchronize ...
    call Par_Barrier
    
    ! ok
    status = 0
    
  end subroutine Par_Broadcast_r3



end module ParTools
