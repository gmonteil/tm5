!###############################################################################
!
! this module holds all specific routines needed to run MPI
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

module mpi_comm

  use mpi_const

  implicit none

  private

  public :: swap_all_mass
  public :: check_mass
  public :: startmpi, stopmpi, abortmpi
  public :: check_domain
  public :: barrier_t, barrier_k, barrier
  public :: scatter_after_read_k_3d, scatter_after_read_k
  public :: scatter_after_read_t
  public :: gather_tracer_k_3d, gather_tracer_k
  public :: gather_tracer_t

contains


  subroutine swap_all_mass(region,iaction)
    !
    !WP! subroutine swaps the arrays m,rm,rxm,rym,and rzm
    !WP! from levels to tracers (iaction=1)
    !WP! or from levels to tracers (iaction=0)
    !WP! for requested regions (region1=>region2)
    !WP! after swapping, old data is set to zero to avoid mistakes
    !WP! with 'old' data
    !WP! 23 december 2002

    use dims,       only : okdebug, im,jm ,lm
    use chem_param, only : ntracet
    use global_data,only : mass_dat
    implicit none

    ! input/output
    integer,intent(in) :: iaction
    integer,intent(in) :: region

    ! local
    integer :: imr,jmr,lmr
    real,dimension(:,:,:),pointer     :: msource
    real,dimension(:,:,:),pointer     :: mtarget
    real,dimension(:,:,:,:),pointer     :: rxmsource
    real,dimension(:,:,:,:),pointer     :: rxmtarget
    integer :: nsend,xx, xx_end

    ! start
    imr=im(region)
    jmr=jm(region)
    lmr=lm(region)

    if ( iaction == 0 ) then

       msource => mass_dat(region)%m_t
       mtarget => mass_dat(region)%m_k
       !WP! the 3d array mass is correct on root_t, and scattered to levels
       call scatter_after_read_k_3d( msource, imr,jmr,lmr, 2,2,0, &
                                     mtarget, root )
       call barrier
       nullify(msource)
       nullify(mtarget)

    else if ( iaction == 1 ) then

       !WP! gather mass to root_k, and broadcast to all
       !WP! only PE's in com_trac need mloc_tracer,
       !WP! but we cannot guarantee that root_k is
       !WP! also in com_trac, whcih would lead to errors.

       msource => mass_dat(region)%m_k
       mtarget => mass_dat(region)%m_t

       call  gather_tracer_k_3d( mtarget, imr,jmr,lmr, 2,2,0, &
                                 msource, .true. ) !WP! gather and broadcast

       nullify(msource)
       nullify(mtarget)

    end if
    if ( okdebug .and. myid == root ) print *,'swap_all_mass: Scattered mass'

#ifdef secmom
    xx_end = 10
#else
    xx_end=4
#endif

    do xx=1,xx_end   !WP! once for each x,y, and z
       select case(xx)
       case(1)
          rxmsource => mass_dat(region)%rm_t
          rxmtarget => mass_dat(region)%rm_k
       case(2)
          rxmsource => mass_dat(region)%rxm_t
          rxmtarget => mass_dat(region)%rxm_k
       case(3)
          rxmsource => mass_dat(region)%rym_t
          rxmtarget => mass_dat(region)%rym_k
       case(4)
          rxmsource => mass_dat(region)%rzm_t
          rxmtarget => mass_dat(region)%rzm_k
#ifdef with_mpi
       case(5)
          rxmsource => mass_dat(region)%rxxm_t
          rxmtarget => mass_dat(region)%rxxm_k
       case(6)
          rxmsource => mass_dat(region)%rxym_t
          rxmtarget => mass_dat(region)%rxym_k
       case(7)
          rxmsource => mass_dat(region)%ryym_t
          rxmtarget => mass_dat(region)%ryym_k
       case(8)
          rxmsource => mass_dat(region)%rxzm_t
          rxmtarget => mass_dat(region)%rxzm_k
       case(9)
          rxmsource => mass_dat(region)%ryzm_t
          rxmtarget => mass_dat(region)%ryzm_k
       case(10)
          rxmsource => mass_dat(region)%rzzm_t
          rxmtarget => mass_dat(region)%rzzm_k
#endif
       end select
       if ( iaction == 0 ) then
          call tracer_to_k(rxmtarget,&
               rxmsource, &
               imr,jmr,lmr,&
               2,2,0,&
               ntracet,ntracetloc,&
               ntracet_ar)

          call barrier
       else if ( iaction == 1 ) then

          call k_to_tracer(rxmtarget,&
               rxmsource, &
               imr,jmr,lmr,&
               2,2,0,&
               ntracet,ntracetloc,&
               ntracet_ar)
       end if
       nullify(rxmsource)
       nullify(rxmtarget)

       if ( okdebug .and. myid == root ) &
            print*,'swap_all_mass: Scattered tracer mass and slopes ',xx

    end do

  end subroutine swap_all_mass



  subroutine tracer_to_k( array_k, array_tracer, im, jm, lm, &
       ih, jh, lh, ntracer, ntracerloc, tracerar )
    !
    ! re-distribute an array which is distributed over the k-index
    ! to an array which is distributed over the tracers
    ! uses MPI_ALLTOALLV and a help array
    ! we have to use MPI_ALLTOALLV: the number of data received from the PEs
    ! can be different from PE to PE
    !
    implicit none

    ! input/output
    integer,intent(in) :: im, jm, lm
    integer,intent(in) :: ih, jh, lh
    integer,intent(in) :: ntracer, ntracerloc
    integer,dimension(0:npes-1),intent(in)  :: tracerar
    real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracerloc) :: array_tracer
    real,dimension(1-ih:im+ih,1-jh:jm+jh,1   :lmloc,ntracer) :: array_k

    ! local
    real,dimension((im+2*ih)*(jm+2*jh)*(lm+2*lh)*ntracerloc) :: array_tracer_help
    integer,dimension(0:npes-1) :: sdispls
    integer,dimension(0:npes-1) :: rdispls
    integer,dimension(0:npes-1) ::  sendcounts
    integer,dimension(0:npes-1) ::  recvcounts
    integer :: recvcount
    integer :: my_vector
    integer :: slab2d
    integer :: itracer, i, j, k, ii, iproc, kadd, kloop
    logical :: okdebug=.false.

    ! start

    slab2d = (im+2*ih)*(jm+2*jh)
    !
    ! fill necessary arrays
    !
    do i=0,npes-1
       recvcounts(i) = lmloc*slab2d*tracerar(i)
    end do

    rdispls(0)       = 0
    !
    do i=1,npes-1
       rdispls   (i) = rdispls(i-1)+recvcounts(i-1)
    end do
    !
    do j=0,npes-1
       sendcounts(j) = lmar(j)*slab2d*ntracerloc
    end do
    !
    sdispls(0)       = 0
    !
    do j=1,npes-1
       sdispls(j)    = sdispls(j-1)+sendcounts(j-1)
    end do
    !
    !
    ii = 0
    array_tracer_help = 0.
    !
    ! data received from npes PE's, numbered 0,... npes-1
    !
    do iproc  =0,npes-1
       !
       ! data received from each processor for ntracerloc tracers
       !
       !
       do itracer=1,ntracerloc
          !
          ! calculate the GLOBAL k value = k + kadd
          ! by adding the k values already treated
          !
          kadd = 0
          !
          ! calculate k-offset
          !
          do kloop=0,iproc-1
             kadd = kadd + lmar(kloop)
          end do
          !
          do kloop=1,lmar(iproc)
             k = kloop + kadd
             do j=1-jh,jm+jh
                do i=1-ih,im+ih
                   ii = ii + 1
                   array_tracer_help(ii) = array_tracer(i,j,k,itracer)
                   !WP! array_tracer set to zero to avoid mistakes!!!
                   array_tracer(i,j,k,itracer)=0.0
                end do
             end do
          end do
       end do
    end do
    !
    !
    !
    ! fill array_tracer_help
    !
    !if(okdebug)then
    !     write(6,*)'COUNTCHECK t_to_k',myid,ntracerloc,lmloc
    !     write(6,*)'COUNTCHECK t_to_k',myid,tracerar
    !     write(6,*)'COUNTCHECK t_to_k',myid,sdispls
    !     write(6,*)'COUNTCHECK t_to_k',myid,sendcounts
    !     write(6,*)'COUNTCHECK t_to_k',myid,recvcounts
    !     write(6,*)'COUNTCHECK t_to_k',myid,'to send ',sum(array_tracer_help)
    !end if

    array_k=0.
    call MPI_ALLTOALLV( array_tracer_help, sendcounts, sdispls, MY_REAL, &
                        array_k          , recvcounts, rdispls, MY_REAL, &
                        MPI_COMM_WORLD, ierr )

    !
    ! array array_tracer_help contains all the data, which are now
    ! distributed properly over the local tracer array
    !
    ! the structure of array array_tracer_help is now as follows:
    !
    !
    !     PE0            itrace1      lmar(0)          2D slabs  ---
    !     PE0            itrace2      lmar(0)          2D slabs    |
    ! itracer=1,ntracerloc from PE0:
    !     PE0            itrace3      lmar(0)          2D slabs    |
    !         ...                                                  |
    !     PE0            ntracerloc   lmar(0)          2D slabs  ---
    !
    !     PE1            itrace1      lmar(1)          2D slabs  ---
    !     PE1            itrace2      lmar(1)          2D slabs    |
    ! itracer=1,ntracerloc from PE1:
    !     PE1            itrace3      lmar(1)          2D slabs    |
    !         ...                                                  |
    !     PE1            ntracerloc   lmar(0)          2D slabs  ---
    !
    !     ....
    !     ....
    !     ....
    !     ....
    !
    !     PEntracerloc   itrace1      lmar(ntracerloc) 2D slabs  ---
    !     PEntracerloc   itrace2      lmar(ntracerloc) 2D slabs    |
    ! itracer=1,ntracerloc from PEnpes-1
    !     PEntracerloc   itrace3      lmar(ntracerloc) 2D slabs    |
    !         ...                                                  |
    !     PEntracerloc   ntracerloc   lmar(ntracerloc) 2D slabs  ---
    !
  end subroutine tracer_to_k



  subroutine k_to_tracer( array_k, array_tracer, im, jm, lm, &
       ih, jh, lh, ntracer, ntracerloc, tracerar )
    !
    ! re-distribute an array which is distributed over the k-index
    ! to an array which is distributed over the tracers
    ! uses MPI_ALLTOALLV and a help array
    ! we have to use MPI_ALLTOALLV: the number of data received from the PEs
    ! can be different from PE to PE

    implicit none

    ! input/output
    integer,intent(in) :: im,jm,lm,ih,jh,lh,ntracer,ntracerloc
    integer,dimension(0:npes-1),intent(in)        :: tracerar
    real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracerloc), intent(out) :: array_tracer
    real,dimension(1-ih:im+ih,1-jh:jm+jh,1   :lmloc,ntracer), intent(inout) :: array_k

    ! local
    real,dimension((im+2*ih)*(jm+2*jh)*(lm+2*lh)*ntracerloc) :: array_tracer_help
    integer sdispls   (0:npes-1)
    integer rdispls   (0:npes-1)
    integer sendcounts(0:npes-1)
    integer recvcounts(0:npes-1)
    integer recvcount
    integer my_vector
    integer slab2d
    integer itracer, i, j, k, ii, iproc, kadd, kloop

    ! start

    slab2d = (im+2*ih)*(jm+2*jh)
    !
    ! fill necessary arrays
    !
    array_tracer=0.0
    do i=0,npes-1
       sendcounts(i) = lmloc*slab2d*tracerar(i)
    end do

    sdispls(0)       = 0
    !
    do i=1,npes-1
       sdispls   (i) = sdispls(i-1)+sendcounts(i-1)
    end do
    !
    do j=0,npes-1
       recvcounts(j) = lmar(j)*slab2d*ntracerloc
    end do
    !
    rdispls(0)       = 0
    !
    do j=1,npes-1
       rdispls(j)    = rdispls(j-1)+recvcounts(j-1)
    end do

    !
    ! fill array_tracer_help
    !
    array_tracer_help = 0.0
    call MPI_ALLTOALLV( array_k          , sendcounts, sdispls, MY_REAL, &
                        array_tracer_help, recvcounts, rdispls, MY_REAL, &
                        MPI_COMM_WORLD, ierr )
    !
    ! array array_tracer_help contains all the data, which are now distributed properly over
    ! the local tracer array
    !
    ! the structure of array array_tracer_help is now as follows:
    !
    !
    !     PE0            itrace1      lmar(0)          2D slabs  ---
    !     PE0            itrace2      lmar(0)          2D slabs    |
    ! itracer=1,ntracerloc from PE0
    !     PE0            itrace3      lmar(0)          2D slabs    |
    !         ...                                                  |
    !     PE0            ntracerloc   lmar(0)          2D slabs  ---
    !     PE1            itrace1      lmar(1)          2D slabs  ---
    !     PE1            itrace2      lmar(1)          2D slabs    |
    ! itracer=1,ntracerloc from PE1
    !     PE1            itrace3      lmar(1)          2D slabs    |
    !         ...                                                  |
    !     PE1            ntracerloc   lmar(0)          2D slabs  ---
    !     ....
    !     ....
    !     ....
    !     ....
    !     PEntracerloc   itrace1      lmar(ntracerloc) 2D slabs  ---
    !     PEntracerloc   itrace2      lmar(ntracerloc) 2D slabs    |
    ! itracer=1,ntracerloc from PEnpes-1
    !     PEntracerloc   itrace3      lmar(ntracerloc) 2D slabs    |
    !         ...                                                  |
    !     PEntracerloc   ntracerloc   lmar(ntracerloc) 2D slabs  ---
    !
    ii = 0
    !
    ! data received from npes PE's, numbered 0,... npes-1
    !
    do iproc  =0,npes-1
       !
       ! data received from each processor for ntracerloc tracers
       !
       !WP! array_tracer is only filled on PEs where ntracetloc is nonzero

       do itracer=1,ntracerloc
          !
          ! calculate the GLOBAL k value = k + kadd
          ! by adding the k values already treated
          !
          kadd = 0
          !
          ! calculate k-offset
          !
          do kloop=0,iproc-1
             kadd = kadd + lmar(kloop)
          end do
          !
          do kloop=1,lmar(iproc)
             k = kloop + kadd
             do j=1-jh,jm+jh
                do i=1-ih,im+ih
                   ii = ii + 1
                   array_tracer(i,j,k,itracer) = array_tracer_help(ii)
                   !AJS! >>> error in subscript 3 :
                   !AJS! lmar(iproc) for some processor iproc might
                   !AJS! exceed the local value 'lmloc' ...
                   !WP! array_k set to zero to avoid mistakes !!!
                   !array_k(i,j,kloop,itracer)=0.0
                   !AJS! <<< see below
                end do
             end do
          end do   !  kloop

          !AJS! now set to zero to avoid mistakes ...
          array_k(:,:,:,itracer) = 0.0

       end do  ! itracer
    end do  ! iproc
    !
    !
  end subroutine k_to_tracer



  subroutine startmpi

    !WP! subroutine initilizes mpi,
    !WP! sets up mpi_comm_world,
    !WP! checks the number of processors, and
    !WP! proceeds to set the paralel mode of the model to 'tracers'

    implicit none

    ! start

    call MPI_INIT(ierr)
    !WP! my_real comes from mpif.h, other declarations may be needed sometimes
    my_real = MPI_DOUBLE_PRECISION
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, npes, ierr )

    ! allocate arrays
    allocate( lmar      (0:npes-1) )
    allocate( ntracet_ar(0:npes-1) )

    call initialize_domains

    which_par='tracer'
    previous_par(:)='tracer'
    root=0 !WP! setting root in mpi_comm_world to PE 0

    if ( myid == root ) then
       print*,'***************************************************************'
       print*,'***************************************************************'
       print*,'**************   MPI initialized, starting         ************'
       print*,'**************                                     ************'
       print*,'**************       data is now paralel       ****************'
       print*,'**************          over ',which_par,&
            '            ****************'
       print*,'***************************************************************'
       print*,'***************************************************************'
    end if

    call barrier

  end subroutine startmpi



  subroutine stopmpi
    implicit none

    ! deallocate arrays
    deallocate( lmar       )
    deallocate( ntracet_ar )

    call barrier
    call mpi_finalize(ierr)
    STOP 'stopmpi: EINDE'

  end subroutine stopmpi



  subroutine abortmpi
    implicit none

    ! deallocate arrays

    call mpi_abort(mpi_comm_world, 999, ierr)
    STOP 'abortmpi: FORCED quit'

  end subroutine abortmpi




  subroutine escape_tm_mpi(msg)
    !--------------------------------------------------------------
    !
    ! abnormal termination of a model run
    !
    ! msg       character string to be printed on unit kmain
    !
    !-------------------------------------------------------------
    use dims,        only : okdebug, kmain, itau
    use datetime,    only : tstamp
    implicit none
    !
    character(*),intent(in) :: msg
    !
    if ( myid == root ) then
       write(kmain,'(///4(1x,39("--")/))')
       call tstamp(kmain,itau,msg)
       write(kmain,'(///4(1x,39("--")/))')
    end if

    call abortmpi

  end subroutine escape_tm_mpi



  subroutine initialize_domains
    !
    !
    !
    use dims
    use chem_param
    implicit none

    !if ( npes > lm(1) .or. npes > ntracet ) then
    !   print *, 'initialize_domains: ERROR ', &
    !        'npes should be less than ntracet and lm!'
    !   call STOPMPI
    !end if

    lmar    = lm(1)
    call determine_lmar(lmar,npes,lm(1))
    call determine_lmar(ntracet_ar,npes,ntracet)

    lmloc   = lmar(myid)

    call determine_first(PE_FIRST_TRACER,ntracet_ar)
    call determine_first(PE_FIRST_L     ,lmar   )

    ntracetloc   = ntracet_ar(myid)

    print*,'initialize_domains: myid, PE_FIRST_TRACER', myid, &
         PE_FIRST_TRACER, lmar
    print*,'initialize_domains: myid, PE_FIRST_L     ', myid, PE_FIRST_L
    call set_domain('tracer')
    call set_domain('levels')
    if ( root_t /= 0 .or. root_k /= 0 .or. root /= 0 ) then
       if ( myid == 0 ) then
          print*,'initialize_domains: ERROR  root_t /= root_k'
          print*,'initialize_domains: problems for budget calculations '
       end if
       call MPI_FINALIZE(ierr)
       stop 'initialize_domains: ERROR'
    end if

    call determine_tracer_active

  contains

    subroutine determine_lmar(xlmar,npesmax,xlm)
      implicit none
      integer resdiv, remainder, i, count
      integer npesmax
      integer xlm
      integer xlmar(0:npesmax-1)
      !
      ! this subroutine fills elements 0,npes-1 of array lmar
      ! with a number equal to lm div npes = resdiv (integer division)
      ! if resdiv*npes not equal to lm, there is a remainder
      ! the remainder is divided over the last elements
      ! (npes-1, npes-remainder)
      !
      resdiv = xlm/npes
      remainder  = xlm - resdiv*npes

      xlmar       = 0
      do i=0,npes-1
         xlmar(i) = resdiv
      end do
      !
      count = remainder
      ! CMK changed!  PE0 should have less work (surface = more iterations)
      i = npes-1
      do
         if(count == 0)exit
         xlmar(i) = xlmar(i) + 1
         count   = count   - 1
         i       = i       - 1   ! CMK changed!
      end do
      do i=0,npes-1
         print '(a,2i4)' ,'distribution over PE s: i, array(i) ', i, xlmar(i)
      end do

    end subroutine determine_lmar

    subroutine determine_first(pe_first,xlmar)
      !
      ! we assume a distribution over PEs with:
      !
      ! lmar(0)    lmar(1)    lmar(2)     ....   lmar(npes-1)
      !
      !
      implicit none
      integer :: pe_first, xlmar(0:npes-1)
      integer :: iproc

      pe_first=0

      !if(xlmar(npes-1) == 0)then
      !   call barrier
      !   call mpi_finalize(ierr)
      !   STOP 'ERROR: NO DATA on LAST PE'
      !end if
      !
      !iproc = npes-1
      !
      !do while (iproc.ge.0)
      !   if(xlmar(iproc) == 0)exit
      !   iproc = iproc - 1
      !end do
      !pe_first = iproc + 1

    end subroutine determine_first

    subroutine determine_tracer_active
      integer :: n,offsetn,offsetnp,ipos
      offsetn = sum(ntracet_ar(0:myid-1)) + 1
      offsetnp = offsetn + ntracetloc - 1

      do n=1,ntracet
         if (n >= offsetn .and. n <= offsetnp) then
            tracer_active(n) = .true.
            tracer_loc(n)    = n-offsetn+1    ! the counter in the
         else
            tracer_active(n) = .false.
            tracer_loc(n)    = -9           !  for checks....
         end if
      end do

      ipos = 1
      do n=0,npes-1
         ! proc_tracer - which processor handles tracer i ?
         proc_tracer(ipos:ipos+ntracet_ar(n)-1) = n
         ipos = ipos + ntracet_ar(n)
      end do
      print '(a30,a30,i4,60l2)' , 'determine_tracer_active: ', &   !CMK corrected
           'Tracers active for processor ', myid, tracer_active
      print '(a30,a30,i4,60i2)' , 'determine_tracer_active: ', &
           'Local tracer number on proc  ', myid, tracer_loc
      print '(a30,a30,i4,60i2)' , 'determine_tracer_active: ', &
           'Tracer numbers processed by  ', myid, proc_tracer

    end subroutine determine_tracer_active

  end subroutine initialize_domains



  subroutine set_domain(switch)
    use dims
    implicit none

    ! input/output
    character(len=6),intent(in) :: switch
    integer :: orig_group2,orig_group,com_trac_group,com_lev_group
    integer :: a,b,i
    integer :: new_rank,lev_rank
    integer,allocatable,dimension(:) :: inc_ranks

    myid_t=999 !WP! only processors in the communicator receive a myid
    myid_k=999 !WP! only processors in the communicator receive a myid

    !WP! set up communication world for tracers
    if ( switch == 'tracer' ) then
       if ( okdebug ) print*,'set_domain: Setting up communicator ', &
            'world for tracer domain ',myid
       call MPI_COMM_GROUP(mpi_comm_world,orig_group,ierr)
       if ( ntracetloc /= 0 ) then   !WP! only relevant processors
          b=npes-PE_FIRST_TRACER !WP! number of involved PE's
          allocate(inc_ranks(b)) !WP! to hold ranks of PE's
          do i=1,b
             inc_ranks(i)=PE_FIRST_TRACER+i-1 !WP! fill ranks array
          end do
          call MPI_GROUP_INCL(orig_group,b,inc_ranks,com_trac_group,ierr)
          deallocate(inc_ranks)
       end if
       call barrier
       call mpi_comm_create(mpi_comm_world,com_trac_group,com_trac,ierr)
       call barrier
       !WP! ranks in new communicator
       if ( ntracetloc /= 0 ) call MPI_COMM_RANK(com_trac,myid_t,ierr)
       root_t=PE_FIRST_TRACER

       print*,'set_domain: In tracer domain myid: ',myid, &
            ' has number ',myid_t,' and uses as root: ',root_t
       !WP! first free old communciator group
       call MPI_GROUP_FREE(com_trac_group,ierr)
       !WP!
       !WP! com_trac is setup now, proceed with com_lev
       !WP!

    end if
    if ( switch == 'levels' ) then

       if ( okdebug ) print*,'set_domain: Setting up communicator ', &
            'world for levels domain '
       call barrier
       call MPI_COMM_GROUP(mpi_comm_world,orig_group2,ierr)
       if ( lmloc /= 0 ) then
          b=npes-PE_FIRST_L
          allocate(inc_ranks(b))
          do i=1,b
             inc_ranks(i)=PE_FIRST_L+i-1
          end do
          call MPI_GROUP_INCL(orig_group2,b,inc_ranks,com_lev_group,ierr)
          deallocate(inc_ranks)
       end if
       call barrier
       call mpi_comm_create(mpi_comm_world,com_lev_group,com_lev,ierr)
       call barrier
       call MPI_COMM_RANK(com_lev,myid_k,ierr)
       !WP! ranks in new communicator
       if ( lmloc /= 0) call MPI_COMM_RANK(com_lev,myid_k,ierr)
       root_k=PE_FIRST_L
       print*,'set_domain: In levels domain myid: ',myid, &
            ' has number ',myid_k,' and uses as root: ',root_k
       !WP! first free old communciator group
       call MPI_GROUP_FREE(com_lev_group,ierr)
    end if

    !call MPI_COMM_FREE(com_trac,ierr)  !WP! first free old communciator group
    !call MPI_COMM_FREE(com_lev,ierr)  !WP! first free old communciator group

  end subroutine set_domain



  subroutine check_domain(region1,others,to_check)

    use mpi_const,   only : previous_par,allocate_mass
    use dims,        only : parent,children, okdebug
    use global_data, only : free_massdat, assign_massdat

    implicit none

    ! input/output
    integer,intent(in)   ::region1
    character(len=*) :: to_check,others

    ! local
    integer :: reg,regstart,regstop,iaction
    integer :: nchild,regiop

    ! start

    select case(others)
    case('p') !WP! also check parents
       regiop=parent(region1)
       if ( regiop == 0 ) then
          regstart=region1
          regstop=region1
       else if ( regiop /= 0 ) then
          regstart=regiop
          regstop=region1
       end if
    case('c') !WP! also check child
       nchild=children(region1,0)
       if ( nchild == 0 ) then
          regstart=region1
          regstop=region1
       else if ( nchild /= 0 ) then
          regstart=region1
          !WP! take number of last child as regstop
          regstop=children(region1,nchild)
       end if
    case('n') !WP! check no parents or children
       regstart=region1
       regstop=region1
    case default
       call escape_tm_mpi('check_domain: illegal value of *others* ')
    end select

    select case(to_check)
    case('levels')
       iaction=0
    case('tracer')
       iaction=1
    case default
       call escape_tm_mpi('check_domain: illegal value of *to_check* ')
    end select

    do reg=regstart,regstop
       ! make sure data is in requested domain
       if ( previous_par(reg) /= to_check ) then
          if ( okdebug .and. myid == root ) &
               print*,'check_domain: data will be swapped'
          previous_par(reg)=to_check   ! reset previous par
          ! assign new memory
          if ( allocate_mass ) call assign_massdat(reg,iaction)
          call barrier
          call swap_all_mass(reg,iaction)      ! swap data for this region
          call barrier
          ! free previously allocated memory
          if ( allocate_mass ) call free_massdat(reg,iaction)
       else
          if ( okdebug .and. myid == root ) then
             print*,'check_domain: data is on proper communicator, nothing to do'
          end if
       end if
    end do

  end subroutine check_domain



  subroutine barrier
    implicit none

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  end subroutine barrier



  subroutine barrier_t
    implicit none

    call MPI_BARRIER(com_trac,ierr)

  end subroutine barrier_t



  subroutine barrier_k
    implicit none

    call MPI_BARRIER(com_lev,ierr)

  end subroutine barrier_k


  !
  ! gather an array
  !    tracer(1-ih,im+ih,1-jh,jm+jh,1-lh,lm+lh,ntracer)
  ! on root_t
  ! from local arrays
  !    tracerloc(1-ih,im+ih,1-jh,jm+jh,1-lh:lm+lh,ntracerloc)
  ! on all PE's where
  ! the array is distributed over the n-index
  !

  subroutine gather_tracer_t(tracer,im,jm,lm,ih,jh,lh,ntracer,tracerloc,broadcast)

    use dims, only : CheckShape

    ! input/output
    logical,intent(in)                    :: broadcast
    integer,intent(in)                    :: im,jm,lm
    integer,intent(in)                    :: ih,jh,lh
    integer,intent(in)                    :: ntracer

    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracer), &
    !     intent(inout)   :: tracer
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracetloc), &
    !     intent(in)      :: tracerloc
    real,dimension(:,:,:,:), intent(inout)   ::  tracer
    real,dimension(:,:,:,:), intent(in)      ::  tracerloc

    ! local
    integer :: slab3d
    integer,dimension(0:npes-1) :: displ
    integer,dimension(0:npes-1) :: recvcounts
    integer :: sendcount, nsend
    integer :: itracer, i , ierr
    logical :: okdebug=.false.

    ! start

    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh,ntracer   /), shape(tracer   ) )
    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh,ntracetloc/), shape(tracerloc) )

    if ( okdebug ) print*,'Gather_tracer_t: total on pe ',myid, &
         ' = ',sum(tracerloc)
    call barrier
    slab3d    = (im+2*ih)*(jm+2*jh)*(lm+2*lh)
    displ(:)=0  ! offset

    do i=1,npes-1
       displ(i) = displ(i-1) + ntracet_ar(i-1)*slab3d ! calculate offset slabs
    end do
    !
    do i=0,npes-1
       recvcounts(i) = ntracet_ar(i)*slab3d ! receive nr of slabs from other PE's
    end do
    !
    sendcount = ntracetloc*slab3d  ! send the local nr of slabs
    !
    !
    call barrier
    ! gather on PE root_t
    !call MPI_GATHERV( &
    !       tracerloc(1-ih,1-jh,1-lh,1), sendcount        , my_real, &
    !       tracer   (1-ih,1-jh,1-lh,1), recvcounts, displ, my_real, &
    !       root_t, mpi_comm_world, ierr )
    call MPI_GATHERV( &
           tracerloc(:,:,:,:), sendcount        , my_real, &
           tracer   (:,:,:,:), recvcounts, displ, my_real, &
           root_t, mpi_comm_world, ierr )
    if ( ierr /= 0 ) call escape_tm_mpi('gather_tracer_t: ERROR')
    call barrier
    !
    if ( myid==root_t .and. okdebug ) then
       print*,'Gather_tracer_t: total on root ',myid,' = ',sum(tracer)
    end if

    if ( broadcast ) then
       nsend=(im+2*ih)*(jm+2*jh)*(lm+2*lh)*ntracer
       call mpi_bcast( tracer, nsend, my_real, root_t, mpi_comm_world, ierr )
    end if

  end subroutine gather_tracer_t



  !
  ! gather   an array
  !    tracer    (1-ih,im+ih,1-jh,jm+jh,1-lh,lm+lh,ntracer)
  ! on root_k
  ! from local arrays
  !    tracerloc(1-ih,im+ih,1-jh,jm+jh,1-lh,lmloc+lh,ntracer)
  ! on all PE's where
  ! the array is distributed over the l-index
  !

  subroutine gather_tracer_k_3d( tracer, im,jm,lm, ih,jh,lh, tracerloc, broadcast )

    use dims, only : CheckShape

    ! input/output
    logical,intent(in)                        :: broadcast
    integer,intent(in)                        :: im,jm,lm
    integer,intent(in)                        :: ih,jh,lh

    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh), intent(inout)   :: tracer
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1   :lmloc), intent(in)      :: tracerloc
    real, dimension(:,:,:), intent(inout)   ::  tracer
    real, dimension(:,:,:), intent(in)      ::  tracerloc

    ! local
    integer :: slab2d
    integer,dimension(0:npes-1) :: displ
    integer,dimension(0:npes-1) :: recvcounts
    integer :: sendcount, nsend
    integer :: i , ierr
    logical :: okdebug=.false.

    ! start

    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh/), shape(tracer   ) )
    call CheckShape( (/im+2*ih,jm+2*jh,lmloc  /), shape(tracerloc) )

    if ( okdebug ) &
         print*,'Gather_tracer_k: total on pe ',myid,' = ',sum(tracerloc)

    slab2d    = (im+2*ih)*(jm+2*jh)
    displ(:)=0  ! offset
    displ(0) = lh*slab2d  ! offset with lh as these do not exist on tracerloc

    do i=1,npes-1
       displ(i) = displ(i-1) + lmar(i-1)*slab2d ! calculate offset slabs
    end do
    !
    do i=0,npes-1
       recvcounts(i) = lmar(i)*slab2d ! receive nr of slabs from other PE's
    end do
    !
    sendcount = lmloc*slab2d  ! send the local nr of slabs
    !
    !
    call barrier
    ! gather on PE root_k
    !call MPI_GATHERV( tracerloc(1-ih,1-jh,   1), sendcount, my_real, &
    !                  tracer   (1-ih,1-jh,1-lh), recvcounts, displ, my_real, &
    !                  root_k, mpi_comm_world, ierr )
    call MPI_GATHERV( tracerloc(:,:,:), sendcount, my_real, &
                      tracer   (:,:,:), recvcounts, displ, my_real, &
                      root_k, mpi_comm_world, ierr )
    if(ierr/=0) call escape_tm_mpi('gather_tracer_k: error ')
    call barrier
    !
    if ( myid==root_k .and. okdebug ) then
       print*,'Gather_tracer_k: total on root ',myid,' = ',sum(tracer)
    end if

    if ( broadcast ) then
       nsend=(im+2*ih)*(jm+2*jh)*(lm+2*lh)
       call mpi_bcast(tracer,nsend,my_real,root_k,mpi_comm_world,ierr)
    end if

  end subroutine gather_tracer_k_3d


  ! ***


  subroutine gather_tracer_k(tracer,im,jm,lm,ih,jh,lh,ntracer,tracerloc,broadcast)

    use dims, only : CheckShape

    ! input/output
    logical,intent(in)                        :: broadcast
    integer,intent(in)                        :: im,jm,lm
    integer,intent(in)                        :: ih,jh,lh
    integer,intent(in)                        :: ntracer

    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracer), &
    !     intent(inout)   :: tracer
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1   :lmloc,ntracer), &
    !     intent(in)      :: tracerloc
    real, dimension(:,:,:,:), intent(inout)   ::  tracer
    real, dimension(:,:,:,:), intent(in)      ::  tracerloc


    ! local
    integer :: slab2d
    integer,dimension(0:npes-1) :: displ
    integer,dimension(0:npes-1) :: recvcounts
    integer :: sendcount, nsend
    integer :: itracer, i , ierr
    logical :: okdebug=.false.

    ! start

    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh,ntracer/), shape(tracer   ) )
    call CheckShape( (/im+2*ih,jm+2*jh,  lmloc,ntracer/), shape(tracerloc) )

    if ( okdebug ) &
         print*,'Gather_tracer_k: total on pe ',myid,' = ',sum(tracerloc)

    slab2d    = (im+2*ih)*(jm+2*jh)
    displ(:)=0  ! offset
    displ(0) = lh*slab2d  ! offset with lh as these do not exist on tracerloc

    do i=1,npes-1
       displ(i) = displ(i-1) + lmar(i-1)*slab2d ! calculate offset slabs
    end do
    !
    do i=0,npes-1
       recvcounts(i) = lmar(i)*slab2d ! receive nr of slabs from other PE's
    end do
    !
    sendcount = lmloc*slab2d  ! send the local nr of slabs
    !
    !
    call barrier
    do itracer = 1, ntracer
       ! gather on PE root_k
       !call MPI_GATHERV( tracerloc(1-ih,1-jh,   1,itracer), sendcount        , my_real, &
       !                  tracer   (1-ih,1-jh,1-lh,itracer), recvcounts, displ, my_real, &
       !                  root_k,mpi_comm_world, ierr)
       call MPI_GATHERV( tracerloc(:,:,:,itracer), sendcount        , my_real, &
                         tracer   (:,:,:,itracer), recvcounts, displ, my_real, &
                         root_k,mpi_comm_world, ierr)
       if(ierr/=0) call escape_tm_mpi('gather_tracer_k: error ')
    end do
    call barrier
    !
    if ( myid==root_k .and. okdebug ) then
       print*,'Gather_tracer_k: total on root ',myid,' = ',sum(tracer)
    end if

    if ( broadcast ) then
       nsend=(im+2*ih)*(jm+2*jh)*(lm+2*lh)*ntracer
       call mpi_bcast(tracer,nsend,my_real,root_k,mpi_comm_world,ierr)
    end if

  end subroutine gather_tracer_k


  ! *********************************************


  !
  ! scatter an array
  !    tracer   (1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracer)
  ! on root (pe=0)
  ! over local arrays
  !    tracerloc(1-ih:im+ih,1-jh:jm+jh,1:lmloc   ,ntracer)
  ! on all PE's so that
  ! the array is distributed over the l-index
  !

  subroutine scatter_after_read_k_3d( tracer_read, im,jm,lm, ih,jh,lh, &
                                      tracerloc, send_id )

    use mpi_const,only : lmar,myid,lmloc
    use dims, only : CheckShape

    ! input/output
    integer, intent(in)              :: im,jm,lm,ih,jh,lh
    integer, intent(in)              :: send_id
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh) :: tracer_read
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1   :lmloc) :: tracerloc
    real, dimension(:,:,:), intent(in)    :: tracer_read
    real, dimension(:,:,:), intent(inout) :: tracerloc

    ! local
    real, dimension(1-ih:im+ih,1-jh:jm+jh) :: sum_one,sum_all
    integer                         :: slab2d
    integer, dimension(0:npes-1)    :: displ
    integer, dimension(0:npes-1)    :: sendcounts
    integer                         :: recvcount
    integer                         :: nsum
    integer                         :: i

    ! start

    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh/), shape(tracer_read) )
    call CheckShape( (/im+2*ih,jm+2*jh,lmloc  /), shape(tracerloc  ) )

    slab2d    = (im+2*ih)*(jm+2*jh)
    !
    ! offset for gather-scatter
    !
    displ(:) = 0
    displ(0) = lh*slab2d
    do i=1,npes-1
       displ(i) = displ(i-1) + lmar(i-1)*slab2d
    end do
    !
    do i=0,npes-1
       sendcounts(i) = lmar(i)*slab2d
    end do
    !
    recvcount = lmloc*slab2d
    !
    ! call MPI_SCATTERV
    !
    !call MPI_SCATTERV( tracer_read(1-ih,1-jh,1-lh), sendcounts, displ, MY_REAL, &
    !                   tracerloc  (1-ih,1-jh,1   ), recvcount , MY_REAL, &
    !                   send_id, MPI_COMM_WORLD, ierr )
    call MPI_SCATTERV( tracer_read(:,:,:), sendcounts, displ, MY_REAL, &
                       tracerloc  (:,:,:), recvcount        , MY_REAL, &
                       send_id, MPI_COMM_WORLD, ierr )

  end subroutine scatter_after_read_k_3d


  ! ***


  subroutine scatter_after_read_k( tracer_read, im,jm,lm, ih,jh,lh, ntracer, &
                                   tracerloc, send_id )

    use mpi_const,only : lmar,myid,lmloc
    use dims, only : CheckShape

    ! input/output
    integer, intent(in)              :: im,jm,lm,ih,jh,lh,ntracer
    integer, intent(in)              :: send_id
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracer) :: tracer_read
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1   :lmloc,ntracer) :: tracerloc
    real, dimension(:,:,:,:), intent(in)    :: tracer_read
    real, dimension(:,:,:,:), intent(inout) :: tracerloc

    ! local
    real,dimension(1-ih:im+ih,1-jh:jm+jh,ntracer) :: sum_one,sum_all
    integer                         :: slab2d
    integer,dimension(0:npes-1)   :: displ
    integer,dimension(0:npes-1)   :: sendcounts
    integer                         :: recvcount
    integer                         ::nsum
    integer                         :: itracer, i

    ! start

    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh,ntracer/), shape(tracer_read) )
    call CheckShape( (/im+2*ih,jm+2*jh,lmloc  ,ntracer/), shape(tracerloc  ) )

    slab2d    = (im+2*ih)*(jm+2*jh)
    !
    ! offset for gather-scatter
    !
    displ(:) = 0
    displ(0) = lh*slab2d
    do i=1,npes-1
       displ(i) = displ(i-1) + lmar(i-1)*slab2d
    end do
    !
    do i=0,npes-1
       sendcounts(i) = lmar(i)*slab2d
    end do
    !
    recvcount = lmloc*slab2d
    !
    ! call MPI_SCATTERV
    !
    do itracer = 1, ntracer
       !call MPI_SCATTERV( tracer_read(1-ih,1-jh,1-lh,itracer), sendcounts, displ, MY_REAL, &
       !                   tracerloc  (1-ih,1-jh,1   ,itracer), recvcount , MY_REAL, &
       !                   send_id, MPI_COMM_WORLD, ierr )
       call MPI_SCATTERV( tracer_read(:,:,:,itracer), sendcounts, displ, MY_REAL, &
                          tracerloc  (:,:,:,itracer), recvcount        , MY_REAL, &
                          send_id, MPI_COMM_WORLD, ierr )
    end do

  end subroutine scatter_after_read_k


  ! **********************************


  subroutine scatter_after_read_t(tracer_read,im,jm,lm,ih,jh,lh,tracerloc,send_id)

    use mpi_const,only : lmar,myid,lmloc,ntracetloc,ntracet_ar
    use dims, only : CheckShape

    !
    ! scatter an array
    !    tracer    (1-ih,im+ih,1-jh,jm+jh,1-lh,lm+lh,ntracer)
    ! on root (pe=0)
    ! over local arrays
    !    tracerloc(1-ih,im+ih,1-jh,jm+jh,1,  ,lmloc,ntracer)
    ! on all PE's so that
    ! the array is distributed over the l-index
    !
    ! input/output
    integer,intent(in)              :: im,jm,lm,ih,jh,lh
    integer,intent(in)              :: send_id
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracet) :: tracer_read
    !real,dimension(1-ih:im+ih,1-jh:jm+jh,1-lh:lm+lh,ntracetloc) :: tracerloc
    real, dimension(:,:,:,:)        :: tracer_read
    real, dimension(:,:,:,:)        :: tracerloc

    ! local
    integer                         :: slab3d
    integer,dimension(0:npes-1)     :: displ
    integer,dimension(0:npes-1)     :: sendcounts
    integer                         :: recvcount,i
    logical                         :: okdebug=.false.

    ! start

    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh,ntracet   /), shape(tracer_read) )
    call CheckShape( (/im+2*ih,jm+2*jh,lm+2*lh,ntracetloc/), shape(tracerloc  ) )

    slab3d    = (im+2*ih)*(jm+2*jh)*(lm+2*lh)
    !
    ! offset for gather-scatter
    !
    displ(:) = 0
    do i=1,npes-1
       displ(i) = displ(i-1) + ntracet_ar(i-1)*slab3d ! calculate offset slabs
    end do
    !
    do i=0,npes-1
       sendcounts(i) = ntracet_ar(i)*slab3d
    end do
    !
    recvcount = ntracetloc*slab3d
    !
    ! call MPI_SCATTERV
    !
    call barrier
    !call MPI_SCATTERV( tracer_read(1-ih,1-jh,1-lh,1), sendcounts, displ, MY_REAL, &
    !                   tracerloc  (1-ih,1-jh,1-lh,1), recvcount        , MY_REAL, &
    !                   send_id, MPI_COMM_WORLD, ierr )
    call MPI_SCATTERV( tracer_read(:,:,:,:), sendcounts, displ, MY_REAL, &
                       tracerloc  (:,:,:,:), recvcount        , MY_REAL, &
                       send_id, MPI_COMM_WORLD, ierr )
    call barrier
    !

  end subroutine scatter_after_read_t



  subroutine dump_field4d(region,fieldloc4d,lmr,ntracel,add,name,filename)
    !
    ! This subroutine is used to collect 3d arrays from all
    ! the processors and write them to a HDF file.
    ! The array add contains the number of halo cells in each dimension (3).
    ! depending on the size of lmr, data is gathered and written
    ! from root_k, or copied and written from root_t.
    ! If lmloc==lm(region), then npes=0 and gathering is also not necessary.
    ! WP january 2003
    !
    use io_hdf,      only : io_write, DFACC_CREATE, DFNT_INT32
    use dims,        only : im, jm, lm
    use chem_param,  only : ntracet
    use dims       , only : CheckShape

    ! input/output
    character(len=*), intent(in)   :: name,filename
    integer, intent(in)            :: region,lmr,ntracel
    integer, dimension(3)          :: add  !WP! to expand arrays, only one or zero
    !real,dimension(&
    !     1-add(1):im(region)+add(1),&
    !     1-add(2):jm(region)+add(2),&
    !     1:lmr,ntracel)   :: fieldloc4d
    real, dimension(:,:,:,:), intent(in)  ::  fieldloc4d

    ! local
    integer :: is,ie,js,je
    integer :: io,istat,sfstart,sfend,sfsnatt,root_id
    real,dimension(&
         1-add(1):im(region)+add(1),&
         1-add(2):jm(region)+add(2),&
         1-add(3):lm(region)+add(3),ntracet) :: fieldglob4d

    ! start

    call CheckShape( (/im(region)+2*add(1),jm(region)+2*add(2),lmr,ntracel/), shape(fieldloc4d) )

    which_par=previous_par(region)
    fieldglob4d=0.0
    if ( which_par == 'levels' ) then ! parallel over levels
       call  gather_tracer_k(fieldglob4d,im(region),jm(region),lm(region), &
            add(1),add(2),add(3),ntracet, fieldloc4d,.false.)
       ! no need to broadcast result
       root_id=root_k
    else if(which_par=='tracer') then  !WP! parallel over tracers
       call  gather_tracer_t(fieldglob4d,im(region),jm(region),lm(region), &
            add(1),add(2),add(3),ntracet, fieldloc4d,.false.)
       ! no need to broadcast result
       root_id=root_t
    end if

    if ( myid == root_id ) then
       io=SFSTART(filename, DFACC_CREATE)
       istat = sfsnatt(io,'im',    DFNT_INT32, 1, im(region))
       istat = sfsnatt(io,'jm',    DFNT_INT32, 1, jm(region))
       istat = sfsnatt(io,'lm',    DFNT_INT32, 1, lm(region))
       istat = sfsnatt(io,'ntracet',    DFNT_INT32, 1, ntracet)
       call io_write( io, im(region), 'lon', jm(region), 'lat', &
            lm(region),'hybrid',ntracet,'ntracet',&
            fieldglob4d(&
            1:im(region),&
            1:jm(region),&
            1:lm(region),1:ntracet) ,name)
       print *,'dump_field4d: closing output file',SFend(io)
    end if

  end subroutine dump_field4d



  subroutine dump_field3d(region,fieldloc3d,lmr,add,name,filename)
    !
    ! This subroutine is used to collect 3d arrays from all the
    ! processors and write them to a HDF file. The array add
    ! contains the number of halo cells in each dimension (3).
    ! depending on the size of lmr, data is gathered and written
    ! from root_k, or copied and written from root_t.
    ! If lmloc==lm(region), then npes=0 and gathering is also not necessary.
    !    WP january 2003
    !
    use io_hdf,      only : io_write, DFACC_CREATE, DFNT_INT32
    use dims,        only : im, jm, lm
    use dims,        only : CheckShape

    ! input/output
    character(len=*), intent(in)         :: name,filename
    integer, intent(in)                  :: region,lmr
    integer, dimension(3), intent(in)    :: add  !WP! to expand arrays, only one or zero
    !real,dimension(&
    !     1-add(1):im(region)+add(1),&
    !     1-add(2):jm(region)+add(2),&
    !     1:lmr)   :: fieldloc3d  !WP! watch out for am!
    real,dimension(:,:,:), intent(in)    :: fieldloc3d  !WP! watch out for am!

    ! local
    integer :: is,ie,js,je
    integer :: io,istat,sfstart,sfend,sfsnatt,root_id
    real,dimension(&
         1-add(1):im(region)+add(1),&
         1-add(2):jm(region)+add(2),&
         1-add(3):lm(region)+add(3)) :: fieldglob3d

    ! start

    call CheckShape( (/im(region)+2*add(1),jm(region)+2*add(2),lmr/), shape(fieldloc3d) )

    fieldglob3d=0.0
    if ( lmr == lmloc ) then
       call gather_tracer_k_3d( &
               fieldglob3d, im(region),jm(region),lm(region), add(1),add(2),add(3), &
               fieldloc3d, .false. )
       ! no need to broadcast result
       root_id=root_k
    else
       fieldglob3d=fieldloc3d
       root_id=root_t
    end if

    if ( myid == root_id ) then
       io=sfstart(filename, DFACC_CREATE)
       istat = sfsnatt(io,'im',    DFNT_INT32, 1, im(region))
       istat = sfsnatt(io,'jm',    DFNT_INT32, 1, jm(region))
       istat = sfsnatt(io,'lm',    DFNT_INT32, 1, lm(region))
       call io_write(io,im(region),'lon',jm(region),'lat', &
            lm(region),'hybrid',&
            fieldglob3d(&
            1:im(region),&
            1:jm(region),&
            1:lm(region)) ,name)
       print *,'dump_field3d: closing output file',sfend(io)
    end if

  end subroutine dump_field3d



  subroutine check_mass(region,text)

    use dims, only : lm
    use global_data, only : mass_dat
    use global_data, only : outdir

    ! --- in/out -----------------------------------

    integer,intent(in)            ::  region
    character(len=*),intent(in)   ::  text

    ! --- local ------------------------------------

    real,dimension(:,:,:),pointer     :: m
    real,dimension(:,:,:,:),pointer   :: rm
    integer :: lmr

    ! --- start ------------------------------------

    call barrier
    if ( previous_par(region) == 'tracer' ) then
       if ( myid == root_t ) print*,'check_mass: checking mass in ',text
       call barrier
       m => mass_dat(region)%m_t
       rm => mass_dat(region)%rm_t
       lmr=lm(region)

       call dump_field3d(region,m,lmr,(/2,2,0/),'m',trim(outdir)//'/mtrac.hdf')
       call dump_field4d(region,rm,lmr,ntracetloc,(/2,2,0/),'rm',trim(outdir)//'/rmtrac.hdf')
       call summasr( region, m, rm )
    else if ( previous_par(region) == 'levels' ) then
       if ( myid == root_k ) print*,'check_mass: checking mass in ',text
       call barrier
       m => mass_dat(region)%m_k
       rm => mass_dat(region)%rm_k
       lmr=lmloc

       call dump_field3d(region,m,lmr,(/2,2,0/),'m',trim(outdir)//'/mlev.hdf')
       call dump_field4d(region,rm,lmr,ntracet,(/2,2,0/),'rm',trim(outdir)//'/rmlev.hdf')
       call summasr( region, m, rm )
    end if
    call barrier
    nullify(m)
    nullify(rm)

  end subroutine check_mass


  ! ***


  subroutine summasr( region, m, rm )

    use dims, only : im, jm, lm
    use chem_param,   only : fscale

    ! --- in/out ---------------------------

    integer, intent(in)                ::  region
    real,dimension(:,:,:), pointer     ::  m
    real,dimension(:,:,:,:), pointer   ::  rm

    ! --- local ---------------------------

    integer i, j, k,n,root_id
    real summas,summin,summax
    real summas_all,summin_all,summax_all
    real,allocatable,dimension(:,:,:) :: mglob
    real,allocatable,dimension(:,:,:,:) :: rmglob

    ! --- begin ------------------------------

    allocate(mglob(-1:im(region)+2,-1:jm(region)+2,1:lm(region)))
    allocate(rmglob(-1:im(region)+2,-1:jm(region)+2,1:lm(region),ntracet))

    mglob=0.0
    rmglob=0.0

    if ( previous_par(region) == 'levels' ) then
       call gather_tracer_k_3d( &
               mglob, im(region),jm(region),lm(region), 2,2,0, &
               m, .false. ) ! no need to broadcast result
       call  gather_tracer_k( &
               rmglob, im(region),jm(region),lm(region), 2,2,0, ntracet, &
               rm, .false. ) ! no need to broadcast result
       root_id = root_k
    else
       mglob=m
       call gather_tracer_t( &
               rmglob, im(region),jm(region),lm(region), 2,2,0, ntracet, &
               rm, .false. ) ! no need to broadcast result
       root_id = root_t
    end if
    summin=1.e20
    summax=0.0
    if ( myid == root_id ) then
       summas=sum(mglob(1:im(region),1:jm(region),1:lm(region)) )
       do k=1,lm(region)
          do j=1,jm(region)
             do i=1,im(region)
                summin = min(summin,mglob(i,j,k))
                summax = max(summax,mglob(i,j,k))
             end do
          end do
       end do
       print *, ' Total mass            ',summas
       print *, ' Maximum value  m      ',summax
       print *, ' Minimum value  m      ',summin
       summin=1.0
       summax=0.0
       do k=1,lm(region)
          do j=1,jm(region)
             do i=1,im(region)
                summin = min(summin,fscale(1)*rmglob(i,j,k,1)/mglob(i,j,k))
                summax = max(summax,fscale(1)*rmglob(i,j,k,1)/mglob(i,j,k))
             end do
          end do
       end do
       print *, ' Maximum value mixing ratio tracer 1 ',summax
       print *, ' Minimum value mixing ratio tracer 1 ',summin
    end if
    deallocate(mglob)
    deallocate(rmglob)
    call barrier

  end subroutine summasr


end module mpi_comm
