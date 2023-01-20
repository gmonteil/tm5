!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!###############################################################################
#include "tm5.inc"

module zoom_tools

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- interface ---

  private

  public :: Zoom_Tools_Init, Zoom_Tools_Done
  public :: mix_edges, coarsen_region, update_parent
  public :: Fill_Interface_Children


  ! --- const ------------------------------

  character(len=*), parameter ::  mname = 'Zoom_Tools'


  ! --- local ------------------------------------

  integer    ::  itim_mix_edges


contains


  ! ====================================================================


  subroutine Zoom_Tools_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Zoom_Tools_Init'

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_mix_edges, 'mix edges', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Zoom_Tools_Init


  ! ***


  subroutine Zoom_Tools_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Zoom_Tools_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine Zoom_Tools_Done


  ! ***


  !
  ! distribute data at the x-edges of the zoom region
  ! over the whole interface cell
  ! written by mike botchev, march-june 1999
  ! adapted for the CRAY by Maarten Krol, 3-2000
  !

  subroutine mix_edges( region, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims       , only : im, jm, lm
    use dims       , only : xref, yref, zref
    use dims       , only : okdebug
    use dims       , only : touch_sp, touch_np
    use dims       , only : parent
    use global_data, only : mass_dat, region_dat
    use MeteoData  , only : m_dat
    use chem_param,  only : ntracet
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_const,   only : myid,lmloc,ntracetloc,which_par,previous_par
#endif

    implicit none

    ! input
    integer,intent(in)            ::  region
    integer, intent(out)             ::  status

    ! const
    character(len=*), parameter  ::  rname = mname//'/mix_edges'

    ! local
    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer    :: m
    integer,dimension(:,:), pointer    :: edge

    integer    :: imr, i, j, l, n, jmr, lmr, ip, jp, imp, jmp, nt
    integer    :: xref_, yref_, zref_, xyzref_, xyref_
    real       :: m_edge, rm_edge, rxm_edge, rym_edge, rzm_edge
    logical    :: xedges , ynp,ysp

    ! relaxation factor for slopes; relax=0.0 means upwinding
    real,parameter:: relax=1.0
    real,dimension(:,:),allocatable   :: to_parent
    real,dimension(:,:),allocatable   :: to_parentx
    real,dimension(:,:),allocatable   :: to_parenty
    real,dimension(:,:),allocatable   :: to_parentz

    ! start

#ifdef MPI
    which_par=previous_par(region)
    if ( which_par == 'tracer' .and. ntracetloc == 0 ) return
    if ( which_par == 'level' .and. lmloc == 0 ) return  !WP!
#endif

    call GO_Timer_Start( itim_mix_edges, status )
    IF_NOTOK_RETURN(status=1)

    xedges = .true.
    ysp = .true.  !mix also NP/SP?
    ynp = .true.

    imr = im(region)
    jmr = jm(region)

#ifdef MPI
    if ( which_par == 'tracer') then
#endif

       m => m_dat(region)%data
       rm => mass_dat(region)%rm_t
       rxm => mass_dat(region)%rxm_t
       rym => mass_dat(region)%rym_t
       rzm => mass_dat(region)%rzm_t

       lmr = lm(region)
#ifdef MPI
       nt=ntracetloc
#else
       nt=ntracet
#endif

#ifdef MPI
    else if ( which_par == 'levels' ) then

       m => mass_dat(region)%m_k
       rm => mass_dat(region)%rm_k
       rxm => mass_dat(region)%rxm_k
       rym => mass_dat(region)%rym_k
       rzm => mass_dat(region)%rzm_k

       lmr = lmloc
       nt=ntracet
    end if
#endif

    edge => region_dat(region)%edge

    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    if ( zref_ /= 1 ) call escape_tm( &
         'mix_edges: zooming in z direction not supported!')

    xyref_ = xref_*yref_

    if ( okdebug ) then
       print *,'mix_edges: region=',region
    end if
    if ( (xref_ == 1) .or. (im(region)/xref(region) == im(1)) ) then
       if ( okdebug ) print *,     &
            'mix_edges: no refinement or periodic bc, skipping x-walls'
       xedges = .false.
       if ( region == 1 ) then
          if ( okdebug ) print *, 'mix_edges: region = 1, returning'
          call GO_Timer_End( itim_mix_edges, status )
          IF_NOTOK_RETURN(status=1)
          return
       end if
    end if

    if ( touch_sp(region) == 1) then
       ysp = .false.
       if ( okdebug ) print *,     &
            'mix_edges: SP touching...skip SP y-walls'
    end if
    if(touch_np(region) == 1) then
       ynp = .false.
       if ( okdebug ) print *,     &
            'mix_edges: NP touching...skip NP y-walls'
    end if
    imp = imr/xref_   !in 'coarse' resolution
    jmp = jmr/yref_

    !$omp parallel default ( shared ) &
    !$omp   private ( i, j, l, ip, jp, to_parent )
    allocate(to_parent(imp,jmp))
    !$omp do schedule(guided)
    do l=1,lmr !WP! depends on tracers or levels
       to_parent = 0.0
       do j=1,jmr
          jp = 1 + (j-1)/yref_
          do i=1,imr
             if ( edge(i,j) /= -1 ) cycle   !do only for edge of child...
             ip = 1 + (i-1)/xref_
             to_parent(ip,jp) = to_parent(ip,jp) + m(i,j,l)/xyref_
          end do
       end do
       if ( xedges ) then     !mix the x-wall of the zoom region.
          do j = 1,jmr
             jp = 1 + (j-1)/yref_
             do i = 1, xref_
                m(i,j,l) = to_parent(1,jp)
             end do
             do i = imr-xref_+1,imr
                m(i,j,l) = to_parent(imp,jp)
             end do
          end do
       end if    !xedges
       if ( ysp ) then
          do j = 1,yref_
             do i = 1, imr
                ip = 1 + (i-1)/xref_
                m(i,j,l) = to_parent(ip,1)
             end do
          end do
       end if
       if ( ynp ) then
          do j = jmr-yref_+1,jmr
             do i = 1, imr
                ip = 1 + (i-1)/xref_
                m(i,j,l) = to_parent(ip,jmp)
             end do
          end do
       end if
    end do  !l
    !$omp   end do
    deallocate( to_parent )
    !$omp end parallel

    do n=1,nt  !WP! depends on tracers or levels
      !$omp parallel default ( shared ) &
      !$omp private ( i, j, l, ip, jp, to_parent, to_parentx, to_parenty, to_parentz)
      allocate(to_parent(imp,jmp))
      allocate(to_parentx(imp,jmp))
      allocate(to_parenty(imp,jmp))
      allocate(to_parentz(imp,jmp))
      !$omp do schedule(guided)
      do l=1,lmr !WP! depends on tracers or levels
        to_parent = 0.0
        to_parentx = 0.0
        to_parenty = 0.0
        to_parentz = 0.0
        do j=1,jmr
           jp = 1 + (j-1)/yref_
           do i=1,imr
              if ( edge(i,j) /= -1 ) cycle   !do only for edge of child...
              ip = 1 + (i-1)/xref_
              to_parent(ip,jp) = to_parent(ip,jp) + rm(i,j,l,n)/xyref_
              to_parentx(ip,jp) = to_parentx(ip,jp) + rxm(i,j,l,n)/xyref_
              to_parenty(ip,jp) = to_parenty(ip,jp) + rym(i,j,l,n)/xyref_
              to_parentz(ip,jp) = to_parentz(ip,jp) + rzm(i,j,l,n)/xyref_
           end do
        end do
        if ( xedges ) then     !mix the x-wall of the zoom region.
           do j = 1,jmr
              jp = 1 + (j-1)/yref_
              do i = 1, xref_
                 rm(i,j,l,n) = to_parent(1,jp)
                 rxm(i,j,l,n) = to_parentx(1,jp)
                 rym(i,j,l,n) = to_parenty(1,jp)
                 rzm(i,j,l,n) = to_parentz(1,jp)
              end do
              do i = imr-xref_+1,imr
                 rm(i,j,l,n) = to_parent(imp,jp)
                 rxm(i,j,l,n) = to_parentx(imp,jp)
                 rym(i,j,l,n) = to_parenty(imp,jp)
                 rzm(i,j,l,n) = to_parentz(imp,jp)
              end do
           end do
        end if !xedges
        if ( ysp ) then
           do j = 1,yref_
              do i = 1, imr
                 ip = 1 + (i-1)/xref_
                 rm(i,j,l,n) = to_parent(ip,1)
                 rxm(i,j,l,n) = to_parentx(ip,1)
                 rym(i,j,l,n) = to_parenty(ip,1)
                 rzm(i,j,l,n) = to_parentz(ip,1)
              end do
           end do
        end if
        if ( ynp ) then
           do j = jmr-yref_+1,jmr
              do i = 1, imr
                 ip = 1 + (i-1)/xref_
                 rm(i,j,l,n) = to_parent(ip,jmp)
                 rxm(i,j,l,n) = to_parentx(ip,jmp)
                 rym(i,j,l,n) = to_parenty(ip,jmp)
                 rzm(i,j,l,n) = to_parentz(ip,jmp)
              end do
           end do
        end if
      end do   !l
      !$omp end do
      deallocate(to_parent)
      deallocate(to_parentx)
      deallocate(to_parenty)
      deallocate(to_parentz)
      !$omp end parallel
    end do   !nt


    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(edge)

    call GO_Timer_End( itim_mix_edges, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine mix_edges



  recursive subroutine coarsen_region(region)
    !
    ! calculates values in the parent identical to the child.
    ! Avoids small numerical differences that are caused by the
    ! limited accuracy of the stored pressure and wind-data.
    ! written by Maarten Krol, december 1999
    !
    use dims
    use global_data,   only : mass_dat, wind_dat
    use MeteoData    , only : m_dat
    use advect_tools,  only : m2phlb, m2phlb1
    use toolbox,       only : escape_tm
#ifdef MPI
    use mpi_comm, only : check_domain,barrier,barrier_t,stopmpi
    use mpi_const,only : ntracetloc,myid,root_t,lmar,root_k,lmloc,ntracet_ar
#endif

    implicit none

    ! in/out:
    integer,intent(in)   :: region

    ! local:
    real,dimension(:,:,:),pointer     :: am,bm,cm,m,amp,bmp,cmp,mp
    real,dimension(:,:,:),allocatable :: top
    integer i,j,l,ip,jp,lp,regiop,imp,jmp,lmp,imr,jmr,lmr
    integer xref_,yref_,zref_,tref_
    real,dimension(0:im(region)+1,0:jm(region)+1,0:lm(region) ) :: fieldglob

    ! start

    regiop = parent(region)
#ifdef MPI
    call check_domain(region,'p','tracer') !WP! check whether region and
#endif
    !parents are on tracer domain
    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)
    if ( regiop > 0 ) then
#ifdef MPI
       if ( ntracetloc /= 0 ) then !only PE's with nonzero ntracetloc continue
#endif
#ifdef MPI
          if ( okdebug .and. myid == root_t ) &
               print *,'coarsen_region: region:',region
#else
          if ( okdebug ) print *,'coarsen_region: region:',region
#endif
          xref_ = xref(region)/xref(regiop)
          yref_ = yref(region)/yref(regiop)
          zref_ = zref(region)/zref(regiop)
          tref_ = tref(region)/tref(regiop)
          imp = imr/xref_
          jmp = jmr/yref_
          lmp = lmr/zref_

          if ( ibeg(region) < iend(region) .and. &
               imp /= iend(region)-ibeg(region)+1 )  &
               call escape_tm('coarsen_region: stop')
          if ( jmp /= jend(region)-jbeg(region)+1 ) &
               call escape_tm('coarsen_region: stop')
          if ( lmp /= lend(region)-lbeg(region)+1 ) &
               call escape_tm('coarsen_region: stop')

          ! info ...
          !print *,'coarsen_region: coarsening parent',regiop, 'with child',region

          allocate(top(0:imp,jmp,lmp))
          amp => wind_dat(regiop)%am_t
          am => wind_dat(region)%am_t
          top = 0.0
          do l=1,lmr
             lp = 1 + (l-1)/zref_
             do j=1,jmr
                jp = 1 + (j-1)/yref_
                do ip=0,imp
                   top(ip,jp,lp) = top(ip,jp,lp) + am(ip*xref_,j,l)*tref_
                end do
             end do
          end do
          !   copy top to parent
          do l=1,lmp
             lp = lbeg(region)-1+l
             do j=1,jmp
                jp = jbeg(region)-1+j
                do i=0,imp
                   ip = mod(ibeg(region)-2+i,im(regiop))+1
                   amp(ip,jp,lp) = top(i,j,l)
                end do
             end do
          end do

          deallocate(top)
          nullify(am)
          nullify(amp)

          bm => wind_dat(region)%bm_t
          bmp => wind_dat(regiop)%bm_t
          allocate(top(imp,jmp+1,lmp))
          top = 0.0
          do l=1,lmr
             lp = 1 + (l-1)/zref_
             do i=1,imr
                ip = 1 + (i-1)/xref_
                do jp = 1,jmp+1
                   top(ip,jp,lp) = top(ip,jp,lp) + bm(i,1 + &
                        (jp-1)*yref_,l)*tref_
                end do
             end do
          end do
          do l=1,lmp
             lp = lbeg(region)-1+l
             do j=1,jmp+1
                jp = jbeg(region)-1+j
                do i=1,imp
                   ip = mod(ibeg(region)-2+i,im(regiop))+1
                   bmp(ip,jp,lp) = top(i,j,l)
                end do
             end do
          end do

          deallocate(top)
          nullify(bm)
          nullify(bmp)

          allocate(top(imp,jmp,lmp))
          cm => wind_dat(region)%cm_t
          cmp => wind_dat(regiop)%cm_t
          top = 0.0
          do i=1,imr
             ip = 1 + (i-1)/xref_
             do j = 1,jmr
                jp = 1 + (j-1)/yref_
                do lp=1,lmp
                   top(ip,jp,lp) = top(ip,jp,lp) + cm(i,j,lp*zref_)*tref_
                end do
             end do
          end do
          do l=1,lmp
             lp = lbeg(region)-1+l
             do j=1,jmp
                jp = jbeg(region)-1+j
                do i=1,imp
                   ip = mod(ibeg(region)-2+i,im(regiop))+1
                   cmp(ip,jp,lp) = top(i,j,l)
                end do
             end do
          end do

          nullify(cm)
          nullify(cmp)

          top = 0.0

          m  => m_dat(region)%data
          mp => m_dat(regiop)%data

          do l=1,lmr
             lp = 1 + (l-1)/zref_
             do j = 1,jmr
                jp = 1 + (j-1)/yref_
                do i=1,imr
                   ip = 1 + (i-1)/xref_
                   top(ip,jp,lp) = top(ip,jp,lp) + m(i,j,l)
                end do
             end do
          end do

          do l=1,lmp
             lp = lbeg(region)-1+l
             do j=1,jmp
                jp = jbeg(region)-1+j
                do i=1,imp
                   ip = mod(ibeg(region)-2+i,im(regiop))+1
                   mp(ip,jp,lp) = top(i,j,l)
                end do
             end do
          end do

          deallocate(top)
          nullify(m)
          nullify(mp)

#ifdef MPI
       end if ! all PE's from here
#endif

       if ( regiop /= 1 ) then
!???          call m2phlb(regiop,parent(regiop))
          call m2phlb(regiop)
       else
          call m2phlb1(regiop)
       end if
       ! call grand-parent recursively...
       call coarsen_region(regiop)
#ifdef MPI
       call barrier
#endif

    end if !regiop>0


  end subroutine coarsen_region




  subroutine update_parent( region )
    !
    ! update rm of its parent in correspondence with itself
    ! written by mike botchev, march-june 1999
    !
    use dims
    use global_data  ,only : mass_dat
    use MeteoData   , only : m_dat
    use budget_global,only : sum_update
    use chem_param,   only : ntracet
    use toolbox,      only : escape_tm
#ifdef MPI
    use mpi_const
    use mpi_comm,     only : check_domain
#endif
    use misctools,    only : collapse_loop

    implicit none

    ! input
    integer,intent(in)   :: region

    ! local
    real,dimension(:,:,:,:),pointer         :: rm,rxm,rym,rzm,rmp,rxmp,rymp,rzmp
    real,dimension(:,:,:),pointer           :: m,mp
    real,dimension(:,:,:),allocatable       :: to_parent

    real,dimension(xref(region)/xref(parent(region)), &
        yref(region)/yref(parent(region)), &
        zref(region)/zref(parent(region)))  :: slope
    real                                    :: dummy
    real                                    :: sum_parent(ntracet), sum_to_parent(ntracet), sum_parent_all(ntracet), sum_to_parent_all(ntracet)
    integer                                 :: i, iend_loop, iip, ip, j, jp, l, lp, n
    integer                                 :: lb, le, ib, ie, jb, je
    integer                                 :: my_parent, xref_, yref_, zref_
    integer                                 :: imp, jmp, lmp, imr, jmr, lmr, nt
    real                                    :: m_value, mxval, sum_mass
    integer,dimension(3)                    :: mxloc
    integer                                 :: communicator, root_id
    ! to speed up nested loops
    integer                                 :: i_iter_3d, N_iter_3d
    integer, allocatable                    :: counter(:,:)

    ! start

    ! global region ? then not necessary to update a parent:
    if (region == 1) return

#ifdef MPI
    which_par=previous_par(region)

    !WP! make sure parents are on same domain
    call check_domain(region,'p',which_par)

    if ( which_par == 'tracer' .and. ntracetloc == 0 ) return
    if ( which_par == 'level' .and. lmloc == 0 ) return  !WP!

#endif

    imr = im(region)
    jmr = jm(region)
    my_parent = parent(region)
    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    zref_ = zref(region)/zref(my_parent)

#ifdef MPI
    if ( which_par == 'tracer' ) then
#endif

       m => m_dat(region)%data
       rm => mass_dat(region)%rm_t
       rxm => mass_dat(region)%rxm_t
       rym => mass_dat(region)%rym_t
       rzm => mass_dat(region)%rzm_t

       mp => m_dat(my_parent)%data
       rmp => mass_dat(my_parent)%rm_t
       rxmp => mass_dat(my_parent)%rxm_t
       rymp => mass_dat(my_parent)%rym_t
       rzmp => mass_dat(my_parent)%rzm_t

       lmr = lm(region)
#ifdef MPI
       nt=ntracetloc
       communicator=com_trac  !WP! assign com_trac as communicator
       root_id=root_t
#else
       nt=ntracet
#endif
#ifdef MPI
    else if ( which_par == 'levels' ) then

       m => mass_dat(region)%data_k
       rm => mass_dat(region)%rm_k
       rxm => mass_dat(region)%rxm_k
       rym => mass_dat(region)%rym_k
       rzm => mass_dat(region)%rzm_k

       mp => m_dat(my_parent)%data_k
       rmp => mass_dat(my_parent)%rm_k
       rxmp => mass_dat(my_parent)%rxm_k
       rymp => mass_dat(my_parent)%rym_k
       rzmp => mass_dat(my_parent)%rzm_k

       lmr = lmloc
       nt=ntracet
       communicator=com_lev  !WP! assign com_lev as communicator
       root_id=root_k
    end if
#endif

    if ( okdebug ) then
       print *,'update_parent: my_parent=',my_parent, &
            ' x-,y-,zref_: ',xref_,yref_,zref_
    end if

    imp = im(region)/xref_
    jmp = jm(region)/yref_
    !lmp = lm(region)/zref_
    lmp = lmr     ! CMK changed: when paralel over levels lmp can differ from lmr

    if ( ibeg(region) < iend(region) .and. imp /= iend(region)-ibeg(region)+1 )  &
         call escape_tm('update_parent: stop')
    if ( jmp /= jend(region)-jbeg(region)+1 ) &
         call escape_tm('update_parent: stop')
    !CMK if( lmp /= lend(region)-lbeg(region)+1 ) &
    !CMK    call escape_tm('update_parent: stop')

    ! storage:
    allocate( to_parent(1:imp,1:jmp,lmp) )

    to_parent = 0.0
    ! Sourish Basu: Speed up triple-nested loop by collapsing it
    N_iter_3d = imr * jmr * lmr
    allocate(counter(N_iter_3d, 3))
    call collapse_loop(1, imr, 1, jmr, 1, lmr, counter)
    !$omp parallel do schedule(guided) private(i_iter_3d, i, j, l, ip, jp, lp) reduction(+:to_parent)
    do i_iter_3d = 1, N_iter_3d
        i = counter(i_iter_3d, 1)
        j = counter(i_iter_3d, 2)
        l = counter(i_iter_3d, 3)
        ip = 1 + (i-1)/xref_
        jp = 1 + (j-1)/yref_
        lp = 1 + (l-1)/zref_
        to_parent(ip,jp,lp) = to_parent(ip,jp,lp) + m(i,j,l)
    end do
    !$omp end parallel do
    deallocate(counter)

    !do i = 1, imr
       !ip = 1 + (i-1)/xref_
       !do j = 1, jmr
          !jp = 1 + (j-1)/yref_
          !do l=1,lmr
             !lp = 1 + (l-1)/zref_
             !to_parent(ip,jp,lp) = to_parent(ip,jp,lp) + m(i,j,l)
          !end do
       !end do
    !end do
    if ( ibeg(region) < iend(region) ) then   !no dateline crossing!
       mp(ibeg(region):iend(region),jbeg(region):jend(region),1:lmr) = &
            to_parent
    else
       mp(ibeg(region):im(region),jbeg(region):jend(region),1:lmr)  = &
            to_parent(1:im(region)-ibeg(region)+1,:,:)
       mp(1:iend(region),jbeg(region):jend(region),1:lmr)  = &
            to_parent(im(region)-ibeg(region)+2:,:,:)
    end if

    ! Sourish Basu: Modifying the following block to extend sum_update to multiple tracers
    sum_parent = 0.0
    sum_to_parent = 0.0
    sum_parent_all = 0.0
    sum_to_parent_all = 0.0
    do n=1,nt
       to_parent = 0.0
       do i = 1, imr
          ip = 1 + (i-1)/xref_
          do j = 1, jmr
             jp = 1 + (j-1)/yref_
             do l=1,lmr
                lp = 1 + (l-1)/zref_
                to_parent(ip,jp,lp) = to_parent(ip,jp,lp) + rm(i,j,l,n)
                !if ( n == 1 ) sum_to_parent = sum_to_parent + rm(i,j,l,1) ! Sourish Basu
                sum_to_parent(n) = sum_to_parent(n) + rm(i,j,l,n)
#ifdef MPI
                if ( which_par == 'tracer' .and. myid/=pe_first_tracer) &
                     sum_to_parent=0.0 !WP! only on one when in tracer
#endif
             end do
          end do
       end do
       if ( ibeg(region) < iend(region) ) then   !no dateline crossing!
          !if ( n == 1 ) then ! Sourish Basu
             sum_parent(n) = sum(rmp(ibeg(region):iend(region), &
                  jbeg(region):jend(region), 1:lmr,n) )
#ifdef MPI
             if ( which_par == 'tracer' .and. myid /= pe_first_tracer ) &
                  sum_parent = 0.0 !WP! only on one when in tracer
#endif
          !end if
          ! CMK this line gave a problem 1:lmr /= to_parent dimension!
          rmp(ibeg(region):iend(region),jbeg(region):jend(region),1:lmr,n) = &
               to_parent
       else
          !if (n == 1) then ! Sourish Basu
             sum_parent(n) = sum(rmp(ibeg(region):im(region), &
                  jbeg(region):jend(region),   &
                  1:lmr,n) ) + &
                  sum(rmp(1:iend(region), &
                  jbeg(region):jend(region),   &
                  1:lmr,n) )
#ifdef MPI
             if ( which_par == 'tracer' .and. myid /= pe_first_tracer ) &
                  sum_parent=0.0 !WP! only on one when in tracer
#endif
          !end if
          rmp(ibeg(region):im(region),jbeg(region):jend(region),1:lmr,n)  = &
               to_parent(1:im(region)-ibeg(region)+1,:,:)
          rmp(1:iend(region),jbeg(region):jend(region),1:lmr,n)  = &
               to_parent(im(region)-ibeg(region)+2:,:,:)
       end if
    end do   ! ntracet

#ifdef MPI
    call mpi_allreduce(sum_parent,sum_parent_all,1,my_real, &
         mpi_sum,communicator,ierr)
    call mpi_allreduce(sum_to_parent,sum_to_parent_all,1,my_real, &
         mpi_sum,communicator,ierr)
#else
    sum_parent_all = sum_parent
    sum_to_parent_all = sum_to_parent
#endif

    sum_update(my_parent,:) = sum_update(my_parent,:) &
         + sum_to_parent_all - sum_parent_all
    ! End modification by Sourish Basu

    deallocate(to_parent)

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

    nullify(mp)
    nullify(rmp)
    nullify(rxmp)
    nullify(rymp)
    nullify(rzmp)

    ! mk still to do : limit the slopes after rm update!
    !if ( okdebug ) then
    !   !call write_rm(rm-rm1,m,'r1.dat')
    !end if

  end subroutine update_parent


  ! *


  subroutine Fill_Interface_Children( region, status )

    use Dims, only : children

    ! --- in/out ----------------------------------------

    integer, intent(in)     ::  region
    integer, intent(out)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Fill_Interface_Children'

    ! --- local ----------------------------------

    integer     ::  nchild, ichild
    integer     ::  cregion

    ! --- begin ----------------------------------

    ! number of children for this region:
    nchild = children(region,0)
    ! loop over children:
    do ichild = 1, nchild
      ! child region:
      cregion = children(region,ichild)
      ! fill:
      call Fill_Interface( cregion, region, status )
      IF_NOTOK_RETURN(status=1)
    end do  ! children

    ! ok
    status = 0

  end subroutine Fill_Interface_Children


  ! *


  subroutine Fill_Interface( region, parent, status )

    !
    ! Fill 'rm' in interface cells of child regions with fractions
    ! of 'rm' values in current region.
    ! Fractions are relative to the air mass in the child cells
    ! compared to the parent cell.
    !

    use Dims        , only : im, jm, lm
    use Dims        , only : isr, ier, jsr, jer     ! range of region cells in core
    use Dims        , only : ibeg, iend, jbeg, jend ! range of parent cells
    use Dims        , only : xref, yref
    use chem_param  , only : ntracet
    use global_data , only : mass_dat
    use MeteoData   , only : m_dat
    use ParTools    , only : which_par

    ! --- in/out ----------------------------------------

    integer, intent(in)     ::  region
    integer, intent(in)     ::  parent
    integer, intent(out)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Fill_Interface'

    ! --- local ----------------------------------

    integer           ::  i, j, l, n
    integer           ::  imr, jmr, lmr
    integer           ::  xref_, yref_
    integer           ::  ip, jp
    real, pointer     ::  rm (:,:,:,:)
    real, pointer     ::  rmp(:,:,:,:)
    real, pointer     ::  m(:,:,:)
    real, pointer     ::  mp(:,:,:)
    real              ::  frac

    ! --- begin ----------------------------------

    ! local grid size:
    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    ! point to mass arrays:
    if ( which_par == 'tracer' ) then
      ! region data:
      m  => m_dat(region)%data
      rm => mass_dat(region)%rm_t
      ! idem for parent:
      mp  => m_dat(parent)%data
      rmp => mass_dat(parent)%rm_t
    else
      write (gol,'("unsupported which_par : ",a)') trim(which_par); call goErr
      TRACEBACK; status=1; return
    end if

    ! refinement factor:
    xref_ = xref(region) / xref(parent)
    yref_ = yref(region) / yref(parent)

    ! loop over y direction:
    do j = 1, jmr
      ! loop over x direction:
      do i = 1, imr
        ! core region ? then skip, only interface:
        if ( (i >= isr(region)) .and. (i <= ier(region)) .and. &
             (j >= jsr(region)) .and. (j <= jer(region))         ) cycle
        ! parent cell:
        ip = ibeg(region) + floor((i-0.5)/xref_)
        jp = jbeg(region) + floor((j-0.5)/yref_)
        ! loop over levels:
        do l = 1, lmr
          ! mass fraction in region cell compared to parent:
          frac = m(i,j,l) / mp(ip,jp,l)
          ! loop over transported tracers:
          do n = 1, ntracet
            ! fill fraction:
            rm(i,j,l,n) = rmp(ip,jp,l,n) * frac
            !! testing ...
            !if (l==1) print *, i, j, ip, jp,  rmp(ip,jp,l,n), frac, rm(i,j,l,n)
          end do ! tracers
        end do ! levels
      end do  ! i
    end do  ! j

    ! ok
    status = 0

  end subroutine Fill_Interface


end module zoom_tools
