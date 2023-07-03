!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################
#include "tm5.inc"

module advect_tools

  implicit none

  ! --- in/out ------------------------------

  private

  public :: advect_m, m2phlb1, m2phlb, coarsen_ambm, mass_to_pressure
  public :: dynam0

#ifdef MPI
  public :: calc_phlb_k
#endif

  character(len=*), parameter   :: mname = 'advect_tools'

contains


  subroutine m2phlb1( region, check_pressure )
    !-----------------------------------------
    ! routine to convert m(im,jm,lm) into pressure layers.
    ! First, calculate the surface pressure in p(im,jm), THEN
    ! using the hybrid coordinate system, calculate phlb(im,jm,lmp1)
    ! special routine for region 1
    !   mk december 1998
    !-----------------------------------------
    use binas,       only : grav
    use dims
    use global_data, only : mass_dat, region_dat
    use MeteoData  , only : sp_dat, phlb_dat, m_dat
#ifdef MPI
    use mpi_const
    use mpi_comm,    only : barrier, stopmpi, check_mass
#endif
    use misctools,   only : collapse_loop

    implicit none

    ! input/output
    integer,intent(in)              ::  region
    logical, intent(in), optional   ::  check_pressure

    ! local
    logical                         :: do_check_pressure
    real,dimension(:,:,:), pointer  :: m, phlb
    real,dimension(:,:,:), pointer  :: p
    real,dimension(:),     pointer  :: dxyp
    integer                         :: i,j,l,imax,xref_,yref_,imr,jmr
    real                            :: maxdiff,m_edge,maxval_one,maxdiff_all,maxval_all
    real                            :: mnew
    real, allocatable               :: pold(:,:), phelp(:,:)
    !real,dimension(-1:im(region)+2,-1:jm(region)+2) :: pold,phelp
    integer                                         :: communicator,ll,root_id,lmr
    ! variables to unroll a loop
    integer                         :: N_iter_3d, i_iter_3d
    integer, allocatable            :: counter(:,:)

    ! start

    ! check flag:
    do_check_pressure = .true.
    if ( present(check_pressure) ) do_check_pressure = check_pressure

#ifdef MPI
    which_par=previous_par(region)
    if ( which_par == 'tracer' .and. ntracetloc == 0 ) return
    if ( which_par == 'levels' .and. lmloc == 0 ) return  !WP!
#endif

#ifdef MPI
    if ( which_par == 'tracer' ) then
#endif
       p    =>   sp_dat(region)%data
       m    =>    m_dat(region)%data
       phlb => phlb_dat(region)%data
       dxyp => region_dat(region)%dxyp

       lmr = lm(region)
#ifdef MPI
       communicator=com_trac  !WP! assign com_trac as communicator
       root_id=root_t
    else if(which_par=='levels') then
       p    =>   sp_dat(region)%data
       m    =>    m_dat(region)%data_k
       phlb => phlb_dat(region)%data_k
       dxyp => region_dat(region)%dxyp

       lmr = lmloc
       communicator=com_lev  !WP! assign com_lev as communicator
       root_id=root_k
    end if
#endif

    allocate(pold(-1:im(region)+2,-1:jm(region)+2))
    allocate(phelp(-1:im(region)+2,-1:jm(region)+2))

    imr=im(region)
    jmr=jm(region)

    pold = p(:,:,1) ! store the old surface pressure
    phelp = zero
    p = zero

    !$omp parallel do schedule(guided) private(j)
    do j = 1, jmr
        phelp(1:imr, j) = sum(m(1:imr, j, 1:lmr), dim=2) * grav/dxyp(j) ! phelp = surface pressure according to air mass
    end do
    !$omp end parallel do

    !do l=1,lmr
       !do j=1,jmr
          !do i=1,imr
             !phelp(i,j) = phelp(i,j) + m(i,j,l)*grav/dxyp(j)
          !end do
       !end do
    !end do

#ifdef MPI
    if ( which_par == 'tracer' ) then
#endif
       p(:,:,1) = phelp ! now p has surface pressure according to air mass
#ifdef MPI
    else if ( which_par == 'levels' ) then
       call mpi_allreduce(phelp,p,(imr+4)*(jmr+4),my_real,mpi_sum,com_lev,ierr)
    end if
#endif

    !$omp parallel do schedule(guided) private(l, ll)
    do  l=1,lmr+1
       ll=l
#ifdef MPI
       ! offset for levels
       if ( myid > 0 .and. which_par == 'levels' ) ll=sum(lmar(0:myid-1))+l
#endif
       phlb(1:imr,1:jmr,l) = at(ll) + bt(ll) * p(1:imr,1:jmr,1)
       !do j=1,jmr
          !do i=1,imr
             !phlb(i,j,l) = at(ll)+bt(ll)*p(i,j,1)
          !end do
       !end do
    end do
    !$omp end parallel do

    ! check pressures ?
    if ( do_check_pressure ) then

    maxdiff = 0.0

    N_iter_3d = imr*jmr*lmr
    allocate(counter(N_iter_3d, 3))
    call collapse_loop(1, lmr, 1, jmr, 1, imr, counter)
    !$omp parallel do schedule(guided) private(i_iter_3d, i, j, l, mnew) reduction(max:maxdiff)
    do i_iter_3d = 1, N_iter_3d
        l = counter(i_iter_3d, 1)
        j = counter(i_iter_3d, 2)
        i = counter(i_iter_3d, 3)
        mnew=(phlb(i,j,l)-phlb(i,j,l+1))*dxyp(j)/grav
        maxdiff = max(maxdiff,abs(mnew-m(i,j,l))/mnew)
    end do
    !$omp end parallel do
    deallocate(counter)

    !do l = 1,lmr
       !do j=1,jmr
          !do  i=1,imr
             !mnew=(phlb(i,j,l)-phlb(i,j,l+1))*dxyp(j)/grav
             !maxdiff = max(maxdiff,abs(mnew-m(i,j,l))/mnew)
          !end do
       !end do
    !end do

      maxval_one=maxval(abs(p(1:im(region),1:jm(region),1)-  &
         pold(1:im(region),1:jm(region)))/ &
         pold(1:im(region),1:jm(region)) )
#ifdef MPI
    call mpi_allreduce(maxdiff,maxdiff_all,1,my_real,mpi_max,communicator,ierr)
    call mpi_allreduce(maxval_one,maxval_all,1,my_real,mpi_max,communicator,ierr)
    if ( okdebug .and. myid == root_id ) &
         print*,'m2phlb1: maxdiff (mnew-m)/mnew on all pes:',maxdiff_all
    if ( okdebug .and. myid == root_id ) &
         print*,'m2phlb1: max pressure difference (%) on all pes :',maxval_all
#else
    if ( okdebug ) print*,'m2phlb1: maxdiff (mnew-m)/mnew :',maxdiff
    if ( okdebug ) print*,'m2phlb1: max pressure difference (%) :',maxval_one
#endif
    end if  ! check pressures

    nullify(p)
    nullify(m)
    nullify(phlb)
    nullify(dxyp)

    deallocate(pold, phelp)

  end subroutine m2phlb1

  subroutine mass_to_pressure(region)
    ! Sourish Basu: Routine to convert mass to pressure levels after advection. One routine for all regions.
    use binas,          only : grav
    use global_data,    only : region_dat
    use MeteoData,      only : sp_dat, phlb_dat, m_dat
    use dims,           only : im, jm, lm, xref, yref, parent, xcyc, touch_np, touch_sp, jbeg, revert, newsrun, at, bt

    implicit none

    integer, intent(in)     :: region

    real, pointer           :: m(:,:,:), phlb(:,:,:), p(:,:,:)
    real, pointer           :: dxyp(:), dxy(:)
    integer                 :: i, j, l, xref_, yref_, imr, jmr, lmr, jpar

    character(len=*), parameter :: rname = mname//'/mass_to_pressure'

    p    =>   sp_dat(region)%data
    m    =>    m_dat(region)%data
    phlb => phlb_dat(region)%data
    dxyp => region_dat(region)%dxyp

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    p = 0.0
    ! construct the surface pressure from air mass
    !$omp parallel do schedule(guided) private(j)
    do j = 1, jmr
        p(1:imr, j, 1) = sum(m(1:imr, j, 1:lmr), dim=2) * grav/dxyp(j)
    end do
    !$omp end parallel do

    ! for nested regions, smooth the pressures at the interface cells
    if (region > 1) then
        dxy  => region_dat(parent(region))%dxyp
        yref_ = yref(region)/yref(parent(region)) ! will give the number of latitudinal interface cells in the region
        xref_ = xref(region)/xref(parent(region)) ! will give the number of longitudinal interface cells in the region
        ! MK: one exception: in adjoint pressure is restored and is already corrected for this if newsrun=.TRUE.
        if (.not. (revert == -1 .and. newsrun)) then
            do j = 1, jmr
                do i = 1, imr
                    if ((i > xref_ .and. i < imr-xref_+1) .and. (j > yref_ .and. j < jmr-yref_+1)) cycle    ! not an interface
                    if ((xcyc(region) == 1) .and. (j > yref_ .and. j < jmr-yref_+1)) cycle                  ! periodic boundary
                    if ((touch_np(region) == 1) .and. (j > jmr-yref_)) cycle                                ! np_touching
                    if ((touch_sp(region) == 1) .and. (j < yref_+1)) cycle                                  ! sp_touching
                    jpar = jbeg(region) + (j-1)/yref_   ! to get grid cell area in the parent region
                    p(i,j,1) = p(i,j,1)*dxyp(j)*xref_*yref_/dxy(jpar)
                end do
            end do
        end if
        nullify(dxy)
    end if

    !$omp parallel do schedule(guided) private(l)
    do l = 1, lmr+1
        !phlb(1:imr, 1:jmr, l) = at(l) + bt(l) * p(1:imr, 1:jmr, 1)
        phlb(:, :, l) = at(l) + bt(l) * p(:, :, 1)
    end do
    !$omp end parallel do

    nullify(p, m, phlb, dxyp)

    ! Debug
    !write(*,'(a, " :: Maximum difference between phlb_dat and sp_dat in region ", i1, " = ", es20.12)') &
        !rname, region, maxval(abs(phlb_dat(region)%data(:,:,1) - sp_dat(region)%data(:,:,1)))
    ! End debug

  end subroutine mass_to_pressure

  subroutine m2phlb(region)
    !-----------------------------------------
    !routine to convert m(im,jm,lm) into pressure layers.
    !First, calculate the surface pressure in p(im,jm), then
    !using the hybrid coordinate system, calculate phlb(im,jm,lmp1)
    ! mk december 1998
    ! special routine for region > 1
    !-----------------------------------------
    use binas,       only : grav
    use dims
!    use io_hdf
    use global_data, only : mass_dat, region_dat
    use MeteoData  , only : sp_dat, phlb_dat, m_dat
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_const
    use mpi_comm , only : barrier_t, barrier_k
#endif

    implicit none

    ! input
    integer,intent(in)   ::region

    ! local
    real,dimension(:,:,:), pointer    :: m, phlb, mpar
    real,dimension(:,:,:), pointer    :: p
    real,dimension(:),     pointer    :: dxyp, dxy
    !real,dimension(-1:im(region)+2,-1:jm(region)+2)  :: pold
    !real,dimension(-1:im(region)+2,-1:jm(region)+2)  :: phelp
    real, allocatable                 :: phelp(:,:)
    integer                           :: i, j, l, imax, xref_, yref_
    integer                           :: imr, jmr, ipar, jpar, lmr
    real                              :: maxdiff, m_edge, maxdiff_all
    real                              :: mnew
    real,dimension(:,:,:),allocatable :: mn, mo
    integer                           :: io, sfstart, sfend
    integer                           :: communicator, ll, root_id

    ! start

#ifdef MPI
    which_par=previous_par(region)
    if ( which_par == 'tracer' .and. ntracetloc == 0) return
    if ( which_par == 'levels' .and. lmloc == 0) return  !WP!
#endif

#ifdef MPI
    if ( which_par == 'tracer' ) then
#endif
       p    =>   sp_dat(region)%data
       m    =>    m_dat(region)%data
       phlb => phlb_dat(region)%data
       mpar =>    m_dat(region)%data
       dxyp => region_dat(region)%dxyp
       dxy  => region_dat(parent(region))%dxyp

       lmr = lm(region)
#ifdef MPI
       communicator=com_trac  !WP! assign com_trac as communicator
       root_id=root_t
    else if ( which_par == 'levels' ) then
       p    =>   sp_dat(region)%data
       m    =>    m_dat(region)%data_k
       phlb => phlb_dat(region)%data_k
       dxyp => region_dat(region)%dxyp
       dxy  => region_dat(parent(region))%dxyp
       mpar => mass_dat(parent(region))%m_k

       lmr = lmloc
       communicator=com_lev  !WP! assign com_lev as communicator
       root_id=root_k
    end if
#endif

    yref_ = yref(region)/yref(parent(region)) ! will give the number of latitudinal interface cells in the region
    xref_ = xref(region)/xref(parent(region)) ! will give the number of longitudinal interface cells in the region
    imr = im(region)
    jmr = jm(region)
    ! check mass in the interface region....
    if ( okdebug ) then
       maxdiff = zero
       do l=1,lmr
          do j=1,jmr-yref_+1,yref_
             do i=1,imr-xref_+1,xref_
                if ((i > xref_ .and. i < imr-xref_+1) .and.  &
                     (j > yref_ .and. j < jmr-yref_+1) ) cycle  ! not on edge...
                if ((xcyc(region) == 1) .and.  &
                     (j > yref_ .and. j < jmr-yref_+1) ) cycle  ! periodic boundary
                if ((touch_np(region) == 1).and.(j > jmr-yref_)) cycle ! np_touching
                if ((touch_sp(region) == 1).and.(j < yref_+1)  ) cycle ! sp_touching
                ipar = ibeg(region) + (i-1)/xref_
                jpar = jbeg(region) + (j-1)/yref_   !to get mass in the parent cell
                m_edge = sum(m(i:i+xref_-1,j:j+yref_-1,l)) !get mass in zoom region
                if ( abs((m_edge-mpar(ipar,jpar,l))/m_edge) > 1e-8) then
                   print*,'m2phlb: Encountered mass-difference at',i,j,l
                   print*,'m2phlb: Parent mass comes from        ',ipar,jpar,l
                   print*,'m2phlb: Parent mass is                ',mpar(ipar,jpar,l)
                   print*,'m2phlb: Child  mass is                ',m_edge
                   print*,'m2phlb: Made out of                   ', &
                        m(i:i+xref_-1,j:j+yref_-1,l)
                end if
                maxdiff = max(maxdiff,100*abs((m_edge-mpar(ipar,jpar,l))/m_edge))
             end do
          end do
       end do
       print*,'m2phlb: maxdiff mass  on edges (%)',maxdIFf, 'region:', region
    end if

    allocate(phelp(-1:im(region)+2,-1:jm(region)+2))

    phelp = zero
    p = zero
    !$omp parallel do schedule(guided) private(j)
    do j = 1, jmr
        phelp(1:imr, j) = sum(m(1:imr, j, 1:lmr), dim=2) * grav/dxyp(j) ! phelp contains the surface pressure computed from air mass
    end do
    !$omp end parallel do
    !do l=1,lmr
       !do j=1,jmr
          !do i=1,imr
             !phelp(i,j) = phelp(i,j) + m(i,j,l)*grav/dxyp(j)
          !end do
       !end do
    !end do

#ifdef MPI
    if ( which_par == 'tracer' ) then
#endif
       p(:,:,1) = phelp ! lowest layer of p contains the newly computed surface pressure
#ifdef MPI
    else if ( which_par == 'levels' ) then
       call mpi_allreduce(phelp,p,(imr+4)*(jmr+4),my_real,mpi_sum,com_lev,ierr)
    end if
#endif
    ! correct pressure edges for coarse grid area...
    ! MKADJ: one exception: in adjoint pressure is restored
    ! and is already corrected for this if newsrun=.TRUE.
       if(revert /= -1 .or. .not. newsrun) then
       do j=1,jmr
          do i=1,imr
             if ((i > xref_ .and. i < imr-xref_+1) .and.  &
                  (j > yref_ .and. j < jmr-yref_+1) ) cycle    !not on edge...
             if ((xcyc(region) == 1) .and.  &
                  (j > yref_ .and. j < jmr-yref_+1) ) cycle    !periodic boundary
             if ((touch_np(region) == 1).and.(j > jmr-yref_)) cycle    !np_touching
             if ((touch_sp(region) == 1).and.(j < yref_+1)  ) cycle    !sp_touching
             jpar = jbeg(region) + (j-1)/yref_   !to get dxyp in the parent cell
             p(i,j,1) = p(i,j,1)*dxyp(j)*xref_*yref_/dxy(jpar)
          end do
       end do
    endif

    !$omp parallel do schedule(guided) private(l, ll)
    do  l=1,lmr+1
       ll=l
#ifdef MPI
       ! offset for levels
       if ( myid > 0 .and. which_par == 'levels' ) ll=sum(lmar(0:myid-1))+l
#endif
       phlb(1:imr, 1:jmr, l) = at(ll) + bt(ll) * p(1:imr, 1:jmr, 1)
       !do  j=1,jmr
          !do  i=1,imr
             !phlb(i,j,l) = at(ll)+bt(ll)*p(i,j,1)
          !end do
       !end do
    end do
    !$omp end parallel do

    if ( okdebug ) then
       allocate(mn(imr,jmr,lmr))
       allocate(mo(imr,jmr,lmr))
       mn = 0.0    !omitting might produce crashes on some systems (PB)
       mo = 0.0
       maxdiff = 0.0
       do l = 1,lmr
          jloop: do  j=1,jmr
             if ( yref_> 1) then
                ! not valid on edge..
                if ( j < yref_+1 .and. touch_sp(region) /= 1) cycle jloop
                if ( j > jmr-yref_ .and. touch_np(region) /= 1) cycle jloop
             endif
             iloop : do  i=1,imr
                if ( xcyc(region) /= 1 .and. (i<xref_+1)) cycle iloop
                if ( xcyc(region) /= 1 .and. (i>imr-xref_)) cycle iloop
                mnew=(phlb(i,j,l)-phlb(i,j,l+1))*dxyp(j)/grav
                mn(i,j,l) = mnew
                mo(i,j,l) = m(i,j,l)
                maxdiff = max(maxdiff,abs(mnew-m(i,j,l))/mnew)
             end do iloop
          end do jloop
       end do

       !WP! execute an allreduce over proper processors (communicator)

#ifdef MPI
       call mpi_allreduce(maxdiff,maxdiff_all,1,my_real,mpi_max,communicator,ierr)

       if ( myid == root_id) &
            print*,'m2phlb: maxdIFf (mnew-m)/mnew on all pes:',maxdiff_all
       if ( maxdiff_all > 1.e-4)then
          call escape_tm('problem with m2phlb: maxdiff > 1e-4')
       endif
#else
       print*,'m2phlb: maxdIFf (mnew-m)/mnew:',maxdiff
       if ( maxdiff > 1.e-4)then
          call escape_tm('problem with m2phlb: maxdiff > 1e-4')
       endif
#endif

       deallocate(mn)
       deallocate(mo)

    end if

    nullify(p)
    nullify(m)
    nullify(phlb)
    nullify(dxyp)
    nullify(dxy)
    nullify(mpar)

    deallocate(phelp)

  end subroutine m2phlb



#ifdef MPI

  subroutine calc_phlb_k(region)
    !
    use dims
    use global_data, only : mass_dat
    use mpi_const,   only : lmar,myid,lmloc

    implicit none

    ! in/out
    integer,intent(in) :: region

    ! local
    real,dimension(:,:,:),pointer  :: phlb
    real,dimension(:,:),pointer    :: p
    integer :: i,j,l,ll

    phlb => mass_dat(region)%phlb_k
    p    => mass_dat(region)%p

    do l=1,lmloc+1
       ll=l
       if ( myid > 0 ) ll = sum(lmar(0:myid-1))+l  !offset for levels
       do  j=1,jm(region)
          do  i=1,im(region)
             phlb(i,j, l) = at(ll)+bt(ll)*p(i,j)
          end do
       end do
    end do

    nullify(phlb)

  end subroutine calc_phlb_k

#endif


  subroutine advect_m(region,ntimes)
    !
    ! subroutine
    !
    use binas,       only : grav
    use dims
    use global_data, only : mass_dat
    use global_data, only : region_dat
    use MeteoData  , only : sp_dat, m_dat
    use global_data, only : wind_dat
#ifdef MPI
    use mpi_const,   only : ntracetloc,myid,root
#endif

    implicit none

    ! input/output
    integer,intent(in) :: region,ntimes

    ! local
    real,dimension(:,:,:),pointer      :: m,am,bm,cm
    real,dimension(:,:,:),pointer      :: p
    real,dimension(:)    ,pointer      :: dxyp

    real,dimension(:,:,:),allocatable                :: mnew, msave
    real,dimension(:,:)  ,allocatable                :: pnew
    integer   :: imr
    integer   :: jmr
    integer   :: lmr
    integer   :: n,i,j,l
    real      :: ss,so

    ! start

#ifdef MPI
    !WP! only processors with nonzero ntracet allowed
    if ( ntracetloc == 0 ) return
    !if(myid==root) print*, 'advect_m: called for region ', region
    !if(myid==root) print*, 'advect_m: ntimes            ', ntimes
#else
    !print*, 'advect_m: called for region ', region
    !print*, 'advect_m: ntimes            ', ntimes
#endif

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    allocate(mnew(imr,jmr,lmr))
    allocate(msave(imr,jmr,lmr))
    allocate(pnew(imr,jmr))

    m =>  m_dat(region)%data
    p => sp_dat(region)%data
    am => wind_dat(region)%am_t
    bm => wind_dat(region)%bm_t
    cm => wind_dat(region)%cm_t
    dxyp => region_dat(region)%dxyp

    msave = m(1:imr,1:jmr,1:lmr)
    mnew  = m(1:imr,1:jmr,1:lmr)
    do i=1,ntimes
       mnew = mnew   + revert*am(0:imr-1,1:jmr  ,1:lmr  ) &   !east   !ADJ re-inserted the revert CMK
                     - revert*am(1:imr  ,1:jmr  ,1:lmr  ) &   !west
                     + revert*bm(1:imr  ,1:jmr  ,1:lmr  ) &   !south
                     - revert*bm(1:imr  ,2:jmr+1,1:lmr  ) &   !north
                     + revert*cm(1:imr  ,1:jmr  ,0:lmr-1) &   !lower
                     - revert*cm(1:imr  ,1:jmr  ,1:lmr  )     !upper
    end do
    pnew = 0.0
    do l=1,lm(region)
       do j=1,jm(region)
          do i=1,im(region)
             pnew(i,j) = pnew(i,j) + mnew(i,j,l)*grav/dxyp(j)
          end do
       end do
    end do

    p(1:imr,1:jmr,1) = pnew   !return new pressure for testing
    m(1:imr,1:jmr,1:lmr) = mnew
    nullify(m)
    nullify(am)
    nullify(bm)
    nullify(cm)
    nullify(p)
    nullify(dxyp)

    deallocate(pnew)
    deallocate(mnew)
    deallocate(msave)

  end subroutine advect_m



  subroutine coarsen_ambm(region)
   !
   ! coarsens boundary values of am,bm of all children regions
   ! in corresponding with region's am,bm
   ! written by mike botchev, march-june 1999
   ! updated by maarten krol, dec 2002
   !
   use dims
   use global_data, only: wind_dat
#ifdef MPI
   use mpi_comm
   use mpi_const
#endif

   implicit none

   ! input
   integer,intent(in) :: region

   ! local
   real,dimension(:,:,:),pointer             :: am,bm,amc,bmc
   real cell1,cellm
   integer child,i,ic,ichild,j,jc,l,lc
   integer xref_,yref_,zref_

   ! start

   if ( okdebug ) print *,'coarsen_ambm:'

   am => wind_dat(region)%am_t
   bm => wind_dat(region)%bm_t

   ichild = 0
   do while(ichild<children(region,0))
      ichild = ichild + 1
      child = children(region,ichild)
      xref_ = xref(child)/xref(region)
      yref_ = yref(child)/yref(region)
      zref_ = zref(child)/zref(region)
      !WP! define these on all PE's as we will scatter them later
      !WP! amck holds am of child over levels

      amc => wind_dat(child)%am_t

#ifdef MPI
      if ( ntracetloc /= 0 ) then !WP! only PE's with nonzero ntracetloc
#endif
         if ( okdebug ) print *,'coarsen_ambm: coarsening am,bm in child ',child

         ! write am: loop through the cells of x-walls of the child
         ! do this only if the region is not 360-degrees wide
         if ( im(child)/xref(child) /= im(1) ) then


            do j=jbeg(child),jend(child)
               do l=lbeg(child),lend(child)

                  jc = (j-jbeg(child))*yref_ + 1    ! jc is i_child
                  lc = (l-lbeg(child))*zref_ + 1    ! lc is l_child

                  cell1 = am(ibeg(child)-1,j,l)/(yref_*zref_)
                  cellm = am(iend(child)  ,j,l)/(yref_*zref_)

                  amc(xref_-1,  jc:jc+yref_-1, lc:lc+zref_-1) = cell1
                  amc(im(child)-xref_+1,  jc:jc+yref_-1, lc:lc+zref_-1) = cellm
               end do
            end do

         end if
#ifdef MPI
      end if
#endif
!PB      nullify(am)
      nullify(amc)

      bmc => wind_dat(child)%bm_t
#ifdef MPI
      if ( ntracetloc /= 0 ) then !WP! only PE's with nonzero ntracetloc
#endif

         ! write bm: loop through the cells of y-walls of the child

         do i=ibeg(child),iend(child)
            do l=lbeg(child),lend(child)

               ic = (i-ibeg(child))*xref_ + 1    ! ic is i_child
               lc = (l-lbeg(child))*zref_ + 1    ! lc is l_child

               cell1 = bm(i,jbeg(child)  ,l)/(xref_*zref_)
               cellm = bm(i,jend(child)+1,l)/(xref_*zref_)
               ! do this only if the region does not touch the SP
               if ( touch_sp(child) /= 1 )   &
                    bmc(ic:ic+xref_-1,yref_,lc:lc+zref_-1)=cell1
               ! idem for NP...
               if ( touch_np(child) /= 1 )  &
                    bmc(ic:ic+xref_-1,jm(child)-yref_+2,lc:lc+zref_-1)=cellm
            end do
         end do
#ifdef MPI
      end if
#endif
!PB       nullify(bm)
      nullify(bmc)
   end do

!PBi
  nullify(am)
  nullify(bm)
!PBe

 end subroutine coarsen_ambm




  ! ***


  subroutine dynam0( region, ndyn, status )
    !-----------------------------------------------------------------------
    !
    !****  dynam0          - calculates vertical massfluxes        v 9.1
    !
    !      programmed by           mh      mpi HH          1-oct-1991
    !
    !      purpose
    !      -------
    !      This sr is CALLed each time after new massfluxes are read in.
    !      Calculate the new air-masses in each grid box from the
    !      surface pressure.
    !      Calculate the amount of air moved in each advection substep,
    !      also in the vertical direction.
    !
    !      interface
    !      ---------
    !      CALL dynam0
    !
    !      method
    !      ------
    !      integrate air conservation equation in the vertical direction
    !      in order to obtain the vertical massfluxes.
    !
    !      externals
    !      ---------
    !      none
    !
    !      reference
    !      ---------
    !      see manual
    !
    ! modified by mk (1999) for zoom version.
    !-----------------------------------------------------------------------

    use dims,        only : im, jm, lm
    use dims,        only : tref, zero, bt, xcyc
    use MeteoData  , only : pu_dat, pv_dat
    use global_data, only : wind_dat
#ifdef MPI
    use mpi_const,   only : ntracetloc
    use mpi_comm,    only : barrier_t
#endif

    implicit none

    ! input
    integer,intent(in)   :: region
    integer,intent(in)   :: ndyn
    integer,intent(out)  :: status

    ! local
    real,dimension(:,:,:),pointer           :: pu,pv,am,bm,cm
    real,dimension(im(region),jm(region))   :: pit
    real,dimension(im(region),jm(region),lm(region)-1)   :: sd
    real,dimension(im(region),jm(region),lm(region))     :: conv_adv

    real    :: dtu, dtv, dtw
    integer :: i, j, l, lmr, imr, jmr

    ! start

#ifdef MPI
    if ( ntracetloc == 0 ) return !WP! only processors with nonzero ntracet
#endif

    pu => pu_dat(region)%data
    pv => pv_dat(region)%data
    am => wind_dat(region)%am_t
    bm => wind_dat(region)%bm_t
    cm => wind_dat(region)%cm_t

    !---- length of advection substeps (in seconds) in the three directions
    dtu=ndyn/(2.*tref(region))  ! again removed the revert! cmk
    dtv=ndyn/(2.*tref(region))
    dtw=ndyn/(2.*tref(region))

    imr = im(region) ; jmr = jm(region) ; lmr = lm(region)

    !$omp parallel
    !$omp sections
    !$omp section
    am(0:imr+1,0:jmr+1,0:lmr+1) = zero
    !$omp section
    bm(0:imr+1,0:jmr+1,0:lmr+1) = zero
    !$omp section
    cm(0:imr+1,0:jmr+1,0:lmr+1) = zero
    !$omp section
    !
    !     compute conv_adv, the horizontal mass convergence
    !      the arrays pu and pv contain the horizontal air mass fluxes crossing
    !      the grid box boundaries. pu(i,j,l) is the eastward mass flux in kg/sec
    !      at the eastern edge of box i,j,l; pv(i,j,l) is the northward
    !      mass flux at the southern edge of box i,j,l
    !
    conv_adv(1:imr,1:jmr,1:lmr) = pu(0:imr-1,1:jmr,1:lmr) - pu(1:imr,1:jmr,1:lmr) + pv(1:imr,1:jmr,1:lmr) - pv(1:imr,2:jmr+1,1:lmr)
    !do l=1,lmr
       !do  j=1,jm(region)
          !do  i=1,im(region)
             !conv_adv(i,j,l)=pu(i-1,j,l)-pu(i,j,l)+pv(i,j,l)-pv(i,j+1,l)
          !end do
       !end do
    !end do
    !$omp end sections
    !$omp end parallel
    !
    !     compute pit, the vertiCALLy integrated conv_adv
    !
    pit = sum(conv_adv, dim=3)
    sd(:,:,lmr-1) = conv_adv(:,:,lmr) - (bt(lmr)-bt(lmr+1)) * pit
    do l = lmr-2,1,-1
        sd(:,:,l) = sd(:,:,l+1) + conv_adv(:,:,l+1) - (bt(l+1)-bt(l+2)) * pit
    end do

    !do j=1,jm(region)
       !do  i=1,im(region)
          !!pit(i,j) = sum(conv_adv(i,j,:))
          !!!pit(i,j)=conv_adv(i,j,1)
          !!!do l=2,lmr
             !!!pit(i,j)=pit(i,j)+conv_adv(i,j,l)
          !!!end do
          !!
          !!     compute the vertical massflux on the box boundaries (in kg/s)
          !!mk  in the hybrid system, the tendency in ps is given by bt. Bt is
          !!mk  already normalized. Now distribute pit such that the new mass
          !!mk  distribution after advection exactly matches the
          !!mk  exactly coincides with the distribution that is obtained
          !!mk  from the new surface pressure.
          !!
          !sd(i,j,lmr-1)=conv_adv(i,j,lmr) - (bt(lmr)-bt(lmr+1))*pit(i,j)
          !do  l=lmr-2,1,-1
             !sd(i,j,l)=sd(i,j,l+1)+conv_adv(i,j,l+1)-(bt(l+1)-bt(l+2))*pit(i,j)
          !end do
       !end do
    !end do
    !
    !     compute amount of air moved each advection substep, am, bm or cm (kg)
    !
    !      am(i,j,l) is the amount of air moved eastward at the eastern
    !      boundary of grid box i,j,l
    !     compute cm
    !      cm(i,j,l) is the amount of air moved in upward direction at the
    !      upper boundary of grid box i,j,l
    !     compute bm
    !      bm(i,j,l) is the amount of air moved northward at the southern
    !      boundary of grid box i,j,l
    !
    !$omp parallel
    !$omp sections
    !$omp section
    am(0:imr,1:jmr,1:lmr) = dtu * pu(0:imr,1:jmr,1:lmr)
    !do l=1,lmr
       !do j=1,jm(region)
          !do i=0,im(region)
             !am(i,j,l)=dtu*pu(i,j,l)
          !end do
       !end do
    !end do
    !$omp section
    bm(1:imr,1:jmr+1,1:lmr) = dtv * pv(1:imr,1:jmr+1,1:lmr)
    !do l=1,lmr
       !do j=1,jm(region)+1
          !do i=1,im(region)
             !bm(i,j,l)=dtv*pv(i,j,l)
          !end do
       !end do
    !end do
    !$omp section
    cm(1:imr,1:jmr,1:lmr-1) = -dtw * sd(1:imr,1:jmr,1:lmr-1)
    !do l=1,lmr-1
       !do j=1,jm(region)
          !do i=1,im(region)
             !cm(i,j,l)=-dtw*sd(i,j,l)
          !end do
       !end do
    !end do
    !$omp end sections
    !$omp end parallel

    if ( xcyc(region) == 1 ) then
       am(0,:,:) = am(imr,:,:)
       am(imr+1,:,:) = am(1,:,:)
    end if

    nullify(pu)
    nullify(pv)
    nullify(am)
    nullify(bm)
    nullify(cm)

#ifdef MPI
    call barrier_t
#endif

    ! ok
    status = 0

  end subroutine dynam0

end module advect_tools
