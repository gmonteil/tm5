!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_zoom_tools

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- interface ---

  private

  public :: adj_mix_edges, adj_update_parent

  ! --- const ---

  character(len=*), parameter  ::  mname = 'adj_zoom_tools'


contains

!===========================================================================================================
!===========================================================================================================

  !=============================================================
  ! - adj_mix_edges is self-adjoint
  ! - however: reconstruction of required varialbles necessary
  !=============================================================

  subroutine adj_mix_edges( region, status )

    use zoom_tools, only:  mix_edges

    !__IO____________________________________________________________________

    integer,intent(in)            ::  region
    integer,intent(out)           ::  status

    !__CONST_________________________________________________________________

    character(len=*), parameter  ::  rname = mname//'adj_mix_edges'

    !__LOCAL_VARIABLES_______________________________________________________


    !__START_SUBROUTINE______________________________________________________

    call mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    call reconstruct_required_var_1(region)

    ! ok
    status = 0

  end subroutine adj_mix_edges

!===========================================================================================================
!===========================================================================================================

  subroutine reconstruct_required_var_1(region)

    use dims
    use global_data, only: wind_dat, mass_dat
    use MeteoData  , only : m_dat

    implicit none

    !__IO____________________________________________________________________

    integer,intent(in)            ::  region

    !__LOCAL_VARIABLES_______________________________________________________

    real,dimension(:,:,:), pointer   :: mc,mp,amc,bmc,cmc,amp,bmp,cmp
    integer                          :: my_parent,i,j,l, xref_, yref_, imc,jmc,imp,jmp,q
    integer                          :: jp,ipw,ipe,step,jps,jpn,ip
    character(len=n_operators)       :: stencil
    logical                          :: xyz
    real                             :: mpw, mpe, mps, mpn   ! parent masses east/west...south/north


    !__START_SUBROUTINE______________________________________________________



! mix_edges nevertheless distroys some information about the air masses.
! mix_edges is called after each advection step. The masses alon the
! x- and y-edges are 'mixed', which involves:
!    o    m = sum over xref_, sum over yref_ (m_edge)   / xref_*yref_
!
! In the adjoint run, the masses on the edges have to be reconstructed and that is
! done here. It involves:
!     o a backintegration of the parent masses
!     o a partial forward integration
! The forward integration depends on
!     o whether the parent has XYZ or ZYX sequence
!     o the position in the child processes.
!

    my_parent = parent(region)
    if(region == 1.or.my_parent ==0) return

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    imp = im(my_parent) ; jmp = jm(my_parent)
    imc = im(region)    ; jmc = jm(region)
    amc => wind_dat(region)%am_t
    bmc => wind_dat(region)%bm_t
    cmc => wind_dat(region)%cm_t
    mc  => m_dat(region)%data
    amp => wind_dat(my_parent)%am_t
    bmp => wind_dat(my_parent)%bm_t
    cmp => wind_dat(my_parent)%cm_t
    mp  => m_dat(my_parent)%data
    q=(status(my_parent)-1)/((nsplitsteps/2))      ! find q - the place in the splitorderzoom of the parent
    do i=1,n_operators
      stencil(i:i)=splitorderzoom(my_parent,q*(nsplitsteps/2)+i)
    enddo
    ! find out if parent has order xyz, or zyx...
    if( scan(stencil,'x') >  scan(stencil,'z')) then
      xyz = .false.
    else
      xyz = .true.
    endif

    ! find out 'step' the position in the xyzzyx / zyxxyz operation sequence (1...6)
    if (xyz) then
      if( mod(status(region)-1,nsplitsteps) < n_operators ) then  ! first step of xyz...zyx
        select case(splitorderzoom(region,status(region)))
        case('x')
          step = 1
        case('y')
          step = 2
        case('z')
          step = 3
        end select
      else
        select case(splitorderzoom(region,status(region)))
        case('x')
          step = 6
        case('y')
          step = 5
        case('z')
          step = 4
        end select
      endif
    else
      if( mod(status(region)-1,nsplitsteps) < n_operators ) then  ! first step of xyz...zyx
        select case(splitorderzoom(region,status(region)))
        case('x')
          step = 3
        case('y')
          step = 2
        case('z')
          step = 1
        end select
      else
        select case(splitorderzoom(region,status(region)))
        case('x')
          step = 4
        case('y')
          step = 5
        case('z')
          step = 6
        end select
      endif
    endif

    if (okdebug) print *, 'in adj_mix_edges (region, parent)', region,my_parent, ' has ',stencil,xyz
    if (okdebug) print *, 'in adj_mix_edges ', splitorderzoom(region,status(region)),mod(status(region)-1,nsplitsteps),step

    if (xyz) then   ! first put_xedges, then put_yedges
      ! x-edges:
      if (xcyc(region)==0) then

        do l=1,lm(region)
          do j=1,jmc
            jp = jbeg(region) + (j-1)/yref_
            ipw = ibeg(region)
            ipe = iend(region)

            select case(step)    ! calculate the parent mass that forms the basis for the child:
            case(1)   ! x: advect parent fully back.
              mpw = mp(ipw,jp,l) + amp(ipw,jp,l) - amp(ipw-1,jp,l) - bmp(ipw,jp,l) + bmp(ipw,jp+1,l) &
                                                                   + cmp(ipw,jp,l) - cmp(ipw,jp,l-1)
              mpw = (mpw/yref_ + amc(xref_-1,j,l) - amc(xref_,j,l))/xref_
              mc(1:xref_,j,l) = mpw
              mpe = mp(ipe,jp,l) + amp(ipe,jp,l) - amp(ipe-1,jp,l) - bmp(ipe,jp,l) + bmp(ipe,jp+1,l) &
                                                                   + cmp(ipe,jp,l) - cmp(ipe,jp,l-1)
              mpe = (mpe/yref_ - amc(imc-xref_+1,j,l) + amc(imc-xref_,j,l))/xref_
              mc(imc-xref_+1:imc,j,l) = mpe
            case(2)   ! y
              mpw = mp(ipw,jp,l) + 0.5*amp(ipw  ,jp,l) - bmp(ipw,jp,l) + bmp(ipw,jp+1,l) + cmp(ipw,jp,l) - cmp(ipw,jp,l-1)
              mpw = mpw/(xref_*yref_)
              mpe = mp(ipe,jp,l) - 0.5*amp(ipe-1,jp,l) - bmp(ipe,jp,l) + bmp(ipe,jp+1,l) + cmp(ipe,jp,l) - cmp(ipe,jp,l-1)
              mpe = mpe/(xref_*yref_)
              do i=1,xref_
                mc(i,j,l) = mpw + bmc(i,j,l)-bmc(i,j+1,l)
                mc(imc-xref_+i,j,l) = mpe + bmc(imc-xref_+i,j,l)-bmc(imc-xref_+i,j+1,l)
              enddo
            case(3)   ! z
              mpw = mp(ipw,jp,l) + 0.5*amp(ipw  ,jp,l) - 0.5*bmp(ipw,jp,l) + 0.5*bmp(ipw,jp+1,l) + cmp(ipw,jp,l) - cmp(ipw,jp,l-1)
              mpw = mpw/(xref_*yref_)
              mpe = mp(ipe,jp,l) - 0.5*amp(ipe-1,jp,l) - 0.5*bmp(ipe,jp,l) + 0.5*bmp(ipe,jp+1,l) + cmp(ipe,jp,l) - cmp(ipe,jp,l-1)
              mpe = mpe/(xref_*yref_)
              do i=1,xref_
                mc(i,j,l) = mpw - cmc(i,j,l) + cmc(i,j,l-1)
                mc(imc-xref_+i,j,l) = mpe - cmc(imc-xref_+i,j,l) + cmc(imc-xref_+i,j,l-1)
              enddo
            case(4)   ! z
              mpw = mp(ipw,jp,l) + 0.5*amp(ipw  ,jp,l) - 0.5*bmp(ipw,jp,l) + 0.5*bmp(ipw,jp+1,l) &
                                                       + 0.5*cmp(ipw,jp,l) - 0.5*cmp(ipw,jp,l-1)
              mpw = mpw/(xref_*yref_)
              mpe = mp(ipe,jp,l) - 0.5*amp(ipe-1,jp,l) - 0.5*bmp(ipe,jp,l) + 0.5*bmp(ipe,jp+1,l) &
                                                       + 0.5*cmp(ipe,jp,l) - 0.5*cmp(ipe,jp,l-1)
              mpe = mpe/(xref_*yref_)
              do i=1,xref_
                mc(i,j,l) = mpw - cmc(i,j,l) + cmc(i,j,l-1)
                mc(imc-xref_+i,j,l) = mpe - cmc(imc-xref_+i,j,l) + cmc(imc-xref_+i,j,l-1)
              enddo
            case(5)   ! y
              mpw = mp(ipw,jp,l) + 0.5*amp(ipw  ,jp,l) - 0.5*bmp(ipw,jp,l) + 0.5*bmp(ipw,jp+1,l)
              mpw = mpw/(xref_*yref_)
              mpe = mp(ipe,jp,l) - 0.5*amp(ipe-1,jp,l) - 0.5*bmp(ipe,jp,l) + 0.5*bmp(ipe,jp+1,l)
              mpe = mpe/(xref_*yref_)
              do i=1,xref_
                mc(i,j,l) = mpw + bmc(i,j,l)-bmc(i,j+1,l)
                mc(imc-xref_+i,j,l) = mpe + bmc(imc-xref_+i,j,l)-bmc(imc-xref_+i,j+1,l)
              enddo
            case(6)
              mpw = mp(ipw,jp,l) + 0.5*amp(ipw  ,jp,l)
              mpw = (mpw/yref_ - amc(xref_,j,l))/xref_
              mc(1:xref_,j,l) = mpw
              mpe = mp(ipe,jp,l) - 0.5*amp(ipe-1,jp,l)
              mpe = (mpe/yref_  + amc(imc-xref_,j,l))/xref_
              mc(imc-xref_+1:imc,j,l) = mpe
            end select
          enddo !j
        enddo !l
      endif !xcyc...

      ! now the y-edges:
      do l=1,lm(region)
        do i=1,imc
          ip = ibeg(region) + (i-1)/xref_
          jps = jbeg(region)
          jpn = jend(region)
          select case(step)
          case(2)   ! y
            if(touch_sp(region)==0) then
              mps = mp(ip,jps,l) - bmp(ip,jps,l) + bmp(ip,jps+1,l) + cmp(ip,jps,l) - cmp(ip,jps,l-1)
              mps = (mps/xref_ + bmc(i,yref_,l) - bmc(i,yref_+1,l) )/yref_
              mc(i,1:yref_,l) = mps
            endif
            if(touch_np(region)==0) then
              mpn = mp(ip,jpn,l) - bmp(ip,jpn,l) + bmp(ip,jpn+1,l) + cmp(ip,jpn,l) - cmp(ip,jpn,l-1)
              mpn = (mpn/xref_ - bmc(i,jmc - yref_ + 2,l) + bmc(i,jmc - yref_ + 1,l) )/yref_
              mc(i,jmc-yref_+1:jmc,l) = mpn
            endif
          case(3)   ! z
            mps = mp(ip,jps,l) + 0.5*bmp(ip,jps+1,l) + cmp(ip,jps,l) - cmp(ip,jps,l-1)
            mps = mps/(xref_*yref_)
            mpn = mp(ip,jpn,l) - 0.5*bmp(ip,jpn  ,l) + cmp(ip,jpn,l) - cmp(ip,jpn,l-1)
            mpn = mpn/(xref_*yref_)
            do j=1,yref_
              if(touch_sp(region)==0) mc(i,j,l) = mps - cmc(i,j,l) + cmc(i,j,l-1)
              if(touch_np(region)==0) mc(i,jmc-yref_+j,l) = mpn - cmc(i,jmc-yref_+j,l) + cmc(i,jmc-yref_+j,l-1)
            enddo
          case(4)   ! z
            mps = mp(ip,jps,l) + 0.5*bmp(ip,jps+1,l) + 0.5*cmp(ip,jps,l) - 0.5*cmp(ip,jps,l-1)
            mps = mps/(xref_*yref_)
            mpn = mp(ip,jpn,l) - 0.5*bmp(ip,jpn  ,l) + 0.5*cmp(ip,jpn,l) - 0.5*cmp(ip,jpn,l-1)
            mpn = mpn/(xref_*yref_)
            do j=1,yref_
              if(touch_sp(region)==0) mc(i,j,l) = mps - cmc(i,j,l) + cmc(i,j,l-1)
              if(touch_np(region)==0) mc(i,jmc-yref_+j,l) = mpn - cmc(i,jmc-yref_+j,l) + cmc(i,jmc-yref_+j,l-1)
            enddo
          case(5)   ! y
            if(touch_sp(region)==0) then
              mps = mp(ip,jps,l) + 0.5*bmp(ip,jps+1,l)
              mps = (mps/xref_ - bmc(i,yref_+1,l) )/yref_
              mc(i,1:yref_,l) = mps
            endif
            if( touch_np(region)==0) then
              mpn = mp(ip,jpn,l) - 0.5*bmp(ip,jpn  ,l)
              mpn = (mpn/xref_ + bmc(i,jmc - yref_ + 1,l) )/yref_
              mc(i,jmc-yref_+1:jmc,l) = mpn
            endif
          case default
          end select
        enddo
      enddo
    else   !zyx
      ! now first the y-edges:
      do l=1,lm(region)
        do i=1,imc
          ip = ibeg(region) + (i-1)/xref_
          jps = jbeg(region)
          jpn = jend(region)
          select case(step)
          case(2)   ! y
            if( touch_sp(region)==0) then
              mps = mp(ip,jps,l) - bmp(ip,jps,l) + bmp(ip,jps+1,l) + amp(ip,jps,l) - amp(ip-1,jps,l)
              mps = (mps/xref_ + bmc(i,yref_,l) - bmc(i,yref_+1,l) )/yref_
              mc(i,1:yref_,l) = mps
            endif
            if( touch_np(region)==0) then
              mpn = mp(ip,jpn,l) - bmp(ip,jpn,l) + bmp(ip,jpn+1,l) + amp(ip,jpn,l) - amp(ip-1,jpn,l)
              mpn = (mpn/xref_ - bmc(i,jmc - yref_ + 2,l) + bmc(i,jmc - yref_ + 1,l) )/yref_
              mc(i,jmc-yref_+1:jmc,l) = mpn
            endif
          case(3)   ! x
            mps = mp(ip,jps,l) + 0.5*bmp(ip,jps+1,l) + amp(ip,jps,l) - amp(ip-1,jps,l)
            mps = mps/(xref_*yref_)
            mpn = mp(ip,jpn,l) - 0.5*bmp(ip,jpn  ,l) + amp(ip,jpn,l) - amp(ip-1,jpn,l)
            mpn = mpn/(xref_*yref_)
            do j=1,yref_
              if(touch_sp(region)==0) mc(i,j,l) = mps - amc(i,j,l) + amc(i-1,j,l)
              if(touch_np(region)==0) mc(i,jmc-yref_+j,l) = mpn - amc(i,jmc-yref_+j,l) + amc(i-1,jmc-yref_+j,l)
            enddo
          case(4)   ! x
            mps = mp(ip,jps,l) + 0.5*bmp(ip,jps+1,l) + 0.5*amp(ip,jps,l) - 0.5*amp(ip-1,jps,l)
            mps = mps/(xref_*yref_)
            mpn = mp(ip,jpn,l) - 0.5*bmp(ip,jpn  ,l) + 0.5*amp(ip,jpn,l) - 0.5*amp(ip-1,jpn,l)
            mpn = mpn/(xref_*yref_)
            do j=1,yref_
              if(touch_sp(region)==0) mc(i,j,l) = mps - amc(i,j,l) + amc(i-1,j,l)
              if(touch_np(region)==0) mc(i,jmc-yref_+j,l) = mpn - amc(i,jmc-yref_+j,l) + amc(i-1,jmc-yref_+j,l)
            enddo
          case(5)   ! y
            if( touch_sp(region)==0) then
              mps = mp(ip,jps,l) + 0.5*bmp(ip,jps+1,l)
              mps = (mps/xref_ - bmc(i,yref_+1,l) )/yref_
              mc(i,1:yref_,l) = mps
            endif
            if( touch_np(region)==0) then
              mpn = mp(ip,jpn,l) - 0.5*bmp(ip,jpn  ,l)
              mpn = (mpn/xref_ + bmc(i,jmc - yref_ + 1,l) )/yref_
              mc(i,jmc-yref_+1:jmc,l) = mpn
            endif
          case default
          end select
        enddo
      enddo
      ! and then the xedges:
      if (xcyc(region)==0) then

        do l=1,lm(region)
          do j=1,jmc
            jp = jbeg(region) + (j-1)/yref_
            ipw = ibeg(region)
            ipe = iend(region)
            select case(step)    ! calculate the parent mass that forms the basis for the child:
            case(3)   ! x
              mpw = mp(ipw,jp,l) + amp(ipw,jp,l) - amp(ipw-1,jp,l)
              mpw = (mpw/yref_ + amc(xref_-1,j,l) - amc(xref_,j,l))/xref_
              mc(1:xref_,j,l) = mpw
              mpe = mp(ipe,jp,l) + amp(ipe,jp,l) - amp(ipe-1,jp,l)
              mpe = (mpe/yref_ - amc(imc-xref_+1,j,l) + amc(imc-xref_,j,l))/xref_
              mc(imc-xref_+1:imc,j,l) = mpe
            case(4)   ! x
              mpw = mp(ipw,jp,l) + 0.5*amp(ipw  ,jp,l)
              mpw = (mpw/yref_ - amc(xref_,j,l))/xref_
              mc(1:xref_,j,l) = mpw
              mpe = mp(ipe,jp,l) - 0.5*amp(ipe-1,jp,l)
              mpe = (mpe/yref_  + amc(imc-xref_,j,l))/xref_
              mc(imc-xref_+1:imc,j,l) = mpe
            case default
            end select
          enddo !j
        enddo !l
      endif !xcyc...

    endif

    nullify(amc)
    nullify(bmc)
    nullify(cmc)
    nullify(mc)
    nullify(amp)
    nullify(bmp)
    nullify(cmp)
    nullify(mp)

  end subroutine reconstruct_required_var_1

!===========================================================================================================
!===========================================================================================================


  subroutine adj_update_parent(region)
    !
    ! update rm of its parent in correspondence with itself
    ! written by mike botchev, march-june 1999
    !
    use dims
    use global_data  ,only : mass_dat
    use MeteoData   , only : m_dat
    use chem_param,   only : ntracet
    use ParTools  ,   only : ntracetloc
    use toolbox,      only : escape_tm

    implicit none

    !__IO____________________________________________________________________

    integer,intent(in)   :: region

    !__LOCAL_VARIABLES_______________________________________________________


    real,dimension(:,:,:,:),pointer          :: rm,rxm,rym,rzm
    real,dimension(:,:,:,:),pointer          :: rmp,rxmp,rymp,rzmp
    real,dimension(:,:,:),pointer            :: m,mp
    real,dimension(:,:,:),allocatable        :: to_parent

    integer :: my_parent
    integer :: xref_, yref_, zref_
    integer :: i, j, l
    integer :: ip, jp, lp
    integer :: n
    integer :: imp, jmp, lmp
    integer :: imr, jmr, lmr
    integer :: nt


    !__START_SUBROUTINE______________________________________________________


    if (region == 1) return

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    my_parent = parent(region)

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    zref_ = zref(region)/zref(my_parent)

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

    nt=ntracetloc


    if ( okdebug ) then
       print *,'adj_update_parent: my_parent=',my_parent, &
            ' x-,y-,zref_: ',xref_,yref_,zref_
    end if

    imp = im(region)/xref_
    jmp = jm(region)/yref_
    lmp = lmr     ! CMK changed: when paralel over levels lmp can differ from lmr

    ! check imp
    if ( ibeg(region) < iend(region) .and. imp /= iend(region)-ibeg(region)+1 )  &
         call escape_tm('update_parent: stop')

    ! check jmp
    if ( jmp /= jend(region)-jbeg(region)+1 ) &
         call escape_tm('update_parent: stop')


    allocate(to_parent(imp,jmp,lmp))

    ! reconstruct masses
    call reconstruct_required_var_2(region)


    do n=1,nt

      do l=1,lmp
        do j=1,jmp
          do i=1,imp

            to_parent(i,j,l)  =  rmp(ibeg(region)-1+i,jbeg(region)-1+j,l,n)
            rmp(ibeg(region)-1+i,jbeg(region)-1+j,l,n) = 0.0
          enddo
        enddo
      enddo

      do l=1,lmr
        do j=1,jmr
          jp = 1 + (j-1)/yref_
          do i=1,imr
            ip = 1 + (i-1)/xref_
            rm(i,j,l,n) = rm(i,j,l,n) + to_parent(ip,jp,l)
          enddo
        enddo
      enddo  !l


    end do ! n

    ! from Maarten,s adjoint version ?????
    ! call mix_edges(region)

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


  end subroutine adj_update_parent

!===========================================================================================================
!===========================================================================================================


  subroutine reconstruct_required_var_2(region)

    ! from Maarten's adjoint version:
    ! provide the correct air masses to the child...
    ! e.g
    ! region 1      |x|yz         |vsc
    ! region 2           xyzvsccsv
    !          in backward processing, the adj_update should provide the correct masses in
    !          the interface cells of the chlid (thus: at the end of the child operations).
    !          the 'internal' cell masses are conforming the parent masses, but
    !          the calls on the edges (0, im+1, jm+1) are sometines 'partly' advected.
    ! The situation is different for the x-interface and the y-interface (no z-zooming).
    ! In the above example, the x cells (0, imr+1) have 'seen' in the forward:
    !           o put_xedges puts the mass of the parent in the x-interface
    !           o edges remain the same during the full advection
    ! The y-cells (0, jmr+1), however.....
    !           o put_ydeges puts FULLY X-advected masses in 0, jmr+1
    ! In summary:
    !           the adj_update_parent now has to put the following masses in the cells (case parent = XYZ...)
    !           (i = 0, imr+1 (j=1,jm)) : no advection
    !           (j = 0, jmr+1 (i=1,im)) : X
    !
    ! To complicate the things even further, the situation is different when the order is different:
    ! region 1     csvz|y|x      |
    ! region 2             csvzyx
    !           the adj_update_parent now has to put the following masses in the cells (case parent = ZYX...)
    !           (i = 0, imr+1 (j=1,jm)) : ZY
    !           (j = 0, jmr+1 (i=1,im)) : Z
    !
    ! Furthermore, the total masses of the child should match the mass of the underlying parent.....
    !

    use dims
    use chem_param,   only : ntracet
    use ParTools  ,   only : ntracetloc
    use global_data,  only : mass_dat, wind_dat
    use MeteoData   , only : m_dat
    use toolbox,      only : escape_tm

    implicit none

    !__IO____________________________________________________________________

    integer,intent(in)            ::  region

    !__LOCAL_VARIABLES_______________________________________________________


    real,dimension(:,:,:,:),pointer          :: rm,rmp
    real,dimension(:,:,:),pointer            :: m,mp,amp,bmp,cmp

    real ::dummy
    integer i,iend_loop,iip,ip,j,jp,l,lp,my_parent,n,xref_,yref_,zref_,lb,le,ib,ie,jb,je
    integer :: imp,jmp,lmp,imr,jmr,lmr,q
    real m_value,mxval,sum_mass
    integer,dimension(3) :: mxloc
    character(len=n_operators) :: stencil
    logical                    :: xyz


    !__START_SUBROUTINE______________________________________________________

    if (region==1) return

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    my_parent = parent(region)

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    zref_ = zref(region)/zref(my_parent)


    if (okdebug) then
       print *,' adj_updating parent=',my_parent,' x-,y-,zref_: ',xref_, yref_, zref_
    endif

    imp = im(region)/xref_
    jmp = jm(region)/yref_
    lmp = lm(region)/zref_

    if(ibeg(region) < iend(region) .and. imp /= iend(region)-ibeg(region)+1)  &
         call escape_tm('stopped in update_parent')
    if(jmp.ne.jend(region)-jbeg(region)+1) call escape_tm('stopped in update_parent')
    if(lmp.ne.lend(region)-lbeg(region)+1) call escape_tm('stopped in update_parent')

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t

    amp => wind_dat(my_parent)%am_t
    bmp => wind_dat(my_parent)%bm_t
    cmp => wind_dat(my_parent)%cm_t
    mp => m_dat(my_parent)%data
    rmp => mass_dat(my_parent)%rm_t

    ! put mass parent in the child cells, 0, imr+1

    q=(status(my_parent)-1)/((nsplitsteps/2))      ! find q - the place in the splitorderzoom of the parent
    do i=1,n_operators
       stencil(i:i)=splitorderzoom(my_parent,q*(nsplitsteps/2)+i)
    enddo

    ! find out if parent has order xyz, or zyx...
    if( scan(stencil,'x') >  scan(stencil,'z')) then
       xyz = .false.
    else
       xyz = .true.
    endif

    if(xcyc(region)==0 ) then
       do l=1,lmr
          do j=1,jmr
             jp = jbeg(region) + (j-1)/yref_
             if(xyz) then ! x-edges: -  (YZX backwards)
                m(0,j,l) = (mp(ibeg(region)-1,jp,l)   &
                     - bmp(ibeg(region)-1,jp,l  ) + bmp(ibeg(region)-1,jp+1,l) &
                     - cmp(ibeg(region)-1,jp,l-1) + cmp(ibeg(region)-1,jp  ,l) &
                     - amp(ibeg(region)-2,jp,l  ) + amp(ibeg(region)-1,jp  ,l)  )/yref_
                m(imr+1,j,l) = (mp(iend(region)+1,jp,l)  &
                     - bmp(iend(region)+1,jp,l  ) + bmp(iend(region)+1,jp+1,l) &
                     - cmp(iend(region)+1,jp,l-1) + cmp(iend(region)+1,jp  ,l) &
                     - amp(iend(region)  ,jp,l  ) + amp(iend(region)+1,jp  ,l)  )/yref_

             else ! x-edges YZ  (X backwards)
                m(0,j,l) = (mp(ibeg(region)-1,jp,l)       + amp(ibeg(region)-1,jp,l) - amp(ibeg(region)-2,jp,l)  )/yref_
                m(imr+1,j,l) = (mp(iend(region)+1,jp,l)   + amp(iend(region)+1,jp,l) - amp(iend(region)  ,jp,l)  )/yref_
             endif
          enddo
       enddo
    endif !if(xcyc(region)==0)

    if(xcyc(region)==0) then
       do l=1,lmr
          do i=1,imr
             ip = ibeg(region) + (i-1)/xref_

             if(xyz) then    ! y-edges: X  (ZY backwards)

                if(touch_sp(region)==0) then
                   m(i,0,l) = (mp(ip, jbeg(region)-1,l) &
                        - cmp(ip, jbeg(region)-1,l-1) + cmp(ip, jbeg(region)-1,l) &
                        + bmp(ip, jbeg(region)  ,l  ) - bmp(ip, jbeg(region)-1,l)  )/xref_
                endif
                if(touch_np(region)==0) then
                   m(i,jmr+1,l) = (mp(ip, jend(region)+1,l)    &
                        - cmp(ip, jend(region)+1,l-1) + cmp(ip, jend(region)+1,l) &
                        - bmp(ip, jend(region)+1,l)   + bmp(ip, jend(region)+2,l)  )/xref_
                endif

             else  ! y-edges Z  (XY backwards)

                if(touch_sp(region)==0) then
                   m(i,0,l) =  (mp(ip, jbeg(region)-1,l) &
                        - amp(ip-1, jbeg(region)-1,l) + amp(ip, jbeg(region)-1,l) &
                        + bmp(ip,   jbeg(region),  l) - bmp(ip, jbeg(region)-1,l)  )/xref_
                endif
                if(touch_np(region)==0) then
                   m(i,jmr+1,l) = (mp(ip, jend(region)+1,l)  &
                        - amp(ip-1, jend(region)+1,l) + amp(ip, jend(region)+1,l) &
                        - bmp(ip,   jend(region)+1,l) + bmp(ip, jend(region)+2,l)  )/xref_
                endif
             endif
          enddo
       enddo
    endif !if(xcyc(region)==0)


  end subroutine reconstruct_required_var_2

!===========================================================================================================
!===========================================================================================================

end module adj_zoom_tools

!===========================================================================================================
!===========================================================================================================
