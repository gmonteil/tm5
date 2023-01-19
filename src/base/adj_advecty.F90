!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_advecty

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------------

  private

  public  ::  adj_AdvectY_Init, adj_AdvectY_Done
  public  ::  adj_advectyzoom


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'adj_AdvectY'


  ! --- local ------------------------------------

  integer    ::  itim_zoom
  integer    ::  itim_dynam
  integer    ::  itim_compress
  integer    ::  itim_uncompress
  integer    ::  itim_put_yedges


contains


  ! ====================================================================


  subroutine adj_AdvectY_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_AdvectY_Init'

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_zoom, 'adj_advecty zoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynam, 'adj_advecty dynam', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_compress, 'adj_advecty compress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_uncompress, 'adj_advecty uncompress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_put_yedges, 'adj_advecty put yedges', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_AdvectY_Init


  ! ***


  subroutine adj_AdvectY_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_AdvectY_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine adj_AdvectY_Done


  ! ***


  subroutine adj_advectyzoom( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only: xref, yref, zref
    use dims,        only: im, jm ,lm
    use dims,        only: tref
    use dims,        only: parent
    use dims,        only: nsplitsteps, n_operators, splitorderzoom, xcyc
    use dims,        only: rstatus => status
    use dims,        only: adv_scheme
    use dims,        only: okdebug
    use dims,        only: touch_sp, touch_np
    use dims,        only: zero, zoom2D
    use global_data, only: wind_dat

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  region
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_advectyzoom'

    ! --- local ----------------------------------

    real,dimension(:,:,:),pointer               :: bm
    real,dimension(:,:),allocatable             :: bm0,bm1

    integer            :: is,ie,ls,le,n,q
    integer            :: imr,jmr,lmr,tref_,xref_,yref_,zref_,my_parent
    logical            :: y_encountered
    character(len=1)   :: dir
    real               :: sum_old,sum_new

    ! --- begin -------------------------

    call GO_Timer_Start( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    allocate(bm0(0:im(region)+1,0:lm(region)+1))  !mkadj: corrected to lm+1
    allocate(bm1(0:im(region)+1,0:lm(region)+1))

    bm => wind_dat(region)%bm_t

    tref_ = tref(region)/tref(parent(region))
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region);  jmr = jm(region);  lmr = lm(region)

     ! determine the scope for advecty:

    if (region==1) then
       xref_ = 0; zref_ = 0  ! to have is/ie and ls/le properly computed
    endif
    my_parent = parent(region)

    if(my_parent == 0) then   ! always full scope: no parent
       is = 1
       ie = im(region)
       ls = 1
       le = lm(region)
    else
      q=(rstatus(my_parent)-1)/((nsplitsteps/2))      ! find q - the place in the splitorderzoom
                                               ! corresponding to the begining of the last
                                               ! processing order of the parent.
      y_encountered=.false.

      do n=1,n_operators                       ! now track n_operators following
                                               ! steps from q !wp! changed to n_op.. steps

        dir=splitorderzoom(my_parent,q*(nsplitsteps/2)+n)

        select case(dir)
        case('x')
           IF  ((.not.y_encountered).or.(xcyc(region).eq.1)) THEN
              is=1                             ! x-substep is before y =>
              ie=im(region)                    ! full i-scope
           ELSE
              is=xref_+1                       ! y-substep is before x =>
              ie=im(region)-xref_              ! restricted i-scope
           ENDIF
        case('y')
           y_encountered=.true.
        case('z')
           IF (.not.y_encountered) THEN
              ls=1                             ! z-substep is before y =>
              le=lm(region)                    ! full l-scope
           ELSE
              ls=zref_+1                       ! y-substep is before z =>
              le=lm(region)-zref_              ! restricted l-scope
           ENDIF
           IF (zoom2D) THEN
             ls=1; le=lm(region)
           ENDIF
        case('c')
        case('v')
        case('d')
        case('s')
        case default
           print *,'strange value in splitorderzoom(',region,',',   &
                q*(nsplitsteps/2)*tref_+n,'): ',q
           TRACEBACK; status=1; return
        end select

     ENDDO
     endif   !my_parent = 0
     IF ((mod(rstatus(region),(nsplitsteps/2)*tref_)>=n_operators).and.(adv_scheme=='slope')) THEN
        ! IF (1) more than fourth substep and (2) slope scheme THEN
        ! at this substep no coarse fluxes will be applied

        ! for slopes only: zero out velocities via the y-edges of the region
        ! to assure that no fluxes will be applied

        IF(okdebug) print *,'        zeroing out outer mass fluxes (yref)',yref_
        IF(touch_sp(region).ne.1) THEN
          bm0 = bm(:,      yref_,:);   bm(:,      yref_,:) = zero
        ELSE
          IF(okdebug) print *,'       skipping zeroing outer mass flux at SP'
        ENDIF
        IF(touch_np(region).ne.1) THEN
          bm1 = bm(:,jmr-yref_+2,:);   bm(:,jmr-yref_+2,:) = zero
        ELSE
          IF(okdebug) print *,'       skipping zeroing outer mass flux at NP'
        ENDIF

     ENDIF

     call adj_dynamv( region, is,ie,ls,le, status )
     IF_NOTOK_RETURN(status=1)

     IF ((mod(rstatus(region),(nsplitsteps/2)*tref_)>=n_operators).and.(adv_scheme=='slope')) THEN

        ! for slopes only: recreate velocities via the x-edges
        IF(touch_sp(region).ne.1) bm(:,      yref_,:) = bm0
        IF(touch_np(region).ne.1) bm(:,jmr-yref_+2,:) = bm1

     ENDIF


    deallocate(bm0)
    deallocate(bm1)
    nullify(bm)

    call adj_put_yedges( region, status )
    IF_NOTOK_RETURN(status=1)

    call GO_Timer_End( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_advectyzoom


  !-----------------------------------------------------------------------
  !
  !**** dynamv      - south-north tracer transport  v 9.1
  !
  ! programmed by       mh  mpi HH      23-feb-1994
  !
  ! purpose
  ! -------
  ! calculate amount of tracer moved in a south-north advection
  ! substep
  !
  ! interface
  ! ---------
  ! call dynamv
  !
  ! method
  ! ------
  ! slopes scheme
  !
  ! externals
  ! ---------
  ! none
  !
  ! reference
  ! ---------
  ! Russel and Lerner, 1979
  !-----------------------------------------------------------------------
  ! fixed bug in calculating maximum courant number
  !             mh, Thu, Feb 24, 1994 12:49:27
  !     included code for limits of slopes to prevent negative tracer
  !     masses                   mh, 20-jun-1994
  !
  !     zoom version written by mike botchev, march-june 1999
  !     adjoint by MK march 2003
  !-----------------------------------------------------------------------

  !***************************************************************
  !***************************************************************
  !** Part of This routine was generated by the                 **
  !** Tangent linear and Adjoint Model Compiler,    TAMC 4.83   **
  !***************************************************************
  !***************************************************************

  subroutine adj_dynamv( region,is,ie,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,           only: okdebug, im,jm,lm, nregions
    use dims,           only: xref, yref, zref, tref, parent
    use dims,           only: touch_sp, touch_np, zero
    use redgridZoom,    only: nred, jred, imredj, clustsize
    use global_data,    only: mass_dat, wind_dat
    use MeteoData     , only : m_dat
    use chem_param,     only: ntracet
    use ParTools  ,     only: ntracetloc
    use adj_zoom_tools, only: adj_mix_edges

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  region
    integer, intent(in)       ::  is,ie,ls,le
    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_dynamv'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer           :: adrm,adrxm,adrym,adrzm
    real,dimension(:,:,:),  pointer           :: m,bm
    real,dimension(:,:),allocatable           :: mnew
    real,dimension(:,:),allocatable           :: adf,adpf,adfx,adfz
    integer                                   :: i, j, l, n
    integer                                   :: je, js, iee, iss
    integer                                   :: imr,imr2,jmr,lmr
    integer                                   :: tref_,xref_,yref_,zref_
    real                                      :: sfs,sfzs,sfn,sfzn
    real                                      :: beta,mxval
    integer,dimension(3)                      :: mxloc
    integer                                   :: lrg, redfact, ixe, ixs
    real                                      :: summ

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    if(okdebug) print *,'    call adj_dynamv, region=',region,' is,ie,ls,le=',is,ie,ls,le
    if ((region<0).or.(region>nregions)) then
      write (gol,'("illegal number of region !!!")'); call goErr
      TRACEBACK; status=1; return
    end if

    call adj_mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    call adj_uncompress_yedges( region, is,ie,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    m => m_dat(region)%data
    adrm => mass_dat(region)%rm_t
    adrxm => mass_dat(region)%rxm_t
    adrym => mass_dat(region)%rym_t
    adrzm => mass_dat(region)%rzm_t
    bm => wind_dat(region)%bm_t

    ! compute refinement factors with respect to the parent

    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    tref_ = tref(region)/tref(parent(region))
    imr=im(region);jmr=jm(region);lmr=lm(region)

    ! check is,ie,ls,le:
    if((is/=xref_+1)  .and.(is/=1))   stop 'Wrong value for IS in dynamv'
    if((ie/=imr-xref_).and.(ie/=imr)) stop 'Wrong value for IE in dynamv'
    if((ls/=zref_+1)  .and.(ls/=1))   stop 'Wrong value for LS in dynamv'
    if((le/=lmr-zref_).and.(le/=lmr)) stop 'Wrong value for LE in dynamv'

    ! compute js/je -- cells is:ie,js:je:ls:le will be updated
    if (region==1) then
       js = 2            !exclode poles
       je = jmr-1
    else
       js = yref_
       je = jmr-yref_+1
       if(touch_sp(region).eq.1) js = 2      !cmk....new option
       if(touch_np(region).eq.1) je = jmr-1  !
    endif

    ! loop over levels:
    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( touch_np, touch_sp ) &
    !$OMP   shared ( imr, jmr ) &
    !$OMP   shared ( is, ie, js, je, le, ls ) &
    !$OMP   shared ( m ) &
    !$OMP   shared ( bm ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( adrm, adrxm, adrym, adrzm ) &
    !$OMP   private ( i, j, l, n ) &
    !$OMP   private ( mnew ) &
    !$OMP   private ( adf, adpf, adfx, adfz ) &
    !$OMP   private ( beta )

    allocate(mnew(imr,jmr))
    allocate( adf(imr,0:jmr))
    allocate(adpf(imr,0:jmr))
    allocate(adfx(imr,0:jmr))
    allocate(adfz(imr,0:jmr))

    ! init fluxes to zero:
    adf  = 0.0
    adpf = 0.0
    adfx = 0.0
    adfz = 0.0

    !$OMP   DO
    do l=le,ls,-1

      ! calculate new air mass distribution
      if (region==1) then
         mnew(1:imr,1:jmr) = m(1:imr,1:jmr,l)
         m(1:imr,1:jmr,l) = mnew(1:imr,1:jmr)- bm(1:imr,1:jmr,l) + bm(1:imr,2:jmr+1,l)
      else if(touch_sp(region).eq.1) then
         mnew(1:imr,1:je) = m(1:imr,1:je,l)
         m(1:imr,1:je,l) = mnew(1:imr,1:je)- bm(1:imr,1:je,l) + bm(1:imr,2:je+1,l)
      else if(touch_np(region).eq.1) then
         mnew(1:imr,js:jmr) = m(1:imr,js:jmr,l)
         m(1:imr,js:jmr,l) = mnew(1:imr,js:jmr)- bm(1:imr,js:jmr,l) + bm(1:imr,js+1:jmr+1,l)
      else
         mnew(is:ie,js:je) = m(is:ie,js:je,l)
         m(is:ie,js:je,l) = mnew(is:ie,js:je) - bm(is:ie,js  :je  ,l)  &
                                                        + bm(is:ie,js+1:je+1,l)
      endif

      ! loop over tracers:
      do n=1,ntracetloc

       if (region==1.or.touch_np(region).eq.1) then
         do i=1,imr
            adf(i,jmr) = adf(i,jmr) + adrm(i,jmr,l,n)
            adfz(i,jmr) = adfz(i,jmr) + adrzm(i,jmr,l,n)
         enddo
       endif
       if (region==1.or.touch_sp(region).eq.1) then
         do i=1,imr
            adf(i,2) = adf(i,2) - adrm(i,1,l,n)
            adfz(i,2) = adfz(i,2) - adrzm(i,1,l,n)
         enddo
       endif
       do j = je,js,-1
          do i = ie,is,-1
             adfz(i,j+1) = adfz(i,j+1)-adrzm(i,j,l,n)
             adfz(i,j) = adfz(i,j)+adrzm(i,j,l,n)
             adfx(i,j+1) = adfx(i,j+1)-adrxm(i,j,l,n)
             adfx(i,j) = adfx(i,j)+adrxm(i,j,l,n)
             adf(i,j+1) = adf(i,j+1)-adrym(i,j,l,n)*(3.*m(i,j,l)/mnew(i,j))
             adf(i,j) = adf(i,j)-adrym(i,j,l,n)*(3.*m(i,j,l)/mnew(i,j))
             adpf(i,j+1) = adpf(i,j+1)-adrym(i,j,l,n)/mnew(i,j)
             adpf(i,j) = adpf(i,j)+adrym(i,j,l,n)/mnew(i,j)
             adrm(i,j,l,n) = adrm(i,j,l,n)+adrym(i,j,l,n)* &
                      (3.*(bm(i,j,l)+bm(i,j+1,l))/mnew(i,j))
             adrym(i,j,l,n) = adrym(i,j,l,n)*(1-(bm(i,j,l)-bm(i,j+1,l))/mnew(i,j))
             adf(i,j+1) = adf(i,j+1)-adrm(i,j,l,n)
             adf(i,j) = adf(i,j)+adrm(i,j,l,n)
          end do
       end do

       if(region.eq.1.or.touch_np(region).eq.1) then       ! north pole
          do i=ie,is,-1
             if (bm(i,jmr,l)>=zero) then
                    beta = bm(i,jmr,l)/m(i,jmr-1,l)
                    adrzm(i,jmr-1,l,n) = adrzm(i,jmr-1,l,n)+adfz(i,jmr)*beta
                    adfz(i,jmr) = 0.
                    adrxm(i,jmr-1,l,n) = adrxm(i,jmr-1,l,n)+adfx(i,jmr)*beta
                    adfx(i,jmr) = 0.
                    adf(i,jmr) = adf(i,jmr)-3*adpf(i,jmr)*bm(i,jmr,l)
                    adrym(i,jmr-1,l,n) = adrym(i,jmr-1,l,n)+ &
                         adpf(i,jmr)*bm(i,jmr,l)*beta*beta
                    adpf(i,jmr) = 0.
                    adrm(i,jmr-1,l,n) = adrm(i,jmr-1,l,n)+adf(i,jmr)*beta
                    adrym(i,jmr-1,l,n) = adrym(i,jmr-1,l,n)+ &
                         adf(i,jmr)*beta*(1.-beta)
                    adf(i,jmr) = 0.
                 else
                    beta = bm(i,jmr,l)/m(i,jmr,l)
                    adrzm(i,jmr,l,n) = adrzm(i,jmr,l,n)+adfz(i,jmr)*beta
                    adfz(i,jmr) = 0.
                    adfx(i,jmr) = 0.
                    adf(i,jmr) = adf(i,jmr)-3*adpf(i,jmr)*bm(i,jmr,l)
                    adpf(i,jmr) = 0.
                    adrm(i,jmr,l,n) = adrm(i,jmr,l,n)+adf(i,jmr)*beta
                    adf(i,jmr) = 0.
                 endif
          enddo
        else   !zoom region not touching north pole
          j = je+1
          do i=ie,is,-1   !no reduced grid allowed
             if (bm(i,j,l)>=zero) then
                    beta=bm(i,j,l)/m(i,j-1,l)
                    adrzm(i,j-1,l,n) = adrzm(i,j-1,l,n)+adfz(i,j)*beta
                    adfz(i,j) = 0.
                    adrxm(i,j-1,l,n) = adrxm(i,j-1,l,n)+adfx(i,j)*beta
                    adfx(i,j) = 0.
                    adf(i,j) = adf(i,j)-3*adpf(i,j)*bm(i,j,l)
                    adrym(i,j-1,l,n) = adrym(i,j-1,l,n)+ &
                         adpf(i,j)*bm(i,j,l)*beta*beta
                    adpf(i,j) = 0.
                    adrm(i,j-1,l,n) = adrm(i,j-1,l,n)+adf(i,j)*beta
                    adrym(i,j-1,l,n) = adrym(i,j-1,l,n)+ adf(i,j)*beta*(1.-beta)
                    adf(i,j) = 0.
             else
                    beta = bm(i,j,l)/m(i,j,l)
                    adrzm(i,j,l,n) = adrzm(i,j,l,n)+adfz(i,j)*beta
                    adfz(i,j) = 0.
                    adrxm(i,j,l,n) = adrxm(i,j,l,n)+adfx(i,j)*beta
                    adfx(i,j) = 0.
                    adf(i,j) = adf(i,j)-3*adpf(i,j)*bm(i,j,l)
                    adrym(i,j,l,n) = adrym(i,j,l,n) + adpf(i,j)*bm(i,j,l)*beta*beta
                    adpf(i,j) = 0.
                    adrm(i,j,l,n) = adrm(i,j,l,n) + adf(i,j)*beta
                    adrym(i,j,l,n) = adrym(i,j,l,n) - adf(i,j)*beta*(1.+beta)
                    adf(i,j) = 0.
             endif
          enddo
       endif ! compute boundary fluxes north pole...

       if(region.eq.1.or.touch_sp(region).eq.1) then          ! south pole or southern edge:
          do i=ie,is,-1
             adfz(i,1)=zero
             adfx(i,1)=zero
             adpf(i,1)=zero
             adf(i,1)=zero
             if (bm(i,2,l)>=zero) then
                beta=bm(i,2,l)/m(i,1,l)
                adrzm(i,1,l,n) = adrzm(i,1,l,n) + beta*adfz(i,2)
                adfz(i,2) = zero
                adfx(i,2) = zero
                adf(i,2) = adf(i,2) - 3.*bm(i,2,l)*adpf(i,2)
                adpf(i,2)= zero
                adrm(i,1,l,n) = adrm(i,1,l,n) + adf(i,2)*beta
                adf(i,2)=zero
             else
                beta=bm(i,2,l)/m(i,2,l)
                adrzm(i,2,l,n) = adrzm(i,2,l,n)+adfz(i,2)*beta
                adfz(i,2) = 0.
                adrxm(i,2,l,n) = adrxm(i,2,l,n)+adfx(i,2)*beta
                adfx(i,2) = 0.
                adf(i,2) = adf(i,2)-3*adpf(i,2)*bm(i,2,l)
                adrym(i,2,l,n) = adrym(i,2,l,n)+ adpf(i,2)*bm(i,2,l)*beta*beta
                adpf(i,2) = 0.
                adrm(i,2,l,n) = adrm(i,2,l,n)+adf(i,2)*beta
                adrym(i,2,l,n) = adrym(i,2,l,n)-adf(i,2)*beta*(1.+beta)
                adf(i,2) = 0.
             endif
          enddo
        else   !zoom region not touching south pole
          j = js
          do i=ie,is,-1   !no reduced grid allowed
             if (bm(i,j,l)>=zero) then
                beta=bm(i,j,l)/m(i,j-1,l)
                adrzm(i,j-1,l,n) = adrzm(i,j-1,l,n) + adfz(i,j)*beta
                adfz(i,j) = zero
                adrxm(i,j-1,l,n) = adrxm(i,j-1,l,n) + adfx(i,j)*beta
                adfx(i,j) = zero
                adf(i,j) = adf(i,j) -3*adpf(i,j)*bm(i,j,l)
                adrym(i,j-1,l,n) = adrym(i,j-1,l,n) + adpf(i,j)*bm(i,j,l)*beta*beta
                adpf(i,j) = 0.
                adrm(i,j-1,l,n) = adrm(i,j-1,l,n)+adf(i,j)*beta
                adrym(i,j-1,l,n) = adrym(i,j-1,l,n)+adf(i,j)*beta*(1.-beta)
                adf(i,j) = 0.
             else
                beta=bm(i,j,l)/m(i,j,l)
                adrzm(i,j,l,n) = adrzm(i,j,l,n) + adfz(i,j)*beta
                adfz(i,j) = zero
                adrxm(i,j,l,n) = adrxm(i,j,l,n) + adfx(i,j)*beta
                adfx(i,j) = zero
                adf(i,j) = adf(i,j) -3*adpf(i,j)*bm(i,j,l)
                adrym(i,j,l,n) = adrym(i,j,l,n) + adpf(i,j)*bm(i,j,l)*beta*beta
                adpf(i,j) = 0.
                adrm(i,j,l,n) = adrm(i,j,l,n)+adf(i,j)*beta
                adrym(i,j,l,n) = adrym(i,j,l,n)-adf(i,j)*beta*(1.+beta)
                adf(i,j) = 0.
             endif
          enddo
       endif ! compute boundary fluxes south pole...

       do j=je,js+1,-1
          do i=ie,is,-1
             if (bm(i,j,l)>=zero) then
                beta=bm(i,j,l)/m(i,j-1,l)
                adrzm(i,j-1,l,n) = adrzm(i,j-1,l,n) + adfz(i,j)*beta
                adfz(i,j) = zero
                adrxm(i,j-1,l,n) = adrxm(i,j-1,l,n) + adfx(i,j)*beta
                adfx(i,j) = zero
                adf(i,j) = adf(i,j) -3*adpf(i,j)*bm(i,j,l)
                adrym(i,j-1,l,n) = adrym(i,j-1,l,n) + adpf(i,j)*bm(i,j,l)*beta*beta
                adpf(i,j) = 0.
                adrm(i,j-1,l,n) = adrm(i,j-1,l,n)+adf(i,j)*beta
                adrym(i,j-1,l,n) = adrym(i,j-1,l,n)+adf(i,j)*beta*(1.-beta)
                adf(i,j) = 0.
             else
                beta=bm(i,j,l)/m(i,j,l)
                adrzm(i,j,l,n) = adrzm(i,j,l,n) + adfz(i,j)*beta
                adfz(i,j) = zero
                adrxm(i,j,l,n) = adrxm(i,j,l,n) + adfx(i,j)*beta
                adfx(i,j) = zero
                adf(i,j) = adf(i,j) -3*adpf(i,j)*bm(i,j,l)
                adrym(i,j,l,n) = adrym(i,j,l,n) + adpf(i,j)*bm(i,j,l)*beta*beta
                adpf(i,j) = 0.
                adrm(i,j,l,n) = adrm(i,j,l,n)+adf(i,j)*beta
                adrym(i,j,l,n) = adrym(i,j,l,n)-adf(i,j)*beta*(1.+beta)
                adf(i,j) = 0.
             endif
          enddo
       enddo

     enddo  ! n

    enddo  ! l
    !$OMP   END DO

    deallocate(mnew)
    deallocate( adf)
    deallocate(adpf)
    deallocate(adfx)
    deallocate(adfz)

    !$OMP END PARALLEL

    nullify(m)
    nullify(adrm)
    nullify(adrxm)
    nullify(adrym)
    nullify(adrzm)
    nullify(bm)

    call adj_compress_yedges( region, is,ie,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    call GO_Timer_End( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)
    ! ok
    status = 0

  end subroutine adj_dynamv


  ! ***


  ! this is for slope only: condense data at the y-edges of the zoom region
  ! to allow for uniform work in dynamv
  ! written by mike botchev, march-june 1999
  ! modified by MK, dec 2002
  ! ADJOINT implementation, mk, frb 2003
  !    forward mode:
  !      xc(xref_) = (xc(1) + xc(2) + ...)
  !      xc(xref_-1) = xc(0)
  !      xc(0:xref_-2) = 0.0
  !      xc(im-xref_+1) = (xc(im-xref_+1) + xc(im-xref_+2) + ...)
  !      xc(im-xref_+2) = xc(im+1)
  !      xc(im-xref_+3:im+1) = 0.0
  !
  !    adjoint modeL
  !      xc(0) = xc(xref_-1)
  !      xc(1:xref_) = xc(xref_)
  !      xc(im+1) = xc(im-xref_+2)
  !      xc(im-xref_+1:im) = xc(im-xref_+1)
  !

  subroutine adj_compress_yedges( region, is,ie,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only: jm,yref,okdebug, touch_sp, touch_np, parent
    use global_data, only: mass_dat
    use MeteoData  , only : m_dat
    use chem_param,  only: ntracet
    use ParTools  ,     only: ntracetloc

    ! --- in/out ---------------------------------

    integer                   ::  region
    integer, intent(in)       ::  is,ie,ls,le
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_compress_yedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer           :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer           :: m
    integer                                   :: jmr,i,l,n,yref_,j
    real                                      :: medge,rmedge, rxmedge, rymedge, rzmedge

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_compress, status )
    IF_NOTOK_RETURN(status=1)

    yref_ = yref(region)/yref(parent(region))
    jmr = jm(region)
    if ((yref_==1).or.(region==1)) then
       if(okdebug) print *,'       adj_compress_yedges: no refinement, nothnig to do'
       call GO_Timer_End( itim_compress, status )
       IF_NOTOK_RETURN(status=1)
       return
    endif
    if(okdebug) then
       print *,'       call adj_compress_yedges, region=',region,' jmr=',jmr
       print *,'       cells to be updated: ',yref_-1,yref_,jmr-yref_+1,jmr-yref_+2
       if(touch_sp(region).eq.1) print *,'     SP will be skipped here'
       if(touch_np(region).eq.1) print *,'     NP will be skipped here'
    endif

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( jmr ) &
    !$OMP   shared ( yref_ ) &
    !$OMP   shared ( is, ie, ls, le ) &
    !$OMP   shared ( touch_sp, touch_np ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m, rm, rxm, rym, rzm ) &
    !$OMP   private ( i, l, n ) &
    !$OMP   private ( medge, rmedge )
    !$OMP   DO
    do l=ls,le
       do i=is,ie

          !for m, calculate the masses as BEFORE compress_yedges was called
          if(touch_sp(region).ne.1) then
             m (i, 0, l) = m (i,yref_-1,l)
             medge = m (i,yref_,  l)/yref_
             do j=1,yref_
                m(i,j,l) = medge
             enddo
          endif
          if(touch_np(region).ne.1) then
             m (i,jmr+1,l) =  m (i,jmr-yref_+2,l)
             medge =    m (i,jmr-yref_+1,l)/yref_
             do j=1,yref_
                m (i,jmr-j+1,l) = medge
             enddo
          endif

          !for rm, rxm, rym ,rzm: adjoint operation (difference is the division by yref_)
          do n=1,ntracetloc

            if(touch_sp(region) /= 1) then
               rm (i, 0,l,n) = rm (i, yref_-1, l,n)
               rm(i,1:yref_,l,n) = rm(i,yref_,l,n)
               rxm(i, 0,l,n) = rxm(i, 0,l,n)  + rxm(i, yref_-1, l,n)
               rxm(i, yref_-1, l,n) = 0.0
               rym(i, 0,l,n) = rym(i, 0,l,n)  + rym(i, yref_-1, l,n)
               rym(i, yref_-1, l,n) = 0.0
               rzm(i, 0,l,n) = rzm(i, 0,l,n)  + rzm(i, yref_-1, l,n)
               rzm(i, yref_-1, l,n) = 0.0
            endif
            if(touch_np(region) /= 1) then
               rm (i, jmr+1,l,n) = rm (i, jmr-yref_+2, l,n)
               rmedge  = rm (i,jmr-yref_+1,l,n)
               do j=1,yref_
                 rm (i,jmr+1-j,l,n) = rmedge
               enddo
               rxm(i, jmr+1,l,n) = rxm(i, jmr+1,l,n) + rxm(i, jmr-yref_+2, l,n)
               rxm(i, jmr-yref_+2, l,n) = 0.0
               rym(i, jmr+1,l,n) = rym(i, jmr+1,l,n) + rym(i, jmr-yref_+2, l,n)
               rym(i, jmr-yref_+2, l,n) = 0.0
               rzm(i, jmr+1,l,n) = rzm(i, jmr+1,l,n) + rzm(i, jmr-yref_+2, l,n)
               rzm(i, jmr-yref_+2, l,n) = 0.0
            endif
          enddo ! n

       enddo  ! i
    enddo   ! l
    !$OMP   END DO
    !$OMP END PARALLEL

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

    call GO_Timer_End( itim_compress, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_compress_yedges


  ! ***


  ! distribute data at the y-edges of the zoom region
  ! over the whole interface cell
  ! written by mike botchev, march-june 1999
  ! ADJOINT implementation, MK, feb 2003.
  !
  !  forward:
  !     xc(1:xref_) = (xc(xref_) )/xref_     !NOTE: the cells 0, and imr+1 are not used any more after this step,
  !                                           because the flux between cells cells xref_-1 and xref will be zero'd
  !     xc(im-xref_+1:im) = (xc(im-xref_+1) )/xref_
  !
  !  adjoint mode:
  !     xc(xref_) = (xc(1) + ..+ xc(xref_))/xref_
  !     xc(1:xref_-1) = 0.0
  !     xc(im+1-xref_) = (xc(im) + ..+ xc(im+1-xref_))/xref_
  !     xc(im+2-xref_......im) = 0.0
  !

  subroutine adj_uncompress_yedges( region, is,ie,ls,le, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims       , only : jm, okdebug, yref, parent, touch_sp, touch_np
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use chem_param , only : ntracet
    use ParTools   , only : ntracetloc

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  region
    integer, intent(in)       ::  is,ie,ls,le
    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_uncompress_yedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer           :: rm,rxm,rym,rzm
    real,dimension(:,:,:),  pointer           :: m
    integer                                   :: jmr,i,l,n,yref_
    real                                      :: m_edge,rm_edge

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_uncompress, status )
    IF_NOTOK_RETURN(status=1)

    yref_ = yref(region)/yref(parent(region))
    jmr = jm(region)
    if ((yref_==1).or.(region==1)) then
      if(okdebug) print *,'       uncompress_yedges: no refinement, nothnig to do'
      call GO_Timer_End( itim_uncompress, status )
      IF_NOTOK_RETURN(status=1)
      return
    end if
    if(okdebug) then
       print *,'       call adj_uncompress_yedges, region=',region,' jmr=',jmr
       print *,'       cells to be updated: ',1,':',yref_,' ',jmr-yref_+1,':',jmr
       if(touch_sp(region).eq.1) print *,'     SP will be skipped here'
       if(touch_np(region).eq.1) print *,'     NP will be skipped here'
    endif

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( jmr ) &
    !$OMP   shared ( yref_ ) &
    !$OMP   shared ( is, ie, ls, le ) &
    !$OMP   shared ( touch_sp, touch_np ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m, rm, rxm, rym, rzm ) &
    !$OMP   private ( i, l, n )
    !$OMP   DO
    do l=ls,le
       do i=is,ie

          ! for m, calculate the masses BEFORE uncompress edges
          if(touch_sp(region) /= 1) then
             m(i,yref_,l) = sum(m(i,1:yref_,l))
             m(i,yref_-1,l) = m(i,0,l)
             m(i,0:yref_-2,l) = 0.0   ! will be overwritten by adj_compress anyhow
          endif
          if(touch_np(region) /= 1) then
             m(i,jmr-yref_+1,l) = sum(m(i,jmr-yref_+1:jmr,l))
             m(i,jmr-yref_+2,l) = m(i,jmr+1,l)
             m(i,jmr-yref_+3:jmr+1,l) = 0.0
          endif

          do n=1,ntracetloc
             if(touch_sp(region) /= 1) then     !rewrite......not correct!

               rm(i,yref_,l,n) = sum(rm(i,1:yref_,l,n))/yref_
               rm(i,1:yref_-1,l,n) = 0.0

               rxm(i,yref_,l,n) = sum(rxm(i,1:yref_,l,n))/yref_
               rxm(i,1:yref_-1,l,n) = 0.0

               rym(i,yref_,l,n) = sum(rym(i,1:yref_,l,n))/yref_
               rym(i,1:yref_-1,l,n) = 0.0

               rzm(i,yref_,l,n) = sum(rzm(i,1:yref_,l,n))/yref_
               rzm(i,1:yref_-1,l,n) = 0.0

             endif

             if(touch_np(region) /= 1) then

                rm(i,jmr-yref_+1,l,n) = sum(rm(i,jmr-yref_+1:jmr,l,n))/yref_
                rm(i,jmr-yref_+2:jmr,l,n) = 0.0

                rxm(i,jmr-yref_+1,l,n) = sum(rxm(i,jmr-yref_+1:jmr,l,n))/yref_
                rxm(i,jmr-yref_+2:jmr,l,n) = 0.0

                rym(i,jmr-yref_+1,l,n) = sum(rym(i,jmr-yref_+1:jmr,l,n))/yref_
                rym(i,jmr-yref_+2:jmr,l,n) = 0.0

                rzm(i,jmr-yref_+1,l,n) = sum(rzm(i,jmr-yref_+1:jmr,l,n))/yref_
                rzm(i,jmr-yref_+2:jmr,l,n) = 0.0

             endif
          enddo ! n

       enddo  ! i
    enddo  ! l
    !$OMP   END DO
    !$OMP END PARALLEL

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

    call GO_Timer_End( itim_uncompress, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_uncompress_yedges


  ! ***


  ! passes values along the x-boundaries to all the children
  ! written by mike botchev, march-june 1999
  ! structure implemented by MK, dec 2002
  ! adjoint routine, MK feb. 2003
  !
  ! in forward mode, in principle:
  !      xc(0) = xp(ibeg-1)
  !      xc(im+1) = xp(iend+1)
  !      xc(1:xref_) = xp(ibeg)/xref_
  !      xc(im-xref_+1: im) = xp(iend)/xref_
  !           or
  !
  ! in adjoint mode:
  !    xp(ibeg-1) = xp(ibeg-1) + xc(0)
  !    xp(ibeg)   = xp(ibeg)   + (xc(1) + xc(2) +...)/xref_
  !    xc(0:xref_) = 0.0
  !    xp(iend)   = xp(iend)   + (xc(im) + xc(im-1) +...)/xref_
  !    xp(iend+1) = xp(iend+1) + xc(im+1)
  !    xc(imr-xref_+1:imr+1) = 0.0
  !

  subroutine adj_put_yedges( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only: splitorderzoom, nsplitsteps, n_operators
    use dims,        only: im, jm, lm
    use dims,        only: touch_sp, touch_np, children
    use dims,        only: xref, yref, zref
    use dims,        only: jbeg, jend, ibeg, lbeg
    use dims,        only: rstatus => status
    use dims,        only: okdebug
    use global_data, only: mass_dat, wind_dat
    use MeteoData  , only : m_dat
    use chem_param,  only: ra, ntracet
    use ParTools  ,     only: ntracetloc

    ! --- in/out ---------------------------------

    integer,intent(in)        ::  region
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_put_yedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer              :: rm,rxm,rym,rzm
    real,dimension(:,:,:,:),pointer              :: rmc,rxmc,rymc,rzmc
    real,dimension(:,:,:),pointer                :: m,mc
    real,dimension(:,:,:),pointer                :: am,bm,cm
    integer            :: child,ichild
    integer            :: i, j, l, n
    integer            :: ip, jps, jpn, lp
    integer            :: imc, jmc, lmc
    integer            :: xref_,yref_,zref_
    real               :: xzref,xyzref
    real               :: mpn,mps
    integer            :: nzne, nzne_v,nznem,imr,q
    character(len=n_operators) :: stencil
    logical                    :: xyz

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_put_yedges, status )
    IF_NOTOK_RETURN(status=1)

    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    m => m_dat(region)%data
    am => wind_dat(region)%am_t
    bm => wind_dat(region)%bm_t
    cm => wind_dat(region)%cm_t

    q=(rstatus(region)-1)/((nsplitsteps/2))      ! find q - the place in the splitorderzoom of the parent
    do i=1,n_operators
       stencil(i:i)=splitorderzoom(region,q*(nsplitsteps/2)+i)
    enddo
    ! find out if parent has order xyz, or zyx...
    if( scan(stencil,'x') >  scan(stencil,'z')) then
       xyz = .false.
    else
       xyz = .true.
    endif


    ichild = children(region,0)    ! start at the back
    do while(ichild>0)             ! as long as children
      child = children(region,ichild)   ! get the number of the child region
      xref_ = xref(child)/xref(region)
      yref_ = yref(child)/yref(region)
      zref_ = zref(child)/zref(region)

      if(touch_sp(child).eq.1) then
         if (okdebug)  &
              print *,'     put_yedges: child ',child,' sp  skipped'
      endif
      if(touch_np(child).eq.1) then
         if (okdebug)  &
              print *,'     put_yedges: child ',child,' np  skipped'
      endif

      imc = im(child)
      jmc = jm(child)
      lmc = lm(child)

      xzref = 1./(xref_*zref_)
      xyzref = 1./(xref_*yref_*zref_)

      jpn = jend(child)
      jps = jbeg(child)

      ! put_yedges puts masses of the parent in the child.
      ! In the case of a XYZ sequence, these masses did not see a full X advection
      ! In the case of the parent sequence ZYX, these masses have seen a full Z advection.
      ! In the adjoint, we have to advect these masses of the child back in time.
      !
      if(okdebug) print *, ' call adj_put_yedges:', region,child,xyz

      mc => m_dat(child)%data
      rmc => mass_dat(child)%rm_t
      rxmc => mass_dat(child)%rxm_t
      rymc => mass_dat(child)%rym_t
      rzmc => mass_dat(child)%rzm_t

      ! loop over layers:
      !$OMP PARALLEL &
      !$OMP   default ( none ) &
      !$OMP   shared ( imc, jmc, lmc ) &
      !$OMP   shared ( child ) &
      !$OMP   shared ( ibeg, jbeg, jend, lbeg ) &
      !$OMP   shared ( jpn, jps ) &
      !$OMP   shared ( xyz ) &
      !$OMP   shared ( touch_np, touch_sp ) &
      !$OMP   shared ( xref_, yref_, zref_ ) &
      !$OMP   shared ( xzref, xyzref ) &
      !$OMP   shared ( m, am, cm ) &
      !$OMP   shared ( rm ) &
      !$OMP   shared ( rxm, rym, rzm ) &
      !$OMP   shared ( mc, rmc ) &
      !$OMP   shared ( rxmc, rymc, rzmc ) &
      !$OMP   private ( i, j, l ) &
      !$OMP   private ( ip, lp ) &
      !$OMP   private ( mpn, mps )
      !$OMP   DO
      do l=1,lmc

        ! parent layer:
        lp = lbeg(child) + (l-1)/zref_

        do i=1,imc
           ip = ibeg(child) + (i-1)/xref_
           if(xyz) then   ! xadvection backwards in interface cells
              if (touch_sp(child)==0) then   !
                 mps = m(ip,jps,l) + am(ip,jps,l) - am(ip-1,jps,l)
                 mc(i,1:yref_,l) = mps/(xref_*yref_)
              endif
              if (touch_np(child)==0) then   !
                 mpn = m(ip,jpn,l) + am(ip,jpn,l) - am(ip-1,jpn,l)
                 mc(i,jmc-yref_+1:jmc,l) = mpn/(xref_*yref_)
              endif
           else          ! z advection in backwards direction  (ZYX)
              if (touch_sp(child)==0) then   !
                 mps = m(ip,jps,l) + cm(ip,jps,l) - cm(ip,jps,l-1)
                 mc(i,1:yref_,l) = mps/(yref_*xref_)
              endif
              if (touch_np(child)==0) then   !
                 mpn = m(ip,jpn,l) + cm(ip,jpn,l) - cm(ip,jpn,l-1)
                 mc(i,jmc-yref_+1:jmc,l) = mpn/(yref_*xref_)
              endif
           endif
        enddo  ! i

        do i=1,imc
           ip = ibeg(child) + (i-1)/xref_
           if(touch_sp(child)==0) then
              rm(ip,jps-1,lp,:) = rm(ip,jps-1,lp,:) + rmc(i,0,l,:)*xzref
              rmc(i,0,l,:) = 0.0
           endif
           if(touch_np(child)==0) then
              rm(ip,jpn+1,lp,:) = rm(ip,jpn+1,lp,:) + rmc(i,jmc+1,l,:)*xzref
              rmc(i,jmc+1,l,:) = 0.0
           endif

           do j=1,yref_
              if(touch_sp(child)==0) then
                 rm(ip,jps,lp,:) = rm(ip,jps,lp,:) + rmc(i,j,l,:)*xyzref
                 rmc(i,j,l,:) = 0.0
              endif
              if(touch_np(child)==0) then
                 rm(ip,jpn,lp,:) = rm(ip,jpn,lp,:) + rmc(i,jmc+1-j,l,:)*xyzref
                 rmc(i,jmc+1-j,l,:) = 0.0
              endif
           enddo
        enddo

        do i=1,imc
           ip = ibeg(child) + (i-1)/xref_
           if(touch_sp(child)==0) then
              rxm(ip,jps-1,lp,:) = rxm(ip,jps-1,lp,:) + rxmc(i,0,l,:)*xzref
              rxmc(i,0,l,:) = 0.0
           endif
           if(touch_np(child)==0) then
              rxm(ip,jpn+1,lp,:) = rxm(ip,jpn+1,lp,:) + rxmc(i,jmc+1,l,:)*xzref
              rxmc(i,jmc+1,l,:) = 0.0
           endif

           do j=1,yref_
              if(touch_sp(child)==0) then
                 rxm(ip,jps,lp,:) = rxm(ip,jps,lp,:) + rxmc(i,j,l,:)*xzref
                 rxmc(i,j,l,:) = 0.0
              endif
              if(touch_np(child)==0) then
                 rxm(ip,jpn,lp,:) = rxm(ip,jpn,lp,:) + rxmc(i,jmc+1-j,l,:)*xzref
                 rxmc(i,jmc+1-j,l,:) = 0.0
              endif
           enddo
        enddo

        do i=1,imc
            ip = ibeg(child) + (i-1)/xref_
            if(touch_sp(child)==0) then
               rym(ip,jps-1,lp,:) = rym(ip,jps-1,lp,:) + rymc(i,0,l,:)*xzref
               rymc(i,0,l,:) = 0.0
            endif
            if(touch_np(child)==0) then
               rym(ip,jpn+1,lp,:) = rym(ip,jpn+1,lp,:) + rymc(i,jmc+1,l,:)*xzref
               rymc(i,jmc+1,l,:) = 0.0
            endif

            do j=1,yref_
               if(touch_sp(child)==0) then
                  rym(ip,jps,lp,:) = rym(ip,jps,lp,:) + rymc(i,j,l,:)*xzref
                  rymc(i,j,l,:) = 0.0
               endif
               if(touch_np(child)==0) then
                  rym(ip,jpn,lp,:) = rym(ip,jpn,lp,:) + rymc(i,jmc+1-j,l,:)*xzref
                  rymc(i,jmc+1-j,l,:) = 0.0
               endif
            enddo
        enddo  ! i

         do i=1,imc
            ip = ibeg(child) + (i-1)/xref_
            if(touch_sp(child)==0) then
               rzm(ip,jps-1,lp,:) = rzm(ip,jps-1,lp,:) + rzmc(i,0,l,:)*xzref
               rzmc(i,0,l,:) = 0.0
            endif
            if(touch_np(child)==0) then
               rzm(ip,jpn+1,lp,:) = rzm(ip,jpn+1,lp,:) + rzmc(i,jmc+1,l,:)*xzref
               rzmc(i,jmc+1,l,:) = 0.0
            endif

            do j=1,yref_
               if(touch_sp(child)==0) then
                  rzm(ip,jps,lp,:) = rzm(ip,jps,lp,:) + rzmc(i,j,l,:)*xzref
                  rzmc(i,j,l,:) = 0.0
               endif
               if(touch_np(child)==0) then
                  rzm(ip,jpn,lp,:) = rzm(ip,jpn,lp,:) + rzmc(i,jmc+1-j,l,:)*xzref
                  rzmc(i,jmc+1-j,l,:) = 0.0
               endif
           enddo  ! j
        enddo  ! i

      enddo  ! layers
      !$OMP   END DO
      !$OMP END PARALLEL

      ! clear:
      nullify(mc)
      nullify(rmc)
      nullify(rxmc)
      nullify(rymc)
      nullify(rzmc)

      ! count down
      ichild = ichild - 1
    end do  ! regions

    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(m)
    nullify(am)
    nullify(bm)
    nullify(cm)

    call GO_Timer_End( itim_put_yedges, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_put_yedges


end module adj_advecty
