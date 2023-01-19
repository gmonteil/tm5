!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_advectx

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public :: adj_AdvectX_Init, adj_AdvectX_Done
  public :: adj_advectxzoom


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'adj_AdvectX'


  ! --- local ------------------------------------

  integer    ::  itim_zoom
  integer    ::  itim_work
  integer    ::  itim_dynam
  integer    ::  itim_compress
  integer    ::  itim_uncompress
  integer    ::  itim_put_xedges


contains


  ! ====================================================================


  subroutine adj_AdvectX_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_AdvectX_Init'

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_zoom, 'adj_advectx zoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_work, 'adj_advectx work', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynam, 'adj_advectx dynam', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_compress, 'adj_advectx compress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_uncompress, 'adj_advectx uncompress', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_put_xedges, 'adj_advectx put xedges', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_AdvectX_Init


  ! ***


  subroutine adj_AdvectX_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_AdvectX_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine adj_AdvectX_Done


  ! ***


  ! set parameters for advectx
  ! written by patrick berkvens and mike botchev, march-june 1999
  ! updated and modified by MK, dec 2002
  ! for adjoint version....feb 2003 MK

  subroutine adj_advectxzoom( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,          only: im,jm,lm
    use dims,          only: tref,xref,yref
    use dims,          only: okdebug
    use dims,          only: zero
    use dims,          only: rstatus => status
    use dims,          only: parent, nsplitsteps
    use dims,          only: zref, n_operators, splitorderzoom
    use dims,          only: touch_sp, touch_np, zoom2d, adv_scheme
    use global_data  , only: mass_dat, wind_dat
    use toolbox,       only: escape_tm

    ! --- in/out ----------------------------------

    integer, intent(in)             ::  region
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_advectxzoom'

    ! --- local ----------------------------------

    real,dimension(:,:,:),  pointer     :: am
    integer                             :: js,je,ls,le,n,q
    integer                             :: imr,jmr,lmr,tref_,xref_,yref_,zref_,my_parent
    real,dimension(:,:),allocatable     :: am0,am1  ! to store mass_flux
    logical                             :: x_encountered
    character(len=1)                    :: dir

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    !cmk changed to lm+1 below
    allocate(am0(0:jm(region)+1, 0:lm(region)+1))
    allocate(am1(0:jm(region)+1, 0:lm(region)+1))

    am => wind_dat(region)%am_t

    tref_ = tref(region)/tref(parent(region))
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region);  jmr = jm(region)

    ! determine the scope for advectx:

    if (region==1) then
        yref_ = 0; zref_ = 0  ! to have js/je and ls/le properly computed
    endif
    my_parent = parent(region)


    if(my_parent == 0) then   ! always full scope: no parent
       js = 1
       je = jm(region)
       ls = 1
       le = lm(region)
    else
      q=(rstatus(my_parent)-1)/((nsplitsteps/2))      ! find q - the place in the splitorderzoom
                                               ! corresponding to the begining of the last
                                               ! processing order of the parent.
      x_encountered=.false.

      do n=1,n_operators                       ! now track n_operators following
                                               ! steps from q !wp! changed to n_op.. steps

       dir=splitorderzoom(my_parent,q*(nsplitsteps/2)+n)
       select case(dir)
       case('x')
          x_encountered=.true.
       case('y')
          IF (.not.x_encountered) THEN
            js=1                             ! y-substep is before x =>
            je=jm(region)                    ! full j-scope
          ELSE
            IF(touch_sp(region).eq.1) THEN
              js=1                             !CMK special case...touching SP
            ELSE
              js=yref_+1                       ! x-substep is before y =>
            ENDIF
            IF(touch_np(region).eq.1) THEN
              je=jm(region)                    ! CMK special case: touching NP
            ELSE
              je=jm(region)-yref_              ! restricted j-scope
            ENDIF
         ENDIF
       case('z')
         IF  (.not.x_encountered) THEN
           ls=1                             ! z-substep is before x =>
           le=lm(region)                    ! full l-scope
         ELSE
           ls=zref_+1                       ! x-substep is before z =>
           le=lm(region)-zref_              ! restricted l-scope
         ENDIF
         IF (zoom2D) THEN
             ls=1; le=lm(region)
         ENDIF
       case ('c')
       case ('v')
       case ('d')
       case ('s')
       case default
          print *,'strange value in splitorderzoom(',region,',',   &
               q*(nsplitsteps/2)*tref_+n,'): ',q
          call escape_tm( ' error in x-advection ')
       end select

     ENDDO
    endif   ! parent = 0
    ! adjoint version remains identical, since we run from back --> front
    IF ((mod(rstatus(region)-1,(nsplitsteps/2)*tref_)>=n_operators).and.       &
        (adv_scheme=='slope').and.(im(region)/xref(region)<im(1))) THEN
        ! IF (1) more than n_operators substep, (2) slope scheme, and
        ! (3) the region is not [0;360] degrees wide THEN
        ! at this substep no coarse fluxes will be applied

        ! for slopes only: zero out velocities via the x-edges of the region
        ! to assure that no fluxes will be applied

        IF(okdebug) print *,'    zeroing out outer mass fluxes'
        am0 = am(    xref_-1,:,:);   am(    xref_-1,:,:) = zero
        am1 = am(imr-xref_+1,:,:);   am(imr-xref_+1,:,:) = zero

    ENDIF

    if(okdebug) print *, 'call adj_advectx_slopes with scope:',js,je,ls,le
    call adj_advectx_slopes( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)
    IF ((mod(rstatus(region)-1,(nsplitsteps/2)*tref_)>=n_operators).and.   &   ! note the -1
        (adv_scheme=='slope').and.(im(region)/xref(region)<im(1))) THEN
        ! for slopes only: recreate velocities via the x-edges
        am(    xref_-1,:,:) = am0
        am(imr-xref_+1,:,:) = am1
    ENDIF

    ! the adjoint of put_xedges!
    call adj_put_xedges( region, status )
    IF_NOTOK_RETURN(status=1)

    nullify(am)
    deallocate(am0)
    deallocate(am1)

    call GO_Timer_End( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_advectxzoom


  ! ***


  ! makes reduced grid pre-/postprocessing and switches between dynamu and dynamu1
  ! written by mike botchev, march-june 1999
  ! adjoint version: assumes use of dynamu

  subroutine adj_advectx_slopes( region, js,je,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,            only: adv_scheme, im, jm, lm
    use redgridZoom,     only: grid_reduced, nred, uni2red_mf
    use adj_redgridZoom, only: adj_uni2red, adj_red2uni_em, adj_red2uni
    use global_data,     only: mass_dat, wind_dat
    use MeteoData      , only : m_dat
    use chem_param,      only: ntracet
    use ParTools  ,      only: ntracetloc
    use toolbox,         only: escape_tm

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  region
    integer, intent(in)       ::  js,je,ls,le
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_advectx_slopes'

    ! --- local ----------------------------------

    real,dimension(:,:,:)  ,pointer   :: m,am
    real,dimension(:,:,:),allocatable :: m_uni, am_uni    ! used for reduced gird...
    real                              :: mxval
    integer,dimension(3)              :: mxloc
    integer                           :: imr,jmr,lmr,j

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_work, status )
    IF_NOTOK_RETURN(status=1)

    if (adv_scheme/='slope') call escape_tm('Wrong advection scheme: adv_scheme/=slope')

    ! cmk increased dimension to lm+1 fro am_uni
    allocate(m_uni(-1:im(region)+2,-1:jm(region)+2,   lm(region)))
    allocate(am_uni(0:im(region)+1, 0:jm(region)+1, 0:lm(region)+1))

    am => wind_dat(region)%am_t
    m  => m_dat(region)%data

    imr = im(region) ; jmr = jm(region) ; lmr = lm(region)
    ! transform to the reduced grid:

    if (grid_reduced.and.(nred(region).ne.0)) then  !check for reduced grid in region

       ! save non-reduced m and am in m_uni and am_uni:
       m_uni = m; am_uni = am

       ! reduce m,rm,rxm,rym,rzm:
       !call adj_red2uni_em(region)
       call adj_red2uni(region)    ! contains uni2red_mf

       ! reduce am:
       !call uni2red_mf(region)
    endif

    call adj_dynamu( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    ! transform from the reduced grid:
    if (grid_reduced.and.nred(region)/=0) then
       ! advection on uniform grid:
       ! redistribute rm,rxm,rym,rzm
       call adj_uni2red(region)
       ! note: now----> m comes back correct!
       ! JFM: recover m from back advection of m_uni
       m(1:imr,1:jmr,1:lmr)=m_uni(1:imr,1:jmr,1:lmr) - am_uni(0:imr-1,1:jmr,1:lmr) &
                                                      +am_uni(1:imr,  1:jmr,1:lmr)
       ! recreate am:
       am = am_uni
    endif

    nullify(am)
    nullify(m)

    deallocate(am_uni)
    deallocate(m_uni)

    call GO_Timer_End( itim_work, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_advectx_slopes


  ! ***


  !***************************************************************
  !***************************************************************
  !** This routine was generated by the                         **
  !** Tangent linear and Adjoint Model Compiler,    TAMC 4.83   **
  !** Adapted for TM5 by Maarten Krol, March 2003               **
  !** Included iteration (MK /04/2005)                          **
  !***************************************************************
  !***************************************************************

  subroutine adj_dynamu( region, js,je,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,           only : okdebug, nregions, parent
    use dims,           only : im, jm, lm, xref, yref, zref, tref
    use dims,           only : zero, one, xi, nxi, nloop_max, limits, xcyc
    use redgridZoom,    only : grid_reduced, nred, imred
    use global_data,    only : mass_dat, wind_dat
    use MeteoData     , only : m_dat
    use AdvectM_CFL   , only : advectx_get_nloop
    use adj_zoom_tools, only : adj_mix_edges
    use toolbox,        only : escape_tm
    use chem_param,     only : ntracet
    use ParTools,       only: ntracetloc

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  region
    integer, intent(in)       ::  js,je,ls,le
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_dynamu'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer    :: adrm,adrxm,adrym,adrzm
    real,dimension(:,:,:),  pointer    :: m,am
    real,dimension(:),allocatable      :: adf,adpf,adfy,adfz
    real,dimension(:),allocatable      :: mnew

    integer                            :: nloop, iloop

    integer                            ::  i, j, l, n
    integer                            ::  ie, is
    integer                            ::  iie
    integer                            ::  imr, jmr, lmr
    integer                            ::  xref_, yref_, zref_
    integer                            ::  ip1
    real                               ::  alpha

    logical            ::  special_grid
    integer            ::  nfail

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    !----------------------------------------------
    ! RESET LOCAL ADJOINT VARIABLES
    !----------------------------------------------

    if(okdebug) print *,'    call adj_dynamu, region=',region,' js,je,ls,le=',js,je,ls,le

    imr = im(region)
    jmr = jm(region)
    lmr = lm(region)

    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))

    ! check js,je,ls,le:

    if((js/=yref_+1)  .and.(js/=1))   stop 'Wrong value for JS in dynamu'
    if((je/=jmr-yref_).and.(je/=jmr)) stop 'Wrong value for JE in dynamu'
    if((ls/=zref_+1)  .and.(ls/=1))   stop 'Wrong value for LS in dynamu'
    if((le/=lmr-zref_).and.(le/=lmr)) stop 'Wrong value for LE in dynamu'

    ! compute is/ie -- cells is:ie,js:je:ls:le will be updated
    is = xref_; ie = imr-xref_+1   ! default scope (for zoom region, global 1:im)
    if ( xcyc(region) == 1 ) then ! periodic boundary condition
      is = 1; ie = imr
    endif

    ! reduced grid ?
    special_grid = grid_reduced .and. (nred(region) /= 0)

    call adj_mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    call adj_uncompress_xedges( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    am => wind_dat(region)%am_t
    m => m_dat(region)%data
    adrm => mass_dat(region)%rm_t
    adrxm => mass_dat(region)%rxm_t
    adrym => mass_dat(region)%rym_t
    adrzm => mass_dat(region)%rzm_t

    !----------------------------------------------
    ! ROUTINE BODY
    !----------------------------------------------

    ! JFM: adjoint of periodic boundary conditions (applied first!)
    if ( xcyc(region) == 1 ) then
       do n=1,ntracetloc
          do l=ls,le
             do j=js,je
                if ( special_grid ) then
                  iie = imred(j,region) ! red grid
                else
                  iie = im(region)
                end if
                adrm (iie,j,l,n) = adrm (iie,j,l,n) + adrm (0,j,l,n)
                adrm (  0,j,l,n) = 0.0
                adrxm(iie,j,l,n) = adrxm(iie,j,l,n) + adrxm(0,j,l,n)
                adrxm(  0,j,l,n) = 0.0
                adrym(iie,j,l,n) = adrym(iie,j,l,n) + adrym(0,j,l,n)
                adrym(  0,j,l,n) = 0.0
                adrzm(iie,j,l,n) = adrzm(iie,j,l,n) + adrzm(0,j,l,n)
                adrzm(  0,j,l,n) = 0.0
                adrm (  1,j,l,n) = adrm (  1,j,l,n) + adrm (iie+1,j,l,n)
                adrm (iie+1,j,l,n) = 0.0
                adrxm(  1,j,l,n) = adrxm(  1,j,l,n) + adrxm(iie+1,j,l,n)
                adrxm(iie+1,j,l,n) = 0.0
                adrym(  1,j,l,n) = adrym(  1,j,l,n) + adrym(iie+1,j,l,n)
                adrym(iie+1,j,l,n) = 0.0
                adrzm(  1,j,l,n) = adrzm(  1,j,l,n) + adrzm(iie+1,j,l,n)
                adrzm(iie+1,j,l,n) = 0.0
             enddo
          enddo
       enddo
    endif

    ! no failues yet:
    nfail = 0

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( region ) &
    !$OMP   shared ( imr ) &
    !$OMP   shared ( xcyc ) &
    !$OMP   shared ( is, ie, js, je, ls, le ) &
    !$OMP   shared ( special_grid, imred ) &
    !$OMP   shared ( m, am ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( adrm ) &
    !$OMP   shared ( adrxm, adrym, adrzm ) &
    !$OMP   reduction ( + : nfail ) &
    !$OMP   private ( status ) &
    !$OMP   private ( i, j, l ) &
    !$OMP   private ( iie ) &
    !$OMP   private ( iloop, nloop ) &
    !$OMP   private ( mnew ) &
    !$OMP   private ( adf, adpf, adfy, adfz ) &
    !$OMP   private ( alpha )

    ! cells:
    allocate( mnew(imr) )
    ! interfaces:
    allocate( adf (0:imr) )
    allocate( adpf(0:imr) )
    allocate( adfy(0:imr) )
    allocate( adfz(0:imr) )

    ! init fluxes to zero:
    adf  = 0.0
    adfy = 0.0
    adfz = 0.0
    adpf = 0.0

    !$OMP   DO
    do l = ls, le
       do j = js, je

          ! number of points in x-direction in case of reduced grid:
          if ( special_grid ) then
            iie = imred(j,region)
          else
            iie = ie
          end if

          do i = iie, is, -1
             mnew(i) = m(i,j,l)   ! mnew is defined as m at t+dt
             m(i,j,l) = mnew(i) + am(i,j,l) - am(i-1,j,l)   ! adjoint transport
          enddo

          ! get number of loops,
          ! divide am by number of loops:
          call advectx_get_nloop( is,iie, xcyc(region)==1, &
                                   m(is-1:iie+1,j,l), &
                                   am(is-1:iie,j,l), &
                                   nloop, status )
          if (status/=0) then
            nfail = nfail + 1
            cycle
          end if

          ! CFL loop
          do iloop = 1, nloop

             do i = iie, is, -1
                m(i,j,l) = mnew(i) + am(i,j,l) - am(i-1,j,l)   ! adjoint transport
             end do
             if ( xcyc(region) == 1 ) then
                m(0    ,j,l) = m(iie,j,l) !periodic BCs...
                m(iie+1,j,l) = m(1,  j,l)
             end if

             ! loop over tracers
             do n = 1, ntracetloc

                do i = iie, is, -1
                   adfz(i-1) = adfz(i-1)+adrzm(i,j,l,n)
                   adfz(i) = adfz(i)-adrzm(i,j,l,n)
                   adfy(i-1) = adfy(i-1)+adrym(i,j,l,n)
                   adfy(i) = adfy(i)-adrym(i,j,l,n)
                   adf(i-1) = adf(i-1)-adrxm(i,j,l,n)*(3.*m(i,j,l)/mnew(i))
                   adf(i) = adf(i)-adrxm(i,j,l,n)*(3.*m(i,j,l)/mnew(i))
                   adpf(i-1) = adpf(i-1)+adrxm(i,j,l,n)/(mnew(i))
                   adpf(i) = adpf(i)-adrxm(i,j,l,n)/(mnew(i))
                   adrm(i,j,l,n) = adrm(i,j,l,n)+adrxm(i,j,l,n)* &
                        (3.*(am(i-1,j,l)+am(i,j,l))/ mnew(i))
                   adrxm(i,j,l,n) = adrxm(i,j,l,n)* &
                        (1-(am(i-1,j,l)-am(i,j,l))/ mnew(i))
                   adf(i-1) = adf(i-1)+adrm(i,j,l,n)
                   adf(i) = adf(i)-adrm(i,j,l,n)
                end do

                if ( xcyc(region) == 1 ) then
                  ! adjoint of redefinition left fluxes using periodicity:
                  !if (grid_reduced.and.nred(region).ne.0) imr = imred(j,region)      ! reduced grid
                  adf (iie) = adf (iie) + adf(0)
                  adf (0  ) = 0.0
                  adpf(iie) = adpf(iie) + adpf(0)
                  adpf(0  ) = 0.0
                  adfy(iie) = adfy(iie) + adfy(0)
                  adfy(0  ) = 0.0
                  adfz(iie) = adfz(iie) + adfz(0)
                  adfz(0  ) = 0.0
                endif

                do i = iie, is-1, -1
                   if (am(i,j,l) >= 0.) then
                      alpha = am(i,j,l)/m(i,j,l)
                      adrzm(i,j,l,n) = adrzm(i,j,l,n)+adfz(i)*alpha
                      adfz(i) = 0.
                      adrym(i,j,l,n) = adrym(i,j,l,n)+adfy(i)*alpha
                      adfy(i) = 0.
                      adf(i) = adf(i)-3*adpf(i)*am(i,j,l)
                      adrxm(i,j,l,n) = adrxm(i,j,l,n)+ &
                           adpf(i)*am(i,j,l)*alpha*alpha
                      adpf(i) = 0.
                      adrm(i,j,l,n) = adrm(i,j,l,n)+adf(i)*alpha
                      adrxm(i,j,l,n) = adrxm(i,j,l,n)+ &
                           adf(i)*alpha*(1.-alpha)
                      adf(i) = 0.
                   else
                      alpha = am(i,j,l)/ m(i+1,j,l)
                      adrzm(i+1,j,l,n) = adrzm(i+1,j,l,n)+adfz(i)*alpha
                      adfz(i) = 0.
                      adrym(i+1,j,l,n) = adrym(i+1,j,l,n)+adfy(i)*alpha
                      adfy(i) = 0.
                      adf(i) = adf(i)-3*adpf(i)*am(i,j,l)
                      adrxm(i+1,j,l,n) = adrxm(i+1,j,l,n)+ &
                           adpf(i)*am(i,j,l)*alpha*alpha
                      adpf(i) = 0.
                      adrm(i+1,j,l,n) = adrm(i+1,j,l,n)+adf(i)*alpha
                      adrxm(i+1,j,l,n) = adrxm(i+1,j,l,n)- &
                           adf(i)*alpha*(1.+alpha)
                      adf(i) = 0.
                   endif
                end do  !i

                ! Circ. Boundary conditions within iloop
                if ( xcyc(region) == 1 ) then
                   ! JFM: added update adrm( iie,: ), etc. July 2006
                   adrm (    1,j,l,n) = adrm (  1,j,l,n) + adrm( iie+1,j,l,n )
                   adrm (iie+1,j,l,n) = 0.0
                   adrxm(    1,j,l,n) = adrxm(  1,j,l,n) + adrxm( iie+1,j,l,n )
                   adrxm(iie+1,j,l,n) = 0.0
                   adrym(    1,j,l,n) = adrym(  1,j,l,n) + adrym( iie+1,j,l,n )
                   adrym(iie+1,j,l,n) = 0.0
                   adrzm(    1,j,l,n) = adrzm(  1,j,l,n) + adrzm( iie+1,j,l,n )
                   adrzm(iie+1,j,l,n) = 0.0
                   adrm (iie  ,j,l,n) = adrm (iie,j,l,n) + adrm( 0,j,l,n )
                   adrm (    0,j,l,n) = 0.0
                   adrxm(iie  ,j,l,n) = adrxm(iie,j,l,n) + adrxm( 0,j,l,n )
                   adrxm(    0,j,l,n) = 0.0
                   adrym(iie  ,j,l,n) = adrym(iie,j,l,n) + adrym( 0,j,l,n )
                   adrym(    0,j,l,n) = 0.0
                   adrzm(iie  ,j,l,n) = adrzm(iie,j,l,n) + adrzm( 0,j,l,n )
                   adrzm(0    ,j,l,n) = 0.0
                endif

             end do  !n
             ! in the next step mnew will be m!
             mnew(is:iie) = m(is:iie,j,l)
          end do ! cfl loop

          ! restore 'old' am
          if ( nloop > 1 ) then
            am(is-1:iie,j,l) = am(is-1:iie,j,l)*nloop
          end if

       end do  !j
    end do  !l
    !$OMP   END DO

    ! clear:
    deallocate( mnew )
    deallocate( adf  )
    deallocate( adpf )
    deallocate( adfy )
    deallocate( adfz )

    !$OMP END PARALLEL

    ! check ...
    if ( nfail > 0 ) then
      write (gol,'("failures from x advection : ",i6)') nfail; call goErr
      TRACEBACK; status=1; return
    end if

    ! clear:
    nullify(am)
    nullify(m)
    nullify(adrm)
    nullify(adrxm)
    nullify(adrym)
    nullify(adrzm)

    call adj_compress_xedges( region, js,je,ls,le, status )
    IF_NOTOK_RETURN(status=1)

    call GO_Timer_End( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_dynamu


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

  subroutine adj_put_xedges( region, status )

    use GO         , only : GO_Timer_Start, GO_Timer_End
    use dims       , only : im, jm, lm
    use dims       , only : rstatus => status
    use dims       , only : nsplitsteps, splitorderzoom, n_operators
    use dims       , only : children
    use dims       , only : xref, yref, zref
    use dims       , only : ibeg, iend, jbeg, jend, lbeg, lend
    use dims       , only : xcyc
    use dims       , only : okdebug
    use global_data, only: mass_dat, wind_dat
    use MeteoData  , only : m_dat
    use chem_param , only: ra, ntracet
    use Partools   , only: ntracetloc

    ! --- in/out ----------------------------------

    integer,intent(in)        ::  region
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_put_xedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer              :: rm, rxm, rym, rzm
    real,dimension(:,:,:,:),pointer              :: rmc, rxmc, rymc, rzmc
    real,dimension(:,:,:),  pointer              :: m
    real,dimension(:,:,:),  pointer              :: mc
    real,dimension(:,:,:),  pointer              :: am,bm,cm
    integer            :: child,ichild
    integer            :: i, j, l, n
    integer            :: ipw, ipe, jp, lp
    integer            :: imc, jmc,lmc
    integer            :: xref_,yref_,zref_
    real               :: yzref,xyzref
    real               :: mpw,mpe
    integer            :: nzne, nzne_v,nznem,imr,q
    logical            :: xyz
    character(len=n_operators) :: stencil

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_put_xedges, status )
    IF_NOTOK_RETURN(status=1)

    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    m => m_dat(region)%data
    am => wind_dat(region)%am_t
    bm => wind_dat(region)%bm_t
    cm => wind_dat(region)%cm_t

    q=(rstatus(region)-1)/((nsplitsteps/2))
    ! find q - the place in the splitorderzoom of the parent
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

      ! get the number of the child region
      child = children(region,ichild)

      ! child is [0,360] degrees wide ? skip it since
      ! periodic boundary conditions will be used
      if ( xcyc(child) == 1 ) then
         if (okdebug)  print *,'put_xedges: child ',child,' skipped'
         cycle
      endif

      xref_ = xref(child)/xref(region)
      yref_ = yref(child)/yref(region)
      zref_ = zref(child)/zref(region)
      yzref = 1./(yref_*zref_)
      xyzref = 1./(xref_*yref_*zref_)

      imc = im(child)
      jmc = jm(child)
      lmc = lm(child)

      ipw = ibeg(child)
      ipe = iend(child)

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
      !$OMP   shared ( ibeg, iend, jbeg, jend, lbeg, lend ) &
      !$OMP   shared ( ipw, ipe ) &
      !$OMP   shared ( xyz ) &
      !$OMP   shared ( xref_, yref_, zref_ ) &
      !$OMP   shared ( yzref, xyzref ) &
      !$OMP   shared ( m, bm, cm ) &
      !$OMP   shared ( rm ) &
      !$OMP   shared ( rxm, rym, rzm ) &
      !$OMP   shared ( mc, rmc ) &
      !$OMP   shared ( rxmc, rymc, rzmc ) &
      !$OMP   private ( i, j, l ) &
      !$OMP   private ( jp, lp ) &
      !$OMP   private ( mpw, mpe )
      !$OMP   DO
      do l=1,lmc
        ! parent layer:
        lp = lbeg(child) + (l-1)/zref_

        if(.not.xyz) then
          ! put_xedges puts masses of the parent in the child.
          ! In the case of a XYZ sequence, these masses did not see any advection.
          ! In the case of the parent sequence ZYX, these masses have seen a full ZY advection.
          ! In the adjoint, we have to advect these masses of the child back in time.
          !
          do j=1,jmc
             jp = jbeg(child) + (j-1)/yref_
             mpw = m(ipw,jp,l) - bm(ipw,jp,l) + bm(ipw,jp+1,l) + cm(ipw,jp,l) - cm(ipw,jp,l-1)
             mpe = m(ipe,jp,l) - bm(ipe,jp,l) + bm(ipe,jp+1,l) + cm(ipe,jp,l) - cm(ipe,jp,l-1)
             mc(1:xref_,j,l) = mpw/(yref_*xref_)
             mc(imc-xref_+1:imc,j,l) = mpe/(yref_*xref_)
           enddo ! j
        endif

        do j=1,jmc
           jp = jbeg(child) + (j-1)/yref_
           rm(ipw-1,jp,lp,:) = rm(ipw-1,jp,lp,:) + rmc(0,j,l,:)*yzref
           rmc(0,j,l,:) = 0.0
           rm(ipe+1,jp,lp,:) = rm(ipe+1,jp,lp,:) + rmc(imc+1,j,l,:)*yzref
           rmc(imc+1,j,l,:) = 0.0
           do i=1,xref_
              rm(ipw,jp,lp,:) = rm(ipw,jp,lp,:) + rmc(i,j,l,:)*xyzref
              rmc(i,j,l,:) = 0.0
              rm(ipe,jp,lp,:) = rm(ipe,jp,lp,:) + rmc(imc+1-i,j,l,:)*xyzref
              rmc(imc+1-i,j,l,:) = 0.0
           enddo
        enddo

        do j=1,jmc
           jp = jbeg(child) + (j-1)/yref_
           rxm(ipw-1,jp,lp,:) = rxm(ipw-1,jp,lp,:) + rxmc(0,j,l,:)*yzref
           rxmc(0,j,l,:) = 0.0
           rxm(ipe+1,jp,lp,:) = rxm(ipe+1,jp,lp,:) + rxmc(imc+1,j,l,:)*yzref
           rxmc(imc+1,j,l,:) = 0.0
           do i=1,xref_
              rxm(ipw,jp,lp,:) = rxm(ipw,jp,lp,:) + rxmc(i,j,l,:)*yzref
              rxmc(i,j,l,:) = 0.0
              rxm(ipe,jp,lp,:) = rxm(ipe,jp,lp,:) + rxmc(imc+1-i,j,l,:)*yzref
              rxmc(imc+1-i,j,l,:) = 0.0
           enddo
        enddo

        do j=1,jmc
           jp = jbeg(child) + (j-1)/yref_
           rym(ipw-1,jp,lp,:) = rym(ipw-1,jp,lp,:) + rymc(0,j,l,:)*yzref
           rymc(0,j,l,:) = 0.0
           rym(ipe+1,jp,lp,:) = rym(ipe+1,jp,lp,:) + rymc(imc+1,j,l,:)*yzref
           rymc(imc+1,j,l,:) = 0.0
           do i=1,xref_
              rym(ipw,jp,lp,:) = rym(ipw,jp,lp,:) + rymc(i,j,l,:)*yzref
              rymc(i,j,l,:) = 0.0
              rym(ipe,jp,lp,:) = rym(ipe,jp,lp,:) + rymc(imc+1-i,j,l,:)*yzref
              rymc(imc+1-i,j,l,:) = 0.0
           end do
        end do

        do j=1,jmc
           jp = jbeg(child) + (j-1)/yref_
           rzm(ipw-1,jp,lp,:) = rzm(ipw-1,jp,lp,:) + rzmc(0,j,l,:)*yzref
           rzmc(0,j,l,:) = 0.0
           rzm(ipe+1,jp,lp,:) = rzm(ipe+1,jp,lp,:) + rzmc(imc+1,j,l,:)*yzref
           rzmc(imc+1,j,l,:) = 0.0
           do i=1,xref_
              rzm(ipw,jp,lp,:) = rzm(ipw,jp,lp,:) + rzmc(i,j,l,:)*yzref
              rzmc(i,j,l,:) = 0.0
              rzm(ipe,jp,lp,:) = rzm(ipe,jp,lp,:) + rzmc(imc+1-i,j,l,:)*yzref
              rzmc(imc+1-i,j,l,:) = 0.0
          end do
        end do

      end do   ! layers
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

    call GO_Timer_End( itim_put_xedges, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_put_xedges


  ! ***


  ! this is for slope only: condense data at the x-edges of the zoom region
  ! to allow for uniform work in dynamu/dynamu1
  ! written by mike botchev, march-june 1999
  !
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

  subroutine adj_compress_xedges( region, js,je,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only: xref,yref,parent, zref, im,jm,okdebug
    use chem_param,  only: ntracet
    use ParTools  ,  only: ntracetloc
    use global_data, only: mass_dat
    use MeteoData  , only : m_dat

    ! --- in/out ----------------------------------

    integer, intent(in)       ::  region
    integer, intent(in)       ::  js,je,ls,le
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_compress_xedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer    :: rm,rxm, rym,rzm
    real,dimension(:,:,:),pointer      :: m
    integer     ::  imr,jmr,lmr
    integer     ::  xref_, i,yref_,zref_,j,l, n
    real        ::  medge1, medge2

    ! --- in/out ---------------------------------

    call GO_Timer_Start( itim_compress, status )
    IF_NOTOK_RETURN(status=1)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region) ; jmr = jm(region)

    if ((xref_==1).or.(im(region)/xref(region)==im(1))) then
       if(okdebug) then
         print *, "       compress_xedges: no refinement or periodic bc's, nothnig to do"
       end if
       call GO_Timer_End( itim_compress, status )
       IF_NOTOK_RETURN(status=1)
       return
    end if

    if(okdebug) then
       print *,'       call adj_compress_xedges, region=',region,' imr=',imr
       print *,'       cells to be updated: ',xref_-1,xref_,imr-xref_+1,imr-xref_+2
    endif

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( imr ) &
    !$OMP   shared ( js, je, ls, le ) &
    !$OMP   shared ( xref_ ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m ) &
    !$OMP   shared ( rm ) &
    !$OMP   shared ( rxm, rym, rzm ) &
    !$OMP   private ( i, j, l, n ) &
    !$OMP   private ( medge1, medge2 )
    !$OMP   DO
    do l = ls, le
      do j = js, je

        ! for m, just calculate the masses BEFORE compress_xedges
        m (0,j,l) = m(xref_ -1,j,l)
        m(imr+1,j,l) = m(imr-xref_+2,j,l)
        medge1 = m(xref_,j,l)/xref_
        medge2 = m(imr-xref_+1,j,l)/xref_
        do i=1,xref_
           m(i,j,l) = medge1
           m(imr-i+1,j,l) = medge2
        end do

        ! for rm, rxm do the adjoint operation
        do n=1,ntracetloc
           rm (0,j,l,n) = rm(xref_ -1,j,l,n)
           rm(imr+1,j,l,n) = rm(imr-xref_+2,j,l,n)
           medge1 = rm(xref_,j,l,n)
           medge2 = rm(imr-xref_+1,j,l,n)
           do i=1,xref_
              rm(i,j,l,n) = medge1
              rm(imr-i+1,j,l,n) = medge2
           end do
           ! for the slopes, the 'mixed' interface cells are not 'treated'  >> should be verified !CMK ADJ
           rxm (0,j,l,n) = rxm(xref_ -1,j,l,n) + rxm(0,j,l,n)
           rxm (xref_-1,j,l,n) = 0.0
           rxm(imr+1,j,l,n) = rxm(imr+1,j,l,n) + rxm(imr-xref_+2,j,l,n)
           rxm(imr-xref_+2,j,l,n) = 0.0
           rym (0,j,l,n) = rym(xref_ -1,j,l,n) + rym(0,j,l,n)
           rym (xref_-1,j,l,n) = 0.0
           rym(imr+1,j,l,n) = rym(imr+1,j,l,n) + rym(imr-xref_+2,j,l,n)
           rym(imr-xref_+2,j,l,n) = 0.0
           rzm (0,j,l,n) = rzm(xref_ -1,j,l,n) + rzm(0,j,l,n)
           rzm (xref_-1,j,l,n) = 0.0
           rzm(imr+1,j,l,n) = rzm(imr+1,j,l,n) + rzm(imr-xref_+2,j,l,n)
           rzm(imr-xref_+2,j,l,n) = 0.0
        end do ! n

      end do ! j
    end do ! l
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

  end subroutine adj_compress_xedges


  ! ***


  ! distribute data at the x-edges of the zoom region
  ! over the whole interface cell
  ! written by mike botchev, march-june 1999
  !
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

  subroutine adj_uncompress_xedges( region, js,je,ls,le, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims,        only: parent, xref, yref, zref, im,jm, okdebug
    use global_data, only: mass_dat
    use MeteoData  , only : m_dat
    use chem_param,  only: ntracet
    use ParTools  ,  only: ntracetloc

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  region
    integer, intent(in)       ::  js,je,ls,le
    integer,intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_uncompress_xedges'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer    :: rm,rxm, rym,rzm
    real,dimension(:,:,:),pointer      :: m
    integer     ::  imr, jmr, lmr
    integer     ::  xref_, yref_, zref_, i,j,l,n
    real        ::  m_edge,rm_edge

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_uncompress, status )
    IF_NOTOK_RETURN(status=1)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    xref_ = xref(region)/xref(parent(region))
    imr = im(region)

    if ((xref_==1).or.(im(region)/xref(region)==im(1))) then
      if(okdebug) then
        print *, "     adj_uncompress_xedges: no refinement or periodic bc's, nothnig to do"
      end if
      call GO_Timer_End( itim_uncompress, status )
      IF_NOTOK_RETURN(status=1)
      return
    end if
    if(okdebug) then
       print *,'   call adj_uncompress_xedges, region=',region,' imr=',imr
       print *,'   cells to be updated: ',1,':',xref_,' ',imr-xref_+1,':',imr
    endif

    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( imr ) &
    !$OMP   shared ( js, je, ls, le ) &
    !$OMP   shared ( xref_ ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( m ) &
    !$OMP   shared ( rm ) &
    !$OMP   shared ( rxm, rym, rzm ) &
    !$OMP   private ( j, l, n )
    !$OMP   DO
    do l=ls,le
       do j=js,je

         ! for m, calculate the masses BEFORE uncompress edges
          m (xref_,j,l) = sum(m(1:xref_,j,l))
          m (xref_-1,j,l) = m(0,j,l)
          m (0:xref_-2,j,l) = 0.0
          m (imr-xref_+1,j,l) = sum(m(imr-xref_+1:imr,j,l))
          m (imr-xref_+2,j,l) = m(imr+1,j,l)
          m (imr-xref_+3:imr+1,j,l) = 0.0

          do n=1,ntracetloc

             rm (xref_,j,l,n) = sum(rm(1:xref_,j,l,n))/xref_
             rm (1:xref_-1,j,l,n) = 0.0
             rm (imr-xref_+1,j,l,n) = sum(rm(imr-xref_+1:imr,j,l,n))/xref_
             rm (imr-xref_+2:imr,j,l,n) = 0.0

             rxm (xref_,j,l,n) = sum(rxm(1:xref_,j,l,n))/xref_
             rxm (1:xref_-1,j,l,n) = 0.0
             rxm (imr-xref_+1,j,l,n) = sum(rxm(imr-xref_+1:imr,j,l,n))/xref_
             rxm (imr-xref_+2:imr,j,l,n) = 0.0

             rym (xref_,j,l,n) = sum(rym(1:xref_,j,l,n))/xref_
             rym (1:xref_-1,j,l,n) = 0.0
             rym (imr-xref_+1,j,l,n) = sum(rym(imr-xref_+1:imr,j,l,n))/xref_
             rym (imr-xref_+2:imr,j,l,n) = 0.0

             rzm (xref_,j,l,n) = sum(rzm(1:xref_,j,l,n))/xref_
             rzm (1:xref_-1,j,l,n) = 0.0
             rzm (imr-xref_+1,j,l,n) = sum(rzm(imr-xref_+1:imr,j,l,n))/xref_
             rzm (imr-xref_+2:imr,j,l,n) = 0.0

          enddo  ! n

        enddo  ! j
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

  end subroutine adj_uncompress_xedges


end module adj_advectx
