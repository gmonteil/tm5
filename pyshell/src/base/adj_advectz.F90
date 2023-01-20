!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_AdvectZ

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------------

  private

  public :: adj_AdvectZ_Init, adj_AdvectZ_Done
  public :: adj_advectzzoom


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'adj_AdvectZ'


  ! --- local ------------------------------------

  integer    ::  itim_zoom
  integer    ::  itim_dynam


contains


  ! ====================================================================


  subroutine adj_AdvectZ_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_AdvectZ_Init'

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_zoom, 'adj_advectz zoom', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_dynam, 'adj_advectz dynam', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_AdvectZ_Init


  ! ***


  subroutine adj_AdvectZ_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_AdvectZ_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine adj_AdvectZ_Done


  ! ***


  subroutine adj_advectzzoom( region, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use dims, only    : im,jm,lm, xref, yref, zref, tref
    use dims, only    : n_operators, splitorderzoom, nsplitsteps, xcyc
    use dims, only    : touch_sp, touch_np, zoom2D, rstatus => status, parent
    use toolbox, only : escape_tm

    ! --- in/out ---------------------------------

    integer,intent(in) :: region
    integer,intent(out) :: status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_advectzzoom'

    ! --- local ----------------------------------

    integer                       :: is,ie,js,je,n,q
    integer                       :: imr,jmr,lmr,tref_,xref_,yref_,zref_,my_parent
    logical                       :: z_encountered
    character(len=1)              :: dir

    ! --- begin ----------------------------------

    call GO_Timer_Start( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    ! write cells along the z-walls of each child
    ! to rm(0,...),rm(imr+1,...) and to interface cells of the child

    ! this option should be there: zooming in z disabled....
    if ( .not. zoom2D ) then
      write (gol,'("zoom2D should be true: zooming along z not allowed")'); call goErr
      TRACEBACK; status=1; return
    end if

    !! for inspiration in case of 3D zooming:
    !CMKCALL put_zedges(region,rm,rxm,rym,rzm,m,  &
    !CMK                        rmgl,rxmgl,rymgl,rzmgl,mgl)

    tref_ = tref(region)/tref(parent(region))
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    zref_ = zref(region)/zref(parent(region))
    imr = im(region);  jmr = jm(region);  lmr = lm(region)

    ! determine the scope for advectz:

    if (region==1) then
       xref_ = 0; yref_ = 0  ! to have is/ie and js/je properly computed
    end if
    my_parent = parent(region)

    if(my_parent == 0) then   ! always full scope: no parent
       is = 1
       ie = im(region)
       js = 1
       je = jm(region)
    else
      q=(rstatus(my_parent)-1)/((nsplitsteps/2))      ! find q - the place in the splitorderzoom
                                               ! corresponding to the begining of the last
                                               ! processing order of the parent.
      z_encountered=.false.

      do n=1,n_operators                                  ! now track following steps from q

        dir=splitorderzoom(my_parent,q*(nsplitsteps/2)+n)

        select case(dir)
        case('x')
           IF  ((.not.z_encountered).or.(xcyc(region).eq.1)) THEN
              is=1                             ! x-substep is before z =>
              ie=im(region)                    ! full i-scope
           ELSE
              is=xref_+1                       ! z-substep is before x =>
              ie=im(region)-xref_              ! restricted i-scope
           ENDIF
        case('y')
           IF (.not.z_encountered) THEN
              js=1                             ! y-substep is before z =>
              je=jm(region)                    ! full j-scope
           ELSE
              IF(touch_sp(region).eq.1) THEN
                 js=1                             ! sp is southern edge--->full scope.
              ELSE
                 js=yref_+1                       ! x-substep is before y =>
              ENDIF
              IF(touch_np(region).eq.1) THEN
                je=jm(region)              ! np is northern edge--->full scope
              ELSE
                je=jm(region)-yref_        ! restricted j-scope
              ENDIF
           ENDIF
        case('z')
           z_encountered=.true.
        case ('c')
        case ('v')
        case ('d')
        case ('s')
        case default
           print *,'strange value in splitorderzoom(',region,',',   &
                q*(nsplitsteps/2)*tref_+n,'): ',q
           call escape_tm('Error in advectz')
        end select

      enddo
    endif  ! my_parent = 0

    call adj_dynamw( region, is,ie,js,je, status )
    IF_NOTOK_RETURN(status=1)

    call GO_Timer_End( itim_zoom, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_advectzzoom

  !-----------------------------------------------------------------------
  !
  !**** dynamw      - vertical tracer transport     v 9.1
  !
  ! programmed by       mh  mpi HH      23-feb-1995
  !
  ! purpose
  ! -------
  ! calculate amount of tracer moved in a vertical advection
  ! substep
  !
  ! interface
  ! ---------
  ! call dynamw
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
  !       changed order of loops for increased performance
  !                               mh, 11-jun-1994
  !     included code for limits of slopes to prevent negative tracer
  !     masses                   mh, 20-jun-1994
  !
  !     zoom version written by mike botchev, march-june 1999
  !
  !     modified for the zoom version by Maarten Krol, March 2003
  !-----------------------------------------------------------------------!

  !***************************************************************
  !***************************************************************
  !** This routine was generated by the                         **
  !** Tangent linear and Adjoint Model Compiler,    TAMC 4.83   **
  !***************************************************************
  !***************************************************************

  subroutine adj_dynamw( region, is,ie,js,je, status )

    use GO            , only : GO_Timer_Start, GO_Timer_End
    use dims,           only : im, jm, lm
    use dims,           only : nregions
    use dims,           only : xref, yref
    use dims,           only : parent
    use global_data,    only : wind_dat, mass_dat
    use MeteoData     , only : m_dat
    use ParTools  ,     only : ntracetloc
    use adj_zoom_tools, only : adj_mix_edges

    ! --- in/out ---------------------------------

    integer,intent(in)  ::  region
    integer,intent(in)  ::  is
    integer,intent(in)  ::  ie
    integer,intent(in)  ::  js
    integer,intent(in)  ::  je
    integer,intent(out) ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/adj_dynamw'

    ! --- local ----------------------------------

    real,dimension(:,:,:,:),pointer   :: adrm,adrxm,adrym,adrzm
    real,dimension(:,:,:),pointer     :: m,cm
    integer                           :: imr, jmr, lmr
    integer                           :: i, j
    integer                           :: n
    integer                           :: xref_, yref_

    ! --- begin ---------------------------------

    call GO_Timer_Start( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    ! region dimensions:
    imr=im(region) ; jmr=jm(region) ; lmr=lm(region)

    call adj_mix_edges( region, status )
    IF_NOTOK_RETURN(status=1)

    adrm => mass_dat(region)%rm_t
    adrxm => mass_dat(region)%rxm_t
    adrym => mass_dat(region)%rym_t
    adrzm => mass_dat(region)%rzm_t
    m => m_dat(region)%data
    cm => wind_dat(region)%cm_t

    ! check region:
    if ((region<0).or.(region>nregions)) then
      write (gol,'("illegal number of region: ",i6)') region
      TRACEBACK; status=1; return
    end if

    ! compute refinement factors with respect to the parent
    xref_ = xref(region)/xref(parent(region))
    yref_ = yref(region)/yref(parent(region))
    ! check is,ie,js,je:
    if ( ( is /= xref_+1 ) .and. ( is /= 1 ) ) then
      write (gol,'("Wrong value for IS in dynamw : ",i6)') is; call goErr
      TRACEBACK; status=1; return
    end if
    if ( ( ie /= imr-xref_ ) .and. ( ie /= imr ) ) then
      write (gol,'("Wrong value for IE in dynamw : ",i6)') is; call goErr
      TRACEBACK; status=1; return
    end if
    if ( ( js /= yref_+1 ) .and. ( js /= 1 ) ) then
      write (gol,'("Wrong value for JS in dynamw : ",i6)') is; call goErr
      TRACEBACK; status=1; return
    end if
    if ( ( je /= jmr-yref_ ) .and. ( je /= jmr ) ) then
      write (gol,'("Wrong value for JE in dynamw : ",i6)') is; call goErr
      TRACEBACK; status=1; return
    end if

    ! check ...
    if ( any(cm(:,:,0) /= 0.0) .or. any(cm(:,:,lmr) /= 0.0) ) then
      write (gol,'("found non-zero halo interfaces in cm")'); call goErr
      TRACEBACK; stop
    end if

    ! loop over cells:
    !$OMP PARALLEL &
    !$OMP   default ( none ) &
    !$OMP   shared ( is, ie, js, je, lmr ) &
    !$OMP   shared ( m, cm ) &
    !$OMP   shared ( ntracetloc ) &
    !$OMP   shared ( adrm ) &
    !$OMP   shared ( adrxm, adrym, adrzm ) &
    !$OMP   private ( i, j )
    !$OMP   DO
    do j = js, je
      do i = is, ie

        ! column advection:
        call adj_dynamw_1d( lmr, cm(i,j,0:lmr), m(i,j,1:lmr), &
                            ntracetloc, &
                            adrm(i,j,1:lmr,:), &
                            adrxm(i,j,1:lmr,:), adrym(i,j,1:lmr,:), adrzm(i,j,1:lmr,:) )

      end do  ! i
    end do   ! j
    !$OMP   END DO
    !$OMP END PARALLEL

    ! clear:
    nullify(adrm)
    nullify(adrxm)
    nullify(adrym)
    nullify(adrzm)
    nullify(m)
    nullify(cm)

    call GO_Timer_End( itim_dynam, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_dynamw


  ! *


  pure subroutine adj_dynamw_1d( lmr, cm, m, &
                        ntr, adrm, adrxm, adrym, adrzm )

    ! --- in/out ---------------------------------

    integer, intent(in)   ::  lmr
    real, intent(in)      ::  cm(0:lmr)
    real, intent(inout)   ::  m(1:lmr)
    integer, intent(in)   ::  ntr
    real, intent(inout)   ::  adrm(1:lmr,1:ntr)
    real, intent(inout)   ::  adrxm(1:lmr,1:ntr)
    real, intent(inout)   ::  adrym(1:lmr,1:ntr)
    real, intent(inout)   ::  adrzm(1:lmr,1:ntr)

    ! --- local ----------------------------------

    integer              ::  l, n
    real                 ::  mnew(1:lmr)
    real                 ::  adf (0:lmr)
    real                 ::  adpf(0:lmr)
    real                 ::  adfx(0:lmr)
    real                 ::  adfy(0:lmr)
    real                 ::  gamma

    ! --- begin ----------------------------------

    ! new mass:
    do l = 1, lmr
      mnew(l) = m(l)
      m(l) =  mnew(l) - cm(l-1) + cm(l)
    end do

    ! init boundary fluxes to zero:
    adf  = 0.0
    adpf = 0.0
    adfx = 0.0
    adfy = 0.0

    ! loop over tracers:
    do n = 1, ntr
      adfy(lmr-1) = adfy(lmr-1) + adrym(lmr,n)
      adfx(lmr-1) = adfx(lmr-1) + adrxm(lmr,n)
      adf (lmr-1) = adf (lmr-1) - adrzm(lmr,n) * (3.0*m(lmr)/mnew(lmr))
      adpf(lmr-1) = adpf(lmr-1) + adrzm(lmr,n) / mnew(lmr)
      adrm(lmr,n) = adrm(lmr,n) + adrzm(lmr,n) * (3.*cm(lmr-1)/mnew(lmr))
      adrzm(lmr,n) = adrzm(lmr,n) * (1-cm(lmr-1)/mnew(lmr))
      adf(lmr-1) = adf(lmr-1)+adrm(lmr,n)
      adfy(1) = adfy(1)-adrym(1,n)
      adfx(1) = adfx(1)-adrxm(1,n)
      adf(1) = adf(1)-adrzm(1,n)*(3.*m(1)/ mnew(1))
      adpf(1) = adpf(1)-adrzm(1,n)/mnew(1)
      adrm(1,n) = adrm(1,n)+adrzm(1,n)*(3.*cm(1)/mnew(1))
      adrzm(1,n) = adrzm(1,n)*(1+cm(1)/mnew(1))
      adf(1) = adf(1)-adrm(1,n)
      do l = lmr-1, 2, -1
        adfy(l-1) = adfy(l-1)+adrym(l,n)
        adfy(l) = adfy(l)-adrym(l,n)
        adfx(l-1) = adfx(l-1)+adrxm(l,n)
        adfx(l) = adfx(l)-adrxm(l,n)
        adf(l-1) = adf(l-1)-adrzm(l,n)*(3.*m(l)/mnew(l))
        adf(l) = adf(l)-adrzm(l,n)*(3.*m(l)/mnew(l))
        adpf(l-1) = adpf(l-1)+adrzm(l,n)/mnew(l)
        adpf(l) = adpf(l)-adrzm(l,n)/mnew(l)
        adrm(l,n) = adrm(l,n)+adrzm(l,n)* (3.*(cm(l-1)+cm(l))/mnew(l))
        adrzm(l,n) = adrzm(l,n)*(1-(cm(l-1)-cm(l))/mnew(l))
        adf(l-1) = adf(l-1)+adrm(l,n)
        adf(l) = adf(l)-adrm(l,n)
      end do  ! l
      do l = lmr-1, 1, -1
        if (cm(l) .ge. 0.) then
          gamma = cm(l)/m(l)
          adrym(l,n) = adrym(l,n)+adfy(l)*gamma
          adfy(l) = 0.
          adrxm(l,n) = adrxm(l,n)+adfx(l)*gamma
          adfx(l) = 0.
          adf(l) = adf(l)-3*adpf(l)*cm(l)
          adrzm(l,n) = adrzm(l,n)+adpf(l)*cm(l)*gamma*gamma
          adpf(l) = 0.
          adrm(l,n) = adrm(l,n)+adf(l)*gamma
          adrzm(l,n) = adrzm(l,n)+adf(l)*gamma*(1.0-gamma)
          adf(l) = 0.
         else
          gamma = cm(l)/m(l+1)
          adrym(l+1,n) = adrym(l+1,n)+adfy(l)*gamma
          adfy(l) = 0.
          adrxm(l+1,n) = adrxm(l+1,n)+adfx(l)*gamma
          adfx(l) = 0.
          adf(l) = adf(l)-3*adpf(l)*cm(l)
          adrzm(l+1,n) = adrzm(l+1,n)+adpf(l)* &
               cm(l)*gamma*gamma
          adpf(l) = 0.
          adrm(l+1,n) = adrm(l+1,n)+adf(l)*gamma
          adrzm(l+1,n) = adrzm(l+1,n)-adf(l)* &
               gamma*(1.+gamma)
          adf(l) = 0.
        endif
      end do   ! l

    end do   ! n

  end subroutine adj_dynamw_1d


end module adj_advectz
