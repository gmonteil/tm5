!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!#################################################################

module TM5_Conv

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  TM5_Conv_Matrix
  public  ::  TM5_Conv_Apply


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'TM5_Conv'


contains


  ! ================================================================


  pure subroutine TM5_Conv_Matrix( &
                        dt, lmx, m, &
                        entu, detu, entd, detd, &
                        ld, li, &
#ifndef without_wet_deposition
                        cp, sceffdeep, &
#endif
                        conv1, lbdcv, lmc, &
                        status )

    ! --- in/out ---------------------------------

    integer, intent(in)         ::  lmx  ! lmax_conv
    real, intent(in)            ::  m(lmx)
    real, intent(in)            ::  entu(lmx)
    real, intent(in)            ::  detu(lmx)
    real, intent(in)            ::  entd(lmx)
    real, intent(in)            ::  detd(lmx)
    integer, intent(in)         ::  ld   ! cloud_lfs
    integer, intent(in)         ::  li   ! cloud_top
#ifndef without_wet_deposition
    real, intent(in)            ::  cp   ! convective precepitation ;
                                         ! actually cp/cp_scale
    real, intent(in)            ::  sceffdeep
#endif
    real, intent(in)            ::  dt
    real, intent(out)           ::  conv1(lmx,lmx)    ! convection matrix
    real, intent(out)           ::  lbdcv(lmx,lmx)    ! wet removal matrix
    integer, intent(out)        ::  lmc               ! convection top layer
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/TM5_Conv_Matrix'

    ! --- local ----------------------------------

    ! Since this is a pure subroutine, we should no longer use static arrays
    ! Instead these arrays will have to be allocatable
    ! tracer fluxes:
    real, allocatable         ::  f(:,:), fu(:,:)
    ! air mass up/downdraft:
    real, allocatable         ::  amu(:), amd(:)

    integer                   ::  k, kk
    integer                   ::  l
    real                      ::  zxi

    ! --- begin ----------------------------------
    allocate(f(0:lmx,lmx), fu(0:lmx,lmx))
    allocate(amu(0:lmx), amd(0:lmx))

    ! calculate convection matrix and directly apply to rm(i,j,*,1:ntracet)
    ! Set all amu's, amd's entu' etc to 0 for gridcells without clouds
    !
    ! updraft   - amu is positive upward
    !
    f (:,:) = 0.0
    fu(:,:) = 0.0
    amu(:)  = 0.0
    do k = 1, li
      amu(k) = amu(k-1) + entu(k) - detu(k)
      if ( amu(k) > 0.0 ) then
        !
        ! we limit zxi from below at 0. in case
        ! there are inconsistent en/detrainment
        ! rates. The case zxi>1 should not occur
        ! since eu and du are >=0.
        !
        zxi = max( 0.0, 1.0-detu(k)/(amu(k-1)+entu(k)) )
      else
        ! set massflux at upper boundary to strictly 0.
        amu(k) = 0.0
        zxi = 0.0
      end if
      ! loop over all previous levels
      do kk = 1, k-1
        fu(k,kk) = fu(k-1,kk) * zxi
      end do
      !cmk: note: the fu(k-1,k) values are zero (see manual Heimann)
      fu(k,k)=entu(k)/m(k)*zxi
    end do !k
    !
    ! downdraft  -  amd is negative downward
    !
    amd(:)=0.
    do  k=ld,2,-1
       amd(k-1)=amd(k)-entd(k)+detd(k)
       if ( amd(k-1) < 0. ) then
          zxi=max(0.,1.+detd(k)/(amd(k)-entd(k)))
       else
          amd(k-1)=0.
          zxi=0.
       end if
       do  kk=k+1,ld   !
          f(k-1,kk)=f(k,kk)*zxi !f is here only fd downdraftmatrix
       end do !kk
       f(k-1,k)=-entd(k)/m(k)*zxi    !note the negatives for j
    end do !k

    !
    ! add coefficients from up and downdraft together
    !
    do k = 1, lmx-1
       do kk = 1, lmx
          f(k,kk) = fu(k,kk) + f(k,kk)
       end do !kk
       ! add coefficient from subsidence
       !CMK SH       f(k,k+1)=f(k,k+1)-(amu(k)+amd(k))/m(k+1)
       ! reported by Sander Houweling.........
       f(k,k+1) = f(k,k+1) - amu(k)/m(k+1)
       f(k,k  ) = f(k,k  ) - amd(k)/m(k  )
    end do

    ! initialize wet-removal
    lbdcv = 0.0

#ifndef without_wet_deposition
    ! fill scaveging stuff if necessary
    !   cp : rain in m/s
    do k = 2, lmx
      do kk = 1, k-1
        ! scale factor between 0 - 1, with ~ 1 for cp >> 0
        lbdcv(k,kk) = sceffdeep * (1.0 - exp(-cp)) * &
                         dt * (fu(k-1,kk) - fu(k,kk))
      end do   ! kk
    end do   ! k
    ! which ensures that the rate of scavenging is equal to dn/dt=-s Mu n
#endif

    ! top of convection, initially zero (no convection):
    lmc = 0
    ! generate forward matrix
    conv1 = 0.0
    do k = 1, lmx
      do kk = 1, lmx
        ! fluxes are different ?
        if ( f(k-1,kk) /= f(k,kk) ) then
          ! change ratio:
          conv1(k,kk) = -dt * ( f(k-1,kk) - f(k,kk) )
          ! update convection top layer:
          lmc = max(max(lmc,k),kk)
        end if
      end do !kk
      conv1(k,k) = conv1(k,k) + 1.0
    end do !k

    !! all layers ...
    !lmc = lmx
    deallocate(f, fu, amu, amd)

    ! ok
    status = 0

  end subroutine TM5_Conv_Matrix


  ! *

#ifdef with_mkl
  pure subroutine TM5_Conv_Apply( trans, &
#else
  subroutine TM5_Conv_Apply( trans, &
#endif
                           lmx, ntr, &
                           conv1, &
                           rm, &
#ifdef slopes
                           rxm, rym, rzm, &
#ifdef secmom
                           rxxm, ryym, rzzm, rxym, rxzm, ryzm
#endif
#endif
#ifndef without_wet_deposition
                           lbdcv, cvsfac, &
                           zwetrm, &
#endif
                           status )

    ! --- use ------------------------------------

#ifdef with_mkl
    use lapack95
#endif

    ! --- in/out ---------------------------------


    character(len=1), intent(in)  ::  trans      ! 'N' for forward and
                                                 ! 'T' for reverse run
    integer, intent(in)         ::  lmx  ! lmax_conv
    integer, intent(in)         ::  ntr
    real, intent(inout)         ::  conv1(lmx,lmx)
#ifndef without_wet_deposition
    real, intent(in)            ::  lbdcv(lmx,lmx)
    real, intent(in)            ::  cvsfac(ntr)   ! cvsfac(nglob)
#endif
    real, intent(inout)         ::  rm(lmx,ntr)
#ifdef slopes
    real, intent(inout)         ::  rxm(lmx,ntr), rym(lmx,ntr), rzm(lmx,ntr)
#ifdef secmom
    real, intent(inout)         ::  rxxm(lmx,ntr), ryym(lmx,ntr), rzzm(lmx,ntr)
    real, intent(inout)         ::  rxym(lmx,ntr), rxzm(lmx,ntr), ryzm(lmx,ntr)
    real, intent(inout)         ::  zwetrm(lmx,ntr)
#endif
#endif
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/TM5_Conv_Apply'

    ! --- local ----------------------------------

    integer, allocatable        ::  ipiv(:)
#ifndef without_wet_deposition
    integer           ::  n
#endif

    ! --- begin ----------------------------------

    allocate(ipiv(lmx))

    ! NOTE: apply same factors to ALL slopes ;
    !   in previous code z-slopes were omitted for not
    !   further described reasons

    ! Solve for all tracers:
    ! If we are using mkl, we will use the F95 interfaces, because they are defined 'pure'
#ifdef with_mkl
    ! LU decomposition (general matrix, triangular factorization)
    call getrf(conv1, ipiv, status)
    if (status/=0) return
    call getrs(conv1, ipiv, rm, trans, status)
    if (status/=0) return
#ifdef slopes
    call getrs(conv1, ipiv, rxm, trans, status)
    if (status/=0) return
    call getrs(conv1, ipiv, rym, trans, status)
    if (status/=0) return
    call getrs(conv1, ipiv, rzm, trans, status)
    if (status/=0) return
#ifdef secmom
    call getrs(conv1, ipiv, rxxm, trans, status)
    if (status/=0) return
    call getrs(conv1, ipiv, ryym, trans, status)
    if (status/=0) return
    call getrs(conv1, ipiv, rzzm, trans, status)
    if (status/=0) return
    call getrs(conv1, ipiv, rxym, trans, status)
    if (status/=0) return
    call getrs(conv1, ipiv, rxzm, trans, status)
    if (status/=0) return
    call getrs(conv1, ipiv, ryzm, trans, status)
    if (status/=0) return
#endif
#endif

#else
    ! LU decomposition (general matrix, triangular factorization)
    call dGeTrf( lmx, lmx, conv1, lmx, ipiv, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rm, lmx, status )
    if (status/=0) return
#ifdef slopes
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rxm, lmx, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rym, lmx, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rzm, lmx, status )
    if (status/=0) return
#ifdef secmom
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rxxm, lmx, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, ryym, lmx, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rzzm, lmx, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rxym, lmx, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, rxzm, lmx, status )
    if (status/=0) return
    call dGeTrs( trans, lmx, ntr, conv1, lmx, ipiv, ryzm, lmx, status )
    if (status/=0) return
#endif
#endif
#endif

#ifndef without_wet_deposition
    ! maximum loss:
    zwetrm(1:lmx,:) = matmul( lbdcv(1:lmx,1:lmx), rm(1:lmx,:) )
    ! actual loss depends on tracer:
    do n = 1, ntr
      zwetrm(1:lmx,n) = zwetrm(1:lmx,n) * cvsfac(n)
    end do
    ! remove from current concentrations:
    rm(1:lmx,:) = rm(1:lmx,:) - zwetrm(1:lmx,:)
    ! ?? this was not applied to slopes, why ??
#endif

    deallocate(ipiv)
    ! ok
    status = 0

  end subroutine TM5_Conv_Apply


end module TM5_Conv
