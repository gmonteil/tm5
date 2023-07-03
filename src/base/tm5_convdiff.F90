!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!#################################################################

module TM5_ConvDiff

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  TM5_ConvDiff_Matrix
  public  ::  TM5_ConvDiff_Apply
  public  ::  TM5_ConvDiff_Apply2


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'TM5_ConvDiff'


contains


  ! ================================================================


  pure subroutine TM5_ConvDiff_Matrix( &
                        dt, lmx, m, &
                        entu, detu, entd, detd, &
                        ld, li, &
#ifndef without_wet_deposition
                        cp, sceffdeep, &
#endif
#ifndef without_diffusion
                        dkg, &
#endif
                        conv1, lbdcv, &
                        status )

    ! --- in/out ---------------------------------

    real, intent(in)            ::  dt
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
#ifndef without_diffusion
    real, intent(in)            ::  dkg(lmx)
#endif
    real, intent(out)           ::  conv1(lmx,lmx)    ! convection matrix
    real, intent(out)           ::  lbdcv(lmx,lmx)    ! wet removal matrix
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/TM5_ConvDiff_Matrix'

    ! --- local ----------------------------------

    ! Since this is a pure subroutine, we should no longer use static arrays
    ! Instead these arrays will have to be allocatable
    ! tracer fluxes:
    real, allocatable         ::  f(:,:)
    real, allocatable         ::  fu(:,:)
    ! air mass up/downdraft:
    real, allocatable         ::  amu(:)
    real, allocatable         ::  amd(:)

    integer                   ::  k, kk
    integer                   ::  l
    real                      ::  zxi

    ! --- begin ----------------------------------
    allocate(f(0:lmx,lmx), fu(0:lmx,lmx))
    allocate(amu(0:lmx), amd(0:lmx))

    ! WARNING: For matrices f(k,kk), fu(k,kk), and fd(k,kk), index k
    !          corresponds to the upper boundary of level k
    !          and index j to level j.
    !          For matrices g(i,j,k,kk) and lbdcv(i,j,k,kk),
    !          both indices k and j
    !          correspond to respective levels and not to boundaries.

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
    do k=1,lmx-1
       do kk=1,lmx
          f(k,kk)=fu(k,kk)+f(k,kk)
       end do !kk
       ! add coefficient from subsidence
       !CMK SH       f(k,k+1)=f(k,k+1)-(amu(k)+amd(k))/m(k+1)
       ! reported by Sander Houweling.........
       f(k,k+1)=f(k,k+1)-amu(k)/m(k+1)
       f(k,k  )=f(k,k  )-amd(k)/m(k  )
    end do

#ifndef without_diffusion
    ! diffusion terms
    ! no diffusion in the stratosphere  k>lmx?
    do k = 1, lmx-1
      f(k,k  ) = f(k,k  ) + dkg(k)/m(k  )
      f(k,k+1) = f(k,k+1) - dkg(k)/m(k+1)
    end do
#endif

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

    ! generate forward matrix
    conv1 = 0.0
    do k=1,lmx
       do kk=1,lmx
          conv1(k,kk) = -dt * (f(k-1,kk)-f(k,kk))
       end do !kk
       conv1(k,k) = conv1(k,k) + 1.0
    end do !k

    deallocate(f, fu, amd, amu)
    ! ok
    status = 0

  end subroutine TM5_ConvDiff_Matrix


  ! *

  subroutine TM5_ConvDiff_Apply( revert, &
                           lmx, ntr, conv1, &
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

    use linalg_interface, only: Inverse
    ! --- in/out ---------------------------------

    integer, intent(in)         ::  revert    ! 1 (forward), -1 (adjoint)
    integer, intent(in)         ::  lmx       ! lmax_conv
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

    character(len=*), parameter ::  rname = mname//'/TM5_ConvDiff_Apply'

    ! --- local ----------------------------------

    real,dimension(lmx)               :: zrm
    real,dimension(lmx)               :: srm
#ifdef slopes
    real,dimension(lmx)               :: zrxm,zrym
    real,dimension(lmx)               :: srxm,srym,srzm
#ifdef secmom
    real,dimension(lmx)               :: zrxxm,zrxym,zrxzm,zryym,zryzm
    real,dimension(lmx)               :: srxxm,srxym,srxzm,sryym,sryzm,srzzm
#endif
#endif
    integer                   ::  k
    integer                   ::  l
    integer                   ::  n
    real                      ::  scavf

    ! --- begin ----------------------------------
    ! matrix inverse:
    call Inverse(lmx, conv1)

    ! adjoint run ? then transpose:
    if (revert==-1) conv1 = transpose(conv1)

    ! loop over local tracers
    do n = 1, ntr

      ! copy tracer concentrations for this column in srm/srxm/etc:
      srm(1:lmx) = rm(1:lmx,n)
#ifdef slopes
      srxm(1:lmx) = rxm(1:lmx,n)
      srym(1:lmx) = rym(1:lmx,n)
      srzm(1:lmx) = rzm(1:lmx,n)
#ifdef secmom
      srxxm(1:lmx) = rxxm(1:lmx,n)
      srxym(1:lmx) = rxym(1:lmx,n)
      srxzm(1:lmx) = rxzm(1:lmx,n)
      sryym(1:lmx) = ryym(1:lmx,n)
      sryzm(1:lmx) = ryzm(1:lmx,n)
      srzzm(1:lmx) = rzzm(1:lmx,n)
#endif
#endif

      ! init budgets to zero:
#ifndef without_wet_deposition
      zwetrm(1:lmx,n) = 0.0
#endif
      !zcnvd (1:lmx,n) = 0.0

      ! init column of concentrations to zero:
      zrm(1:lmx )= 0.0
#ifdef slopes
      zrxm(1:lmx) = 0.0
      zrym(1:lmx) = 0.0
#ifdef secmom
      zrxxm(1:lmx) = 0.0
      zrxym(1:lmx) = 0.0
      zrxzm(1:lmx) = 0.0
      zryym(1:lmx) = 0.0
      zryzm(1:lmx) = 0.0
#endif
#endif

      ! apply convection matrix:
      do l = 1, lmx
#ifdef slopes
        srzm(l) = conv1(l,l)*srzm(l)   !rzm is directly calculated
#ifdef secmom
        srzzm(l) = conv1(l,l)*srzzm(l)   !rzzm is directly calculated
#endif
#endif
        do k = 1, lmx

          ! species dependent removal....
#ifndef without_wet_deposition
          scavf = min( lbdcv(l,k)*cvsfac(n), conv1(l,k) )
#else
          scavf = 0.0
#endif

          ! for budgets:
          !zcnvd (l,n) = zcnvd (l,n) +  conv1(l,k)        * srm(k)
#ifndef without_wet_deposition
          zwetrm(l,n) = zwetrm(l,n) +             scavf  * srm(k)  ! loss is positive
#endif

          ! store new concentrations in zrm/zrxm/etc
          !  using old concentrations in srm/srxm/etc
          zrm  (l) = zrm  (l) + (conv1(l,k)-scavf) * srm(k)
#ifdef slopes
          zrxm (l) = zrxm (l) + (conv1(l,k)-scavf) * srxm(k)
          zrym (l) = zrym (l) + (conv1(l,k)-scavf) * srym(k)
#ifdef secmom
          zrxxm(l) = zrxxm(l) + (conv1(l,k)-scavf) * srxxm(k)
          zrxym(l) = zrxym(l) + (conv1(l,k)-scavf) * srxym(k)
          zrxzm(l) = zrxzm(l) + (conv1(l,k)-scavf) * srxzm(k)
          zryym(l) = zryym(l) + (conv1(l,k)-scavf) * sryym(k)
          zryzm(l) = zryzm(l) + (conv1(l,k)-scavf) * sryzm(k)
#endif
#endif
        end do    !k

        ! substract original tracer mass
        ! VH, TvN (2007-11-16)
        !zcnvd(l,n) = zcnvd(l,n) - srm(l)

      end do      !l

      ! restore concentrations:
      rm(1:lmx,n) = zrm(1:lmx)
#ifdef slopes
      rxm(1:lmx,n) = zrxm(1:lmx)
      rym(1:lmx,n) = zrym(1:lmx)
      rzm(1:lmx,n) = srzm(1:lmx)
#ifdef secmom
      rxxm(1:lmx,n) = zrxxm(1:lmx)
      rxym(1:lmx,n) = zrxym(1:lmx)
      rxzm(1:lmx,n) = srxzm(1:lmx)
      ryym(1:lmx,n) = zryym(1:lmx)
      ryzm(1:lmx,n) = zryzm(1:lmx)
      rzzm(1:lmx,n) = srzzm(1:lmx)
#endif
#endif

    end do  ! tracers
    ! ok
    status = 0

  end subroutine TM5_ConvDiff_Apply


  subroutine TM5_ConvDiff_Apply2( trans, &
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

!#ifdef with_mkl
!    ! SKIP: something wrong with interfaces;
!    ! seem to require 'END SUBROUTINE CBDSQR' etc.
!    !include "mkl_lapack.inc"
!#endif

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

    character(len=*), parameter ::  rname = mname//'/TM5_ConvDiff_Apply2'

    ! --- local ----------------------------------

    integer           ::  ipiv(lmx)
#ifndef without_wet_deposition
    integer           ::  n
#endif

    ! --- begin ----------------------------------
    ! LU decomposition (general matrix, triangular factorization)
    call dGeTrf( lmx, lmx, conv1, lmx, ipiv, status )
    if (status/=0) return

    ! NOTE: apply same factors to ALL slopes ;
    !   in previous code z-slopes were omitted for not
    !   further described reasons

    ! Solve for all tracers:

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
    ! ok
    status = 0

  end subroutine TM5_ConvDiff_Apply2


end module TM5_ConvDiff
