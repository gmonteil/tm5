!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!#################################################################

module TM5_Diff

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

#ifndef without_diffusion
  public  ::  TM5_Diff_Matrix
  public  ::  TM5_Diff_Apply
#endif


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'TM5_Diff'


#ifndef without_diffusion
contains


  ! ================================================================


  !pure
  pure subroutine TM5_Diff_Matrix( lmx, m, dkg, dt, eps, &
                                  dl, dm, du, lmd, status )

    ! --- in/out ---------------------------------

    integer, intent(in)         ::  lmx  ! lmax_conv
    real, intent(in)            ::  m(lmx)
    real, intent(in)            ::  dkg(lmx)
    real, intent(in)            ::  dt
    real, intent(in)            ::  eps          ! significance treshold
    real, intent(out)           ::  dl(lmx-1)    ! lower diagonal
    real, intent(out)           ::  dm(lmx)      ! main diagonal
    real, intent(out)           ::  du(lmx-1)    ! upper diagonal
    integer, intent(out)        ::  lmd          ! layers with significant diffusion
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/TM5_Diff_Matrix'

    ! --- local ----------------------------------

    real                      ::  fd(0:lmx)
    real                      ::  fu(0:lmx)
    integer                   ::  k

    ! --- begin ----------------------------------
    ! input 'dkg(0:lmx)' : air mass flux through interface k
    ! between layer k and k+1 due to diffusion

    ! cut-off highest rows/columns of diffusion matrix
    ! if diffusion flux is small ; init top level to zero:
    lmd = 0

    ! fu(k) : fraction of tracer mass in layer k
    !  transported upwards to layer k+1 through interface k
    ! no diffusion through bottom:
    fu(0) = 0.0
    ! loop over internal interfaces:
    do k = 1, lmx-1
      ! fraction of tracer mass k transported upwards to layer k+1:
      fu(k) = dt * dkg(k) / m(k)
      ! include next layer if flux is significant:
      if ( abs(fu(k)) >= eps ) lmd = k+1
    end do
    ! no diffusion through top:
    fu(lmx) = 0.0

    ! fd(k) : fraction of tracer mass in layer k+1
    !  transported downwards to layer k through interface k
    ! no diffusion through bottom:
    fd(0) = 0.0
    ! loop over internal interfaces:
    do k = 1, lmx-1
      ! fraction of tracer mass in k+1 transported downwards to layer k:
      fd(k) = dt * dkg(k) / m(k+1)
      ! include next layer if flux is significant:
      if ( abs(fd(k)) >= eps ) lmd = k+1
    end do
    fd(lmx) = 0.0

    ! generate backwards discretization matrix:
    !   A c(t+1) = c(t)
    ! lmd could be zero (no signficant diffusion)
    ! or >= 2 (exchange between at least 2 layers):
    if ( lmd > 0 ) then
      dm(1:lmd) = 1.0
      do k = 1, lmd-1
        dl(k)   =       - fu(k)
        dm(k)   = dm(k) + fu(k)
      end do
      do k = 2, lmd
        du(k-1) =       - fd(k-1)
        dm(k  ) = dm(k) + fd(k-1)
      end do
    end if

    !! check ...
    !if ( lmd > 0 ) then
    !  if ( abs(sum(dl(1:lmd-1))+sum(dm(1:lmd))+sum(du(1:lmd-1)) - lmd) > 1e-9 ) then
    !    print *, 'diffusion matrix not ok; column sum should be 1.0 :'
    !    print *, 0.0, dm(1), dl(1), dm(1)+dl(1)
    !    do k = 2, lmd-1
    !      print *, du(k-1), dm(k), dl(k), du(k-1)+dm(k)+dl(k)
    !    end do
    !    print *, du(lmd-1), dm(lmd), 0.0, du(lmd-1)+dm(lmd)
    !    print *, 'total sum : ', abs(sum(dl(1:lmd-1))+sum(dm(1:lmd))+sum(du(1:lmd-1)) - lmd)
    !    stop
    !  end if
    !end if
    ! ok
    status = 0

  end subroutine TM5_Diff_Matrix


  ! *


  subroutine TM5_Diff_Apply( trans, &
                                lmx, ntr, &
                                dl, dm, du, &
                                rm, &
#ifdef slopes
                                rxm, rym, rzm, &
#ifdef secmom
                                rxxm, ryym, rzzm, rxym, rxzm, ryzm
#endif
#endif
                                status )

    !use omp_lib

    ! --- in/out ---------------------------------

    character(len=1), intent(in)  ::  trans ! 'N' for forward, and
                                            ! 'T' for reverse run
    integer, intent(in)           ::  lmx   ! lmax_conv
    integer, intent(in)           ::  ntr
    real, intent(inout)           ::  dl(lmx-1)
    real, intent(inout)           ::  dm(lmx  )
    real, intent(inout)           ::  du(lmx-1)
    real, intent(inout)           ::  rm(lmx,ntr)
#ifdef slopes
    real, intent(inout)           ::  rxm(lmx,ntr), rym(lmx,ntr), rzm(lmx,ntr)
#ifdef secmom
    real, intent(inout)           ::  rxxm(lm,ntr), ryym(lm,ntr), rzzm(lm,ntr)
    real, intent(inout)           ::  rxym(lm,ntr), rxzm(lm,ntr), ryzm(lm,ntr)
#endif
#endif
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/TM5_Diff_Apply'

    ! --- local ----------------------------------

    real              :: du2(lmx-2)
    integer           :: ipiv(lmx), itr, info
    real              :: c(lmx), d(lmx), e(lmx), f(lmx)
    logical           :: had_error

    ! --- begin ----------------------------------

    !write(*, '(a, ": thread=", i2, ", du2=", i15, ", ipiv=", i15, ", c=", i15, ", d=", i15, ", e=", i15, ", f=", i15)') &
        !rname, omp_get_thread_num(), loc(du2), loc(ipiv), loc(c), loc(d), loc(e), loc(f)

    d = dm
    ! SB: use ESSL for solving instead of LAPACK
#ifdef with_essl
    if (trans == 'N') then
        c(2:lmx) = dl
        e(1:lmx-1) = du
    else if (trans == 'T') then
        c(2:lmx) = du
        e(1:lmx-1) = dl
    end if

    call DGTF(lmx, c, d, e, f, ipiv)
    do itr = 1, ntr
        call DGTS(lmx, c, d, e, f, ipiv, rm(:,itr))
#ifdef slopes
        call DGTS(lmx, c, d, e, f, ipiv, rxm(:,itr))
        call DGTS(lmx, c, d, e, f, ipiv, rym(:,itr))
        call DGTS(lmx, c, d, e, f, ipiv, rzm(:,itr))
#endif
#ifdef secmom
        call DGTS(lmx, c, d, e, f, ipiv, rxxm(:,itr))
        call DGTS(lmx, c, d, e, f, ipiv, ryym(:,itr))
        call DGTS(lmx, c, d, e, f, ipiv, rzzm(:,itr))
        call DGTS(lmx, c, d, e, f, ipiv, rxym(:,itr))
        call DGTS(lmx, c, d, e, f, ipiv, rxzm(:,itr))
        call DGTS(lmx, c, d, e, f, ipiv, ryzm(:,itr))
#endif
    end do

#endif
#ifdef with_mkl

    !write(0, *) "XXXX DESTRUCTION AND HORRIBLE DEATH FOLLOWS"
    call dgttrf(lmx, dl, d, du, du2, ipiv, info)

    if (info /= 0) then
       write(0,'("tridiag factor failed. dgttrf info = ",i3)') info
    end if

    had_error = .false.
    do itr = 1, ntr
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rm(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
#ifdef slopes
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rxm(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rym(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rzm(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
#endif
#ifdef secmom
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rxxm(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, ryym(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rzzm(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rxym(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, rxzm(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
        call dgttrs(trans, lmx, 1, dl, d, du, du2, ipiv, ryzm(:,itr), lmx, info)
        had_error = had_error .or. (info /= 0)
#endif
    end do
    if (had_error) then
       write(0,*) "tridiag solve had an error."
    end if
#endif

    ! LU decomposition (general tri-diagnonal matrix, triangular factorization)
    ! NOTE: apply same factors to ALL slopes ;
    !   in previous code z-slopes were omitted for not
    !   further described reasons
!    call dGtTrf( lmx, dl, dm, du, du2, ipiv, status )
!    if (status/=0) return

!    ! Solve for all tracers:
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rm, lmx, status )
!    if (status/=0) return
!#ifdef slopes
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rxm, lmx, status )
!    if (status/=0) return
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rym, lmx, status )
!    if (status/=0) return
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rzm, lmx, status )
!    if (status/=0) return
!#ifdef secmom
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rxxm, lmx, status )
!    if (status/=0) return
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, ryym, lmx, status )
!    if (status/=0) return
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rzzm, lmx, status )
!    if (status/=0) return
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rxym, lmx, status )
!    if (status/=0) return
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, rxzm, lmx, status )
!    if (status/=0) return
!    call dGtTrs( trans, lmx, ntr, dl, dm, du, du2, ipiv, ryzm, lmx, status )
!    if (status/=0) return
!#endif
!#endif
    status = 0

  end subroutine TM5_Diff_Apply

#endif  // diffusion


end module TM5_Diff
