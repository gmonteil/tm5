!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!#################################################################
#include "tm5.inc"

!------------------------------------------------------------------------------
!                    TM5                                                      !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE:  convection
!
! !DESCRIPTION: methods to init and apply convection
!\\
!\\
! !INTERFACE:
!
module convection
  !
  ! !USES:
  !
  use GO, only : gol, goPr, goErr
  use misctools, only : collapse_loop, T_region_iter

  implicit none

  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public  ::  Convection_Init, Convection_Done
  public  ::  Convec
  !
  ! !PRIVATE DATA MEMBERS:
  !
  character(len=*), parameter  ::  mname = 'convection'

  type(T_region_iter), allocatable :: region_ij_counter(:) ! should be dimension(nregions)
  !
  ! no diffusion in top layers if diagonal elemens are below:
  real                               ::  eps_d
  !
  ! timers:
  integer             ::  itim_convdiff
  integer             ::  itim_conv
  integer             ::  itim_diff
  !
  !! histogram with conv/diff top:
  !integer, allocatable   ::  lmc_hist(:)
  !integer, allocatable   ::  lmd_hist(:)
  !
  ! !REVISION HISTORY:
  !   30 Mar 2010 - P. Le Sager - re-instate test on wet dep flag in convec
  !
  ! !REMARKS:
  !
  !EOP
  !------------------------------------------------------------------------------

contains



  ! ================================================================


  subroutine Convection_Init( status )

    use GO                 , only : GO_Timer_Def
    use GO                 , only : ReadRc
    use global_data        , only : rcF
    use dims               , only : nregions, isr, ier, jsr, jer
    use MeteoData          , only : mdat_set
    use MeteoData          , only : entu_dat, entd_dat, detu_dat, detd_dat

    ! --- in/out --------------------------------

    integer, intent(out)           ::  status

    ! --- const ------------------------------

    character(len=*), parameter ::  rname = mname//'/Convection_Init'

    ! --- local -------------------------------------
    integer :: i, j, Ni, Nj, region, idx, jdx

    ! --- begin --------------------------------
    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! define timers:
    call GO_Timer_Def( itim_convdiff, 'convdiff', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_conv, 'conv', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_diff, 'diff', status )
    IF_NOTOK_RETURN(status=1)

    ! loop over regions:
    do region = 1, nregions
      ! enable meteo:
      call mdat_set( entu_dat(region), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( entd_dat(region), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( detu_dat(region), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable meteo:
      call mdat_set( detd_dat(region), status, used=.true. )
      IF_NOTOK_RETURN(status=1)
    end do

    !! histograms:
    !allocate( lmc_hist(0:lmax_conv) ) ; lmc_hist = 0
    !allocate( lmd_hist(0:lmax_conv) ) ; lmd_hist = 0

    ! settings:
    call ReadRc( rcF, 'proces.diffusion.eps_d', eps_d, status )
    IF_NOTOK_RETURN(status=1)

    ! Fill the region counter
    allocate(region_ij_counter(nregions))
    do region = 1, nregions
        Ni = ier(region) - isr(region) + 1
        Nj = jer(region) - jsr(region) + 1
        allocate(region_ij_counter(region)%counter_2d(Ni*Nj, 2))
        region_ij_counter(region)%n_iter = Ni*Nj
        !call collapse_loop((/ isr(region):ier(region):1 /), (/ jsr(region):jer(region):1 /), region_ij_counter(region)%counter_2d)
        call collapse_loop( (/ (idx, idx=isr(region),ier(region),1) /), (/ (jdx, jdx=jsr(region),jer(region),1) /), region_ij_counter(region)%counter_2d )
    end do

    write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Convection_Init


  ! ***


  subroutine Convection_Done( status )

    !use TM5_Geometry, only : levi

    ! --- in/out --------------------------------

    integer, intent(out)           ::  status

    ! --- const ------------------------------

    character(len=*), parameter ::  rname = mname//'/Convection_Done'

    ! --- local --------------------------------

    !integer           :: k

    ! --- begin --------------------------------
    write (gol,'(a," : entering")') trim(rname) ; call goPr
    !! show histograms:
    !write (gol,'("lmc histogram (",i8," total)")') sum(lmc_hist); call goPr
    !do k = lbound(lmc_hist,1), ubound(lmc_hist,1)
    !  write (gol,'(i6,f8.2,i8)') k, levi%p0(k)/1e2, lmc_hist(k); call goPr
    !end do
    !write (gol,'("lmd histogram (",i8," total; eps_d = ",e8.1,")")') &
    !                sum(lmd_hist), eps_d; call goPr
    !do k = lbound(lmd_hist,1), ubound(lmd_hist,1)
    !  write (gol,'(i6,f8.2,i8)') k, levi%p0(k)/1e2, lmd_hist(k); call goPr
    !end do

    deallocate(region_ij_counter)

    write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Convection_Done



  ! ================================================================


  !-----------------------------------------------------------------------
  !
  !****   convec          - mix tracers by convection     v 8.5
  !
  !       programmed by           mh      mpi HH          1-oct-1991
  !       modified by             mh      mpi HH         11-jun-1994
  !
  !       purpose
  !       -------
  !       mix tracers vertically by convection
  !
  !       interface
  !       ---------
  !       call convec
  !
  !       method
  !       ------
  !       the matrix conv is applied simultaneously to rm, rxm and rym.
  !       The diagonal elements of conv are applied to rzm.
  !
  !       externals
  !       ---------
  !       subroutines:    tstamp
  !
  !       reference
  !       ---------
  !       see manual
  !-----------------------------------------------------------------------
  !       Following the method by Walter Guelle (1997) convective removal of
  !        tracers is added with some little changes:
  !       a vector cvsfac is added containing scavenging efficiencies (0-1)
  !       per tracer.
  !                                   aj, knmi may 1998
  !      Maarten Krol, july 2000, implemented zoom version
  !      lmax_conv is now set as the maximum layer to which convection
  !      is taken into account.
  !
  !       Second moments have been added: Bram Bregman, August 2004
  !--------------------------------------------------------------

  subroutine convec( region, step, tr, status )

    use GO           , only : GO_Timer_Start, GO_Timer_End
    use GO           , only : TDate, operator(-), rTotal
    use dims,          only : dx, dy, xref, yref, tref
    use dims,          only : isr, ier, jsr, jer
    use dims,          only : im, jm, lm
    use dims,          only : itau, okdebug, kdebug
    use dims,          only : lmax_conv
    use dims,          only : nregions
    use dims,          only : revert
    use global_data,   only : mass_dat,region_dat
    use MeteoData    , only : m_dat
    use global_data,   only : emis_data
    use MeteoData    , only : entu_dat, entd_dat, detu_dat, detd_dat
    use global_data,   only : conv_dat
    use chem_param,    only : ra, ntracet
#ifdef with_budgets
    use budget_global, only : nzon_vg, apply_budget_global
#ifndef without_wet_deposition
    use wet_deposition,only : sum_wetdep => sum_wet, buddep_dat
#endif
#endif
    use toolbox,       only : lvlpress
    use TM5_ConvDiff , only : TM5_ConvDiff_Matrix
    use TM5_ConvDiff , only : TM5_ConvDiff_Apply
    use TM5_ConvDiff , only : TM5_ConvDiff_Apply2
    use TM5_Conv     , only : TM5_Conv_Matrix, TM5_Conv_Apply
    use TM5_Diff     , only : TM5_Diff_Matrix, TM5_Diff_Apply
#ifndef without_wet_deposition
    use surface,       only : cp
    use wet_deposition,only : cvsfac
#endif
    use datetime,      only : tstamp
#ifdef MPI
    use mpi_const, only: my_real, &
         mpi_min,mpi_max,com_trac,ierr
    use mpi_const, only: pe_first_tracer
#endif
    use ParTools, only : myid, root_t
    use ParTools, only : ntracetloc, ntracet_ar
    use omp_lib

    ! input/output
    integer, intent(in)             ::  region
    character(len=1), intent(in)    ::  step
    type(TDate), intent(in)         ::  tr(2)
    integer, intent(out)            ::  status

    ! const
    character(len=*), parameter ::  rname = mname//'/Convec'

    ! local
    real,dimension(:,:,:,:), pointer  ::  rm
#ifdef slopes
    real,dimension(:,:,:,:), pointer  ::  rxm,rym,rzm
#ifdef secmom
    real,dimension(:,:,:,:), pointer  ::  rxxm,rxym,rxzm,ryym,ryzm,rzzm
#endif
#endif
    real,dimension(:,:,:)  , pointer  ::  m
    real,dimension(:,:,:)  , pointer  ::  entu, entd, detu, detd
#ifndef without_diffusion
    real,dimension(:,:,:)  , pointer  ::  dkg
#endif
    integer,dimension(:,:) , pointer  ::  cloud_base, cloud_top, cloud_lfs, zoomed
    real, allocatable                 ::  conv1(:,:), lbdcv(:,:)
    integer                           ::  imr, jmr, lmr
    integer                           ::  i, j, counter
    real                              ::  dt
    character(len=1)                  ::  trans
#ifndef without_wet_deposition
    real                              ::  sceffdeep
    real, allocatable                 ::  zwetrm2(:,:)
#ifdef with_budgets
    integer                           ::  nglob
    real                              ::  g_sum_wet(ntracet)
#endif
#endif
    integer                           ::  nfail
    ! conv+diff
    integer                           ::  lmc
    integer                           ::  lmd
    real, allocatable                 ::  dl(:), dm(:), du(:)
    real, allocatable                 ::  fd(:,:), fu(:,:), amu(:), amd(:)
    integer, allocatable              ::  ipiv(:)

    ! start
    !write (gol,'(a," : entering")') trim(rname) ; call goPr

    call GO_Timer_Start( itim_convdiff, status )
    IF_NOTOK_RETURN(status=1)

    imr=im(region) ; jmr=jm(region) ; lmr=lm(region)

    !WP! only PE's with nonzero ntracet
    if ( ntracetloc == 0 ) then
      call Go_Timer_End(itim_convdiff, status)
      return
    endif

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
#ifdef slopes
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
#ifdef secmom
    rxxm => mass_dat(region)%rxxm_t
    rxym => mass_dat(region)%rxym_t
    rxzm => mass_dat(region)%rxzm_t
    ryym => mass_dat(region)%ryym_t
    ryzm => mass_dat(region)%ryzm_t
    rzzm => mass_dat(region)%rzzm_t
#endif
#endif
    zoomed => region_dat(region)%zoomed
    entu => entu_dat(region)%data
    entd => entd_dat(region)%data
    detu => detu_dat(region)%data
    detd => detd_dat(region)%data
    cloud_top  => conv_dat(region)%cloud_top
    cloud_base => conv_dat(region)%cloud_base
    cloud_lfs  => conv_dat(region)%cloud_lfs
#ifndef without_diffusion
    dkg => conv_dat(region)%dkg
#endif

    if ( okdebug ) call tstamp(kdebug,itau,'convec ')
    if ( okdebug ) print *,'convec: Convection called (region)', region

    !timestep and position in packed arrays
    dt = abs(rTotal( tr(2) - tr(1), 'sec' ))

    ! reverse run ?
    if ( revert == -1 ) then
      trans = 'T'   ! adjoint run (transpose)
    else
      trans = 'N'   ! forward run
    end if

    ! *

    ! choose step:
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    select case ( step )
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'v' )   ! convection (+diffusion if solved together)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call GO_Timer_Start( itim_conv, status )
        IF_NOTOK_RETURN(status=1)

#ifndef without_wet_deposition
        ! some tuning:
        sceffdeep = 1.0
        ! wet removal budget:
        allocate( zwetrm2(lmax_conv,ntracetloc) )
#ifdef with_budgets
        ! initialize local variables:
        g_sum_wet = 0.0
#endif
#endif

        ! start main loop over surface cells
        ! no failures yet:
        nfail = 0
        ! PLS : re-incorporate the test on wet dep, from older cy3base/src.
        !       It hides cvsfac and cp_scale when wet dep is off.
        !$OMP PARALLEL &
        !$OMP   default (shared) &
#ifndef without_wet_deposition
        !$OMP   private ( zwetrm2, nglob ) &
#ifdef with_budgets
        !$OMP   reduction( + : g_sum_wet ) &
#endif
#endif
        !$OMP   private ( i, j, lmc, counter ) &
        !$OMP   private ( conv1, lbdcv, fd, fu, amu, amd, ipiv, status ) &
        !$OMP   reduction ( + : nfail )
        allocate(conv1(lmax_conv,lmax_conv), lbdcv(lmax_conv,lmax_conv))
        allocate(fd(0:lmax_conv,lmax_conv), fu(0:lmax_conv,lmax_conv), amu(0:lmax_conv), amd(0:lmax_conv))
        allocate(ipiv(lmax_conv))
        !$OMP   DO schedule(dynamic, 10)
        do counter = 1, region_ij_counter(region)%n_iter
            i = region_ij_counter(region)%counter_2d(counter, 1)
            j = region_ij_counter(region)%counter_2d(counter, 2)

        !do  j=jsr(region),jer(region)
           !do  i=isr(region),ier(region)

#ifdef with_zoom
              if (zoomed(i,j) /= region) cycle
#endif

#ifdef with_convdiff

              ! compute convection+diffusion matrix:
              !if (counter == 1) write(*, '(a, " : calling TM5_ConvDiff_Matrix")') rname
              call TM5_ConvDiff_Matrix( dt, lmax_conv, m(i,j,1:lmax_conv), &
                                         entu(i,j,:), detu(i,j,:), &
                                         entd(i,j,:), detd(i,j,:), &
                                         cloud_lfs(i,j), cloud_top(i,j), &
#ifndef without_wet_deposition
                                         cp_dat(region)%data(i,j,1)/cp_scale, &
                                         sceffdeep, &
#endif
#ifndef without_diffusion
                                         dkg(i,j,:), &
#endif
                                         conv1, lbdcv, status )
              if (status/=0) then
                nfail = nfail + 1
                cycle
              end if

              ! apply to tracers:
              !if (counter == 1) write(*, '(a, " : calling TM5_ConvDiff_Apply")') rname
              call TM5_ConvDiff_Apply( revert, &
                               lmax_conv, ntracetloc, conv1, &
                               rm(i,j,1:lmax_conv,:), &
#ifdef slopes
                               rxm(i,j,1:lmax_conv,:), rym(i,j,1:lmax_conv,:), rzm(i,j,1:lmax_conv,:), &
#ifdef secmom
                               rxxm(i,j,1:lmax_conv,:), ryym(i,j,1:lmax_conv,:), rzzm(i,j,1:lmax_conv,:), &
                               rxym(i,j,1:lmax_conv,:), rxzm(i,j,1:lmax_conv,:), ryzm(i,j,1:lmax_conv,:)
#endif
#endif
#ifndef without_wet_deposition
                               lbdcv, cvsfac(offsetn+1:offsetn+ntracetloc), &
                               zwetrm2, &
#endif
                               status )
              if (status/=0) then
                nfail = nfail + 1
                cycle
              end if

#else

              ! compute convection matrix:
              !if (counter == 1) write(*, '(a, " : calling TM5_Conv_Matrix")') rname
              call TM5_Conv_Matrix( dt, lmax_conv, m(i,j,1:lmax_conv), &
                                         entu(i,j,:), detu(i,j,:), &
                                         entd(i,j,:), detd(i,j,:), &
                                         cloud_lfs(i,j), cloud_top(i,j), &
#ifndef without_wet_deposition
                                         cp_dat(region)%data(i,j,1)/cp_scale, &
                                         sceffdeep, &
#endif
                                         conv1, lbdcv, lmc, &
                                         fd, fu, amu, amd, status )
              if (status/=0) then
                nfail = nfail + 1
                cycle
              end if

              !! add to histogram:
              !lmc_hist(lmd) = lmc_hist(lmc) + 1

              ! some layers with convection ?
              if ( lmc > 0 ) then

                ! apply to tracers:
                !if (counter == 1) write(*, '(a, " : calling TM5_Conv_Apply")') rname
                call TM5_Conv_Apply( trans, &
                                 lmc, ntracetloc, conv1(1:lmc,1:lmc), &
                                 rm(i,j,1:lmc,:), &
#ifdef slopes
                                 rxm(i,j,1:lmc,:), rym(i,j,1:lmc,:), rzm(i,j,1:lmc,:), &
#ifdef secmom
                                 rxxm(i,j,1:lmc,:), ryym(i,j,1:lmc,:), rzzm(i,j,1:lmc,:), &
                                 rxym(i,j,1:lmc,:), rxzm(i,j,1:lmc,:), ryzm(i,j,1:lmc,:)
#endif
#endif
#ifndef without_wet_deposition
                                 lbdcv, cvsfac(offsetn+1:offsetn+ntracetloc), &
                                 zwetrm2, &
#endif
                                 ipiv, status )
                status = 0
                if (status/=0) then
                  nfail = nfail + 1
                  cycle
                end if

              end if

#endif  ! without convdiff

#ifndef without_wet_deposition
#ifdef with_budgets
              if (apply_budget_global) then
                ! loop over local tracers
                do n = 1, ntracetloc
                  ! global tracer index:
                  nglob=n+offsetn
                  if (cvsfac(nglob).gt.0.) then
                    do l=1,lmax_conv
                      nzv=nzon_vg(l)
                      buddep_dat(region)%cp(i,j,nzv,nglob) = &
                           buddep_dat(region)%cp(i,j,nzv,nglob) + &
                           zwetrm2(l,n)/ra(nglob)*1.e3       !kg=> moles
                      g_sum_wet(nglob) = g_sum_wet(nglob) + zwetrm(l,n)
                    end do !l
                  end if  ! cvsfac > 0
                end do  ! tracers
              end if
#endif
#endif

           !end do   ! i
         !end do     ! j
         end do ! counter
        !$OMP   END DO
        deallocate(ipiv)
        deallocate(fd, fu, amu, amd)
        deallocate(conv1, lbdcv)
        !$OMP END PARALLEL

        ! failures ?
        if ( nfail > 0 ) then
          write (gol,'("matrix inversion faild for ",i6," locations")') nfail; call goErr
          TRACEBACK; status=1; return
        end if

#ifndef without_wet_deposition
#if defined (with_budgets)
        sum_wetdep(region,:) = sum_wetdep(region,:) + g_sum_wet
#endif
#endif

        call GO_Timer_End( itim_conv, status )
        IF_NOTOK_RETURN(status=1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'd' )   ! diffusion (solved seperately)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef with_convdiff
#ifndef without_diffusion

        call GO_Timer_Start( itim_diff, status )
        IF_NOTOK_RETURN(status=1)

        ! start main loop over surface cells
        ! no failures yet:
        nfail = 0
        !$OMP PARALLEL &
        !$OMP   default (none) &
        !$OMP   shared ( region ) &
        !$OMP   shared ( zoomed, isr, ier, jsr, jer, region_ij_counter ) &
        !$OMP   shared ( ntracetloc ) &
        !$OMP   shared ( dt ) &
        !$OMP   shared ( trans ) &
        !$OMP   shared ( rm ) &
#ifdef slopes
        !$OMP   shared ( rxm, rym, rzm ) &
#ifdef secmom
        !$OMP   shared ( rxxm, ryym, rzzm ) &
        !$OMP   shared ( rxym, ryzm, rxzm ) &
#endif
#endif
        !$OMP   shared ( m, dkg ) &
        !$OMP   shared ( eps_d ) &
        !$OMP   private ( i, j, lmd, counter ) &
        !$OMP   private ( dl, dm, du ) &
        !$OMP   private ( status ) &
        !$OMP   reduction ( + : nfail )
        allocate(dl(lmax_conv-1), dm(lmax_conv), du(lmax_conv-1))
        !$OMP   DO schedule(dynamic, 50)
        do counter = 1, region_ij_counter(region)%n_iter
            i = region_ij_counter(region)%counter_2d(counter, 1)
            j = region_ij_counter(region)%counter_2d(counter, 2)
        !do  j=jsr(region),jer(region)
           !do  i=isr(region),ier(region)

#ifdef with_zoom
              if (zoomed(i,j) /= region) cycle
#endif

              ! compute diffustion matrix, stored as tri-diagonal:
              !if (counter == 1) write(*, '(a, " : calling TM5_Diff_Matrix")') rname
              call TM5_Diff_Matrix( lmax_conv, m(i,j,1:lmax_conv), &
                                     dkg(i,j,:), dt, eps_d, &
                                     dl, dm, du, lmd, &
                                     status )
              if (status/=0) then
                print *, 'xxx diff matrix failed for ', i, j
                nfail = nfail + 1
                cycle
              end if

              !! add to histogram:
              !lmd_hist(lmd) = lmd_hist(lmd) + 1

              ! some layers with significant diffusion ?
              if ( lmd > 0 ) then
                ! apply to tracers:
                !if (counter == 1) write(*, '(a, " : calling TM5_Diff_Apply")') rname
                call TM5_Diff_Apply( trans, &
                             lmd, ntracetloc, &
                             dl, dm, du, &
                             rm(i,j,1:lmd,:), &
#ifdef slopes
                             rxm(i,j,1:lmd,:), rym(i,j,1:lmd,:), rzm(i,j,1:lmd,:), &
#ifdef secmom
                             rxxm(i,j,1:lmd,:), ryym(i,j,1:lmd,:), rzzm(i,j,1:lmd,:), &
                             rxym(i,j,1:lmd,:), rxzm(i,j,1:lmd,:), ryzm(i,j,1:lmd,:)
#endif
#endif
                             status )
                if (status/=0) then
                  print *, 'xxx solving diffusion failed for ', i, j
                  nfail = nfail + 1
                  cycle
                end if
              end if

           !end do   ! i
         !end do     ! j
         end do ! counter
        !$OMP   END DO
        deallocate(dl, dm, du)
        !$OMP END PARALLEL

        ! failures ?
        if ( nfail > 0 ) then
          write (gol,'("matrix inversion faild for ",i6," locations")') nfail; call goErr
          TRACEBACK; status=1; return
        end if

        call GO_Timer_End( itim_diff, status )
        IF_NOTOK_RETURN(status=1)

#endif  !  with_diffusion
#endif  !  without convdiff

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        write (gol,'("unsupported step : ",a)') trim(step); call goErr
        TRACEBACK; status=1; return

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end select
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! *


    !if (okdebug)  then
    !
    !   call mpi_allreduce(minvalue,minvalue_all,1,my_real, &
    !        mpi_min,com_trac,ierr)
    !   call mpi_allreduce(mxval,mxval_all,1,my_real,mpi_max,com_trac,ierr)
    !   call mpi_allreduce(mnval,mnval_all,1,my_real,mpi_min,com_trac,ierr)
    !   call mpi_allreduce(mxvald,mxvald_all,1,my_real,mpi_max,com_trac,ierr)
    !   call barrier_t
    !   if (myid==root_t) then
    !    print* ,'-----convection matrix information----------------------'
    !    print* ,'Minimum value of convection matrix :',minvalue_all
    !    print *,'Max/Min along the columns  (?/=1)  :',mxval_all,mnval_all
    !     !these should be one/zero to ensure mass-conservation...
    !    print *,'Maximum fractional deviation from uniform mixing ', &
    !       'ratio column:',mxvald_all
    !    print* ,'------convection matrix information---------------------'
    !   end if
    !   call barrier_t
    !   if (mxvald_all > 0.1) &
    !      call exitus(1, ' Convection matrix appears wrong')
    !end if

    nullify(m)
    nullify(rm)
#ifdef slopes
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
#ifdef secmom
    nullify(rxxm)
    nullify(rxym)
    nullify(rxzm)
    nullify(ryym)
    nullify(ryzm)
    nullify(rzzm)
#endif
#endif
    nullify(zoomed)
    nullify(entu)
    nullify(detu)
    nullify(entd)
    nullify(detd)
    nullify(cloud_top)
    nullify(cloud_base)
    nullify(cloud_lfs)
#ifndef without_diffusion
    nullify(dkg)
#endif

    call GO_Timer_End( itim_convdiff, status )
    IF_NOTOK_RETURN(status=1)

    !write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine convec


end module convection
