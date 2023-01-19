!###############################################################################
!
! version March 2003, adapted for TM5 MK
! This module contains all subroutines dealing with wet removal 
! by large scale precipitation
! adapted from the KNMI version.
! 
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################

module wet_deposition

!  use dims,      only: nregions
!#ifdef MPI  
!  use mpi_const, only: root, myid
!#endif
!  use chem_param, only: ih2o2,ihno3,ich3o2h,ich2o,irooh, &
!       iorgntr,iso2,inh3,iso4,inh4,imsa,&
!       ino3_a,ipb210,d3_data,names
  use chem_param, only : ntracet

  implicit none

  ! --- in/out --------------------------------

  private

  public :: cvsfac
!  public :: declare_rloss, free_rloss 
!  public :: calc_cvsfac, calculate_lsp_scav, lspscav

  ! --- var ------------------------------------

  !
  ! scavenging efficiencies, used in convection
  !
  real, dimension(ntracet)   :: cvsfac=0.0

!  ! private data
!  !
!  ! large scale scavenging coefficients [s-1]
!  ! 1: incloud completely soluble gas 
!  ! 2: below cloud completely soluble gas
!  ! 3: below cloud accumulation range aerosol
!  !
!  type(d3_data),dimension(nregions) :: rloss1
!  type(d3_data),dimension(nregions) :: rloss2
!  type(d3_data),dimension(nregions) :: rloss3
!
!  !
!  ! rain-out can not be higher than maximum level of convection
!  ! thus maximum level of convection lmax_conv(=>ntot_ed) is used
!  !
!  ! nscav       : selected species for scavenging
!  ! nscav_index : index for scavenging:
!  ! nscav_type  : type of scavenging:
!  !               0 no scavenging
!  !               1 scavenging 100 % solubility assumed
!  !               2 scavenging henry solubility assumed
!  !               3 scavenging, aerosol removal assumed
!  !
!  integer,parameter                    :: nscav=13
!
!  integer,parameter,dimension(nscav)   :: &
!       nscav_index  = (/ih2o2,   ihno3,  ich3o2h, ich2o, &
!       irooh, iorgntr,     iso2,  inh3, &
!       iso4,    inh4,     imsa        ,&
!       ino3_a,  ipb210                /)
!
!  ! note CMK CFD: wetS accounts for removal of SO2, so nscav_type = 0
!  !               nh3 is taken out anyhow, since most will be scavenged
!  !               by acidic falling frops below cloud
!  integer,parameter,dimension(nscav)   :: &
!       nscav_type   = (/    2,       1,        2,     2, &
!       2,       2,        0,     1, &    !corrected CMK CFD
!       3,       3,        3,        &
!       3,       3                 /)
!  !
!  !----------------------------------------------
!  ! acidity needed for explicit calculation of reactive removal of SO2.
!  ! Parameterisation calculates reaction of SO2 with H2O2 and O3.
!  ! Not yet implemented: for information Frank Dentener    ! see routine wetS
!  ! nacid       : selected species for acidity
!  ! nacid_index : indices of species for acidity : iso4,imsa,ihno3,inh3,inh4
!  ! nacid_eq    : equivalents acid per mole
!  ! integer,parameter                    :: nacid=5
!  ! integer,parameter,dimension(nacid)   :: &
!  !                          nacid_index    = (/iso4,imsa,ihno3,ino3_a,inh3,inh4/)
!  ! integer,parameter,dimension(nacid)   :: &
!  !                             nacid_eq    = (/   2,   1,    1,  1   ,  -1,  -1/)
!  !----------------------------------------------
!  !
!
!contains
!  !
!  ! subroutine lspscav, do_lspscav, declare_rloss, free_rloss
!  ! subroutine calculate_lsp_scav,do_calculate_lsp_scav
!  ! subroutine calfk(prf)
!  !
!
!  subroutine declare_rloss
!    !
!    use dims, only: im,jm,lmax_conv
!    implicit none
!    ! local
!    integer      :: region, imr,jmr,lmr
!
!    do region = 1, nregions
!       imr = im(region) ; jmr = jm(region) ; lmr = lmax_conv
!       allocate(rloss1(region)%d3(imr,jmr,lmr))
!       allocate(rloss2(region)%d3(imr,jmr,lmr))
!       allocate(rloss3(region)%d3(imr,jmr,lmr))
!    end do
!
!  end subroutine declare_rloss
!
!
!
!  subroutine free_rloss
!    implicit none
!    ! local
!    integer :: region
!
!    do region = 1, nregions
!       deallocate(rloss1(region)%d3)
!       deallocate(rloss2(region)%d3)
!       deallocate(rloss3(region)%d3)
!    end do
!
!  end subroutine free_rloss
!
!
!
!  subroutine calc_cvsfac
!    !
!    use chem_param, only: henry, ntlow 
!    ! lookup tables, calculated before by routine rates...
!    implicit none
!    ! local
!    integer :: iscav,n,k
!    real    :: rtl
!    !
!    ! calculate scale factor relative to 100% scavenging (of HNO3)
!    ! assume constant temperature in convective updraft of 273K
!    ! factors for different scavenging types are:
!    ! 0) no/low solubility  cvsfac=0
!    ! 1) high solubility    cvsfac=1
!    ! 2) henry solubility   cvsfac=variable
!    ! 3) aerosol            cvsfac=1
!    !
!    cvsfac=0.0
!    do iscav=1,nscav
!       n=nscav_index(iscav)
!       select case(nscav_type(iscav))
!       case(0)
!          cvsfac(n) = 0.0
!       case(1,3)
!          cvsfac(n) = 1.0
!       case(2)
!          rtl=8.3148e-8*273.
!          k=nint(273.-float(ntlow))
!          if ( henry(n,k) > 10. ) then
!             cvsfac(n) =  rtl*henry(n,k)/(1.+rtl*henry(n,k))
!          else
!             cvsfac(n) = 0.0
!          end if
!       case default
!          cvsfac(n) = 0.0
!       end select
!    end do
!
!#ifdef MPI
!    if(myid == root) then
!#endif
!       do n=1,ntracet
!          print *, 'calc_cvsfac: Scavenging factor species ', &
!               n, names(n), cvsfac(n)
!       end do
!#ifdef MPI
!    end if
!#endif
!    !
!  end subroutine calc_cvsfac
!
!
!
!  subroutine lspscav(region)
!    !----------------------------------------------------------------------
!    ! Calculation of wet removal by large scale precipitation of soluble tracers
!    !
!    ! Remove tracers with previously calculated rainout rate [s-1]
!    ! separately for in- and below cloud scavenging and only for the
!    ! cloud covered fraction of the gridcell
!    !      
!    ! Reference:
!    !      Langner and Rodhe (1990)
!    !      Junge (1959)
!    !
!    ! Programmed by Frank Dentener, August 1995;
!    ! Ad Jeuken, KNMI, November 1998
!    ! Modifications Bas Henzing, KNMI, 2002
!    ! Adapted for TM5, Frank Dentener, JRC, 2002
!    ! Paralel, Maarten Krol, Jul 2003
!    !----------------------------------------------------------------------------
!    use dims,          only: im, jm, lm, isr, ier, jsr, jer, nchem
!    use dims,          only: tref, lmax_conv, adv_scheme
!    use binas,         only: rgas
!    use chem_param,    only: ntrace, ntracet, henry, ntlow, ra
!    use budget_global, only: buddep_dat, nzon_vg, sum_wet
!    use global_data,   only: mass_dat, region_dat, cloud_dat, meteo_dat
!#ifdef MPI  
!    use mpi_const, only: tracer_active, tracer_loc
!#endif
!    implicit none
!
!    ! input/output
!    integer,intent(in)  :: region
!
!    ! local
!    real,dimension(:,:,:,:), pointer            :: rm, rxm, rym, rzm
!    integer,dimension(:,:), pointer             :: zoomed
!    real,dimension(:,:,:), pointer              :: t,cc
!
!    real,parameter     ::  rtl_fac=rgas/1e2*1e-6 
!    !  rgas (8.314 J/mol/K) ---> 0.08314 atm/(mol/l)/K
!    !  (thesis Frank Dentener, p. 31)
!    !  1e-6 corresponds to 1 g/m3 dimensionless 
!    !  liquid water content of raining cloud
!
!    real               ::  rtl,f,f1,f2,vol1,vol2,vol3,ahelp1,ahelp2
!    real               ::  incloud,belowcloud,newcov,c_old,corr_diff,fnchem
!    integer            ::  n,iscav,i,j,k,itemp,nzone_v, nloc
!
!    ! oldcov: cloud cover in layer above, to calculate below-cloud scaveging.
!    real,dimension(:,:),allocatable    :: oldcov
!    !
!    rm => mass_dat(region)%rm_t   ! paralel over tracers
!    rxm => mass_dat(region)%rxm_t
!    rym => mass_dat(region)%rym_t
!    rzm => mass_dat(region)%rzm_t
!    t => meteo_dat(region)%t
!    cc => cloud_dat(region)%cc
!    zoomed => region_dat(region)%zoomed
!
!    fnchem=real(nchem/(2*tref(region)))
!    !
!    allocate(oldcov(im(region),jm(region)))
!
!    do iscav=1,nscav
!       !
!       n=nscav_index(iscav)                ! tracer number in global count
!#ifdef MPI     
!       if(.not.tracer_active(n)) cycle     ! process only trcaers on this PE
!       if(.not.tracer_active(n)) cycle     ! process only trcaers on this PE
!       nloc = tracer_loc(n)                ! tracer number in local count  
!#else
!       nloc = n                            ! tracer number in local count  
!#endif
!       oldcov=0.
!       !
!       ! assumption no stratiform precipitation above the maximum 
!       ! level of convection
!       !
!       do k=lmax_conv,1,-1   ! top to bottom
!          do j=jsr(region),jer(region)
!             do i=isr(region),ier(region)
!                if(zoomed(i,j)/=region) cycle
!                !
!                ! calculate, depending on solubility, scale factor relative 
!                ! to 100% scavenging (of HNO3)
!                !
!                ! rtl - composite factor of liquid water content, rgas 
!                !       and liquid water content
!                rtl=rtl_fac*t(i,j,k)
!                ! multiplaction with Henry constant gives phase factor
!                itemp=nint(t(i,j,k)-float(ntlow))
!                rtl=rtl*henry(n,itemp)/ (1.+rtl*henry(n,itemp))
!                !
!                corr_diff=sqrt(ra(ihno3)/ra(n))   
!                ! 
!                select case (nscav_type(iscav))
!                case(0)
!                   belowcloud=0.0   ! corrected CMK CFD
!                   incloud = 0.0
!                case(1)
!                   belowcloud=rloss2(region)%d3(i,j,k)*corr_diff
!                   incloud=rloss1(region)%d3(i,j,k)*rtl
!                case(2) 
!                   belowcloud=rloss2(region)%d3(i,j,k)*rtl*corr_diff
!                   incloud=rloss1(region)%d3(i,j,k)*rtl
!                case(3) 
!                   belowcloud=rloss3(region)%d3(i,j,k)
!                   incloud=rloss1(region)%d3(i,j,k)*rtl
!                case default
!                   belowcloud=0.0   ! corrected CMK CFD
!                   incloud = 0.0
!                end select
!
!                !if(incloud > 1e-4) then
!                !print *, i,j,k,names(n),rtl, rloss1(region)%d3(i,j,k), rtl_fac
!                !end if
!                !
!                ! f1, f2 are the fractions of the tracermass that remain in the 
!                ! gridbox after scavenging.
!                !
!                f1=exp(-fnchem*incloud)
!                f2=exp(-fnchem*belowcloud)
!                !
!                ! A grid box can be divided into three volumes
!                ! 1) incloud scavenging       (newcov)
!                ! 2) below cloud scavenging   (oldcov-newcov)
!                ! 3) no in-cloud scavenging and no below cloud 
!                !    scavenging by precipitation (no removal)
!                !
!                newcov=cc(i,j,k)
!                vol1 = newcov
!                vol2 = max(0.,oldcov(i,j)-newcov)
!                vol3 = 1.-vol1-vol2
!                f=f1*vol1+f2*vol2+vol3
!                c_old=rm(i,j,k,nloc)
!                rm(i,j,k,nloc)=c_old*f
!                if (adv_scheme=='slope') then
!                   rxm(i,j,k,nloc)=rxm(i,j,k,nloc)*f
!                   rym(i,j,k,nloc)=rym(i,j,k,nloc)*f
!                   rzm(i,j,k,nloc)=rzm(i,j,k,nloc)*f
!                end if
!                ! _____budget
!                nzone_v = nzon_vg(k)
!                buddep_dat(region)%lsp(i,j,nzone_v,n)= &
!                     buddep_dat(region)%lsp(i,j,nzone_v,n)+ &
!                     (c_old-rm(i,j,k,nloc))/ra(n)*1000.   ! in moles
!                if ( n == 1 ) sum_wet(region) = sum_wet(region) + &
!                     (c_old-rm(i,j,k,nloc))   ! in kg
!                if ( region == nregions ) then
!                   nzone_v = nzon_v(k)
!                   depwet_lsp(i,j,nzone_v,n) = depwet_lsp(i,j,nzone_v,n) + &
!                        (c_old-rm(i,j,k,nloc))/ra(n)*1000. 
!                end if
!                ! _____budget
!                if ( rloss1(region)%d3(i,j,k) > 0.0 )  &
!                     oldcov(i,j)=max(oldcov(i,j),cc(i,j,k)) 
!             end do !i
!          end do !j
!       end do !k
!    end do !iscav
!
!    deallocate(oldcov)
!
!    nullify(rm)
!    nullify(rxm)
!    nullify(rym)
!    nullify(rzm)
!    nullify(t)
!    nullify(cc)
!    nullify(zoomed)
!
!    !
!  end subroutine lspscav
!
!
!
!  subroutine calculate_lsp_scav(region) 
!    !----------------------------------------------------------------------------
!    ! Calculate wet removal rates rloss1,rloss2,rloss3 (s-1) for 
!    ! stratiform precipitation from cloud-top to ground, 
!    ! distinguishing between below-cloud and in-cloud scavenging. 
!    !
!    ! Method: 
!    !        adapted from GJ Roelofs and Lelieveld, JGR, October 1995
!    !
!    ! fills array "rloss" should be called once new precipitation fields 
!    ! are available (MK: in trace_after_read)
!    !
!    ! Programmed by:  Frank Dentener, IMAU, 1996
!    ! Ad Jeuken, KNMI, 1998
!    ! Modification, Bas Henzing, KNMI, 2002.
!    ! Adapted for TM5, August 2002, Frank Dentener, JRC.
!    ! And finally for the new version, MK, IMAU, March 2003
!    ! Parallel, MK Jul 2003
!    !----------------------------------------------------------------------------
!    use binas,       only: grav, rgas, xmair
!    use dims,        only: im,jm,lm,lmax_conv,isr,jsr,ier,jer
!    use surface,     only: lspr=>lsp
!    use global_data, only: emis_data
!    use global_data, only: mass_dat, meteo_dat, cloud_dat, use_clouds
!    use toolbox,      only : escape_tm
!#ifdef MPI
!    use mpi_const, only: ntracetloc, myid, root
!#endif
!    implicit none
!
!    ! input
!    integer, intent(in)              :: region
!
!    ! local
!    real,dimension(:,:)  ,pointer    :: lsp
!    real,dimension(:,:,:),pointer    :: t, lwc, iwc, cc, phlb
!    real,parameter                   :: max_lwc=2.e-3        ! kg/m3 
!    !
!    ! Microphysical parameters:
!    !
!    !  rdrad2:  square of raindroplet radius (20 microns)
!    !  dghno3: in [cm2/s] 
!    !  dgair:  in [cm2/s] 
!    !  rol:  density of water in [kg/m3]
!    !  roi:  density of ice in [kg/m3]
!    !
!    real,parameter :: rdrad2 = (2E-5)**2
!    real,parameter :: dghno3 = 0.136
!    real,parameter :: dgair  = 0.133
!    real,parameter :: rol    = 1000.
!    real,parameter :: roi    = 917.
!    !
!    ! quantity used in calculation of Sherwood number
!    !
!    ! bas henzing: bug repair  real,parameter :: znsc=dgair/dghno3**(-3)  
!    ! bas henzing: should be: znsc=(dgair/dghno3)**(1./3.); 
!    !    znsc is now defined a real
!    !
!    real,dimension(:),allocatable   :: dzk
!    real :: rflx,beta,fac,beta1,beta2,beta3,rlwc,rdrad,rn,ru
!    real :: press,aird,dgairx,dghno3x
!    real :: znre,znsh,znsc,zkg,totw,sfz,sfz_no
!    !
!    integer :: nfz,i,j,k,l,no_fz
!    integer,dimension(:,:),allocatable :: kcltop
!    real,dimension(:,:),allocatable    :: oldcov,fz
!
!    real,dimension(:,:,:),allocatable :: precip ! precipitation per level (kg/m2/s)
!    real,dimension(:,:,:),allocatable :: prf ! precipitation formation rate.
!    ! evaporation NOT YET AVAILABLE
!    !
!    ! how much less efficient is tracer scavenged from ice
!    ! cloud droplet compared to water cloud droplet. 
!    ! This should be tracer dependent. 
!    !
!    real,parameter  :: ice_eff=0.2
!    real            :: inc_rdf
!    real,parameter  :: scale_heigth= rgas/grav/xmair*1e3 ! about 29.3*T (m)
!
!    ! start
!
!#ifdef MPI
!    if(ntracetloc == 0) return      !if no tracer on this processor do nothing...
!#endif
!
!    lsp => lspr(region)%surf
!    t => meteo_dat(region)%t
!    if( .not. use_clouds ) then
!       call escape_tm('calculate_lsp_scav: use_clouds = .false.')
!    end if
!    lwc => cloud_dat(region)%lwc   !mk: lm, and not lmax_conv
!    iwc => cloud_dat(region)%iwc   !mk: lm, and not lmax_conv
!    cc => cloud_dat(region)%cc     !mk: lm, and not lmax_conv
!    phlb => mass_dat(region)%phlb_t  !mk: 1:lm+1
!
!    allocate(oldcov(im(region),jm(region)))
!    allocate(fz(im(region),jm(region)))
!    allocate(precip(im(region),jm(region),lmax_conv))
!    allocate(dzk(lmax_conv))
!    allocate(kcltop(im(region),jm(region)))
!    allocate(prf(im(region),jm(region),lmax_conv))
!    !
!    ! calculate precipitation formation rate  prf
!    !
!    call calfk
!    !
!    ! initialize cloud top
!    !
!    kcltop(:,:)=lmax_conv
!    !
!    ! calculate and rescale precip
!    !
!    sfz=0.
!    nfz=0     
!    precip(:,:,:)=0.
!    !
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          !
!          ! Calculate precipitation intensity at the bottom of each layer
!          !
!          do k=1,lmax_conv-1
!             ! thickness of layer, only correct in troposphere
!             dzk(k)=scale_heigth*t(i,j,k)*alog(phlb(i,j,k)/(1.+phlb(i,j,k+1)))
!          end do
!          !
!          do k=lmax_conv-1,1,-1
!             precip(i,j,k)=precip(i,j,k+1)+prf(i,j,k)*dzk(k)    !precip: kg/m2/s
!          end do
!          !
!          ! Rescale prf and precip such that these are consistent with ground lsp 
!          !
!          no_fz = 0   ! for statistics   !CMK was not initialised!
!          sfz_no = 0.0
!
!          fz(i,j)=0.
!          !cmk if (lsp(i,j) > 1.e-5) then  old data came in mm/day
!          if (lsp(i,j)*1e3*3600.*24. > 1.e-5) then   !new data are in m/s
!             if (precip(i,j,1) > 0.) then
!                fz(i,j)=lsp(i,j)*1e3/precip(i,j,1) ! m/s ->kg/m2/s   !new unit...
!                !cmk fz(i,j)=lsp(i,j)/3600./24./precip(i,j,1) ! mm/day->kg/m2/s
!                nfz=nfz+1
!                ! precipitation statistics 
!                ! (avoid 'strange' statistics by only few values)
!                sfz=sfz+min(10.,fz(i,j))
!             else
!                ! no precipitation calculated, but ECMWF fields contain rain value
!                no_fz=no_fz+1
!                sfz_no=sfz_no+lsp(i,j)*1e3*86400.   !  (in mm/day)
!             end if
!          else
!             precip(i,j,:)=0.
!          end if
!          do k=1,lmax_conv
!             precip(i,j,k)=precip(i,j,k)*fz(i,j)
!             prf(i,j,k)=prf(i,j,k)*fz(i,j)
!          end do !k
!       end do ! i
!    end do ! j
!    !
!#ifdef MPI  
!    if(myid == root) then
!#endif
!       write(6,*) 'calculate_lsp_scav: region',region
!       write(6,*) '   rainout: average scale factor for precipitation   = ',sfz/nfz
!       write(6,*) '   rainout: percentage of columns with precipitation = ', &
!            100.*nfz/real(im(region)*jm(region)),' %'
!       if ( no_fz > 0 ) write(6,*) 'rainout: lsp and prf not consistent ', &
!            no_fz,'average rainfall [mm/day]',sfz_no/no_fz
!#ifdef MPI   
!    end if
!#endif
!    !
!    ! initialise
!    !
!    oldcov=0.
!    ! evap=0.
!    rloss1(region)%d3=0.
!    rloss2(region)%d3=0.
!    rloss3(region)%d3=0.
!    !
!    do k=lmax_conv-1,1,-1  ! top-down
!       do j=jsr(region),jer(region)
!          do i=isr(region),ier(region)
!             !
!             ! pressure correction for diffusion coefficient
!             !
!             press   = (phlb(i,j,k)+phlb(i,j,k+1))/2.
!             dgairx  = dgair*1e5/press ! dgair at 1 atmosphere
!             dghno3x = dghno3*1e5/press
!             beta1=0.
!             beta2=0.
!             beta3=0.
!             !
!             ! total cloudwater  (kg H2O / kg air)
!             !
!             totw=lwc(i,j,k)+iwc(i,j,k)
!             !
!             ! no influx set kcltop to low number  
!             !
!             if (precip(i,j,k+1)<=1.e-15) kcltop(i,j)=0 
!             !
!             !==========================================
!             !  below-cloud scavenging (beta2 and beta3)
!             !==========================================
!             !
!             ! with evaporation do:
!             !    precip(i,j,k+1)=precip(i,j,k+1)-evap(i,j,k))>1E-15
!             !
!             if( (precip(i,j,k+1)>1e-15) .and. (k<kcltop(i,j)) &
!                  .and. (oldcov(i,j)>0.) )then
!                !
!                ! rflx  in [mm/hr]   !MK?   thought precip was in kg/m2/s ??
!                ! rlwc  in [vol mixing ratio]
!                ! rdrad in [cm]
!                ! ru    in [cm/s] (terminal velocity)
!                ! znre = Reynolds number
!                ! znsc = (Sherwood number)**(1./3.)
!                ! zkg in [cm/s] = dghno3[cm^2/s]/rdrad[cm]
!                !
!                rflx  = precip(i,j,k+1)/oldcov(i,j)*3600.
!                rlwc  = 72.*rflx**0.88*1.e-9
!                rdrad = 0.1*0.3659*rflx**0.21
!                ru    = 9.58*(1.-exp(-(rdrad*10./0.885)**1.147)) 
!                znre  = 20.*rdrad*ru/dgair
!                znsc  = (dgair/dghno3)**(1./3.)
!                znsh  = 1.+0.3*(znre**0.5)*znsc
!                zkg   = dghno3/rdrad*znsh
!                beta2 = 3.*zkg*rlwc/rdrad
!                !
!                ! beta3 for below cloud scavenging of accumulation range aerosol 
!                ! (Dana and Hales, Atmos. Env. 1976, pp. 45-50
!                ! assuming a Whitby aerosol distrbution r=0.034 sigma=2; 
!                ! mass-mean-diameter r=0.14 microm; 
!                ! figure 2 gives a wash-out coefficient of 0.05 mm^-1 (raindepth)
!                ! for other sigma and r look-up tables can be generated
!                !
!                ! mm-1*[mm hr-1] * [hr/s] => [s-1] 
!                !
!                beta3= 0.05*rflx/3600.                          
!                !
!             end if
!             !
!             !  revaporation ( not implemented yet!, evap put to 0.)
!             !
!             !  IF ((1.-cc(i,j,k))>1E-20.AND.precip(i,j,k+1)>1E-20) THEN
!             !    rev1=max(0.,EVAP(i,j,k)/precip(i,j,k+1))
!             !      ! evaporation fraction    
!             !    rev(i,j,k)=MIN(1.,rev1)     
!             !  END IF
!             !
!             !==============================
!             !  in cloud scavenging (beta1)
!             !==============================
!             !
!             if (totw>1.e-9.and.prf(i,j,k)>0.and.cc(i,j,k)>0.05) then     
!                !
!                if(k>kcltop(i,j)) kcltop(i,j)=k     !set new cloud top
!                !
!                ! rlwc: convert mass mixing ratios to volume mixing ratios in cloud
!                ! aird: air density in kg air / m^3
!                ! max_lwc: in kg H2O / m^3
!                ! prf: in kg H2O / m3 / s
!                ! beta: in [s-1] = [vol mixing ratio]*[cm2/s]*1e-4/[m2]
!                ! fac, beta1: in [s-1]
!                !
!                aird=press/t(i,j,k)/rgas*xmair*1.e-3
!                rlwc=(lwc(i,j,k)/rol+iwc(i,j,k)/roi)*aird/cc(i,j,k)
!                rlwc=min(max_lwc/rol,rlwc)
!                !
!                !bas henzing: bug repair: beta=3.*rlwc*(dghno3*1e-4)/rdrad2
!                !bas henzing: should be:  
!                !             beta=(3.*rlwc*(dghno3*1e-4)/(2.*rdrad2))*znsh
!                !bas henzing: reference:  (Roelofs and Lelieveld, 1995)
!                !fd               beta=(3.*rlwc*(dghno3*1e-4)/(2.*rdrad2))*znsh
!                ! fd 13/08/2002 use old defenition again (pers. comm Henzing)
!                beta=3.*rlwc*(dghno3*1e-4)/rdrad2
!                !
!                inc_rdf=(iwc(i,j,k)/totw)*ice_eff+lwc(i,j,k)/totw
!                fac=prf(i,j,k)*inc_rdf/(totw*aird)
!                !
!                !resistance analogy: the slowest determines the overall resistance
!                !
!                beta1=1./(1./beta+1./fac)    
!                !
!                !if no precipitation formation oldcov remains old value
!                !
!                oldcov(i,j)=max(oldcov(i,j),cc(i,j,k))           
!                !
!             end if ! in cloud scavenging
!             !
!             rloss1(region)%d3(i,j,k)=beta1 ! in cloud completely soluble
!             rloss2(region)%d3(i,j,k)=beta2 ! below cloud gas
!             rloss3(region)%d3(i,j,k)=beta3 ! below cloud aerosol
!          end do !i
!       end do !j
!    end do !k 
!#ifdef MPI  
!    if(myid == root) then
!#endif
!       write(*,*) 'calculate_lsp_scav: average rain-out loss rate 1 region', &
!            region,sum(rloss1(region)%d3)/im(region)/jm(region)/lmax_conv
!       write(*,*) 'calculate_lsp_scav: average rain-out loss rate 2 region', &
!            region,sum(rloss2(region)%d3)/im(region)/jm(region)/lmax_conv
!       write(*,*) 'calculate_lsp_scav: average rain-out loss rate 3 region', &
!            region,sum(rloss3(region)%d3)/im(region)/jm(region)/lmax_conv
!#ifdef MPI  
!    end if
!#endif
!
!    deallocate(oldcov)
!    deallocate(fz)
!    deallocate(precip)
!    deallocate(prf)
!    deallocate(dzk)
!    deallocate(kcltop)
!
!    nullify(lsp)
!    nullify(t)
!    nullify(lwc)
!    nullify(iwc)
!    nullify(cc)
!    nullify(phlb)
!
!  contains
!
!    subroutine calfk
!      !--------------------------------------------------------------
!      ! calculate vertical distribution of large scale precipitation
!      !               
!      ! OUTPUT:  prf = precipitation formation rate (kg m-3 s-1)
!      !
!      ! Written by Ad Jeuken, KNMI, May 1998
!      ! Adapted for TM5, MK, march 2003
!      !--------------------------------------------------------------
!      use toolbox, only: lvlpress
!      !
!      ! microphysical constants
!      real,parameter :: cr1=1.e-4    ! s-1
!      real,parameter :: cr2=1.0      ! m2 kg-1      
!      real,parameter :: qcrs=3.e-4   ! kg/kg
!      !bas henzing: replace alfa  real,parameter :: alfa=1.77, beta=0.16
!      !bas henzing: new value alfa=3.29 (Heymsfield and Donner, 1990)
!      real,parameter :: alfa=3.29, beta=0.16
!      !bas henzing: replace b1   real,parameter :: b1=100., b2=0.5, Tberg=268.
!      !bas henzing: new value b1=300. (Sunquist et al., 1989)
!      real,parameter :: b1=300., b2=0.5, Tberg=268.
!      !
!      real                 :: plocal,f1,f2,delta_prec,ahelp,rain_frac,c0
!      real                 :: pr0,qcr,qcl,qci,dzk,aird,press
!      real                 :: qup,qdo,rup,rdo,vtup,vtdo,zcc,ccp,ccm
!      integer              :: iclb,iclt,icldtop,i,j,k,l,l1,l2
!      real,dimension(:),allocatable   :: zrho,pcl,pci
!      !
!      allocate(zrho(lmax_conv))
!      allocate(pci(lmax_conv))
!      allocate(pcl(lmax_conv))
!
!      prf=0.
!      do j=jsr(region),jer(region)
!         do i=isr(region),ier(region)
!            !
!            ! calculate airdensity "zrho" in kg/m3
!            !
!            do l=1,lmax_conv
!               press=(phlb(i,j,l)+phlb(i,j,l+1))/2. 
!               zrho(l)=press/t(i,j,l)/rgas*xmair*1.e-3
!            end do
!            !        
!            iclb=0
!            iclt=0
!            !
!            ! do not consider columns with lsp less than 1.e-5 mm/day
!            !
!            ! if (lsp(i,j)>1.e-5) then 
!            if ( lsp(i,j)*1e3*3600.*24. > 1.e-5 ) then   !new data are in m/s
!               k=1
!               ! determine cloud base
!               do while ((cc(i,j,k)<0.01).and.(k /= lmax_conv)) 
!                  k=k+1
!               end do
!               iclb=k
!               k=lmax_conv
!               ! determine cloud top
!               do while ((cc(i,j,k)<0.01).and.(k /= 1))
!                  k=k-1
!               end do
!               iclt=k
!               ! check for inconsistency in cloud cover fields
!               if ( iclb >= iclt ) iclb=iclt 
!               if ( iclb < 1 ) iclb=1
!               !mvw: replace fixed iclb/t=14 by 120 hPa level (about 15km)
!               !           if (iclb>14) iclb=14  
!               !           if (iclt>14) iclt=14
!               icldtop=lvlpress(region,12000., 98400.)
!               if ( iclb > icldtop ) iclb=icldtop
!               if ( iclt > icldtop ) iclt=icldtop
!               !
!               pr0=0.
!               pcl=0.
!               pci=0.
!               rain_frac=0.00001
!               !
!               ! start at cloudtop
!               do k=iclt,iclb,-1
!                  zcc=max(0.01,cc(i,j,k))
!                  !
!                  ! precipitation formation from ice clouds
!                  !
!                  ! pci in kg_H2O / (kg_air sec)
!                  !
!                  if ( ( t(i,j,k) < 258.0 ) .and. ( k > 1 ) ) then
!                     ccp=max(0.01,cc(i,j,k+1))
!                     ccm=max(0.01,cc(i,j,k-1))
!                     qup=(iwc(i,j,k)/zcc+iwc(i,j,k+1)/ccp)*0.5
!                     qdo=(iwc(i,j,k)/zcc+iwc(i,j,k-1)/ccm)*0.5
!                     rup=(zrho(k)+zrho(k+1))*0.5
!                     rdo=(zrho(k)+zrho(k-1))*0.5
!                     vtup=alfa*(rup*qup)**beta
!                     vtdo=alfa*(rdo*qdo)**beta
!                     pci(k)=grav*(vtup*qup*rup-vtdo*qdo*rdo)/ &
!                          (phlb(i,j,k+1)-phlb(i,j,k))
!                     pci(k)=max(pci(k),0.)
!                  end if
!                  !
!                  ! precipitation formation from liquid clouds
!                  ! formulation ECMWF
!                  !
!                  qcl=lwc(i,j,k)/zcc
!                  qcl=min(max_lwc/zrho(k),qcl)
!                  qcl=max(0.001e-3/zrho(k),qcl)
!                  !
!                  ! pcl in kg_H2O / (kg_air sec)
!                  !
!                  plocal=pr0/rain_frac
!                  f1=1.+b1*sqrt(plocal)
!                  f2=1.+b2*sqrt(max(0.,Tberg-t(i,j,k)))
!                  c0=cr1*f1*f2
!                  qcr=qcrs/(f1*f2)
!                  pcl(k)=zcc*c0*qcl*(1.-exp(-(qcl/qcr)**2))
!                  !
!                  ! prec at top next layer in kg m-2 s-1
!                  !
!                  dzk=29.3*t(i,j,k)*alog(phlb(i,j,k)/(1.+phlb(i,j,k+1)))
!                  delta_prec=(pcl(k)+pci(k))*zrho(k)*dzk
!                  ahelp=(zcc*delta_prec+rain_frac*pr0)/(delta_prec+pr0) 
!                  rain_frac=max(rain_frac,ahelp)
!                  pr0=pr0+(pcl(k)+pci(k))*zrho(k)*dzk
!                  !
!                  ! liquid+ice precipitation formation rates in  kg m-3 s-1
!                  !
!                  prf(i,j,k)= (pcl(k)+pci(k))*zrho(k)
!                  !
!               end do ! k
!            end if ! lsp gt 1.e-5
!         end do  ! i
!      end do ! j
!
!      deallocate(zrho)
!      deallocate(pci)
!      deallocate(pcl)
!      !
!    end subroutine calfk
!
!  end subroutine calculate_lsp_scav
!
!
end module wet_deposition
