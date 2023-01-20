!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module chem_rates
!
!  private
!
!  public :: rates, calrates, calchetnew1, calchetnew2
!
!
!contains
!
!
!  subroutine rates
!    !----------------------------------------------------------------------
!    !     
!    !****  rates calculation of look up tables for rate constants and
!    !      Henry's law constants
!    !
!    !      purpose
!    !      -------
!    !      calculation of look up tables rate constants/henry coefficients
!    !
!    !      interface
!    !      ---------
!    !      call rates
!    !
!    !      method
!    !      ------
!    !      use known array of temperatures (integers) to calculate rate constants
!    !      3 body reactions are explicitly calculated in chemistry
!    !      
!    !      external
!    !      --------
!    !      none
!    !      
!    !      reference
!    !      ---------
!    !      none
!    !
!    !------------------------------------------------------------------
!    use toolbox, only: zfarr
!    use chem_param
!
!    implicit none
!
!    ! local
!    integer ::  itemp,k,i
!    real    :: ztrec,zt3rec,temp
!
!    ! start
!    do k=1,ntemp
!       itemp=k+ntlow
!       temp=float(itemp)
!       ztrec=1./float(itemp)
!       zt3rec=300./float(itemp)     
!       rates_lut(knoo3,k)=zfarr(2.e-12,-1400.,ztrec)
!       rates_lut(kho2no,k)=zfarr(3.7e-12,250.,ztrec)
!       rates_lut(kno2oha,k)=2.47e-30*zt3rec**2.97          !!!wp!!! new ravi
!       rates_lut(kno2ohb,k)=1.45e-11*zt3rec**2.77          !!!wp!!! new ravi
!       rates_lut(kohhno3a,k)=zfarr(2.41e-14,460.,ztrec)    !!!wp!!! new ravi
!       rates_lut(kohhno3b,k)=zfarr(6.51e-34,1335.,ztrec)   !!!wp!!! new ravi
!       rates_lut(kohhno3c,k)=zfarr(2.69e-17,2199.,ztrec)   !!!wp!!! new ravi
!       rates_lut(kno2o3,k)=zfarr(1.2e-13,-2450.,ztrec)
!       rates_lut(knono3,k)=zfarr(1.5e-11,170.,ztrec)
!       rates_lut(kno2no3a,k)=2.2e-30*zt3rec**3.9
!       rates_lut(kno2no3b,k)=1.5e-12*zt3rec**0.7
!       rates_lut(kn2o5,k)=zfarr(2.7e-27,11000.,ztrec)
!       rates_lut(khno4oh,k)=zfarr(1.3e-12,380.,ztrec)
!       rates_lut(kno2ho2a,k)=1.8e-31*zt3rec**3.2
!       rates_lut(kno2ho2b,k)=4.7e-12*zt3rec**1.4
!       rates_lut(khno4m,k)=zfarr(2.1e-27,10900.,ztrec)
!       rates_lut(kodm,k)=.2095*zfarr(3.2e-11,70.,ztrec)+ &
!            .7808*zfarr(1.8e-11,110.,ztrec)
!       rates_lut(kh2ood,k)=2.2e-10
!       rates_lut(ko3ho2,k)=zfarr(1.1e-14,-500.,ztrec)
!       rates_lut(ko3oh,k)=zfarr(1.6e-12,-940.,ztrec)
!       rates_lut(khpoh,k)=zfarr(2.9e-12,-160.,ztrec)
!       rates_lut(kfrmoh,k)=1.e-11
!       rates_lut(kch4oh,k)=zfarr(2.65e-12,-1800.,ztrec)
!       rates_lut(kohmper,k)=zfarr(3.8e-12,200.,ztrec)
!       rates_lut(kohrooh,k)=3.e-12
!       rates_lut(kho2oh,k)=zfarr(4.8e-11,250.,ztrec)
!       rates_lut(kmo2ho2,k)=zfarr(3.8e-13,800.,ztrec)
!       rates_lut(kmo2no,k)=zfarr(4.2e-12,180.,ztrec)
!       rates_lut(kmo2mo2,k)=zfarr(2.5e-13,190.,ztrec)
!       rates_lut(kho2ho2a,k)=zfarr(2.3e-13,600.,ztrec)
!       rates_lut(kho2ho2b,k)=zfarr(1.7e-33,1000.,ztrec)
!       rates_lut(kho2ho2c,k)=zfarr(1.4e-21,2200.,ztrec)
!       rates_lut(kc41,k)=5.8e-16
!       rates_lut(kc43,k)=zfarr(7.0e-12,250.,ztrec)
!       rates_lut(kc44,k)=2.5e-15
!       rates_lut(kc46,k)=zfarr(3.5e-11,-180.,ztrec)
!       rates_lut(kc47a,k)=2.7e-28*zt3rec**(-7.1)
!       rates_lut(kc47b,k)=1.2e-11*zt3rec**0.9
!       rates_lut(kc48,k)=zfarr(2.0e16,-13500.,ztrec)
!       rates_lut(kc49,k)=2.e-12
!       rates_lut(kc50,k)=6.5e-12
!       rates_lut(kc52,k)=8.1e-13
!       rates_lut(kc53,k)=zfarr(1.e15,-8000.,ztrec)
!       rates_lut(kc54,k)=1.6e3
!       rates_lut(kc57,k)=zfarr(5.2e-12,504.,ztrec)
!       rates_lut(kc58,k)=zfarr(4.33e-15,-1800.,ztrec)
!       rates_lut(kc59,k)=7.7e-15
!       rates_lut(kc61a,k)=1.e-28*zt3rec**0.8
!       rates_lut(kc61b,k)=8.8e-12
!       rates_lut(kc62,k)=zfarr(9.14e-15,-2580.,ztrec)
!       rates_lut(kc73,k)=1.7e-11
!       rates_lut(kc76,k)=zfarr(2.54e-11,410.,ztrec)
!       rates_lut(kc77,k)=zfarr(12.3e-15,-2013.,ztrec)
!       rates_lut(kc78,k)=7.8e-13
!       rates_lut(kc79,k)=zfarr(4.2e-12,180.,ztrec)
!       rates_lut(kc80,k)=zfarr(1.7e-14,1300.,ztrec)
!       rates_lut(kc81,k)=6.8e-13
!       rates_lut(kc82,k)=zfarr(3.5e-13,1000.,ztrec)
!       rates_lut(kc83,k)=8.e-11
!       rates_lut(kc84,k)=1.78e-12
!
!       ! sulfur and ammonia gas phase reactions
!
!       ! the hynes et al. parameterisation is replaced by chin et al. 1996
!
!       rates_lut(kdmsoha,k)= 9.6e-12*exp(-234./temp)
!       rates_lut(kdmsohb,k)=1.7e-42*exp(7810./temp)
!       rates_lut(kdmsohc,k)=5.5e-31*exp(7460./temp)
!       rates_lut(kdmsno3,k)=zfarr(1.9e-13,500.,ztrec)!at 298 1.01e-12
!       rates_lut(kso2oha,k)=3.0e-31*(temp/300.)**(-3.3)
!       rates_lut(kso2ohb,k)= 1.5e-12*(temp/300.)
!
!       !
!       ! fate of ammonia
!       !
!       rates_lut(knh3oh,k)=  zfarr(1.7e-12,-710.,ztrec)!1.56e-13 at 298k
!       rates_lut(knh2no,k)=  zfarr(3.8e-12,+450.,ztrec)!1.72e-11
!       rates_lut(knh2no2,k)= zfarr(2.1e-12,650.,ztrec)!1.86e-11
!       rates_lut(knh2ho2,k)= 3.4e-11
!       rates_lut(knh2o2,k)= 6.0e-21
!       rates_lut(knh2o3,k)= zfarr(4.3e-12,-930.,ztrec)!1.89e-13 at 298k
!
!       !
!       ! **** solubility Henry equilibrium
!       !    HNO3/so4/nh4 just a very high number to take H and 
!       !      dissociation into account
!       !
!       henry(:,k)=0.
!       henry(ih2o2,k)=zfarr(1.67e-5,6621.,ZTREC)  
!       henry(ihno3,k)=1e7 
!       henry(ich3o2h,k)=zfarr(1.5e-6,5607.,ZTREC) 
!       henry(ich2o,k)=zfarr(2.7e-6,6425.,ZTREC)
!       henry(irooh,k)=zfarr(1.5e-6,5607.,ZTREC)!(as CH3O2H)
!       henry(iorgntr,k)=zfarr(1.8e-6,6000.,ZTREC) 
!       henry(iso4,k)=1.e7
!       henry(inh4,k)=1.e7
!       henry(imsa,k)=1.e7 
!       henry(iso2,k) =1.2*exp(3120.*ZTREC)*3.41e-5 !correction for the 1/298. part
!       henry(inh3,k) =75.0*exp(3400.*ZTREC)*1.10e-5
!       henry(io3,k)=1.1e-2*exp(2300.*ZTREC)*4.45e-4
!    end do !k temperature loop
!    !
!    ! marked tracers:
!    !
!    henry(io3s,:)      = henry(io3,:)
!
!
!  end subroutine rates
!
!
!
!  subroutine calrates(region,rjx,rr,ye)
!    !**************************************************************
!    !
!    ! CALRATES 	calculate rate constants using lookup table rates_lut
!    !		calculate third bodies 
!    !		calculate heterogeneous removal on aerosols
!    !
!    ! External: CALCHET
!    !
!    !************************************************************
!    !debug use hdf
!
!    use global_data, only: region_dat
!    use chem_param
!    use dims, only: isr,ier, jsr,jer, im, jm
!
!    implicit none
!
!    ! input/output
!    integer, intent(in)                                            :: region
!    real,dimension(im(region),jm(region),njnum)                    :: rjx
!    ! output 
!    ! rr: reaction rates...
!    real,dimension(im(region),jm(region),nreac),intent(out)        :: rr
!    ! ye: extra 2D fields
!    real,dimension(im(region),jm(region),n_extra),intent(inout)    :: ye
!
!    ! local:
!    ! heterogeneous removal rates
!    real,dimension(:,:),allocatable  :: het_nh3, het_n2o5
!
!    ! help variables
!    integer, dimension(:,:), pointer :: zoomed
!    integer                          :: itemp,i,j
!    real                             :: tr, temp, wv, airp, rx1, rx2, rx3
!    real                             :: dum, h2ox, aird, o2
!    real                             :: x1, x2, xice, xliq
!    !
!    ! cloud chemistry of n2o5
!    real                             :: dg, kt_liq, kt_ice
!    real, parameter                  :: r_liq=1e-3, r_ice=5e-3  ! cm
!    real, parameter                  :: g_n2o5_i=0.01, g_n2o5_l=0.04
!
!    !debug integer             :: io,sfstart,sfend
!
!    ! start
!    zoomed => region_dat(region)%zoomed
!    allocate(het_nh3(im(region),jm(region)))
!    allocate(het_n2o5(im(region),jm(region)))
!
!    rx3=0.6
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          temp=ye(i,j,i_temp)
!          itemp=nint(temp-float(ntlow))
!          itemp=min(max(itemp,1),ntemp) !limit temperatures in loop-up table range
!          airp=ye(i,j,i_pres)
!          !
!          ! Calculation of relative humidity [%]
!          ! Richardson's approximation for water vapor pressure
!          ! should be calculated in subroutine rates
!          !
!          h2ox = ye(i,j,ih2on)
!          aird = ye(i,j,iairn)
!          tr=1.-373.15/temp
!          wv=exp((((-.1299*tr-.6445)*tr-1.976)*tr+13.3185)*tr)
!          ye(i,j,irh)=h2ox*temp/(1013.25*wv*7.24e16)
!          ye(i,j,irh) = max(min(ye(i,j,irh),100.0),0.0)   !limit rh between 0-100%
!
!          o2=0.209476*aird
!          !
!          !**** calculate three body reaction rates
!          !
!          rx1=rates_lut(KNO2OHA,itemp)*aird
!          rx2=rates_lut(KNO2OHB,itemp)
!          rr(i,j,kno2oh)=RX1/(1.+RX1/RX2)*RX3**(1./(1.+LOG10(RX1/RX2)**2))
!          rx1=rates_lut(KOHHNO3C,itemp)
!          rx2=rates_lut(KOHHNO3B,itemp)*aird
!          rr(i,j,kohhno3)=rates_lut(KOHHNO3A,itemp)+rx1*rx2/(rx1+rx2)
!          rx1=rates_lut(KNO2NO3A,itemp)*aird
!          rx2=rates_lut(KNO2NO3B,itemp)
!          rr(i,j,kno2no3)=RX1/(1.+RX1/RX2)*RX3**(1./(1.+LOG10(RX1/RX2)**2))
!          rx1=rates_lut(KNO2HO2A,itemp)*aird
!          rx2=rates_lut(KNO2HO2B,itemp)
!          rr(i,j,kno2ho2)=RX1/(1.+RX1/RX2)*RX3**(1./(1.+LOG10(RX1/RX2)**2))
!          rx1=rates_lut(KC61A,itemp)*aird
!          rx2=rates_lut(KC61B,itemp)
!          rr(i,j,kc61)=RX1/(1.+RX1/RX2)*RX3**(1./(1.+LOG10(RX1/RX2)**2))
!          !
!          ! photolysis rates and 2 body reactions    
!          !
!          rr(i,j,knoo3)=rates_lut(KNOO3,itemp)  
!          rr(i,j,kho2no)=rates_lut(KHO2NO,itemp)
!          rr(i,j,kmo2no)=rates_lut(Kmo2NO,itemp)
!          rr(i,j,kno2o3)=rates_lut(KNO2O3,itemp)
!          rr(i,j,knono3)=rates_lut(KNONO3,itemp)
!          rr(i,j,kn2o5)=rr(i,j,kno2no3)/rates_lut(KN2O5,itemp)
!          rr(i,j,khno4oh)=rates_lut(KHNO4OH,itemp) 
!          rr(i,j,khno4m)=rr(i,j,kno2ho2)/rates_lut(KHNO4M,itemp)
!          rr(i,j,kodm)=rates_lut(KODM,itemp)
!          rr(i,j,kh2ood)=rates_lut(KH2OOD,itemp)
!          rr(i,j,ko3ho2)=rates_lut(KO3HO2,itemp)
!          rr(i,j,kcooh)=1.5E-13+9E-14*airp/101325. 
!          rr(i,j,ko3oh)=rates_lut(KO3OH,itemp)
!          rr(i,j,khpoh)=rates_lut(KHPOH,itemp)
!          rr(i,j,kfrmoh)=rates_lut(KFRMOH,itemp)
!          rr(i,j,kch4oh)=rates_lut(KCH4OH,itemp)
!          rr(i,j,kh2oh)=1.06*rates_lut(KCH4OH,itemp)*550.e-9*aird !H2=550ppbv
!          rr(i,j,kohmper)=rates_lut(KOHMPER,itemp)
!          rr(i,j,kohrooh)=rates_lut(KOHROOH,itemp)
!          rr(i,j,kmo2ho2)=rates_lut(KMO2HO2,itemp)
!          rr(i,j,kmo2mo2)=rates_lut(KMO2MO2,itemp)
!          rr(i,j,kho2oh)=rates_lut(KHO2OH,itemp)
!          rr(i,j,kho2ho2)=(rates_lut(KHO2HO2A,itemp) + &
!               rates_lut(KHO2HO2B,itemp)*aird)*  &
!               (1.+rates_lut(KHO2HO2C,itemp)*h2ox)
!          rr(i,j,kc41)=rates_lut(KC41,itemp)
!          rr(i,j,kc43)=rates_lut(KC43,itemp)
!          rr(i,j,kc44)=rates_lut(KC44,itemp)
!          rr(i,j,kc46)=rates_lut(KC46,itemp)
!          rx1=rates_lut(KC47A,itemp)*aird*0.7808
!          rx2=rates_lut(KC47B,itemp)
!          rr(i,j,kc47)=0.96*RX1/(1.+RX1/RX2)*0.3**(1./(1.+LOG10(RX1/RX2)**2))
!          rr(i,j,kc48)=rates_lut(KC48,itemp)
!          rr(i,j,kc49)=rates_lut(KC49,itemp)
!          rr(i,j,kc50)=rates_lut(KC50,itemp)
!          rr(i,j,kc52)=rates_lut(KC52,itemp)    
!          rr(i,j,kc53)=rates_lut(KC53,itemp)
!          rr(i,j,kc54)=rates_lut(KC54,itemp)
!          rr(i,j,kc57)=rates_lut(KC57,itemp)
!          rr(i,j,kc58)=rates_lut(KC58,itemp)
!          rr(i,j,kc59)=rates_lut(KC59,itemp)
!          rr(i,j,kc62)=rates_lut(KC62,itemp)
!          rr(i,j,kc73)=rates_lut(KC73,itemp)
!          rr(i,j,kc76)=rates_lut(KC76,itemp)
!          rr(i,j,kc77)=rates_lut(KC77,itemp)
!          rr(i,j,kc78)=rates_lut(KC78,itemp)
!          rr(i,j,kc79)=rates_lut(KC79,itemp)
!          rr(i,j,kc80)=rates_lut(KC80,itemp)
!          rr(i,j,kc81)=rates_lut(KC81,itemp)
!          rr(i,j,kc82)=rates_lut(KC82,itemp)
!          rr(i,j,kc83)=rates_lut(KC83,itemp)
!          rr(i,j,kc84)=rates_lut(KC84,itemp)
!          rr(i,j,kc85)=rr(i,j,kc81)*rr(i,j,kc82)/rr(i,j,kc79)
!          rjx(i,j,jo3d)=rjx(i,j,jo3d)*rr(i,j,kh2ood)*h2ox / &
!               ( rr(i,j,kodm)*aird + rr(i,j,kh2ood)*h2ox ) 
!          ! this is fraction rjo3d=>OH 
!
!          RX1=rates_lut(kso2oha,itemp)*aird
!          RX2=rates_lut(kso2ohb,itemp)
!          rr(i,j,kso2oh)=RX1/(1.+RX1/RX2)*0.6**(1./(1.+LOG10(RX1/RX2)**2))
!          !
!          ! dmsoha => so2
!          ! dmsohb => 0.75 SO2 + 0.25 MSA
!          !
!          rr(i,j,kdmsoha)=rates_lut(kdmsoha,itemp)
!          rr(i,j,kdmsohb)=rates_lut(kdmsohb,itemp)*o2/ &
!               (1.+rates_lut(kdmsohc,itemp)*o2)
!          rr(i,j,kdmsno3)=rates_lut(kdmsno3,itemp)
!
!          ! ammonia chemistry
!          rr(i,j,knh3oh)=rates_lut(knh3oh,itemp)
!          rr(i,j,knh2no)=rates_lut(knh2no,itemp)
!          rr(i,j,knh2no2)=rates_lut(knh2no2,itemp)
!          rr(i,j,knh2ho2)=rates_lut(knh2ho2,itemp)
!          rr(i,j,knh2o2)=rates_lut(knh2o2,itemp)*o2
!          rr(i,j,knh2o3)=rates_lut(knh2o3,itemp)  
!          rr(i,j,krn222)=2.11e-6    ! s-1 decay time Rn222 to Pb210
!       end do
!    end do   !_lat/lon loop
!
!    ! calculate heterogeneous removal constants of n2o5
!
!    call calchetnew2(region,ye,het_n2o5,1) !n2o5
!    call calchetnew2(region,ye,het_nh3,2) !nh3
!
!    !
!    !  heterogeneous reaction of N2O5 and H2O -> 2 HNO3 on cloud and aerosol
!    !  included in gas phase chemistry
!    !
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          !
!          ! kt= (r2/3Dg + 4*r/3vgamma)^-1
!          ! ice   r=50 micrometer gamma = 0.01 
!          ! water r=10 micrometer  gamma = 0.05
!          ! v is 4e5 cm/s and Dg is 0.1 cm2/s at standard press
!          airp=ye(i,j,i_pres)
!          dg=0.1*1e5/airp      !simple approximation for diffusion coeff. [cm2/s]
!          kt_liq=1./(r_liq*r_liq/3./dg+4.*r_liq/3./4e5/g_n2o5_l)
!          kt_ice=1./(r_ice*r_ice/3./dg+4.*r_ice/3./4e5/g_n2o5_i)
!          aird = ye(i,j,iairn)
!          xliq = ye(i,j,ilwc)
!          xice = ye(i,j,iiwc)
!          rr(i,j,kn2o5l)=(kt_liq*xliq+kt_ice*xice)           ! cloud
!          !
!          ! kn2o5aq and nh3so4 can be done implicitly, 
!          ! it has occurred that these rates have
!          ! become negative over antarctica (aug 1993), 
!          ! therefore put minimum value of 0. (AJ jul1999)
!          !
!          !cmk    rr(i,j,kn2o5aq)=max(0.,het_n2o5(i,j)/
!          !       1e-9*(y(jl,iso4)+y(jl,imsa))/aird
!          !cmk multiplication moved to EBI
!          rr(i,j,kn2o5aq)=max(0.,het_n2o5(i,j))/1e-9/aird
!          !
!          ! knh3so4 is uptake coefficient on H2SO4. 
!          ! 1 uptake of NH3 consumes 1 acid molecule.  
!          ! 
!          rr(i,j,knh3so4)=max(0.,het_nh3(i,j))/1e-9/aird
!          !  rr(i,j,knh3so4)=0.0   !CMK sep2003: why was this here?
!          !
!       end do
!    end do
!
!    deallocate(het_n2o5)
!    deallocate(het_nh3)
!    nullify(zoomed)
!    !debug if(level==1) then
!    !debug io = sfstart('rr.hdf',dfacc_create)
!    !debug call io_write(io,im(region),'im',jm(region),'jm',lm(region),'lm',t,'t')
!    !debug call io_write(io,im(region),'im',jm(region),'jm',lm(region),'lm',q,'q')
!    !debug call io_write(io,im(region),'im',jm(region),'jm',lm(region),'lm',clwc,'lwc')
!    !debug call io_write(io,im(region),'im',jm(region),'jm',lm(region),'lm',ciwc,'iwc')
!    !debug call io_write(io,im(region),'im',jm(region),'jm',nrat,'nrat',rr,'rr')
!    !debug call io_write(io,nrat,'nrat',ntemp,'ntemp',rates_lut,'rates_lut')
!    !debug call io_write(io,ntrace,'ntrace',ntemp,'ntemp',henry,'henry')
!    !debug io = sfend(io)
!    !debug endif
!  end subroutine calrates
!
!
!
!  subroutine calchet1(gamma,xmw,icomp)
!    !----------------------------------------------------------------------
!    !     
!    !****  calchet1 - pre-calculate heterogeneous removal rates on sulfate aerosol
!    !
!    !      programmed by frank dentener 01.04.96
!    !      modified by Maarten krol
!    !
!    !      purpose
!    !      -------
!    !      calculate heterogeneous removal constants for specified species
!    !      on sulfate aerosol as function of concentration, 
!    !      relative humidity and pressure
!    !    
!    !      interface
!    !      ---------
!    !      calchet1(gamma,xmw,icomp)
!    !
!    !            gamma dimensionless accomodation coefficient
!    !            xmw molar weight [g/mol]
!    !            icomp (compound number)
!    !
!    !      method
!    !      ------
!    !      use Whitby sulfate distribution, and Fuchs' rate expression
!    !      to integrate rate coefficient on aerosol distribution
!    !      
!    !      external
!    !      --------
!    !      none
!    !      
!    !      reference
!    !      ---------
!    !      Dentener (1993) Ph.D. thesis
!    !
!    !------------------------------------------------------------------
!    use binas, xgamma=> gamma
!    use chem_param
!
!    implicit none
!
!    ! input
!    real,intent(in)   :: xmw,gamma
!    integer,intent(in):: icomp
!
!    ! local
!    integer           :: ip,isat,i
!    real              :: press,temp,dxm,dn2o5,vsp,xl,aird,aervol
!    real              :: hsat1,hsat2,raer,rx1,zlogs,rx2
!    real              :: FN1,FN2,FR1,FR2,FA1,FA2,FV1,FV2,rmean,qi30,xkn,xlab,xfac
!    real,parameter    :: RG=8.314E3,VENT=1.0, FLN10=2.302585,W2PI=2.506638
!    real,parameter    :: xmnso4=xmso4+xmh+xmnh4,p1=1.,t1=288.,g=1.40,conc=1e-9
!    real,dimension(3),parameter :: apar=(/1.,3.4e-8,0.301/)
!    ! quantities of integration (e.g.number surface, volume and rate coefficient  
!    integer,parameter :: nt=4
!    integer,parameter :: nint1=2000    ! number of integration intervals
!    ! rint: integration stepsizes [m],0-1 um,1-100 um
!    real,dimension(nint1) ::  rint = &
!         (/ (.001E-6, i=1,1000), (0.1E-6, i=1,1000) /)
!    real,dimension(nt) ::  qi
!
!    !
!    !     1 particle (unity)/cm3, radius and log(sigma) measurements from Whitby
!    !     the radius is assumed to be 'dry radius
!    !     We take aerosol size from Whitby accumulation mode (1978)
!    !      Numbermean radius: 0.034um, sigma=2, 1 (unity) particles cm-3
!    !     Molecular weight NH4HSO4 111 g/mol
!    !     aerosol density of dry NH4HSO4 1.8 E3 kg/m3= 1.8 gcm-3
!    !     temperature is not a determining factor is implicitly accounted for
!    !     as a function of pressure.
!    !     temperature is assumed to follow an adiabatic lapse rate:
!    !    (T2/T1)=(P2/P1)^{(g-1)/g} function of pressure with g=Cp/Cv=ca. 1.40
!    !
!
!    ! start
!
!    print *,'calchet1: initialize heterogeneous rem. rates'
!
!    ! pressure from 1000 to 0 mbar     
!    do ip=1,11
!       press=max(0.001,1.1-ip*0.1)         !atmosphere (minimum is 1 hPa)
!       temp=max(210.,t1*(press/p1)**((g-1)/g)) 
!       !this estimate of temp is not very accurate
!       DN2O5=4.6e-6*TEMP**1.75/PRESS*1E-4  ! diffusion coefficient for n2o5 [m2/s] 
!       dxm=dn2o5*xmw/xmn2o5             ! diffusion coefficient for other component
!       VSP=SQRT(8.*RG*TEMP/PI/XMW)      ! Molecular speed [m/s]
!       XL=3.*DXM/VSP                    ! free molecular path length [m]
!
!       aird=press/(rg*temp)*1e2         ! (mole/cm3) 
!       aervol=conc*aird*xmnso4/aerdens  ! (mole/cm3) * (g/mole)/(g/cm3)=>[cm3/cm3]
!       ! aervol is the volume of 1 pbbv dry nh4hso4 at temp and press
!       do isat=1,11
!          hsat1=-10.+10.*isat
!          HSAT1=AMIN1(HSAT1,95.)            !max RH of 95 %
!          HSAT2=HSAT1*HSAT1
!          RAER=AMAX1(1.,1.0129231-0.0041328044*hsat1+0.00070336143*hsat2&
!               -1.4388956e-05*hsat1*hsat2+9.1359802e-08*hsat2*hsat2)
!          ! growth of aerosol radius due to deliquescence for NH4HSO4
!
!          ! actual integration
!          RX1=0.0
!          QI(:)=0. 
!          ZLOGS=APAR(3)
!          rmean=apar(2)*raer
!          do I=1,NINT1-1
!             RX1=RX1+RINT(I)  
!             RX2=RX1+RINT(I+1)                 !RX2-RX1 is integration size interval
!             XKN= 3.*DXM/RX1/VSP               !Knudsen number
!             XLAB=(XKN*4./3.+0.71)/(XKN+1.) 
!             FN1=(LOG10(RX1/rmean))**2  
!             FN2=(LOG10(RX2/rmean))**2
!             FN1= EXP(-FN1/2./ZLOGS/ZLOGS)/RX1 !number integration
!             FN2= EXP(-FN2/2./ZLOGS/ZLOGS)/RX2 !number integration      
!             FR1=1./(1.+XKN*(XLAB+(4.*(1.-GAMMA)/3./GAMMA)))*RX1*FN1 
!             !reactivity integration 
!             FR2=1./(1.+XKN*(XLAB+(4.*(1.-GAMMA)/3./GAMMA)))*RX2*FN2
!             FA1= RX1*RX1*FN1                  !surface integration
!             FA2= RX2*RX2*FN2                  !surface integration
!             FV1= RX1*RX1*RX1*FN1              !volume integration
!             FV2= RX2*RX2*RX2*FN2              !volume integration
!             QI(1)=QI(1)+RINT(I)/2.*(FN1+FN2)  !EULER INTEGRATION
!             QI(2)=QI(2)+RINT(I)/2.*(FA1+FA2)            
!             QI(3)=QI(3)+RINT(I)/2.*(FV1+FV2)       
!             QI(4)=QI(4)+RINT(I)/2.*(FR1+FR2)   
!          end do  !I=1,NINT1-1 
!          xfac=APAR(1)*1.e6/FLN10/W2PI/zlogs   ! constant integration factor     
!          QI(1)=QI(1)*xfac*1.E-6               ! conversion cm3=>m3 number  
!          QI(2)=QI(2)*4.*PI*xfac*1.e-2         ! m=>cm surface
!          QI(3)=QI(3)*4./3.*PI*xfac            !volume
!          QI(4)=QI(4)*4.*PI*DXM*xfac*vent      !reactivity
!          if (isat == 1) qi30=qi(3)            ! dry volume
!          hetrem(isat,ip,icomp)=aervol/qi30*qi(4)! removal coefficient
!       end do !isat
!    end do !ip
!
!  end subroutine calchet1
!
!
!
!  subroutine calchet2(region,ye,het,icomp)
!    !----------------------------------------------------------------------
!    !     
!    !****  calchet2 - calculate heterogeneous removal rates on sulfate aerosol
!    !
!    !      programmed by frank dentener 01.04.96
!    !      modified for TM5 by maarten krol jan 2002
!    !
!    !      purpose
!    !      -------
!    !      calculate heterogeneous removal constants for specified species
!    !      on sulfate aerosol as function of concentration, 
!    !      relative humidity and pressure
!    !    
!    !      interface
!    !      ---------
!    !      calchet2(region,ye,het,icomp)
!    !
!    !            region to indicate for which zoom region
!    !            ye  : extra fields (for rh,pressure).
!    !            het heterogeneous removal constant [s-1/ppbv]
!    !            icomp (1,2)   compound (N2O5, NH3)
!    !
!    !      method
!    !      ------
!    !      use Whitby sulfate distribution, and Fuchs' rate expression
!    !      to integrate rate coefficient on aerosol distribution
!    !      
!    !      external
!    !      --------
!    !      none
!    !      
!    !      reference
!    !      ---------
!    !      Dentener (1993) Ph.D. thesis
!    !
!    use global_data, only: region_dat
!    use binas, xgamma=>gamma
!    use dims, only: isr,ier, jsr,jer, im, jm
!    use chem_param
!
!    implicit none
!
!    ! input
!    integer,intent(in)                            :: region,icomp
!    real,dimension(im(region),jm(region))         :: het      ! result
!    real,dimension(im(region),jm(region),n_extra) :: ye
!    ! ye: extra fields (rh, pressure)
!
!    ! local
!    integer,dimension(:,:),pointer                :: zoomed
!    real    :: pres,px,hx,hp1,hp2
!    integer :: np1,np2,nh1,nh2,i,j
!
!    ! relative humidity should be between 0-100 % and 
!    ! pressure between 105000 and 0 Pa
!    ! actual interpolation....
!    ! hetrem field in in module.....
!
!    zoomed => region_dat(region)%zoomed
!
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          pres = ye(i,j,i_pres)
!          np1=min(11,1+nint(10.-pres/10000.))  !pressure
!          np1=max(1,np1)
!          np2=min(11,np1+1)            
!          nh1=max(1,nint(ye(i,j,irh)/10.+0.5))    !relative humidity
!          nh1=min(11,nh1)    
!          nh2=min(11,nh1+1)
!          px=((11-np1)*10000.-pres)/10000.
!          hx=(ye(i,j,irh)-(nh1-1)*10.)/10.
!          hp1=px*hetrem(nh1,np2,icomp)+(1.-px)*hetrem(nh1,np1,icomp)
!          hp2=px*hetrem(nh2,np2,icomp)+(1.-px)*hetrem(nh2,np1,icomp)
!          het(i,j)=hx*hp2+(1.-hx)*hp1   
!       end do
!    end do
!
!    nullify(zoomed)
!
!  end subroutine calchet2
!
!
!
!  subroutine calchetnew1(gamma,xmw,icomp)
!    !----------------------------------------------------------------------
!    !
!    !****  calchetnew1 - pre- calculate heterogeneous removal rates 
!    !                    on sulfate aerosol
!    !
!    !      programmed by frank dentener 01.04.96
!    !      modified MK oct 2003: splitted in two. 
!    !
!    !      purpose
!    !      -------
!    !      calculate heterogeneous removal constants for specified species
!    !      on sulfate aerosol as function of concentration, 
!    !      relative humidity and pressure
!    !
!    !      interface
!    !      ---------
!    !      calchetnew1(gamma,xmw,icomp)
!    !
!    !            gamma dimensionless accomodation coefficient
!    !            xmw molar weight [g/mol]
!    !            icomp component number: 1: n2o5 2:nh3
!    !
!    !      method
!    !      ------
!    !      use Whitby sulfate distribution, and Fuchs' rate expression
!    !      to integrate rate coefficient on aerosol distribution
!    !
!    !      external
!    !      --------
!    !      none
!    !
!    !      reference
!    !      ---------
!    !      Dentener (1993) Ph.D. thesis
!    !
!    !------------------------------------------------------------------
!    use toolbox,       only : escape_tm 
!    use reaction_data, only : ncomponent, hetrem, aerdens, rra
!    use reaction_data, only : nr_interval, np_interval
!
!    implicit none
!
!    ! input/output
!    real, intent(in)    :: gamma
!    real, intent(in)    :: xmw
!    integer, intent(in) :: icomp
!
!    !local
!    integer                     :: ip,i,iaero
!    ! quantities of integration e.g. number surface, volume and rate coefficient
!    integer,parameter           :: nt=4
!    integer,parameter           :: nint1=2000    !number of integration intervals
!    real :: rx1,rx2,fn1,fn2,fr1,fr2,fa1,fa2,fv1,fv2
!    real :: zlogs,rmean,qi(nt),qi30
!    real,dimension(nint1) ::  rint = &
!         (/ (.001E-6, i=1,1000), (0.1E-6, i=1,1000) /)
!    !
!    ! needed for calculation of reaction constants
!    real,parameter :: xmn2o5 = 108.
!    real,parameter :: rg     = 8.314e3
!    real,parameter :: vent   = 1.0
!    real,parameter :: pi     = 3.14159
!    real,parameter :: fln10  = 2.302585
!    real,parameter :: w2pi   = 2.506638
!    real,parameter :: avo    = 6.0e23
!    real,parameter :: xmnso4 = 111.
!    real,parameter :: p1     = 1.
!    real,parameter :: t1     = 288.
!    real,parameter :: g      = 1.40
!    real,parameter :: conc   = 1e-9
!    !
!    real :: temp, press, dn2o5, vsp, xl, xfac, xkn, xlab, dxm
!    real :: raer,aervol,aird
!    real,dimension(3) :: apar=(/1.,3.4e-8,0.301/)
!
!    !     1 particle (unity)/cm3, radius and log(sigma) measurements from Whitby
!    !     the radius is assumed to be 'dry radius
!    !     We take aerosol size from Whitby accumulation mode (1978)
!    !      Numbermean radius: 0.034um, sigma=2, 1 (unity) particles cm-3
!    !     Molecular weight NH4HSO4 115 g/mol
!    !     aerosol density of dry NH4HSO4 1.8 E3 kg/m3= 1.8 gcm-3
!    !     temperature is not a determining factor is implicitly accounted for
!    !     as a function of pressure.
!    !     temperature is assumed to follow an adiabatic lapse rate:
!    !    (T2/T1)=(P2/P1)^{(g-1)/g} function of pressure with g=Cp/Cv=ca. 1.40
!    !
!
!    if ( icomp > ncomponent ) then
!       call escape_tm('calchetnew1: Check component in calchetnew1')
!    end if
!    print *,'calchetnew1: initialize heterogeneous rem. rates ', icomp
!
!    do ip=1,np_interval !pressure from 1000 to 0 mbar
!
!       press=max(0.001,1.1-ip*0.1)         !atmosphere (minimum is 1 hPa)
!       temp=max(210.,t1*(press/p1)**((g-1)/g))
!       !this estimate of temp is not very accurate
!       DN2O5=4.6e-6*TEMP**1.75/PRESS*1E-4  ! diffusion coefficient for n2o5 [m2/s]
!       dxm=dn2o5*xmw/xmn2o5              ! diffusion coefficient for other component
!       VSP=SQRT(8.*RG*TEMP/PI/XMW)       ! Molecular speed [m/s]
!       XL=3.*DXM/VSP                     ! free molecular path length [m]
!
!       aird=press/(rg*temp)*1e2          ! (mole/cm3)
!       aervol=conc*aird*xmnso4/aerdens   ! (mole/cm3) * (g/mole)/(g/cm3)=>[cm3/cm3]
!       ! aervol is the volume of 1 pbbv dry nh4hso4 at temp and press
!
!       do iaero=1,nr_interval ! aerosol increased radius loop
!
!          ! RAER: growth of aerosol radius due to deliquescence for NH4HSO4
!          RAER=rra(iaero)
!
!          ! actual integration
!          RX1=0.0
!          QI(:)=0.
!          ZLOGS=APAR(3)
!          rmean=apar(2)*raer
!
!          do I=1,NINT1-1
!
!             RX1=RX1+RINT(I)
!             RX2=RX1+RINT(I+1)                 !RX2-RX1 is integration size interval
!             XKN= 3.*DXM/RX1/VSP               !Knudsen number
!             XLAB=(XKN*4./3.+0.71)/(XKN+1.)
!             FN1=(LOG10(RX1/rmean))**2
!             FN2=(LOG10(RX2/rmean))**2
!             FN1= EXP(-FN1/2./ZLOGS/ZLOGS)/RX1 !number integration
!             FN2= EXP(-FN2/2./ZLOGS/ZLOGS)/RX2 !number integration
!             FR1=1./(1.+XKN*(XLAB+(4.*(1.-GAMMA)/3./GAMMA)))*RX1*FN1
!             !reactivity integration
!             FR2=1./(1.+XKN*(XLAB+(4.*(1.-GAMMA)/3./GAMMA)))*RX2*FN2
!             FA1= RX1*RX1*FN1                  !surface integration
!             FA2= RX2*RX2*FN2                  !surface integration
!             FV1= RX1*RX1*RX1*FN1              !volume integration
!             FV2= RX2*RX2*RX2*FN2              !volume integration
!             QI(1)=QI(1)+RINT(I)/2.*(FN1+FN2)  !EULER INTEGRATION
!             QI(2)=QI(2)+RINT(I)/2.*(FA1+FA2)
!             QI(3)=QI(3)+RINT(I)/2.*(FV1+FV2)
!             QI(4)=QI(4)+RINT(I)/2.*(FR1+FR2)
!
!          end do  !I=1,NINT1-1
!
!          xfac=APAR(1)*1.e6/FLN10/W2PI/zlogs   ! constant integration factor
!          QI(1)=QI(1)*xfac*1.E-6               ! conversion cm3=>m3 number
!          QI(2)=QI(2)*4.*PI*xfac*1.e-2         ! m=>cm surface
!          QI(3)=QI(3)*4./3.*PI*xfac            !volume
!          QI(4)=QI(4)*4.*PI*DXM*xfac*vent      !reactivity
!          if (iaero == 1) qi30=qi(3)            ! dry volume
!          hetrem(iaero,ip,icomp)=aervol/qi30*qi(4)! removal coefficient
!
!       end do !iaero
!       write(*,'(a,i3,19(1x,1pe9.1))') ' calchetnew1: ',ip,hetrem(:,ip,icomp)
!
!    end do !ip
!    
!  end subroutine calchetnew1
!
!
!
!  subroutine calchetnew2(region,ye,het,icomp)
!    !----------------------------------------------------------------------
!    !     
!    !****  calchetnew2 - calculate heterogeneous removal rates on sulfate aerosol
!    !
!    !      programmed by frank dentener 01.04.96
!    !      modified for TM5 by maarten krol jan 2002
!    !
!    !      purpose
!    !      -------
!    !      calculate heterogeneous removal constants for specified species
!    !      on sulfate aerosol as function of concentration, 
!    !      relative humidity and pressure
!    !    
!    !      interface
!    !      ---------
!    !      calchetnew2(region,ye,het,icomp)
!    !
!    !            region to indicate for which zoom region
!    !            ye  : extra fields (for rh,pressure).
!    !            het heterogeneous removal constant [s-1/ppbv]
!    !            icomp (1,2)   compound (N2O5, NH3)
!    !            
!    !
!    !      method
!    !      ------
!    !      use Whitby sulfate distribution, and Fuchs' rate expression
!    !      to integrate rate coefficient on aerosol distribution
!    !      1 particle (unity)/cm3, radius and log(sigma) measurements from Whitby
!    !      the radius is assumed to be 'dry radius
!    !      We take aerosol size from Whitby accumulation mode (1978)
!    !        Numbermean radius: 0.034um, sigma=2, 1 (unity) particles cm-3
!    !      Molecular weight NH4HSO4 115 g/mol
!    !      aerosol density of dry NH4HSO4 1.8 E3 kg/m3= 1.8 gcm-3
!    !      temperature is not a determining factor is implicitly accounted for
!    !      as a function of pressure.
!    !      temperature is assumed to follow an adiabatic lapse rate:
!    !      (T2/T1)=(P2/P1)^{(g-1)/g} function of pressure with g=Cp/Cv=ca. 1.40
!    !      
!    !      external
!    !      --------
!    !      none
!    !      
!    !      reference
!    !      ---------
!    !      Dentener (1993) Ph.D. thesis
!    !----------------------------------------------------------------------
!
!    use global_data,   only : region_dat
!    use binas,         only : xgamma=>gamma
!    use dims,          only : isr,ier, jsr,jer, im, jm
!    use reaction_data, only : ncomponent, hetrem, nr_interval
!    use reaction_data, only : np_interval, aerdens, rra
!    use chem_param,    only : n_extra, irinc, i_pres
!
!    implicit none
!
!    ! input
!    integer, intent(in)                           :: region
!    integer, intent(in)                           :: icomp
!    real,dimension(im(region),jm(region))         :: het      !result (time scale)
!    real,dimension(im(region),jm(region),n_extra) :: ye 
!    ! ye: extra fields ( pressure, rinc)
!
!    ! local
!    integer,dimension(:,:),pointer                :: zoomed
!
!    real    :: pres,px,hx,hp1,hp2
!    integer :: np1,np2,nh1,nh2,i,j,jr, nr1, nr2
!
!    ! start
!
!    zoomed => region_dat(region)%zoomed
!
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          pres = ye(i,j,i_pres)
!          np1=min(np_interval,1+nint(10.-pres/10000.))
!          np1=max(1,np1)
!          np2=min(np_interval,np1+1)
!          nr1=1
!          do jr=1,nr_interval
!             if(ye(i,j,irinc).ge.rra(jr)) nr1=jr  ! lower bound of rinc array
!          end do
!          nr2=min(nr1+1,nr_interval)  ! upper bound of rinc
!          px=((np_interval-np1)*10000.-pres)/10000.
!          hx=1.
!          if (nr1.ne.nr2) hx=(ye(i,j,irinc)-rra(nr1))/(rra(nr2)-rra(nr1))
!          hp1=px*hetrem(nr1,np2,icomp)+(1.-px)*hetrem(nr1,np1,icomp)
!          hp2=px*hetrem(nr2,np2,icomp)+(1.-px)*hetrem(nr2,np1,icomp)
!          het(i,j)=hx*hp2+(1.-hx)*hp1
!       end do
!    end do
!    nullify(zoomed)
!
!  end subroutine calchetnew2
!
!
end module chem_rates
