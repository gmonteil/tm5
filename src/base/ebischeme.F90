!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module ebischeme
!
!  !--------------------------------------------------------------------------
!  ! Eulerian backward Iteration
!  ! Chemistry solver for the CBM4 scheme 
!  !--------------------------------------------------------------------------
!
!  implicit none
!
!  private
!
!  public :: ebi
!
!
!contains
!
!
!  subroutine ebi(region,level,rj,rr,y,ye)
!    !
!    ! perform Eulerian backwards chemistry at one model layer level in region 
!    ! input: 
!    ! region : region 
!    ! level:   model layer...
!    ! rj    : array of photolysis rates
!    ! rr    : array of reaction rates
!    ! y     : vector of concentrations (to be returned...)
!    !
!#ifdef MPI
!    use mpi_const,only : lmar,lmloc,myid
!#endif
!    use chem_param     
!    use dims, only: isr, ier, jsr,jer,im ,jm, nchem, tref, okdebug
!    use global_data, only: region_dat
!    use toolbox, only: dumpfield
!    use dry_deposition, only: vd  
!    !type(emis_data),dimension(nregions,ntrace) :: vd
!
!    implicit none
!
!    ! input/output
!    integer, intent(in)                                   :: region,level
!    real,dimension(im(region),jm(region),nreac),intent(in):: rr
!    real,dimension(im(region),jm(region),njnum),intent(in):: rj
!    real,dimension(im(region),jm(region),maxtrace)        :: y
!    real,dimension(im(region),jm(region),n_extra)         :: ye
!
!    ! local
!    integer,dimension(:,:),pointer        :: zoomed
!    real,dimension(:,:,:),allocatable     :: cr2,cr3,cr4   !reaction budgets...
!    real,dimension(:,:,:),allocatable     :: y0
!
!    real    :: dtime
!    integer :: iterebi,i,j,ib,maxit,iter
!    integer :: offsetl
!
!    real :: R57,R89,P1,R12,R21,XL1,P2,XL2,P3,XL3,X1,X2,X3
!    real :: C1,C2,C3,Y2,XJT,R21T,R12T,R12TC,R21TC,XJTC,ACUB,BCUB,CCUB
!    real :: CUBDET,DNO2,R56,R65,R75,P5,XL5,R66,X5,P6,XL6,X6,C6,XL7
!    real :: Y1,C7,R98,P8,XL8,X4,C5,XL9,R1920,R1919,P19,XL19,R2019,XL20
!    real :: XLHNO3,PH2O2,XLH2O2,PCH2O,PCO,PHNO3,XLCH2O,PCH3O2,XLCH3O2
!    real :: PCH3O2H,XLCH3O2H,PALD2,XLALD2,PMGLY,XLMGLY,POLE,XLETH
!    real :: XLOLE,XLISOP,PRXPAR,XLRXPAR,PPAR,XLPAR,PROR,XLROR,PXO2
!    real :: XLXO2,PXO2N,XLXO2N,PROOH,XLROOH,PORGNTR,XLORGNTR,XLCO
!
!    real :: dt,dt2
!    real :: qdms,pso2,qso2,qso2d,qnh3,pnh2,qnh2,qdms1,qdms2,pmsa,fnoy
!
!    integer :: io,sfstart,sfend
!
!    ! for vectorization/blocking ....
!    ! npts can be varied to optimize cache memory management.
!    integer,parameter             :: npts=15 
!    integer,dimension(npts)       :: ipts,jpts
!    real,dimension(npts,nreac)    :: rrv
!    real,dimension(npts,njnum)    :: rjv
!    real,dimension(npts,ntrace)   :: vdv  ! deposition velocities
!    real,dimension(npts)          :: emino
!    real,dimension(npts,maxtrace) :: yv,y0v
!    integer                       :: iv,itrace,ivt,n
!
!    ! start
!
!    zoomed => region_dat(region)%zoomed
!
!    offsetl=0
!#ifdef MPI
!    if(myid>0) offsetl=sum(lmar(0:myid-1))
!#endif
!    dtime=nchem/(2*tref(region))  
!    !CMK  iterebi=max(1,nint(dtime/2400))  !needed if nchem <2400
!    iterebi=max(1,nint(dtime/1350))  !needed if nchem <2400
!    dt=dtime/iterebi
!
!    if ( okdebug ) print *, 'EBI: called with:',region,offsetl+level, &
!         'with dt:',dt,iterebi
!
!    !cmkif(level==1) then
!    !cmkio = sfstart('ebi.hdf',dfacc_create)
!    !cmkcall io_write3d_32(io,im(region),'im',jm(region),'jm', &
!    !  nreac,'nreac',rr,'rr')
!    !cmkcall io_write3d_32(io,im(region),'im',jm(region),'jm', &
!    !  njnum,'njnum',rj,'rj')
!    !cmkcall io_write3d_32(io,im(region),'im',jm(region),'jm', &
!    !  maxtrace,'maxtrace',y,'y')
!    !cmkcall io_write3d_32(io,im(region),'im',jm(region),'jm', &
!    !  n_extra,'n_extra',ye,'ye')
!    !cmkio = sfend(io)
!    !cmkend if
!
!    allocate(y0(im(region),jm(region),maxtrace))
!    allocate(cr2(im(region),jm(region),njnum))
!    allocate(cr3(im(region),jm(region),nreac))
!    allocate(cr4(im(region),jm(region),nreacw))
!
!    !*** SCALING OF NOx, which has changed due to transport/deposition
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          y(i,j,ino)  =max(1e-3,y(i,j,ino))
!          y(i,j,ino2) =max(1e-3,y(i,j,ino2))
!          y(i,j,ino3) =max(1e-3,y(i,j,ino3))
!          y(i,j,in2o5)=max(1e-3,y(i,j,in2o5))
!          y(i,j,ihno4)=max(1e-3,y(i,j,ihno4))   
!          fnoy=y(i,j,ino)+y(i,j,ino2)+y(i,j,ino3)+2.*y(i,j,in2o5)+y(i,j,ihno4)   
!          fnoy=y(i,j,inox)/fnoy
!          y(i,j,ino)  =fnoy*y(i,j,ino)
!          y(i,j,ino2) =fnoy*y(i,j,ino2)
!          y(i,j,ino3) =fnoy*y(i,j,ino3)
!          y(i,j,in2o5)=fnoy*y(i,j,in2o5)
!          y(i,j,ihno4)=fnoy*y(i,j,ihno4) 
!       end do
!    end do
!    !
!    ! set budget accumulators to zero
!    !
!    cr2=0.
!    cr3=0.
!    cr4=0.
!    !===========================================================
!    ! ** Start iterating over CHEMISTRY
!    !===========================================================
!    do ib=1,iterebi  
!       maxit=8   !CMKTEMP
!       if(offsetl+level<=3) maxit = maxit*2   ! lowest layers more iterations
!
!       y0 = y
!       !-------------------------------
!       ! wet sulphur/ammonia chemistry
!       !------------------------------
!       call wetS(region,level,y0,dt,y,ye,cr4)
!       !-------------------------------------
!       ! gasphase chemistry using EBI solver
!       !-------------------------------------
!
!       !cmk NOTE this statement was  there before but someone (me?) removed it.
!       y0 = y   
!
!       !cmk ______do EBI solver_______
!
!       !if (level  ==  19) then
!       !call dumpfield(0,'dumpb.hdf',rr,'rr')
!       !call dumpfield(1,'dumpb.hdf',rj,'rj')
!       !call dumpfield(1,'dumpb.hdf',y,'y')
!       !call dumpfield(1,'dumpb.hdf',ye,'ye')
!       !end if
!
!       dt2 = dt*dt
!       ! block the input for EBI for efficiency
!       ! copy values with faster running index in inside loop
!       iv = 0
!       do j=jsr(region),jer(region)
!          do i=isr(region),ier(region)
!             if(zoomed(i,j)/=region) cycle
!             iv = iv+1
!             ipts(iv) = i
!             jpts(iv) = j
!             if(iv==npts) then
!                ! copy reaction rates...
!                do itrace=1,nreac
!                   do ivt=1,npts
!                      rrv(ivt,itrace) = rr(ipts(ivt),jpts(ivt),itrace)
!                   end do
!                end do
!                ! copy photolysis rates....
!                do itrace=1,njnum
!                   do ivt=1,npts
!                      rjv(ivt,itrace) = rj(ipts(ivt),jpts(ivt),itrace)
!                   end do
!                end do
!                ! copy tracer concentrations ....
!                do itrace=1,maxtrace
!                   do ivt=1,npts
!                      yv(ivt,itrace) = y(ipts(ivt),jpts(ivt),itrace)
!                      y0v(ivt,itrace) = y0(ipts(ivt),jpts(ivt),itrace)
!                   end do
!                end do
!                ! deposition ....
!                if(offsetl + level == 1) then 
!                   do itrace=1,ntrace
!                      do ivt=1,npts
!                         vdv(ivt,itrace) = &
!                              vd(region,itrace)%surf(ipts(ivt),jpts(ivt)) &
!                              / ye(ipts(ivt),jpts(ivt),idz)   !1/s 
!                      end do
!                   end do
!                else
!                   vdv(:,:) = 0.0
!                end if
!                ! copy nox emissions....
!                do ivt=1,npts
!                   emino(ivt) = ye(ipts(ivt),jpts(ivt),ieno)
!                end do
!                ! do ebi solver....
!
!                call do_ebi(npts)    !the actual EBI solver on vectorized arrays  
!
!                do itrace=1,maxtrace
!                   do ivt=1,npts
!                      y(ipts(ivt),jpts(ivt),itrace)=yv(ivt,itrace)
!                   end do
!                end do
!                iv=0
!             end if
!          end do
!       end do
!       ! do the 'remaining' points...
!       if(iv > 0) then
!          do itrace=1,nreac
!             do ivt=1,iv
!                rrv(ivt,itrace) = rr(ipts(ivt),jpts(ivt),itrace)
!             end do
!          end do
!          do itrace=1,njnum
!             do ivt=1,iv
!                rjv(ivt,itrace) = rj(ipts(ivt),jpts(ivt),itrace)
!             end do
!          end do
!          do itrace=1,maxtrace
!             do ivt=1,iv
!                yv(ivt,itrace) = y(ipts(ivt),jpts(ivt),itrace)
!                y0v(ivt,itrace) = y0(ipts(ivt),jpts(ivt),itrace)
!             end do
!          end do
!          ! deposition ....
!          if(offsetl + level == 1) then 
!             do itrace=1,ntrace
!                do ivt=1,iv
!                   vdv(ivt,itrace) = &
!                        vd(region,itrace)%surf(ipts(ivt),jpts(ivt)) &
!                        / ye(ipts(ivt),jpts(ivt),idz)   !1/s 
!                end do
!             end do
!          else
!             vdv(:,:) = 0.0
!          end if
!          do ivt=1,iv
!             emino(ivt) = ye(ipts(ivt),jpts(ivt),ieno)
!          end do
!
!          call do_ebi(iv)    !the actual EBI solver on remaining cells  
!
!          do itrace=1,maxtrace
!             do ivt=1,iv
!                y(ipts(ivt),jpts(ivt),itrace)=yv(ivt,itrace)
!             end do
!          end do
!       end if
!
!       call NOymass
!
!       !-------------------------------------
!       ! marked tracers
!       ! apply after correction of nitrogen components
!       !----------------------------------------------
!
!       call mark_trac(region,level,y,rr,rj,dt,ye)
!
!       !------------------------------------------------------------
!       ! increase budget accumulators cr2 and cr3 (cr4 is done in wetS)
!       !------------------------------------------------------------
!
!       call incc2c3
!       !===========================================================
!       ! ** END iterating over CHEMISTRY
!       !===========================================================
!    end do   !iterebi
!
!    call reacbud   !add budgets for this timestep
!
!    deallocate(y0)
!    deallocate(cr2)
!    deallocate(cr3)
!    deallocate(cr4)
!
!    nullify(zoomed)
!
!  contains 
!
!    subroutine do_ebi(lvec)
!
!      integer, intent(in) :: lvec
!      integer             :: ivec
!
!      do ITER=1,MAXIT
!
!
!         do ivec=1,lvec
!            ! --- Short living compounds & groups
!            ! --- First group: NO NO2 O3
!            P1=rjv(ivec,jbno3)*yv(ivec,ino3)+emino(ivec)
!            R12=0.
!            R21=rrv(ivec,kho2no)*yv(ivec,iho2)+rrv(ivec,kmo2no)*yv(ivec,ich3o2)&
!                 +rrv(ivec,kc79)*yv(ivec,ixo2)+rrv(ivec,kc46)*yv(ivec,ic2o3)
!            XL1=rrv(ivec,knono3)*yv(ivec,ino3)+rrv(ivec,kc81)*yv(ivec,ixo2n)
!            XL1 = XL1 + vdv(ivec,ino)
!            P2=rjv(ivec,jhno3)*yv(ivec,ihno3)+rjv(ivec,jn2o5)*yv(ivec,in2o5)&
!                 +rrv(ivec,kn2o5)*yv(ivec,in2o5)+rjv(ivec,jano3)*yv(ivec,ino3)&
!                 +yv(ivec,ihno4)*(rjv(ivec,jhno4)+rrv(ivec,khno4m)+rrv(ivec,khno4oh)&
!                 *yv(ivec,ioh))+2.*rrv(ivec,knono3)*yv(ivec,ino3)*yv(ivec,ino)&
!                 +rrv(ivec,kc48)*yv(ivec,ipan)+rrv(ivec,kc59)*yv(ivec,iole)*yv(ivec,ino3)&
!                 +rrv(ivec,kc84)*yv(ivec,iorgntr)*yv(ivec,ioh)+rjv(ivec,jorgn)&
!                 *yv(ivec,iorgntr)+rjv(ivec,jpan)*yv(ivec,ipan)+0.1*rrv(ivec,kc78) *yv(ivec,iisop)*yv(ivec,ino3)
!            XL2=rrv(ivec,kno2oh)*yv(ivec,ioh)+rrv(ivec,kno2no3)*yv(ivec,ino3)&
!                 +rrv(ivec,kno2ho2)*yv(ivec,iho2)+rrv(ivec,kno2o3)*yv(ivec,io3)+rrv(ivec,kc47)*yv(ivec,ic2o3)
!            XL2 = XL2 + vdv(ivec,ino2)
!            P3=rjv(ivec,jano3)*yv(ivec,ino3)
!            XL3=rrv(ivec,ko3ho2)*yv(ivec,iho2)+rrv(ivec,ko3oh)*yv(ivec,ioh)+rrv(ivec,kno2o3)&
!                 *yv(ivec,ino2)+rjv(ivec,jo3d)+rrv(ivec,kc58)*yv(ivec,iole)&
!                 +rrv(ivec,kc62)*yv(ivec,ieth)+rrv(ivec,kc77)*yv(ivec,iisop)
!            XL3 = XL3 + vdv(ivec,io3)
!
!            X1=y0v(ivec,ino)+P1*DT
!            X2=y0v(ivec,ino2)+P2*DT
!            X3=y0v(ivec,io3)+P3*DT
!            C1=1.+XL1*DT
!            C2=1.+XL2*DT
!            C3=1.+XL3*DT
!            Y1=rrv(ivec,knoo3)*DT
!            R21T=R21*DT
!            R12T=R12*DT
!            XJT=rjv(ivec,jno2)*DT
!            ! --- Oplossen onbekende x
!            ACUB=-2.*Y1*(C2+R12T+C2*R21T/C1)
!            BCUB=2.*C1*C2*C3+2.*C1*C3*(R12T+XJT)+2.*C2*C3*R21T+&
!                 2.*Y1*(R12T*(X1-X2)+2.*C2*R21T*X1/C1+C2*(X1+X3))
!            CCUB=2.*C1*C3*X2*(R12T+XJT)-2.*C2*C3*X1*R21T+2.*Y1*X1*&
!                 (X2*R12T-C2*X3-C2*R21T*X1/C1)
!            CUBDET=BCUB*BCUB-4.*ACUB*CCUB
!            DNO2=(-1.*BCUB+SQRT(CUBDET))/(2.*ACUB)
!            dno2=min(x1,dno2)
!            yv(ivec,ino2)=(X2+DNO2)/C2
!            yv(ivec,ino)=(X1-DNO2)/C1
!            yv(ivec,io3)=(X3+XJT*yv(ivec,ino2))/(C3+Y1*yv(ivec,ino))
!            ! --- Second group: yv(ivec,iho2) yv(ivec,ioh) yv(ivec,ihno4)
!            R57=rjv(ivec,jhno4)+rrv(ivec,khno4m)
!            R56=rrv(ivec,kcooh)*yv(ivec,ico)+rrv(ivec,ko3oh)*yv(ivec,io3)+rrv(ivec,khpoh)&
!                 *yv(ivec,ih2o2)+rrv(ivec,kfrmoh)*yv(ivec,ich2o)+rrv(ivec,kh2oh)
!            P5=2.*rjv(ivec,jbch2o)*yv(ivec,ich2o)+rrv(ivec,kc46)*yv(ivec,ic2o3)*yv(ivec,ino)&
!                 +rrv(ivec,kmo2no)*yv(ivec,ich3o2)*yv(ivec,ino)+0.66*rrv(ivec,kmo2mo2)&
!                 *yv(ivec,ich3o2)*yv(ivec,ich3o2)+rjv(ivec,jmepe)*yv(ivec,ich3o2h)&
!                 +2.*rrv(ivec,kc49)*yv(ivec,ic2o3)*yv(ivec,ic2o3)+2.*rjv(ivec,j45)&
!                 *yv(ivec,iald2)+rjv(ivec,j74)*yv(ivec,imgly)+0.11*rrv(ivec,kc52)*yv(ivec,ipar)&
!                 *yv(ivec,ioh)+0.94*rrv(ivec,kc53)*yv(ivec,iror)+rrv(ivec,kc54)*yv(ivec,iror)&
!                 +rrv(ivec,kc57)*yv(ivec,iole)*yv(ivec,ioh)+0.25*rrv(ivec,kc58)*yv(ivec,io3)&
!                 *yv(ivec,iole)+rrv(ivec,kc61)*yv(ivec,ieth)*yv(ivec,ioh)+0.26*rrv(ivec,kc62)&
!                 *yv(ivec,ieth)*yv(ivec,io3)+0.85*rrv(ivec,kc76)*yv(ivec,iisop)*yv(ivec,ioh)&
!                 +0.3*rrv(ivec,kc77)*yv(ivec,iisop)*yv(ivec,io3)+rrv(ivec,kc41)*yv(ivec,ich2o)&
!                 *yv(ivec,ino3)+rjv(ivec,jorgn)*yv(ivec,iorgntr)+0.9*rrv(ivec,kc78)*yv(ivec,iisop)*yv(ivec,ino3)
!            XL5=rrv(ivec,kho2no)*yv(ivec,ino)+rrv(ivec,kno2ho2)*yv(ivec,ino2)&
!                 +rrv(ivec,ko3ho2)*yv(ivec,io3)+rrv(ivec,kmo2ho2)*yv(ivec,ich3o2)&
!                 +rrv(ivec,kho2oh)*yv(ivec,ioh)+rrv(ivec,kc82)*yv(ivec,ixo2)+rrv(ivec,kc85)*yv(ivec,ixo2n)
!            R66=2.*rrv(ivec,kho2ho2)
!            X5=y0v(ivec,iho2)+P5*DT
!            R65=rrv(ivec,kho2no)*yv(ivec,ino)+rrv(ivec,ko3ho2)*yv(ivec,io3)
!            P6=rjv(ivec,jhno3)*yv(ivec,ihno3)+2.*rjv(ivec,jo3d)*yv(ivec,io3)&
!                 +2.*rjv(ivec,jh2o2)*yv(ivec,ih2o2)+rjv(ivec,jmepe)*yv(ivec,ich3o2h)&
!                 +0.79*rrv(ivec,kc50)*yv(ivec,ic2o3)*yv(ivec,iho2)+0.4*rrv(ivec,kc58)&
!                 *yv(ivec,io3)*yv(ivec,iole)+0.28*rrv(ivec,kc77)*yv(ivec,iisop)*yv(ivec,io3)&
!                 +rjv(ivec,jrooh)*yv(ivec,irooh)+0.12*rrv(ivec,kc62)*yv(ivec,ieth)*yv(ivec,io3)
!            XL6=rrv(ivec,khno4oh)*yv(ivec,ihno4)+rrv(ivec,kho2oh)*yv(ivec,iho2)&
!                 +rrv(ivec,kno2oh)*yv(ivec,ino2)+rrv(ivec,kohhno3)*yv(ivec,ihno3)&
!                 +rrv(ivec,kcooh)*yv(ivec,ico)+rrv(ivec,ko3oh)*yv(ivec,io3)+rrv(ivec,khpoh)&
!                 *yv(ivec,ih2o2)+rrv(ivec,kfrmoh)*yv(ivec,ich2o)+rrv(ivec,kch4oh)&
!                 *yv(ivec,ich4)+0.7*rrv(ivec,kohmper)*yv(ivec,ich3o2h)+rrv(ivec,kc43)&
!                 *yv(ivec,iald2)+rrv(ivec,kc73)*yv(ivec,imgly)+rrv(ivec,kc52)&
!                 *yv(ivec,ipar)+rrv(ivec,kc57)*yv(ivec,iole)+rrv(ivec,kc61)*yv(ivec,ieth)&
!                 +rrv(ivec,kc76)*yv(ivec,iisop)+0.7*rrv(ivec,kohrooh)*yv(ivec,irooh)&
!                 +rrv(ivec,kc84)*yv(ivec,iorgntr)+rrv(ivec,kh2oh)&
!                 +(rrv(ivec,kdmsoha)+rrv(ivec,kdmsohb)) *yv(ivec,idms) +rrv(ivec,knh3oh)*yv(ivec,inh3)  !sulfur
!
!            X6=y0v(ivec,ioh)+P6*DT
!            C6=1.+XL6*DT
!            R75=rrv(ivec,kno2ho2)*yv(ivec,ino2)
!            XL7=rjv(ivec,jhno4)+rrv(ivec,khno4oh)*yv(ivec,ioh)+rrv(ivec,khno4m)
!            XL7 = XL7 + vdv(ivec,ihno4)
!            C7=1.+XL7*DT
!            Y1=R57/C7
!            Y2=R56/C6
!            ACUB=R66*DT
!            BCUB=1.+XL5*DT-DT2*(Y1*R75+Y2*R65)
!            CCUB=-1.*X5-DT*(Y1*y0v(ivec,ihno4)+Y2*X6)
!            CUBDET=BCUB*BCUB-4.*ACUB*CCUB
!            CUBDET=MAX(CUBDET,1.E-20)
!            yv(ivec,iho2)=max(0.1,(-1.*BCUB+SQRT(CUBDET))/(2.*ACUB))
!            yv(ivec,ioh)=(X6+R65*yv(ivec,iho2)*DT)/C6
!            yv(ivec,ihno4)=(y0v(ivec,ihno4)+R75*DT*yv(ivec,iho2))/C7
!            ! --- Third group: NO3 N2O5
!            R89=rjv(ivec,jn2o5)+rrv(ivec,kn2o5)
!            P8=rrv(ivec,kohhno3)*yv(ivec,ihno3)*yv(ivec,ioh)+rrv(ivec,kno2o3)*yv(ivec,ino2)*yv(ivec,io3)
!            XL8=rjv(ivec,jbno3)+rjv(ivec,jano3)+rrv(ivec,knono3)*yv(ivec,ino)&
!                 +rrv(ivec,kno2no3)*yv(ivec,ino2)+rrv(ivec,kc44)*yv(ivec,iald2)+rrv(ivec,kc59)&
!                 *yv(ivec,iole)+rrv(ivec,kc78)*yv(ivec,iisop)+rrv(ivec,kc41)*yv(ivec,ich2o)+rrv(ivec,kdmsno3)*yv(ivec,idms)
!            XL8 = XL8 + vdv(ivec,ino3) 
!            X4=y0v(ivec,ino3)+P8*DT
!            C5=1.+XL8*DT
!            R98=rrv(ivec,kno2no3)*yv(ivec,ino2)
!            XL9=rjv(ivec,jn2o5)+rrv(ivec,kn2o5)+rrv(ivec,kn2o5aq)+rrv(ivec,kn2o5l)   !cmk rates now idependent from y
!            XL9 = XL9 + vdv(ivec,in2o5)
!            C6=1.+XL9*DT
!            C7=(C5*C6-R89*R98*DT2)
!            yv(ivec,in2o5)=(C5*y0v(ivec,in2o5)+R98*DT*X4)/C7
!            yv(ivec,ino3)=(C6*X4+R89*DT*y0v(ivec,in2o5))/C7
!            ! --- Fourth group: C2O3 PAN
!            R1920=rrv(ivec,kc48)+rjv(ivec,jpan)
!            R1919=rrv(ivec,kc49)
!            P19=rrv(ivec,kc43)*yv(ivec,iald2)*yv(ivec,ioh)+rrv(ivec,kc44)*yv(ivec,iald2)&
!                 *yv(ivec,ino3)+rrv(ivec,kc73)*yv(ivec,imgly)*yv(ivec,ioh)+rjv(ivec,j74)&
!                 *yv(ivec,imgly)+0.15*rrv(ivec,kc77)*yv(ivec,iisop)*yv(ivec,io3)
!            XL19=rrv(ivec,kc46)*yv(ivec,ino)+rrv(ivec,kc50)*yv(ivec,iho2)+rrv(ivec,kc47)*yv(ivec,ino2)
!            XL19 = XL19 + vdv(ivec,ic2o3)
!            R2019=rrv(ivec,kc47)*yv(ivec,ino2)
!            XL20=rrv(ivec,kc48)+rjv(ivec,jpan)
!            XL20 = XL20 + vdv(ivec,ipan)
!            ACUB=2*R1919*DT*(1+XL20*DT)
!            BCUB=(1+XL20*DT)*(1+XL19*DT)-R1920*DT*R2019*DT
!            CCUB=(1+XL20*DT)*(y0v(ivec,ic2o3)+P19*DT)+R1920*DT*y0v(ivec,ipan)
!            CUBDET=BCUB*BCUB+4.*ACUB*CCUB
!            yv(ivec,ic2o3)=max(1e-8,(-1.*BCUB+SQRT(CUBDET))/(2.*ACUB))    !cmk  put max here....
!            yv(ivec,ipan)=(y0v(ivec,ipan)+R2019*yv(ivec,ic2o3)*DT)/(1.+XL20*DT)
!            ! --- CH4 chemistry (short living radicals)
!            PCH3O2=rrv(ivec,kch4oh)*yv(ivec,ich4)*yv(ivec,ioh)+0.7*rrv(ivec,kohmper)*yv(ivec,ioh)*yv(ivec,ich3o2h)
!            XLCH3O2=rrv(ivec,kmo2no)*yv(ivec,ino)+rrv(ivec,kmo2ho2)*yv(ivec,iho2)&
!                 +2*rrv(ivec,kmo2mo2)*yv(ivec,ich3o2)
!            yv(ivec,ich3o2)=(y0v(ivec,ich3o2)+PCH3O2*DT)/(1.+XLCH3O2*DT)
!            ! --- CBM4 chem.(short living compounds & operators)
!            PRXPAR=0.11*rrv(ivec,kc52)*yv(ivec,ioh)*yv(ivec,ipar)+2.1*rrv(ivec,kc53)&
!                 *yv(ivec,iror)+rrv(ivec,kc57)*yv(ivec,iole)*yv(ivec,ioh)+0.9*rrv(ivec,kc58)&
!                 *yv(ivec,io3)*yv(ivec,iole)+rrv(ivec,kc59)*yv(ivec,iole)*yv(ivec,ino3)
!            XLRXPAR=rrv(ivec,kc83)*yv(ivec,ipar)
!            yv(ivec,irxpar)=(y0v(ivec,irxpar)+PRXPAR*DT)/(1.+XLRXPAR*DT)
!            XLISOP=rrv(ivec,kc76)*yv(ivec,ioh)+rrv(ivec,kc77)*yv(ivec,io3)+rrv(ivec,kc78)*yv(ivec,ino3)
!            yv(ivec,iisop)=y0v(ivec,iisop)/(1.+XLISOP*DT)
!            PROR=0.76*rrv(ivec,kc52)*yv(ivec,ipar)*yv(ivec,ioh)+0.02*rrv(ivec,kc53)*yv(ivec,iror)
!            XLROR=rrv(ivec,kc53)+rrv(ivec,kc54)
!            yv(ivec,iror)=(y0v(ivec,iror)+PROR*DT)/(1.+XLROR*DT)
!            PXO2=rrv(ivec,kc46)*yv(ivec,ic2o3)*yv(ivec,ino)+2.*rrv(ivec,kc49)&
!                 *yv(ivec,ic2o3)*yv(ivec,ic2o3)+rrv(ivec,kc50)*yv(ivec,ic2o3)&
!                 *yv(ivec,iho2)+rjv(ivec,j45)*yv(ivec,iald2)+rrv(ivec,kc73)&
!                 *yv(ivec,imgly)*yv(ivec,ioh)+0.87*rrv(ivec,kc52)*yv(ivec,ipar)*yv(ivec,ioh)&
!                 +0.96*rrv(ivec,kc53)*yv(ivec,iror)+rrv(ivec,kc57)*yv(ivec,iole)*yv(ivec,ioh)&
!                 +0.29*rrv(ivec,kc58)*yv(ivec,io3)*yv(ivec,iole)+0.91*rrv(ivec,kc59)&
!                 *yv(ivec,iole)*yv(ivec,ino3)+rrv(ivec,kc61)*yv(ivec,ieth)*yv(ivec,ioh)&
!                 +0.85*rrv(ivec,kc76)*yv(ivec,iisop)*yv(ivec,ioh)+0.7*rrv(ivec,kohrooh)&
!                 *yv(ivec,irooh)*yv(ivec,ioh)+rrv(ivec,kc84)*yv(ivec,ioh)*yv(ivec,iorgntr)&
!                 +0.18*rrv(ivec,kc77)*yv(ivec,iisop)*yv(ivec,io3)
!            XLXO2=rrv(ivec,kc79)*yv(ivec,ino)+2.*rrv(ivec,kc80)*yv(ivec,ixo2)+rrv(ivec,kc82)*yv(ivec,iho2)
!            yv(ivec,ixo2)=(y0v(ivec,ixo2)+PXO2*DT)/(1.+XLXO2*DT)
!            PXO2N=0.13*rrv(ivec,kc52)*yv(ivec,ipar)*yv(ivec,ioh)+0.04*rrv(ivec,kc53)&
!                 *yv(ivec,iror)+0.09*rrv(ivec,kc59)*yv(ivec,iole)*yv(ivec,ino3)+0.15*rrv(ivec,kc76)*yv(ivec,iisop)*yv(ivec,ioh)
!            XLXO2N=rrv(ivec,kc81)*yv(ivec,ino)+rrv(ivec,kc85)*yv(ivec,iho2)
!            yv(ivec,ixo2n)=(y0v(ivec,ixo2n)+PXO2N*DT)/(1.+XLXO2N*DT)
!         end do !ivec
!
!         if ( mod(iter,2) == 0 ) then
!
!            do ivec=1,lvec
!               ! --- Species with intermediate lifetimes
!               ! --- Inorganic compounds (HNO3 H2O2)
!               ! 
!               PHNO3=rrv(ivec,kno2oh)*yv(ivec,ino2)*yv(ivec,ioh)+2.*(rrv(ivec,kn2o5aq)+rrv(ivec,kn2o5l))&
!                    *yv(ivec,in2o5)+rrv(ivec,kc44)*yv(ivec,iald2)*yv(ivec,ino3)+rrv(ivec,kc41)*yv(ivec,ich2o)*yv(ivec,ino3)
!               XLHNO3=rjv(ivec,jhno3)+rrv(ivec,kohhno3)*yv(ivec,ioh)
!               XLHNO3=XLHNO3 + vdv(ivec,ihno3)
!               yv(ivec,ihno3)=(y0v(ivec,ihno3)+PHNO3*DT)/(1.+XLHNO3*DT)
!               PH2O2=rrv(ivec,kho2ho2)*yv(ivec,iho2)*yv(ivec,iho2)
!               XLH2O2=rjv(ivec,jh2o2)+rrv(ivec,khpoh)*yv(ivec,ioh)
!               XLH2O2=XLH2O2 + vdv(ivec,ih2o2)
!               yv(ivec,ih2o2)=(y0v(ivec,ih2o2)+PH2O2*DT)/(1.+XLH2O2*DT)
!               ! --- CH4-chemistry (methyl peroxide formaldehyde)
!               PCH3O2H=rrv(ivec,kmo2ho2)*yv(ivec,ich3o2)*yv(ivec,iho2)
!               XLCH3O2H=rrv(ivec,kohmper)*yv(ivec,ioh)+rjv(ivec,jmepe)
!               XLCH3O2H=XLCH3O2H + vdv(ivec,ich3o2h)
!               yv(ivec,ich3o2h)=(y0v(ivec,ich3o2h)+PCH3O2H*DT)/(1.+XLCH3O2H*DT)
!               PCH2O=rrv(ivec,kc46)*yv(ivec,ic2o3)*yv(ivec,ino)+0.3*rrv(ivec,kohmper)&
!                    *yv(ivec,ich3o2h)*yv(ivec,ioh)+rrv(ivec,kmo2no)*yv(ivec,ich3o2)*yv(ivec,ino)&
!                    +1.33*rrv(ivec,kmo2mo2)*yv(ivec,ich3o2)*yv(ivec,ich3o2)+rjv(ivec,jmepe)&
!                    *yv(ivec,ich3o2h)+2.*rrv(ivec,kc49)*yv(ivec,ic2o3)*yv(ivec,ic2o3)&
!                    +rrv(ivec,kc50)*yv(ivec,ic2o3)*yv(ivec,iho2)+rjv(ivec,j45)*yv(ivec,iald2)&
!                    +rrv(ivec,kc57)*yv(ivec,iole)*yv(ivec,ioh)+0.64*rrv(ivec,kc58)*yv(ivec,iole)&
!                    *yv(ivec,io3)+rrv(ivec,kc59)*yv(ivec,iole)*yv(ivec,ino3)+1.56*rrv(ivec,kc61)&
!                    *yv(ivec,ieth)*yv(ivec,ioh)+rrv(ivec,kc62)*yv(ivec,ieth)*yv(ivec,io3)&
!                    +0.61*rrv(ivec,kc76)*yv(ivec,iisop)*yv(ivec,ioh)+0.9*rrv(ivec,kc77)&
!                    *yv(ivec,iisop)*yv(ivec,io3)+0.03*rrv(ivec,kc78)*yv(ivec,iisop)*yv(ivec,ino3)
!               XLCH2O=rjv(ivec,jach2o)+rjv(ivec,jbch2o)+yv(ivec,ioh)*rrv(ivec,kfrmoh)+rrv(ivec,kc41)*yv(ivec,ino3)
!               XLCH2O=XLCH2O + vdv(ivec,ich2o)
!               yv(ivec,ich2o)=(y0v(ivec,ich2o)+PCH2O*DT)/(1.+XLCH2O*DT)
!               ! --- CBIV-elements for higher HC-chemistry: ALD2 MGLY
!               ! --- ETH OLE ISOP ROOH ORGNTR
!               PALD2=0.11*rrv(ivec,kc52)*yv(ivec,ipar)*yv(ivec,ioh)+1.1*rrv(ivec,kc53)&
!                    *yv(ivec,iror)+rrv(ivec,kc57)*yv(ivec,iole)*yv(ivec,ioh)+0.44*rrv(ivec,kc58)&
!                    *yv(ivec,iole)*yv(ivec,io3)+rrv(ivec,kc59)*yv(ivec,iole)*yv(ivec,ino3)&
!                    +0.22*rrv(ivec,kc61)*yv(ivec,ieth)*yv(ivec,ioh)+0.12*rrv(ivec,kc78)*yv(ivec,iisop)*yv(ivec,ino3)
!               XLALD2=rrv(ivec,kc43)*yv(ivec,ioh)+rrv(ivec,kc44)*yv(ivec,ino3)+rjv(ivec,j45)
!               XLALD2=XLALD2 + vdv(ivec,iald2)
!               yv(ivec,iald2)=(y0v(ivec,iald2)+PALD2*DT)/(1.+XLALD2*DT)
!               PMGLY=0.03*rrv(ivec,kc76)*yv(ivec,iisop)*yv(ivec,ioh)+0.03*rrv(ivec,kc77)&
!                    *yv(ivec,iisop)*yv(ivec,io3)+0.08*rrv(ivec,kc78)*yv(ivec,iisop)*yv(ivec,ino3)
!               XLMGLY=rrv(ivec,kc73)*yv(ivec,ioh)+rjv(ivec,j74)
!               yv(ivec,imgly)=(y0v(ivec,imgly)+PMGLY*DT)/(1.+XLMGLY*DT)
!               XLETH=rrv(ivec,kc61)*yv(ivec,ioh)+rrv(ivec,kc62)*yv(ivec,io3)
!               yv(ivec,ieth)=y0v(ivec,ieth)/(1.+XLETH*DT)
!               POLE=0.58*rrv(ivec,kc76)*yv(ivec,iisop)*yv(ivec,ioh)+0.55*rrv(ivec,kc77)&
!                    *yv(ivec,iisop)*yv(ivec,io3)+0.45*rrv(ivec,kc78)*yv(ivec,iisop) *yv(ivec,ino3)
!               XLOLE=rrv(ivec,kc57)*yv(ivec,ioh)+rrv(ivec,kc58)*yv(ivec,io3)+rrv(ivec,kc59)*yv(ivec,ino3)
!               yv(ivec,iole)=(y0v(ivec,iole)+POLE*DT)/(1.+XLOLE*DT)
!               PROOH=rrv(ivec,kc82)*yv(ivec,ixo2)*yv(ivec,iho2)+0.21*rrv(ivec,kc50)&
!                    *yv(ivec,ic2o3)*yv(ivec,iho2)+rrv(ivec,kc85)*yv(ivec,iho2)*yv(ivec,ixo2n)
!               XLROOH=rjv(ivec,jrooh)+rrv(ivec,kohrooh)*yv(ivec,ioh)
!               XLROOH = XLROOH + vdv(ivec,irooh)
!               yv(ivec,irooh)=(y0v(ivec,irooh)+PROOH*DT)/(1.+XLROOH*DT)
!
!               PORGNTR=rrv(ivec,kc81)*yv(ivec,ino)*yv(ivec,ixo2n)+0.9*rrv(ivec,kc78)*yv(ivec,iisop)*yv(ivec,ino3)
!               XLORGNTR=rrv(ivec,kc84)*yv(ivec,ioh)+rjv(ivec,jorgn)
!               XLORGNTR=XLORGNTR+vdv(ivec,iorgntr)
!
!               yv(ivec,iorgntr)=(y0v(ivec,iorgntr)+PORGNTR*DT)/(1.+XLORGNTR*DT)
!
!               ! gas phase sulfur   & ammonia
!
!               qdms1=rrv(ivec,kdmsoha)*yv(ivec,ioh)+rrv(ivec,kdmsno3)*yv(ivec,ino3)
!               qdms2=rrv(ivec,kdmsohb)*yv(ivec,ioh)
!               qdms=qdms1+qdms2
!               yv(ivec,idms)=y0v(ivec,idms)/(1.+qdms*DT)
!               pso2=yv(ivec,idms)*(qdms1+0.75*qdms2)
!               pmsa=yv(ivec,idms)*0.25*qdms2
!               qso2=rrv(ivec,kso2oh)*yv(ivec,ioh)
!               qso2d=qso2 + vdv(ivec,iso2)
!               yv(ivec,iso2)=(y0v(ivec,iso2)+pso2*DT) /(1.+qso2d*DT)  !qso2d includes deposition
!               yv(ivec,imsa)=(y0v(ivec,imsa)+pmsa*DT) /(1.+vdv(ivec,imsa)*DT)  
!               yv(ivec,iso4)=(y0v(ivec,iso4)+qso2*yv(ivec,iso2)*DT) /(1. + vdv(ivec,iso4)*DT)   !corredted CMK qso2/qso2d
!               yv(ivec,iacid)=(y0v(ivec,iacid)+(pmsa+2.*qso2*yv(ivec,iso2))*DT)/&
!                    (1.+rrv(ivec,knh3SO4)*yv(ivec,inh3)*DT)
!               yv(ivec,inh4)=(y0v(ivec,inh4)+yv(ivec,iacid)*rrv(ivec,knh3SO4)*yv(ivec,inh3)*DT)/(1.+vdv(ivec,inh4)*DT)
!               pnh2=yv(ivec,ioh)*rrv(ivec,knh3oh)
!               qnh3=yv(ivec,iacid)*rrv(ivec,knh3SO4)+pnh2
!               qnh3 = qnh3 + vdv(ivec,inh3)
!               yv(ivec,inh3)=y0v(ivec,inh3)/(1.+qnh3*DT)
!               qnh2= rrv(ivec,knh2no)*yv(ivec,ino)+rrv(ivec,knh2no2)*yv(ivec,ino2)&
!                    +rrv(ivec,knh2ho2)*yv(ivec,iho2) +rrv(ivec,knh2o2)+rrv(ivec,knh2o3)*yv(ivec,io3)
!               yv(ivec,inh2)=(y0v(ivec,inh2)+yv(ivec,inh3)*pnh2*dt)/(1.+qnh2*dt)
!            end do  !ivec
!
!         end if
!
!         if ( mod(iter,maxit) == 0 ) then
!
!            ! --- Long living compounds
!            do ivec=1,lvec
!               yv(ivec,ich4)=y0v(ivec,ich4)/(1.+rrv(ivec,kch4oh)*yv(ivec,ioh)*DT)
!               PCO=yv(ivec,ich2o)*(rjv(ivec,jach2o)+rjv(ivec,jbch2o)+yv(ivec,ioh)&
!                    *rrv(ivec,kfrmoh))+rjv(ivec,j45)*yv(ivec,iald2)+rjv(ivec,j74)*yv(ivec,imgly)&
!                    +0.37*rrv(ivec,kc58)*yv(ivec,iole)*yv(ivec,io3)+0.43*rrv(ivec,kc62)&
!                    *yv(ivec,ieth)*yv(ivec,io3)+0.36*rrv(ivec,kc77)*yv(ivec,iisop)*yv(ivec,io3)&
!                    +rrv(ivec,kc41)*yv(ivec,ich2o)*yv(ivec,ino3)
!               XLCO = rrv(ivec,kcooh)*yv(ivec,ioh)
!               XLCO = XLCO + vdv(ivec,ico)
!               yv(ivec,ico)=(y0v(ivec,ico)+PCO*DT)/(1.+XLCO*DT)
!               PPAR=0.63*rrv(ivec,kc77)*yv(ivec,iisop)*yv(ivec,io3)+0.63*rrv(ivec,kc76)*yv(ivec,iisop)*yv(ivec,ioh)
!               XLPAR=rrv(ivec,kc52)*yv(ivec,ioh)+rrv(ivec,kc83)*yv(ivec,irxpar)
!               yv(ivec,ipar)=(y0v(ivec,ipar)+PPAR*DT)/(1.+XLPAR*DT)
!               !cmk ____added rn222 chemistry in EBI language
!               yv(ivec,irn222) = y0v(ivec,irn222)/(1.+rrv(ivec,krn222)*dt)
!               yv(ivec,ipb210) = y0v(ivec,ipb210)+y0v(ivec,irn222)-yv(ivec,irn222)
!               !if(yv(ivec,ipb210) < 0.0 .or. yv(ivec,irn222) < 0.0) then
!               !   print *, 'Negatives .....rn222, pb210', y0v(ivec,irn222), yv(ivec,irn222) , y0v(ivec,ipb210),  yv(ivec,ipb210)
!               !end if
!            end do   !ivec
!
!         end if
!
!      end do !ITER
!
!    end subroutine do_ebi
!
!
!    subroutine NOYmass
!#ifdef MPI    
!      use mpi_comm,only: stopmpi
!#endif
!      implicit none
!      integer i,j,imax, offsetl
!      real :: ncormax,ncorav,totn,totn0,fnoy,fnoy1
!      real :: ncorr,ncorr1,ncorr2,ncorr3, totdep
!      logical :: nerr
!
!      offsetl=0
!#ifdef MPI
!      if(myid>0) offsetl=sum(lmar(0:myid-1))
!#endif
!
!      ncormax=0.
!      ncorav=0.
!      nerr=.false.
!      imax = 0
!      do j=jsr(region),jer(region)
!         do i=isr(region),ier(region)
!            if(zoomed(i,j)/=region) cycle
!            imax = imax + 1 
!            !
!            !** Guarantee exact mass conservation of NOY 
!            !   (this may matter a few percent)
!            ! 
!            fnoy=y(i,j,ino)+y(i,j,ino2)+y(i,j,ino3)+2.*y(i,j,in2o5)+y(i,j,ihno4)
!            if (level+offsetl == 1) then
!               totdep = (y(i,j,ino) *vd(region,ino )%surf(i,j)  + &
!                    y(i,j,ino2)*vd(region,ino2)%surf(i,j)  + &
!                    y(i,j,ino3)*vd(region,ino3)%surf(i,j)  + &
!                    y(i,j,ihno3)*vd(region,ihno3)%surf(i,j)  + &
!                    y(i,j,ipan)*vd(region,ipan)%surf(i,j)  + &
!                    y(i,j,iorgntr)*vd(region,iorgntr)%surf(i,j)  + &
!                    2*y(i,j,in2o5)*vd(region,in2o5)%surf(i,j)  + &
!                    y(i,j,ihno4)*vd(region,ihno4)%surf(i,j) )*dt/ye(i,j,idz)  
!            else
!               totdep  = 0.0 
!            end if
!            totn0=y0(i,j,inox)+y0(i,j,ihno3)+y0(i,j,ipan)+ &
!                 y0(i,j,iorgntr) + ye(i,j,ieno)*dt -  totdep  
!            ! note that emino is added here and the total deposition is subtracted
!            !
!            ! totn0 contains all nitrogen at beginning of timestep + nox emissions
!            !
!            !
!            ! totn contains all nitrogen at end of timestep
!            !
!            totn=fnoy+y(i,j,ihno3)+y(i,j,ipan)+y(i,j,iorgntr)
!            ! correction factor for all nitrogen compounds
!            ncorr=totn-totn0
!
!            if ( totn < tiny(totn) ) cycle
!
!            if ( (abs(ncorr)/totn) > 0.05 ) then    !CMK changed from 0.1 to 0.05
!
!               nerr=.true.
!               print *,'NOYmass: N-error....',region,offsetl+level,i,j,totn0,totn
!               print *,'NOYmass: emino    ',ye(i,j,ieno)*dt/y0(i,j,iair)*1e9
!               print *,'NOYmass: NO(0)    ', &
!                    y0(i,j,ino)/y0(i,j,iair)*1e9,y(i,j,ino)/y(i,j,iair)*1e9
!               print *,'NOYmass: NO2(0)   ', &
!                    y0(i,j,ino2)/y0(i,j,iair)*1e9,y(i,j,ino2)/y(i,j,iair)*1e9
!               print *,'NOYmass: O3(0)    ', &
!                    y0(i,j,io3)/y0(i,j,iair)*1e9,y(i,j,io3)/y(i,j,iair)*1e9
!               print *,'NOYmass: NO3(0)   ', &
!                    y0(i,j,ino3)/y0(i,j,iair)*1e9,y(i,j,ino3)/y(i,j,iair)*1e9
!               print *,'NOYmass: N2O5(0)  ', &
!                    y0(i,j,in2o5)/y0(i,j,iair)*1e9,y(i,j,in2o5)/y(i,j,iair)*1e9
!               print *,'NOYmass: HNO4(0)  ', &
!                    y0(i,j,ihno4)/y0(i,j,iair)*1e9,y(i,j,ihno4)/y(i,j,iair)*1e9
!               print *,'NOYmass: HNO3(0)  ', &
!                    y0(i,j,ihno3)/y0(i,j,iair)*1e9,y(i,j,ihno3)/y(i,j,iair)*1e9
!               print *,'NOYmass: PAN(0)   ', &
!                    y0(i,j,ipan)/y0(i,j,iair)*1e9,y(i,j,ipan)/y(i,j,iair)*1e9
!               print *,'NOYmass: ORGNT(0) ', &
!                    y0(i,j,iorgntr)/y0(i,j,iair)*1e9,y(i,j,iorgntr)/y(i,j,iair)*1e9
!               print *,'NOYmass: NOx(0)   ', &
!                    y0(i,j,inox)/y0(i,j,iair)*1e9,y(i,j,inox)/y(i,j,iair)*1e9
!               print *,'NOYmass: ',rj(i,j,jhno3),rr(i,j,kohhno3)*y(i,j,ioh), &
!                    y(i,j,ioh)/y(i,j,iair)*1e9
!
!            end if
!            ! maximum and average correction factor in this loop
!            ncormax=max(abs(ncormax),abs(ncorr/totn)) 
!            ncorav=ncorav+abs(ncorr/totn)
!            !
!            ! first correct hno3, pan and organic nitrates 
!            ! (as a group of reservoir tracers)
!            ! 
!            totn=y(i,j,ihno3)+y(i,j,ipan)+y(i,j,iorgntr)
!            if ( totn < tiny(totn) ) cycle
!            ncorr1=y(i,j,ihno3)  *(1.-ncorr/totn)
!            ncorr2=y(i,j,ipan)   *(1.-ncorr/totn)
!            ncorr3=y(i,j,iorgntr)*(1.-ncorr/totn)
!            y(i,j,ihno3)  =max(0.,ncorr1)
!            y(i,j,ipan)   =max(0.,ncorr2)
!            y(i,j,iorgntr)=max(0.,ncorr3)
!            ncorr=ncorr1+ncorr2+ncorr3-y(i,j,ihno3)-y(i,j,ipan)- y(i,j,iorgntr)
!            !
!            ! the remainder is used to scale the noy components
!            !
!            fnoy1=(fnoy+ncorr)/fnoy
!            y(i,j,ino)  =fnoy1*y(i,j,ino)
!            y(i,j,ino2) =fnoy1*y(i,j,ino2)
!            y(i,j,ino3) =fnoy1*y(i,j,ino3)
!            y(i,j,in2o5)=fnoy1*y(i,j,in2o5)
!            y(i,j,ihno4)=fnoy1*y(i,j,ihno4)
!            y(i,j,inox)=y(i,j,ino)+y(i,j,ino2)+y(i,j,ino3)+ &
!                 2.*y(i,j,in2o5)+y(i,j,ihno4)
!
!         end do
!
!      end do
!
!      if ( nerr ) print*,'NOYmass: N-mass balance error, ncorr>5% ',&
!           'Maximum correction:',  ncormax,& 
!           'Average correction in this loop:',  ncorav/imax, imax, '(imax)'
!
!    end subroutine NOYmass
!
!
!
!    subroutine incc2c3
!      !
!      use budget_global,only : buddep_dat, sum_deposition, sum_chemistry
!      use budget_fine,only: depdry
!      implicit none
!      integer :: i1,n1,n2,jl,i,j,offsetl
!      ! nrj and nrr used for reaction budget calculations
!      integer,dimension(njnum),parameter    ::  nrj=(/io3,ino2,ih2o2,ihno3,ihno4,in2o5,ich2o,ich2o, &
!           ich3o2h,ino3,ino3,ipan,iorgntr,iald2,imgly,irooh/)
!      integer,dimension(nreac,2),parameter  ::  nrr = reshape((/ &
!           ino,iho2,ich3o2,ino2,ioh,ino2,ino,ino2,in2o5,ihno4,&
!           ino2,ihno4,iair,ih2o,io3,ico,io3,ih2o2,ich2o,ich4, &
!           ioh,ioh,ich3o2, ich3o2, iho2,iho2,in2o5,in2o5,ioh,ich2o,&
!           iald2,iald2,ic2o3,ic2o3,ipan,ic2o3,ic2o3,ipar,iror,iror,&
!           ioh,io3,ino3,ioh,io3,ioh,ioh,io3,ino3,ixo2,&
!           ixo2,ixo2n,ixo2,irxpar,iorgntr,ixo2n,idms,idms,idms,iso2,&
!           inh3,inh3,inh2,inh2,inh2,inh2,inh2,irn222, &
!           !second reaction partner (if monmolecular = 0)
!      io3,ino,ino,ioh,ihno3,io3,ino3,ino3, 0, ioh, &
!           iho2,0,0,0,iho2,ioh,ioh,ioh,ioh,ioh, &
!           ich3o2h,irooh,iho2, ich3o2, ioh,iho2,0,0,0,ino3,&
!           ioh,ino3,ino,ino2,0,ic2o3,iho2,ioh, 0, 0,&
!           iole,iole,iole,ieth,ieth,imgly,iisop,iisop,iisop,ino,&
!           ixo2,ino,iho2,ipar,ioh,iho2,ioh,ioh,ino3,ioh,&
!           iacid,ioh,ino,ino2,iho2,0,io3,0/),(/nreac,2/))
!
!      real :: c1,xdep
!
!      c1=dt*1000./xmair   !conversion to moles...
!
!      offsetl=0
!#ifdef MPI
!      if(myid>0) offsetl=sum(lmar(0:myid-1))
!#endif
!
!      ! reaction budgets
!
!      do i1=1,njnum !photolysis rates
!         n1=nrj(i1)
!         do j=jsr(region),jer(region)
!            do i=isr(region),ier(region)
!               if(zoomed(i,j)/=region) cycle
!               if(n1 > 0) cr2(i,j,i1)=cr2(i,j,i1)+rj(i,j,i1)*y(i,j,n1)
!            end do
!         end do
!      end do!i1=1,njnum
!      !
!      do i1=1,nreac !reactions
!         n1=nrr(i1,1) !make sure n1 > 0
!         n2=nrr(i1,2)
!         if (n2 > 0.) then
!            do j=jsr(region),jer(region)
!               do i=isr(region),ier(region)
!                  if(zoomed(i,j)/=region) cycle
!                  cr3(i,j,i1)= cr3(i,j,i1)+y(i,j,n1)*y(i,j,n2)*rr(i,j,i1)
!               end do
!            end do
!         else
!            do j=jsr(region),jer(region)
!               do i=isr(region),ier(region)
!                  if(zoomed(i,j)/=region) cycle
!                  cr3(i,j,i1)= cr3(i,j,i1)+y(i,j,n1)*rr(i,j,i1)
!               end do
!            end do
!         end if
!      end do  !i1=1,nreac
!
!      if ( offsetl + level == 1 ) then   ! deposition budget
!
!         do i1=1,ntrace
!            do j=jsr(region),jer(region)
!               do i=isr(region),ier(region)
!                  if(zoomed(i,j)/=region) cycle
!                  xdep = y(i,j,i1)*vd(region,i1)%surf(i,j)/ye(i,j,idz)* &
!                       c1*ye(i,j,iairm)/y(i,j,iair)  !from updated concentrations
!                  buddep_dat(region)%dry(i,j,i1) = &
!                       buddep_dat(region)%dry(i,j,i1) + xdep
!                  if ( region == nregions ) then
!                     depdry(i,j,i1) = depdry(i,j,i1) + xdep   !in mole
!                  end if
!                  if ( i1 == 1 ) then   !seperate deposition from 'other' chemistry
!                     sum_deposition(region) = sum_deposition(region) - &
!                          xdep*ra(1)*1e-3   ! in kg
!                     sum_chemistry(region) = sum_chemistry(region) +  &
!                          (y(i,j,1)-y0(i,j,1))/y(i,j,iair)* &
!                          ye(i,j,iairm)/xmair*ra(1) + xdep*ra(1)*1e-3  
!                  end if
!               end do
!            end do
!         end do   !i1
!      else   ! other layers
!         do j=jsr(region),jer(region)
!            do i=isr(region),ier(region)
!               if(zoomed(i,j)/=region) cycle
!               sum_chemistry(region) = sum_chemistry(region) + &
!                    (y(i,j,1)-y0(i,j,1))/y(i,j,iair)*ye(i,j,iairm)/xmair*ra(1)
!            end do
!         end do
!      end if   !level ==1
!    end subroutine incc2c3
!
!
!
!    subroutine reacbud
!      !------------------------------------------------------------------------
!      !
!      ! REACBUD increase reaction budgets for each reaction
!      !         arrays nrr and nrj determine which species are
!      !         involved in a reaction
!      !
!      !
!      !------------------------------------------------------------------------
!      use budget_global,only : budrjg,budrrg,budrwg,budg_dat,nzon_vg
!      use budget_fine,only: budrj,budrr,budrw,nzon,nzon_v
!      implicit none
!      integer :: i1,i,j,nzone,nzone_v
!      real    :: c1
!      !
!      ! attribute to regions
!      !
!      c1=dt*1000./xmair   !conversion to moles...
!      do j=jsr(region),jer(region)
!         do i=isr(region),ier(region)
!            if(zoomed(i,j)/=region) cycle
!            if(region==nregions) then    !finest region budget....
!               nzone=nzon(i,j)
!               nzone_v=nzon_v(offsetl+level)    !level is passed to ebi...
!               do i1=1,njnum
!                  budrj(nzone,nzone_v,i1)=budrj(nzone,nzone_v,i1)+ &
!                       cr2(i,j,i1)*c1*ye(i,j,iairm)/y(i,j,iair)   !units mole
!               end do !njnum
!               do i1=1,nreac
!                  budrr(nzone,nzone_v,i1)=budrr(nzone,nzone_v,i1)+ &
!                       cr3(i,j,i1)*c1*ye(i,j,iairm)/y(i,j,iair)   !units mole
!               end do
!               do i1=1,nreacw
!                  budrw(nzone,nzone_v,i1)=budrw(nzone,nzone_v,i1)- &
!                       cr4(i,j,i1)*(1000./xmair)*ye(i,j,iairm)/y(i,j,iair)
!                  !note: changed sign to get 'positive' budget, just a 
!                  !      matter of definition, !CMK
!               end do
!            end if
!            nzone=budg_dat(region)%nzong(i,j)    !global budget
!            nzone_v=nzon_vg(offsetl+level)    !level is passed to ebi...
!            do i1=1,njnum
!               budrjg(nzone,nzone_v,i1)=budrjg(nzone,nzone_v,i1)+ &
!                    cr2(i,j,i1)*c1*ye(i,j,iairm)/y(i,j,iair)   !units mole
!            end do !njnum
!            do i1=1,nreac
!               budrrg(nzone,nzone_v,i1)=budrrg(nzone,nzone_v,i1)+ &
!                    cr3(i,j,i1)*c1*ye(i,j,iairm)/y(i,j,iair)   !units mole
!            end do
!            do i1=1,nreacw
!               budrwg(nzone,nzone_v,i1)=budrwg(nzone,nzone_v,i1)- &
!                    cr4(i,j,i1)*(1000./xmair)*ye(i,j,iairm)/y(i,j,iair)   ! mole
!               !note: changed sign to get 'positive' budget, just a
!               !      matter of definition, !CMK
!            end do
!         end do
!      end do
!
!    end subroutine REACBUD
!
!
!  end subroutine ebi
!
!
!
!  subroutine wetS(region,level,y0,dt,y,ye,c4)
!    !**********************************************************************
!    !     
!    !wetS - aqueous phase chemistry of sulfur  (and other)
!    !programmed by Ad Jeuken (KNMI) and Frank Dentener (IMAU)
!    !adapted for TM5 by Maarten Krol (IMAU) 1-2002
!    !
!    !purpose
!    !-------
!    !oxidation of SO2 and uptake of other gases in the aqueous phase
!    !
!    !interface
!    !---------
!    !call wetS(region,level,zoomed,y0,dt,y,ye,c4)
!    !region  region under operation (provides im,jm,lm via chemistry.mod)
!    !level   vertical level 
!    !zoomed       which cells to skip (chemistrty is done by child)
!    !y0    initial concentration
!    !dt    chemistry timestep
!    !yt    concentrations at time is t
!    !ye    extra fields (temperature, clouds, ....)
!    !c4    budget accumulatior
!    !   
!    !method
!    !------
!    !implicit solution of oxidation of SO2
!    !
!    !external
!    !--------
!    !none
!    !
!    !reference
!    !---------
!    !-
!    !**********************************************************************
!
!    use global_data,  only: region_dat
!    use reaction_data,only: nreacw,ntlow,kso2hp,kso2o3
!    use chem_param
!    use budget_global, only: sum_wet
!    use dims, only: isr, jsr, ier,jer, im, jm
!    use Binas, only: Avog
!
!    !    use toolbox, only: dumpfield, escape_tm
!    !    use dims, only: lm
!    implicit none
!
!    ! input/output
!
!    integer,intent(in)                                        :: region
!    integer,intent(in)                                        :: level
!    real,intent(in),dimension(im(region),jm(region),maxtrace) :: y0
!    real,intent(in)                                           :: dt 
!    real,intent(out),dimension(im(region),jm(region),maxtrace):: y 
!    real,dimension(im(region),jm(region),n_extra)             :: ye   !extra fields (temp, cc, pH)
!    real,dimension(im(region),jm(region),nreacw),intent(inout):: c4
!
!    ! local
!
!    integer,dimension(:,:),pointer                            :: zoomed
!    integer n,i,j,l,itemp,iter
!    real x1,x2,x3,b1,b2,so2x,dh2o2,dso2,disc,dnh3,dn2o5,xso2o3a,xso2o3b
!    real,parameter :: co2=3.20e-4,rg=0.08314
!    real,dimension(:,:),allocatable :: hkso2  ! Henry's constant for sulfur dioxide
!    real,dimension(:,:),allocatable :: hkh2o2 ! Henry's constant for hydroperoxide
!    real,dimension(:,:),allocatable :: hko3       ! Henry's constant for ozone
!    real,dimension(:,:),allocatable :: dkso2      ! Dissociation constant for SO2
!    real,dimension(:,:),allocatable :: dkhso3     ! Dissociation constant for HSO3-
!    real,dimension(:,:),allocatable :: dkh2o      ! dissociation constant water
!    real,dimension(:,:),allocatable :: dknh3      ! dissociation constant ammonia
!    real,dimension(:,:),allocatable :: hknh3      ! Henry's constant ammonia
!    real,dimension(:,:),allocatable :: hkco2      ! Henry's constant CO2
!    real,dimension(:,:),allocatable :: dkco2      ! Dissociation constant CO2
!    real phs4                       ! effective dissolvation of S(IV)
!    real phso2                      ! effective dissolvation of SO2
!    real phh2o2                     ! effective dissolvation of H2O2
!    real phozone                    ! effective dissolvation of O3
!    real,dimension(:,:),allocatable :: hplus      !concen§tration h+
!    real a1,a2,a,b,c,z              ! help variables
!    real xcov,xliq,xl,temp,rt,ztr,h2o,air,press ! meteo
!    real,dimension(:,:,:),allocatable   :: rw ! reaction rates
!    logical,dimension(:,:),allocatable ::  cloudy
!    !    character(len=2) :: levelc
!
!    ! start
!    
!    zoomed => region_dat(region)%zoomed
!
!    allocate(hkso2      (im(region),jm(region)))
!    allocate(hkh2o2     (im(region),jm(region)))
!    allocate(hko3       (im(region),jm(region)))
!    allocate(dkso2      (im(region),jm(region)))
!    allocate(dkhso3     (im(region),jm(region)))
!    allocate(dkh2o      (im(region),jm(region)))
!    allocate(dknh3      (im(region),jm(region)))
!    allocate(hknh3      (im(region),jm(region)))
!    allocate(hkco2      (im(region),jm(region)))
!    allocate(dkco2      (im(region),jm(region)))
!    allocate(hplus      (im(region),jm(region)))
!    allocate(rw (im(region),jm(region),nreacw))
!    allocate(cloudy (im(region),jm(region)))
!
!    !-----------------------------
!    ! wet phase reactions
!    !-----------------------------
!    rw   =0.0 
!    hplus=0.0
!
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          cloudy(i,j)=.false.
!
!          ! lwc is dimensionless 
!          if ((ye(i,j,ilwc) > 1e-10).and.(ye(i,j,icc) > 0.01)) then 
!             cloudy(i,j)=.true. 
!             TEMP=ye(i,j,i_temp)
!             ZTR=(1./TEMP-1./298)
!             RT=TEMP*rg
!             ITEMP=nint(TEMP-float(ntlow))
!             !
!             ! Henry and dissociation equilibria
!             !
!             dkh2o(i,j) =1.01e-14*exp(-6706.0 *ztr)   !h2o<=>hplus+so3--
!             hkco2(i,j) =3.4e-2*(2420.*ztr)           ! is already dimensionless
!             dkco2(i,j) =4.5E-7*exp(-1000.*ztr)       !co2aq<=>hco3- + hplus
!             hkso2(i,j) =henry(iso2,itemp)*rt         !dimensionless
!             dknh3(i,j) =1.8e-5*exp(-450.*ztr)        !nh3<=>nh4+ + OH-
!             hknh3(i,j) =henry(inh3,itemp)*rt         !dimensionless
!             hkh2o2(i,j)=henry(ih2o2,itemp)*rt        !dimensionless
!             hko3(i,j)  =henry(io3,itemp)*rt          !dimensionless
!             dkso2(i,j) =1.7e-2*exp(2090.*ztr)        !so2<=>hso3m+hplus
!             dkhso3(i,j)=6.6e-8*exp(1510.*ztr)        !hso3m<=>so3-- + hplus
!             !
!             ! calculate H+ from initial sulfate, ammonium, hno3, and nh3
!             ! if solution is strongly acidic no further calculations are performed
!             !
!
!             xl=ye(i,j,ilwc)*Avog*1e-3/ye(i,j,icc)
!             !x1 is initial strong acidity from SO4 and NO3
!             !
!             !acidity from strong acids alone
!             !  
!             hplus(i,j)=(2.*y0(i,j,iso4)+y0(i,j,imsa)-y0(i,j,inh4)+ &
!                  y0(i,j,ihno3)+y0(i,j,ino3_a))/xl
!          end if
!       end do
!    end do
!    do iter=1,10
!       do j=jsr(region),jer(region)
!          do i=isr(region),ier(region)
!             if ( zoomed(i,j) /= region ) cycle
!             ! only if solution pH>4.5
!             if ( cloudy(i,j) .and. hplus(i,j) < 3e-5 ) then
!                xl=ye(i,j,ilwc)*Avog*1e-3/ye(i,j,icc)
!                x1=(2.*y0(i,j,iso4)+y0(i,j,imsa)+y0(i,j,ihno3)+ &
!                     y0(i,j,ino3_a))/xl    
!                !x2 is initial total NHx
!                x2=(y0(i,j,inh3)+y0(i,j,inh4))/xl
!                !x3 is combined dissolution and solubility const for CO2
!                x3=dkco2(i,j)*hkco2(i,j)*co2 
!                a1=dkh2o(i,j)/dknh3(i,j)*(1.+1./hknh3(i,j)) ! integration constant
!                a2=y0(i,j,iso2)/xl      !initial SO2
!                z=a2/(hplus(i,j)/dkso2(i,j)*(1.+1./hkso2(i,j))+ &
!                     dkhso3(i,j)/hplus(i,j)+1.)
!                a=1.+x2/(a1+hplus(i,j))
!                b=-x1-z
!                c=-x3-2.*dkhso3(i,j)*z
!                z=max(0.,(b*b-4.*a*c))
!                hplus(i,j)=max(1.e-10,(-b+sqrt(z))/(2.*a))
!             end if
!          end do !
!       end do ! i,j loop§
!    end do   !iter
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          if (cloudy(i,j)) then
!             temp=ye(i,j,i_temp)
!             ZTR=(1./TEMP-1./298)
!             xliq=ye(i,j,ilwc)/ye(i,j,icc)
!             xl=ye(i,j,ilwc)*Avog*1e-3/ye(i,j,icc)
!             ye(i,j,iph)=-log10(hplus(i,j))     ! pH for diagnostics 
!
!             ! phase factor ratio of aqueous phase to gas phase concentration
!
!             phs4   =hkso2(i,j) *(1.+dkso2(i,j)/hplus(i,j)+ &
!                  dkhso3(i,j)/hplus(i,j)/hplus(i,j))*xliq
!             phso2  =hkso2(i,j) *xliq
!             phh2o2 =hkh2o2(i,j)*xliq
!             phozone=hko3(i,j)  *xliq
!
!             ! the original rate equations could be partly in look-up table
!
!             rw(i,j,KSO2HP) =8e4*exp(-3560.*ztr)/(0.1+hplus(i,j))
!             XSO2O3A=4.39e11*exp(-4131./temp)+2.56e3*exp(-966./temp)  !S(IV)
!             XSO2O3B=2.56e3*exp(-966./temp)/hplus(i,j)  
!             !divide by [H+]!S(IV)
!
!             !  make rate constants dimensionless by multiplying 
!             !  by (1./xliq/avo=6e20)  
!             !  multiply with the fractions of concentrations residing 
!             !  in the aqueous phase
!
!             rw(i,j,KSO2HP)=rw(i,j,KSO2HP)/xl*phso2/(1.+phs4)*phh2o2/(1.+phh2o2)
!             rw(i,j,KSO2O3)=(XSO2O3A+XSO2O3B)/xl*phs4/(1.+phs4)*phozone/ &
!                  (1.+phozone)
!          end if !cloudy
!       end do ! 
!    end do ! I,J, LOOP
!!    write(levelc,'(i2.2)') level
!!    if(level == 1) then 
!!       call dumpfield(0,'rw.hdf',im(region),jm(region),nreacw,1,rw,'rw'//levelc)
!!       call dumpfield(1,'rw.hdf',im(region),jm(region),n_extra,1,ye,'ye'//levelc)
!!    else
!!       call dumpfield(1,'rw.hdf',im(region),jm(region),nreacw,1,rw,'rw'//levelc)
!!       call dumpfield(1,'rw.hdf',im(region),jm(region),n_extra,1,ye,'ye'//levelc)
!!    end if
!!    if(level == lm(1) ) call escape_tm(' forced stop ')
!
!    ! Start main loop
!    do j=jsr(region),jer(region)
!       do i=isr(region),ier(region)
!          if(zoomed(i,j)/=region) cycle
!          !
!          ! only cloud chemistry if substantial amount of clouds are present 
!          !
!          if (cloudy(i,j)) then
!             !
!             ! oxidation of S(IV) by O3
!             !
!             so2x=y0(i,j,iso2)
!             xcov=ye(i,j,icc)
!             x1=min(100.,rw(i,j,kso2o3)*y0(i,j,io3)*dt)
!             dso2=y0(i,j,iso2)*xcov*(exp(-x1)-1.)
!             !only applied to xcov part of cloud
!    !CMK         print *, i,j, xcov, x1, y0(i,j,iso2), dso2
!             dso2=max(-y0(i,j,io3)*xcov,dso2)! limit to o3 availability
!             y(i,j,iso2)=y0(i,j,iso2)+dso2 
!             !NOTE CMK: paralel MPI should take care here!
!             y(i,j,iso4)=y0(i,j,iso4)-dso2
!             y(i,j,io3)=y0(i,j,io3)+dso2
!             if ( io3 == 1 ) sum_wet(region) = sum_wet(region)- &
!                  dso2 *ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!             if ( iso2 == 1 ) sum_wet(region) = sum_wet(region)+ &
!                  dso2   *ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!             if ( iso4 == 1 ) sum_wet(region) = sum_wet(region)- &
!                  dso2  *ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!             c4(i,j,1)=c4(i,j,1)+dso2
!             xliq=ye(i,j,ilwc)/ye(i,j,icc)
!             !
!             ! oxidation of S(IV) by H2O2
!             !
!             !*** here we explicitly solve the dv: 
!             !    y'=P-Q*y-R*y*y (P and Q are 0=>b3=0.)
!             !
!             so2x=y(i,j,iso2)
!             if ( so2x > tiny(so2x) ) then
!               b1=rw(i,j,kso2hp)
!               b2=b1*(y0(i,j,ih2o2)-so2x)
!               disc=min(100.,sqrt(b2*b2))                ! disc is b2 for b3=0.0
!               x1=(b2-disc)/(-2.*b1)                     ! in this case x1 =0.
!               x2=(b2+disc)/(-2.*b1)
!               x3=(so2x-x1)/(so2x-x2)*exp(-disc*dt)
!               so2x=(x1-x2*x3)/(1.-x3)
!               dso2=(so2x-y(i,j,iso2))*xcov
!               dso2=max(dso2,-y0(i,j,ih2o2)*xcov)
!               y(i,j,iso2) =y(i,j,iso2)+dso2       ! dso2 is loss of SO2 and H2O2
!               y(i,j,iso4)=y(i,j,iso4)-dso2
!               y(i,j,ih2o2) =y0(i,j,ih2o2)+dso2 
!               if ( ih2o2 == 1 ) sum_wet(region) = sum_wet(region)- &
!                    dso2  *ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!               if ( iso2 == 1 ) sum_wet(region) = sum_wet(region)- &
!                    dso2   *ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!               if ( iso4 == 1 ) sum_wet(region) = sum_wet(region)+ &
!                    dso2  *ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!               c4(i,j,2)=c4(i,j,2)+dso2
!             end if
!
!             !
!             ! NH3 uptake in cloud droplets is limited by H2SO4 availability
!             ! no HNO3 is considered at this point
!             ! assume instantaneous uptake of NH3 incloud  only in cloudy part
!             !
!             dnh3=max((2.*y(i,j,iso4)+y0(i,j,imsa)-y0(i,j,inh4))*xcov,0.)
!             dnh3=max(-y0(i,j,inh3)*xcov,-dnh3)
!             y(i,j,inh3)=y0(i,j,inh3)+dnh3                 ! dnh3 is loss of NH3  
!             y(i,j,inh4)=y0(i,j,inh4)-dnh3
!             if ( inh3 == 1 ) sum_wet(region) = sum_wet(region) - &
!                  dnh3*ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!             if ( inh4 == 1 ) sum_wet(region) = sum_wet(region) + &
!                  dnh3*ye(i,j,iairm)/ (fscale(1)*y(i,j,iair))
!             c4(i,j,3)=c4(i,j,3)+dnh3
!          end if     !cloudy
!       end do ! i,j,loop
!    end do ! 
!
!    !free memory
!    deallocate(hkso2      )
!    deallocate(hkh2o2     )
!    deallocate(hko3       )
!    deallocate(dkso2      )
!    deallocate(dkhso3     )
!    deallocate(dkh2o      )
!    deallocate(dknh3      )
!    deallocate(hknh3      )
!    deallocate(hkco2      )
!    deallocate(dkco2      )
!    deallocate(hplus      )
!    deallocate(rw )
!    deallocate(cloudy )
!
!    nullify(zoomed)
!
!  !  write(levelc, '(i2.2)') level
!  !  call dumpfield(1,'wetS.hdf', im(region), jm(region), maxtrace, 1, y0, 'y0'//levelc)
!  !  call dumpfield(1,'wetS.hdf', im(region), jm(region), n_extra , 1, ye, 'ye'//levelc)
!  !  call dumpfield(1,'wetS.hdf', ntracet, ntemp, 1 , 1, henry, 'henry'//levelc)
!  !  call dumpfield(1,'wetS.hdf', im(region), jm(region), maxtrace, 1, y, 'y'//levelc)
!
!  end subroutine wets
!
!
!
!  subroutine mark_trac(region,level,y,rr,rj,dt,ye)
!    !      ---------
!    !   call subroutine mark_trac(region,level,zoomed,y,rr,rj,dt,ye)
!    !      region,level               :: where and which cells
!    !      zoomed                     :: which cells are done by child?
!    !      y    :: concentrations in layer
!    !      rr   :: reaction rates
!    !      rj   :: photolysis rates
!    !      dt   :: time step
!    !      ye   :: help fields ( air masses )
!    !         
!    !      method
!    !      ------
!    !      calculate nox/pan/orgn/hno3 analogous to ebi scheme
!    !      ozone production from marked nox
!    !      simple nhx chemistry, scaled to real nhx
!    !
!    !      fjd Mon Aug 10 16:55:35 MET 1998/Fri Jan  1 17:08:37 MET 1999
!    !      mk   adapted for TM5 jan/2002
!    !-------------------------------------------------------------
!    use global_data,   only : region_dat
!    use budget_global, only : budmarkg,budg_dat,nzon_vg
!    use budget_fine,   only : budmark,nzon,nzon_v
!    use chem_param
!    use dims, only: isr, ier, jsr, jer, at, bt, im, jm
!
!    implicit none
!
!    ! input/output
!    integer, intent(in)                            :: region,level
!    real,dimension(im(region),jm(region),maxtrace) :: y
!    real,dimension(im(region),jm(region),nreac),intent(in):: rr
!    real,dimension(im(region),jm(region),njnum),intent(in):: rj
!    real                                           :: dt
!    real, dimension(im(region),jm(region),n_extra) :: ye
!
!    ! local
!    integer, dimension(:,:),pointer                :: zoomed
!
!    ! start
!
!    zoomed => region_dat(region)%zoomed
!
!    call mark_o3s
!    !
!    ! more marked tracers possible here
!    !
!    nullify(zoomed)
!
!  contains
!
!
!    subroutine mark_o3s
!      !---------------------------------------------------
!      ! marked tracer  O3S stratospheric ozone
!      !---------------------------------------------------
!
!#ifdef MPI
!      use mpi_const,only: myid,lmar
!#endif
!      use dry_deposition, only: vd
!      implicit none
!      integer :: i,j,nzone,nzone_v
!      real                           :: p3,xl3,o3old
!      integer :: offsetl
!
!      offsetl=0
!#ifdef MPI
!      if(myid>0) offsetl=sum(lmar(0:myid-1))
!#endif
!
!      do j=jsr(region),jer(region)
!         do i=isr(region),ier(region)
!            if(zoomed(i,j)/=region) cycle
!            if (at(offsetl+level+1)+bt(offsetl+level+1)*1e5<= 14000) then ! 
!               ! well, you want to count all layers below 140 hPa 
!               ! (given surface pressure of 1e5 Pa)
!               ! in the current model setup (25 layers) this means  
!               ! 12077 + 1e5*0.00181 = 12258 Pa and above...
!
!               ! p3: production of o3 in stratosphere
!               P3 = rj(i,j,jano3)*y(i,j,ino3)+ &
!                    rj(i,j,jno2)*y(i,j,ino2)
!               XL3= rr(i,j,ko3ho2)*y(i,j,iho2)+&
!                    rr(i,j,ko3oh)*y(i,j,ioh)+ &
!                    rr(i,j,kno2o3)*y(i,j,ino2)+&
!                    rj(i,j,jo3d)+&
!                    rr(i,j,knoo3)*y(i,j,ino)+&
!                    rr(i,j,kc62)*y(i,j,ieth)+&
!                    rr(i,j,kc58)*y(i,j,iole)+&
!                    rr(i,j,kc77)*y(i,j,iisop)
!            else
!               !
!               ! these are only the net destruction reactions
!               !
!               P3 = 0.
!               XL3= rr(i,j,ko3ho2)*y(i,j,iho2)+&
!                    rr(i,j,ko3oh)*y(i,j,ioh)+&
!                    rj(i,j,jo3d)+&
!                    rr(i,j,kc62)*y(i,j,ieth)+&
!                    rr(i,j,kc58)*y(i,j,iole)+&
!                    rr(i,j,kc77)*y(i,j,iisop)
!               ! add up deposition....
!               if ( offsetl + level == 1 ) &
!                    XL3 = XL3 + vd(region,io3)%surf(i,j)/ye(i,j,idz)
!            end if
!            o3old=y(i,j,io3s)
!            y(i,j,io3s)=(o3old+p3*dt)/(1.+xl3*dt)
!            if ( region == nregions ) then
!               nzone=nzon(i,j)
!               nzone_v=nzon_v(level+offsetl) 
!               ! budget in mole
!               budmark(nzone,nzone_v,1)=budmark(nzone,nzone_v,1)+ &
!                    (y(i,j,io3s)-o3old)*ye(i,j,iairm)*1000./xmair/y(i,j,iair)
!            end if
!            nzone=budg_dat(region)%nzong(i,j)    ! global budget
!            nzone_v=nzon_vg(level+offsetl) 
!            budmarkg(nzone,nzone_v,1)=budmarkg(nzone,nzone_v,1)+ &
!                 (y(i,j,io3s)-o3old)*ye(i,j,iairm)*1000./xmair/y(i,j,iair)
!         end do
!
!      end do !i,j, l
!
!    end subroutine mark_o3s
!
!  end subroutine mark_trac
!
!
!
end module ebischeme

