!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################

module reaction_data

  implicit none

  ! --- in/out ---------------------------------------

  public


  ! --- const ----------------------------------------

  integer, parameter :: nreac=68       !number of thermal reactions
  integer, parameter :: nreacw=4       !number of aqueous phase reactions
!  integer, parameter :: nthigh=320
!  integer, parameter :: ntlow=165
!  integer, parameter :: ntemp=nthigh-ntlow
!  integer, parameter :: nrat=nreac+19
  character(len=8), dimension(nreac),parameter :: ratnam=(/&
       'RNOO3   ','RHO2NO  ','RMO2NO  ','RNO2OH  ','ROHHNO3 ',&
       'RNO2O3  ','RNONO3  ','RNO2NO3 ','RN2O5   ','RHNO4OH ',&
       'RNO2HO2 ','RHNO4M  ','RODM    ','RH2OOD  ','RO3HO2  ',&
       'RCOOH   ','RO3OH   ','RHPOH   ','RFRMOH  ','RCH4OH  ',&
       'ROHMPER ','ROHROOH ','RMO2HO2 ','RMO2MO2 ',&
       'RHO2OH  ','RHO2HO2 ','RN2O5AQ ','RN2O5L  ','RH2OH   ',&
       'RC41    ','RC43    ','RC44    ','RC46    ','RC47    ','RC48    ',&
       'RC49    ','RC50    ','RC52    ','RC53    ','RC54    ','RC57    ',&
       'RC58    ','RC59    ','RC61    ','RC62    ','RC73    ','RC76    ',&
       'RC77    ','RC78    ','RC79    ','RC80    ','RC81    ','RC82    ',&
       'RC83    ','RC84    ','RC85    ',&
       'RCdmsoha','RCdmsohb','RCdmsno3','RCso2oh ','RCnh3so4',&
       'RCnh3oh ','RCnh2no ','RCnh2no2','RCnh2ho2','RCnh2o2 ',&
       'RCnh2o3 ','DRN222  '/)
  character(len=8), dimension(nreacw),parameter :: rwnam = &
       (/'kso2o3  ','kso2hp  ','knh3so4 ','void    '/)
  !
  ! reaction rates
  !
!  integer, parameter :: KNOO3=1
!  integer, parameter :: KHO2NO=2
!  integer, parameter :: KMO2NO=3
!  integer, parameter :: KNO2OH=4
!  integer, parameter :: KOHHNO3=5
!  integer, parameter :: KNO2O3=6
!  integer, parameter :: KNONO3=7
!  integer, parameter :: KNO2NO3=8
!  integer, parameter :: KN2O5=9
!  integer, parameter :: KHNO4OH=10
!  integer, parameter :: KNO2HO2=11
!  integer, parameter :: KHNO4M=12
!  integer, parameter :: KODM=13
!  integer, parameter :: KH2OOD=14
!  integer, parameter :: KO3HO2=15
!  integer, parameter :: KCOOH=16
!  integer, parameter :: KO3OH=17
!  integer, parameter :: KHPOH=18
!  integer, parameter :: KFRMOH=19
!  integer, parameter :: KCH4OH=20
!  integer, parameter :: KOHMPER=21
!  integer, parameter :: KOHROOH=22
!  integer, parameter :: KMO2HO2=23
!  integer, parameter :: KMO2MO2=24
!  integer, parameter :: KHO2OH=25
!  integer, parameter :: KHO2HO2=26
!  integer, parameter :: KN2O5AQ=27
!  integer, parameter :: KN2O5L=28
!  integer, parameter :: KH2OH=29
!  integer, parameter :: KC41=30
!  integer, parameter :: KC43=31
!  integer, parameter :: KC44=32
!  integer, parameter :: KC46=33
!  integer, parameter :: KC47=34
!  integer, parameter :: KC48=35
!  integer, parameter :: KC49=36
!  integer, parameter :: KC50=37
!  integer, parameter :: KC52=38
!  integer, parameter :: KC53=39
!  integer, parameter :: KC54=40
!  integer, parameter :: KC57=41
!  integer, parameter :: KC58=42
!  integer, parameter :: KC59=43
!  integer, parameter :: KC61=44
!  integer, parameter :: KC62=45
!  integer, parameter :: KC73=46
!  integer, parameter :: KC76=47
!  integer, parameter :: KC77=48
!  integer, parameter :: KC78=49
!  integer, parameter :: KC79=50
!  integer, parameter :: KC80=51
!  integer, parameter :: KC81=52
!  integer, parameter :: KC82=53
!  integer, parameter :: KC83=54
!  integer, parameter :: KC84=55
!  integer, parameter :: KC85=56
!  integer, parameter :: Kdmsoha=57
!  integer, parameter :: Kdmsohb=58
!  integer, parameter :: Kdmsno3=59
!  integer, parameter :: Kso2oh=60
!  integer, parameter :: Knh3so4=61
!  integer, parameter :: Knh3oh=62
!  integer, parameter :: Knh2no=63
!  integer, parameter :: Knh2no2=64
!  integer, parameter :: Knh2ho2=65
!  integer, parameter :: Knh2o2=66
!  integer, parameter :: Knh2o3=67
!  integer, parameter :: krn222=68
!  !
!  ! additional parameters needed for 3-body and other reactions
!  !
!  integer, parameter :: KNO2OHA=nreac+1
!  integer, parameter :: KNO2OHB=nreac+2
!  integer, parameter :: KOHHNO3A=nreac+3
!  integer, parameter :: KOHHNO3B=nreac+4
!  integer, parameter :: KOHHNO3C=nreac+5
!  integer, parameter :: KNO2NO3A=nreac+6
!  integer, parameter :: KNO2NO3B=nreac+7
!  integer, parameter :: KNO2HO2A=nreac+8
!  integer, parameter :: KNO2HO2B=nreac+9
!  integer, parameter :: KHO2HO2A=nreac+10
!  integer, parameter :: KHO2HO2B=nreac+11
!  integer, parameter :: KHO2HO2C=nreac+12
!  integer, parameter :: KC47A=nreac+13
!  integer, parameter :: KC47B=nreac+14
!  integer, parameter :: KC61A=nreac+15
!  integer, parameter :: KC61B=nreac+16
!  integer, parameter :: kso2oha=nreac+17
!  integer, parameter :: kso2ohb =nreac+18
!  integer, parameter :: kdmsohc =nreac+19
!
!  integer, parameter :: kso2hp=1
!  integer, parameter :: kso2o3=2
!
!  ! rates_lut : for reaction rate look up table !WP!
!  real,dimension(nrat,ntemp) :: rates_lut
!
!  ! for N2O5 removal on aerosol parameterisation:
!  ! nr_interval : number of interval for intergration over aerosol size
!  integer,parameter          :: nr_interval=13
!  ! np_interval : number of interval for integration over pressure
!  integer,parameter          :: np_interval=11
!  integer,parameter          :: ncomponent=2
!  ! aerdens : density of aerosols(g/cm3) (water 1.0 )
!  real,parameter             :: aerdens=1.8
!  ! rra     : lookup table radius..
!  real,dimension(nr_interval):: rra = &
!       (/1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.,4.,5.,10.,15./)
!  ! hetrem  : heterogeneous removal coefficient [s-1]
!  !           for 1 ppbv NH4HSO4 at given aerosol radius and pressure
!  real,dimension(nr_interval,np_interval,ncomponent) :: hetrem
!
!
end module reaction_data
