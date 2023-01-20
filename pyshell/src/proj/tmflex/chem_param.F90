!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################

module chem_param

  use binas, only : xmair, Avog, grav, amu

  implicit none

  ! define specific data structures

  type emis_data
     ! emis_data%surf  : any type of surface data
     real,dimension(:,:),pointer     :: surf
  end type emis_data

  type d3_data
     ! d3_data%d3   : 3D data, e.g. nox emissions
     real,dimension(:,:,:), allocatable :: d3
  end type d3_data

  type d4_data
    type(d3_data), allocatable  :: field_t(:) ! time dimension added to d3_data
  end type d4_data

  type d23_data
     ! d23_data%d23  : lat/pres fields, e.g. o3 climatology in ppmv
     real,dimension(:,:),pointer     :: d23
  end type d23_data

  type d2_data
     ! d2_data%d2  : lat fields, e.g. hno3/o3 ratios at 10 hPa
     real,dimension(:),pointer       :: d2
  end type d2_data

  type isop_data
    ! isop_data%scalef_isop  : (jm,ntim) scalefactor isoprene emissions
     real,dimension(:,:),pointer   :: scalef_isop
  end type isop_data

  type chem_data
     ! chem_data%rm_k  : 'chemistry' tracers are parallel over layers...
     real,dimension(:,:,:,:),pointer :: rm_k
  end type chem_data

  !
  ! definition of the chemistry: #reactions, order of species, etc.
  ! parameters needed for chemistry and rate constants
  !

  ! nmark: number of 'marked' tracers
  integer, parameter :: nmark = 2
  character(len=8), dimension(nmark) :: marknam = (/'CO2     ','CO2fg   '/)

  ! ntrace: number of tracers for chemistry
  integer, parameter :: ntrace  = 2
  ! ntracet: number of transported tracers
  integer, parameter :: ntracet = 2    !
  ! ntrace_chem: non-transported tracers
  integer, parameter :: ntrace_chem  = ntrace-ntracet
  ! maxtrace: total number of tracers
  integer, parameter :: maxtrace = ntrace + 0
  ! max length of a tracer name
  integer, parameter :: tracer_name_len = 8
  !
  ! components numbers
  !
  integer, parameter :: ico2=1
  !
  ! additional fields used in chemistry routine alone
  ! (more meteo-like files in units different from #/cm3)
  !
  integer,parameter                   :: n_extra = 0
  integer,parameter                   :: nstd    = 1
  integer, dimension(nstd), parameter :: istd    = (/ico2/)
  !
  ! species name
  !
  character*8, dimension(maxtrace) :: names = (/'CO2     ','CO2fg   '/)
  !
  ! some stuff for the budgets....
  !
  integer, parameter           :: njnum = 1
  integer, parameter           :: nrdep = 0
  integer, parameter           :: nreac = 1
  integer, parameter           :: nreacw = 1
  character(len=tracer_name_len), dimension(njnum)  :: jnam   = (/'dummy   '/)
  character(len=tracer_name_len), dimension(nreacw) :: rwnam  = (/'dummy   '/)
  character(len=tracer_name_len), dimension(nreac)  :: ratnam = (/'dummy   '/)

  ! molar weights of components

  real,    parameter :: xmh=1.0079
  real,    parameter :: xmn=14.0067
  real,    parameter :: xmc=12.01115
  real,    parameter :: xms=32.064
  real,    parameter :: xmo=15.9994
  real,    parameter :: xmf=18.9984
  real,    parameter :: xmcl=35.453
  real,    parameter :: xmsf6 = xms + 6.0 * xmf
  real,    parameter :: xmo3=xmo*3
  real,    parameter :: xmnox=xmn
  real,    parameter :: xmh2o2=xmo*2.+xmh*2.
  real,    parameter :: xmch4=xmc+xmh*4.
  real,    parameter :: xmco=xmc+xmo
  real,    parameter :: xmco2 = xmc + 2.0*xmo
  real,    parameter :: xmhno3=xmh+xmn+xmo*3.
  real,    parameter :: xmmepe=xmc+xmh*4.+xmo*2.
  real,    parameter :: xmch2o=xmc+xmh*2.+xmo
  real,    parameter :: xmno=xmn+xmo
  real,    parameter :: xmho2=xmh+xmo*2.
  real,    parameter :: xmch3o2=xmc+2.*xmo+3.*xmh
  real,    parameter :: xmoh=xmo+xmh
  real,    parameter :: xmno2=xmn+2.*xmo
  real,    parameter :: xmno3=xmn+3.*xmo
  real,    parameter :: xmn2o5=2.*xmn+5.*xmo
  real,    parameter :: xmhno4=xmn+4.*xmo+xmh
  real,    parameter :: xmn2o=2.*xmn + xmo
  !FD real,    parameter :: xmair=28.94
  real,    parameter :: xmh2o=xmh*2+xmo
  ! xmpar is the results of the CBM4 implementation...calculate as C
  real,    parameter :: xmpar=xmc
  real,    parameter :: xmeth=2.*xmc
  real,    parameter :: xmole=2.*xmc
  real,    parameter :: xmisop=5.*xmc+8.*xmh
  real,    parameter :: xmgly=3.*xmc
  real,    parameter :: xmald2=2.*xmc
  real,    parameter :: xmc2o3=2.*xmc+3.*xmo+3.*xmh
  ! xmpan: ch3co-o2-no2
  real,    parameter :: xmpan=2.*xmc+3.*xmh+3.*xmo+xmn+2.*xmo
  real,    parameter :: xmror=2.*xmc+4.*xmh
  real,    parameter :: xmrxpar=xmc
  real,    parameter :: xmrooh=xmc+3.*xmh+2*xmo
  real,    parameter :: xmorgntr=xmn
  real,    parameter :: xmxo2=2.*xmo+xmc
  real,    parameter :: xmxo2n=2.*xmo+xmc
  ! attention xmso2: conversion emissions done when added...
  real,    parameter :: xmso2=xms+2.*xmo
  real,    parameter :: xmdms=xms+2*xmc+6*xmh
  ! attention xmnh3: conversion emissions when added...
  real,    parameter :: xmnh3=xmn+3.*xmh
  ! attention xmnh4: conversion emissions when added...
  real,    parameter :: xmnh4=xmn+4.*xmh
  real,    parameter :: xmmsa=xms+xmc+3*xmo+4*xmh
  real,    parameter :: xmnh2=xmn+xmh*2.
  real,    parameter :: xmso4=xms+4.*xmo
  real,    parameter :: xmno3_a=xmn+xmo*3
  real,    parameter :: xmrn222=222.
  real,    parameter :: xmpb210=210.
  ! xmnmv: this is a dummy molecular mass
  real,    parameter :: xmnmv=999.


  real,dimension(maxtrace),parameter :: ra = (/ xmco2, xmco2 /)
  real,dimension(maxtrace),parameter :: mixrat_unit = (/ 1.e6, 1.e6 /)
  character(len=3), dimension(maxtrace), parameter :: mixrat_unit_name = (/'ppm', 'ppm'/)
  ! some tracer emissions are best expressed in Tg/yr, some others in Kg/yr
  real,dimension(maxtrace),parameter :: emis_unit = (/ 1.0e-12, 1.0e-12 /)
  character(len=2), dimension(maxtrace), parameter :: emis_unit_name = (/ 'Pg','Pg' /)

  !   fscale(ntrace): scaling factor for conversion of mixing ratios
  !   in kg tracer per kg air to practical mixing ratio units (e.g. ppm)
  !   In this version: ratio of molecular weight  of tracer to that of air

  real,dimension(maxtrace),parameter :: fscale = xmair/ra(:)*mixrat_unit(:)


  real, parameter                   :: factor_airdens = 1.0D-6/(xmair * amu * grav) ! without the 'D-6', we would get molec/m^3

  ! tracer indices on which dry depositions should be applied:
  ! deposition parameters:
  integer, parameter   :: ndep = 0

end module chem_param
