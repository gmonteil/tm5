!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module dry_deposition
!  !
!  ! Purpose:
!  ! --------
!  !
!  ! This module  contains all subroutines needed for dry deposition calculations.
!  ! The mean purpose is to perform on a high resolution of  1x1 degree:
!  !  
!  ! 0. allocate the vd on the model resolution
!  !
!  ! 1. read fields needed for further calculations in subroutines:
!  !    
!  ! 2. calculate the tracer independent friction velocity, aerodynamic 
!  !    and sub-laminar resistance 
!  !
!  ! 3. calculate fields needed for resistance calculations 
!  !
!  ! 4. the tracer dependent vd
!  !
!  ! 5. coarsen the vd
!  !
!  ! 6. apply the coarsened deposition velocities to concentration field rm
!  !
!  ! 7. deallocate vd
!  !
!  ! Reference:
!  ! ---------
!  ! general documentationon ECMWF fields can be found on:
!  ! http://www.ecmwf.int/research/ifsdocs/PHYSICS/
!  ! ---------
!  ! Authors:
!  ! -------
!  ! original code by Laurens Ganzeveld and Ad Jeuken (1997)
!  ! revised and adapted for use in TM5 and use of ECMWF data 
!  ! by F. Dentener and M. Krol (2003)
!  !------------------------------------
!  !
!  use chem_param, only: io3,ihno3,ino,ino2,iso2,iso4,ih2o2,iald2,ich2o,ich3o2h,&
!       inh3,ipan,ico,irooh,ino3,ihno4,iorgntr,in2o5,ihno3,inh4, &
!       imsa,ino3_a,io3s,names, ntrace, emis_data
!  use dims     , only: nlon360,nlat180,nregions,okdebug
!
!  implicit none
!
!  ! interface with other modules
!
!  private
!
!  public :: declare_vd, free_vd, dd_surface_fields
!  public :: vd
!
!  ! module variables
!
!  integer, parameter                     :: ndep = 14
!  integer, parameter, dimension(ndep)    :: &
!       idep = (/ io3, ihno3, ino, ino2, iso2, iso4, ih2o2, iald2, &
!                 ich2o, ich3o2h, inh3, ipan, ico, ihno4 /)
!
!  ! vd      : final deposition velocities for all species on model resolution.
!  ! vd_temp : to store the vd temporarily for one specie
!  type(emis_data),dimension(nregions,ntrace) :: vd
!  type(emis_data),dimension(nregions)        :: vd_temp  
!
!
!contains 
!
!
!
!  subroutine declare_vd
!    !
!    ! Allocate emission data structure "vd",
!    ! final deposition velocities for all species on model resolution
!    !
!    use dims , only: im,jm,nregions
!
!    implicit none
!
!    ! local
!    integer  :: imr,jmr,region,n
!
!    ! start
!    do region = 1, nregions
!       imr = im(region) ; jmr = jm(region)
!       allocate(vd_temp(region)%surf(imr,jmr))
!       do n=1,ntrace
!          allocate(vd(region,n)%surf(imr,jmr))
!          vd(region,n)%surf = 0.0 
!       end do
!    end do
!
!  end subroutine declare_vd
!
!
!
!  subroutine free_vd
!    !
!    ! De-allocate emission data structure "vd",
!    ! final deposition velocities for all species on model resolution
!    !
!    use dims , only: nregions
!
!    implicit none
!
!    ! local
!    integer  :: region, n
!
!    ! start
!    do region = 1, nregions
!       deallocate(vd_temp(region)%surf)
!       do n=1,ntrace
!          deallocate(vd(region,n)%surf)
!       end do
!    end do
!
!  end subroutine free_vd
!
!
!  subroutine dd_surface_fields
!
!    !------------------------------------
!    !
!    ! Purpose:
!    ! -------
!    ! this subroutine reads and prepares all surface fields
!    !
!    ! External
!    ! --------
!    ! dd_get_3_hourly_surface_fields
!    ! dd_calc_ustar_raero_rb
!    ! dd_calc_rstom_rahcan
!    ! dd_calc_inisurf
!    ! dd_calc_rs
!    ! dd_coarsen_vd
!    !
!    ! Reference
!    ! ---------
!    ! Ganzeveld and Lelieveld (1996) and references therein.
!    !
!    !------------------------------------
!    use io_save,         only : readtm3hdf
!    use binas,           only : vkarman
!    use dims,            only : im, jm, idate
!    use photolysis_data, only : phot_dat
!#ifdef MPI
!    use mpi_const,only : myid, root_k, ierr, com_lev, my_real
!#endif
!
!    implicit none
!    !
!    ! input/output
!    ! fsurf: the already opened surface file (can be changed!)
!
!    ! local
!    character(len=10)                            :: c_time
!    real,dimension(:,:,:),allocatable            :: soilph
!    real,dimension(:,:),allocatable              :: lai, lai1
!    real,dimension(:,:),allocatable              :: vgrat_low, vgrat_high
!    integer(kind=1),dimension(:,:), allocatable  :: lsm                     
!    real,dimension(:,:),allocatable              :: sr_ecm, sr_ols, ci, sd, swvl1
!    real,dimension(:,:),allocatable              :: sstr, wind10m, slhf, sshf
!    real,dimension(:,:),allocatable              :: src,d2m,t2m,ssr
!    real,dimension(:,:),allocatable              :: ustar_loc,rb,raero
!    real,dimension(:,:,:),allocatable            :: vd11  
!    real,dimension(:,:),allocatable              :: rahcan
!    real,dimension(:,:),allocatable              :: rstom
!    real,dimension(:,:),allocatable              :: snow_cover
!    real,dimension(:,:),allocatable              :: wet_skin
!    real,dimension(:,:),allocatable              :: fws
!    real,dimension(:,:),allocatable              :: rh2m
!    real,dimension(:,:),allocatable              :: ags
!
!    real,dimension(ntrace)     :: rsoil,rws,rwat,rsnow,diffcf,diffrb,rmes,rcut  
!
!    integer :: i, nsend, region,itrace
!
!#ifdef MPI
!    if(myid == root_k) then    ! reading /calculations by root!
!#endif
!
!       allocate(soilph(nlon360,nlat180,5))
!       allocate(lai(nlon360,nlat180))
!       allocate(lai1(nlon360,nlat180))
!       allocate(vgrat_low(nlon360,nlat180))
!       allocate(vgrat_high(nlon360,nlat180))
!       allocate(lsm(nlon360,nlat180))
!       allocate(sr_ecm(nlon360,nlat180))
!       allocate(sr_ols(nlon360,nlat180))
!       allocate(ci(nlon360,nlat180))
!       allocate(sd(nlon360,nlat180))
!       allocate(swvl1(nlon360,nlat180))
!       allocate(sstr(nlon360,nlat180))
!       allocate(wind10m(nlon360,nlat180))
!       allocate(ssr(nlon360,nlat180))
!       allocate(slhf(nlon360,nlat180))
!       allocate(sshf(nlon360,nlat180))
!       allocate(src(nlon360,nlat180))
!       allocate(d2m(nlon360,nlat180))
!       allocate(t2m(nlon360,nlat180))
!       allocate(ustar_loc(nlon360,nlat180))
!       allocate(raero(nlon360,nlat180))
!       allocate(rb(nlon360,nlat180))
!       allocate(rahcan(nlon360,nlat180))
!       allocate(rstom(nlon360,nlat180))
!       allocate(snow_cover(nlon360,nlat180))
!       allocate(wet_skin(nlon360,nlat180))
!       allocate(fws(nlon360,nlat180))
!       allocate(rh2m(nlon360,nlat180))
!       allocate(ags(nlon360,nlat180))
!       allocate(vd11(nlon360,nlat180,ndep))
!
!       ! -- for output of time header on debug fields
!
!       write(c_time,'(i4,3(i2))') idate(1:4)   !CMK corrected
!       do i=1,10
!          if (c_time(i:i)==' ') c_time(i:i)='0'
!       end do
!
!       !
!       ! surface data sets may have the following characteristics
!       ! 1. instanteneous. 
!       !      The time written in the file is valid from t-1.5 to t+1.5
!       ! 2. Accumulated. 
!       !      The time written in the file is valid from t to t+3
!       ! 3. Daily averaged. 
!       !      The time on the file is always 00 hours, and valid until t+24
!       !
!       ! ***   It is essential to understand this timing when 
!       !       reading and applying these sets.
!       !
!       ! fd2mk it has to be decided to 'save' fields like 
!       ! vgrat_low and daily average surface fields
!       ! for later use, or to read it every 3 hours time step
!       !
!
!       call dd_get_soilph_lai         ! non-ECMWF provided data: soilph and lai
!
!       call dd_get_surface_vegetation      ! vgrat_low, vgrat_high and lsm
!
!       call dd_get_daily_av_surf_fields    ! sr_ecm, sr_ols, ci, sd, swvl1
!
!       call dd_get_3_hourly_surface_fields ! sstr,wind10m,slhf,sshf,src,d2m,t2m
!
!       call dd_calc_ustar_raero_rb         ! raero, ustar_loc, rb
!
!       call dd_calc_rstom_rahcan           ! rstom, rahcan
!
!       call dd_calc_inisurf                ! snow_cover, wet_skin, fws, rh2m, ags
!
!       !okdebug = .true.
!       call dd_calc_rs                     ! vd
!
!       call dd_coarsen_vd                  ! coarsen vd to model resolution
!       !okdebug = .false.
!
!       deallocate(soilph)
!       deallocate(lai)
!       deallocate(lai1)
!       deallocate(vgrat_low)
!       deallocate(vgrat_high)
!       deallocate(sr_ecm)
!       deallocate(sr_ols)
!       deallocate(ci)
!       deallocate(sd)
!       deallocate(swvl1)
!       deallocate(sstr)
!       deallocate(wind10m)
!       deallocate(ssr)
!       deallocate(slhf)
!       deallocate(sshf)
!       deallocate(src)
!       deallocate(d2m)
!       deallocate(t2m)
!       deallocate(ustar_loc)
!       deallocate(raero)
!       deallocate(rb)
!       deallocate(rahcan)
!       deallocate(rstom)
!       deallocate(snow_cover)
!       deallocate(wet_skin)
!       deallocate(fws)
!       deallocate(rh2m)
!       deallocate(ags)
!       deallocate(vd11)
!
!#ifdef MPI
!    end if   !myid == root_k
!#endif
!
!    !ALL PEs from here:
!
!#ifdef MPI
!    do region = 1, nregions
!       nsend = im(region)*jm(region) 
!       do itrace = 1, ntrace
!          call mpi_bcast(vd(region,itrace)%surf, nsend, my_real, root_k, &
!               com_lev, ierr)
!       end do
!       call mpi_bcast(phot_dat(region)%albedo ,nsend, my_real, root_k, &
!            com_lev, ierr)
!    end do
!    if ( myid == root_k .and. okdebug ) &
!         print *, 'dd_surface_fields: Deposition velocities broadcasted' 
!#endif
!
!  contains
!    !
!    ! 
!    subroutine dd_get_soilph_lai
!      !---------------------------
!      !
!      ! Purpose
!      ! -------
!      ! Reads two independent datasets on soils and vegetation not provide by ECMWF
!      ! 
!      ! Reference:
!      ! ---------
!      ! Data provided by Laurens Ganzeveld
!      !
!      !--------------------------
!      use dims, only: datadir
!
!      integer,parameter              :: rank2=2,level1=1,annual_field=0
!      character(len=5),dimension(5)  :: &
!           name_soil = (/'spec1','spec2','spec3','spec4','spec5'/)
!      integer :: i,mo,n
!
!      !
!      !  -- Reading of input file with LAI data which is derived
!      ! from the Olson ecosystem database (0.5 X 0.5 degrees), which
!      ! discerns 46 ecosystems and their characteristics, e.g. LAI. 
!      ! Fields are monthly averaged.
!      ! 
!      i = len_trim(datadir)
!      mo=idate(2)
!
!      call readtm3hdf(datadir(1:i)//'lsmlai.hdf',rank2, &
!           nlon360,nlat180,level1,mo-1,lai,'spec1')  !lai index
!
!      if ( okdebug ) call dd_field_statistics('lai',lai,nlon360,nlat180)
!      !cmk if (okdebug) call dumpfield(0,'lai.hdf',lai,'lai')
!
!      !
!      ! -- Reading of input file with soil pH data,which are derived from
!      !    a 0.5 x 0.5 degree soil pH database, Batjes, 1995  
!      !    over land soilph 1 to 5 should add up to 1.
!      !    SOILPH(I,J,1) - soil pH <5.5
!      !    SOILPH(I,J,2) - soil 5.5 <pH <7.3
!      !    SOILPH(I,J,3) - soil 7.3< pH <8.5
!      !    SOILPH(I,J,4) - soil 8.5 <pH
!      !    SOILPH(I,J,5) - soil 4 < pH <8.5
!      !
!
!      do n=1,5 
!         call readtm3hdf(datadir(1:i)//'soilph.hdf',rank2, &
!              nlon360,nlat180,level1,annual_field,soilph(:,:,n),name_soil(n))
!      end do
!      if ( okdebug ) write(*,*) 'dd_get_soilph_lai: end'
!
!    end subroutine dd_get_soilph_lai
!    !
!    !
!    subroutine dd_get_surface_vegetation
!      !-----------------------------------
!      ! Purpose
!      ! -------
!      ! read ECMWF dataset for vegetation and landmask
!      ! 
!      ! External
!      ! --------
!      ! read_meteo
!      ! dd_field_statistics and dd_field_statistics_i8
!      !
!      ! Reference
!      ! ----------
!      ! None
!      !------------------------
!
!      !
!      ! routine to read fields that need to be initialised only once
!      !
!      use MeteoFiles, only : ReadSurfaceField
!      use datetime,   only : date2tau, tau2date
!      use toolbox,    only : escape_tm
!
!      implicit none
!      integer,parameter :: nveg=20            ! number of surface types in ECMWF
!      !
!      ! following parameters are directly from Table 7.1 
!      ! www.ecmwf.int/research/ifsdocs/Physics
!      !
!      ! cveg is fractional coverage given a certain type of vegetation present
!      real,dimension(nveg),parameter       :: cveg=& 
!           (/0.9, 0.85,  0.9,  0.9,  0.9, &
!             0.99, 0.7,  0. ,  0.5,  0.9, &
!             0.1,   9.,  0.6,  9.0,  9.0, &
!             0.5,  0.5,  0.9,  0.9,  0.6 /)
!      !high_low: e.g. high are trees, low schrub
!      character(len=1),dimension(nveg),parameter :: high_low=&
!           (/'L','L','H','H','H', &
!             'H','L','O','L','L', &
!             'L','O','L','O','O', &
!             'L','L','H','H','O' /)
!      !lai_veg assigned lai per type of vegetation. No seasonal cycle
!      real,dimension(nveg),parameter       :: lai_veg=&
!           (/3.0, 2.0, 5.0, 5.0, 5.0, &
!             6.0, 2.0, 0.5, 1.0, 3.0, &
!             0.5,-9.0, 4.0,-9.0,-9.0, &
!             3.0, 1.5, 5.0, 2.5, 4.0 /)
!      !evergreen is the vegetation seasonal or evergren
!      integer,dimension(nveg),parameter :: evergreen =& 
!           (/0,0,1,0,1, &
!             0,0,0,0,0, &
!             0,0,0,0,0, &
!             1,0,0,0,0 /)
!      character(len=30),dimension(nveg),parameter :: vegtype=(/&
!           'Crops and mixed farming       ','Short grass                   ',&
!           'Evergreen needle leaf trees   ','Deciduous needle leaf trees   ',&
!           'Evergreen broadleaf trees     ','Deciduous broad leaf trees    ',&
!           'Tall grass                    ','Desert                        ',&
!           'Tundra                        ','Irrigated crops               ',&
!           'Semidesert                    ','Ice caps and glaciers         ',&
!           'Bogs and marshes              ','Inland water                  ',&
!           'Ocean                         ','Evergreen shrubs              ',&
!           'Deciduous shrubs              ','Mixed forest/woodland         ',&
!           'Interrupted forest            ','Water and land mixtures       '/)
!      !
!      ! --
!      !
!      integer(kind=1),dimension(:,:,:),allocatable :: tv
!      real, dimension(:,:), allocatable            :: cvl,cvh
!      integer,dimension(6)                         :: idater_acc,idater
!
!      character(len=2),dimension(nveg),parameter:: tvname=&
!           (/'01','02','03','04','05', &
!             '06','07','08','09','10', &
!             '11','12','13','14','15', &
!             '16','17','18','19','20' /)
!
!      integer    :: i,j,nv,itaux
!
!      integer            ::  status
!      real, allocatable  ::  field_r(:,:)
!
!      ! start
!
!      allocate( field_r(nlon360,nlat180) )
!
!      call date2tau(idate,itaux)
!      call tau2date(itaux-21*3600,idater)
!      idater(4) = 21   ! that is the date stored in the file: 21 previous day
!      call ReadSurfaceField( field_r, 'lsm', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_surface_vegetation: error from ReadSurfaceField' )
!      lsm = int(field_r,kind=1)
!
!      if ( okdebug ) call dd_field_statistics_i8('lsm',lsm,nlon360,nlat180)
!
!      allocate(tv(nlon360,nlat180,nveg))
!      allocate(cvl(nlon360,nlat180))
!      allocate(cvh(nlon360,nlat180))
!
!      do nv=1,nveg
!         if (high_low(nv)/='O') then
!
!            call ReadSurfaceField( field_r, 'tv'//tvname(nv), idater, 24, status )
!            if (status/=0) call escape_tm( 'dd_get_surface_vegetation: error from ReadSurfaceField' )
!            tv(:,:,nv) = int(field_r,kind=1)
!            if ( okdebug )  &
!              call dd_field_statistics_i8('tv'//tvname(nv),tv(:,:,nv),nlon360,nlat180)
!
!         end if
!      end do
!
!      ! fraction low vegation
!      call ReadSurfaceField( cvl, 'cvl', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_surface_vegetation: error from ReadSurfaceField' )
!      if ( okdebug ) call dd_field_statistics('cvl',cvl,nlon360,nlat180)
!
!      ! fraction high vegetation
!      call ReadSurfaceField( cvh, 'cvh', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_surface_vegetation: error from ReadSurfaceField' )
!      if ( okdebug ) call dd_field_statistics('cvh',cvh,nlon360,nlat180)
!
!      !
!      ! -- vgrat_low :  fractional coverage of low vegetation
!      ! -- vgrat_high : fractional coverage of high vegetation
!      ! -- lai1      : lai derived following ECMWF, 
!      !                unfortunately there is no yearly cycle on this product
!      ! -- fd Note:
!      !    ewsmx : maximum field capacity could possibly be approximated 
!      !            with the root depth
!      !  
!      !   ewsmx(:,:)=ewsmx(:,:)+cvl(:,:)*tv(:,:,nv)/100.*
!      !              root_depth(nv)/depth_tessel
!      !   ewsmx(:,:)=ewsmx(:,:)+cvh(:,:)*tv(:,:,nv)/100.*
!      !              root_depth(nv)/depth_tessel
!      !   Note that in the ECMWF model a constant value of 0.353 m3/m3 is used
!      !   depth_tessel = 2.89  ! tessel soil model has four layers, 
!      !   the lowest of four ends at 2.89 m 
!      !   root depth is derived from ECMWF root distribution, 
!      !   and is used to calculate wilting point
!      !   real,dimension(nveg),parameter :: root_depth=(/& 
!      !   0.4017,0.3455,0.4223,0.4391,0.4521,0.5479,0.4401,0.0700,0.1850,0.4017,&
!      !   0.6737,0.0000,0.4912,0.0000,0.0000,0.5156,0.5156,0.5350,0.5350,0.0000/)
!      ! -- fd not implemented
!
!      vgrat_high=0.
!      vgrat_low=0.
!      lai1=0.
!      do nv=1,20 
!         if (high_low(nv)=='L') then
!            vgrat_low(:,:)=vgrat_low(:,:)+cvl(:,:)*tv(:,:,nv)/100.*cveg(nv)
!            lai1(:,:)=lai1(:,:)+cvl(:,:)*tv(:,:,nv)/100.*cveg(nv)*lai_veg(nv) 
!         end if
!         if (high_low(nv)=='H') then
!            vgrat_high(:,:)=vgrat_high(:,:)+cvh(:,:)*tv(:,:,nv)/100.*cveg(nv)
!            lai1(:,:)=lai1(:,:)+cvh(:,:)*tv(:,:,nv)/100.*cveg(nv)*lai_veg(nv)
!         end if
!      end do !nv
!
!      if ( okdebug ) call dd_field_statistics('vgrat_low',vgrat_low, &
!           nlon360,nlat180)
!      !cmk if ( okdebug ) call dumpfield(0,'vgrat_low.hdf',vgrat_low,'vgrat_low')
!      if ( okdebug ) call dd_field_statistics('vgrat_high',vgrat_high, &
!           nlon360,nlat180)
!      !cmk if ( okdebug ) &
!      !cmk   call dumpfield(0,'vgrat_high.hdf',vgrat_high,'vgrat_high')
!      if ( okdebug ) call dd_field_statistics('lai1',lai1,nlon360,nlat180)
!      !cmk   if ( okdebug ) call dumpfield(0,'lai1.hdf',lai1,'lai1')
!
!      deallocate(cvl)
!      deallocate(cvh)
!      deallocate(tv)
!      if ( okdebug ) write(*,*) 'dd_get_surface_vegetation: end'
!
!      deallocate( field_r )
!
!    end subroutine dd_get_surface_vegetation
!    !
!    !
!    subroutine dd_get_daily_av_surf_fields
!      !---------------------------------------------
!      !
!      ! Purpose
!      ! -------
!      ! read daily averaged ECMWF datasets
!      ! 
!      ! External
!      ! --------
!      ! read_meteo
!      ! dd_field_statistics and dd_field_statistics_i8
!      !
!      ! Reference
!      ! ----------
!      ! None
!      !------------------------
!      use MeteoFiles, only : ReadSurfaceField
!      use datetime,   only : date2tau, tau2date
!      use toolbox,     only : escape_tm
!
!      integer :: itaux
!      integer,dimension(6) :: idater
!
!      integer            ::  status
!
!      call date2tau(idate,itaux)
!      call tau2date(itaux-21*3600,idater)
!
!      ! large scale surface roughness for momentum [m]
!      idater(4) = 21
!      call ReadSurfaceField( sr_ecm, 'sr_ecm', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_daily_av_surf_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('sr_ecm',sr_ecm,nlon360,nlat180)
!
!      !surface roughness from Olson data base [m]
!      call ReadSurfaceField( sr_ols, 'sr_ols', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_daily_av_surf_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('sr_ols',sr_ols,nlon360,nlat180)
!
!      ! sea ice cover [fraction,0-1]
!      call ReadSurfaceField( ci, 'ci', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_daily_av_surf_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('ci',ci,nlon360,nlat180)
!
!      ! snow depth [m]
!      call ReadSurfaceField( sd, 'sd', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_daily_av_surf_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('sd',sd,nlon360,nlat180)
!
!      !soil water volume [m3/m3]
!      call ReadSurfaceField( swvl1, 'swvl1', idater, 24, status )
!      if (status/=0) call escape_tm( 'dd_get_daily_av_surf_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('swvl1',swvl1,nlon360,nlat180)
!
!      if ( okdebug ) write(*,*) 'end of dd_get_daily_average surface_fields'
!
!    end subroutine  dd_get_daily_av_surf_fields
!
!    !
!    !
!
!    subroutine dd_get_3_hourly_surface_fields
!      !----------------------------------------
!      !
!      ! Purpose
!      ! -------
!      ! read 3-hourly ECMWF datasets, can be instanteneous or acccumulated
!      ! 
!      ! External
!      ! --------
!      ! date2tau
!      ! tau2date
!      ! read_meteo
!      ! dd_field_statistics and dd_field_statistics_i8
!      !
!      ! Reference
!      ! ----------
!      ! None
!      !------------------------
!      use MeteoFiles, only : ReadSurfaceField
!      use datetime,   only : date2tau, tau2date
!      use toolbox,     only : escape_tm
!
!      integer,dimension(6):: idater_acc,idater
!      real,dimension(:,:),allocatable :: u10m,v10m
!      integer :: itaux ! auxiliary variable
!
!      integer            ::  status
!
!      allocate(u10m(nlon360,nlat180))
!      allocate(v10m(nlon360,nlat180))
!
!      call date2tau(idate,itaux)
!      call tau2date(itaux,idater_acc)
!      call tau2date(itaux+3600*3,idater)
!      if ( okdebug ) write(*,*) 'dd_get_3_hourly_surface_fields read time',&
!           ' on instanteneous field',idater
!      if ( okdebug ) write(*,*) 'dd_get_3_hourly_surface_fields read time',&
!           ' on accumulated field',idater_acc
!      !
!      ! read instanteneous fields
!      !
!
!      ! temperature at 2 m [K]
!      call ReadSurfaceField( t2m, 't2m', idater, 0, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('t2m',t2m,nlon360,nlat180)
!
!      ! dewpoint at 2 m [K]
!      call ReadSurfaceField( d2m, 'd2m', idater, 0, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('d2m',d2m,nlon360,nlat180)
!
!      ! u-component [East West] component of diagnosed 10m wind [m/s]
!      call ReadSurfaceField( u10m, 'u10m', idater, 0, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('u10m',u10m,nlon360,nlat180)
!
!      ! v-component [North-South] component of diagnosed 10 m wind [m/s]
!      call ReadSurfaceField( v10m, 'v10m', idater, 0, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('v10m',v10m,nlon360,nlat180)
!
!      ! wind at 10 meters is diagnosed from wind at higher levels
!      ! documentation can be found on 
!      ! http:www.ecmwf.int/ifsdocs/research/PHYSICS   chapter 3
!      !
!      wind10m=sqrt(u10m*u10m+v10m*v10m) 
!
!
!      ! skin reservoir content [m]
!      call ReadSurfaceField( src, 'src', idater, 0, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('ci',src,nlon360,nlat180)
!
!      !
!      ! read accumulated fields
!      !
!
!      ! surface solar radiation
!      call ReadSurfaceField( ssr, 'ssr', idater_acc, 3, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('ssr',ssr,nlon360,nlat180)
!
!      ! surface sensible heat flux [W m-2]
!      call ReadSurfaceField( sshf, 'sshf', idater_acc, 3, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('sshf',sshf,nlon360,nlat180)
!
!      ! surface latent heat flux [W m-2]
!      call ReadSurfaceField( slhf, 'slhf', idater_acc, 3, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('slhf',slhf,nlon360,nlat180)
!
!      ! surface stress [Nm-2]
!      call ReadSurfaceField( sstr, 'sstr', idater_acc, 3, status )
!      if (status/=0) call escape_tm( 'dd_get_3_hourly_surface_fields: error from ReadSurfaceField' )
!
!      if ( okdebug ) call dd_field_statistics('sstr',sstr,nlon360,nlat180)
!
!      deallocate(u10m)
!      deallocate(v10m)
!
!      if ( okdebug ) write(*,*) 'end of dd_get_3_hourly_surface_fields'
!
!    end subroutine dd_get_3_hourly_surface_fields
!
!
!    !--------------------------------
!    subroutine dd_calc_ustar_raero_rb
!      !--------------------------------
!      !
!      ! Purpose
!      ! -------
!      ! Calculate friction velocity (ustar), aerodynamic resistance (raero)
!      ! and quasi laminar surface resistance (rb), 
!      !
!      ! Method
!      ! ------
!      ! uses the log normal wind profile information stored by ECMWF 
!      ! in 10 meter wind
!      ! to estimate a local ustar over land
!      ! uses Charnock equation over sea to estimate ustar
!      ! aerodynamic resistance is calculated from heat fluxes and ustar 
!      ! using a constant reference height
!      ! sub laminar resistance from ustar
!      ! 
!      ! External
!      ! --------
!      ! dd_field_statistics 
!      ! dumpfield
!      !
!      ! Reference
!      ! ----------
!      ! Ge Verver (personal communication, 2003)
!      !------------------------
!      use binas, only : grav, vKarman
!      implicit none
!      ! Href: constant reference height for calculations of raero
!      real, parameter :: Href=30.     
!      ! some constants specific to calculation of ustar and raero
!      real,parameter  :: alfa_charnock1 = 0.11
!      real,parameter  :: rhoCp          = 1231.0
!      real,parameter  :: rhoLv          = 3013000.0
!      real,parameter  :: alfa_charnock2 = 0.018
!      real,parameter  :: v_charnock     = 1.5e-5 !(m2/s)
!      real,parameter  :: rz0            = 2.0
!      real,dimension(:,:),allocatable :: sr_mix
!      integer :: i, j
!
!      real :: buoy, tstv, obuk, ra, y0, yra
!      real :: ustar_sea, ustar_land, xland, sr_sea, sr_land, sr_help
!
!      allocate(sr_mix(nlon360,nlat180))
!      do j=1,nlat180
!         do i=1,nlon360
!            xland=lsm(i,j)/100.
!
!            ! SEA:
!            ! surface roughness  from Charnock equation
!            ! friction velocity from surface stress
!            !
!            ustar_sea=sqrt(sstr(i,j))
!            sr_sea = alfa_charnock1*v_charnock/ustar_sea + &
!                     alfa_charnock2*sstr(i,j)/grav 
!
!            !
!            ! LAND
!            ! calculate the 'local' surface roughness for momentum that 
!            ! is consistent with 10 m wind and the Olsson data base, 
!            ! we assume that the windspeed at 75 m is independent of
!            ! surface roughness
!            !
!            if ( xland > 0. ) then
!               sr_land=max(1e-2,sr_ols(i,j))  ! occurs at Islands, etc.
!               sr_help=min(sr_ecm(i,j),0.03)
!               !fd ustar_land=vKarman*wind10m(i,j)/alog(10./sr_help)
!               ! !ustar consistent with 'clipped' large scale roughness
!               ustar_land=vKarman*wind10m(i,j)/alog(10./sr_help)*&
!                    alog(75./sr_help)/alog(10./sr_land)
!            end if
!
!            ustar_loc(i,j)=xland*ustar_land+(1-xland)*ustar_sea
!            sr_mix(i,j) =xland*sr_land  +(1-xland)*sr_sea
!
!         end do !i
!      end do !j
!
!      !
!      ! Calculation of quasi laminar resistance of the layer above the surface.
!      ! E.g. Walcek, Atmos. Environ, 20, pp. 949-964, 1986, and earlier refs.
!      !     
!      do j=1,nlat180
!         do i=1,nlon360
!            rb(i,j)=(rz0/(ustar_loc(i,j)*vkarman)) 
!         end do
!      end do
!
!      !
!      ! calculate the aerodynamic resistance
!      !
!      do j=1,nlat180
!         do i=1,nlon360
!            buoy=-sshf(i,j)/rhoCp-0.61*t2m(i,j)*slhf(i,j)/rhoLv
!            tstv=-buoy/ustar_loc(i,j) 
!            obuk=1e-6
!            if (tstv/=0) then 
!               obuk=ustar_loc(i,j)*ustar_loc(i,j)*t2m(i,j)/(tstv*grav*vKarman)
!            end if
!
!            if ( obuk> 0. ) then !   stable conditions
!               ra = 0.74*( alog(Href/sr_mix(i,j)) + 6.4*(Href-sr_mix(i,j))/obuk )&
!                    / (vKarman*ustar_loc(i,j)) 
!            else              !   unstable
!               y0 = sqrt(1.-9.*sr_mix(i,j)/obuk)+1.       
!               yra = sqrt(1.-9.*Href/obuk)+1. 
!               ra = 0.74*(alog(Href/sr_mix(i,j))+2.*(alog(y0)-alog(yra)))/ &
!                    (vKarman*ustar_loc(i,j)) 
!            end if
!            raero(i,j)=max(10.,min(ra,1e5))  !
!         end do !i
!      end do !j
!
!      if ( okdebug ) &
!           call dd_field_statistics('ustar_loc',ustar_loc,nlon360,nlat180)
!      !cmk if ( okdebug ) &
!      !cmk  call dumpfield(0,'ustar_loc'//c_time//'.hdf',ustar_loc,'ustar_loc')
!      if ( okdebug ) call dd_field_statistics('raero',raero,nlon360,nlat180)
!      !cmk  if ( okdebug ) call dumpfield(0,'raero'//c_time//'.hdf',raero,'raero')
!      if ( okdebug ) call dd_field_statistics('rb',rb,nlon360,nlat180)
!      !cmk  if ( okdebug ) call dumpfield(0,'rb'//c_time//'.hdf',rb,'rb')
!      if ( okdebug ) write(*,*) 'end calc_aero_ustar'
!
!      deallocate(sr_mix)
!
!    end subroutine dd_calc_ustar_raero_rb
!    !
!    ! 
!    subroutine dd_calc_rstom_rahcan
!      ! -----------------------------
!      ! Purpose
!      ! ---------
!      ! Calculate water vapour stomatal resistance from the PAR
!      ! (Photosythetically Active Radiation) which is 0.55 times
!      ! the net short wave radiation (ssr) at the surface.
!      ! Calculate rahcan from ustar and lai
!      ! 
!      ! External
!      ! --------
!      ! none
!      !
!      ! Reference
!      ! ----------
!      ! none
!      !--------------------------------------------------
!      !
!      !  -- constants used for the calculation of the stomatal resistance 
!      !     (see ECHAM3 manual)
!      implicit none
!      real,parameter  :: a=5000.
!      real,parameter  :: b=10.
!      real,parameter  :: c=100.
!      real,parameter  :: vk=0.9
!      real,parameter  :: vlt=1.    
!      !
!      real, parameter :: foresth=20. 
!      integer         :: j, i
!      real            :: vpar, d
!      real            :: canht,laihelp
!
!      !
!      !  -- recalculated rstom for a LAI of 1 (see ECHAM3 manual)
!      !
!      do j=1,nlat180
!         do i=1,nlon360
!            !
!            !    --  calculation of PAR from net short wave radiation
!            !    --  radiation is corrected for daily cycle
!            !
!            rstom(i,j)=1e5    !high resistance during the night
!            if (ssr(i,j) > 0.) then
!               vpar=max(1.,0.55*ssr(i,j))
!               d=(a+b*c)/(c*vpar) 
!               rstom(i,j)=(vk*c)/((b/(d*vpar))*alog((d*exp(vk*vlt)+1.)/(d+1.))-&
!                    alog((d+exp(-vk*vlt))/(d+1.)))
!            end if !ssr(i,j) > 0.
!         end do !i
!      end do !j
!      !
!      ! -- computation of in-canopy aerodynamic resistance from canopy height
!      !    and the friction velocity and the Leaf Area Index 
!      !    (see Erisman & Van Pul)
!      !
!      do j=1,nlat180 
!         do i=1,nlon360
!            !
!            !
!            ! -- calculation of canopy height CANHT from effective fraction 
!            !    of high vegetation assuming that forest has a canopy height 
!            !    of 20 m (other vegetation 0 m)
!            !    
!            canht= vgrat_high(i,j)*foresth
!            ! laihelp: the maximum values derived from ECMWF are more realistic    
!            laihelp=min(lai1(i,j),lai(i,j))
!            rahcan(i,j)=max(1.e-5,14.*laihelp*canht/ustar_loc(i,j))
!         end do !j
!      end do !i
!
!      !cmk   if ( okdebug ) &
!      !cmk    call dumpfield(0,'rahcan'//c_time//'.hdf',rahcan,'rahcan')
!      !cmk   if ( okdebug ) call dumpfield(0,'rstom'//c_time//'.hdf',rstom,'rstom')
!      if ( okdebug ) write(*,*) 'dd_calc_rstom_rahcan: end'
!      !
!    end subroutine dd_calc_rstom_rahcan
!    !
!    !
!    subroutine dd_calc_inisurf
!      ! ------------------------
!      ! Purpose
!      !---------
!      ! Calculate some surface fields later needed in the calculation
!      ! 
!      ! External
!      ! --------
!      ! none
!      !
!      ! Reference
!      ! ----------
!      ! none
!      !--------------------------------------------------
!      implicit none
!      real,parameter    :: sncr     = 0.015       ! critical snow cover depth
!      real,parameter    :: tstar    = 273.
!      real,parameter    :: wlmax    = 2.e-4              
!      real,parameter    :: ewsmx    = 0.323
!      real,parameter    :: ewsmx_sat= 0.472
!      real,parameter    :: wpwp     = 0.171
!      !
!      ! for albedo calculations
!      real              :: snow, freesea, bare, desert
!      real              :: e,es,wcr,wlmx,xland,ahelp,vgr,laihelp
!      integer           :: i,j    ! auxiliary variables
!      !
!      ! start
!      ! 
!      rsoil=0.0
!      rws=0.0
!      rwat=0.0
!      rsnow=0.0
!      diffcf=0.0
!      diffrb=0.0
!      rmes=0.0
!      rcut=0.0
!      !
!      !  -- definition of terms which correct for the diffusivity 
!      ! for the computation of the boundary layer resistance, (v/D)**2/3
!      ! v = 0.189 cm2 s-1 (heat) and D from DH2O/DX where DH2O=0.212
!      !
!      diffrb(io3)     = 1.2
!      diffrb(ihno3)   = 1.4
!      diffrb(ino)     = 1.1
!      diffrb(ino2)    = 1.2
!      diffrb(iso2)    = 1.4
!      diffrb(iso4)    = 999.
!      diffrb(ih2o2)   = 1.2
!      diffrb(iald2)   = 1.4
!      diffrb(ich2o)   = 1.1
!      diffrb(ich3o2h) = 1.3
!      diffrb(inh3)    = 0.9
!      diffrb(ipan)    = 1.7
!      diffrb(ico)     = 1.2 
!      !  DIFFRB(ich3coo2h)=1.5,DIFFRB(ihcooh)=1.3,DIFFRB(ihno2)=1.3
!      ! not needed in this version
!
!      !  -- parallel soil resistance, the resistances of the trace gases
!      ! of the Wesely scheme are calculated assuming a O3 soil resistance
!      ! of 400 s m-1, for SO2 rsoil=100 s m-1 and the wet skin resistance of
!      ! O3 is 2000 s m-1 whereas rws SO2 is 100 s m-1.
!
!      rsoil(io3)      = 400.
!      rsoil(ihno3)    = 1.
!      rsoil(ino)      = 1.e+5
!      rsoil(ino2)     = 600.
!      rsoil(iso2)     = 100.
!      rsoil(iso4)     = 999.
!      rsoil(ih2o2)    = 80.
!      rsoil(iald2)    = 1.e+5
!      rsoil(ich2o)    = 1666.
!      rsoil(ich3o2h)  = 3650.
!      rsoil(inh3)     = 100.
!      rsoil(ipan)     = 3994.
!      rsoil(ico)      = 5000.
!      !      RSOIL(ich3coo2h)=3290./,RSOIL(ihcooh)=1./,RSOIL(ihno2)=97.
!      !
!      !  -- sea water resistance, which is generally similar to the wet skin
!      ! resistance except of that for SO2 due to a different pH of sea water
!      ! compared to the pH of a wet canopy 
!      ! note for ammonia ocean emissions we use 2-way exchange approach FD
!      !
!      rwat(io3)       = 2000.
!      rwat(ihno3)     = 1.
!      rwat(ino)       = 1.e+5
!      rwat(ino2)      = 3000.
!      rwat(iso2)      = 1.
!      rwat(iso4)      = 999.
!      rwat(ih2o2)     = 72.
!      rwat(iald2)     = 300.
!      rwat(ich2o)     = 254.
!      rwat(ich3o2h)   = 293.
!      rwat(inh3)      = 1.
!      rwat(ipan)      = 295.
!      rwat(ico)       = 1e5
!      !      RWAT(ich3coo2h)=291.,RWAT(ihcooh)=1.,RWAT(ihno2)=75.
!      !
!      !  -- wet skin reservoir resistance
!      !
!      rws(io3)        = 2000.
!      rws(ihno3)      = 1.
!      rws(ino)        = 1.e+5
!      rws(ino2)       = 3000.
!      rws(iso2)       = 100.
!      rws(iso4)       = 999.
!      rws(ih2o2)      = 72.
!      rws(iald2)      = 300.
!      rws(ich2o)      = 254.
!      rws(ich3o2h)    = 293.
!      rws(inh3)       = 1.
!      rws(ipan)       = 295.
!      rws(ico)        = 1e5
!      !        RWS(ich3coo2h)=291., RWS(ihno2)=75.,RWS(ihcooh)=1.
!      !
!      !  -- snow resistance, the snow/ice resistances of the trace gases
!      ! of Wesely's scheme are taken same as that of the rsoil and corrected
!      ! for temperatures smaller than 271 K
!      !
!      rsnow(io3)      = 2000.
!      rsnow(ihno3)    = 1.
!      rsnow(ino)      = 1.e+5
!      rsnow(ino2)     = 3000.
!      rsnow(iso2)     = 1.
!      rsnow(iso4)     = 999.
!      rsnow(ih2o2)    = 80.
!      rsnow(iald2)    = 1.e+5
!      rsnow(ich2o)    = 1666.
!      rsnow(ich3o2h)  = 3650.
!      rsnow(inh3)     = 1e5
!      rsnow(ipan)     = 3394.
!      rsnow(ico)      = 1e5
!      !       RSNOW(ich3coo2h)=3290.,RSNOW(ihcooh)=1.E+5,RSNOW(ihno2)=97.
!      !
!      !  -- mesophyll resistance
!      !
!      rmes(io3)       = 1.
!      rmes(ihno3)     = 1.
!      rmes(ino)       = 500.
!      rmes(ino2)      = 1.
!      rmes(iso2)      = 1.
!      rmes(iso4)      = 999.
!      rmes(ih2o2)     = 1.
!      rmes(iald2)     = 200.
!      rmes(ich2o)     = 1.
!      rmes(ich3o2h)   = 1.
!      rmes(inh3)      = 1.
!      rmes(ipan)      = 1.
!      rmes(ico)       = 5000.
!      !      RMES(ich3coo2h)=1.,RMES(ihcooh)=1.,RMES(ihno2)=1.
!      !
!      !  -- cuticle resistance
!      !
!      rcut(io3)       = 1.e+5
!      rcut(ihno3)     = 1.
!      rcut(ino)       = 1.e+5
!      rcut(ino2)      = 1.e+5
!      rcut(iso2)      = 1.e+5
!      rcut(iso4)      = 999.
!      rcut(ih2o2)     = 1.e+5
!      rcut(iald2)     = 1.e+5
!      rcut(ich2o)     = 1.e+5
!      rcut(ich3o2h)   = 1.e+5
!      rcut(inh3)      = 1e5
!      rcut(ipan)      = 1.e+5
!      rcut(ico)       = 1e5
!      !     RCUT(ich3coo2h)=1.E+5, RCUT(ihcooh)=2500.,RCUT(ihno2)=1.E+5
!      ! 
!      !  -- Diffusivity coefficent, to correct stomatal resistance for differences
!      ! in diffusivity between water vapour and the specific trace gas
!      !
!      diffcf(io3)     = 1.6
!      diffcf(ihno3)   = 1.9
!      diffcf(ino)     = 1.3
!      diffcf(ino2)    = 1.6
!      diffcf(iso2)    = 1.9
!      diffcf(iso4)    = 999.
!      diffcf(ih2o2)   = 1.4
!      diffcf(iald2)   = 1.6
!      diffcf(ich2o)   = 1.3
!      diffcf(ich3o2h) = 1.6
!      diffcf(inh3)    = 1.0
!      diffcf(ipan)    = 2.6
!      diffcf(ico)     = 1.2
!      !
!      !    DIFF(ich3coo2h)=2.0)=,DIFF(ihcooh)=1.6,DIFF(ihno2)=1.6
!      !
!      do j=1,nlat180
!         do i=1,nlon360
!            xland=lsm(i,j)/100.       ! land-sea fraction directly from ECMWF
!            ! the maximum values derived from ECMWF are more realistic
!            laihelp=min(lai1(i,j),lai(i,j))
!
!            !
!            ! -- calculation of monthly snow cover fraction sd using 
!            !    constant critical snow depth
!            !    alternatively this critical snow depth may be chosen as 
!            !    a linear function of ln(sr_ecm)
!            !    B. v.d. Hurk (2002) suggests 0.015 m (=sncr) for z0<0.25
!            !    and 1 m for z0>5m and log-linear in between in between.
!            !    factor 0.9 is introduced to account for the fact that 
!            !    it is unlikely 
!            !    that 'high' vegetation is completely snow covered
!            ! 
!            ! FD25.03.2003
!            snow_cover(i,j)=min(xland,sd(i,j)/sncr)
!
!            !
!            ! -- calculation of wet skin fraction from the skin reservoir 
!            !    content src (prognostic variable in ECMWF)
!            !    the vegetation fractions and their attributed LAI, 
!            !    Wlmax which is the maximum amount of water that can be held 
!            !    on one layer of leaves or bare ground, Wlmx is the 
!            !    maximum skin reservoir content
!            !   
!            if ( xland > 0.0 ) then
!               ! bare soil fraction small discrepancy when lakes are present
!               bare=max(0.,1.-vgrat_high(i,j)-vgrat_low(i,j))
!               wlmx=wlmax*(bare+laihelp)             
!               wet_skin(i,j)=amin1(1.,src(i,j)/wlmx)
!            else
!               wet_skin(i,j)=1.0   !for sea  CMK bug 01-07-2003
!            end if
!            !
!            ! -- calculation of water stress factor FWS with Wsmx is the 
!            !    field capacity.
!            !
!            !    Field capacity is defined as the maximum amount of water 
!            !    that an column of soil can hold
!            !    against gravity 24-48 hours after wetting of the soil
!            !    Wcr is a critical value, Wpwp is the permanent wilting point 
!            !    and Ws is the total amount of water available in the 
!            !    root zone (ECHAM3 manual)
!            !
!            fws(i,j)=1e-5
!            if (xland > 0) then
!               wcr=0.5*ewsmx_sat
!               ! used ewsmx_sat instead of ewsmx, is more consistent 
!               ! with values of swvl1    
!               if ( swvl1(i,j) > wpwp .and. swvl1(i,j) < wcr ) &
!                    fws(i,j)=(swvl1(i,j)-wpwp)/(wcr-wpwp)
!               if ( swvl1(i,j) > wcr ) fws(i,j)=1.
!            end if
!            !
!            ! -- Computation of relative humidity (2 m) out of the air and dew 
!            !    temperature at 2 m  which is used 
!            !    
!            es=exp(19.59*(t2m(i,j)-tstar)/t2m(i,j))
!            e=exp(19.59*(1.-tstar/d2m(i,j)))
!            rh2m(i,j)=e/es
!            !
!            ! calculate albedo (needed for photolysis rates) from 
!            ! different fractions
!            ! 
!            ! vgr: coverage by low and high vegetation
!            vgr= vgrat_low(i,j)+vgrat_high(i,j)
!            !
!            ! --  if high vegetation is present and snow this may alter 
!            !     the effective snow albedo: let high vegetation prevail
!            snow = snow_cover(i,j)- max(0.0,snow_cover(i,j)+vgrat_high(i,j)-1.0) 
!
!            bare = max(0.,1.0 - vgr - snow)
!            ! soils with pH values larger than 7.3 
!            desert = soilph(i,j,3) + soilph(i,j,4)
!            ! when more desert is present than bare land...
!            desert = desert - max(0.0,desert-bare)  
!            bare = bare -desert 
!            freesea = max(0.,1.-xland-ci(i,j))
!            !
!            ! albedo array ags  (for photolysis routine) 
!            ! store 1x1 field in ags1
!            ags(i,j) = freesea*0.05 + ci(i,j)*0.70 + &
!                 xland*(desert*0.10+bare*0.05+vgr*0.01+snow*0.70)   
!         end do
!      end do
!      !cmk  if ( okdebug ) call dumpfield(0,'ags'//c_time//'.hdf',ags,'ags')
!      !cmk  if ( okdebug ) call dumpfield(0,'rh2m'//c_time//'.hdf',rh2m,'rh2m')
!      if ( okdebug ) write(*,*) 'after dd_cal_inisurf'
!    end subroutine dd_calc_inisurf
!    !
!    !
!    subroutine dd_calc_rs
!      !--------------------
!      !
!      ! This subroutine calculates the surface resistance as
!      ! a function of a series of resistances
!      !
!      !  purpose
!      !  -------
!      !  determine surface resistance "rsurf" for dry deposition 
!      !  calculations
!      !
!      !  external
!      !  --------
!      !  
!      !  reference
!      !  ---------
!      !  Ganzeveld and Lelieveld (1996) and references therein
!      !
!      !------------------------------------------------------------------
!      !
!      !  -- resistances and auxiliary variables
!      !
!      implicit none
!      real                             :: xland1,xland2,xland3,xland4
!      real                             :: xwat1,xwat2,xsum
!      real                             :: rstomx,rsveg,ra1
!      real                             :: vdland1,vdland2,vdland3,vdland4
!      real                             :: rssoil,rsws,vdwat1,vdwat2
!      real                             :: rssnow,rveg,rleaf,rcanopy
!      real                             :: w10,zrsnhno3,ustarh
!      real                             :: cvsh,cvwh,evgrath,tsoilph
!      real                             :: xland,ts1,laihelp
!      real,dimension(:,:), allocatable :: rwatso4
!      real,dimension(:,:), allocatable :: rsoilnh3,rsoilso2,rsoilso4
!      real,dimension(:,:), allocatable :: rsnowhno3,rsnowso2
!      !
!      integer :: jt,j,i,jdep
!      character(len=8) :: adum
!
!      allocate(rwatso4(nlon360,nlat180))
!      allocate(rsoilnh3(nlon360,nlat180))
!      allocate(rsoilso2(nlon360,nlat180))
!      allocate(rsoilso4(nlon360,nlat180))
!      allocate(rsnowhno3(nlon360,nlat180))
!      allocate(rsnowso2(nlon360,nlat180))
!      !
!      ! -- extension with the Wesely scheme, 1989 
!      !    (Atmospheric Environ.,1989)
!      !    in which the resistances of trace gases are estimated from 
!      !    the reactivity relative to 
!      !    ozone and the dissolvation relativeto SO2. The deposition process of 
!      !    these two trace is used as a reference for the estimation.
!      !
!      ! calculate some tracer specific resistances
!      !
!      do j=1,nlat180
!         do i=1,nlon360
!            ustarh=ustar_loc(i,j)
!
!            ! -- calculation of the integrated resistance from the mass size 
!            !    distribution for rural continental aerosol with an unimodal distr.
!            !    with the mean mass radius at about 0.4 um. 
!            !    (observations by Mehlmann (1986)
!            !    as referenced in Warneck, 1988) and the friction velocity. 
!            !    Polynominal fit!
!            !    Over land the surface resistance of bare soil and snow 
!            !    covered surfaces is calculated
!            !
!            ! -- following distributions are used:
!            ! 1. deposition over land using  a rural continental size distribution
!            !    i.e. mostly accumulation range aerosol
!            ! 2. deposition over water using a marine size distribution 
!            !    i.e. bimodal distribution with a 30% coarse fraction
!            !    fjd this is most relevant when explicit chemistry in seasalt 
!            !    is considered
!            ! 3. deposition over water using rural continental size distribution
!            !    the deposition of 'seasalt' sulfate may be calculated using
!            !    the relation vd2=0.7*vd1+0.3*vdsssulfate
!            !
!            ! 1st distribution:
!            rsoilso4(i,j) = max( 1., &
!                 100./(0.08-0.57*ustarh+1.66*ustarh**2-0.36*ustarh**3) )
!            ! 2nd distribution:
!            ! rwatso4(i,j) =
!            ! max(1.,100./(0.12-0.24*ustarh+5.25*ustarh**2 -1.84*ustarh**3))
!            ! 3rd distribution: 
!            rwatso4(i,j)  = max( 1., &
!                 100./(0.07-0.1*ustarh+2.1*ustarh**2 -0.20*ustarh**3) )
!            !   
!            ! -- parameterization of snow/ice SO2 resistance as a function of the 
!            !    snow/ice temperature, based on observations by Dasch and Cadle
!            !
!            ts1=max(200.,t2m(i,j))
!            !
!            ! FD25032003 the measurements are not completely consistent; 
!            ! put a lower
!            ! limit of vd=0.10 cm/s for SO2 to avoid excessive accumulation
!            ! if snow is melting high uptake
!            !
!            rsnowso2(i,j)=amin1(max(10.,10.**(-0.09*(ts1-273.)+2.4)),1.e+3)
!            ! snow is melting changed (advise Wouter Greull) FD072003
!            if ( ts1 > (273.+3.) .and. sshf(i,j) > 10. ) rsnowso2(i,j)=50.
!
!            !
!            ! arbitrary attribution of soil resistance to NH3 deposition
!            ! acid soils have lowest deposition
!            !
!            tsoilph=sum(soilph(i,j,1:5))
!            rsoilnh3(i,j)=1e5
!            if (tsoilph > 0.) &
!                 rsoilnh3(i,j)=  max(25.,soilph(i,j,1)*25. +soilph(i,j,2)*25.&
!                 +soilph(i,j,3)*200.+soilph(i,j,4)*200.+soilph(i,j,5)*100.)
!            !
!            ! -- Computation of rsoil for SO2 from soil pH three different 
!            !    rh classes
!            !    The rsoil values for each class are determined out of 
!            !    regression and the average value of pH of the class. 
!            !    Concerning rh, a wet, dry and very dry class 
!            !    is discerned with the threshold value of 60% and 40%
!            ! 
!            rsoilso2(i,j)=1e5
!            !
!            ! -- for rh > 60 %
!            !
!            if ( tsoilph > 0. ) rsoilso2(i,j)=max(25.,(soilph(i,j,1)*115.+&
!                 soilph(i,j,2)*65+soilph(i,j,3)*25.+ &
!                 soilph(i,j,4)*25.+soilph(i,j,5)*70.)/tsoilph+&
!                 amin1(max(0.,1000.*exp(269.-ts1)),1.e+5))
!
!            ! -- for rh less than 60% and larger than 40 %
!
!            if ( rh2m(i,j) < 0.6 .and. rh2m(i,j) > 0.4) &
!                 !cmk    rsoilso2(i,j)=max(50.,(rsoilso2(i,j)*3.41-85.)+
!                 !cmk amin1(max(0.,1000.*exp(269.-ts1)),1.e+5))
!                 rsoilso2(i,j)=max(50.,rsoilso2(i,j)*3.41-85.)+ &
!                 amin1(max(0.,1000.*exp(269.-ts1)),1.e+5)
!
!            ! -- for rh less than 40 %
!
!            if (rh2m(i,j).le.0.4) &
!                 !cmk    rsoilso2(i,j)=max(50.,amin1(1000.,(rsoilso2(i,j)*3.41-&
!                 !cmk    85.+((0.4-rh2m(i,j))/0.4)*1.e+5)+&
!                 rsoilso2(i,j)=max(50.,amin1(1000.,(rsoilso2(i,j)*3.41-85.+ &
!                 ((0.4-rh2m(i,j))/0.4)*1.e+3)+& 
!                 ! 1e3, see paper Laurens, JGR
!                 amin1(max(0.,1000.*exp(269.-ts1)),1.e+5)))
!
!            ! -- Temperature correction term for HNO3 above snow surface and ice
!            !    (Wesely, 1989), recomputation from K to 0C
!
!            rsnowhno3(i,j)=amin1(max(10.,1000.*exp(269.-ts1)),1.e+5)   
!            ! snow is melting changed (advise Wouter Greull) FD072003
!            if ( ts1 > (273.+3.) .and. sshf(i,j) > 10. ) rsnowhno3(i,j)=50. 
!            !
!            !
!         end do !nlon
!      end do !nlat
!
!
!      !cmk if ( okdebug ) 
!      !cmk   call dumpfield(0,'rsoilso4'//c_time//'.hdf',rsoilso4,'rsoilso4')
!      !cmk if ( okdebug ) 
!      !cmk   call dumpfield(0,'rsoilso2'//c_time//'.hdf',rsoilso2,'rsoilso2')
!      !cmk if ( okdebug ) 
!      !cmk   call dumpfield(0,'rsoilnh3'//c_time//'.hdf',rsoilnh3,'rsoilnh3')
!      !cmk if ( okdebug ) 
!      !cmk   call dumpfield(0,'rwatso4'//c_time//'.hdf',rwatso4,'rwatso4')
!      !cmk if ( okdebug ) 
!      !cmk   call dumpfield(0,'rsnowhno3'//c_time//'.hdf',rsnowhno3,'rsnowhno3')
!      !cmk if ( okdebug ) 
!      !cmk   call dumpfield(0,'rsnowso2'//c_time//'.hdf',rsnowso2,'rsnowso2')
!      !
!      !
!      vd11(:,:,:) = 1e-5
!      !CMK rsurf(:,:,:)=1e+5 
!      !
!      !  -- Loop for tracers, all timesteps
!      !
!      do jdep=1,ndep
!         jt = idep(jdep)
!         !
!         !  -- Only if a value for DIFFCF (diffusivity) is defined, the program 
!         !   calculates a surface resistance
!         !
!         if ( diffcf(jt) > 1e-5 ) then
!            !  
!            do j=1,nlat180
!               do i=1,nlon360
!
!                  xland=lsm(i,j)/100.
!                  ! laihelp: the max values derived from ECMWF are more realistic
!                  laihelp=min(lai(i,j),lai1(i,j))
!                  !   -- Assigning of values calculated in previous routine
!                  !   rsnow/rsoil of other components in data statement
!                  if (jt == iso2) then
!                     rsnow(jt)=rsnowso2(i,j)
!                     rsoil(jt)=rsoilso2(i,j)
!                  end if
!                  if (jt == inh3) rsoil(jt)=rsoilnh3(i,j)
!                  if (jt == iso4) then
!                     !   set snow resistance of so4 equal to soil resistance
!                     rsnow(jt)=rsoilso4(i,j) 
!                     rsoil(jt)=rsoilso4(i,j)
!                     rwat(jt)=rwatso4(i,j)
!                  end if
!                  !
!                  if (jt == ihno3) rsnow(jt)=rsnowhno3(i,j)
!
!                  !   -- Computation of stomatal resistance for component X, 
!                  !      correction based on differences in molecular
!                  !      diffusion of H2O and X and the water
!                  !      stress is also considered, FWS
!
!                  rstomx=diffcf(jt)*rstom(i,j)/fws(i,j)
!
!                  !   -- Calculation of mesophyll resistance of NO 
!                  !      and NO2 as function of
!                  !      the stomatal resistance of ozone
!
!                  if (jt == ino) rmes(jt)= 5.*diffcf(io3)*rstom(i,j)/fws(i,j)
!                  if (jt == ino2) rmes(jt)=0.5*diffcf(io3)*rstom(i,j)/fws(i,j)
!
!                  rleaf=(1./((1./rcut(jt))+(1./(rstomx+rmes(jt)))))
!                  ! linear upscaling from leaf scale to canopy scale 
!                  ! applying the LAI derived from Olson 
!
!                  rcanopy=rleaf/max(1e-5,laihelp)
!                  rveg=(1./((1./(rahcan(i,j)+rb(i,j)*diffrb(jt)+rsoil(jt)))+ &
!                       (1./rcanopy)))
!
!                  !-- It can be assumed that HNO3 deposition only depends on ra
!                  !   so that surface resistance = 0, however definition of 
!                  !   minimal surface resistance in order to avoid 
!                  !   extreme Vd values.
!
!                  if (jt == ihno3) rveg=1.
!
!                  !-- Computation of surface resistance Rs (s m-1) above land
!                  !   Vd for surface with vegetation and already incorporated 
!                  !   is the vegetation boudanry layer resistance
!
!                  rsveg=rveg+rb(i,j)*diffrb(jt)
!
!                  !-- Sulfate deposition calculation over vegetation according 
!                  !   to the observations by Wesely (1985), 
!                  !   see the Dry Deposition Inferential Model
!
!                  if (jt == iso4) then 
!                     !
!                     ! fd removed the factor stheta 
!                     ! (standard deviation of the wind direction)
!
!                     if ( ssr(i,j) > 250 ) then
!                        rsveg=1./(0.01*ustar_loc(i,j))
!                     else
!                        rsveg=1./(0.002*ustar_loc(i,j))
!                     end if
!                     !
!                     !   -- for particle deposition the rb is already 
!                     !      implicitly included
!                     !
!                     rssoil=rsoil(jt)
!                     !
!                     !   assumed that the wet skin fraction consists of vegetation
!                     !
!                     rsws=rsveg 
!                     rssnow=rsnow(jt)
!                     !
!                  else !jt neq iso4
!                     ! 
!                     !-- Rs for surface with bare soil
!                     !
!                     rssoil=rsoil(jt)+rb(i,j)*diffrb(jt)
!                     !
!                     !-- Wet skin fraction, 
!                     !   It is assumed that the wet skin fraction 
!                     !   mainly consists of vegetation thus laminar 
!                     !   resistance for vegetation may be applied
!                     !
!                     rsws=rws(jt)+rb(i,j)*diffrb(jt)
!                     !
!                     !-- Snow coverage
!                     !
!                     rssnow=rsnow(jt)+rb(i,j)*diffrb(jt)
!                     !
!                  end if
!                  !
!                  !-- land surface deposition velocity
!                  !
!
!                  cvsh=snow_cover(i,j) - &
!                       max(0.0,snow_cover(i,j)+vgrat_high(i,j)-1.0)
!                  ! fd same assumption as photolysis
!
!                  !cmk: need to correct the rsoil resistance here.
!                  !     fraction will become bare soil, which in the 
!                  !     vegetation deposition calculation
!                  !     will result in a too low resistance, 
!                  !     since the 'soil' part will be snow covered!
!
!                  cvwh=min(xland,wet_skin(i,j) )
!                  evgrath=min(xland,vgrat_low(i,j)+vgrat_high(i,j))
!
!                  ra1=raero(i,j)
!                  ! snow covered land fraction
!                  xland1=cvsh
!                  vdland1=1./(rssnow+ra1)
!                  ! wet skin covered land fraction (not covered by snow)
!                  xland2=max(0.,cvwh-xland1)
!                  vdland2=1./(rsws+ra1)
!                  ! vegetation covered land fraction 
!                  ! (not covered by snow and wet skin)
!                  xland3=max(0.,evgrath-xland1-xland2)
!                  vdland3=1./(rsveg+ra1)
!                  ! bare soil covered land fraction   
!                  ! (landfraction not covered by snow, wet skin, or vegetation)
!                  xland4=max(0.,xland-xland1-xland2-xland3)
!                  vdland4=1./(rssoil+ra1)  
!
!                  xwat1=min(1-xland,ci(i,j))
!                  vdwat1=1./(rsnow(jt)+rb(i,j)*diffrb(jt)+ra1)
!                  xwat2 = max(0.,1.-xland-ci(i,j))
!                  vdwat2=1./(rwat(jt)+rb(i,j)*diffrb(jt)+ra1)
!                  if (jt == iso4) then
!                     vdwat1=1./(rsnow(jt)+ra1)
!                     vdwat2=1./(rwat(jt)+ra1)
!                  end if
!                  xsum=xland1+xland2+xland3+xland4+xwat1+xwat2
!                  !fd    if (xsum.ne.1.) then
!                  !fd    print *,'sum',i,j,xsum,xland,xland1,xland2,'xl3',&
!                  !fd    xland3,xland4,xwat1,xwat2,cvsh,cvwh,evgrath,ci(i,j)
!                  !fd    end if
!                  !
!                  !-- Computation of deposition velocity 
!                  !
!                  vd11(i,j,jdep)=xland1*vdland1+xland2*vdland2+xland3*vdland3+&
!                       xland4*vdland4+xwat1*vdwat1+xwat2*vdwat2
!               end do !End of loop longitude
!            end do  !End of loop latitude
!
!            adum=names(jt)
!            i=len_trim(adum)
!
!            !cmk   if ( okdebug ) call dumpfield(0,'vd'//&
!            !cmk      adum(1:i)//c_time//'.hdf',vd11(:,:,jdep),'vd'//adum(1:i))
!            if ( okdebug ) &
!                 call dd_field_statistics('vd'//adum(1:i),vd11(:,:,jdep), &
!                 nlon360,nlat180)
!
!         end if !  (diffcf(jt) > 1e-5) 
!         !
!      end do   !End of loop tracer
!      deallocate(rwatso4)
!      deallocate(rsoilnh3)
!      deallocate(rsoilso2)
!      deallocate(rsoilso4)
!      deallocate(rsnowhno3)
!      deallocate(rsnowso2)
!      ! 
!      if ( okdebug ) write(*,*) 'dd_calc_rs: end'
!    end subroutine dd_calc_rs
!    !
!    !
!    subroutine dd_coarsen_vd 
!      !
!      ! 
!      ! 
!      use toolbox, only : coarsen_emission
!      implicit none
!      integer i,j,jt,jdep,idno2,idh2o2,idhno4,region
!      integer,parameter  :: avg_field = 1
!      character(len=8) :: adum
!      character(len=2) :: bdum
!      real, dimension(:,:),allocatable :: vdhelp
!
!      allocate(vdhelp(nlon360,nlat180))
!      !
!      ! --
!      ! Scale the surface resistance of a number of trace gases and aerosols
!      ! NOx deposition for all NOx components separately and scaling 
!      ! of NOx afterwards
!      ! 
!      do jdep = 1,ndep   
!         select case(idep(jdep))
!         case(ihno4)
!            idhno4 = jdep
!         case(ih2o2)
!            idh2o2 = jdep
!         case(ino2)
!            idno2 = jdep
!         case default
!         end select
!      end do
!
!      vd11(:,:,idhno4)=(vd11(:,:,idno2)+vd11(:,:,idh2o2))/2.
!      do jdep=1,ndep
!         jt = idep(jdep)
!         do j=1,nlat180
!            do i=1,nlon360
!               vdhelp(i,j)=vd11(i,j,jdep)
!            end do
!         end do
!         call coarsen_emission('vd',nlon360,nlat180,vdhelp,vd_temp,avg_field)
!         do region = 1,nregions
!            vd(region,jt)%surf(:,:) = vd_temp(region)%surf(:,:)
!
!            if ( okdebug ) then 
!               adum=names(jt)
!               write(bdum,'(i2.2)') region
!               i=len_trim(adum)
!               !cmk  call dumpfield(0,'vd.hdf',&
!               !cmk     vd(region,jt)%surf(:,:),'vd_'//adum(1:i)//bdum)
!            end if
!         end do
!
!
!      end do !jdep
!      deallocate(vdhelp)
!
!      do region = 1,nregions
!         vd(region,irooh  )%surf(:,:) = vd(region,ich3o2h)%surf(:,:)
!         vd(region,iorgntr)%surf(:,:) = vd(region,ipan   )%surf(:,:)
!         vd(region,in2o5  )%surf(:,:) = vd(region,ihno3  )%surf(:,:)
!         vd(region,ino3   )%surf(:,:) = vd(region,ino2   )%surf(:,:)
!         vd(region,inh4   )%surf(:,:) = vd(region,iso4   )%surf(:,:)
!         vd(region,imsa   )%surf(:,:) = vd(region,iso4   )%surf(:,:)
!         vd(region,ino3_a )%surf(:,:) = vd(region,iso4   )%surf(:,:)
!         vd(region,io3s   )%surf(:,:) = vd(region,io3    )%surf(:,:)
!      end do
!
!      ! coarsen ags for photolysis...(temp in vd_temp)
!      call coarsen_emission('ags',nlon360,nlat180,ags,vd_temp,avg_field)
!      do region = 1,nregions
!         phot_dat(region)%albedo(:,:) = vd_temp(region)%surf(:,:)
!         phot_dat(region)%ags_av = phot_dat(region)%ags_av + &
!              phot_dat(region)%albedo   !statistics...
!         phot_dat(region)%nalb_av = phot_dat(region)%nalb_av + 1
!      end do
!
!    end subroutine dd_coarsen_vd
!    !
!  end subroutine dd_surface_fields
!
!
!
!  subroutine dd_field_statistics(name,field,ix,jx)
!    !-----------------------------------------------
!    !
!    ! Purpose:
!    ! -------
!    ! this subroutine calculate the min,mean,max of a real field
!    !
!    ! External
!    ! --------
!    ! none
!    !
!    ! Reference
!    ! ---------
!    ! none
!    !------------------------------------
!    !
!    implicit none
!    character(len=*),intent(in) :: name
!    integer,intent(in) :: ix,jx
!    real,dimension(ix,jx),intent(in) :: field
!    integer :: i,j,ntel_non_zero
!    real :: maxf,minf,meanf,mean_non_zero
!    maxf=-1e20
!    minf=1e12
!    meanf=0.
!    mean_non_zero=0.
!    ntel_non_zero=0
!    do i=1,ix
!       do j=1,jx
!          meanf=meanf+field(i,j)
!          maxf=max(maxf,field(i,j))
!          minf=min(minf,field(i,j))
!          if (field(i,j).ne.0) then
!             ntel_non_zero=ntel_non_zero+1
!             mean_non_zero=mean_non_zero+field(i,j)
!          end if
!       end do
!    end do
!    meanf=meanf/ix/jx
!    if (ntel_non_zero > 0) mean_non_zero= mean_non_zero/ntel_non_zero
!    write(6,'(a10,4(a3,1x,1pe10.3))') &
!         name,'max',maxf,'min',minf,'mean',meanf,'mn0',mean_non_zero
!
!  end subroutine dd_field_statistics
!
!
!
!  subroutine dd_field_statistics_i8(name,field,ix,jx)
!    !------------------------------------
!    !
!    ! Purpose:
!    ! -------
!    ! this subroutine calculate the min,mean,max of a field of int8
!    !
!    ! External
!    ! --------
!    ! none
!    !
!    ! Reference
!    ! ---------
!    ! none
!    !------------------------------------
!
!    implicit none
!    character(len=*),intent(in) :: name
!    integer,intent(in) :: ix,jx
!    integer(kind=1),dimension(ix,jx),intent(in) :: field
!    integer :: i,j,ntel_non_zero
!    real :: maxf,minf,meanf,mean_non_zero
!    maxf=-1e20
!    minf=1e12
!    meanf=0.
!    mean_non_zero=0.
!    ntel_non_zero=0
!    do i=1,ix
!       do j=1,jx
!          meanf=meanf+field(i,j)
!          maxf=max(maxf,real(field(i,j)))
!          minf=min(minf,real(field(i,j)))
!          if (field(i,j).ne.0) then
!             ntel_non_zero=ntel_non_zero+1
!             mean_non_zero=mean_non_zero+field(i,j)
!          end if
!       end do
!    end do
!    meanf=meanf/ix/jx
!    if (ntel_non_zero > 0) mean_non_zero= mean_non_zero/ntel_non_zero
!    write(6,'(a10,4(a3,1x,1pe10.3))') &
!         name,'max',maxf,'min',minf,'mean',meanf,'mn0',mean_non_zero
!
!  end subroutine dd_field_statistics_i8
!
!
end module dry_deposition
