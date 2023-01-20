!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#define CHECK_NCSTAT(nc_ret_code) if (nc_ret_code /= NF90_NOERR) then; status=1; IF_NOTOK_RETURN(status=1); end if
!
#include "tm5.inc"
!
!###############################################################################

module diffusion

  use GO, only : gol, goPr, goErr

  use dims, only                  : nregions, revert
  use binas, only                 : xmair
  use chem_param, only            : emis_data
  use Go, only                    : TDate
  use os_specs, only              : MAX_FILENAME_LEN

  implicit none

  ! --- in/out ------------------------------------

  private
  public   ::  Diffusion_Init, Diffusion_Done
  public   ::  Diffusion_Setup

  ! --- custom types -----------------------------
  type T_diffcoeff
    real, allocatable           :: dkg(:,:,:,:)
    real, allocatable           :: blh(:,:,:)
    integer                     :: nsteps
    type(TDate), allocatable    :: t1(:), t2(:)
    integer                     :: cur_date(3)
  end type T_diffcoeff


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'diffusion'


  ! --- local ------------------------------------
  type(T_diffcoeff), allocatable    :: dkg_ncfile(:)

  ! diffusion files:
  character(len=MAX_FILENAME_LEN) :: diffusion_dir
  character(len=20)             :: mlevs
  logical, allocatable          :: file_complete(:)
  integer                       :: deflate_lvl ! compression level
  logical                       :: write_kvh   ! whether to store the diffusion coefficients or not

  type(emis_data), target       :: ss(nregions)
  type(emis_data), target       :: ustar_loc(nregions)
  type(emis_data), target       :: sr_mix(nregions)

  !--------------------------------------------------------------------------------
  ! global input fields:
  ! phlb(:,:,lm+1)          ! half level pressure [Pa]
  ! gph(:,:,lm+1)           ! height of half-levels [m]    !CMK note gph starts NOW at 1!
  ! t                       ! temperature [K]
  ! q                       ! specific humidity [kg/kg]
  ! uwind                   ! mean wind W-E [m/s]
  ! vwind                   ! mean wind S-N [m/s]
  ! sshf                    ! surface sensible heat flux [W / m^2]
  ! slhf                    ! surface latent heat flux [W / m^2]
  ! ustar (ustar_loc)       ! velocity scale
  !
  ! global output fields:
  ! kvh                     ! eddy diffusivity for heat at top [m^2/s]
  ! pblh                    ! boundary layer height[m]
  ! dkg                                     ! vertical diffusion coefficient [ks s-1]
  ! Adjusted for TM5 by Bram Bregman
  ! April 2005/ adjoint version
  !--------------------------------------------------------------------------------


contains


  ! =====================================================


  subroutine Diffusion_Init( status )

    use GO                 , only : ReadRc
    use global_data        , only : rcF
    use dims               , only : nregions, im, jm, lmax_conv
    use MeteoData          , only : lsmask_dat
    use MeteoData          , only : sr_ols_dat, sr_ecm_dat
    use MeteoData          , only : u10m_dat, v10m_dat
    use MeteoData          , only : slhf_dat, sshf_dat
    use MeteoData          , only : ewss_dat, nsss_dat
    use Meteo              , only : Set

    implicit none

    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Diffusion_Init'

    ! --- local ---------------------------------

    integer           ::  n

    ! --- begin ---------------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! directory for temporary diffusion files:
    call ReadRc(rcF, 'diffusion.dir', diffusion_dir, status)
    IF_NOTOK_RETURN(status=1)
    ! should we store the diffusion coefficients?
    call ReadRc(rcF, 'diffusion.coefficients.store', write_kvh, status, default=.false.)
    IF_ERROR_RETURN(status=1)
    ! read the string 'tropo25' or its equivalent somewhere
    call ReadRc(rcF, 'my.levs', mlevs, status)
    IF_NOTOK_RETURN(status=1)
    ! Should we compress the netcdf dkg files?
    call ReadRc(rcF, 'diffusion.files.deflate.level', deflate_lvl, status, default=1)
    IF_ERROR_RETURN(status=1)

    allocate(dkg_ncfile(nregions))
    do n = 1, nregions
        dkg_ncfile(n)%cur_date(:) = -9999

    end do
    ! assume that the first day's dkg file is complete, possibly to be proved wrong later
    allocate(file_complete(nregions))
    file_complete = .true.

    !
    ! enable meteo
    !

    ! only forward model, adjoint should read pre-computed files:
    if ( revert > 0 ) then

      ! loop over zoom regions:
      do n = 1, nregions

        ! enable meteo:
        call Set( lsmask_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)

        ! enable meteo:
        call Set( sr_ecm_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)
        ! enable meteo:
        call Set( sr_ols_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)

        ! enable meteo:
        call Set( u10m_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)
        ! enable meteo:
        call Set( v10m_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)

        ! enable meteo:
        call Set( slhf_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)
        ! enable meteo:
        call Set( sshf_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)

        ! enable meteo:
        call Set( ewss_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)
        ! enable meteo:
        call Set( nsss_dat(n), status, used=.true. )
        IF_NOTOK_RETURN(status=1)

      end do  ! regions
    end if  ! forward run

    !
    ! done
    !
    write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Diffusion_Init


  ! ***


  subroutine Diffusion_Done( status )

    implicit none
    ! --- in/out --------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Diffusion_Done'

    ! --- begin ---------------------------------

    deallocate(dkg_ncfile)
    deallocate(file_complete)

    ! ok
    status = 0

  end subroutine Diffusion_Done


  ! ***


  ! called every timestep, setup OH fields etc

  subroutine Diffusion_Setup( t1, t2, status )

    use GO,     only : TDate
    use dims,   only : nregions

    implicit none

    ! --- in/out --------------------------------

    type(TDate), intent(in)   ::  t1, t2
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Diffusion_Setup'

    integer     :: region

    ! --- begin ---------------------------------
    !write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! ...
    ! copied from 'sources_sinks/trace_after_read
    ! ...
    do region = 1, nregions

        call read_diffusion(region, t1, t2, status)

        if ( status > 0 ) then
            ! error reading
            TRACEBACK; status=1; return
        else if ( status < 0 ) then
            ! dkg file could not be found --> calculate
            call calc_kzz(region, status)
            IF_NOTOK_RETURN(status=1)
            ! write to file
            call write_diffusion(region, t1, t2, status)
            IF_NOTOK_RETURN(status=1)
        end if

    end do ! region

    !write (gol,'(a," : done")') trim(rname) ; call goPr

    ! ok
    status = 0

  end subroutine Diffusion_Setup


  subroutine calc_kzz(region, status)
    !------------------------------------
    !
    ! Purpose:
    ! -------
    ! this subroutine reads and prepares all fields needed for the calculation of kzz
    !
    ! External
    ! --------
    ! dd_get_3_hourly_surface_fields
    ! dd_calc_ustar_raero_rb
    ! dd_coarsen_fields
    !
    ! Reference
    ! ---------
    ! Ganzeveld and Lelieveld (1996) and references therein.
    ! Adjusted for TM5, Bram Bregman, August 2003
    !
    !

    use binas,       only: vkarman, grav
    use dims,        only: im, jm, lmax_conv, idate, lm, okdebug
    use dims,        only: revert
    use dims,        only: ndyn, tref, at, bt
    use dims,        only: isr,jsr,ier,jer
    use global_data, only : mass_dat, conv_dat, wind_dat, region_dat
    use MeteoData  , only : m_dat, phlb_dat
    use MeteoData  , only : lsmask_dat
!    use MeteoData  , only : gph_dat, humid_dat, temper_dat, pu_dat, pv_dat
!    use MeteoData  , only : sr_ols_dat, sr_ecm_dat
!    use Meteo      , only : u10m_dat, v10m_dat
!    use Meteo      , only : slhf_dat, sshf_dat
!    use Meteo      , only : ewss_dat, nsss_dat
#ifdef MPI
       use mpi_const,only: myid, root_k, ierr, com_lev, my_real, lmloc
#endif

    implicit none

    ! --- in/out ----------------------------

    integer, intent(in)     :: region
    integer, intent(out)    :: status

    ! --- const -------------------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_kzz'

    ! --- local ------------------------------

    character(len=10)                            :: c_time
    real,dimension(:,:,:),pointer                :: phlb, m
!    real, dimension(:,:,:), pointer              :: pu, pv, t, q, gph, dkg
    real,dimension(:,:,:),pointer                :: am,bm,cm ! fluxes in kg/timestep(region)
    real,dimension(:),pointer                    :: dxyp
    integer,parameter                            :: avg_field=1,nlon360=360,nlat180=180
    integer                                      :: i, nsend, ntimes
    real, dimension(:,:,:),allocatable,target    :: m_adj, phlb_adj
    real, dimension(:,:),allocatable             :: p_adj
    integer                                      :: imr, jmr, lmr, j,l

    ! --- begin ----------------------------------

    if (okdebug) write(*,*) 'start of calc_kzz'

    ! safety check ...
    if ( revert < 0 ) then
      write (gol,'("reverse model should read Kzz, implementation not trusted yet ...")'); call goErr
      TRACEBACK; status=1; return
    end if

#ifdef MPI
    if(myid == root_k) then    ! reading /calculations by root!
#endif

    allocate(ss(region)%surf(im(region),jm(region)))
    allocate(ustar_loc(region)%surf(im(region),jm(region)))
    allocate(sr_mix(region)%surf(im(region),jm(region)))

    ! -- for output of time header on debug fields

    write(c_time,'(i4,3(i2))') idate(1:4)   !CMK corrected

    do i=1,10
    if (c_time(i:i)==' ') c_time(i:i)='0'
    end do

    !
    ! surface data sets may have the following characteristics
    ! 1. instanteneous. The time written in the file is valid from t-1.5 to t+1.5
    ! 2. Accumulated. The time written in the file is valid from t to t+3
    ! 3. Daily averaged. The time on the file is always 00 hours, and valid until t+24
    !
    ! ***   It is essential to understand this timing when reading and applying these sets.
    !
    !
    ! fd2mk it has to be decided to 'save' fields like vgrat_low and daily average surface fields
    ! for later use, or to read it every 3 hours time step

    ! wind10m,slhf,sshf
    call dd_get_3_hourly_surface_fields(region, status)
    IF_NOTOK_RETURN(status=1)

    if (revert == -1) then
       !mkadj
       !mass and pressure are at t=t+dt in adjoint
       !these should be brought back to t to get
       !exactly the same diffusion as in forward run
       imr = im(region) ; jmr = jm(region) ; lmr = lm(region)
       allocate ( m_adj(1:imr,1:jmr,1:lmr) )
       allocate ( phlb_adj(imr,jmr,lmr+1))
       allocate ( p_adj(imr,jmr))
       am => wind_dat(region)%am_t
       bm => wind_dat(region)%bm_t
       cm => wind_dat(region)%cm_t
       dxyp => region_dat(region)%dxyp

       ! this point should not been reached ...
       write (gol,'("expecting adjoint run to read kzz from files saved in apri run")'); call goErr
       TRACEBACK; status=1; return
       ! no 'nreadX' anymore, so stop here ...
       !ntimes = (nreadX/(ndyn))*tref(region)

       m_adj = m_dat(region)%data(1:imr,1:jmr,1:lmr)

       do i=1,ntimes
          m_adj = m_adj + revert*am(0:imr-1,1:jmr  ,1:lmr  ) &   !east
                        - revert*am(1:imr  ,1:jmr  ,1:lmr  ) &   !west
                        + revert*bm(1:imr  ,1:jmr  ,1:lmr  ) &   !south
                        - revert*bm(1:imr  ,2:jmr+1,1:lmr  ) &   !north
                        + revert*cm(1:imr  ,1:jmr  ,0:lmr-1) &   !lower
                        - revert*cm(1:imr  ,1:jmr  ,1:lmr  )     !upper
       enddo
       p_adj(:,:) = 0.0
       do l=1,lmr
          ! note that on the EDGES masses are coarse>
          ! diffusion iis applied only on core zoom>
          do j=jsr(region), jer(region)
             do i=isr(region), ier(region)
                p_adj(i,j) = p_adj(i,j) + m_adj(i,j,l)*grav/dxyp(j)
             end do
          end do
       end do
       do l=1,lmr+1
          do j=jsr(region), jer(region)
             do i=isr(region), ier(region)
                phlb_adj(i,j,l) = at(l)+bt(l)*p_adj(i,j)
             end do
          end do
       end do
       deallocate(p_adj)
       nullify(am)
       nullify(bm)
       nullify(cm)
       nullify(dxyp)

       phlb => phlb_adj
       m => m_adj
    else
       phlb => phlb_dat(region)%data
       m    =>    m_dat(region)%data
    end if  ! revert is 1 or -1

!    pu => pu_dat(region)%data
!    pv => pv_dat(region)%data
!    t   => temper_dat(region)%data
!    q   =>  humid_dat(region)%data
!    gph =>    gph_dat(region)%data
!    dkg => conv_dat(region)%dkg

    ! ustar_loc
    call dd_calc_ustar(region)

    ! calculate kzz
    call bldiff(region, phlb, m)

    nullify(phlb)
    nullify(m)
!    nullify(pu)
!    nullify(pv)
!    nullify(dkg)
!    nullify(t)
!    nullify(q)
!    nullify(gph)
    if(revert == -1) then
       deallocate(m_adj)
       deallocate(phlb_adj)
    endif

#ifdef MPI
    end if   !myid == root_k
#endif

#ifdef MPI
    nsend = im(region)*jm(region)*lmax_conv
    if(lmloc.ne.0) call mpi_bcast(conv_dat(region)%dkg ,nsend, my_real, root_k, com_lev, ierr)

    if ( myid == root_k ) then
#endif
      if (okdebug) print *, 'vertical diffusion kzz (dkg) broadcasted'

    deallocate(ss(region)%surf)
    deallocate(ustar_loc(region)%surf)
    deallocate(sr_mix(region)%surf)

#ifdef MPI
    endif
#endif

    if (okdebug) write(*,*) 'end of calc_kzz'

    ! ok
    status = 0

    !------------------------------------------------------------------
    ! Procedure:
    ! (1) calc_Kv: Diffusion coefficients full atmosphere from shear
    ! (2) pblhght: Determine the boundary layer height
    ! (3) difcoef: Replace diffusion coefficients in the boundary layer
    !------------------------------------------------------------------
    ! subroutine adjusted for TM5, Bram Bregman, August 2003.
    !-------------------------IO---------------------------------------


  end subroutine calc_kzz

  subroutine dd_get_3_hourly_surface_fields(region, status)

    !----------------------------------------
    !
    ! Purpose
    ! -------
    ! read 3-hourly ECMWF datasets, can be instanteneous or acccumulated
    !
    ! External
    ! --------
    ! date2tau
    ! tau2date
    ! read_meteo
    ! dd_field_statistics and dd_field_statistics_i8
    !
    ! Reference
    ! ----------
    ! None
    !------------------------

      use datetime,   only: date2tau, tau2date
      use toolbox,    only: escape_tm, coarsen_emission
      use Dims      , only : nregions, idate, okdebug
      use MeteoData , only : ewss_dat, nsss_dat

      implicit none

      ! --- in/out -----------------------------------------

      integer, intent(in)       :: region
      integer, intent(out)      :: status

      ! --- const -------------------------------------

      character(len=*), parameter  ::  rname = mname//'/dd_get_3_hourly_surface_fields'

      integer,parameter                    :: avg_field=1,nlon360=360,nlat180=180

      ! --- local --------------------------------------------

      integer,dimension(6)                 :: idater_acc,idater
      integer                              :: itaux ! auxiliary variable
      real,dimension(:,:), allocatable     :: field
      real,dimension(:,:), allocatable     :: field2

      ! --- begin ---------------------------------------------

      allocate(field(nlon360,nlat180))
      allocate(field2(nlon360,nlat180))

      !
      !     [t   ,        t+3]
      !
      !    idate
      !    idater_acc     idater
      !

      call date2tau(idate,itaux)
      call tau2date(itaux,idater_acc)
      if (revert == 1) then
         call tau2date(itaux+3*3600,idater) !CMK use 00 for 21-00, 03 for 00-03 etc... (more handy!)
         call tau2date(itaux,idater_acc)
      else
         call tau2date(itaux,idater)
         call tau2date(itaux-3*3600,idater_acc)        !CMKADJ
      endif

      if (okdebug) write(*,*) 'dd_get_3_hourly_surface_fields read time on instanteneous field',idater
      if (okdebug) write(*,*) 'dd_get_3_hourly_surface_fields read time on accumulated field',idater_acc

      !
      ! read accumulated fields
      !

      ss(region)%surf = sqrt( ewss_dat(region)%data(:,:,1)**2 + nsss_dat(region)%data(:,:,1)**2 )
      if (okdebug) write(*,*) 'end of dd_get_3_hourly_surface_fields'

!PBi
!     missing deallocate caused major memory leak !!!!!
!     360 x 180 x 8 Byte (real(8)) x 8 (1/day) = 4 MByte / day = 120 MByte / month
      deallocate(field)
      deallocate(field2)
!PBf
      ! ok
      status = 0

  end subroutine dd_get_3_hourly_surface_fields

  subroutine dd_calc_ustar(region)

    !--------------------------------
    !
    ! Purpose
    ! -------
    ! Calculate friction velocity (ustar)
    ! Method
    ! ------
    ! uses the log normal wind profile information stored by ECMWF in 10 meter wind
    ! to estimate a local ustar over land
    ! uses Charnock equation over sea to estimate ustar
    ! aerodynamic resistance is calculated from heat fluxes and ustar using a constant reference height
    ! sub laminar resistance from ustar
    !
    ! External
    ! --------
    ! dd_field_statistics
    ! dumpfield
    !
    ! Reference
    ! ----------
    ! Ge Verver (personal communication, 2003)
    !------------------------

      use binas,     only : grav, vKarman
      use dims,      only : im, jm, okdebug
      use MeteoData, only : lsmask_dat
      use MeteoData, only : sr_ecm_dat, sr_ols_dat
      use MeteoData, only : u10m_dat, v10m_dat

      implicit none

      !--- in/out ------------------------------------

      integer,intent(in)              :: region

      ! --- const -------------------------------------

      ! Garret relation:
      real,parameter  :: alfa_garret    = 0.11
      ! air density at seal level and 293 K :
      real, parameter :: rho_air = 1.2  ! kg/m3

      ! Charnock relation:
      real,parameter  :: alfa_charnock  = 0.018
      ! kinematic viscosity of air at about 300 K :
      real,parameter  :: nu_air         = 1.5e-5    ! m2/s

      real, parameter :: Href=30.     ! constant reference height for calculations of raero
                                      ! some constants specific to calculation of ustar and raero
      real,parameter  :: rhoCp          = 1231.0
      real,parameter  :: rhoLv          = 3013000.0

      !real, parameter :: rz0 =2.0

      !real,parameter  ::alfa_charnock1=0.11,alfa_charnock2=0.018,v_charnock=1.5e-5 !(m2/s)

      !--- local -------------------------------------

      integer ::i,j

      real :: ustar_sea,ustar_land,xland,sr_sea,sr_land,sr_help,wind10m
      real,dimension(:,:),allocatable :: field

      ! --- begin ----------------------------------

      allocate(field(im(region),jm(region)))

      do j=1,jm(region)
        do i=1,im(region)
          xland = lsmask_dat(region)%data(i,j,1) / 100.0

          !
          ! From M.Z. Jacobson, "Fundamentals of Atmospheric Modeling",
          !   section "8.3 Friction velocity" :
          ! The friction velocity :
          !   u* = sqrt( |tau_z|/rho_a )   [m/s]
          ! where:
          !   tau_z = (ewss,nsss)       surface stress [N/m2=kg/m/s2]
          !   rho_air                   air density [kg/m3] ; about 1.2  at sea level and 293 K
          !
          ustar_sea = sqrt( ss(region)%surf(i,j) / rho_air )  ! m/s
          !
          ! limitation to avoid division by zero:
          ustar_sea = max( 0.001, ustar_sea )    ! m/s ;  minium of 1 mm/s
          !
          ! From M.Z. Jacobson, "Fundamentals of Atmospheric Modeling",
          !   section "8.4 Surface roughness lengths" :
          ! "For smooth surfaces, such as over a smooth ocean with low wind speeds (Garret):
          !   z0m = 0.11 nu_a / u*
          ! where nu_a is the kinematic viscosity of air.
          ! Over a rough ocean with high wind speeds:
          !   z0m = alpha_c (u*)**2 / grav
          ! which is the 'Charnock relation',
          ! where alpha_c ~ 0.016 is the 'Charnock constant'."
          !
          ! Here the sum of these two parameterizations is used,
          ! where the first part is dominant for u*<0.1 and the second for u*>0.1 ;
          ! minimum value for u* prevents division by zero:
          !
          ! surface roughness at sea:
          !    m                 m^2/s    m/s                         (m/s)^2        m/s^2
          sr_sea = alfa_garret * nu_air / ustar_sea + alfa_charnock * ustar_sea**2 / grav

          !
          ! LAND
          ! calculate the 'local' surface roughness for momentum that is consistent with 10 m wind
          ! and the Olsson data base, we assume that the windspeed at 75 m is independent of
          ! surface roughness
          !
          wind10m = sqrt( u10m_dat(region)%data(i,j,1)**2 + &
                          v10m_dat(region)%data(i,j,1)**2     )

          if (xland>0.) then
            !>>> TvN
            ! Over land, the friction velocity ustar
            ! is calculated from the 10 m wind speed, u10.
            ! In IFS u10 is a diagnostic variable, calculated to be
            ! compatible with "open-terrain" wind speed observations.
            ! Over rough or inhomogeneous terrain (z0M > 0.03 m),
            ! it is calculated using an aerodynamic roughness
            ! length typical for open terrain with grassland (0.03 m)
            ! (see IFS cy31r1 documentation, p. 46).

            ! In the expressions applied in TM5,
            ! the stability profile functions Psi_M have disappeared,
            ! and only the logarithmic terms are kept.
            ! Apart from this, the expression for ustar was incorrect:
            ! 10./sr_land should be 75./sr_land.
            ! This has now been corrected.

            ! Moreover, over islands and near coast lines,
            ! zeros can occur in sr_ols_dat,
            ! which have to be removed.
            ! However, a lower bound of 1e-2 m = 1 cm
            ! seems too high for smooth surfaces,
            ! like deserts and ice caps.
            ! The minimum positive value in sr_ols_dat
            ! is 2.5e-4 m = 0.025 cm,
            ! which is obtained in parts of Antarctica.
            ! This is likely the result of regridding of
            ! a field with minimum value of 0.1 cm
            ! from 0.5x0.5 to 1x1 degrees.

            ! Please note that with the introduction of CY31R1,
            ! the orographic contribution to the aerodynamic roughness length
            ! in IFS has been replaced by an explicit specification of stress
            ! on model levels due to turbulent orographic form drag,
            ! and the climatological variable previously used
            ! has been replaced by a prognostic variable.
            ! In TM5, the climatological variable is still used,
            ! but it would be better to use the prognostic variable instead.

            ! Note also that the same calculation
            ! is done in dry_deposition.F90.

            ! occurs at Islands, etc.
            sr_land = max( 1e-3, sr_ols_dat(region)%data(i,j,1)       )
            sr_help = min(       sr_ecm_dat(region)%data(i,j,1), 0.03 )

            ! ustar consistent with 'clipped' large scale roughness
            !fd ustar_land=vKarman*wind10m(i,j)/alog(10./sr_help)

            !ustar_land=vKarman*wind10m/&
            !                        alog(10./sr_help)*alog(75./sr_help)/alog(10./sr_land)
            ustar_land=vKarman*wind10m/alog(10./sr_help)*alog(75./sr_help)/alog(75./sr_land)
            ! <<< TvN
          else
            sr_land = 0.0
            ustar_land = 0.0
          endif

                  ustar_loc(region)%surf(i,j)=xland*ustar_land+(1-xland)*ustar_sea
                  sr_mix(region)%surf(i,j) =xland*sr_land  +(1-xland)*sr_sea
        enddo !i
      enddo !j
      field = ustar_loc(region)%surf

      if (okdebug) call dd_field_statistics('ustar_loc',field,im(region),jm(region))
      !if (okdebug) call dumpfield(0,'ustar_loc'//c_time//'.hdf',field,'ustar_loc')

      deallocate(field)

      if (okdebug) write(*,*) 'end calc_ustar'

  end subroutine dd_calc_ustar

  subroutine bldiff(region, phlb, m)

      use binas,        only : ae, cp_air, Rgas, grav, Lvap, vkarman
      use dims,         only : gtor, dx, dy, ybeg, xref, yref, okdebug
      use dims,         only : isr, ier, jsr, jer, idate, im, jm, lm, lmax_conv
      use toolbox,      only : dumpfield
      use MeteoData,    only : slhf_dat, sshf_dat, temper_dat, humid_dat, gph_dat, pu_dat, pv_dat
      use global_data,  only : conv_dat

      implicit none

      ! --- in/out ----------------------------------------------

      integer, intent(in)           :: region
      real, pointer, intent(in)     :: phlb(:,:,:), m(:,:,:)

      ! --- const ----------------------------------

      ! minimum value diffusion coefficient
      real, parameter          :: ckmin = 1.e-15

      real, parameter          :: sffrac=  0.1
      real, parameter          :: onet   = 1./3.
      real, parameter          :: rair   = Rgas*1.e3/xmair
      real, parameter          :: ccpq   = 1869.46/cp_air -1.
      real, parameter          :: epsilo = rair/461.51
      real, parameter          :: apzero = 100000. ! reference pressure (usually 1000 hPa)
      real, parameter          :: ccon   = sffrac*vkarman
      real, parameter          :: betam=15.,betah=15.
      real, parameter          :: binm   = betam*sffrac
      real, parameter          :: binh   = betah*sffrac
      real, parameter          :: fr_excess = 1.
      real, parameter          :: fak = fr_excess * 8.5
      real, parameter          :: fakn = fr_excess * 7.2

      real, parameter          :: zkappa = rair / cp_air
      real, parameter          :: zcrdq = 1.0/epsilo - 1.0

      real, parameter          :: tiny     = 1.e-9   ! bound wind**2 to avoid dividing by 0
      real, parameter          :: zacb     = 100.    !factor in ustr-shearproduction term
      real, parameter          :: ricr     = 0.3     ! critical richardson number

      real, parameter          :: ssfrac=0.1

      ! --- local ------------------------------------

      real, dimension(:,:,:), allocatable:: wkvf,&     ! outer layer Kv (at the top of each layer)
                                            zdup2,&    ! vertical shear squared
                                            zrinub,&   ! Richardson number
                                            wthm,&     ! potential temperature
                                            zthvk,&    ! virtual potential temperature
                                            pf,&       ! full level pressure
                                            whm, &     ! height of layer centers
                                            uwind,vwind,& ! wind velocities [m s-1]
                                            kvh                     !

      real, dimension(:,:), allocatable  :: ustr,&     ! velocity scale
                                            wheat,&    ! surface sensible heat flux [K m/s]
                                            wqflx,&    ! surface latent heat flux [kg m/ (kg s)]
                                            wobkl,&    ! Monin Obukhov length
                                            wheatv,&   ! surface virtual heat flux [K m/s]
                                            wtseff   ! theta_v at lowest level [K]

      integer             :: i, j, l
      real                :: intfac    !cmk BUG was integer!

      ! variables for the calculation of free atmosphere vertical diffusion coefficients wkvf
      real                        :: cml2,zlambdac,arg,zdz,zthva,zfunst, zfstabh,zsstab,zkvn,dyy
      real,dimension(jm(region))  :: dxx,lat
      integer                     :: jq! jqif function

      ! variabales for the calculation of the boundary layer height
      real                        :: wrino      ! Richardson number
      real                        :: thvparc    ! parcel virtual potential temperature
      real                        :: fmt,wsc,tkv,zrino,vvk,dtkv,zexcess
      integer                     :: jwork
      real,dimension(:,:),pointer   :: pblh         ! planetary boundary layer height [m]

      ! variables to calculate difusion coefficient in the boundary layer
      !real                        :: wsc,cml2
      real                        :: z,zwstr,zkvh,kvhmin,zm,zp,zslask,zsl1, &
                                     zh,zzh,zl,pr,zpblk,zfstabh1,zfstabh2,&
                                     term1,term2,term3,obukov
      integer                     :: jq1,jq2,jq3,jq4,jq5,jck,jqq, jq6, jq7

      real                      :: yres
      real                      :: xres

      real                      :: thvgrad, kvhentr
      real, dimension(:,:,:), pointer   :: t, q, gph, pu, pv, dkg


      ! --- begin ------------------------------------------

      if (okdebug) write(*,*) 'start of bldiff'

      allocate(wkvf(im(region),jm(region),lm(region)))
      allocate(zdup2(im(region),jm(region),lm(region)))
      allocate(zrinub(im(region),jm(region),lm(region)))
      allocate(wthm(im(region),jm(region),lm(region)))
      allocate(zthvk(im(region),jm(region),lm(region)))
      allocate(pf(im(region),jm(region),lm(region)))
      allocate(whm(im(region),jm(region),lm(region)))
      allocate(uwind(im(region),jm(region),lm(region)))
      allocate(vwind(im(region),jm(region),lm(region)))
      allocate(kvh(im(region),jm(region),lm(region)))
      allocate(ustr(im(region),jm(region)))
      allocate(wheat(im(region),jm(region)))
      allocate(wqflx(im(region),jm(region)))
      allocate(wobkl(im(region),jm(region)))
      allocate(wheatv(im(region),jm(region)))
      allocate(wtseff(im(region),jm(region)))

      pblh  => conv_dat(region)%blh
      t     => temper_dat(region)%data
      q     => humid_dat(region)%data
      gph   => gph_dat(region)%data
      pu    => pu_dat(region)%data
      pv    => pv_dat(region)%data
      dkg   => conv_dat(region)%dkg

      ! mid-layer pressure levels
      ! potential temperature field
      ! virtual potential temperature field
      !
      do l=1,lm(region)
         do j=jsr(region), jer(region)
            do i=isr(region), ier(region)
               pf(i,j,l)    = (phlb(i,j,l) + phlb(i,j,l+1) )*0.5
               wthm(i,j,l)  = t(i,j,l) * ( apzero/pf(i,j,l) )**zkappa
               zthvk(i,j,l) = wthm(i,j,l) * (1.0 + zcrdq*q(i,j,l))
            end do
         end do
      end do

      ! mid-layer heights (note: gph array is for levels 1:lm+1)
      !
      do l=1,lm(region)
         do j=jsr(region), jer(region)
            do i=isr(region), ier(region)
               intfac = (phlb(i,j,l)-pf(i,j,l)) / (phlb(i,j,l)-phlb(i,j,l+1))
               whm(i,j,l) = (gph(i,j,l)-gph(i,j,1)) * (1.0-intfac) + (gph(i,j,l+1)  -gph(i,j,1)) * intfac
               ! cmk bug corrected: was (gph(i,j,l+1)-gph(i,j,1))
            end do
         end do
      end do

      ! convert units of sensible and latent heat fluxes
      ! sshf: W/m2 / (J/kg/K) / (kg/m3) =  J/s/m2 /J kg K /kg m3 =  K m/s
      ! slhf: W/m2 / (J/kg) / (kg/m3) = J/s/m2 /J kg /kg m3 = m/s
      ! ustr is local velocity scale,  NOT from ECMWF surface stresses
      !
      do j=jsr(region), jer(region)
         do i=isr(region), ier(region)
            wheat(i,j) = sshf_dat(region)%data(i,j,1) * (gph(i,j,2)-gph(i,j,1)) * grav / &
                                  (phlb(i,j,1) - phlb(i,j,2)) / cp_air
            wqflx(i,j) = slhf_dat(region)%data(i,j,1) * (gph(i,j,2)-gph(i,j,1)) * grav / &
                                  (phlb(i,j,1) - phlb(i,j,2)) / Lvap
            ! if from surface stress:   ustr(:,:) = sqrt(sstr(:,:))
            ustr(i,j) = max(ustar_loc(region)%surf(i,j), 0.01)
            wheat(i,j)=-1*wheat(i,j)   !WP! bug fix
            wqflx(i,j)=-1*wqflx(i,j)   !WP! bug fix
         enddo
      enddo
      !
      !-----------------------------------------------------------------
      ! compute the free atmosphere vertical diffusion coefficients wkvf
      ! height-dependent asymptotic mixing length implemented
      ! Holtslag and Boville, 1993, J. Climate, 6, 1825-1842.
      !------------------------------------------------------------------

      ! first calculate wind velocities [m s-1] from the mass fluxes

      yres = dy/yref(region)
      xres = dx/xref(region)

      do j = 1,jm(region)
          lat(j) = ybeg(region) + 0.5 * yres + yres * (j-1)
      enddo

      dxx(:) = ae * xres * gtor * cos(lat(:)*gtor)
      dyy = ae * yres * gtor

      do l=1,lm(region)
        do j=jsr(region), jer(region)
          do i=isr(region), ier(region)
             uwind(i,j,l) = dxx(j)*(pu(i,j,l) + pu(i-1,j,l))*0.5 / m(i,j,l)
             ! cmk Bug repared, 06/2005 yres --> dyy
             vwind(i,j,l) = dyy*  (pv(i,j,l) + pv(i,j+1,l))*0.5 / m(i,j,l)
          enddo
        enddo
      enddo

      !
      ! Initialize gradient Richardson number and free tropospheric k values
      !

      if (okdebug) write(*,*) '  Richardson number'

      zrinub(:,:,:) = 0.0
      wkvf(:,:,:)   = 0.0
      zdup2(:,:,:)  = 0.0
      !
      do l=1,lm(region)-1
        do j=jsr(region), jer(region)
          do i=isr(region), ier(region)
            !
            ! if (delta gph < 1000) then jq=1
            jq = jqif( 1000.0, (gph(i,j,l+1)-gph(i,j,1)) )
            zlambdac=jq*300.+(1-jq)*(30. + 270. * exp (1.-(gph(i,j,l+1)-gph(i,j,1))/1000.))
            cml2=(1./( 1./zlambdac + 1./(vkarman*(gph(i,j,l+1)-gph(i,j,1))) ))**2.
            !
            !  vertical shear squared,
            !  min value of (delta v)**2 prevents zero shear
            !
            zdup2(i,j,l) = (uwind(i,j,l+1)-uwind(i,j,l))**2 + &
                                               (vwind(i,j,l+1)-vwind(i,j,l))**2
            zdup2(i,j,l) = max(zdup2(i,j,l),1.e-10)
            zdz   = whm(i,j,l+1) - whm(i,j,l)
            zdup2(i,j,l) = zdup2(i,j,l)/(zdz**2)
            !
            ! static stability (use virtual potential temperature)
            ! dv now calculated at interface height (lineair in pressure)
            !
            arg=(log(phlb(i,j,l+1))-log(pf(i,j,l)))/(log(pf(i,j,l+1))-log(pf(i,j,l)))
            zthva = zthvk(i,j,l) + arg * (zthvk(i,j,l+1)-zthvk(i,j,l))
            zsstab = grav/zthva*(zthvk(i,j,l+1) - zthvk(i,j,l))/zdz
            !
            ! gradient Richardson number
            !
            zrinub(i,j,l)  = zsstab/zdup2(i,j,l)
            !
            ! stability functions
            !
            zfunst = max(1. - 18.*zrinub(i,j,l),0.)
            zfunst = sqrt(zfunst)
            !
            zfstabh = 1.0 / (1.0 + 10.0*zrinub(i,j,l)*(1.0+8.0*zrinub(i,j,l)))
            !
            ! neutral diffusion coefficient
            !
            zkvn = cml2 * sqrt(zdup2(i,j,l))
            !
            ! correction with stability functions
            !
            jq = jqif( zrinub(i,j,l),0.0 )  ! thus: if (zrinub>0) then jq=1
            wkvf(i,j,l)=jq*(max(ckmin,zkvn*zfstabh))+(1-jq)*(max(ckmin,zkvn*zfunst))
          end do
        end do
      end do
      !
      ! Determine the boundary layer height
      !--------------------------------------------------------------------------------
      ! based on: Vogelezang and Holtslag, 1996, Bound. Layer Meteorol., 81, 245-269.
      !--------------------------------------------------------------------------------

      if (okdebug) write(*,*) '  boundary layer height'


      pblh(:,:) = 0.0  ! avoid non-initialised on edges.


      do j=jsr(region),jer(region)
        do i=isr(region),ier(region)
          ! compute bottom level virtual temperature and heat flux and Monin-Obukhov Length
          wtseff(i,j)  = wthm(i,j,1)*(1. + zcrdq*q(i,j,1))
          wheatv(i,j)  = wheat(i,j) + zcrdq*wthm(i,j,1)*wqflx(i,j)
          wobkl(i,j)   = -wtseff(i,j)*ustr(i,j)**3 / (grav*vkarman*(wheatv(i,j)+sign(1.e-10,wheatv(i,j))) )
          ! initialize
          jwork  = 1
          wrino  = 0.0
          do l=1,lm(region)

            vvk   = (uwind(i,j,l)-uwind(i,j,1))**2 + (vwind(i,j,l)-vwind(i,j,1))**2 + &
                                     zacb*ustr(i,j)*ustr(i,j)
            vvk   = max(vvk,tiny)

            tkv   = wthm(i,j,l)*(1. + zcrdq*q(i,j,l))
            dtkv  = tkv-wtseff(i,j)
            !
            ! Bulk Richardson number
            !
            zrino = grav*dtkv * (whm(i,j,l)-whm(i,j,1)) / (wtseff(i,j)*vvk)
            zrino = zrino+sign(1.e-10,zrino) ! prevent zrino becoming zero
            zrino = jwork*zrino

            jq    = jqif(ricr,zrino)  ! thus: if (zrino < ricr) then jq=1
            !
            ! calculate pblh from linear interpolation in Ri(bulk)
            !
            if ( l == 1 ) then
              pblh(i,j) = jq*pblh(i,j) + (1-jq) *  &
                    (                (ricr-wrino)/(zrino-wrino)*(whm(i,j,l)             ) )
            else
              pblh(i,j) = jq*pblh(i,j) + (1-jq) *  &
                    ( whm(i,j,l-1) + (ricr-wrino)/(zrino-wrino)*(whm(i,j,l)-whm(i,j,l-1)) )
            end if
            !
            ! first time zrino > ricr we set jwork to zero, thus all further times zrino=0
            ! and pblh does not change anymore
            !
            jwork = jq*jwork
            !
            ! if pblh is found already avoid dividing by zero next level by a faked value
            ! of wrino (jwork=0 already)
            !
            wrino = zrino+(1-jwork)*0.1

          end do

          pblh(i,j)  = (1-jwork)*pblh(i,j) + whm(i,j,lm(region))*jwork
          !
          ! second calculation of pblh including excess surface temperature using
          ! a convective velocity scale wsc (wm, not equal to w*!!!), using result
          ! of the first calculation of pblh
          !
          jq         = jqif(wheatv(i,j),0.) ! thus: if (wheatv(i,j) > 0.) then jq=1
          jwork = jq
          fmt = ( jq*(1.0 - binm*pblh(i,j)/wobkl(i,j)) + (1.0-jq) )**onet
          wsc     =  ustr(i,j)*fmt
          zexcess =  wheatv(i,j)*fak/wsc
          !
          ! improve pblh estimate under convective conditions using
          ! convective temperature excess (zexcess)
          !
          thvparc = wtseff(i,j)+jq*zexcess
          !
          jwork  = 1
          wrino  = 0.0
          !
          do l=2,lm(region)
            !
            ! Bulk Richardson number:
            !  if stable corrected with extra shear term
            !  if unstable now also corrected for temperature excess
            !
            vvk   = (uwind(i,j,l)-uwind(i,j,1))**2 + (vwind(i,j,l)-vwind(i,j,1))**2 + &
                                      zacb*ustr(i,j)*ustr(i,j)
            vvk   = max(vvk,tiny)
            tkv   = wthm(i,j,l)*(1. + zcrdq*q(i,j,l))
            zrino = grav*(tkv-thvparc) * (whm(i,j,l)-whm(i,j,1))/(vvk*wtseff(i,j))
            zrino = zrino+sign(1.e-10,zrino) ! prevent zrino becoming zero
            zrino = zrino*jwork

            jq    = jqif(ricr,zrino)  ! thus: if (zrino < ricr) then jq=0
            !
            ! calculate pblh from linear interpolation in Ri(bulk)
            !
            pblh(i,j)  = jq*pblh(i,j)+(1-jq)* (whm(i,j,l-1)+(ricr-wrino)/ &
                                                (zrino-wrino)*(whm(i,j,l)-whm(i,j,l-1)))
            !
            ! first time zrino > ricr we set jwork to zero, thus all further times zrino=0
            ! and pblh does not change anymore
            !
            jwork = jq*jwork
            !
            ! if pblh is found already avoid dividing by zero next level by a faked value
            ! of wrino (jwork=0 already)
            !
            wrino = zrino+(1-jwork)*0.1

          end do

          pblh(i,j) = (1-jwork)*pblh(i,j) + jwork*whm(i,j,lm(region))

        end do
      end do

      pblh = max(pblh, 100.0) ! set minimum value for pblh

      !
      ! Diffusion coefficients in the surface layer
      !
      !--------------------------------------------------------------------------------
      ! the revised LTG scheme (Beljaars and Viterbo, 1998)
      ! Modification - April 1, 1999: prevent zpblk to become too small in very stable conditions
      !--------------------------------------------------------------------------------

      if (okdebug) write(*,*) '  diffusion coeff'

      ! initialize output array with minimum value
      kvhmin = 0.1
      !
      ! Loop to calculate the diffusivity
      !
      do l=1,lm(region)-1
        do j=jsr(region), jer(region)
          do i=isr(region), ier(region)


            jq  = jqif(pblh(i,j),(gph(i,j,l+1)-gph(i,j,1))) ! if top of level hh < pblh jq=1 inside BL
            z   = gph(i,j,l+1)-gph(i,j,1)
            zh  = jq*z/pblh(i,j) ! only defined in BL (jq=1)

            ! zl must not be zero
            zl  = jq*z/wobkl(i,j)+(1-jq) ! z/L outside BL = 1
            zzh = jq*(1.-zh)**2           ! (1-z/h)^2 only in BL

            jq1 = jqif(wheatv(i,j),0.)    ! jq1 is 1 in unstable case
            jq2 = jq*(1-jq1)              ! jq2 is 1 in stable BL
            jq3 = jq*jq1                  ! jq3 is 1 in unstable BL

            jq6 = jqif(pblh(i,j), whm(i,j,l))
            jq7 = jqif(whm(i,j,l+1), pblh(i,j))
            jq6 = jq1*jq6*jq7  ! jq6 is 1 at at the kvh-level that is located between the
                               ! mid-levels that surround pblh

            jq1 = jqif(sffrac,zh)         ! jq1 is 1 in surface layer
            jq4 = jq3*jq1                 ! jq4 is 1 in unstable surface layer
            jq5 = jq3*(1-jq1)             ! jq5 is 1 outside unstable surface layer
            !
            ! calculate coefficients for momentum, the values for heat are obtained by dividing by Prantl number
            !
            ! stable and neutral:
            !
            cml2    = (1./( 1./(vkarman*z) + 1./450. ))**2.
            jqq     = jqif(zrinub(i,j,l),0.)
            !
            ! stability function for momentum
            !
            zfstabh  = jqq/(1.+10.*zrinub(i,j,l)*sqrt(1.+jqq*zrinub(i,j,l))) + (1-jqq)

            zpblk   = cml2 * sqrt(zdup2(i,j,l))*zfstabh
            zkvh    = jq2*max(zpblk,kvhmin)+(1-jq2)*wkvf(i,j,l) ! take free troposp. values outside BL
            !
            ! unstable case in surface layer
            !
            zslask = 1.-betah*zl
            zslask = jq4*zslask+(1-jq4)
            zpblk  = (1-jq4)*zkvh+jq4*(ustr(i,j)*vkarman*z*zzh*zslask**(1./3.))
            !
            ! unstable case above surface layer
            !
            ! evaluate stability function for momentum at height z = 0.1*pblh (top of surface layer)
            !
            zslask = 1.-ssfrac*betah*pblh(i,j)/wobkl(i,j)
            zslask = jq5*zslask+(1-jq5)
            fmt    = zslask**(1./3.)
            !
            ! use surface layer definition for w_m
            wsc    = ustr(i,j)*fmt
            zpblk  = (1-jq5)*zpblk + jq5*wsc*vkarman*z*zzh
            !
            ! Determine Prandt Number
            !
            ! Pr-number is assumed to be constant throughout the CBL, i.e. in surface and mixed layer
            ! NOTE: this is different from the original formulation
            ! zwstr not used if wheatv<0.
            !
            zwstr  = (abs(wheatv(i,j))*grav*pblh(i,j)/wtseff(i,j))**(1./3.)         ! bh

            obukov = jq3*wobkl(i,j)-(1-jq3)
            term1 = (1.-ssfrac*betah*pblh(i,j)/obukov)**(1./3.)
            term2 = sqrt(1.-ssfrac*betah*pblh(i,j)/obukov)
            term3 = ccon*fakn*zwstr/(ustr(i,j)*(1.-ssfrac*betah*pblh(i,j)/obukov)**(1./3.))

            pr = jq3*(term1/term2+term3)+(1-jq3)
            !
            ! NOTE that kvh applies to the top of the level
            !
            kvh(i,j,l) = jq3*max(zpblk/pr,kvhmin)+(1-jq3)*zkvh

            ! if in the entrainment zone, override with prescribed entrainment
            thvgrad = (zthvk(i,j,l+1)-zthvk(i,j,l))/(whm(i,j,l+1)-whm(i,j,l))
            kvhentr = 0.2*wheatv(i,j)/thvgrad
            kvh(i,j,l) = jq6*kvhentr + (1-jq6)*kvh(i,j,l)

          end do
        end do
      end do

      ! calculate vertical diffusion mass flux:
      ! air mass exchanged between layer l and l+1 through interface l
      ! due to vertical diffusion in [(kg air)/s]
      dkg = 0.0
      do  l=1,lmax_conv-1
        do  j=jsr(region), jer(region)
          do  i=isr(region), ier(region)
            ! dh between layer centers:
            !   dh(l+1/2) = [h(l+1)+h(l+2)]/2 - [h(l)+h(l+1)]/2 = [h(l+2)-h(l)]/2
            ! average mass:
            !   mh(l+1/2) = [m(l)+m(l+1)]/2
            ! conversion:
            !   dkg = kvh mh / dh**2
            ! if original 'kvh' has unit [m2/s], then new unit is [kg/s]
            dkg(i,j,l)=max(0.,kvh(i,j,l)) * 2. * (m(i,j,l+1)+m(i,j,l)) / (gph(i,j,l+2)-gph(i,j,l))**2
            ! CMK dec 2002 ----> gph changed to 1--->lm+1 boundaries
          end do
        end do
      end do
      dkg(:,:,lmax_conv) = 0.0

      ! Also store the diffusion coefficients
      conv_dat(region)%kvh = kvh

      deallocate(wkvf)
      deallocate(zdup2)
      deallocate(zrinub)
      deallocate(wthm)
      deallocate(zthvk)
      deallocate(pf)
      deallocate(whm)
      deallocate(ustr)
      deallocate(wheat)
      deallocate(wqflx)
      deallocate(wobkl)
      deallocate(wheatv)
      deallocate(wtseff)
      deallocate(uwind)
      deallocate(vwind)

      nullify(pblh, t, q, gph, dkg)
      deallocate(kvh)

      if (okdebug) write(*,*) 'start of bldiff'

  end subroutine bldiff


  ! ***


  !-----------------------------------------------
  !
  ! Purpose:
  ! -------
  ! this subroutine calculate the min,mean,max of a real field
  !
  ! External
  ! --------
  ! none
  !
  ! Reference
  ! ---------
  ! none
  !------------------------------------

  subroutine dd_field_statistics(name,field,ix,jx)

    implicit none

    ! --- in/out -------------------------------

    character(len=*),intent(in)        :: name
    integer,intent(in)                 :: ix,jx
    real,dimension(ix,jx),intent(in)   :: field

    ! --- local --------------------------------

    integer :: i,j,ntel_non_zero
    real :: maxf,minf,meanf,mean_non_zero

    ! --- begin -------------------------------

    maxf=-1e20
    minf=1e12
    meanf=0.
    mean_non_zero=0.
    ntel_non_zero=0
    do i=1,ix
    do j=1,jx
      meanf=meanf+field(i,j)
      maxf=max(maxf,field(i,j))
      minf=min(minf,field(i,j))
      if (field(i,j).ne.0) then
        ntel_non_zero=ntel_non_zero+1
        mean_non_zero=mean_non_zero+field(i,j)
      endif
    enddo
    enddo
    meanf=meanf/ix/jx
    if (ntel_non_zero.gt.0) mean_non_zero= mean_non_zero/ntel_non_zero
    write(6,'(a10,4(a3,1x,1pe10.3))') name,'max',maxf,'min',minf,'mean',meanf,'mn0',mean_non_zero

  end subroutine dd_field_statistics


  ! ****


  !------------------------------------
  !
  ! Purpose:
  ! -------
  ! this subroutine calculate the min,mean,max of a field of int8
  !
  ! External
  ! --------
  ! none
  !
  ! Reference
  ! ---------
  ! none
  !------------------------------------

  subroutine dd_field_statistics_i8(name,field,ix,jx)

    implicit none

    ! --- in/out ----------------------------

    character(len=*),intent(in)                   :: name
    integer,intent(in)                            :: ix,jx
    integer(kind=1),dimension(ix,jx),intent(in)   :: field

    ! --- local ----------------------------

    integer :: i,j,ntel_non_zero
    real :: maxf,minf,meanf,mean_non_zero

    ! --- begin ---------------------------

    maxf=-1e20
    minf=1e12
    meanf=0.
    mean_non_zero=0.
    ntel_non_zero=0
    do i=1,ix
      do j=1,jx
        meanf=meanf+field(i,j)
        maxf=max(maxf,real(field(i,j)))
        minf=min(minf,real(field(i,j)))
        if (field(i,j).ne.0) then
          ntel_non_zero=ntel_non_zero+1
          mean_non_zero=mean_non_zero+field(i,j)
        endif
      enddo
    enddo
    meanf=meanf/ix/jx
    if (ntel_non_zero.gt.0) mean_non_zero= mean_non_zero/ntel_non_zero
    write(6,'(a10,4(a3,1x,1pe10.3))') name,'max',maxf,'min',minf,'mean',meanf,'mn0',mean_non_zero

  end subroutine dd_field_statistics_i8


  ! ***


  !
  ! define vectorizable 'if' function
  ! if xxxz>yyyz jqif=1 else jqif=0
  !

  integer function jqif(xxxz,yyyz)

    implicit none

    ! --- in/out ---------------------------

    real, intent(in)    ::  xxxz, yyyz

    ! --- begin ----------------------------

    !
    !jqif = nint( 0.5 - sign( 0.5, -dim(xxxz,yyyz)) )
    !
    jqif = nint( 0.5 + sign( 0.5, xxxz-yyyz ) )

  end function jqif


  !=====================================================================================================
  !=====================================================================================================

  !subroutine diffusion_filename( region_name, t1, t2, filename, status )

    !use GO  , only : TDate
    !use GO  , only : pathsep
    !use GO  , only : Get
    !use Dims, only : revert

    !implicit none

    !! --- in/out ----------------------------------------------

    !character(len=*), intent(in)     ::  region_name
    !type(TDate), intent(in)          ::  t1, t2
    !character(len=*), intent(out)    ::  filename
    !integer, intent(out)             ::  status

    !! --- const -----------------------------------------------

    !character(len=*), parameter      :: rname = mname//'diffusion_filename'

    !! --- local -----------------------------------------------

!!    integer       :: idater1(6)
!!    integer       :: idater2(6)
    !integer     :: idater(6)

    !! --- begin -----------------------------------------------

    !! time range in 'forward' order, used in filename:
    !if ( revert == 1 ) then
      !! forward order:
      !call Get( t1, time6=idater )
!!      call Get( t2, time6=idater2 )
    !else
      !! reverse order:
      !call Get( t2, time6=idater )
!!      call Get( t1, time6=idater2 )
    !end if

    !! filename:
    !!    <dir>/region_name/YYYY/MM/DD/dkg_YYYYMMDDHHMM.hdf
    !!
    !write(filename, '(4a, i4.4, a1, i2.2, a1, i2.2, a1, "dkg_", i4, 4i2.2, ".hdf")') &
            !trim(diffusion_dir), pathsep, trim(region_name), pathsep, &
            !idater(1), pathsep, idater(2), pathsep, idater(3), pathsep, idater(1:5)

    !! ok
    !status = 0

  !end subroutine diffusion_filename


  ! ***


  !
  ! Read vertical diffusion from file and store in conv_dat%dkg
  !
  ! Output:
  !  o status      0 = ok
  !               -1 = file not available for at least one region
  !                1 = read error
  !

  subroutine read_diffusion(region, t1, t2, status)

    ! --- modules ------------------------------

    use file_netcdf
    use Go,             only : TDate, operator(==), assignment(=), Pretty
    use Go,             only : Get, NewDate
    use dims,           only : revert, region_name, newday, newsrun, lmax_conv
    use global_data,    only : conv_dat

    implicit none

    ! --- in/out ----------------------------------------------

    integer, intent(in)             :: region
    type(TDate), intent(in)         ::  t1, t2
    integer, intent(out)            ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/read_diffusion_nc'

    ! --- local -----------------------------------------------

    character(len=MAX_FILENAME_LEN) :: fname
    integer, dimension(6)           :: file_date
    logical                         :: new_file, file_exist, index_found
    integer                         :: tidx, nc_id, i, nsteps
    type(TDate)                     :: tplus, tminus
    integer, allocatable            :: t_temp(:,:)

    ! --- begin -----------------------------------------------

    if (newday .or. newsrun) then
        ! It's a new day! Rejoice, be optimistic! Assume the dkg file is complete!
        file_complete(region) = .true.
    end if

    ! Does a new file need to be read?
    if ( revert == 1 ) then
        ! forward order:
        call Get( t1, time6=file_date )
    else
        ! reverse order:
        call Get( t2, time6=file_date )
    end if

    ! Do we need to read a new file? Yes if file_date is not the same as cur_date.
    new_file = any(file_date(1:3) /= dkg_ncfile(region)%cur_date)

    !write(gol, '(a, " called for region ", i1, " with t1 = ", a, " and t2 = ", a, ", need to read new file = ", l1)') &
        !rname, region, trim(Pretty(t1)), trim(Pretty(t2)), new_file
    !call goPr

    if (new_file) then
        ! Is the file potentially complete? If it's not, no need to look further. Calculate dkg instead.
        ! On a new day, file_complete will be true. It only gets set to false if a previous attempt to read it
        ! has resulted in a 'file not found'.
        if (.not. file_complete(region)) then
            status = -1
            return
        end if

        ! get file name:
        call diffusion_filename( region_name(region), t1, t2, fname, status )
        IF_NOTOK_RETURN(status=1)

        ! If the file is not there, return with status = -1, so that calc_kzz can be called.
        ! However, once a file for a particular day has not been found, always return status = -1
        ! for that day. This is to prevent reading of incomplete files.
        inquire( file=trim(fname), exist=file_exist )
        if ( .not. file_exist ) then
            ! the file doesn't exist for a day, so doesn't exist for all time steps in that day
            file_complete(region) = .false. ! this will be reset when a day boundary is crossed
            write (gol,'(a,": dkg file not found for region ", a, " for date ", i4.4, "-", i2.2, "-", i2.2)') &
                rname, region_name(region), file_date(1:3); call goPr
            status = -1; return
        end if

        ! deallocate the old arrays
        if (allocated(dkg_ncfile(region)%dkg)) deallocate(dkg_ncfile(region)%dkg)
        if (allocated(dkg_ncfile(region)%blh)) deallocate(dkg_ncfile(region)%blh)
        if (allocated(dkg_ncfile(region)%t1))  deallocate(dkg_ncfile(region)%t1)
        if (allocated(dkg_ncfile(region)%t2))  deallocate(dkg_ncfile(region)%t2)

        nc_id = nc_open(fname, 'r', status)
        IF_NOTOK_RETURN(status=1)

        nsteps = nc_get_dim(nc_id, 'steps')

        ! check if the dkg coefficients were calculated with the same lmax_conv
        if (nc_get_dim(nc_id, 'lmax_conv') /= lmax_conv) then
            write(gol, '("lmax_conv = ", i2, " in ", a, " does not match ", i2, " in model")') &
                nc_get_dim(nc_id, 'lmax_conv'), trim(fname), lmax_conv
            call goErr
            status = 1
            IF_NOTOK_RETURN(status=1)
        end if

        dkg_ncfile(region)%nsteps = nsteps
        dkg_ncfile(region)%blh = nc_read_var(nc_id, 'blh')
        dkg_ncfile(region)%dkg = nc_read_var(nc_id, 'dkg')

        ! the t1 and t2 arrays require a bit more work
        t_temp = nc_read_var(nc_id, 't1')
        allocate(dkg_ncfile(region)%t1(nsteps))
        do i = 1, nsteps
            dkg_ncfile(region)%t1(i) = NewDate(time6=t_temp(:,i))
        end do
        deallocate(t_temp)

        t_temp = nc_read_var(nc_id, 't2')
        allocate(dkg_ncfile(region)%t2(nsteps))
        do i = 1, nsteps
            dkg_ncfile(region)%t2(i) = NewDate(time6=t_temp(:,i))
        end do
        deallocate(t_temp)

        dkg_ncfile(region)%cur_date = file_date(1:3)

        call nc_close(nc_id)

    else

        nsteps = dkg_ncfile(region)%nsteps

    end if

    ! Now find the correct time index
    if ( revert == 1 ) then
        ! forward order:
        tminus = t1
        tplus = t2
    else
        tminus = t2
        tplus = t1
    end if

    ! we now need to search for tminus in t1 and tplus in t2, and see if the indices match
    index_found = .false.
    do i = 1, nsteps
        if ((dkg_ncfile(region)%t1(i) == tminus) .and. (dkg_ncfile(region)%t2(i) == tplus)) then
            index_found = .true.
            tidx = i
            exit
        end if
    end do

    if (.not. index_found) then
        ! A file for this day and zoom region exists, but does not have the diffusion coefficients for the correct
        ! time step. This can happen, for example, if the model is run with a different (base) dynamic timestep
        ! (ndyn in the code) with the same zoom config. In that case, call calc_kzz and store the new values in the
        ! file. We do not, however, want to set file_complete(region) to .false., because on this same day the dkg
        ! and blh values may already exist in the file for a different time step. If we just return status = -1, then
        ! during the next time step the file will be read again. Will be a bit long during a first run with different
        ! ndyn, but eventually the file will contain dkg and blh for all time steps.
        write(gol,'(a, ": diffusion coeff between ", a, " and ", a, " for region ", a, " does not exist in file, calculating")') &
            rname, trim(Pretty(t1)), trim(Pretty(t2)), region_name(region)
        call goPr
        status = -1
        return
    end if

    !write(gol, '(a, ": found tidx = ", i2)') rname, tidx
    !call goPr

    conv_dat(region)%dkg = dkg_ncfile(region)%dkg(:,:,:,tidx)
    conv_dat(region)%blh = dkg_ncfile(region)%blh(:,:,tidx)

    status = 0

  end subroutine read_diffusion

  !subroutine read_diffusion(region, t1, t2, status)

    !! --- modules ------------------------------

    !use GO                     , only : TDate
    !use dims,                    only : region_name
    !use global_data,             only : conv_dat
    !use file_hdf,                only : Init, ReadData, Done, THdfFile, TSds
!!    use Go,                      only : Pretty
!!    use dims,                    only : newday, newsrun

    !implicit none

    !! --- in/out ----------------------------------------------

    !integer, intent(in)             :: region
    !type(TDate), intent(in)         ::  t1, t2
    !integer, intent(out)            ::  status

    !! --- const -----------------------------------------------

    !character(len=*), parameter      :: rname = mname//'/read_diffusion'

    !! --- local -----------------------------------------------

    !character(len=MAX_FILENAME_LEN)               :: fname
    !logical                          :: exist
    !type(THdfFile)                   :: hdf
    !type(TSds)                       :: sds

    !! --- begin -----------------------------------------------

!!    write(*,'(a, " called with t1 = ", a, ", t2 = ", a, ", newday = ", l1, ", newsrun = ", l1)') &
!!        rname, trim(Pretty(t1)), trim(Pretty(t2)), newday, newsrun

    !! get file name:
    !call diffusion_filename( region_name(region), t1, t2, fname, status )
    !IF_NOTOK_RETURN(status=1)

    !! present ?
    !inquire( file=trim(fname), exist=exist )
    !! leave with warning status if not found:
    !if ( .not. exist ) then
      !write (gol,'(a,": dkg file not found for region",i2)') rname, region; call goPr
      !write (gol,'(a,":   ",a)') rname, trim(fname); call goPr
      !!call GO_Timer_End( itim_read_diffusion, status )
      !!IF_NOTOK_RETURN(status=1)
      !status = -1; return
    !end if

    !! info ...
    !!write(*,'(a,": reading ", a)') rname, trim(fname)

    !! Open file
    !call Init( hdf, fname, 'read', status )
    !IF_NOTOK_RETURN(status=1)

    !! Read dkg
    !call Init( sds, hdf, 'dkg', status )
    !IF_NOTOK_RETURN(status=1)

    !call ReadData( sds, conv_dat(region)%dkg, status )
    !IF_NOTOK_RETURN(status=1)

    !call Done( sds, status )
    !IF_NOTOK_RETURN(status=1)

    !! Read blh
    !call Init( sds, hdf, 'blh', status )
    !IF_NOTOK_RETURN(status=1)

    !call ReadData( sds, conv_dat(region)%blh, status )
    !IF_NOTOK_RETURN(status=1)

    !call Done( sds, status )
    !IF_NOTOK_RETURN(status=1)

    !! Close file
    !call Done( hdf, status )
    !IF_NOTOK_RETURN(status=1)


    !! ok
    !status = 0

  !end subroutine read_diffusion


  ! ***


  subroutine write_diffusion(region, t1, t2, status)

    use netcdf
    use file_netcdf,    only : nc_open, nc_close, nc_get_dim
    use dims,           only : region_name, im, jm, lm, lmax_conv, revert
    use global_data,    only : conv_dat
    use misctools,      only : check_dir
    use Go,             only : Get

    implicit none

    ! --- in/out ----------------------------------------------

    integer, intent(in)             :: region
    type(TDate), intent(in)         ::  t1, t2
    integer, intent(out)            ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'/write_diffusion_nc'

    ! --- local -----------------------------------------------

    character(len=MAX_FILENAME_LEN) :: fname
    logical                         :: old_file
    integer                         :: output_fid, tidx, io_stat
    integer                         :: dim_im, dim_jm, dim_lm, dim_lmax, dim_t, dim_d
    integer                         :: var_dkg, var_kvh, var_t1, var_t2, var_blh, deflate, shuffle
    integer                         :: temp_date(6)

    ! --- begin -----------------------------------------------

    ! This should never be called from the adjoint, because diffusion files should be written during a forward run
    if (revert /= 1) then
        write(gol, '(a, " should only be called from a forward run")') rname; call goErr
        status = 1
        IF_NOTOK_RETURN(status=1)
    end if

    ! get file name:
    call diffusion_filename( region_name(region), t1, t2, fname, status )
    IF_NOTOK_RETURN(status=1)

    ! create the folder if it doesn't already exist
    call check_dir(fname)

    ! is the file already there?
    inquire(file=trim(fname), exist=old_file)

    ! if the file is not there, create it and write the dimensions, etc.
    if (.not. old_file) then
        output_fid = nc_open(fname, 'c', status)
        IF_NOTOK_RETURN(status=1)
        
        io_stat = nf90_def_dim(output_fid, 'im', im(region), dim_im) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_def_dim(output_fid, 'jm', jm(region), dim_jm) ; CHECK_NCSTAT(io_stat)
        if (write_kvh) then
            io_stat = nf90_def_dim(output_fid, 'lm', lm(region), dim_lm) ; CHECK_NCSTAT(io_stat)
        end if
        io_stat = nf90_def_dim(output_fid, 'lmax_conv', lmax_conv, dim_lmax) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_def_dim(output_fid, 'steps', NF90_UNLIMITED, dim_t) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_def_dim(output_fid, 'idate', 6, dim_d) ; CHECK_NCSTAT(io_stat)

        deflate = 0
        shuffle = 0
        if (deflate_lvl > 0) then
            deflate = 1
            shuffle = 1
        end if

        io_stat = nf90_def_var(output_fid, 'dkg', NF90_DOUBLE, (/dim_im, dim_jm, dim_lmax, dim_t/), var_dkg) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_def_var_deflate(output_fid, var_dkg, shuffle, deflate, deflate_lvl) ; CHECK_NCSTAT(io_stat)

        if (write_kvh) then
            io_stat = nf90_def_var(output_fid, 'kvh', NF90_DOUBLE, (/dim_im, dim_jm, dim_lm, dim_t/), var_kvh) ; CHECK_NCSTAT(io_stat)
            io_stat = nf90_def_var_deflate(output_fid, var_kvh, shuffle, deflate, deflate_lvl) ; CHECK_NCSTAT(io_stat)
        end if

        io_stat = nf90_def_var(output_fid, 'blh', NF90_DOUBLE, (/dim_im, dim_jm, dim_t/), var_blh) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_def_var_deflate(output_fid, var_blh, shuffle, deflate, deflate_lvl) ; CHECK_NCSTAT(io_stat)

        io_stat = nf90_def_var(output_fid, 't1',  NF90_INT,    (/dim_d, dim_t/), var_t1) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_def_var(output_fid, 't2',  NF90_INT,    (/dim_d, dim_t/), var_t2) ; CHECK_NCSTAT(io_stat)

        ! put comments to say that t1 and t2 are not sorted, and put units of blh
        io_stat = nf90_put_att(output_fid, var_t1,  'description',  'starting point of time step over which dkg and blh are valid') ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_put_att(output_fid, var_t1,  'comment',      't1 may not be sorted') ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_put_att(output_fid, var_t2,  'description',  'ending point of time step over which dkg and blh are valid') ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_put_att(output_fid, var_t2,  'comment',      't2 may not be sorted') ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_put_att(output_fid, var_blh, 'description',  'TM5 diagnosed boundary layer height') ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_put_att(output_fid, var_blh, 'unit',         'meters above ground level') ; CHECK_NCSTAT(io_stat)
        if (write_kvh) then
            io_stat = nf90_put_att(output_fid, var_kvh, 'description',  'Diffusion coefficient') ; CHECK_NCSTAT(io_stat)
            io_stat = nf90_put_att(output_fid, var_kvh, 'unit',         'meter^2/second') ; CHECK_NCSTAT(io_stat)
        end if

        io_stat = nf90_enddef(output_fid) ; CHECK_NCSTAT(io_stat)

        tidx = 1

    else
        output_fid = nc_open(fname, 'a', status)
        IF_NOTOK_RETURN(status=1)
        ! read the variable and dimension IDs
        io_stat = nf90_inq_dimid(output_fid, 'im', dim_im) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_inq_dimid(output_fid, 'jm', dim_jm) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_inq_dimid(output_fid, 'lmax_conv', dim_lmax) ; CHECK_NCSTAT(io_stat)
        if (write_kvh) then
            io_stat = nf90_inq_dimid(output_fid, 'lm', dim_lm) ; CHECK_NCSTAT(io_stat)
        end if
        io_stat = nf90_inq_dimid(output_fid, 'steps', dim_t) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_inq_dimid(output_fid, 'idate', dim_d) ; CHECK_NCSTAT(io_stat)

        io_stat = nf90_inq_varid(output_fid, 'dkg', var_dkg) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_inq_varid(output_fid, 'blh', var_blh) ; CHECK_NCSTAT(io_stat)
        if (write_kvh) then
            io_stat = nf90_inq_varid(output_fid, 'kvh', var_kvh) ; CHECK_NCSTAT(io_stat)
        end if
        io_stat = nf90_inq_varid(output_fid, 't1', var_t1) ; CHECK_NCSTAT(io_stat)
        io_stat = nf90_inq_varid(output_fid, 't2', var_t2) ; CHECK_NCSTAT(io_stat)

        tidx = nc_get_dim(output_fid, 'steps') + 1

    end if

    io_stat = nf90_put_var(output_fid, var_dkg, conv_dat(region)%dkg(1:im(region), 1:jm(region), 1:lmax_conv), (/1, 1, 1, tidx/)) ; CHECK_NCSTAT(io_stat)
    io_stat = nf90_put_var(output_fid, var_blh, conv_dat(region)%blh(1:im(region), 1:jm(region)), (/1, 1, tidx/)) ; CHECK_NCSTAT(io_stat)
    if (write_kvh) then
        io_stat = nf90_put_var(output_fid, var_kvh, conv_dat(region)%kvh(1:im(region), 1:jm(region), 1:lm(region)), (/1, 1, 1, tidx/)) ; CHECK_NCSTAT(io_stat)
    end if

    ! convert t1 and t2 to integers before writing them out
    call Get( t1, time6=temp_date )
    io_stat = nf90_put_var(output_fid, var_t1, temp_date, (/1, tidx/)) ; CHECK_NCSTAT(io_stat)

    call Get( t2, time6=temp_date )
    io_stat = nf90_put_var(output_fid, var_t2, temp_date, (/1, tidx/)) ; CHECK_NCSTAT(io_stat)

    call nc_close(output_fid)

    status = 0

  end subroutine write_diffusion

  !subroutine write_diffusion(region, t1, t2, status)

    !!
    !! Write vertical diffusion to files (for all regions)
    !!
    !! Output:
    !!  o status      0 = ok, else not ok
    !!

    !! --- modules ------------------------------

    !use GO                     , only : TDate
    !use dims,                    only : region_name
    !use global_data,             only : conv_dat
    !use file_hdf,                only : Init, WriteData, Done, THdfFile, TSds
    !use misctools,               only : check_dir

    !implicit none

    !! --- in/out ----------------------------------------------

    !integer, intent(in)             :: region
    !type(TDate), intent(in)         ::  t1, t2
    !integer, intent(out)            ::  status

    !! --- const -----------------------------------------------

    !character(len=*), parameter      :: rname = mname//', write_diffusion'

    !! --- local -----------------------------------------------

    !character(len=MAX_FILENAME_LEN)               :: fname
    !type(THdfFile)                   :: hdf
    !type(TSds)                       :: sds

    !! --- begin -----------------------------------------------

    !! get file name:
    !call diffusion_filename( region_name(region), t1, t2, fname, status )
    !IF_NOTOK_RETURN(status=1)

    !! info ...
    !write (gol,'(a,": creating ", a)') rname, trim(fname); call goPr
    !call check_dir(fname)

    !! Create file
    !call Init( hdf, fname, 'create', status )
    !IF_NOTOK_RETURN(status=1)

    !! Write dkg
    !call Init( sds, hdf, 'dkg', shape(conv_dat(region)%dkg), 'real(8)', status )
    !IF_NOTOK_RETURN(status=1)

    !call WriteData( sds, conv_dat(region)%dkg, status )
    !IF_NOTOK_RETURN(status=1)

    !call Done( sds, status )
    !IF_NOTOK_RETURN(status=1)


    !! Write blh
    !call Init( sds, hdf, 'blh', shape(conv_dat(region)%blh), 'real(8)', status )
    !IF_NOTOK_RETURN(status=1)

    !call WriteData( sds, conv_dat(region)%blh, status )
    !IF_NOTOK_RETURN(status=1)

    !call Done( sds, status )
    !IF_NOTOK_RETURN(status=1)

    !! Close file
    !call Done( hdf, status )
    !IF_NOTOK_RETURN(status=1)

    !status = 0

  !end subroutine write_diffusion

  subroutine diffusion_filename( region_name, t1, t2, filename, status )

    use GO  , only : TDate
    use GO  , only : pathsep
    use GO  , only : Get
    use Dims, only : revert

    implicit none

    ! --- in/out ----------------------------------------------

    character(len=*), intent(in)     ::  region_name
    type(TDate), intent(in)          ::  t1, t2
    character(len=*), intent(out)    ::  filename
    integer, intent(out)             ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter      :: rname = mname//'diffusion_filename'

    ! --- local -----------------------------------------------

    integer     :: idater(6)

    ! --- begin -----------------------------------------------

    ! time range in 'forward' order, used in filename:
    if ( revert == 1 ) then
      ! forward order:
      call Get( t1, time6=idater )
    else
      ! reverse order:
      call Get( t2, time6=idater )
    end if

    write(filename, '(4a, i4.4, a1, i2.2, a1, "dkg_", i4, 2i2.2, ".nc")') &
            trim(diffusion_dir), pathsep, trim(region_name), pathsep, &
            idater(1), pathsep, idater(2), pathsep, idater(1:3)

    ! ok
    status = 0

  end subroutine diffusion_filename


end module diffusion
