!
! Subgrid stuff from TMPP.
!
! Copied from:
!   TMPP/SOURCE/tmpp_sub_subg.f90
!   TMPP/SOURCE/tmpp_sub_work.f90
!   TMPP/SOURCE/tmpp_geometry.f90
!   TMPP/SOURCE/phys_geopot.f90
!

module phys_convec_tmpp

  implicit none

  ! --- in/out -------------------------------

  private

  public   ::  subscal
  public   ::  subscal_2d

  public   ::  mid2bound_uv
  public   ::  mid2bound_w
  public   ::  mid2bound_t
  public   ::  mid2bound_q
  public   ::  mid2bound_z
  public   ::  mid2bound_p

  public   ::  phlev_ec_rev

  public   ::  geopot

contains


  ! ==========================================================
  ! ===
  ! ===  TMPP/SOURCE/tmpp_sub_subg.f90
  ! ===
  ! ==========================================================


  !-----------------------------------------------------------------------
  ! calculate subscal parameters at a given horizontal location.
  ! Three different pressure levels are used in this subroutine:
  ! ppress (boundaries ECMWF levels)
  ! ppresg (boundaries TM3)
  ! zpl (centers ECMWF levels)
  ! zplg (centers TM3)
  !
  !-----------------------------------------------------------------------

  !subroutine subscal(pz,ppress,pu,pv,pw,pt,pq,pqac,pqam,pevap, &
  !                   pentug,pdetug,pentdg,pdetdg,pdkg,pzg)

  subroutine subscal_2d( np, lme, at, bt, &
                     pz,ppress,pw,pt,pq,pqac,pqam,pevap, &
                     pentug,pdetug,pentdg,pdetdg)

    ! --- in/out -----------------------------------

    integer, intent(in)  ::  np
    integer, intent(in)  ::  lme
    real, intent(in)     ::  at(0:lme), bt(0:lme)
    real, intent(in)     ::  pz(np,0:lme)
    real, intent(in)     ::  ppress(np,0:lme)
    real, intent(in)     ::  pw(np,0:lme)
    real, intent(in)     ::  pt(np,0:lme)
    real, intent(in)     ::  pq(np,0:lme)
    real, intent(in)     ::  pqac(np,0:lme)
    real, intent(in)     ::  pqam(np,0:lme)
    real, intent(in)     ::  pevap(np)
    real, intent(out)    ::  pentug(np,lme)
    real, intent(out)    ::  pdetug(np,lme)
    real, intent(out)    ::  pentdg(np,lme)
    real, intent(out)    ::  pdetdg(np,lme)

    ! --- local -------------------------------------

    integer     :: i

    ! --- begin ------------------------------------

    print *, '                subscal_2d'
    !$OMP PARALLEL &
    !$OMP       default ( none ) &
    !$OMP       shared  ( np, lme, at, bt ) &
    !$OMP       shared  ( pz, ppress, pw, pt, pq ) &
    !$OMP       shared  ( pqac, pqam, pevap ) &
    !$OMP       shared  ( pentug, pdetug, pentdg, pdetdg ) &
    !$OMP       private ( i )
    !$OMP   DO
    do i = 1, np
      call subscal( lme, at, bt, &
                     pz(i,:), ppress(i,:), pw(i,:), pt(i,:), pq(i,:), &
                     pqac(i,:), pqam(i,:), pevap(i), &
                     pentug(i,:), pdetug(i,:), pentdg(i,:), pdetdg(i,:) )
    end do
    !$OMP   END DO
    !$OMP END PARALLEL
    print *, '                end subscal_2d'

  end subroutine subscal_2d

  ! *

  subroutine subscal( lme, at, bt, &
                     pz,ppress,pw,pt,pq,pqac,pqam,pevap, &
                     pentug,pdetug,pentdg,pdetdg) !,pdkg,pzg)

    ! >>> adhoc: not from ECMWF to TM levels >>>>>>>>>>>>>>>>>>
    !     (set tm stuff to ec stuff)
    !use Geometry, only : lm
    !use Geometry, only : at => a_tm, bt => b_tm
!    use Geometry, only : lm => lme
!    use Geometry, only : at => a_ec, bt => b_ec
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!    use Geometry, only : lme

    use Num, only : IntervalQuad_Lin
    use Num, only : interp_muherm

    ! --- const ------------------------------------

!    integer, parameter  ::  jplm = lme

    !integer, parameter  ::  jplmg = lm
    !integer, parameter  ::  jplmgp1 = lm+1
!    integer, parameter  ::  jplmg = lme
!    integer, parameter  ::  jplmgp1 = lme+1

    real, parameter   ::  ppg = 9.80665

    ! --- in/out ------------------------------------

    integer              ::  status

    integer, intent(in)  ::  lme
    real, intent(in)     ::  at(0:lme), bt(0:lme)

    real, intent(in)   :: pz(0:lme)
    real, intent(in)   :: ppress(0:lme)
    real, intent(in)   :: pw(0:lme)
    real, intent(in)   :: pt(0:lme)
    real, intent(in)   :: pq(0:lme)
    real, intent(in)   :: pqac(0:lme)
    real, intent(in)   :: pqam(0:lme)
    real, intent(in)   :: pevap

    !real, intent(in)   :: pu(0:lme)
    !real, intent(in)   :: pv(0:lme)

    ! entrainment, detrainment rates, diffusion coefficient

    real, intent(out)  :: pentug(lme)
    real, intent(out)  :: pdetug(lme)
    real, intent(out)  :: pentdg(lme)
    real, intent(out)  :: pdetdg(lme)

!    real, intent(out)  :: pdkg(lme)

!    real, intent(out)  :: pzg(0:lme)

    ! --- local ------------------------------------

    integer :: jl,ilast

    real :: sumup, sumdown


    ! output variables on TM vertical grid

    real ::  zamu(0:lme),zptc(0:lme),zpqc(0:lme),zplc(0:lme), &
         zpqp(0:lme),ppresg(0:lme)
!    real ::  zpdk(0:lme)

    real ::  zpam(lme),zped(lme),zpdd(lme),zpeu(lme),zpdu(lme)
    real ::  zpl(lme),zplg(lme)
    real ::  zam

    integer ::  lct

    ! --- begin -------------------------------------------

    ! calculate convection by clouds

    call cloud( lme, pz,ppress,pt,pq,pw,pqac,pqam,pevap, &
               zamu,zptc,zpqc,zplc,zpqp, &
               zpam,zpeu,zpdu,zped,zpdd)


!    ! calculate turbulent diffusion coefficient
!
!    call louis(pz,pt,pu,pv,zpdk)

    ! generate pressures on TM3 level boundaries and centers respectively

    do jl=0,lme
      ppresg(jl)=at(jl)+ppress(lme)*bt(jl)
    end do

    do jl=1,lme
      zplg(jl)=0.5*(ppresg(jl)+ppresg(jl-1))
    end do

    ! --- Vertical averaging of variables defined on level centers

    ! * interpolate entrainment/detrainment rates on TM3 vertical coordinate:

    ! Dirk Olivie, 12 May 2004
    ! No interpolation needed if lme = lme

    if ( lme .ne. lme ) then

      do jl=1,lme
        zpl(jl)=0.5*(ppress(jl)+ppress(jl-1))
        zpeu(jl)=zpeu(jl)/(ppress(jl)-ppress(jl-1))
        zpdu(jl)=zpdu(jl)/(ppress(jl)-ppress(jl-1))
        zped(jl)=zped(jl)/(ppress(jl)-ppress(jl-1))
        zpdd(jl)=zpdd(jl)/(ppress(jl)-ppress(jl-1))
      end do
      ilast=1
      do jl=1,lme
        call IntervalQuad_Lin(zpl,zpeu,ppresg(jl-1),ppresg(jl),pentug(jl), ilast, status )
        if (status/=0) stop 'ERROR in subscal'
        call IntervalQuad_Lin(zpl,zpdu,ppresg(jl-1),ppresg(jl),pdetug(jl), ilast, status )
        if (status/=0) stop 'ERROR in subscal'
        call IntervalQuad_Lin(zpl,zped,ppresg(jl-1),ppresg(jl),pentdg(jl), ilast, status )
        if (status/=0) stop 'ERROR in subscal'
        call IntervalQuad_Lin(zpl,zpdd,ppresg(jl-1),ppresg(jl),pdetdg(jl), ilast, status )
        if (status/=0) stop 'ERROR in subscal'
      end do

    else

      pentug = zpeu
      pdetug = zpdu
      pentdg = zped
      pdetdg = zpdd

    endif

!    ! * Interpolate diffusion coefficients from upper boundaries of ECMWF layers to
!    !   upper boundaries of TM3 layers and
!    !   reorganize the vertical index of the diffusion coefficient
!    !   pdkg(1) is diffusion coefficient at the top of layer 1 (TOA)
!    !   pdkg(lme) is diffusion coefficient at the top of layer lme
!
!    call interp_muherm( ppress, zpdk, ppresg, pdkg )
!    pdkg(1)=0.
!
!    ! * Interpolate geopotential height from ECMWF layer boundaries to TM3 layer
!    !  boundaries
!
!    call interp_muherm( ppress, pz, ppresg, pzg )
!    pzg(0)=pz(0)
!    pzg(lme)=pz(lme)

    !--------------------------------------------------------------------------
    ! correct massfluxes on TM3 grid to conserve mass
    !
    ! updraft

    ! correction level is uppermost level with nonzero detrainment
    lct=lme
    do jl=1,lme
      if ( pdetug(jl) > 0.0 ) then
        lct=jl
        exit
      endif
    end do

    zam=0.
    do jl=lme,1,-1
      zam=zam+pentug(jl)-pdetug(jl)
    end do

    pdetug(lct)=pdetug(lct)+zam
    ! downdraft

    ! correction level is lowermost level

    lct=lme
    zam=0.
    do jl=lme,1,-1
      zam=zam+pentdg(jl)-pdetdg(jl)
    end do


    pdetdg(lct)=pdetdg(lct)+zam
    !
    ! check conservation
    !
    sumup=0.
    sumdown=0.
    do jl=1,lme
      sumup=sumup+pentug(jl)-pdetug(jl)
      sumdown=sumdown+pentdg(jl)-pdetdg(jl)
    enddo
    if ( abs(sumup) > 1.0e-5 ) then
      write(*,*)' ERROR no massconserv in updraft'
      stop
    endif
    if ( abs(sumdown) > 1.0e-5 ) then
      write(*,*)' ERROR no massconserv in downdraft'
      stop
    endif

!!pentug = 0.0
!!pdetug = 0.0
!!pentdg = 0.0
!!pdetdg = 0.0

  end subroutine subscal


  ! =================================================


      !***********************************************************************
      !****     cloud   - routine to calculate cloud properties
      !
      ! m.heimann       mpi HH  Nov-12-1990
      !
      ! Purpose
      ! -------
      ! cloud calculates all properties associated with subgridscale
      ! cloud mixing, i.e. massflux in updraft and downdraft, entrainment
      ! and detrainment rates per level, and cloud properties: temperature,
      ! moisture and liquid water and precipitation rate on each level
      ! boundary.
      !
      !**       interface
      ! ---------
      !
      !       call cloud(pz,ppress,pt,pq,pw,pqac,pqam,
      !     .            pqtur,pamu,ptc,pqc,plc,pgp,
      !     .            pam,peu,pdu,ped,pdd,
      !     .            klc,klt)
      !
      ! input: (all variables are defined on level boundaries)
      !         pz      geopotential height (m)
      !         ppress  pressure        (Pa)
      !         pt      temperature     (K)
      !         pq      moisture        (kg water/kg air)
      !         pw      vertcal velocity (negative upward)   (Pa s**-1)
      !         pqac    div(q v)        (kg water/kg air s**-1)
      !         pqam    div(v)          (s**-1)
      !         pqtur   Fq surf-air     (kg water/m**2 s**-1)
      !
      ! output: variables defined on level boundaries:
      !         pamu    massflux                (kg m**-2 s**-1)
      !         ptc     updraft temperature     (K)
      !         pqc     updraft moisture        (kg water/kg air)
      !         plc     updraft liquid water    (kg water/kg air)
      !         pgp     updraft precipitation   (kg water m**-2 s**-1)
      !
      !         variables defined on each model level:
      !         pam     air mass                        (kg m**-2)
      !         peu     entrainment in updraft          (kg m**-2 s**-1)
      !         pdu     detrainment in updraft          (same)
      !         ped     entrainment in downdraft        (kg m**-2 s**-1)
      !         pdd     detrainment in downdraft        (same)
      !         klc     lowest level in cloud
      !         klt     first level above cloud
      !
      ! method
      ! ------
      !
      ! updraft cloud properties are calculated according to Tiedke scheme
      !
      ! externals
      ! ---------
      !
      ! needs functions qsat and dqsat_dt
      !
      ! references
      ! ----------
      !
      ! Tiedke, Mon. Wea. Rev., 117, 1779-1800, 1989.
      !
      ! revisions
      ! ---------
      !   26-jan-1995 , sr
      !
      !   including now cumulus downdraft and midlevel convection.
      !   Calculates corresponding entrainment and detrainment rates
      !
      !-------------------------------------------------------------------------

      subroutine cloud(lme, pz,ppress,pt,pq,pw,pqac,pqam, pqtur,&
                       pamu,ptc,pqc,plc,pgp, &
                       pam,peu,pdu,ped,pdd)

        use phys_humidity, only : QSat, dQSat_dt

        ! --- flags -----------------------------------------

        ! parameter llcudo=.true.  : cumulus downdraft turned on
        ! parameter llcudo=.false. : cumulus downdraft turned off
        ! parameter llmilc=.true.  : midlevel convection turned on
        ! parameter llmilc=.false. : midlevel convection turned off

        logical       :: llcudo = .true.
        logical       :: llmilc = .true.

        ! --- const (dims) -------------------------------------

        !! vertical resolution (no of model layers)
        !integer   ::  jplm = lme
        !integer   ::  jplmm1 = jplm-1

        ! --- in/out ----------------------------------------

        integer, intent(in) ::  lme

        real, intent(in)  :: pz(0:lme)
        real, intent(in)  :: pt(0:lme)
        real, intent(in)  :: pq(0:lme)
        real, intent(in)  :: pw(0:lme)
        real, intent(in)  :: ppress(0:lme)
        real, intent(in)  :: pqac(0:lme),pqam(0:lme)
        real, intent(in)  :: pqtur

        real, intent(out) :: pam(lme),peu(lme),pdu(lme),ped(lme),pdd(lme)
        real, intent(out) :: pamu(0:lme),ptc(0:lme),pqc(0:lme),plc(0:lme),pgp(0:lme)

        ! --- const (phys) -------------------------------------

        ! physical constants

        real, parameter  ::  pprgasd = 287.05
        real, parameter  ::  pprgasv = 461.51
        real, parameter  ::  ppeps = pprgasd/pprgasv

        real, parameter  ::  ppg = 9.80665
        real, parameter  ::  ppcpd = 1005.46
        real, parameter  ::  ppalv = 2.5008e6
        real, parameter  ::  ppzeta = ppalv/ppcpd
        real, parameter  ::  ppvtcf = (1.0-ppeps)/ppeps

        ! * constants for turbulent entrainment and detrainment rates

        ! penetrative and midlevel convection
        real, parameter   ::  ppepsu = 1.e-4
        real, parameter   ::  ppdeltu = 1.e-4

        ! shallow convection
        real, parameter   ::  ppepsus = 3.e-4
        real, parameter   ::  ppdeltus = 3.e-4

        ! downdraft
        real, parameter   ::  ppepsd = 2.e-4
        real, parameter   ::  ppdeltd = 2.e-4

        ! * constants for precipitation parametrization
        real, parameter   ::  ppkmin = 1500.0
        real, parameter   ::  ppkval = 2.e-3

        ! * parameter beta determines the overshoot of the detrainment plumes
        ! above level of no buoyancy (beta=0 : no overshoot)

        ! penetrative and midlevel convection
        real, parameter   ::  ppbeta = 0.0
        ! shallow convection
        real, parameter   ::  ppbetas = 0.15

        ! * parameter gamma determines downward massflux at the level of free
        ! sinking (LFS)

        real, parameter   ::  ppgamma = 0.3

        ! parameter alpha determines specific humidity of air parcel at surface
        ! before starting the dry adiabatic ascent
        ! (alpha = 0 : air parcel has the specific humidity of the environment,
        !  alpha = 1 : air parcel has saturation specific humidity of the
        !              environment)
        real, parameter   ::  ppalpha = 0.2

        ! no of iterations for saturation calculation
        integer, parameter   ::   jpitermax = 5

        ! --- local --------------------------------------

        integer :: klc,klt,klfs

        real :: zfck,zamub,zamdlfs
        real :: zamu,zamd
        real :: zlc,zqc,ztc
        real :: zlcklc,zqcklc,ztcklc
        real :: ztd,zqd

        real :: zdq1,zdq2,zsv,zscv,zgp
        real :: zpgp(lme)

        real :: zepsu,zdeltu
        real :: zbeta
        integer :: jl,jjl,jiter
        integer :: imlc

        real  :: zdqcmin,zdqdmin
        real  :: ztenwb,zqenwb
        real  :: zttest,zqtest,zstv

        ! extra
        character(len=10) :: convection_type

        ! --- begin ----------------------------------------------

        convection_type = 'none'

        ! default cloud properties on each layer boundary

        do jl = 0, lme
          pamu(jl) = 0.0
          ptc(jl) = pt(jl)
          pqc(jl) = pq(jl)
          plc(jl) = 0.0
          pgp(jl) = 0.0
        end do

        ! air-masses (kg/m**2) in each layer, default entrainment/detrainment
        ! and precipitation rates

        do jl = 1, lme
          pam(jl)=(ppress(jl)-ppress(jl-1))/ppg
          peu(jl)=0.
          pdu(jl)=0.
          ped(jl)=0.
          pdd(jl)=0.
          zpgp(jl)=0.
        end do

        ! find condensation level by lifting an air parcel until supersaturation occurs

        ! we start the ascent with moisture and temperature of the upper boundary of
        ! the surface layer

        ztc = pt(lme-1)
        zqc = pq(lme-1) + ppalpha * ( qsat(pt(lme-1),ppress(lme-1)) - pq(lme-1) )

        do jl = lme-1-1, 1, -1

          ! preliminary value of parcel temperature (dry adiabatic ascent)
          ztc = ztc - ppg*(pz(jl)-pz(jl+1))/ppcpd

          ! check for supersaturation
          if ( zqc > qsat(ztc,ppress(jl)) ) then

            ! if supersaturation is detected we adjust moisture and
            ! temperature by condensation
            ! and set liquid water content equal to the condensate

            zdq2 = 0.0
            do jiter=1,jpitermax
              zdq1=(zqc-QSat(ztc,ppress(jl))) &
                   /(1.+ppzeta*dQSat_dt(ztc,ppress(jl)))
              zqc=zqc-zdq1
              ztc=ztc+zdq1*ppzeta
              zdq2=zdq2+zdq1
            end do
            zlc = zdq2

            ! check if parcel is buoyant:
            !   virtual temperature of parcel
            zscv = ztc*( 1.0 + ppvtcf*zqc - zlc )

            ! virtual temperature of environment
            zsv = pt(jl) * ( 1.0 + ppvtcf*pq(jl) )

            if ( zscv > zsv ) then
              ! if parcel is still buoyant then we have detected the cloud base
              klc=jl
              goto 20
            else
              !  if parcel is not buoyant then there is no penetrative or shallow
              !  convection
              if (llmilc) then
                goto 700
              else
                goto 3000
              endif
            endif

          endif

        end do

        ! no condensation level found
        goto 3000

20      continue

        ! if we arrive here a cloud base is detected:
        ! klc is cloud base level, ztc is cloud temperature, zqc moisture (at
        ! saturation)
        ! and zlc the liquid water content in the updraft at base level

        ztcklc = ztc
        zqcklc = zqc
        zlcklc = zlc

        ! calculate large scale moisture convergence below cloud base
        ! (use trapezoidal integration formula)

        !zfck=
        !    ((pq(klc)*pqam(klc)-pqac(klc))*pam(klc)+
        ! Correction Zoe Stockwell - Peter van Velthoven, 5 January 2000 &
        zfck = ((pq(klc)*pqam(klc) -pqac(klc) )*pam(klc+1)+ &
                (pq(klc)*pqam(lme)-pqac(lme))*pam(lme)      )

        do jl=klc+1,lme-1
          zfck=zfck+(pq(klc)*pqam(jl)-pqac(jl))*(pam(jl)+pam(jl+1))
        end do

        ! check for shallow or penetrative convection, set proper parameter
        ! values

        if (zfck.gt.0.) then
          ! penetrative and midlevel convection
          zepsu=ppepsu
          zdeltu=ppdeltu
          zbeta=ppbeta
          convection_type = 'deep'
        else
          ! shallow convection
          zepsu=ppepsus
          zdeltu=ppdeltus
          zbeta=ppbetas
          convection_type = 'shallow'
        endif

        zfck=pqtur+0.5*zfck


        ! if no moisture convergence then no penetrative or shallow
        ! convection is present

        if (zfck <= 0.0 ) then
          if (llmilc) then
             goto 700
          else
             goto 3000
          endif
        else
          goto 900
        endif

 700    continue

        ! check possibility for midlevel convection

        imlc = 0

        ! Check atmosphere from 2 layers above the surface to the middle of
        ! the atmosphere to see if midlevel convection might occur

        do jl=lme-1-1,lme-int(lme/2.),-1

          ! large scale ascent and an environmental relative humidity of more than
          ! 90% are needed for midlevel convection to occur

          if ( (pq(jl).gt.(0.9*qsat(pt(jl),ppress(jl))) ).and. &
              (pw(jl).lt.0.)) then

            if (imlc.eq.0) then
              ztc = pt(jl+1)
              zqc = pq(jl+1)
              zlc = 0.
              imlc = jl
            else if (imlc.gt.0) then
              if ((imlc-jl).eq.1) then
                imlc = jl
                goto 720
              else
                ztc = pt(jl+1)
                zqc = pq(jl+1)
                zlc = 0.
                imlc = jl
              end if
            end if

 720        continue

            ! do adiabatic ascent

            ztc = ztc - ppg*(pz(jl)-pz(jl+1))/ppcpd

            ! check for supersaturation

            if (zqc.gt.qsat(ztc,ppress(jl))) then

              ! if supersaturation is detected we adjust moisture and
              ! temperature by condensation
              ! and set liquid water content equal to the condensate

              zdq2=0.
              do jiter=1,jpitermax
                zdq1 = (zqc-qsat(ztc,ppress(jl)))/(1.+ppzeta*dqsat_dt(ztc,ppress(jl)))
                zqc=zqc-zdq1
                ztc=ztc+zdq1*ppzeta
                zdq2=zdq2+zdq1
              end do
              zlc=zdq2

            endif

            ! check if parcel is buoyant

            ! virtual temperature of parcel

            zscv = ztc*(1.+ppvtcf*zqc-zlc)

            ! virtual temperature of environment

            zsv = pt(jl)*(1.+ppvtcf*pq(jl))

            if (zscv.gt.zsv) then

              ! parcel is buoyant and we have detected the cloud base of midlevel
              ! convection

              klc = jl
              zamub = -pw(klc)/ppg
              zepsu = ppepsu
              zdeltu = ppdeltu
              zbeta = ppbeta
              llcudo = .false.
              convection_type = 'mid-level'
              goto 1000

            endif

          endif

        end do

        goto 3000

 900    continue

        ! massflux at base of cloud

        ! limit specific humidity difference between cloud and environment at
        ! cloud base

        zdqcmin = max(0.01*pq(klc),1.e-10)
        zdqcmin = max(zdqcmin,zqc+zlc-pq(klc))

        zamub=zfck/zdqcmin

 1000   continue

        ! limit mass flux at cloud base
        zamub=min(zamub,1.0)

        ! set updraft entrainment rates below cloud base proportional
        ! to layer air masses
        ! set updraft detrainment rates below cloud base to zero

        do jl=lme,klc+1,-1
          peu(jl) = zamub*pam(jl)*ppg/(ppress(lme)-ppress(klc))
          pdu(jl) = 0.0
        end do

        ! calculate now parcel ascent within cloud updraft
        !    cloud mass flux               zamu,
        !    cloud moisture                zqc,
        !    cloud temperature             ztc,
        !    cloud liquid water            zlc

        zamu = zamub

        do jl = klc,2,-1

          ! mass entrainment and detrainment
          peu(jl)=zepsu*zamu*(pz(jl-1)-pz(jl))
              pdu(jl)=zdeltu*zamu*(pz(jl-1)-pz(jl))

          ! preliminary values of cloud temperature, moisture and cloud liquid water
          ztc=ztc &
                -ppg*(pz(jl-1)-pz(jl))/ppcpd &
                +zepsu*(pz(jl-1)-pz(jl))*(pt(jl)-ztc)

              zqc=zqc &
                +zepsu*(pz(jl-1)-pz(jl))*(pq(jl)-zqc)

          zlc=zlc &
                +zepsu*(pz(jl-1)-pz(jl))*(-zlc)

          ! adjust moisture and temperature by condensation

          zdq2=0.
          do jiter=1,jpitermax
            zdq1=(zqc-qsat(ztc,ppress(jl))) &
                   /(1.+ppzeta*dqsat_dt(ztc,ppress(jl)))
            zqc=zqc-zdq1
            ztc=ztc+zdq1*ppzeta
            zdq2=zdq2+zdq1
          end do

          ! precipitation rate out of layer jl

          if ((pz(jl)+pz(jl-1))*0.5-pz(klc) .gt. ppkmin) then
            zgp=pam(jl)*ppkval/zamu
          else
            zgp=0.
          endif

          ! adjust liquid cloud water in updraft (use implicit scheme to prevent
          ! instability)

          zgp = zgp*zlc/(1+zgp)

          zpgp(jl) = zgp*zamu

          zlc = zlc-zgp+zdq2

          ! check for level of free sinking (LFS) where cumulus downdraft starts

          if (.not.llcudo) then

            ! downdraft calculation already done or turned off

            goto 800

          end if

          if ( zpgp(jl) == 0.0 ) then

            ! no downdraft exists since downdrafts are associated with convective
            ! precipitation from the updrafts

            goto 800

          end if

          if (jl.lt.3) then

            ! no downdraft since maximum possible cloud height is reached

            goto 800

          end if

          ! The LFS is assumed to be the highest model level where a mixture of equal
          ! parts of cloud air and environmental air (at wet bulb temperature) becomes
          ! negative buoyant with respect to the environmental air

          ! first step :
          ! calculate wet bulb temperature and moisture of the environmental air

          ztenwb = pt(jl-1)
          zqenwb = pq(jl-1)

          ! adjust temperature and moisture by evaporation
          ! zdq1 must be less equal 0 (zdq1=0 : environmental air is saturated,
          ! i.e. zqenwb = pq)

          do jiter = 1,jpitermax
            zdq1 = (zqenwb-qsat(ztenwb,ppress(jl-1)))/ &
                   (1.+ppzeta*dqsat_dt(ztenwb,ppress(jl-1)))
            zqenwb = zqenwb-zdq1
            ztenwb = ztenwb+zdq1*ppzeta
          end do

          ! second step :
          ! do mixing with cloud air

          zttest = 0.5*(ztc+ztenwb)
          zqtest = 0.5*(zqc+zqenwb)

          ! third step :
          ! check for negative buoyancy

          ! virtual temperature of the air mixture

          zstv = zttest*(1.+ppvtcf*zqtest)

          ! virtual temperature of the environmental air

          zsv = pt(jl-1)*(1.+ppvtcf*pq(jl-1))

          if (zstv.lt.zsv) then

            ! level of free sinking (LFS) is found, start downdraft calculation

            ! massflux at LFS is assumed to be directly proportional to the (preliminary)
            ! upward massflux at cloud base

            klfs = jl
            zamdlfs = -ppgamma*zamub
            zamd = zamdlfs
            ztd = zttest
            zqd = zqtest

            ped(klfs) = (-zamd)
            pdd(klfs) = 0.

            if (klfs.eq.klc) goto 45

            do jjl = klfs+1,klc,1

              ! mass entrainment and detrainment

              ped(jjl) = ppepsd*zamd*(pz(jjl)-pz(jjl-1))
              pdd(jjl) = ppdeltd*zamd*(pz(jjl)-pz(jjl-1))

              ! preliminary values of cloud temperature and moisture in downdraft

              ztd = ztd &
                      -ppg*(pz(jjl)-pz(jjl-1))/ppcpd &
                      +ppepsd*(pz(jjl)-pz(jjl-1))*(ztd-pt(jjl-1))

              zqd = zqd &
                      +ppepsd*(pz(jjl)-pz(jjl-1))*(zqd-pq(jjl-1))

              ! adjust moisture and temperature by evaporation

              do jiter=1,jpitermax
                zdq1 = (zqd-qsat(ztd,ppress(jjl-1)))/ &
                       (1.+ppzeta*dqsat_dt(ztd,ppress(jjl-1)))
                zqd = zqd-zdq1
                ztd = ztd+zdq1*ppzeta
              end do

              ! downdraft massflux at lower layer boundary

              zamd = zamd - ped(jjl) + pdd(jjl)

            end do
 45         continue

            zamd = min(0.,zamd)

            ! set downdraft detrainment rates below cloud base proportional to layer
            ! air masses
            ! set downdraft entrainment rates below cloud base to zero

            do jjl = lme,klc+1,-1
              ped(jjl) = 0.
              pdd(jjl) = zamd*pam(jjl)*ppg/(ppress(klc)-ppress(lme))
            end do

            ! recalculate updraft massflux at cloud base,
            ! now allowing for the downdraft massflux

            if (zamd.lt.0.) then
              zdqdmin = zqd-pq(klc)
              zamub = (zfck-zamd*zdqdmin)/zdqcmin
              if (zamub.le.0.) then
                do jjl=1,lme
                  peu(jjl)=0.
                  pdu(jjl)=0.
                  ped(jjl)=0.
                  pdd(jjl)=0.
                end do
                goto 3000
              endif

              ! go back to cloud base and start updraft calculation again

              ztc = ztcklc
              zqc = zqcklc
              zlc = zlcklc
              llcudo = .false.
              goto 1000

            else
              goto 800
            endif

          else

            goto 800

          endif

 800      continue

          ! check for buoyancy (in updraft)
          !   virtual temperature in updraft at upper boundary of layer jl:
          zscv=ztc*(1.+ppvtcf*zqc-zlc)

          !   virtual temperature outside of cloud
          zsv=pt(jl-1)*(1.+ppvtcf*pq(jl-1))

          if ( zscv <= zsv ) then
            klt=jl
            goto 400
          endif

          ! updraft massflux at upper layer boundary
          zamu=zamu+peu(jl)-pdu(jl)

          ! store cloud properties on upper layer boundary
          ptc(jl-1)=ztc
          pqc(jl-1)=zqc
          plc(jl-1)=zlc

        end do

        klt=2

400     continue

        ! set detrainment in two layers above cloud

        pdu(klt-1)=zbeta*zamu
        peu(klt-1)=0.
        pdu(klt)=(1-zbeta)*zamu
        peu(klt)=0.

        ! add up rainrate on each of the layer boundaries

        do jl=klt+1,lme
          pgp(jl)=pgp(jl-1)+zpgp(jl)
        end do

        ! determine net mass flux on each of the layer boundaries

        do jl=lme,1,-1
          pamu(jl-1)=pamu(jl)+(peu(jl)-pdu(jl))-(ped(jl)-pdd(jl))
        end do

        llcudo = .true.
        llmilc = .true.

        return

        ! no cloud present, set cloud base and top to 0 and return

3000    continue
        klc=0
        klt=0
        llcudo = .true.
        llmilc = .true.
        convection_type = 'none'

      end subroutine cloud


  ! ==========================================================
  ! ===
  ! ===  TMPP/SOURCE/tmpp_sub_work.f90
  ! ===
  ! ==========================================================


  !---------------------------------------------------------------------
  ! calculate en/detrainment rates and diffusion coefficient on TM grid
  !---------------------------------------------------------------------
  ! History:
  ! Increased vertical dimension of z,t,q,u,v,w from lme to lme + 1
  ! in order to be able to use the same memory location in worksub
  ! for u and wu, for t and wt, etc.
  ! Added subroutine cen2bound for the same reason
  ! Removed dummy fields for geopotential height and zonal means
  ! Program just fits into memory of SGI machines (max stacksize = 524288) if
  !  TM and ECMWF both have 1x1 degree resolution and 60 levels
  !                                                         pvv, 5-2-2000
  !---------------------------------------------------------------------

  ! =========================================================

  ! interpolate variables from the center of parent-model layers to the
  ! boundaries of parent-model layers and save result in same memory location


  subroutine mid2bound_uv( lme, npe, u, v, ps, mask, a, b )

    ! --- in/out ----------------------------------

    integer, intent(in)             ::  lme, npe
    real, intent(inout)             ::  u(npe,0:lme)
    real, intent(inout)             ::  v(npe,0:lme)
    real, intent(in)                ::  ps(npe)
    logical, intent(in)             ::  mask(npe)
    real, intent(in)                ::  a(0:lme)
    real, intent(in)                ::  b(0:lme)

    ! --- local ---------------------------------

    integer              ::  i
    real                 ::  wcol(0:lme)

    ! --- begin -------------------------

    do i = 1, npe
      if ( mask(i) ) then

        call cen2bound_col( lme, u(i,1:lme), ps(i), 1, wcol, a, b )
        u(i,:) = wcol

        call cen2bound_col( lme, v(i,1:lme), ps(i), 1, wcol, a, b )
        v(i,:) = wcol

      end if
    end do

  end subroutine mid2bound_uv


  ! ===


  subroutine mid2bound_w( lme, npe, w, ps, mask, a, b )

    ! --- in/out ----------------------------------

    integer, intent(in)             ::  lme, npe
    real, intent(inout)             ::  w(npe,0:lme)
    real, intent(in)                ::  ps(npe)
    logical, intent(in)             ::  mask(npe)
    real, intent(in)                ::  a(0:lme)
    real, intent(in)                ::  b(0:lme)

    ! --- local ---------------------------------

    integer              ::  i
    real                 ::  wcol(0:lme)

    ! --- begin -------------------------

    do i = 1, npe
      if ( mask(i) ) then

        call cen2bound_col( lme, w(i,1:lme), ps(i), 1, wcol, a, b )
        w(i,:) = wcol

        ! set to zero at top
        w(i,0) = 0.0

      end if
    end do

  end subroutine mid2bound_w


  ! ===


  subroutine mid2bound_t( lme, npe, t, ps, mask, a, b )

    ! --- in/out ----------------------------------

    integer, intent(in)             ::  lme, npe
    real, intent(inout)             ::  t(npe,0:lme)
    real, intent(in)                ::  ps(npe)
    logical, intent(in)             ::  mask(npe)
    real, intent(in)                ::  a(0:lme)
    real, intent(in)                ::  b(0:lme)

    ! --- local ---------------------------------

    integer              ::  i
    real                 ::  wcol(0:lme)

    ! --- begin -------------------------

    do i = 1, npe
      if ( mask(i) ) then

        call cen2bound_col( lme, t(i,1:lme), ps(i), 1, wcol, a, b )
        t(i,:) = wcol

      end if
    end do

  end subroutine mid2bound_t


  ! ===


  subroutine mid2bound_q( lme, npe, q, ps, mask, a, b, t )

    use phys_humidity, only : qsat

    ! --- in/out ----------------------------------

    integer, intent(in)             ::  lme, npe
    real, intent(inout)             ::  q(npe,0:lme)
    real, intent(in)                ::  ps(npe)
    logical, intent(in)             ::  mask(npe)
    real, intent(in)                ::  a(0:lme)
    real, intent(in)                ::  b(0:lme)
    real, intent(in)                ::  t(npe,0:lme)

    ! --- local ---------------------------------

    integer              ::  i, l
    real                 ::  wcol(0:lme)
    real            :: tmpress(0:lme)

    ! --- begin -------------------------

    do i = 1, npe
      if ( mask(i) ) then

        call cen2bound_col( lme, q(i,1:lme), ps(i), 1, wcol, a, b )
        q(i,:) = wcol

        ! limit specific humidity at 0 and qsat(t,p)
        ! first establish hybrid vertical coordinate at i,j ;
        ! note that ps is expressed in Pa
        do l = 0, lme
          tmpress(l) = a(l) + ps(i)*b(l)
        end do
        do l = 0, lme
          q(i,l) = min( qsat(t(i,l),tmpress(l)), max(0.0,q(i,l)) )
        end do

      end if
    end do

  end subroutine mid2bound_q


  ! ===


  subroutine mid2bound_z( lme, npe, z, ps, mask, a, b, zsurf )

    use Binas, only : g => grav

    ! --- in/out ----------------------------------

    integer, intent(in)             ::  lme, npe
    real, intent(inout)             ::  z(npe,0:lme)
    real, intent(in)                ::  ps(npe)
    logical, intent(in)             ::  mask(npe)
    real, intent(in)                ::  a(0:lme)
    real, intent(in)                ::  b(0:lme)
    real, intent(in)                ::  zsurf(npe)

    ! --- local ---------------------------------

    integer              ::  i
    real                 ::  wcol(0:lme)

    ! --- begin -------------------------

    do i = 1, npe
      if ( mask(i) ) then

        call cen2bound_col( lme, z(i,1:lme), ps(i), 1, wcol, a, b )
        z(i,:) = wcol

        ! set to known value at surface:
        z(i,lme) = zsurf(i)/g

      end if
    end do

  end subroutine mid2bound_z


  ! ===


  subroutine mid2bound_p( lme, npe, p, ps, mask, a, b )

    ! --- in/out ----------------------------------

    integer, intent(in)             ::  lme, npe
    real, intent(inout)             ::  p(npe,0:lme)
    real, intent(in)                ::  ps(npe)
    logical, intent(in)             ::  mask(npe)
    real, intent(in)                ::  a(0:lme)
    real, intent(in)                ::  b(0:lme)

    ! --- local ---------------------------------

    integer              ::  i, j
    real                 ::  wcol(0:lme)

    ! --- begin -------------------------

    do i = 1, npe
      if ( mask(i) ) then

        call cen2bound_col( lme, p(i,:), ps(i), 0, wcol, a, b )
        p(i,:) = wcol

      end if
    end do

  end subroutine mid2bound_p



  ! interpolate from mid-levels of parent model to
  ! level boundaries
  !
  ! Peter van Velthoven  - 4 January 2000
  ! This subroutine is included to save memory
  ! field and field2 use the same space in memory
  !
  ! iopt = 0    : no field as input : fill wfield with pressure
  ! iopt = other: use field as input
  !

  subroutine cen2bound_col( lme, field, ps, iopt, wfield, a, b )

    use Num, only : interp_muherm

    ! --- in/out -------------------------------

    integer, intent(in)  :: lme
    real, intent(in)     :: field(lme)    ! input on (mid-)levels
    real, intent(in)     :: ps            ! surface pressure
    integer, intent(in)  :: iopt
    real, intent(out)    :: wfield(0:lme) ! output on vertical level boundaries
    real, intent(in)     ::  a(0:lme)
    real, intent(in)     ::  b(0:lme)

    ! --- begin -------------------------------

    integer       ::  status
    real          ::  ztemp(lme)
    real          ::  tmpress(0:lme) ! pressure at ECMWF vertical level boundaries
    real          ::  tcmpress(lme)  ! pressure at ECMWF (mid-)levels
    integer       ::  l

    ! --- begin --------------------------------

    ! establish hybrid vertical coordinate at i,j
    ! note that ps is expressed in Pa
    tmpress = a + ps * b

    ! calculate pressure at model layer center
    do l=1,lme
      tcmpress(l) = (tmpress(l-1)+tmpress(l))/2.
    end do

    if ( iopt == 0 ) then
      wfield = tmpress
    else
      call interp_muherm( tcmpress, field, tmpress, wfield, status )
      if (status/=0) stop 'ERROR in cen2bound_col'
    end if


  end subroutine cen2bound_col



  ! ==========================================================
  ! ===
  ! ===  TMPP/SOURCE/tmpp_geometry.f90
  ! ===
  ! ==========================================================


  ! pressure at half leves from bottom to top

  subroutine phlev_ec_rev( lme, a_ec, b_ec, ps, pb )

    ! --- in/out --------------------------

    integer, intent(in)        ::  lme
    real, intent(in)           ::  a_ec(0:lme)
    real, intent(in)           ::  b_ec(0:lme)
    real, intent(in)           ::  ps
    real, intent(out)          ::  pb(0:lme)

    ! --- local --------------------------

    integer          ::  l

    ! --- in/out -------------------------

    do l = 0, lme
      pb(lme-l) = a_ec(l) + b_ec(l) * ps
    end do

  end subroutine phlev_ec_rev


  ! ==========================================================
  ! ===
  ! ===  TMPP/SOURCE/phys_geopot.f90
  ! ===
  ! ==========================================================


  !
  ! NAME
  !   GeoPot_col  -  calculate geopotential height
  !
  ! DESCRIPTION
  !   Calculate geopotential height from halflevel pressures
  !   and full level virtual temperature.
  !
  ! USAGE
  !
  !   call GeoPot( z, zsurf, pt, pq, pb, lm )
  !
  !     integer, intent(in)  ::  lm        ! number of levels
  !
  !       (levels numbered downwards (top -> down) )
  !
  !     real, intent(out)    ::  z(lm)     ! geopotential height (m ?).
  !     real, intent(in)     ::  zsurf     ! orography (m ?)
  !     real, intent(in)     ::  pt(lm)    ! temperature (K ?)
  !     real, intent(in)     ::  pq(lm)    ! specific humidity (??)
  !
  !       (levels numbered upwards (bottom -> up) )
  !
  !     real, intent(in)     ::  pb(0:lm)  ! pressure at half levels
  !
  ! HISTORY
  !
  !   06-11-2001, Arjo Segers
  !     Extracted from original routines 'geopot' and 'auxhyb'
  !     by Ad Jeuken
  !

  subroutine GeoPot( lm, zsurf, pt, pq, pb, z )

    ! --- in/out -------------------------

    integer, intent(in)  ::  lm
    real, intent(out)    ::  z(lm)
    real, intent(in)     ::  zsurf
    real, intent(in)     ::  pt(lm)
    real, intent(in)     ::  pq(lm)
    real, intent(in)     ::  pb(0:lm)

    ! --- const ------------------------

    real, parameter      ::  rd = 287.05
    real, parameter      ::  g0 = 9.801

    ! --- local ------------------------------

    integer  ::  linv
    real     ::  pdelp, prdelp

    real     ::  palfa(lm)
    real     ::  plnr(lm)

    real     ::  tv(lm)

    integer  ::  l, lp1

    ! --- begin ---------------------------------

    ! >>> former routine 'auxhyb' >>>>>>>>>>>>>>>>>>>>>>

    ! loop from top to bottom:
    do l = 1, lm
      linv = lm - l   ! lm-1, 0
      pdelp = pb(linv) - pb(linv+1)
      prdelp = 1.0 / pdelp
      if ( l == 1 ) then
        plnr(l) = rd * 1.3862944
      else
        plnr(l) = rd * log( pb(linv)/pb(linv+1) )
      end if
      palfa(l) = rd - pb(linv+1) * plnr(l) * prdelp
    end do

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! loop from bottom to top:
    do l = lm, 1, -1
      tv(l) = pt(l) * ( 1.0 + 0.608*pq(l) )
      if ( l == lm ) then
        z(l) = palfa(l)*tv(l)/g0 + zsurf/g0
      else
        lp1=l+1
        z(l) = z(lp1) + ( palfa(l)*tv(l) + (plnr(lp1)-palfa(lp1))*tv(lp1) )/g0
      end if
    end do

  end subroutine GeoPot



end module phys_convec_tmpp
