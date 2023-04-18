!###############################################################################
!
! open and read all HDF data
! perform some meteo_dependend calculations
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module meteo

  use GO        , only : gol, goErr, goPr
  !
  use TMM       , only : TtmMeteo
  !
  use dims, only : nregions
  !
  use meteodata  , only : mdat_set, mdat_check
  use meteodata  , only : mdat_timeinterpolation
  !
  use meteodata  , only : sp1_dat, sp2_dat, sp_dat, spm_dat
  use meteodata  , only : tsp_dat
  use meteodata  , only : phlb_dat, m_dat
  use meteodata  , only : mfu_dat, mfv_dat, mfw_dat
  use meteodata  , only : pu_dat,  pv_dat,  pw_dat
  use meteodata  , only : temper_dat, humid_dat
  use meteodata  , only : gph_dat
  use meteodata  , only : omega_dat
  use meteodata  , only : lwc_dat, iwc_dat, cc_dat, cco_dat, ccu_dat
  use meteodata  , only : entu_dat, entd_dat, detu_dat, detd_dat
  use meteodata  , only : oro_dat
  use meteodata  , only : lsmask_dat
  use meteodata  , only : pclim_dat
  use meteodata  , only : sr_ecm_dat, sr_ols_dat
  use meteodata  , only : u10m_dat, v10m_dat
  use meteodata  , only : sshf_dat, slhf_dat
  use meteodata  , only : ewss_dat, nsss_dat
  use meteodata  , only : blh_dat, ci_dat
  ! new for dry_deposition
#ifndef without_dry_deposition
  use meteodata  , only : sd_dat, sf_dat, swvl1_dat, src_dat
  use meteodata  , only : d2m_dat, t2m_dat, ssr_dat, tv_dat, cvl_dat, cvh_dat, nveg
#endif

  implicit none

  ! --- local ----------------------------

  ! Interface to TM meteo data
  type(TtmMeteo), save         ::  tmmd

  ! single cell global surface pressure (region 0)
  real                         ::  sp_region0(1,1)


  ! --- in/out -------------------------------------

  private

  public  ::  Meteo_Def_Timers
  public  ::  Meteo_Init, Meteo_Done, Meteo_Alloc
  public  ::  Meteo_Setup_Mass
  public  ::  Meteo_Setup_Other
  public  ::  mdat_Set, mdat_Check
  public  ::  mdat_TimeInterpolation

  !public  ::  sp1_dat, sp2_dat, sp_dat, spm_dat
  !public  ::  tsp_dat
  !public  ::  phlb_dat, m_dat
  !public  ::  mfu_dat, mfv_dat, mfw_dat
  !public  ::   pu_dat,  pv_dat,  pw_dat
  !public  ::  temper_dat, humid_dat
  !public  ::  gph_dat
  !public  ::  omega_dat
  !public  ::  lwc_dat, iwc_dat, cc_dat, cco_dat, ccu_dat
  !public  ::  entu_dat, entd_dat, detu_dat, detd_dat
  !public  ::  oro_dat
  !public  ::  lsmask_dat
  !public  ::  pclim_dat
  !public  ::  sr_ecm_dat, sr_ols_dat
  !public  ::  u10m_dat, v10m_dat
  !public  ::  sshf_dat, slhf_dat
  !public  ::  ewss_dat, nsss_dat
  !public  ::  blh_dat


  ! --- const ----------------------------------

  ! module name
  character(len=*), parameter  ::  mname = 'Meteo'


  ! --- interfaces -------------------------------

  interface Setup
    module procedure Setup_2d
    module procedure Setup_3d
  end interface

  ! --- local -----------------------------------

  !! timers:
  !integer     ::  itim_read_wind, itim_read_tempq, itim_read_cloud, itim_read_conv


contains


  ! ==========================================================


  subroutine Meteo_Def_Timers( status )

!    use TMM, only : TMM_Def_Timers

    ! --- in/out -------------------------------

    integer, intent(out)         ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Meteo_Def_Timers'

    ! --- begin ----------------------------

    ! define timers:
    !call GO_Timer_Def( itim_read_wind , 'meteo read wind' )
    !call GO_Timer_Def( itim_read_tempq, 'meteo read temperature/humidity' )
    !call GO_Timer_Def( itim_read_cloud, 'meteo read clouds' )
    !call GO_Timer_Def( itim_read_conv , 'meteo read convection' )

!    ! define TMM timers:
!    call TMM_Def_Timers( status )
!    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Meteo_Def_Timers


  ! ==========================================================


  subroutine Meteo_Init( status )

    use Binas      , only : p_global
    use TMM        , only : Init
    use global_data, only : rcF
    use dims       , only : nregions_all, iglbsfc
    use dims       , only : nregions, region_name
    use dims       , only : im, jm
    use dims       , only : lm, lmax_conv
    use meteodata  , only : mdat_init
    use global_data, only : rcfile

    ! --- in/out -------------------------------

    integer, intent(out)         ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Meteo_Init'

    ! --- local -----------------------------

    integer             ::  region, n, iveg
    integer             ::  imr, jmr, lmr
    integer             ::  halo
    character(len=4)    ::  sveg

    ! --- begin ----------------------------

    write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! initialize timers:
    call Meteo_Def_Timers( status )
    IF_NOTOK_RETURN(status=1)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! meteo database
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! setup interface to TM meteo:
    call Init( tmmd, rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! define meteo data
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! global mean surface pressure
    sp_region0 = p_global

    ! setup meteo fields: not in use, not allocated:
    do region = 1, nregions_all

      ! local grid sizes
      if ( (1 <= region) .and. (region <= nregions) ) then
        imr = im(region)
        lmr = lm(region)
        jmr = jm(region)
      else if ( region == iglbsfc ) then
        imr = 360
        jmr = 180
        lmr = lm(1)
      else
        write (gol,'("could not set imr/jmr for region ",i0)') region; call goErr
        TRACEBACK; status=1; return
      end if

      !
      ! ** surface pressure *************************************
      !

      ! two extra horizontal cells
      halo = 2

      ! end of interval; also reads for sp1 and spm :
      call Init_MeteoData( sp2_dat(region), 'sp', 'Pa', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'* ','ml','sp'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! check time interpolation:
      if ( sp2_dat(region)%tinterp(1:6) /= 'interp' ) then
        write (gol,'("surface pressure should be interpolated:")'); call goErr
        write (gol,'("  requested tinterp : ",a)') trim(sp2_dat(region)%tinterp); call goErr
        TRACEBACK; status=1; return
      end if

      ! start of interval (copied from sp2_dat):
      call mdat_init( sp1_dat(region), 'sp', 'Pa', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)

      ! current pressure:
      call mdat_init( sp_dat(region), 'sp', 'Pa', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)

      ! surface pressure at mid of dynamic time interval:
      call mdat_init( spm_dat(region), 'sp', 'Pa', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)


      !
      ! ** 3D pressure and mass **************************
      !

      ! two extra horizontal cells (same as surface pressures)
      halo = 2

      ! pressure at half levels (lm+1):
      call mdat_init( phlb_dat(region), 'phlb', 'Pa', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr+1/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)

      ! air mass:
      call mdat_init( m_dat(region), 'm', 'kg', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)


      !
      ! ** massfluxes *************************************
      !

      ! ~~ vertical

      ! no extra cells
      halo = 0

      ! vertical flux (kg/s)
      call Init_MeteoData( mfw_dat(region), 'mfw', 'kg/s', &
                     (/1,imr/), (/1,jmr/), halo, (/0,lmr/), &
                     rcF, (/'*      ','ml     ','mflux_w'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! vertical flux (kg/s) : BALANCED
      ! NOTE: data is copied from mfw, thus use same tinterp
      !       for correct allocation of data arrays
      call mdat_init( pw_dat(region), 'pw', 'kg/s', mfw_dat(region)%tinterp, &
                    (/1,imr/), (/1,jmr/), halo, (/0,lmr/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)

      ! tendency of surface pressure:
      call Init_MeteoData( tsp_dat(region), 'tsp', 'Pa/s', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*      ','ml     ','mflux_w'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! ~~ horizontal

      ! NOTE: strange old indexing:
      !   pu_tmpp  -->  pu(0:imr,1:jmr  ,1:lmr)  in  pu_t(0:imr+1,0:jmr+1,0:lmr)
      !                                                         ^ ^     ^ ^       too large !
      !   pv_tmpp  -->  pv(1:imr,1:jmr+1,1:lmr)  in  pv_t(0:imr+1,0:jmr+1,0:lmr)
      !                                                   ^     ^ ^       ^       too large !
      ! The extra cells are implemented below as halo cells.

      ! one extra cell
      halo = 1

      !! east/west flux (kg/s)
      !call Init( mfu_dat(region), 'mfu', 'kg/s', tinterp_curr, &
      !               (/0,imr/), (/1,jmr/), halo, (/1,lmr/), &
      !               sourcekey_curr, write_meteo, status, default=destkey_curr )
      !IF_NOTOK_RETURN(status=1)

      !! south/north flux (kg/s)
      !call Init( mfv_dat(region), 'mfv', 'kg/s', tinterp_curr, &
      !               (/1,imr/), (/0,jmr/), halo, (/1,lmr/), &
      !               sourcekey_curr, write_meteo, status, default=destkey_curr )
      !IF_NOTOK_RETURN(status=1)


      ! east/west flux (kg/s)
      call Init_MeteoData( mfu_dat(region), 'mfu', 'kg/s', &
                     (/1,imr/), (/1,jmr/), halo, (/0,lmr/), &
                     rcF, (/'*       ','ml      ','mflux_uv'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! south/north flux (kg/s)
      call Init_MeteoData( mfv_dat(region), 'mfv', 'kg/s', &
                     (/1,imr/), (/1,jmr/), halo, (/0,lmr/), &
                     rcF, (/'*       ','ml      ','mflux_uv'/), region, status )
      IF_NOTOK_RETURN(status=1)


      !! east/west flux (kg/s) : BALANCED
      !call Init( pu_dat(region), 'pu', 'kg/s', 'computed', &
      !               (/0,imr/), (/1,jmr/), halo, (/1,lmr/), &
      !               'no-sourcekey', .false., 'no-destkey', status )
      !IF_NOTOK_RETURN(status=1)
      !
      !! south/north flux (kg/s) : BALANCED
      !call Init( pv_dat(region), 'pv', 'kg/s', 'computed', &
      !               (/1,imr/), (/0,jmr/), halo, (/1,lmr/), &
      !               'no-sourcekey', .false., 'no-destkey', status )
      !IF_NOTOK_RETURN(status=1)

      halo = 1

      ! east/west flux (kg/s) : BALANCED
      ! NOTE: data is copied from mfu, thus use same tinterp
      !       for correct allocation of data arrays
      call mdat_init( pu_dat(region), 'pu', 'kg/s', mfu_dat(region)%tinterp, &
                     (/1,imr/), (/1,jmr/), halo, (/0,lmr/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)

      ! south/north flux (kg/s) : BALANCED
      ! NOTE: data is copied from mfv, thus use same tinterp
      !       for correct allocation of data arrays
      call mdat_init( pv_dat(region), 'pv', 'kg/s', mfv_dat(region)%tinterp, &
                     (/1,imr/), (/1,jmr/), halo, (/0,lmr/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)


      !
      ! ** temperature *************************************
      !

      ! no extra cells
      halo = 0

      ! temperature (K)   (halo=0)
      call Init_MeteoData( temper_dat(region), 'T', 'K', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     rcF, (/'*     ','ml    ','temper'/), region, status )
      IF_NOTOK_RETURN(status=1)


      !
      ! ** humidity *************************************
      !

      ! no extra cells
      halo = 0

      ! humidity (kg/kg)   (halo = 0)
      call Init_MeteoData( humid_dat(region), 'Q', 'kg/kg', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     rcF, (/'*    ','ml   ','humid'/), region, status )
      IF_NOTOK_RETURN(status=1)


      !
      ! ** computed *************************************
      !

      ! no extra cells
      halo = 0

      ! geopotential height(m)  (lm+1, halo=0)
      call mdat_init( gph_dat(region), 'gph', 'm', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr+1/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)

      ! vertical velocity (Pa/s)  (lm+1, halo=0)
      call mdat_init( omega_dat(region), 'omega', 'Pa/s', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr+1/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)


      !
      ! ** clouds *************************************
      !

      ! no extra cells
      halo = 0

      ! lwc    liquid water content (kg/kg) (halo=0)
      call Init_MeteoData( lwc_dat(region), 'CLWC', 'kg/kg', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     rcF, (/'*    ','ml   ','cloud'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! iwc    ice water content (kg/kg) (halo=0)
      call Init_MeteoData( iwc_dat(region), 'CIWC', 'kg/kg', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     rcF, (/'*    ','ml   ','cloud'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! cc     cloud cover (fraction) (halo=0)
      call Init_MeteoData( cc_dat(region), 'CC', '1', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     rcF, (/'*    ','ml   ','cloud'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! cco    cloud cover overhead = above bottom of box (fraction) (halo=0)
      call Init_MeteoData( cco_dat(region), 'CCO', '1', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     rcF, (/'*    ','ml   ','cloud'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! ccu    cloud cover underfeet = below top of box (fraction) (halo=0)
      call Init_MeteoData( ccu_dat(region), 'CCU', '1', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmr/), &
                     rcF, (/'*    ','ml   ','cloud'/), region, status )
      IF_NOTOK_RETURN(status=1)


      !
      ! ** convection *************************************
      !

      ! no extra cells
      halo = 0

      ! entu        entrainement updraft
      call Init_MeteoData( entu_dat(region), 'eu', 'kg/m2/s', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmax_conv/), &
                     rcF, (/'*     ','ml    ','convec'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! entd        entrainement downdraft (im,jm,lmax_conv)
      call Init_MeteoData( entd_dat(region), 'ed', 'kg/m2/s', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmax_conv/), &
                     rcF, (/'*     ','ml    ','convec'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! detu        detrainement updraft
      call Init_MeteoData( detu_dat(region), 'du', 'kg/m2/s', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmax_conv/), &
                     rcF, (/'*     ','ml    ','convec'/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! detd        detrainement downdraft
      call Init_MeteoData( detd_dat(region), 'dd', 'kg/m2/s', &
                     (/1,imr/), (/1,jmr/), halo, (/1,lmax_conv/), &
                     rcF, (/'*     ','ml    ','convec'/), region, status )
      IF_NOTOK_RETURN(status=1)

      !
      ! *** surface fields
      !

      ! no extra cells
      halo = 0

      ! orography (m*[g])
      call Init_MeteoData( oro_dat(region), 'oro', 'm m/s2', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*        ','sfc      ','sfc.const','oro      '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! land/sea mask (%)
      call Init_MeteoData( lsmask_dat(region), 'lsm', '%', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*        ','sfc      ','sfc.const','lsm      '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! climatological surface pressure (Pa) ;
      ! derived from orography:
      call mdat_init( pclim_dat(region), 'pclim', 'Pa', 'computed', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     'no-sourcekey', .false., 'no-destkey', status )
      IF_NOTOK_RETURN(status=1)

      ! ~~~ instantaneous fields

      !! sea surface temperature:
      !call Init_MeteoData( sst_dat(region), 'sst', 'K', &
                     !(/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     !rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','sst     '/), region, status )
      !IF_NOTOK_RETURN(status=1)

      ! 10m u wind (m/s)
      call Init_MeteoData( u10m_dat(region), 'u10m', 'm/s', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','u10m    '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! 10m v wind (m/s)
      call Init_MeteoData( v10m_dat(region), 'v10m', 'm/s', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','v10m    '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! boundary layer height (m)  ; instant
      call Init_MeteoData( blh_dat(region), 'blh', 'm', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','blh     '/), region, status )
      IF_NOTOK_RETURN(status=1)

#ifndef without_dry_deposition
      ! skin reservoir content (m water) ; instant
      call Init_MeteoData( src_dat(region), 'src', 'm', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','src     '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! 2 meter dewpoint temperature (K)  ; instant
      call Init_MeteoData( d2m_dat(region), 'd2m', 'K', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','d2m     '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! 2 meter temperature (K)  ; instant
      call Init_MeteoData( t2m_dat(region), 't2m', 'K', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','t2m     '/), region, status )
      IF_NOTOK_RETURN(status=1)
#endif

      !! skin temperature (K)  ; instant
      !call Init_MeteoData( skt_dat(region), 'skt', 'K', &
                     !(/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     !rcF, (/'*       ','sfc     ','sfc.inst','sfc.fc  ','skt     '/), region, status )
      !IF_NOTOK_RETURN(status=1)

      ! ~~~ average field (accumulated)

      ! east-west surface stress (N/m2); time aver
      call Init_MeteoData( ewss_dat(region), 'ewss', 'N/m2', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.aver','sfc.fc  ','ewss    '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! north-south surface stress (N/m2) ; time aver
      call Init_MeteoData( nsss_dat(region), 'nsss', 'N/m2', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.aver','sfc.fc  ','nsss    '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! surface roughness (ecmwf,ncep)
      call Init_MeteoData( sr_ecm_dat(region), 'sr', 'm', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.an  ','sr      '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! surface sensible heat flux (W/m2) ; time aver
      call Init_MeteoData( sshf_dat(region), 'sshf', 'W/m2', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.aver','sfc.fc  ','sshf    '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! surface latent heat flux (W/m2) ; time aver
      call Init_MeteoData( slhf_dat(region), 'slhf', 'W/m2', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.aver','sfc.fc  ','slhf    '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! sea ice is also used for C14 disequilibrium flux, not just dry deposition
      call Init_MeteoData( ci_dat(region), 'ci', '1', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.fc  ','ci      '/), region, status )
      IF_NOTOK_RETURN(status=1)

#ifndef without_dry_deposition

      ! snow fall (m water eqv); time aver
      call Init_MeteoData( sf_dat(region), 'sf', 'm', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.aver','sfc.fc  ','sf      '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! snow depth (m water eqv); day aver ?
      call Init_MeteoData( sd_dat(region), 'sd', 'm', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.fc  ','sd      '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! volumetric soil water layer 1 ( m3 water / m3 soil) ; day aver ?
      call Init_MeteoData( swvl1_dat(region), 'swvl1', '1', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.fc  ','swvl1   '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! vegetation type (%) ; day aver
      do iveg = 1, nveg
        write (sveg,'("tv",i2.2)') iveg
        call Init_MeteoData( tv_dat(region,iveg), sveg, '%', &
                       (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                       rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.an  ','veg     '/), region, status )
        IF_NOTOK_RETURN(status=1)
      end do

      ! low vegetation cover (0-1) ; day aver
      call Init_MeteoData( cvl_dat(region), 'cvl', '1', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.an  ','veg     '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! high vegetation cover (0-1) ; day aver
      call Init_MeteoData( cvh_dat(region), 'cvh', '1', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.an  ','veg     '/), region, status )
      IF_NOTOK_RETURN(status=1)

      ! surface solar radiation ( W/m2 ) ; time aver
      call Init_MeteoData( ssr_dat(region), 'ssr', 'W/m2', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.aver','sfc.fc  ','ssr     '/), region, status )
      IF_NOTOK_RETURN(status=1)
#endif

      ! ~~~ monthly data

      ! surface roughness (olsson) ; monthly
      call Init_MeteoData( sr_ols_dat(region), 'srols', 'm', &
                     (/1,imr/), (/1,jmr/), halo, (/1,1/), &
                     rcF, (/'*       ','sfc     ','sfc.inst','sfc.day ','sfc.an  ','srols   '/), region, status )
      IF_NOTOK_RETURN(status=1)

    end do  ! regions

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! done
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Meteo_Init


  ! ***


  !
  ! Read multiple keys in rcfile to setup meteodata structure.
  ! The following keys are read:
  !
  !   meteo.tinterp.<param>              <-- time interpolation
  !   tmm.sourcekey.<grid>.<param>       <-- input file name description
  !   tmm.output.<grid>.<param>          <-- write meteo ?
  !   tmm.destkey.<grid>.<param>         <-- output file name description
  !
  ! where <grid> is first '*' and then set to the grid name,
  ! and <param> is set to each of the provided keys.
  !
  subroutine Init_MeteoData( md, name, unit, is, js, halo, ls, &
                                 rcF, rcs, region, status )

    use GO          , only : TRcFile, ReadRc
    use Dims        , only : nregions, nregions_max, okdebug_tmm
    use TM5_Geometry, only : lli
    use MeteoData   , only : TMeteoData, mdat_init, mdat_set
    use os_specs    , only : MAX_RCKEY_LEN

    ! --- in/out -------------------------------------

    type(TMeteoData), intent(out)         ::  md
    !! SB: Sometimes we need to read in met fields for emission routines, such as u10m for C14 disequilibrium, or ssr for the
    !! Olsen & Randerson formula. However, emission_init is called before this. So even if we set u10m_dat(:)%used = .true. in
    !! emission_init, it gets reset to .false. in this routine. A workaround -- and Arjo can say if it will have other side
    !! effects -- is to make md intent(inout) instead of intent(out), and not modify the 'used' component if it's already set
    !! to .true.
    !type(TMeteoData), intent(inout)       ::  md
    character(len=*), intent(in)          ::  name, unit
    integer, intent(in)                   ::  is(2), js(2)
    integer, intent(in)                   ::  halo
    integer, intent(in)                   ::  ls(2)
    type(TRcFile), intent(inout)          ::  rcF
    character(len=*), intent(in)          ::  rcs(:)
    integer, intent(in)                   ::  region
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Init_MeteoData'

    ! --- local -------------------------------------

    character(len=10)   ::  tinterp
    character(len=MAX_RCKEY_LEN)  ::  sourcekey
    logical             ::  write_meteo
    character(len=MAX_RCKEY_LEN)  ::  destkey
    logical             ::  used
    integer             ::  mi

    ! --- begin -------------------------------------

    ! time interpolation :
    call ReadRc( rcF, 'meteo.tinterp', rcs, tinterp, status )
    IF_NOTOK_RETURN(status=1)

    ! source filenames:
    call ReadRc( rcF, 'tmm.sourcekey.*', rcs, sourcekey, status, default='no-sourcekey', match_index=mi )
    IF_ERROR_RETURN(status=1)
    if (mi > 0 .and. okdebug_tmm) then
        write(*,'(a, " :: for field ", a, ", ", a, " = ", a)') rname, trim(name), &
        'tmm.sourcekey.*'//'.'//trim(rcs(mi)), trim(sourcekey)
    end if

    call ReadRc( rcF, 'tmm.sourcekey.'//trim(lli(region)%name), rcs, sourcekey, status, default=sourcekey, match_index=mi )
    IF_ERROR_RETURN(status=1)
    if (mi > 0 .and. okdebug_tmm) then
        write(*,'(a, " :: for field ", a, ", ", a, " = ", a)') rname, trim(name), &
        'tmm.sourcekey.'//trim(lli(region)%name)//'.'//trim(rcs(mi)), trim(sourcekey)
    end if

    ! write flag:
    call ReadRc( rcF, 'tmm.output.*', rcs, write_meteo, status, default=.false., match_index=mi )
    IF_ERROR_RETURN(status=1)
    if (mi > 0 .and. okdebug_tmm) then
        write(*,'(a, " :: for field ", a, ", ", a, " = ", l1)') rname, trim(name), &
        'tmm.output.*'//'.'//trim(rcs(mi)), write_meteo
    end if
    call ReadRc( rcF, 'tmm.output.'//trim(lli(region)%name), rcs, write_meteo, status, default=write_meteo, match_index=mi )
    IF_ERROR_RETURN(status=1)
    if (mi > 0 .and. okdebug_tmm) then
        write(*,'(a, " :: for field ", a, ", ", a, " = ", l1)') rname, trim(name), &
        'tmm.output.'//trim(lli(region)%name)//'.'//trim(rcs(mi)), write_meteo
    end if

    ! destination filenames:
    if ( write_meteo ) then
      call ReadRc( rcF, 'tmm.destkey.*', rcs, destkey, status, default='no-destkey', match_index=mi )
      IF_ERROR_RETURN(status=1)
      if (mi > 0 .and. okdebug_tmm) then
        write(*,'(a, " :: for field ", a, ", ", a, " = ", a)') rname, trim(name), &
        'tmm.destkey.*'//'.'//trim(rcs(mi)), trim(destkey)
      end if

      call ReadRc( rcF, 'tmm.destkey.'//trim(lli(region)%name), rcs, destkey, status, default=destkey, match_index=mi )
      IF_ERROR_RETURN(status=1)
      if (mi > 0 .and. okdebug_tmm) then
        write(*,'(a, " :: for field ", a, ", ", a, " = ", a)') rname, trim(name), &
        'tmm.destkey.'//trim(lli(region)%name)//'.'//trim(rcs(mi)), trim(destkey)
      end if

      !if (okdebug_tmm) write(*,'("For field ", a, ", region ", a, ", destkey = ", a)') trim(name), trim(lli(region)%name), trim(destkey)
    else
      destkey = 'no-destkey'
    end if

    ! define meteo data,
    ! but should be marked as 'used' to be allocated and filled:
    call mdat_init( md, name, unit, tinterp, is, js, halo, ls, &
                   sourcekey, write_meteo, destkey, status )
    IF_NOTOK_RETURN(status=1)

    ! read this type of meteo ?
    ! only regions 1..nregions or the extra fiels above nregions_max
    ! could be in use:
    if ( (region <= nregions) .or. (region > nregions_max) ) then
      call ReadRc( rcF, 'meteo.read.*', rcs, used, status, default=.false. )
      IF_ERROR_RETURN(status=1)
      call ReadRc( rcF, 'meteo.read.'//trim(lli(region)%name), rcs, used, status, default=used )
      IF_ERROR_RETURN(status=1)
    else
      used = .false.
    end if

    ! in use ?
    call mdat_set( md, status, used=used )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Init_MeteoData


  ! ***


  subroutine Meteo_Done( status )

    use TMM        , only : Done
    use Dims       , only : nregions_all
    use meteodata  , only : mdat_Done

    ! --- in/out -------------------------------

    integer, intent(out)         ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/Meteo_Done'

    ! --- local -----------------------------

    integer   ::  n, iveg

    ! --- begin --------------------------------

    ! interface to TM meteo:
    call Done( tmmd, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done meteo data
    !

    ! destroy meteo fields:
    do n = 1, nregions_all

      ! ***

      call mdat_Done( sp1_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( sp2_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( sp_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( spm_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Done( phlb_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( m_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Done( mfu_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( mfv_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( mfw_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( tsp_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( pu_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( pv_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( pw_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Done( temper_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( humid_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( gph_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( omega_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Done( lwc_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( iwc_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done(  cc_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( cco_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( ccu_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Done( entu_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( entd_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( detu_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( detd_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Done( oro_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( lsmask_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( pclim_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( sr_ecm_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( sr_ols_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( u10m_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( v10m_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( blh_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( sshf_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( slhf_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( ewss_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( nsss_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( ci_dat(n), status )
      IF_NOTOK_RETURN(status=1)

#ifndef without_dry_deposition

      call mdat_Done( sf_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( sd_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( swvl1_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( src_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( d2m_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( t2m_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( ssr_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      do iveg = 1, nveg
        call mdat_Done( tv_dat(n,iveg), status )
        IF_NOTOK_RETURN(status=1)
      end do

      call mdat_Done( cvl_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Done( cvh_dat(n), status )
      IF_NOTOK_RETURN(status=1)
#endif

    end do   ! regions

    ! ok
    status = 0

  end subroutine Meteo_Done


  ! ***


  subroutine Meteo_Alloc( status )

    use dims       , only : nregions_all
    use meteodata  , only : mdat_Alloc

    ! --- in/out -------------------------------

    integer, intent(out)         ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Meteo_Alloc'

    ! --- local -----------------------------

    integer   ::  region
    integer   ::  iveg

    ! --- begin --------------------------------

    ! allocate meteo fields if in use:
    do region = 1, nregions_all

      ! ***

      call mdat_Alloc( sp1_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( sp2_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( sp_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( spm_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Alloc( phlb_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( m_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Alloc( mfu_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( mfv_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( mfw_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( tsp_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( pu_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( pv_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( pw_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Alloc( temper_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( humid_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( gph_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( omega_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Alloc( lwc_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( iwc_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc(  cc_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( cco_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( ccu_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Alloc( entu_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( entd_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( detu_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( detd_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      ! ***

      call mdat_Alloc( oro_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( lsmask_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( pclim_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( sr_ecm_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( sr_ols_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( u10m_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( v10m_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( blh_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( sshf_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( slhf_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( ewss_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( nsss_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( ci_dat(region), status )
      IF_NOTOK_RETURN(status=1)

#ifndef without_dry_deposition
      call mdat_Alloc( sf_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( sd_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( swvl1_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( src_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( d2m_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( t2m_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( ssr_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( cvl_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      call mdat_Alloc( cvh_dat(region), status )
      IF_NOTOK_RETURN(status=1)

      do iveg = 1, nveg
        call mdat_Alloc( tv_dat(region,iveg), status )
        IF_NOTOK_RETURN(status=1)
      end do
#endif

    end do  ! regions

    ! ok
    status = 0

  end subroutine Meteo_Alloc


  ! **************************************************************
  ! ***
  ! *** setup
  ! ***
  ! **************************************************************


  subroutine Meteo_Setup_Mass( tr1, tr2, status, &
                                 isfirst, check_pressure, restore_airmass )

    ! isfirst is .true. when called from initexit/start_TM5, not present otherwise
    ! restore_airmass is .true. when called from initexit/start_TM5 at the beginning of an adjoint run
    ! !USES: -------------------------------------
    !
    use go           , only : TDate, rTotal, operator(-), wrtgol
    use go           , only : IncrDate, operator(+), operator(<), Get
    use grid         , only : Match
    use Grid         , only : FillMassChange, BalanceMassFluxes, CheckMassBalance
    use dims         , only : nregions, nregions_all, im, jm, lm, parent
    use dims         , only : xcyc
    use TM5_Geometry , only : lli, levi
    use meteodata    , only : mdat_SetData  ! to copy %data and %tr from one MD to another
!    use restart      , only : Restart_Read
    use Var4D_IO_Mass, only : restore_masses

    ! !INPUT PARAMETERS: -------------------------
    !
    type(TDate), intent(in)        ::  tr1, tr2
    logical, intent(in), optional  ::  check_pressure
    logical, intent(in), optional  ::  isfirst
    logical, intent(in), optional  ::  restore_airmass

    ! !OUTPUT PARAMETERS: ------------------------
    !
    integer, intent(out)           ::  status

    !EOP

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Meteo_Setup_Mass'

    ! --- local -----------------------------

    logical                 ::  do_check_pressure
    logical                 ::  do_isfirst
    logical                 ::  do_restore_airmass
    integer                 ::  n, p
    integer                 ::  idater(6)
    real, allocatable       ::  dm_dt(:,:,:)
    real                    ::  dt_sec
    integer                 ::  l
    real                    ::  tol_rms, tol_diff

    ! --- begin --------------------------------
    !write (gol,'(a," : entering")') trim(rname) ; call goPr

    ! check pressure ? true if called from tracer.F90, false if called from initexit.F90
    if ( present(check_pressure) ) then
      do_check_pressure = check_pressure
    else
      do_check_pressure = .false.
    end if

    ! initial call ?
    if ( present(isfirst) ) then
      do_isfirst = isfirst
    else
      do_isfirst = .false.
    end if

    ! restore airmass ? true at the start of an adjoint run
    if ( present(restore_airmass) ) then
      do_restore_airmass = restore_airmass
    else
      do_restore_airmass = .false.
    end if

    ! info ..
!    call wrtgol( 'setup mass meteo from: ', tr1, ' to: ', tr2 ); call goPr

    !
    ! ** mass fluxes *************************************
    !

    ! loop over regions:
    do n = 1, nregions_all

       !if (mfu_dat(n)%used) then
       !   write (gol,'("mfu,mfv ",a)') trim(lli(n)%name); call goPr
       !end if

       ! update horizontal u flux (unbalanced)
       call Setup_MFUV( mfu_dat(n), mfv_dat(n), (/tr1,tr2/), lli(n), levi, status )
       IF_NOTOK_RETURN(status=1)

    end do  ! regions

    ! **

    ! loop over regions:
    do n = 1, nregions_all

       !if (mfw_dat(n)%used) then
       !   write (gol,'("mfw ",a)') trim(lli(n)%name); call goPr
       !end if

       ! update vertical flux;
       ! tendency of surface pressure is by-product of vertical flux from spectral fields
       ! or filled with zero's
       call Setup_MFW( mfw_dat(n), tsp_dat(n), (/tr1,tr2/), lli(n), 'n', levi, 'w', status )
       IF_NOTOK_RETURN(status=1)

    end do  ! regions


    !
    ! ** surface pressures *****************************
    !

    ! loop over regions:
    do n = 1, nregions_all

      ! skip ?
      if ( .not. sp1_dat(n)%used ) cycle

      !write (gol,'("sp1 ",a)') trim(lli(n)%name); call goPr

      ! Advance 'next' surface pressure (a/k/a sp2%data) to start of
      ! new interval tr1. If start of a new meteo interval, then data
      ! is automatically read from file, or recieved from coupler
      ! with OASIS/prism
      call Setup( sp2_dat(n), (/tr1,tr1/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! copy SP2 into SP1 (%data and %tr)...
      call mdat_SetData( sp1_dat(n), sp2_dat(n), status )
      IF_NOTOK_RETURN(status=1)

      !! testing ...
      !call wrtgol( ' PPP NEW: fill sp1 for tr1 = ', tr1 ); call goPr
      !write (gol,*) 'PPP NEW: sp1 region ', n, ' mean = ', sum(sp1_dat(n)%data)/size(sp1_dat(n)%data); call goPr

      ! global field (first region) ?
      ! then match with average global surface pressure to ensure
      ! global mass balance;
      ! otherwise, match with parent grid:
      if ( n == 1 ) then
        call Match( 'area-aver', 'n', lli(0), sp_region0, &
                                      lli(n), sp1_dat(n)%data(1:im(n),1:jm(n),1), status )
        IF_NOTOK_RETURN(status=1)
      else
        p = parent(n)
        call Match( 'area-aver', 'n', lli(p), sp1_dat(p)%data(1:im(p),1:jm(p),1), &
                                      lli(n), sp1_dat(n)%data(1:im(n),1:jm(n),1), status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! Initial call ? then set current surface pressure to just
      ! read/advanced sp1.
      ! otherwise, sp remains filled with the advected pressure.
      if ( do_isfirst ) then
        ! info ...
        write (gol,'("  copy SP1 to SP ...")'); call goPr

        ! copy sp1 into sp :
        call mdat_SetData( sp_dat(n), sp1_dat(n), status )
        IF_NOTOK_RETURN(status=1)

        ! fill pressure and mass from sp:
        call Pressure_to_Mass( n, status )
        IF_NOTOK_RETURN(status=1)

        ! eventually replace by fields in restart file, since meteo
        ! from hdf meteo files is in real(4) while computed
        ! pressures and mass are probably in real(8) ;

        !!write (gol,'("  Reading initial sp from restart file")'); call goPr
        !call Restart_Read( status, region=n, surface_pressure=.true., pressure=.true., air_mass=.true. )
        !IF_NOTOK_RETURN(status=1)

        ! used by adjoint run to restore initial air mass:
        if ( do_restore_airmass ) then
          ! read final mass array from forward run:
          call restore_masses( status, only_region=n )
          IF_NOTOK_RETURN(status=1)
          ! fill phlb and sp (surface phlb) from m:
          call Mass_to_Pressure( n, status )
          IF_NOTOK_RETURN(status=1)
        end if

        !AJS>>> don't do this! sp1 contains data interpolated between
        ! fields received from the archive or the coupled model,
        ! while sp contains the actual pressure after advection.
        !! copy sp into sp1 (PLS, 29-3-2010)
        !call mdat_setdata( sp1_dat(n), sp_dat(n), status )
        !IF_NOTOK_RETURN(status=1)
        !<<<

        ! fill halo cells:
        call FillHalo_pm( n, status )
        IF_NOTOK_RETURN(status=1)

      end if ! end first

      !! fill initial pressure and mass arrays,
      !! eventually apply cyclic boundaries to mass
      !call Meteo_SetupMass( n, status )
      !IF_NOTOK_RETURN(status=1)

      ! check 'advected' pressure ?
      if ( do_check_pressure) then
        !! testing ...
        !write (gol,'("WARNING - skipped advected pressure check")'); call goPr
        ! compare 'advected' pressure still in sp with just read
        ! pressure sp1 : diff b/w sp%data and sp1%data
        call Meteo_CheckPressure( n, status )
        IF_NOTOK_RETURN(status=1)
      end if

    end do  ! regions

    ! **

    ! loop over regions:
    do n = 1, nregions_all

      ! skip ?
      if ( .not. sp2_dat(n)%used ) cycle

      !write (gol,'("sp2 ",a)') trim(lli(n)%name); call goPr

      ! advance 'next' surface pressure to end of interval:
      call Setup( sp2_dat(n), (/tr2,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      !! testing ...
      !call wrtgol( ' PPP NEW: read sp2 for tr2 = ', tr2 ); call goPr
      !write (gol,*) 'PPP NEW: sp2 region ', n, ' mean = ', sum(sp2_dat(n)%data)/size(sp2_dat(n)%data); call goPr

      ! global field (first region) ?
      ! then match with average global surface pressure to ensure global mass balance;
      ! otherwise, match with parent grid:
      if ( n == 1 ) then
        call Match( 'area-aver', 'n', lli(0), sp_region0, &
                                      lli(n), sp2_dat(n)%data(1:im(n),1:jm(n),1), status )
        IF_NOTOK_RETURN(status=1)
      else
        p = parent(n)
        call Match( 'area-aver', 'n', lli(p), sp2_dat(p)%data(1:im(p),1:jm(p),1), &
                                      lli(n), sp2_dat(n)%data(1:im(n),1:jm(n),1), status )
        IF_NOTOK_RETURN(status=1)
      end if

    end do  ! regions


#ifndef without_advection

    !
    ! ** mass balance *****************************
    !

    ! NOTE: since only the surface pressure gradient is used,
    ! it is not necessary to use the data1 and data2 arrays

    ! loop over regions:
    do n = 1, nregions_all

      ! skip ?
      if ( .not. pu_dat(n)%used ) cycle
      if ( .not. pv_dat(n)%used ) cycle
      if ( .not. pw_dat(n)%used ) cycle

      !write (gol,'("balance ",a)') trim(lli(n)%name); call goPr

      ! length of time step between sp1 and sp2:
      dt_sec = rTotal( sp2_dat(n)%tr(1) - sp1_dat(n)%tr(1), 'sec' )

      ! allocate temporary array:
      allocate( dm_dt(im(n),jm(n),lm(n)) )

      ! mass change (kg) :
      call FillMassChange( dm_dt, lli(n), levi, &
                           sp1_dat(n)%data(1:im(n),1:jm(n),1), &
                           sp2_dat(n)%data(1:im(n),1:jm(n),1), &
                           status )
      IF_NOTOK_RETURN(status=1)

      ! mass tendency (kg/s) :
      dm_dt = dm_dt / dt_sec   ! kg/s

      ! >>> data1 >>>

      ! initial guess for balanced fluxes are unbalanced fluxes:
      pu_dat(n)%data1   = mfu_dat(n)%data1
      pu_dat(n)%filled1 = mfu_dat(n)%filled1
      pu_dat(n)%tr1     = mfu_dat(n)%tr1
      pv_dat(n)%data1   = mfv_dat(n)%data1
      pv_dat(n)%filled1 = mfv_dat(n)%filled1
      pv_dat(n)%tr1     = mfv_dat(n)%tr1
      pw_dat(n)%data1   = mfw_dat(n)%data1
      pw_dat(n)%filled1 = mfw_dat(n)%filled1
      pw_dat(n)%tr1     = mfw_dat(n)%tr1

      ! match with parent grid if necessary; note strange indexing:
      !   pu_dat(n)%data1( 0:im(n), 1:jm(n)  , 1:lm(n) )
      !   pv_dat(n)%data1( 1:im(n), 1:jm(n)+1, 1:lm(n) )
      if ( n > 1 ) then
        p = parent(n)
        !
        do l = 1, lm(n)
          call Match( 'sum', 'u', lli(p), pu_dat(p)%data1(0:im(p),1:jm(p),l), &
                                  lli(n), pu_dat(n)%data1(0:im(n),1:jm(n),l), status )
          IF_NOTOK_RETURN(status=1)
        end do
        !
        do l = 1, lm(n)
          call Match( 'sum', 'v', lli(p), pv_dat(p)%data1(1:im(p),1:jm(p)+1,l), &
                                  lli(n), pv_dat(n)%data1(1:im(n),1:jm(n)+1,l), status )
          IF_NOTOK_RETURN(status=1)
        end do
        !
        do l = 0, lm(n)
          call Match( 'sum', 'n', lli(p), pw_dat(p)%data1(1:im(p),1:jm(p),l), &
                                  lli(n), pw_dat(n)%data1(1:im(n),1:jm(n),l), status )
          IF_NOTOK_RETURN(status=1)
        end do
        !
      end if


!#ifdef with_prism
      ! skip initial mass balance; relative large differences might exist
      ! between pressure imposed by mass fluxes and pressure according to
      ! surface pressure tendencies since the later is based on:
      !
      !    sp(t-1)+tsp(t-1)      _ *
      !                      _ -   o-------*     sp(t), sp(t)+tsp(t)
      !    sp(t-1)         o
      !
      ! PLS : I do not understand that diagram... tsp is for an
      !       interval, and sp for a point in time. This may be
      !       wrong then. What we had at the first time step was:
      !
      !    sp(t+1)+tsp(t:t+1)      _ *
      !                        _ -       =>  sp(t) to sp(t)+tsp(t:t+1)
      !    sp(t+1)           o
      !
      ! AJS : This describes what the CTM received before the above
      !   described update. The 'tsp' was *not* for an interval but
      !   an instantaneous field describing the 'direction' of the surface
      !   pressure in time (you might call this 'tendency', but that is a
      !   dangerous word in GEMS IFS-CTM coupling context).
      !   Thus, at time 't-1' the only estimate of 'sp(t)' we could make was:
      !     sp(t-1)+tsp(t-1)
      !   At time 't' we received the actual 'sp(t)' and this was of course
      !   different from the initial guess.
      !
      ! PLS : Just need to be sure that we have the correct sp to start
      ! with. Code above has been modified, so that we have:
      !
      !    sp(t)+tsp(t:t+1)      _ *
      !                      _ -       =>  sp(t) to sp(t)+tsp(t:t+1)
      !    sp(t)           o
      !
!#else
      ! check initial mass balance:
      ! NOTE: strange old indexing:
      !   pu_tmpp  -->  pu(0:im(n),1:jm(n)  ,1:lm(n))  in  pu_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
      !   pv_tmpp  -->  pv(1:im(n),1:jm(n)+1,1:lm(n))  in  pv_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
      ! tolerance for difference between sp from mass fluxes and sp from tendency:
      tol_rms  = 1.0e-4    ! max rms
      tol_diff = 1.0e-3    ! max absolute difference
      call CheckMassBalance( lli(n), &
                               pu_dat(n)%data1(0:im(n),1:jm(n)  ,1:lm(n)), &
                               pv_dat(n)%data1(1:im(n),1:jm(n)+1,1:lm(n)), &
                               sp1_dat(n)%data(1:im(n),1:jm(n),1), &
                               sp2_dat(n)%data(1:im(n),1:jm(n),1), &
                               dt_sec, tol_rms, tol_diff, status )
      if (status/=0) then
        write (gol,'("initial mass imbalance too large for region ",i2)') n; call goErr
        TRACEBACK; status=1; return
      end if
!#endif

      ! balance horizontal fluxes:
      ! NOTE: strange old indexing:
      !   pu_tmpp  -->  pu(0:im(n),1:jm(n)  ,1:lm(n))  in  pu_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
      !   pv_tmpp  -->  pv(1:im(n),1:jm(n)+1,1:lm(n))  in  pv_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
      call BalanceMassFluxes( lli(n), &
                               pu_dat(n)%data1(0:im(n),1:jm(n)  ,1:lm(n)), &
                               pv_dat(n)%data1(1:im(n),1:jm(n)+1,1:lm(n)), &
                               pw_dat(n)%data1, dm_dt, lli(parent(n)), dt_sec, status )
      IF_NOTOK_RETURN(status=1)

      ! check final mass balance:
      ! NOTE: strange old indexing:
      !   pu_tmpp  -->  pu(0:im(n),1:jm(n)  ,1:lm(n))  in  pu_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
      !   pv_tmpp  -->  pv(1:im(n),1:jm(n)+1,1:lm(n))  in  pv_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
      ! tolerance for difference between sp from mass fluxes and sp from tendency:
      tol_rms  = 1.0e-7    ! max rms
      tol_diff = 1.0e-6    ! max absolute difference
      call CheckMassBalance( lli(n), &
                               pu_dat(n)%data1(0:im(n),1:jm(n)  ,1:lm(n)), &
                               pv_dat(n)%data1(1:im(n),1:jm(n)+1,1:lm(n)), &
                               sp1_dat(n)%data(1:im(n),1:jm(n),1), &
                               sp2_dat(n)%data(1:im(n),1:jm(n),1), &
                               dt_sec, tol_rms, tol_diff, status )
      if (status/=0) then
        write (gol,'("final mass imbalance too large for region ",i2)') n; call goErr
        TRACEBACK; status=1; return
      end if

      ! periodic boundary
      if ( xcyc(n) == 1 ) then
        pu_dat(n)%data1(0      ,:,:) = pu_dat(n)%data1(im(n),:,:)
        pu_dat(n)%data1(im(n)+1,:,:) = pu_dat(n)%data1(1    ,:,:)
      end if

      ! >>> data2 >>>

      if ( any((/mfu_dat%filled2,mfv_dat%filled2,mfw_dat%filled2/)) ) then

        ! check ...
        if ( .not. all((/mfu_dat(n)%filled2,mfv_dat(n)%filled2,mfw_dat(n)%filled2/)) ) then
          write (gol,'("either none or all secondary data should be in use:")'); call goErr
          write (gol,'("  mfu_dat%filled2  : ",l1)') mfu_dat(n)%filled2; call goErr
          write (gol,'("  mfv_dat%filled2  : ",l1)') mfv_dat(n)%filled2; call goErr
          write (gol,'("  mfw_dat%filled2  : ",l1)') mfw_dat(n)%filled2; call goErr
          TRACEBACK; status=1; return
        end if

        ! initial guess for balanced fluxes are unbalanced fluxes:
        pu_dat(n)%data2   = mfu_dat(n)%data2
        pu_dat(n)%filled2 = .true.
        pu_dat(n)%tr2     = mfu_dat(n)%tr2
        pv_dat(n)%data2   = mfv_dat(n)%data2
        pv_dat(n)%filled2 = .true.
        pv_dat(n)%tr2     = mfv_dat(n)%tr2
        pw_dat(n)%data2   = mfw_dat(n)%data2
        pw_dat(n)%filled2 = .true.
        pw_dat(n)%tr2     = mfw_dat(n)%tr2

        ! match with parent grid if necessary; note strange indexing:
        !   pu_dat(n)%data2( 0:im(n), 1:jm(n)  , 1:lm(n) )
        !   pv_dat(n)%data2( 1:im(n), 1:jm(n)+1, 1:lm(n) )
        if ( n > 1 ) then
          p = parent(n)
          !
          do l = 1, lm(n)
            call Match( 'sum', 'u', lli(p), pu_dat(p)%data2(0:im(p),1:jm(p),l), &
                                    lli(n), pu_dat(n)%data2(0:im(n),1:jm(n),l), status )
            IF_NOTOK_RETURN(status=1)
          end do
          !
          do l = 1, lm(n)
            call Match( 'sum', 'v', lli(p), pv_dat(p)%data2(1:im(p),1:jm(p)+1,l), &
                                    lli(n), pv_dat(n)%data2(1:im(n),1:jm(n)+1,l), status )
            IF_NOTOK_RETURN(status=1)
          end do
          !
          do l = 0, lm(n)
            call Match( 'sum', 'n', lli(p), pw_dat(p)%data2(1:im(p),1:jm(p),l), &
                                    lli(n), pw_dat(n)%data2(1:im(n),1:jm(n),l), status )
            IF_NOTOK_RETURN(status=1)
          end do
          !
        end if

        ! check initial mass balance:
        ! NOTE: strange old indexing:
        !   pu_tmpp  -->  pu(0:im(n),1:jm(n)  ,1:lm(n))  in  pu_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
        !   pv_tmpp  -->  pv(1:im(n),1:jm(n)+1,1:lm(n))  in  pv_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
        call CheckMassBalance( lli(n), &
                                 pu_dat(n)%data2(0:im(n),1:jm(n)  ,1:lm(n)), &
                                 pv_dat(n)%data2(1:im(n),1:jm(n)+1,1:lm(n)), &
                                 sp1_dat(n)%data(1:im(n),1:jm(n),1), &
                                 sp2_dat(n)%data(1:im(n),1:jm(n),1), &
                                 dt_sec, 1.0e-4, 1.0e-3, status )
        if (status/=0) then
          write (gol,'("initial mass imbalance too large for region ",i2)') n; call goErr
          TRACEBACK; status=1; return
        end if

        ! balance horizontal fluxes:
        ! NOTE: strange old indexing:
        !   pu_tmpp  -->  pu(0:im(n),1:jm(n)  ,1:lm(n))  in  pu_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
        !   pv_tmpp  -->  pv(1:im(n),1:jm(n)+1,1:lm(n))  in  pv_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
        call BalanceMassFluxes( lli(n), &
                                 pu_dat(n)%data2(0:im(n),1:jm(n)  ,1:lm(n)), &
                                 pv_dat(n)%data2(1:im(n),1:jm(n)+1,1:lm(n)), &
                                 pw_dat(n)%data2, dm_dt, lli(parent(n)), dt_sec, status )
        IF_NOTOK_RETURN(status=1)

        ! check final mass balance:
        ! NOTE: strange old indexing:
        !   pu_tmpp  -->  pu(0:im(n),1:jm(n)  ,1:lm(n))  in  pu_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
        !   pv_tmpp  -->  pv(1:im(n),1:jm(n)+1,1:lm(n))  in  pv_t(0:im(n)+1,0:jm(n)+1,0:lm(n))
        call CheckMassBalance( lli(n), &
                                 pu_dat(n)%data2(0:im(n),1:jm(n)  ,1:lm(n)), &
                                 pv_dat(n)%data2(1:im(n),1:jm(n)+1,1:lm(n)), &
                                 sp1_dat(n)%data(1:im(n),1:jm(n),1), &
                                 sp2_dat(n)%data(1:im(n),1:jm(n),1), &
                                 dt_sec, 1.0e-7, 1.0e-6, status )
        if (status/=0) then
          write (gol,'("final mass imbalance too large for region ",i2)') n; call goErr
          TRACEBACK; status=1; return
        end if

        ! periodic boundary
        if ( xcyc(n) == 1 ) then
          pu_dat(n)%data2(0      ,:,:) = pu_dat(n)%data2(im(n),:,:)
          pu_dat(n)%data2(im(n)+1,:,:) = pu_dat(n)%data2(1    ,:,:)
        end if

      end if  ! filled2

      ! >>>

      ! clear
      deallocate( dm_dt )

    end do  ! regions

#endif


    !
    ! ** done **************************************************
    !
    !write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Meteo_Setup_Mass


  ! ***


  subroutine Meteo_Setup_Other( tr1, tr2, status, isfirst )

    use GO          , only : TDate, wrtgol, operator(<)
    use GO          , only : operator(-), operator(+), operator(/)
    use GO          , only : InterpolFractions
    use Dims        , only : nregions_all
    use dims        , only : im, jm
    use dims        , only : lmax_conv
    use Dims        , only : czeta
    use TM5_Geometry, only : lli, levi
    use global_data , only : region_dat
#ifndef without_convection
    use global_data , only : conv_dat
#endif
    use Phys       , only : ConvCloudDim

    ! --- in/out -------------------------------

    type(TDate), intent(in)        ::  tr1, tr2
    integer, intent(out)           ::  status
    logical, intent(in), optional  ::  isfirst

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Meteo_Setup_Other'

    ! --- local -----------------------------

    !logical                 ::  do_isfirst

    integer                 ::  n, iveg

    integer                 ::  i, j, l
    integer                 ::  lsave
    real                    ::  tote, totd, maxe
    real, pointer           ::  dxyp(:)

    type(TDate)             ::  tmid
    real                    ::  alfa1, alfa2

    ! --- begin --------------------------------
    !write (gol,'(a," : entering")') trim(rname) ; call goPr

    !! initial call ?
    !if ( present(isfirst) ) then
    !  do_isfirst = isfirst
    !else
    !  do_isfirst = .false.
    !end if

    ! info ...
!    call wrtgol( 'setup other meteo from: ', tr1, ' to: ', tr2 ); call goPr

    !
    ! ** orography *****************************
    !

    ! read orographies (if necessary):
    do n = 1, nregions_all
      !if (oro_dat(n)%used) then; write (gol,'("oro ",a)') trim(lli(n)%name); call goPr; end if
      call Setup( oro_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! climatological pressure (derived from oro):
    do n = 1, nregions_all
      !if (pclim_dat(n)%used) then; write (gol,'("pclim ",a)') trim(lli(n)%name); call goPr; end if
      call Compute_PClim( n, status )
      IF_NOTOK_RETURN(status=1)
    end do

    !
    ! ** spm **************************************
    !

    ! loop over regions:
    do n = 1, nregions

      ! skip ?
      if ( .not. spm_dat(n)%used ) cycle

      !write (gol,'("spm ",a)') trim(lli(n)%name); call goPr

      ! mid time:
      tmid = tr1 + ( tr2 - tr1 )/2

      ! deterimine weights to sp1 and sp2 :
      call InterpolFractions( tmid, sp1_dat(n)%tr(1), sp2_dat(n)%tr(1), alfa1, alfa2, status )
      IF_NOTOK_RETURN(status=1)

      ! interpolate:
      spm_dat(n)%data(1:im(n),1:jm(n),1) = &
          alfa1 * sp1_dat(n)%data(1:im(n),1:jm(n),1) + &
          alfa2 * sp2_dat(n)%data(1:im(n),1:jm(n),1)

      !! testing ...
      !if ( n == 1 ) then
      !  call wrtgol( '  spm: interpolate between ', sp1_dat(n)%tr(1), ' and ', sp2_dat(n)%tr(1) ); call goPr
      !  call wrtgol( '  spm: to tmid ', tmid ); call goPr
      !  write (gol,'("  spm: weights ",2f10.4)') alfa1, alfa2; call goPr
      !  write (gol,'("  spm: pressures ",3f16.4)') sp1_dat(n)%data(10,10,1), sp2_dat(n)%data(10,10,1), spm_dat(n)%data(10,10,1); call goPr
      !end if

      ! store time:
      spm_dat(n)%tr = (/tr1,tr2/)
      ! reset flag:
      spm_dat(n)%changed = .true.


    end do  ! regions

    !
    ! ** omega **************************************
    !

    ! loop over regions:
    do n = 1, nregions_all

      !if (omega_dat(n)%used) then; write (gol,'("omega ",a)') trim(lli(n)%name); call goPr; end if

      ! re-compute omega from vertical mass flux:
      call Compute_Omega( omega_dat(n), lli(n), mfw_dat(n), status )
      IF_NOTOK_RETURN(status=1)

    end do  ! regions


    !
    ! ** temperature and humid **************************************
    !

    ! loop over regions:
    do n = 1, nregions

      ! read temperature (if necessary):
      call Setup( temper_dat(n), (/tr1,tr2/), lli(n), 'n', levi, 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! read humidity (if necessary):
      call Setup( humid_dat(n), (/tr1,tr2/), lli(n), 'n', levi, 'n', status )
      IF_NOTOK_RETURN(status=1)

    end do  ! regions


    !
    ! ** gph **************************************
    !

    ! loop over regions:
    do n = 1, nregions_all

       !write (gol,'("gph for region ",a," (#",i1,") :",l2)') trim(lli(n)%name), n, gph_dat(n)%used; call goPr
       !if (gph_dat(n)%used) then; write (gol,'("gph ",a)') trim(lli(n)%name); call goPr; end if

       ! re-compute gph from pressure, temperature, and humidity:
       call compute_gph( n, status )
       IF_NOTOK_RETURN(status=1)

    end do  ! regions


    !
    ! ** clouds **************************************
    !

    ! loop over regions:
    do n = 1, nregions

      !if (any((/lwc_dat(n)%used,iwc_dat(n)%used,cc_dat(n)%used,cco_dat(n)%used,ccu_dat(n)%used/))) then
      !  write (gol,'("clouds ",a)') trim(lli(n)%name); call goPr
      !end if

      call Setup( lwc_dat(n), (/tr1,tr2/), lli(n), 'n', levi, 'n', status )
      IF_NOTOK_RETURN(status=1)

      call Setup( iwc_dat(n), (/tr1,tr2/), lli(n), 'n', levi, 'n', status )
      IF_NOTOK_RETURN(status=1)

      call Setup_CloudCovers( cc_dat(n), cco_dat(n), ccu_dat(n), (/tr1,tr2/), lli(n), levi, status )
      IF_NOTOK_RETURN(status=1)

    end do


    !
    ! ** convection **************************************
    !

    ! loop over regions:
    do n = 1, nregions

      !if (entu_dat(n)%used) then; write (gol,'("convection ",a)') trim(lli(n)%name); call goPr; end if

      ! read (if necessary):
      call Setup_Convec( entu_dat(n), entd_dat(n), detu_dat(n), detd_dat(n), &
                            omega_dat(n), gph_dat(n), (/tr1,tr2/), lli(n), levi, status )
      IF_NOTOK_RETURN(status=1)

    end do

#ifndef without_convection
    ! ~~ convective clouds

    ! loop over regions:
    do n = 1, nregions

      ! skip ?
      if ( .not. entu_dat(n)%used ) cycle
      if ( .not. entd_dat(n)%used ) cycle

      ! update necessary ?
      if ( any((/entu_dat(n)%changed,entd_dat(n)%changed/)) ) then

        ! loop over grid cells
        do j = 1, jm(n)
          do i = 1, im(n)

            ! compute convective cloud dimensions for this column:
            call ConvCloudDim( 'u', size(detu_dat(n)%data,3), &
                     detu_dat(n)%data(i,j,:), entd_dat(n)%data(i,j,:),&
                     conv_dat(n)%cloud_base(i,j), &
                     conv_dat(n)%cloud_top (i,j), &
                     conv_dat(n)%cloud_lfs (i,j), &
                     status )
            IF_NOTOK_RETURN(status=1)

          end do    ! i
        end do   ! j

      end if  ! changed ?

    end do  ! regions
#endif

    ! ~~ unit conversion

    ! loop over regions:
    do n = 1, nregions

      ! skip ?
      if ( .not. entu_dat(n)%used ) cycle
      if ( .not. entd_dat(n)%used ) cycle
      if ( .not. detu_dat(n)%used ) cycle
      if ( .not. detd_dat(n)%used ) cycle

      ! update necessary ?
      if ( any((/entu_dat(n)%changed,entd_dat(n)%changed,&
                 detu_dat(n)%changed,detd_dat(n)%changed/)) ) then

        !cmk calculate the rates in kg/gridbox and scale with czeta

        dxyp => region_dat(n)%dxyp

        do j = 1, jm(n)
          do i = 1, im(n)

            ! kg/m2/s -> kg/gridbox/s * scale_factor
            entu_dat(n)%data(i,j,:) = entu_dat(n)%data(i,j,:)*dxyp(j)*czeta
            detu_dat(n)%data(i,j,:) = detu_dat(n)%data(i,j,:)*dxyp(j)*czeta

            ! ensure netto zero tracer transport by updraught in column
            ! (add difference between total entrement and detrement
            ! to level where entrement reaches maximum):
            tote = sum( entu_dat(n)%data(i,j,:) )
            totd = sum( detu_dat(n)%data(i,j,:) )
            maxe = entu_dat(n)%data(i,j,1)  ! changed: reported by PB feb 2003
            lsave = 1
            do l = 2, lmax_conv
              if ( entu_dat(n)%data(i,j,l) > maxe ) then
                maxe = entu_dat(n)%data(i,j,l)
                lsave = l
              end if
            end do
            entu_dat(n)%data(i,j,lsave) = entu_dat(n)%data(i,j,lsave) - tote + totd

            ! kg/m2/s -> kg/gridbox/s * scale_factor
            entd_dat(n)%data(i,j,:) = entd_dat(n)%data(i,j,:)*dxyp(j)*czeta
            detd_dat(n)%data(i,j,:) = detd_dat(n)%data(i,j,:)*dxyp(j)*czeta

            ! ensure netto zero tracer transport by downdraught in column
            ! (add difference between total entrement and detrement
            ! to level where entrement reaches maximum):
            tote = sum( entd_dat(n)%data(i,j,:) )   ! total entrainement
            totd = sum( detd_dat(n)%data(i,j,:) )   ! total detrainement
            maxe = 0.0
            lsave = lmax_conv
            do l = 1, lmax_conv
              if ( entd_dat(n)%data(i,j,l) > maxe ) then
                maxe = entd_dat(n)%data(i,j,l)
                lsave = l
              end if
            end do
            entd_dat(n)%data(i,j,lsave) = entd_dat(n)%data(i,j,lsave) - tote + totd

          end do
        end do

      end if  ! changed ?

    end do  ! regions


    !
    ! ** surface fields *****************************
    !

    do n = 1, nregions_all

      !write (gol,'("surface fields ",a)') trim(lli(n)%name); call goPr

      ! * lsmask

      call Setup( lsmask_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * sr_ecm

      call Setup( sr_ecm_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * sr_ols

      call Setup( sr_ols_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * u10m

      call Setup( u10m_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * v10m

      call Setup( v10m_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * boundary layer height

      call Setup( blh_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * slhf

      call Setup( slhf_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * sshf

      call Setup( sshf_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * surface stress

      call Setup( ewss_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      call Setup( nsss_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! sea ice is also required for C14 isoflux, not just dry deposition
      call Setup( ci_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

#ifndef without_dry_deposition

      ! * 2m dewpoint temperature

      call Setup( d2m_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * 2m temperature

      call Setup( t2m_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * surface solar radiation

      call Setup( ssr_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * snow fall and depth

      call Setup( sf_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      call Setup( sd_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * soil water level 1

      call Setup( swvl1_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * vegetation types

      do iveg = 1, nveg
        select case ( iveg )
          case ( 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 16, 17, 18, 19 )
            call Setup( tv_dat(n,iveg), (/tr1,tr2/), lli(n), 'n', status )
            IF_NOTOK_RETURN(status=1)
          case ( 8, 12, 14, 15, 20 )
            if ( tv_dat(n,iveg)%used ) tv_dat(n,iveg)%data = 0.0
          case default
            write (gol,'("do not know how to setup vegetation type ",i2)') iveg
            call goErr; status=1; return
        end select
      end do

      ! * low vegetation cover

      call Setup( cvl_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! * high vegetation cover

      call Setup( cvh_dat(n), (/tr1,tr2/), lli(n), 'n', status )
      IF_NOTOK_RETURN(status=1)

#endif

    end do   ! regions


    !
    ! ** done ********************************************
    !
    !write (gol,'(a," : done")') trim(rname) ; call goPr
    ! ok
    status = 0

  end subroutine Meteo_Setup_Other


  ! ***

  subroutine SetupSetup( md, tr, &
                          data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
                          data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
                          status )

    use GO         , only : TDate, NewDate, IncrDate, AnyDate, IsAnyDate, Get, Set, wrtgol
    use GO         , only : rTotal, iTotal
    use GO         , only : operator(+), operator(-), operator(/)
    use GO         , only : operator(==), operator(/=), operator(<), operator(>), operator(<=)
    use meteodata  , only : TMeteoData
    !use global_data, only : fcmode, tfcday0

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  md
    type(TDate), intent(in)               ::  tr(2)
    logical, intent(out)                  ::  data1_read, data1_copy
    type(TDate), intent(out)              ::  data1_tref, data1_t1, data1_t2
    logical, intent(out)                  ::  data2_read, data2_copy
    type(TDate), intent(out)              ::  data2_tref, data2_t1, data2_t2
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SetupSetup'

    ! --- local ----------------------------------

    integer                         ::  direction
    integer                         ::  dth, baseh
    integer                         ::  year, month, day, hour, minu
    type(TDate)                     ::  trm
    type(TDate)                     ::  tmid
    type(TDate)                     ::  tc(2)
    real                            ::  frac
    integer                         ::  ifrac
    integer                         ::  dth_int
    type(TDate)                     ::  tprev, tnext
    real                            ::  dhr

    ! --- begin -----------------------------

    !call goLabel(rname)

    ! default output:
    data1_read    = .false.
    data1_copy    = .false.
    data2_read    = .false.
    data2_copy    = .false.

    ! time direction:
    if ( tr(1) <= tr(2) ) then
      direction =  1
    else
      direction = -1
    end if

    !
    ! trap constant fields ...
    !

    ! constant and already filled ? then leave
    if ( (md%tinterp == 'const') .and. md%filled1 ) then
      !call goLabel()
      status = 0; return
    end if

    !
    ! fc stuff
    !

    !! 3 hourly data only available up to 72h , then 6 hourly
    !if ( fcmode ) then
    !  ! number of hours from fcday 00:00 to end of requested interval:
    !  dhr = rTotal( tr(2) - tfcday0, 'hour' )
    !  ! lower time resolution after a while ...
    !  if ( tfcday0 < NewDate(year=2006,month=03,day=14) ) then
    !    ! after 12+72 hour ?
    !    if ( dhr > 12.0 + 72.0 ) then
    !      ! convert time interpolation:
    !      select case ( md%tinterp )
    !        case ( 'aver3'   )
    !          write (gol,'("WARNING - convert time interpolation from `aver3` to `aver6`")'); call goPr
    !          md%tinterp = 'aver6'
    !        case ( 'interp3' )
    !          write (gol,'("WARNING - convert time interpolation from `interp3` to `interp6`")'); call goPr
    !          md%tinterp = 'interp6'
    !      end select
    !    end if  ! > 72 hour
    !  else
    !    ! after 12+96 hour ?
    !    if ( dhr > 12.0 + 96.0 ) then
    !      ! convert time interpolation:
    !      select case ( md%tinterp )
    !        case ( 'aver3'   )
    !          write (gol,'("WARNING - convert time interpolation from `aver3` to `aver6`")'); call goPr
    !          md%tinterp = 'aver6'
    !        case ( 'interp3' )
    !          write (gol,'("WARNING - convert time interpolation from `interp3` to `interp6`")'); call goPr
    !          md%tinterp = 'interp6'
    !      end select
    !    end if  ! > 96 hour
    !  end if  ! change in fc resolution
    !end if  ! fcmode


    !
    ! time stuff
    !

    ! basic time resolution in hours
    select case ( md%tinterp )
      case ( 'const', 'month' )
        ! nothing to be set here ...
      case ( 'aver24' )
        ! constant fields produced valid for [00,24]
        dth   = 24
        baseh = 00
      case ( 'aver24_3' )
        ! constant fields produced by tmpp valid for [21,21] = [09-12,09+12]
        dth   = 24
        baseh = -3
      case ( 'const3', 'interp3', 'aver3', 'cpl3' )
        dth = 3
        baseh = 0
      case ( 'interp2', 'cpl2' )
        dth = 2
        baseh = 0
      case ( 'const1', 'interp1', 'aver1', 'cpl1' )
        dth = 1
        baseh = 0
      case ( 'const6', 'interp6', 'aver6', 'cpl6' )
        dth = 6
        baseh = 0
      case ( 'interp6_3' )
        dth = 6
        baseh = 3
      case default
        write (gol,'("unsupported time interpolation : ",a)') md%tinterp; call goErr
        TRACEBACK; status=1; return
    end select

    ! set time parameters for field to be read:
    select case ( md%tinterp )
      !
      ! ** constant fields
      !
      case ( 'const' )
        ! read main field ?
        data1_read = .not. md%filled1
        ! read or leave ?
        if ( data1_read ) then
          data1_tref = tr(1)        !  <--- used for file names
          data1_t1   = AnyDate()
          data1_t2   = AnyDate()
        else
          ! field valid around requested interval, thus leave:
          !call goLabel();
          status=0; return
        end if
      !
      ! ** constant fields, valid for complete month
      !
      case ( 'month' )
        ! mid of requested interval:
        trm = tr(1) + (tr(2)-tr(1))/2
        ! extract time values:
        call Get( trm, year=year, month=month )
        ! interval for this month:
        if ( direction > 0 ) then
        tc(1) = NewDate( year=year, month=month, day=01, hour=00 )
        else
          tc(2) = NewDate( year=year, month=month, day=01, hour=00 )
        end if
        month = month + 1
        if ( month > 12 ) then
          month = 1
          year = year + 1
        end if
        if ( direction > 0 ) then
        tc(2) = NewDate( year=year, month=month, day=01, hour=00 )
        else
          tc(1) = NewDate( year=year, month=month, day=01, hour=00 )
        end if
        ! check for strange values:
        if ( ((direction > 0) .and. ((tr(1) < tc(1)) .or. (tc(2) < tr(2)))) .or. &
             ((direction < 0) .and. ((tr(2) < tc(2)) .or. (tc(1) < tr(1))))      ) then
          write (gol,'("determined invalid constant interval:")'); call goErr
          call wrtgol( '  requested   : ', tr(1), ' - ', tr(2) ); call goErr
          call wrtgol( '  guessed     : ', tc(1), ' - ', tc(2) ); call goErr
          write (gol,'("  for tinterp : ",a)') md%tinterp; call goErr
          TRACEBACK; status=1; return
          !write (gol,'("          WARNING - requested interval exceeds meteo interval; should be improved")')
        end if
        ! read main field ?
        if ( md%filled1 ) then
          data1_read = md%tr1(1) /= tc(1)
        else
          data1_read = .true.
        end if
        ! read or leave ?
        if ( data1_read ) then
          data1_tref = tr(1)
          data1_t1   = tc(1)
          data1_t2   = tc(2)
        else
          ! field valid around requested interval, thus leave:
          !call goLabel();
          status=0; return
        end if
      !
      ! ** constant fields, valid for 24hr intervals  [21:00,21:00]
      !    constant fields, valid for  6hr intervals  [21:00,03:00]  etc
      !    constant fields, valid for  3hr intervals  [22:30,01:30]  etc
      !
      !  Assignment of instant time (tr(1)==tr(2)) to intervals:
      !    [21,03) [03,09) [09,12) ...
      !
      !
      case ( 'const6', 'const3' )
        ! mid of requested interval:
        trm = tr(1) + (tr(2)-tr(1))/2
        ! extract time values:
        call Get( trm, year, month, day, hour, minu )
        ! round hour to 00/06/12/18 or 00/03/06/09/12/15/18/21 or 09 ;
        ! do not use 'nint', this gives unpredictable rounding of 0.5
        ! with some compilers ...
        frac = real(hour+minu/60.0-baseh) / real(dth)
        if ( frac - int(frac) < 0.5 ) then
          ifrac = int(frac)
        else
          ifrac = int(frac) + 1
        end if
        hour = baseh + dth * ifrac
        ! set mid of 3 or 6 hour interval:
        tmid = NewDate( year, month, day, hour )
        ! interval with constant field
        tc(1) = tmid - IncrDate(hour=direction*dth)/2
        tc(2) = tmid + IncrDate(hour=direction*dth)/2
        ! check for strange values:
        if ( ((direction > 0) .and. ((tr(1) < tc(1)) .or. (tc(2) < tr(2)))) .or. &
             ((direction < 0) .and. ((tr(2) < tc(2)) .or. (tc(1) < tr(1))))      ) then
          write (gol,'("determined invalid constant interval:")'); call goErr
          call wrtgol( '  requested   : ', tr(1), ' - ', tr(2) ); call goErr
          call wrtgol( '  guessed     : ', tc(1), ' - ', tc(2) ); call goErr
          write (gol,'("  for tinterp : ",a)') md%tinterp; call goErr
          TRACEBACK; status=1; return
        end if
        ! read main field ?
        if ( md%filled1 ) then
          data1_read = md%tr1(1) /= tmid
        else
          data1_read = .true.
        end if
        ! read or leave ?
        if ( data1_read ) then
          data1_tref = tmid
          data1_t1   = tmid
          data1_t2   = tmid
        else
          ! field valid around requested interval, thus leave:
          !call goLabel();
          status=0; return
        end if
      !
      ! ** couple fields , valid for  3hr intervals  [00:00,03:00]  etc
      !    input filed valid for BEGIN of interval !
      !
      case ( 'cpl6', 'cpl3', 'cpl2', 'cpl1' )

        ! safety ...
        if ( direction < 0 ) then
          write (gol,'("check implementation for reverse tinterp ",a)') trim(md%tinterp); call goErr
          TRACEBACK; status=1; return
        end if

        ! extract time values for begin of current interval:
        call Get( tr(1), year, month, day, hour, minu )

        ! round hour to previous  baseh + 00/03/06/09/12/15/18/21
        hour = dth * floor(real(hour-baseh)/real(dth)) + baseh

        ! interval with constant field
        tc(1) = NewDate( year, month, day, hour )
        tc(2) = tc(1) + IncrDate(hour=dth)

        ! check for strange values:
        if ( (tr(1) < tc(1)) .or. (tc(2) < tr(1)) ) then
          write (gol,'("determined invalid first interval:")'); call goErr
          call wrtgol( '  requested   : ', tr(1), ' - ', tr(2) ); call goErr
          call wrtgol( '  guessed     : ', tc(1), ' - ', tc(2) ); call goErr
          write (gol,'("  for tinterp : ",a)') md%tinterp; call goErr
          TRACEBACK; status=1; return
       end if

        ! read primary field ?
        if ( md%filled1 ) then
          ! read new field if times are different:
          data1_read = (md%tr1(1) /= tc(1)) .or. (md%tr1(2) /= tc(1))
        else
          ! not filled yet, thus must read:
          data1_read = .true.
        end if

        ! read or leave ?
        if ( data1_read ) then
          data1_tref = tc(1)     ! begin of time interval
          data1_t1   = tc(1)
          data1_t2   = tc(1)
        end if
      !
      ! ** average fields , valid for  3hr intervals  [00:00,03:00]  etc
      !    average fields , valid for  3hr intervals  [00:00,06:00]  etc
      !
      case ( 'aver1', 'aver3', 'aver6', 'aver24', 'aver24_3' )

        ! extract time values for begin of current interval:
        call Get( tr(1), year, month, day, hour, minu )
        ! round hour to previous (next)  baseh + 00/03/06/09/12/15/18/21
        if ( direction > 0 ) then
          hour = dth *   floor(real(hour+minu/60.0-baseh)/real(dth)) + baseh
        else
          hour = dth * ceiling(real(hour+minu/60.0-baseh)/real(dth)) + baseh
        end if
        ! interval with constant field
        tc(1) = NewDate( year, month, day, hour )
        tc(2) = tc(1) + IncrDate(hour=direction*dth)
        ! check for strange values:
        if ( ((direction > 0) .and. ((tr(1) < tc(1)) .or. (tc(2) < tr(2)))) .or. &
             ((direction < 0) .and. ((tr(2) < tc(2)) .or. (tc(1) < tr(1))))      ) then
          write (gol,'("determined invalid first interval:")'); call goErr
          call wrtgol( '  requested   : ', tr(1), ' - ', tr(2) ); call goErr
          call wrtgol( '  guessed     : ', tc(1), ' - ', tc(2) ); call goErr
          write (gol,'("  for tinterp : ",a)') md%tinterp; call goErr
          TRACEBACK; status=1; return
        end if
        ! read primary field ?
        if ( md%filled1 ) then
          ! read new field if times are different:
          data1_read = (md%tr1(1) /= tc(1)) .or.  (md%tr1(2) /= tc(2))
        else
          ! not filled yet, thus must read:
          data1_read = .true.
        end if
        if ( data1_read ) then
          data1_tref = tc(1)
          data1_t1   = tc(1)
          data1_t2   = tc(2)
        end if

        ! setup reading of secondary data only if end of requested
        !  interval is later than primary interval:
        if ( (direction > 0) .and. (tc(2) < tr(2)) .or. &
             (direction < 0) .and. (tc(2) > tr(2))      ) then
          ! extract time values for end of requested interval:
          call Get( tr(2), year, month, day, hour, minu )
          ! round hour to next (previous)  baseh + 00/03/06/09/12/15/18/21
          if ( direction > 0 ) then
          hour = dth * floor(real(hour+minu/60.0-baseh)/real(dth)) + baseh
          else
            hour = dth * ceiling(real(hour+minu/60.0-baseh)/real(dth)) + baseh
          end if
          ! interval with constant field
          tc(1) = NewDate( year, month, day ) + IncrDate(hour=hour)
          tc(2) = tc(1) + IncrDate(hour=direction*dth)
          ! check for strange values:
          if ( ((direction > 0) .and. ((tr(1) < tc(1)) .or. (tc(2) < tr(2)))) .or. &
               ((direction < 0) .and. ((tr(2) < tc(2)) .or. (tc(1) < tr(1))))      ) then
            write (gol,'("determined invalid second interval:")'); call goErr
            call wrtgol( '  requested   : ', tr(1), ' - ', tr(2) ); call goErr
            call wrtgol( '  guessed     : ', tc(1), ' - ', tc(2) ); call goErr
            write (gol,'("  for tinterp : ",a)') md%tinterp; call goErr
            TRACEBACK; status=1; return
          end if
          ! read secondary field ?
          if ( md%filled2 ) then
            ! read new field if times are different;
            data2_read = (md%tr2(1) /= tc(1)) .or. (md%tr2(2) /= tc(2))
          else
            ! not filled yet, thus must read:
            data2_read = .true.
          end if
          if ( data2_read ) then
            data2_tref = tc(1)
            data2_t1   = tc(1)
            data2_t2   = tc(2)
          end if
        end if   ! tr partly after primary interval

      !
      ! ** interpolated between 6 hourly times 00/06/12/18
      !    interpolated between 6 hourly times 03/09/15/21
      !    interpolated between 3 hourly times 00/03/06/09/12/15/18/21
      !
      case ( 'interp6', 'interp6_3', 'interp3', 'interp2', 'interp1' )

        ! extract time values for begin of current interval:
        call Get( tr(1), year, month, day, hour, minu )
        ! truncate hour to previous 00/06/12/18, 03/09/15/21,
        !    or 00/03/06/09/12/15/18/21
        if ( direction > 0 ) then
          hour = dth *   floor(real(hour+minu/60.0-baseh)/real(dth)) + baseh
        else
          hour = dth * ceiling(real(hour+minu/60.0-baseh)/real(dth)) + baseh
        end if
        ! set begin of 3 or 6 hour interval:
        tprev = NewDate( year, month, day, hour )

        ! extract time values for end of current interval:
        call Get( tr(2), year, month, day, hour, minu )
        ! truncate hour to previous 00/06/12/18
        if ( direction > 0 ) then
          hour = dth * ceiling(real(hour+minu/60.0-baseh)/real(dth)) + baseh
        else
          hour = dth *   floor(real(hour+minu/60.0-baseh)/real(dth)) + baseh
        end if
        ! set end of 3 or 6 hour interval:
        tnext = NewDate( year, month, day, hour )

        !! testing ...
        !write (gol,'(a,": name ",a,"; tinterp ",a)') rname, trim(md%name), trim(md%tinterp); call goPr
        !call wrtgol( rname//':   tr(1),tr(2) ', tr(1), ', ', tr(2) ); call goPr
        !call wrtgol( rname//':   tprev,tnext ', tprev, ', ', tnext ); call goPr

        ! checks:
        !   [tprev,tnext] should be dth hours
        !   [tprev,tnext] should contain [tr(1),tr(2)]
        dth_int = abs( iTotal(tnext-tprev,'hour') )
        if ( ( (direction > 0) .and. ((tr(1) < tprev) .or. (tnext < tr(2))) ) .or. &
             ( (direction < 0) .and. ((tr(1) > tprev) .or. (tnext > tr(2))) ) .or. &
             ( (dth_int /= 0) .and. (dth_int /= dth) ) ) then
          write (gol,'("determined invalid interpolation interval:")'); call goErr
          call wrtgol( '  requested interval     : ', tr(1), ' - ', tr(2) ); call goErr
          call wrtgol( '  guesses input interval : ', tprev, ' - ', tnext ); call goErr
          write (gol,'("  dth, dth_int           : ",2f12.4)') dth, dth_int; call goErr
          write (gol,'("  temporal interpolation : ",a)') md%tinterp; call goErr
          TRACEBACK; status=1; return
        end if

        !
        !   .                         <-- previous field at dth hours
        !        o                    <-- latest interpolated field
        !             x               <-- target
        !                  o          <-- next field at dth hours
        !       tr1  tr   tr2
        ! --+--------------+------
        !  tprev          tnext
        !

        ! read main field ?
        if ( md%filled1 ) then
          ! md%data should be defined in [tprev,tr]
          if ( direction > 0 ) then
            data1_read = (md%tr1(1) < tprev) .or. (tr(2) < md%tr1(1))
          else
            data1_read = (md%tr1(1) > tprev) .or. (tr(2) > md%tr1(1))
          end if
        else
          data1_read = .true.
        end if
        if ( data1_read ) then
          data1_tref  = tprev
          data1_t1    = tprev
          data1_t2    = tprev
        end if

        !! testing ...
        !write (gol,'(a,":   data1_read ",l1)') rname, data1_read; call goPr
        !call wrtgol( rname//':   data1_t1,data_t2 ',data1_t1, ',', data1_t2 ); call goPr

        ! read second field ?
        if ( md%filled2 ) then
          ! md%data should be defined for tnext
          data2_read = md%tr2(1) /= tnext
        else
          data2_read = .true.
        end if
        if ( data2_read ) then
          data2_tref = tnext
          data2_t1   = tnext
          data2_t2   = tnext
        end if

        !! testing ...
        !write (gol,'(a,":   data2_read ",l1)') rname, data2_read; call goPr
        !call wrtgol( rname//':   data2_t1,data_t2 ',data2_t1, ',', data2_t2 ); call goPr

      !
      ! ** error ...
      !
      case default
        write (gol,'("unsupported time interpolation : ",a)') md%tinterp ; call goErr
        TRACEBACK; status=1; return
    end select


    !
    ! set ref times
    !

    !if ( fcmode ) then
    !
    !  ! in forecast mode, tfcday0 is 00:00 at the day the forecast starts;
    !  data1_tref = tfcday0
    !  data2_tref = tfcday0
    !
    !else

      ! dummy tref's : begin of day in which [data?_t1,data?_t2] starts:

      data1_tref = data1_t1
      if ( IsAnyDate(data1_tref) ) data1_tref = tr(1)
      call set( data1_tref, hour=0, min=0, sec=0, mili=0 )

      data2_tref = data2_t1
      if ( IsAnyDate(data2_tref) ) data2_tref = tr(1)
      call set( data2_tref, hour=0, min=0, sec=0, mili=0 )

    !end if


    !
    ! trap double reading
    !

    ! data1 already in data2 ?
    if ( data1_read .and. md%filled2 ) then
      if ( (data1_t1 == md%tr2(1)) .and. (data1_t2 == md%tr2(2)) ) then
        data1_read = .false.
        data1_copy = .true.
      end if
    end if

    ! data2 just read ?
    if ( data2_read .and. data1_read ) then
      ! data2 is same as data ?
      if ( (data2_tref == data1_tref) .and. &
           (data2_t1 == data1_t1) .and. (data2_t2 == data1_t2) ) then
        data2_read = .false.
        data2_copy = .true.
      end if
    end if

    !write (gol,'("SetupSetup:")'); call goPr
    !write (gol,'("  fcmode      : ",l1)') fcmode; call goPr
    !call wrtgol( '  tfcday0     : ', tfcday0 ); call goPr
    !write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goPr
    !call wrtgol( '  tr(1)       : ', tr(1) ); call goPr
    !call wrtgol( '  tr(2)       : ', tr(2) ); call goPr
    !write (gol,'("  1 read,copy : ",2l2)') data1_read, data1_copy; call goPr
    !call wrtgol( '  1 tref      : ', data1_tref ); call goPr
    !call wrtgol( '  1 t1        : ', data1_t1 ); call goPr
    !call wrtgol( '  1 t2        : ', data1_t2 ); call goPr
    !write (gol,'("  2 read,copy : ",2l2)') data2_read, data2_copy; call goPr
    !call wrtgol( '  2 tref      : ', data2_tref ); call goPr
    !call wrtgol( '  2 t1        : ', data2_t1 ); call goPr
    !call wrtgol( '  2 t2        : ', data2_t2 ); call goPr

    ! ok
    status = 0

    !call goLabel()

  end subroutine SetupSetup


  ! ***

  subroutine Setup_2d( md, tr, lli, nuv, status )

    use GO       , only : TDate, wrtgol
    use Grid     , only : TllGridInfo
    use TMM      , only : ReadField, Read_SP, Read_SR_OLS, WriteField
    use meteodata, only : TMeteoData, mdat_TimeInterpolation
    use ParTools , only : myid, root, Par_Broadcast

    ! --- in/out ------------------------------------

    type(TMeteoData), intent(inout)       ::  md
    type(TDate), intent(in)               ::  tr(2)
    type(TllGridInfo), intent(in)         ::  lli
    character(len=1), intent(in)          ::  nuv
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Setup_2d'

    ! --- local ----------------------------------
    logical                         ::  data1_read, data1_copy
    type(TDate)                     ::  data1_tref, data1_t1, data1_t2
    logical                         ::  data2_read, data2_copy
    type(TDate)                     ::  data2_tref, data2_t1, data2_t2

    ! --- begin -----------------------------

    !call goLabel(rname)

    ! leave if not in use:
    if ( .not. md%used ) then
      !call goLabel();
      status=0; return
    end if

    ! not changed by default
    md%changed = .false.

    !
    ! Get time interval valid for that met field, and read/copy flags
    !
    call SetupSetup( md, tr, &
                      data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
                      data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
                      status )
    IF_NOTOK_RETURN(status=1)

    !
    ! read/write field
    !
    ! read or copy primairy field ?
    if ( data1_read ) then

      if (myid==root) then ! only root does IO

          ! fill data:
          select case ( md%name )
             !
             ! special routine for surface pressure
             !
          case ( 'sp', 'sps' )

             call Read_SP( tmmd, md%sourcekey, trim(md%name), trim(md%unit), &
                  data1_tref, data1_t1, data1_t2, &
                  lli, md%data1(md%is(1):md%is(2),md%js(1):md%js(2),1), &
                  md%tmi1, status )
             IF_NOTOK_RETURN(status=1)

          !
          ! special routine for Olsson surface roughness:
          !
          case ( 'srols' )

            call Read_SR_OLS( tmmd, md%sourcekey, &
                            data1_tref, data1_t1, data1_t2, &
                            lli, md%data1(md%is(1):md%is(2),md%js(1):md%js(2),1), &
                            md%tmi1, status )
            IF_NOTOK_RETURN(status=1)

          !
          ! general field
          !
          case default

            call ReadField( tmmd, md%sourcekey, trim(md%name), trim(md%unit), &
                              data1_tref, data1_t1, data1_t2, &
                              lli, nuv, md%data1(md%is(1):md%is(2),md%js(1):md%js(2),1), &
                              md%tmi1, status )
            IF_NOTOK_RETURN(status=1)

        end select

        ! write meteofiles ?
        if ( md%putout ) then
          call WriteField( tmmd, md%destkey, &
                   md%tmi1, trim(md%name), trim(md%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, nuv, md%data1(md%is(1):md%is(2),md%js(1):md%js(2),1), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md%data1, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md%filled1  = .true.
      md%tr1(1)   = data1_t1
      md%tr1(2)   = data1_t2
      md%changed  = .true.

    else if ( data1_copy ) then

      ! copy data from secondary array:
      md%data1 = md%data2

      ! data array is filled now:
      md%filled1  = .true.
      md%tr1(1)   = data1_t1
      md%tr1(2)   = data1_t2
      md%changed  = .true.

    end if

    ! read or copy secondary field ?
    if ( data2_read ) then

      if (myid==root) then ! only root does IO

        ! fill data:
        select case ( md%name )
          !
          ! special routine for surface pressure
          !
          case ( 'sp', 'sps' )

            call Read_SP( tmmd, md%sourcekey, trim(md%name), trim(md%unit), &
                            data2_tref, data2_t1, data2_t2, &
                            lli, md%data2(md%is(1):md%is(2),md%js(1):md%js(2),1), &
                            md%tmi2, status )
            IF_NOTOK_RETURN(status=1)

          !
          ! general field
          !
          case default

            call ReadField( tmmd, md%sourcekey, trim(md%name), trim(md%unit), &
                              data2_tref, data2_t1, data2_t2, &
                              lli, nuv, md%data2(md%is(1):md%is(2),md%js(1):md%js(2),1), &
                              md%tmi2, status )
            IF_NOTOK_RETURN(status=1)

        end select

        ! write meteofiles ?
        if ( md%putout ) then
          call WriteField( tmmd, md%destkey, &
                   md%tmi2, trim(md%name), trim(md%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, nuv, md%data2(md%is(1):md%is(2),md%js(1):md%js(2),1), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md%data2, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md%filled2 = .true.
      md%tr2(1)  = data2_t1
      md%tr2(2)  = data2_t2

    else if ( data2_copy ) then

      ! copy data from secondary array:
      md%data2 = md%data1

      ! data array is filled now:
      md%filled2 = .true.
      md%tr2(1)  = data2_t1
      md%tr2(2)  = data2_t2

    end if

    !
    ! time interpolation
    !

    ! apply time interpolation:
    call mdat_TimeInterpolation( md, tr, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0

    !call goLabel()

  end subroutine Setup_2d


  ! ***


  subroutine Setup_3d( md, tr, lli, nuv, levi, nw, status )

    use GO       , only : TDate, wrtgol, operator(/=)
    use Grid     , only : TllGridInfo, TLevelInfo
    use TMM      , only : TMeteoInfo, ReadField, WriteField
    use meteodata, only : TMeteoData, mdat_TimeInterpolation
    use ParTools , only : myid, root, Par_Broadcast

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  md
    type(TDate), intent(in)               ::  tr(2)
    type(TllGridInfo), intent(in)         ::  lli
    character(len=1), intent(in)          ::  nuv
    type(TLevelInfo), intent(in)          ::  levi
    character(len=1), intent(in)          ::  nw
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Setup_3d'

    ! --- local ----------------------------------

    logical                         ::  data1_read, data1_copy
    type(TDate)                     ::  data1_tref, data1_t1, data1_t2
    logical                         ::  data2_read, data2_copy
    type(TDate)                     ::  data2_tref, data2_t1, data2_t2

    real, allocatable               ::  tmp_sp(:,:)

    ! --- begin -----------------------------

    !call goLabel(rname)

    !if (okdebug) then
    !   write(gol,'("    ",a,": ",a,l2)') rname, trim(md%name), md%used; call goPr
    !end if

    ! leave if not in use:
    if ( .not. md%used ) then
      !call goLabel()
      status=0; return
    end if

    ! not changed by default
    md%changed = .false.

    !
    ! time stuff
    !

    call SetupSetup( md, tr, &
                      data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
                      data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
                      status )
    IF_NOTOK_RETURN(status=1)

    !
    ! read/write field
    !

    ! primairy field ?
    if ( data1_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data1_t2 /= data1_t1 ) then
          !write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          !call wrtgol( '  data1_t1  : ', data1_t1  ); call goErr
          !call wrtgol( '  data1_t2  : ', data1_t2  ); call goErr
          !write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          !TRACEBACK; status=1; return
          write (gol,'("WARNING - using instant surface pressure for regridding temporal averaged 3D field ...")'); call goPr
        end if

        ! surface pressure field:
        allocate( tmp_sp(md%is(1):md%is(2),md%js(1):md%js(2)) )

        ! fill data:
        call ReadField( tmmd, md%sourcekey, trim(md%name), trim(md%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, nuv, levi, nw, &
                   tmp_sp, md%data1(md%is(1):md%is(2),md%js(1):md%js(2),md%ls(1):md%ls(2)), &
                   md%tmi1, status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( md%putout ) then
          call WriteField( tmmd, md%destkey, &
                   md%tmi1, 'sp', trim(md%name), trim(md%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, nuv, levi, nw, &
                   tmp_sp, md%data1(md%is(1):md%is(2),md%js(1):md%js(2),md%ls(1):md%ls(2)), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md%data1, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md%filled1  = .true.
      md%tr1(1)   = data1_t1
      md%tr1(2)   = data1_t2
      md%changed  = .true.

    else if ( data1_copy ) then

      ! copy data from secondary array:
      md%data1 = md%data2

      ! data array is filled now:
      md%filled1  = .true.
      md%tr1(1)   = data1_t1
      md%tr1(2)   = data1_t2
      md%changed  = .true.

    end if

    ! secondary field ?
    if ( data2_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data2_t2 /= data2_t1 ) then
          write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          call wrtgol( '  data2_t1  : ', data2_t1  ); call goErr
          call wrtgol( '  data2_t2  : ', data2_t2  ); call goErr
          write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure field:
        allocate( tmp_sp(md%is(1):md%is(2),md%js(1):md%js(2)) )

        ! fill data:
        call ReadField( tmmd, md%sourcekey, trim(md%name), trim(md%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, nuv, levi, nw, &
                   tmp_sp, md%data2(md%is(1):md%is(2),md%js(1):md%js(2),md%ls(1):md%ls(2)), &
                   md%tmi2, status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( md%putout ) then
          call WriteField( tmmd, md%destkey, &
                   md%tmi2, 'sp', trim(md%name), trim(md%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, nuv, levi, nw, &
                   tmp_sp, md%data2(md%is(1):md%is(2),md%js(1):md%js(2),md%ls(1):md%ls(2)), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md%data2, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md%filled2 = .true.
      md%tr2(1)  = data2_t1
      md%tr2(2)  = data2_t2

    else if ( data2_copy ) then

      ! copy data from secondary array:
      md%data2 = md%data1

      ! data array is filled now:
      md%filled2 = .true.
      md%tr2(1)  = data2_t1
      md%tr2(2)  = data2_t2

    end if


    !
    ! time interpolation
    !

    ! apply time interpolation:
    call mdat_TimeInterpolation( md, tr, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0

    !call goLabel()

  end subroutine Setup_3d


  ! **************************************************************
  ! ***
  ! *** mass fluxes
  ! ***
  ! **************************************************************


  subroutine Setup_MFUV( md_mfu, md_mfv, tr, lli, levi, status )

    use GO       , only : TDate, wrtgol, operator(/=)
    use Grid     , only : TllGridInfo, TLevelInfo
    use TMM      , only : TMeteoInfo, Read_MFUV, WriteField
    use meteodata, only : TMeteoData, mdat_TimeInterpolation
    use ParTools , only : myid, root, Par_Broadcast
    use Dims     , only : okdebug

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  md_mfu
    type(TMeteoData), intent(inout)       ::  md_mfv
    type(TDate), intent(in)               ::  tr(2)
    type(TllGridInfo), intent(in)         ::  lli
    type(TLevelInfo), intent(in)          ::  levi
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Setup_MFUV'

    ! --- local ----------------------------------

    logical                         ::  data1_read, data1_copy
    type(TDate)                     ::  data1_tref, data1_t1, data1_t2
    logical                         ::  data2_read, data2_copy
    type(TDate)                     ::  data2_tref, data2_t1, data2_t2

    real, allocatable               ::  tmp_spu(:,:)
    real, allocatable               ::  tmp_spv(:,:)

    ! --- begin -----------------------------

    !call goLabel(rname)

    ! testing ...
    if (okdebug) then
       write(gol,'("    ",a,": ",a,l2)') rname, "MFUV", md_mfu%used; call goPr
    end if

    ! leave if not in use:
    if ( md_mfu%used .neqv. md_mfv%used ) then
      write (gol,'("either none or both mfu and mfv should be in use")'); call goErr
      TRACEBACK; status=1; return
    end if
    if ( .not. md_mfu%used ) then
      !call goLabel();
      status=0; return
    end if

    ! not changed by default
    md_mfu%changed = .false.
    md_mfv%changed = .false.

    !
    ! time stuff
    !

    ! (sufficient to setup from mfu only)
    call SetupSetup( md_mfu, tr, &
                      data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
                      data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
                      status )
    IF_NOTOK_RETURN(status=1)

    !! testing ...
    !write (gol,'(a,": data1 read = ",l1,", copy = ",l2)') rname, data1_read, data1_copy; call goPr
    !call wrtgol( rname//':   t1,t2 = ', data1_t1, ',', data1_t2 ); call goPr
    !write (gol,'(a,": data2 read = ",l1,", copy = ",l2)') rname, data2_read, data2_copy; call goPr
    !call wrtgol( rname//':   t1,t2 = ', data2_t1, ',', data2_t2 ); call goPr

    !
    ! read/write field
    !

    ! primairy field ?
    if ( data1_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data1_t2 /= data1_t1 ) then
          write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          call wrtgol( '  data1_t1  : ', data1_t1  ); call goErr
          call wrtgol( '  data1_t2  : ', data1_t2  ); call goErr
          write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure field:
        allocate( tmp_spu(md_mfu%is(1)-1:md_mfu%is(2),md_mfu%js(1):md_mfu%js(2)  ) )
        allocate( tmp_spv(md_mfv%is(1)  :md_mfv%is(2),md_mfv%js(1):md_mfv%js(2)+1) )

        ! NOTE: strange old indexing:
        !   pu_tmpp  -->  pu(0:imr,1:jmr  ,1:lmr)  in  pu_t(0:imr+1,0:jmr+1,0:lmr)
        !   pv_tmpp  -->  pv(1:imr,1:jmr+1,1:lmr)  in  pv_t(0:imr+1,0:jmr+1,0:lmr)

        ! fill data:
        call Read_MFUV( tmmd, md_mfu%sourcekey, &
                 data1_tref, data1_t1, data1_t2, lli, levi, &
                 tmp_spu, &
                   md_mfu%data1(md_mfu%is(1)-1:md_mfu%is(2),&
                                md_mfu%js(1)  :md_mfu%js(2),&
                                md_mfu%ls(1)+1:md_mfu%ls(2)), &
                   md_mfu%tmi1, &
                 tmp_spv, &
                   md_mfv%data1(md_mfv%is(1)  :md_mfv%is(2)  , &
                                md_mfv%js(1)  :md_mfv%js(2)+1, &
                                md_mfv%ls(1)+1:md_mfv%ls(2)   ), &
                   md_mfv%tmi1, &
                 status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( md_mfu%putout ) then
          call WriteField( tmmd, md_mfu%destkey, &
                   md_mfu%tmi1, 'spu', trim(md_mfu%name), trim(md_mfu%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'u', levi, 'n', &
                   tmp_spu, md_mfu%data1(md_mfu%is(1)-1:md_mfu%is(2),&
                                         md_mfu%js(1)  :md_mfu%js(2),&
                                         md_mfu%ls(1)+1:md_mfu%ls(2)), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( md_mfv%putout ) then
          call WriteField( tmmd, md_mfv%destkey, &
                   md_mfv%tmi1, 'spv', trim(md_mfv%name), trim(md_mfv%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'v', levi, 'n', &
                   tmp_spv, md_mfv%data1(md_mfv%is(1)  :md_mfv%is(2)  , &
                                         md_mfv%js(1)  :md_mfv%js(2)+1, &
                                         md_mfv%ls(1)+1:md_mfv%ls(2)   ), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_spu )
        deallocate( tmp_spv )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md_mfu%data1, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( md_mfv%data1, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md_mfu%filled1  = .true.
      md_mfu%tr1(1)   = data1_t1
      md_mfu%tr1(2)   = data1_t2
      md_mfu%changed  = .true.
      md_mfv%filled1  = .true.
      md_mfv%tr1(1)   = data1_t1
      md_mfv%tr1(2)   = data1_t2
      md_mfv%changed  = .true.

    else if ( data1_copy ) then

      ! copy data from secondary array:
      md_mfu%data1 = md_mfu%data2
      md_mfv%data1 = md_mfv%data2

      ! data array is filled now:
      md_mfu%filled1  = .true.
      md_mfu%tr1(1)   = data1_t1
      md_mfu%tr1(2)   = data1_t2
      md_mfu%changed  = .true.
      md_mfv%filled1  = .true.
      md_mfv%tr1(1)   = data1_t1
      md_mfv%tr1(2)   = data1_t2
      md_mfv%changed  = .true.

    end if

    ! secondary field ?
    if ( data2_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data2_t2 /= data2_t1 ) then
          write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          call wrtgol( '  data2_t1  : ', data2_t1  ); call goErr
          call wrtgol( '  data2_t2  : ', data2_t2  ); call goErr
          write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure field:
        allocate( tmp_spu(md_mfu%is(1)-1:md_mfu%is(2),md_mfu%js(1):md_mfu%js(2)  ) )
        allocate( tmp_spv(md_mfv%is(1)  :md_mfv%is(2),md_mfv%js(1):md_mfv%js(2)+1) )

        ! NOTE: strange old indexing:
        !   pu_tmpp  -->  pu(0:imr,1:jmr  ,1:lmr)  in  pu_t(0:imr+1,0:jmr+1,0:lmr)
        !   pv_tmpp  -->  pv(1:imr,1:jmr+1,1:lmr)  in  pv_t(0:imr+1,0:jmr+1,0:lmr)

        ! fill data:
        call Read_MFUV( tmmd, md_mfu%sourcekey, &
                 data2_tref, data2_t1, data2_t2, lli, levi, &
                 tmp_spu, &
                   md_mfu%data2(md_mfu%is(1)-1:md_mfu%is(2),&
                                md_mfu%js(1)  :md_mfu%js(2),&
                                md_mfu%ls(1)+1:md_mfu%ls(2)), &
                   md_mfu%tmi2, &
                 tmp_spv, &
                   md_mfv%data2(md_mfv%is(1)  :md_mfv%is(2)  , &
                                md_mfv%js(1)  :md_mfv%js(2)+1, &
                                md_mfv%ls(1)+1:md_mfv%ls(2)   ), &
                   md_mfv%tmi2, &
                 status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( md_mfu%putout ) then
          call WriteField( tmmd, md_mfu%destkey, &
                   md_mfu%tmi2, 'spu', trim(md_mfu%name), trim(md_mfu%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'u', levi, 'n', &
                   tmp_spu, md_mfu%data2(md_mfu%is(1)-1:md_mfu%is(2),&
                                         md_mfu%js(1)  :md_mfu%js(2),&
                                         md_mfu%ls(1)+1:md_mfu%ls(2)), &
                   status )
          IF_NOTOK_RETURN(status=1)
        endif
        if ( md_mfv%putout ) then
          call WriteField( tmmd, md_mfv%destkey, &
                   md_mfv%tmi2, 'spv', trim(md_mfv%name), trim(md_mfv%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'v', levi, 'n', &
                   tmp_spv, md_mfv%data2(md_mfv%is(1)  :md_mfv%is(2)  , &
                                         md_mfv%js(1)  :md_mfv%js(2)+1, &
                                         md_mfv%ls(1)+1:md_mfv%ls(2)   ), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_spu )
        deallocate( tmp_spv )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md_mfu%data2, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( md_mfv%data2, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md_mfu%filled2 = .true.
      md_mfu%tr2(1)  = data2_t1
      md_mfu%tr2(2)  = data2_t2
      md_mfv%filled2 = .true.
      md_mfv%tr2(1)  = data2_t1
      md_mfv%tr2(2)  = data2_t2

    else if ( data2_copy ) then

      ! copy data from primary array:
      md_mfu%data2 = md_mfu%data
      md_mfv%data2 = md_mfv%data

      ! data array is filled now:
      md_mfu%filled2 = .true.
      md_mfu%tr2(1)  = data2_t1
      md_mfu%tr2(2)  = data2_t2
      md_mfv%filled2 = .true.
      md_mfv%tr2(1)  = data2_t1
      md_mfv%tr2(2)  = data2_t2

    end if


    !
    ! time interpolation
    !

    ! apply time interpolation:
    call mdat_TimeInterpolation( md_mfu, tr, status )
    IF_NOTOK_RETURN(status=1)
    call mdat_TimeInterpolation( md_mfv, tr, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0
    !call goLabel()

  end subroutine Setup_MFUV


  ! ***


  subroutine Setup_MFW( md_mfw, md_tsp, tr, lli, nuv, levi, nw, status )

    use GO       , only : TDate, wrtgol, operator(/=)
    use Grid     , only : TllGridInfo, TLevelInfo
    use TMM      , only : TMeteoInfo, ReadField, Read_MFW, WriteField
    use meteodata, only : TMeteoData, mdat_TimeInterpolation
    use ParTools , only : myid, root, Par_Broadcast
    use Dims     , only : okdebug

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  md_mfw
    type(TMeteoData), intent(inout)       ::  md_tsp
    type(TDate), intent(in)               ::  tr(2)
    type(TllGridInfo), intent(in)         ::  lli
    character(len=1), intent(in)          ::  nuv
    type(TLevelInfo), intent(in)          ::  levi
    character(len=1), intent(in)          ::  nw
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Setup_MFW'

    ! --- local ----------------------------------

    logical                         ::  data1_read, data1_copy
    type(TDate)                     ::  data1_tref, data1_t1, data1_t2
    logical                         ::  data2_read, data2_copy
    type(TDate)                     ::  data2_tref, data2_t1, data2_t2

    real, allocatable               ::  tmp_sp(:,:)

    ! --- begin -----------------------------

    !call goLabel(rname)

    if (okdebug) then
       write(gol,'("    ",a,": ",a,l2)') rname, trim(md_mfw%name), md_mfw%used; call goPr
    end if

    ! leave if not in use:
    if ( .not. md_mfw%used ) then
      !call goLabel();
      status=0; return
    end if

    ! error if tsp is not in use ...
    if ( .not. md_tsp%used ) then
      write (gol,'("mfw is in use but tsp not ..")'); call goErr
      !call goLabel();
      status=1; return
    end if

    ! not changed by default
    md_mfw%changed = .false.
    md_tsp%changed = .false.

    !
    ! time stuff
    !

    call SetupSetup( md_mfw, tr, &
                      data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
                      data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
                      status )
    IF_NOTOK_RETURN(status=1)

    !
    ! read/write field
    !

    ! primairy field ?
    if ( data1_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data1_t2 /= data1_t1 ) then
          write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          call wrtgol( '  data1_t1  : ', data1_t1  ); call goErr
          call wrtgol( '  data1_t2  : ', data1_t2  ); call goErr
          write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure field:
        allocate( tmp_sp(md_mfw%is(1):md_mfw%is(2),md_mfw%js(1):md_mfw%js(2)) )

        ! fill data:
        call Read_MFW( tmmd, md_mfw%sourcekey, &
                            data1_tref, data1_t1, data1_t2, &
                            lli, levi, &
                            tmp_sp, md_mfw%data1(md_mfw%is(1):md_mfw%is(2),md_mfw%js(1):md_mfw%js(2),md_mfw%ls(1):md_mfw%ls(2)), &
                            md_tsp%data1(md_tsp%is(1):md_tsp%is(2),md_tsp%js(1):md_tsp%js(2),1), &
                            md_mfw%tmi1, status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( md_mfw%putout ) then
          call WriteField( tmmd, md_mfw%destkey, &
                   md_mfw%tmi1, 'sp', trim(md_mfw%name), trim(md_mfw%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, nuv, levi, nw, &
                   tmp_sp, md_mfw%data1(md_mfw%is(1):md_mfw%is(2),md_mfw%js(1):md_mfw%js(2),md_mfw%ls(1):md_mfw%ls(2)), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( md_tsp%putout ) then
          ! use history from mfw ...
          call WriteField( tmmd, md_tsp%destkey, &
                   md_mfw%tmi1, trim(md_tsp%name), trim(md_tsp%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, nuv, md_tsp%data1(md_tsp%is(1):md_tsp%is(2),md_tsp%js(1):md_tsp%js(2),1), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md_mfw%data1, root, status )
      IF_NOTOK_RETURN(status=1)
      !
      call Par_Broadcast( md_tsp%data1, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md_mfw%filled1  = .true.
      md_mfw%tr1(1)   = data1_t1
      md_mfw%tr1(2)   = data1_t2
      md_mfw%changed  = .true.
      !
      md_tsp%filled1  = .true.
      md_tsp%tr1(1)   = data1_t1
      md_tsp%tr1(2)   = data1_t2
      md_tsp%changed  = .true.

    else if ( data1_copy ) then

      ! copy data from secondary array:
      md_mfw%data1 = md_mfw%data2

      ! data array is filled now:
      md_mfw%filled1  = .true.
      md_mfw%tr1(1)   = data1_t1
      md_mfw%tr1(2)   = data1_t2
      md_mfw%changed  = .true.
      !
      md_tsp%filled1  = .true.
      md_tsp%tr1(1)   = data1_t1
      md_tsp%tr1(2)   = data1_t2
      md_tsp%changed  = .true.

    end if

    ! secondary field ?
    if ( data2_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data2_t2 /= data2_t1 ) then
          write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          call wrtgol( '  data2_t1  : ', data2_t1  ); call goErr
          call wrtgol( '  data2_t2  : ', data2_t2  ); call goErr
          write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure field:
        allocate( tmp_sp(md_mfw%is(1):md_mfw%is(2),md_mfw%js(1):md_mfw%js(2)) )

        ! fill data:
        call Read_MFW( tmmd, md_mfw%sourcekey, &
                   data2_tref, data2_t1, data2_t2, &
                   lli, levi, &
                   tmp_sp, md_mfw%data2(md_mfw%is(1):md_mfw%is(2),md_mfw%js(1):md_mfw%js(2),md_mfw%ls(1):md_mfw%ls(2)), &
                   md_tsp%data2(md_tsp%is(1):md_tsp%is(2),md_tsp%js(1):md_tsp%js(2),1), &
                   md_mfw%tmi2, status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( md_mfw%putout ) then
          call WriteField( tmmd, md_mfw%destkey, &
                   md_mfw%tmi2, 'sp', trim(md_mfw%name), trim(md_mfw%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, nuv, levi, nw, &
                   tmp_sp, md_mfw%data2(md_mfw%is(1):md_mfw%is(2),md_mfw%js(1):md_mfw%js(2),md_mfw%ls(1):md_mfw%ls(2)), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( md_tsp%putout ) then
          ! use history from mfw ...
          call WriteField( tmmd, md_tsp%destkey, &
                   md_mfw%tmi2, trim(md_tsp%name), trim(md_tsp%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, nuv, &
                   md_tsp%data2(md_tsp%is(1):md_tsp%is(2),md_tsp%js(1):md_tsp%js(2),1), &
                   status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( md_mfw%data2, root, status )
      IF_NOTOK_RETURN(status=1)
      !
      call Par_Broadcast( md_tsp%data2, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      md_mfw%filled2 = .true.
      md_mfw%tr2(1)  = data2_t1
      md_mfw%tr2(2)  = data2_t2
      !
      md_tsp%filled2 = .true.
      md_tsp%tr2(1)  = data2_t1
      md_tsp%tr2(2)  = data2_t2

    else if ( data2_copy ) then

      ! copy data from secondary array:
      md_mfw%data2 = md_mfw%data1

      ! data array is filled now:
      md_mfw%filled2 = .true.
      md_mfw%tr2(1)  = data2_t1
      md_mfw%tr2(2)  = data2_t2
      !
      md_tsp%filled2 = .true.
      md_tsp%tr2(1)  = data2_t1
      md_tsp%tr2(2)  = data2_t2

    end if


    !
    ! time interpolation
    !

    ! apply time interpolation:
    call mdat_TimeInterpolation( md_mfw, tr, status )
    IF_NOTOK_RETURN(status=1)
    !
    call mdat_TimeInterpolation( md_tsp, tr, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0
    !call goLabel()

  end subroutine Setup_MFW


  ! ***


  subroutine Meteo_CheckPressure( n, status )

    use ParTools   , only : myid, root, ntracetloc
    !use ParTools   , only : Par_AllReduce
    use dims       , only : idate, newsrun
    use dims       , only : xcyc, im, jm
    use redgridZoom, only : calc_pdiff
    use io_hdf     , only : io_write2d_32d, DFACC_CREATE

    ! --- in/out -----------------------------

    integer, intent(in)       ::  n              ! region
    integer, intent(out)      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Meteo_CheckPressure'

    ! maximum accepted pressure difference:
    real, parameter  ::  pdiffmax_treshold = 1.0e2   ! Pa

    ! --- external -------------------------

    integer(4), external    ::  sfStart, sfEnd

    ! --- local -----------------------------

    real                    ::  pdiffmax, pdiffmax_l
    integer(4)              ::  io

    ! --- begin ------------------------------

    !call goLabel(rname)

    ! compare 'advected' pressure with read pressure
    if ( .not. newsrun ) then

      ! compute difference between 'advected' pressure sp and read pressure sp1;
      ! if no tracer have to be transported on this pe, use dummy:
      if ( ntracetloc == 0 ) then
        pdiffmax_l = 0.0
      else
        call calc_pdiff( n, sp1_dat(n)%data(:,:,1), &
                             sp_dat(n)%data(:,:,1), pdiffmax_l )
      end if

      ! compute maximum over all pe's:
      if (myid/=root) then
        write (gol,'("NO MPI YET!")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      pdiffmax = pdiffmax_l
      ! not yet ...
      !call Par_AllReduce( pdiffmax_l, pdiffmax, 'max', status )
      !IF_NOTOK_RETURN(status=1)

      ! check ...
      if ( pdiffmax > pdiffmax_treshold ) then
        write (gol,'("difference between advected and read-in pressure exceeds treshold :")'); call goErr
        write (gol,'("  max diff. : ",es12.4," [Pa]")') pdiffmax; call goErr
        write (gol,'("  treshold  : ",es12.4," [Pa]")') pdiffmax_treshold; call goErr
        write (gol,'("pressure arrays saved to local `pressure.hdf`")'); call goErr
        if ( myid == root ) then
          io = sfStart( 'pressure.hdf', DFACC_CREATE )
          if ( io > 0 ) then
            call io_write2d_32d( io, im(n)+4, 'LON', jm(n)+4, 'LAT', sp1_dat(n)%data(:,:,1), 'p'   , idate )
            call io_write2d_32d( io, im(n)+4, 'LON', jm(n)+4, 'LAT', sp_dat(n)%data(:,:,1), 'pold', idate )
            status = sfend(io)
          else
            write (gol,'("writing pressures")'); call goErr
          end if
        end if   ! root
        TRACEBACK; status=1; return
      end if   ! max diff

    end if  ! no newsrun

    ! ok
    status = 0
    !call goLabel()

  end subroutine Meteo_CheckPressure



  ! **************************************************************
  ! ***
  ! *** vertical velocity
  ! ***
  ! **************************************************************

  subroutine Compute_Omega( omega, lli, mfw, status )

    use binas    , only : grav
    use grid     , only : TllGridInfo, AreaOper
    use meteodata, only : TMeteoData
    use tmm      , only : SetHistory, AddHistory

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  omega          ! Pa/s downward
    type(TllGridInfo), intent(in)         ::  lli
    type(TMeteoData), intent(in)          ::  mfw            ! kg/s upward
    integer, intent(out)                  ::  status

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Compute_Omega'

    ! --- local ----------------------------------

    integer         ::  l

    ! --- begin ----------------------------------

    ! not in use ?
    if ( .not. omega%used ) return

    ! leave if not in use:
    if ( .not. mfw%used ) then
      write (gol,'("omega (Pa/s) requires mfw (kg/s)")'); call goErr
      TRACEBACK; status=1; return
    end if

    !call goLabel(rname)

    ! Pa/s  =  kg/s / m2 * g

    ! init with mass flux; revert sign from upward to downard, devide by gravity accelaration:
    omega%data = - mfw%data * grav   !  Pa/s m2

    ! loop over levels:
    do l = 1, size(omega%data,3)
      ! devide by cell area (m2) :
      call AreaOper( lli, omega%data(:,:,l), '/', 'm2', status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! info ..
    !call SetHistory( omega%tmi, mfw%tmi, status )
    !call AddHistory( omega%tmi, 'convert to Pa/s', status )

    ! ok
    status = 0
    !call goLabel()

  end subroutine Compute_Omega


  ! **************************************************************
  ! ***
  ! *** convective fluxes
  ! ***
  ! **************************************************************


  subroutine Setup_Convec( entu, entd, detu, detd, omega, gph, &
                                tr, lli, levi, status )

    use GO       , only : TDate, wrtgol, operator(/=), operator(>)
    use Grid     , only : TllGridInfo, TLevelInfo
    use TMM      , only : TMeteoInfo, Read_Convec, WriteField
    use meteodata, only : TMeteoData, mdat_TimeInterpolation
    use ParTools , only : myid, root, Par_Broadcast
    use Dims     , only : okdebug

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  entu, entd, detu, detd
    type(TMeteoData), intent(in)          ::  omega, gph
    type(TDate), intent(in)               ::  tr(2)
    type(TllGridInfo), intent(in)         ::  lli
    type(TLevelInfo), intent(in)          ::  levi
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Setup_Convec'

    ! --- local ----------------------------------

    logical                         ::  data1_read, data1_copy
    type(TDate)                     ::  data1_tref, data1_t1, data1_t2
    logical                         ::  data2_read, data2_copy
    type(TDate)                     ::  data2_tref, data2_t1, data2_t2
    type(TDate)                     ::  ttmp

    real, allocatable               ::  tmp_sp(:,:)

    ! --- begin -----------------------------

    !call goLabel(rname)
    if (okdebug) then
       write(gol,'("    ",a,": ",a,l2)') rname, "Conv.", entu%used; call goPr
    end if

    ! leave if not in use:
    if ( (.not. all((/entu%used,entd%used,detu%used,detd%used/)) ) &
         .and.  any((/entu%used,entd%used,detu%used,detd%used/))   ) then
      write (gol,'("either none or all of entu/entd/detu/detd should be in use")'); call goErr
      TRACEBACK; status=1; return
    end if
    if ( .not. entu%used ) then
      !call goLabel();
      status=0; return
    end if

    ! gph is required as input:
    if ( .not. gph%used ) then
      write (gol,'("gph should be in use to compute convective stuff from EC convective fluxes")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! omega is required as input:
    if ( .not. omega%used ) then
      write (gol,'("omega should be in use to compute convective stuff")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! not changed by default
    entu%changed = .false.
    entd%changed = .false.
    detu%changed = .false.
    detd%changed = .false.

    !
    ! time stuff
    !

    ! (sufficient to setup from entu only)
    call SetupSetup( entu, tr, &
                      data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
                      data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
                      status )
    IF_NOTOK_RETURN(status=1)

    !
    ! read/write field
    !

    ! primairy field ?
    if ( data1_read ) then

      ! swap reverse intervals, read/write in positive order:
      if ( data1_t1 > data1_t2 ) then
        ttmp = data1_t1
        data1_t1 = data1_t2
        data1_t2 = ttmp
      end if

      if (myid==root) then ! only root does IO

        !! safety check ...
        !if ( data1_t2 /= data1_t1 ) then
        !  !write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
        !  !call wrtgol( '  data1_t1  : ', data1_t1  ); call goErr
        !  !call wrtgol( '  data1_t2  : ', data1_t2  ); call goErr
        !  !write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
        !  !TRACEBACK; status=1; return
        !  write (gol,'("WARNING - convec for interval, but pressure/gph/etc instant ...")'); call goPr
        !end if

        ! surface pressure field:
        allocate( tmp_sp(entu%is(1):entu%is(2),entu%js(1):entu%js(2)) )

        ! fill data:
        call Read_Convec( tmmd, entu%sourcekey, &
                 data1_tref, data1_t1, data1_t2, lli, levi, &
                 omega%data, omega%tmi, &
                 gph%data, gph%tmi, &
                 tmp_sp, &
                 entu%data1, entu%tmi1, entd%data1, entd%tmi1, &
                 detu%data1, detu%tmi1, detd%data1, detd%tmi1, &
                 status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( entu%putout ) then
          call WriteField( tmmd, entu%destkey, &
                   entu%tmi1, 'sp', trim(entu%name), trim(entu%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, entu%data1, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( entd%putout ) then
          call WriteField( tmmd, entd%destkey, &
                   entd%tmi1, 'sp', trim(entd%name), trim(entd%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, entd%data1, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( detu%putout ) then
          call WriteField( tmmd, detu%destkey, &
                   detu%tmi1, 'sp', trim(detu%name), trim(detu%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, detu%data1, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( detd%putout ) then
          call WriteField( tmmd, detd%destkey, &
                   detd%tmi1, 'sp', trim(detd%name), trim(detd%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, detd%data1, status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( entu%data1, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( entd%data1, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( detu%data1, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( detd%data1, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
      entu%filled1  = .true.
      entu%tr1(1)   = data1_t1
      entu%tr1(2)   = data1_t2
      entu%changed  = .true.
      entd%filled1  = .true.
      entd%tr1(1)   = data1_t1
      entd%tr1(2)   = data1_t2
      entd%changed  = .true.
      detu%filled1  = .true.
      detu%tr1(1)   = data1_t1
      detu%tr1(2)   = data1_t2
      detu%changed  = .true.
      detd%filled1  = .true.
      detd%tr1(1)   = data1_t1
      detd%tr1(2)   = data1_t2
      detd%changed  = .true.

    else if ( data1_copy ) then

      ! copy data from secondary array:
      entu%data1 = entu%data2
      entd%data1 = entd%data2
      detu%data1 = detu%data2
      detd%data1 = detd%data2

      ! data array is filled now:
      entu%filled1  = .true.
      entu%tr1(1)   = data1_t1
      entu%tr1(2)   = data1_t2
      entu%changed  = .true.
      entd%filled1  = .true.
      entd%tr1(1)   = data1_t1
      entd%tr1(2)   = data1_t2
      entd%changed  = .true.
      detu%filled1  = .true.
      detu%tr1(1)   = data1_t1
      detu%tr1(2)   = data1_t2
      detu%changed  = .true.
      detd%filled1  = .true.
      detd%tr1(1)   = data1_t1
      detd%tr1(2)   = data1_t2
      detd%changed  = .true.

    end if

    ! secondary field ?
    if ( data2_read ) then

      ! swap reverse intervals, read/write in positive order:
      if ( data2_t1 > data2_t2 ) then
        ttmp = data2_t1
        data2_t1 = data2_t2
        data2_t2 = ttmp
      end if

      if (myid==root) then ! only root does IO

        !! safety check ...
        !if ( data2_t2 /= data2_t1 ) then
        !  !write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
        !  !call wrtgol( '  data2_t1  : ', data2_t1  ); call goErr
        !  !call wrtgol( '  data2_t2  : ', data2_t2  ); call goErr
        !  !write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
        !  !TRACEBACK; status=1; return
        !  write (gol,'("WARNING - convec for interval, but pressure/gph/etc instant ...")'); call goPr
        !end if

        ! surface pressure field:
        allocate( tmp_sp(entu%is(1):entu%is(2),entu%js(1):entu%js(2)) )

        ! fill data2:
        call Read_Convec( tmmd, entu%sourcekey, &
                 data2_tref, data2_t1, data2_t2, lli, levi, &
                 omega%data, omega%tmi, &
                 gph%data, gph%tmi, &
                 tmp_sp, &
                 entu%data2, entu%tmi2, entd%data2, entd%tmi2, &
                 detu%data2, detu%tmi2, detd%data2, detd%tmi2, &
                 status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( entu%putout ) then
          call WriteField( tmmd, entu%destkey, &
                   entu%tmi2, 'sp', trim(entu%name), trim(entu%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, entu%data2, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( entd%putout ) then
          call WriteField( tmmd, entd%destkey, &
                   entd%tmi2, 'sp', trim(entd%name), trim(entd%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, entd%data2, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( detu%putout ) then
          call WriteField( tmmd, detu%destkey, &
                   detu%tmi2, 'sp', trim(detu%name), trim(detu%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, detu%data2, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( detd%putout ) then
          call WriteField( tmmd, detd%destkey, &
                   detd%tmi2, 'sp', trim(detd%name), trim(detd%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'n', levi, '*', &
                   tmp_sp, detd%data2, status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( entu%data2, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( entd%data2, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( detu%data2, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( detd%data2, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data2 array is filled now:
      entu%filled2  = .true.
      entu%tr2(1)   = data2_t1
      entu%tr2(2)   = data2_t2
      entd%filled2  = .true.
      entd%tr2(1)   = data2_t1
      entd%tr2(2)   = data2_t2
      detu%filled2  = .true.
      detu%tr2(1)   = data2_t1
      detu%tr2(2)   = data2_t2
      detd%filled2  = .true.
      detd%tr2(1)   = data2_t1
      detd%tr2(2)   = data2_t2

    else if ( data2_copy ) then

      ! copy data2 from primary array:
      entu%data2 = entu%data1
      entd%data2 = entd%data1
      detu%data2 = detu%data1
      detd%data2 = detd%data1

      ! data2 array is filled now:
      entu%filled2  = .true.
      entu%tr2(1)   = data2_t1
      entu%tr2(2)   = data2_t2
      entd%filled2  = .true.
      entd%tr2(1)   = data2_t1
      entd%tr2(2)   = data2_t2
      detu%filled2  = .true.
      detu%tr2(1)   = data2_t1
      detu%tr2(2)   = data2_t2
      detd%filled2  = .true.
      detd%tr2(1)   = data2_t1
      detd%tr2(2)   = data2_t2

    end if


    !
    ! time interpolation
    !

    ! apply time interpolation:
    call mdat_TimeInterpolation( entu, tr, status )
    IF_NOTOK_RETURN(status=1)
    call mdat_TimeInterpolation( entd, tr, status )
    IF_NOTOK_RETURN(status=1)
    call mdat_TimeInterpolation( detu, tr, status )
    IF_NOTOK_RETURN(status=1)
    call mdat_TimeInterpolation( detd, tr, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0
    !call goLabel()

  end subroutine Setup_Convec


!  ! **************************************************************
!  ! ***
!  ! *** diffusive fluxes
!  ! ***
!  ! **************************************************************
!
!
!  subroutine Setup_Diffus( Kzz, tr, lli, levi, status )
!
!    use GO       , only : TDate, wrtgol, operator(/=)
!    use Grid     , only : TllGridInfo, TLevelInfo
!    use TMM      , only : TMeteoInfo, Read_Diffus, WriteField
!    use meteodata, only : TMeteoData, mdat_TimeInterpolation
!    use ParTools , only : myid, root, Par_Broadcast
!    use Dims     , only : okdebug
!
!    ! --- in/out ----------------------------------
!
!    type(TMeteoData), intent(inout)       ::  Kzz
!    type(TDate), intent(in)               ::  tr(2)
!    type(TllGridInfo), intent(in)         ::  lli
!    type(TLevelInfo), intent(in)          ::  levi
!    integer, intent(out)                  ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/Setup_Diffus'
!
!    ! --- local ----------------------------------
!
!    logical                         ::  data1_read, data1_copy
!    type(TDate)                     ::  data1_tref, data1_t1, data1_t2
!    logical                         ::  data2_read, data2_copy
!    type(TDate)                     ::  data2_tref, data2_t1, data2_t2
!
!    real, allocatable               ::  tmp_sp(:,:)
!
!    ! --- begin -----------------------------
!
!    ! not in use ?
!    if ( .not. Kzz%used ) then
!      status=0; return
!    end if
!
!    !call goLabel(rname)
!    if (okdebug) then
!       write(gol,'("    ",a,": ",a,l2)') rname, "Diffus", Kzz%used; call goPr
!    end if
!
!    ! not changed by default
!    Kzz%changed = .false.
!
!    !
!    ! time stuff
!    !
!
!    ! (sufficient to setup from entu only)
!    call SetupSetup( Kzz, tr, &
!                      data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
!                      data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
!                      status )
!    IF_NOTOK_RETURN(status=1)
!
!    !
!    ! read/write field
!    !
!
!    ! primairy field ?
!    if ( data1_read ) then
!
!      if (myid==root) then ! only root does IO
!
!        !! safety check ...
!        !if ( data1_t2 /= data1_t1 ) then
!        !  !write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
!        !  !call wrtgol( '  data1_t1  : ', data1_t1  ); call goErr
!        !  !call wrtgol( '  data1_t2  : ', data1_t2  ); call goErr
!        !  !write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
!        !  !TRACEBACK; status=1; return
!        !  write (gol,'("WARNING - convec for interval, but pressure/gph/etc instant ...")'); call goPr
!        !end if
!
!        ! surface pressure field:
!        allocate( tmp_sp(Kzz%is(1):Kzz%is(2),Kzz%js(1):Kzz%js(2)) )
!
!        ! fill data:
!        call Read_Diffus( tmmd, Kzz%sourcekey, &
!                 data1_tref, data1_t1, data1_t2, lli, levi, &
!                 tmp_sp, &
!                 Kzz%data1, Kzz%tmi1, &
!                 status )
!        IF_NOTOK_RETURN(status=1)
!
!        ! write meteofiles ?
!        if ( Kzz%putout ) then
!          call WriteField( tmmd, Kzz%destkey, &
!                   Kzz%tmi1, 'sp', trim(Kzz%name), trim(Kzz%unit), &
!                   data1_tref, data1_t1, data1_t2, &
!                   lli, 'n', levi, '*', &
!                   tmp_sp, Kzz%data1, status )
!          IF_NOTOK_RETURN(status=1)
!        end if
!
!        ! clear
!        deallocate( tmp_sp )
!
!      end if  ! root ?
!
!      ! send to other processors if necessary:
!      call Par_Broadcast( Kzz%data1, root, status )
!      IF_NOTOK_RETURN(status=1)
!
!      ! data array is filled now:
!      Kzz%filled1  = .true.
!      Kzz%tr1(1)   = data1_t1
!      Kzz%tr1(2)   = data1_t2
!      Kzz%changed  = .true.
!
!    else if ( data1_copy ) then
!
!      ! copy data from secondary array:
!      Kzz%data1 = Kzz%data2
!
!      ! data array is filled now:
!      Kzz%filled1  = .true.
!      Kzz%tr1(1)   = data1_t1
!      Kzz%tr1(2)   = data1_t2
!      Kzz%changed  = .true.
!
!    end if
!
!    ! secondary field ?
!    if ( data2_read ) then
!
!      if (myid==root) then ! only root does IO
!
!        !! safety check ...
!        !if ( data2_t2 /= data2_t1 ) then
!        !  !write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
!        !  !call wrtgol( '  data2_t1  : ', data2_t1  ); call goErr
!        !  !call wrtgol( '  data2_t2  : ', data2_t2  ); call goErr
!        !  !write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
!        !  !TRACEBACK; status=1; return
!        !  write (gol,'("WARNING - convec for interval, but pressure/gph/etc instant ...")'); call goPr
!        !end if
!
!        ! surface pressure field:
!        allocate( tmp_sp(Kzz%is(1):Kzz%is(2),Kzz%js(1):Kzz%js(2)) )
!
!        ! fill data2:
!        call Read_Diffus( tmmd, Kzz%sourcekey, &
!                 data2_tref, data2_t1, data2_t2, lli, levi, &
!                 tmp_sp, &
!                 Kzz%data2, Kzz%tmi2, &
!                 status )
!        IF_NOTOK_RETURN(status=1)
!
!        ! write meteofiles ?
!        if ( Kzz%putout ) then
!          call WriteField( tmmd, Kzz%destkey, &
!                   Kzz%tmi2, 'sp', trim(Kzz%name), trim(Kzz%unit), &
!                   data2_tref, data2_t1, data2_t2, &
!                   lli, 'n', levi, '*', &
!                   tmp_sp, Kzz%data2, status )
!          IF_NOTOK_RETURN(status=1)
!        end if
!
!        ! clear
!        deallocate( tmp_sp )
!
!      end if  ! root ?
!
!      ! send to other processors if necessary:
!      call Par_Broadcast( Kzz%data2, root, status )
!      IF_NOTOK_RETURN(status=1)
!
!      ! data2 array is filled now:
!      Kzz%filled2  = .true.
!      Kzz%tr2(1)   = data2_t1
!      Kzz%tr2(2)   = data2_t2
!
!    else if ( data2_copy ) then
!
!      ! copy data2 from primary array:
!      Kzz%data2 = Kzz%data1
!
!      ! data2 array is filled now:
!      Kzz%filled2  = .true.
!      Kzz%tr2(1)   = data2_t1
!      Kzz%tr2(2)   = data2_t2
!
!    end if
!
!
!    !
!    ! time interpolation
!    !
!
!    ! apply time interpolation:
!    call mdat_TimeInterpolation( Kzz, tr, status )
!    IF_NOTOK_RETURN(status=1)
!
!    !
!    ! done
!    !
!
!    ! ok
!    status = 0
!    !call goLabel()
!
!  end subroutine Setup_Diffus
  ! **************************************************************
  ! ***
  ! *** cloud cover
  ! ***
  ! **************************************************************


  subroutine Setup_CloudCovers( cc, cco, ccu, tr, lli, levi, status )

    use GO       , only : TDate, wrtgol, operator(/=)
    use Grid     , only : TllGridInfo, TLevelInfo
    use TMM      , only : TMeteoInfo, Read_CloudCovers, WriteField
    use meteodata, only : TMeteoData, mdat_TimeInterpolation
    use ParTools , only : myid, root, Par_Broadcast
    use Dims     , only : okdebug

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  cc, cco, ccu
    type(TDate), intent(in)               ::  tr(2)
    type(TllGridInfo), intent(in)         ::  lli
    type(TLevelInfo), intent(in)          ::  levi
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Setup_CloudCovers'

    ! --- local ----------------------------------

    logical                         ::  data1_read, data1_copy
    type(TDate)                     ::  data1_tref, data1_t1, data1_t2
    logical                         ::  data2_read, data2_copy
    type(TDate)                     ::  data2_tref, data2_t1, data2_t2

    real, allocatable               ::  tmp_sp(:,:)

    ! --- begin -----------------------------

    !call goLabel(rname)
    if (okdebug) then
       write(gol,'("    ",a,": ",a,l2)') rname, "Cloud Cover", cc%used; call goPr
    end if

    ! leave if not in use:
    if ( (.not. all((/cc%used,cco%used,ccu%used/)) ) .and. any((/cc%used,cco%used,ccu%used/)) ) then
      write (gol,'("either none or all of cc/cco/ccu should be in use")'); call goErr
      TRACEBACK; status=1; return
    end if
    if ( .not. cc%used ) then
      !call goLabel();
      status=0; return
    end if

    ! not changed by default
     cc%changed = .false.
    cco%changed = .false.
    ccu%changed = .false.

    !
    ! time stuff
    !

    ! (sufficient to setup from cc only)
    call SetupSetup( cc, tr, &
                      data1_read, data1_copy, data1_tref, data1_t1, data1_t2, &
                      data2_read, data2_copy, data2_tref, data2_t1, data2_t2, &
                      status )
    IF_NOTOK_RETURN(status=1)

    !
    ! read/write field
    !

    ! primairy field ?
    if ( data1_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data1_t2 /= data1_t1 ) then
          write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          call wrtgol( '  data1_t1  : ', data1_t1  ); call goErr
          call wrtgol( '  data1_t2  : ', data1_t2  ); call goErr
          write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure field:
        allocate( tmp_sp(cc%is(1):cc%is(2),cc%js(1):cc%js(2)) )

        ! fill data:
        call Read_CloudCovers( tmmd, cc%sourcekey, &
                 data1_tref, data1_t1, data1_t2, lli, levi, &
                 tmp_sp, cc%data1, cc%tmi1, &
                    cco%data1, cco%tmi1, ccu%data1, ccu%tmi1, &
                 status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if ( cc%putout ) then
          call WriteField( tmmd, cc%destkey, &
                   cc%tmi1, 'sp', trim(cc%name), trim(cc%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'n', levi, 'n', &
                   tmp_sp, cc%data1, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( cco%putout ) then
          call WriteField( tmmd, cco%destkey, &
                   cco%tmi1, 'sp', trim(cco%name), trim(cco%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'n', levi, 'n', &
                   tmp_sp, cco%data1, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( ccu%putout ) then
          call WriteField( tmmd, ccu%destkey, &
                   ccu%tmi1, 'sp', trim(ccu%name), trim(ccu%unit), &
                   data1_tref, data1_t1, data1_t2, &
                   lli, 'n', levi, 'n', &
                   tmp_sp, ccu%data1, status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast( cc%data1, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( cco%data1, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( ccu%data1, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data array is filled now:
       cc%filled1  = .true.
       cc%tr1(1)   = data1_t1
       cc%tr1(2)   = data1_t2
       cc%changed  = .true.
      cco%filled1  = .true.
      cco%tr1(1)   = data1_t1
      cco%tr1(2)   = data1_t2
      cco%changed  = .true.
      ccu%filled1  = .true.
      ccu%tr1(1)   = data1_t1
      ccu%tr1(2)   = data1_t2
      ccu%changed  = .true.

    else if ( data1_copy ) then

      ! copy data from secondary array:
       cc%data1 =  cc%data2
      cco%data1 = cco%data2
      ccu%data1 = ccu%data2

      ! data array is filled now:
       cc%filled1  = .true.
       cc%tr1(1)   = data1_t1
       cc%tr1(2)   = data1_t2
       cc%changed  = .true.
      cco%filled1  = .true.
      cco%tr1(1)   = data1_t1
      cco%tr1(2)   = data1_t2
      cco%changed  = .true.
      ccu%filled1  = .true.
      ccu%tr1(1)   = data1_t1
      ccu%tr1(2)   = data1_t2
      ccu%changed  = .true.

    end if

    ! secondary field ?
    if ( data2_read ) then

      if (myid==root) then ! only root does IO

        ! safety check ...
        if ( data2_t2 /= data2_t1 ) then
          write (gol,'("not sure that this routine is correct for time intervals:")'); call goErr
          call wrtgol( '  data2_t1  : ', data2_t1  ); call goErr
          call wrtgol( '  data2_t2  : ', data2_t2  ); call goErr
          write (gol,'("please deceide what to do with surface pressures ... ")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure field:
        allocate( tmp_sp(cc%is(1):cc%is(2),cc%js(1):cc%js(2)) )

        ! fill data2:
        call Read_CloudCovers( tmmd,  cc%sourcekey, &
                 data2_tref, data2_t1, data2_t2, lli, levi, &
                 tmp_sp, cc%data2,  cc%tmi2, &
                    cco%data2, cco%tmi2, ccu%data2, ccu%tmi2, &
                 status )
        IF_NOTOK_RETURN(status=1)

        ! write meteofiles ?
        if (  cc%putout ) then
          call WriteField( tmmd,  cc%destkey, &
                    cc%tmi2, 'sp', trim( cc%name), trim( cc%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'n', levi, 'n', &
                   tmp_sp,  cc%data2, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( cco%putout ) then
          call WriteField( tmmd, cco%destkey, &
                   cco%tmi2, 'sp', trim(cco%name), trim(cco%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'n', levi, 'n', &
                   tmp_sp, cco%data2, status )
          IF_NOTOK_RETURN(status=1)
        end if
        if ( ccu%putout ) then
          call WriteField( tmmd, ccu%destkey, &
                   ccu%tmi2, 'sp', trim(ccu%name), trim(ccu%unit), &
                   data2_tref, data2_t1, data2_t2, &
                   lli, 'n', levi, 'n', &
                   tmp_sp, ccu%data2, status )
          IF_NOTOK_RETURN(status=1)
        end if

        ! clear
        deallocate( tmp_sp )

      end if  ! root ?

      ! send to other processors if necessary:
      call Par_Broadcast(  cc%data2, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( cco%data2, root, status )
      IF_NOTOK_RETURN(status=1)
      call Par_Broadcast( ccu%data2, root, status )
      IF_NOTOK_RETURN(status=1)

      ! data2 array is filled now:
       cc%filled2  = .true.
       cc%tr2(1)   = data2_t1
       cc%tr2(2)   = data2_t2
      cco%filled2  = .true.
      cco%tr2(1)   = data2_t1
      cco%tr2(2)   = data2_t2
      ccu%filled2  = .true.
      ccu%tr2(1)   = data2_t1
      ccu%tr2(2)   = data2_t2

    else if ( data2_copy ) then

      ! copy data2 from primary array:
       cc%data2 =  cc%data1
      cco%data2 = cco%data1
      ccu%data2 = ccu%data1

      ! data2 array is filled now:
       cc%filled2  = .true.
       cc%tr2(1)   = data2_t1
       cc%tr2(2)   = data2_t2
      cco%filled2  = .true.
      cco%tr2(1)   = data2_t1
      cco%tr2(2)   = data2_t2
      ccu%filled2  = .true.
      ccu%tr2(1)   = data2_t1
      ccu%tr2(2)   = data2_t2

    end if


    !
    ! time interpolation
    !

    ! apply time interpolation:
    call mdat_TimeInterpolation( cc, tr, status )
    IF_NOTOK_RETURN(status=1)
    call mdat_TimeInterpolation( cco, tr, status )
    IF_NOTOK_RETURN(status=1)
    call mdat_TimeInterpolation( ccu, tr, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0
    !call goLabel()

  end subroutine Setup_CloudCovers


  ! **************************************************************
  ! ***
  ! *** compute
  ! ***
  ! **************************************************************


  !
  ! sp -> phlb -> mass
  !
  ! NOTE: assume that halo cells in sp have been filled correctly ...
  !

  subroutine Pressure_to_Mass( region, status )

    use Binas, only : grav
    use Grid , only : HPressure
    !use Grid , only : FillMass
    use Grid , only : AreaOper
    use dims , only : im, jm, lm
    use dims , only : xcyc
    use TM5_Geometry, only : lli, levi

    ! --- in/out ----------------------------------

    integer, intent(in)                   ::  region
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Pressure_to_Mass'

    ! --- local ----------------------------------

    integer         ::  imr, jmr, lmr
    integer         ::  l

    ! --- begin ----------------------------------

    ! local grid size:
    imr = im(region) ;  jmr = jm(region) ;  lmr = lm(region)

    ! fill pressure boundaries (Pa):
    if ( phlb_dat(region)%used ) then
      call HPressure( levi, sp_dat(region)%data(1:imr,1:jmr,1), &
                              phlb_dat(region)%data(1:imr,1:jmr,:), status )
      IF_NOTOK_RETURN(status=0)
    end if

    ! fill air mass (kg):
    if ( m_dat(region)%used ) then
      !call FillMass( m_dat(region)%data(1:imr,1:jmr,:), lli(region), levi, &
      !                  sp_dat(region)%data(1:imr,1:jmr,1), status )
      !IF_NOTOK_RETURN(status=0)
      ! pressure difference between top and bottom of layer:
      do l = 1, lmr
        m_dat(region)%data(:,:,l) = phlb_dat(region)%data(:,:,l) - phlb_dat(region)%data(:,:,l+1)  ! Pa
      end do
      ! convert to kg/m2 :
      m_dat(region)%data = m_dat(region)%data / grav  ! Pa/g = kg/m2
      ! convert to kg :
      call AreaOper( lli(region), m_dat(region)%data(1:imr,1:jmr,:), '*', 'm2', status )   ! kg
      IF_NOTOK_RETURN(status=0)
    end if

    ! ok
    status = 0

  end subroutine Pressure_to_Mass


  ! ***


  subroutine FillHalo_pm( region, status )

    use dims , only : im, jm
    use dims , only : xcyc

    ! --- in/out ----------------------------------

    integer, intent(in)                   ::  region
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/FillHalo_pm'

    ! --- local ----------------------------------

    integer         ::  imr, jmr

    ! --- begin ----------------------------------

    ! local grid size:
    imr = im(region) ;  jmr = jm(region)

    ! cyclic boundary conditions ?
    if ( xcyc(region) == 1 ) then
      ! fill western halo cells with east-most cells:
        sp_dat(region)%data(   -1:0    ,:,1) =   sp_dat(region)%data(imr-1:imr,:,1)
      phlb_dat(region)%data(   -1:0    ,:,:) = phlb_dat(region)%data(imr-1:imr,:,:)
         m_dat(region)%data(   -1:0    ,:,:) =    m_dat(region)%data(imr-1:imr,:,:)
      ! fill eastern halo cells with west-most cells:
        sp_dat(region)%data(imr+1:imr+2,:,1) =   sp_dat(region)%data(    1:2  ,:,1)
      phlb_dat(region)%data(imr+1:imr+2,:,:) = phlb_dat(region)%data(    1:2  ,:,:)
         m_dat(region)%data(imr+1:imr+2,:,:) =    m_dat(region)%data(    1:2  ,:,:)
    end if

    ! ok
    status = 0

  end subroutine FillHalo_pm


  ! ***


  ! mass -> phlb -> sp

  subroutine Mass_to_Pressure( region, status )

    use Binas       , only : grav
    use Grid        , only : AreaOper
    use dims        , only : im, jm, lm
    use TM5_Geometry, only : lli

    ! --- in/out ----------------------------------

    integer, intent(in)                   ::  region
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Mass_to_Pressure'

    ! --- local ----------------------------------

    integer         ::  imr, jmr, lmr
    integer         ::  l

    ! --- begin ----------------------------------

    ! local grid size:
    imr = im(region) ;  jmr = jm(region);  lmr = lm(region)

    ! other meteo required:
    if ( (.not. m_dat(region)%used) .or. (.not. phlb_dat(region)%used) &
         .or. (.not. sp_dat(region)%used) ) then
      write (gol,'("mass_to_pressure requires m, phlb, and sp")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! fill pressure at half level boundaries:
    !  o zero in space:
    phlb_dat(region)%data(:,:,lmr+1) = 0.0  ! kg m/s2 = Pa m2
    !  o add for each level pressure gradient:
    do l = lmr, 1, -1
      phlb_dat(region)%data(1:imr,1:jmr,l) = phlb_dat(region)%data(1:imr,1:jmr,l+1) &
                               + m_dat(region)%data(1:imr,1:jmr,l)*grav  ! kg m/s2 = Pa m2
    end do
    ! devide by grid cell area:
    call AreaOper( lli(region), phlb_dat(region)%data(1:imr,1:jmr,:), '/', 'm2', status )  ! Pa
    IF_NOTOK_RETURN(status=0)

    ! copy surface pressure:
    sp_dat(region)%data(1:imr,1:jmr,1) = phlb_dat(region)%data(1:imr,1:jmr,1)  ! Pa

    ! ok
    status = 0

  end subroutine Mass_to_Pressure


  ! ***

  subroutine Compute_PClim( region, status )

    use binas,       only: T0, p0, xm_air, grav, Rgas

    ! --- in/out ----------------------------------

    integer, intent(in)                   ::  region
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Compute_PClim'

    ! inv. scale height
    real, parameter                :: c = xm_air/(Rgas*T0)

    ! --- local ----------------------------------

    ! --- begin -----------------------------

    ! leave if not in use:
    if ( .not. pclim_dat(region)%used ) then
      !write (*,'(" pclim_dat not in use")')
      status=0; return
    end if

    ! only if input changed:
    if ( oro_dat(region)%changed ) then

      ! derive 'climatological' surface pressure
      ! Note: oro is h*grav --> scale height does not contain grav
      pclim_dat(region)%data = p0 * exp( - c * oro_dat(region)%data )

    end if  ! changed

    ! ok
    status = 0

  end subroutine Compute_PClim


  ! ***


  subroutine compute_gph( region, status )

    use Dims,        only : cdebug, kdebug, okdebug, kmain
    use Dims,        only : idate, itau, revert, newsrun
    use Dims,        only : im, jm, lm
    use Dims,        only : at, bt
    use binas,       only : grav
    use global_data, only : mass_dat
    use datetime,    only : tstamp
    use ParTools,    only : myid, root

    ! --- in/out ----------------------------------

    integer, intent(in)                   ::  region
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/compute_gph'

    ! --- local ----------------------------------

    real,dimension(:,:,:),pointer        :: gph, t, q
    real,dimension(:,:,:),pointer        :: ps

    integer                              :: i,j,l
    real                                 :: tv,pdown,pup

    ! --- begin -----------------------------

    ! leave if not in use:
    if ( .not. gph_dat(region)%used ) then
      !write (*,'(" gph_dat not in use")')
      status=0; return
    end if

    ! other meteo required:
    if ( (.not. temper_dat(region)%used) .or. (.not. humid_dat(region)%used) &
         .or. (.not. spm_dat(region)%used) .or. (.not. oro_dat(region)%used)) then
      write (gol,'("computation of gph for region ",i0," requires:")') region; call goErr
      write (gol,'("  temper (",l1,"), humid (",l1,"), spm (",l1,"), and oro (",l1,")")') &
          temper_dat(region)%used, humid_dat(region)%used, spm_dat(region)%used, oro_dat(region)%used; call goErr
      TRACEBACK; status=1; return
    end if

    ! leave if input did not change:
    if ( (.not.    spm_dat(region)%changed) .and. &
         (.not. temper_dat(region)%changed) .and. &
         (.not.  humid_dat(region)%changed) ) then
      !write (gol,'(a,": because ",3(l2,x))') rname, temper_dat(region)%changed, humid_dat(region)%changed, sp_dat(region)%changed; call goErr
      return
    end if

    ! field will be changed ...
    gph_dat(region)%changed = .true.

    ! pointers to meteo field
    ps  =>    spm_dat(region)%data
    t   => temper_dat(region)%data
    q   =>  humid_dat(region)%data
    gph =>    gph_dat(region)%data

    if ( myid==root .and. cdebug  ) call tstamp(kdebug,itau,rname)
    if ( myid==root .and. okdebug ) call tstamp(kmain ,itau,rname)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! compute geo potential height
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! surface altitude, oro is stored in g*m
    gph(:,:,1) = oro_dat(region)%data(:,:,1)/grav
    ! loop over cells:
    do l = 1, lm(region)-1
       do j = 1, jm(region)
          do i = 1, im(region)
            ! virtual temperature:
            tv = t(i,j,l) * ( 1.0 + 0.608*q(i,j,l) )
            ! pressure levels:
            pdown = at(l  ) + bt(l  )*ps(i,j,1)
            pup   = at(l+1) + bt(l+1)*ps(i,j,1)
             ! rgas in different units!
             gph(i,j,l+1) = gph(i,j,l) + tv * 287.05 * alog(pdown/pup) / grav
             ! note dec 2002 (MK) gph now from 1--->lm+1
          end do
       end do
    end do
    !set top of atmosphere at 200 km !note_mk: now in lm+1
    gph(:,:,lm(region)+1) = 200000.0

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! done
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! done
    nullify( ps )
    nullify( t )
    nullify( q )
    nullify( gph )

    ! ok
    status = 0

  end subroutine compute_gph



end module meteo
