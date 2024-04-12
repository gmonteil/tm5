!#################################################################
!
! Input/output of meteofiles : prism version.
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_PRISM_ERROR_RETURN(action) if (status/=0) then; call prism_error(status,gol); call goErr; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module tmm_mf_prism

  use GO       , only : gol, goPr, goErr, goBug, goLabel
  use GO       , only : TDate
  use TM5_Prism, only : TshRemap
  
  implicit none
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  TMeteoFile_prism
  public  ::  mfPrism_Init, mfPrism_Done
  public  ::  Init, Done
  public  ::  ReadRecord

!  public  ::  ifs_tmin, ifs_tmid, ifs_tplus
!  public  ::  ifs_T, ifs_shn, ifs_lm
!  public  ::  ifs_lsp_tmin, ifs_lsp_tmid, ifs_lsp_tplus
!  public  ::  ifs_vor_tmid, ifs_div_tmid
!  public  ::  ifs_temper_tmid

  
  ! --- const ------------------------------
  
  character(len=*), parameter  ::  mname = 'module tmm_mf_prism'
  
  !--- type ---------------------------------
  
  type TMeteoFile_prism
    ! reference date:
    type(TDate)            ::  tday
  end type TMeteoFile_prism
  
  
  ! --- interfaces -------------------------
  
  interface Init
    module procedure mf_Init
  end interface

  interface Done
    module procedure mf_Done
  end interface

  interface ReadRecord
    module procedure mf_ReadRecord
  end interface


  ! --- local -----------------------------------
  
  type TCache
    type(TDate)              ::  tmid
    real, pointer            ::  data(:,:,:)
  end type TCache
  
  integer, parameter         ::  ncache = 4
  integer, parameter         ::  icache_sp   = 1
  integer, parameter         ::  icache_lnsp = 2
  integer, parameter         ::  icache_vor  = 3
  integer, parameter         ::  icache_div  = 4
  type(TCache), save         ::  cache(ncache)

  type(TshRemap), save       ::  shRemap2d, shRemap3d



contains


  ! =============================================================
  
  !
  ! Init prism stuff.
  !
  
  subroutine mfPrism_Init( status )
  
    use GO       , only : NewDate
    use TM5_Prism, only : ifs_nlev, ifs_shn
    use TM5_Prism, only : ifs_nlon, ifs_nlat
    use TM5_Prism, only : Init
  
    ! --- in/out -----------------------------------------
    
    integer, intent(out)         ::  status
    
    ! --- const -----------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/mfPrism_Init'
    
    ! --- local -----------------------------------------
    
    integer             ::  icache
    
    ! --- begin -----------------------------------------
    
    !
    ! cache
    !
    
    icache = icache_sp
    cache(icache)%tmid = NewDate(year=9999)
    allocate( cache(icache)%data(ifs_nlon,ifs_nlat,1) )
    
    icache = icache_lnsp
    cache(icache)%tmid = NewDate(year=9999)
    allocate( cache(icache)%data(2,ifs_shn,1) )
    
    icache = icache_vor
    cache(icache)%tmid = NewDate(year=9999)
    allocate( cache(icache)%data(2,ifs_shn,ifs_nlev) )
    
    icache = icache_div
    cache(icache)%tmid = NewDate(year=9999)
    allocate( cache(icache)%data(2,ifs_shn,ifs_nlev) )
    
    !
    ! spectral remapping
    !
    
    call Init( shRemap2d, status )
    IF_NOTOK_RETURN(status=1)  

    call Init( shRemap3d, status )
    IF_NOTOK_RETURN(status=1)  

    ! ok
    status = 0

  end subroutine mfPrism_Init
    
    
  ! ***
  
  
  subroutine mfPrism_Done( status )
  
    use TM5_Prism, only : Done
  
    ! --- in/out -----------------------------------------
    
    integer, intent(out)         ::  status
    
    ! --- const -----------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/mfPrism_Done'
    
    ! --- begin -----------------------------------------
  
    ! clear:
    deallocate( cache(icache_sp  )%data )
    deallocate( cache(icache_lnsp)%data )
    deallocate( cache(icache_vor )%data )
    deallocate( cache(icache_div )%data )
    
    ! done:
    call Done( shRemap2d, status )
    IF_NOTOK_RETURN(status=1)  
    call Done( shRemap3d, status )
    IF_NOTOK_RETURN(status=1)  

    ! ok
    status = 0

  end subroutine mfPrism_Done


  ! ********************************************************************************
  ! ***
  ! *** definitions
  ! ***
  ! ********************************************************************************
  
  
  
!  subroutine mfprism_GetInVarID( fieldname, id, status )
!  
!    ! --- in/out -----------------------------------------
!    
!    character(len=*), intent(in)   ::  fieldname
!    integer, intent(out)           ::  id
!    integer, intent(out)           ::  status
!    
!    ! --- const -----------------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/mfprism_GetInVarID'
!    
!    ! --- begin -----------------------------------------
!    
!    ! select var_id for this input field:
!    select case ( fieldname )
!      case ( 'sst' ); id = 1
!      case default
!        write (gol,'("unsupported field name `",a,"`")') fieldname; call goErr
!        write (gol,'("in ",a)') rname; call goErr; status=1; return
!    end select
!    
!    ! ok
!    status = 0
!  
!  end subroutine mfprism_GetInVarID
!
!
!  ! ***
!  
!  
!  subroutine mfprism_GetInVar( var_id, status, ggN )
!  
!    ! --- in/out -----------------------------------------
!    
!    integer, intent(in)              ::  var_id
!    integer, intent(out)             ::  status
!    
!    integer, intent(out), optional   ::  ggN
!    
!    ! --- const -----------------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/mfprism_GetInVar'
!    
!    ! --- begin -----------------------------------------
!    
!    ! return gg size ?
!    if ( present(ggN) ) then
!      write (gol,'("ggN not implemented yet ...")'); call goErr
!      write (gol,'("in ",a)') rname; call goErr; status=1; return
!      ggN = 0
!    end if
!    
!    ! ok
!    status = 0
!  
!  end subroutine mfprism_GetInVar
!
!
!  ! ***
!  
!
!  subroutine mfprism_GetField( fieldname, itap_sec, gg, status )
!  
!    ! --- in/out -----------------------------------------
!    
!    character(len=*), intent(in)   ::  fieldname
!    integer, intent(in)            ::  itap_sec
!    real, intent(out)              ::  gg(:)
!    integer, intent(out)           ::  status
!    
!    ! --- const -----------------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/mfprism_GetField'
!
!    ! --- begin -----------------------------------------
!    
!!    ! select var_id for this input field:
!!    select case ( fieldname )
!!
!!      case ( 'sst' )
!!      
!!        ! receive field ...
!!        call prism_get_proto (il_var_id_in(var_id), itap_sec, sst, status)
!!        write (il_mparout,fmt=*) 'itap_sec, statussst=',itap_sec, status
!!        if ( (status /= PRISM_ok) .and. (status < PRISM_recvd) ) then
!!          write (il_mparout,'("Pb in reading ",A8)') cl_read(var_id)
!!          write (il_mparout,'("Time is ",I8)') itap_sec
!!          write (il_mparout,'("Error code is ",I2)') status
!!          call prism_abort_proto(il_comp_id, 'tm5.f90','abort3') 
!!        end if
!!
!!        ! type conversion ...
!!        gg = sst
!!
!!      case default
!!
!!        write (il_mparout,'("ERROR - unsupported field name `",a,"`")') fieldname
!!        call prism_abort_proto(il_comp_id, 'tm5.f90','abort3') 
!!
!!    end select
!    gg = 0.0
!  
!    ! ok
!    status = 0
!  
!  end subroutine mfprism_GetField



  ! ==============================================================
  

  subroutine mf_Init( mf, tday, status )
  
    use GO, only : TDate
  
    ! --- in/out ----------------------------
    
    type(TMeteoFile_prism), intent(out)   ::  mf
    type(TDate), intent(in)               ::  tday
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/mf_Init'
    
    ! --- begin --------------------------------
    
    ! store reference date:
    mf%tday = tday
    
!    ! select var id in coupler for requested field:
!    call PRISM_GetInVarID( paramkey, mf%id, status )
!    if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if

    ! ok
    status = 0
 
  end subroutine mf_Init
  
  
  ! ***
  
  
  subroutine mf_Done( mf, status )
  
    ! --- in/out ------------------------------------

    type(TMeteoFile_prism), intent(inout) ::  mf
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/mf_Done'
    
    ! --- begin -------------------------------------
    
    ! ok
    status = 0
    
  end subroutine mf_Done
  


  ! ***
  

!  subroutine mf_ReadRecord( mf, paramkey, t1, t2, nuv, nw, &
!                                gridtype, levi, &
!                                lli, ll, sp_ll, &
!                                ggi, gg, sp_gg, &
!                                shi, sh, lnsp_sh, &
!                                tmi, status )
!
!    use parray    , only : pa_Init, pa_Done
!    use GO        , only : TDate
!    use Grid      , only : TLevelInfo
!    use Grid      , only : TllGridInfo, TggGridInfo, TshGridInfo
!    use tmm_info  , only : TMeteoInfo, AddHistory
!
!#ifdef pgf90
!    use parray_r4
!    use parray_r8
!    use parray_c4
!    use parray_c8
!#endif
!
!    ! --- in/out -------------------------------
!    
!    type(TMeteoFile_prism), intent(inout)  ::  mf
!    character(len=*), intent(in)           ::  paramkey
!    type(TDate), intent(in)                ::  t1, t2
!    character(len=1), intent(in)           ::  nuv
!    character(len=1), intent(in)           ::  nw
!
!    character(len=2), intent(out)          ::  gridtype
!    type(TLevelInfo), intent(out)          ::  levi
!
!    type(TllGridInfo), intent(inout)       ::  lli
!    real, pointer                          ::  ll(:,:,:)
!    real, pointer                          ::  sp_ll(:,:)
!
!    type(TggGridInfo), intent(inout)       ::  ggi
!    real, pointer                          ::  gg(:,:)
!    real, pointer                          ::  sp_gg(:)
!
!    type(TshGridInfo), intent(inout)       ::  shi
!    complex, pointer                       ::  sh(:,:)
!    complex, pointer                       ::  lnsp_sh(:)
!
!    type(TMeteoInfo), intent(out)          ::  tmi
!    integer, intent(out)                   ::  status
!    
!    ! --- const --------------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord'
!    
!    ! --- local -----------------------------------
!    
!    real, pointer         ::  ll2(:,:,:)
!    real, pointer         ::  gg2(:,:)
!    complex, pointer      ::  sh2(:,:)
!    
!    ! --- begin ---------------------------------
!
!    ! combined field ?
!    select case ( paramkey )
!    
!      ! *** surface stress
!      
!      case ( 'sstr' )
!      
!        ! read first field:
!        call mf_ReadRecord_1( mf, 'ewss', t1, t2, nuv, nw, &
!                                gridtype, levi, &
!                                lli, ll, sp_ll, &
!                                ggi, gg, sp_gg, &
!                                shi, sh, lnsp_sh, &
!                                tmi, status )
!        IF_NOTOK_RETURN(status=1)  
!
!        ! init pointer:
!        call pa_Init( ll2 )
!        call pa_Init( gg2 )
!        call pa_Init( sh2 )
!
!        ! read second field:
!        call mf_ReadRecord_1( mf, 'nsss', t1, t2, nuv, nw, &
!                                gridtype, levi, &
!                                lli, ll2, sp_ll, &
!                                ggi, gg2, sp_gg, &
!                                shi, sh2, lnsp_sh, &
!                                tmi, status )
!         IF_NOTOK_RETURN(status=1)      
!
!        ! process:
!        select case ( gridtype )
!          case ( 'll' ) ; ll = sqrt( ll**2 + ll2**2 )
!          case ( 'gg' ) ; gg = sqrt( gg**2 + gg2**2 )
!          case default
!            write (gol,'("unsupported gridtype for surface stress :",a)') gridtype; call goErr
!            write (gol,'("in ",a)') rname; call goErr; status=1; return
!        end select   
!        call AddHistory( tmi, 'sstr=sqrt(ewss**2+nsss**2)', status )
!
!        ! clear pointers:
!        call pa_Done( ll2 )
!        call pa_Done( gg2 )
!        call pa_Done( sh2 )
!
!      ! *** default
!      
!      case default
!      
!        call mf_ReadRecord_1( mf, paramkey, t1, t2, nuv, nw, &
!                                gridtype, levi, &
!                                lli, ll, sp_ll, &
!                                ggi, gg, sp_gg, &
!                                shi, sh, lnsp_sh, &
!                                tmi, status )
!         IF_NOTOK_RETURN(status=1)
!      
!    end select
!    
!    ! ok
!    status = 0
!
!  end subroutine mf_ReadRecord
  

  ! ***


  subroutine mf_ReadRecord( mf, paramkey, t1, t2, nuv, nw, &
                                  gridtype, levi, &
                                  lli, ll, sp_ll, &
                                  ggi, gg, sp_gg, &
                                  shi, sh, lnsp_sh, &
                                  tmi, status )

    use parray  , only : pa_SetShape
#ifdef pgf90
    use parray_r4
    use parray_r8
    use parray_c4
    use parray_c8
#endif
    use GO      , only : TDate, wrtgol, IsAnyDate
    use GO      , only : operator(+), operator(-), operator(/), operator(==)
    use Grid    , only : TllGridInfo, TggGridInfo, TshGridInfo, TLevelInfo
    use Grid    , only : Init, Set
    use tmm_info, only : TMeteoInfo, Init, AddHistory
    use binas   , only : grav

    use PRISM    , only : PRISM_Time_Struct
    use PRISM    , only : PRISM_Jobstart_date
    use PRISM    , only : prism_calc_newdate
    !use PRISM    , only : prism_get      ! compile error : not public entity of module
    use PRISM    , only : PRISM_CPL, PRISM_CPLIO

    use TM5_Prism, only : exchange_period

    use TM5_Prism, only : PrsmGrid_2d_sfc, PrsmGrid_2d
    use TM5_Prism, only : ifs_nlon, ifs_nlat, ifs_nlev
    use TM5_Prism, only : ifs_nlon_sfc, ifs_nlat_sfc
    use TM5_Prism, only : ifs_shT, ifs_shn, ifs_shn_recv  
    use TM5_Prism, only : SetPrismTime, wrtgol
    use TM5_Prism, only : var_ids, var_names
    use TM5_Prism, only : ivar_ctm_spinf3d, ivar_ctm_spvor, ivar_ctm_spdiv
    use TM5_Prism, only : ivar_ctm_spinf2d, ivar_ctm_spsp
!    use TM5_Prism, only : ivar_ctm_slon, ivar_ctm_lat, ivar_ctm_lev
    use TM5_Prism, only : ivar_ctm_sp, ivar_ctm_t, ivar_ctm_q
    use TM5_Prism, only : ivar_ctm_sshf, ivar_ctm_slhf
    use TM5_Prism, only : ivar_ctm_ewss, ivar_ctm_nsss
    use TM5_Prism, only : ivar_ctm_oro, ivar_ctm_lsm
    use TM5_Prism, only : Setup, Remap
    
!use binas   , only : p_global
!use tm5_prism, only : comp_id
!use prism, only : prism_abort
!use file_hdf
!type(TTimeSeriesHDF)   :: F

    ! --- in/out -------------------------------
    
    type(TMeteoFile_prism), intent(inout)  ::  mf
    character(len=*), intent(in)           ::  paramkey
    type(TDate), intent(in)                ::  t1, t2
    character(len=1), intent(in)           ::  nuv
    character(len=1), intent(in)           ::  nw

    character(len=2), intent(out)          ::  gridtype
    type(TLevelInfo), intent(out)          ::  levi

    type(TllGridInfo), intent(inout)       ::  lli
    real, pointer                          ::  ll(:,:,:)
    real, pointer                          ::  sp_ll(:,:)

    type(TggGridInfo), intent(inout)       ::  ggi
    real, pointer                          ::  gg(:,:)
    real, pointer                          ::  sp_gg(:)

    type(TshGridInfo), intent(inout)       ::  shi
    complex, pointer                       ::  sh(:,:)
    complex, pointer                       ::  lnsp_sh(:)

    type(TMeteoInfo), intent(out)          ::  tmi
    integer, intent(out)                   ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord'
    
    ! value filled into spinf arrays to define not a number:    
    real, parameter           ::  spinf_nan = -1.2345

    ! --- local -----------------------------------
    
    type(TDate)               ::  tmid
    type(PRISM_Time_Struct)   ::  model_time
    type(PRISM_Time_Struct)   ::  model_time_bounds(2)
    real                      ::  tm5_dt_sec
    
    real, allocatable         ::  sh_raw(:,:,:)
    real, allocatable         ::  sh_zcl(:,:,:)

    integer                   ::  ivar
    integer                   ::  info
    
    real                      ::  tb_sec
    
    integer                   ::  icache
    
    ! --- begin ---------------------------------
    
    call goLabel(rname)
    
    ! no times defined in t1 and t2 ?
    if ( IsAnyDate(t1) .and. IsAnyDate(t2) ) then
      ! for constant fields (orography), t1 and t2 are any date;
      ! use the tday stored in mf structure for the orography time:
      tmid = mf%tday
    else
      ! mid time of [t1,t2]
      tmid = t1 + (t2-t1)/2
    end if

    ! convert from tm5 time structure to prism time structure:
    call SetPrismTime( model_time, tmid, status )
    IF_ERROR_RETURN(status=1)
    
    ! time bounds : half of exchange period
    tb_sec = exchange_period * 3600.0 / 2.0   ! sec
    !
    model_time_bounds(1) = model_time
    call prism_calc_newdate ( model_time_bounds(1), -1.0*tb_sec, status )
    IF_PRISM_ERROR_RETURN(status=1)
    !
    model_time_bounds(2) = model_time
    call prism_calc_newdate ( model_time_bounds(2), +1.0*tb_sec, status )
    IF_PRISM_ERROR_RETURN(status=1)
    
    ! info ...
    write (gol,'("  paramkey          : ",a)') trim(paramkey); call goPr
    call wrtgol( '  t2                : ', t2   ); call goPr
    call wrtgol( '  tmid              : ', tmid ); call goPr
    call wrtgol( '  model_time        : ', model_time        ); call goPr
    call wrtgol( '  model_time_bounds : ', model_time_bounds ); call goPr
    
    !
    ! return spectral fields:
    !
    
    select case ( paramkey )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 3D lat/lon fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'T', 'Q' )
      
!call Init( F, '/var/tmp/nas/run.TM5/tmm-out/test-'//trim(paramkey)//'.hdf', 40, status )

        ! return lat/lon field:
        gridtype = 'll'

        ! intialize lat/lon info:
        !call Init( lli, -180.0+0.5*(360.0/ifs_nlon), (360.0/ifs_nlon), ifs_nlon, &
        !                 -90.0+0.5*(180.0/ifs_nlat), (180.0/ifs_nlat), ifs_nlat, status )
        call Init( lli, PrsmGrid_2d(1)%lli_west , PrsmGrid_2d(1)%lli_dlon, PrsmGrid_2d(1)%lli_nlon, &
                        PrsmGrid_2d(1)%lli_south, PrsmGrid_2d(1)%lli_dlat, PrsmGrid_2d(1)%lli_nlat, status )
        IF_NOTOK_RETURN(status=1)
        
        ! allocate arrays if necessary:
        call pa_SetShape( sp_ll, ifs_nlon, ifs_nlat )
        call pa_SetShape( ll, ifs_nlon, ifs_nlat, ifs_nlev )
        
        ! ~~ 3d field ~~
        
        ! receive field:
        select case ( paramkey )
          case ( 'T' ) ; ivar = ivar_ctm_t
          case ( 'Q' ) ; ivar = ivar_ctm_q
          case default
            write (gol,'("ivar not implemented for paramkey `",a,"`")') trim(paramkey); call goErr
            write (gol,'("in ",a)') rname; status=1; return; call goErr
        end select
        !
        call prism_get ( var_ids(ivar), model_time, model_time_bounds, ll, info, status )
        IF_PRISM_ERROR_RETURN(status=1)
        select case ( info )
          !case ( -1 )
          !  write (gol,'("prism_get : no ",a," for this time")') trim(var_names(ivar)); call goErr
          !  write (gol,'("in ",a)') rname; call goErr; status=1; return
          case ( PRISM_CPL, PRISM_CPLIO )
            write (gol,'("    prism_get : received ",a)') trim(var_names(ivar)); call goPr
            write (gol,*) '    min,max : ', minval(ll), maxval(ll); call goPr
          case default
            write (gol,'("unknown info from prism_get : ",i6)') info; call goPr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
        
!call AddRecord( F, 't_raw', 'prism_get', 'K', 'real(4)', (/60,45,19/), tmp_ll, status )

        ! ~~ surface pressure ~~
        
        ! id's        
        icache = icache_sp

        ! in cache ?
        if ( cache(icache)%tmid == tmid ) then
        
          ! copy from cache:
          sp_ll = cache(icache)%data(:,:,1)
          
        else
        
          ! receive field:
          ivar = ivar_ctm_sp
          call prism_get ( var_ids(ivar), model_time, model_time_bounds, sp_ll, info, status )
          IF_PRISM_ERROR_RETURN(status=1)
          select case ( info )
            !case ( -1 )
            !  write (gol,'("    prism_get : no ",a," for this time")') trim(var_names(ivar)); call goPr
            case ( PRISM_CPL, PRISM_CPLIO )
              write (gol,'("    prism_get : received ",a)') trim(var_names(ivar)); call goPr
              write (gol,*) '    min,max : ', minval(sp_ll), maxval(sp_ll); call goPr
            case default
              write (gol,'("unknown info from prism_get : ",i6)') info; call goPr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select

          ! store in chache:
          cache(icache)%tmid = tmid
          cache(icache)%data(:,:,1) = sp_ll
          
        end if

!call AddRecord( F, var_names(ivar), 'prism_get', 'unit', 'real(4)', (/60,45/), tmp_ll(:,:,1), status )


        ! ~~ levels ~~

        ! level info
        select case ( ifs_nlev )
          case ( 19 )
            call Init( levi, 'ec19', status )
            IF_NOTOK_RETURN(status=1)
          case ( 60 )
            call Init( levi, 'ec60', status )
            IF_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported ifs nlev : ",i4)') ifs_nlev; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
        
!call Done( F, status )


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2D surface fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'sshf', 'slhf', 'ewss', 'nsss', 'oro', 'lsm' )
      
        ! return lat/lon field:
        gridtype = 'll'

        ! intialize lat/lon info:
        !call Init( lli, -180.0+0.5*(360.0/ifs_nlon_sfc), (360.0/ifs_nlon_sfc), ifs_nlon_sfc, &
        !                 -90.0+0.5*(180.0/ifs_nlat_sfc), (180.0/ifs_nlat_sfc), ifs_nlat_sfc, status )
        call Init( lli, PrsmGrid_2d_sfc%lli_west , PrsmGrid_2d_sfc%lli_dlon, PrsmGrid_2d_sfc%lli_nlon, &
                        PrsmGrid_2d_sfc%lli_south, PrsmGrid_2d_sfc%lli_dlat, PrsmGrid_2d_sfc%lli_nlat, status )
        IF_NOTOK_RETURN(status=1)
        
        ! allocate arrays if necessary:
        call pa_SetShape( ll, ifs_nlon_sfc, ifs_nlat_sfc, 1 )
        
        ! receive field:
        select case ( paramkey )
          case ( 'sshf' ) ; ivar = ivar_ctm_sshf
          case ( 'slhf' ) ; ivar = ivar_ctm_slhf
          case ( 'ewss' ) ; ivar = ivar_ctm_ewss
          case ( 'nsss' ) ; ivar = ivar_ctm_nsss
          case ( 'oro'  ) ; ivar = ivar_ctm_oro
          case ( 'lsm'  ) ; ivar = ivar_ctm_lsm
          case default
            write (gol,'("ivar not implemented for paramkey `",a,"`")') trim(paramkey); call goErr
            write (gol,'("in ",a)') rname; status=1; return; call goErr
        end select
        !
        call prism_get ( var_ids(ivar), model_time, model_time_bounds, ll(:,:,1), info, status )
        IF_PRISM_ERROR_RETURN(status=1)
        select case ( info )
          !case ( -1 )
          !  write (gol,'("prism_get : no ",a," for this time")') trim(var_names(ivar)); call goErr
          !  write (gol,'("in ",a)') rname; call goErr; status=1; return
          case ( PRISM_CPL, PRISM_CPLIO )
            write (gol,'("    prism_get : received ",a)') trim(var_names(ivar)); call goPr
            write (gol,*) '    min,max : ', minval(ll(:,:,1)), maxval(ll(:,:,1)); call goPr
          case default
            write (gol,'("unknown info from prism_get : ",i6)') info; call goPr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
        
        ! convert ...
        select case ( paramkey )
          !case ( 'oro'  ) ; ll = ll/grav      ! m m/s2  ->  m   <-- oro is in m*m/s2
          case ( 'lsm'  ) ; ll = ll * 100.0   ! 0-1     ->  %
        end select
        
        ! dummy levels
        call Init( levi, 1, (/0.0,0.0/), (/1.0,0.0/), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2D spectral fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'LNSP' )
      
        ! return spectral field:
        gridtype = 'sh'

        ! intialize spherical harmonic field info:
        call Init( shi, ifs_shT, status )
        IF_NOTOK_RETURN(status=1)
        
        ! allocate array if necessary:
        call pa_SetShape( sh, ifs_shn, 1 )
        
        ! raw spectral field from coupler:
        allocate( sh_raw(ifs_shn_recv*2,1,1) )
        allocate( sh_zcl(2,ifs_shn,1) )

        ! setup remapping if not done yet ...
        if ( .not. shRemap2d%filled ) then

          ! fill with nan; values unchanged if receiving array is larger than necessary:
          sh_raw = spinf_nan
    
          ! receive field:
          ivar = ivar_ctm_spinf2d
          call prism_get ( var_ids(ivar), model_time, model_time_bounds, sh_raw, info, status )
          IF_PRISM_ERROR_RETURN(status=1)
          select case ( info )
            !case ( -1 )
            !  write (gol,'("    prism_get : no ",a," for this time")') trim(var_names(ivar)); call goPr
            case ( PRISM_CPL, PRISM_CPLIO )
              write (gol,'("    prism_get : received ",a)') trim(var_names(ivar)); call goPr
            case default
              write (gol,'("unknown info from prism_get : ",i6)') info; call goPr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select

          ! setup remapping:      
          call Setup( shRemap2d, sh_raw, spinf_nan, status )
          IF_NOTOK_RETURN(status=1)
          
        end if

        ! id's        
        icache = icache_lnsp

        ! in cache ?
        if ( cache(icache)%tmid == tmid ) then

          ! copy from cache:
          sh_zcl = cache(icache)%data
          
        else
        
          ! receive field:
          ivar = ivar_ctm_spsp
          call prism_get ( var_ids(ivar), model_time, model_time_bounds, sh_raw, info, status )
          IF_PRISM_ERROR_RETURN(status=1)
          select case ( info )
            !case ( -1 )
            !  write (gol,'("    prism_get : no ",a," for this time")') trim(var_names(ivar)); call goPr
            case ( PRISM_CPL, PRISM_CPLIO )
              write (gol,'("    prism_get : received ",a)') trim(var_names(ivar)); call goPr
            case default
              write (gol,'("unknown info from prism_get : ",i6)') info; call goPr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select

          ! remap ...
          call Remap( shRemap2d, sh_raw, shi, sh_zcl, status )
          IF_NOTOK_RETURN(status=1)
        
          ! store in chache:
          cache(icache)%tmid = tmid
          cache(icache)%data = sh_zcl
          
        end if

        ! convert to complex:
        sh = cmplx(sh_zcl(1,:,:),sh_zcl(2,:,:))
        
        ! clear
        deallocate( sh_raw )
        deallocate( sh_zcl )
        
        ! dummy levels
        call Init( levi, 1, (/0.0,0.0/), (/1.0,0.0/), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 3D spectral fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      ! NOTE: no lnsp provided, we should get rid of this ...

      case ( 'VO', 'D' )
      
        ! return spectral field:
        gridtype = 'sh'

        ! intialize spherical harmonic field info:
        call Init( shi, ifs_shT, status )
        IF_NOTOK_RETURN(status=1)
        
        ! allocate arrays if necessary:
        call pa_SetShape( sh, ifs_shn, ifs_nlev )
        
        ! raw spectral field from coupler:
        allocate( sh_raw(ifs_shn_recv*2,ifs_nlev,1) )
        allocate( sh_zcl(2,ifs_shn,ifs_nlev) )

        ! setup remapping if not done yet ...
        if ( .not. shRemap3d%filled ) then

          ! fill with nan; values unchanged if receiving array is larger than necessary:
          sh_raw = spinf_nan
    
          ! receive field:
          ivar = ivar_ctm_spinf3d
          call prism_get ( var_ids(ivar), model_time, model_time_bounds, sh_raw, info, status )
          IF_PRISM_ERROR_RETURN(status=1)
          select case ( info )
            !case ( -1 )
            !  write (gol,'("    prism_get : no ",a," for this time")') trim(var_names(ivar)); call goPr
            case ( PRISM_CPL, PRISM_CPLIO )
              write (gol,'("    prism_get : received ",a)') trim(var_names(ivar)); call goPr
            case default
              write (gol,'("unknown info from prism_get : ",i6)') info; call goPr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select

          ! setup remapping:      
          call Setup( shRemap3d, sh_raw, spinf_nan, status )
          IF_NOTOK_RETURN(status=1)
          
        end if

        ! id's        
        select case ( paramkey )
          case ( 'VO' ) ; ivar = ivar_ctm_spvor ; icache = icache_vor
          case ( 'D'  ) ; ivar = ivar_ctm_spdiv ; icache = icache_div
          case default
            write (gol,'("unsupported 3d sh : ",a)') trim(paramkey); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

        ! in cache ?
        if ( cache(icache)%tmid == tmid ) then
        
          ! copy from cache:
          sh_zcl = cache(icache)%data
          
        else
        
          ! receive field:
          call prism_get ( var_ids(ivar), model_time, model_time_bounds, sh_raw, info, status )
          IF_PRISM_ERROR_RETURN(status=1)
          select case ( info )
            !case ( -1 )
            !  write (gol,'("    prism_get : no ",a," for this time")') trim(var_names(ivar)); call goPr
            case ( PRISM_CPL, PRISM_CPLIO )
              write (gol,'("    prism_get : received ",a)') trim(var_names(ivar)); call goPr
            case default
              write (gol,'("unknown info from prism_get : ",i6)') info; call goPr
              write (gol,'("in ",a)') rname; call goErr; status=1; return
          end select

          ! remap ...
          call Remap( shRemap3d, sh_raw, shi, sh_zcl, status )
          IF_NOTOK_RETURN(status=1)

          ! store in chache:
          cache(icache)%tmid = tmid
          cache(icache)%data = sh_zcl
          
        end if

        ! convert to complex:
        sh = cmplx(sh_zcl(1,:,:),sh_zcl(2,:,:))
        
        ! clear
        deallocate( sh_raw )
        deallocate( sh_zcl )
        
        ! level info
        select case ( ifs_nlev )
          case ( 19 )
            call Init( levi, 'ec19', status )
            IF_NOTOK_RETURN(status=1)
          case ( 60 )
            call Init( levi, 'ec60', status )
            IF_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported ifs nlev : ",i4)') ifs_nlev; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! not supported yet ...
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported paramkey `",a,"`")') trim(paramkey); call goErr
        write (gol,'("in ",a)') rname; status=1; return; call goErr

    end select
    
    ! fill some info values:
    call Init( tmi, paramkey, 'unkown', status )
    call AddHistory( tmi, 'model==oasis_coupler', status )

    call goLabel()

    ! ok
    status = 0

  end subroutine mf_ReadRecord
  

  
end module tmm_mf_prism
