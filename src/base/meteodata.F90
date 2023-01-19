!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module MeteoData

  use GO         , only : gol, goErr, goPr
  use go         , only : TDate
  use tmm        , only : TMeteoInfo
  use Dims       , only : nregions, nregions_all
  use os_specs   , only : MAX_RCKEY_LEN

  implicit none

  ! --- in/out -----------------------------------

  private

  public  :: TMeteoData

  public  :: Init, Done
  public  :: Set, Check
  public  :: Alloc
  public  :: SetData
  public  :: TimeInterpolation

  public  :: sp1_dat, sp2_dat, sp_dat, spm_dat
  public  :: tsp_dat
  public  :: phlb_dat, m_dat
  public  :: mfu_dat, mfv_dat, mfw_dat
  public  :: pu_dat,  pv_dat,  pw_dat
  public  :: temper_dat, humid_dat
  public  :: gph_dat
  public  :: omega_dat
  public  :: lwc_dat, iwc_dat, cc_dat, cco_dat, ccu_dat
  public  :: entu_dat, entd_dat, detu_dat, detd_dat
  public  :: oro_dat
  public  :: pclim_dat
  public  :: lsmask_dat
  public  :: sr_ecm_dat, sr_ols_dat
  public  :: u10m_dat, v10m_dat
  public  :: ewss_dat, nsss_dat
  public  :: blh_dat
  ! New for dry deposition
  public  :: ci_dat, sd_dat, swvl1_dat, slhf_dat, sshf_dat, src_dat, d2m_dat, t2m_dat
  public  :: ssr_dat, tv_dat, cvl_dat, cvh_dat, nveg, sf_dat


  ! --- const ----------------------------------

  ! module name
  character(len=*), parameter  ::  mname = 'MeteoData'

  ! number of surface types in ECMWF
  integer, parameter  ::  nveg = 20


  ! --- types -----------------------------------

  ! storage for single meteo field:

  type TMeteoData
    ! in use ?
    logical              ::  used
    ! changed ?
    logical              ::  changed
    ! description:
    character(len=16)    ::  name                  ! field name
    character(len=16)    ::  unit                  ! kg, K, ...
    ! time interpolation:
    character(len=10)    ::  tinterp               ! const6, interp3, ...
    ! shapes:
    integer              ::  is(2), js(2), ls(2)
    integer              ::  halo
    ! target data:
    real, pointer        ::  data(:,:,:)
    type(TDate)          ::  tr(2)                 ! timerange
    type(TMeteoInfo)     ::  tmi                   ! history info
    ! primary data:
    real, pointer        ::  data1(:,:,:)
    logical              ::  filled1
    type(TDate)          ::  tr1(2)                ! timerange
    type(TMeteoInfo)     ::  tmi1                  ! history info
    ! secondary data:
    real, pointer        ::  data2(:,:,:)
    logical              ::  filled2
    type(TDate)          ::  tr2(2)                ! timerange
    type(TMeteoInfo)     ::  tmi2                  ! history info
    ! input:
    character(len=MAX_RCKEY_LEN)   ::  sourcekey
    ! output
    logical              ::  putout
    character(len=MAX_RCKEY_LEN)   ::  destkey
  end type TMeteoData


  ! --- interfaces -----------------------------------

  interface Init
    module procedure mdat_Init
  end interface

  interface Done
    module procedure mdat_Done
  end interface

  interface Set
    module procedure mdat_Set
  end interface

  interface Check
    module procedure mdat_Check
  end interface

  interface Alloc
    module procedure mdat_Alloc
  end interface

  interface TimeInterpolation
    module procedure mdat_TimeInterpolation
  end interface

  interface SetData
    module procedure mdat_SetData
  end interface


  ! --- var ---------------------------------------------

  ! meteo fields:

  type(TMeteoData), save, target    ::      sp1_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      sp2_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       sp_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      spm_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      tsp_dat(1:nregions_all)

  type(TMeteoData), save, target    ::     phlb_dat(1:nregions_all)
  type(TMeteoData), save, target    ::        m_dat(1:nregions_all)

  type(TMeteoData), save, target    ::      mfu_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      mfv_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      mfw_dat(1:nregions_all)

  type(TMeteoData), save, target    ::       pu_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       pv_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       pw_dat(1:nregions_all)
  type(TMeteoData), save, target    ::    omega_dat(1:nregions_all)

  type(TMeteoData), save, target    ::   temper_dat(1:nregions_all)
  type(TMeteoData), save, target    ::    humid_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      gph_dat(1:nregions_all)

  type(TMeteoData), save, target    ::      lwc_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      iwc_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       cc_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      cco_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      ccu_dat(1:nregions_all)

  type(TMeteoData), save, target    ::     entu_dat(1:nregions_all)
  type(TMeteoData), save, target    ::     entd_dat(1:nregions_all)
  type(TMeteoData), save, target    ::     detu_dat(1:nregions_all)
  type(TMeteoData), save, target    ::     detd_dat(1:nregions_all)

  type(TMeteoData), save, target    ::        oro_dat(1:nregions_all)
  type(TMeteoData), save, target    ::     lsmask_dat(1:nregions_all)
  ! climatological surface pressure (derived from oro by uniform scale height)
  type(TMeteoData), save, target    ::      pclim_dat(1:nregions_all)
  !type(TMeteoData), save, target    ::     albedo_dat(1:nregions_all)
  type(TMeteoData), save, target    ::     sr_ecm_dat(1:nregions_all)
  type(TMeteoData), save, target    ::     sr_ols_dat(1:nregions_all)
  !type(TMeteoData), save, target    ::        sst_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       u10m_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       v10m_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       ewss_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       nsss_dat(1:nregions_all)
  !type(TMeteoData), save, target    ::        lsp_dat(1:nregions_all)
  !type(TMeteoData), save, target    ::         cp_dat(1:nregions_all)
  type(TMeteoData), save, target    ::        blh_dat(1:nregions_all)
!  type(TMeteoData), save, target    ::        blh_sfc
  ! New for dry deposition
  type(TMeteoData), save, target    ::         ci_dat(1:nregions_all)
  type(TMeteoData), save, target    ::         sf_dat(1:nregions_all)
  type(TMeteoData), save, target    ::         sd_dat(1:nregions_all)
  type(TMeteoData), save, target    ::      swvl1_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       sshf_dat(1:nregions_all)
  type(TMeteoData), save, target    ::       slhf_dat(1:nregions_all)
  type(TMeteoData), save, target    ::        src_dat(1:nregions_all)
  type(TMeteoData), save, target    ::        d2m_dat(1:nregions_all)
  type(TMeteoData), save, target    ::        t2m_dat(1:nregions_all)
  type(TMeteoData), save, target    ::        ssr_dat(1:nregions_all)
  type(TMeteoData), save, target    ::         tv_dat(1:nregions_all,nveg)
  type(TMeteoData), save, target    ::        cvl_dat(1:nregions_all)
  type(TMeteoData), save, target    ::        cvh_dat(1:nregions_all)


contains


  ! ==========================================================


  subroutine mdat_Init( md, name, unit, tinterp, is, js, halo, ls, &
                          sourcekey, putout, destkey, status )

    ! --- in/out -----------------------------------

    type(TMeteoData), intent(out)     ::  md
    character(len=*), intent(in)      ::  name, unit
    character(len=*), intent(in)      ::  tinterp
    integer, intent(in)               ::  is(2), js(2)
    integer, intent(in)               ::  halo
    integer, intent(in)               ::  ls(2)
    character(len=*), intent(in)      ::  sourcekey
    logical, intent(in)               ::  putout
    character(len=*), intent(in)      ::  destkey
    integer, intent(out)              ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/mdat_Init'

    ! --- begin --------------------------------

    ! not in use:
    md%used = .false.

    ! not changed yet
    md%changed = .false.

    ! store description:
    md%name = name
    md%unit = unit

    ! store time info:
    md%tinterp = tinterp

    ! store data shape:
    md%is   = is
    md%js   = js
    md%halo = halo
    md%ls   = ls

    ! no data allocated yet:
    nullify( md%data  )

    ! no primary data allocated yet:
    nullify( md%data1  )
    md%filled1 = .false.

    ! no secondary data allocated yet:
    nullify( md%data2 )
    md%filled2 = .false.

    ! store input info:
    md%sourcekey = sourcekey

    ! store output info:
    md%putout  = putout
    md%destkey = destkey

    ! ok
    status = 0

  end subroutine mdat_Init


  ! ***


  subroutine mdat_Done( md, status )

    ! --- in/out -----------------------------------

    type(TMeteoData), intent(inout)   ::  md
    integer, intent(out)              ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/mdat_Done'

    ! --- begin --------------------------------

    ! deallocate target data if neccesary:
    if ( associated(md%data) ) then
      ! target data points to data1 ?
      if ( associated(md%data,md%data1) ) then
        ! data points to data1 :
        nullify( md%data )
      else
        ! data is allocated; clear:
        deallocate( md%data  )
      end if
    end if
    ! deallocate primary and seondary data:
    if ( associated(md%data1) ) deallocate( md%data1 )
    if ( associated(md%data2) ) deallocate( md%data2 )

    ! for safety ...
    md%used = .false.
    md%name = 'none'
    md%unit = 'none'

    ! ok
    status = 0

  end subroutine mdat_Done


  ! ***


  subroutine mdat_Set( md, status, used )

    ! --- in/out -----------------------------------

    type(TMeteoData), intent(inout)   ::  md
    integer, intent(out)              ::  status

    logical, intent(in), optional     ::  used

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/mdat_Set'

    ! --- begin --------------------------------

    if ( present(used) ) md%used = used

    ! ok
    status = 0

  end subroutine mdat_Set


  ! ***


  subroutine mdat_Check( md, status )

    ! --- in/out -----------------------------------

    type(TMeteoData), intent(inout)   ::  md
    integer, intent(out)              ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/mdat_Check'

    ! --- begin --------------------------------

    if ( .not. md%used ) then
      write (gol,'("meteo data `",a,"` not in use ...")') trim(md%name); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine mdat_Check


  ! ***


  subroutine mdat_Alloc( md, status )

    ! --- in/out -----------------------------------

    type(TMeteoData), intent(inout)   ::  md
    integer, intent(out)              ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/mdat_Alloc'

    ! --- begin --------------------------------

    ! allocate if field is in use:
    if ( md%used ) then

      ! allocate target, primary, and/or secondary data array;
      ! set all to zero to avoid fpe in halo cells during array operations:
      select case ( md%tinterp )

        ! computed data is stored in target data array:
        case ( 'computed' )

          ! allocate target data:
          allocate( md%data( md%is(1)-md%halo:md%is(2)+md%halo, &
                             md%js(1)-md%halo:md%js(2)+md%halo, &
                             md%ls(1)        :md%ls(2)         ) )
          md%data = 0.0

        ! constant data is storred in primary data;
        ! target data points to data1 :
        case ( 'const', 'month', 'const6', 'const3', 'cpl6', 'cpl3', 'cpl2', 'cpl1' )

          ! allocate primary data:
          allocate( md%data1( md%is(1)-md%halo:md%is(2)+md%halo, &
                              md%js(1)-md%halo:md%js(2)+md%halo, &
                              md%ls(1)        :md%ls(2)         ) )
          md%data1 = 0.0

          ! point target data:
          md%data => md%data1


       case ( 'interp6_3', 'interp6', 'interp3', 'interp2', 'interp1', &
                     'aver1', 'aver3', 'aver6', 'aver24', 'aver24_3'  )

          ! allocate target data:
          allocate( md%data( md%is(1)-md%halo:md%is(2)+md%halo, &
                             md%js(1)-md%halo:md%js(2)+md%halo, &
                             md%ls(1)        :md%ls(2)         ) )
          md%data = 0.0

          ! allocate primary data:
          allocate( md%data1( md%is(1)-md%halo:md%is(2)+md%halo, &
                              md%js(1)-md%halo:md%js(2)+md%halo, &
                              md%ls(1)        :md%ls(2)         ) )
          md%data1 = 0.0

          ! allocate secondary data:
          allocate( md%data2( md%is(1)-md%halo:md%is(2)+md%halo, &
                              md%js(1)-md%halo:md%js(2)+md%halo, &
                              md%ls(1)        :md%ls(2)         ) )
          md%data2 = 0.0

        case default

          write (gol,'("unsupported time interpolation:")'); call goErr
          write (gol,'("  md%tinterp : ",a)') trim(md%tinterp); call goErr
          write (gol,'("  md%name    : ",a)') trim(md%name); call goErr
          TRACEBACK; status=1; return

      end select

    end if

    ! ok
    status = 0

  end subroutine mdat_Alloc


  ! ***

  ! Fill data in md from md2.
  ! If optional time range 'tr' is provided,
  ! the data in mdin might be interpolated to the requested interval.

  subroutine mdat_SetData( md, mdin, status )

    use GO, only : TDate

    ! --- in/out -----------------------------------

    type(TMeteoData), intent(inout)     ::  md
    type(TMeteoData), intent(in)        ::  mdin
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/mdat_SetData'

    ! --- begin --------------------------------

    ! skip rest ?
    if ( .not. md%used ) then
      write (gol,'("WARNING - destination meteo data not in use ...")'); call goPr
      write (*,'("WARNING in ",a)') rname; status=1; return
      status=0; return
    end if

    ! check shapes
    if ( any(md%is  /=mdin%is  ) .or. &
         any(md%js  /=mdin%js  ) .or. &
             md%halo/=mdin%halo  .or. &
         any(md%ls  /=mdin%ls  )  ) then
      write (gol,'("destination and source shapes should be the same:")'); call goErr
      write (gol,'("  is   :  ",2i4," , ",2i4)') md%is  , mdin%is; call goErr
      write (gol,'("  js   :  ",2i4," , ",2i4)') md%js  , mdin%js; call goErr
      write (gol,'("  halo :  ", i8," , ", i8)') md%halo, mdin%halo; call goErr
      write (gol,'("  ls   :  ",2i4," , ",2i4)') md%ls  , mdin%ls; call goErr
      TRACEBACK; status=1; return
    end if

    ! check source data:
    if ( .not. associated(mdin%data) ) then
      write (gol,'("source data not allocated ...")'); call goErr
      TRACEBACK; status=1; return
    end if
    !if ( .not. mdin%filled ) then
    !  write (gol,'("source data not filled")'); call goErr
    !  TRACEBACK; status=1; return
    !end if

    ! check target data:
    if ( md%tinterp /= 'computed' ) then
      write (gol,'("destination data has wrong tinterp:")'); call goErr
      write (gol,'("  expected : ",a)') 'computed'; call goErr
      write (gol,'("  found    : ",a)') trim(md%tinterp); call goErr
      TRACEBACK; status=1; return
    end if
    if ( .not. associated(md%data) ) then
      write (gol,'("destination data not allocated ...")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! check shapes
    if ( any( shape(md%data) /= shape(mdin%data) ) ) then
      write (gol,'("shapes are not the same:")'); call goErr
      write (gol,'("  md   : ",3i5)') shape(  md%data); call goErr
      write (gol,'("  mdin : ",3i5)') shape(mdin%data); call goErr
      TRACEBACK; status=1; return
    end if

    ! copy data:
    md%data = mdin%data
    md%tr   = mdin%tr

    ! ok
    status = 0

  end subroutine mdat_SetData


  ! ***


  subroutine mdat_TimeInterpolation( md, tr, status )

    use go , only : TDate, NewDate, IncrDate, Get
    use go , only : wrtgol, InterpolFractions, rTotal
    use go , only : operator(/=), operator(<), operator(<=)
    use go , only : operator(+), operator(-), operator(/)
    use tmm, only : SetHistory, AddHistory

    ! --- in/out ----------------------------------

    type(TMeteoData), intent(inout)       ::  md
    type(TDate), intent(in)               ::  tr(2)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/mdat_TimeInterpolation'

    ! --- local ----------------------------------

    integer           ::  direction
    integer           ::  dth, baseh
    integer           ::  year, month, day, hour, minu
    type(TDate)       ::  tmid, tc(2)
    real              ::  alfa1, alfa2

    ! --- begin -----------------------------

    ! not used ? error
    if ( .not. md%used ) then
      write (gol,'("meteo data `",a,"` not used")') trim(md%name); call goErr
      TRACEBACK; status=1; return
    end if

    ! time direction:
    if ( tr(1) <= tr(2) ) then
      direction =  1
    else
      direction = -1
    end if

    ! different actions based on time interpolation type:
    select case ( md%tinterp )

      !
      ! constant data
      !

      case ( 'const' )

        ! md%data points to md%data1, so nothing to be done

      !
      ! data constant during interval
      !

      case ( 'month' )

        ! check time: tr should be in md%tr1
        if ( ((direction > 0) .and. ((tr(2) < md%tr1(1)) .or. (md%tr1(2) < tr(1)))) .or. &
             ((direction < 0) .and. ((tr(1) < md%tr1(2)) .or. (md%tr1(1) < tr(2))))      ) then
          write (gol,'("model data does not include requested interval:")'); call goErr
          write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
          call wrtgol( '  md%tr1      : ', md%tr1(1), ' - ', md%tr1(2) ); call goErr
          call wrtgol( '  tr          : ', tr(1), ' - ', tr(2) ); call goErr
          TRACEBACK; status=1; return
        end if

        ! md%data points to md%data1, so no changes

      case ( 'const6', 'const3' )

        select case ( md%tinterp )
          case ( 'const3' ) ; baseh = 0 ; dth = 3
          case ( 'const6' ) ; baseh = 0 ; dth = 6
        end select

        ! extract time values for begin of current interval:
        call Get( tr(1), year, month, day, hour, minu )
        ! round hour to 00/06/12/18 or 00/03/06/09/12/15/18/21 or 09

        !WP! changed line to NOT use the nint function, instead just use int+0.5 which is the same
        !WP! for positive numbers, see http://h21007.www2.hp.com/portal/download/files/unprot/Fortran/docs/lrm/lrm0299.htm
        !WP! For negative numbers, the lines below need to be added, this can happen only if baseh!=0
        !WP! This fixes a problem on the NOAA machines that cause unpredictable crashes when using nint()
        !WP!
        !dummy=real(hour+minu/60.0-baseh)/real(dth)
        !if(dummy>0.0) hour=dth*int(dummy+0.5)+baseh
        !if(dummy<0.0) hour=dth*int(dummy-0.5)+baseh
        !endif
        !WP!
        !WP!
        ! ported from cy2, 24 Jan 2008, ARJ
        hour = dth * int(real(hour+minu/60.0-baseh)/real(dth)+0.5) + baseh


        ! set mid of 3 or 6 hour interval:
        tmid = NewDate( year, month, day, hour )
        ! interval with constant field
        tc(1) = tmid - IncrDate(hour=dth)/2
        tc(2) = tmid + IncrDate(hour=dth)/2

        ! check interval
        if ( (tr(1) < tc(1)) .or. (tc(2) < tr(2)) ) then
          write (gol,'("time intervals do not match:")'); call goErr
          call wrtgol( '  requested     : ', tr(1), ' - ', tr(2) ); call goErr
          call wrtgol( '  mdin valid for : ', tc(1), ' - ', tc(2) ) ; call goErr
          write (gol,'(" mdin%tinterp : ",a)') trim(md%tinterp); call goErr
          TRACEBACK; status=1; return
        end if

        ! md%data points to md%data1, so no changes

      !
      ! coupling: get field for t, use it for [t,t+dt]
      !
      case ( 'cpl6', 'cpl3', 'cpl2', 'cpl1' )

        select case ( md%tinterp )
          case ( 'cpl6' ) ; baseh = 0 ; dth = 6
          case ( 'cpl3' ) ; baseh = 0 ; dth = 3
          case ( 'cpl2' ) ; baseh = 0 ; dth = 2
          case ( 'cpl1' ) ; baseh = 0 ; dth = 1
        end select

        ! extract time values for begin of current interval:
        call Get( tr(1), year, month, day, hour, minu )
        ! round hour to 00/06/12/18 or 00/03/06/09/12/15/18/21
        hour = dth * floor(real(hour+minu/60.0-baseh)/real(dth)) + baseh
        ! interval with constant field
        tc(1) = NewDate( year, month, day, hour )
        tc(2) = tc(1) + IncrDate(hour=dth)

        ! check interval
        if ( (tr(1) < tc(1)) .or. (tc(2) < tr(2)) ) then
          write (gol,'("time intervals do not match:")'); call goErr
          call wrtgol( '  requested     : ', tr(1), ' - ', tr(2) ); call goErr
          call wrtgol( '  mdin valid for : ', tc(1), ' - ', tc(2) ) ; call goErr
          write (gol,'(" mdin%tinterp : ",a)') trim(md%tinterp); call goErr
          TRACEBACK; status=1; return
        end if

        ! md%data points to md%data1, so no changes

      !
      ! linear interpolation between instant times
      !
      case ( 'interp6', 'interp6_3', 'interp3', 'interp2', 'interp1' )

        ! not filled ? error
        if ( (.not. md%filled1) .or. (.not. md%filled2) ) then
          write (gol,'("meteo data not filled:")'); call goErr
          write (gol,'("  name        : ",a)') trim(md%name); call goErr
          write (gol,'("  filled      : ",2l2)') md%filled1, md%filled2; call goErr
          write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
          TRACEBACK; status=1; return
        end if

        ! interpolation between instant times, not between intervals ...
        if ( (md%tr1(1) /= md%tr1(2)) .or. (md%tr2(1) /= md%tr2(2)) ) then
          write (gol,'("time interpolation not between intervals:")'); call goErr
          write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
          call wrtgol( '  tr1         : ', md%tr1(1), ' - ', md%tr1(2) ); call goErr
          call wrtgol( '  tr2         : ', md%tr2(1), ' - ', md%tr2(2) ); call goErr
          TRACEBACK; status=1; return
        end if

        ! interpolate to mid of interval:
        tmid = tr(1) + (tr(2)-tr(1))/2

        ! deterimine weights to data and data2 :
        call InterpolFractions( tmid, md%tr1(1), md%tr2(1), alfa1, alfa2, status )
        if (status/=0) then; TRACEBACK; return; end if

        ! Since ifort 15.0, simple array additions such as the one below *should* be parallelized by workshare
        ! https://software.intel.com/en-us/articles/openmp-workshare-constructs-now-parallelize-with-intel-fortran-compiler-150
        ! However, this results in a segfault with ifort 15.0.3. So I'll get rid of the openmp directives for now.
        !!$omp parallel
        !!$omp workshare
        md%data = alfa1 * md%data1 + alfa2 * md%data2
        !!$omp end workshare
        !!$omp end parallel

        ! data is changed ...
        md%changed = .true.

      !
      ! fractions of time average fields:
      !   data1  :   [tr1(1),tr1(2)]
      !   data2  :                  [tr1(1),tr1(2)]
      !   tr     :               [tr(1),tr(2)]
      !

      case ( 'aver1', 'aver3', 'aver6', 'aver24', 'aver24_3' )

        ! primary data not filled ? error
        if ( .not. md%filled1 ) then
          write (gol,'("meteo data1 not filled:")'); call goErr
          write (gol,'("  name        : ",a)') trim(md%name); call goErr
          write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
          TRACEBACK; status=1; return
        end if

        ! tr earlier than tr1 ? error ...
        if ( tr(1) < md%tr1(1) ) then
          write (gol,'("requested time interval earlier than data:")'); call goErr
          write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
          call wrtgol( '  md%tr1      : ', md%tr1(1), ' - ', md%tr1(2) ); call goErr
          call wrtgol( '  tr          : ', tr(1), ' - ', tr(2) ); call goErr
          TRACEBACK; status=1; return
        end if

        ! tr complete in tr1 ? simple ...
        if ( tr(2) <= md%tr1(2) ) then

          ! just copy ...
          md%data = md%data1

          ! data is changed ...
          md%changed = .true.

        else

          ! fractions of data1 and data2

          ! secondary data not filled ? error
          if ( .not. md%filled2 ) then
            write (gol,'("meteo data2 not filled:")'); call goErr
            write (gol,'("  name        : ",a)') trim(md%name); call goErr
            write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
            TRACEBACK; status=1; return
          end if

          ! time ranges for data1 and data2 should be connected:
          if ( md%tr1(2) /= md%tr2(1) ) then
            write (gol,'("time intervals not connected:")'); call goErr
            call wrtgol( '  md%tr1      : ', md%tr1(1), ' - ', md%tr1(2) ); call goErr
            call wrtgol( '  md%tr2      : ', md%tr2(1), ' - ', md%tr2(2) ); call goErr
            write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
            TRACEBACK; status=1; return
          end if

          ! check requested time range:
          if ( (tr(1) < md%tr1(1)) .or. (md%tr2(2) < tr(2)) ) then
            write (gol,'("requested time interval not covered by data :")'); call goErr
            call wrtgol( '  md%tr1      : ', md%tr1(1), ' - ', md%tr1(2) ); call goErr
            call wrtgol( '  md%tr2      : ', md%tr2(1), ' - ', md%tr2(2) ); call goErr
            call wrtgol( '  tr          : ', tr(1), ' - ', tr(2) ); call goErr
            write (gol,'("  md%tinterp  : ",a)') trim(md%tinterp); call goErr
            TRACEBACK; status=1; return
          end if

          ! first fraction:
          if ( tr(1) < md%tr1(2) ) then
            if ( tr(2) < md%tr1(2) ) then
              ! tr complete in tr1 ...
              alfa1 = 1.0
            else
              ! fraction of tr inside tr1 :
              alfa1 = rTotal(md%tr1(2)-tr(1),'sec') / rTotal(tr(2)-tr(1),'sec')
            end if
          else
            ! tr later than tr1, thus not covered:
            alfa1 = 0.0
          end if

          ! second fraction:
          if ( md%tr2(1) < tr(2) ) then
            if ( md%tr2(1) < tr(1) ) then
              ! tr complete in tr2 ...
              alfa2 = 1.0
            else
              ! fraction of tr inside tr2 :
              alfa2 = rTotal(tr(2)-md%tr2(1),'sec') / rTotal(tr(2)-tr(1),'sec')
            end if
          else
            ! tr before tr2, thus not covered:
            alfa2 = 0.0
          end if

          ! replace data array:
          md%data = alfa1 * md%data1 + alfa2 * md%data2

          ! data is changed ...
          md%changed = .true.

        end if

      !
      ! unknown ...
      !

      case default

        write (gol,'("unsupported time interpolation:")'); call goErr
        write (gol,'("  md%tinterp : ",a)') trim(md%tinterp); call goErr
        write (gol,'("  md%name    : ",a)') trim(md%name); call goErr
        TRACEBACK; status=1; return

    end select

    ! store new time:
    md%tr = tr

    ! copy history:
    call SetHistory( md%tmi, md%tmi1, status )

    ! ok
    status = 0

  end subroutine mdat_TimeInterpolation

end module MeteoData
