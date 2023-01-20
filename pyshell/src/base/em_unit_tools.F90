!----------------------------------------------------------
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!----------------------------------------------------------

module em_unit_tools

  ! --- modules ------------------------------

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------

  private

  public :: Convert_Units
  public :: UnitStringToType
  public :: UnitTypeToString

  ! --- const ------------------------------

  character(len=*), parameter  :: mname = 'em_unit_tools'

contains

  subroutine Convert_Units( em, unit, status )

    ! --------------------------------------------
    ! Convert units in 'em' to 'unit'
    ! Exceptions:
    !  1. Units yr-1 are changed to se-1 or mo-1 for seasonal
    !     or monthly fields, respectively
    !  2. If species full molecular mass is undefined then
    !     output is in element mass (e.g. for NOx)
    ! --------------------------------------------

    ! --- modules --------------------------------

    use tipp_base, only : TEmis
    use tipp_base, only : TEmUnit

    ! --- in/out ---------------------------------
    
    type(TEmis), intent(inout)    :: em
    character(len=*), intent(in)  :: unit
    integer, intent(out)          :: status

    ! --- const --------------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Convert_Units'

    ! --- local --------------------------------------

    integer                       :: ix, iy, iz, ip
    integer                       :: ival
    type(TEmUnit)                 :: unit_in_tp, unit_out_tp
    real, dimension(em%ny,em%np)  :: cf_src, cf_tgt

    ! --- begin ---------------------------------

    ! Put source units in type
    call UnitStringToType( em%unit, unit_in_tp, status )
    if ( status /= 0 ) then
       write (gol,'("failed to put source units in type")'); call goErr
       TRACEBACK; status=1; return
    end if

    ! Put target units in type
    call UnitStringToType( unit, unit_out_tp, status )
    if ( status /= 0 ) then
       write (gol,'("failed to put target units in type")'); call goErr
       TRACEBACK; status=1; return
    end if

    ! If full molecular mass undefined (as for NOx) convert to element mass
    if ( em%xm_fm == -1. ) unit_out_tp%amount_type = 'el'

    ! Check whether timeframe must be changed
    ! NOTE: this is set such that sum over all periods gives you the yearly total ...
    if ( unit_out_tp%timeframe == 'yr-1' ) then
      select case ( em%np )
        case ( 4          ) ; unit_out_tp%timeframe = 'se-1'
        case ( 12         ) ; unit_out_tp%timeframe = 'mo-1'
        case ( 46         ) ; unit_out_tp%timeframe = '8d-1'
        case ( 365, 366   ) ; unit_out_tp%timeframe = 'da-1'  ! daily
        case ( 8760, 8784 ) ; unit_out_tp%timeframe = 'h-1'   ! hourly
      end select
    end if

    ! Get conversion factor from source to standard
    call Convert_Emission( unit_in_tp, em%xm_fm, em%xm_el, &
         em%nx, em%ny, em%np, em%year, cf_src, status )
    if ( status /= 0 ) then
       write (gol,'("failed conversion source units to standard units")'); call goErr
       TRACEBACK; status=1; return
    end if

    ! Get conversion factor from target to standard
    call Convert_Emission( unit_out_tp, em%xm_fm, em%xm_el, &
         em%nx, em%ny, em%np, em%year, cf_tgt, status )
    if ( status /= 0 ) then
       write (gol,'("failed conversion target units to standard units")'); call goErr
       TRACEBACK; status=1; return
    end if

    ! Apply conversion to field
    select case ( trim(em%storage) )
      !~~ full 4D field
      case ( 'full' )
        ! loop over grid cells:
        do ix = 1, em%nx
           do iz = 1, em%nz
              em%field(ix,:,iz,:) = em%field(ix,:,iz,:) * cf_src / cf_tgt
           end do
        end do
      !~~ sparse field
      case ( 'sparse' )
        ! loop over non-zero values:
        do ival = 1, em%sparse_n
          ! current y and p indices:
          iy = em%sparse_index(2,ival)
          ip = em%sparse_index(4,ival)
          ! convert value:
          em%sparse_field(ival) = em%sparse_field(ival) * cf_src(iy,ip) / cf_tgt(iy,ip)
        end do
      !~~
      case default
        write (gol,'("unsupported storage type : ",a)') trim(em%storage); call goErr
        TRACEBACK; status=1; return
    end select

    ! Rebuild target-unit string and put in em
    call UnitTypeToString( unit_out_tp, em%unit, status )
    if ( status /= 0 ) then
       write (gol,'("failed to put target unit type in string")'); call goErr
       TRACEBACK; status=1; return
    end if

    status = 0

  end subroutine Convert_Units
  
  
  ! ***
  
  

  subroutine UnitStringToType( unit_str, unit_tp, status )
    !----------------------------------------------------------
    ! Convert unit string to an emission unit type
    ! Emission string should (in random order) contain:
    !  1. amount of species, specified as
    !        mass plus specification of whether this mass refers
    !        to full molecule or main element (e.g. C in CH4)
    !                  [g fmm, g el, kg fmm, kg el, etc.]
    !     or number of molecules   [mlcs]
    !  2. spatial domain           [m-2, km-2, cl-1, ..]
    !  3. temporal domain          [s-1, da-1, mo-1, yr-1]
    ! These components should be separated by one or more spaces.
    !
    ! INPUT:
    !  o unit_str     - units as string
    !
    ! OUTPUT:
    !  o unit         - units as type
    !  o status
    !----------------------------------------------------------

    ! --- modules ---------------------------------------------

    use tipp_base, only : amount_units, amount_type_units, area_units, timeframe_units
    use tipp_base, only : TEmUnit

    ! --- in/out ----------------------------------------------

    character(len=*), intent(in)  :: unit_str
    type(TEmUnit), intent(out)    :: unit_tp
    integer, intent(out)          :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter   :: rname = mname//'/UnitStringToType'

    ! --- local -----------------------------------------------

    character(len=len(unit_str))  :: un
    character(len=8)              :: tag
    integer                       :: i1, i2
    logical                       :: amount_present = .false.
    logical                       :: amount_type_present = .false.
    logical                       :: area_present = .false.
    logical                       :: timeframe_present = .false.

    ! --- begin -----------------------------------------------

    un = unit_str
    do
       if ( len_trim(un) == 0 ) exit
       i1 = scan( un, ' ' )
       select case ( i1 )
       case ( 0 )
          exit
       case ( 1 )
          un = un(2:)
          cycle
       case default
          tag = un( 1:i1-1 )
          if ( any( amount_units == tag ) ) then
             unit_tp%amount = tag
             amount_present = .true.
             if ( tag == 'mlcs' .or. tag == 'mol' .or. tag == 'mmol' ) then
                unit_tp%amount_type = 'el'
                amount_type_present = .true.
             end if
          else if ( any( amount_type_units == tag ) ) then
             unit_tp%amount_type = tag
             amount_type_present = .true.
          else if ( any( area_units == tag ) ) then
             unit_tp%area = tag
             area_present = .true.
          else if ( any( timeframe_units == tag ) ) then
             unit_tp%timeframe = tag
             timeframe_present = .true.
          else
             write (gol,'("unknown substring `",a,"` in unit `",a,"`")') &
                  trim(tag), trim(unit_str); call goErr
             TRACEBACK; status=1; return
          end if
          un = un(i1+1:)
       end select
    end do

    if ( .not. amount_present ) then
       write (gol,'("no amount unit found in `",a,"`")') unit_str; call goErr
       TRACEBACK; status=1; return
    end if
    if ( .not. amount_type_present ) then
       write (gol,'("no amount type unit found in `",a,"`")') unit_str; call goErr
       TRACEBACK; status=1; return
    end if
    if ( .not. area_present ) then
       write (gol,'("no area unit found in `",a,"`")') unit_str; call goErr
       TRACEBACK; status=1; return
    end if
    if ( .not. timeframe_present ) then
       write (gol,'("no time unit found in `",a,"`")') unit_str; call goErr
       TRACEBACK; status=1; return
    end if

    status = 0

  end subroutine UnitStringToType


  !----------------------------------------------------------
  ! Convert unit type to string
  !
  ! INPUT:
  !  o unit         - units as type
  !
  ! OUTPUT:
  !  o unit_str     - units as string
  !  o status
  !----------------------------------------------------------

  subroutine UnitTypeToString( unit_tp, unit_str, status )

    use tipp_base, only : TEmUnit

    ! --- in/out ----------------------------------------------

    type(TEmUnit), intent(in)      :: unit_tp
    character(len=50), intent(out) :: unit_str
    integer, intent(out)           :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter   :: rname = mname//'/UnitTypeToString'

    ! --- begin -----------------------------------------------

    unit_str = trim(unit_tp%amount)//' '//trim(unit_tp%amount_type)//' '// &
         trim(unit_tp%area)//' '//trim(unit_tp%timeframe)

    status = 0

  end subroutine UnitTypeToString


  subroutine Convert_Emission( un, xm_fm, xm_el, nx, ny, np, year, &
       ConversionFactor, status )
    !----------------------------------------------------------
    ! Calculate the multiplicitave factor for conversion from
    ! specified units to 'standard' units
    ! 
    ! Standard units for emissions: g el cl-1 yr-1
    !
    ! INPUT:
    !  o un           - source units
    !  o xm_fm, xm_el - molecular mass of full molecule / element
    !  o nx, ny, np   - nr of longitudes, latitudes, periods
    !  o year         - year of data [0 is arbitrary non-leap year]
    !
    ! OUTPUT:
    !  o ConversionFactor   - array of size `(ny,np)`
    !  o status
    !
    !----------------------------------------------------------

    use tipp_base, only : TEmUnit

    ! --- in/out ----------------------------------------------

    type(TEmUnit), intent(in)     :: un
    real, intent(in)              :: xm_fm, xm_el
    integer, intent(in)           :: nx, ny, np, year
    real, dimension(ny,np), intent(out) :: ConversionFactor
    integer, intent(out)          :: status
    
    ! --- const -----------------------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Convert_Emission'
    
    ! --- local -----------------------------------------------

    integer                       :: i, j
    real                          :: cf_am
    real, dimension(ny)           :: cf_ar
    real, dimension(np)           :: cf_tf

    ! --- begin -----------------------------------------------

    ! Get respective conversion factors
    call Convert_Amount( un, xm_fm, xm_el, cf_am, status )
    if ( status /= 0 ) then
       TRACEBACK; status=1; return
    end if

    call Convert_Area( un, nx, ny, cf_ar, status )
    if ( status /= 0 ) then
       TRACEBACK; status=1; return
    end if

    call Convert_Timeframe( un, year, np, cf_tf, status )
    if ( status /= 0 ) then
       TRACEBACK; status=1; return
    end if

    ! Combine into total conversion factor
    do i = 1, np
       do j = 1, ny
          ConversionFactor(j,i) = cf_am*cf_ar(j)*cf_tf(i)
       end do
    end do

    status = 0

  end subroutine Convert_Emission


  subroutine Convert_Amount( un, xm_fm, xm_el, ConversionFactor, status )
    !----------------------------------------------------------
    ! Calculate the multiplicitave factor for conversion from
    ! specfified units to 'standard' units
    ! 
    ! Standard amount units for emissions: g el
    !
    ! INPUT:
    !  o un           - source units
    !  o xm_fm, xm_el - molecular mass of full molecule / element
    !
    ! OUTPUT:
    !  o ConversionFactor
    !  o status
    !
    !----------------------------------------------------------
    
    ! --- modules ---------------------------------------------

    use binas
    use tipp_base, only            : namount, amount_units, &
         namount_type, amount_type_units
    use tipp_base, only : TEmUnit

    ! --- in/out ----------------------------------------------

    type(TEmUnit), intent(in)     :: un
    real, intent(in)              :: xm_fm, xm_el
    real, intent(out)             :: ConversionFactor
    integer, intent(out)          :: status
    
    ! --- const -----------------------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Convert_Amount'
    
    ! --- local -----------------------------------------------

    integer                       :: i, iam, iamt
    real, dimension(namount)      :: am_cf
    real, dimension(namount_type) :: amt_cf

    ! --- begin -----------------------------------------------

    ! Amount conversion to g
    am_cf = (/ 1.e-6, 1.e-3, 1., 1.e3, 1.e6, 1.e12, &
         xm_el/Avog, xm_el, xm_el*1.e-3 /)
          
    ! Amount type conversion to el
    amt_cf = (/ xm_el/xm_fm, 1. /)

    ! Check amount units
    iam = 0
    do i = 1, namount
       if ( un%amount == amount_units(i) ) iam = i
    end do
    if ( iam == 0 ) then
       write (gol,'("unknown amount units `",a,"`")') un%amount; call goErr
       TRACEBACK; status=1; return
    end if

    ! Check amount type units
    iamt = 0
    do i = 1, namount_type
       if ( un%amount_type == amount_type_units(i) ) iamt = i
    end do
    if ( iamt == 0 ) then
       write (gol,'("unknown amount type units `",a,"`")') un%amount_type; call goErr
       TRACEBACK; status=1; return
    end if

    ! Conversion factor
    ConversionFactor = am_cf(iam)*amt_cf(iamt)

    if ( ConversionFactor <= 0. ) then
       write (gol,'("cannot perform mass unit conversion")'); call goErr
       write (gol,'("xm_fm       : ",f5.1)') xm_fm; call goErr
       write (gol,'("xm_el       : ",f5.1)') xm_el; call goErr
       write (gol,'("source units: ",a," ",a)') trim( un%amount ), &
            trim( un%amount_type ); call goErr
       write (gol,'("target units: g el")'); call goErr
       TRACEBACK; status=1; return
    end if

    status = 0

  end subroutine Convert_Amount


  subroutine Convert_Area( un, nx, ny, ConversionFactor, status )
    !----------------------------------------------------------
    ! Calculate the multiplicitave factor for conversion from
    ! specfified units to 'standard' units
    ! 
    ! Standard units for emissions: cl-1 (standard grid = 1x1 deg)
    ! A regular lon-lat grid is assumed.
    !
    ! INPUT:
    !  o un          - source units
    !  o nx, ny      - nr of longitudes, latitudes
    !
    ! OUTPUT:
    !  o ConversionFactor   - array of size `ny`
    !  o status
    !
    !----------------------------------------------------------
    
    ! --- modules ---------------------------------------------

    use binas
    use Grid_Tools,                only : ll_area
    use tipp_base,                 only : narea, area_units
    use tipp_base, only : TEmUnit

    ! --- in/out ----------------------------------------------

    type(TEmUnit), intent(in)        :: un
    integer, intent(in)              :: nx, ny
    real, dimension(ny), intent(out) :: ConversionFactor
    integer, intent(out)             :: status
    
    ! --- const -----------------------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Convert_Area'
    
    ! --- local -----------------------------------------------

    integer                       :: i, j, iar
    real, dimension(narea,ny)     :: area_cf
    real, dimension(narea-1)      :: to_m2

    ! --- begin -----------------------------------------------

    ! Area conversion to cl-1
    to_m2 = (/ 1.e4, 1., 1.e-6 /)
    do j = 1, ny
       area_cf(1:narea-1,j) = to_m2 * ae**2 * &
            ll_area( 0., 360./nx*deg2rad, &
            (-90. + (j-1)*180./ny)*deg2rad, (-90. + j*180./ny)*deg2rad)
    end do
    area_cf(narea,:) = 1.

    ! Check area units
    iar = 0
    do i = 1, narea
       if ( un%area == area_units(i) ) iar = i
    end do
    if ( iar == 0 ) then
       write (gol,'("unknown area units `",a,"`")') un%area; call goErr
       TRACEBACK; status=1; return
    end if

    ! Conversion factor
    ConversionFactor = area_cf(iar,:)

    status = 0

  end subroutine Convert_Area


  subroutine Convert_Timeframe( un, year, np, ConversionFactor, status )
    !----------------------------------------------------------
    ! Calculate the multiplicitave factor for conversion from
    ! specfified units to 'standard' units
    ! 
    ! Standard timeframe units for emissions: yr-1
    !
    ! INPUT:
    !  o un          - source units
    !  o year        - year of data [0 is default non-leap year]
    !  o np          - number of periods
    !
    ! OUTPUT:
    !  o ConversionFactor   - array of size `np`
    !  o status
    !
    !----------------------------------------------------------
    
    ! --- modules ---------------------------------------------

    use go_date
    use tipp_base, only            : ntimeframe, timeframe_units
    use tipp_base, only : TEmUnit

    ! --- in/out ----------------------------------------------

    type(TEmUnit), intent(in)     :: un
    integer, intent(in)           :: year
    integer, intent(in)           :: np
    real, dimension(np), intent(out)  :: ConversionFactor
    integer, intent(out)          :: status
    
    ! --- const -----------------------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Convert_Timeframe'
    
    ! --- local -----------------------------------------------

    integer                        :: i, j, imo, ise, itf, ip
    integer                        :: ye
    real                           :: diy, dis, dim, dip
    type(TDate)                    :: t
    real, dimension(ntimeframe,np) :: tf_cf

    ! --- begin -----------------------------------------------

    ! Timeframe conversion to yr-1
    
    ! set year to be used:
    ye = year
    ! if year is '0000' this actually means 'constant' ;
    ! a unit 'kg/year' should then be interpreted as 'kg/(year of 365 days)'
    if ( year == 0000 ) ye = 1999   ! arbitrary non-leap year of 365 days
    
    ! januari 1 of (current) or dummy year:
    t = NewDate( year=ye, month=1, day=1 )
    
    ! days in this year:
    diy = real( Days_in_Year( t ) )
    
    !
    ! conversion factors for each (timeframe,period)
    !
    ! ~ undefined by default:
    tf_cf = -1
    !
    ! ~ standard conversions:
    !   (names and meaning of time frames defined in tipp_base)
    tf_cf(1,:) = 60*60*24*diy        ! seconds to year
    tf_cf(2,:) = 60*24*diy           ! minutes to year
    tf_cf(3,:) = 24*diy              ! hours   to year
    tf_cf(4,:) = diy                 ! days    to year
    tf_cf(5,:) = -1                  ! month   to year ; set below
    tf_cf(6,:) = -1                  ! season  to year ; set below
    tf_cf(7,:) = 1.0                 ! year    to year
    tf_cf(8,:) = -1                  ! 8day    to year ; set below
    !
    ! ~ specials:
    select case ( np )
      ! ~~ yearly data
      case( 1 )
        ! no extra conversion factors needed
      ! ~~ seasonal data (4 times 3 months)
      case( 4 )                     ! seasonal data
        ! loop over seasons:
        do ise = 1, 4
          ! days in season:
          dis = 0.
          do j = 1, 3
             imo = modulo( (ise - 1)*3 - 2 + j, 12 ) + 1
             !t = NewDate( year=ye, month=imo, day=1 )
             t = NewDate( year=year, month=imo, day=1 )
             dis = dis + real( Days_in_Month( t ) )
          end do
          ! conversion to yearly number:
          tf_cf(6,ise) = diy/dis       ! season to year
        end do
      ! ~~ monthly data
      case( 12 )
        ! loop over months:
        do imo = 1, 12
          ! days in month:
          !t = NewDate( year=ye, month=imo, day=1 )
          t = NewDate( year=year, month=imo, day=1 )
          dim = real( Days_in_Month( t ) )
          ! conversion to yearly number:
          tf_cf(5,imo) = diy/dim     ! month to year
        end do
      ! ~~ 8 daily data
      case ( 46 )
        ! periods valid for 8 days:
        tf_cf(8,:) = diy/8.0   ! 8-days to year
        ! ... except the last one which has less days:
        tf_cf(8,np) = diy/(diy-(np-1)*8.0)
      ! ~~ daily data
      case ( 365, 366 )
        ! no extra conversion factors needed
      ! ~~ hourly data
      case ( 8760, 8784 )
        ! no extra conversion factors needed
      ! ~~ error
      case default
        write (gol,'("unknown number of periods ",i5)') np; call goErr
        TRACEBACK; status=1; return
    end select

    ! Check timeframe units
    itf = 0
    do i = 1, ntimeframe
       if ( un%timeframe == timeframe_units(i) ) itf = i
    end do
    if ( itf == 0 ) then
       write (gol,'("unknown timeframe units `",a,"`")') un%timeframe; call goErr
       TRACEBACK; status=1; return
    end if

    ! Conversion factor
    ConversionFactor = tf_cf(itf,:)

    if ( any( ConversionFactor < 0. ) ) then
       write (gol,'("cannot convert to yr-1")'); call goErr
       write (gol,'("timeframe units `",a,"`")') un%timeframe; call goErr
       write (gol,'("number of periods of data ",i5)') np; call goErr
       TRACEBACK; status=1; return
    end if

    status = 0

  end subroutine Convert_Timeframe

end module em_unit_tools
