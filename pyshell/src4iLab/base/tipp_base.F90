!----------------------------------------------------------
!
!  Definitions, types and parameters related to
!  emissions and boundary conditions
!
!----------------------------------------------------------
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!----------------------------------------------------------

module tipp_base

  ! --- modules ------------------------------

  use GO, only                     : TDate
  use grid, only                   : TLevelInfo, TllGridInfo, Init, Done
  use os_specs, only               : MAX_FILENAME_LEN, DUMMY_STR_LEN

  implicit none


  ! --- in/out -----------------------------

  private

  public   :: c_in_nmv
  public   :: namount, amount_units, namount_type, amount_type_units
  public   :: narea, area_units, ntimeframe, timeframe_units
  public   :: TEmis, TEmUnit, TBouCon
  public   :: Init
  public   :: Emis_Sparse_Init, Emis_Sparse_Done, Emis_Sparse_Add


  ! --- const ------------------------------

  character(len=*), parameter     :: mname = 'module tipp_base'

  ! mole masses
  real, parameter               :: xmch4 = 16.
  real, parameter               :: xmco2 = 44.
  real, parameter               :: xmco = 28.
  real, parameter               :: xm14co = 30.
  real, parameter               :: xmno = 30.
  real, parameter               :: xmno2 = 46.
  real, parameter               :: xmnh3 = 17.
  real, parameter               :: xmso2 = 64.
  real, parameter               :: xmdms = 62.   ! (CH3)2S
  real, parameter               :: xmh2s = 34.
  real, parameter               :: xmnmv = 55.   ! 'Indicative value for NMV mixture'
  real, parameter               :: xmc2h4 = 28.  ! ethene
  real, parameter               :: xmc2h6 = 30.  ! ethane
  real, parameter               :: xmc3h6 = 42.  ! propene
  real, parameter               :: xmc3h8 = 44.  ! propane
  real, parameter               :: xmbute = 56.  ! butene
  real, parameter               :: xmbuta = 58.  ! butane
  real, parameter               :: xmch2o = 30.  ! formaldehyde
  real, parameter               :: xmalde = 44.  ! higher aldehydes (take mass of CH3CHO)
  real, parameter               :: xmch4o = 32.  ! methanol
  real, parameter               :: xmalco = 46   ! higher alcohols (take mass of C2H5OH)
  real, parameter               :: xmacet = 58.  ! aceton = CH3COCH3
  real, parameter               :: xmmek = 72.   ! methyl-ethyl ketone and higher ketones
  real, parameter               :: xmtolu = 92.  ! toluene C7H8 and all other aromatic species
  real, parameter               :: xmisop = 68.  ! isoprene C5H8
  real, parameter               :: xmterp = 136  ! terpenes (take mass of monoterpene (C5H8)2)
  real, parameter               :: xmc = 12.
  real, parameter               :: xm14c = 14.
  real, parameter               :: xmn = 14.
  real, parameter               :: xms = 32.

  ! kg C / kg NMV in typical NMV mixture
  real, parameter               :: c_in_nmv = 0.8

  ! number of species
  integer, parameter            :: nspec = 27
  ! names
  character(len=4), parameter, dimension(nspec)  :: spec_names = &
       (/ 'CH4 ', 'CO2 ', 'CO  ', 'NO  ', 'NO2 ', 'NOx ', &
       'NH3 ', 'SO2 ', 'DMS ', 'H2S ', 'NMV ', 'C2H4', &
       'C2H6', 'C3H6', 'C3H8', 'BUTE', 'BUTA', 'CH2O', &
       'ALDE', 'CH4O', 'ALCO', 'ACET', 'MEK ', 'TOLU', &
       'ISOP', 'TERP', '14CO' /)
  ! full molecular masses
  real, parameter, dimension(nspec) :: spec_xm_fm = &
       (/ xmch4, xmco2, xmco, xmno, xmno2, -1., &
       xmnh3, xmso2, xmdms, xmh2s, xmnmv, xmc2h4, &
       xmc2h6, xmc3h6, xmc3h8, xmbute, xmbuta, xmch2o, &
       xmalde, xmch4o, xmalco, xmacet, xmmek, xmtolu, &
       xmisop, xmterp, xm14co /)
  ! mass of main element
  real, parameter, dimension(nspec) :: spec_xm_el = &
       (/ xmc, xmc, xmc, xmn, xmn, xmn, &
       xmn, xms, xms, xms, xmnmv*c_in_nmv, 2*xmc, &
       2*xmc, 3*xmc, 3*xmc, 4*xmc, 4*xmc, xmc, &
       2*xmc, xmc, 2*xmc, 3*xmc, 4*xmc, 7*xmc, &
       5*xmc, 10*xmc, xm14c /)

  ! Allowed units
  ! NOTE: make sure all elements of a character array have equal length
  integer, parameter           :: namount = 9
  character(len=8), dimension(namount), parameter    :: amount_units = &
       (/ 'mug ', 'mg  ', 'g   ', 'kg  ', 'Mg  ', 'Tg  ', 'mlcs', 'mol ', 'mmol' /)
  integer, parameter           :: namount_type = 2
  character(len=8), dimension(namount_type), parameter :: amount_type_units = &
       (/ 'fm', 'el' /)
  integer, parameter           :: narea = 4
  character(len=8), dimension(narea), parameter      :: area_units = &
       (/ 'cm-2', 'm-2 ', 'km-2', 'cl-1' /)

  ! allowed time frame units:
  integer, parameter           :: ntimeframe = 8
  character(len=8), dimension(ntimeframe), parameter :: timeframe_units = &
       (/ 's-1 ', 'mi-1', 'h-1 ', 'da-1', 'mo-1', 'se-1', 'yr-1', '8d-1' /)

  ! --- type -------------------------------

  ! General structure for emission field incl metadata
  type TEmis
    ! * filename:
    character(len=DUMMY_STR_LEN)          ::  name              ! short name
    character(len=MAX_FILENAME_LEN)         ::  filename          ! full name of HDF file
    ! * file format number:
    !     1 = classic
    !     2 = option to store sparse array
    integer                     ::  format_nr
    ! * how stored ?
    character(len=32)           ::  storage           ! 'full' | 'sparse'
    ! * 3D grid:
    integer                     ::  nx                ! longitudinal dimension
    integer                     ::  ny                ! latitudinal dimension
    integer                     ::  nz                ! vertical dimension
    integer                     ::  np                ! time dimension (nr. of periods)
    character(len=80)           ::  x_def             ! definition of long. dimension
    character(len=80)           ::  y_def             ! definition of lat. dimension
    character(len=80)           ::  z_def             ! definition of vert. dimension
    real, pointer               ::  x_grid(:)         ! long. grid
    real, pointer               ::  y_grid(:)         ! lat. grid
    real, pointer               ::  z_a_grid(:)       ! hybrid grid, a coeff
    real, pointer               ::  z_b_grid(:)       ! hybrid grid, b coeff
    character(len=20)           ::  region_name       ! glb100x100, eur050x050, ...
    character(len=20)           ::  hybrid_name       ! sfc, ml60, ...
    type(TllGridInfo)           ::  lli               ! lat-lon grid structure
    type(TLevelInfo)            ::  levi              ! hybrid level grid structure
    ! * temporal ax:
    integer                     ::  year              ! year (in TIPP, thus 0000 for year-independent)
    character(len=80)           ::  p_def             ! definition of time dimension
    character(len=5), pointer   ::  p_grid(:)         ! time grid
    type(TDate), pointer        ::  t1(:)             ! start of period (np)
    type(TDate), pointer        ::  t2(:)             ! end of period (np)
    integer                     ::  valid_period_range(2)  ! periods within this range have been read
    ! * original data
    character(len=MAX_FILENAME_LEN)          ::  original_filename ! original (ascii) name
    integer                     ::  nx_org            ! original longitudinal dimension
    integer                     ::  ny_org            ! original latitudinal dimension
    integer                     ::  nz_org            ! original vertical dimension
    integer                     ::  np_org            ! original time dimension (nr. of periods)
    character(len=32)           ::  filetype          ! file type ('E', 'G', ..)
    character(len=50)           ::  source            ! data source
    character(len=50)           ::  species           ! species
    character(len=50)           ::  category          ! emission category
    character(len=DUMMY_STR_LEN)          ::  category_long     ! description of emission category
    character(len=50)           ::  original_unit     ! unit in original file
    ! * full storage
    real, pointer               ::  field(:,:,:,:)    ! emission field
    logical, pointer            ::  valid(:,:,:,:)    ! fortran mask: true if data is valid
    ! * sparse storage:
    integer                     ::  sparse_n          ! number of values
    real, pointer               ::  sparse_field(:)   ! actual values (sparse_n)
    integer, pointer            ::  sparse_index(:,:) ! indices in x,y,z,t arrays (4,sparse_n)
    ! * units
    character(len=50)           ::  unit              ! unit in output file
    real                        ::  xm_fm             ! mole mass full molecule
    real                        ::  xm_el             ! mole mass main element
    character(len=50)           ::  xm_unit           ! unit of mole mass
    real                        ::  xm_undefined      ! value of xm in case undefined
    ! * scaling factor
    real                        ::  scaling_factor    ! scaling factor applied to original field
    ! * diagnostic: total
    real                        ::  total             ! total of field
    character(len=50)           ::  total_unit        ! units of total
    ! * other
    character(len=DUMMY_STR_LEN)          ::  additional_comments  ! general comments
  end type TEmis

  ! General emission unit is amount(amount_type)/area/timeframe
  type TEmUnit
     character(len=8)              :: amount
     character(len=8)              :: amount_type
     character(len=8)              :: area
     character(len=8)              :: timeframe
  end type TEmUnit

  ! General structure for boundary condition field incl metadata
  type TBouCon
     real, pointer       :: field(:,:,:,:)  ! emission field
     integer             :: nx      ! longitudinal dimension
     integer             :: ny      ! latitudinal dimension
     integer             :: nz      ! vertical dimension
     integer             :: np      ! time dimension (nr. of periods)
     character(len=80)   :: x_def   ! definition of long. dimension
     character(len=80)   :: y_def   ! definition of lat. dimension
     character(len=80)   :: z_def   ! definition of vert. dimension
     character(len=80)   :: p_def   ! definition of time dimension
     real, pointer       :: x_grid(:)  ! long. grid
     real, pointer       :: y_grid(:)  ! lat. grid
     real, pointer       :: z_a_grid(:)  ! hybrid grid, a coeff
     real, pointer       :: z_b_grid(:)  ! hybrid grid, b coeff
     character(len=4), pointer  :: p_grid(:)  ! time grid
     character(len=20)   :: region_name       ! glb100x100, eur050x050, ...
     character(len=20)   :: hybrid_name       ! sfc, ml60, ...
     type(TllGridInfo)   :: lli               ! lat-lon grid structure
     type(TLevelInfo)    :: levi              ! hybrid level grid structure
     character(len=50)   :: name              ! file name
     character(len=50)   :: original_filename ! original (ascii) name
     integer             :: nx_org      ! original longitudinal dimension
     integer             :: ny_org      ! original latitudinal dimension
     integer             :: nz_org      ! original vertical dimension
     integer             :: np_org      ! original time dimension (nr. of periods)
     character(len=50)   :: description        ! description of data set
     character(len=DUMMY_STR_LEN)  :: description_long   ! long description
     character(len=50)   :: source            ! data source
     character(len=50)   :: unit              ! units of field
     real                :: scaling_factor    ! scaling factor applied to original field
     integer             :: year              ! year
     character(len=DUMMY_STR_LEN)  :: additional_comments  ! general comments
  end type TBouCon

  ! --- interfaces -------------------------------------

  interface Init
     module procedure Init_ll
     module procedure Init_hyb
     module procedure Init_xm
  end interface


contains


  !----------------------------------------------------------
  ! Define latlon grid based on region_name (e.g. glb100x100)
  !
  ! Optionally make 1D grid
  !  o no_lon  :  no longitudinal dimension (nx=1)
  !  o no_lat  :  no latitudinal dimension (ny=1)
  !----------------------------------------------------------

  subroutine Init_ll( region_name, no_lon, no_lat, lli, status )

    ! --- modules --------------------------------

    ! --- in/out ----------------------------------------------

    character(len=*), intent(in)   :: region_name
    logical, intent(in)            :: no_lon, no_lat
    type(TllGridInfo), intent(out) :: lli
    integer, intent(out)           :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//', Init_ll'

    ! --- local -----------------------------------------------

    integer                        :: nx, ny
    real                           :: xbeg, xend, ybeg, yend, dx, dy

    ! --- begin -----------------------------------------------

    select case ( region_name )
    case ( 'glb600x400' )
       nx = 60
       ny = 45
       xbeg = -180.
       xend = 180.
       ybeg = -90.
       yend = 90.
    case ( 'glb100x100' )
       nx = 360
       ny = 180
       xbeg = -180.
       xend = 180.
       ybeg = -90.
       yend = 90.
    case ( 'eur100x100' )
       nx = 60
       ny = 54
       xbeg = -21.
       xend = 39.
       ybeg = 12.
       yend = 66.
    case ( 'eur050x050' )
       nx = 120
       ny = 108
       xbeg = -21.
       xend = 39.
       ybeg = 12.
       yend = 66.
    case default
       write (*,'("ERROR - unknown region name: ",a)') trim(region_name)
       write (*,'("ERROR in ",a)') name; status=1; return
    end select

    ! Check whether 1D grid
    if ( no_lon ) nx = 1
    if ( no_lat ) ny = 1

    dx = (xend - xbeg) / nx
    dy = (yend - ybeg) / ny

    call Init( lli, xbeg + dx/2., dx, nx, ybeg + dy/2., dy, ny, status )
    if ( status /= 0. ) then
       write (*,'("ERROR in ",a)') name; status=1; return
    end if

    status = 0

  end subroutine Init_ll


  subroutine Init_hyb( hybrid_name, levi, status )

    !----------------------------------------------------------
    ! Define hybrid grid based on name (e.g. 'ml60')
    !----------------------------------------------------------

    ! --- modules --------------------------------

    ! --- in/out ----------------------------------------------

    character(len=*), intent(in)   :: hybrid_name
    type(TLevelInfo), intent(out)  :: levi
    integer, intent(out)           :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//', Init_hyb'

    ! --- local -----------------------------------------------

    character(len=4)               :: levi_key

    ! --- begin -----------------------------------------------

    select case ( hybrid_name )
    case ( 'sfc' )
       levi_key = 'tm60'   ! Initialize 60-level grid as dummy
    case ( 'ml60' )
       levi_key = 'tm60'
    case ( 'ml91' )
       levi_key = 'tm91'
    case default
       write (*,'("ERROR - unknown hybrid grid name: ",a)') trim(hybrid_name)
       write (*,'("ERROR in ",a)') name; status=1; return
    end select

    call Init( levi, levi_key, status )
    if ( status /= 0. ) then
       write (*,'("ERROR in ",a)') name; status=1; return
    end if

    status = 0

  end subroutine Init_hyb


  subroutine Init_xm( species, xm_fm, xm_el, status )

    !----------------------------------------------------------
    ! Get mole mass of species (g mol-1)
    !----------------------------------------------------------

    ! --- modules --------------------------------

    ! --- in/out ----------------------------------------------

    character(len=*), intent(in)   :: species
    real, intent(out)              :: xm_fm, xm_el
    integer, intent(out)           :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//', Init_xm'

    ! --- local -----------------------------------------------

    integer                        :: i, ispec

    ! --- begin -----------------------------------------------

    ispec = 0
    do i = 1, nspec
       if ( species == spec_names(i) ) ispec = i
    end do
    if ( ispec == 0 ) then
       write (*,'("ERROR - unknown species `",a,"`")') species
       write (*,'("ERROR in ",a)') name; status=1; return
    end if

    xm_fm = spec_xm_fm(ispec)
    xm_el = spec_xm_el(ispec)

    status = 0

  end subroutine Init_xm


  ! *************************************************************************


  subroutine Emis_Sparse_Init( em, status )

    ! --- in/out ---------------------------------

    type(TEmis), intent(inout)      ::  em
    integer, intent(out)            ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//'/Emis_Sparse_Init'

    ! --- local -----------------------------------------------

    ! --- begin -----------------------------------------------

    ! nothing stored yet:
    em%sparse_n = 0
    nullify( em%sparse_field )
    nullify( em%sparse_index )

    ! ok
    status = 0

  end subroutine Emis_Sparse_Init


  ! ***


  subroutine Emis_Sparse_Done( em, status )

    ! --- in/out ---------------------------------

    type(TEmis), intent(inout)      ::  em
    integer, intent(out)            ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//'/Emis_Sparse_Done'

    ! --- local -----------------------------------------------

    ! --- begin -----------------------------------------------

    ! clear:
    if ( associated(em%sparse_field) ) nullify( em%sparse_field )
    if ( associated(em%sparse_index) ) nullify( em%sparse_index )
    ! nothing stored anymore:
    em%sparse_n = 0

    ! ok
    status = 0

  end subroutine Emis_Sparse_Done


  ! ***


  subroutine Emis_Sparse_Add( em, field, valid, ip, status )

    ! --- in/out ---------------------------------

    type(TEmis), intent(inout)      ::  em
    real, intent(in)                ::  field(:,:,:)   ! (nx,ny,nz) values
    logical, intent(in)             ::  valid(:,:,:)   ! (nx,ny,nz) valid data ?
    integer, intent(in)             ::  ip             ! period index
    integer, intent(out)            ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//'/Emis_Sparse_Add'

    ! --- local -----------------------------------------------

    integer             ::  n
    integer             ::  ix, iy, iz
    real, pointer       ::  new_sparse_field(:)
    integer, pointer    ::  new_sparse_index(:,:)

    ! --- begin -----------------------------------------------

    ! any positive values ?
    if ( any(valid) ) then

      ! number of positive values:
      n = count( valid )

      ! extend arrays:
      allocate( new_sparse_field(  em%sparse_n+n) )
      allocate( new_sparse_index(4,em%sparse_n+n) )
      ! old data present ?
      if ( em%sparse_n > 0 ) then
        ! copy old data:
        new_sparse_field(  1:em%sparse_n) = em%sparse_field(  1:em%sparse_n)
        new_sparse_index(:,1:em%sparse_n) = em%sparse_index(:,1:em%sparse_n)
        ! clear old arrays:
        deallocate( em%sparse_field )
        deallocate( em%sparse_index )
        ! reset for safety:
        nullify( em%sparse_field )
        nullify( em%sparse_index )
      end if
      ! point to new arrays:
      em%sparse_field => new_sparse_field
      em%sparse_index => new_sparse_index
      ! for safety ...
      nullify( new_sparse_field )
      nullify( new_sparse_index )

      ! loop over all cells:
      do iz = 1, size(field,3)
        do iy = 1, size(field,2)
          do ix = 1, size(field,1)
            ! valid ?
            if ( valid(ix,iy,iz) ) then
              ! increase counter:
              em%sparse_n = em%sparse_n + 1
              ! add value:
              em%sparse_field(em%sparse_n) = field(ix,iy,iz)
              ! add index:
              em%sparse_index(:,em%sparse_n) = (/ix,iy,iz,ip/)
            end if  ! non zero
            !
          end do  ! iz
        end do  ! iy
      end do  ! ix

    end if  ! any non zero

    ! update ip range:
    !  o first value:
    if ( em%valid_period_range(1) < 0 ) em%valid_period_range(1) = ip
    !  o latest vauue:
    em%valid_period_range(2) = max( em%valid_period_range(2), ip )

    ! ok
    status = 0

  end subroutine Emis_Sparse_Add



end module tipp_base
