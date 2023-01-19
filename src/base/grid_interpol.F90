!
! NAME
!   grid_interpol  -  from sh/gg/ll to gg/ll
!
! HISTORY
!   v1.?
!     Original
!   v1.1
!     Bug fixed: southpole in gg -> ll
!

module Grid_Interpol

  implicit none

  ! --- in/out ------------------------------

  private

  public  ::  Interpol
  public  ::  NewInterpol
  public  ::  InterpolMask
  public  ::  Aver
  public  ::  IntArea
  public  ::  IntLat, IntLon

  !public  ::  ShTruncation, ShRefinement

  public  ::  Tgg2llFracs, Tll2ggFracs, Init, Done, FracSum


  ! --- const ----------------------------------

  character(len=*), parameter  ::  mname = 'grid_interpol'


  ! --- types -------------------------------

  type Tgg2llFracs
    integer           ::  nlon, nlat, np
    integer, pointer  ::  ncov(:,:)
    integer, pointer  ::  indx(:,:,:)
    real, pointer     ::  frac(:,:,:)
    real, pointer     ::  A_gg(:)          ! m2
    real, pointer     ::  A_ll(:,:)        ! m2
  end type Tgg2llFracs


  type Tll2ggFracs
    integer           ::  nlon, nlat, np
    integer, pointer  ::  ncov(:)          ! (1:np)
    integer, pointer  ::  ii(:,:)          ! (1:np,1:max(ncov(k)))
    integer, pointer  ::  jj(:,:)          ! (1:np,1:max(ncov(k)))
    real, pointer     ::  ff(:,:)          ! (1:np,1:max(ncov(k)))
    real, pointer     ::  A_gg(:)          ! m2
    real, pointer     ::  A_ll(:,:)        ! m2
  end type Tll2ggFracs


  ! --- interfaces --------------------------

  interface Init
    module procedure gg2ll_Init
    module procedure ll2gg_Init
  end interface

  interface Done
    module procedure gg2ll_Done
    module procedure ll2gg_Done
  end interface

  interface FracSum
    module procedure gg2ll_FracSum
    module procedure ll2gg_FracSum
  end interface

  interface Interpol
    module procedure Interpol_sh_gg
    module procedure Interpol_shi_gg
    module procedure Interpol_sh_ll
    module procedure Interpol_shi_ll
    module procedure Interpol_gg_ll
    module procedure Interpol_ll_gg
  end interface

  interface NewInterpol
    module procedure NewInterpol_shi_gg
  end interface

  interface Aver
    module procedure Aver_gg_ll
    module procedure Aver_sh_ll
  end interface

  interface IntArea
    module procedure IntArea_sh_ll_f
    module procedure IntArea_shi_ll_f
    module procedure IntArea_sh_ll_fgh
    module procedure IntArea_shi_ll_fgh
    module procedure IntArea_sh_ll_fh
    module procedure IntArea_shi_ll_fh
  end interface

  interface IntLat
    module procedure IntLat_sh_ll
  end interface

  interface IntLon
    module procedure IntLon_sh_ll
  end interface



contains



  ! =========================================================
  ! ===
  ! === evaluate spectral fields
  ! ===
  ! =========================================================


  ! from spectral to reduced gausian grid


  subroutine Interpol_sh_gg( sh, ggi, gg, status )

    use Grid_Type_sh, only : TshGrid, Eval_Lons, Check
    use Grid_Type_gg, only : TggGridInfo, Check

    ! --- in/out ----------------------------------

    type(TshGrid), intent(in)       ::  sh
    type(TggGridInfo), intent(in)   ::  ggi
    real, intent(out)               ::  gg(ggi%np)
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/Interpol_sh_gg'

    ! --- local -----------------------------------

    !real, allocatable       ::  llgrid(:,:)
    real, allocatable       ::  llgrid(:)
    integer                 ::  nlon
    integer                 ::  jn !, js

    ! --- begin -----------------------------------

    call Check( sh )
    call Check( ggi, gg )

    !allocate( llgrid(maxval(ggi%nlon),2) )
    allocate( llgrid(maxval(ggi%nlon)) )

    ! northern rows:
    !do jn = 1, ggi%nlat/2
    do jn = 1, ggi%nlat

      ! southern row:
      !js = ggi%nlat + 1 - jn

      ! only if one of the rows is marked:
      !if ( ggi%latflag(jn) .or. ggi%latflag(js) ) then
      if ( ggi%latflag(jn) ) then

        nlon = ggi%nlon(jn)

        !call Eval_Lons( llgrid(1:nlon,1:2), sh, ggi%lat(jn), nlon, 0.0, nlon )
        !gg(ggi%i1(jn):ggi%im(jn)) = llgrid(1:nlon,1)
        !gg(ggi%i1(js):ggi%im(js)) = llgrid(1:nlon,2)

        call Eval_Lons( llgrid(1:nlon), sh, ggi%lat(jn), nlon, 0.0, nlon, status )
        gg(ggi%i1(jn):ggi%im(jn)) = llgrid(1:nlon)

      end if

    end do

    deallocate( llgrid )

  end subroutine Interpol_sh_gg


  ! *

  subroutine NewInterpol_shi_gg( shi, shc, ggi, gg, status )

    use Grid_Type_sh, only : TshGridInfo, Check, Eval_Lons
    use Grid_Type_sh, only : TshGrid, Init, Done, Set
    use Grid_Type_gg, only : TggGridInfo, Check

    ! --- in/out ----------------------------------

    type(TshGridInfo), intent(in)   ::  shi
    complex, intent(in)             ::  shc(shi%np)
    type(TggGridInfo), intent(in)   ::  ggi
    real, intent(out)               ::  gg(ggi%np)
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/NewInterpol_shi_gg'

    ! --- local -----------------------------------

    real, pointer           ::  llgrid(:)
    integer                 ::  nlon
    integer                 ::  jn !, js

    ! --- begin -----------------------------------

    allocate( llgrid(maxval(ggi%nlon)) )

    ! northern rows:
    !do jn = 1, ggi%nlat/2
    do jn = 1, ggi%nlat

      ! southern row:
      !js = ggi%nlat + 1 - jn

      ! only if one of the rows is marked:
      !if ( ggi%latflag(jn) .or. ggi%latflag(js) ) then
      if ( ggi%latflag(jn) ) then

        nlon = ggi%nlon(jn)

        !call Eval_Lons( llgrid(1:nlon,1:2), sh, ggi%lat(jn), nlon, 0.0, nlon )
        !gg(ggi%i1(jn):ggi%im(jn)) = llgrid(1:nlon,1)
        !gg(ggi%i1(js):ggi%im(js)) = llgrid(1:nlon,2)

        call Eval_Lons( shi, shc, ggi%lat(jn), nlon, 0.0, nlon, llgrid(1:nlon), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        gg(ggi%i1(jn):ggi%im(jn)) = llgrid(1:nlon)

      end if

    end do

    ! done
    deallocate( llgrid )

    ! ok
    status = 0

  end subroutine NewInterpol_shi_gg

  ! *


  subroutine Interpol_shi_gg( shi, shc, ggi, gg, status )

    use Grid_Type_sh, only : TshGridInfo, Check, Eval_Lons
    use Grid_Type_sh, only : TshGrid, Init, Done, Set
    use Grid_Type_gg, only : TggGridInfo, Check

    ! --- in/out ----------------------------------

    type(TshGridInfo), intent(in)   ::  shi
    complex, intent(in)             ::  shc(shi%np)
    type(TggGridInfo), intent(in)   ::  ggi
    real, intent(out)               ::  gg(ggi%np)
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/Interpol_shi_gg'

    ! --- local -----------------------------------

    type(TshGrid)           ::  sh

    !real, allocatable       ::  llgrid(:,:)
    real, allocatable       ::  llgrid(:)
    integer                 ::  nlon
    integer                 ::  jn !, js

    ! --- begin -----------------------------------

    ! store input in old type grid:
    call Init( sh )
    call Set( sh, shi%T, shc )

    !allocate( llgrid(maxval(ggi%nlon),2) )
    allocate( llgrid(maxval(ggi%nlon)) )

    ! northern rows:
    !do jn = 1, ggi%nlat/2
    do jn = 1, ggi%nlat

      ! southern row:
      !js = ggi%nlat + 1 - jn

      ! only if one of the rows is marked:
      !if ( ggi%latflag(jn) .or. ggi%latflag(js) ) then
      if ( ggi%latflag(jn) ) then

        nlon = ggi%nlon(jn)

        !call Eval_Lons( llgrid(1:nlon,1:2), sh, ggi%lat(jn), nlon, 0.0, nlon )
        !gg(ggi%i1(jn):ggi%im(jn)) = llgrid(1:nlon,1)
        !gg(ggi%i1(js):ggi%im(js)) = llgrid(1:nlon,2)

        call Eval_Lons( llgrid(1:nlon), sh, ggi%lat(jn), nlon, 0.0, nlon, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        gg(ggi%i1(jn):ggi%im(jn)) = llgrid(1:nlon)

      end if

    end do

    ! done
    deallocate( llgrid )
    call Done( sh )

    ! ok
    status = 0

  end subroutine Interpol_shi_gg


  ! ===


  ! from sh to ll

  subroutine Interpol_sh_ll( sh, lli, ll, status )

    use Grid_Type_sh, only : TshGrid, Eval_Lons, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -------------------------------

    type(TshGrid), intent(in)        ::  sh
    type(TllGridInfo), intent(in)    ::  lli
    real, intent(out)                ::  ll(lli%im,lli%jm)
    integer, intent(out)             ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/Interpol_sh_ll'

    ! --- local --------------------------------

    integer              ::  j

    ! --- begin --------------------------------

    call Check( sh )

    call Check( lli, 'n', ll, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! loop over latitudes
    do j = 1, lli%jm

      ! evaluate at oposite latidutes:
      call Eval_Lons( ll(:,j), sh, lli%lat(j), int(360.0/lli%dlon_deg), &
                      lli%lon(1), lli%im, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    end do

    ! ok
    status = 0

  end subroutine Interpol_sh_ll


  ! *


  ! from shi/sh to ll

  subroutine Interpol_shi_ll( shi, shc, lli, ll, status )

    use Grid_Type_sh, only : TshGridInfo, Check, Eval_Lons
    use Grid_Type_sh, only : TshGrid, Init, Done, Set
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -------------------------------

    type(TshGridInfo), intent(in)    ::  shi
    complex, intent(in)              ::  shc(shi%np)
    type(TllGridInfo), intent(in)    ::  lli
    real, intent(out)                ::  ll(lli%im,lli%jm)
    integer, intent(out)             ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/Interpol_shi_ll'

    ! --- local --------------------------------

    type(TshGrid)        ::  sh
    integer              ::  j

    ! --- begin --------------------------------

    ! store input in old type grid:
    call Init( sh )
    call Set( sh, shi%T, shc )

    ! check lat/lon arguments:
    call Check( lli, 'n', ll, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! loop over latitudes
    do j = 1, lli%jm

      ! evaluate spectral field:
      call Eval_Lons( ll(:,j), sh, lli%lat(j), int(360.0/lli%dlon_deg), &
                      lli%lon(1), lli%im, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    end do

    ! done
    call Done( sh )

    ! ok
    status = 0

  end subroutine Interpol_shi_ll



  ! =========================================================
  ! ===
  ! === interpolate from gg
  ! ===
  ! =========================================================

  ! Return a logical array with size gg;
  ! true if a point is required for interpolation from gg to ll
  !
  ! Include 'depth' extra rows.
  !
  ! Also set flags in ggi for each row.
  !

  subroutine InterpolMask( ggi, gg, lli, depth )

    use Grid_Type_gg, only : TggGridInfo, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -------------------------------

    type(TggGridInfo), intent(inout) ::  ggi
    logical, intent(out)             ::  gg(ggi%np)
    type(TllGridInfo), intent(in)    ::  lli
    integer, intent(in)              ::  depth

    ! --- local --------------------------------

    integer        ::  j, jn, js

    ! --- begin ------------------------------

    ! note: gg lats from north to south!

    !print *, 'llilat:',lli%blat_deg(0), lli%blat_deg(lli%nlat)

    ! search south most gg lat outside lli
    js = 1
    do
      if ( js == ggi%nlat ) exit
      if ( ggi%lat(js) < lli%blat(0) ) exit
      js = js + 1
    end do
    js = min( ggi%nlat, js+depth )

    ! search north-most gg lat row outside lli
    jn = ggi%nlat
    do
      if ( jn == 1 ) exit
      if ( ggi%lat(jn) > lli%blat(lli%nlat) ) exit
      jn = jn - 1
    end do
    jn = max( 1, jn-depth )

    ! set all rows between js and jn to true:
    gg = .false.
    ggi%latflag = .false.
    do j = jn, js
      ggi%latflag(j) = .true.
      gg(ggi%i1(j):ggi%im(j)) = .true.
    end do

  end subroutine InterpolMask


  ! *************************************************************

  !
  ! Determine fraction of gg cell that covers a ll cell.
  !
  ! Returns three arrays:
  !
  !   integer  ::  ncov(im,jm)
  !   integer  ::  indx(im,jm,max_nconv)
  !   real     ::  frac(im,jm,max_nconv)
  !
  ! For an ll cell (i,j), 'ncov' gives the number of gg cells covering it.
  ! Their gg indices are stored in : indx(i,j,1:ncov(i,j)),
  ! corresponding fractions are in : frac(i,j,1:ncov(i,j)) .
  !

  subroutine gg2ll_Init( gg2ll, ggi, lli, status )

    use num, only : Interval
    use grid_tools, only : pi, ll_area_frac
    use grid_type_gg, only : TggGridInfo, AreaOper
    use grid_type_ll, only : TllGridInfo, AreaOper

    ! --- in/out ---------------------------

    type(Tgg2llFracs), intent(out)     ::  gg2ll
    type(TggGridInfo), intent(in)      ::  ggi
    type(TllGridInfo), intent(in)      ::  lli
    integer, intent(out)               ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/gg2ll_Init'

    ! --- local ----------------------------

    integer    ::  max_ncov_lon, max_ncov_lat, max_ncov
    integer    ::  i, j
    integer    ::  iflag
    integer    ::  gi, gi1, gi2
    integer    ::  gj, gj1, gj2
    integer    ::  nlon
    real       ::  gblon(0:ggi%nlon_reg)
    real       ::  west1, east1, south1, north1
    real       ::  west2, east2, south2, north2
    real       ::  frac
    integer    ::  ncov

    ! --- begin ----------------------------

    ! save dimension
    gg2ll%nlon = lli%nlon
    gg2ll%nlat = lli%nlat
    gg2ll%np   = ggi%np

    ! estimate maximum number of gg cells covering an ll cell:
    max_ncov_lon = ceiling(lli%dlon/minval(ggi%dlon))
    max_ncov_lat = ceiling(lli%dlat/minval(ggi%dlat))
    max_ncov = (max_ncov_lon+1) * (max_ncov_lat+1)

    ! allocate arrays
    allocate( gg2ll%ncov(lli%nlon,lli%nlat) )
    allocate( gg2ll%indx(lli%nlon,lli%nlat,max_ncov) )
    allocate( gg2ll%frac(lli%nlon,lli%nlat,max_ncov) )

    ! zero by default:
    gg2ll%ncov = 0
    gg2ll%indx = 0
    gg2ll%frac = 0.0

    gj1 = 0
    gj2 = 0

    ! loop over ll cells from north to south (gg direction):
    do j = lli%nlat, 1, -1
      do i = 1, lli%nlon

        ! extract boundaries of ll cell:
        west2  = lli%blon(i-1)
        east2  = lli%blon(i)
        if ( east2 < 0.0 ) then
          west2 = west2 + 2*pi    ! (0,2pi)
          east2 = east2 + 2*pi    ! (0,2pi)
        end if
        south2 = lli%blat(j-1)
        north2 = lli%blat(j)

        ! search gg rows including north and south ll lat;
        ! negative lats to get increasing values ...
        call Interval( -ggi%blat, -north2  , gj1, iflag )
        if ( iflag /= 0 ) stop 'BUG IN gg2ll_Init : wrong iflag gj1'
        call Interval( -ggi%blat, -south2, gj2, iflag )
        if ( iflag /= 0 ) stop 'BUG IN gg2ll_Init : wrong iflag gj2'

        gi1 = 0
        gi2 = 0
        ! loop over gg lat rows:
        do gj = gj1, gj2

          ! boundary lons
          nlon = ggi%nlon(gj)
          do gi = 0, nlon
            gblon(gi) = (gi-0.5)*ggi%dlon(gj)
          end do

          ! search cells including west and east bound of ll cell
          if ( west2 < gblon(0) ) then
            call Interval( gblon(0:nlon), west2+2*pi, gi1, iflag )
          else if ( west2 > gblon(nlon) ) then
            call Interval( gblon(0:nlon), west2-2*pi, gi1, iflag )
          else
            call Interval( gblon(0:nlon), west2     , gi1, iflag )
          end if
          if ( iflag /= 0 ) then
            print *, 'gblon=', gblon(0:nlon)
            print *, 'west2=',west2
            print *, 'iflag=',iflag
            stop 'BUG IN gg2ll_Init : wrong iflag gi1'
          end if
          if ( east2 < gblon(0) ) then
            call Interval( gblon(0:nlon), east2+2*pi, gi2, iflag )
          else if ( east2 > gblon(nlon) ) then
            call Interval( gblon(0:nlon), east2-2*pi, gi2, iflag )
          else
            call Interval( gblon(0:nlon), east2     , gi2, iflag )
          end if
          if ( iflag /= 0 ) then
            print *, 'gblon=', gblon(0:nlon)
            print *, 'east2=',east2
            stop 'BUG IN gg2ll_Init : wrong iflag gi2'
          end if

          ! loop over all gg cells in current row:
          gi = gi1
          do

            ! extract boundaries of gg cell:
            west1  = (gi-1.5)*ggi%dlon(gj)    ! (0,2pi)
            east1  = (gi-0.5)*ggi%dlon(gj)    ! (0,2pi)
            south1 = ggi%blat(gj)
            north1 = ggi%blat(gj-1)

            ! shift if gg cell is right from [0,2pi]
            if ( west1 > east2 ) then
              west1 = west1 - 2*pi    ! (0,2pi)
              east1 = east1 - 2*pi    ! (0,2pi)
            end if

            ! determine covarage fraction:
            frac = ll_area_frac( west1, east1, south1, north1, &
                                 west2, east2, south2, north2 )

            ! fill fraction:
            if ( frac > 0.0 .and. frac <= 1.0 ) then
              ncov = gg2ll%ncov(i,j) + 1
              gg2ll%ncov(i,j) = ncov
              gg2ll%indx(i,j,ncov) = ggi%i1(gj)-1 + gi
              gg2ll%frac(i,j,ncov) = frac
            else if ( abs(frac) < 1.0e-4 ) then
              ! almost no coverage ...
            else if ( abs(frac-1.0) < 1.0e-4 ) then
              !print *, 'WARNING in module grid_interpol, gg2ll_Init'
              !print *, 'frac=',frac
              !print *, '  1:', west1, east1, south1, north1
              !print *, '  2:', west2, east2, south2, north2
              ncov = gg2ll%ncov(i,j) + 1
              gg2ll%ncov(i,j) = ncov
              gg2ll%indx(i,j,ncov) = ggi%i1(gj)-1 + gi
              gg2ll%frac(i,j,ncov) = 1.0
            else
              print *, 'ERROR in module grid_interpol, gg2ll_Init'
              print *, 'frac=',frac
              print *, '  1:', west1, east1, south1, north1
              print *, '  2:', west2, east2, south2, north2
              stop
            end if

            if ( gi == gi2 ) exit
            gi = gi + 1
            if ( gi == nlon+1 ) gi = 1
          end do

        end do

      end do
    end do

    ! store cell area's
    allocate( gg2ll%A_gg(ggi%np) )
    call AreaOper( ggi, gg2ll%A_gg, '=', 'm2' )

    allocate( gg2ll%A_ll(lli%nlon,lli%nlat) )
    call AreaOper( lli, gg2ll%A_ll, '=', 'm2', status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0

  end subroutine gg2ll_Init


  ! ***


  subroutine gg2ll_Done( gg2ll )

    ! --- in/out ---------------------------

    type(Tgg2llFracs), intent(inout)     ::  gg2ll

    ! --- begin ----------------------------

    ! deallocate arrays
    deallocate( gg2ll%ncov )
    deallocate( gg2ll%indx )
    deallocate( gg2ll%frac )
    deallocate( gg2ll%A_gg )
    deallocate( gg2ll%A_ll )

  end subroutine gg2ll_Done


  ! ***


  !
  !           ncov(i,j)
  ! ll(i,j) =   sum      gg(indx(i,j,k)) * frac(i,j,k)
  !             k=1
  !

  subroutine gg2ll_FracSum( gg2ll, gg, ll, status, key )

    ! --- in/out ---------------------------

    type(Tgg2llFracs), intent(in)      ::  gg2ll
    real, intent(in)                   ::  gg(:)
    real, intent(out)                  ::  ll(:,:)
    integer, intent(out)               ::  status

    character(len=*), intent(in), optional ::  key

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/gg2ll_FracSum'

    ! --- local ----------------------------

    character(len=10)  ::  the_key
    integer            ::  i, j, k
    integer            ::  ncov, indx
    real               ::  fac

    ! --- begin ----------------------------

    the_key = 'none'
    if ( present(key) ) the_key = key

    if ( any(shape(gg2ll%ncov)/=shape(ll)) ) then
      write (*,'("ERROR - shapes of ll and gg2ll do not match:")')
      write (*,'("ERROR -   ll    : ",2i4)') shape(ll)
      write (*,'("ERROR -   gg2ll : ",2i4)') shape(gg2ll%ncov)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    if ( minval(gg2ll%indx)<0 .or. maxval(gg2ll%indx)>size(gg) ) then
      write (*,'("ERROR - indices of gg array in gg2ll out of range:")')
      write (*,'("ERROR -   gg    :    1 - ",i4)') size(gg)
      write (*,'("ERROR -   gg2ll : ",i4," - ",i4)') minval(gg2ll%indx), maxval(gg2ll%indx)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    if ( minval(gg2ll%frac)<0.0 .or. maxval(gg2ll%frac)>1.0 ) then
      write (*,'("ERROR - fraction in gg2ll out of range:")')
      write (*,'("ERROR -   gg2ll: ",i4," - ",i4)') minval(gg2ll%frac), maxval(gg2ll%frac)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ll = 0.0
    do j = 1, gg2ll%nlat
      do i = 1, gg2ll%nlon

        ncov = gg2ll%ncov(i,j)
        if ( ncov > 0 ) then

          do k = 1, ncov
            indx = gg2ll%indx(i,j,k)

            select case ( the_key )
              case ( 'none' )
                fac = 1.0
              case ( 'area-aver', 'area-sum' )
                fac = gg2ll%A_gg(indx)
              case default
                write (*,'("ERROR - key `",a,"` not supported")') trim(the_key)
                write (*,'("ERROR in ",a)') rname; status=1; return
            end select

            ll(i,j) = ll(i,j) + gg(indx) * fac * gg2ll%frac(i,j,k)

          end do

        end if

        select case ( the_key )
          case ( 'none', 'area-sum' )
            ! nothing to be done
          case ( 'area-aver' )
            ll(i,j) = ll(i,j) / gg2ll%A_ll(i,j)
          case default
            write (*,'("ERROR - key `",a,"` not supported")') trim(the_key)
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

      end do
    end do

    ! ok
    status = 0

  end subroutine gg2ll_FracSum




  ! *************************************************************

  !
  ! Determine fraction of gg cell that covers a ll cell.
  !
  ! Returns three arrays:
  !
  !   integer  ::  ncov(im,jm)
  !   integer  ::  indx(im,jm,max_nconv)
  !   real     ::  frac(im,jm,max_nconv)
  !
  ! For a gg cell 'k', 'ncov(k)' gives the number of ll cells covering it.
  ! Their ll indices are stored in : ii(k,1:ncov(k)), jj(k,1:ncov(k)) ;
  ! corresponding fractions are in : frac(k,1:ncov(k)) .
  !

  subroutine ll2gg_Init( ll2gg, lli, ggi, status )

    use num, only : Interval
    use grid_tools, only : pi, ll_area_frac
    use grid_type_gg, only : TggGridInfo, AreaOper
    use grid_type_ll, only : TllGridInfo, AreaOper

    ! --- in/out ---------------------------

    type(Tll2ggFracs), intent(out)     ::  ll2gg
    type(TllGridInfo), intent(in)      ::  lli
    type(TggGridInfo), intent(in)      ::  ggi
    integer, intent(out)               ::  status

    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ll2gg_Init'

    ! --- local ----------------------------

    integer    ::  max_ncov_lon, max_ncov_lat, max_ncov
    integer    ::  gi, gi1, gi2
    integer    ::  gj, gj1, gj2
    integer    ::  li, li1, li2
    integer    ::  lj, lj1, lj2
    real       ::  west2, east2, south2, north2
    real       ::  west1, east1, south1, north1
    integer    ::  iflag
    integer    ::  gp
    real       ::  frac
    integer    ::  ncov

    ! --- begin ----------------------------

    ! save dimension
    ll2gg%nlon = lli%nlon
    ll2gg%nlat = lli%nlat
    ll2gg%np   = ggi%np

    ! estimate maximum number of ll cells covering a gg cell:
    max_ncov_lon = ceiling(maxval(ggi%dlon)/lli%dlon)
    max_ncov_lat = ceiling(maxval(ggi%dlat)/lli%dlat)
    max_ncov = (max_ncov_lon+1) * (max_ncov_lat+1)

    ! allocate arrays
    allocate( ll2gg%ncov(ggi%np) )
    allocate( ll2gg%ii(ggi%np,max_ncov) )
    allocate( ll2gg%jj(ggi%np,max_ncov) )
    allocate( ll2gg%ff(ggi%np,max_ncov) )

    ! zero by default:
    ll2gg%ncov = 0
    ll2gg%ii   = 0
    ll2gg%jj   = 0
    ll2gg%ff   = 0.0

    ! index in gg arrays:
    gp = 0

    ! row range in ll grid:
    lj1 = 0
    lj2 = 0

    ! loop over gg rows:
    do gj = 1, ggi%nlat

      ! extract north/south boundaries of gg cell:
      north2 = ggi%blat(gj-1)  ! [-pi,pi]
      south2 = ggi%blat(gj)    ! [-pi,pi]

      ! search ll rows including north and south gg lat;
      call Interval( lli%blat, south2, lj1, iflag )
      if ( iflag /= 0 ) stop 'BUG IN ll2gg_Init : wrong iflag lj1'
      call Interval( lli%blat, north2, lj2, iflag )
      if ( iflag /= 0 ) stop 'BUG IN ll2gg_Init : wrong iflag lj2'

      ! loop over cells in row:
      do gi = 1, ggi%nlon(gj)

        ! next index in 1d row ...
        gp = gp + 1

        ! set east/west boundaries of gg cell:
        west2 = (gi-1.5) * ggi%dlon(gj)  ! [0,2pi]
        east2 = (gi-0.5) * ggi%dlon(gj)  ! [0,2pi]

        ! loop over ll lat rows:
        do lj = lj1, lj2

          ! search ll cells including west and east bound of gg cell
          if ( west2 > lli%blon(lli%nlon) ) then
            call Interval( lli%blon, west2-2*pi, li1, iflag )
          else
            call Interval( lli%blon, west2     , li1, iflag )
          end if
          if ( iflag /= 0 ) then
            print *, 'lli%blon=', lli%blon
            print *, 'west2=',west2
            print *, 'iflag=',iflag
            stop 'BUG IN ll2gg_Init : wrong iflag li1'
          end if
          if ( east2 > lli%blon(lli%nlon) ) then
            call Interval( lli%blon, east2-2*pi, li2, iflag )
          else
            call Interval( lli%blon, east2     , li2, iflag )
          end if
          if ( iflag /= 0 ) then
            print *, 'lli%blon=', lli%blon
            print *, 'east2=',east2
            print *, 'east2-2pi=',east2-2*pi
            stop 'BUG IN ll2gg_Init : wrong iflag li2'
          end if

          ! loop over ll lon cels:
          li = li1
          do

            ! extract boundaries of ll cell:
            west1  = lli%blon(li-1)     ! [-pi,pi]
            east1  = lli%blon(li  )     ! [-pi,pi]
            south1 = lli%blat(lj-1)     ! [-pi,pi]
            north1 = lli%blat(lj  )     ! [-pi,pi]

            ! shift if completely left from gg cell:
            if ( east1 < west2 ) then
              west1 = west1 + 2*pi
              east1 = east1 + 2*pi
            end if

            ! determine covarage fraction:
            frac = ll_area_frac( west1, east1, south1, north1, &
                                 west2, east2, south2, north2 )

            ! fill fraction:
            if ( frac > 0.0 .and. frac <= 1.0 ) then
              ncov = ll2gg%ncov(gp) + 1
              ll2gg%ncov(gp)    = ncov
              ll2gg%ii(gp,ncov) = li
              ll2gg%jj(gp,ncov) = lj
              ll2gg%ff(gp,ncov) = frac
            else if ( abs(frac) < 1.0e-4 ) then
              ! almost no coverage ...
            else if ( abs(frac-1.0) < 1.0e-4 ) then
              !print *, 'WARNING in module grid_interpol, ll2gg_Init'
              !print *, 'frac=',frac
              !print *, '  1:', west1, east1, south1, north1
              !print *, '  2:', west2, east2, south2, north2
              ncov = ll2gg%ncov(gp) + 1
              ll2gg%ncov(gp)    = ncov
              ll2gg%ii(gp,ncov) = li
              ll2gg%jj(gp,ncov) = lj
              ll2gg%ff(gp,ncov) = 1.0
            else
              print *, 'ERROR in module grid_interpol, ll2gg_Init'
              print *, 'frac=',frac
              print *, '  1:', west1, east1, south1, north1
              print *, '  2:', west2, east2, south2, north2
              stop
            end if

            if ( li == li2 ) exit
            li = li + 1
            if ( li == lli%nlon+1 ) li = 1

          end do  ! ll i
        end do    ! ll j

      end do  ! gg i
    end do    ! gg j

    ! store cell area's
    allocate( ll2gg%A_gg(ggi%np) )
    call AreaOper( ggi, ll2gg%A_gg, '=', 'm2' )
    !
    allocate( ll2gg%A_ll(lli%nlon,lli%nlat) )
    call AreaOper( lli, ll2gg%A_ll, '=', 'm2', status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0

  end subroutine ll2gg_Init


  ! ***


  subroutine ll2gg_Done( ll2gg )

    ! --- in/out ---------------------------

    type(Tll2ggFracs), intent(inout)     ::  ll2gg

    ! --- begin ----------------------------

    ! deallocate arrays
    deallocate( ll2gg%ncov )
    deallocate( ll2gg%ii   )
    deallocate( ll2gg%jj   )
    deallocate( ll2gg%ff   )
    deallocate( ll2gg%A_gg )
    deallocate( ll2gg%A_ll )

  end subroutine ll2gg_Done


  ! ***


  !
  !          ncov(p)
  ! gg(p) =   sum      ll(ii(p,k),jj(p,k)) * ff(p,k)
  !           k=1
  !

  subroutine ll2gg_FracSum( ll2gg, ll, gg, key )

    ! --- in/out ---------------------------

    type(Tll2ggFracs), intent(in)      ::  ll2gg
    real, intent(in)                   ::  ll(:,:)
    real, intent(out)                  ::  gg(:)

    character(len=*), intent(in), optional ::  key

    ! --- local ----------------------------

    character(len=10)  ::  the_key
    integer            ::  gp, i, j, k
    integer            ::  ncov
    real               ::  fac

    ! --- begin ----------------------------

    the_key = 'none'
    if ( present(key) ) the_key = key

    if ( size(ll2gg%ncov) /= size(gg) ) then
      print *, 'shapes of gg and ll2gg do not match:'
      print *, '  ll   : ', size(gg)
      print *, '  ll2gg: ', size(ll2gg%ncov)
      stop 'FATAL ERROR IN ll2gg_FracSum'
    end if
    if ( minval(ll2gg%ii)<0 .or. maxval(ll2gg%ii)>size(ll,1) ) then
      print *, 'indices of ll array in ll2gg out of range:'
      print *, '  ll i : 1 - ', size(ll,1)
      print *, '  ll2gg: ', minval(ll2gg%ii),'-',maxval(ll2gg%ii)
      stop 'FATAL ERROR IN ll2gg_FracSum'
    end if
    if ( minval(ll2gg%jj)<0 .or. maxval(ll2gg%jj)>size(ll,2) ) then
      print *, 'indices of ll array in ll2gg out of range:'
      print *, '  ll j : 1 - ', size(ll,2)
      print *, '  ll2gg: ', minval(ll2gg%jj),'-',maxval(ll2gg%jj)
      stop 'FATAL ERROR IN ll2gg_FracSum'
    end if
    if ( minval(ll2gg%ff)<0.0 .or. maxval(ll2gg%ff)>1.0 ) then
      print *, 'fraction in ll2gg out of range:'
      print *, '  ll2gg: ', minval(ll2gg%ff),'-',maxval(ll2gg%ff)
      stop 'FATAL ERROR IN ll2gg_FracSum'
    end if

    ! init to zero:
    gg = 0.0

    ! loop over gg
    do gp = 1, ll2gg%np

      ncov = ll2gg%ncov(gp)
      if ( ncov > 0 ) then

        do k = 1, ncov
          i = ll2gg%ii(gp,k)
          j = ll2gg%jj(gp,k)

          select case ( the_key )
            case ( 'none' )
              fac = 1.0
            case ( 'area-aver', 'area-sum' )
              fac = ll2gg%A_ll(i,j)
            case default
              print *, 'Sorry, key "'//trim(the_key)//'" not supported.'
              stop 'BUG IN ll2gg_FracSum'
          end select
!if (abs(ll(i,j))>0.0) then
!print *, gp,':',gg(gp) ,'+', ll(i,j) ,'*', fac ,'*', ll2gg%ff(gp,k),'=',gg(gp) + ll(i,j) * fac * ll2gg%ff(gp,k)
!endif
          gg(gp) = gg(gp) + ll(i,j) * fac * ll2gg%ff(gp,k)

        end do

      end if

      select case ( the_key )
        case ( 'none', 'area-sum' )
          ! nothing to be done
!if (abs(gg(gp))>0.0) then
!print *, ' '
!endif
        case ( 'area-aver' )
          gg(gp) = gg(gp) / ll2gg%A_gg(gp)
!if (abs(gg(gp))>0.0) then
!print *, gp,':',gg(gp),'/',ll2gg%A_gg(gp),'=',gg(gp) / ll2gg%A_gg(gp)
!print *, ' '
!endif
        case default
          print *, 'Sorry, key "'//trim(the_key)//'" not supported.'
          stop 'BUG IN ll2gg_FracSum'
      end select

    end do   ! gg cells

  end subroutine ll2gg_FracSum




  ! *************************************************************


  ! gg to ll

  subroutine Interpol_gg_ll( ggi, gg, lli, ll, status )

    use Binas, only : deg2rad

    use Num, only : Interp_Lin, CircInterp_Lin

    use Grid_Type_gg, only : TggGridInfo, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -------------------------------

    type(TggGridInfo), intent(in)    ::  ggi
    real, intent(in)                 ::  gg(ggi%np)
    type(TllGridInfo), intent(in)    ::  lli
    real, intent(out)                ::  ll(lli%im,lli%jm)
    integer, intent(out)             ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Interpol_gg_ll'

    ! --- local --------------------------------

    integer         ::  nlon, nlon_max
    integer         ::  nlat
    integer         ::  i1, im

    real, allocatable     ::  lons(:)
    real, allocatable     ::  row(:)

    real, allocatable     ::  gl(:,:)
    real, allocatable     ::  gl_lat(:)

    integer               ::  i, j
    integer               ::  j1, jm
    integer               ::  ilast

    ! --- begin --------------------------------

    call Check( ggi, gg )

    call Check( lli, 'n', ll, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    nlat = ggi%nlat
    nlon_max = maxval(ggi%nlon)

    ! ll in lon, gg in lat
    allocate( gl(size(ll,1),0:nlat+1) ); gl = 0.0

    ! latitudes from south->north
    allocate( gl_lat(0:nlat+1) )
    gl_lat(0) = -90.0 * deg2rad       ! south pole (rad)
    do j = 1, nlat
      gl_lat(nlat+1-j) = ggi%lat(j)   ! rad
    end do
    gl_lat(nlat+1) = 90.0 * deg2rad   ! north pole (rad)

    ! row in gg grid; doubled from -360.0 to 360.0
    allocate( lons(nlon_max) ); lons = 0.0
    allocate( row(nlon_max) ) ; row  = 0.0

    ! select first and last Gaussian lat:
    j1 = nlat
    do
      if ( (j1 == 1) .or. (ggi%lat(j1) > maxval(lli%lat)) ) exit
      j1 = j1 - 1
    end do
    jm = 1
    do
      if ( (jm == nlat) .or. (ggi%lat(jm) < minval(lli%lat)) ) exit
      jm = jm + 1
    end do

    ! loop over Gaussian latitudes, from north to south!
    do j = j1, jm

      ! number of lon points at this latitude:
      nlon = ggi%nlon(j)

      ! start and end indices in grid point array
      i1 = ggi%i1(j)
      im = ggi%im(j)

      ! lons in [0,2pi)
      do i = 1, nlon
        lons(i) = (i-1)*360.0/nlon
      end do
      lons(1:nlon) = lons(1:nlon) * deg2rad

      ! grid values (doubled)
      row(1:nlon) = gg(i1:im)

      ! set north pole  (j=1 in gg, j=nlat+1 in gl)
      if ( j == 1 ) then
        gl(:,nlat+1) = sum(row(1:nlon))/nlon
      end if

      ! set south pole
      if ( j == nlat ) then
        gl(:,0) = sum(row(1:nlon))/nlon
      end if

      ! Interpolate over lon (circular arrays);
      ! swap lats to ensure south -> north:
      ilast = 1
      do i = 1, lli%im
        call CircInterp_Lin( lons(1:nlon), 360.0*deg2rad, row(1:nlon), &
                             lli%lon(i), gl(i,nlat+1-j), ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      end do

    end do

    ! Linear interpolation over lat:
    ilast = 1
    do j = 1, lli%jm
      do i = 1, lli%im
        call Interp_Lin( gl_lat, gl(i,:), lli%lat(j), ll(i,j), ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      end do
    end do

    ! free memory
    deallocate( gl )
    deallocate( gl_lat )
    deallocate( lons )
    deallocate( row )

    ! ok
    status = 0

  end subroutine Interpol_gg_ll


  ! =========================================================
  ! ===
  ! === interpolate from ll
  ! ===
  ! =========================================================


  ! ll to gg
  !
  ! First interpol in lat direction to Gaussian lats,
  ! then interpolate in lon direction.

  ! NOTE: MIPSpro compiler gives errors if argument names are the
  !       same as used for Interpol_gg_ll ...

  subroutine Interpol_ll_gg( lli, ll, xggi, xgg, status )

    use Binas, only : deg2rad

    use Num, only : Interp_Lin, CircInterp_Lin

    use Grid_Type_gg, only : TggGridInfo, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -------------------------------

    type(TllGridInfo), intent(in)    ::  lli
    real, intent(in)                 ::  ll(lli%im,lli%jm)
    type(TggGridInfo), intent(in)    ::  xggi
    real, intent(out)                ::  xgg(1:xggi%np)
    integer, intent(out)             ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Interpol_ll_gg'

    ! --- local --------------------------------

    integer               ::  nlon
    integer               ::  nlat

    integer               ::  i, j
    integer               ::  i0
    integer               ::  ilast

    real                  ::  dlon_deg
    real                  ::  glat, glon
    real                  ::  period
    real, allocatable     ::  row(:)

    ! --- begin --------------------------------

    call Check( xggi, xgg )

    call Check( lli, 'n', ll, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    nlat = xggi%nlat

    ! ll grid interpolated to gaussian lat:
    allocate( row(1:lli%im) )

    ! loop over Gaussian latitudes:
    do j = 1, nlat

      ! number of lon points at this latitude:
      nlon = xggi%nlon(j)

      ! current Gaussian lat
      glat = xggi%lat(j)

      !! error if ll grid does not cover this gaussian lat
      !if ( glat < minval(lli%lat) .or. glat > maxval(lli%lat) ) then
      !  write (*,'("ERROR - ll grid does not cover current gg latitude:")')
      !  write (*,'("ERROR -   ll lat range : ",2f12.4)') minval(lli%lat),  maxval(lli%lat)
      !  write (*,'("ERROR -   gg latitude  : ",f12.4)') glat
      !  write (*,'("ERROR in ",a)') rname; status=1
      !end if

      ! * interpolate ll grid to this gaussian lat:
      ! check for direction of ll lats:
      if ( lli%lat(1) < lli%lat(lli%jm) ) then
        ! 'normal' : grid stored from south to north
        ilast = 1
        do i = 1, lli%im
          call Interp_Lin( lli%lat, ll(i,:), glat, row(i), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do
      else
        ! grid stored from north to south; fake with negative lats:
        ilast = 1
        do i = 1, lli%im
          call Interp_Lin( -lli%lat, ll(i,:), -glat, row(i), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do
      end if

      ! base index of current row in xgg array:
      i0 = xggi%i1(j)-1

      ! interpolate from ll lons to lons in this xgg row:
      period = 360.0*deg2rad
      ilast = 1
      do i = 1, nlon
        ! curren lon in xgg grid:
        glon = (i-1)*xggi%dlon(j)
        ! periodic interpolation:
        call CircInterp_Lin( lli%lon, period, row, glon, xgg(i0+i), ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      end do

    end do

    ! free memory
    deallocate( row )

    ! ok
    status = 0

  end subroutine Interpol_ll_gg




  ! =========================================================
  ! ===
  ! === area average gg
  ! ===
  ! =========================================================


  subroutine Aver_gg_ll( ggi, gg, lli, ll )

    use Binas, only : deg2rad

    use Num, only : IntervalQuad_Lin, IntervalQuad_Cos_Lin

    use Grid_Type_gg, only : TggGridInfo, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -------------------------------

    type(TggGridInfo), intent(in)    ::  ggi
    real, intent(in)                 ::  gg(ggi%np)
    type(TllGridInfo), intent(in)    ::  lli
    real, intent(out)                ::  ll(lli%im,lli%jm)

    integer  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ Aver_gg_ll'

    ! --- local --------------------------------

    integer         ::  nlon, nlon_max
    integer         ::  nlat
    integer         ::  i1, im

    real, allocatable            ::  lons(:)
    real, allocatable            ::  row(:)

    real, allocatable            ::  gl(:,:)
    real, allocatable            ::  gl_lat(:)
    real, allocatable            ::  gl_dim2(:)

    integer         ::  i, j
    integer         ::  ilast

    ! --- begin --------------------------------

    call Check( ggi, gg )

    call Check( lli, 'n', ll, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if

    nlat = ggi%nlat
    nlon_max = maxval(ggi%nlon)

    ! ll in lon, gg in lat
    allocate( gl(size(ll,1),0:nlat+1) )
    allocate( gl_dim2(0:nlat+1) )

    ! latitudes from south->north
    allocate( gl_lat(0:nlat+1) )

    ! row in gg grid; doubled from -360.0 to 360.0
    allocate( lons(2*nlon_max) )
    allocate( row(2*nlon_max) )

    ! loop over Gaussian latitudes, from north to south!
    do j = 1, nlat

      ! number of lon points at this latitude:
      nlon = ggi%nlon(j)

      ! start and end indices in grid point array
      i1 = ggi%i1(j)
      im = ggi%im(j)

      ! lons from -360 to 360
      do i = 1, nlon
        lons(i) = -360.0 + (i-1)*360.0/nlon
        lons(nlon+i) = (i-1)*360.0/nlon
      end do
      lons(1:2*nlon) = lons(1:2*nlon) * deg2rad

      ! grid values (doubled)
      row(1:2*nlon) = (/ gg(i1:im), gg(i1:im) /)

      ! integrate over dlon assuming linear interpolation;
      ! swap lats to ensure south -> north;
      ! result in   [g] rad
      ilast = 1
      do i = 1, lli%im
        call IntervalQuad_Lin( lons(1:2*nlon), row(1:2*nlon), &
                               lli%blon(i-1), lli%blon(i), gl(i,nlat+1-j), &
                               ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if
      end do
      gl_lat(nlat+1-j) = ggi%lat(j)   ! rad

      ! take care of poles (reverse lats);
      ! set polar value to average of surrounding points:
      if ( j == 1 ) then
        ! north pole
        gl_lat(nlat+1) = 90.0 * deg2rad
        gl(:,nlat+1) = lli%dlon * sum(gg(i1:im))/(im-i1+1)
      end if
      if ( j == nlat ) then
        ! south pole
        gl_lat(0) = -90.0 * deg2rad
        gl(:,0) = lli%dlon * sum(gg(i1:im))/(im-i1+1)
      end if

    end do

    ! integrate over dlat assuming linear interpolation;
    ! weight with cos(lat) to account for smaller cells near the poles,
    ! result in  [g] rad^2
    ilast = 1
    do j = 1, lli%jm
      do i = 1, lli%im
        !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !call IntervalQuad_Cos_Lin( gl_lat, gl(i,:), lli%blat(j-1), lli%blat(j), ll(i,j), ilast )
        gl_dim2 = gl(i,:)
        call IntervalQuad_Cos_Lin( gl_lat, gl_dim2, lli%blat(j-1), lli%blat(j), &
                                         ll(i,j), ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      end do
    end do

    ! area average; lli%area(j) for cells in row j in rad^2
    do j = 1, lli%jm
      ll(:,j) = ll(:,j) / lli%area(j)
    end do

    ! free memory
    deallocate( gl )
    deallocate( gl_dim2 )
    deallocate( gl_lat )
    deallocate( lons )
    deallocate( row )

  end subroutine Aver_gg_ll


  ! =========================================================
  ! ===
  ! === area average sh
  ! ===
  ! =========================================================


  ! Deterimine spectral truncation for area integration.

  integer function ShTruncation( T, dlon, dlat )

    ! --- in/out --------------------------

    integer, intent(in)          ::  T
    real, intent(in)             ::  dlon, dlat   ! deg

    ! --- begin --------------------------

print *, 'WARNING - spectral fields are not truncated ...'
    ! o no truncation:
    ShTruncation = T

    ! o choose minium T based on grid resolution:
!    ShTruncation = min( T, ceiling( (360.0/min(dlon,dlat))/2 - 1 ) )

  end function ShTruncation


  ! Deterimine refinement for averaging spectral fields
  ! over distance 'cellspacing' in degrees.
  ! Cell is divided in number of intervals returned by this function.

  integer function ShRefinement( T, cellspacing )

    ! --- in/out ----------------------------------------

    integer, intent(in)     ::  T
    real, intent(in)        ::  cellspacing

    ! --- local -----------------------------------------

    real         ::  shres
    integer      ::  nstep

    ! --- beging -----------------------------------------

    ! o  fixec number of intervals per cell
    !ShRefinement = 20
    !ShRefinement = 1
    !ShRefinement = nint( cellspacing * 2 )

          !write (*,'("        WARNING: ShRefinement = ",i4)') ShRefinement

    ! o resultion specified by spectral truncation
    !     truncation T  <--> resolution  360.0/(2(T+1))
    !   nstep points within resolution implied by T
    shres = 360.0 / (2*(T+1))
    nstep = 2
    ShRefinement = max( 1, ceiling(nstep*cellspacing/shres) )

  end function ShRefinement


  ! ***********************************


  subroutine Aver_sh_ll( sh, lli, ll, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin, IntervalQuad_Cos_Lin

    use Grid_Type_sh, only : TshGrid, Check, Eval_Lons
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -------------------------------

    type(TshGrid), intent(in)        ::  sh
    type(TllGridInfo), intent(in)    ::  lli
    real, intent(out)                ::  ll(lli%im,lli%jm)
    integer, intent(out)             ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Aver_sh_ll'

    ! --- local --------------------------------

    integer               ::  i, j, jf
    integer               ::  ilast
    integer               ::  nlon_fine

    integer               ::  T
    integer               ::  refinement_i, refinement_j

    real, allocatable     ::  llgrid(:)
    real, allocatable     ::  lons_fine(:), row_fine(:)
    real, allocatable     ::  llf(:,:), lat(:)

    ! --- begin --------------------------------

    call Check( sh )
    T = sh%T

    call Check( lli, 'n', ll, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; stop; end if

    ! determine refinement (5 points per sh resolution)
    refinement_i = ShRefinement( T, lli%dlon_deg )
    refinement_j = ShRefinement( T, lli%dlon_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine = 360.0/lli%dlon_deg * refinement_i

    ! store evaluation of spectral field:
    allocate( llgrid(nlon_fine) )

    ! lons on complete circle from westb+[0,2pi)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine
    end do
    lons_fine = lli%blon(0) + lons_fine

    ! ll in lon, fine in lat
    allocate( llf(lli%im,0:refinement_j) )
    allocate( lat(0:refinement_j) )

    ! loop over latitudes in ll grid:
    do j = 1, lli%jm

      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        ! evaluate row:
        ! (oposite latitudes, but use only one)
        call Eval_Lons( llgrid, sh, lat(jf), &
                        nlon_fine, lons_fine(0), nlon_fine, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! copy result, cyclic:
        row_fine = (/ llgrid(:), llgrid(1) /)

        ! integral in lon direction assuming linear interpolation,
        ! result in  [sh] rad
        ilast = 1
        do i = 1, lli%im
          call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), llf(i,jf), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do

      end do

      ! integral in lat direction;
      ! weight with cos(lat) to account for smaller cells near the poles,
      ! result in  [sh] rad^2
      do i = 1, lli%im
        ilast = 1
        call IntervalQuad_Cos_Lin( lat, llf(i,:), lli%blat(j-1), lli%blat(j), ll(i,j), ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      end do

      ! area average; lli%area(j) for cells in row j in rad^2
      ll(:,j) = ll(:,j) / lli%area(j)

    end do

    ! free memory
    deallocate( llgrid )
    deallocate( lons_fine )
    deallocate( row_fine )
    deallocate( llf )
    deallocate( lat )

  end subroutine Aver_sh_ll



  ! ===========================================================

  !
  ! key == 'aver'
  !
  !   int int   F(k)  dA    /   A
  !
  ! key == 'exp,aver'
  !
  !   int int   exp(F(k))  dA    /   A
  !


  subroutine IntArea_shi_ll_f( key, shi, shc, lli, ll, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin, IntervalQuad_Cos_Lin

    use Grid_Type_sh, only : TshGrid, TshGridInfo
    use Grid_Type_sh, only : Init, Done, Check, Set, Truncate, SpN
    use Grid_Type_sh, only : sh_Pnm, Eval_Lons
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGridInfo), intent(in)      ::  shi
    complex, intent(in)                ::  shc(:)
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(out)                  ::  ll(lli%im,lli%jm)
    integer, intent(out)               ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntArea_shi_ll_f'

    ! --- local ------------------------------------------

    integer               ::  i, j, jf
    integer               ::  ilast

    type(TshGrid)         ::  sh
    integer               ::  Ttr
    real, allocatable     ::  Pnm(:)

    integer               ::  refinement_i, refinement_j

    integer               ::  nlon_fine_360, nlon_fine
    real, allocatable     ::  ff(:)
    real, allocatable     ::  lons_fine(:), row_fine(:)
    real, allocatable     ::  llf(:,:), lat(:)

    real, allocatable     ::  llf_dim2(:)

    !logical               ::  aver_to_prev

    ! --- begin ------------------------------------------

    ! store input in old type grid:
    call Init( sh )
    call Set( sh, shi%T, shc )

    ! use truncation up to grid resolution:
    Ttr = ShTruncation( shi%T, lli%dlon_deg, lli%dlat_deg )
    call Truncate( sh, Ttr )

    ! allocate space for legendre coeff:
    allocate( Pnm(SpN(Ttr)) )

    ! determine refinement (5 points per sh resolution)
    refinement_i = ShRefinement( Ttr, lli%dlon_deg )
    refinement_j = ShRefinement( Ttr, lli%dlat_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg * refinement_i
    nlon_fine = lli%im * refinement_i

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )

    ! lons on arc  westb+[0,..)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine_360
    end do
    lons_fine = lli%blon(0) + lons_fine

    ! ll in lon, fine in lat
    allocate( llf(lli%im,0:refinement_j) )
    allocate( llf_dim2(0:refinement_j) )
    allocate( lat(0:refinement_j) )

    ! loop over latitudes in ll grid:
    !aver_to_prev = .false.
    do j = 1, lli%jm

      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        !! southpole ?
        !if ( abs(lat(jf) - (-pi/2)) < 1.0e-4 ) then
        !  ! fill with average of next row
        !  aver_to_prev = .true.
        !  cycle
        !end if
        !
        !! northpole ?
        !if ( abs(lat(jf) - (pi/2)) < 1.0e-4 ) then
        !  ! fill with average of previous row:
        !  llf(:,jf) = sum(llf(:,jf-1)) / size(llf,1)
        !  exit
        !end if

        ! evaluate Legendre functions:
        call sh_Pnm( Pnm, Ttr, lat(jf), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! evaluate rows:
        call Eval_Lons( ff, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! combine fields:
        select case ( key )

          !
          !   int int   F(k)  dA    /   A
          !
          case ( 'aver' )

            row_fine = ff / lli%area(j)

          !
          !   int int   exp(F(k))  dA    /   A
          !
          case ( 'exp,aver' )

            row_fine = exp(ff) / lli%area(j)

          !
          ! error ...
          !
          case default
            write (*,'("ERROR - unknown key `",a,"`")') trim(key)
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select

        ! integral in lon direction assuming linear interpolation,
        ! result in  [sh] rad
        ilast = 1
        do i = 1, lli%im
          call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), llf(i,jf), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do

        !! copy to southpole ..
        !if ( (jf==1) .and. aver_to_prev ) then
        !  llf(:,0) = sum(llf(:,1)) / size(llf,1)
        !  aver_to_prev = .false.
        !end if

      end do ! loop over rows in fine grid

      ! integral in lat direction;
      ! weight with cos(lat) to account for smaller grid cells:
      do i = 1, lli%im
        ilast = 1
        !call IntervalQuad_Cos_Lin( lat, llf(i,:,l), lli%blat(j-1), lli%blat(j), ll(i,j,l), ilast )
        !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        llf_dim2 = llf(i,:)
        call IntervalQuad_Cos_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), ll(i,j), ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      end do

    end do ! loop over rows

    ! free memory
    deallocate( Pnm )
    deallocate( ff )
    deallocate( lons_fine )
    deallocate( row_fine )
    deallocate( llf )
    deallocate( llf_dim2 )
    deallocate( lat )
    call Done( sh )

  end subroutine IntArea_shi_ll_f


  ! ***


  subroutine IntArea_sh_ll_f( key, F, lli, ll, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin, IntervalQuad_Cos_Lin

    use Grid_Type_sh, only : TshGrid, Init, Done, SpN, sh_Pnm, Eval_Lons, Set, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGrid), intent(in)          ::  F(:)
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(out)                  ::  ll(lli%im,lli%jm,size(F))
    integer, intent(out)               ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/IntArea_sh_ll_f'


    ! --- local ------------------------------------------

    integer               ::  i, j, jf
    integer               ::  l, lm
    integer               ::  ilast

    integer               ::  T
    type(TshGrid)         ::  sh
    real, allocatable     ::  Pnm(:)

    integer               ::  refinement_i, refinement_j

    integer               ::  nlon_fine_360, nlon_fine
    real, allocatable     ::  ff(:)
    real, allocatable     ::  lons_fine(:), row_fine(:)
    real, allocatable     ::  llf(:,:,:), lat(:)

    real, allocatable     ::  llf_dim2(:)

    ! --- begin ------------------------------------------

    ! number of levels:
    lm = size(F)

    ! check arguments:
    T = F(1)%T
    do l = 2, lm
      call Check( F(l), T )
    end do

    ! use truncation up to grid resolution:
    T = ShTruncation( T, lli%dlon_deg, lli%dlat_deg )
    call Init( sh, T )

    ! allocate space for legendre coeff:
    allocate( Pnm(SpN(T)) )

    ! determine refinement (5 points per sh resolution)
    refinement_i = ShRefinement( T, lli%dlon_deg )
    refinement_j = ShRefinement( T, lli%dlat_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg * refinement_i
    nlon_fine = lli%im * refinement_i

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )

    ! lons on arc  westb+[0,..)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine_360
    end do
    lons_fine = lli%blon(0) + lons_fine

    ! ll in lon, fine in lat
    allocate( llf(lli%im,0:refinement_j,lm) )
    allocate( llf_dim2(0:refinement_j) )
    allocate( lat(0:refinement_j) )

    ! loop over latitudes in ll grid:
    !aver_to_prev = .false.
    do j = 1, lli%jm

      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        !! southpole ?
        !if ( abs(lat(jf) - (-pi/2)) < 1.0e-4 ) then
        !  ! fill with average of next row
        !  aver_to_prev = .true.
        !  cycle
        !end if
        !
        !! northpole ?
        !if ( abs(lat(jf) - (pi/2)) < 1.0e-4 ) then
        !  ! fill with average of previous row:
        !  do l = 1, lm
        !    llf(:,jf,l) = sum(llf(:,jf-1,l)) / size(llf,1)
        !  end do
        !  exit
        !end if

        ! evaluate Legendre functions:
        call sh_Pnm( Pnm, T, lat(jf), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! evaluate rows:
        do l = 1, lm

          call Set( sh, T, F(l) )
          call Eval_Lons( ff, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

          ! combine fields:
          select case ( key )

            !
            !   int int   F(k)  dA    /   A
            !
            case ( 'aver' )

              row_fine = ff / lli%area(j)

            !
            !   int int   exp(F(k))  dA    /   A
            !
            case ( 'exp,aver' )

              row_fine = exp(ff) / lli%area(j)

            !
            ! error ...
            !
            case default
              write (*,'("ERROR - unsupported integrand key : ",a)') trim(key)
              write (*,'("ERROR in ",a)') rname; status=1; return
          end select

          ! integral in lon direction assuming linear interpolation,
          ! result in  [sh] rad
          ilast = 1
          do i = 1, lli%im
            call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), llf(i,jf,l), ilast, status )
            if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          end do

        end do  ! loop over levels

        !! copy to southpole ..
        !if ( (jf==1) .and. aver_to_prev ) then
        !  do l = 1, lm
        !    llf(:,0,l) = sum(llf(:,1,l)) / size(llf,1)
        !  end do
        !  aver_to_prev = .false.
        !end if

      end do ! loop over rows in fine grid

      ! integral in lat direction;
      ! weight with cos(lat) to account for smaller grid cells:
      do l = 1, lm
        do i = 1, lli%im
          ilast = 1
          !call IntervalQuad_Cos_Lin( lat, llf(i,:,l), lli%blat(j-1), lli%blat(j), ll(i,j,l), ilast )
          !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          llf_dim2 = llf(i,:,l)
          call IntervalQuad_Cos_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), ll(i,j,l), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        end do
      end do

    end do ! loop over rows

    ! free memory
    deallocate( Pnm )
    deallocate( ff )
    deallocate( lons_fine )
    deallocate( row_fine )
    deallocate( llf )
    deallocate( llf_dim2 )
    deallocate( lat )
    call Done( sh )

    ! ok
    status = 0

  end subroutine IntArea_sh_ll_f


  !
  ! key == 'F*(da+db*exp(H))*cos'
  !
  !   int int   F(k)    (da+db*exp(H)) cos(lat)  dA
  !
  ! key == 'F*G*(db*exp(H))/cos'
  !
  !   int int   F(k)  G     db*exp(H) / cos(lat)  dA
  !
  ! (for these keys, results are added to ll)
  !
  !

  subroutine IntArea_sh_ll_fgh( key, F, G, H, da, db, lli, ll, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin

    use Grid_Type_sh, only : TshGrid, Init, Done, SpN, sh_Pnm, Eval_Lons, Set, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGrid), intent(in)          ::  F(:)
    type(TshGrid), intent(in)          ::  G, H
    real, intent(in)                   ::  da(size(F)), db(size(F))
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(inout)                ::  ll(lli%im,lli%jm,size(F))
    integer, intent(out)               ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntArea_sh_ll_fgh'

    ! --- local ------------------------------------------

    integer               ::  i, j, jf
    integer               ::  l, lm
    integer               ::  ilast

    integer               ::  T
    type(TshGrid)         ::  sh
    real, allocatable     ::  Pnm(:)

    integer               ::  refinement_i, refinement_j

    integer               ::  nlon_fine_360, nlon_fine
    real, allocatable     ::  ff(:)
    real, allocatable     ::  gg(:)
    real, allocatable     ::  exp_hh(:)
    real, allocatable     ::  lons_fine(:), row_fine(:)
    real, allocatable     ::  llf(:,:,:), lat(:)

    real, allocatable     ::  llf_dim2(:)

    logical               ::  aver_to_prev
    real                  ::  res

    ! --- begin ------------------------------------------

    ! number of levels:
    lm = size(F)

    ! check arguments:
    call Check( H )
    T = H%T
    do l = 1, lm
      call Check( F(l), T )
    end do
    call Check( G, T )

    ! use truncation up to grid resolution:
    T = ShTruncation( T, lli%dlon_deg, lli%dlat_deg )
    call Init( sh, T )

    ! allocate space for legendre coeff:
    allocate( Pnm(SpN(T)) )

    ! determine refinement (5 points per sh resolution)
    refinement_i = ShRefinement( T, lli%dlon_deg )
    refinement_j = ShRefinement( T, lli%dlat_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg * refinement_i
    nlon_fine = lli%im * refinement_i

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )
    allocate( gg(0:nlon_fine) )
    allocate( exp_hh(0:nlon_fine) )

    ! lons on arc  westb+[0,..)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine_360
    end do
    lons_fine = lli%blon(0) + lons_fine

    ! ll in lon, fine in lat
    allocate( llf(lli%im,0:refinement_j,lm) )
    allocate( lat(0:refinement_j) )
    allocate( llf_dim2(0:refinement_j) )


    ! *** integrals in lon direction

    ! loop over latitudes in ll grid:
    aver_to_prev = .false.
    do j = 1, lli%jm

      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        ! southpole ?
        if ( abs(lat(jf) - (-pi/2)) < 1.0e-4 ) then
          ! fill with average of next row
          aver_to_prev = .true.
          cycle
        end if

        ! northpole ?
        if ( abs(lat(jf) - (pi/2)) < 1.0e-4 ) then
          ! fill with average of previous row:
          do l = 1, lm
            llf(:,jf,l) = sum(llf(:,jf-1,l)) / size(llf,1)
          end do
          exit
        end if

        ! evaluate Legendre functions:
        call sh_Pnm( Pnm, T, lat(jf), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! evaluate rows:

        call Set( sh, T, G )
        call Eval_Lons( gg, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        call Set( sh, T, H )
        call Eval_Lons( exp_hh, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        exp_hh = exp(exp_hh)

        do l = 1, lm

          ! combine fields:
          select case ( key )

            !
            !   int int   F(k)    (da+db*exp(H)) cos(lat)  dA
            !
            case ( 'F*(da+db*exp(H))*cos' )

              call Set( sh, T, F(l) )
              call Eval_Lons( ff, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
              if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

              row_fine = ff * ( da(l) + exp_hh*db(l) ) * cos(lat(jf))

            !
            !   int int   F(k)  G     db*exp(H) / cos(lat)  dA
            !
            case ( 'F*G*(db*exp(H))/cos' )

              if ( db(l) == 0.0 ) then
                row_fine = 0.0
              else
                call Set( sh, T, F(l) )
                call Eval_Lons( ff, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
                if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                row_fine = ff * gg * ( exp_hh*db(l) ) / cos(lat(jf))
              end if

            !
            ! error ...
            !
            case default
              print *, 'IntArea_sh_ll_fgh - unknown key "'//trim(key)//'"'
              stop
          end select

          ! integral in lon direction assuming linear interpolation,
          ! result in  [sh] rad
          ilast = 1
          do i = 1, lli%im
            call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), llf(i,jf,l), ilast, status )
            if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          end do

        end do  ! loop over levels

        ! copy to southpole ..
        if ( (jf==1) .and. aver_to_prev ) then
          do l = 1, lm
            llf(:,0,l) = sum(llf(:,1,l)) / size(llf,1)
          end do
          aver_to_prev = .false.
        end if

      end do ! loop over rows in fine grid


      ! *** integral in lat direction

      ! 3D field:
      do l = 1, lm
        ilast = 1
        do i = 1, lli%im
          !call lquad( lat, llf(i,:,l), lli%blat(j-1), lli%blat(j), res, ilast )
          !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          llf_dim2 = llf(i,:,l)
          call IntervalQuad_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), res, ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ! add result
          ll(i,j,l) = ll(i,j,l) + res
        end do
      end do

    end do ! loop over rows

    ! free memory
    deallocate( Pnm )
    deallocate( ff )
    deallocate( gg )
    deallocate( exp_hh )
    deallocate( lons_fine )
    deallocate( row_fine )
    deallocate( llf )
    deallocate( llf_dim2 )
    deallocate( lat )
    call Done( sh )

    ! ok
    status = 0

  end subroutine IntArea_sh_ll_fgh


  ! ***


  subroutine IntArea_shi_ll_fgh( key, shi_in, nlev, F, G, H, da, db, lli, ll, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin

    use Grid_Type_sh, only : TshGridInfo, Init, Done, sh_Pnm, Eval_Lons, Set, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGridInfo), intent(in)      ::  shi_in
    integer, intent(in)                ::  nlev
    complex, intent(in)                ::  F(shi_in%np,nlev)
    complex, intent(in)                ::  G(shi_in%np)
    complex, intent(in)                ::  H(shi_in%np)
    real, intent(in)                   ::  da(nlev), db(nlev)
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(inout)                ::  ll(lli%im,lli%jm,nlev)
    integer, intent(out)               ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntArea_sh_ll_fgh'

    ! --- local ------------------------------------------

    integer               ::  i, j, jf
    integer               ::  l
    integer               ::  ilast

    type(TshGridInfo)     ::  shi
    integer               ::  T
    real, pointer         ::  Pnm(:)

    integer               ::  refinement_i, refinement_j

    integer               ::  nlon_fine_360, nlon_fine
    real, pointer         ::  ff(:)
    real, pointer         ::  gg(:)
    real, pointer         ::  exp_hh(:)
    real, pointer         ::  lons_fine(:), row_fine(:)
    real, pointer         ::  llf(:,:,:), lat(:)

    real, pointer         ::  llf_dim2(:)

    logical               ::  aver_to_prev
    real                  ::  res

    ! --- begin ------------------------------------------

    ! use truncation up to grid resolution:
    T = ShTruncation( shi_in%T, lli%dlon_deg, lli%dlat_deg )
!T = ShTruncation( 159, lli%dlon_deg, lli%dlat_deg )
!print *, 'aaa1 ', shi_in%T, shi_in%np, T

    ! spectral info:
    call Init( shi, T, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
!print *, 'aaa2 ', shi%T, shi%np

    ! allocate space for legendre coeff:
    allocate( Pnm(shi%np) )
!print *, 'aaa3 '

    ! determine refinement (5 points per sh resolution)
    refinement_i = ShRefinement( shi%T, lli%dlon_deg )
    refinement_j = ShRefinement( shi%T, lli%dlat_deg )
!print *, 'aaa4 ', refinement_i, refinement_j

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg * refinement_i
    nlon_fine = lli%im * refinement_i
!print *, 'aaa5 ', nlon_fine_360, nlon_fine

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )
    allocate( gg(0:nlon_fine) )
    allocate( exp_hh(0:nlon_fine) )
!print *, 'aaa6 '

    ! lons on arc  westb+[0,..)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine_360
    end do
    lons_fine = lli%blon(0) + lons_fine
!print *, 'aaa7 '

    ! ll in lon, fine in lat
    allocate( llf(lli%im,0:refinement_j,nlev) )
    allocate( lat(0:refinement_j) )
    allocate( llf_dim2(0:refinement_j) )
!print *, 'aaa8 '


    ! *** integrals in lon direction

    ! loop over latitudes in ll grid:
    aver_to_prev = .false.
    do j = 1, lli%jm
!print *, 'aaa9 ', j, lli%jm

      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j
!print *, 'aaa10 ', jf, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        ! southpole ?
        if ( abs(lat(jf) - (-pi/2)) < 1.0e-4 ) then
          ! fill with average of next row
          aver_to_prev = .true.
          cycle
        end if

        ! northpole ?
        if ( abs(lat(jf) - (pi/2)) < 1.0e-4 ) then
          ! fill with average of previous row:
          do l = 1, nlev
            llf(:,jf,l) = sum(llf(:,jf-1,l)) / size(llf,1)
          end do
          exit
        end if

!print *, 'aaa10a ', shi%T, shi%np, lat(jf)
        ! evaluate Legendre functions:
        call sh_Pnm( Pnm, shi%T, lat(jf), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! evaluate rows:

!print *, 'aaa10b shi_in ', shi_in%T, shi_in%np
!print *, 'aaa10b G size ', size(G)
!print *, 'aaa10b G val ', G(1:10)
!print *, 'aaa10b shi ', shi%T, shi%np
!print *, 'aaa10b Pnm size ', size(Pnm)
!print *, 'aaa10b Pnm val ', Pnm(1:10)
        call Eval_Lons( shi_in, G, shi, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, gg, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
!print *, 'aaa10c '

        call Eval_Lons( shi_in, H, shi, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, exp_hh, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
!print *, 'aaa10d '
        exp_hh = exp(exp_hh)
!print *, 'aaa10e '

        ! loop over levels:
        !xOMP PARALLEL &
        !xOMP       default ( none ) &
        !xOMP       shared  ( nlev, key, T, F, Pnm, nlon_fine_360, nlon_fine, lons_fine ) &
        !xOMP       shared  ( ff, gg, da, exp_hh, db, lat, jf, lli, llf ) &
        !xOMP       private ( l, sh, row_fine, ilast, status )
        !xOMP   DO
        do l = 1, nlev

          ! combine fields:
          select case ( key )

            !
            !   int int   F(k)    (da+db*exp(H)) cos(lat)  dA
            !
            case ( 'F*(da+db*exp(H))*cos' )

!print *, 'aaa11 ', l
              call Eval_Lons( shi_in, F(:,l), shi, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, ff, status )
              !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
!print *, 'aaa12 '

              row_fine = ff * ( da(l) + exp_hh*db(l) ) * cos(lat(jf))

            !
            !   int int   F(k)  G     db*exp(H) / cos(lat)  dA
            !
            case ( 'F*G*(db*exp(H))/cos' )

!print *, 'aaa13 ', l
              if ( db(l) == 0.0 ) then
                row_fine = 0.0
              else
                call Eval_Lons( shi_in, F(:,l), shi, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, ff, status )
                !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
                row_fine = ff * gg * ( exp_hh*db(l) ) / cos(lat(jf))
              end if
!print *, 'aaa14 '

            !
            ! error ...
            !
            case default
              print *, 'IntArea_sh_ll_fgh - unknown key "'//trim(key)//'"'
              !stop
          end select

          ! integral in lon direction assuming linear interpolation,
          ! result in  [sh] rad
          ilast = 1
          do i = 1, lli%im
            call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), llf(i,jf,l), ilast, status )
            !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          end do
!print *, 'aaa15 '

        end do  ! levels
        !xOMP   END DO
        !xOMP END PARALLEL

        ! copy to southpole ..
        if ( (jf==1) .and. aver_to_prev ) then
          do l = 1, nlev
            llf(:,0,l) = sum(llf(:,1,l)) / size(llf,1)
          end do
          aver_to_prev = .false.
        end if

      end do ! loop over rows in fine grid
!print *, 'aaa16 '


      ! *** integral in lat direction

      ! loop over levels:
      !xOMP PARALLEL DO
      do l = 1, nlev
        ilast = 1
        do i = 1, lli%im
          !call lquad( lat, llf(i,:,l), lli%blat(j-1), lli%blat(j), res, ilast )
          !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          llf_dim2 = llf(i,:,l)
          call IntervalQuad_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), res, ilast, status )
          !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ! add result
          ll(i,j,l) = ll(i,j,l) + res
        end do
      end do
      !xOMP END PARALLEL DO
!print *, 'aaa17 '

    end do ! loop over rows

    ! free memory
    deallocate( Pnm )
    deallocate( ff )
    deallocate( gg )
    deallocate( exp_hh )
    deallocate( lons_fine )
    deallocate( row_fine )
    deallocate( llf )
    deallocate( llf_dim2 )
    deallocate( lat )
!print *, 'aaa18 '

    call Done( shi )!, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
!print *, 'aaa19 '

    ! ok
    status = 0

  end subroutine IntArea_shi_ll_fgh


  ! ***


  !
  ! key == '[F*(da+db*exp(H))*cos]/[()*cos]'
  !
  !   int int F(k) (da+db*exp(H)) cos(lat) dA  /  int int (da+db*exp(H)) cos(lat) dA
  !
  ! Result is put in ll, not added.
  ! Uses cos integration in lat direction.
  !

  subroutine IntArea_sh_ll_fh( key, F, H, da, db, lli, ll, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin, IntervalQuad_Cos_Lin

    use Grid_Type_sh, only : TshGrid, Init, Done, SpN, sh_Pnm, Eval_Lons, Set, Check
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGrid), intent(in)          ::  F(:)
    type(TshGrid), intent(in)          ::  H
    real, intent(in)                   ::  da(size(F)), db(size(F))
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(inout)                ::  ll(lli%im,lli%jm,size(F))
    integer, intent(out)               ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntArea_sh_ll_fh'

    ! --- local ------------------------------------------

    integer               ::  i, j, jf
    integer               ::  l, lm
    integer               ::  ilast

    integer               ::  T
    type(TshGrid)         ::  sh
    real, allocatable     ::  Pnm(:)

    integer               ::  refinement_i, refinement_j

    integer               ::  nlon_fine_360, nlon_fine
    real, allocatable     ::  ff(:)
    real, allocatable     ::  exp_hh(:)
    real, allocatable     ::  lons_fine(:), row_fine(:)
    real, allocatable     ::  llf(:,:,:), lat(:)
    real, allocatable     ::  llf_dim2(:)

    real                  ::  res

    real, allocatable     ::  llf_expH(:,:)

    ! --- begin ------------------------------------------

    ! number of levels:
    lm = size(F)

    ! check arguments:
    call Check( H )
    T = H%T
    do l = 1, lm
      call Check( F(L), T )
    end do

    T = ShTruncation( T, lli%dlon_deg, lli%dlat_deg )
    call Init( sh, T )

    ! allocate space for legendre coeff:
    allocate( Pnm(SpN(T)) )

    ! determine refinement (5 points per sh resolution)
    refinement_i = ShRefinement( T, lli%dlon_deg )
    refinement_j = ShRefinement( T, lli%dlat_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg * refinement_i
    nlon_fine = lli%im * refinement_i

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )
    allocate( exp_hh(0:nlon_fine) )

    ! lons on arc  westb+[0,..)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine_360
    end do
    lons_fine = lli%blon(0) + lons_fine

    ! ll in lon, fine in lat
    allocate( llf(lli%im,0:refinement_j,lm) )
    allocate( llf_expH(lli%im,0:refinement_j) )
    allocate( lat(0:refinement_j) )
    allocate( llf_dim2(0:refinement_j) )


    ! *** integrals in lon direction

    ! loop over latitudes in ll grid:
    do j = 1, lli%jm
      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        ! evaluate Legendre functions:
        call sh_Pnm( Pnm, T, lat(jf), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! evaluate rows:
        call Set( sh, T, H )
        call Eval_Lons( exp_hh, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! exponent:
        exp_hh = exp(exp_hh)

        ! compute lon integral over exp(H), as a part of integral over second term
        ! (use linear interpolation):
        ilast = 1
        do i = 1, lli%im
          call IntervalQuad_Lin( lons_fine, exp_hh, lli%blon(i-1), lli%blon(i), llf_expH(i,jf), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do

        do l = 1, lm

          call Set( sh, T, F(l) )
          call Eval_Lons( ff, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

          ! combine fields:
          select case ( key )

            !
            !   int int F(k) (da+db*exp(H)) cos(lat) dA / int int (da+db*exp(H)) cos(lat) dA
            !
            case ( '[F*(da+db*exp(H))*cos]/[()*cos]' )

              ! lon integral of first term:
              row_fine = ff * ( da(l) + exp_hh*db(l) )

            !
            ! error ...
            !
            case default
              print *, 'IntArea_sh_ll_fh - unknown key "'//trim(key)//'"'
              stop
          end select

          ! integral in lon direction assuming linear interpolation,
          ! result in  [sh] rad
          ilast = 1
          do i = 1, lli%im
            call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), llf(i,jf,l), ilast, status )
            if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          end do

        end do  ! loop over levels

      end do ! loop over rows in fine grid


      ! *** integral in lat direction

      ! 3D field:
      do l = 1, lm
        ilast = 1
        do i = 1, lli%im
          !call IntervalQuad_Cos_Lin( lat, llf(i,:,l), lli%blat(j-1), lli%blat(j), res, ilast )
          !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          llf_dim2 = llf(i,:,l)
          call IntervalQuad_Cos_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), res, ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ! set result
          ll(i,j,l) = res
        end do
      end do

      ! weight with integral over exp_hh :
      ilast = 1
      do i = 1, lli%im
        !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !call IntervalQuad_Cos_Lin( lat, llf_expH(i,:), lli%blat(j-1), lli%blat(j), res, ilast )
        llf_dim2 = llf_expH(i,:)
        call IntervalQuad_Cos_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), res, ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ! weight with  int int (da+db*exp(H))*cos(lat) dA
        do l = 1, lm
          ll(i,j,l) = ll(i,j,l) / (da(l)*lli%area(j) + db(l)*res)
        end do
      end do

    end do ! loop over rows

    ! free memory
    deallocate( Pnm )
    deallocate( ff )
    deallocate( exp_hh )
    deallocate( lons_fine )
    deallocate( row_fine )
    deallocate( llf )
    deallocate( llf_dim2 )
    deallocate( llf_expH )
    deallocate( lat )

    call Done( sh )

    ! ok
    status = 0

  end subroutine IntArea_sh_ll_fh


  ! ***


  subroutine IntArea_shi_ll_fh( key, shi_in, nlev, F, H, da, db, lli, ll, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin, IntervalQuad_Cos_Lin

    use Grid_Type_sh, only : TshGrid, Init, Done
    use Grid_Type_sh, only : SpN, sh_Pnm, Eval_Lons, Set, Check, SpN, Truncate
    use Grid_Type_sh, only : TshGridInfo
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGridInfo), intent(in)      ::  shi_in
    integer, intent(in)                ::  nlev
    complex, intent(in)                ::  F(shi_in%np,nlev)
    complex, intent(in)                ::  H(shi_in%np)
    real, intent(in)                   ::  da(nlev), db(nlev)
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(inout)                ::  ll(lli%im,lli%jm,nlev)
    integer, intent(out)               ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntArea_shi_ll_fh'

    ! --- local ------------------------------------------

    integer               ::  i, j, jf
    integer               ::  l !, lm
    integer               ::  ilast

    type(TshGridInfo)     ::  shi
!    type(TshGrid)         ::  sh
    integer               ::  Ttr
    real, allocatable     ::  Pnm(:)

    integer               ::  refinement_i, refinement_j

    integer               ::  nlon_fine_360, nlon_fine
    real, allocatable     ::  ff(:)
    real, allocatable     ::  exp_hh(:)
    real, allocatable     ::  lons_fine(:), row_fine(:)
    real, allocatable     ::  llf(:,:,:), lat(:)
    real, allocatable     ::  llf_dim2(:)

    real                  ::  res

    real, allocatable     ::  llf_expH(:,:)

    ! --- begin ------------------------------------------

!    ! number of levels:
!    lm = size(F,2)
!
!    ! check arguments:
!    if ( (size(F,1) /= shi%np) .or. (size(H) /= shi%np) ) then
!      write (*,'("ERROR - number of complex coeff does not match with sh grid definition:")')
!      write (*,'("ERROR -   shi%np    : ",i6)') shi%np
!      write (*,'("ERROR -   size(F,1) : ",i6)') size(F,1)
!      write (*,'("ERROR -   size(H)   : ",i6)') size(H)
!      write (*,'("ERROR in ",a)') rname; status=1; return
!    end if
!
!    ! input temporary stored in old type grid:
!    call Init( sh, shi%T )

    ! use truncation up to grid resolution:
    Ttr = ShTruncation( shi%T, lli%dlon_deg, lli%dlat_deg )
!    call Truncate( sh, Ttr )
    call Init( shi, Ttr, status )

    ! allocate space for legendre coeff:
!    allocate( Pnm(SpN(Ttr)) )
    allocate( Pnm(shi%np) )

    ! determine refinement (5 points per sh resolution)
    refinement_i = ShRefinement( shi%T, lli%dlon_deg )
    refinement_j = ShRefinement( shi%T, lli%dlat_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg * refinement_i
    nlon_fine = lli%im * refinement_i

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )
    allocate( exp_hh(0:nlon_fine) )

    ! lons on arc  westb+[0,..)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine_360
    end do
    lons_fine = lli%blon(0) + lons_fine

    ! ll in lon, fine in lat
    !allocate( llf(lli%im,0:refinement_j,lm) )
    allocate( llf(lli%im,0:refinement_j,nlev) )
    allocate( llf_expH(lli%im,0:refinement_j) )
    allocate( lat(0:refinement_j) )
    allocate( llf_dim2(0:refinement_j) )


    ! *** integrals in lon direction

    ! loop over latitudes in ll grid:
    do j = 1, lli%jm
      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        ! evaluate Legendre functions:
        !call sh_Pnm( Pnm, Ttr, lat(jf), status )
        call sh_Pnm( Pnm, shi%T, lat(jf), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! evaluate rows:

        !call Set( sh, shi%T, H )
        !call Truncate( sh, Ttr )

        !call Eval_Lons( exp_hh, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
        !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        call Eval_Lons( shi_in, H, shi, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, exp_hh, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        exp_hh = exp(exp_hh)

        ! compute lon integral over exp(H), as a part of integral over second term
        ! (use linear interpolation):
        ilast = 1
        do i = 1, lli%im
          call IntervalQuad_Lin( lons_fine, exp_hh, lli%blon(i-1), lli%blon(i), llf_expH(i,jf), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do

        ! combine fields:
        select case ( key )

          !
          !   int int F(k) (da+db*exp(H)) cos(lat) dA / int int (da+db*exp(H)) cos(lat) dA
          !
          case ( '[F*(da+db*exp(H))*cos]/[()*cos]' )

            !$OMP PARALLEL &
            !$OMP       default ( none ) &
            !$OMP       shared  ( shi_in, F, shi, Pnm, nlon_fine_360, lons_fine, nlon_fine, lli, llf ) &
            !$OMP       shared  ( exp_hh, da, db, jf ) &
            !$OMP       private ( l, nlev, ff, status, row_fine, ilast )
            !$OMP   DO
            do l = 1, nlev
            !do l = 1, lm

              !call Set( sh, shi%T, F(:,l) )
              !call Truncate( sh, Ttr )

              call Eval_Lons( shi_in, F(:,l), shi, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, ff, status )
              !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

              ! lon integral of first term:
              row_fine = ff * ( da(l) + exp_hh*db(l) )

              ! integral in lon direction assuming linear interpolation,
              ! result in  [sh] rad
              ilast = 1
              do i = 1, lli%im
                call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), llf(i,jf,l), ilast, status )
                !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
              end do

            end do  ! loop over levels
            !$OMP   END DO
            !$OMP END PARALLEL

          !
          ! error ...
          !
          case default
            print *, 'IntArea_shi_ll_fh - unknown key "'//trim(key)//'"'
            stop
        end select

      end do ! loop over rows in fine grid


      ! *** integral in lat direction

      ! 3D field:
      do l = 1, nlev
      !do l = 1, lm
        ilast = 1
        do i = 1, lli%im
          !call IntervalQuad_Cos_Lin( lat, llf(i,:,l), lli%blat(j-1), lli%blat(j), res, ilast )
          !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          llf_dim2 = llf(i,:,l)
          call IntervalQuad_Cos_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), res, ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ! set result
          ll(i,j,l) = res
        end do
      end do

      ! weight with integral over exp_hh :
      ilast = 1
      do i = 1, lli%im
        !>>> bsf15k bug fix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !call IntervalQuad_Cos_Lin( lat, llf_expH(i,:), lli%blat(j-1), lli%blat(j), res, ilast )
        llf_dim2 = llf_expH(i,:)
        call IntervalQuad_Cos_Lin( lat, llf_dim2, lli%blat(j-1), lli%blat(j), res, ilast, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ! weight with  int int (da+db*exp(H))*cos(lat) dA
        do l = 1, nlev
        !do l = 1, lm
          ll(i,j,l) = ll(i,j,l) / (da(l)*lli%area(j) + db(l)*res)
        end do
      end do

    end do ! loop over rows

    ! free memory
    deallocate( Pnm )
    deallocate( ff )
    deallocate( exp_hh )
    deallocate( lons_fine )
    deallocate( row_fine )
    deallocate( llf )
    deallocate( llf_dim2 )
    deallocate( llf_expH )
    deallocate( lat )

    !call Done( sh )

    ! ok
    status = 0

  end subroutine IntArea_shi_ll_fh


  ! ***


  subroutine IntLat_sh_ll( key, F, H, da, db, lli, ll_u, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin

    use Grid_Type_sh, only : TshGrid, Init, Done, Set, Check
    use Grid_Type_sh, only : sh_Pnm, Eval_Lons, SpN
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGrid), intent(in)          ::  F(:)
    type(TshGrid), intent(in)          ::  H
    real, intent(in)                   ::  da(size(F))
    real, intent(in)                   ::  db(size(F))
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(inout)                ::  ll_u(0:lli%im,lli%jm,size(F))
    integer, intent(out)               ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntLat_sh_ll'

    ! --- local ------------------------------------------

    integer               ::  i, j, jf
    integer               ::  l, lm
    integer               ::  ilast

    integer               ::  T
    type(TshGrid)         ::  sh
    real, allocatable     ::  Pnm(:)

    integer               ::  refinement_j

    integer               ::  nlon_fine_360, nlon_fine
    real, allocatable     ::  ff(:)
    real, allocatable     ::  exp_hh(:)
    real, allocatable     ::  llf(:,:,:), lat(:), llf_row(:)

    logical               ::  aver_to_prev

    ! --- begin ------------------------------------------

    ! number of levels:
    lm = size(F)

    ! check arguments:
    call Check( H )
    T = H%T
    do l = 1, lm
      call Check( F(l), T )
    end do

    T = ShTruncation( T, lli%dlon_deg, lli%dlat_deg )
    call Init( sh, T )

    ! allocate space for legendre coeff:
    allocate( Pnm(SpN(T)) )

    ! determine refinement based on spectral resolution
    ! and length of required integration interval:
    refinement_j = ShRefinement( T, lli%dlat_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg
    nlon_fine = lli%im

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )
    allocate( exp_hh(0:nlon_fine) )

    ! ll in lon, fine in lat
    allocate( llf(0:lli%im,0:refinement_j,lm) )
    allocate( llf_row(0:refinement_j) )
    allocate( lat(0:refinement_j) )

    ! loop over latitudes in ll grid:
    aver_to_prev = .false.
    do j = 1, lli%jm
      ! loop over latitudes in fine grid:
      do jf = 0, refinement_j

        ! latitude in fine grid:
        lat(jf) = lli%blat(j-1) + jf*lli%dlat/refinement_j

        ! southpole ?
        if ( abs(lat(jf) - (-pi/2)) < 1.0e-4 ) then
          ! fill with average of next row
          aver_to_prev = .true.
          cycle
        end if

        ! northpole ?
        if ( abs(lat(jf) - (pi/2)) < 1.0e-4 ) then
          ! fill with average of previous row:
          do l = 1, lm
            llf(:,jf,l) = sum(llf(:,jf-1,l)) / size(llf,1)
          end do
          exit
        end if

        ! evaluate Legendre functions:
        call sh_Pnm( Pnm, T, lat(jf), status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! evaluate rows:

        call Set( sh, T, H )
        call Eval_Lons( exp_hh, sh, Pnm, nlon_fine_360, lli%blon(0), nlon_fine+1, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        exp_hh = exp(exp_hh)

        do l = 1, lm
          call Set( sh, T, F(l) )
          call Eval_Lons( ff, sh, Pnm, nlon_fine_360, lli%blon(0), nlon_fine+1, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

          ! combine fields:
          select case ( key )

            !
            !  [ int exp(H) dlat ] / deltalat
            !
            case ( 'exp(H),aver' )

              llf(:,jf,l) = exp_hh / lli%dlat


            !
            !  int  F (da+db*exp(H)) / cos(lat)  dlat
            !
            case ( '(da+exp*db)/cos' )

              llf(:,jf,l) = ff * ( da(l) + exp_hh*db(l) ) / cos(lat(jf))

            !
            ! error ...
            !
            case default
              print *, 'IntLat_sh_ll - unknown key "'//trim(key)//'"'
              stop
          end select

        end do  ! loop over levels

        ! copy to southpole ..
        if ( (jf==1) .and. aver_to_prev ) then
          do l = 1, lm
            llf(:,0,l) = sum(llf(:,1,l)) / size(llf,1)
          end do
          aver_to_prev = .false.
        end if

      end do ! loop over rows in fine grid

      ! integral in lat direction:
      do l = 1, lm
        do i = 0, lli%im
          ilast = 1
          llf_row = llf(i,:,l)
          call IntervalQuad_Lin( lat, llf_row, lli%blat(j-1), lli%blat(j), ll_u(i,j,l), ilast, status )
          if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do
      end do

    end do ! loop over rows

    ! free memory
    call Done( sh )
    deallocate( Pnm )
    deallocate( ff )
    deallocate( exp_hh )
    deallocate( llf )
    deallocate( llf_row )
    deallocate( lat )

    ! ok
    status = 0

  end subroutine IntLat_sh_ll


  ! ***


  subroutine IntLon_sh_ll( key, F, H, da, db, lli, ll_v, status )

    use grid_tools, only : deg2rad, pi

    use Num, only : IntervalQuad_Lin

    use Grid_Type_sh, only : TshGrid, Init, Done, SpN, Set, Check
    use Grid_Type_sh, only : sh_Pnm, Eval_Lons
    use Grid_Type_ll, only : TllGridInfo, Check

    ! --- in/out -----------------------------------------

    character(len=*), intent(in)       ::  key
    type(TshGrid), intent(in)          ::  F(:)
    type(TshGrid), intent(in)          ::  H
    real, intent(in)                   ::  da(size(F))
    real, intent(in)                   ::  db(size(F))
    type(TllGridInfo), intent(in)      ::  lli
    real, intent(inout)                ::  ll_v(lli%im,0:lli%jm,size(F))
    integer, intent(out)               ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntLon_sh_ll'

    ! --- local ------------------------------------------

    integer               ::  i, j
    integer               ::  l, lm
    integer               ::  ilast

    integer               ::  T
    type(TshGrid)         ::  sh
    real, allocatable     ::  Pnm(:)

    integer               ::  refinement_i

    integer               ::  nlon_fine_360, nlon_fine
    real, allocatable     ::  ff(:)
    real, allocatable     ::  exp_hh(:)
    real, allocatable     ::  lons_fine(:), row_fine(:)

    logical               ::  pole

    ! --- begin ------------------------------------------

    ! number of levels:
    lm = size(F)

    ! check arguments:
    call Check( H )
    T = H%T
    do l = 1, lm
      call Check( F(l), T )
    end do

    ! use truncated harmonics ?
    T = ShTruncation( T, lli%dlon_deg, lli%dlat_deg )
    call Init( sh, T )

    ! allocate space for legendre coeff:
    allocate( Pnm(SpN(T)) )

    ! determine refinement based on spectral resolution
    ! and length of required integration interval:
    refinement_i = ShRefinement( T, lli%dlon_deg )

    ! number of lons in fine grid on complete circle:
    nlon_fine_360 = 360.0/lli%dlon_deg * refinement_i
    nlon_fine = lli%im * refinement_i

    ! store evaluation of spectral field:
    allocate( ff(0:nlon_fine) )
    allocate( exp_hh(0:nlon_fine) )

    ! lons on arc  westb+[0,..)
    allocate( row_fine(0:nlon_fine) )
    allocate( lons_fine(0:nlon_fine) )
    do i = 0, nlon_fine
      lons_fine(i) = i*2*pi/nlon_fine_360
    end do
    lons_fine = lli%blon(0) + lons_fine

    ! loop over boundary latitudes in ll grid:
    do j = 0, lli%jm

      ! pole ? then  int f(x) dx = f
      pole = ( abs(lli%blat(j) - (-pi/2)) < 1.0e-4 ) .or. &
             ( abs(lli%blat(j) - ( pi/2)) < 1.0e-4 )

      ! evaluate Legendre functions:
      call sh_Pnm( Pnm, T, lli%blat(j), status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

      ! evaluate rows:

      call Set( sh, T, H )
      call Eval_Lons( exp_hh, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      exp_hh = exp(exp_hh)

      ! loop over levels
      do l = 1, lm
        call Set( sh, T, F(l) )
        call Eval_Lons( ff, sh, Pnm, nlon_fine_360, lons_fine(0), nlon_fine+1, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

        ! combine fields:
        select case ( key )

          !
          !  [ int exp(H) dlon ] / deltalon
          !
          case ( 'exp(H),aver' )

            if ( pole ) then
              row_fine = exp_hh
            else
              row_fine = exp_hh / lli%dlon
            end if


          !
          !   int  F (da+db*exp(H)) dlon
          !
          case ( '(da+exp*db)' )

            if ( pole ) then
              row_fine = 0.0
            else
              row_fine = ff * ( da(l) + exp_hh*db(l) )
            end if

          !
          ! error ...
          !
          case default
            print *, 'IntLon_sh_ll - unknown key "'//trim(key)//'"'
            stop
        end select

        ! special treatment of poles:
        if ( pole ) then
          ! same value at al 'longitudes', thus for example the average ...
          ll_v(:,j,l) = sum(row_fine)/size(row_fine)
        else
          ! integral in lon direction assuming linear interpolation,
          ! result in  [sh] rad
          ilast = 1
          do i = 1, lli%im
            call IntervalQuad_Lin( lons_fine, row_fine, lli%blon(i-1), lli%blon(i), ll_v(i,j,l), ilast, status )
            if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
          end do
        end if

      end do  ! loop over levels

    end do ! loop over rows

    ! free memory
    call Done( sh )
    deallocate( Pnm )
    deallocate( ff )
    deallocate( exp_hh )
    deallocate( lons_fine )
    deallocate( row_fine )

    ! ok
    status = 0

  end subroutine IntLon_sh_ll


end module Grid_Interpol
