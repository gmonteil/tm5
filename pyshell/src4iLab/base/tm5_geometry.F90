!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module TM5_Geometry

  use GO  , only : gol, goPr, goErr
  use Dims, only : nregions, nregions_all
  use Dims, only : lmax_conv
  use Grid, only : TllGridInfo, TLevelInfo

  implicit none

  ! --- in/out -----------------------------------

  private

  public  ::  TM5_Geometry_Init, TM5_Geometry_Done

  public  ::  lli, lli_1x1
  public  ::  levi, levi_ec, levi_tm60
  public  ::  lmax_conv

  ! --- const ----------------------------------

  ! module name
  character(len=*), parameter  ::  mname = 'TM5_Geometry'


  ! --- var ---------------------------------------------

  ! horizontal grid definitions:
  type(TllGridInfo)            ::  lli(0:nregions_all)
  type(TllGridInfo)            ::  lli_1x1

  ! vertical grid definition:
  type(TLevelInfo)             ::  levi
  type(TLevelInfo)             ::  levi_ec
  type(TLevelInfo)             ::  levi_tm60



contains


  ! ==========================================================


  subroutine TM5_Geometry_Init( status )

    use Grid         , only : Init
    use dims,          only : nregions
    use dims,          only : im, jm
    use dims,          only : dxy11, nlat180
    use Dims         , only : region_name
    use Dims         , only : iglbsfc
    use Dims         , only : xbeg, dx, xref, im
    use Dims         , only : ybeg, dy, yref, jm
    use Dims         , only : echlevs, lme, a_ec, b_ec
    use redgridZoom,   only : grid_reduced, initredgrid
    use geometry,      only : geomtryh, calc_dxy

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TM5_Geometry_Init'

    ! --- local ----------------------------------

    integer             :: region

    ! --- begin ----------------------------------

    write (gol,'(a,":     reduced grid ...")') rname; call goPr

    ! set up model geometry...
    if ( grid_reduced ) then
      ! loop over regions:
      do region = 1, nregions
        ! setup reduced grid for this region:
        call initredgrid( region, im(region), jm(region), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    write (gol,'(a,":     horizontal geometry ...")') rname; call goPr

    ! Area of grid boxes on 1x1 degree
    call calc_dxy(dxy11,nlat180)

    do region = 1, nregions

       ! horizontal geometry
       call geomtryh(region)

    enddo

    write (gol,'(a,":     horizontal grids ...")') rname; call goPr

    ! global grid of single cell:
    call Init( lli(0), 0.0, 360.0, 1, 0.0, 180.0, 1, status, &
                 name='globe' )
    IF_NOTOK_RETURN(status=1)

    ! define mid of first west/south cells, spacing, and size:
    do region = 1, nregions_all
      ! switch ...
      if ( (1 <= region) .and. (region <= nregions) ) then
      call Init( lli(region), xbeg(region)+dx/xref(region)/2.0, real(dx)/xref(region), im(region), &
                                ybeg(region)+dy/yref(region)/2.0, real(dy)/yref(region), jm(region), status, &
                                name=trim(region_name(region)) )
      IF_NOTOK_RETURN(status=1)
      else if ( region == iglbsfc ) then
        call Init( lli(region), -179.5, 1.0, 360, -89.5, 1.0, 180, status, &
                 name='glb100x100' )
        IF_NOTOK_RETURN(status=1)
      else
        write (gol,'("could not set grid definition for region ",i0)') region; call goErr
        TRACEBACK; status=1; return
      end if
    end do

    ! old, should be replaced by lli(iglbsfc)
    call Init( lli_1x1, -179.5, 1.0, 360, -89.5, 1.0, 180, status, &
                 name='glb100x100' )
    IF_NOTOK_RETURN(status=1)

    write (gol,'(a,":     vertical layers ...")') rname; call goPr

    ! setup parent level definition:
    call Init( levi_ec, lme, a_ec, b_ec, status )          ! ecmwf levels
    IF_NOTOK_RETURN(status=1)

    ! setup level definition:
    call Init( levi, levi_ec, echlevs, status )   ! tm half level selection
    IF_NOTOK_RETURN(status=1)

    ! all levels ...
    call Init( levi_tm60, 'tm60', status )          ! TM 60 levels

    ! ok
    status = 0

  end subroutine TM5_Geometry_Init


  ! ***


  subroutine TM5_Geometry_Done( status )

    use Dims, only : nregions
    use Grid, only : Done

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'TM5_Geometry_Done'

    ! --- local ----------------------------------

    integer           ::  region

    ! --- begin ----------------------------------

    ! horizontal grids:
    do region = 1, nregions
      call Done( lli(region), status )
      IF_NOTOK_RETURN(status=1)
    end do
    !
    call Done( lli_1x1, status )
    IF_NOTOK_RETURN(status=1)

    ! done 60level definition:
    call Done( levi_tm60, status )
    IF_NOTOK_RETURN(status=1)

    ! done parent level definition:
    call Done( levi_ec, status )
    IF_NOTOK_RETURN(status=1)

    ! level definition:
    call Done( levi, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TM5_Geometry_Done



end module TM5_Geometry
