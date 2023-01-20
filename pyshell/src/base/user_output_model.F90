!###############################################################################
!
! Put out information on model:
!  o regions
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

module User_Output_Model

  use GO, only : gol, goPr, goErr
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

  implicit none


  ! --- in/out -----------------------------------

  private

  public :: User_Output_Model_Init, User_Output_Model_Done


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'User_Output_Model'


  ! --- var ------------------------------------

  ! base path:
  character(len=MAX_FILENAME_LEN) ::  model_output_dir


contains


  ! ====================================================================


  subroutine User_Output_Model_Init( status )

    use GO, only : ReadRc
    use GO, only : pathsep

    use MDF, only : MDF_Init, MDF_Done

    use global_data, only : rcF
    use global_data, only : outdir

    ! --- in/out ---------------------------------

    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/User_Output_Model_Init'

    ! --- local ----------------------------------

    character(len=MAX_RCKEY_LEN) :: subdir

    ! --- begin ----------------------------------

    ! read output subdirectory from settings:
    call ReadRc( rcF, 'model.output.subdir', subdir, status)
    IF_NOTOK_RETURN(status=1)
    ! base path:
    write (model_output_dir,'(3a)') trim(outdir), pathsep, trim(subdir)

    ! setup MDF interface to HDF/NetCDF :
    call MDF_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! write file with region defintions:
    call User_Output_Model_Regions( status )
    IF_NOTOK_RETURN(status=1)

    ! done with MDF interface:
    call MDF_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine User_Output_Model_Init


  ! ***


  subroutine User_Output_Model_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/User_Output_Model_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! nothing to be done ...

    ! ok
    status = 0

  end subroutine User_Output_Model_Done


  ! ***


  subroutine User_Output_Model_Regions( status )

    use GO       , only : pathsep
    use MDF      , only : MDF_Create, MDF_Close, MDF_EndDef
    use MDF      , only : MDF_HDF4, MDF_REPLACE, MDF_GLOBAL, MDF_CHAR, MDF_INT, MDF_FLOAT
    use MDF      , only : MDF_Put_Att
    use MDF      , only : MDF_Def_Dim
    use MDF      , only : MDF_Def_Var, MDF_Put_Var
    use dims     , only : nregions
    use dims     , only : len_region_name, region_name
    use dims     , only : xbeg, xend, im
    use dims     , only : ybeg, yend, jm
    use dims     , only : parent
    use TM5_Geometry,   only : lli
    use RedgridZoom,    only : nred, jred, clustsize

    ! --- in/out ---------------------------------

    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/User_Output_Model_Regions'

    ! --- local ----------------------------------

    character(len=MAX_FILENAME_LEN) ::  fname
    integer                   ::  hid
    integer                   ::  dimid_region, dimid_len_region_name
    integer                   ::  varid_region_name
    integer                   ::  varid_xbeg, varid_xend, varid_nx, varid_dx
    integer                   ::  varid_ybeg, varid_yend, varid_ny, varid_dy
    integer                   ::  varid_parent
    integer                   ::  dimid_lon, dimid_blon
    integer                   ::  dimid_lat, dimid_blat
    integer                   ::  varid_lon, varid_blon
    integer                   ::  varid_lat, varid_blat
    integer                   ::  varid_rg_clustsize

    integer                   ::  region
    integer                   ::  imr, jmr
    integer, allocatable      ::  rg_clustsize(:)
    integer                   ::  ired

    ! --- begin ----------------------------------

    ! * overview file

    ! compose filename:
    write (fname,'(a,a,"regions.hdf")') trim(model_output_dir), pathsep

    ! new file:
    call MDF_Create( trim(fname), MDF_HDF4, MDF_REPLACE, hid, status )
    IF_NOTOK_RETURN(status=1)

    ! define dimensions:
    call MDF_Def_Dim( hid, 'region', nregions, dimid_region, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Dim( hid, 'len_region_name', len_region_name, dimid_len_region_name, status )
    IF_NOTOK_RETURN(status=1)

    ! variables:
    call MDF_Def_Var( hid, 'region_name', MDF_CHAR, (/dimid_len_region_name,dimid_region/), varid_region_name, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'xbeg', MDF_FLOAT, (/dimid_region/), varid_xbeg, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'xend', MDF_FLOAT, (/dimid_region/), varid_xend, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'ybeg', MDF_FLOAT, (/dimid_region/), varid_ybeg, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'yend', MDF_FLOAT, (/dimid_region/), varid_yend, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'nx', MDF_INT, (/dimid_region/), varid_nx, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'ny', MDF_INT, (/dimid_region/), varid_ny, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'dx', MDF_FLOAT, (/dimid_region/), varid_dx, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'dy', MDF_FLOAT, (/dimid_region/), varid_dy, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Def_Var( hid, 'parent', MDF_INT, (/dimid_region/), varid_parent, status )
    IF_NOTOK_RETURN(status=1)

    ! finished definition:
    call MDF_EndDef( hid, status )
    IF_NOTOK_RETURN(status=1)

    ! fill:
    call MDF_Put_Var( hid, varid_region_name, region_name, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_xbeg, xbeg, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_xend, xend, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_ybeg, ybeg, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_yend, yend, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_nx, im, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_ny, jm, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_dx, (xend-xbeg)/float(im), status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_dy, (yend-ybeg)/float(jm), status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Put_Var( hid, varid_parent, parent, status )
    IF_NOTOK_RETURN(status=1)

    ! close file:
    call MDF_Close( hid, status )
    IF_NOTOK_RETURN(status=1)

    ! * region files

    ! loop over regions:
    do region = 1, nregions

      ! local dimensions:
      imr = im(region)
      jmr = jm(region)

      ! compose filename:
      write (fname,'(a,a,"region_",a,".hdf")') trim(model_output_dir), pathsep, trim(region_name(region))

      ! new file:
      call MDF_Create( trim(fname), MDF_HDF4, MDF_REPLACE, hid, status )
      IF_NOTOK_RETURN(status=1)

      ! global attributes:
      call MDF_Put_Att( hid, MDF_GLOBAL, 'region_name', trim(region_name(region)), status )
      IF_NOTOK_RETURN(status=1)
      if ( parent(region) == 0 ) then
        call MDF_Put_Att( hid, MDF_GLOBAL, 'parent', 'globe', status )
        IF_NOTOK_RETURN(status=1)
      else
        call MDF_Put_Att( hid, MDF_GLOBAL, 'parent', trim(region_name(parent(region))), status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! define dimensions:
      call MDF_Def_Dim( hid, 'lon' , imr  , dimid_lon , status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Def_Dim( hid, 'blon', imr+1, dimid_blon, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Def_Dim( hid, 'lat' , jmr  , dimid_lat , status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Def_Dim( hid, 'blat', jmr+1, dimid_blat, status )
      IF_NOTOK_RETURN(status=1)

      ! grid variables:
      call MDF_Def_Var( hid, 'lon' , MDF_FLOAT, (/dimid_lon /), varid_lon , status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Def_Var( hid, 'blon', MDF_FLOAT, (/dimid_blon/), varid_blon, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Def_Var( hid, 'lat' , MDF_FLOAT, (/dimid_lat /), varid_lat , status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Def_Var( hid, 'blat', MDF_FLOAT, (/dimid_blat/), varid_blat, status )
      IF_NOTOK_RETURN(status=1)
      ! reduced grid variables:
      call MDF_Def_Var( hid, 'rg_clustsize', MDF_INT, (/dimid_lat/), varid_rg_clustsize, status )
      IF_NOTOK_RETURN(status=1)

      ! finished definition:
      call MDF_EndDef( hid, status )
      IF_NOTOK_RETURN(status=1)

      ! write grid variables:
      call MDF_Put_Var( hid, varid_lon , lli(region)%lon_deg , status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Put_Var( hid, varid_blon, lli(region)%blon_deg, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Put_Var( hid, varid_lat , lli(region)%lat_deg , status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Put_Var( hid, varid_blat, lli(region)%blat_deg, status )
      IF_NOTOK_RETURN(status=1)

      ! write reduced grid clust size:
      allocate( rg_clustsize(jmr) )
      rg_clustsize = 1
      do ired = 1, nred(region)
        rg_clustsize(jred(ired,region)) = clustsize(ired,region)
      end do
      call MDF_Put_Var( hid, varid_rg_clustsize, rg_clustsize, status )
      IF_NOTOK_RETURN(status=1)
      deallocate( rg_clustsize )

      ! close file:
      call MDF_Close( hid, status )
      IF_NOTOK_RETURN(status=1)

    end do  ! regions

    ! ok
    status = 0

  end subroutine User_Output_Model_Regions


end module User_Output_Model
