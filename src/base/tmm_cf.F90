!###############################################################################
!
! Converstion from TM5 names/units to CF standard names/units
!
! USAGE
!
!   use TMM_CF
!
!   ! initialize module:
!   !   o read name of CF standard table from rcfile key 'cf-standard-name-table'
!   call TMM_CF_Init( rcf, status )
!     type(TRcFile), intent(in)           ::  rcf
!     integer, intent(out)                ::  status
!
!   ! get standard units for standard name:
!   call TMM_CF_Standard_Units( cf_standard_name, cf_units, status )
!     character(len=*), intent(in)                ::  cf_standard_name
!     character(len=*), intent(out)               ::  cf_units
!     integer, intent(out)                        ::  status
!
!   ! convert from TM5 variable name to CF standard name:
!   call TMM_CF_Convert_Name( tm5_name, cf_standard_name, status )
!     character(len=*), intent(in)                ::  tm5_name
!     character(len=*), intent(out)               ::  cf_standard_name
!     integer, intent(out)                        ::  status
!
!   ! get conversion factor from TM5 units to CF standard unit:
!   call subroutine TMM_CF_Convert_Units( tm5_units, cf_units, tm5_to_cf_factor, status )
!     character(len=*), intent(in)                ::  tm5_units
!     character(len=*), intent(in)                ::  cf_units
!     real, intent(out)                           ::  tm5_to_cf_factor
!     integer, intent(out)                        ::  status
!
!   ! done with module:
!   call TMM_CF_Done( status )
!     integer, intent(out)                        ::  status
!
!
! EXTERNAL LIBRARIES
!
!   Module uses the 'UDUnits' library (version 1), and the Fortran90 interface module 'UDUnits' .
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_UDUNITS_NOTOK_RETURN(action) if (status/=UDUNITS_NOERR) then; gol=trim(UDUnits_StrError(status)); call goErr; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module TMM_CF

  use GO     , only : gol, goErr, goPr
#ifdef with_udunits
  use UDUnits, only : UDUNITS_NOERR, UDUnits_StrError
#endif

  use Standard_Name_Table, only : T_Standard_Name_Table

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  TMM_CF_Init, TMM_CF_Done
  public  ::  TMM_CF_Standard_Units
  public  ::  TMM_CF_Convert_Name
  public  ::  TMM_CF_Convert_Units


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'TMM_CF'

  ! rcfile key with name of cf table:
  character(len=*), parameter  ::  rckey_cf_table = 'cf-standard-name-table'

  ! --- local ------------------------------------

  logical                               ::  with_cf_table
  type(T_Standard_Name_Table), save     ::  cf_table



contains


  ! ==============================================================


  subroutine TMM_CF_Init( rcf, status )

    use GO                 , only : TrcFile, ReadRc
    use Standard_Name_Table, only : Standard_Name_Table_Init
#ifdef with_udunits
    use UDUnits            , only : UDUnits_Init
#endif
    use os_specs           , only : WRITE_STR_LEN

    ! --- in/out ---------------------------------

    type(TRcFile), intent(in)           ::  rcf
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TMM_CF_Init'

    ! --- local ----------------------------------

    character(len=WRITE_STR_LEN)      ::  cf_standard_name_table

    ! --- begin ----------------------------------

    ! read full path of cf table:
    call ReadRc( rcF, rckey_cf_table, cf_standard_name_table, status, default='None' )
    IF_ERROR_RETURN(status=1)
    ! defined ?
    with_cf_table = trim(cf_standard_name_table) /= 'None'
    ! continue with cf table ?
    if ( with_cf_table ) then
      ! really exists ?
      inquire( file=trim(cf_standard_name_table), exist=with_cf_table )
    end if

    ! standard names and units:
    if ( with_cf_table ) then
      call Standard_Name_Table_Init( cf_table, trim(cf_standard_name_table), status )
      if ( status /= 0 ) then
        write (gol,'("opening standard name table `",a,"`")') trim(cf_standard_name_table); call goErr
        TRACEBACK; status=1; return
      end if
    end if

#ifdef with_udunits
    ! initialize UDUnits module:
    call UDUnits_Init( status )
    IF_UDUNITS_NOTOK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine TMM_CF_Init


  ! ***


  subroutine TMM_CF_Done( status )

    use Standard_Name_Table, only : Standard_Name_Table_Done
#ifdef with_udunits
    use UDUnits            , only : UDUnits_Done
#endif

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TMM_CF_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! done with standard names and units:
    if ( with_cf_table ) then
      call Standard_Name_Table_Done( cf_table, status )
      IF_NOTOK_RETURN(status=1)
    end if

#ifdef with_udunits
    ! done with UDUnits module:
    call UDUnits_Done( status )
    IF_UDUNITS_NOTOK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine TMM_CF_Done


  ! ***


  subroutine TMM_CF_Standard_Units( cf_standard_name, cf_units, status )

    ! --- in/out ---------------------------------

    character(len=*), intent(in)                ::  cf_standard_name
    character(len=*), intent(out)               ::  cf_units
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TMM_CF_Standard_Units'

    ! --- local ----------------------------------

    integer                         ::  i

    ! --- begin ----------------------------------

    ! no unit yet ...
    cf_units = ''
    ! check for specials ...
    select case ( trim(cf_standard_name) )
      ! TM5 special fields
      case ( 'tm5_integrated_eastward_mass_flux_of_air'  ) ; cf_units = 'kg s-1'
      case ( 'tm5_integrated_northward_mass_flux_of_air' ) ; cf_units = 'kg s-1'
      case ( 'tm5_integrated_upward_mass_flux_of_air'    ) ; cf_units = 'kg s-1'
      ! CF standard fields
      case default
        ! check ..
        if ( .not. with_cf_table ) then
          write (gol,'("no CF table: either rcfile key `",a,"` not defined,",&
                  &" or specified table not found")') trim(rckey_cf_table); call goErr
          TRACEBACK; status=1; return
        end if
        ! loop over all entries:
        do i = 1, size(cf_table%entry)
          ! compare ...
          if ( trim(cf_table%entry(i)%id) == trim(cf_standard_name) ) then
            ! copy values:
            cf_units = trim(cf_table%entry(i)%canonical_units)
            ! found!
            exit
          end if
        end do   ! CF table entries
    end select
    ! not found ?
    if ( len_trim(cf_units) == 0 ) then
      write (gol,'("id not found in cf standard name table : ",a)') trim(cf_standard_name); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine TMM_CF_Standard_Units


  ! ***


  subroutine TMM_CF_Convert_Name( tm5_name, cf_standard_name, status )

    ! --- in/out ---------------------------------

    character(len=*), intent(in)                ::  tm5_name
    character(len=*), intent(out)               ::  cf_standard_name
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TMM_CF_Convert_Name'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! convert from TM5 internal name to CF standard name:
    select case ( trim(tm5_name) )
      case ( 'lon'        ) ; cf_standard_name = 'longitude'
      case ( 'lat'        ) ; cf_standard_name = 'latitude'
      case ( 'sp'         ) ; cf_standard_name = 'surface_air_pressure'
      case ( 'mfu'        ) ; cf_standard_name = 'tm5_integrated_eastward_mass_flux_of_air'
      case ( 'mfv'        ) ; cf_standard_name = 'tm5_integrated_northward_mass_flux_of_air'
      case ( 'mfw'        ) ; cf_standard_name = 'tm5_integrated_upward_mass_flux_of_air'
      case ( 'tsp'        ) ; cf_standard_name = 'tendency_of_surface_air_pressure'
      case ( 'oro'        ) ; cf_standard_name = 'geopotential'
      case ( 'T'          ) ; cf_standard_name = 'air_temperature'
      case ( 'Q'          ) ; cf_standard_name = 'specific_humidity'
      case ( 'CLWC'       ) ; cf_standard_name = 'mass_fraction_of_cloud_liquid_water_in_air'
      case ( 'CIWC'       ) ; cf_standard_name = 'mass_fraction_of_cloud_ice_in_air'
      case ( 'CC'         ) ; cf_standard_name = 'cloud_area_fraction_in_atmosphere_layer'
      case ( 'CCO'        ) ; cf_standard_name = 'cloud_area_fraction_in_atmosphere_layer'       ! dummy for overhead cloud cover
      case ( 'CCU'        ) ; cf_standard_name = 'cloud_area_fraction_in_atmosphere_layer'       ! dummy for underfeet cloud cover
      case ( 'eu'         ) ; cf_standard_name = 'atmosphere_updraft_convective_mass_flux'       ! dummy for entrainement updraft
      case ( 'ed'         ) ; cf_standard_name = 'atmosphere_downdraft_convective_mass_flux'     ! dummy for entrainement downdraft
      case ( 'du'         ) ; cf_standard_name = 'atmosphere_updraft_convective_mass_flux'       ! dummy for detrainement updraft
      case ( 'dd'         ) ; cf_standard_name = 'atmosphere_downdraft_convective_mass_flux'     ! dummy for detrainement downdraft
      case ( 'cloud_base' ) ; cf_standard_name = 'model_level_number_at_convective_cloud_base'
      case ( 'cloud_top'  ) ; cf_standard_name = 'model_level_number_at_convective_cloud_top'
      case ( 'cloud_lfs'  ) ; cf_standard_name = 'model_level_number_at_convective_cloud_top'    ! dummy for level-of-free-sinking
      case ( 'ssr'        ) ; cf_standard_name = 'surface_net_upward_shortwave_flux'
      case ( 'ssrd'       ) ; cf_standard_name = 'surface_net_downward_shortwave_flux'
      case ( 'str'        ) ; cf_standard_name = 'surface_net_upward_longwave_flux'
      case ( 'strd'       ) ; cf_standard_name = 'surface_net_downward_longwave_flux'
      case ( 'lsm'        ) ; cf_standard_name = 'land_area_fraction'
      case ( 'albedo'     ) ; cf_standard_name = 'surface_albedo'
      case ( 'sr'         ) ; cf_standard_name = 'surface_roughness_length'
      case ( 'srols'      ) ; cf_standard_name = 'surface_roughness_length'
      case ( 'ci'         ) ; cf_standard_name = 'sea_ice_area_fraction'
      case ( 'sst'        ) ; cf_standard_name = 'sea_surface_temperature'
      case ( 'u10m'       ) ; cf_standard_name = 'eastward_wind'
      case ( 'v10m'       ) ; cf_standard_name = 'northward_wind'
      case ( 'g10m'       ) ; cf_standard_name = 'wind_speed_of_gust'
      case ( 'd2m'        ) ; cf_standard_name = 'dew_point_temperature'
      case ( 't2m'        ) ; cf_standard_name = 'air_temperature'
      case ( 'skt'        ) ; cf_standard_name = 'canopy_temperature'
      case ( 'blh'        ) ; cf_standard_name = 'atmosphere_boundary_layer_thickness'
      case ( 'sshf'       ) ; cf_standard_name = 'surface_downward_sensible_heat_flux'
      case ( 'slhf'       ) ; cf_standard_name = 'surface_downward_latent_heat_flux'
      case ( 'ewss'       ) ; cf_standard_name = 'surface_downward_eastward_stress'
      case ( 'nsss'       ) ; cf_standard_name = 'surface_downward_northward_stress'
      case ( 'cp'         ) ; cf_standard_name = 'lwe_convective_precipitation_rate'
      case ( 'lsp'        ) ; cf_standard_name = 'lwe_large_scale_precipitation_rate'
      case ( 'sf'         ) ; cf_standard_name = 'lwe_thickness_of_snowfall_amount'
      case ( 'sd'         ) ; cf_standard_name = 'lwe_thickness_of_surface_snow_amount'
      case ( 'src'        ) ; cf_standard_name = 'lwe_thickness_of_canopy_water_amount'
      case ( 'swvl1'      ) ; cf_standard_name = 'volume_fraction_of_condensed_water_in_soil'
      case ( 'tv01', 'tv02', 'tv03', 'tv04', 'tv05', &
             'tv06', 'tv07', 'tv08', 'tv09', 'tv10', &
             'tv11', 'tv12', 'tv13', 'tv14', 'tv15', &
             'tv16', 'tv17', 'tv18', 'tv19', 'tv20' )
                              cf_standard_name = 'vegetation_area_fraction'
      case ( 'cvl', 'cvh' ) ; cf_standard_name = 'vegetation_area_fraction'
      case default
        !write (gol,'("do not know cf standard name for tm5 name : ",a)') trim(tm5_name); call goErr
        !TRACEBACK; status=1; return
        ! assume name follows CF already ..
        cf_standard_name = trim(tm5_name)
    end select

    ! ok
    status = 0

  end subroutine TMM_CF_Convert_Name


  ! ***


  subroutine TMM_CF_Convert_Units( tm5_units, cf_units, tm5_to_cf_factor, status )

#ifdef with_udunits
    use UDUnits, only : UDUnits_ConversionFactor
#endif

    ! --- in/out ---------------------------------

    character(len=*), intent(in)                ::  tm5_units
    character(len=*), intent(in)                ::  cf_units
    real, intent(out)                           ::  tm5_to_cf_factor
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TMM_CF_Convert'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef with_udunits
    ! unit conversion factor:
    call UDUnits_ConversionFactor( trim(tm5_units), trim(cf_units), tm5_to_cf_factor, status )
    if ( status /= UDUNITS_NOERR ) then
      write (gol,'("from conversion of TM5 units to CF units:")'); call goErr
      write (gol,'("  TM5 units    : ",a)') trim(tm5_units); call goErr
      write (gol,'("  CF units     : ",a)') trim(cf_units); call goErr
      TRACEBACK; status=1; return
    end if
#else
    ! dummy assignment to avoid compiler warnings about undefined output arguments:
    tm5_to_cf_factor = 1.0
    ! error ...
    write (gol,'("need UDUnits library but macro `with_udunits` not defined")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! ok
    status = 0

  end subroutine TMM_CF_Convert_Units


end module TMM_CF



