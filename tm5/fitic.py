#!/usr/bin/env python

import xarray as xr
from pathlib import Path
from omegaconf import DictConfig
from pandas import DataFrame, DatetimeIndex
import numpy as np
from loguru import logger
from collections import OrderedDict
from types import SimpleNamespace
#-- library packages
from tm5.gridtools import TM5Grids

#
#-- table of FIT-IC regions
#   - their will be one (global) instance for all queries
#   TODO: should not be hard-coded here!
#
fitic_region_table = OrderedDict()
def _fitic_region_table_init():
    global fitic_region_table
    if len(fitic_region_table)>0:
        return #-- already initialised
    #
    #-- create grid instances
    #
    glb_grid = TM5Grids.from_corners(west=-180, east=180, south=-90, north=90, dlon=6, dlat=4)
    eur_grid = TM5Grids.from_corners(west=-36, east=54, south=22, north=74, dlon=3, dlat=2)
    gns_grid = TM5Grids.from_corners(west=0, east=18, south=42, north=58, dlon=1, dlat=1)
    #
    #-- initial attributes (per region)
    #
    fitic_region_table['glb600x400'] = SimpleNamespace(
        grid=glb_grid,
        rlat=4, rlon=6,
        child='eur300x200', parent=None)
    fitic_region_table['eur300x200'] = SimpleNamespace(
        grid=eur_grid,
        rlat=2, rlon=3,
        child='gns100x100', parent='glb600x400')
    fitic_region_table['gns100x100'] = SimpleNamespace(
        grid=gns_grid,
        rlat=1, rlon=1,
        child=None, parent='eur300x200')
    #
    #-- extending attributes
    #   - 1D flattened arrays of coordinates
    #   - prepare for filtering out child domain within in parent
    #
    for reg,reg_info in fitic_region_table.items():
        grid = reg_info.grid
        area1D = grid.area.ravel()
        #
        lonmesh,latmesh = np.meshgrid(grid.lonc,grid.latc)
        fitic_region_table[reg].lonmesh = lonmesh
        fitic_region_table[reg].latmesh = latmesh
        #
        #--
        #
        grid_mask = xr.DataArray(
            np.ones((grid.nlat, grid.nlon), dtype='i1'),
            dims = ('lat','lon'),
            coords = { 'lon' : grid.lonc,
                       'lat' : grid.latc },
            name = 'mask'
        )
        if reg=='glb600x400':
            #-- drop the HALO corrected child domain
            child_info = fitic_region_table[reg_info.child]
            child_grid = child_info.grid
            #
            lon_min = child_grid.west + grid.dlon
            lon_max = child_grid.east - grid.dlon
            lat_min = child_grid.south + grid.dlat
            lat_max = child_grid.north - grid.dlat
            drop_mask = (
                (grid_mask.lat >= lat_min) &
                (grid_mask.lat <= lat_max) &
                (grid_mask.lon >= lon_min) &
                (grid_mask.lon <= lon_max)
                )
        elif reg=='eur300x200':
            #-- drop the HALO corrected child domain
            child_grid = fitic_region_table[reg_info.child].grid
            #
            lon_min = child_grid.west + grid.dlon
            lon_max = child_grid.east - grid.dlon
            lat_min = child_grid.south + grid.dlat
            lat_max = child_grid.north - grid.dlat
            drop_mask1 = (
                (grid_mask.lat >= lat_min) &
                (grid_mask.lat <= lat_max) &
                (grid_mask.lon >= lon_min) &
                (grid_mask.lon <= lon_max)
                )
            #-- drop HALO part of domain itself
            parent_grid = fitic_region_table[reg_info.parent].grid
            drop_mask2 = (
                (grid_mask.lat<=grid.south + parent_grid.dlat) |
                (grid_mask.lat>=grid.north - parent_grid.dlat) |
                (grid_mask.lon<=grid.west + parent_grid.dlon) |
                (grid_mask.lon>=grid.east - parent_grid.dlon)
            )
            drop_mask = drop_mask1 | drop_mask2
        elif reg=='gns100x100':
            #-- drop HALO part of domain itself
            parent_grid = fitic_region_table[reg_info.parent].grid
            drop_mask = (
                (grid_mask.lat<=grid.south + parent_grid.dlat) |
                (grid_mask.lat>=grid.north - parent_grid.dlat) |
                (grid_mask.lon<=grid.west + parent_grid.dlon) |
                (grid_mask.lon>=grid.east - parent_grid.dlon)
            )
        #
        #-- set drop_mask
        #
        fitic_region_table[reg].drop_mask = drop_mask


def get_fitic_region_table():
    if len(fitic_region_table)==0:
        msg = f"...initialise region table"
        logger.info(msg)
        _fitic_region_table_init()
    return fitic_region_table


def tm5emisdir_load_emissions2D( emisdir : str | Path, emis_prefix : str, day_range : DatetimeIndex, regions : list, drop : bool = False ) -> SimpleNamespace:
    """Read in daily emissions as prepared for TM5 for the selected temporal range
    and regions.
    The emissions array will be 2D with only one single dimension in the spatial domain,
    and concatenating the contributions from each region in the spatial domain as well.
    """
    nday = len(day_range)
    #
    #-- get spatial information as 1D vector
    #
    region_table = get_fitic_region_table()
    ng = 0
    regionid_1D = []
    lon_1D = None
    lat_1D = None
    area_1D = None
    for region,region_info in region_table.items():
        if drop:
            keep_mask = ~region_info.drop_mask
            ng_reg = np.count_nonzero(keep_mask)
            lon_reg = region_info.lonmesh[keep_mask]
            lat_reg = region_info.latmesh[keep_mask]
            area_reg = region_info.grid.area[keep_mask]
        else:
            ng_reg = region_info.grid.nlat*region_info.grid.nlon
            lon_reg = region_info.lonmesh.ravel()
            lat_reg = region_info.latmesh.ravel()
            area_reg = region_info.grid.area.ravel()
        ng += ng_reg
        regionid_1D += [region,]*ng_reg
        if lon_1D is None:
            lon_1D = lon_reg
        else:
            lon_1D= np.hstack((lon_1D,lon_reg))
        if lat_1D is None:
            lat_1D = lat_reg
        else:
            lat_1D= np.hstack((lat_1D,lat_reg))
        if area_1D is None:
            area_1D = area_reg
        else:
            area_1D= np.hstack((area_1D,area_reg))
    msg = f"...preparing emissions for nday={nday} and ng={ng}"
    logger.info(msg)
    #
    #-- prepare emissions field
    #
    missval = -99999.
    emissions2D = np.full((nday,ng), missval)
    for iday,day in enumerate(day_range):
        emis_list = []
        for reg in regions:
            #-- TODO:: 'CH4' is still hard-coded here!
            fpath = Path(emisdir) / day.strftime(f'{emis_prefix}.CH4.{reg}.%Y%m%d.nc')
            if not fpath.exists():
                msg = f"expected emissions file ***{str(fpath)}*** not found on system."
                raise FileNotFoundError(msg)
            #
            #-- open emissions file for current region
            #
            em = xr.open_dataset(fpath)
            #
            #-- sum-up total emissions
            #
            emtot = em.to_array().sum('variable').values
            #
            #-- turn lat/lon into 1D vector (lat major ordering)
            #
            if drop:
                drop_mask = region_table[reg].drop_mask
                keep_mask = ~drop_mask
                emtot = emtot[keep_mask]
            else:
                emtot = emtot.ravel()
            emis_list.append(emtot)
        #-- add concatenated emissions
        emissions2D[iday,:] = np.hstack(emis_list)
    #-- consistency: emissions should be filled completely!
    assert np.count_nonzero(emissions2D==missval)==0
    
    return SimpleNamespace(emis2D=emissions2D,
                           reg1D=np.array(regionid_1D),
                           lonc1D=lon_1D, latc1D=lat_1D,
                           area1D=area_1D,
                           region_table=region_table,
                           emisdir=emisdir)


def read_obs_table(filename: Path | str, drop_missing_value : bool = False) -> DataFrame:
    """
    """
    obs_ds = xr.open_dataset(filename)
    #-- 'mixing_ratio' is the variable providing the observed concentrations
    obs_missval = None
    try:
        if obs_ds.mixing_ratio.comment.startswith('missing value set to'):
            obs_missval = float(obs_ds.mixing_ratio.comment.split()[-1])
            # print(f"detected observation missing-value -->{obs_missval}<--")
    except AttributeError:
        pass
    #
    #-- prefer dataframe upstream
    #
    obs_table = obs_ds.to_dataframe()
    obs_ds.close()
    #
    #-- remove missing values (if any)
    #
    if drop_missing_value and (not obs_missval is None):
        #-- detect number of missing values
        nmiss = obs_table['mixing_ratio'].value_counts().get(obs_missval, 0)
        if nmiss>0:
            msg = f"...detected {nmiss} missing values in observations file -->{filename}<--"
            logger.debug(msg)
            cnd_nomiss = obs_table['mixing_ratio']!=obs_missval
            obs_table = obs_table.loc[cnd_nomiss,:]
    #
    #-- adjust type sampling_strategy/time_window_length
    #
    convert_dict = {
        'sampling_strategy': np.int32,
        'time_window_length': np.int32,
        'station_id': np.int32
        }
    obs_table = obs_table.astype(convert_dict)
    
    # Add a "obsid" column to identify each "tracer", if not already present in the obs table
    if 'obsid' not in obs_table:
        for tracer in set(obs_table.tracer):
            df = obs_table[obs_table.tracer == tracer].reset_index()
            obsid = df.tracer + '_' + df.index.astype(str)
            obs_table.loc[obs_table.tracer == tracer, 'obsid'] = obsid.values
    #
    return obs_table

    
def create_departure_files(dconf: DictConfig):
    """
    The departure files have the following columns:
    - sampling_strategy, lat, lon, alt, forcing, time_window_length, nsamples, total_weight, date_components
    
    File structure:
      /region
        /tracer
          /variables
            
    Origin of the data:
    - point_output.nc4: total_weight, nsamples, sampling_strategy
    - point_input.nc4: lat, lon, alt, time_window_length, date_components
    - imposed: forcing
    """
    reglist = dconf.run.regions
    trlist = dconf.run.tracers   # This refers to the *forward* tracers!!!
    # Remove any pre-existing file
    dep_file = Path(dconf.run.paths.output) / 'point/point_departures.nc4'
    dep_file.unlink(missing_ok=True)
    
    for tracer in trlist:
        inp = xr.open_dataset(Path(dconf.output.point.input_dir) / 'point_input.nc4', group=tracer)
        for region in reglist:
            try:
                outp = xr.open_dataset(Path(dconf.run.paths.output) / 'point/point_output.nc4', group=f'{region}/{tracer}')

                # Select the slice of the input file corresponding to that region
                reg_inp = inp.isel(id=inp.id.isin(outp.id))
            
                # Copy the data:
                for varname in ['lat', 'lon', 'alt', 'time_window_length', 'obsid']:
                    outp[varname] = ('samples', reg_inp[varname].values)
                
                outp['date_components'] = (('samples', 'idate'), reg_inp.date_components.values)
                
                # Set the "forcing" to 1:
                outp['forcing'] = ('samples', [1.] * outp.sizes['samples'])
            
                # Write into the dep file:
                for iobs in range(outp.sizes['samples']):
                
                    # Select the subset containing just one obs
                    data_out = outp.isel(samples=[iobs])
                
                    # Write it to its own netCDF group, following the "obsid" value
                    data_out.to_netcdf(dep_file, group=f'{region}/{outp.obsid.values[iobs]}', mode='a')
            except OSError:
                # This can happen when there is no obs in that region. In this case, just cycle to the next ...
                logger.info(f"No observations found for region {region} in point_output file")
                pass
