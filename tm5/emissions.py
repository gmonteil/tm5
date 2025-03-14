#!/usr/bin/env python

"""
Prepare daily emissions from the requested grid
"""

from omegaconf import DictConfig
from tm5.gridtools import TM5Grids
import xarray as xr
from pandas import date_range, Timestamp
from pathlib import Path
from loguru import logger
from tqdm.auto import tqdm


def coarsen_file(path_or_pattern : str, reg : TM5Grids, start : Timestamp, end : Timestamp, use_esmf_regridder : bool = False) -> xr.Dataset:
    # Open the source file(s). It should be on a global 1x1 resolution (for now ...)
    ds = xr.open_mfdataset(path_or_pattern)

    #-- expected grid of global emission inputs
    # !!!TODO:should add consistency check
    glb1x1 = TM5Grids.global1x1()

    #
    #-- MVO 2025-03-14::
    #   - re-gridding results seem not appropriate
    #     with ESMF regridder used underneath,
    #     by default apply again approach as used in old implementation
    #
    # Regrid
    if use_esmf_regridder:
        # Define a regridder
        regridder = glb1x1.regridder(reg)
        coarse = regridder(ds)
    else:
        #
        #-- iterate over (emission) datasets and
        #   (1) crop to target region
        #   (2) re-grid to target resolution
        data_coarse_dict = {}
        for da_name in ds.data_vars:
            if da_name=='area':
                continue
            da = ds[da_name]
            #-- make sure units is as expected
            if not da.units=='kg/cl/s':
                msg = f"emission {da_name} with unexpected units ==>{da.units}<=="
                raise RuntimeError(msg)
            #-- (1) crop (result is numpy array)
            data_out = da.sel(lat=slice(reg.south, reg.north),
                              lon=slice(reg.west, reg.east)).values
            #-- (2) coarsen (only if needed):
            if (reg.dlat!=glb1x1.dlat) or (reg.dlon!=glb1x1.dlon):
                assert (reg.dlat - int(reg.dlat) == 0) & (reg.dlon - int(reg.dlon) == 0), "lat and lon steps must be multiple of the base (1°) step"
                nt, nlat, nlon = data_out.shape
                data_out = data_out.reshape(nt, -1, int(reg.dlat), nlon).sum(2)
                nlat = data_out.shape[1]
                data_out = data_out.reshape(nt, nlat, -1, int(reg.dlon)).sum(3)
            data_coarse_dict[da_name] = xr.DataArray(
                data_out,
                dims = ('time', 'lat', 'lon'),
                coords = {
                    'time': ds.time,
                    'lat': reg.latc,
                    'lon': reg.lonc
                }
            )
        #
        #-- create xarray dataset
        #
        coarse = xr.Dataset(
            coords=dict(lon=('lon', reg.lonc),
                        lat=('lat', reg.latc),
                        time=ds['time']),
            data_vars=data_coarse_dict
            )

    #
    # ensure result is on a daily time step (at least for  now ...)
    #
    coarse = coarse.reindex(time=date_range(start, end, freq='D')).ffill('time')
    #
    # ensure float64 precision (as TM5 I/O expects)
    #
    coarse = coarse.astype('f8')

    return coarse


def prepare_emissions(conf: DictConfig) -> None:
    
    start = Timestamp(conf.run.start)
    end = Timestamp(conf.run.end)
    
    for region in conf.run.regions:
        
        # Create the region object:
        reg = conf.regions[region]
        reg = TM5Grids.from_corners(west=reg.lons[0], east=reg.lons[1], south=reg.lats[0], north=reg.lats[1], dlon=reg.lons[2], dlat=reg.lats[2])
        
        # Regrid the emissions for that region and tracer
        for trname, tracer in conf.emissions.items():

            if 'categories' not in tracer:
                continue
            
            datasets = {}
            
            # Do the actual coarsening.
            for catname, cat in tracer.categories.items():
                datasets[catname] = coarsen_file(cat.path, reg, start, end)[cat.field]
            # Ensure that the dest path exists!
            Path(tracer.prefix).parent.mkdir(exist_ok=True, parents=True)

            datasets = xr.Dataset(datasets)
            
            for day in tqdm(date_range(start, end, freq='D'), desc=f'Writing emission files for {trname} in region {region}'):
                destfile = f'{tracer.prefix}.{trname}.{region}.{day:%Y%m%d}.nc'
                datasets.sel(time=day).to_netcdf(destfile)
                
#            # group the datasets in daily emission files for that region and tracer:
#            for day in tqdm(date_range(start, end, freq='D'), desc=f'Writing emission files for {trname} in region {region}'):
#                destfile = f'{tracer.prefix}.{trname}.{region}.{day:%Y%m%d}.nc'
#                xr.Dataset({k : datasets[k].sel(time=day) for k in datasets}).to_netcdf(destfile)
