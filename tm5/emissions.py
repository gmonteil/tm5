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


def coarsen_file(path_or_pattern : str, reg : TM5Grids, start : Timestamp, end : Timestamp) -> xr.Dataset:
    # Open the source file(s). It should be on a global 1x1 resolution (for now ...)
    ds = xr.open_mfdataset(path_or_pattern)

    # Define a regridder
    glb1x1 = TM5Grids.global1x1()
    regridder = glb1x1.regridder(reg)

    # Regrid and ensure the result is on a daily time step (for now ...)
    coarse = regridder(ds)
    coarse = coarse.reindex(time=date_range(start, end, freq='D')).ffill('time')

    # TM5 expects float64 precision 
    coarse = coarse.astype('f8')
    # Return ordered the way TM5 wants it!
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
