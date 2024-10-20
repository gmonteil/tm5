#!/usr/bin/env python

# Prepare daily emission files for the requested grid

import xarray as xr
from tm5.gridtools import TM5Grids
from pandas import Timestamp, date_range
from tm5.settings import load_config
import os
from pathlib import Path

def coarsen_file(path_or_pattern : str, reg : TM5Grids, start : Timestamp, end : Timestamp) -> xr.Dataset:
    # Open the source file(s). It should be on a global 1x1 resolution (for now ...)
    print(path_or_pattern)
    ds = xr.open_mfdataset(path_or_pattern)

    # Define a regridder
    glb1x1 = TM5Grids.global1x1()
    regridder = glb1x1.regridder(reg)

    # Regrid and ensure the result is on a daily time step (for now ...)
    coarse = regridder(ds)
    coarse.reindex(time=date_range(start, end, freq='D')).ffill('time')

    # Return ordered the way TM5 wants it!
    return coarse.transpose('lon', 'lat', 'time')

# Load config file
conf = load_config('forward.yaml', host=os.environ['TM5_HOST'])

# Read start and end
start = Timestamp(conf.run.start)
end = Timestamp(conf.run.end)

# This should write one emission file per day, region and tracer, containing all the emissions of that tracer
for region in conf.run.regions :

    # Create the region object
    reg = conf.regions[region]
    reg = TM5Grids.from_corners(west=reg.lons[0], east=reg.lons[1], south=reg.lats[0], north=reg.lats[1], dlon=reg.lons[2], dlat=reg.lats[2])

    # regrid the emissions for that region and tracer
    for trname, tracer in conf.emissions.items():
        datasets = {}

        # Do the actual coarsening.
        for catname, cat in tracer.categories.items():
            datasets[catname] = coarsen_file(cat.path, reg, start, end)[cat.field]

        # Ensure that the dest path exists!
        Path(tracer.prefix).parent.mkdir(exist_ok=True, parents=True)

        # group the datasets in daily emission files for that region and tracer:
        for day in date_range(start, end, freq='D'):
            destfile = f'{tracer.prefix}.{trname}.{region}.{day:%Y%m%d}.nc'
            xr.Dataset({k:datasets[k].sel(time=day).transpose() for k in datasets}).to_netcdf(destfile)