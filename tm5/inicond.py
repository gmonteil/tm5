#!/usr/bin/env python
import os
import urllib.request
from pandas import Timestamp
from omegaconf import DictConfig
from tm5 import gridtools as gt
import xarray as xr
from pathlib import Path
from loguru import logger


def get_iniconc_carbontracker(url_pattern : str, date : Timestamp | str, regions : dict | DictConfig, filename : str | Path):

    # Download the file
    date = Timestamp(date)
    url = date.strftime(url_pattern)
    ct_filename = url.split('/')[-1]

    logger.info(f"Reading initial condition from file {ct_filename}")
    if not Path(ct_filename).exists():
        logger.info("Downloading the file ...")
        urllib.request.urlretrieve(url, url.split('/')[-1])

    co2_3x2 = xr.open_dataset(ct_filename)
    glb3x2 = gt.TM5Grids.global3x2()

    # Regrid it onto TM5 resolution
    file_mode = 'w'
    for regname, region in regions.items():

        destgrid = gt.TM5Grids.from_corners(
            west=region['lons'][0], east=region['lons'][1], dlon=region['lons'][2],
            south=region['lats'][0], north=region['lats'][1], dlat=region['lats'][2]
        )

        co2 = gt.SpatialData(co2_3x2.co2.values[0, :, :, :], glb3x2, lon_axis=2, lat_axis=1, density=True).regrid(destgrid).data

        co2 = xr.Dataset(
            data_vars= {'mixing_ratio': (('level', 'latitude', 'longitude'), co2)},
            coords = {
                'level': co2_3x2.level.values,
                'latitude': destgrid.latc,
                'longitude': destgrid.lonc
            }
        )

        logger.info(f'Writing CO2 initial condition to group CO2/{regname} of file {filename}')
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
        co2.to_netcdf(filename, group=f'{regname}/CO2', mode=file_mode)
        file_mode = 'a'
