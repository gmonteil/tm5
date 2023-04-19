#!/usr/bin/env python

import xarray as xr
from pandas import Timestamp, date_range, DatetimeIndex
from pandas.tseries.frequencies import to_offset
from omegaconf import DictConfig
from tm5.gridtools import TM5Grids
from typing import Tuple, Dict
from netCDF4 import Dataset
from loguru import logger
from numpy import int16
from pathlib import Path
from tm5.units import units_registry as ureg


def crop_and_coarsen(glo1x1: xr.DataArray, dest: TM5Grids) -> xr.DataArray:
    """
    Crop and coarsen a global 1x1 emission field to a subregion at equal or lower resolution.
    """

    # First, crop:
    data = glo1x1.sel(
        latitude=slice(dest.south, dest.north),
        longitude=slice(dest.west, dest.east)).values

    # Then, coarsen:
    assert (dest.dlat - int(dest.dlat) == 0) & (dest.dlon - int(dest.dlon) == 0), "lat and lon steps must be multiple of the base (1°) step"
    nt, nlat, nlon = data.shape
    data = data.reshape(nt, -1, int(dest.dlat), nlon).sum(2)
    nlat = data.shape[1]
    data = data.reshape(nt, nlat, -1, int(dest.dlon)).sum(3)

    return xr.DataArray(
        data,
        dims = ('time', 'latitude', 'longitude'),
        coords = {
            'time': glo1x1.time,
            'latitude': dest.latc,
            'longitude': dest.lonc
        })


def load_preprocessed_emis(categ: DictConfig, tracer: DictConfig, destreg: TM5Grids, period) -> Tuple[xr.DataArray, Dict]:
    """
    This reads pre-processed emissions (global, 1x1˚, mol/m2/s) and convert them to emissions in kg[tracer]/s limited to the requested TM5 region.
    """

    file_pattern = categ.pattern
    field = categ.get('field', 'emis')
    area = categ.get('area_field', None)
    data = xr.open_mfdataset(file_pattern)

    # Ensure we have an area array in memory (should be for global1x1)
    emis_glo1x1 = data[field].sel(time=(data.time >= period.start) & (data.time < period.stop))

    # Convert from mol/m2/s to kg[tracer]/s:
    if area is None :
        area = TM5Grids.global1x1().area
    else :
        area = data[area].values
    emis_glo1x1 *= area * ureg.Quantity('mol').to(f'kg{tracer.species}').m

    # Regrid:
    emcoarse = crop_and_coarsen(emis_glo1x1, dest=destreg)

    # Split between baselines and anomalies (daily cycle files):
    return split_baseline_anomalies(emcoarse, period.start, period.stop, categ.optim_freq)


def split_baseline_anomalies(emreg: xr.DataArray, start: Timestamp, end: Timestamp, freq: str) -> Tuple[xr.DataArray, Dict]:
    baseline = emreg.resample(time=freq).mean()
    anom = {}

    assert freq in ['D', 'MS']
    baseline['time'] = date_range(start, end, freq=freq, inclusive='left')

    # Calculate daily cycles:
    for day in date_range(start, end, freq='D', inclusive='left'):
        emis = emreg.sel(time=DatetimeIndex(emreg.time.dt.date) == day).values

        # Determine the corresponding period:
        bline = baseline.sel(time = baseline.time <= day).values[0, :, :] # 1st interval whose start time is at max equal to the current day start time

        anom[day] = xr.Dataset(
            data_vars = {'emission_anomaly': (('timesteps', 'latitude', 'longitude'), emis - bline)}
        )

    return baseline, anom


def prepare_emissions(dconf : DictConfig, filename : Path) -> Path:

    with Dataset(str(filename), 'w') as ds :
        ds.createDimension('itime', 6)
        for regname, region in dconf.regions.items():

            # Create region object
            reg = TM5Grids.from_corners(
                east = region.lons[1],
                west = region.lons[0],
                south = region.lats[0],
                north = region.lats[1],
                dlon = region.lons[2],
                dlat = region.lats[2]
            )

            # Create netCDF group
            gid = ds.createGroup(regname)
            gid.createDimension('latitude', reg.nlat)
            gid.createDimension('longitude', reg.nlon)
            gid.latmin = reg.south
            gid.latmax = reg.north
            gid.lonmin = reg.west
            gid.lonmax = reg.east

            # Add the grid cell area
            gid.createVariable('area', 'd', ('latitude', 'longitude'))
            gid['area'][:] = reg.area

            # Create tracer subgroup(s):
            for trname, tracer in dconf.tracers.items():

                # Ensure that the folder for dailycycle files exists:
                Path(tracer.dailycycle_filename_format).parent.mkdir(parents=True, exist_ok=True)

                # Create the tracer group
                gtr = gid.createGroup(trname)

                # Create the categories subgroups
                firstcat = True
                for catname, cat in dconf[regname][trname].items():

                    logger.info(f'{regname}; {trname}; {catname}')

                    grp = gtr.createGroup(catname)

                    # Load the emissions
                    emcat, anomcat = load_preprocessed_emis(
                        cat, tracer, reg, slice(Timestamp(dconf.start), Timestamp(dconf.end))
                    )

                    # Write the anomalies for the category
                    for day, anom in anomcat.items():
                        fname = day.strftime(tracer.dailycycle_filename_format)
                        Path(fname).parent.mkdir(parents=True, exist_ok=True)
                        anom.to_netcdf(
                            fname,
                            group=f'{regname}/{catname}',
                            mode = 'w' if firstcat else 'a'
                        )

                    # Ensure that the next cats don't overwrite the dailycycle files
                    firstcat = False

                    # Store the emission themselves:
                    grp.createDimension('nt', emcat.shape[0])
                    grp.createVariable('emission', 'd', ('nt', 'latitude', 'longitude'))
                    grp.createVariable('time_start', 'i2', ('nt', 'itime'))
                    grp.createVariable('time_end', 'i2', ('nt', 'itime'))
                    grp.time_resolution = cat.optim_freq
                    grp.treskey = {'MS': 'monthly', 'D': 'daily'}[cat.optim_freq]
                    grp.optimize = int(cat.get('optim', 0))

                    grp['emission'][:] = emcat.values
                    grp['time_start'][:] = [t.timetuple()[:6] for t in DatetimeIndex(emcat.time)]
                    grp['time_end'][:] = [t.timetuple()[:6] for t in DatetimeIndex(emcat.time) + to_offset(cat.optim_freq)]

    return filename
