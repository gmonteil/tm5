#!/usr/bin/env python

from pathlib import Path
import xarray as xr
from pandas import DatetimeIndex, Timestamp
from loguru import logger
from numpy import int16, float64, int32


def prepare_point_obs(dconf) -> Path:

    # Determine the filename and ensure that the parent exist
    filename = Path(dconf.input_dir) / 'point_input.nc4'
    filename.parent.mkdir(parents=True, exist_ok=True)

    # Load the observations:
    logger.info(f"Loading observations from {dconf.filename}")
    ds = xr.open_mfdataset(dconf.filename)

    # Replace the time dimension by a regular index
    ds = ds.to_dataframe().reset_index().to_xarray()
    logger.info(f"Done ...")

    # Filter the time interval:
    ds = ds.sel(index = (ds.time >= Timestamp(dconf.start)) & (ds.time < Timestamp(dconf.end)))

    # Rename variables:
    ds = ds.rename_dims(dict(index='id')).rename_vars({
        'latitude': 'lat',
        'longitude': 'lon',
        'sampling_altitude': 'alt',
        'obs': 'mixing_ratio',
        'err_obs': 'mixing_ratio_err',
        'index': 'id',
        'TM5_station_ID': 'station_id'
    })

    # Add a "date_components" variable that TM5 can read:
    ds['date_components'] = xr.DataArray(
        [t.timetuple()[:6] for t in DatetimeIndex(ds.time)],
        dims=('id', 'idate')
    )

    # Adjust the file types:
    ds['lat'] = ds.lat.astype(float64)
    ds['lon'] = ds.lon.astype(float64)
    ds['alt'] = ds.alt.astype(float64)
    ds['mixing_ratio'] = ds.mixing_ratio.astype(float64)
    ds['mixing_ratio_err'] = ds.mixing_ratio_err.astype(float64)
    ds['id'] = ds['id'].astype(int32)
    ds['time_window_length'] = ds.time_window_length.astype(int32)
    ds['sampling_strategy'] = ds.sampling_strategy.astype(int16)
    ds['date_components'] = ds.date_components.astype(int16)
    ds['station_id'] = ds.station_id.astype(int32)

    for itrac, tracer in enumerate(dconf.tracers):
        mode = 'a' if itrac > 0 else 'w'
        ds.sel(id = ds.tracer.str.lower() == tracer.lower()).to_netcdf(filename, group=tracer, mode=mode)

    logger.info(f"Wrote {filename}")

    return filename
