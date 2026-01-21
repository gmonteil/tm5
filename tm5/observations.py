#!/usr/bin/env python

from pathlib import Path
import xarray as xr
from pandas import DatetimeIndex, Timestamp, concat
from loguru import logger
from numpy import int16, float64, int32
from dataclasses import dataclass
from omegaconf import DictConfig
from pandas import DataFrame, concat, TimedeltaIndex
from typing import Tuple, List, Callable
from tqdm import tqdm
from glob import glob
from functools import partial
from tm5.debug import trace_args
from types import SimpleNamespace

def read_obspack_file(filename: str,
                      start: Timestamp | str | None = None,
                      end:   Timestamp | str | None = None ) -> SimpleNamespace: #Tuple[DataFrame, DataFrame]:
    ds = xr.open_dataset(filename)
    #
    #-- apply temporal restriction
    #
    if start!=None and end!=None:
        cnd_time = (ds.time>=Timestamp(start)) & (ds.time <= Timestamp(end))
        ds = ds.sel(obs=cnd_time)
    elif start!=None:
        cnd_time = (ds.time>=Timestamp(start))
        ds = ds.sel(obs=cnd_time)
    elif end!=None:
        cnd_time = (ds.time<=Timestamp(end))
        ds = ds.sel(obs=cnd_time)
    #
    obspack_variables = ['longitude', 'latitude', 'altitude', 'time', 'start_time', 'value',]
    # obspack_variables = ['longitude', 'latitude', 'altitude', 'elevation', 'intake_height', 'time', 'start_time', 'value',]
    #
    #-- obspack uncertainties are provides as variable
    #   'value_unc' or 'value_std_dev'
    if  'value_unc' in ds.variables:
        uncvar = 'value_unc'
        obspack_variables.append('value_unc')
    elif 'value_std_dev' in ds.variables:
        obspack_variables.append('value_std_dev')
        uncvar = 'value_std_dev'
    else:
        #-- not clear whether uncertainty might be a "must" downstream,
        #   better raise Exception (or discard this obspack file)
        uncvar = None #-- not clear whether uncertainty m
    #
    #-- some obspack files have elevation/intake_height as variables
    #
    #   NOTE: measurement height can always be accessed from variable
    #         'altitude'
    #
    for var in ['elevation', 'intake_height',]:
        if var in ds.variables:
            obspack_variables.append(var)
    data = ds[obspack_variables].to_dataframe()

    metadata = {k:ds.attrs.get(k) for k in ['site_code', 'site_name', 'site_latitude', 'site_elevation', 'site_elevation_unit', 'site_utc2lst', 'dataset_name', 'dataset_globalview_prefix', 'dataset_parameter', 'dataset_project', 'dataset_platform', 'dataset_selection', 'dataset_selection_tag', 'dataset_calibration_scale', 'dataset_start_date', 'dataset_stop_date', 'dataset_data_frequency', 'dataset_data_frequency_unit', 'dataset_intake_ht', 'dataset_intake_ht_unit', 'dataset_usage_url', 'dataset_usage_description', 'dataset_contribution', 'obspack_name']}
    metadata['filename'] = Path(filename).name
    return SimpleNamespace(data=data, metadata=metadata)

@dataclass
class PointObs:
    dconf : DictConfig
    observations : DataFrame | None = None
    metadata : DataFrame | None = None

    @classmethod
    def from_obspack(self) -> "PointObs":
        data, metadata = [], []
        for file in Path(self.dconf.observations.obspack.folder).glob('*.nc'):
            dat, mdat = read_obspack_file(file)
            data.append(dat)
            metadata.append(mdat)
        self.observations = concat(data)
        self.metadata = concat(metadata)

    def write(self, filename:str) -> str:
        raise NotImplementedError
    
    
def read_obspack_file_(fname: str | Path, start: Timestamp | str, end: Timestamp | str) -> DataFrame | None :
    """
    This should return a DataFrame with at least the following columns: 
    lat, lon, alt, time, obs, err, tracer, station_id
    """
    ds = xr.open_dataset(fname)
    df = ds.sel(obs=(ds.time > Timestamp(start)) & (ds.time <= Timestamp(end))).to_dataframe().reset_index()
    if len(df) == 0:
        msg = f"...reading obspack@{fname} for {start} -- {end} " \
            f"yielded empty frame. Location will be discarded."
        logger.trace(msg)
        return
    #
    #-- there are obspack files for which 'value_std_dev' is missing
    #   (e.g. ch4_wbi_aircraft-pfp_1_allvalid.nc, which has 'value_unc' instead)
    #
    if 'value_std_dev' in df.columns:
        unc_column = 'value_std_dev'
        df = df.loc[:, ['time', 'start_time', 'longitude', 'latitude', 'altitude', 'elevation', 'intake_height', 'value', unc_column]]
    elif 'value_unc' in df.columns:
         unc_column = 'value_unc'
         df = df.loc[:, ['time', 'start_time', 'longitude', 'latitude', 'altitude', 'elevation', 'intake_height', 'value', unc_column]]
    else:
        msg = f"...reading obspack@{fname}, " \
            f"neither 'value_std_dev' nor 'value_unc' found, " \
            f"location will be discarded."
        logger.warning(msg)
        return
    #
    #-- normalise variable names
    #
    df = df.rename(columns={
        'longitude': 'lon', 
        'latitude': 'lat', 
        'altitude': 'alt', 
        'value': 'mixing_ratio', 
        unc_column: 'mixing_ratio_err', 
    })
    df.loc[:, 'station_id'] = ds.attrs['dataset_num']
    df.loc[:, 'station_code'] = ds.attrs['site_code']
    df.loc[:, 'tracer'] = ds.attrs['dataset_parameter'].upper()
    df.loc[:, 'mixing_ratio'] = df.mixing_ratio * 1.e-9
    df.loc[:, 'mixing_ratio_err'] = df.mixing_ratio_err * 1.e-9
    df.loc[:, 'time_window_length'] = TimedeltaIndex(df.time - df.start_time) * 2
    df.loc[:, 'sampling_strategy'] = 2
    return df
    
    
@trace_args('file_list_or_pattern', 'start', 'end')
def read_obspack(file_list_or_pattern: str | List[str], start: Timestamp | str, end: Timestamp | str) -> DataFrame:
    if isinstance(file_list_or_pattern, str):
        file_list_or_pattern = glob(file_list_or_pattern)

    import_fn = partial(read_obspack_file_, start=start, end=end)
    
    # Import all the obs files
    observations = [import_fn(f) for f in tqdm(file_list_or_pattern)]
    
    # Concat, filtering out the None
    observations = concat([o for o in observations if o is not None])
    
    # Return, adding an "id" variable in the process
    return observations.reset_index(drop=True).reset_index().rename(columns={'index':'id'})
    
    
@trace_args('dconf')
def write_point_obs(df: DataFrame, filename: Path | str) -> Path: # dconf: DictConfig, reader: Callable = read_obspack) -> Path:
#    df = reader(
#        dconf.filename,
#        start = Timestamp(dconf.start),
#        end = Timestamp(dconf.end)
#    )
    
    # The part after this should be independent from the reader
    ds = df.set_index('id').to_xarray()
    
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
    ds['time_window_length'] = (ds.time_window_length).astype(int32)
    ds['sampling_strategy'] = ds.sampling_strategy.astype(int16)
    ds['date_components'] = ds.date_components.astype(int16)
    ds['station_id'] = ds.station_id.astype(int32)
    
    # Determine the filename and ensure that the parent exist
    # filename = Path(dconf.input_dir) / 'point_input.nc4'
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    
    for itrac, tracer in enumerate(set(df.tracer)):
        mode = 'a' if itrac > 0 else 'w'
        ds.sel(id = ds.tracer.str.lower() == tracer.lower()).to_netcdf(filename, group=tracer, mode=mode)
    
    logger.info(f"Wrote {filename}")
        
    return filename
    
        
def prepare_point_obs_(dconf : DictConfig) -> Path:

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
    ds['time_window_length'] = (ds.time_window_length * 1.e-9).astype(int32)
    ds['sampling_strategy'] = ds.sampling_strategy.astype(int16)
    ds['date_components'] = ds.date_components.astype(int16)
    ds['station_id'] = ds.station_id.astype(int32)

    for itrac, tracer in enumerate(dconf.tracers):
        mode = 'a' if itrac > 0 else 'w'
        ds.sel(id = ds.tracer.str.lower() == tracer.lower()).to_netcdf(filename, group=tracer, mode=mode)

    logger.info(f"Wrote {filename}")

    return filename
