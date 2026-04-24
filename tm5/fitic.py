#!/usr/bin/env python

import xarray as xr
from pathlib import Path
from omegaconf import DictConfig
from pandas import DataFrame
import numpy as np
from loguru import logger


def read_obs_table(filename: Path | str) -> DataFrame:
    obs_table = xr.open_dataset(filename).to_dataframe()
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
                
