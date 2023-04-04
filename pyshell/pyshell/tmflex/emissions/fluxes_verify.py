#!/bin/env python

from pyshell.tmflex.emissions.emissions import Emissions
import os
from pyshell.base.helper.Utilities import checkDir
from datetime import timedelta
import xarray as xr
from pyshell.gridtools import TM5Grids
from pandas import DatetimeIndex, Timedelta
from pandas.tseries.frequencies import to_offset
import logging


logger = logging.getLogger(__name__)


def crop_and_coarsen(glo1x1, latb, lonb):
    """
    Crop and coarsen a global 1x1 emission field to a subregion/lower resolution
    """

    # First, crop:
    data = glo1x1.sel(latitude=slice(latb[0], latb[-1]), longitude=slice(lonb[0], lonb[-1])).values
    
    # Then coarsen:
    dlat = latb[1] - latb[0]
    dlon = lonb[1] - lonb[0]
    assert dlat - int(dlat) == 0
    assert dlon - int(dlon) == 0
    nt, nlat, nlon = data.shape
    data = data.reshape(nt, -1, int(dlat), nlon).sum(2)
    nlat = data.shape[1]
    data = data.reshape(nt, nlat, -1, int(dlon)).sum(3)
    return data
    

class PreprocessedEmissions(Emissions):
    def __init__(self, rcf, dconf, tracer='CO2', *args, **kwargs):
        Emissions.__init__(self, rcf, *args, **kwargs)
        self.dconf = dconf
        self.tracer = tracer
        self.writeCycle = self.rcf.get('CO2.emission.dailycycle', 'bool')
        self.MolarMass = {'CO2':44.00995e-3}[tracer]
        self.granularity = timedelta(hours=3)
        self.output_dir = os.path.dirname(self.emission_file_name)
        
    def read_preprocessed_emis(self, category, region, period):
        """
        category : name of the category to import
        region   : region definition
        period   : slice(start, end, timestep)
        """
        file_pattern = self.dconf[category].pattern
        field = self.dconf[category].get('field', 'emis')
        area = self.dconf[category].get('area_field', None)
        data = xr.open_mfdataset(file_pattern)
        
        if area is None :
            area = TM5Grids.global1x1().area
        else :
            area = data[area].values
            
        # Load the emissions. They should be on a global 1x1, 3-hourly grid and in mol/m2/s
        # The following will convert them in kg[tracer]/s
        idx = DatetimeIndex(data.time.values)
        emis_glo1x1 = data[field].sel(time=(idx >= period.start) & (idx < period.stop)) * area * self.MolarMass
        
        # destination region
        destreg = TM5Grids.from_corners(latb = self.lat_grid[region], lonb = self.lon_grid[region])
        
        # regrid:
        emis_coarsened = crop_and_coarsen(emis_glo1x1, latb=destreg.latb, lonb=destreg.lonb)
        
        # Calculate granularity:
        # Granularity (time step) is read, by order of priority 1) in the netcdf file attributes, 2) in the yaml file
        granularity = data.attrs.get('granularity', self.dconf[category].get('granularity'))
        timestep = DatetimeIndex(emis_glo1x1.time.values) + to_offset(granularity) - DatetimeIndex(emis_glo1x1.time.values)
        
        # Return a DataArray:
        return xr.Dataset(
            data_vars = {
                'emis': (('time', 'latitude', 'longitude'), emis_coarsened),
                'timestep': (('time'), [t.total_seconds() for t in timestep])
            },
            coords={
                'time': DatetimeIndex(emis_glo1x1.time.values),
                'latitude': destreg.latc,
                'longitude': destreg.lonc,
            },
            attrs = {
                'freq': to_offset(granularity),
                'dailycycle': timestep[0] < Timedelta(days=1)
            }
        )
        
    def LoopThroughPeriods(self):
        """
        A single call to self.read*Flux gives the flux (in Kg CO2/cell/sec) during a three hour time window for a single
        category. This routine aggregates those fluxes into total flux over a given period, and writes out the daily
        cycle files. This is complicated if the granularity of emissions is different for different categories. For each
        category, we will wrxr.Dataset(ite a different daily cycle file, and at the end of the subroutine we will assemble them
        together into one daily cycle file per day.
        """
        for ireg, region in enumerate(self.zoom_regions_names):
            categories = self.Emission[region][self.tracer]['categories']
            for icat, cat in enumerate(categories):
                for time_index, _ in enumerate(self.Emission[region][self.tracer][cat]['time_interval']['time_mid']):

                    # Get time interval
                    time_start = self.Emission[region][self.tracer][cat]['time_interval']['time_start'][time_index]
                    time_end = self.Emission[region][self.tracer][cat]['time_interval']['time_end'][time_index]
                    sec_period = (time_end-time_start).total_seconds()
                    
                    emission_fine = self.read_preprocessed_emis(cat, region, slice(time_start, time_end))
                    emission_coarse = (emission_fine.emis * emission_fine.timestep).sum('time').values
                    
                    #for emis in emission_fine:
                    logger.info("%20s %s-%s %10s : %20.6f Tg" % (cat, time_start.strftime("%Y/%m/%d"), time_end.strftime("%Y/%m/%d"), region, emission_coarse.sum()*1.0e-9))
                    emission_coarse = emission_coarse / sec_period  # bring back to per second
                    self.Emission[region][self.tracer][cat]['emission_data'][time_index, :, :] = emission_coarse

                    if emission_fine.dailycycle:
                        emission_anomaly = emission_fine.emis - emission_fine.emis.mean('time')
                        cur_day = time_start
                        while cur_day < time_end :
                            anomaly_day = emission_anomaly.sel(time=
                                (emission_anomaly.time.dt.year == cur_day.year) &
                                (emission_anomaly.time.dt.month == cur_day.month) &
                                (emission_anomaly.time.dt.day == cur_day.day)
                            )
                            filename = cur_day.strftime(self.dconf.dailycycle_filename_format)
                            checkDir(filename)
                            ds = xr.Dataset(data_vars = {'emission_anomaly': (('timesteps', 'latitude', 'longitude'), anomaly_day.values)})
                            write_mode = 'w' if icat == 0 else 'a'
                            ds.to_netcdf(filename, group=region + '/' + cat, mode=write_mode)
                            cur_day += timedelta(days=1)

                            