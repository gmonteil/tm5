#!/bin/env python
import subprocess

from pyshell.tmflex.emissions.emissions import Emissions
from netCDF4 import Dataset, num2date
import os
from pyshell.tm5_utils import redistrib_flux
from pyshell.base.helper.Utilities import checkDir
from dateutil.relativedelta import relativedelta
from datetime import datetime, timedelta
from numpy import array, arange, amax, amin, where, linspace, zeros, float64, loadtxt, average
import xarray as xr
from pyshell.gridtools import TM5Grids
from pandas import DatetimeIndex, Timedelta, PeriodIndex, date_range
from pandas.tseries.frequencies import to_offset
from glob import glob


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
                    print("%20s %s-%s %10s : %20.6f Tg" % (cat, time_start.strftime("%Y/%m/%d"), time_end.strftime("%Y/%m/%d"), region, emission_coarse.sum()*1.0e-9))
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

                            
class CO2_Emissions(Emissions):
    def __init__(self, rcf, tracer='CO2', *args, **kwargs):
        Emissions.__init__(self, rcf, *args, **kwargs)
        self.tracer = tracer

        self.writeCycle = self.rcf.get('CO2.emission.dailycycle', 'bool')
        self.MolarMass = 44.00995 * 1.e-3  # to convert from moles of CO2 to kg CO2

        self.granularity = timedelta(hours=3)
        self.fossil_src_granularity = timedelta(hours=3)
        self.secs_per_granule = self.granularity.total_seconds()
        self.bio_src_granularity = timedelta(minutes=90)

#        self.consolidated_dir = self.rcf.get('consolidated.flux.output.dir')
#        self.consolidated_bio_dir = self.rcf.get('consolidated.bio.flux.dir', default=self.consolidated_dir)
#        self.consolidated_oce_dir = self.rcf.get('consolidated.oce.flux.dir', default=self.consolidated_dir)
#        self.consolidated_fire_dir = self.rcf.get('consolidated.fire.flux.dir', default=self.consolidated_dir)
#        self.consolidated_fos_dir = self.rcf.get('consolidated.fos.flux.dir', default=self.consolidated_dir)

#        self.sibcasa_dir = self.rcf.get('emission.CO2_nat.sibcasa.dir')
        self.output_dir = os.path.dirname(self.emission_file_name)
        self.fossil_sep_error = self.rcf.get('CO2.fossil fuel.sep_error', 'bool', default=False)
        self.bio_sep_error = self.rcf.get('CO2.terrestrial flux.sep_error', 'bool', default=True)

        # Some variables to prevent reading in the same file over and over
        self.current_files = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}
        self.current_emis = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}
        self.current_dates = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}

        self.warnings = []

    def readCT2019w_ocean(self, start_time, region, categ='ocean'):
        file_path = self.rcf.get('emissions.input.path')
        filename = os.path.join(file_path, 'ocean/CT2019.prior_oc.nc')
        if filename != self.current_files[categ]:
            if not os.path.exists(filename):
                raise RuntimeError('%s does not exist' % filename)
            with Dataset(filename, 'r') as fid :
                self.current_emis[categ] = fid.variables['oc'][:]
                self.current_dates[categ] = array([datetime(2000, 1, 1)+timedelta(x)-timedelta(7)/2 for x in fid.variables['time'][:]])
                self.current_files[categ] = filename

        if start_time > amax(self.current_dates[categ]):
            warn = '[WARNING] Ocean flux of 2018 will be used for %i' % start_time.year
            self.current_dates[categ] = self.current_dates[categ] + relativedelta(years=1)
            if warn not in self.warnings :
                print(warn)
                self.warnings.append(warn)
        return self.redistrib(categ, start_time, region, 'mol/m2/s')

    def readGFED(self, start_time, region, categ=None):
        file_path = self.rcf.get('emissions.input.path')
        filename = os.path.join(file_path, start_time.strftime('GFED/GFED_CO2_%Y%m%d_3hrly.nc'))
        with Dataset(filename, 'r') as ds :
            emis = ds['co2_emissions'][start_time.hour/3, :, :]
        lats_out = self.lat_grid[region]
        lons_out = self.lon_grid[region]
        rlat = int((lats_out[1]-lats_out[0])/0.25)
        rlon = int((lons_out[1]-lons_out[0])/0.25)
        nlat, nlon = emis.shape
        # Coarsen:
        emis = emis.reshape(nlat/rlat, rlat, nlon/rlon, rlon).sum(3).sum(1)
        # Crop (compare thesouth/east cell boundaries):
        lons_temp = arange(-180, 180, lons_out[1]-lons_out[0])
        lats_temp = arange(-90, 90, lats_out[1]-lats_out[0])
        emis = emis[[l in lats_out[:-1] for l in lats_temp], :]
        emis = emis[:, [l in lons_out[:-1] for l in lons_temp]]
        self.current_dates[categ] = array([start_time])
        return emis

    def readEDGAR(self, start_time, region, categ='fossil'):
        file_path = self.rcf.get('emissions.input.path')
        filename = os.path.join(file_path, start_time.strftime('fossil/EDGARv4.3_BP2019/EDGARv4.3_BP2019_emissions.co2.global.0.5x0.5.1hr.%Y%m.nc'))
        if start_time.year >= 2020:
            filename = os.path.join(file_path, start_time.strftime('fossil/EDGARv4.3_BP2019/EDGARv4.3_BP2019_emissions.co2.global.0.5x0.5.1hr.2019%m.nc'))
        
        # Load the file in memory, if it's not already there
        if filename != self.current_files[categ]:
            self.current_files[categ] = filename
            with Dataset(filename, 'r') as ds:
                self.current_emis[categ] = ds['emission'][:] * 1.e-6  # we want the emission in mol/m2.s
                self.current_dates[categ] = num2date(ds['time'][:], ds['time'].units)
                if start_time >= datetime(2020, 1, 1):
                    self.current_dates[categ] = self.current_dates[categ] + relativedelta(year=start_time.year)
                    
        # Calculate scaling factors based on carbon monitor, for dates > 2019
        ratio = 1.
        if start_time.year > 2019 :
            ratios = loadtxt(os.path.join(file_path, 'ratios_edgar.csv'), delimiter=',')
            times = [datetime(*t) for t in ratios[:, :3].astype(int)]
            ratio = ratios[times.index(datetime(start_time.year, start_time.month, start_time.day)), 4]

        return self.redistrib('fossil', start_time, region, 'mol/m2/s') * ratio

    def readSiBCASA2018_bio(self, start_time, region, categ='bio'):
        file_path = self.rcf.get('emissions.input.path')
        filename = os.path.join(file_path, start_time.strftime('biosphere/SiBCASA-eurocom2018/global_1x1/biofireparams_sibcasa_ei_%Y%m.nc'))
        field = 'nep'
        if start_time >= datetime(2018, 1, 1) :
            filename = os.path.join(file_path, start_time.strftime('biosphere/SiBCASA/global_1x1/biofireparams_sibcasa_ei_%m_clim.nc'))
            warn = '[WARNING] Using SiBCASA climatological prior after 2018!'
            if warn not in self.warnings :
                print(warn)
                self.warnings.append(warn)
            field = 'bio'
        if filename != self.current_files[categ]:
            if not os.path.exists(filename):
                raise RuntimeError('%s does not exist' % filename)
            with Dataset(filename, 'r') as ds:
                self.current_emis[categ] = ds[field][:]  # we want the emission in mol/m2.s
                self.current_dates[categ] = array([datetime(start_time.year, start_time.month, 1) + timedelta(hours=3*x) for x in ds['time']])
            self.current_files[categ] = filename
        return self.redistrib(categ, start_time, region, 'mol/m2/s')

    def readSiBCASA2018_bmb(self, start_time, region, categ='fire'):
        file_path = self.rcf.get('emissions.input.path')
        filename = os.path.join(file_path, start_time.strftime('fire/SiBCASA-eurocom2018/global_1x1/biofireparams_sibcasa_ei_%Y%m.nc'))
        if start_time >= datetime(2018, 1, 1) :
            filename = os.path.join(file_path, start_time.strftime('fire/SiBCASA/global_1x1/biofireparams_sibcasa_ei_%m_clim.nc'))
            warn = '[WARNING] Using SiBCASA cl imatological prior after 2018!'
            if warn not in self.warnings :
                print(warn)
                self.warnings.append(warn)
        if filename != self.current_files[categ]:
            if not os.path.exists(filename): 
                raise RuntimeError('%s does not exist'%filename)
            with Dataset(filename, 'r') as ds:
                self.current_emis[categ] = ds['bb'][:]  # we want the emission in mol/m2.s
                self.current_dates[categ] = array([datetime(start_time.year, start_time.month, 1) + timedelta(hours=3*x) for x in ds['time']])
            self.current_files[categ] = filename
        return self.redistrib('fire', start_time, region, 'mol/m2/s')

    def redistrib(self, cat, start_time, region, unit):
        try :
            idx = where(self.current_dates[cat] >= start_time)[0][0]-1
        except :
            try :
                if start_time.month == 2 and start_time.day == 29 :
                    # 29 february
                    idx = where(self.current_dates[cat] >= start_time-timedelta(1))[0][0]-1
                elif start_time > amax(self.current_dates[cat]) or start_time < amin(self.current_dates[cat]) :
                    msg = '%s %s' % (start_time.strftime("Requested date %d/%m/%Y not found in emission file"), self.current_files[cat])
                    raise RuntimeError(msg)
            except TypeError :
                import pdb; pdb.set_trace()
        if unit == 'mol/m2/s':
            emis = self.MolarMass*self.current_emis[cat][idx]  # kgCO2/m2/s
        else :
            raise
        nlat, nlon = emis.shape
        ip_lats = linspace(-90., 90., nlat+1)
        ip_lons = linspace(-180., 180., nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return redistrib_flux.regrid_flux_submul(ip_lats, ip_lons, op_lats, op_lons, emis, False)  # in Kg CO2/gridbox/sec

    def LoopThroughPeriods(self):
        """
        A single call to self.read*Flux gives the flux (in Kg CO2/cell/sec) during a three hour time window for a single
        category. This routine aggregates those fluxes into total flux over a given period, and writes out the daily
        cycle files. This is complicated if the granularity of emissions is different for different categories. For each
        category, we will write a different daily cycle file, and at the end of the subroutine we will assemble them
        together into one daily cycle file per day.
        """
        self.dailycycle_dates = set()
        self.Emission[self.tracer]['tf_bb_diurnal'] = None
        self.Emission[self.tracer]['emi_class'] = self.__class__.__name__
        print(self.tracer)
        emfine = {}
        for ireg, region in enumerate(self.zoom_regions_names):
            categories = self.Emission[region][self.tracer]['categories']
            emfine[region] = {}
            for icat, cat in enumerate(categories):
                emfine[region][cat] = {'time': [], 'emis': []}
                # What is the routine which is supposed to read fluxes for this category?
                flux_routine = self.rcf.get('%s.%s.routine' % (self.tracer, cat))
                readFlux = getattr(self, flux_routine)

                for time_index, time in enumerate(self.Emission[region][self.tracer][cat]['time_interval']['time_mid']):
                    time_start = self.Emission[region][self.tracer][cat]['time_interval']['time_start'][time_index]
                    time_end = self.Emission[region][self.tracer][cat]['time_interval']['time_end'][time_index]
                    sec_period = (time_end-time_start).total_seconds()
                    # How many three hour periods are there in this time window?
                    del_period = int(sec_period/self.secs_per_granule)
                    nlat, nlon = self.zoom_regions_lat[ireg][2], self.zoom_regions_lon[ireg][2]
                    emission_fine = zeros((del_period, nlat, nlon), float64)
                    emission_coarse = zeros((nlat, nlon), float64)  # this one will hold the total emission over this time period

                    i_time = 0
                    cur_time = time_start
                    while cur_time < time_end:
                        emission_fine[i_time] = readFlux(cur_time, region)
                        i_time += 1
                        cur_time += self.granularity

                    # what is the total emission?
                    emfine[region][cat]['time'].append(arange(time_start, time_end, self.granularity))
                    emfine[region][cat]['emis'].append(emission_fine+0.)
                    for emis in emission_fine:
                        emission_coarse += emis * self.secs_per_granule
                    # print total, remembering that emission_coarse is now in Kg CO2/gridbox/period
                    print("%20s %s-%s %10s : %20.6f Tg" % (cat, time_start.strftime("%Y/%m/%d"), time_end.strftime("%Y/%m/%d"), region, emission_coarse.sum()*1.0e-9))
                    emission_coarse = emission_coarse/sec_period  # bring back to per second
                    self.Emission[region][self.tracer][cat]['emission_data'][time_index, :, :] = emission_coarse

                    # Fill up the emission error structure, if required
                    if self.Emission[region][self.tracer][cat]['optimize']:
                        if cat == 'terrestrial flux' and self.bio_sep_error:
                            readGrossFlux = getattr(self, self.rcf.get('%s.terrestrial flux.gross.routine' % self.tracer))
                            self.CO2_nat_err_proxy = self.rcf.get('emission.CO2_nat.error.proxy')  # should be GPP or resp
                            self.CO2_nat_err_factor = self.rcf.get('emission.CO2_nat.error.proxy.factor', 'float')
                            self.Emission[region][self.tracer][cat]['emission_error'][time_index, :, :] = self.CO2_nat_err_factor * abs(readGrossFlux(time_start, time_end, region, self.CO2_nat_err_proxy))
                    
                        elif cat == 'fossil fuel' and self.fossil_sep_error:
                            readFossilSpread = getattr(self, self.rcf.get('%s.fossil fuel.spread.routine'))
                            self.fossil_spread_prefactor = self.rcf.get('emission.CO2.interprior.spread', 'float')
                            self.Emission[region][self.tracer][cat]['emission_error'][time_index, :, :] = readFossilSpread(time_start, time_end, region)

                    if self.writeCycle:
                        # take the average out of emission_fine
                        emission_fine = emission_fine - average(emission_fine, axis=0)
                        cur_time = time_start
                        i_time = 0
                        # How many self.granularity are there in one day?
                        time_steps_per_day = int(86400.0/self.secs_per_granule)
                        while cur_time < time_end:
                            t_idx_1 = i_time * time_steps_per_day
                            t_idx_2 = (i_time+1) * time_steps_per_day
                            file_name = os.path.join(self.output_dir, 'dailycycle_temp_%s_cat%i_%s.nc' % (region, icat, cur_time.strftime("%Y%m%d")))
                            self.dailycycle_dates.add(cur_time)
                            if os.path.exists(file_name):
                                os.remove(file_name)
                            checkDir(file_name)
                            fid = Dataset(file_name, 'w')
                            fid.createDimension('nt', time_steps_per_day)
                            fid.createDimension('nlat', emission_fine.shape[1])
                            fid.createDimension('nlon', emission_fine.shape[2])
                            v = fid.createVariable('emission_anomaly', float64, ('nt', 'nlat', 'nlon'))
                            v[:] = emission_fine[t_idx_1:t_idx_2]
                            fid.close()
                            i_time += 1
                            cur_time += timedelta(days=1)
        if self.writeCycle:
            self.consolidateDailyCycles()
        # Free up memory
        del self.current_files, self.current_dates, self.current_emis
