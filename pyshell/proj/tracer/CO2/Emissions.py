#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True

from pyshell.tmflex import rc
import re, os.path, shutil, subprocess

from numpy import *
from scipy import interpolate
from pyshell.base.helper.Utilities import checkDir
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4 import OrderedDict, date2num
from calendar import monthrange, isleap
from collections import defaultdict
from contextlib import contextmanager
from pyshell.tm5_utils import redistrib_flux

from pyshell.base.main.Emissions import Emissions, del_time, my_Dataset

@contextmanager
def nested_break():
    class NestedBreakException(Exception):
        pass
    try:
        yield NestedBreakException
    except NestedBreakException:
        pass

class CO2_Emission_Components(Emissions):
    """
    This class basically provides a bunch of routines to read CO2 emissions from different sources, such as Ivar's
    SiB-CASA GFED4, CT2013, etc. Previously all these routines used to be part of the CO2_Emissions class, but I realized
    that due to many tracer combinations needing these routines (radio_co2, CO2, COCO2), there was a lot of replication
    going on. So I've separated these routines from the emission assembler classes.
    """
    def __init__(self, *args, **kwargs):
        super(CO2_Emission_Components, self).__init__(*args, **kwargs)

    def readBioFlux_SiBCASA(self, start_time, region):
        # Read Ivar's SiB-CASA files for flux
        # The files do not have an extra day in February for leap years, so if start_time is Feb 29 in some leap
        # year, read the flux for Feb 28 by setting back start_time by one day.
        if isleap(start_time.year) and start_time.month == 2 and start_time.day == 29:
            start_time = start_time - timedelta(days=1)

        file_name = os.path.join(self.sibcasa_dir, start_time.strftime("biofireparams_sibcasa_ei_%Y%m.nc"))
        orig_date = datetime(start_time.year, start_time.month, 1, 0, 0, 0)
        if file_name != self.current_files['bio']:
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                self.current_emis['bio'] = fid.variables['nep'][:]
                self.current_dates['bio'] = [orig_date + timedelta(hours=3*i) for i in fid.variables['time'][:]]
            self.current_files['bio'] = file_name
        # find the index of the time dimension which would yield start_time
        time_in_month = start_time - orig_date
        idx = del_time(time_in_month)/del_time(self.granularity)
        if start_time != self.current_dates['bio'][idx]:
            sys.stderr.write('Index inside month does not match date/time\n')
            raise
        # Emission is in moles CO2/m^2/sec, convert to Kg CO2/gridbox/sec
        emis = self.MolarMass * self.current_emis['bio'][idx] # Kg CO2/m^2/sec
        nlat, nlon = emis.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,emis,False) # in Kg CO2/gridbox/sec
        return return_arr

    def readGrossFlux_SibCASA(self, start_time, end_time, region, var_name):
        """
        Reads the GPP or RH from Ivar's SiB-CASA monthly files, for use as a proxy for NEE error later
        """
        # what are the year/month tuples within this time window?
        time_windows = []
        cur_date = start_time
        t1 = start_time
        while cur_date < end_time:
            if cur_date.month != t1.month: # has crossed over a month boundary
                t2 = cur_date
                time_windows.append((t1,t2))
                t1 = t2
            cur_date += timedelta(days=1)
        time_windows.append((t1, end_time))
        # I already know that the emission is on a 180x360 grid, at 3 hourly resolution
        total_emis = zeros((180,360), float64)
        total_time = 0.0
        dt_per_step = 180.0 * 60.0 # number of seconds in 3 hours
        for t1, t2 in time_windows:
            file_name = os.path.join(self.sibcasa_dir, t1.strftime("biofireparams_sibcasa_ei_%Y%m.nc"))
            orig_date = datetime(t1.year, t1.month, 1, 0, 0, 0)
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                all_dates = [orig_date + timedelta(hours=3*i+1.5) for i in fid.variables['time'][:]]
                valid_idx = array([i for i,d in enumerate(all_dates) if t1 < d < t2])
                gross_flux = sum(dt_per_step * fid.variables[var_name][valid_idx, :, :], axis=0) # in mol C/m^2
            total_emis += self.MolarMass * gross_flux # in Kg CO2/m^2
            total_time += del_time(t2-t1).to_seconds()
        # Coarsen for the region
        ip_lats = linspace(-90., 90., 181)
        ip_lons = linspace(-180., 180., 361)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,total_emis,False) # in Kg CO2/gridbox
        return return_arr/total_time # in Kg CO2/gridbox/sec

    def readFireFlux_GFED4(self, start_time, region):
        # Read Ivar's GFED4 flux
        # The files do not have an extra day in February for leap years, so if start_time is Feb 29 in some leap
        # year, read the flux for Feb 28 by setting back start_time by one day.
        if isleap(start_time.year) and start_time.month == 2 and start_time.day == 29:
            start_time = start_time - timedelta(days=1)

        file_name = os.path.join(self.sibcasa_dir, start_time.strftime("biofireparams_sibcasa_ei_%Y%m.nc"))
        orig_date = datetime(start_time.year, start_time.month, 1, 0, 0, 0)
        if file_name != self.current_files['fire']:
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                self.current_emis['fire'] = fid.variables['bb'][:]
                self.current_dates['fire'] = [orig_date + timedelta(hours=3*i) for i in fid.variables['time'][:]]
            self.current_files['fire'] = file_name
        # find the index of the time dimension which would yield start_time
        time_in_month = start_time - orig_date
        idx = del_time(time_in_month)/del_time(self.granularity)
        if start_time != self.current_dates['fire'][idx]:
            sys.stderr.write('Index inside month does not match date/time\n')
            raise
        # Emission is in moles CO2/m^2/sec, convert to Kg CO2/gridbox/sec
        emis = self.MolarMass * self.current_emis['fire'][idx] # Kg CO2/m^2/sec
        nlat, nlon = emis.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,emis,False) # in Kg CO2/gridbox/sec
        return return_arr

    def readOceFlux(self, start_time, region):
        file_name = os.path.join(self.consolidated_oce_dir, start_time.strftime("flux_%Y%m.nc"))
        if file_name != self.current_files['ocean']:
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                self.current_emis['ocean'] = fid.variables['ocn_flux_opt'][:]
                self.current_dates['ocean'] = [datetime(*d) for d in fid.variables['date_components'][:]]
            self.current_files['ocean'] = file_name
        # find the index of the time dimension which would yield start_time
        time_in_month = start_time - datetime(start_time.year, start_time.month, 1, 0, 0, 0)
        idx = del_time(time_in_month)/del_time(self.granularity)
        if start_time + self.granularity/2 != self.current_dates['ocean'][idx]:
            raise RuntimeError('Index inside month does not match date/time')
        emis = self.MolarMass * self.current_emis['ocean'][idx] # Kg CO2/m^2/sec
        # now convert to Kg CO2/gridbox/sec
        nlat, nlon = emis.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,emis,False) # in Kg CO2/gridbox/sec
        return return_arr

    def readBioFlux_CASA(self, start_time, region):
        file_name = os.path.join(self.consolidated_bio_dir, start_time.strftime("flux_%Y%m.nc"))
        if file_name != self.current_files['bio']:
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                self.current_emis['bio'] = fid.variables['bio_flux_opt'][:]
                self.current_dates['bio'] = [datetime(*d) for d in fid.variables['date_components'][:]]
            self.current_files['bio'] = file_name
        # find the index of the time dimension which would yield start_time
        time_in_month = start_time - datetime(start_time.year, start_time.month, 1, 0, 0, 0)
        idx = del_time(time_in_month)/del_time(self.granularity)
        if start_time + self.granularity/2 != self.current_dates['bio'][idx]:
            sys.stderr.write('Index inside month does not match date/time\n')
            raise
        emis = self.MolarMass * self.current_emis['bio'][idx] # Kg CO2/m^2/sec
        # now convert to Kg CO2/gridbox/sec
        nlat, nlon = emis.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,emis,False) # in Kg CO2/gridbox/sec
        return return_arr

    def readGrossFlux_CASA(self, start_time, end_time, region, var_name):
        """
        Reads the GPP or RH from CASA monthly files, for use as a proxy for NEE error later
        """
        # what are the year/month tuples within this time window?
        time_windows = []
        cur_date = start_time
        t1 = start_time
        while cur_date < end_time:
            if cur_date.month != t1.month: # has crossed over a month boundary
                t2 = cur_date
                time_windows.append((t1,t2))
                t1 = t2
            cur_date += timedelta(days=1)
        time_windows.append((t1, end_time))
        # I already know that the emission is on a 180x360 grid, at 90 minute resolution
        total_emis = zeros((180,360), float64)
        total_time = 0.0
        dt_per_step = 90.0 * 60.0 # number of seconds in 90 minutes
        for t1, t2 in time_windows:
            file_name = t1.strftime(self.bio_file_pattern)
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                all_dates = [datetime(1970,1,1) + timedelta(days=d) for d in fid.variables['date'][:]]
                valid_idx = array([i for i,d in enumerate(all_dates) if t1 < d < t2])
                gross_flux = sum(dt_per_step * fid.variables[var_name][valid_idx, :, :], axis=0) # in mol C/m^2
            total_emis += self.MolarMass * gross_flux # in Kg CO2/m^2
            total_time += del_time(t2-t1).to_seconds()
        # Coarsen for the region
        ip_lats = linspace(-90., 90., 181)
        ip_lons = linspace(-180., 180., 361)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,total_emis,False) # in Kg CO2/gridbox
        return return_arr/total_time # in Kg CO2/gridbox/sec

    def readFireFlux_GFED3(self, start_time, region):
        # GFED3.1 fire fluxes are given in 3-hour windows, and our granularity is also 3 hours,
        # so let's cut the charade of taking time intervals, splicing, etc.
        monthly_file = os.path.join(self.fire_flux_input_dir, start_time.strftime(self.fire_file_monthly))
        daily_file   = os.path.join(self.fire_flux_input_dir, start_time.strftime(self.fire_file_daily))
        hourly_file  = os.path.join(self.fire_flux_input_dir, start_time.strftime(self.fire_file_hourly))
        if monthly_file != self.current_files['fire']['monthly']:
            self.current_emis['fire']['monthly'] = flipud(loadtxt(monthly_file))/1000.0 # Kg CO2/m^2/month
            self.current_files['fire']['monthly'] = monthly_file
        if daily_file != self.current_files['fire']['daily']:
            with my_Dataset(daily_file, 'r') as fid:
                self.current_emis['fire']['daily'] = fid.variables['Fraction_of_Emissions'][:]
            self.current_files['fire']['daily'] = daily_file
        if hourly_file != self.current_files['fire']['hourly']:
            with my_Dataset(hourly_file, 'r') as fid:
                self.current_emis['fire']['hourly'] = fid.variables['Fraction_of_Emissions'][:]
            self.current_files['fire']['hourly'] = hourly_file
        # what is the time index within the day?
        day_idx = del_time(start_time - datetime(start_time.year, start_time.month, start_time.day, 0, 0, 0))\
            /del_time(self.granularity)
        fire_flux = self.current_emis['fire']['monthly'] * self.current_emis['fire']['daily'] \
            * self.current_emis['fire']['hourly'][day_idx] # Kg CO2/m^2/3 hours
        fire_flux /= self.secs_per_granule # Kg CO2/m^2/sec
        # now coarsen over regions
        nlat, nlon = fire_flux.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,fire_flux,False) # in Kg CO2/gridbox/sec
        return return_arr

    def readFireFlux_CT2013(self, start_time, region):
        """
        CT2013 uses GFED 3.1 fires. However, it uses the GFED model on CASA carbon pools. Which means that even though the
        burnt areas are going to be identical with the official GFED 3.1 release, the actual CO2 emissions might be
        different. This routine reads GFED 3.1 fire fluxes "as interpreted by CASA GFED". Note that CASA GFED only gives
        monthly fires, so there is no need to read in daily and three-hourly factors.
        """
        file_name = os.path.join(self.consolidated_fire_dir, start_time.strftime("flux_%Y%m.nc"))
        if file_name != self.current_files['fire']['monthly']:
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                self.current_emis['fire']['monthly'] = fid.variables['fire_flux_imp'][:]
                self.current_dates['fire']['monthly'] = [datetime(*d) for d in fid.variables['date_components'][:]]
            self.current_files['fire']['monthly'] = file_name
        # find the index of the time dimension which would yield start_time
        time_in_month = start_time - datetime(start_time.year, start_time.month, 1, 0, 0, 0)
        idx = del_time(time_in_month)/del_time(self.granularity)
        if start_time + self.granularity/2 != self.current_dates['fire']['monthly'][idx]:
            raise RuntimeError('Index inside month does not match date/time')
        emis = self.MolarMass * self.current_emis['fire']['monthly'][idx] # Kg CO2/m^2/sec
        # now convert to Kg CO2/gridbox/sec
        nlat, nlon = emis.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,emis,False) # in Kg CO2/gridbox/sec
        return return_arr

    def readFossil_CT(self, start_time, region):
        # on Zeus, self.consolidated_dir is /scratch2/portfolios/BMC/co2/Sourish.Basu/tm5_inputs/fluxes/CO2/CT2013/unopt/ei-b3oifmu
        file_name = os.path.join(self.consolidated_fos_dir, start_time.strftime("flux_%Y%m.nc"))
        if file_name != self.current_files['fossil']:
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                self.current_emis['fossil'] = fid.variables['fossil_flux_imp'][:]
                self.current_dates['fossil'] = [datetime(*d) for d in fid.variables['date_components'][:]]
            self.current_files['fossil'] = file_name
        # find the index of the time dimension which would yield start_time
        time_in_month = start_time - datetime(start_time.year, start_time.month, 1, 0, 0, 0)
        idx = del_time(time_in_month)/del_time(self.fossil_src_granularity)
        if start_time + self.fossil_src_granularity/2 != self.current_dates['fossil'][idx]:
            sys.stderr.write('Index inside month does not match date/time\n')
            raise
        emis = self.MolarMass * self.current_emis['fossil'][idx] # Kg CO2/m^2/sec
        # now convert to Kg CO2/gridbox/sec
        nlat, nlon = emis.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,emis,False) # in Kg CO2/gridbox/sec
        if self.writeCycle:
            hr_idx = del_time(start_time - datetime(*start_time.timetuple()[:3]))/del_time(self.fossil_src_granularity)
            if self.apply_weekly_fac:
                day_idx = start_time.timetuple().tm_wday
                return_arr = return_arr * self.daily_factors[region][hr_idx] * self.weekly_factors[region][day_idx]
            else:
                return_arr = return_arr * self.daily_factors[region][hr_idx]
        return return_arr

    def readFossil_Miller(self, start_time, region):
        fossil_fuel_dir = self.rcf.get('emissions.miller.output.dir')
        fossil_fuel_var = self.rcf.get('emissions.miller.fossil.var_name')

        # On zeus, fossil_fuel_dir is /scratch1/portfolios/BMC/co2/Sourish.Basu/tm5_inputs/fluxes/CO2/Miller/ct13_c14basev7c
        file_name = os.path.join(fossil_fuel_dir, start_time.strftime('%Y/%m/fossil_fluxes_1x1.nc'))
        if file_name != self.current_files['fossil']:
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                self.current_emis['fossil'] = fid.variables[fossil_fuel_var][:]
            self.current_files['fossil'] = file_name
        # Which day does this start_time correspond to?
        i_day = start_time.day - 1
        emis = self.MolarMass * self.current_emis['fossil'][i_day] # Kg CO2/m^2/sec
        # now convert to Kg CO2/gridbox/sec
        nlat, nlon = emis.shape
        ip_lats = linspace(-90.,90.,nlat+1)
        ip_lons = linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,emis,False) # in Kg CO2/gridbox/sec
        if self.writeCycle:
            hr_idx = del_time(start_time - datetime(*start_time.timetuple()[:3]))/del_time(self.fossil_src_granularity)
            if self.apply_weekly_fac:
                day_idx = start_time.timetuple().tm_wday
                return_arr = return_arr * self.daily_factors[region][hr_idx] * self.weekly_factors[region][day_idx]
            else:
                return_arr = return_arr * self.daily_factors[region][hr_idx]
        return return_arr

    def readFossilSpread(self, start_time, end_time, region):
        """
        For a given zoom region, read the spread in fossil fuel fluxes between Miller FF, Odiac, and Miller FF with
        Vulcan pattern. This will serve as a proxy for grid-scale FF flux error.
        """
        input_dir = self.rcf.get('emissions.miller.output.dir') # this is where the monthly fluxes are
        # On zeus, this is /scratch2/portfolios/BMC/co2/Sourish.Basu/tm5_inputs/emissions_CO2/Miller
        # create the output array, of dimension 3 x nlat x nlon
        nlat = self.lat_grid[region].shape[0] - 1
        nlon = self.lon_grid[region].shape[0] - 1
        flux_arr = zeros((3, nlat, nlon), float64)
        # fluxes are stored in one file per month, with one flux array per day
        year_month_dict = defaultdict(list)
        cur_date = start_time
        while cur_date < end_time:
            year, month = cur_date.year, cur_date.month
            year_month_dict[(year,month)].append(cur_date.day-1)
            cur_date += timedelta(days=1)
        sec_per_day = 86400.0
        for (year,month), indices in year_month_dict.items():
            file_name = os.path.join(input_dir, "%04i"%year, "%02i"%month, "fossil_fluxes_1x1.nc")
            if not os.path.exists(file_name):
                raise RuntimeError('%s does not exist'%file_name)
            with my_Dataset(file_name, 'r') as fid:
                nlat = len(fid.dimensions['lat'])
                nlon = len(fid.dimensions['lon'])
                flux_a = fid.variables['flux_a'][:]
                flux_b = fid.variables['flux_b'][:]
                flux_c = fid.variables['flux_c'][:]
            # calculate total flux over the overlap period between this month and (start_time, end_time)
            tot_flux_overlap = zeros((3, nlat, nlon), float64)
            for i in indices:
                tot_flux_overlap[0] += (flux_a[i] * sec_per_day)
                tot_flux_overlap[1] += (flux_b[i] * sec_per_day)
                tot_flux_overlap[2] += (flux_c[i] * sec_per_day)
            # now tot_flux_overlap is in gC/m^2, globally
            tot_flux_overlap = self.MolarMass * tot_flux_overlap # now in Kg CO2/m^2
            ip_lats = linspace(-90.,90.,nlat+1)
            ip_lons = linspace(-180.,180.,nlon+1)
            op_lats = self.lat_grid[region]
            op_lons = self.lon_grid[region]
            for i in range(3):
                temp_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,tot_flux_overlap[i],False) # in Kg CO2/gridbox
                flux_arr[i] += temp_arr
            del flux_a, flux_b, flux_c, tot_flux_overlap, temp_arr
        # flux_arr contains the three fossil fuel fluxes over (start_time, end_time), in Kg CO2/gridbox
        total_seconds = del_time(end_time-start_time).to_seconds()
        flux_arr = flux_arr/total_seconds # now in Kg CO2/gridbox/sec
        flux_sigma = self.fossil_spread_prefactor * ptp(flux_arr, axis=0)
        del flux_arr
        return flux_sigma

    def readFFScaling(self):
        daily_scale_file = self.rcf.get('CO2_ff.daily.scaling.factor')
        weekly_scale_file = self.rcf.get('CO2_ff.weekly.scaling.factor')
        with my_Dataset(daily_scale_file, 'r') as fid:
            daily_factors = fid.variables['diurnal_scale_factors'][:]
        with my_Dataset(weekly_scale_file, 'r') as fid:
            weekly_factors = fid.variables['weekly_scale_factors'][:]
        # need to coarsen them down to zoom regions, first the daily factors
        nlat, nlon, nhours = daily_factors.shape
        src_lat_grid = linspace(-90., 90., nlat+1)
        src_lon_grid = linspace(-180., 180., nlon+1)
        dlat_src = 180.0/nlat
        dlon_src = 360.0/nlon
        self.daily_factors = {}
        for region in self.zoom_regions_names:
            dlat_dst = (self.zoom_info[region]['yend']-self.zoom_info[region]['ybeg'])/self.zoom_info[region]['jm']
            dlon_dst = (self.zoom_info[region]['xend']-self.zoom_info[region]['xbeg'])/self.zoom_info[region]['im']
            lat_mul = int(dlat_dst/dlat_src)
            lon_mul = int(dlon_dst/dlon_src)
            i_beg = src_lon_grid.searchsorted(self.zoom_info[region]['xbeg'])
            i_end = src_lon_grid.searchsorted(self.zoom_info[region]['xend'])
            j_beg = src_lat_grid.searchsorted(self.zoom_info[region]['ybeg'])
            j_end = src_lat_grid.searchsorted(self.zoom_info[region]['yend'])
            daily_fac_subarr = daily_factors[j_beg:j_end, i_beg:i_end, :]
            daily_fac_reg = zeros((nhours, self.zoom_info[region]['jm'], self.zoom_info[region]['im']), float64)
            for hr in range(nhours):
                for j in range(self.zoom_info[region]['jm']):
                    for i in range(self.zoom_info[region]['im']):
                        daily_fac_reg[hr,j,i] = \
                            sum(daily_fac_subarr[j*lat_mul:(j+1)*lat_mul,i*lon_mul:(i+1)*lon_mul,hr])
            # The daily factors are still over 24 hours, need to reduce them to self.fossil_src_granularity
            daily_fac_res = del_time(timedelta(hours=24)/nhours)
            hr_mul = del_time(self.fossil_src_granularity)/daily_fac_res
            num_steps_per_day = del_time(timedelta(hours=24))/del_time(self.fossil_src_granularity)
            self.daily_factors[region] = zeros((num_steps_per_day,daily_fac_reg.shape[1],daily_fac_reg.shape[2]), float64)
            for hr in range(num_steps_per_day):
                self.daily_factors[region][hr,:,:] = sum(daily_fac_reg[hr*hr_mul:(hr+1)*hr_mul,:,:], axis=0)
            # Now normalize so that if the emission at each time step is the same, the factors are all one
            sum_of_factors = sum(self.daily_factors[region], axis=0)
            for hr in range(num_steps_per_day):
                self.daily_factors[region][hr,:,:] = num_steps_per_day * self.daily_factors[region][hr,:,:]/sum_of_factors
        # now the weekly factors
        # WARNING : Ray Nassar told me that the time dimension of the weekly scaling factors does not reflect
        # WARNING : the UTC weekday. So, the first temporal element for each pixel is the scaling to be applied
        # WARNING : on Monday for that pixel. Currently, I am checking the day of week of the UTC date, and
        # WARNING : applying the scaling factor for that day. This needs to be changed. But there is no simple
        # WARNING : way of doing that. The problem is that when it's 00:30 UTC Monday, over N America the scaling
        # WARNING : factor for Sunday should be applied. However, when it's 22:00 UTC Monday, over N America the
        # WARNING : scaling factor for Monday should be applied. So there needs to be one scaling factor array
        # WARNING : per time step within a day. This is too much effort for something that Ray says doesn't matter
        # WARNING : very much, i.e., the day of week variation is much smaller than the hour of day variation.
        # WARNING : So for now, I'll just not use the day of week scaling.
        nlat, nlon, ndays = weekly_factors.shape
        src_lat_grid = linspace(-90., 90., nlat+1)
        src_lon_grid = linspace(-180., 180., nlon+1)
        dlat_src = 180.0/nlat
        dlon_src = 360.0/nlon
        self.weekly_factors = {}
        for region in self.zoom_regions_names:
            dlat_dst = (self.zoom_info[region]['yend']-self.zoom_info[region]['ybeg'])/self.zoom_info[region]['jm']
            dlon_dst = (self.zoom_info[region]['xend']-self.zoom_info[region]['xbeg'])/self.zoom_info[region]['im']
            lat_mul = int(dlat_dst/dlat_src)
            lon_mul = int(dlon_dst/dlon_src)
            i_beg = src_lon_grid.searchsorted(self.zoom_info[region]['xbeg'])
            i_end = src_lon_grid.searchsorted(self.zoom_info[region]['xend'])
            j_beg = src_lat_grid.searchsorted(self.zoom_info[region]['ybeg'])
            j_end = src_lat_grid.searchsorted(self.zoom_info[region]['yend'])
            weekly_fac_subarr = weekly_factors[j_beg:j_end, i_beg:i_end, :]
            weekly_fac_reg = zeros((ndays, self.zoom_info[region]['jm'], self.zoom_info[region]['im']), float64)
            # Both Ray Nassar and Python consider Monday to be the first day of the week, thank God!
            for day in range(ndays):
                for j in range(self.zoom_info[region]['jm']):
                    for i in range(self.zoom_info[region]['im']):
                        weekly_fac_reg[day,j,i] = \
                            sum(weekly_fac_subarr[j*lat_mul:(j+1)*lat_mul,i*lon_mul:(i+1)*lon_mul,day])
            self.weekly_factors[region] = weekly_fac_reg
            # Now normalize so that if the emission in each day of week is the same, the factors are all one
            sum_of_factors = sum(self.weekly_factors[region], axis=0)
            for day in range(ndays):
                self.weekly_factors[region][day,:,:] = ndays * self.weekly_factors[region][day,:,:]/sum_of_factors

class CO2_Emissions(CO2_Emission_Components):
    """
    This class assembles the biospheric and oceanic fluxes that go into making CO2_nat. The biospheric fluxes
    come from CASA GFED 3.1, and the fire fluxes from GFED 3.1, by default.
    """
    def __init__(self, StartTime, EndTime, subdir_tag=''):
        super(CO2_Emissions, self).__init__(StartTime, EndTime, subdir_tag)
        self.tracer = 'CO2'

        self.bio_flux_input_dir = self.rcf.get('emission.CO2_nat.bio.input.dir')
        self.bio_flux_input_file = self.rcf.get('emission.CO2_nat.bio.input.file')
        self.bio_file_pattern = os.path.join(self.bio_flux_input_dir, self.bio_flux_input_file)

        self.fire_flux_input_dir = self.rcf.get('emission.CO2_nat.fire.input.dir')
        self.fire_file_monthly = self.rcf.get('emission.CO2_nat.fire.monthly.file')
        self.fire_file_daily = self.rcf.get('emission.CO2_nat.fire.daily.file')
        self.fire_file_hourly = self.rcf.get('emission.CO2_nat.fire.hourly.file')

        self.ocean_flux_input_dir = self.rcf.get('emission.CO2_nat.ocean.input.dir')
        self.fossil_flux_input_dir = self.rcf.get('emission.CO2_ff.fossil.input.dir')

        self.writeCycle = self.rcf.get('CO2.emission.dailycycle', 'bool')
        self.MolarMass = 44.00995 * 1.e-3 # to convert from moles of CO2 to kg CO2

        self.granularity = timedelta(hours=3)
        self.fossil_src_granularity = timedelta(hours=3)
        self.secs_per_granule = del_time(self.granularity).to_seconds()
        self.bio_src_granularity = timedelta(minutes=90)

        self.consolidated_dir = self.rcf.get('consolidated.flux.output.dir')
        self.consolidated_bio_dir = self.rcf.get('consolidated.bio.flux.dir', default=self.consolidated_dir)
        self.consolidated_oce_dir = self.rcf.get('consolidated.oce.flux.dir', default=self.consolidated_dir)
        self.consolidated_fire_dir = self.rcf.get('consolidated.fire.flux.dir', default=self.consolidated_dir)
        self.consolidated_fos_dir = self.rcf.get('consolidated.fos.flux.dir', default=self.consolidated_dir)

        self.sibcasa_dir = self.rcf.get('emission.CO2_nat.sibcasa.dir')
        self.output_dir = os.path.dirname(self.emission_file_name)
        self.fossil_sep_error = self.rcf.get('CO2.fossil fuel.sep_error', 'bool', default=False)
        self.bio_sep_error = self.rcf.get('CO2.terrestrial flux.sep_error', 'bool', default=True)

        # Some variables to prevent reading in the same file over and over
        self.current_files = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}
        self.current_emis = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}
        self.current_dates = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}

        # We always apply Ray Nassar's daily scaling factors, but perhaps not the weekly ones
        self.apply_weekly_fac = self.rcf.get('CO2_ff.emission.weeklycycle', 'bool', default=False)
        self.readFFScaling()

        # After a certain year, the CT2013 fluxes do not exist any more, and we need to read in climatologically extended fluxes
        self.last_valid_ct13_year = 2012

    def consolidateFluxFiles(self, start_year=2000, end_year=2012):
        """
        Andy has supplied three hourly flux files, such as
        $SCRATCH/tm5_inputs/emissions_CO2/CT2013/unopt/ei-b3oifmu/2012/11/flux1x1_201211122100_201211130000.nc
        and reading them while assembling fluxes is very slow (think of the number of files that need to be
        opened and closed). So this routine assembles all those three hourly flux files into monthly flux files
        using ncrcat. The idea is to generate a list of relevant files in the correct order and feed them to
        ncrcat with an output file name.

        The following months have problems that Andy talked about:

        09/2003
        10/2003
        11/2003
        12/2003
        01/2004

        Update 14/07/2014: He has fixed those problems, and I have updated my links and consolidated those months.
        """
        for year in range(start_year, end_year+1):
            for month in range(1,13):
                self.consolidateMonth(year, month)

    def consolidateMonth(self, year, month):
        # Helper function to consolidateFluxFiles
        dst_dir = self.rcf.get('consolidated.flux.output.dir')
        fmt = '%Y%m%d%H%M'
        dt = timedelta(hours=3)
        num_days = monthrange(year, month)[1]
        file_list = []
        for day in range(1, num_days+1):
            for hr in range(8):
                cur_date = datetime(year, month, day) + hr*dt
                next_date = cur_date + dt
                src_dir = cur_date.strftime(self.ocean_flux_input_dir)
                file_name = '_'.join(['flux1x1', cur_date.strftime(fmt), next_date.strftime(fmt)+'.nc'])
                file_name = os.path.join(src_dir, file_name)
                if not os.path.exists(file_name):
                    sys.stderr.write('File %s does not exist\n'%os.path.basename(file_name))
                    return
                file_list.append(file_name)
        dst_file = os.path.join(dst_dir, 'flux_%04i%02i.nc'%(year, month))
        if os.path.exists(dst_file):
            os.remove(dst_file)
        sbp_comm = ['ncrcat', '-h'] + file_list + [dst_file]
        print 'Merging %04i/%02i'%(year, month)
        subprocess.check_call(sbp_comm)

    def readBioFireFlux(self, start_time, region):
        bio_flux_routine = getattr(self, self.rcf.get('%s.biosphere flux.routine'%self.tracer))
        bio_arr = bio_flux_routine(start_time, region)

        fire_flux_routine = getattr(self, self.rcf.get('%s.fire flux.routine'%self.tracer))
        fire_arr = fire_flux_routine(start_time, region)

        return bio_arr+fire_arr

    def LoopThroughPeriods(self):
        """
        A single call to self.read*Flux gives the flux (in Kg CO2/cell/sec) during a three hour time window for a single
        category. This routine aggregates those fluxes into total flux over a given period, and writes out the daily
        cycle files. This is complicated if the granularity of emissions is different for different categories. For each
        category, we will write a different daily cycle file, and at the end of the subroutine we will assemble them
        together into one daily cycle file per day.
        """
        self.dailycycle_dates = set()
        self.Emission[self.tracer]['tf_bb_diurnal'] =  None
        for ireg, region in enumerate(self.zoom_regions_names):
            categories = self.Emission[region][self.tracer]['categories']
            for icat, cat in enumerate(categories):
                # What is the routine which is supposed to read fluxes for this category?
                flux_routine = self.rcf.get('%s.%s.routine'%(self.tracer, cat))
                readFlux = getattr(self, flux_routine)

                for time_index, time in enumerate(self.Emission[region][self.tracer][cat]['time_interval']['time_mid']):
                    time_start = self.Emission[region][self.tracer][cat]['time_interval']['time_start'][time_index]
                    time_end   = self.Emission[region][self.tracer][cat]['time_interval']['time_end'][time_index]
                    sec_period = del_time(time_end - time_start).to_seconds()
                    # How many three hour periods are there in this time window?
                    del_period = int(sec_period/self.secs_per_granule)
                    nlat, nlon = self.zoom_regions_lat[ireg][2], self.zoom_regions_lon[ireg][2]
                    emission_fine = zeros((del_period, nlat, nlon), float64)
                    emission_coarse = zeros((nlat, nlon), float64) # this one will hold the total emission over this time period

                    i_time = 0
                    cur_time = time_start
                    while cur_time < time_end:
                        emission_fine[i_time] = readFlux(cur_time, region)
                        i_time += 1
                        cur_time += self.granularity

                    # what is the total emission?
                    for emis in emission_fine:
                        emission_coarse += emis * self.secs_per_granule
                    # print total, remembering that emission_coarse is now in Kg CO2/gridbox/period
                    print "%20s %s-%s %10s : %20.6f Tg"%(cat, time_start.strftime("%Y/%m/%d"), \
                        time_end.strftime("%Y/%m/%d"), region, emission_coarse.sum()*1.0e-9)
                    emission_coarse = emission_coarse/sec_period # bring back to per second
                    self.Emission[region][self.tracer][cat]['emission_data'][time_index,:,:] = emission_coarse

                    # Fill up the emission error structure, if required
                    if self.Emission[region][self.tracer][cat]['optimize']:
                        if cat == 'terrestrial flux' and self.bio_sep_error:
                            readGrossFlux = getattr(self, self.rcf.get('%s.terrestrial flux.gross.routine'%self.tracer))
                            self.CO2_nat_err_proxy = self.rcf.get('emission.CO2_nat.error.proxy') # should be be GPP or resp
                            self.CO2_nat_err_factor = self.rcf.get('emission.CO2_nat.error.proxy.factor', 'float')
                            self.Emission[region][self.tracer][cat]['emission_error'][time_index,:,:] = self.CO2_nat_err_factor * \
                                abs(readGrossFlux(time_start, time_end, region, self.CO2_nat_err_proxy))

                        elif cat == 'fossil fuel' and self.fossil_sep_error:
                            readFossilSpread = getattr(self, self.rcf.get('%s.fossil fuel.spread.routine'))
                            self.fossil_spread_prefactor = self.rcf.get('emission.CO2.interprior.spread', 'float')
                            self.Emission[region][self.tracer][cat]['emission_error'][time_index,:,:] = \
                                readFossilSpread(time_start, time_end, region)

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
                            file_name = os.path.join(self.output_dir, \
                                'dailycycle_temp_%s_cat%i_%s.nc'%(region,icat,cur_time.strftime("%Y%m%d")))
                            self.dailycycle_dates.add(cur_time)
                            if os.path.exists(file_name):
                                os.remove(file_name)
                            checkDir(file_name)
                            fid = my_Dataset(file_name, 'w')
                            fid.createDimension('nt', time_steps_per_day)
                            fid.createDimension('nlat', emission_fine.shape[1])
                            fid.createDimension('nlon', emission_fine.shape[2])
                            v = fid.createVariable('emission_anomaly', float64, ('nt','nlat','nlon'))
                            v[:] = emission_fine[t_idx_1:t_idx_2]
                            fid.close()
                            i_time += 1
                            cur_time += timedelta(days=1)
        if self.writeCycle:
            self.consolidateDailyCycles()
        # Free up memory
        del self.current_files, self.current_dates, self.current_emis

    #def consolidateDailyCycles(self):
        #"""
        #Takes the individual dailycycle_temp_* files written out by LoopThroughPeriods and creates one
        #daily cycle file per day
        #"""
        #time_steps_per_day = int(86400.0/self.secs_per_granule)
        #for cur_time in self.dailycycle_dates:
            #output_filename = os.path.join(self.rcf.get('dailycycle.folder'), cur_time.strftime("%Y/%m"), \
                #self.rcf.get('%s.dailycycle.prefix'%self.tracer)+cur_time.strftime("%Y%m%d.nc4"))
            #if os.path.exists(output_filename):
                ## just delete the temp files
                #for ireg, region in enumerate(self.zoom_regions_names):
                    #categories = self.Emission[region][self.tracer]['categories']
                    #for icat, cat in enumerate(categories):
                        #input_filename = os.path.join(self.output_dir, 'dailycycle_temp_%s_cat%i_%s.nc'%\
                            #(region,icat,cur_time.strftime("%Y%m%d")))
                        #os.remove(input_filename)
            #else:
                #checkDir(output_filename)
                #ofid = my_Dataset(output_filename, 'w')
                #for ireg, region in enumerate(self.zoom_regions_names):
                    #ogid_reg = ofid.createGroup(region)
                    #ogid_reg.createDimension('latitude', self.zoom_regions_lat[ireg][2])
                    #ogid_reg.createDimension('longitude', self.zoom_regions_lon[ireg][2])
                    #categories = self.Emission[region][self.tracer]['categories']
                    #for icat, cat in enumerate(categories):
                        #ogid_cat = ogid_reg.createGroup(cat)
                        #ogid_cat.createDimension('timesteps', time_steps_per_day)
                        #input_filename = os.path.join(self.output_dir, 'dailycycle_temp_%s_cat%i_%s.nc'%\
                            #(region,icat,cur_time.strftime("%Y%m%d")))
                        #with my_Dataset(input_filename, 'r') as ifid:
                            #emis = ifid.variables['emission_anomaly'][:]
                        #v = ogid_cat.createVariable('emission_anomaly', float64, ('timesteps','latitude','longitude'))
                        #v[:] = emis
                        #os.remove(input_filename)
                #ofid.close()

class CO2_daynight(CO2_Emissions, Emissions):
    """
    The fact that the "true" daily cycle is usually different from the prior daily cycle has been bothering me for a while.
    One idea is to optimize respiration and photosynthesis separately, and impose the daily cycle as a set of scaling factors
    on the respective average fluxes. So we will now have four categories of CO2 fluxes, oceanic, terrestrial respiration (which
    should include fires), terrestrial productivity, and fossil fuels.
    """
    def __init__(self, StartTime, EndTime):
        Emissions.__init__(self, StartTime, EndTime)
        self.tracer = 'CO2'

        self.bio_flux_input_dir = self.rcf.get('emission.CO2_nat.bio.input.dir')
        self.bio_flux_input_file = self.rcf.get('emission.CO2_nat.bio.input.file')
        self.bio_file_pattern = os.path.join(self.bio_flux_input_dir, self.bio_flux_input_file)

        # The following block is only relevant for reading the gross biosphere flux, in case we
        # use CASA GFED 3.1 as our prior. For SibCASA, self.sibcasa_dir is used.
        self.prior_flux_input_dir = self.rcf.get('emission.CO2_nat.prior.input.dir')
        self.prior_flux_input_file = self.rcf.get('emission.CO2_nat.prior.input.file')

        self.fire_flux_input_dir = self.rcf.get('emission.CO2_nat.fire.input.dir')
        self.fire_file_monthly = self.rcf.get('emission.CO2_nat.fire.monthly.file')
        self.fire_file_daily = self.rcf.get('emission.CO2_nat.fire.daily.file')
        self.fire_file_hourly = self.rcf.get('emission.CO2_nat.fire.hourly.file')

        self.ocean_flux_input_dir = self.rcf.get('emission.CO2_nat.ocean.input.dir')
        self.fossil_flux_input_dir = self.rcf.get('emission.CO2_ff.fossil.input.dir')

        self.writeCycle = self.rcf.get('CO2.emission.dailycycle', 'bool', default=True)
        self.MolarMass = 44.00995 * 1.e-3 # to convert from moles of CO2 to kg CO2

        self.granularity = timedelta(hours=3)
        self.fossil_src_granularity = timedelta(hours=3)
        self.secs_per_granule = del_time(self.granularity).to_seconds()
        self.bio_src_granularity = timedelta(minutes=90)

        self.consolidated_dir = self.rcf.get('consolidated.flux.output.dir')
        self.sibcasa_dir = self.rcf.get('emission.CO2_nat.sibcasa.dir')
        self.output_dir = os.path.dirname(self.emission_file_name)

        # Some variables to prevent reading in the same file over and over
        self.current_files = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}
        self.current_emis = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}
        self.current_dates = {'ocean': None, 'bio': None, 'fire': {'monthly': None, 'daily': None, 'hourly': None}, 'fossil': None}

        # We always apply Ray Nassar's daily scaling factors, but perhaps not the weekly ones
        self.apply_weekly_fac = self.rcf.get('CO2_ff.emission.weeklycycle', 'bool', default=False)
        self.readFFScaling()

        # After a certain year, the CT2013 fluxes do not exist any more, and we need to read in climatologically extended fluxes
        self.last_valid_ct13_year = 2012

class ExtendCTEmissions(Emissions):
    """
    Getting prior fluxes for recent years is a pain. CT2013 fluxes go till the end of 2012, while SiB-CASA goes till July
    2012. So we need to climatologically extend them to 2013 and 2014. For fossil fuel and oceanic fluxes, we will linearly
    extrapolate data from ten years (2003 to 2012). For fire and biosphere, we will take the average of the previous ten years.
    """
    def __init__(self, StartTime, EndTime):
        Emissions.__init__(self, StartTime, EndTime)
        # Which folder are the monthly flux files in?
        self.input_folder = self.rcf.get('flux.extrapolate.input.dir')
        self.output_folder = os.path.join(self.input_folder, 'extended')
        self.flux_file = self.rcf.get('flux.extrapolate.input.file')
        # Which variables should be linearly extrapolated?
        self.extrapolate_vars = self.rcf.get('flux.extrapolate.vars').split()
        self.average_vars = self.rcf.get('flux.average.vars').split()
        start_year = self.rcf.get('flux.existing.year.start', 'int')
        end_year = self.rcf.get('flux.existing.year.end', 'int')
        self.prev_years = arange(start_year, end_year+1)
        self.future_years = [int(y) for y in self.rcf.get('flux.extend.years').split()]
        # Are the fluxes three hourly or every 90 minutes?
        self.steps_per_day = self.rcf.get('flux.steps.per.day', 'int')
        self.flux_time_step = timedelta(hours=24)/self.steps_per_day

    def loadMonth(self, month):
        """
        Load the fluxes for a given month for all previous years. February will pose a special challenge. Luckily,
        we do not want to extend to cover a leap year (we will stop at 2015), so we can get away with throwing
        away the flux for Feb 29, should it exist on any given year.
        """
        # First, get the shape
        junk_dt = datetime(self.prev_years[0], month, 1)
        sample_filename = os.path.join(self.input_folder, junk_dt.strftime(self.flux_file))
        # The first year with existing flux could be a leap year, so get the number of days from a fixed non-leap year
        n_days = monthrange(1999,month)[1]

        with my_Dataset(sample_filename, 'r') as fid:
            try:
                n_lat = len(fid.dimensions['lat'])
                n_lon = len(fid.dimensions['lon'])
            except:
                n_lat = len(fid.dimensions['latitude'])
                n_lon = len(fid.dimensions['longitude'])
            n_tim = n_days * self.steps_per_day

        ret_dict = dict.fromkeys(self.extrapolate_vars + self.average_vars)

        num_prev_years = len(self.prev_years)
        for var_name in ret_dict.keys():
            ret_dict[var_name] = zeros((num_prev_years, n_tim, n_lat, n_lon), float64)

        # Now read all the previous years, one by one
        for i, year in enumerate(self.prev_years):
            junk_dt = datetime(year, month, 1)
            input_filename = os.path.join(self.input_folder, junk_dt.strftime(self.flux_file))
            with my_Dataset(input_filename, 'r') as fid:
                for var_name in ret_dict.keys():
                    ret_dict[var_name][i] = fid.variables[var_name][:n_tim,:,:]
        # Return those values
        return ret_dict

    def extrapolateMonth(self, month):
        """
        Create the fluxes for a given month for the next self.num_future_years years
        """
        ret_dict = dict.fromkeys(self.extrapolate_vars + self.average_vars)
        # First read the fluxes for the previous years
        data_dict = self.loadMonth(month)

        # Now create the extrapolated fluxes
        future_years = array([y-self.prev_years[0]+1 for y in self.future_years])
        for var_name in self.extrapolate_vars:
            ret_dict[var_name] = redistrib_flux.extrapolate3D(data_dict[var_name], future_years)

        for var_name in self.average_vars:
            _, n_time, n_lat, n_lon = data_dict[var_name].shape
            num_future_years = len(self.future_years)
            ret_dict[var_name] = zeros((num_future_years, n_time, n_lat, n_lon), float64)
            for i_year in range(num_future_years):
                ret_dict[var_name][i_year] = average(data_dict[var_name], axis=0)
        # Save some memory
        del data_dict
        return ret_dict

    def writeNewMonth(self, month):
        """
        Write out the new (extrapolated) fluxes to a file
        """
        data_dict = self.extrapolateMonth(month)
        for i_year, year in enumerate(self.future_years):
            junk_dt = datetime(year, month, 1)
            output_file = os.path.join(self.output_folder, junk_dt.strftime(self.flux_file))
            checkDir(output_file)

            # Do not overwrite existing files. This is needed because for some flux products, part years exist. E.g.,
            # SiB-CASA exists for Jan and Feb 2015, but not for the rest of the year.
            if os.path.exists(output_file):
                continue

            with my_Dataset(output_file, 'w') as fid:
                self.make_dimensions(fid)
                self.write_latlon(fid)
                self.write_datetimes(fid, year, month)

                for var_name in self.extrapolate_vars:
                    v = fid.createVariable(var_name, float64, ('date', 'lat', 'lon'))
                    v[:] = data_dict[var_name][i_year]
                    v.units = "mol m-2 s-1"
                    v.comment = 'Linear extrapolation of fluxes from years %i-%i'%(self.prev_years[0], self.prev_years[-1])
                for var_name in self.average_vars:
                    v = fid.createVariable(var_name, float64, ('date', 'lat', 'lon'))
                    v[:] = data_dict[var_name][i_year]
                    v.units = "mol m-2 s-1"
                    v.comment = 'Climatological average of fluxes from years %i-%i'%(self.prev_years[0], self.prev_years[-1])
                self.write_attributes(fid)

            print junk_dt.strftime("Flux for %b %Y written")

    def __call__(self):
        for month in range(1, 13):
            self.writeNewMonth(month)

    def make_dimensions(self, fid):
        fid.createDimension('lon', 360)
        fid.createDimension('lat', 180)
        fid.createDimension('date', 0)
        fid.createDimension('date_components', 6)

    def write_latlon(self, fid):
        v = fid.createVariable('lon', float64, ('lon',))
        v[:] = linspace(-179.5, 179.5, 360)
        v.units = "degrees_east"
        v.actual_range = array([-180., 180.])

        v = fid.createVariable('lat', float64, ('lat',))
        v[:] = linspace(-89.5, 89.5, 180)
        v.units = "degrees_north"
        v.actual_range = array([-90., 90.])

    def write_datetimes(self, fid, year, month):
        # First create the datetime objects
        num_days = monthrange(year, month)[1]
        num_time = num_days * self.steps_per_day
        all_datetimes = [datetime(year,month,1) + self.flux_time_step/2 + i*self.flux_time_step for i in range(num_time)]

        days_since = [d-datetime(2000,1,1) for d in all_datetimes]
        days_since = array([d.days + d.seconds/86400.0 for d in days_since], float64)
        v = fid.createVariable('date', float64, ('date',))
        v[:] = days_since
        v.units = "days since 2000-01-01 00:00:00 UTC"

        decimal_date = self.decimal_date(all_datetimes)
        v = fid.createVariable('decimal_date', float64, ('date',))
        v[:] = decimal_date
        v.units = "years"

        date_comps = array([d.timetuple()[:6] for d in all_datetimes], int32)
        v = fid.createVariable('date_components', int32, ('date', 'date_components'))
        v[:] = date_comps
        v.units = "NA"
        v.comment = "Calendar date components as integers.  Times and dates are UTC."
        v.order = "year, month, day, hour, minute, second"

        # SiB-CASA GFED4 files have a variable called 'time' which goes from 1 to the total number of timesteps in the file
        v = fid.createVariable('time', int16, ('date',))
        v[:] = arange(num_time, dtype=int16)
        v.units = "Time step of %i seconds"%(timedelta(hours=24)/self.steps_per_day).total_seconds()

    def decimal_date(self, datetime_array):
        # Given a datetime array, convert it into an array of floating point years
        ret_arr = zeros(len(datetime_array), float64)
        for i, d in enumerate(datetime_array):
            num_seconds_in_year = 86400.0 * (365.0 + isleap(d.year))
            dt = d - datetime(d.year, 1, 1, 0, 0, 0)
            num_seconds = dt.seconds + 86400.0 * dt.days
            ret_arr[i] = float64(d.year) + num_seconds/num_seconds_in_year
        return ret_arr

    def write_attributes(self, fid):
        fid.Notes = 'This file contains fluxes extrapolated from years %i-%i'%(self.prev_years[0], self.prev_years[-1])
        fid.email = 'Sourish.Basu@noaa.gov'
        fid.institution = "NOAA Earth System Research Laboratory"

class ReadFLT:
    """
    This is not a class to read and/or assemble emissions. This is a helper class to read/convert CABLE NEE
    estimates, which are in ArcGIS flt/hdr format. The code has been adapted from
    http://pydoc.net/Python/PyTOPKAPI/0.2.0/pytopkapi.arcfltgrid/
    """
    def __init__(self):
        self.dlat = 0.05
        self.dlon = 0.05
        self.rcf = rc.RcFile(os.environ['pyshell.rc'])
        self.CABLE_dir = self.rcf.get('emission.CO2.CABLE.dir')
        # filenames are YYYY/mth_FCNEP_YYYYMMMXX.(flt,hdr), where XX = # of days in month
        self.dS = None
        self.start_year = 1911
        self.end_year = 2011

    def createAreaGrid(self, (lat_specs), (lon_specs)):
        """
        Returns the area of a rectangular lat/lon grid in square meters. 'lat_specs' is a tuple
        (lat_beg, lat_end, lat_divs) and 'lon_end' is a similar tuple for longitudes.
        """
        lat_beg, lat_end, lats = lat_specs
        lon_beg, lon_end, lons = lon_specs
        EarthRad = 6.371e6 # meters
        dLon = (pi/180.) * (lon_end-lon_beg)/lons
        dS = zeros((lats+1, lons), float64)
        Lats = (pi/180.) * linspace(lat_beg, lat_end, lats+1)
        for i, lat in enumerate(Lats):
            dS[i] = EarthRad * EarthRad * dLon * sin(lat)
        dS = diff(dS, axis=0)
        return dS

    def read_bin(self, fname):
        """
        Read data from a ArcGIS binary file into an array. Read the data from a binary file created by ArcGIS into
        a Numpy array. The file is expected to be in binary format with floating point precision. e.g. ".flt" extension.
        """

        with open(fname, "rb") as f:
            data = fromstring(f.read(), 'f')

        return data

    def read(self, bingrid_name):
        """
        Read the data field and headers from an ArcGIS binary grid. This function reads the header and data from the
        ArcGIS binary data files produced by the "Raster to Float" tool in ArcGIS 9.1
        """

        if bingrid_name[-4:] == '.flt':
            hdr_name = bingrid_name[:-4]
            bin_name = bingrid_name
        else:
            hdr_name = bingrid_name
            bin_name = bingrid_name + '.flt'

        li_headers = self.read_headers(hdr_name)

        rows = li_headers['nrows']
        cols = li_headers['ncols']

        a = self.read_bin(bin_name)

        if (li_headers['byteorder'].lower() == 'lsbfirst' and sys.byteorder == 'big') or \
            (li_headers['byteorder'].lower() == 'msbfirst' and sys.byteorder == 'little'):
                a = a.byteswap()

        a = a.reshape(rows, cols)

        a = ma.masked_values(a, li_headers['nodata_value'])

        return a, li_headers

    def read_headers(self, bingrid_name):
        """
        Read the ascii headers of the ArcGIS binary grid file. The headers have the following format:

        ncols         62
        nrows         121
        xllcorner     -288595.47161281
        yllcorner     -3158065.5722693
        cellsize      1000
        nodata_value  -9999
        byteorder     LSBFIRST
        """

        hdr_name = bingrid_name + '.hdr'
        with open(hdr_name,'r') as f:
            tab_read=f.readlines()

        li_headers={}
        i=-1
        for line in tab_read:
            i=i+1
            donnees=line.split()
            if i in [0,1]:
                li_headers[donnees[0].strip()] = int(donnees[1])
            elif i in [2,3,4,5]:
                li_headers[donnees[0].strip()] = float(donnees[1])
            else:
                li_headers[donnees[0].strip()] = donnees[1].strip()

        return li_headers

    def total(self, year, month):
        num_days = calendar.monthrange(year, month)[1]
        file_base = os.path.join(self.CABLE_dir, "%4.4i"%year, "mth_FCNEP_%4.4i%2.2i%2.2i"%(year, month, num_days))
        data, header = self.read(file_base)
        if self.dS == None:
            lat_beg = header['xllcorner']
            lon_beg = header['yllcorner']
            lats = header['nrows']
            lons = header['ncols']
            lat_end = lat_beg + self.dlat * lats
            lon_end = lon_beg + self.dlon * lons
            dS = self.createAreaGrid((lat_beg,lat_end,lats), (lon_beg,lon_end,lons))
            self.dS = ma.masked_array(dS, mask=data.mask)
        # unit of 'data' is gC/m^2/day, we need to convert to TgC/month
        per_grid = num_days * data * self.dS # now in gC/cell/month
        return per_grid.sum() * 1.0E-12

    def extract_monthly_totals(self, outfile):
        date_array = []
        emis_array = []
        for year in range(self.start_year, self.end_year+1):
            print 'Calculating year %i totals... '%year,
            for month in range(1,13):
                em = self.total(year,month)
                date_array.append([year,month,1])
                emis_array.append(em)
            print 'done'
        if os.path.exists(outfile):
            os.remove(outfile)
        n_obs = len(emis_array)
        date_array = array(date_array, int16)
        emis_array = array(emis_array, float64)
        with my_Dataset(outfile, 'w') as fid:
            fid.createDimension('n_month', n_obs)
            fid.createDimension('n_date', 3)
            var = fid.createVariable('date', int16, ('n_month', 'n_date'))
            var[:] = date_array
            var = fid.createVariable('emis', float64, ('n_month',))
            var[:] = emis_array
            var.unit = 'TgC/month'
