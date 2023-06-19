#!/usr/bin/env python

import logging
from numpy import linspace, pi, float64, zeros, sin, diff, arange, random, int32, array, int16
import os
from collections import defaultdict
import sys
from dateutil.relativedelta import relativedelta
from datetime import timedelta, datetime
from pyshell.utilities import checkDir, my_Dataset
from pyshell.gridtools import TM5Grids
from pandas import DatetimeIndex, Timedelta
from pandas.tseries.frequencies import to_offset
from pandas import Timestamp
import xarray as xr


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


class Emissions:
    def __init__(self, rcf):
        self.StartDate = rcf.ti
        self.EndDate = rcf.tf
        self.subdir_tag = ''
        self.sec_day = 86400.
        self.rcf = rcf
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.GetZoomRegions()
        # Monte Carlo estimation of posterior covariance or not?
        self.optim_mc = self.rcf.get('my.optimizer.class') in ['conGrad_MC', 'm1qn3_MC']
        emission_file_name = self.rcf.get('PyShell.em.filename')
        self.emission_file_name = self.putDateString(emission_file_name, True)
        self.createDifferentialArea()
        # Define a destination horizontal grid for later
        self.lat_grid = {}
        self.lon_grid = {}
        for region_name, lat_data in zip(self.zoom_regions_names, self.zoom_regions_lat):
            ybeg, yend, jm = lat_data
            self.lat_grid[region_name] = linspace(ybeg, yend, jm+1)
        for region_name, lon_data in zip(self.zoom_regions_names, self.zoom_regions_lon):
            xbeg, xend, im = lon_data
            self.lon_grid[region_name] = linspace(xbeg, xend, im+1)

    def putDateString(self, input_string, file=True):
        output_string = input_string
        for atom in ['Y', 'y', 'm', 'M', 'd', 'H', 'b', 'B', 'j', 'p', 'S']:
            for i, t in zip([1,2], [self.StartDate, self.EndDate]):
                rep_string = '<%1s%1i>' % (atom, i)
                form_string = '%%%1s' % atom
                output_string = output_string.replace(rep_string, t.strftime(form_string))
        # Now add the subdir_tag
        if file:
            dirname, filename = os.path.split(output_string)
            output_string = os.path.join(dirname, self.subdir_tag, filename)
        else:
            output_string = os.path.join(output_string, self.subdir_tag)
        return output_string

    def GetZoomRegions(self):
        """
        Reads a config file to get zoom region info
        """
        self.zoom_regions_names = self.rcf.get('regions').split()
        xbeg = []
        xend = []
        im = []
        ybeg = []
        yend = []
        jm = []
        self.zoom_info = defaultdict(dict)
        for region in self.zoom_regions_names:
            xbeg.append(self.rcf.get('region.'+region+'.xbeg', 'float'))
            self.zoom_info[region]['xbeg'] = self.rcf.get('region.'+region+'.xbeg', 'float')
            ybeg.append(self.rcf.get('region.'+region+'.ybeg', 'float'))
            self.zoom_info[region]['ybeg'] = self.rcf.get('region.'+region+'.ybeg', 'float')
            xend.append(self.rcf.get('region.'+region+'.xend', 'float'))
            self.zoom_info[region]['xend'] = self.rcf.get('region.'+region+'.xend', 'float')
            yend.append(self.rcf.get('region.'+region+'.yend', 'float'))
            self.zoom_info[region]['yend'] = self.rcf.get('region.'+region+'.yend', 'float')
            im.append(self.rcf.get('region.'+region+'.im', 'int'))
            self.zoom_info[region]['im'] = self.rcf.get('region.'+region+'.im', 'int')
            jm.append(self.rcf.get('region.'+region+'.jm', 'int'))
            self.zoom_info[region]['jm'] = self.rcf.get('region.'+region+'.jm', 'int')
        self.zoom_regions_lon = zip(xbeg, xend, im)
        self.zoom_regions_lat = zip(ybeg, yend, jm)

    def createDifferentialArea(self):
        """
        Create a bunch of arrays containing the areas of regular gridboxes in square meters
        """
        self.dS_m2 = {}
        for lats, lons in zip([45,90,180,360,1800], [60,120,360,720,3600]):
            self.dS_m2[(lats,lons)] = self.diffArea(-90., 90., lats, -180., 180., lons)

    def diffArea(self, lat_min, lat_max, lats, lon_min, lon_max, lons):
        EarthRad = 6.371e6 # meters
        dLon = (pi/180.) * ((lon_max-lon_min)/lons)
        dS = zeros((lats+1, lons), float64)
        Lats = (pi/180.) * linspace(lat_min, lat_max, lats+1)
        for i, lat in enumerate(Lats):
            dS[i] = EarthRad * EarthRad * dLon * sin(lat)
        dS = diff(dS, axis=0)
        return dS

    def create_emission_structure(self):
        # create dictionary to contain emission info to be read from rc file:
        self.Emission = dict.fromkeys(self.zoom_regions_names)
        # self.Emission = TM5_emission.fromkeys(self.zoom_regions_names)
        self.Emission['regions'] = self.zoom_regions_names
        for region in self.zoom_regions_names:
            self.Emission[region] = dict.fromkeys(self.species)
            self.Emission[region]['tracers'] = self.species
        for tracer in self.species:
            cat_opt = defaultdict(list)
            cat_list = []
            optim_cat_names = []
            self.Emission[tracer] = {}
            # self.Emission[tracer]['emi_class'] = self.rcf.get(tracer+'.emission.class')
            for region, xlims, ylims in zip(self.zoom_regions_names, self.zoom_regions_lon, self.zoom_regions_lat):
                self.Emission[region][tracer] = {}
                ncat = self.rcf.get('emission.'+tracer+'.'+region+'.categories','int')
                categories = []
                for icat in arange(ncat):
                    infoline = self.rcf.get('emission.'+tracer+'.'+region+'.category%1i'%(icat+1))
                    categ = infoline[:infoline.find(';')].strip()
                    categories.append(categ)
                    error = float64(infoline.split(';')[1])
                    corr = infoline.split(';')[2]
                    tcorr = infoline.split(';')[3]
                    time_res = infoline.split(';')[3]
                    res_key = time_res.split('-')[2]
                    optimize = int(infoline.split(';')[4])
                    remarks = infoline.split(';')[5].strip()
                    time_interval = self.split_time_interval(self.StartDate, self.EndDate, res_key)
                    nt = len(time_interval['time_start'])
                    self.Emission[region][tracer][categ] = {
                        'emission_data'  : zeros((nt, int(ylims[2]), int(xlims[2])), float64),
                        'time_resolution': res_key,
                        'time_interval'  : time_interval,
                        'error'          : error,
                        'corr'           : corr,
                        'tcorr'          : tcorr,
                        'optimize'       : optimize,
                        'remarks'        : remarks
                    }
                    # sometimes we want to specify the errors in emissions separately, instead of as a single fraction of the emission
                    sep_err = self.rcf.get('%s.%s.sep_error'%(tracer, categ), 'bool', default=False)
                    if sep_err:
                        self.Emission[region][tracer][categ]['emission_error'] = zeros((nt, int(ylims[2]), int(xlims[2])), float64)
                    # create library with tuple as key:
                    if optimize == 1:
                        tup = (categ, corr, tcorr)
                        cat_opt[tup].append({'region': region, 'error': error})
                        if tup not in cat_list:
                            cat_list.append(tup)
                        if categ not in optim_cat_names:
                            optim_cat_names.append(categ)
                self.Emission[region][tracer]['categories'] = categories
            self.Emission[tracer]['cat_list'] = cat_list
            self.Emission[tracer]['cat_opt'] = cat_opt
            self.Emission[tracer]['optim_cat_names'] = optim_cat_names

    def readOptimFromFile(self, file_name, add_perturb = False):
        """
        Sometimes we want to read the emissions from a single file instead of assembling it. Two most common
        examples of this need are (a) doing a forward run with optimized emissions, (b) doing an inversion with
        an unbiased prior, where prior = posterior + Gaussian noise. It assumes that self.Emission has already
        been initialized, and only reads in the emission data, and the emission errors where relevant.
        """
        fid = my_Dataset(file_name, 'r')
        for region in self.zoom_regions_names:
            rgid = fid.groups[region]
            for tracer in self.species:
                tgid = rgid.groups[tracer]
                categories = self.Emission[region][tracer]['categories']
                for cat in categories:
                    cgid = tgid.groups[cat]
                    self.Emission[region][tracer][cat]['emission_data'] = cgid.variables['emission'][:]
                    if 'emission_error' in cgid.variables:
                        self.Emission[region][tracer][cat]['emission_error'] = cgid.variables['emission_error'][:]
                    else:
                        error = self.Emission[region][tracer][cat]['error']
                        self.Emission[region][tracer][cat]['emission_error'] = 0.01 * error * \
                            abs(self.Emission[region][tracer][cat]['emission_data'])
                    if add_perturb:
                        perturb = random.standard_normal(self.Emission[region][tracer][cat]['emission_data'].shape)
                        self.Emission[region][tracer][cat]['emission_data'] += \
                            (perturb * self.Emission[region][tracer][cat]['emission_error'])
        fid.close()

    def split_time_interval(self, StartDate, EndDate, res_key):
        if res_key.find('+') == -1:
            nres = 1
        else:
            nres = int(res_key.split('+')[1])
        if res_key.find('monthly') != -1:
            delta = relativedelta(months=+nres)
            if nres == 1:
                delta2 = timedelta(days=15)
            else:
                logger.error('optimization for periods longer than month not yet implemented')
                sys.exit()
            # set to start of month: user should start optimization at 1st day month, and end last day month
            # SB: Try to start from the correct start date and not the first of the month
            # TheStartDate = datetime(StartDate.year,StartDate.month,1,0,0,0)
        elif res_key.find('daily') != -1:
            # TheStartDate = datetime(StartDate.year,StartDate.month,StartDate.day,0,0,0)
            delta = timedelta(days=+nres)
            delta2 = delta/2
        else:
            logger.error("invalid time resolution:", res_key)
            sys.exit()
        # d = TheStartDate
        d = datetime(StartDate.year,StartDate.month,StartDate.day,0,0,0)
        time_start = []
        time_end = []
        time_mid = []
        while d < EndDate:
            time_start.append(d)
            # what is the ending time for this window?
            if res_key.find('monthly') != -1:
                d_end = self.nextMonth(d)
            elif res_key.find('daily') != -1:
                d_end = d + delta

            # Need to check here whether d+delta exceeds the end date or not
            if d_end <= EndDate:
                time_end.append(d_end)
                time_mid.append(d + (d_end-d)/2)
            else:
                time_end.append(EndDate)
                time_mid.append(d + (EndDate - d)/2)

            d = d_end
        time_interval = {'time_start': time_start, 'time_end': time_end, 'time_mid': time_mid, 'dt': delta}
        return time_interval

    def nextMonth(self, cur_date):
        # return the beginning of next month
        year, month = cur_date.timetuple()[:2]
        if month == 12:
            next_year = year + 1
            next_month = 1
        else:
            next_year = year
            next_month = month + 1
        next_date = datetime(next_year, next_month, 1, 0, 0, 0)
        return next_date

    def __call__(self, randomize=False, zero=False):
        if randomize or zero:
            raise NotImplementedError
            # raise RuntimeError('You can choose to either randomize emissions or zero them out, but not both')
        if randomize:
            raise NotImplementedError
            # self.Emission.setRandom()
        elif zero:
            raise NotImplementedError
            # self.Emission.setZero()
        self.WriteEmissions()
        return self.Emission

    def WriteEmissions(self):
        outFile = self.emission_file_name
        checkDir(outFile)
        if os.path.exists(outFile):
            os.remove(outFile)
        fid = my_Dataset(outFile, 'w')
        fid.createDimension('itime', 6)
        for region, xlims, ylims in zip(self.zoom_regions_names, self.zoom_regions_lon, self.zoom_regions_lat):
            group = fid.createGroup(region)
            group.createDimension('latitude', ylims[2])
            group.createDimension('longitude', xlims[2])
            group.latmin = ylims[0]
            group.latmax = ylims[1]
            group.lonmin = xlims[0]
            group.lonmax = xlims[1]
            for tracer in self.species:
                tgroup = group.createGroup(tracer)
                categories = self.Emission[region][tracer]['categories']
                for cat in categories:
                    cgroup = tgroup.createGroup(cat)
                    cgroup.createDimension('nt', self.Emission[region][tracer][cat]['emission_data'].shape[0])
                    cgroup.time_resolution = self.Emission[region][tracer][cat]['time_resolution']
                    cgroup.optimize = int32(self.Emission[region][tracer][cat]['optimize'])
                    # This is Maarten's modification for CO biomass burning, applicable to other tracers as well
                    if self.Emission[tracer].get('tf_bb_diurnal'):
                        cgroup.createDimension('hourly', 24)
                        var = cgroup.createVariable('tf_diurnal', float64, ('hourly',))
                        var[:] = self.Emission[tracer]['tf_bb_diurnal']

                    var = cgroup.createVariable('time_start', int16, ('nt', 'itime'))
                    var[:] = array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_start']], int16)

                    var = cgroup.createVariable('time_end', int16, ('nt', 'itime'))
                    var[:] = array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_end']], int16)

                    var = cgroup.createVariable('time_mid', int16, ('nt', 'itime'))
                    var[:] = array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_mid']], int16)

                    var = cgroup.createVariable('emission', self.Emission[region][tracer][cat]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                    var[:] = self.Emission[region][tracer][cat]['emission_data']

                    write_err = self.rcf.get('%s.%s.sep_error' % (tracer, cat), 'bool', default=False)
                    if write_err:
                        var = cgroup.createVariable('emission_error', self.Emission[region][tracer][cat]['emission_error'].dtype, ('nt', 'latitude', 'longitude'))
                        var[:] = self.Emission[region][tracer][cat]['emission_error']

                    # Write the total emissions per category per zoom region per timestep, for possible debugging
                    # The emissions are in Kg tracer/grid cell/sec
                    nt = self.Emission[region][tracer][cat]['emission_data'].shape[0]
                    emis_totals = zeros(nt, float64)
                    time_intervals = zeros(nt, float64)
                    for i in range(nt):
                        t1 = self.Emission[region][tracer][cat]['time_interval']['time_start'][i]
                        t2 = self.Emission[region][tracer][cat]['time_interval']['time_end'][i]
                        dt = (t2-t1).total_seconds()
                        emis_totals[i] = self.Emission[region][tracer][cat]['emission_data'][i].sum() * dt
                        time_intervals[i] = dt

                    var = cgroup.createVariable('emission_total', emis_totals.dtype, ('nt',))
                    var[:] = emis_totals
                    var.comment = 'Total emission for region %s and category %s' % (region, cat)
                    var.unit = 'Kg %s/time step' % tracer

                    var = cgroup.createVariable('time_step_length', time_intervals.dtype, ('nt',))
                    var[:] = time_intervals
                    var.comment = 'Length of time period (seconds)'

                    del emis_totals, time_intervals

        fid.close()


class PreprocessedEmissions(Emissions):
    def __init__(self, rcf, dconf, tracer='CO2', *args, **kwargs):
        Emissions.__init__(self, rcf) #, *args, **kwargs)
        self.dconf = dconf
        self.tracer = tracer
        self.writeCycle = self.rcf.get('CO2.emission.dailycycle', 'bool')
        self.MolarMass = {'CO2': 44.00995e-3}[tracer]
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
        try :
            data = xr.open_mfdataset(file_pattern)
        except IOError as e :
            logger.error('No file matching pattern ' + file_pattern + ' found')
            logger.exception(e)
            raise e

        if area is None :
            area = TM5Grids.global1x1().area
        else :
            area = data[area].values

        # Calculate granularity:
        # Granularity (time step) is read, by order of priority 1) in the netcdf file attributes, 2) in the yaml file
        granularity = data.attrs.get('granularity', self.dconf[category].get('granularity'))

        # Load the emissions. They should be on a global 1x1, 3-hourly grid and in mol/m2/s
        # The following will convert them in kg[tracer]/s
        idx = DatetimeIndex(data.time.values)
        emis_glo1x1 = data[field].sel(time=(idx >= period.start) & (idx < period.stop)) * area * self.MolarMass

        # Check that the simulation period is fully covered by the preprocessed emissions:
        tmin = Timestamp(data.time.values.min())
        tmax = Timestamp(data.time.values.max())
        assert (tmin <= period.start) & (tmax + to_offset(granularity) >= period.stop), logger.error(
             "Pre-processed emission files for category %s (%s) don't cover the full simulation period. " % (category, file_pattern) +
             "Emissions in the files range from %s to %s" % (tmin.strftime("%-d %B %Y"), tmax.strftime("%-d %B %Y")))

        # destination region
        destreg = TM5Grids.from_corners(latb = self.lat_grid[region], lonb = self.lon_grid[region])

        # regrid:
        try :
            emis_coarsened = crop_and_coarsen(emis_glo1x1, latb=destreg.latb, lonb=destreg.lonb)
        except ValueError as e:
            logger.error("Couldn't import emissions. Maybe some files are missing?")
            raise e

        timestep = DatetimeIndex(emis_glo1x1.time.values) + to_offset(granularity) - DatetimeIndex(emis_glo1x1.time.values)

        # Return a DataArray:
        return xr.Dataset(
            data_vars = {
                'emis': (('time', 'latitude', 'longitude'), emis_coarsened),
                'timestep': (('time', ), [t.total_seconds() for t in timestep])
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

                    # for emis in emission_fine:
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
                            ds.to_netcdf(filename, group = region + '/' + cat, mode = write_mode)
                            cur_day += timedelta(days=1)
