#!/usr/bin/env python2.7

from dateutil.relativedelta import relativedelta
from datetime import timedelta
from netCDF4 import Dataset
from numpy import zeros, int32, float64, sqrt, arange, ones, datetime64, array, int16
from pandas import DatetimeIndex
from collections import OrderedDict
import os
# from pyshell.base.helper.Utilities import checkDir
from pyshell.utilities import checkDir
import xarray as xr
import logging


logger = logging.getLogger(__name__)


class Observations:
    def __init__(self, ti, tf, rcf):
        self.ti = ti
        self.tf = tf
        self.rcf = rcf
        self.species = [s.strip() for s in rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.point_split_period = rcf.get('output.point.split.period', default='a')
        self.PointObservation = {}
        for tracer in self.species:
            self.PointObservation[tracer] = {}

    def SetupPointObs(self, tracer, ds):
        assimilated = self.rcf.get('%s.point.assimilate' % tracer, 'bool', default=True)
        if not assimilated:
            unassim_mdm = self.rcf.get('%s.point.unassimilate.mdm' % tracer, 'float')
        # obs = obsobj.get(tracer, self.ti, self.tf)

        ds = ds.sel(index = (ds.time > datetime64(self.ti)) & (ds.time <= datetime64(self.tf)) & (ds.tracer == tracer.lower()))

        #self.PointObservation[tracer]['dimensions'] = {'id': obs['n']}
        self.PointObservation[tracer]['dimensions'] = {'id': len(ds.obs)}
        self.PointObservation[tracer]['variable_shapes'] = OrderedDict()
        self.PointObservation[tracer]['attr_dict'] = OrderedDict()
        self.PointObservation[tracer]['variable_attrs'] = OrderedDict()
        self.PointObservation[tracer]['id'] = arange(1, len(ds.obs) + 1, dtype=int32)
        self.PointObservation[tracer]['variable_shapes']['id'] = ('id',)
        self.PointObservation[tracer]['lat'] = ds.latitude.values.astype(float64)
        self.PointObservation[tracer]['variable_shapes']['lat'] = ('id',)
        self.PointObservation[tracer]['lon'] = ds.longitude.values.astype(float64)
        self.PointObservation[tracer]['variable_shapes']['lon'] = ('id',)
        self.PointObservation[tracer]['alt'] = ds.sampling_altitude.values.astype(float64)
        self.PointObservation[tracer]['variable_shapes']['alt'] = ('id',)
        self.PointObservation[tracer]['mixing_ratio'] = ds.obs.values.astype(float64)
        self.PointObservation[tracer]['variable_shapes']['mixing_ratio'] = ('id',)
        if assimilated:
            self.PointObservation[tracer]['mixing_ratio_error'] = ds.err_obs.values.astype(float64)
        else:
            self.PointObservation[tracer]['mixing_ratio_error'] = unassim_mdm * ones(len(ds.obs), float64)
            self.PointObservation[tracer]['variable_attrs']['mixing_ratio_error'] = [('comment', 'These data will not be assimilated')]
        self.PointObservation[tracer]['variable_shapes']['mixing_ratio_error'] = ('id',)
        self.PointObservation[tracer]['date_components'] = array([t.timetuple()[:6] for t in DatetimeIndex(ds.time.values)]).astype(int16)
        self.PointObservation[tracer]['variable_shapes']['date_components'] = ('id', 'idate')
        self.PointObservation[tracer]['sampling_strategy'] = ds.sampling_strategy.values.astype(int16)
        self.PointObservation[tracer]['variable_shapes']['sampling_strategy'] = ('id',)
        self.PointObservation[tracer]['station_id'] = ds.TM5_station_ID.values.astype(int32)
        self.PointObservation[tracer]['variable_shapes']['station_id'] = ('id',)
        self.PointObservation[tracer]['time_window_length'] = (ds.time_window_length.values * 1.e-9).astype(int32)
        self.PointObservation[tracer]['variable_shapes']['time_window_length'] = ('id',)

        # Also write the obs file to the output folder, for reference:
        dir_name = self.rcf.get('output.point.input.dir')
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        if len(ds.obs) > 0 :
            ds.to_netcdf(os.path.join(dir_name, 'observations.%s.nc' % tracer))

    def writePointFile(self):
        dir_name = self.rcf.get('output.point.input.dir')
        if not os.path.exists(dir_name): os.makedirs(dir_name)
        if self.point_split_period == 'a':
            file_name = 'point_input.nc4'
        elif self.point_split_period == 'm':
            file_name = self.StartDate.strftime("point_input_%Y%m.nc4")
        elif self.point_split_period == 'd':
            file_name = self.StartDate.strftime("point_input_%Y%m%d.nc4")
        point_file = os.path.join(dir_name, file_name)
        if os.path.exists(point_file): os.remove(point_file)

        # only create a file if there are observations to write
        total_obs = 0
        for tracer in self.species:
            total_obs += self.PointObservation[tracer]['dimensions']['id']
        if total_obs == 0: return

        file_id = Dataset(point_file, 'w')
        file_id.createDimension('idate', 6)
        for tracer in self.species:
            assim_tracer = self.rcf.get('%s.point.assimilate' % tracer, 'bool', default=True)
            if not assim_tracer:
                n_obs = self.PointObservation[tracer]['dimensions']['id']
                self.PointObservation[tracer]['mixing_ratio_error'] = 1.0E10 * ones(n_obs, dtype=float64)

            fid = file_id.createGroup(tracer)
            for dim_name, dim_len in self.PointObservation[tracer]['dimensions'].items():
                fid.createDimension(dim_name, dim_len)
            for var_name, var_shape in self.PointObservation[tracer]['variable_shapes'].items():
                var_value = self.PointObservation[tracer][var_name]
                v = fid.createVariable(var_name, var_value.dtype, var_shape)
                v[:] = var_value
                if var_name in self.PointObservation[tracer]['variable_attrs']:
                    for attr_name, attr_value in self.PointObservation[tracer]['variable_attrs'][var_name]:
                        setattr(v, attr_name, attr_value)
            for attr_name, attr_value in self.PointObservation[tracer]['attr_dict'].items():
                setattr(fid, attr_name, attr_value)
        file_id.close()
        logger.info("Wrote %s" % point_file)


class InputObservations:
    def __init__(self):
        self.n = 0
        self.tracers = []

    @classmethod
    def from_xr_nc(cls, filename_or_pattern):
        obs = cls()
        data = xr.open_mfdataset(filename_or_pattern)
        obs.data2 = data
        return obs


class PointObs:
    def __init__(self, model):
        self.StartTime = model.StartTime
        self.EndTime = model.EndTime
        self.rcf = model.rcf
        self.Optim_dict = model.Optim_dict
        self.output_dir = model.output_dir

        self.runid = self.rcf.get('runid')
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        # self.track_file = self.getTrackFileName()
        self.region_names = self.rcf.get('regions').split()
        # Monte Carlo estimation of posterior covariance or not?
        self.optim_mc = self.rcf.get('my.optimizer.class') in ['conGrad_MC', 'm1qn3_MC']
        self.months = []
        self.split_period = self.rcf.get('output.point.split.period',
                                         default='a')  # 'm' if files are split monthly, 'd' if split daily, 'a' if no split
        if self.split_period == 'm':
            d = self.StartTime.replace(day=1).replace(hour=0).replace(minute=10)
            while d <= self.EndTime:
                self.months.append((d.year, d.month))
                d += relativedelta(months=1)
        elif self.split_period == 'd':
            d = self.StartTime.replace(hour=0).replace(minute=10)
            while d <= self.EndTime:
                self.months.append((d.year, d.month, d.day))
                d += timedelta(days=1)
        elif self.split_period == 'a':
            self.months.append(None)

    def PointCost(self):
        J_obs = 0.0
        for month_tuple in self.months:
            track_file = self.getMismatchFileName(month_tuple)
            if os.path.exists(track_file):
                track_fid = Dataset(track_file, 'r')
                for region, region_data in track_fid.groups.items():
                    for tracer, tracer_data in region_data.groups.items():
                        J_obs += 0.5 * sum((tracer_data.variables['mismatch'][:]/tracer_data.variables['error'][:])**2)
                track_fid.close()
        return J_obs

    def create_pointdeparture_structure(self):
        self.PointDepartures = dict.fromkeys(self.region_names)
        for region in self.region_names:
            self.PointDepartures[region] = {}
            for tracer in self.species:
                self.PointDepartures[region][tracer] = {}
        for tracer in self.species:
            self.PointDepartures[tracer] = {}
            self.PointDepartures[tracer]['departure_class'] = self.rcf.get(tracer+'.departures.point.class')

    def applyPointObs(self, tracer, month_tuple = None):
        mdm_obs_only = self.rcf.get('point.%s.only.obs.error' % tracer, 'bool', default=False)
        track_file = self.getTrackFileName(month_tuple)
        if not os.path.exists(track_file):
            return
        track_fid = Dataset(track_file, 'r')
        for region in self.region_names:
            group_data = track_fid.groups[region]
            if tracer in group_data.groups and len(group_data.groups[tracer].dimensions['samples']) > 0:
                date_components, lat, lon, alt, mixing_ratio, obs_error, tw_length = self.getInputVars(
                    group_data.groups[tracer], tracer, month_tuple)
            mod_mixing = group_data.groups[tracer].variables['mixing_ratio'][:]
            mod_error = group_data.groups[tracer].variables['mixing_ratio_sigma'][:]

            self.PointDepartures[region][tracer]['dimensions'] = {'samples': len(mod_mixing)}
            self.PointDepartures[region][tracer]['variable_shapes'] = OrderedDict()
            self.PointDepartures[region][tracer]['date_components'] = int32(date_components)
            self.PointDepartures[region][tracer]['variable_shapes']['date_components'] = ('samples', 'idate')
            self.PointDepartures[region][tracer]['lat'] = float64(lat)
            self.PointDepartures[region][tracer]['variable_shapes']['lat'] = ('samples',)
            self.PointDepartures[region][tracer]['lon'] = float64(lon)
            self.PointDepartures[region][tracer]['variable_shapes']['lon'] = ('samples',)
            self.PointDepartures[region][tracer]['alt'] = float64(alt)
            self.PointDepartures[region][tracer]['variable_shapes']['alt'] = ('samples',)
            self.PointDepartures[region][tracer]['nsamples'] = group_data.groups[tracer].variables['nsamples'][:]
            self.PointDepartures[region][tracer]['variable_shapes']['nsamples'] = ('samples',)
            self.PointDepartures[region][tracer]['total_weight'] = group_data.groups[tracer].variables['total_weight'][
                                                                   :]
            self.PointDepartures[region][tracer]['variable_shapes']['total_weight'] = ('samples',)
            self.PointDepartures[region][tracer]['sampling_strategy'] = group_data.groups[tracer].variables[
                                                                            'sampling_strategy'][:]
            self.PointDepartures[region][tracer]['variable_shapes']['sampling_strategy'] = ('samples',)
            self.PointDepartures[region][tracer]['time_window_length'] = int32(tw_length)
            self.PointDepartures[region][tracer]['variable_shapes']['time_window_length'] = ('samples',)
            self.PointDepartures[region][tracer]['model_error'] = float64(mod_error)
            self.PointDepartures[region][tracer]['variable_shapes']['model_error'] = ('samples',)
            self.PointDepartures[region][tracer]['variable_shapes']['model_error'] = ('samples',)
            self.PointDepartures[region][tracer]['obs_error'] = float64(obs_error)
            self.PointDepartures[region][tracer]['variable_shapes']['obs_error'] = ('samples',)
            # Should the total error be only the measurement error or the total error?
            if mdm_obs_only:
                total_error = obs_error
            else:
                total_error = sqrt(float64(obs_error) ** 2 + float64(mod_error) ** 2)
            self.PointDepartures[region][tracer]['error'] = float64(total_error)
            self.PointDepartures[region][tracer]['variable_shapes']['error'] = ('samples',)
            self.PointDepartures[region][tracer]['mismatch'] = float64(mod_mixing - mixing_ratio)
            self.PointDepartures[region][tracer]['variable_shapes']['mismatch'] = ('samples',)
            self.PointDepartures[region][tracer]['forcing'] = float64(mod_mixing - mixing_ratio) / total_error ** 2
            self.PointDepartures[region][tracer]['variable_shapes']['forcing'] = ('samples',)
        track_fid.close()

    def getMismatchFileName(self, month_tuple):
        if month_tuple is None:
            return os.path.join(self.output_dir, 'point', 'point_departures.nc4')
        elif len(month_tuple) == 2:
            return os.path.join(self.output_dir, 'point', 'point_departures_%04i%02i.nc4'%month_tuple)
        elif len(month_tuple) == 3:
            return os.path.join(self.output_dir, 'point', 'point_departures_%04i%02i%02i.nc4'%month_tuple)

    def writePointDepartureFile(self, month_tuple):
        self.mismatch_file = self.getMismatchFileName(month_tuple)
        if os.path.exists(self.mismatch_file):
            os.remove(self.mismatch_file)
        checkDir(self.mismatch_file)
        # do we need to write this file?
        total_nobs = 0
        for region in self.region_names:
            for tracer in self.species:
                total_nobs += len(self.PointDepartures[region][tracer].keys()) # not really the number of observations ...
        if total_nobs == 0:
            return
        fid = Dataset(self.mismatch_file, 'w')
        fid.createDimension('tracer', self.ntracer)
        fid.createDimension('idate', 6)
        for region in self.region_names:
            # do we need to write this group?
            total_nobs = 0
            for tracer in self.species:
                total_nobs += len(self.PointDepartures[region][tracer].keys())
            if total_nobs == 0:
                continue # with the next region
            rgrp = fid.createGroup(region)
            for tracer in self.species:
                if len(self.PointDepartures[region][tracer].keys()) == 0:
                    # there are no observations for this tracer in this region
                    continue # with the next tracer
                tgrp = fid.groups[region].createGroup(tracer)
                for dim_name, dim_len in self.PointDepartures[region][tracer]['dimensions'].items():
                    tgrp.createDimension(dim_name, dim_len)
                for var_name, var_shape in self.PointDepartures[region][tracer]['variable_shapes'].items():
                    var_value = self.PointDepartures[region][tracer][var_name]
                    v = tgrp.createVariable(var_name, var_value.dtype, var_shape)
                    v[:] = var_value
        fid.close()

    def getTrackFileName(self, month_tuple):
        if month_tuple is None:
            return os.path.join(self.output_dir, 'point', 'point_output.nc4')
        elif len(month_tuple) == 2:
            return os.path.join(self.output_dir, 'point', 'point_output_%04i%02i.nc4'%month_tuple)
        elif len(month_tuple) == 3:
            return os.path.join(self.output_dir, 'point', 'point_output_%04i%02i%02i.nc4'%month_tuple)

    def getInputVars(self, gid, tracer, month_tuple):
        sample_id = gid.variables['id'][:] - 1
        with Dataset(self.getInputFileName(month_tuple), 'r') as ifid:
            date_components = ifid.groups[tracer].variables['date_components'][:][sample_id]
            mixing_ratios = ifid.groups[tracer].variables['mixing_ratio'][:][sample_id]
            mixing_ratio_errs = ifid.groups[tracer].variables['mixing_ratio_error'][:][sample_id]
            lat = ifid.groups[tracer].variables['lat'][:][sample_id]
            lon = ifid.groups[tracer].variables['lon'][:][sample_id]
            alt = ifid.groups[tracer].variables['alt'][:][sample_id]
            # The variable 'time_window_length' may or may not exist
            if 'time_window_length' in ifid.groups[tracer].variables:
                time_window_length = ifid.groups[tracer].variables['time_window_length'][:][sample_id]
            else:
                time_window_length = zeros(len(alt), int32)
        return date_components, lat, lon, alt, mixing_ratios, mixing_ratio_errs, time_window_length

    def getInputFileName(self, month_tuple):
        indir_point = self.rcf.get('output.point.input.dir')
        if month_tuple is None:
            return os.path.join(indir_point, 'point_input.nc4')
        elif len(month_tuple) == 2:
            return os.path.join(indir_point, 'point_input_%04i%02i.nc4'%month_tuple)
        elif len(month_tuple) == 3:
            return os.path.join(indir_point, 'point_input_%04i%02i%02i.nc4'%month_tuple)