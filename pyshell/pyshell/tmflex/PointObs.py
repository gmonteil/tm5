#!/bin/env python
from pyshell.base.main.PointObs import PointObs as PointObs_base
from netCDF4 import Dataset
from numpy import *
from collections import OrderedDict
import os

class PointObs(PointObs_base):
    def __init__(self, tm5Obj):
        self.StartTime = tm5Obj.StartTime
        self.EndTime = tm5Obj.EndTime
        self.rcf = tm5Obj.rcf
        self.runid = self.rcf.get('runid')
        self.Optim_dict = tm5Obj.Optim_dict
        self.output_dir = tm5Obj.output_dir
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
        # PointObs_base.__init__(self, tm5Obj)

    def applyPointObs(self, tracer, month_tuple=None):
        print(tracer)
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
