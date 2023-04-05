#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True

from pyshell.tmflex import rc
import re, os
from numpy import *
from pyshell.base.helper.Utilities import *
from netCDF4 import Dataset
from dateutil.relativedelta import relativedelta

class PointObs(object):

    def __init__(self, tm5Obj):
        self.StartTime = tm5Obj.StartTime
        self.EndTime   = tm5Obj.EndTime
        self.rcf = tm5Obj.rcf # rc.RcFile(os.environ['pyshell.rc'])
        self.runid = self.rcf.get('runid')
        self.Optim_dict = tm5Obj.Optim_dict
        self.output_dir = tm5Obj.output_dir
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        #self.track_file = self.getTrackFileName()
        self.region_names = self.rcf.get('regions').split()
        # Monte Carlo estimation of posterior covariance or not?
        self.optim_mc = self.rcf.get('my.optimizer.class') in ['conGrad_MC', 'm1qn3_MC']
        self.months = []
        self.split_period = self.rcf.get('output.point.split.period', default='a') # 'm' if files are split monthly, 'd' if split daily, 'a' if no split
        if self.split_period == 'm':
            d = self.StartTime.replace(day=1).replace(hour=0).replace(minute=10)
            while d <= self.EndTime:
                self.months.append((d.year,d.month))
                d += relativedelta(months=1)
        elif self.split_period == 'd':
            d = self.StartTime.replace(hour=0).replace(minute=10)
            while d <= self.EndTime:
                self.months.append((d.year,d.month,d.day))
                d += timedelta(days=1)
        elif self.split_period == 'a':
            self.months.append(None)

    def create_pointdeparture_structure( self):
        self.PointDepartures = dict.fromkeys(self.region_names)
        for region in self.region_names:
            self.PointDepartures[region] = {}
            for tracer in self.species:
                self.PointDepartures[region][tracer] = {}
        for tracer in self.species:
            self.PointDepartures[tracer] = {}
            self.PointDepartures[tracer]['departure_class'] = self.rcf.get(tracer+'.departures.point.class')

    def get_class_from_name( self, class_name):
        _temp = __import__('PointObs', fromlist=[class_name])
        class_from_name = _temp.__dict__[class_name]
        return class_from_name
        #try:
        #    class_from_name = _temp.__dict__[class_name]
        #    return class_from_name
        #except KeyError:
        #    sys.stderr.write("Class %s not defined in %s.\n"%(class_name,'PointObs'))
        #    sys.exit()
        #except:
        #    sys.stderr.write("Unknown error importing %s\n"%class_name)
        #    sys.exit()

    def getTrackFileName(self, month_tuple=None):
        if month_tuple == None:
            fname = os.path.join(self.output_dir, 'point', 'point_output.nc4')
        elif len(month_tuple) == 2:
            fname = os.path.join(self.output_dir, 'point', 'point_output_%04i%02i.nc4'%month_tuple)
        elif len(month_tuple) == 3:
            fname = os.path.join(self.output_dir, 'point', 'point_output_%04i%02i%02i.nc4'%month_tuple)
        return fname

    def getMismatchFileName(self, month_tuple):
        if month_tuple == None:
            fname = os.path.join(self.output_dir, 'point', 'point_departures.nc4')
        elif len(month_tuple) == 2:
            fname = os.path.join(self.output_dir, 'point', 'point_departures_%04i%02i.nc4'%month_tuple)
        elif len(month_tuple) == 3:
            fname = os.path.join(self.output_dir, 'point', 'point_departures_%04i%02i%02i.nc4'%month_tuple)
        return fname

    def getInputFileName(self, month_tuple):
        indir_point = self.rcf.get('output.point.input.dir')
        if month_tuple == None:
            fname = os.path.join(indir_point, 'point_input.nc4')
        elif len(month_tuple) == 2:
            fname = os.path.join(indir_point, 'point_input_%04i%02i.nc4'%month_tuple)
        elif len(month_tuple) == 3:
            fname = os.path.join(indir_point, 'point_input_%04i%02i%02i.nc4'%month_tuple)
        return fname

    def getInputVars(self, gid, month_tuple):
        # needs to return the date components and observed mixing ratios from the input file
        sample_id = gid.variables['id'][:] - 1
        with Dataset(self.getInputFileName(month_tuple), 'r') as ifid:
            date_components = ifid.groups[self.tracer].variables['date_components'][:][sample_id]
            mixing_ratios = ifid.groups[self.tracer].variables['mixing_ratio'][:][sample_id]
            mixing_ratio_errs = ifid.groups[self.tracer].variables['mixing_ratio_error'][:][sample_id]
            lat = ifid.groups[self.tracer].variables['lat'][:][sample_id]
            lon = ifid.groups[self.tracer].variables['lon'][:][sample_id]
            alt = ifid.groups[self.tracer].variables['alt'][:][sample_id]
            # The variable 'time_window_length' may or may not exist
            if 'time_window_length' in ifid.groups[self.tracer].variables:
                time_window_length = ifid.groups[self.tracer].variables['time_window_length'][:][sample_id]
            else:
                time_window_length = zeros(len(alt), int32)

        return date_components, lat, lon, alt, mixing_ratios, mixing_ratio_errs, time_window_length

    def readPerturbations(self, region, tracer=None):
        """
        Reads the perturbations corresponding to a particular file from the file measurement_perturbations.nc4, and return
        them as a single array. If that particular file does not have any perturbations recorded, return an array of zeros.
        The only argument is the group id, from which the region name and track file name will be derived.
        """
        if tracer == None:
            tracer = self.tracer
        # what is the file name?
        file_name = os.path.basename(self.track_file)
        # where are the perturbations?
        perturb_file = os.path.join(self.output_dir, 'measurement_perturbations.nc4')
        with Dataset(perturb_file, 'r') as fid:
            pert = fid.groups['point_perturbations'].groups[file_name].groups[region].groups[tracer].variables['random_vector'][:]
        return pert

    def applyPointObs(self, month_tuple):
        mdm_obs_only = self.rcf.get('point.%s.only.obs.error'%self.tracer, 'bool', default=False)

        self.track_file = self.getTrackFileName(month_tuple)
        if not os.path.exists(self.track_file):
            return
        track_fid = Dataset(self.track_file, 'r')
        for region in self.region_names:
            group_data = track_fid.groups[region]
            if self.tracer in group_data.groups and len(group_data.groups[self.tracer].dimensions['samples']) > 0:
                date_components, lat, lon, alt, mixing_ratio, obs_error, tw_length = \
                    self.getInputVars(group_data.groups[self.tracer], month_tuple)

                mod_mixing = group_data.groups[self.tracer].variables['mixing_ratio'][:]
                mod_error = group_data.groups[self.tracer].variables['mixing_ratio_sigma'][:]

                self.PointDepartures[region][self.tracer]['dimensions'] = {'samples': len(mod_mixing)}
                self.PointDepartures[region][self.tracer]['variable_shapes'] = OrderedDict()

                self.PointDepartures[region][self.tracer]['date_components'] = int32(date_components)
                self.PointDepartures[region][self.tracer]['variable_shapes']['date_components'] = ('samples', 'idate')

                self.PointDepartures[region][self.tracer]['lat'] = float64(lat)
                self.PointDepartures[region][self.tracer]['variable_shapes']['lat'] = ('samples',)

                self.PointDepartures[region][self.tracer]['lon'] = float64(lon)
                self.PointDepartures[region][self.tracer]['variable_shapes']['lon'] = ('samples',)

                self.PointDepartures[region][self.tracer]['alt'] = float64(alt)
                self.PointDepartures[region][self.tracer]['variable_shapes']['alt'] = ('samples',)

                self.PointDepartures[region][self.tracer]['nsamples'] = group_data.groups[self.tracer].variables['nsamples'][:]
                self.PointDepartures[region][self.tracer]['variable_shapes']['nsamples'] = ('samples',)

                self.PointDepartures[region][self.tracer]['sampling_strategy'] = group_data.groups[self.tracer].variables['sampling_strategy'][:]
                self.PointDepartures[region][self.tracer]['variable_shapes']['sampling_strategy'] = ('samples',)

                self.PointDepartures[region][self.tracer]['time_window_length'] = int32(tw_length)
                self.PointDepartures[region][self.tracer]['variable_shapes']['time_window_length'] = ('samples',)

                self.PointDepartures[region][self.tracer]['model_error'] = float64(mod_error)
                self.PointDepartures[region][self.tracer]['variable_shapes']['model_error'] = ('samples',)

                self.PointDepartures[region][self.tracer]['obs_error'] = float64(obs_error)
                self.PointDepartures[region][self.tracer]['variable_shapes']['obs_error'] = ('samples',)

                # Should the total error be only the measurement error or the total error?
                if mdm_obs_only:
                    total_error = obs_error
                else:
                    total_error = sqrt(float64(obs_error)**2 + float64(mod_error)**2)

                self.PointDepartures[region][self.tracer]['error'] = float64(total_error)
                self.PointDepartures[region][self.tracer]['variable_shapes']['error'] = ('samples',)

                self.PointDepartures[region][self.tracer]['mismatch'] = float64(mod_mixing - mixing_ratio)
                self.PointDepartures[region][self.tracer]['variable_shapes']['mismatch'] = ('samples',)

                self.PointDepartures[region][self.tracer]['forcing'] = float64(mod_mixing - mixing_ratio)/total_error**2
                self.PointDepartures[region][self.tracer]['variable_shapes']['forcing'] = ('samples',)

        track_fid.close()

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

    def writePointDepartureFile(self, month_tuple=None):
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

#class Adjoint_dep(PointObs):
    #"""
    #For doing an adjoint test, it is not necessary to compute the adjoint forcings. An adjoint test verifies the equality
    #(y, Hx) = (x, H^T y), where H and H^T are the forward and adjoint transport operators. To do this, we start from a random
    #emission vector x and a random set of observations points where the model is to be sampled. Then in the track file, the
    #set of points given by 'mixing_ratio' constitutes the vector Hx. Then, we construct a set of random measurements y, of the
    #same length and structure as Hx, call them y. The inner product of y with Hx is the LHS. To get the RHS, we feed those
    #random measurements y as forcings (the variable 'forcing' in the point departures file) to the adjoint transport model
    #H^T. This results in an adjoint flux vector H^T y. The inner product of this with the previous random flux vector x gives
    #us the RHS.
    #"""

    #def applyPointObs(self, month_tuple=None):
        #"""
        #For every track file, create a departure/mismatch file, with bogus forcings
        #"""
        #self.track_file = self.getTrackFileName(month_tuple)
        #if not os.path.exists(self.track_file):
            #return
        #track_fid = Dataset(self.track_file, 'r')
        #for region in self.region_names:
            #group_data = track_fid.groups[region]
            #if self.tracer in group_data.groups and len(group_data.groups[self.tracer].dimensions['samples']) > 0:
                #date_components, lat, lon, alt, mixing_ratio, obs_error, tw_length = \
                    #self.getInputVars(group_data.groups[self.tracer], month_tuple)

                #n_obs = len(group_data.groups[self.tracer].dimensions['samples'])

                #self.PointDepartures[region][self.tracer]['dimensions'] = {'samples': n_obs}
                #self.PointDepartures[region][self.tracer]['variable_shapes'] = OrderedDict()

                #self.PointDepartures[region][self.tracer]['date_components'] = int32(date_components)
                #self.PointDepartures[region][self.tracer]['variable_shapes']['date_components'] = ('samples', 'idate')

                #self.PointDepartures[region][self.tracer]['lat'] = float64(lat)
                #self.PointDepartures[region][self.tracer]['variable_shapes']['lat'] = ('samples',)

                #self.PointDepartures[region][self.tracer]['lon'] = float64(lon)
                #self.PointDepartures[region][self.tracer]['variable_shapes']['lon'] = ('samples',)

                #self.PointDepartures[region][self.tracer]['alt'] = float64(alt)
                #self.PointDepartures[region][self.tracer]['variable_shapes']['alt'] = ('samples',)

                #self.PointDepartures[region][self.tracer]['nsamples'] = group_data.groups[self.tracer].variables['nsamples'][:]
                #self.PointDepartures[region][self.tracer]['variable_shapes']['nsamples'] = ('samples',)

                #self.PointDepartures[region][self.tracer]['sampling_strategy'] = group_data.groups[self.tracer].variables['sampling_strategy'][:]
                #self.PointDepartures[region][self.tracer]['variable_shapes']['sampling_strategy'] = ('samples',)

                #self.PointDepartures[region][self.tracer]['time_window_length'] = int32(tw_length)
                #self.PointDepartures[region][self.tracer]['variable_shapes']['time_window_length'] = ('samples',)

                #self.PointDepartures[region][self.tracer]['forcing'] = zeros(n_obs, float64)
                #self.PointDepartures[region][self.tracer]['variable_shapes']['forcing'] = ('samples',)

        #track_fid.close()

