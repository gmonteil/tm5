#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
sys.dont_write_bytecode = True

from pyshell.tmflex import rc
from numpy import *
from datetime import datetime, timedelta
from pyshell.base.helper.Utilities import checkDir
from netCDF4 import Dataset, OrderedDict

from pyshell.base.main.PointObs import PointObs #, Adjoint_dep

#class Adjoint_CO2_dep(Adjoint_dep):

    #def __init__(self, tm5Obj):
        #Adjoint_dep.__init__(self, tm5Obj)
        #self.tracer = 'CO2'

class AMDAR_OSSE_assim(PointObs):

    def __init__(self, tm5Obj):
        PointObs.__init__(self, tm5Obj)
        self.tracer = 'CO2'

class CT2013_CO2_dep(PointObs):

    def __init__(self, tm5Obj):
        PointObs.__init__(self, tm5Obj)
        self.tracer = 'CO2'
        self.only_obs_error = self.rcf.get('point.%s.only.obs.error'%self.tracer, 'bool', default=False)

    def applyPointObs(self, month_tuple=None):
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

                # This is a modification for weighting the samples by variable time steps, do not use this for running code!
                try:
                    self.PointDepartures[region][self.tracer]['total_weight'] = group_data.groups[self.tracer].variables['total_weight'][:]
                    self.PointDepartures[region][self.tracer]['variable_shapes']['total_weight'] = ('samples',)
                except KeyError:
                    # Old code (pre-cycle 3) will not generate a 'total_weight' variable in the track file
                    pass

                self.PointDepartures[region][self.tracer]['sampling_strategy'] = group_data.groups[self.tracer].variables['sampling_strategy'][:]
                self.PointDepartures[region][self.tracer]['variable_shapes']['sampling_strategy'] = ('samples',)

                self.PointDepartures[region][self.tracer]['time_window_length'] = int32(tw_length)
                self.PointDepartures[region][self.tracer]['variable_shapes']['time_window_length'] = ('samples',)

                self.PointDepartures[region][self.tracer]['model_error'] = float64(mod_error)
                self.PointDepartures[region][self.tracer]['variable_shapes']['model_error'] = ('samples',)

                self.PointDepartures[region][self.tracer]['obs_error'] = float64(obs_error)
                self.PointDepartures[region][self.tracer]['variable_shapes']['obs_error'] = ('samples',)

                if self.only_obs_error:
                    total_error = obs_error
                else:
                    total_error = sqrt(obs_error**2 + mod_error**2)

                self.PointDepartures[region][self.tracer]['error'] = float64(total_error)
                self.PointDepartures[region][self.tracer]['variable_shapes']['error'] = ('samples',)

                self.PointDepartures[region][self.tracer]['mismatch'] = float64(mod_mixing - mixing_ratio)
                self.PointDepartures[region][self.tracer]['variable_shapes']['mismatch'] = ('samples',)

                self.PointDepartures[region][self.tracer]['forcing'] = float64(mod_mixing - mixing_ratio)/total_error**2
                self.PointDepartures[region][self.tracer]['variable_shapes']['forcing'] = ('samples',)

        track_fid.close()
