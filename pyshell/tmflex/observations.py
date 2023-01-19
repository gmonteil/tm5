#!/bin/env python

import os
from collections import OrderedDict

import numpy as np

class Observations:
    def __init__(self, ti, tf, rcf):
        self.ti = ti
        self.tf = tf
        self.rcf = rcf
        self.species = [s.strip() for s in rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        # self.station_name_list = []
        # self.obspack_cats = {}
        # self.sat_split_period = rcf.get('output.satellite.split.period', default='d')
        self.point_split_period = rcf.get('output.point.split.period', default='a')
        # self.subdir_tag = ''
        self.PointObservation = {}
        for tracer in self.species:
            self.PointObservation[tracer] = {}

    def SetupPointObs(self, tracer, obsobj):
        """ obsobj must be a python object with a method 'get(tracer, startTime, endTime)', which returns a dictionary with the following keys:
            - n : the number of observations
            - latitudes
            - longitudes
            - altitudes
            - mixing_ratio
            - mixing_ratio_err (only if the observations are to be assimilated)
            - sampling_strategy
            - date_components
            - station_id
        """
        assimilated = self.rcf.get('%s.point.assimilate' % tracer, 'bool', default=True)
        if not assimilated:
            unassim_mdm = self.rcf.get('%s.point.unassimilate.mdm' % tracer, 'float')
        obs = obsobj.get(tracer, self.ti, self.tf)
        self.PointObservation[tracer]['dimensions'] = {'id': obs['n']}
        self.PointObservation[tracer]['variable_shapes'] = OrderedDict()
        self.PointObservation[tracer]['attr_dict'] = OrderedDict()
        self.PointObservation[tracer]['variable_attrs'] = OrderedDict()
        self.PointObservation[tracer]['id'] = np.arange(1, obs['n'] + 1, dtype=np.int32)
        self.PointObservation[tracer]['variable_shapes']['id'] = ('id',)
        self.PointObservation[tracer]['lat'] = np.float64(obs['latitudes'])
        self.PointObservation[tracer]['variable_shapes']['lat'] = ('id',)
        self.PointObservation[tracer]['lon'] = np.float64(obs['longitudes'])
        self.PointObservation[tracer]['variable_shapes']['lon'] = ('id',)
        self.PointObservation[tracer]['alt'] = obs['altitudes']
        self.PointObservation[tracer]['variable_shapes']['alt'] = ('id',)
        self.PointObservation[tracer]['mixing_ratio'] = obs['mixing_ratio']
        self.PointObservation[tracer]['variable_shapes']['mixing_ratio'] = ('id',)
        if assimilated:
            self.PointObservation[tracer]['mixing_ratio_error'] = obs['mixing_ratio_err']
        else:
            self.PointObservation[tracer]['mixing_ratio_error'] = unassim_mdm * np.ones(n_obs, np.float64)
            self.PointObservation[tracer]['variable_attrs']['mixing_ratio_error'] = [
                ('comment', 'These data will not be assimilated')]
        self.PointObservation[tracer]['variable_shapes']['mixing_ratio_error'] = ('id',)
        self.PointObservation[tracer]['date_components'] = np.array([x.timetuple()[:6] for x in obs['dates']],
                                                                    dtype=np.int16)
        self.PointObservation[tracer]['variable_shapes']['date_components'] = ('id', 'idate')
        self.PointObservation[tracer]['sampling_strategy'] = np.array(obs['sampling_strategy'], dtype=np.int16)
        self.PointObservation[tracer]['variable_shapes']['sampling_strategy'] = ('id',)
        self.PointObservation[tracer]['station_id'] = np.array(obs['station_id'], dtype=np.int32)
        self.PointObservation[tracer]['variable_shapes']['station_id'] = ('id',)
        #self.PointObservation[tracer]['obspack_num'] = obs['obspack_id']
        #self.PointObservation[tracer]['variable_shapes']['obspack_num'] = ('id',)

    def writePointFile(self):
        from netCDF4 import Dataset
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
        print "Wrote %s" % point_file
