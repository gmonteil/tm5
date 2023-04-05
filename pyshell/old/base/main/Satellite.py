#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True

import re, os, warnings, shutil
from pyshell.tmflex import rc
from pyshell.base.helper.Utilities import *
from dateutil.relativedelta import relativedelta
import numpy as np
from datetime import datetime
from pyshell.base.main.helper_functions import my_Dataset
from collections import defaultdict, OrderedDict
#from multiprocessing import Pool

from pyshell.base.main.meteo import LevelDefinitions

class ApplySatelliteObs(LevelDefinitions):

    def __init__(self, tm5):
        LevelDefinitions.__init__(self) # Get the pressure level boundaries
        # rcf from LevelDefinitions is overwritten by rcf from tm5
        self.rcf = tm5.rcf
        self.StartTime = tm5.StartTime
        self.EndTime   = tm5.EndTime
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.split_period = self.rcf.get('output.satellite.split.period', default='m') # 'm' if files are split monthly, 'd' if split daily
        self.output_dir = tm5.output_dir
        self.Optim_dict = tm5.Optim_dict
        self.months = []
        self.region_names = self.rcf.get('regions').split()
        # modification for multithreading
        #self.num_procs = self.rcf.get('par.maxthreads', 'int')
        # end modifications

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

        self.monthFileExists = dict(zip(self.months, np.zeros(len(self.months), bool)))
        for month in self.months:
            fname = self.getTrackfileName(month)
            if os.path.isfile(fname):
                self.monthFileExists[month] = True

        self.add_model_error = {}
        for tracer in self.species:
            self.add_model_error[tracer] = not self.rcf.get('sat.%s.only.obs.error'%tracer, 'bool', default=False)

    def getInputFileName(self, month_tuple):
        # month_tuple is a (year,month) tuple, such as (2009,12), or a (year,month,day) tuple, such as (2009,4,24)
        if self.split_period == 'm':
            fileName = os.path.join(self.rcf.get('output.satellite.output.directory'), self.rcf.get('output.satellite.ipfile.prefix')+"%04i%02i.nc4"%month_tuple)
        elif self.split_period == 'd':
            fileName = os.path.join(self.rcf.get('output.satellite.output.directory'), self.rcf.get('output.satellite.ipfile.prefix')+"%04i%02i%02i.nc4"%month_tuple)
        return fileName

    def getTrackfileName(self, month):
        # month is a (year,month) tuple, such as (2009,12), or a (year,month,day) tuple, such as (2009,4,24)
        if self.split_period == 'm':
            dummy_date = datetime(month[0], month[1], 1)
            fname   = os.path.join(self.output_dir, 'satellite', dummy_date.strftime("sat-track_%Y%m.nc4"))
        elif self.split_period == 'd':
            dummy_date = datetime(month[0], month[1], month[2])
            fname   = os.path.join(self.output_dir, 'satellite', dummy_date.strftime("sat-track_%Y%m%d.nc4"))
        return fname

    def getDepartureFileName(self, month):
        # month is a (year,month) tuple, such as (2009,12), or including the day, such as (2009,12,10)
        if self.split_period == 'm':
            dummy_date = datetime(month[0], month[1], 1)
            fname   = os.path.join(self.output_dir, 'satellite', dummy_date.strftime("sat-track_departures_%Y%m.nc4"))
        elif self.split_period == 'd':
            dummy_date = datetime(month[0], month[1], month[2])
            fname   = os.path.join(self.output_dir, 'satellite', dummy_date.strftime("sat-track_departures_%Y%m%d.nc4"))
        return fname

    def create_satdeparture_structure(self):
        self.SatDepartures = dict.fromkeys(self.region_names)
        for region in self.region_names:
            self.SatDepartures[region] = {}
            for tracer in self.species:
                self.SatDepartures[region][tracer] = {}

        optim_bias_param = False
        n_bias = 0
        for tracer in self.species:
            self.SatDepartures[tracer] = {}
            self.SatDepartures[tracer]['departure_class'] = self.rcf.get(tracer+'.departures.sat.class')
            # in case we're optimizing bias parameters, create further entries for them
            self.SatDepartures[tracer]['optimize bias'] = self.rcf.get("%s.optimize.sat.bias.parameter"%tracer, 'bool', default=False)
            optim_bias_param = optim_bias_param or self.SatDepartures[tracer]['optimize bias']
            if self.SatDepartures[tracer]['optimize bias']:
                self.SatDepartures[tracer]['bias'] = {'n_param': self.rcf.get("%s.sat.bias.num_params"%tracer, 'int')}
                n_bias += self.SatDepartures[tracer]['bias']['n_param']
                self.SatDepartures[tracer]['bias']['param'] = np.zeros(self.SatDepartures[tracer]['bias']['n_param'], np.float64)
                self.SatDepartures[tracer]['bias']['adj_param'] = np.zeros(self.SatDepartures[tracer]['bias']['n_param'], np.float64)
                self.SatDepartures[tracer]['bias']['param_err'] = np.zeros(self.SatDepartures[tracer]['bias']['n_param'], np.float64)
        self.optim_bias_param = optim_bias_param

        # At this point, the bias parameters, their adjoints, and their errors are set to zeros for all tracers

        if optim_bias_param:
            parameter_scratch_file = os.path.join(self.output_dir, self.rcf.get('parameter.scratch.filename'))
            with my_Dataset(parameter_scratch_file, 'r') as fid:
                for tracer in self.species:
                    if self.SatDepartures[tracer]['optimize bias']:
                        # assume for the moment that all parameters for a tracer are satellite biases
                        self.SatDepartures[tracer]['bias']['param'][:] = fid.groups[tracer].variables['params'][:]

            apri_file_name = self.rcf.get('parameter.apri.filename')
            with my_Dataset(apri_file_name, 'r') as fid:
                for tracer in self.species:
                    if self.SatDepartures[tracer]['optimize bias']:
                        self.SatDepartures[tracer]['bias']['param_err'][:] = fid.groups[tracer].variables['param_errs'][:]

            #self.bias_parameter = np.loadtxt(parameter_scratch_file, dtype=np.float64)
            #dapri_file_name = self.rcf.get('parameter.prior.error.filename')
            #self.bias_error = np.loadtxt(dapri_file_name, dtype=np.float64)
            #self.AdjBiasParameters = zeros_like(self.bias_error)
            #if len(self.bias_parameter) != n_bias:
                #raise ValueError("Number of bias parameters does not match size of vector in %s"%parameter_scratch_file)
            #if len(self.bias_error) != n_bias:
                #raise ValueError("Number of bias parameters does not match size of vector in %s"%dapri_file_name)

        ## Now, self.bias_parameter is set to the current value of the bias parameters, self.bias_error is set to the
        ## apri errors of those parameters, and self.AdjBiasParameters is set to a vector of zeros

    def get_class_from_name( self, class_name):
        _temp = __import__('Satellite', fromlist=[class_name])
        try:
            class_from_name = _temp.__dict__[class_name]
            return class_from_name
        except KeyError:
            sys.stderr.write("Class %s not defined in %s.\n"%(class_name,'Satellite'))
            sys.exit()
        except:
            sys.stderr.write("Unknown error importing %s\n"%class_name)
            sys.exit()

    def ApplyAveragingKernel(self, i_obs):
        """
        This ApplyAveragingKernel is for a very simple satellite instrument, which measures the pressure-weighted total
        column average from the ground to the TOA. This should be superceded by a satellite class' own ApplyAveragingKernel
        method. Alternatively, a satellite class can also overwrite the MapAveragingKernel method, which calls this method.
        For this method to function, self.tracer, self.only_obs_error and self.instr_num must be set.
        """
        model_profile = self.profiles[i_obs,:]
        std_model_profile = self.std_profiles[i_obs,:]
        model_psurf = self.psurf_model[i_obs]
        index = self.ind[i_obs] # this points to common variables in the input file
        instr_index = self.instr_ind[i_obs] # this points to variables specific to each instrument

        # What is the instrument?
        instrument = self.instrument[index]
        if instrument != self.instr_num:
            raise ValueError('Observations are from instrument %i, not from this instrument (%i)'%(instrument, self.instr_num))

        # What are the pressure levels?
        model_pres_levels = (self.AT + self.BT*model_psurf).astype(np.float64)
        pres_weights = -np.diff(model_pres_levels)/model_psurf

        meas_tc = self.input_vars[instrument]['column_mixing'][instr_index]
        meas_dtc = self.input_vars[instrument]['sigma_column_mixing'][instr_index]**2 # the variance, not the standard deviation

        model_tc = np.sum(model_profile * pres_weights)
        mod_dtc = self.modeled_column_avg_err[i_obs]**2
        if self.only_obs_error:
            tot_dtc = meas_dtc
        else:
            tot_dtc = mod_dtc+meas_dtc # change this, possibly

        self.mismatches['input_index'].append(index)
        self.mismatches['modeled_column'].append(model_tc)
        self.mismatches['measured_column'].append(meas_tc)
        self.mismatches['sigma_column'].append(tot_dtc)
        self.mismatches['track_index'].append(i_obs)
        self.mismatches['instrument_index'].append(instr_index)

        deps = (model_tc - meas_tc)/tot_dtc # scalar
        J = 0.5 * deps * (model_tc - meas_tc)
        self.J_obs_satellite += J
        self.departures[i_obs,:] = pres_weights*deps
        return True

    def OpenInputFile(self, month_tuple):
        # The argument month_tuple can be either a (year,month) tuple or a (year,month,day) tuple.
        if self.split_period == 'm':
            dummyDate = datetime(month_tuple[0],month_tuple[1],1)
            fileName = os.path.join(self.rcf.get('output.satellite.output.directory'), self.rcf.get('output.satellite.ipfile.prefix')+dummyDate.strftime("%Y%m.nc4"))
        elif self.split_period == 'd':
            dummyDate = datetime(month_tuple[0],month_tuple[1],month_tuple[2])
            fileName = os.path.join(self.rcf.get('output.satellite.output.directory'), self.rcf.get('output.satellite.ipfile.prefix')+dummyDate.strftime("%Y%m%d.nc4"))

        # For some special purposes, it is necessary set the measured mole fractions to zero
        zero_mix = self.rcf.get('output.satellite.mole_frac.zero', 'bool', default=False)
        zero_mix = self.rcf.get('output.satellite.%s.mole_frac.zero'%self.tracer, 'bool', default=zero_mix)

        if os.path.isfile(fileName):
            with my_Dataset(fileName, 'r') as fid:
                gid = fid.groups[self.tracer]
                self.latitude = gid.variables['latitude'][:]
                self.longitude = gid.variables['longitude'][:]
                self.sample_times = gid.variables['cdate'][:]
                self.instrument = gid.variables['instrument'][:]

                self.common_input_vars = ['latitude', 'longitude', 'sample_times', 'instrument']

                # Input variables are specific to each instrument, so loop over instruments
                self.input_vars = {}
                try:
                    instrument_nums = gid.variables['instrument_nums'][:]
                except KeyError:
                    instrument_nums = np.unique(gid.variables['instrument'][:])

                for instr in instrument_nums:
                    igid = gid.groups['instrument_%02i'%instr]
                    self.input_vars[instr] = {}

                    for var_name, var_value in igid.variables.items():
                        self.input_vars[instr][var_name] = var_value[:]

                    if zero_mix:
                        self.input_vars[instr]['column_mixing'] = np.zeros_like(self.input_vars[instr]['column_mixing'])
                        print 'Set observed satellite mole fractions for %s to zeros'%self.tracer

    def cleanInputData(self):
        try:
            del self.sample_times, self.latitude, self.longitude, self.instrument, self.input_vars
        except AttributeError:
            print 'Input variables not read in, so nothing to clean'
        try:
            del self.profiles, self.std_profiles, self.psurf_model, self.ind, self.instr_ind
            del self.nsamples, self.sampling_strategy
        except AttributeError:
            print 'Track variables not read in, so nothing to clean'

    def checkDepartures(self, region, month):
        # Check whether there are NaNs or Infs in the departure-related fields
        if self.split_period == 'm':
            dummy_date = datetime(month[0], month[1], 1)
            date_string = dummy_date.strftime("%b %Y")
        elif self.split_period == 'd':
            dummy_date = datetime(month[0], month[1], month[2])
            date_string = dummy_date.strftime("%d %b %Y")

        if np.any(np.isnan(self.departures)) or np.any(np.isinf(self.departures)):
            # find which index is causing trouble
            nan_values = np.any(np.isnan(self.departures), axis=1)
            inf_values = np.any(np.isinf(self.departures), axis=1)
            nan_values = np.where(nan_values)[0]
            inf_values = np.where(inf_values)[0]
            if len(nan_values) > 0:
                nan_string = ', '.join(["%i"%i for i in nan_values])
                warnings.warn("The following indices have NaN values: %s"%nan_string, RuntimeWarning)
                for i in nan_values:
                    write_line = "departure(%i) = ["%i + '  '.join(["%f"%f for f in self.departures[i]]) + "]\n"
                    sys.stderr.write(write_line)
            if len(inf_values) > 0:
                inf_string = ', '.join(["%i"%i for i in inf_values])
                warnings.warn("The following indices have Inf values: %s"%inf_string, RuntimeWarning)
                for i in inf_values:
                    write_line = "departure(%i) = ["%i + '  '.join(["%f"%f for f in self.departures[i]]) + "]\n"
                    sys.stderr.write(write_line)
            raise RuntimeError("Invalid departure value for tracer %s in region %s for %s"%(self.tracer, region, date_string))

    def applySatelliteObs(self, month):
        """
        Although this routine is in the parent class, it will only be called by tracer-specific classes. As long
        as the following syntax for generating departures is valid,

        self.valid_deps = self.MapAveragingKernel()

        which in turn contains the single line

        return np.array(map(self.ApplyAveragingKernel, range(self.n_obs)))

        If tracer-specific classes write their ApplyAveragingKernel routines to conform to this standard, then
        this should work. If a tracer-specific class wants to have a different interface to ApplyAveragingKernel,
        e.g., the IASI CO class, then it needs to have their own MapAveragingKernel or applySatelliteObs routine.
        """
        # the bias parameters should not be modified here, so we can create a local copy from self.SatDepartures
        #if self.SatDepartures[self.tracer]['optimize bias']:
            #self.bias_parameter = self.SatDepartures[self.tracer]['bias']['param']
            #self.bias_error = self.SatDepartures[self.tracer]['bias']['param_err']

        if not self.monthFileExists[month]:
            if self.split_period == 'm':
                dummy_date = datetime(month[0], month[1], 1)
                print dummy_date.strftime("    No track file for %b %Y")
            elif self.split_period == 'd':
                dummy_date = datetime(month[0], month[1], month[2])
                print dummy_date.strftime("    No track file for %e %b %Y")
            # Set the number of observations to zero
            for region in self.region_names:
                self.SatDepartures[region][self.tracer]['n_obs'] = 0
                self.SatDepartures[region][self.tracer]['J_obs_satellite'] = np.zeros(1, np.float64)
            return

        self.OpenInputFile(month)
        fname = self.getTrackfileName(month)
        for region in self.region_names:
            self.OpenSatelliteTrackFile(fname,region)
            self.SatDepartures[region][self.tracer]['n_obs'] = self.n_obs
            if self.n_obs == 0:
                continue # go to next region
                # print 'No observations for region ', region, ' for tracer ', self.tracer

            self.departures = np.zeros((self.n_obs,self.n_lev))
            self.mismatches = defaultdict(list)
            self.J_obs_satellite = 0.0
            # modifications for multiprocessing
            ## save OMP_NUM_THREADS
            #orig_omp_threads = os.environ['OMP_NUM_THREADS']
            ## set OMP_NUM_THREADS to 1
            #os.environ['OMP_NUM_THREADS'] = '1'
            #pool = Pool(processes=self.num_procs)
            #valid_deps = pool.map_async(self.ApplyAveragingKernel, range(self.n_obs))
            #valid_deps = valid_deps.get()
            ## set OMP_NUM_THREADS back to its original value
            #os.environ['OMP_NUM_THREADS'] = orig_omp_threads
            #self.valid_deps = np.array(valid_deps)
            # end modifications
            self.valid_deps = self.MapAveragingKernel()

            self.checkDepartures(region, month)

            # Write departures for this region to structure
            self.SatDepartures[region][self.tracer]['dimensions'] = {'n_obs': self.n_obs, 'n_lev': self.n_lev}
            self.SatDepartures[region][self.tracer]['variable_shapes'] = OrderedDict()

            self.SatDepartures[region][self.tracer]['J_obs_satellite'] = np.array([self.J_obs_satellite], np.float64)
            self.SatDepartures[region][self.tracer]['variable_shapes']['J_obs_satellite'] = ('nJ',)

            self.SatDepartures[region][self.tracer]['valid_departure'] = self.valid_deps.astype(np.int8)
            self.SatDepartures[region][self.tracer]['variable_shapes']['valid_departure'] = ('n_obs',)

            sort_order = np.argsort(self.mismatches['track_index'])
            for k,v in self.mismatches.items():
                self.mismatches[k] = np.array(v)[sort_order]

            # Are we assimilating profiles?
            is_profile = False
            if self.mismatches['modeled_column'].ndim == 2:
                n_channel = self.mismatches['modeled_column'].shape[1]
                is_profile = True

            if is_profile:
                self.SatDepartures[region][self.tracer]['dimensions']['n_channel'] = n_channel
                col_dims = ('n_obs', 'n_channel')
            else:
                col_dims = ('n_obs',)

            self.SatDepartures[region][self.tracer]['input_index'] = self.mismatches['input_index'].astype(np.int32)
            self.SatDepartures[region][self.tracer]['variable_shapes']['input_index'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['modeled_column'] = self.mismatches['modeled_column'].astype(np.float64)
            self.SatDepartures[region][self.tracer]['variable_shapes']['modeled_column'] = col_dims

            self.SatDepartures[region][self.tracer]['measured_column'] = self.mismatches['measured_column'].astype(np.float64)
            self.SatDepartures[region][self.tracer]['variable_shapes']['measured_column'] = col_dims

            self.SatDepartures[region][self.tracer]['sigma_column'] = np.sqrt(self.mismatches['sigma_column'].astype(np.float64))
            self.SatDepartures[region][self.tracer]['variable_shapes']['sigma_column'] = col_dims

            self.SatDepartures[region][self.tracer]['instrument_index'] = self.mismatches['instrument_index'].astype(np.int32)
            self.SatDepartures[region][self.tracer]['variable_shapes']['instrument_index'] = ('n_obs',)

            if 'bias_correction' in self.mismatches:
                self.SatDepartures[region][self.tracer]['bias_correction'] = self.mismatches['bias_correction'].astype(np.float64)
                self.SatDepartures[region][self.tracer]['variable_shapes']['bias_correction'] = col_dims

            # write variables from the input file
            relevant_indices = self.mismatches['input_index'].astype(np.int32)
            self.SatDepartures[region][self.tracer]['idate'] = self.sample_times[relevant_indices].astype(np.int16)
            self.SatDepartures[region][self.tracer]['variable_shapes']['idate'] = ('n_obs', 'date_components')

            self.SatDepartures[region][self.tracer]['lat'] = self.latitude[relevant_indices].astype(np.float64)
            self.SatDepartures[region][self.tracer]['variable_shapes']['lat'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['lon'] = self.longitude[relevant_indices].astype(np.float64)
            self.SatDepartures[region][self.tracer]['variable_shapes']['lon'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['instrument'] = self.instrument[relevant_indices].astype(np.int8)
            self.SatDepartures[region][self.tracer]['variable_shapes']['instrument'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['departures'] = self.departures
            self.SatDepartures[region][self.tracer]['variable_shapes']['departures'] = ('n_obs', 'n_lev')

            self.SatDepartures[region][self.tracer]['model_minus_obs'] = (self.mismatches['modeled_column'] - self.mismatches['measured_column']).astype(np.float64)
            self.SatDepartures[region][self.tracer]['variable_shapes']['model_minus_obs'] = col_dims

            self.SatDepartures[region][self.tracer]['nsamples'] = self.nsamples
            self.SatDepartures[region][self.tracer]['variable_shapes']['nsamples'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['total_weight'] = self.total_weight
            self.SatDepartures[region][self.tracer]['variable_shapes']['total_weight'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['sampling_strategy'] = self.sampling_strategy
            self.SatDepartures[region][self.tracer]['variable_shapes']['sampling_strategy'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['raw_modeled_column'] = self.modeled_column_avg
            self.SatDepartures[region][self.tracer]['variable_shapes']['raw_modeled_column'] = ('n_obs',)

            self.SatDepartures[region][self.tracer]['raw_modeled_column_err'] = self.modeled_column_avg_err
            self.SatDepartures[region][self.tracer]['variable_shapes']['raw_modeled_column_err'] = ('n_obs',)

        # we need to clear the memory of all data for this month/period
        self.cleanInputData()

    def MapAveragingKernel(self):
        return np.array(map(self.ApplyAveragingKernel, range(self.n_obs)))

    def OpenSatelliteTrackFile(self,inputFileName,region):
        # print 'Track file %s opened for region %s'%(os.path.basename(inputFileName), region)
        self.n_obs = 0 # default
        if os.path.isfile(inputFileName):
            with my_Dataset(inputFileName,'r') as f:
                if region in f.groups:
                    # Group Dimensions
                    g = f.groups[region]
                    self.n_lev = len(g.dimensions['n_lev'])
                    if self.tracer in g.groups:
                        t = g.groups[self.tracer]
                        self.n_obs = len(t.dimensions['n_obs'])
                        # Variables
                        self.profiles  = t.variables['profiles'][:] if self.n_obs > 0 else None
                        self.std_profiles = t.variables['std_profiles'][:] if self.n_obs > 0 else None
                        self.psurf_model = t.variables['psurf'][:] if self.n_obs > 0 else None
                        self.ind = t.variables['input_positions'][:] - 1 if self.n_obs > 0 else None # indices are hereby converted to the pythonic numbering system
                        self.instr_ind = t.variables['instrument_positions'][:] - 1 if self.n_obs > 0 else None # indices are hereby converted to the pythonic numbering system
                        #self.measured_total_col = t.variables['column_mixing'][:] if self.n_obs > 0 else None
                        #self.sigma_measured_total_col = t.variables['sigma_column_mixing'][:] if self.n_obs > 0 else None
                        self.nsamples = t.variables['nsamples'][:] if self.n_obs > 0 else None
                        self.sampling_strategy = t.variables['sampling_strategy'][:] if self.n_obs > 0 else None
                        self.total_weight = t.variables['total_weight'][:] if self.n_obs > 0 else None
                        self.modeled_column_avg = t.variables['model_column'][:] if self.n_obs > 0 else None
                        self.modeled_column_avg_err = t.variables['sigma_model_column'][:] if self.n_obs > 0 else None

                        # What are track variables common to all instruments?
                        self.common_track_vars = ['n_obs', 'profiles', 'std_profiles', 'psurf_model', 'ind', 'instr_ind',
                            'nsamples', 'sampling_strategy', 'total_weight', 'modeled_column_avg', 'modeled_column_avg_err']

                        if 'meteo' in t.groups:
                            self.mole_frac_water = (28.94/18.0152) * t.groups['meteo'].variables['specific_humidity'][:] if self.n_obs > 0 else None
                            self.gph = t.groups['meteo'].variables['geopotential_height'][:] if self.n_obs > 0 else None
                            self.temperature = t.groups['meteo'].variables['temperature'][:] if self.n_obs > 0 else None

                            self.common_track_vars.extend(['mole_frac_water', 'gph', 'temperature'])

    def writeDeparturesFile(self, month):
        DepartureFile = self.getDepartureFileName(month)
        if os.path.exists(DepartureFile):
            os.remove(DepartureFile)
        total_n_obs = 0
        for region in self.region_names:
            for tracer in self.species:
                total_n_obs += self.SatDepartures[region][tracer]['n_obs']

        if total_n_obs == 0:
            return

        with my_Dataset(DepartureFile, 'w') as fid:
            fid.createDimension('date_components', 6)
            fid.createDimension('nJ', 1)
            J_obs_satellite = 0.
            for region in self.region_names:
                region_n_obs = 0
                for tracer in self.species:
                    region_n_obs += self.SatDepartures[region][tracer]['n_obs']
                if region_n_obs != 0:
                    rgrp = fid.createGroup(region)
                    rgrp.createDimension('tracer', self.ntracer)
                    for tracer in self.species:
                        tracer_n_obs = self.SatDepartures[region][tracer]['n_obs']
                        if tracer_n_obs != 0:
                            J_obs_satellite += self.SatDepartures[region][tracer]['J_obs_satellite'][0]
                            tgrp = rgrp.createGroup(tracer)
                            for dim_name, dim_len in self.SatDepartures[region][tracer]['dimensions'].items():
                                tgrp.createDimension(dim_name, dim_len)
                            for var_name, var_shape in self.SatDepartures[region][tracer]['variable_shapes'].items():
                                try:
                                    var_value = self.SatDepartures[region][tracer][var_name]
                                    v = tgrp.createVariable(var_name, var_value.dtype, var_shape)
                                    v[:] = var_value
                                    # delete the data from memory
                                    del self.SatDepartures[region][tracer][var_name]
                                except:
                                    print 'var_name = ', var_name
                                    print 'array shape = ', var_value.shape
                                    print 'shape in file = ', v.shape
                                    print 'dimensions in file = ', v.dimensions
                                    raise

            v = fid.createVariable('J_obs_satellite', np.float64, ('nJ',))
            v[:] = J_obs_satellite

            if self.optim_bias_param:
                bias_grp = fid.createGroup('bias_params')
                for tracer in self.species:
                    if self.SatDepartures[tracer]['optimize bias']:
                        bias_tr_grp = bias_grp.createGroup(tracer)
                        n_param = self.SatDepartures[tracer]['bias']['n_param']
                        bias_tr_grp.createDimension('n_param', n_param)
                        v = bias_tr_grp.createVariable('params', np.float64, ('n_param',))
                        v[:] = self.SatDepartures[tracer]['bias']['param']
                        v = bias_tr_grp.createVariable('adj_params', np.float64, ('n_param',))
                        v[:] = self.SatDepartures[tracer]['bias']['adj_param']
                        v = bias_tr_grp.createVariable('param_err', np.float64, ('n_param',))
                        v[:] = self.SatDepartures[tracer]['bias']['param_err']

    def WriteAdjointBiasParameters(self):
        paramFile = os.path.join(self.output_dir, self.rcf.get('adjoint.parameter.filename'))
        if os.path.isfile(paramFile):
            os.remove(paramFile)

        #adj_bias_params = {}
        #for tracer in self.species:
            #if self.SatDepartures[tracer]['optimize bias']:
                #adj_bias_params[tracer] = self.SatDepartures[tracer]['bias']['adj_param']

        # Need to sum up the adjoint forcing from all departure files
        adj_bias_params = {}
        for tracer in self.species:
            if self.SatDepartures[tracer]['optimize bias']:
                n_param = self.SatDepartures[tracer]['bias']['n_param']
                adj_bias_params[tracer] = np.zeros(n_param, np.float64)

        for month_tuple in self.months:
            dep_file = self.getDepartureFileName(month_tuple)
            if not os.path.exists(dep_file):
                continue

            with my_Dataset(dep_file, 'r') as fid:
                gid = fid.groups['bias_params']

                for tracer in adj_bias_params.keys():
                    if tracer not in gid.groups.keys():
                        continue

                    adj_bias_params[tracer][:] = adj_bias_params[tracer][:] + gid.groups[tracer].variables['adj_params'][:]

        with my_Dataset(paramFile, 'w') as fid:
            # for now, all bias parameters are satellite related
            for tracer, adj_params in adj_bias_params.items():
                gid = fid.createGroup(tracer)
                gid.createDimension('n_param', len(adj_params))
                v = gid.createVariable('adj_params', np.float64, ('n_param',))
                v[:] = adj_params

    def SatelliteCost(self):
        J_obs = 0.0
        for month in self.months:
            DepartureFile = self.getDepartureFileName(month)
            if os.path.isfile(DepartureFile):
                with my_Dataset(DepartureFile,'r') as fid:
                    J_obs += fid.variables['J_obs_satellite'][0]
        return J_obs

class Multi_Sat_per_Tracer(ApplySatelliteObs):

    def __init__(self, *args, **kwargs):
        super(Multi_Sat_per_Tracer, self).__init__(*args, **kwargs)

    def OpenInputFile(self, month):
        # This should call the generic OpenInputFile, then create pointers from self.object_dict to the stuff read in
        super(Multi_Sat_per_Tracer, self).OpenInputFile(month)
        for instr in self.instruments:
            for var_name in self.common_input_vars:
                setattr(self.object_dict[instr], var_name, getattr(self, var_name))
            self.object_dict[instr].input_vars = self.input_vars

    def OpenSatelliteTrackFile(self, fname, region):
        # This should call the generic OpenSatelliteTrackFile, then create pointers from self.object_dict to the stuff read in
        super(Multi_Sat_per_Tracer, self).OpenSatelliteTrackFile(fname, region)
        for instr in self.instruments:
            for var_name in self.common_track_vars:
                setattr(self.object_dict[instr], var_name, getattr(self, var_name))

    def MapAveragingKernel(self):
        track_indices = np.arange(self.n_obs) # 0, 1, 2, ... total number of obs across all satellites from the track file, for a tracer and region
        input_idx = self.ind[track_indices] # the position in the input file, across all instruments for this tracer
        instrument = self.instrument[input_idx] # self.instrument has been read from the input file
        # 'instrument' is an array inditicating which obs in the track file belongs to which instrument

        ret_arr = np.zeros(self.n_obs, bool)

        for instr in self.instruments:
            # print 'Calculating departures for ', instr
            # point from self.object_dict to self, so that mismatches and departures could be added to the global variables
            self.object_dict[instr].departures = self.departures
            self.object_dict[instr].mismatches = self.mismatches
            self.object_dict[instr].J_obs_satellite = 0.0

            optim_bias_param = self.rcf.get('%s.%s.optim.bias_param'%(self.tracer, instr), 'bool', default=False)
            if optim_bias_param:
                self.object_dict[instr].SatDepartures[self.tracer]['bias'] = self.SatDepartures[self.tracer]['bias']

            # which obs belong to this instrument?
            instr_num = self.object_dict[instr].instr_num
            # are there any obs in the current track file/region/tracer that belongs to this instrument?
            instr_idx = (instrument == instr_num)
            if instr_idx.sum() > 0:
                ret_arr[instr_idx] = self.object_dict[instr].ApplyAllAveragingKernel(track_indices[instr_idx])
                # track_indices[instr_idx] are the indices in the track file that belong to this instrument

            # add to the satellite cost function
            self.J_obs_satellite += self.object_dict[instr].J_obs_satellite

        #print
        #for instr in self.instruments:
            #deps = self.object_dict[instr].departures
            #mism = self.object_dict[instr].mismatches['modeled_column']
            #J = self.object_dict[instr].J_obs_satellite
            #print 'For instrument ', instr
            #print 'departures shape = ', deps.shape, ', max val = ', np.abs(deps).max()
            #print 'mismatches shape = ', len(mism), ', max val = ', np.abs(mism).max()
            #print 'J satellite = ', J
            #print 'self.departures.shape = ', self.departures.shape
            #print 'self.mismatches.shape = ', len(self.mismatches['modeled_column'])
            #print

        return ret_arr
