#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True

import logging
logging.basicConfig(datefmt='%d %b %Y %H:%M:%S', format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

import re, os, shutil, hashlib, glob
from random import choice
import string
from pyshell.tmflex import rc
from numpy import *
#from Utilities import *
from netCDF4 import Dataserfrom datetime import datetime
import subprocess as sbp
import cPickle as pickle
# import dill as pickle

try :
    rcf = rc.RcFile(os.environ['pyshell.rc'])
    # Are we going to monitor all instruments?
    #monitor_instruments = rcf.get('optimize.monitor.outputs', 'bool', default=False)
    # The postprocess module
    if rcf.get('optimize.postprocess', 'bool', default=False):
        pst_module = rcf.get('my.postprocess.class', default='Postprocess')
        _temp = __import__('Postprocess', fromlist=[pst_module])
        try:
            Postprocess = _temp.__dict__[pst_module]
        except KeyError:
            sys.stderr.write("Class %s not defined in %s\n"%(pst_module,'Postprocess'))
            sys.exit()
        except:
            sys.stderr.write("Unknown error importing %s\n"%pst_module)
            sys.exit()
except :
    sys.stderr.write("Warning: postprocess module not loaded\n")
    

class Optimizer(object):

    def __init__(self, tm5Obj, resume=False, cleanup=False):
        self.output_dir = tm5Obj.output_dir
        self.tm5 = tm5Obj

        if resume:
            self.load(cleanup=cleanup)
            logging.debug("Loading previously saved Optimizer data")
        else:
            self.rcf = tm5Obj.rcf
            self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
            self.ntracer = len(self.species)
            if tm5Obj.Optim_dict['emission']:
                self.EmisFile = tm5Obj.emission_file_name
            else:
                self.EmisFile = None
            self.ParamFile = os.path.join(self.output_dir, self.rcf.get('parameter.scratch.filename')) if tm5Obj.Optim_dict['parameters'] else None
            self.IniFile = 'NotDoneYet.txt' if tm5Obj.Optim_dict['iniconc'] else None
            self.max_iter = self.rcf.get('optimize.maximum.iterations', 'int', default=300)
            self.iter = 0
            self.grad_norm_reduc_factor = self.tm5.dconf.optim.get('gradient_reduc', 1.e12)
            self.initial_grad_norm_preco = None
            self.converged_eigvals = None
            self.converged_eigvecs = None
            self.postprocess_flag = self.rcf.get('optimize.postprocess', 'bool', default=False)
            self.monitor_flag = self.rcf.get('optimize.monitor.outputs', 'bool', default=False)
            if self.monitor_flag:
                outputs_to_monitor = self.rcf.get('optimize.monitor.output.list')
                outputs_to_monitor = outputs_to_monitor.split()
                # Which ones have adjoint inputs as well as outputs?
                self.outputs_to_monitor = []
                for key in outputs_to_monitor:
                    self.outputs_to_monitor.append('output.' + key)
                    if self.rcf.has_key('adjoint.input.'+key):
                        self.outputs_to_monitor.append('adjoint.input.'+key)
            self.store_steps = self.rcf.get('optimize.store.indiv.steps', 'bool', default=False)
            if self.store_steps:
                self.step_num = 0

    def StoreStep(self):
        # Store the point outputs and departures, satellite outputs and departures, the emission and the adjoint emission files
        # in a sub-folder of self.output_dir. This is for debugging. I want to know why my Radiocarbon inversion with 2010
        # coverage and old convection gives NaNs.
        dir_name = os.path.join(self.output_dir, 'record_iter_%03i'%self.step_num)
        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)

        # Make a list of all files to be copied
        file_list = [os.path.join(self.output_dir, 'emission.nc4'), os.path.join(self.output_dir, 'adj_emissions.nc4')]

        point_outs = glob.glob(os.path.join(self.output_dir, 'point', 'point_output*.nc4'))
        point_deps = glob.glob(os.path.join(self.output_dir, 'point', 'point_departures*.nc4'))
        for file_name in point_outs + point_deps:
            if os.path.exists(file_name):
                file_list.append(file_name)

        sat_outs = glob.glob(os.path.join(self.output_dir, 'satellite', 'sat-track*.nc4'))
        for file_name in sat_outs:
            if os.path.exists(file_name):
                file_list.append(file_name)

        for file_name in file_list:
            shutil.copy(file_name, dir_name)

        self.step_num += 1

    def read_eigsys(self):
        # This dummy routine has to be here to provide compatibility with congrad
        pass

    def moveFiles(self, prefix, postfix):
        # Move {prefix}*.nc4 to {prefix}_{postfix}*.nc4
        if prefix == 'stations':
            split_period = self.rcf.get('output.point.split.period', default='a')
        elif prefix == 'tccon':
            split_period = self.rcf.get('output.satellite.split.period', default='a')
        else:
            sys.stderr.write('Split period not defined for output type %s\n'%prefix)
            raise
        if split_period == 'a':
	    src_files = [os.path.join(self.output_dir, prefix, '%s.nc4'%prefix)]
	    dst_files = [os.path.join(self.output_dir, prefix, '%s_%s.nc4'%(prefix, postfix))]
        elif split_period in ['m', 'd']:
            if split_period == 'd':
                file_patt = os.path.join(self.output_dir, prefix, '%s_????????.nc4'%prefix)
            elif split_period == 'm':
                file_patt = os.path.join(self.output_dir, prefix, '%s_??????.nc4'%prefix)
            src_files = glob.glob(file_patt)
            dst_files = []
            for src_file in src_files:
                dir_name = os.path.dirname(src_file)
                base_name = os.path.basename(src_file)
                pre, post = base_name.split('_')
                base_name = '_'.join([pre, postfix, post])
                dst_file = os.path.join(dir_name, base_name)
                dst_files.append(dst_file)
        for src_file, dst_file in zip(src_files, dst_files):
            try:
                shutil.move(src_file, dst_file)
            except IOError:
                sys.stderr.write("Failed to move %s to %s\n"%(src_file, dst_file))

    def optim_state_file_name(self):
        optim_output_file = os.path.join(self.output_dir, self.rcf.get('optimize.optimizedstate'))
        ## add optim-specific part to output file
        #filename, extension = os.path.splitext(optim_output_file)
        ## self.optim_type should be set from the derived class
        #optim_output_file = filename+'-'+self.optim_type+extension
        return optim_output_file

    def SetupOptimizer(self, debug=False, restart=False, optimized_prior_state=False, emclasses=None):
        # Clean up prior adjoint states
        self.tm5.cleanUpOldAdjoint()
        # Now we need to write the apri emissions
        if optimized_prior_state:
            state_file_name = self.optim_state_file_name()
            self.tm5.SetupEmissions(True, state_file_name)
        else:
            self.tm5.SetupEmissions(emclasses) # This writes the prior flux, i.e., the flux corresponding to x_c=0
        # Create the preconditioner
        self.tm5.makePrecon()
        # write the temporal correlation to a file
        self.tm5.writeTemporalCorrelation()
        # read the necessary file names
        emis_init = self.EmisFile
        iniconc_init = 'NotDoneYet.txt' if self.tm5.Optim_dict['iniconc'] else None
        param_init = self.rcf.get('parameter.apri.filename') if self.tm5.Optim_dict['parameters'] else None
        # Set up the prior x_c to be 0, and of the same shape as the prior emission
        self.tm5.ReadPriorState(iniconc_init,param_init)

        if self.optim_type in ['m1qn3','conGrad']:
            exec_name = self.rcf.get('optimize.' + self.optim_type + '.exec')
            comm_file = os.path.join(self.output_dir, os.path.basename(self.rcf.get('optimize.communication.file')))
            if os.path.exists(comm_file) and not restart:
                os.remove(comm_file)
            self.comm_file = comm_file
            self.OptimCommand = [exec_name, '--state-file', comm_file]
            #if self.optim_type == 'conGrad':
                #self.OptimCommand = [exec_name, '--state-file', comm_file, "--write-traject"]
            #else:
                #self.OptimCommand = [exec_name, '--state-file', comm_file]
        self.Restart = restart
        if restart:
            sys.stderr.write("Restarting optimization, check state_trajectory.nc4\n")
        # If we're restarting, 'touch' all files in the date-specific output folder, to prevent their deletion by the time the optimization converges
        for fileName in glob.glob(os.path.join(self.output_dir, '*')):
            if os.path.isfile(fileName):
                os.utime(fileName, None)

    def EvalFunc_pre(self, state_vector):
        self.iter += 1
        self.first_run = self.true_iter == 0
        # The idea is to create a function that returns a (f, g) tuple where f is the functional value and g is the gradient
        # First, the state vector has to be converted to model space and stored in the input file
        self.tm5.state_preco = state_vector
        self.tm5.StoreToStateFile(self.EmisFile,self.IniFile,self.ParamFile)
        # If we're running forward for the first time and this is not a restart, turn on all outputs to generate apri mismatches
        # If we're running forward for the last time after convergence, turn on all outputs to generate apos mismatches
        # However, we need to store the original values in a dictionary, to restore after the first iteration
        if (self.first_run or self.converged) and self.monitor_flag:
            self.original_outputs = {}
            for key_to_set in self.outputs_to_monitor:
                self.original_outputs[key_to_set] = self.rcf.get(key_to_set)
                self.tm5.rcf.replace(key_to_set, 'T')

    def EvalFunc_fwd(self):
        # Now run the model forward to generate the cost function
        self.tm5.RunForward()

    def EvalFunc_mid(self):
        if self.converged:
            # If this run was after convergence, save stations.nc4 to stations_apos.nc4, and tccon.nc4 to tccon_apos.nc4
            if self.tm5.rcf.get('output.station.timeseries', 'bool'):
                self.moveFiles('stations', 'apos')
            if self.tm5.rcf.get('output.tccon', 'bool'):
                self.moveFiles('tccon', 'apos')
        elif self.first_run:
            # For the first run, rename stations.nc4 to stations_apri.nc4, and tccon.nc4 to tccon_apri.nc4
            if self.tm5.rcf.get('output.station.timeseries', 'bool'):
                self.moveFiles('stations', 'apri')
            if self.tm5.rcf.get('output.tccon', 'bool'):
                self.moveFiles('tccon', 'apri')
        # Store the model representation errors from the first run
        if self.rcf.get('adjoint.input.satellite', 'bool') or self.rcf.get('adjoint.input.point', 'bool'):
            if self.first_run:
                self.tm5.StoreModelReprErrors()
            else:
                self.tm5.RestoreModelReprErrors()
        if self.tm5.rcf.get('adjoint.input.satellite','bool'):
            self.tm5.SatDepartures()
        if self.tm5.rcf.get('adjoint.input.point', 'bool'):
            self.tm5.PointDepartures()
            if self.first_run :
                if os.path.isdir(os.path.join(self.output_dir, 'point_apri')): 
                    shutil.rmtree(os.path.join(self.output_dir, 'point_apri'))
                shutil.copytree(os.path.join(self.output_dir, 'point'), os.path.join(self.output_dir, 'point_apri'))
        if (self.first_run or self.converged) and self.postprocess_flag:
            pst = Postprocess(self.tm5)
            if self.first_run:
                pst.summarizeOutputs('first')
                # store the apri emission
                shutil.copy(self.EmisFile, os.path.join(self.output_dir, 'emission_apri.nc4'))
            elif self.converged:
                pst.summarizeOutputs('last')
                # store the apos emission
                shutil.copy(self.EmisFile, os.path.join(self.output_dir, 'emission_apos.nc4'))
        # If we had earlier set some outputs to true during the first run to generate apri/apos mismatches, restore them to their proper values
        if (self.first_run or self.converged) and self.monitor_flag:
            for key, value in self.original_outputs.items():
                self.tm5.rcf.replace(key, value)

    def EvalFunc_adj(self):
        # Run the model backward to generate the gradient
        self.tm5.RunBackward()

    def EvalFunc_post(self):
        # Read in the cost function
        self.CalculateCosts()
        sys.stderr.write('Cost function calculated\n')
        # Read in the gradient
        self.tm5.CalculateGradient()
        sys.stderr.write('Gradient calculated\n')
        # return the cost function and its gradient
        self.cost_fn = self.J_tot
        self.gradient_preco = self.tm5.adjoint_state_preco
        # Store the step in a folder, if required
        if self.store_steps:
            self.StoreStep()
            sys.stderr.write('Step stored\n')

    def OptimFunc(self, state_vector):
        # OptimFunc should set the values for self.cost_fn and self.gradient_preco
        self.EvalFunc_pre(state_vector) # before running TM5 forward
        self.EvalFunc_fwd() # the TM5 forward run
        self.EvalFunc_mid() # stuff to do between TM5 forward and adjoint runs
        self.EvalFunc_adj() # the TM5 adjoint run
        self.EvalFunc_post() # after running TM5 adjoint

    def writeOptimStateFromEmis(self):
        """
        Sometimes we want to take emission_apri.nc4 or emission_apos.nc4 or just plain emission.nc4 and write it out to
        an optimized state file. This is useful if we want to, say, use VPP to aggregate prior emissions without doing
        an optimization.
        """
        optim_state_file = self.optim_state_file_name()
        if os.path.exists(optim_state_file):
            os.remove(optim_state_file)
        # Check which of emission_apri, emission_apos and emission files exist
        apri_file = os.path.join(self.output_dir, 'emission_apri.nc4')
        if not os.path.exists(apri_file):
            apri_file = os.path.join(self.output_dir, 'emission.nc4')
        apos_file = os.path.join(self.output_dir, 'emission_apos.nc4')
        read_apri = os.path.isfile(apri_file)
        read_apos = os.path.isfile(apos_file)

        if read_apri:
            apri_fid = Dataset(apri_file, 'r')
            ref_fid = apri_fid
        if read_apos:
            apos_fid = Dataset(apos_file, 'r')
            ref_fid = apos_fid

        if not (read_apri or read_apos):
            raise ValueError('At least one emission file must exist')

        ofid = Dataset(optim_state_file, 'w')
        ofid.start_time = self.tm5.StartTime.strftime("%Y-%m-%d %H:%M:%S")
        ofid.end_time = self.tm5.EndTime.strftime("%Y-%m-%d %H:%M:%S")
        ofid.createDimension('itime', 6)
        for region in self.tm5.region_names:
            reg_gid = ofid.createGroup(region)
            reg_gid.createDimension('latitude', len(ref_fid.groups[region].dimensions['latitude']))
            reg_gid.createDimension('longitude', len(ref_fid.groups[region].dimensions['longitude']))
            for tracer in self.species:
                tra_gid = reg_gid.createGroup(tracer)
                # What are the categories for this tracer?
                categories = ref_fid.groups[region].groups[tracer].groups.keys()
                for category in categories:
                    cat_gid = tra_gid.createGroup(category)
                    ref_gid = ref_fid.groups[region].groups[tracer].groups[category]

                    nt = len(ref_gid.dimensions['nt'])
                    cat_gid.createDimension('nt', nt)

                    var = cat_gid.createVariable('time_start', int16, ('nt', 'itime'))
                    var[:] = ref_gid.variables['time_start'][:]

                    var = cat_gid.createVariable('time_end', int16, ('nt', 'itime'))
                    var[:] = ref_gid.variables['time_end'][:]

                    var = cat_gid.createVariable('time_mid', int16, ('nt', 'itime'))
                    var[:] = ref_gid.variables['time_mid'][:]

                    cat_gid.time_resolution = ref_gid.time_resolution
                    cat_gid.optimize = ref_gid.optimize

                    if read_apri:
                        var = cat_gid.createVariable('prior_emission', float64, ('nt', 'latitude', 'longitude'))
                        var[:] = apri_fid.groups[region].groups[tracer].groups[category].variables['emission'][:]
                        var.unit = 'Kg %s/grid cell/second'%(tracer)

                    if read_apos:
                        var = cat_gid.createVariable('poste_emission', float64, ('nt', 'latitude', 'longitude'))
                        var[:] = apos_fid.groups[region].groups[tracer].groups[category].variables['emission'][:]
                        var.unit = 'Kg %s/grid cell/second'%(tracer)

        ofid.close()

        if read_apri:
            apri_fid.close()
        if read_apos:
            apos_fid.close()

    def writeOnlyPrior(self, optim_output_file):
        """
        Sometimes we need to aggregate the prior errors with different correlation lengths, times and magnitudes.
        It's not necessary to run a 4DVAR for that, but we do need to write out the prior emissions and errors
        for the VPP routines to read. Hence this routine writes out an optimized state file a la self.FinishUp,
        but without any information about posterior emissions.
        """
        if os.path.exists(optim_output_file):
            os.remove(optim_output_file)
        # Fill self.tm5 with zoom region info
        self.tm5.determine_children_etc()
        emission_model_prior, iniconc_model_prior, param_model_prior = self.tm5.preco.state2struct(self.tm5.state_prior_model)
        apri_std_emis, apri_std_iniconc, apri_std_param = self.tm5.preco.state2struct(self.tm5.preco.G_state, replace_field=self.tm5.preco.emis_dapri)
        outFid = Dataset(optim_output_file, 'w')
        outFid.start_time = self.tm5.StartTime.strftime("%Y-%m-%d %H:%M:%S")
        outFid.end_time = self.tm5.EndTime.strftime("%Y-%m-%d %H:%M:%S")
        outFid.createDimension('itime', 6)
        for reg_num, region in enumerate(self.tm5.region_names):
            reg_gid = outFid.createGroup(region)
            reg_gid.createDimension('latitude', self.tm5.region_jm[reg_num])
            reg_gid.createDimension('longitude', self.tm5.region_im[reg_num])
            for tracer in self.species:
                tra_gid = reg_gid.createGroup(tracer)
                for category in self.tm5.Emission[region][tracer]['categories']:
                    cat_gid = tra_gid.createGroup(category)
                    cat_data = self.tm5.Emission[region][tracer][category]
                    # We need to write the timing information from cat_data['time_interval']
                    n_time = len(cat_data['time_interval']['time_mid'])
                    cat_gid.createDimension('nt', n_time)
                    var = cat_gid.createVariable('time_start', int16, ('nt', 'itime'))
                    var[:] = array([d.timetuple()[:6] for d in cat_data['time_interval']['time_start']], int16)
                    var = cat_gid.createVariable('time_end', int16, ('nt', 'itime'))
                    var[:] = array([d.timetuple()[:6] for d in cat_data['time_interval']['time_end']], int16)
                    cat_gid.time_resolution = cat_data['time_resolution'].strip()
                    var = cat_gid.createVariable('prior_emission', emission_model_prior[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                    var[:] = emission_model_prior[region][tracer][category]['emission_data']
                    var = cat_gid.createVariable('prior_emission_std', apri_std_emis[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                    var[:] = apri_std_emis[region][tracer][category]['emission_data']
        outFid.close()

    def FinishUp(self, optim_output_file):
        if os.path.exists(optim_output_file):
            os.remove(optim_output_file)
        # Fill self.tm5 with zoom region info
        sys.stdout.write('Optimizer.FinishUp :: determining children... ') ; sys.stdout.flush()
        self.tm5.determine_children_etc()
        sys.stdout.write('done!\n') ; sys.stdout.flush()
        # Let's decide on what to write out in the output file:
        # The optimized state vector in the preco space, optimized emissions in the model space
        # The gradient in the preco space, whether iterations have converged or not
        # Number of iterations, the final value of the cost function
        # Convert the state vector to an emission structure
        sys.stdout.write('Optimizer.FinishUp :: converting final state vector to model space... ') ; sys.stdout.flush()
        # At this point, self.state_vec_model is identical to self.tm5.state_model at the last iteration
        emission_model, iniconc_model, param_model = self.tm5.preco.state2struct(self.state_vec_model)

        # At this point, emission_model is identical to emissions in RunTM5_base/StoreToStateFile
        sys.stdout.write('done!\n') ; sys.stdout.flush()
        sys.stdout.write('Optimizer.FinishUp :: converting initial state vector to model space... ') ; sys.stdout.flush()
        emission_model_prior, iniconc_model_prior, param_model_prior = self.tm5.preco.state2struct(self.tm5.state_prior_model)
        sys.stdout.write('done!\n') ; sys.stdout.flush()

        outFid = Dataset(optim_output_file, 'w')
        outFid.J_bg = self.J_bg
        outFid.J_obs = self.J_obs
        outFid.J_tot = self.J_tot
        outFid.start_time = self.tm5.StartTime.strftime("%Y-%m-%d %H:%M:%S")
        outFid.end_time = self.tm5.EndTime.strftime("%Y-%m-%d %H:%M:%S")
        if self.tm5.Optim_dict['emission']:
            outFid.createDimension('itime', 6)
            for reg_num, region in enumerate(self.tm5.region_names):
                reg_gid = outFid.createGroup(region)
                reg_gid.createDimension('latitude', self.tm5.region_jm[reg_num])
                reg_gid.createDimension('longitude', self.tm5.region_im[reg_num])
                for tracer in self.species:
                    t_gid = reg_gid.createGroup(tracer)
                    for category in self.tm5.Emission[region][tracer]['categories']:

                        cat_gid = t_gid.createGroup(category)
                        cat_data = self.tm5.Emission[region][tracer][category]
                        # We need to write the timing information from cat_data['time_interval']
                        n_time = len(cat_data['time_interval']['time_mid'])
                        cat_gid.createDimension('nt', n_time)

                        var = cat_gid.createVariable('time_start', int16, ('nt', 'itime'))
                        var[:] = array([d.timetuple()[:6] for d in cat_data['time_interval']['time_start']], int16)
                        var = cat_gid.createVariable('time_end', int16, ('nt', 'itime'))
                        var[:] = array([d.timetuple()[:6] for d in cat_data['time_interval']['time_end']], int16)
                        cat_gid.time_resolution = cat_data['time_resolution'].strip()

                        # Now write the prior and posterior emissions. The code below has to be a duplicate of
                        # RunTM5_base:StoreToStateFile, but I can't put the common part in a single subroutine
                        # because the file structure is different. If a category is not being optimized, then we
                        # do not need to perform this gymnastics with transformations, so we can just add that
                        # to the optim_type=1 category.
                        if not self.tm5.Emission[region][tracer][category]['optimize']:
                            optim_type = 1
                        else:
                            optim_type = self.tm5.preco.optim_type[tracer][category]

                        apos_var = cat_gid.createVariable('poste_emission', emission_model[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                        apri_var = cat_gid.createVariable('prior_emission', emission_model_prior[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))

                        if optim_type == 1:
                            apos_var[:] = emission_model[region][tracer][category]['emission_data']
                            apri_var[:] = emission_model_prior[region][tracer][category]['emission_data']

                        elif optim_type == 2: # SB's transformation
                            apos_var[:] = emission_model[region][tracer][category]['emission_data']
                            apri_var[:] = emission_model_prior[region][tracer][category]['emission_data']

                        elif optim_type == 3: # MK's transformation
                            emisx = emission_model[region][tracer][category]['emission_data']
                            positive_part = where(emisx >= 0.0, emisx, 0.0)
                            negative_part = where(emisx < 0.0, emisx, 0.0)
                            emisx = self.tm5.Prior_Emission[region][tracer][category]['emission_data'] * \
                                where(emisx < 0.0, exp(negative_part), 1.0+positive_part)
                            del positive_part, negative_part
                            apos_var[:] = emisx

                            emisx = emission_model_prior[region][tracer][category]['emission_data']
                            positive_part = where(emisx >= 0.0, emisx, 0.0)
                            negative_part = where(emisx < 0.0, emisx, 0.0)
                            emisx = self.tm5.Prior_Emission[region][tracer][category]['emission_data'] * \
                                where(emisx < 0.0, exp(negative_part), 1.0+positive_part)
                            del positive_part, negative_part
                            apri_var[:] = emisx

        if self.tm5.Optim_dict['parameters']:
            outGid = outFid.createGroup('bias_parameter')
            outGid.createDimension('parameters', self.tm5.preco.n_param)
            var = outGid.createVariable('poste_param', param_model.dtype, ('parameters',))
            var[:] = param_model
            var = outGid.createVariable('prior_param', param_model_prior.dtype, ('parameters',))
            var[:] = param_model_prior
        outFid.func_calls = int32(self.true_iter)
        outFid.converged = int32(self.converged)
        outFid.close()
        del emission_model, emission_model_prior, iniconc_model, iniconc_model_prior, param_model, param_model_prior

    def repairOptimizedState(self):
        """
        It sometimes happens that the damn routine FinishUp crashes, resulting in a corrupt optimized_state.nc4. In that case,
        we can at least fix the prior and posterior emissions by copying over the appropriate fields from emission_apri.nc4
        and emission_apos.nc4, although to fix the prior and posterior errors the inversion will have to be restarted.
        """
        emis_apos_file = os.path.join(self.output_dir, 'emission_apos.nc4')
        emis_apri_file = os.path.join(self.output_dir, 'emission_apri.nc4')
        optim_output_file = self.optim_state_file_name()
        fid_apri = Dataset(emis_apri_file,'r')
        fid_apos = Dataset(emis_apos_file,'r')
        fid_optim = Dataset(optim_output_file, 'a')
        for region_name, region_data in fid_apri.groups.items():
            for tracer, tracer_data in region_data.groups.items():
                for cat_name in tracer_data.groups.keys():
                    fid_optim.groups[region_name].groups[tracer].groups[cat_name].variables['prior_emission'][:] = \
                        fid_apri.groups[region_name].groups[tracer].groups[cat_name].variables['emission'][:]
                    fid_optim.groups[region_name].groups[tracer].groups[cat_name].variables['poste_emission'][:] = \
                        fid_apos.groups[region_name].groups[tracer].groups[cat_name].variables['emission'][:]
        fid_optim.close()
        fid_apos.close()
        fid_apri.close()

    def CalculateCosts(self):
        self.J_tot = self.tm5.CalculateCosts()
        self.J_bg = self.tm5.J_bg
        self.J_obs = self.tm5.J_obs
        line = "Iteration %i: J_tot=%.10f; J_obs=%.10f; J_bg=%.10f\n"%(self.iter, self.J_tot, self.J_obs, self.J_bg)
        with open(os.path.join(self.output_dir, 'costFunction.out'), 'a') as fid :
            fid.write(line)
            fid.close()
        print(line)

    def save(self):
        t1 = datetime.now()
        state_dir = os.path.join(self.output_dir, 'state')
        # if the 'state' folder does not exist, create it
        if not os.path.isdir(state_dir):
            os.makedirs(state_dir)

        file_name = os.path.join(state_dir, 'Optimizer.state')
        if os.path.exists(file_name):
            os.remove(file_name)

        with open(file_name, 'wb') as fid:
            # Do not dump the entire self.__dict__ at once, that creates problems loading when it's really large
            for k, v in self.__dict__.items():

                # If v is an instance of the RunTM5 class, invoke its own save() method
                if k == 'tm5':
                    v.save()

                # If v is a numpy array, save it as an npy file
                elif isinstance(v, ndarray):
                    fname = os.path.join(state_dir, 'optimizer_%s.npy'%k)
                    save(fname, v)

                # Otherwise, pickle it and save it
                else:
                    pickle.dump((k, v), fid, pickle.HIGHEST_PROTOCOL)

        t2 = datetime.now()
        print("Optimizer state saved in %s"%str(t2-t1))

    def load(self, cleanup=True):
        t1 = datetime.now()
        state_dir = os.path.join(self.output_dir, 'state')

        file_name = os.path.join(state_dir, 'Optimizer.state')
        if not os.path.exists(file_name):
            msg = '%s not found'%file_name
            logging.error(msg)
            raise RuntimeError(msg)

        # Load the non-RunTM5 and non-ndarray things first
        with open(file_name, 'rb') as fid:
            try:
                while True:
                    k,v = pickle.load(fid)
                    setattr(self, k, v)
            except EOFError:
                pass
        # Delete the state file
        if cleanup:
            os.remove(file_name)

        # Load the npy files
        npy_files = glob.glob(os.path.join(state_dir, 'optimizer_*.npy'))
        for file_name in npy_files:
            base_name = os.path.basename(file_name)
            obj_name = os.path.splitext(base_name)[0].split('_', 1)[1]
            setattr(self, obj_name, load(file_name))
            # Delete the npy file
            if cleanup:
                os.remove(file_name)

        # Load the RunTM5 object
        self.tm5.load(cleanup=cleanup)

        t2 = datetime.now()
        print "Optimizer state loaded in %s"%str(t2-t1)

    def save_convergence(self):
        # Save the convergence status to a file in the output dir
        converged_rc = os.path.join(self.output_dir, 'convergence.state')
        with file(converged_rc, 'a'):
            os.utime(converged_rc, None)
        conv_rc = rc.RcFile(converged_rc)
        conv_rc.replace_add('converged', self.converged)
        conv_rc.replace_add('finished', self.finished)
        conv_rc.WriteFile(converged_rc)

    def load_convergence(self):
        # Load the convergence status from a previously saved file
        converged_rc = os.path.join(self.output_dir, 'convergence.state')
        with file(converged_rc, 'a'):
            os.utime(converged_rc, None)
        conv_rc = rc.RcFile(converged_rc)
        self.converged = conv_rc.get('converged', 'bool', default=False)
        self.finished = conv_rc.get('finished', 'bool', default=False)

class conGrad(Optimizer):

    def __init__(self, tm5Obj, resume=False, cleanup=False):
        self.output_dir = tm5Obj.output_dir
        self.tm5 = tm5Obj

        if resume:
            self.load(cleanup=cleanup)
            logging.debug("Loading previously saved Optimizer data")
        else:
            Optimizer.__init__(self, tm5Obj)
            self.optim_type = 'conGrad'
            self.true_iter = 0
            self.fixed_iterations = self.tm5.dconf.optim.get('n_iter', 1000)# rcf.get('optimize.fixed.iterations', 'int', default=1000)

    def Var4Dsetup_fwd(self):
        # first forward run during Var4Dsetup
        self.state_vec = self.tm5.state_preco
        self.EvalFunc_pre(self.state_vec)
        self.EvalFunc_fwd()

    def Var4Dsetup_adj(self):
        # first adjoint run during Var4Dsetup
        self.EvalFunc_mid()
        self.EvalFunc_adj()
        self.EvalFunc_post()
        fid = Dataset(self.comm_file, 'w')
        fid.createDimension('n_state', len(self.state_vec))
        fid.createDimension('dim_x', 0)
        fid.createDimension('dim_g', 0)
        fid.iter_max = self.max_iter
        fid.iter_convergence = self.fixed_iterations
        fid.preduc = 1./self.grad_norm_reduc_factor
        fid.congrad_finished = 0
        fid.J_tot = self.cost_fn
        var = fid.createVariable('x_c', self.state_vec.dtype, ('n_state', 'dim_x'))
        var[:,0] = self.state_vec
        var = fid.createVariable('g_c', self.state_vec.dtype, ('n_state', 'dim_g'))
        var[:,0] = self.gradient_preco
        fid.close()
        self.initial_grad_norm_preco = linalg.norm(self.gradient_preco)

    def Var4Dsetup(self):
        self.converged = False
        self.finished = False
        if self.Restart:
            fid = Dataset(self.comm_file, 'r')
            congrad_finished = fid.congrad_finished
            fid.close()
            self.restoreFile(congrad_finished)
        else:
            # Should run self.OptimFunc(self.tm5.state_preco), but split it up into fwd and adj runs
            self.Var4Dsetup_fwd()
            self.Var4Dsetup_adj()

    def Var4Ddone(self):
        optim_output_file = self.optim_state_file_name()
        self.converged = int(self.converged)
        if self.converged == 1:
            # convert the optimized state to model space
            self.state_vec_model = self.tm5.preco.xc_to_x(self.state_vec, self.tm5.state_prior_model)
            print 'Calling FinishUp with ', optim_output_file
            self.FinishUp(optim_output_file, True)
            sys.stderr.write("VAR4D :: Optimization finished successfully\n")
        else:
            if self.finished:
                sys.stderr.write("VAR4D :: Optimization failed because max_iter was reached\n")
            else:
                sys.stderr.write("VAR4D :: Optimization failed for unknown reason, please investigate\n")
        sys.stderr.flush()

    def Var4Dstep_check(self):
        # decide where we are in the otimization process, and what to do
        # when the optimizer is called, there should be n gradients and n state vectors in the file
        self.checkStateFile(self.comm_file)
        sbp.check_call(self.OptimCommand)
        # After calling the optimizer, the number of state vectors in the file should have increased by 1,
        # whereas the number of gradients should not have changed. Check that.
        with Dataset(self.comm_file, 'r') as fid:
            len_x = len(fid.dimensions['dim_x'])
            len_g = len(fid.dimensions['dim_g'])

            if len_x != self.state_file_dims['x'] + 1:
                raise RuntimeError('The number of state vectors should have increased from %i to %i, instead it is %i'\
                    %(self.state_file_dims['x'], self.state_file_dims['x']+1, len_x))

            if len_g != self.state_file_dims['g']:
                raise RuntimeError('The number of gradients should have remained unchanged at %i, instead it is %i'\
                    %(self.state_file_dims['g'], len_g))

            self.state_file_dims['x'] = len_x
            self.state_file_dims['g'] = len_g

            # after the call to the optimizer, there are n gradients and n+1 state vectors in the file
            # the (n+1)th state vector is where the optimizer is asking for the gradient now
            self.true_iter = len_x - 1
            self.congrad_finished = fid.congrad_finished
        self.optim_finished = self.congrad_finished

    def read_xc(self):
        with Dataset(self.comm_file, 'r') as fid:
            len_x = len(fid.dimensions['dim_x'])
            self.state_vec = fid.variables['x_c'][:,len_x-1]
            # Has the last run also been done?
            len_g = len(fid.dimensions['dim_g'])
            n_eigen = len(fid.dimensions['n_eigen']) if 'n_eigen' in fid.dimensions else 0
            self.last_run_done = (len_x == len_g) and (len_x == n_eigen + 2)

    def write_gc(self):
        with Dataset(self.comm_file, 'a') as fid:
            len_g = len(fid.dimensions['dim_g'])
            fid.variables['g_c'][:,len_g] = self.gradient_preco
            fid.J_tot = self.cost_fn

    def read_eigsys(self):
        with Dataset(self.comm_file, 'r') as fid:
            self.converged_eigvals = fid.variables['eigenvalues'][:]
            self.converged_eigvecs = fid.variables['eigenvectors'][:]

    def Var4Dstep(self):
        self.Var4Dstep_check()
        if self.congrad_finished == 0: # not yet converged and max_iter not reached
            self.read_xc()
            self.OptimFunc(self.state_vec)
            self.write_gc()
        elif self.congrad_finished == 1: # not converged but max_iter reached
            self.finished = True
        elif self.congrad_finished == 2: # converged, do the last run, then read the eigenvectors and eigenvalues
            self.converged = True
            self.finished = True
            self.read_xc()
            if not self.last_run_done:
                sys.stderr.write("Optimization complete, running forward and adjoint with final state vector\n")
                self.OptimFunc(self.state_vec)
                self.write_gc()
            else:
                print 'Last run need not be done, simply reading the current state vector for finishing up'
            self.read_eigsys()

    def Var4D(self):
        self.Var4Dsetup()
        while (not self.converged) and (not self.finished):
            self.Var4Dstep()
        self.Var4Ddone()

    def FinishUp(self, optim_output_file, use_eigen_file=False):
        if os.path.exists(optim_output_file):
            os.remove(optim_output_file)
        # Fill self.tm5 with zoom region info
        sys.stdout.write('Optimizer.FinishUp :: determining children... ') ; sys.stdout.flush()
        self.tm5.determine_children_etc()
        sys.stdout.write('done!\n') ; sys.stdout.flush()
        # Let's decide on what to write out in the output file:
        # The optimized state vector in the preco space, optimized emissions in the model space
        # The gradient in the preco space, whether iterations have converged or not
        # Number of iterations, the final value of the cost function
        # Convert the state vector to an emission structure
        sys.stdout.write('Optimizer.FinishUp :: converting final state vector to model space... ') ; sys.stdout.flush()
        emission_model, iniconc_model, param_model = self.tm5.preco.state2struct(self.state_vec_model)
        sys.stdout.write('done!\n') ; sys.stdout.flush()
        sys.stdout.write('Optimizer.FinishUp :: converting initial state vector to model space... ') ; sys.stdout.flush()
        emission_model_prior, iniconc_model_prior, param_model_prior = self.tm5.preco.state2struct(self.tm5.state_prior_model)
        sys.stdout.write('done!\n') ; sys.stdout.flush()
        # If using conjugate gradients, estimate the posterior covariance matrix
        if self.converged_eigvals is not None:
            # create an emission structure filled with zeros
            sys.stdout.write('Optimizer.FinishUp :: creating null_emission... ') ; sys.stdout.flush()
            null_emission = dict.fromkeys(self.tm5.region_names)
            for region in self.tm5.region_names:
                null_emission[region] = dict.fromkeys(self.species)
                for tracer in self.species:
                    null_emission[region][tracer] = dict.fromkeys(self.tm5.Emission[region][tracer]['categories'])
                    for category in self.tm5.Emission[region][tracer]['categories']:
                        null_emission[region][tracer][category] = {'emission_data': zeros_like(self.tm5.Emission[region][tracer][category]['emission_data'])}
            LE = zeros_like(self.converged_eigvecs)
            sys.stdout.write('done!\n') ; sys.stdout.flush()
            eigenvec_dict = {}
            param_eigenvec_list = []
            for i in range(len(self.converged_eigvals)):
                # If use_eigen_file is true, for each eigenvector, check first whether the file eigenvector_modelspace_%03i.nc
                # exists or not. If it does, read the eigenvector in model space from that file. Else, calculate the
                # eigenvector in model space and store it in such a file.
                t1 = datetime.now()
                sys.stdout.write('Optimizer.FinishUp :: converting eigenvector %3i to model space... '%(i+1)) ; sys.stdout.flush()
                if use_eigen_file:
                    file_name = os.path.join(self.output_dir, 'eigenvector_modelspace_%03i.nc'%(i+1))
                    if os.path.exists(file_name):
                        # read the eigenvector from the file
                        with Dataset(file_name, 'r') as fid:
                            LE[:,i] = fid.variables['eigenvector'][:]
                    else:
                        LE[:,i] = self.tm5.preco.xc_to_x(self.converged_eigvecs[:,i], zeros_like(self.tm5.state_prior_model))
                        # store the eigenvector to the file
                        with Dataset(file_name, 'w') as fid:
                            fid.createDimension('n_state', LE.shape[0])
                            v = fid.createVariable('eigenvector', LE.dtype, ('n_state',))
                            v[:] = LE[:,i]
                else:
                    LE[:,i] = self.tm5.preco.xc_to_x(self.converged_eigvecs[:,i], zeros_like(self.tm5.state_prior_model))
                eigenvec_dict[i], junk_iniconc, param = self.tm5.preco.state2struct(LE[:,i], replace_field=null_emission)
                param_eigenvec_list.append(param)
                t2 = datetime.now()
                sys.stdout.write('done in %s\n'%str(t2-t1)) ; sys.stdout.flush()
                # eigenvec_dict[i] is now a dictionary with the same structure as emission_model and emission_model_prior

            Mat2 = 1.0/self.converged_eigvals - 1.0
            sys.stdout.write('Optimizer.FinishUp :: calculating apos errors... ') ; sys.stdout.flush()
            apos_std_emis, apos_std_iniconc, apos_std_param = self.tm5.preco.state2struct(nan_to_num(sqrt(self.tm5.preco.G_state**2 + inner(LE**2, Mat2))), replace_field=self.tm5.preco.emis_dapri)
            sys.stdout.write('done!\n') ; sys.stdout.flush()
            sys.stdout.write('Optimizer.FinishUp :: calculating apri errors... ') ; sys.stdout.flush()
            apri_std_emis, apri_std_iniconc, apri_std_param = self.tm5.preco.state2struct(self.tm5.preco.G_state, replace_field=self.tm5.preco.emis_dapri)
            sys.stdout.write('done!\n') ; sys.stdout.flush()
            # apos_std_emis and apri_std_emis now have the same structure as emission_model
            del LE, Mat2, null_emission

            # Delete the temporary eigenvector files
            if use_eigen_file:
                for i in range(len(self.converged_eigvals)):
                    file_name = os.path.join(self.output_dir, 'eigenvector_modelspace_%03i.nc'%(i+1))
                    if os.path.exists(file_name):
                        os.remove(file_name)

        outFid = Dataset(optim_output_file, 'w')
        outFid.J_bg = self.J_bg
        outFid.J_obs = self.J_obs
        outFid.J_tot = self.J_tot
        outFid.start_time = self.tm5.StartTime.strftime("%Y-%m-%d %H:%M:%S")
        outFid.end_time = self.tm5.EndTime.strftime("%Y-%m-%d %H:%M:%S")
        if self.tm5.Optim_dict['emission']:
            outFid.createDimension('itime', 6)
            if self.converged_eigvals is not None:
                outFid.createDimension('eigenvalues', len(self.converged_eigvals))
                var = outFid.createVariable('eigenvalues', self.converged_eigvals.dtype, ('eigenvalues',))
                var[:] = self.converged_eigvals
            for reg_num, region in enumerate(self.tm5.region_names):
                reg_gid = outFid.createGroup(region)
                reg_gid.createDimension('latitude', self.tm5.region_jm[reg_num])
                reg_gid.createDimension('longitude', self.tm5.region_im[reg_num])
                for tracer in self.species:
                    t_gid = reg_gid.createGroup(tracer)
                    for category in self.tm5.Emission[region][tracer]['categories']:
                        cat_gid = t_gid.createGroup(category)
                        cat_data = self.tm5.Emission[region][tracer][category]
                        # We need to write the timing information from cat_data['time_interval']
                        n_time = len(cat_data['time_interval']['time_mid'])
                        cat_gid.createDimension('nt', n_time)
                        var = cat_gid.createVariable('time_start', int16, ('nt', 'itime'))
                        var[:] = array([d.timetuple()[:6] for d in cat_data['time_interval']['time_start']], int16)
                        var = cat_gid.createVariable('time_end', int16, ('nt', 'itime'))
                        var[:] = array([d.timetuple()[:6] for d in cat_data['time_interval']['time_end']], int16)
                        cat_gid.time_resolution = cat_data['time_resolution'].strip()
                        # Now write the prior and posterior emissions, their errors, and the Lanczos eigenvectors
                        var = cat_gid.createVariable('poste_emission', emission_model[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                        var[:] = emission_model[region][tracer][category]['emission_data']
                        var = cat_gid.createVariable('prior_emission', emission_model_prior[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                        var[:] = emission_model_prior[region][tracer][category]['emission_data']
                        if self.converged_eigvals is not None:
                            print "Writing eigenvectors and eigenvalues for region %s tracer %s category %s"%(region,tracer,category)
                            var = cat_gid.createVariable('poste_emission_std', apos_std_emis[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                            var[:] = apos_std_emis[region][tracer][category]['emission_data']
                            var = cat_gid.createVariable('prior_emission_std', apri_std_emis[region][tracer][category]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                            var[:] = apri_std_emis[region][tracer][category]['emission_data']
                            var = cat_gid.createVariable('eigenvectors', eigenvec_dict[0][region][tracer][category]['emission_data'].dtype, ('eigenvalues', 'nt', 'latitude', 'longitude'))
                            for j in range(len(self.converged_eigvals)):
                                var[j] = eigenvec_dict[j][region][tracer][category]['emission_data']
        if self.tm5.Optim_dict['parameters']:
            outGid = outFid.createGroup('bias_parameter')
            outGid.createDimension('parameters', self.tm5.preco.n_param)
            var = outGid.createVariable('poste_param', param_model.dtype, ('parameters',))
            var[:] = param_model
            var = outGid.createVariable('prior_param', param_model_prior.dtype, ('parameters',))
            var[:] = param_model_prior
            if self.converged_eigvals is not None:
                var = outGid.createVariable('poste_param_std', apos_std_param.dtype, ('parameters',))
                var[:] = apos_std_param
                var = outGid.createVariable('prior_param_std', apri_std_param.dtype, ('parameters',))
                var[:] = apri_std_param
            # write out the eigenvectors
            param_eigenvec_list = array(param_eigenvec_list, float64)
            var = outGid.createVariable('eigenvectors', float64, ('eigenvalues','parameters'))
            var[:] = param_eigenvec_list
        outFid.func_calls = int32(self.true_iter)
        outFid.converged = int32(self.converged)
        outFid.close()
        del emission_model, emission_model_prior, iniconc_model, iniconc_model_prior, param_model, param_model_prior
        if self.converged_eigvals is not None:
            del apos_std_emis, apri_std_emis, apos_std_iniconc, apri_std_iniconc, apos_std_param, apri_std_param, eigenvec_dict, param_eigenvec_list

    def checkStateFile(self, fileName):
        # The communication file, for example state_trajectory.nc4, should contain the same number of
        # gradients as state vectors when it is fed to the optimizer. Here we check if that is the
        # case or not. If that is not the case, we shave off some iterations, and issue a warning
        # message to the standard error.
        with Dataset(fileName, 'r') as fid:
            len_x = len(fid.dimensions['dim_x'])
            len_g = len(fid.dimensions['dim_g'])
        self.state_file_dims = {'x': len_x, 'g': len_g}

        if min(len_x, len_g) < 2 and self.Restart:
            sys.stderr.write("The state trajectory file %s seems to have at least one vector with no or very little data!\n"%fileName)
            sys.stderr.write("This is normal for a new run, but not for a restart\n")
            raise RuntimeError

        if len_x != len_g:
            min_len = min(len_x, len_g)
            sys.stderr.write("The state trajectory file %s did not have the same dimension for the state vector and the gradient\n"%fileName)
            sys.stderr.write("There were %i state vectors and %i gradients\n"%(len_x,len_g))
            sys.stderr.write("The file has been truncated to %i state vectors and gradients\n"%min_len)
            self.shaveOffIterations(min_len, min_len)
            self.state_file_dims = {'x': min_len, 'g': min_len}

    def restoreFile(self, mode):
        print "restoreFile called with mode %i"%mode
        # mode=0 for interrupted run, mode=1 for finished but non-converged run, mode=2 for restart with stricter convergence criterion
        if mode == 0:
            # the run was interrupted
            # if the run was interrupted while running TM5, there will be one more state vector than gradient
            # we need to knock out that last state vector
            randString = ''.join(choice(string.ascii_letters + string.digits) for x in range(10))
            randFile = os.path.join(os.path.dirname(self.comm_file), randString+'.nc4')
            inFid = Dataset(self.comm_file, 'r')
            outFid = Dataset(randFile, 'w')
            # copy the attributes
            outFid.congrad_finished = 0
            # Sometimes we want to restart with either a stricter/looser convergence criterion, or a higher/lower
            # maximum number of iterations. So here we check if either is true, and update the state file accordingly.
            # iter_max is the number of iterations to do before giving up
            if self.max_iter != inFid.iter_max:
                outFid.iter_max = self.max_iter
            else:
                outFid.iter_max = inFid.iter_max
            # iter_convergence is the number of iterations at which to declare that we've converged
            if self.fixed_iterations != inFid.iter_convergence:
                outFid.iter_convergence = self.fixed_iterations
            else:
                outFid.iter_convergence = inFid.iter_convergence
            # grad norm reduction will be a bit tougher to match, because it's floating point
            tolerance = 1.0E-10 * 1./self.grad_norm_reduc_factor
            if abs(1./self.grad_norm_reduc_factor - inFid.preduc) > tolerance:
                outFid.preduc = 1./self.grad_norm_reduc_factor
            else:
                outFid.preduc = inFid.preduc
            outFid.J_tot = inFid.J_tot
            # copy the dimensions
            outFid.createDimension('n_state', len(inFid.dimensions['n_state']))
            outFid.createDimension('dim_x', 0)
            outFid.createDimension('dim_g', 0)
            # copy the variables
            dim_g = len(inFid.dimensions['dim_g'])
            v = outFid.createVariable('x_c', inFid.variables['x_c'].dtype, ('n_state','dim_x'))
            v[:] = inFid.variables['x_c'][:,:dim_g]
            v = outFid.createVariable('g_c', inFid.variables['g_c'].dtype, ('n_state','dim_g'))
            v[:] = inFid.variables['g_c'][:,:dim_g]
            outFid.close()
            inFid.close()
            print "Both x_c and g_c reduced to %i vectors"%dim_g
            shutil.move(randFile, self.comm_file)
        elif mode == 1:
            # max_iter was reached
            # we keep all the states and gradients (there are n each), but modify the necessary file attributes
            fid = Dataset(self.comm_file, 'a')
            fid.iter_max = self.max_iter
            fid.preduc = 1./self.grad_norm_reduc_factor
            fid.congrad_finished = int32(0)
            fid.close()
        elif mode == 2:
            # previous run converged
            # so now we either want to run with a stricter convergence
            # or we've stupidly deleted the optimized state file and want to regenerate it
            randString = ''.join(choice(string.ascii_letters + string.digits) for x in range(10))
            randFile = os.path.join(os.path.dirname(self.comm_file), randString+'.nc4')
            inFid = Dataset(self.comm_file, 'r')
            outFid = Dataset(randFile, 'w')
            # copy the attributes
            outFid.preduc = 1./self.grad_norm_reduc_factor
            outFid.iter_max = self.max_iter
            outFid.iter_convergence = self.fixed_iterations
            outFid.congrad_finished = int32(0)
            outFid.J_tot = inFid.J_tot
            # copy the dimensions
            outFid.createDimension('n_state', len(inFid.dimensions['n_state']))
            outFid.createDimension('dim_x', 0)
            outFid.createDimension('dim_g', 0)
            # and now the variables
            stricter_convergence = False # means we just want to run the optimizer through
            if inFid.preduc > outFid.preduc :
                stricter_convergence = True
            if inFid.iter_convergence < outFid.iter_convergence:
                stricter_convergence = True
            # are we looking for stricter convergence?
            if stricter_convergence:
                dim_g = len(inFid.dimensions['dim_g']) - 1 # skip the last state
            # or do we want to run the optimizer through?
            else:
                dim_g = len(inFid.dimensions['dim_g'])
                self.converged = True
                self.finished = True
                self.converged_eigvals = inFid.variables['eigenvalues'][:]
                self.converged_eigvecs = inFid.variables['eigenvectors'][:]
                self.cost_fn = inFid.J_tot
                self.state_vec = inFid.variables['x_c'][:,dim_g-1]
            v = outFid.createVariable('x_c', inFid.variables['x_c'].dtype, ('n_state','dim_x'))
            v[:] = inFid.variables['x_c'][:,:dim_g]
            v = outFid.createVariable('g_c', inFid.variables['g_c'].dtype, ('n_state','dim_g'))
            v[:] = inFid.variables['g_c'][:,:dim_g]
            outFid.close()
            inFid.close()
            shutil.move(randFile, self.comm_file)

    def shaveOffIterations(self, num_iters_x, num_iters_g=None):
        """
        Sometimes a 4DVAR run has gone awry, and we need to shave off a few iterations from the state trajectory.
        In that case, supply the number of iterations to keep for the state vector, and number of iterations to keep for the gradient.
        Remember that to restart, we need one more x_c than g_c.
        """
        if num_iters_g == None:
            num_iters_g = num_iters_x - 1
        comm_file = os.path.join(self.output_dir, os.path.basename(self.rcf.get('optimize.communication.file')))
        randString = ''.join(choice(string.ascii_letters + string.digits) for x in range(10))
        randFile = os.path.join(self.output_dir, randString+'.nc4')
        inFid = Dataset(comm_file,'r')
        outFid = Dataset(randFile, 'w')
        # copy the attributes
        outFid.congrad_finished = int32(0)
        outFid.iter_max = inFid.iter_max
        outFid.iter_convergence = inFid.iter_convergence
        outFid.preduc = inFid.preduc
        outFid.J_tot = inFid.J_tot
        # copy the dimensions
        outFid.createDimension('n_state', len(inFid.dimensions['n_state']))
        outFid.createDimension('dim_x', 0)
        outFid.createDimension('dim_g', 0)
        # copy the variables
        v = outFid.createVariable('x_c', inFid.variables['x_c'].dtype, ('n_state','dim_x'))
        v[:] = inFid.variables['x_c'][:,:num_iters_x]
        v = outFid.createVariable('g_c', inFid.variables['g_c'].dtype, ('n_state','dim_g'))
        v[:] = inFid.variables['g_c'][:,:num_iters_g]
        outFid.close()
        inFid.close()
        shutil.move(randFile, comm_file)

