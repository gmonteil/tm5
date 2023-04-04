#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True


import logging
logging.basicConfig(datefmt='%d %b %Y %H:%M:%S', format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


import re, os, shutil, tempfile, hashlib, glob, subprocess
from dateutil.relativedelta import relativedelta
from dateutil import rrule
from numpy import *
from pyshell.base.helper.Utilities import *
from datetime import datetime, timedelta
import cPickle as pickle
# import dill as pickle
from copy import deepcopy

from pyshell.tmflex import rc
from pyshell.proj.tracer.CO2.Observations import Observations
from pyshell.proj.tracer.CO2.PointObs import PointObs
from pyshell.proj.tracer.CO2.Emissions import Emissions
from pyshell.proj.tracer.CO2.IniConc import InitialConcentration
from pyshell.base.main.Precon import Precon
from pyshell.base.main.Emissions import TM5_emission, my_Dataset
from pyshell.base.main.Satellite import ApplySatelliteObs


class ExecEnvironment(object):
    # A class to make sure that foreground runs come back to the current folder
    def __init__(self, cur_dir, exec_dir):
        self.cur_dir = cur_dir
        self.exec_dir = exec_dir

    def __enter__(self):
        os.chdir(self.exec_dir)

    def __exit__(self, type, value, tb):
        os.chdir(self.cur_dir)


class RunTM5(object):

    def __init__(self, StartTime, EndTime, subdir_tag=''):
        """
        The times are given as tuples (year,month,day,hour) with hour optional. The subdir_tag is meant to create a
        subdir within which the output will go. For example, if normally the output would have gone to

        /scratch/shared/sbasu/var4d/CO2/glb/ml60/tropo25/carbonsat_unbiased/output/2008100100-2008101100/

        then with a subdir_tag of 'no_slopes', they end up in

        /scratch/shared/sbasu/var4d/CO2/glb/ml60/tropo25/carbonsat_unbiased/output/2008100100-2008101100/no_slopes/

        This allows multiple forward/adjoint runs over the same period and using the same input data, but perhaps with
        subtle changes in parameters.
        """
        super(RunTM5, self).__init__()
        self.rcf = rc.RcFile(os.environ['pyshell.rc'])
        # The key 'my.tracer' is a comma-separated list of tracer, whereas the key 'my.tracer.name' is a name for the
        # combination. For example, my.tracer could be 'CO,CO2' and my.tracer.name could be 'COCO2'.
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.runid = self.rcf.get('runid')
        self.StartTime = datetime(*StartTime)
        self.EndTime = datetime(*EndTime)
        self.subdir_tag = subdir_tag
        # SBi : Check whether all the following keys are necessary
        self.rcf.add('timerange.start', self.StartTime.strftime("%Y-%m-%d %H:%M:%S"))
        self.rcf.add('timerange.end', self.EndTime.strftime("%Y-%m-%d %H:%M:%S"))
        # SBf
        self.J_bg = 0.0
        self.J_obs = 0.0
        self.J_tot = self.J_bg + self.J_obs
        self.Optim_dict = {'emission': 'optimize.emission', 'iniconc': 'optimize.initialconcentration', 'parameters': 'optimize.parameter'}
        for k,v in self.Optim_dict.items():
            self.Optim_dict[k] = self.rcf.get(v, 'bool')
        self.months = []
        d = self.StartTime.replace(day=1).replace(hour=0).replace(minute=10)
        while d <= self.EndTime:
            self.months.append(d.strftime('%Y%m'))
            d += relativedelta(months=1)
        # Substitute all date-dependent keys with the right dates
        self.substituteTimes()
        # Add the subdir tag to necessary keys
        self.output_dir = self.putTag(self.rcf.get('output.dir'), False)
        self.rcf.replace('output.dir', self.output_dir)
        self.region_names = self.rcf.get('regions').split()
        self.nregions = len(self.region_names)
        # get the emission file name and add the subdir tag
        self.emission_file_name = self.putTag(self.rcf.get('PyShell.em.filename'), True)
        self.rcf.replace('PyShell.em.filename', self.emission_file_name)
        # delete all the tm5-pyshell.*.rc files in the output directory
        files_to_delete = glob.glob(os.path.join(self.rcf.get('my.run.dir'), self.rcf.get('my.basename')+'.*.rc'))
        for fileName in files_to_delete:
            try:
                os.remove(fileName)
            except:
                print 'File ', fileName, ' has been already deleted!'
        self.GetZoomRegions()

    def splitRunTime(self, interval='Y', count=1):
        """
        Split the total model run period into monthly ('M') or yearly ('Y') intervals. Also create the 'save'
        file names for each interval.
        """
        if interval in ['M', 'm']:
            period = rrule.MONTHLY
        elif interval in ['Y', 'y']:
            period = rrule.YEARLY
        elif interval in ['D', 'd']:
            period = rrule.DAILY
        elif interval in ['W', 'w']:
            period = rrule.WEEKLY
        else:
            sys.stderr.write("Invalid interval %s specified for RunTM5.splitRunTime\n"%interval)
            raise

        boundary_times = list(rrule.rrule(period, dtstart=self.StartTime, until=self.EndTime, interval=count))
        # if the ending time is not the last element in boundary times, include it
        if boundary_times[-1] != self.EndTime:
            boundary_times.append(self.EndTime)

        # now create pairs
        ret_dict = {'time intervals': [], 'save files': [], 'adj save files': []}
        n_time = len(boundary_times)
        for i in range(n_time-1):
            start_time = boundary_times[i]
            end_time = boundary_times[i+1]
            save_file = os.path.join(self.output_dir, 'save', end_time.strftime("save_%Y%m%d%H%M.nc4"))
            adj_save_file = os.path.join(self.output_dir, 'save', start_time.strftime("adj_save_%Y%m%d%H%M.nc4"))
            ret_dict['time intervals'].append((start_time, end_time))
            ret_dict['save files'].append(save_file)
            ret_dict['adj save files'].append(adj_save_file)
        ret_dict['n_time'] = n_time-1

        self.split_run_periods = ret_dict

    def substituteTimes(self):
        """
        The rc file is allowed to have character combinations such as <Y1> and <m2>, which mean,
        respectively, the starting year and the ending month. This routine goes through all the
        rc keys and makes the necessary substitutions.
        """
        match_string = []
        for atom in ['Y', 'y', 'm', 'M', 'd', 'H', 'b', 'B', 'j', 'p', 'S']:
            for i in range(1,3):
                rep_string = '<%1s%1i>'%(atom, i)
                match_string.append(rep_string)
        match_string = '|'.join(match_string)
        match_pattern = re.compile(match_string)
        for key in self.rcf.keys():
            value = self.rcf.get(key)
            if match_pattern.search(value):
                value = self.putDateString(value)
                self.rcf.replace(key, value)

    def putDateString(self, input_string):
        output_string = input_string
        for atom in ['Y', 'y', 'm', 'M', 'd', 'H', 'b', 'B', 'j', 'p', 'S']:
            for i, t in zip([1,2], [self.StartTime, self.EndTime]):
                rep_string = '<%1s%1i>'%(atom, i)
                form_string = '%%%1s'%atom
                output_string = output_string.replace(rep_string, t.strftime(form_string))
        return output_string

    def putTag(self, input_path, is_file=True):
        # add self.subdir_tag to input_path, which could be a file or a folder
        if is_file:
            dirname, filename = os.path.split(input_path)
            output_string = os.path.join(dirname, self.subdir_tag, filename)
        else:
            output_string = os.path.join(input_path, self.subdir_tag)
        return output_string

    def cleanUpOldAdjoint(self):
        """
        Before starting a 4DVAR run or a gradient test, it's a good idea to clean up the old adjoints. This
        includes the old adjoint emissions, as well as the old adjoint bias parameters. The latter is particularly
        important. If I'm doing a point-only inversion inside a project that has the potential to optimize
        satellite bias parameters, then an old adjoint bias parameter file can be read in by self.ReadAdjointState
        despite it not being relevant, and could result in an incorrect gradient.
        """
        if self.Optim_dict['emission']:
            adj_emis_file = os.path.join(self.output_dir, 'adj_emissions.nc4')
            if os.path.exists(adj_emis_file):
                os.remove(adj_emis_file)
        if self.Optim_dict['parameters']:
            adj_param_file = os.path.join(self.output_dir, self.rcf.get('adjoint.parameter.filename'))
            if os.path.exists(adj_param_file):
                os.remove(adj_param_file)

    def SetupEmissions(self, randomize=False, zero=False):
        """
        This will either assemble the emissions (default), or read them from a file, depending on
        an rc key, emission.read.optimized.
        """
        t1 = datetime.now()

        emisFile = self.emission_file_name
        if os.path.exists(emisFile):
            os.remove(emisFile)
        checkDir(emisFile)

        # First, process the emissions and write the emissions file
        StartTuple = self.StartTime.timetuple()[:5]
        EndTuple = self.EndTime.timetuple()[:5]
        em = Emissions(StartTuple, EndTuple, self.subdir_tag)
        em.create_emission_structure()

        read_optim_emis = self.rcf.get('emission.read.optimized', 'bool', default=False)
        if read_optim_emis:
            perturb_optim_emis = self.rcf.get('emission.read.optimized.perturb', 'bool', default=False)
            optim_emis_file = self.rcf.get('emission.read.optimized.filename')
            em.readOptimFromFile(optim_emis_file, perturb_optim_emis)
        else:
            for tracer in em.species:
                emfill=em.get_class_from_name(em.Emission[tracer]['emi_class'])(StartTuple, EndTuple, self.subdir_tag)
                emfill.Emission=em.Emission
                emfill.LoopThroughPeriods()

        self.Emission = em(randomize, zero)   # Write emissions and keep them in memory

        for region in self.region_names:
            for tracer in self.species:
                categories = self.Emission[region][tracer]['categories']
                for cat in categories:
                    print region, tracer, cat, self.Emission[region][tracer][cat]['time_resolution'], ' Optimize?:', \
                    self.Emission[region][tracer][cat]['optimize'],self.Emission[region][tracer][cat]['emission_data'].shape
                for cat_key in self.Emission[tracer]['cat_list']:
                    re = self.Emission[tracer]['cat_opt'][cat_key]
                    print '========================================'
                    print 'Tracer: ',tracer,'Category marked for optimization:', cat_key
                    for rei in re:
                        print '----->region:', rei['region'], 'error:', rei['error']
                    print '========================================'
        print '========================================'
        print 'end info parsed from emission routine'
        print '========================================'

        t2 = datetime.now()
        sys.stderr.write("Emission between %s and %s assembled in %s\n"%(self.StartTime.strftime("%c"), self.EndTime.strftime("%c"), str(t2-t1)))

    def makePrecon(self):
        # For coding ease, it's better to create the Precon object here instead of in the Optimizer routines
        self.preco = Precon(self)

    def SetupPointObservations(self):

        point_split_period = self.rcf.get('output.point.split.period', default='a')
        if point_split_period == 'a':
            startTime = self.StartTime
            d_time = self.EndTime - self.StartTime
        elif point_split_period == 'm':
            startTime = datetime(self.StartTime.year, self.StartTime.month, 1, 0, 0, 0)
            d_time = relativedelta(months=1)
        elif point_split_period == 'd':
            startTime = datetime(self.StartTime.year, self.StartTime.month, self.StartTime.day, 0, 0, 0)
            d_time = timedelta(days=1)
        else:
            raise ValueError("Invalid split period %s specified"%point_split_period)

        curTime = startTime
        while curTime < self.EndTime:
            StartTuple = curTime.timetuple()[:6]
            EndTuple = (curTime+d_time).timetuple()[:6]

            obs = Observations(StartTuple, EndTuple)
            obs.create_point_observation_structure()
            for tracer in obs.species:
                obsfill=obs.get_class_from_name(obs.PointObservation[tracer]['obs_class'])(StartTuple, EndTuple)
                obsfill.PointObservation=obs.PointObservation
                obsfill()
            obs.writePointFile()
            curTime = curTime + d_time

    def SetupSatObservations(self):
        sat_split_period = self.rcf.get('output.satellite.split.period')
        if sat_split_period == 'm':
            startTime = datetime(self.StartTime.year, self.StartTime.month, 1, 0, 0, 0)
            d_time = relativedelta(months=1)
        elif sat_split_period == 'd':
            startTime = datetime(self.StartTime.year, self.StartTime.month, self.StartTime.day, 0, 0, 0)
            d_time = timedelta(days=1)
        else:
            raise ValueError("Invalid split period %s specified"%sat_split_period)

        curTime = startTime
        while curTime < self.EndTime:
            StartTuple = curTime.timetuple()[:6]
            EndTuple = (curTime + d_time).timetuple()[:6]
            obs = Observations(StartTuple, EndTuple)
            obs.create_sat_observation_structure()
            for tracer in obs.species:
                obsfill=obs.get_class_from_name(obs.SatObservation[tracer]['obs_class'])(StartTuple, EndTuple)
                obsfill.SatObservation=obs.SatObservation
                obsfill()
            obs.writeSatFile()
            curTime = curTime + d_time

    def SetupIniConc(self, IniDate):
        """
        WARNING: The IniConc classes have not yet been modified to deal with this sytax. Fix that.
        I am not really using this right now. That may change in the future.
        """
        inic = InitialConcentration(IniDate)
        inic.create_iniconc_structure()

        for tracer in inic.species:
            inicfill = inic.get_class_from_name(inic.IniConc[tracer]['iniconc_class'])(IniDate)
            inicfill.IniConc = inic.IniConc
            inicfill.AssembleIniconc()
        inic.WriteIniConc()

    def PointDepartures(self):

        dp = PointObs(self)

        for month_tuple in dp.months:
            dp.create_pointdeparture_structure()
            for tracer in dp.species:
                dpfill = dp.get_class_from_name(dp.PointDepartures[tracer]['departure_class'])(self)
                dpfill.PointDepartures = dp.PointDepartures
                dpfill.applyPointObs(month_tuple)

            dp.writePointDepartureFile(month_tuple)

    def SatDepartures(self):

        ds = ApplySatelliteObs(self)
        ds.create_satdeparture_structure()

        """
        SPECIAL NOTE FOR OPTIMIZING BIAS PARAMETERS

        At this point, ds.SatDepartures[tracer]['optimize bias'] is a boolean denoting whether a specific tracer
        includes bias corrections. If that flag is True, then ds.SatDepartures[tracer]['bias'] is a dictionary with
        'n_param' being the number of bias parameters, and 'param', 'adj_param' and 'param_err' being all 1D float64
        arrays of length 'n_param'. Those three arrays are right now filled with zeros. However, ds.bias_parameter,
        a 1D array, now contains the current/latest values of the bias parameters for all tracers, ds.bias_error
        contains the errors in those parameters, and ds.AdjBiasParameters is a vector of zeros of the same length as
        ds.bias_parameter, which is the total number of bias parameters across all tracers. The tracer-specific
        departure-creating routines should not touch any ds.* variables, but should only modify them through
        ds.SatDepartures[tracer]. Therefore, before calling any tracer-specific routine, ds.bias_parameter, ds.bias_error
        and ds.AdjBiasParameters should be disaggregated into values for respective tracers.

        """

        if ds.optim_bias_param:
            i = 0
            for tracer in ds.species:
                if ds.SatDepartures[tracer]['optimize bias']:
                    n_param = ds.SatDepartures[tracer]['bias']['n_param']
                    # ds.bias_parameter and ds.bias_error contain the current bias parameters and their apri
                    # errors respectively, so at this point disaggregate those into tracer-specific bias
                    # parameters and errors that can be accessed by tracer-specific routines to apply the
                    # averaging kernel
                    ds.SatDepartures[tracer]['bias']['param'] = ds.bias_parameter[i:i+n_param]
                    ds.SatDepartures[tracer]['bias']['param_err'] = ds.bias_error[i:i+n_param]
                    i += n_param

        ds.J_obs_satellite = 0.0
        for month_tuple in ds.months:
            for tracer in ds.species:
                dsfill = ds.get_class_from_name(ds.SatDepartures[tracer]['departure_class'])(self)
                dsfill.SatDepartures = ds.SatDepartures
                # dsfill.SatDepartures now contains information about all previous tracers
                dsfill.applySatelliteObs(month_tuple)
                # ds.SatDepartures[tracer]['bias']['adj_param'] has now been augmented with information
                # from this particular period
            ds.writeDeparturesFile(month_tuple)

        if ds.optim_bias_param:
            # We need to aggregate bias information from all tracers
            i = 0
            for tracer in ds.species:
                if ds.SatDepartures[tracer]['optimize bias']:
                    n_param = ds.SatDepartures[tracer]['bias']['n_param']
                    # the bias parameters and their errors should not have been modified by applySatelliteObs
                    # only adjoint bias parameters should have been updated
                    ds.AdjBiasParameters[i:i+n_param] = ds.SatDepartures[tracer]['bias']['adj_param']
                    i += n_param
            ds.WriteAdjointBiasParameters()

    def writeTemporalCorrelation(self):
        """
        Write the temporal correlation matrices to the emission file. This is only called from within Optimizer.SetupOptimizer,
        since this is only needed in case we want to do a 4DVAR run. Precon() must be called before calling this. After Precon()
        is called, self.Emission['cat_list'] contains n_cat elements, each of which is a category to be optimized. Each element
        looks like:

        [('vary', ' 1000.0-g ', ' 1.00-e-daily+2 ')]

        Also, self.Emission['cat_Temp_L'], self.Emission['cat_Temp_Lt'] and self.Emission['cat_nt'] are all lists. The first two
        contain temporal correlation matrices for all the categories to be optimized, and the third contains the number of time
        steps for each of the optimizable categories. We need to write this info into a file within self.output_dir.

        After a call to SetupOptimizer, self.Emission[tracer] has the following keys:
        ['cat_list', 'cat_Hor_L', 'emi_class', 'cat_Bh_file', 'cat_opt', 'cat_nt', 'cat_n_hor',
        'cat_Hor_Lt', 'tf_bb_diurnal', 'cat_vec2ll_file', 'cat_Temp_L', 'cat_Temp_Lt']
        """
        corrFile = os.path.join(self.output_dir, 'temporal_correlation.nc4')
        if os.path.exists(corrFile):
            os.remove(corrFile)
        checkDir(corrFile)
        fid = my_Dataset(corrFile, 'w')
        for tracer in self.species:
            tid = fid.createGroup(tracer)
            for category_name, category_dim, category_matrix in zip(self.Emission[tracer]['cat_list'], self.Emission[tracer]['cat_nt'], self.Emission[tracer]['cat_Temp_L']):
                gid = tid.createGroup(category_name[0])
                gid.createDimension('nt', category_dim)
                var = gid.createVariable('bt', float64, ('nt','nt'))
                var[:] = dot(category_matrix, category_matrix.T)
        fid.close()

    def Compile(self, new_build=False):
        """
        Compile TM5 via pycasso scripting: build, make
        """
        # rcfile for TM5 following pycasso requirements:
        if not os.path.isdir(self.rcf.get('my.run.dir')):
            os.makedirs(self.rcf.get('my.run.dir'))
        rcfile = os.path.join(self.rcf.get('my.run.dir'),'compile.rc')
        self.rcf.WriteFile(rcfile)
        # command to setup and submit a run:
        if new_build:
            command = [os.path.join(os.curdir,'setup_tm5'), '--new', rcfile]
        else:
            command = [os.path.join(os.curdir,'setup_tm5'), rcfile]
        # run, check status:
        subprocess.check_call( command )

    def RunForward(self):
        """
        Run TM5 via pycasso scripting: build, make, submit.
        """
        ## Save and delete the precon object, because it's huge
        #if 'preco' in self.__dict__.keys():
            #self.preco.save()
            #del self.preco
            #had_precon = True
        #else:
            #had_precon = False

        t1 = datetime.now()
        # rcfile for TM5 following pycasso requirements:
        checkDir(self.output_dir, True)
        curdir = os.getcwd()
        exec_dir = self.rcf.get('my.run.dir')
        self.rcf.replace_add('my.runmode','1')
        rcfile = os.path.join(self.output_dir, 'forward.rc')
        self.rcf.WriteFile(rcfile)
        # Delete the old tm5.ok file
        ok_file = os.path.join(self.output_dir, 'tm5.ok')
        if os.path.exists(ok_file):
            os.remove(ok_file)
        print('Nthreads = ', os.getenv('OMP_NUM_THREADS'))
        # command to setup and submit a run:
        # command = ['submit_tm5', self.rcf.get('job.step.run.exe'), rcfile]
        command = [os.path.join(self.rcf.get('pyshell2.build_directory'), 'tm5.x'), rcfile]
        # run, check status:
        with ExecEnvironment(curdir, exec_dir):
            logger.info(command)
            subprocess.check_call(command)
        t2 = datetime.now()
        # Check if run completed successfully
        if not os.path.exists(ok_file):
            sys.stderr.write("Forward run did not complete successfully\n")
            sys.exit(41)
        sys.stderr.write("Forward run from %s to %s took %s on %s threads\n"%(self.StartTime.strftime("%c"), self.EndTime.strftime("%c"), str(t2-t1), os.getenv('OMP_NUM_THREADS')))

        ## Re-load the precon object
        #if had_precon:
            ## Load the Precon object
            #self.preco = Precon(self, resume=True, cleanup=True)

    def RunForward_step(self, step_num):
        """
        Same as RunForward, except that the total run period is split into several smaller ones, and
        only one of them is run. The split period are taken from self.split_run_periods
        """
        if step_num >= self.split_run_periods['n_time']:
            sys.stderr.write("You are asking me to run period %i, but there are only %i periods\n"%\
                (step_num+1, self.split_run_periods['n_time']))
            raise
        start_time, end_time = self.split_run_periods['time intervals'][step_num]
        self.rcf.replace('timerange.start', start_time.strftime("%Y-%m-%d %H:%M:%S"))
        self.rcf.replace('timerange.end', end_time.strftime("%Y-%m-%d %H:%M:%S"))
        self.rcf.replace_add('restart.write', 'T')
        if step_num > 0:
            self.rcf.replace('istart', 33)

        t1 = datetime.now()
        # rcfile for TM5 following pycasso requirements:
        checkDir(self.output_dir, True)
        curdir = os.getcwd()
        exec_dir = self.rcf.get('my.run.dir')
        self.rcf.replace_add('my.runmode','1')
        rcfile = os.path.join(self.output_dir, 'forward.rc')
        self.rcf.WriteFile(rcfile)
        # Delete the old tm5.ok file
        ok_file = os.path.join(self.output_dir, 'tm5.ok')
        if os.path.exists(ok_file):
            os.remove(ok_file)
        # command to setup and submit a run:
        #command = [os.path.join(os.curdir,'submit_tm5'), self.rcf.get('job.step.run.exe'), rcfile]
        command = [os.path.join(self.rcf.get('pyshell2.build_directory'), 'tm5.x'), rcfile]
        # run, check status:
        print('Nthreads = ', os.getenv('OMP_NUM_THREADS'))
        with ExecEnvironment(curdir, exec_dir):
            logger.info(command)
            subprocess.check_call(command)
        t2 = datetime.now()
        # Check if run completed successfully
        if not os.path.exists(ok_file):
            sys.stderr.write("Forward run did not complete successfully\n")
            sys.exit(41)
        sys.stderr.write("Forward run from %s to %s took %s on %s threads\n"%(start_time.strftime("%c"), end_time.strftime("%c"), str(t2-t1), os.getenv('OMP_NUM_THREADS')))

    def RunBackward(self):
        """
        Run adjoint TM5 via pycasso scripting: build, make, submit.
        """
        ## Save and delete the precon object, because it's huge
        #if 'preco' in self.__dict__.keys():
            #self.preco.save()
            #del self.preco
            #had_precon = True
        #else:
            #had_precon = False

        t1 = datetime.now()
        checkDir(self.output_dir, True)
        curdir = os.getcwd()
        exec_dir = self.rcf.get('my.run.dir')
        self.rcf.replace_add('my.runmode','2')
        rcfile = os.path.join(self.output_dir, 'backward.rc')
        self.rcf.WriteFile(rcfile)
        # Delete the old tm5.ok file
        ok_file = os.path.join(self.output_dir, 'tm5.ok')
        if os.path.exists(ok_file):
            os.remove(ok_file)
        # command to setup and submit a run:
        # command = ['submit_tm5', self.rcf.get('job.step.run.exe'), rcfile]
        command = [os.path.join(self.rcf.get('pyshell2.build_directory'), 'tm5.x'), rcfile]
        # run, check status:
        with ExecEnvironment(curdir, exec_dir):
            print('Nthreads = ', os.getenv('OMP_NUM_THREADS'))
            logger.info(command)
            subprocess.check_call(command)
        t2 = datetime.now()
        # Check if run completed successfully
        if not os.path.exists(ok_file):
            sys.stderr.write("Adjoint run did not complete successfully\n")
            sys.exit(42)
        sys.stderr.write("Adjoint run from %s to %s took %s on %s threads\n"%(self.StartTime.strftime("%c"), self.EndTime.strftime("%c"), str(t2-t1), os.getenv('OMP_NUM_THREADS')))

        ## Re-load the precon object
        #if had_precon:
            ## Load the Precon object
            #self.preco = Precon(self, resume=True, cleanup=True)

    def CalculateCosts(self):
        self.J_bg = self.backgroundCost()
        sys.stderr.write('Background cost calculated\n')
        self.J_obs = self.ObservationCost()
        sys.stderr.write('Observation cost calculated\n')
        self.J_tot = self.J_bg + self.J_obs
        return self.J_tot

    def CalculateGradient(self, add_bg=True):
        adjoint_emission_model, adjoint_iniconc_model, adjoint_parameter_model = self.ReadAdjointState()
        sys.stderr.write('Adjoint state read\n')
        self.adjoint_state_model = self.preco.struct2state(adjoint_emission_model,adjoint_iniconc_model,adjoint_parameter_model)
        sys.stderr.write('Adjoint state calculated in model space\n')
        self.adj_state_norm = linalg.norm(self.adjoint_state_model) # Just the observations
        self.adjoint_state_preco = self.preco.g_to_gc(self.adjoint_state_model) # Convert g to preconditioned space (g_c)
        sys.stderr.write('Adjoint state calculated in preco space\n')
        # This gradient is just the gradient of J_obs. To get the total gradient we have to add x_c
        # However, for the adjoint test we might not want to do that, so add a flag
        if add_bg:
            self.adjoint_state_preco = self.adjoint_state_preco + self.backgroundGradient()
            sys.stderr.write('Background gradient added\n')
        self.adj_state_norm_preco = linalg.norm(self.adjoint_state_preco)

    def backgroundGradient(self):
        # The portion of the gradient of the cost function from the background term, in preco space
        return self.state_preco - self.state_prior_preco

    def backgroundCost(self):
        # The cost function from the background term
        state_vec_diff = self.state_preco - self.state_prior_preco
        return 0.5*sum(state_vec_diff*state_vec_diff)

    def ReadAdjointState(self):
        AdjEmission = None
        AdjParameters = None
        AdjIniconc = None
        if self.Optim_dict['emission']:
            AdjEmisFile = os.path.join(self.rcf.get('output.dir'), 'adj_emissions.nc4')
            fid = my_Dataset(AdjEmisFile, 'r')
            # put info in the default structure:
            AdjEmission= dict.fromkeys(self.region_names)
            for region in (self.region_names):
                rgroup = fid.groups[region]
                AdjEmission[region] = dict.fromkeys(self.species)
                for tracer in self.species:
                    tgroup = rgroup.groups[tracer]
                    categories = self.Emission[region][tracer]['categories']
                    AdjEmission[region][tracer] = dict.fromkeys(categories)
                    for cat in categories:
                        cgroup = tgroup.groups[cat]

                        if self.Emission[region][tracer][cat]['optimize']:

                            if self.preco.optim_type[tracer][cat] == 1:
                                AdjEmission[region][tracer][cat] = {'emission_data': cgroup.variables['adj_emis'][:]}

                            elif self.preco.optim_type[tracer][cat] == 2:
                                AdjEmission[region][tracer][cat] = {'emission_data': cgroup.variables['adj_emis'][:]}

                            elif self.preco.optim_type[tracer][cat] == 3:
                                emisx = self.Emission[region][tracer][cat]['emission_data']
                                emisp = self.Prior_Emission[region][tracer][cat]['emission_data']
                                Adex = cgroup.variables['adj_emis'][:]
                                emisx = where(emisx < emisp, Adex*emisx, Adex*emisp)
                                AdjEmission[region][tracer][cat] = {'emission_data': emisx}
                                del emisx, emisp, Adex

                        else: # category not optimized
                            AdjEmission[region][tracer][cat] = {'emission_data': cgroup.variables['adj_emis'][:]}

            fid.close()

        if self.Optim_dict['iniconc']:
            print 'Not implemented yet, we should open the file BACKWARD_SAVE for this run... ReadAdjointState, class RunTM5...'
            sys.exit(0)

        if self.Optim_dict['parameters']:
            AdjParamFile = os.path.join(self.output_dir, self.rcf.get('adjoint.parameter.filename'))
            if os.path.isfile(AdjParamFile):
                AdjParameters = loadtxt(AdjParamFile, dtype=float64)
            else:
                sys.stderr.write("WARNING :: I am supposed to optimize bias parameters, but the adjoint bias parameter file\n")
                sys.stderr.write("WARNING :: %s\n"%AdjParamFile)
                sys.stderr.write("WARNING :: does not exist or is not a regular file\n")
                # We still need to supply some dummy adjoint bias parameters, for which we need the number of parameters
                dapri_file_name = self.rcf.get('parameter.prior.error.filename')
                parameters_dapri = loadtxt(dapri_file_name, dtype=float64)
                n_param = len(parameters_dapri)
                AdjParameters = zeros(n_param, float64)

        return AdjEmission, AdjIniconc, AdjParameters

    #def ReadAdjointState(self):
        #AdjEmission = None
        #AdjParameters = None
        #AdjIniconc = None
        #if self.Optim_dict['emission']:
            #AdjEmisFile = os.path.join(self.rcf.get('output.dir'), 'adj_emissions.nc4')
            #fid = my_Dataset(AdjEmisFile, 'r')
            ## put info in the default structure:
            #AdjEmission= dict.fromkeys(self.region_names)
            #for region in (self.region_names):
                #rgroup = fid.groups[region]
                #AdjEmission[region] = dict.fromkeys(self.species)
                #for tracer in self.species:
                    #tgroup = rgroup.groups[tracer]
                    #categories = self.Emission[region][tracer]['categories']
                    #AdjEmission[region][tracer] = dict.fromkeys(categories)
                    #for cat in categories:
                        #cgroup = tgroup.groups[cat]
                        #AdjEmission[region][tracer][cat] = {'emission_data': cgroup.variables['adj_emis'][:]}

            #fid.close()

        #if self.Optim_dict['iniconc']:
            #print 'Not implemented yet, we should open the file BACKWARD_SAVE for this run... ReadAdjointState, class RunTM5...'
            #sys.exit(0)

        #if self.Optim_dict['parameters']:
            #AdjParamFile = os.path.join(self.output_dir, self.rcf.get('adjoint.parameter.filename'))
            #AdjParameters = loadtxt(AdjParamFile, dtype=float64)

        #return AdjEmission, AdjIniconc, AdjParameters

    def ObservationCost(self):
        J_obs = 0.0
        if self.rcf.get('adjoint.input.satellite', 'bool'):
            satobs = ApplySatelliteObs(self)
            sys.stderr.write('Satellite departures calculated\n')
            J_obs += satobs.SatelliteCost()
            sys.stderr.write('Satellite cost function calculated\n')
        if self.rcf.get('adjoint.input.point', 'bool'):
            pointobs = PointObs(self)
            sys.stderr.write('Point departures calculated\n')
            J_obs += pointobs.PointCost()
            sys.stderr.write('Point cost function calculated\n')
        return J_obs

    def CreateMeasurementPerturbations(self):
        """
        In order to use Chevallier's MC method for estimating posterior errors, we need to perturb the measurements in accordance
        with their error statistics. This is done by (a) creating a vector of normally distributed elements during the first loop
        of 4DVAR, (b) multiplying that with the total error and adding that to the measurements when creating departures. This
        routine creates the standard normal perturbations -- one set per track file -- and stores them in a single netcdf file.

        Update 11-04-2015: It is useful to not perturb the measurements in some cases, to figure out the spread in the posterior
        just due to the spread in the prior. We implement that by perturbing the measurements by a vector of zeros.
        """
        random_error_file = os.path.join(self.output_dir, 'measurement_perturbations.nc4')
        if os.path.exists(random_error_file):
            os.remove(random_error_file)
        error_scale = self.rcf.get('mc.4dvar.measurement.scale.factor', 'float', default=1.0)

        fid = my_Dataset(random_error_file, 'w')

        point_list = glob.glob(os.path.join(self.output_dir, 'point', 'point_output*.nc4'))
        if self.rcf.get('adjoint.input.point', 'bool') and len(point_list) > 0:
            grp_id = fid.createGroup('point_perturbations')
            for point_track_file in point_list:
                track_fid = my_Dataset(point_track_file, 'r')
                grp_month = grp_id.createGroup(os.path.basename(point_track_file))
                for region, region_data in track_fid.groups.items():
                    grp_region = grp_month.createGroup(region)
                    for tracer, tracer_data in region_data.groups.items():
                        if len(tracer_data.dimensions['samples']) > 0:
                            grp_tracer = grp_region.createGroup(tracer)
                            grp_tracer.createDimension('samples', len(tracer_data.dimensions['samples']))
                            var = grp_tracer.createVariable('random_vector', 'd', ('samples'))
                            var[:] = error_scale * random.standard_normal(len(tracer_data.dimensions['samples']))
                track_fid.close()

        sat_list = glob.glob(os.path.join(self.output_dir, 'satellite', 'sat-track_[0-9]*.nc4'))
        if self.rcf.get('adjoint.input.satellite', 'bool') and len(sat_list) > 0:
            grp_id = fid.createGroup('satellite_perturbations')
            for dep_file in sat_list:
                fid_i = my_Dataset(dep_file, 'r')
                grp_month = grp_id.createGroup(os.path.basename(dep_file))
                for region, region_data in fid_i.groups.items():
                    grp_reg = grp_month.createGroup(region)
                    for tracer, tracer_data in region_data.groups.items():
                        grp_tracer = grp_reg.createGroup(tracer)
                        grp_tracer.createDimension('n_obs', len(tracer_data.dimensions['n_obs']))
                        var = grp_tracer.createVariable('random_vector', 'd', ('n_obs',))
                        var[:] = error_scale * random.standard_normal(len(tracer_data.dimensions['n_obs']))
                fid_i.close()
        fid.close()

    def StoreModelReprErrors(self):
        # store the model representation error in a file
        repr_error_file = os.path.join(self.output_dir, 'representation_errors.nc4')
        if os.path.exists(repr_error_file):
            os.remove(repr_error_file)
        fid = my_Dataset(repr_error_file, 'w')

        point_split_period = self.rcf.get('output.point.split.period', default='a')
        if point_split_period == 'a':
            point_patt = 'point_output.nc4'
        elif point_split_period == 'm':
            point_patt = 'point_output_??????.nc4'
        elif point_split_period == 'd':
            point_patt = 'point_output_????????.nc4'

        point_list = glob.glob(os.path.join(self.output_dir, 'point', point_patt))
        if self.rcf.get('adjoint.input.point', 'bool') and len(point_list) > 0:
            grp_id = fid.createGroup('point_errors')
            for point_track_file in point_list:
                track_fid = my_Dataset(point_track_file, 'r')
                grp_month = grp_id.createGroup(os.path.basename(point_track_file))
                for region, region_data in track_fid.groups.items():
                    grp_reg = grp_month.createGroup(region)
                    for tracer, tracer_data in region_data.groups.items():
                        if len(tracer_data.dimensions['samples']) > 0:
                            out_grp_tracer = grp_reg.createGroup(tracer)
                            out_grp_tracer.createDimension('samples', len(tracer_data.dimensions['samples']))
                            var = out_grp_tracer.createVariable('mixing_ratio_sigma', 'd', ('samples'))
                            var[:] = tracer_data.variables['mixing_ratio_sigma'][:]
                track_fid.close()

        sat_split_period = self.rcf.get('output.satellite.split.period', default='m')
        if sat_split_period == 'm':
            sat_patt = 'sat-track_??????.nc4'
        elif sat_split_period == 'd':
            sat_patt = 'sat-track_????????.nc4'

        sat_list = glob.glob(os.path.join(self.output_dir, 'satellite', sat_patt))
        if self.rcf.get('adjoint.input.satellite', 'bool') and len(sat_list) > 0:
            grp_id = fid.createGroup('satellite_errors')
            for dep_file in sat_list:
                fid_i = my_Dataset(dep_file, 'r')
                grp_month = grp_id.createGroup(os.path.basename(dep_file))
                for region, region_data in fid_i.groups.items():
                    grp_reg = grp_month.createGroup(region)
                    grp_reg.createDimension('n_lev', len(region_data.dimensions['n_lev']))
                    for tracer, tracer_data in region_data.groups.items():
                        grp_tra = grp_reg.createGroup(tracer)
                        grp_tra.createDimension('n_obs', len(tracer_data.dimensions['n_obs']))
                        var = grp_tra.createVariable('std_profiles', 'd', ('n_obs','n_lev'))
                        var[:] = tracer_data.variables['std_profiles'][:]
                        var = grp_tra.createVariable('sigma_model_column', 'd', ('n_obs',))
                        var[:] = tracer_data.variables['sigma_model_column'][:]
                fid_i.close()
        fid.close()

    def RestoreModelReprErrors(self):
        # restore the model representation errors from a file
        repr_error_file = os.path.join(self.output_dir, 'representation_errors.nc4')
        fid_i = my_Dataset(repr_error_file, 'r')
        if self.rcf.get('adjoint.input.point', 'bool') and 'point_errors' in fid_i.groups:
            gid_i = fid_i.groups['point_errors']
            for fileName, fileData in gid_i.groups.items(): # point_output_*.nc4
                fid_o = my_Dataset(os.path.join(self.output_dir, 'point', fileName), 'a')
                for region, region_data in fileData.groups.items():
                    for tracer, tracer_data in region_data.groups.items():
                        fid_o.groups[region].groups[tracer].variables['mixing_ratio_sigma'][:] = tracer_data.variables['mixing_ratio_sigma'][:]
                fid_o.close()
        if self.rcf.get('adjoint.input.satellite', 'bool') and 'satellite_errors' in fid_i.groups:
            gid_i = fid_i.groups['satellite_errors']
            for fileName, fileData in gid_i.groups.items():
                fid_o = my_Dataset(os.path.join(self.output_dir, 'satellite', fileName), 'a')
                for region, region_data in fileData.groups.items():
                    for tracer, tracer_data in region_data.groups.items():
                        fid_o.groups[region].groups[tracer].variables['std_profiles'][:] = tracer_data.variables['std_profiles'][:]
                        fid_o.groups[region].groups[tracer].variables['sigma_model_column'][:] = tracer_data.variables['sigma_model_column'][:]
                fid_o.close()
        fid_i.close()

    def StoreToStateFile(self, EmisFile=None, iniFile=None, paramFile=None):
        """
        In the following, we are going to implement nonlinear emissions, i.e., emissions which are strictly positive. At the
        beginning of the routine, self.preco.xc_to_x will return emission values in self.state_model which will be both positive
        and negative. So after self.preco.state2struct, the variable 'emissions' will contain all sorts of values. We will
        transform those to positive-only values using the following transform. If x is a Gaussian random variable with mean mu
        and standard deviation sigma, then

        y = <y> sqrt(2*pi)/F(mu/sigma) * log(exp(x/sigma) + 1)

        where <y> is the desired mean of the transformed emission and

        F(z) = \int_{-\infty}^\infty du log(exp(u) + 1) exp(-(u-z)^2/2)

        is a function which we will get from a lookup table, results in a random variable y which is strictly positive. This
        transformation will correspond to optim_type = 2. We can implement other transformations for optim_type > 2. The prior
        emissions, which will serve as mu, are in self.Prior_Emission, whereas the prior emission errors, which will serve as
        sigma, are in self.preco.emis_dapri. This still leaves <y> to be fixed. Since y is in emission space, it makes sense to
        make <y> = mu, i.e., self.Prior_Emission. This leads to a problem, however. If the prior emission is zero over a cell,
        it will never deviate from that value! So for each region/tracer/category combo, let's make <y> equal to the average
        prior emission.

        Maarten's old transformation is coded up as optim_type = 3.
        """
        self.state_model = self.preco.xc_to_x(self.state_preco, self.state_prior_model)
        # Convert to emissions, iniconc and parameters
        emissions,iniconc,parameters = self.preco.state2struct(self.state_model)

        if self.Optim_dict['emission']:
            if EmisFile == None: raise NameError('RunTM5.StoreToStateFile :: no EmisFile supplied')
            if os.path.exists(EmisFile):
               os.remove(EmisFile)
            fid = my_Dataset(EmisFile, 'w')
            fid.createDimension('itime',6)
            for region, xlims, ylims in zip(self.region_names, self.region_limits_lon, self.region_limits_lat):
                group = fid.createGroup(region)
                group.createDimension('latitude', ylims[2])
                group.createDimension('longitude', xlims[2])
                group.latmin = ylims[0]
                group.latmax = ylims[1]
                group.lonmin = xlims[0]
                group.lonmax = xlims[1]
                for tracer in self.species:
                    categories = self.Emission[region][tracer]['categories']
                    tgroup = group.createGroup(tracer)
                    for cat in categories:
                        cgroup = tgroup.createGroup(cat)
                        cgroup.createDimension('nt', emissions[region][tracer][cat]['emission_data'].shape[0])
                        cgroup.time_resolution =  self.Emission[region][tracer][cat]['time_resolution']
                        cgroup.optimize =  int32(self.Emission[region][tracer][cat]['optimize'])

                        var = cgroup.createVariable('time_start', int16, ('nt', 'itime'))
                        var[:] = array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_start']], int16)

                        var = cgroup.createVariable('time_end', int16, ('nt', 'itime'))
                        var[:] = array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_end']], int16)

                        var = cgroup.createVariable('emission', emissions[region][tracer][cat]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))

                        if not self.Emission[region][tracer][cat]['optimize']:
                            var[:] = emissions[region][tracer][cat]['emission_data']

                        else:
                            if self.preco.optim_type[tracer][cat] == 1:
                                var[:] = emissions[region][tracer][cat]['emission_data']

                            elif self.preco.optim_type[tracer][cat] == 2: # SB's transformation
                                var[:] = emissions[region][tracer][cat]['emission_data']

                            elif self.preco.optim_type[tracer][cat] == 3: # MK's transformation
                                emisx = emissions[region][tracer][cat]['emission_data']

                                # I used to have a simple construct, where(emisx < 0.0, exp(emisx), 1.0 + emisx), but that
                                # evaluated the exponential for all values of emisx, even large positive ones, which were
                                # unnecessary, but generated overflow errors. Therefore, I first need to create two temporary
                                # arrays separating the positive from the negative elements of emisx, and only take the
                                # exponential of the negative part. What a pain! But there is currently no way to do lazy
                                # evaluation of numpy.where, i.e., only evaluate the necessary elements.
                                positive_part = where(emisx >= 0.0, emisx, 0.0)
                                negative_part = where(emisx < 0.0, emisx, 0.0)
                                emisx = self.Prior_Emission[region][tracer][cat]['emission_data'] * \
                                    where(emisx < 0.0, exp(negative_part), 1.0+positive_part)
                                del positive_part, negative_part

                                var[:] = emisx
                                # keep this non-negative emission in the self.Emission structure
                                self.Emission[region][tracer][cat]['emission_data'] = emisx

                        # If we need to write out the emission errors, do that
                        write_err = self.rcf.get('%s.%s.sep_error'%(tracer, cat), 'bool', default=False)
                        if write_err:
                            var = cgroup.createVariable('emission_error', self.preco.emis_dapri[region][tracer][cat]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                            var[:] = self.preco.emis_dapri[region][tracer][cat]['emission_data']

            fid.close()

        if self.Optim_dict['iniconc']:
            if iniFile == None:
                raise NameError('RunTM5.StoreToStateFile :: no iniFile supplied')
            else:
                raise KeyError('Optimization of initial concentration not yet implemented')

        if self.Optim_dict['parameters']:
            if paramFile == None:
                raise NameError('RunTM5.StoreToStateFile :: no paramFile supplied')
            if os.path.isfile(paramFile):
                os.remove(paramFile)
            checkDir(paramFile)
            savetxt(paramFile, parameters)

    def ReadPriorState(self, iniFile=None, paramFile=None):
        emissions = None
        iniconc   = None
        biasParameters = None

        if self.Optim_dict['emission']:
            # Emissions are already in struct: self.Emission[region][tracer][cat]['emission_data']. However, for simplicity of
            # coding, we will read in the prior emissions in self.Prior_Emission. This will help in Precon.state2struct,
            # although strictly speaking we only need this for nonlinear emissions.
            self.Prior_Emission = dict.fromkeys(self.region_names)
            for region in self.region_names:
                self.Prior_Emission[region] = dict()
                for tracer in self.species:
                    self.Prior_Emission[region][tracer] = dict()
                    for cat_name in self.Emission[region][tracer]['categories']:
                        self.Prior_Emission[region][tracer][cat_name] = dict()
                        self.Prior_Emission[region][tracer][cat_name]['emission_data'] = \
                            zeros_like(self.Emission[region][tracer][cat_name]['emission_data'])
                        self.Prior_Emission[region][tracer][cat_name]['emission_data'][:] = \
                            self.Emission[region][tracer][cat_name]['emission_data']

                        if self.Emission[region][tracer][cat_name]['optimize']:
                            optim_type = self.preco.optim_type[tracer][cat_name]
                            if optim_type == 3: # For MK's transformation
                                # Maarten says that the following line is needed
                                self.Emission[region][tracer][cat_name]['emission_data'][:] = 0.0

        if self.Optim_dict['iniconc']:
            if iniFile == None:
                raise NameError('RunTM5.ReadPriorState :: no iniFile supplied')
            else:
                raise KeyError('Optimization of initial concentration not yet implemented')

        if self.Optim_dict['parameters']:
            if paramFile == None:
                raise NameError('RunTM5.ReadPriorState :: no paramFile supplied')
            biasParameters = loadtxt(paramFile, dtype=float64)

        # Apply struct2state
        self.state_prior_model = self.preco.struct2state(self.Emission,iniconc,biasParameters)
        self.state_preco = zeros_like(self.state_prior_model)
        # We will store the prior in preco space for possible later use. I know that normally this is a string of zeros,
        # but for MC congrad it is not!
        self.state_prior_preco = zeros_like(self.state_prior_model)

    def GetZoomRegions(self):
        """
        Reads a config file to get zoom region info
        """
        self.region_names = self.rcf.get('regions').split()
        self.nregions = len(self.region_names)
        self.xbeg = []
        self.xend = []
        self.im = []
        self.ybeg = []
        self.yend = []
        self.jm = []
        for region in self.region_names:
            self.xbeg.append(self.rcf.get('region.'+region+'.xbeg','float'))
            self.ybeg.append(self.rcf.get('region.'+region+'.ybeg','float'))
            self.xend.append(self.rcf.get('region.'+region+'.xend','float'))
            self.yend.append(self.rcf.get('region.'+region+'.yend','float'))
            self.im.append(self.rcf.get('region.'+region+'.im','int'))
            self.jm.append(self.rcf.get('region.'+region+'.jm','int'))
        self.region_limits_lon = array(zip(self.xbeg,self.xend,self.im))
        self.region_limits_lat = array(zip(self.ybeg,self.yend,self.jm))

    def determine_children_etc(self):
        """
        This routine calculates the geometry for the TM5 model. How many
        zoom regions, patents, children, where do they start & end.
        """
        zoomfile = os.path.join(self.rcf.get('my.run.dir'), 'Zoomed.nc4')
        region_names = self.rcf.get('regions').split()
        self.nregions = len(region_names)
        self.region_names = []
        for i in range(self.nregions):
            self.region_names.append(region_names[i].strip())
        if os.path.isfile(zoomfile):
            f = my_Dataset(zoomfile,'r')
            self.xbeg = f.variables['xbeg'][:]
            self.xend = f.variables['xend'][:]
            self.ybeg = f.variables['ybeg'][:]
            self.yend = f.variables['yend'][:]
            self.xref = f.variables['xref'][:]
            self.yref = f.variables['yref'][:]
            self.im   = f.variables['im'][:]
            self.isr  = f.variables['isr'][:]
            self.ier  = f.variables['ier'][:]
            self.jm   = f.variables['jm'][:]
            self.jsr  = f.variables['jsr'][:]
            self.jer  = f.variables['jer'][:]
            self.ibeg = f.variables['ibeg'][:]
            self.iend = f.variables['iend'][:]
            self.jbeg = f.variables['jbeg'][:]
            self.jend = f.variables['jend'][:]
            self.parents  = f.variables['parents'][:]
            self.children = []
            self.zoomed   = {}
            self.edge     = {}
            self.region_limits_lon = array([l for l in zip(self.xbeg,self.xend,self.im)])
            self.region_limits_lat = array([l for l in zip(self.ybeg,self.yend,self.jm)])
            self.region_xref = self.xref
            self.region_yref = self.yref
            self.region_im = self.im
            self.region_jm = self.jm
            for i in range(self.nregions):
                self.children.append(f.variables['children_%2.2i'%(i+1)][:])
                grp = f.groups[self.region_names[i]]
                self.zoomed[self.region_names[i]] = grp.variables['zoomed'][:]
                self.edge[self.region_names[i]] = grp.variables['edge'][:]
            f.close()
        else:
            sys.stderr.write('Zoomed file is created in forward run, do this first\n')
            sys.exit(2)

    def save(self):
        state_dir = os.path.join(self.output_dir, 'state')
        # if the 'state' folder does not exist, create it
        if not os.path.isdir(state_dir):
            os.makedirs(state_dir)

        file_name = os.path.join(state_dir, 'RunTM5.state')
        if os.path.exists(file_name):
            os.remove(file_name)

        # Do not dump the entire self.__dict__ at once, that creates problems loading when it's really large
        with open(file_name, 'wb') as fid:
            for k, v in self.__dict__.items():

                # If v is an instance of Precon, invoke its own save() method
                if isinstance(v, Precon):
                    v.save()

                # If v is a TM5_Emission object, call its own save() method
                elif k == 'Emission':
                    state_file = os.path.join(state_dir, 'emission_state.h5')
                    v.save(state_file)

                # If v is a numpy array, save it as an npy file
                elif isinstance(v, ndarray):
                    fname = os.path.join(state_dir, 'runtm5_%s.npy'%k)
                    save(fname, v)

                # Otherwise, pickle it and save it
                else:
                    pickle.dump((k, v), fid, pickle.HIGHEST_PROTOCOL)

    def load(self, cleanup=True):
        state_dir = os.path.join(self.output_dir, 'state')

        file_name = os.path.join(state_dir, 'RunTM5.state')
        if not os.path.exists(file_name):
            logging.error('%s not found'%file_name)
            sys.exit()

        # Load the simple, pickled things first
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
        npy_files = glob.glob(os.path.join(state_dir, 'runtm5_*.npy'))
        for file_name in npy_files:
            base_name = os.path.basename(file_name)
            obj_name = os.path.splitext(base_name)[0].split('_', 1)[1]
            setattr(self, obj_name, load(file_name))
            # Delete the npy file
            if cleanup:
                os.remove(file_name)

        # Load self.Emission
        self.Emission = TM5_emission()
        state_file = os.path.join(state_dir, 'emission_state.h5')
        self.Emission.load(state_file, cleanup=cleanup)

        # Load the Precon object
        self.preco = Precon(self, resume=True, cleanup=cleanup)
