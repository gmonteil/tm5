#!/usr/bin/env python2.7

import os
from dateutil.relativedelta import relativedelta
import glob
from numpy import array, dot, float64, zeros_like, loadtxt, int16, where, exp, int32, savetxt, linalg, zeros
from datetime import datetime
import sys
import logging
import subprocess
from pyshell.emissions import Emissions
from pyshell.observations import PointObs, Observations
from pyshell.utilities import checkDir, my_Dataset
from pyshell.precon import Precon
logging.basicConfig(datefmt='%d %b %Y %H:%M:%S', format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


class ExecEnvironment(object):
    # A class to make sure that foreground runs come back to the current folder
    def __init__(self, cur_dir, exec_dir):
        self.cur_dir = cur_dir
        self.exec_dir = exec_dir

    def __enter__(self):
        os.chdir(self.exec_dir)

    def __exit__(self, type, value, tb):
        os.chdir(self.cur_dir)


class RunTM5:
    def __init__(self, rcf, dconf):
        self.rcf = rcf
        self.dconf = dconf
        # The key 'my.tracer' is a comma-separated list of tracer, whereas the key 'my.tracer.name' is a name for the
        # combination. For example, my.tracer could be 'CO,CO2' and my.tracer.name could be 'COCO2'.
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.runid = self.rcf.get('runid')
        self.StartTime = rcf.ti
        self.EndTime = rcf.tf
        self.subdir_tag = ''
        # SBi : Check whether all the following keys are necessary
        self.rcf.add('timerange.start', self.StartTime.strftime("%Y-%m-%d %H:%M:%S"))
        self.rcf.add('timerange.end', self.EndTime.strftime("%Y-%m-%d %H:%M:%S"))
        # SBf
        self.J_bg = 0.0
        self.J_obs = 0.0
        self.J_tot = self.J_bg + self.J_obs
        self.Optim_dict = {'emission': 'optimize.emission', 'iniconc': 'optimize.initialconcentration',
                           'parameters': 'optimize.parameter'}
        for k, v in self.Optim_dict.items():
            self.Optim_dict[k] = self.rcf.get(v, 'bool')
        self.months = []
        d = self.StartTime.replace(day=1).replace(hour=0).replace(minute=10)
        while d <= self.EndTime:
            self.months.append(d.strftime('%Y%m'))
            d += relativedelta(months=1)
        self.region_names = self.rcf.get('regions').split()
        self.nregions = len(self.region_names)
        self.output_dir = self.rcf.get('output.dir')
        self.emission_file_name = self.rcf.get('PyShell.em.filename')
        # delete all the tm5-pyshell.*.rc files in the output directory
        files_to_delete = glob.glob(os.path.join(self.rcf.get('my.run.dir'), self.rcf.get('my.basename') + '.*.rc'))
        for fileName in files_to_delete:
            try:
                os.remove(fileName)
            except:
                print('File ', fileName, ' has been already deleted!')
        self.GetZoomRegions()

    def SetupObservations(self, obsobj):
        ti = self.StartTime
        dt = self.EndTime - ti
        while ti < self.EndTime:
            obs = Observations(ti, ti + dt, self.rcf)
            for tracer in obs.species:
                obs.SetupPointObs(tracer, obsobj)
            obs.writePointFile()
            ti = ti + dt

    def SetupEmissions(self, emclasses, randomize=False, zero=False, step=None):
        """
        This will either assemble the emissions (default), or read them from a file, depending on
        an rc key, emission.read.optimized.
        """

        t1 = datetime.now()

        emisFile = self.emission_file_name

        # First, process the emissions and write the emissions file
        # StartTuple = self.StartTime.timetuple()[:5]
        # EndTuple = self.EndTime.timetuple()[:5]
        em = Emissions(self.rcf)
        em.create_emission_structure()

        read_optim_emis = self.rcf.get('emission.read.optimized', 'bool', default=False)
        if read_optim_emis :
            optim_emis_file = self.rcf.get('emission.read.optimized.filename')
            if os.path.basename(optim_emis_file) == optim_emis_file :
                optim_emis_file = os.path.join(self.rcf.get('output.dir'), optim_emis_file)
            em.readOptimFromFile(optim_emis_file, False)
        else :
            for tracer in em.species:
                emfill = emclasses[tracer](self.rcf, self.dconf.emissions[tracer], step=step, tracer=tracer)
                emfill.Emission = em.Emission
                emfill.LoopThroughPeriods()

        # Remove emission file that may exist AFTER it has been read
        if os.path.exists(emisFile):
            os.remove(emisFile)
        checkDir(emisFile)

        self.Emission = em(randomize, zero)  # Write emissions and keep them in memory

        for region in self.region_names:
            for tracer in self.species:
                categories = self.Emission[region][tracer]['categories']
                for cat in categories:
                    print(region, tracer, cat, self.Emission[region][tracer][cat]['time_resolution'], ' Optimize?:', self.Emission[region][tracer][cat]['optimize'], self.Emission[region][tracer][cat]['emission_data'].shape)
                for cat_key in self.Emission[tracer]['cat_list']:
                    re = self.Emission[tracer]['cat_opt'][cat_key]
                    print('========================================')
                    print('Tracer: ', tracer, 'Category marked for optimization:', cat_key)
                    for rei in re:
                        print('----->region:', rei['region'], 'error:', rei['error'])
                    print('========================================')
        print('========================================')
        print('end info parsed from emission routine')
        print('========================================')

        t2 = datetime.now()
        sys.stderr.write("Emission between %s and %s assembled in %s\n" % (self.StartTime.strftime("%c"), self.EndTime.strftime("%c"), str(t2 - t1)))

    def RunForward(self):
        """
        Run TM5 via pycasso scripting: build, make, submit.
        """

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

    def RunBackward(self):
        """
        Run adjoint TM5 via pycasso scripting: build, make, submit.
        """
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
        logger.info("Adjoint run from %s to %s took %s on %s threads\n" % (self.StartTime.strftime("%c"), self.EndTime.strftime("%c"), str(t2-t1), os.getenv('OMP_NUM_THREADS')))

    def CalculateGradient(self, add_bg=True):
        adjoint_emission_model, adjoint_iniconc_model, adjoint_parameter_model = self.ReadAdjointState()
        logger.info('Adjoint state read\n')
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

    def CalculateCosts(self):
        self.J_bg = self.backgroundCost()
        logger.info('Background cost calculated\n')
        self.J_obs = self.ObservationCost()
        logger.info('Observation cost calculated\n')
        self.J_tot = self.J_bg + self.J_obs
        return self.J_tot

    def backgroundGradient(self):
        # The portion of the gradient of the cost function from the background term, in preco space
        return self.state_preco - self.state_prior_preco

    def backgroundCost(self):
        # The cost function from the background term
        state_vec_diff = self.state_preco - self.state_prior_preco
        return 0.5*sum(state_vec_diff*state_vec_diff)

    def ObservationCost(self):
        J_obs = 0.0
        if self.rcf.get('adjoint.input.satellite', 'bool'):
            raise NotImplementedError
            # satobs = ApplySatelliteObs(self)
            # logger.info('Satellite departures calculated\n')
            # J_obs += satobs.SatelliteCost()
            # logger.info('Satellite cost function calculated\n')
        if self.rcf.get('adjoint.input.point', 'bool'):
            pointobs = PointObs(self)
            logger.info('Point departures calculated\n')
            J_obs += pointobs.PointCost()
            logger.info('Point cost function calculated\n')
        return J_obs

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
            logger.error('Not implemented yet, we should open the file BACKWARD_SAVE for this run... ReadAdjointState, class RunTM5...')
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

    def makePrecon(self):
        # For coding ease, it's better to create the Precon object here instead of in the Optimizer routines
        self.preco = Precon(self)

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

    def ReadPriorState(self, iniFile=None, paramFile=None):
        emissions = None
        iniconc = None
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
            if iniFile is None:
                raise NameError('RunTM5.ReadPriorState :: no iniFile supplied')
            else:
                raise KeyError('Optimization of initial concentration not yet implemented')

        if self.Optim_dict['parameters']:
            if paramFile is None:
                raise NameError('RunTM5.ReadPriorState :: no paramFile supplied')
            biasParameters = loadtxt(paramFile, dtype=float64)

        # Apply struct2state
        self.state_prior_model = self.preco.struct2state(self.Emission, iniconc, biasParameters)
        self.state_preco = zeros_like(self.state_prior_model)
        # We will store the prior in preco space for possible later use. I know that normally this is a string of zeros,
        # but for MC congrad it is not!
        self.state_prior_preco = zeros_like(self.state_prior_model)

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
            if EmisFile is None:
                raise NameError('RunTM5.StoreToStateFile :: no EmisFile supplied')
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
                        cgroup.time_resolution = self.Emission[region][tracer][cat]['time_resolution']
                        cgroup.optimize = int32(self.Emission[region][tracer][cat]['optimize'])

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
            if iniFile is None:
                raise NameError('RunTM5.StoreToStateFile :: no iniFile supplied')
            else:
                raise KeyError('Optimization of initial concentration not yet implemented')

        if self.Optim_dict['parameters']:
            if paramFile is None:
                raise NameError('RunTM5.StoreToStateFile :: no paramFile supplied')
            if os.path.isfile(paramFile):
                os.remove(paramFile)
            checkDir(paramFile)
            savetxt(paramFile, parameters)

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

    def PointDepartures(self):
        dp = PointObs(self)
        for month_tuple in dp.months:
            dp.create_pointdeparture_structure()
            for tracer in dp.species:
                dp.applyPointObs(tracer, month_tuple)
            dp.writePointDepartureFile(month_tuple)

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
            logger.error('Zoomed file is created in forward run, do this first\n')
            sys.exit(2)
