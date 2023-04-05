#!/usr/bin/env python2.7

import os
import sys
from netCDF4 import Dataset
from numpy import zeros, where, float64, dot, loadtxt, transpose, ones, average, maximum, ones_like, zeros_like, eye, meshgrid, arange, linalg, argsort, diag, sqrt, array, load, linspace, exp, shape
from tqdm import tqdm
from collections import defaultdict
import logging
import pickle
import glob
from datetime import datetime
from copy import deepcopy


logger = logging.getLogger(__name__)


def xc_to_x(G_state, Temp_L, Hor_L, x_c, ipos):
    n_state = len(G_state)
    nt = shape(Temp_L)[0]
    nhor = shape(Hor_L)[0]
    x = zeros(n_state)
    for i in tqdm(range(nt), desc='xc_to_x', leave=True):
        for j in tqdm(range(nt), desc='step %i/%i'%(i, nt), leave=False):
            x[ipos+i*nhor:ipos+(i+1)*nhor] += G_state[ipos+i*nhor:ipos+(i+1)*nhor]* dot(Temp_L[i,j]*Hor_L, x_c[ipos+j*nhor:ipos+(j+1)*nhor])
    return x


def g_to_gc(G_state, Temp_Lt, Hor_Lt, g, ipos):
    n_state = len(G_state)
    nt = shape(Temp_Lt)[0]
    nhor = shape(Hor_Lt)[0]
    g_c = zeros([n_state])
    for i in tqdm(range(nt), desc='preconditioning gradient', leave=False):
        for j in range(nt):
            g_c[ipos+i*nhor:ipos+(i+1)*nhor] += dot(Temp_Lt[i,j]*Hor_Lt, G_state[ipos+j*nhor:ipos+(j+1)*nhor] * g[ipos+j*nhor:ipos+(j+1)*nhor])
    return g_c


class Precon:

    def __init__(self, tm5Obj, resume=False, cleanup=True):

        self.output_dir = tm5Obj.output_dir
        self.tm5 = tm5Obj

        if resume:
            self.load(cleanup=cleanup)

        else:
            self.copyRunObjects()
            self.hor_corr_files = defaultdict(dict)
            self.hor_corr_choices = defaultdict(dict)
            # Read variables needed for zoom...
            self.GetZoomRegions()
            self.ReadRegionDefinitions()
            # Read the topo file for separating between land and sea regions
            topo_file = self.rcf.get('topo.database')
            self.createLandSeaMask(topo_file)
            #  Note MK: here we read a new parameter to allow for non-linear optimizations (optim_emis.type > 1)
            #  Note MK: all changes to the code will contain this optim_type.
            #  Note SB: Modified so that some (but not all) tracers can have nonlinear emissions
            #  Note SB: Actually, we want to be able to specify that only some categories have nonlinear emissions
            self.optim_type = dict.fromkeys(self.species)
            for tracer in self.species:
                self.optim_type[tracer] = dict()
                for cat_key in self.Emission[tracer]['cat_list']: # only categories that will be optimized
                    cat_name = cat_key[0]
                    self.optim_type[tracer][cat_name] =  \
                        self.rcf.get('var4d.%s.%s.optim_emis.type'%(tracer,cat_name), 'int', default=1)
            self.setup_corr()
            # What is the maximum number of OpenMP threads?
            try:
                self.max_threads = self.rcf.get('par.maxthreads', 'int')
            except :
                self.max_threads = int(os.environ['OMP_NUM_THREADS'])

    def xc_to_x(self, x_c, x_prior):
        # xc_to_c converts the preconditioned space to normal space
        x = zeros_like(x_prior)
        x += x_prior
        if self.Optim_dict['emission']:
            t1 = datetime.now()
            ipos = 0
            for tracer in self.species:
                for c,cat_key in enumerate(self.Emission[tracer]['cat_list']):
                    nt = self.Emission[tracer]['cat_nt'][c]
                    n_hor = self.Emission[tracer]['cat_n_hor'][c]
                    Hor_L = self.Emission[tracer]['cat_Hor_L'][c]
                    Temp_L = self.Emission[tracer]['cat_Temp_L'][c]
                    cat_blk = nt*n_hor
                    x += xc_to_x(self.G_state, Temp_L, Hor_L, x_c, ipos)
                    ipos += cat_blk
            t2 = datetime.now()
            sys.stderr.write("Emission preconditioned in %s\n"%str(t2-t1))

        if self.Optim_dict['iniconc']:
            logger.error('Iniconc preconditioning is not implemented yet in xc_to_x in class Precon...')
            sys.exit(0)

        if self.Optim_dict['parameters']:
            x[-self.n_param:] += self.G_state[-self.n_param:]*dot(self.Bias_L,x_c[-self.n_param:])

        return x

    def g_to_gc(self, g): # Pim: add stuff for adj_iniconc
        g_c = zeros_like(g)
        if self.Optim_dict['emission']:
            t1 = datetime.now()
            ipos = 0
            for tracer in self.species:
                for c,cat_key in enumerate(self.Emission[tracer]['cat_list']):
                    nt =    self.Emission[tracer]['cat_nt'][c]
                    n_hor = self.Emission[tracer]['cat_n_hor'][c]
                    Hor_Lt = transpose(self.Emission[tracer]['cat_Hor_L'][c])
                    Temp_Lt = transpose(self.Emission[tracer]['cat_Temp_L'][c])
                    cat_blk = nt*n_hor
                    g_c += g_to_gc(self.G_state, Temp_Lt, Hor_Lt, g, ipos)
                    ipos += cat_blk
                    del Hor_Lt # save memory

            t2 = datetime.now()
            sys.stderr.write("Adjoint emission preconditioned in %s\n"%str(t2-t1))

        if self.Optim_dict['iniconc']:
            logger.error('Iniconc preconditioning is not implemented yet in g_to_gc in class Precon...')
            sys.exit(0)

        if self.Optim_dict['parameters']:
            g_c[-self.n_param:] = dot(self.Bias_Lt,self.G_state[-self.n_param:]*g[-self.n_param:])

        return g_c

    def ReadRegionDefinitions(self):
        # This routine will read Zoomed.nc4 (in output folder) only after a forward run...
        zoomfile = os.path.join(self.rcf.get('my.zoom.dir'), 'Zoomed_%s.nc4'%self.rcf.get('my.zoom'))
        if os.path.isfile(zoomfile):
            f = Dataset(zoomfile,'r')
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
            for i in range(self.nregions):
                self.children.append(f.variables['children_%2.2i'%(i+1)][:])
                grp = f.groups[self.region_names[i]]
                self.zoomed[self.region_names[i]] = grp.variables['zoomed'][:]
                self.edge[self.region_names[i]] = grp.variables['edge'][:]
            f.close()
        else:
            raise RuntimeError('Zoomed file %s not yet created, do a forward run first'%zoomfile)

    def copyRunObjects(self):
        # Copy over RunTM5 attributes onto self
        self.rcf = self.tm5.rcf
        self.Emission = self.tm5.Emission
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.Optim_dict = self.tm5.Optim_dict
        self.StartDate = self.tm5.StartTime
        self.EndDate = self.tm5.EndTime
        self.emis_dapri_file = self.tm5.emission_file_name

    def read_latlon(self, file_name):
        if not os.path.exists(file_name):
            raise RuntimeError("%s does not exist"%file_name)
        f = Dataset(file_name)
        P = f.variables['P'][:]
        D = f.variables['sqrt_lam'][:]
        f.close()
        # Thanks to Fortran's idiosyncracy, the indices of arrays in netcdf files are swapped
        # So P is actually transpose of the eigenvector matrix, and we need to transpose it back
        return transpose(P), D

    def load(self, cleanup=True):
        state_dir = os.path.join(self.output_dir, 'state')

        file_name = os.path.join(state_dir, 'Precon.state')
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

        # Resume link to RunTM5 attributes
        self.copyRunObjects()

        # Load the npy files
        npy_files = glob.glob(os.path.join(state_dir, 'precon_*.npy'))
        for file_name in npy_files:
            base_name = os.path.basename(file_name)
            obj_name = os.path.splitext(base_name)[0].split('_', 1)[1]
            setattr(self, obj_name, load(file_name))
            # Delete the npy file
            if cleanup:
                os.remove(file_name)

    def createLandSeaMask(self, topo_file):
        """
        When we optimize CO emissions we have a problem with the minimum error per grid cell. We do not want to
        set it to zero, since that will preclude biomass burning emissions from cells where they do not exist.
        At the same time, if we set it to something positive then ocean pixels end up getting biomass burning
        emissions. So what we are going to do is read in a land sea mask, and then set a non-zero minimum emission
        error only over the land/sea pixels, depending on the category. The category will be determined by the
        last keyword (currently unused) in the emission specification line, such as

        CO.region3.category2    : biomass burning ; 250.0  ;  200.0-e ; 0.03-e-daily   ; 1 ; def-def-0
                                                                                             ^^^^^^^^^

        Instead of 'def-def-0', we will put 'land' or 'sea' to specify whether a certain category is terrestrial
        or oceanic. For land categories we will enforce zero prior emission error over the ocean, and vice versa.
        Reading in of the 'def-def-0' keyword will have to be done by the __init__ method of Emissions_base.py.
        This is done, and the variable to look at is self.Emission[region][cat]['remarks']. Further, a land-sea
        mask will have to be read in.
        """
        with Dataset(topo_file, 'r') as fid:
            topo = fid.variables['TerrainHeight'][:]
            ocean_value = fid.Ocean_value
        # the ocean pixels are ocean_value, everything else is land
        ocean_mask = zeros(topo.shape, float64)
        ocean_mask[where(topo == ocean_value)] = 1.0
        # we want masks for each zoom region
        self.ocean_masks = {}
        self.land_masks = {}
        for region, lat_spec, lon_spec in zip(self.zoom_regions_names, self.zoom_regions_lat, self.zoom_regions_lon):
            # what is the resolution of this zoom region?
            dlat = (lat_spec[1]-lat_spec[0])/lat_spec[2]
            dlon = (lon_spec[1]-lon_spec[0])/lon_spec[2]
            # what is the resolution of the topo file?
            mask_dlat = 180.0/ocean_mask.shape[0]
            mask_dlon = 360.0/ocean_mask.shape[1]
            # therefore, what is the coarsening factor?
            coar_fac_lat = int(dlat/mask_dlat)
            coar_fac_lon = int(dlon/mask_dlon)
            self.ocean_masks[region] = zeros((lat_spec[2], lon_spec[2]), float64)
            # select the correct portion of the ocean mask
            mask_lats = linspace(-90.,90.,ocean_mask.shape[0]+1)
            mask_lons = linspace(-180.,180.,ocean_mask.shape[1]+1)
            lat_start_index = mask_lats.searchsorted(lat_spec[0])
            lat_end_index = mask_lats.searchsorted(lat_spec[1])
            lon_start_index = mask_lons.searchsorted(lon_spec[0])
            lon_end_index = mask_lons.searchsorted(lon_spec[1])
            ocean_chunk = ocean_mask[lat_start_index:lat_end_index, lon_start_index:lon_end_index]
            # now create the region-specific ocean mask
            for j in range(lat_spec[2]):
                for i in range(lon_spec[2]):
                    chunk = ocean_chunk[coar_fac_lat*j:coar_fac_lat*(j+1), coar_fac_lon*i:coar_fac_lon*(i+1)]
                    self.ocean_masks[region][j,i] = chunk.mean()
            self.land_masks[region] = 1.0 - self.ocean_masks[region]

    def setup_corr(self):
        # Read the a priori standard deviation
        self.G_state = self.read_apri_std(self.emis_dapri_file)
        # Sourish, April 2011: Add funtionality for the different categories to have different correlation lengths
        # Read horizontal correlation matrix for current region setup
        if self.Optim_dict['emission'] or self.Optim_dict['iniconc']:
            for tracer in self.species:
                # Get the necessary correlation lengths
                self.Emission[tracer]['cat_Hor_L'] = []
                for hor_cor_file in self.Emission[tracer]['cat_Bh_file']:
                    P_h,D_h = self.read_latlon(hor_cor_file)
                    Hor_L = P_h * D_h
                    self.Emission[tracer]['cat_Hor_L'].append(Hor_L)
                    del P_h, D_h # save memory

        # Construct temporal correlation matrices with different temporal correlation length per category
        if self.Optim_dict['emission']:
            for tracer in self.species:
                cat_Temp_L = []
                cat_nt = []
                for cat_key in self.Emission[tracer]['cat_list']:
                    tcorr = cat_key[2]  # time-correlation key 9.50-e-monthly
                    cat = cat_key[0]
                    re = self.Emission[tracer]['cat_opt'][cat_key]  # info about regions and errors
                    region = re[0]['region']   # this info is same for possible other regions...
                    time_interval = self.Emission[region][tracer][cat]['time_interval']
                    res_key = self.Emission[region][tracer][cat]['time_resolution']
                    temp_corlen = float(tcorr[:4].strip())    # actually unit is month:
                    if res_key.find('monthly') != -1:
                        dt = time_interval['dt'].months
                    elif res_key.find('daily') != -1:
                        dt = time_interval['dt'].days/30.0   # assumes 30 day in a month
                    else:
                        print("invalid time resolution:", res_key)
                        sys.exit()
                    # note: routine calc_temp_corr assumes time-steps of one month
                    # now we allow for longer/shorter steps also. Expressed in months units it is now
                    # in variable dt. t-corr is calculated as corr = exp(-(t2-t1)/temp_corlen)
                    nt = len(time_interval['time_start'])
                    Temp_L  = zeros((nt,nt))
                    P_t, D_t = self.calc_temp_corr(temp_corlen,dt,nt)
                    Temp_L = dot(P_t, D_t)
                    cat_Temp_L.append(Temp_L)
                    cat_nt.append(nt)

                self.Emission[tracer]['cat_Temp_L'] = cat_Temp_L
                self.Emission[tracer]['cat_nt'] = cat_nt

        # Construct spatial correlation matrix for bias parameters
        if self.Optim_dict['parameters']:
            # Given a bias covariance matrix B_b, we need the LL^T decomposition of B_b
            corr_file_name = self.rcf.get('parameter.correlation.filename')
            corr_matrix = loadtxt(corr_file_name, dtype=float64)
            self.n_param = corr_matrix.shape[0]
            P_bias, D_bias = self.matrix_square_root(corr_matrix)
            self.Bias_L = dot(P_bias,D_bias)
            self.Bias_Lt = transpose(self.Bias_L)

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
        for region in self.zoom_regions_names:
            xbeg.append(self.rcf.get('region.'+region+'.xbeg', 'float'))
            ybeg.append(self.rcf.get('region.'+region+'.ybeg', 'float'))
            xend.append(self.rcf.get('region.'+region+'.xend', 'float'))
            yend.append(self.rcf.get('region.'+region+'.yend', 'float'))
            im.append(self.rcf.get('region.'+region+'.im', 'int'))
            jm.append(self.rcf.get('region.'+region+'.jm', 'int'))
        self.zoom_regions_lon = zip(xbeg, xend, im)
        self.zoom_regions_lat = zip(ybeg, yend, jm)
        self.region_names = self.zoom_regions_names
        self.nregions = len(self.region_names)

    def read_apri_std(self, apri_file):
        if self.Optim_dict['emission']:
            self.emis_dapri = dict.fromkeys(self.region_names)
            for region in self.region_names:
                self.emis_dapri[region] = dict.fromkeys(self.species)
                for tracer in self.species:
                    categories = self.Emission[region][tracer]['categories']
                    self.emis_dapri[region][tracer] = dict.fromkeys(categories)
                    for cat in categories:
                        error = self.Emission[region][tracer][cat]['error']
                        emission = self.Emission[region][tracer][cat]['emission_data']
                        separate_error = self.rcf.get('%s.%s.sep_error'%(tracer,cat), 'bool', default=False)
                        if separate_error:
                            grid_error = self.Emission[region][tracer][cat]['emission_error']
                        else:
                            grid_error = abs(emission) * 0.01 * error
                        if self.Emission[region][tracer][cat]['optimize'] == 1:
                            # is this a land or a sea category? mask the minimum error accordingly
                            err_mask = ones(emission.shape, float64)
                            if self.Emission[region][tracer][cat]['remarks'] == 'land':
                                for i in range(err_mask.shape[0]):
                                    err_mask[i] = self.land_masks[region]
                            elif self.Emission[region][tracer][cat]['remarks'] == 'ocean':
                                for i in range(err_mask.shape[0]):
                                    err_mask[i] = self.ocean_masks[region]

                            if self.optim_type[tracer][cat] == 1:
                                # There might be some categories that want to specify whether minimum error should be
                                # enforced or not. E.g., maybe we want to have a minimum error for the biosphere but
                                # not for fossil fuel emissions.
                                key_name = '%s.%s.enforce.minimum.error'%(tracer, cat)
                                enforce_minerror = self.rcf.get(key_name, 'bool', default=False)
                                if enforce_minerror:
                                    # We might want to have different minimum errors for different categories, so check that
                                    key_name = '%s.%s.minimum.error'%(tracer, cat)
                                    min_std = self.rcf.get(key_name, 'float')
                                    key_name = '%s.%s.minimum.error.depends.average'%(tracer, cat)
                                    min_std_average_flag = self.rcf.get(key_name, 'bool')
                                    if min_std_average_flag:
                                        self.emis_dapri[region][tracer][cat] = {'emission_data': err_mask * maximum(grid_error, min_std * average(abs(emission)) * ones_like(emission))}
                                    else:
                                        self.emis_dapri[region][tracer][cat] = {'emission_data': err_mask * maximum(grid_error, min_std * ones_like(grid_error))}
                                else:
                                    self.emis_dapri[region][tracer][cat] = {'emission_data': err_mask * grid_error}

                            elif self.optim_type[tracer][cat] == 2: # SB's transformation
                                self.emis_dapri[region][tracer][cat] = {'emission_data': err_mask * grid_error}

                            elif self.optim_type[tracer][cat] == 3: # MK/PB's transformation
                                #   the problem is the following:  The emissions are being transformed and we start from a
                                # state = 0 that is transformed to emissions using Emission = Prior_Emission * exp(state) if
                                # state < 0, and Emission = Prior_Emission (1 + state) for state > 0. So we optimize state using
                                # a start value = 0 and an error. Now the error does not make sense when Prior_Emission = 0,
                                # best thing is to set the dapri to error everywhere, e.g. 250*0.01 ...this optimizes state in
                                # -2.5 <----0----> 2.5 and the emissions in Prior_emission*exp(-2.5) <---- Prior_emission --->
                                # Prior_Emission*(1 + 2.5)
                                self.emis_dapri[region][tracer][cat] = {'emission_data': err_mask * ones_like(grid_error)*error*0.01}

                            else: # optim_type is not 1, 2 or 3, but code below is identical to 3
                                self.emis_dapri[region][tracer][cat] = {'emission_data': err_mask * ones_like(emission)*error*0.01}

                        else: # category not optimized
                            self.emis_dapri[region][tracer][cat] = {'emission_data': zeros_like(emission)}
        else:
            self.emis_dapri = None

        if self.Optim_dict['parameters']:
            dapri_file_name = self.rcf.get('parameter.prior.error.filename')
            parameters_dapri = loadtxt(dapri_file_name, dtype=float64)
            self.n_param = len(parameters_dapri)
        else:
            parameters_dapri = None

        if self.Optim_dict['iniconc']:
            logger.error('Not implemented yet, quit program in read_apri_std, class Precon...')
            sys.exit(0)
        else:
            iniconc_dapri = None

        return self.struct2state(self.emis_dapri,iniconc_dapri,parameters_dapri)

    def calc_temp_corr(self, corlen, dt, n):
        A = zeros((n,n))
        P = zeros((n,n))
        D = zeros((n,n))
        lam = zeros(n)
        if corlen<1.e-20:
            A = eye(n)
            P = eye(n)
            D = eye(n)
            lam = ones(n)
        else:
            dummy_X, dummy_Y = meshgrid(arange(n), arange(n))
            A = exp(-abs(dummy_X-dummy_Y)*dt/corlen)
            P,D = self.matrix_square_root(A)
        return P, D

    def matrix_square_root(self, B):
        # Given a real symmetric matrix B, calculate L such that LL^T = B
        # Actually, calculate L as the product of a unitary matrix and a diagonal matrix
        lam, P = linalg.eigh(B)
        # Sort eigenvalues and switch columns of P accordingly
        sort_order = argsort(lam)
        lam = lam[sort_order]
        lam[lam<0.0] = 0.0
        P = P[:, sort_order]
        D = diag(sqrt(lam))
        # Make sure that the elements in the top row of P are non-negative
        col_sign = where(P[0]<0.0, -1.0, 1.0)
        P = P*col_sign
        return P, D

    def struct2state(self,emission=None,iniconc=None,parameters=None):
        state_x = []   # list to map the state
        if emission is not None:
            for tracer in self.species:
                # now only include the info for categories marked for opimization:
                cat_Bh_file = []
                cat_n_hor = []
                cat_vec2ll_file = []
                for cat_key in self.Emission[tracer]['cat_list']:
                    cat = cat_key[0]   # name cat, e.g. anthropogenic
                    corr = cat_key[1]  # horizontal correlation key, e.g. 1000.0-g
                    re = self.Emission[tracer]['cat_opt'][cat_key]  # info about regions and errors

                    # get the number of time intervals for this category:
                    region = re[0]['region']
                    time_interval = self.Emission[region][tracer][cat]['time_interval']
                    nt = len(time_interval['time_start'])

                    fname = 'Bh:'+self.rcf.get('my.zoom')+':'
                    for rei in re:
                        fname += rei['region'].strip()+'_'
                    cl = int(float(corr[:corr.index('-')]))
                    ct = corr[corr.index('-')+1:]
                    fname += '%5.5i'%(cl)+'_'+ct.strip()+'.nc'
                    corr_file = os.path.join(self.rcf.get('correlation.inputdir'),fname)
                    cat_Bh_file.append(corr_file)
                    self.hor_corr_files[tracer][cat] = corr_file
                    self.hor_corr_choices[tracer][cat] = ct.strip()

                    # get the vec2ll mapping for this category:

                    vec2ll_file = os.path.join(self.rcf.get('my.zoom.dir'), 'vec2ll_')
                    vec2ll_file += '_'.join([rei['region'].strip() for rei in re]) + '.nc'
                    cat_vec2ll_file.append(vec2ll_file)

                    if os.path.isfile(vec2ll_file):
                        with Dataset(vec2ll_file, 'r') as f:
                            f = Dataset(vec2ll_file,'r')
                            n_hor = len(f.dimensions['n_hor'])
                            cat_n_hor.append(n_hor)
                            # note: switch here to python count....
                            vec2ll_i = f.variables['vec2ll_i'][:] - 1
                            vec2ll_j = f.variables['vec2ll_j'][:] - 1
                            vec2ll_region = f.variables['vec2ll_region'][:] - 1
                            # note: switch here to python count....
                            vec2ll_lon = f.variables['vec2ll_lon'][:]
                            vec2ll_lat = f.variables['vec2ll_lat'][:]
                            vec2ll_enter_region = f.variables['vec2ll_enter_region'][:]
                            vec2ll_leave_region = f.variables['vec2ll_leave_region'][:]
                    else:
                        sys.stderr.write('vec2ll file not yet created. Done in forward run => do this first\n')
                        sys.exit(2)
                    # loop over time
                    for itime, time in enumerate(time_interval['time_start']):
                        for i_hor in range(n_hor):
                            if vec2ll_enter_region[i_hor] == 1:
                                iregion = vec2ll_region[i_hor]
                                region = self.region_names[iregion]
                                em = emission[region][tracer][cat]['emission_data']
                            i = vec2ll_i[i_hor]
                            j = vec2ll_j[i_hor]
                            state_x.append(em[itime,j,i])
                # save file names in Emission dictionary:
                self.Emission[tracer]['cat_Bh_file'] = cat_Bh_file
                self.Emission[tracer]['cat_n_hor']   = cat_n_hor
                self.Emission[tracer]['cat_vec2ll_file'] = cat_vec2ll_file

        if iniconc is not None:
           logger.error('Iniconc optimization not implemented')
           sys.exit(3)

        if parameters is not None:
            state_x.extend(parameters)

        return array(state_x)

    def state2struct(self,state, replace_field = None):
        emis = None
        iniconc = None
        parameters = None
        ipos = 0

        if self.Optim_dict['emission']:
            emis = dict.fromkeys(self.region_names)
            # start with the base emissions:
            for region in self.region_names:
                emis[region] = dict.fromkeys(self.species)
                for tracer in self.species:
                    categories = self.Emission[region][tracer]['categories']
                    emis[region][tracer] = dict.fromkeys(categories)
                    for cat in categories:
                        if replace_field is None:
                            emstore = deepcopy(self.Emission[region][tracer][cat]['emission_data'])
                        else:
                            emstore = deepcopy(replace_field[region][tracer][cat]['emission_data'])
                        emis[region][tracer][cat] = {'emission_data': emstore}

            # replace the part that has been optimized:
            for tracer in self.species:
                for c, cat_key in enumerate(self.Emission[tracer]['cat_list']):
                    cat = cat_key[0]
                    re = self.Emission[tracer]['cat_opt'][cat_key]
                    corr_file = self.Emission[tracer]['cat_Bh_file'][c]
                    nt = self.Emission[tracer]['cat_nt'][c]
                    vec2ll_file = self.Emission[tracer]['cat_vec2ll_file'][c]

                    if os.path.isfile(vec2ll_file):
                        with Dataset(vec2ll_file,'r') as f:
                            n_hor = len(f.dimensions['n_hor'])
                            # note: switch here to python count....
                            vec2ll_i = f.variables['vec2ll_i'][:] - 1
                            vec2ll_j = f.variables['vec2ll_j'][:] - 1
                            vec2ll_region = f.variables['vec2ll_region'][:] - 1
                            # note: switch here to python count....
                            vec2ll_lon = f.variables['vec2ll_lon'][:]
                            vec2ll_lat = f.variables['vec2ll_lat'][:]
                            vec2ll_enter_region = f.variables['vec2ll_enter_region'][:]
                            vec2ll_leave_region = f.variables['vec2ll_leave_region'][:]
                    else:
                        logger.error('vec2ll file not yet created. Done in forward run => do this first')
                        sys.exit(2)

                    for it in range(nt):
                        for i_hor in range(n_hor):
                            if vec2ll_enter_region[i_hor] == 1:
                                iregion = vec2ll_region[i_hor]
                                region = self.region_names[iregion]
                            i = vec2ll_i[i_hor]
                            j = vec2ll_j[i_hor]
                            emis[region][tracer][cat]['emission_data'][it,j,i] = state[ipos]
                            ipos += 1

            del emstore

        if self.Optim_dict['iniconc']:
            pass

        if self.Optim_dict['parameters']:
            parameters = zeros(self.n_param)
            for i in range(self.n_param):
                parameters[i] = state[ipos]
                ipos += 1

        return emis, iniconc, parameters