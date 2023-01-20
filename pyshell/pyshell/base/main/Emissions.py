#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True

import logging
logging.basicConfig(datefmt='%d %b %Y %H:%M:%S', format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

#import re, os.path, shutil, tempfile, hashlib, rc
import os.path
from pyshell.tmflex import rc
from dateutil.relativedelta import relativedelta
#import calendar, string, progressbar
#import random as noncon_random
import h5py
import cPickle as pickle
#from time import mktime, sleep
import numpy as np
from pyshell.base.helper.Utilities import *
from datetime import datetime, timedelta
#from gzip import GzipFile
#from tarfile import open as TarFile
#from scipy import optimize
from netCDF4 import Dataset, OrderedDict
from collections import defaultdict
from copy import deepcopy
from pyshell.tm5_utils import redistrib_flux

class del_time(object):

    def __init__(self, del_obj):
        super(del_time, self).__init__()
        # del_obj should be a timedelta object
        self.days = del_obj.days
        self.seconds = del_obj.seconds
        self.microseconds = del_obj.microseconds
        self.resolution = del_obj.resolution
        self.total_seconds = del_obj.total_seconds

    def to_seconds(self):
        # I'm sick of timedelta not implementing a to_seconds() method
	#
	# ... what about:
        return self.total_seconds()

    def __div__(self, denom):
        self_sec = self.total_seconds()
        denom_sec = denom.total_seconds()
        if self_sec%denom_sec == 0.0:
            return int(self_sec/denom_sec)
        else:
            return self_sec/denom_sec

class my_Dataset(Dataset):

    def __init__(self, file_name, *args, **kwargs):

        # We need to check if the mode is 'r' or 'a'. Basically, if it's not 'w', check for file existence.
        check_exist = True
        if 'mode' in kwargs:
            check_exist = not kwargs['mode'].startswith('w')
        if len(args) > 0:
            check_exist = not args[0].startswith('w')
        if check_exist:
            if not os.path.exists(file_name):
                raise RuntimeError('File %s does not exist'%file_name)
        Dataset.__init__(self, file_name, *args, **kwargs)

class TM5_emission(dict):
    """
    TM5 emissions are stored in a dictionary structure within python. This does not allow for simple operations such as
    'set all emissions to zero', 'halve all emissions', 'duplicate emissions', etc. This derived class will contain those
    methods.
    """
    def setZero(self):
        """
        Set all the emissions to zero
        """
        for region in self['regions']:
            for tracer in self[region]['tracers']:
                for category in self[region][tracer]['categories']:
                    self[region][tracer][category]['emission_data'][:] = 0.0

    def setRandom(self):
        """
        Set emissions to random values, of roughly the same magnitude as actual emissions
        """
        for region in self['regions']:
            for tracer in self[region]['tracers']:
                for category in self[region][tracer]['categories']:
                    orig_emis = self[region][tracer][category]['emission_data']
                    mean_emis = orig_emis.mean()
                    std_emis = orig_emis.std()
                    self[region][tracer][category]['emission_data'][:] = mean_emis + std_emis * random.standard_normal(orig_emis.shape)

    def zeros_like(self):
        """
        Return a copy of self with zero emissions
        """
        dummy_struct = TM5_emission.fromkeys(self['regions'])
        for region in self['regions']:
            dummy_struct[region] = deepcopy(self[region])
        dummy_struct['regions'] = deepcopy(self['regions'])
        dummy_struct.setZero()
        return dummy_struct

    def __sub__(self, base_emis):
        # first make a copy
        ret_copy = deepcopy(self)
        for region in self['regions']:
            for tracer in self[region]['tracers']:
                for category in self[region][tracer]['categories']:
                    ret_copy[region][tracer][category]['emission_data'] = \
                        self[region][tracer][category]['emission_data'] - base_emis[region][tracer][category]['emission_data']
        return ret_copy

    def __add__(self, base_emis):
        # first make a copy
        ret_copy = deepcopy(self)
        for region in self['regions']:
            for tracer in self[region]['tracers']:
                for category in self[region][tracer]['categories']:
                    ret_copy[region][tracer][category]['emission_data'] = \
                        self[region][tracer][category]['emission_data'] + base_emis[region][tracer][category]['emission_data']
        return ret_copy

    def __mul__(self, factor):
        """
        Overloaded multiplication operator, to multiply all emissions by some factor
        """
        for region in self['regions']:
            for tracer in self[region]['tracers']:
                for category in self[region][tracer]['categories']:
                    self[region][tracer][category]['emission_data'][:] = factor * self[region][tracer][category]['emission_data'][:]

    def __div__(self, divisor):
        """
        Overloaded division operator, to multiply all emissions by some factor
        """
        for region in self['regions']:
            for tracer in self[region]['tracers']:
                for category in self[region][tracer]['categories']:
                    self[region][tracer][category]['emission_data'][:] = self[region][tracer][category]['emission_data'][:] / divisor

    def save(self, file_name, destruct=True):
        """
        Save self in a bunch of HDF5 and state files. The file_name refers to the netcdf file, add '.state' to it for
        the pickle file. Save is by default destructive, as in 'self' is changed.
        """
        pickle_file = file_name + '.state'

        with h5py.File(file_name, 'w') as h5file, open(pickle_file, 'wb') as pfid:

            for region in self['regions']:
                tracer_list = self[region]['tracers']
                r_gid = h5file.create_group(region)
                for tracer in tracer_list:
                    t_gid = r_gid.create_group(tracer)
                    for category in self[region][tracer]['categories']:
                        c_gid = t_gid.create_group(category)
                        _ = c_gid.create_dataset('emission_data', data=self[region][tracer][category]['emission_data'])
                        if destruct:
                            del self[region][tracer][category]['emission_data']

            for tracer in tracer_list:
                # Also store the cat_Hor_L and cat_Temp_L arrays if present
                if 'cat_Hor_L' in self[tracer] and 'cat_Temp_L' in self[tracer]:
                    t_gid = h5file.create_group(tracer)
                    num_cats = len(self[tracer]['cat_Hor_L'])

                    for i_cat in range(num_cats):
                        _ = t_gid.create_dataset('hor_L_%02i'%i_cat, data=self[tracer]['cat_Hor_L'][i_cat])
                        _ = t_gid.create_dataset('temp_L_%02i'%i_cat, data=self[tracer]['cat_Temp_L'][i_cat])

                    if destruct:
                        del self[tracer]['cat_Hor_L']
                        del self[tracer]['cat_Temp_L']

            # now store the rest in the pickle file
            for k, v in self.items():
                pickle.dump((k,v), pfid, pickle.HIGHEST_PROTOCOL)

    def load(self, file_name, cleanup=True):
        """
        Load the data saved by self.save
        """
        pickle_file = file_name + '.state'

        if not os.path.exists(file_name):
            logging.error('%s not found'%file_name)
            sys.exit(1)

        if not os.path.exists(pickle_file):
            logging.error('%s not found'%pickle_file)
            sys.exit(1)

        # First, load the pickled items
        with open(pickle_file, 'rb') as fid:
            try:
                while True:
                    k,v = pickle.load(fid)
                    self[k] = v
            except EOFError:
                pass
        # Delete the state file
        if cleanup:
            os.remove(pickle_file)

        # Now load all the large matrices
        with h5py.File(file_name, 'r') as h5file:
            for region in self['regions']:
                tracer_list = self[region]['tracers']
                for tracer in tracer_list:
                    for category in self[region][tracer]['categories']:
                        self[region][tracer][category]['emission_data'] = h5file[region][tracer][category]['emission_data'][:]

            for tracer in tracer_list:
                num_cats = len(self[tracer]['optim_cat_names'])
                if num_cats > 0:
                    self[tracer]['cat_Hor_L'] = []
                    self[tracer]['cat_Temp_L'] = []

                    for i_cat in range(num_cats):
                        self[tracer]['cat_Hor_L'].append(h5file[tracer]['hor_L_%02i'%i_cat][:])
                        self[tracer]['cat_Temp_L'].append(h5file[tracer]['temp_L_%02i'%i_cat][:])

        # Delete the data file
        if cleanup:
            os.remove(file_name)

class Emissions(object):
    """
    The main purpose of this class is to define an emission structure and some general routines for coarsening/refining emissions.
    Emission data for TM5 is stored in Emissions.Emission, which is a dictionary with the following structure:

    {
    'region1' :
        {
        'tracer1' :
            {
            'category1' :
                {
                'emission_data'   : n_time x n_lat x n_lon floating point array with tracer emission in kg/sec/gridbox
                'time_resolution' : a string of the form 'monthly', 'daily', 'daily+2', etc. specifying the time resolution
                'time_interval'   :
                    {
                    'time_start': [t0, t1, t2, ... t(n_time-1)] -- list of datetime objects specifying the starting times for each of the n_time elements of 'emission_data'
                    'time_end'  : [t1, t2, t3, ... t(n_time)] -- list of datetime objects specifying the ending times for each of the n_time elements of 'emission_data'
                    'time_mid'  : [(t0+t1)/2, (t1+t2)/2, ... ] -- list of datetime objects specifying the midpoints for each of the n_time elements of 'emission_data'
                    'dt'        : (ti-t(i-1)) -- single timedelta (or relativedelta for monthly emissions) object specifying the size of the time step
                    }
                'error'           : prior error as percentage of emission at the grid scale, e.g., 100.0 (error as large as absolute emission) or 250.0
                'corr'            : spatial correlation information specified as a string, such as '1000.0-g' or '500.0-e'
                'tcorr'           : temporal correlation information specified as a string, such as '9.50-e-monthly' or '1.00-e-daily+2'
                'optimize'        : 0 or 1, to specify whether this category within this region is to be optimized or not
                }
            'category2' :
            }
        'tracer2' :
            {
            'category1' :
                {
                ...
                }
            'category2' :
                {
                ...
                }
            }
        }
    'region2' :
        {
        'tracer1' :
            {
            'category_21' :
                {
                ...
                }
            'category_22' :
            ...
            }
        'tracer2' :
            {
            ...
            }
        }
    ...
    }

    Any species-specific derived class must have one routine called LoppThroughPeriods(), which fills this structure.
    Of all the elements in this structure, LoopThroughPeriods only needs to fill in 'emission_data'. All the others
    are automatically filled in by Emission.__init__().
    """
    def __init__(self, StartTime, EndTime, subdir_tag='', **kwargs):
        """
        The times are given as tuples (year,month,day,hour) with hour optional. The subdir_tag is meant to create a
        subdir within which the emission file will be placed. For example, if normally the emission file would
        have been placed in

        /scratch/shared/sbasu/var4d/CO2/glb/ml60/tropo25/carbonsat_unbiased/output/2008100100-2008101100/emission.nc4

        then with a subdir_tag of 'no_slopes', it ends up in

        /scratch/shared/sbasu/var4d/CO2/glb/ml60/tropo25/carbonsat_unbiased/output/2008100100-2008101100/no_slopes/emission.nc4

        This allows multiple forward/adjoint runs over the same period and using the same input data, but perhaps with
        subtle changes in parameters.
        """
        super(Emissions, self).__init__()
        self.StartDate = datetime(*StartTime)
        self.EndDate = datetime(*EndTime)
        self.subdir_tag = subdir_tag
        self.sec_day = 86400.
        self.rcf = rc.RcFile(os.environ['pyshell.rc'])
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.GetZoomRegions()
        # Monte Carlo estimation of posterior covariance or not?
        self.optim_mc = self.rcf.get('my.optimizer.class') in ['conGrad_MC', 'm1qn3_MC']

        emission_file_name = self.rcf.get('PyShell.em.filename')
        self.emission_file_name = self.putDateString(emission_file_name, True)
        self.createDifferentialArea()
        # Define a destination horizontal grid for later
        self.lat_grid = {}
        self.lon_grid = {}
        for region_name, lat_data in zip(self.zoom_regions_names, self.zoom_regions_lat):
            ybeg, yend, jm = lat_data
            self.lat_grid[region_name] = np.linspace(ybeg, yend, jm+1)
        for region_name, lon_data in zip(self.zoom_regions_names, self.zoom_regions_lon):
            xbeg, xend, im = lon_data
            self.lon_grid[region_name] = np.linspace(xbeg, xend, im+1)

    def create_emission_structure( self):
        # create dictionary to contain emission info to be read from rc file:
        self.Emission = TM5_emission.fromkeys(self.zoom_regions_names)
        self.Emission['regions'] = self.zoom_regions_names
        for region in self.zoom_regions_names:
            self.Emission[region] = dict.fromkeys(self.species)
            self.Emission[region]['tracers'] = self.species
        for tracer in self.species:
            cat_opt = defaultdict(list)
            cat_list  = []
            optim_cat_names = []
            self.Emission[tracer] = {}
            #self.Emission[tracer]['emi_class'] = self.rcf.get(tracer+'.emission.class')
            for region, xlims, ylims in zip(self.zoom_regions_names, self.zoom_regions_lon, self.zoom_regions_lat):
                self.Emission[region][tracer] = {}
                ncat = self.rcf.get('emission.'+tracer+'.'+region+'.categories','int')
                categories = []
                for icat in np.arange(ncat):
                    infoline = self.rcf.get('emission.'+tracer+'.'+region+'.category%1i'%(icat+1))
                    categ = infoline[:infoline.find(';')].strip()
                    categories.append(categ)
                    error = np.float64(infoline.split(';')[1])
                    corr  = infoline.split(';')[2]
                    tcorr = infoline.split(';')[3]
                    time_res = infoline.split(';')[3]
                    res_key = time_res.split('-')[2]
                    optimize = int(infoline.split(';')[4])
                    remarks = infoline.split(';')[5].strip()
                    time_interval = self.split_time_interval( self.StartDate, self.EndDate, res_key)
                    nt = len(time_interval['time_start'])
                    self.Emission[region][tracer][categ] = {'emission_data'  : np.zeros((nt,int(ylims[2]),int(xlims[2])), np.float64), \
                        'time_resolution': res_key,      \
                        'time_interval'  : time_interval, \
                        'error'          : error, \
                        'corr'           : corr, \
                        'tcorr'          : tcorr, \
                        'optimize'       : optimize, \
                        'remarks'        : remarks}
                    # sometimes we want to specify the errors in emissions separately, instead of as a single fraction of the emission
                    sep_err = self.rcf.get('%s.%s.sep_error'%(tracer, categ), 'bool', default=False)
                    if sep_err:
                        self.Emission[region][tracer][categ]['emission_error'] = np.zeros((nt,int(ylims[2]),int(xlims[2])), np.float64)
                    # create library with tuple as key:
                    if optimize == 1:
                        tup = (categ,corr,tcorr)
                        cat_opt[tup].append({'region':region, 'error':error})
                        if tup not in cat_list:
                            cat_list.append(tup)
                        if categ not in optim_cat_names:
                            optim_cat_names.append(categ)
                self.Emission[region][tracer]['categories'] = categories
            self.Emission[tracer]['cat_list'] = cat_list
            self.Emission[tracer]['cat_opt'] = cat_opt
            self.Emission[tracer]['optim_cat_names'] = optim_cat_names

    def putDateString(self, input_string, file=True):
        output_string = input_string
        for atom in ['Y', 'y', 'm', 'M', 'd', 'H', 'b', 'B', 'j', 'p', 'S']:
            for i, t in zip([1,2], [self.StartDate, self.EndDate]):
                rep_string = '<%1s%1i>'%(atom, i)
                form_string = '%%%1s'%atom
                output_string = output_string.replace(rep_string, t.strftime(form_string))
        # Now add the subdir_tag
        if file:
            dirname, filename = os.path.split(output_string)
            output_string = os.path.join(dirname, self.subdir_tag, filename)
        else:
            output_string = os.path.join(output_string, self.subdir_tag)
        return output_string

    def split_time_interval( self, StartDate, EndDate, res_key):
        if res_key.find('+') == -1:
            nres = 1
        else:
            nres = int(res_key.split('+')[1])
        if res_key.find('monthly') != -1:
            delta = relativedelta(months=+nres)
            if nres == 1:
                delta2 = timedelta(days=15)
            else:
                print 'optimization for periods longer than month not yet implemented'
                sys.exit()
            # set to start of month: user should start optimization at 1st day month, and end last day month
            # SB: Try to start from the correct start date and not the first of the month
            #TheStartDate = datetime(StartDate.year,StartDate.month,1,0,0,0)
        elif res_key.find('daily') != -1:
            #TheStartDate = datetime(StartDate.year,StartDate.month,StartDate.day,0,0,0)
            delta = timedelta(days=+nres)
            delta2 = delta/2
        else:
            print "invalid time resolution:", res_key
            sys.exit()
        #d = TheStartDate
        d = datetime(StartDate.year,StartDate.month,StartDate.day,0,0,0)
        time_start = []
        time_end = []
        time_mid = []
        while d < EndDate:
            time_start.append(d)
            # what is the ending time for this window?
            if res_key.find('monthly') != -1:
                d_end = self.nextMonth(d)
            elif res_key.find('daily') != -1:
                d_end = d + delta

            # Need to check here whether d+delta exceeds the end date or not
            if d_end <= EndDate:
                time_end.append(d_end)
                time_mid.append(d + (d_end-d)/2)
            else:
                time_end.append(EndDate)
                time_mid.append(d + (EndDate - d)/2)

            d = d_end
        time_interval = {'time_start': time_start, 'time_end': time_end, 'time_mid': time_mid, 'dt': delta}
        return time_interval

    def nextMonth(self, cur_date):
        # return the beginning of next month
        year, month = cur_date.timetuple()[:2]
        if month == 12:
            next_year = year + 1
            next_month = 1
        else:
            next_year = year
            next_month = month + 1
        next_date = datetime(next_year, next_month, 1, 0, 0, 0)
        return next_date

    def makeOverlap(self, ts, te, period):
        """
        Given a time period (ts,te), calculate the fractional coverage per period, which can be 'year' or 'month'
        """
        overlap_dict = OrderedDict()
        t1 = ts
        if period == 'month':
            dt = relativedelta(months=1)
            t2 = datetime(t1.year, t1.month, 1) + dt
        elif period == 'year':
            dt = relativedelta(years=1)
            t2 = datetime(t1.year, 1, 1) + dt
        t2 = min(t2, te)

        year, month = t1.year, t1.month
        num_secs = del_time(t2-t1).to_seconds()
        if period == 'month':
            overlap_dict[(year, month)] = num_secs
        elif period == 'year':
            overlap_dict[year] = num_secs

        while t2 < te:
            t1 = t2
            t2 = t2 + dt
            t2 = min(t2, te)
            year, month = t1.year, t1.month
            num_secs = del_time(t2-t1).to_seconds()
            if period == 'month':
                overlap_dict[(year, month)] = num_secs
            elif period == 'year':
                overlap_dict[year] = num_secs

        total_seconds = 0.0
        for v in overlap_dict.values():
            total_seconds += v

        for k,v in overlap_dict.items():
            overlap_dict[k] = v/total_seconds

        return overlap_dict

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
        self.zoom_info = defaultdict(dict)
        for region in self.zoom_regions_names:
            xbeg.append(self.rcf.get('region.'+region+'.xbeg','float'))
            self.zoom_info[region]['xbeg'] = self.rcf.get('region.'+region+'.xbeg','float')
            ybeg.append(self.rcf.get('region.'+region+'.ybeg','float'))
            self.zoom_info[region]['ybeg'] = self.rcf.get('region.'+region+'.ybeg','float')
            xend.append(self.rcf.get('region.'+region+'.xend','float'))
            self.zoom_info[region]['xend'] = self.rcf.get('region.'+region+'.xend','float')
            yend.append(self.rcf.get('region.'+region+'.yend','float'))
            self.zoom_info[region]['yend'] = self.rcf.get('region.'+region+'.yend','float')
            im.append(self.rcf.get('region.'+region+'.im','int'))
            self.zoom_info[region]['im'] = self.rcf.get('region.'+region+'.im','int')
            jm.append(self.rcf.get('region.'+region+'.jm','int'))
            self.zoom_info[region]['jm'] = self.rcf.get('region.'+region+'.jm','int')
        self.zoom_regions_lon = zip(xbeg,xend,im)
        self.zoom_regions_lat = zip(ybeg,yend,jm)

    def multiplyEmission(self, emis, region, cat_name):
        """
        Sometimes, I need to artificially multiply the flux of one specific category over a specific region. An example
        is multiplying the terrestrial disequilibrium flux by 3 over North America for a real data radiocarbon inversion
        experiment.
        """
        # get the scales and masks
        scale_dict = self.readEmisScales(cat_name)
        # To start off, we create an emission identical to the input emission
        base_emis = np.zeros_like(emis)
        base_emis[:] = emis
        for i_scale, scale_info in scale_dict.items():
            # For each mask, we need to coarsen it to a rectangle fitting this zoom region.
            coarse_mask = self.coarsenMask(np.float64(scale_info['mask']), region)
            # Now coarse_mask is a mask of zeros and ones. Where it is zero, we want to keep the original emissions (or
            # the emissions from the previous iteration). Where it is one, we want to replace it completely by
            # scale_info['factor'] * emis, where emis is the original emission. Where it is in between, we want to linearly
            # interpolate between the base_emis and scale_info['factor'] * emis.
            base_emis = (1.0-coarse_mask) * base_emis + coarse_mask * scale_info['factor'] * emis
        return base_emis

    def multiplyError(self, emis, region, cat_name):
        """
        In my standard inversions, the error on the terrestrial disequilibrium flux is set to be a constant
        fraction of the gross flux itself. However, Scott and John think that we should allow more freedom
        to the flux from the boreal North American region in order to fit the bump at WLEF. Quite frankly,
        to me it seems equally likely that the bump in summer is from southerly winds, in which case the
        cause would be enhanced respiration in the SE US. However, it's still worthwhile to be able to specify
        different prior error percentages over different geographical regions. This routine multiplies the
        gross flux by a heterogeneous mask to give a heterogeneous flux error.

        This routine takes a (net) flux, a zoom region name, and a category name, and returns an array of
        flux errors.
        """
        # First, get the gross flux over this region
        gross_flux = abs(emis)
        # What is the standard percentage to be applied globally?
        cats_opt = self.Emission[self.tracer]['cat_opt'].keys()
        # We need to select the key corresponding to correct category
        found_cat = False
        for cat_opt in cats_opt:
            if cat_opt[0] == cat_name:
                found_cat = True
                break
        if not found_cat:
            raise RuntimeError("You want to scale errors for %s/%s, but looks like that category is not being optimized"%(self.tracer, cat_name))
        # Element 0 is selected because this is a list with one element for each zoom region. However, since the error fraction
        # does not change with zoom region, it is sufficient to read that for just one region.
        err_dict = self.Emission[self.tracer]['cat_opt'][cat_opt][0]
        err_frac = err_dict['error']/100.0
        scale_dict = self.readErrorScales(cat_name)
        # To start off, we create an error with the scaling factor err_frac everywhere
        base_error = err_frac * gross_flux
        for i_scale, scale_info in scale_dict.items():
            # For each mask, we need to coarsen it to a rectangle fitting this zoom region.
            coarse_mask = self.coarsenMask(np.float64(scale_info['mask']), region)
            # Now coarse_mask is a mask of zeros and ones. Where it is zero, we want to keep the error from the previous
            # iteration, i.e., the previous base_error. Where it is one, we want to replace it completely by
            # scale_info['factor'] * gross_flux. Where it is in between, we want to linearly interpolate between the
            # previous base_error and scale_info['factor'] * gross_flux.
            base_error = (1.0 - coarse_mask)*base_error + coarse_mask * scale_info['factor'] * gross_flux
        return base_error

    def coarsenMask(self, msk, region):
        """
        Sometimes we have a mask (msk) of zeros and ones, indicating which pixels belong to a particular region (say
        Canada). This mask may be on a 1x1 or finer grid. This routine takes such a mask and returns a coarsened
        mask for a specific zoom region.
        """
        nlat, nlon = msk.shape
        ip_lats = np.linspace(-90.,90.,nlat+1)
        ip_lons = np.linspace(-180.,180.,nlon+1)
        op_lats = self.lat_grid[region]
        op_lons = self.lon_grid[region]
        # We have a mask of ones and zeros that needs to be converted into a new mask of ones and zeros after coarsening.
        # The operation is identical to coarsening an intensive flux field to another intensive flux field.
        return_arr = redistrib_flux.regrid_flux_submul(ip_lats,ip_lons,op_lats,op_lons,msk,True)
        # Make sure that return_arr is between 0 and 1
        if return_arr.min() < -1.0E-10 or return_arr.max() > 1.0+1.0E-10:
            sys.stderr.write("After coarsening mask, elements are larger than 1 or smaller than 0\n")
            sys.stderr.write("Maximum of mask = 1 + %20.12e\n"%(return_arr.max()-1))
            sys.stderr.write("Minimum of mask = %20.12e\n"%return_arr.min())
            raise RuntimeError("Mask must be between 0 and 1")
        # If it's roughly right, make it exactly right
        return_arr = np.where(return_arr > 1.0, 1.0, return_arr)
        return_arr = np.where(return_arr < 0.0, 0.0, return_arr)
        return return_arr

    def readEmisScales(self, cat_name):
        # For a given tracer/category, read the (optional) scales to apply
        key = 'emission.%s.%s.num_scales'%(self.tracer, cat_name)
        num_scales = self.rcf.get(key, 'int', default=0)
        scale_dict = OrderedDict.fromkeys(range(num_scales))
        for i_scale in range(num_scales):
            file_name = self.rcf.get('emission.%s.%s.%02i.maskfile'%(self.tracer, cat_name, i_scale+1))
            region    = self.rcf.get('emission.%s.%s.%02i.region'%(self.tracer, cat_name, i_scale+1))
            varname   = self.rcf.get('emission.%s.%s.%02i.varname'%(self.tracer, cat_name, i_scale+1))
            reg_list  = self.rcf.get('emission.%s.%s.%02i.regionlist'%(self.tracer, cat_name, i_scale+1))
            scale_fac = self.rcf.get('emission.%s.%s.%02i.scale'%(self.tracer, cat_name, i_scale+1), 'float') # not a percentage
            with my_Dataset(file_name, 'r') as fid:
                reg_names = [s.strip() for s in getattr(fid, reg_list).split(',')]
                reg_idx = reg_names.index(region)
                msk = fid.variables[varname][reg_idx, :, :]
                scale_dict[i_scale] = {'mask': msk[:], 'factor': scale_fac}
        # We do not check for overlapping masks, because it is in fact an advantage not to check. We could, for example,
        # assign a certain factor to N America and a different factor to a subregion N American temperate, simply by specifying
        # the N American temperate factor after the N American factor in the rc file.
        return scale_dict

    def readErrorScales(self, cat_name):
        """
        For a given category, read the scaling and masks for the error specification for a category
        """
        key = 'emission.error.%s.%s.num_scales'%(self.tracer, cat_name)
        num_scales = self.rcf.get(key, 'int')
        scale_dict = OrderedDict.fromkeys(range(num_scales))
        for i_scale in range(num_scales):
            file_name = self.rcf.get('emission.error.%s.%s.%02i.maskfile'%(self.tracer, cat_name, i_scale+1))
            region    = self.rcf.get('emission.error.%s.%s.%02i.region'%(self.tracer, cat_name, i_scale+1))
            varname   = self.rcf.get('emission.error.%s.%s.%02i.varname'%(self.tracer, cat_name, i_scale+1))
            reg_list  = self.rcf.get('emission.error.%s.%s.%02i.regionlist'%(self.tracer, cat_name, i_scale+1))
            scale_fac = self.rcf.get('emission.error.%s.%s.%02i.scale'%(self.tracer, cat_name, i_scale+1), 'float') # not a percentage
            with my_Dataset(file_name, 'r') as fid:
                reg_names = [s.strip() for s in getattr(fid, reg_list).split(',')]
                reg_idx = reg_names.index(region)
                msk = fid.variables[varname][reg_idx, :, :]
                scale_dict[i_scale] = {'mask': msk[:], 'factor': scale_fac}
        # We do not check for overlapping masks, because it is in fact an advantage not to check. We could, for example,
        # assign a certain factor to N America and a different factor to a subregion N American temperate, simply by specifying
        # the N American temperate factor after the N American factor in the rc file.
        return scale_dict

    def diffArea(self, lat_min, lat_max, lats, lon_min, lon_max, lons):
        EarthRad = 6.371e6 # meters
        dLon = (np.pi/180.) * ((lon_max-lon_min)/lons)
        dS = np.zeros((lats+1, lons), np.float64)
        Lats = (np.pi/180.) * np.linspace(lat_min, lat_max, lats+1)
        for i, lat in enumerate(Lats):
            dS[i] = EarthRad * EarthRad * dLon * np.sin(lat)
        dS = np.diff(dS, axis=0)
        return dS

    def createDifferentialArea(self):
        """
        Create a bunch of arrays containing the areas of regular gridboxes in square meters
        """
        self.dS_m2 = {}
        for lats, lons in zip([45,90,180,360,1800], [60,120,360,720,3600]):
            self.dS_m2[(lats,lons)] = self.diffArea(-90., 90., lats, -180., 180., lons)

    def CoarsenedEmis(self, in_array, l_x, l_y):
        """
        in_array has emissions, always in a 1x1 grid, because that's the resolution of the priors
        """
        lump_x = in_array.shape[0]/l_x
        lump_y = in_array.shape[1]/l_y
        dummyVar = np.zeros((l_x, l_y), in_array.dtype)
        for i in range(l_x):
            for j in range(l_y):
                dummyVar[i, j] = np.sum(in_array[i*lump_x:(i+1)*lump_x, j*lump_y:(j+1)*lump_y])
        return dummyVar

    def SpreadEmissionsToFineGrid(self, in_array, dy_min, dx_min):
        dx_orig = 360./in_array.shape[1] # longitude bin of the input array
        dy_orig = 180./in_array.shape[0] # latitude bin of the input array
        if fmod(dx_orig,dx_min) != 0.0 or fmod(dy_orig,dy_min) != 0.0 or fmod(180.,dy_orig) != 0.0 or fmod(360.,dx_orig) != 0.0:
            raise RuntimeError('Coarse array cannot be redistributed')
        out_lats = int(round(180./dy_min))
        out_lons = int(round(360./dx_min))
        out_array = np.zeros((out_lats,out_lons), in_array.dtype)
        lump_x = out_array.shape[1]/in_array.shape[1]
        lump_y = out_array.shape[0]/in_array.shape[0]
        for i in range(in_array.shape[0]):
            for j in range(in_array.shape[1]):
                out_array[i*lump_y:(i+1)*lump_y,j*lump_x:(j+1)*lump_x] = in_array[i,j]/(lump_x*lump_y)
        return out_array

    def CoarsenEmissionsByZoom(self, in_array):
        dx_orig = 360./in_array.shape[1] # longitude bin of the input array
        dy_orig = 180./in_array.shape[0] # latitude bin of the input array
        # We assume that the input array covers the entire globe, but not that it's on a 1x1 grid
        # For this routine to work, we need the input array to be on at least as fine a resolution
        # as the finest zoom region. Check if that is the case:
        dx_regs = []
        dy_regs = []
        for lon, lat in zip(self.zoom_regions_lon, self.zoom_regions_lat):
           dx_regs.append((lon[1] - lon[0])/lon[2])
           dy_regs.append((lat[1] - lat[0])/lat[2])
        if dx_orig > np.array(dx_regs).min() or dy_orig > np.array(dy_regs).min():
           # That is not the case, so we need to redistribute the input array
           in_array = self.SpreadEmissionsToFineGrid(in_array, np.array(dy_regs).min(), np.array(dx_regs).min())
           dx_orig = 360./in_array.shape[1]
           dy_orig = 180./in_array.shape[0]
        return_dict = dict.fromkeys(self.zoom_regions_names)
        for region, xlims, ylims in zip(self.zoom_regions_names, self.zoom_regions_lon, self.zoom_regions_lat):
           # region will be a name, such as nam1x1
           # xlims will be the starting point, ending point and divisions of longitude, such as [-132,-66,66]
           # ylims will be the starting point, ending point and divisions of latitude, such as [20,64,44]
           xs_id = int((xlims[0]+180.0)/dx_orig)
           ys_id = int((ylims[0]+90.0)/dy_orig)
           xe_id = xs_id + int((xlims[1]-xlims[0])/dx_orig)
           ye_id = ys_id + int((ylims[1]-ylims[0])/dy_orig)
           relevant_in_array = in_array[ys_id:ye_id, xs_id:xe_id] # lat x lon
           out_array = np.zeros((int(ylims[2]), int(xlims[2])), in_array.dtype)
           lump_x = relevant_in_array.shape[1]/int(xlims[2]) # longitude
           lump_y = relevant_in_array.shape[0]/int(ylims[2]) # latitude
           for i in range(int(ylims[2])): # latitude
               for j in range(int(xlims[2])): # longitude
                   out_array[i, j] = (relevant_in_array[i*lump_y:(i+1)*lump_y, j*lump_x:(j+1)*lump_x]).sum()
           return_dict[region] = out_array
           del out_array
        return return_dict

    def readOptimFromFile(self, file_name, add_perturb=False):
        """
        Sometimes we want to read the emissions from a single file instead of assembling it. Two most common
        examples of this need are (a) doing a forward run with optimized emissions, (b) doing an inversion with
        an unbiased prior, where prior = posterior + Gaussian noise. It assumes that self.Emission has already
        been initialized, and only reads in the emission data, and the emission errors where relevant.
        """
        fid = my_Dataset(file_name, 'r')
        for region in self.zoom_regions_names:
            rgid = fid.groups[region]
            for tracer in self.species:
                tgid = rgid.groups[tracer]
                categories = self.Emission[region][tracer]['categories']
                for cat in categories:
                    cgid = tgid.groups[cat]
                    self.Emission[region][tracer][cat]['emission_data'] = cgid.variables['emission'][:]
                    if 'emission_error' in cgid.variables:
                        self.Emission[region][tracer][cat]['emission_error'] = cgid.variables['emission_error'][:]
                    else:
                        error = self.Emission[region][tracer][cat]['error']
                        self.Emission[region][tracer][cat]['emission_error'] = 0.01 * error * \
                            abs(self.Emission[region][tracer][cat]['emission_data'])
                    if add_perturb:
                        perturb = random.standard_normal(self.Emission[region][tracer][cat]['emission_data'].shape)
                        self.Emission[region][tracer][cat]['emission_data'] += \
                            (perturb * self.Emission[region][tracer][cat]['emission_error'])
        fid.close()

    def WriteEmissions(self):
        outFile = self.emission_file_name
        checkDir(outFile)
        if os.path.exists(outFile):
            os.remove(outFile)
        fid = my_Dataset(outFile, 'w')
        fid.createDimension('itime',6)
        for region, xlims, ylims in zip(self.zoom_regions_names, self.zoom_regions_lon, self.zoom_regions_lat):
            group = fid.createGroup(region)
            group.createDimension('latitude', ylims[2])
            group.createDimension('longitude', xlims[2])
            group.latmin = ylims[0]
            group.latmax = ylims[1]
            group.lonmin = xlims[0]
            group.lonmax = xlims[1]
            for tracer in self.species:
                tgroup = group.createGroup(tracer)
                categories = self.Emission[region][tracer]['categories']
                for cat in categories:
                    cgroup = tgroup.createGroup(cat)
                    cgroup.createDimension('nt', self.Emission[region][tracer][cat]['emission_data'].shape[0])
                    cgroup.time_resolution =  self.Emission[region][tracer][cat]['time_resolution']
                    cgroup.optimize =  np.int32(self.Emission[region][tracer][cat]['optimize'])
                    # This is Maarten's modification for CO biomass burning, applicable to other tracers as well
                    if self.Emission[tracer]['tf_bb_diurnal'] != None:
                        cgroup.createDimension('hourly',24)
                        var = cgroup.createVariable('tf_diurnal', np.float64, ('hourly',))
                        var[:] = self.Emission[tracer]['tf_bb_diurnal']

                    var = cgroup.createVariable('time_start', np.int16, ('nt', 'itime'))
                    var[:] = np.array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_start']], np.int16)

                    var = cgroup.createVariable('time_end', np.int16, ('nt', 'itime'))
                    var[:] = np.array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_end']], np.int16)

                    var = cgroup.createVariable('time_mid', np.int16, ('nt', 'itime'))
                    var[:] = np.array([d.timetuple()[:6] for d in self.Emission[region][tracer][cat]['time_interval']['time_mid']], np.int16)

                    var = cgroup.createVariable('emission', self.Emission[region][tracer][cat]['emission_data'].dtype, ('nt', 'latitude', 'longitude'))
                    var[:] = self.Emission[region][tracer][cat]['emission_data']

                    write_err = self.rcf.get('%s.%s.sep_error'%(tracer, cat), 'bool', default=False)
                    if write_err:
                        var = cgroup.createVariable('emission_error', self.Emission[region][tracer][cat]['emission_error'].dtype, ('nt', 'latitude', 'longitude'))
                        var[:] = self.Emission[region][tracer][cat]['emission_error']

                    # Write the total emissions per category per zoom region per timestep, for possible debugging
                    # The emissions are in Kg tracer/grid cell/sec
                    nt = self.Emission[region][tracer][cat]['emission_data'].shape[0]
                    emis_totals = np.zeros(nt, np.float64)
                    time_intervals = np.zeros(nt, np.float64)
                    for i in range(nt):
                        t1 = self.Emission[region][tracer][cat]['time_interval']['time_start'][i]
                        t2 = self.Emission[region][tracer][cat]['time_interval']['time_end'][i]
                        dt = (t2-t1).total_seconds()
                        emis_totals[i] = self.Emission[region][tracer][cat]['emission_data'][i].sum() * dt
                        time_intervals[i] = dt

                    var = cgroup.createVariable('emission_total', emis_totals.dtype, ('nt',))
                    var[:] = emis_totals
                    var.comment = 'Total emission for region %s and category %s'%(region, cat)
                    var.unit = 'Kg %s/time step'%tracer

                    var = cgroup.createVariable('time_step_length', time_intervals.dtype, ('nt',))
                    var[:] = time_intervals
                    var.comment = 'Length of time period (seconds)'

                    del emis_totals, time_intervals

        fid.close()

        # print emission totals for debugging
        self.print_emission_totals()

    def print_emission_totals(self):
        # For debugging purposes, print the emission totals through the entire simulation period
        for tracer in self.species:
            print
            print '******** Emission totals for %s ********'%tracer
            for region in self.zoom_regions_names:
                print '    ---- region %s ----'%region
                for category in self.Emission[region][tracer]['categories']:
                    nt = self.Emission[region][tracer][category]['emission_data'].shape[0]
                    emis_total = 0.0
                    for i in range(nt):
                        t1 = self.Emission[region][tracer][category]['time_interval']['time_start'][i]
                        t2 = self.Emission[region][tracer][category]['time_interval']['time_end'][i]
                        dt = del_time(t2-t1).to_seconds()
                        # emission data are in Kg tracer/grid cell/second
                        emis = self.Emission[region][tracer][category]['emission_data'][i]
                        emis_total = emis_total + 1.0E-9 * np.sum(emis*dt)
                    print '%25s : %11.4g Tg %8s'%(category.rjust(25), emis_total, tracer.ljust(8))
        print

    def get_class_from_name( self, class_name):
        _temp = __import__('Emissions', fromlist=[class_name])
        try:
            class_from_name = _temp.__dict__[class_name]
            return class_from_name
        except KeyError:
            sys.stderr.write("Class %s not defined in %s.\n"%(class_name,'Emission'))
            sys.exit()
        except:
            sys.stderr.write("Unknown error importing %s\n"%class_name)
            sys.exit()

    def __call__(self, randomize=False, zero=False):
        if randomize and zero:
            raise RuntimeError('You can choose to either randomize emissions or zero them out, but not both')
        if randomize:
            self.Emission.setRandom()
        elif zero:
            self.Emission.setZero()
        self.WriteEmissions()
        return self.Emission

    def readEmission(self, input_file, emis_type=None):
        """
        Sometimes we want to start from emissions that have already been optimized. In that case, we need to know
        which file to read it from (should be some optimized_state*.nc4). We will fill the self.Emission structure
        from that file.
        """
        if emis_type == 'apos':
            var_name = 'poste_emission'
        elif emis_type == 'apri':
            var_name = 'prior_emission'
        else:
            var_name = 'emission'
        fid = my_Dataset(input_file, 'r')
        for region_name in self.zoom_regions_names:
            reg_gid = fid.groups[region_name]
            for tracer in self.species:
                tra_gid = reg_gid.groups[tracer]
                for cat_name in self.Emission[region_name][tracer]['categories']:
                    cat_gid = tra_gid.groups[cat_name]
                    self.Emission[region_name][tracer][cat_name]['emission_data'] = cat_gid.variables[var_name][:]
        fid.close()

    def consolidateDailyCycles(self, overwrite=True):
        """
        Takes the individual dailycycle_temp_* files written out by LoopThroughPeriods and creates one
        daily cycle file per day. In some cases we may not want to overwrite existing daily cycle files,
        such as for MC inversions. In those cases, setting overwrite=False will just delete the temp
        files (which are in the output folder) without touching the daily cycle files in the input folder.
        """
        time_steps_per_day = int(86400.0/self.secs_per_granule)
        for cur_time in self.dailycycle_dates:
            output_filename = os.path.join(self.rcf.get('dailycycle.folder'), cur_time.strftime("%Y/%m"), \
                self.rcf.get('%s.dailycycle.prefix'%self.tracer)+cur_time.strftime("%Y%m%d.nc4"))
            if os.path.exists(output_filename) and (not overwrite):
                # just delete the temp files
                for ireg, region in enumerate(self.zoom_regions_names):
                    categories = self.Emission[region][self.tracer]['categories']
                    for icat, cat in enumerate(categories):
                        input_filename = os.path.join(self.output_dir, 'dailycycle_temp_%s_cat%i_%s.nc'%\
                            (region,icat,cur_time.strftime("%Y%m%d")))
                        os.remove(input_filename)
            else:
                checkDir(output_filename)
                ofid = my_Dataset(output_filename, 'w')
                for ireg, region in enumerate(self.zoom_regions_names):
                    ogid_reg = ofid.createGroup(region)
                    ogid_reg.createDimension('latitude', self.zoom_regions_lat[ireg][2])
                    ogid_reg.createDimension('longitude', self.zoom_regions_lon[ireg][2])
                    categories = self.Emission[region][self.tracer]['categories']
                    for icat, cat in enumerate(categories):
                        ogid_cat = ogid_reg.createGroup(cat)
                        ogid_cat.createDimension('timesteps', time_steps_per_day)
                        input_filename = os.path.join(self.output_dir, 'dailycycle_temp_%s_cat%i_%s.nc'%\
                            (region,icat,cur_time.strftime("%Y%m%d")))
                        with my_Dataset(input_filename, 'r') as ifid:
                            emis = ifid.variables['emission_anomaly'][:]
                        v = ogid_cat.createVariable('emission_anomaly', np.float64, ('timesteps','latitude','longitude'))
                        v[:] = emis
                        os.remove(input_filename)
                ofid.close()
