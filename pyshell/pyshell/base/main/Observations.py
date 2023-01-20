#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
For lack of a place to put, I'm putting the gory details of the datetime.strftime()
accepted formats here:

    %a - abbreviated weekday name
    %A - full weekday name
    %b - abbreviated month name
    %B - full month name
    %c - preferred date and time representation
    %C - century number (the year divided by 100, range 00 to 99)
    %d - day of the month (01 to 31)
    %D - same as %m/%d/%y
    %e - day of the month (1 to 31)
    %g - like %G, but without the century
    %G - 4-digit year corresponding to the ISO week number (see %V).
    %h - same as %b
    %H - hour, using a 24-hour clock (00 to 23)
    %I - hour, using a 12-hour clock (01 to 12)
    %j - day of the year (001 to 366)
    %m - month (01 to 12)
    %M - minute
    %n - newline character
    %p - either am or pm according to the given time value
    %r - time in a.m. and p.m. notation
    %R - time in 24 hour notation
    %S - second
    %t - tab character
    %T - current time, equal to %H:%M:%S
    %u - weekday as a number (1 to 7), Monday=1. Warning: In Sun Solaris Sunday=1
    %U - week number of the current year, starting with the first Sunday as the first day of the first week
    %V - The ISO 8601 week number of the current year (01 to 53), where week 1 is the first week that has at least 4 days in the current year, and with Monday as the first day of the week
    %W - week number of the current year, starting with the first Monday as the first day of the first week
    %w - day of the week as a decimal, Sunday=0
    %x - preferred date representation without the time
    %X - preferred time representation without the date
    %y - year without a century (range 00 to 99)
    %Y - year including the century
    %Z or %z - time zone or name or abbreviation
    %% - a literal % character

"""

import sys
sys.dont_write_bytecode = True

from pyshell.tmflex import rc
import calendar, itertools, progressbar, glob
import re, os, shutil, tempfile, hashlib
from dateutil.relativedelta import relativedelta
from random import sample
from time import mktime
from numpy import *
from pyshell.base.helper.Utilities import *
from datetime import datetime, timedelta, time
from netCDF4 import Dataset, OrderedDict, date2num, num2date
# from pyshell.tm5_utils import interpolate_fields
from collections import defaultdict
import cPickle as pickle
from warnings import warn

from pyshell.base.main.Emissions import del_time

class Observations(object):

    def __init__(self, StartTime, EndTime):
        super(Observations, self).__init__()
        if type(StartTime) == datetime:
            self.StartDate = StartTime
        else:
            self.StartDate = datetime(*StartTime)
        if type(EndTime) == datetime:
            self.EndDate = EndTime
        else:
            self.EndDate = datetime(*EndTime)
        # Define strings for debug printing
        self.starttime_string = self.StartDate.strftime("%Y-%m-%d %H:%M:%S")
        self.endtime_string = self.EndDate.strftime("%Y-%m-%d %H:%M:%S")
        self.rcf = rc.RcFile(os.environ['pyshell.rc'])
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        #self.point_file = self.rcf.get('output.point.infile')
        self.GetZoomRegions()
        self.StartTime = self.StartDate
        self.EndTime = self.EndDate
        # and where is the input list?
        self.station_name_list = []
        self.obspack_cats = {}
        #self.stationlist_file = self.rcf.get('input.station.filename')
        #self.getStationList(self.stationlist_file)
        #if self.rcf.get('output.satellite', 'bool'):
        self.sat_split_period = self.rcf.get('output.satellite.split.period', default='d')
        #if self.rcf.get('output.point', 'bool'):
        self.point_split_period = self.rcf.get('output.point.split.period', default='a')
        # hack to shut up things which need subdir_tag
        self.subdir_tag = ''

    def reverseLookup(self):
        ret_dict = OrderedDict()

        for site_cat, site_cat_data in self.obspack_cats.items():
            for site_code, site_data in site_cat_data.items():
                if site_data['n_obs'] > 0:
                    site_id = site_data['site_id']
                    ret_dict[site_id] = {'site_code': site_code, 'site_category': site_cat}
                    for k,v in site_data.items():
                        if k != 'site_id':
                            ret_dict[site_id][k] = v

        return ret_dict

    def writeStationFile(self):
        # Create the lines first
        write_lines = [' ID     LAT     LON     ALT TP STATION NAME\n']
        for site_cat, site_cat_data in self.obspack_cats.items():
            for site_code, site_data in site_cat_data.items():
                if site_data['n_obs'] > 0 and site_data['static']:
                    line = '%3s %7.2f %7.2f %7.1f %2s %s\n'%\
                        (site_code.upper(), site_data['latitude'], site_data['longitude'], site_data['altitude'], \
                        site_data['type'], site_data['name'])
                    write_lines.append(line)

        file_name = self.rcf.get('output.station.timeseries.filename')
        with open(file_name, 'w') as fid:
            fid.writelines(write_lines)

    def updateSiteIDFile(self, out_file, site_id_dict):
        # When we are writing out point observations one file per day, we do not want to overwrite the dictionary
        # from the previous day. Therefore, we read it in, update it with self.site_id_dict, and write it back.
        checkDir(out_file)
        existing_dict = {}
        if os.path.exists(out_file):
            with open(out_file, 'rb') as fid:
                existing_dict = pickle.load(fid)
            os.remove(out_file)
        existing_dict.update(site_id_dict)
        with open(out_file, 'wb') as fid:
            pickle.dump(existing_dict, fid, pickle.HIGHEST_PROTOCOL)

    def decimal_date_to_datetime(self, decimal_date):
        """
        Given a single decimal date such as 2003.4932, convert it to a datetime object datetime(2003, 6, 30, 0, 25, 55, 200000)
        """
        year = int(decimal_date)
        days = (365.0+calendar.isleap(year)) * (decimal_date - year)
        return datetime(year,1,1,0) + timedelta(days=days)

    def filter_by_local_time(self, datetime_array, lon_array, site_type):
        """
        Given a list of datetimes and longitudes, return a list of indices of the datetimes corresponding
        to either mid-afternoon or mid-night samples, as required.
        """
        local_times = [d+timedelta(hours=lon/15.0) for d,lon in zip(datetime_array, lon_array)]
        if site_type == 'afternoon':
            allowed_hours = [11, 12, 13, 14, 15, 16]
        elif site_type == 'night':
            allowed_hours = [0, 1, 2, 3, 4, 5, 6]
        else:
            allowed_hours = []
        valid_indices = array([i for i,d in enumerate(local_times) if d.hour in allowed_hours])
        return valid_indices

    def decimal_date(self, date_array):
        """
        Given an array or list of datetime objects, convert them to decimal dates. For example, datetime(2009,3,15,9,12,30)
        would be converted to 2009.2010512.
        """
        # check if date_array is a single object
        if isinstance(date_array, datetime):
            return self.decimal_date_atomic(date_array)
        else:
            ret_array = array([self.decimal_date_atomic(d) for d in date_array], float64)
        return ret_array

    def decimal_date_atomic(self, date_obj):
        # This is called by decimal_date to convert a single datetime object into a decimal date
        year = date_obj.year
        dt = date_obj - datetime(year,1,1)
        secs_in_year = dt.seconds + 86400.0 * dt.days
        total_secs = 86400.0 * (int(calendar.isleap(year)) + 365)
        return float(year) + secs_in_year/total_secs

    def shorten_path(self, file_name, keep_N):
        # Given a really long file name with path, shorten to keep the last keep_N elements, including the base name
        atoms = file_name.split(os.sep)
        if len(atoms) <= keep_N:
            return file_name
        else:
            return os.sep.join(atoms[-keep_N:])

    def get_YM_tuples(self):
        # What are the year/month tuples for the current run period?
        cur_date = datetime(self.StartDate.year, self.StartDate.month, 1, 0, 0, 0)
        ym_list = []
        while cur_date < self.EndDate:
            ym_list.append((cur_date.year, cur_date.month))
            cur_date = cur_date + relativedelta(months=1)
        return ym_list

    def get_YMD_tuples(self):
        # What are the year/month/day tuples for the current run period?
        cur_date = datetime(self.StartDate.year, self.StartDate.month, self.StartDate.day, 0, 0, 0)
        ymd_list = []
        while cur_date < self.EndDate:
            ymd_list.append(cur_date.timetuple()[:3])
            cur_date += timedelta(days=1)
        return ymd_list

    def getStationList(self, fileName):
        self.station_name_dict = OrderedDict()
        self.station_name_list = []
        self.station_coords = {}

        if not os.path.exists(fileName):
            warn("Station list %s does not exist"%fileName, RuntimeWarning, stacklevel=2)
            lines = []
        else:
            with open(fileName, 'r') as fid:
                lines = fid.readlines()
            # ignore header line
            lines = lines[1:]

        for line in lines:
            stat_ID = line.split()[0].upper()
            stat_type = line.split()[4].upper()
            stat_name = ' '.join(line.split()[5:])
            stat_lat = float(line.split()[1])
            stat_lon = float(line.split()[2])
            stat_alt = float(line.split()[3])
            self.station_coords[(stat_ID, stat_type)] = {'lat': stat_lat, 'lon': stat_lon, 'alt': stat_alt}
            self.station_name_dict[(stat_ID, stat_type)] = stat_name
            self.station_name_list.append((stat_ID, stat_type))

    def get_class_from_name(self, class_name):
        _temp = __import__('Observations', fromlist=[class_name])
        try:
           class_from_name = _temp.__dict__[class_name]
           return class_from_name
        except KeyError:
           sys.stderr.write("Class %s not defined in %s.\n"%(class_name, 'Observations'))
           sys.exit()
        except:
           sys.stderr.write("Unknown error importing %s\n"%class_name)
           sys.exit()

    def create_point_observation_structure( self):
        # create dictionary to contain observation info to be read from rc file:
        self.PointObservation = {}
        for tracer in self.species:
            self.PointObservation[tracer] = {}
            self.PointObservation[tracer]['obs_class'] = self.rcf.get(tracer+'.obs.point.class')

    def create_sat_observation_structure( self):
        # create dictionary to contain observation info to be read from rc file:
        self.SatObservation = {}
        for tracer in self.species:
            self.SatObservation[tracer] = {}
            self.SatObservation[tracer]['obs_class'] = self.rcf.get(tracer+'.obs.sat.class')

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
            # Store the dx and dy values as well
            self.zoom_info[region]['dx'] = (self.zoom_info[region]['xend'] - self.zoom_info[region]['xbeg'])/self.zoom_info[region]['im']
            self.zoom_info[region]['dy'] = (self.zoom_info[region]['yend'] - self.zoom_info[region]['ybeg'])/self.zoom_info[region]['jm']
        self.zoom_regions_lon = zip(xbeg,xend,im)
        self.zoom_regions_lat = zip(ybeg,yend,jm)

    def GetZoomRegion(self,loc_array):
        """
        Given an array of (lat, lon) pairs, returns an array of region indices where they belong
        """
        refinement_factors = array(self.zoom_region_xref[1:]) * array(self.zoom_region_yref[1:])
        region_order = argsort(refinement_factors)[::-1]
        indices = zeros(len(loc_array), int16) # the output array
        region_indices = range(1,len(self.zoom_regions_names)+1) # 1,2,3,4,... the index of the zoom region
        for i,location in enumerate(loc_array):
            for idx in region_order:
                if self.zoom_regions_lat[idx,0] <= location[0] <= self.zoom_regions_lat[idx,1] and self.zoom_regions_lon[idx,0] <= location[1] <= self.zoom_regions_lon[idx,1]:
                    break
            indices[i] = region_indices[idx]
        return indices

    #def InflateErrors(self):
        #"""
        #For inflating errors in total column measurements, based on some binning length and binning time
        #"""
        #if len(self.CO2_sigma_column_mixing) > 0:
            ## convert sample times to floating point seconds
            #sample_seconds = [t-self.sample_times[0] for t in self.sample_times]
            #sample_seconds = [dt.days * 86400 + dt.seconds for dt in sample_seconds]
            #locations = zeros((len(self.Latitudes),2), self.Latitudes.dtype)
            #locations[:,0] = self.Latitudes
            #locations[:,1] = self.Longitudes
            #self.CO2_sigma_column_mixing, num_samples = interpolate_fields.movingAverage(float64(sample_seconds), self.CO2_sigma_column_mixing, locations, self.binning_time, self.binning_length)
            #del sample_seconds, locations, num_samples
        #else:
            #raise ValueError('No samples in this month, please investigate!')

    def putDateString(self, input_string, file=True):
        output_string = input_string
        for atom in ['Y', 'y', 'm', 'M', 'd', 'H', 'b', 'B', 'j', 'p', 'S']:
            for i, t in zip([1,2], [self.StartTime, self.EndTime]):
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

    def copyOutputToInput(self):
        """
        For an OSSE, we create a bunch of input files with sampling locations prescribed by the virtual network,
        with dummy concentrations. Then we run TM5 forward with 'true' fluxes to produce pseudo-obs. This routine
        copies over the pseudo-obs to the input files. Also, to be compatible
        """
        # First, the point files
        if self.rcf.get('output.point', 'bool'):
            track_dir = self.putDateString(self.rcf.get('output.dir'))
            track_dir = os.path.join(track_dir, 'point')
            all_track_files = glob.glob(os.path.join(track_dir, 'point_output*.nc4'))
            input_dir = self.rcf.get('output.point.input.dir')
            for track_file in sorted(all_track_files):
                bname = os.path.basename(track_file)
                bname = bname.replace('output', 'input')
                input_file = os.path.join(input_dir, bname)
                fid_i = Dataset(input_file, 'a')
                fid_o = Dataset(track_file, 'r')
                for tracer in self.species:
                    in_mix = fid_i.groups[tracer].variables['mixing_ratio'][:]
                    n_obs = 0
                    for region in fid_o.groups:
                        ipos = fid_o.groups[region].groups[tracer].variables['id'][:] - 1
                        out_mix = fid_o.groups[region].groups[tracer].variables['mixing_ratio'][:]
                        in_mix[ipos] = out_mix
                        n_obs += len(fid_o.groups[region].groups[tracer].dimensions['samples'])
                    fid_i.groups[tracer].variables['mixing_ratio'][:] = in_mix
                    # replace the instantaneous sampling strategy with a symmetric one
                    n_in = len(fid_i.groups[tracer].dimensions['id'])
                    samp_strat = 3 * ones(n_in, int16)
                    fid_i.groups[tracer].variables['sampling_strategy'][:] = samp_strat
                    # check if total number of obs match
                    if n_obs != len(fid_i.groups[tracer].dimensions['id']):
                        print "========================="
                        print "Number of observations in %s = %i"%(os.path.basename(track_file), n_obs)
                        print "Number of observations in %s = %i"%(os.path.basename(input_file), len(fid_i.groups[tracer].dimensions['id']))
                        print "========================="
                setattr(fid_i, 'obs_simulated', 1)
                fid_i.close()
                fid_o.close()
                print "%s done"%os.path.basename(input_file)
        # Now, satellite files
        if self.rcf.get('output.satellite', 'bool'):
            track_dir = self.putDateString(self.rcf.get('output.dir'))
            track_dir = os.path.join(track_dir, 'satellite')
            all_dep_files = glob.glob(os.path.join(track_dir, 'sat-track_departures*.nc4'))
            input_dir = self.rcf.get('output.satellite.output.directory')
            for dep_file in sorted(all_dep_files):
                bname = os.path.basename(dep_file)
                last_part = bname.split('_')[2] # something like '20100222.nc4'
                track_file = os.path.join(track_dir, 'sat-track_%s'%last_part)
                input_file = os.path.join(input_dir, 'inputfile_%s'%last_part)
                fid_i = Dataset(input_file, 'a')
                fid_d = Dataset(dep_file, 'r')
                fid_t = Dataset(track_file, 'r')
                for tracer in self.species:
                    in_mix = fid_i.groups[tracer].variables['column_mixing'][:]
                    n_obs = 0
                    for region in fid_d.groups:
                        ipos = fid_t.groups[region].groups[tracer].variables['input_positions'][:] - 1
                        out_mix = fid_d.groups[region].groups[tracer].variables['modeled_column'][:]
                        in_mix[ipos] = out_mix
                        n_obs += len(fid_d.groups[region].groups[tracer].dimensions['n_obs'])
                    fid_i.groups[self.tracer].variables['column_mixing'][:] = in_mix
                    # for satellites, the sampling strategy was symmetric to begin with, so change nothing
                    # check if the total number of obs match
                    if n_obs != len(fid_i.groups[tracer].dimensions['n_obs']):
                        print "========================="
                        print "Number of observations in %s = %i"%(os.path.basename(dep_file), n_obs)
                        print "Number of observations in %s = %i"%(os.path.basename(input_file), len(fid_i.groups[tracer].dimensions['n_obs']))
                        print "========================="
                setattr(fid_i, 'obs_simulated', 1)
                fid_i.close()
                fid_d.close()
                fid_t.close()
                print "%s done"%os.path.basename(input_file)

    def writePointFile(self):
        """
        This overwrites the basic writePointFile, so that the point file gets a YYYYMM[DD] tag.
        """
        dir_name = self.rcf.get('output.point.input.dir')
        if self.point_split_period == 'a':
            file_name = 'point_input.nc4'
        elif self.point_split_period == 'm':
            file_name = self.StartDate.strftime("point_input_%Y%m.nc4")
        elif self.point_split_period == 'd':
            file_name = self.StartDate.strftime("point_input_%Y%m%d.nc4")
        point_file = os.path.join(dir_name, file_name)

        if os.path.exists(point_file):
            os.remove(point_file)
        # only create a file if there are observations to write
        total_obs = 0
        for tracer in self.species:
            total_obs += self.PointObservation[tracer]['dimensions']['id']
        if total_obs == 0:
            return

        checkDir(point_file)
        file_id = Dataset(point_file, 'w')
        file_id.createDimension('idate', 6)
        for tracer in self.species:
            #print 'Writing observations for tracer ', tracer
            # If a certain tracer is not to be assimilated, set its error very high
            assim_tracer = self.rcf.get('%s.point.assimilate'%tracer, 'bool', default=True)
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

        print "Wrote %s"%point_file

    def writeSatFile(self):
        output_dir = self.rcf.get('output.satellite.output.directory')
        output_pfx = self.rcf.get('output.satellite.ipfile.prefix')
        if self.sat_split_period == 'm':
            fileName = os.path.join(output_dir, output_pfx + self.StartDate.strftime("%Y%m") + ".nc4")
        elif self.sat_split_period == 'd':
            fileName = os.path.join(output_dir, output_pfx + self.StartDate.strftime("%Y%m%d") + ".nc4")
        if os.path.exists(fileName):
            os.remove(fileName)

        # only create a file if there are observations to write
        total_obs = 0
        for tracer in self.species:
            total_obs += self.SatObservation[tracer]['dimensions']['n_obs']
        if total_obs == 0:
            return

        checkDir(fileName)
        fid = Dataset(fileName, 'w')
        fid.createDimension('idate', 6)
        for tracer in self.species:
            if self.SatObservation[tracer]['dimensions']['n_obs'] > 0:
                gid = fid.createGroup(tracer)
                for dim_name, dim_len in self.SatObservation[tracer]['dimensions'].items():
                    gid.createDimension(dim_name, dim_len)
                for attr_name, attr_value in self.SatObservation[tracer]['file_attrs'].items():
                    setattr(gid, attr_name, attr_value)
                for var_name, var_shape in self.SatObservation[tracer]['variable_shapes'].items():
                    var_value = self.SatObservation[tracer][var_name]
                    v = gid.createVariable(var_name, var_value.dtype, var_shape)
                    v[:] = var_value
                    if var_name in self.SatObservation[tracer]['variable_attrs']:
                        for attr_name, attr_value in self.SatObservation[tracer]['variable_attrs'][var_name]:
                            setattr(v, attr_name, attr_value)
        fid.close()

class MakeSyntheticObs_point(Observations):
    """
    Sometimes it is necessary to make synthetic observations -- for all tracers -- that lie within a specific zoom region.
    For example, I may want to test the 4DVAR system with N observations within the innermost zoom region, M within the
    first layer of halo cells, etc.
    """
    def __init__(self, *args, **kwargs):
        super(MakeSyntheticObs_point, self).__init__(*args, **kwargs)

    def readRcKeys(self):
        # For N zoom regions, there are N-1 halo cell regions
        self.obs_num_dict = {}
        for region in self.zoom_regions_names:
            self.obs_num_dict[region] = self.rcf.get('point.%s.synthetic.num.obs.%s'%(self.tracer,region), 'int', default=0)
        if len(self.zoom_regions_names) > 1:
            for region in self.zoom_regions_names[1:]:
                self.obs_num_dict['%s.halo'%region] = self.rcf.get('point.%s.synthetic.num.obs.%s.halo'%(self.tracer,region), 'int', default=0)
        self.mean_obs = self.rcf.get('point.%s.synthetic.mean.obs'%self.tracer, 'float')
        self.std_obs = self.rcf.get('point.%s.synthetic.std.obs'%self.tracer, 'float')
        self.mdm = self.rcf.get('point.%s.synthetic.mdm'%self.tracer, 'float')
        self.sampling_strategies = [int(s) for s in self.rcf.get('point.%s.synthetic.sampl_strategies'%self.tracer).split()]
        # To create realistic altitudes, we need a topo map
        self.topo_file = self.rcf.get('topo.database')
        self.createTopoMap()
        # Make a dictionary of all children regions for every parent region
        self.child_dict = defaultdict(list)
        for region in self.zoom_regions_names:
            parent = self.rcf.get('region.%s.parent'%region)
            self.child_dict[parent].append(region)

    def in_region(self, lat, lon, region):
        # Checks whether (lat,lon) is in a rectangular box defined by the region boundaries
        # This does not account for halo cells!
        ret_val = (self.zoom_info[region]['xbeg'] < lon < self.zoom_info[region]['xend'])
        ret_val = ret_val and (self.zoom_info[region]['ybeg'] < lat < self.zoom_info[region]['yend'])
        return ret_val

    def createTopoMap(self):
        """
        Create a lat/lon array of the topography in meters, which will be used for creating the altitudes of synthetic observations
        """
        with Dataset(self.topo_file, 'r') as fid:
            topo = fid.variables['TerrainHeight'][:]
            ocean_value = fid.Ocean_value
        self.topo_alt = float64(topo)
        self.topo_alt[where(topo == ocean_value)] = 1.0
        # Now add 100 meters as buffer
        self.topo_alt += 100.0

    def createObs(self, region, num_obs, exclude_parent, halo=False):
        if not halo:
            # Create observations for a region excluding the halo cells, and excluding all child regions
            if exclude_parent:
                parent = self.rcf.get('region.%s.parent'%region)
                xbeg = self.zoom_info[region]['xbeg'] + self.zoom_info[parent]['dx'] + 0.00001
                xend = self.zoom_info[region]['xend'] - self.zoom_info[parent]['dx']
                ybeg = self.zoom_info[region]['ybeg'] + self.zoom_info[parent]['dy'] + 0.00001
                yend = self.zoom_info[region]['yend'] - self.zoom_info[parent]['dy']
            else:
                xbeg = self.zoom_info[region]['xbeg'] + 0.00001
                xend = self.zoom_info[region]['xend']
                ybeg = self.zoom_info[region]['ybeg'] + 0.00001
                yend = self.zoom_info[region]['yend']
            # now make the obs, excluding child regions
            random_lats = zeros(num_obs, float64)
            random_lons = zeros(num_obs, float64)
            valid_obs = 0
            children = self.child_dict[region]
            while valid_obs < num_obs:
                lat = (yend-ybeg) * random.random_sample() + ybeg
                lon = (xend-xbeg) * random.random_sample() + xbeg
                in_child = False
                for child in children:
                    in_child = in_child or self.in_region(lat, lon, child)
                if not in_child:
                    random_lats[valid_obs] = lat
                    random_lons[valid_obs] = lon
                    valid_obs += 1
        else:
            # Create observations only within the halo region
            parent = self.rcf.get('region.%s.parent'%region)
            xbeg_out = self.zoom_info[region]['xbeg'] + 0.00001
            xbeg_in  = xbeg_out + self.zoom_info[parent]['dx']
            xend_out = self.zoom_info[region]['xend']
            xend_in  = xend_out - self.zoom_info[parent]['dx']
            ybeg_out = self.zoom_info[region]['ybeg'] + 0.00001
            ybeg_in  = ybeg_out + self.zoom_info[parent]['dy']
            yend_out = self.zoom_info[region]['yend']
            yend_in  = yend_out - self.zoom_info[parent]['dy']
            valid_obs = 0
            random_lats = zeros(num_obs, float64)
            random_lons = zeros(num_obs, float64)
            while valid_obs < num_obs:
                lat = (yend_out-ybeg_out) * random.random_sample() + ybeg_out
                lon = (xend_out-xbeg_out) * random.random_sample() + xbeg_out
                if lat < ybeg_in or lat > yend_in or lon < xbeg_in or lon > xend_in:
                    random_lats[valid_obs] = lat
                    random_lons[valid_obs] = lon
                    valid_obs += 1
        return random_lats, random_lons

    def createObservations(self):
        # Generate the lats and lons first
        lats = []
        lons = []
        for region in self.zoom_regions_names:
            globe = self.rcf.get('region.%s.parent'%region) == 'globe'
            n_obs = self.obs_num_dict[region]
            if globe:
                r_lats, r_lons = self.createObs(region, n_obs, False, False)
            else:
                r_lats, r_lons = self.createObs(region, n_obs, True, False)
            lats.extend(list(r_lats))
            lons.extend(list(r_lons))
        if len(self.zoom_regions_names) > 1:
            for region in self.zoom_regions_names[1:]:
                n_obs = self.obs_num_dict['%s.halo'%region]
                r_lats, r_lons = self.createObs(region, n_obs, None, True)
                lats.extend(list(r_lats))
                lons.extend(list(r_lons))
        # How many observations in total?
        n_obs = len(lats)
        # Create the times
        time_start = date2num(self.StartDate, units='seconds since 1900-01-01 00:00:00-0:00', calendar='standard') + 0.0001
        time_end = date2num(self.EndDate, units='seconds since 1900-01-01 00:00:00-0:00', calendar='standard')
        r_times = (time_end-time_start) * random.random_sample(n_obs) + time_start
        r_times.sort()
        r_times = [num2date(d, units='seconds since 1900-01-01 00:00:00-0:00', calendar='standard') for d in r_times]
        # Store them in a dictionary
        self.write_dict = {}
        self.write_dict['times'] = array([d.timetuple()[:6] for d in r_times], int16)
        self.write_dict['lat'] = array(lats, float64)
        self.write_dict['lon'] = array(lons, float64)
        nlat, nlon = self.topo_alt.shape
        lat_array = linspace(-90.,90.,nlat+1)
        lon_array = linspace(-180.,180.,nlon+1)
        self.write_dict['alt'] = zeros(n_obs, float64)
        for i, (lat, lon) in enumerate(zip(lats,lons)):
            lat_idx = lat_array.searchsorted(lat) - 1
            lon_idx = lon_array.searchsorted(lon) - 1
            self.write_dict['alt'][i] = self.topo_alt[lat_idx, lon_idx]
        # Now create the mixing ratios
        self.write_dict['mix'] = self.std_obs * random.standard_normal(n_obs) + self.mean_obs
        self.write_dict['mix_std'] = self.mdm * ones(n_obs, float64)
        # Sampling strategy should be chosen from the list self.sampling_strategies
        self.write_dict['samp_strat'] = int16(random.choice(self.sampling_strategies, size=n_obs))
        # For sampling strategy 4, we need a specification of the time window lengths
        random_window_lengths = random.random_integers(3600, 604800, size=n_obs) # between 1 hour and 1 week
        self.write_dict['tw_length'] = where(self.write_dict['samp_strat'] == 4, random_window_lengths, 0)

    def __call__(self):
        self.createObservations()
        n_obs = len(self.write_dict['mix'])
        self.PointObservation[self.tracer]['dimensions'] = {'id': n_obs}
        self.PointObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.PointObservation[self.tracer]['variable_attrs'] = OrderedDict()
        self.PointObservation[self.tracer]['attr_dict'] = OrderedDict()

        self.PointObservation[self.tracer]['id'] = arange(1, n_obs+1, dtype=int32)
        self.PointObservation[self.tracer]['variable_shapes']['id'] = ('id',)

        self.PointObservation[self.tracer]['lat'] = self.write_dict['lat']
        self.PointObservation[self.tracer]['variable_shapes']['lat'] = ('id',)

        self.PointObservation[self.tracer]['lon'] = self.write_dict['lon']
        self.PointObservation[self.tracer]['variable_shapes']['lon'] = ('id',)

        self.PointObservation[self.tracer]['alt'] = self.write_dict['alt']
        self.PointObservation[self.tracer]['variable_shapes']['alt'] = ('id',)

        self.PointObservation[self.tracer]['mixing_ratio'] = float64(self.write_dict['mix'])
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio'] = ('id',)

        self.PointObservation[self.tracer]['mixing_ratio_error'] = self.write_dict['mix_std']
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio_error'] = ('id',)

        self.PointObservation[self.tracer]['date_components'] = self.write_dict['times']
        self.PointObservation[self.tracer]['variable_shapes']['date_components'] = ('id', 'idate')

        self.PointObservation[self.tracer]['sampling_strategy'] = self.write_dict['samp_strat']
        self.PointObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('id',)

        self.PointObservation[self.tracer]['time_window_length'] = int32(self.write_dict['tw_length'])
        self.PointObservation[self.tracer]['variable_shapes']['time_window_length'] = ('id',)

        self.PointObservation[self.tracer]['station_id'] = arange(1, n_obs+1, dtype=int32)
        self.PointObservation[self.tracer]['variable_shapes']['station_id'] = ('id',)

class MakePerfectObs(Observations):
    """
    We would like to test a CO+CO2 inversion system against individual CO and CO2 systems. The problem is that since we
    are not iterating up to perfect convergence, the intermediate state of the joint tracer system will be different
    from the intermediate states of the single tracer systems. So Maarten suggested that I could manufacture observations
    for one of the tracers in a joint tracer system that are perfect, so only the second tracer contributes to the cost
    function. To take an example:

    (1) Run CO forward with prior fluxes and create modeled total columns and point observations.
    (2) Copy over those modeled total columns and point observations into CO observations assimilated by a CO+CO2
        system.
    (3) Run a CO+CO2 system for N iterations. Since the CO observations match the prior CO fluxes, the CO flux
        should not be changed. Only the CO2 flux should be changed.
    (4) The posterior CO2 flux from the aforementioned CO+CO2 optimization should be identical to the posterior
        CO2 flux from a CO2-only optimization after N iterations.

    This class contains methods to do step (2), i.e., copy over modeled observations of one tracer into the observations
    to be assimilated in a joint tracer inversion.
    """
    def __init__(self, *args, **kwargs):
        super(MakePerfectObs, self).__init__(*args, **kwargs)
        self.tracer_name = self.rcf.get('my.tracer.name')
        # hack to shut up things which need subdir_tag
        self.subdir_tag = ''
        # if we want to ignore an obs during assimilation, we set the error for
        # that particular obs to this factor times the maximum observed value
        self.inflate_error_factor = self.rcf.get('unassimilated.obs.error.factor', 'float', default=1.0E15)

    def readPointObs(self, tracer, base='var4d'):
        # read modeled point obs from 1-tracer output folders
        joint_op_dir = self.putDateString(self.rcf.get('output.dir'))
        # both the base and the tracer are potentially different
        tracer_op_dir = joint_op_dir.replace(self.tracer_name, tracer)
        tracer_op_dir = tracer_op_dir.replace(self.rcf.get('my.branch'), base)
        id_list = []
        meas_list = []
        ip_file = os.path.join(tracer_op_dir, 'point_output.nc4')
        print "Reading point observations from %s for tracer %s"%(ip_file, tracer)
        with Dataset(ip_file, 'r') as fid:
            for region in fid.groups.keys():
                if len(fid.groups[region].dimensions['samples']) > 0:
                    id_list.extend(fid.groups[region].variables['id'][:] - 1)
                    meas_list.extend(fid.groups[region].variables['mixing_ratio'][:])
        return {'id': array(id_list), 'model': array(meas_list)}

    def readSatObs(self, tracer, base='var4d'):
        # read modeled total columns from 1-tracer output folders
        joint_op_dir = self.putDateString(self.rcf.get('output.dir'))
        # both the base and the tracer are potentially different
        tracer_op_dir = joint_op_dir.replace(self.tracer_name, tracer)
        tracer_op_dir = tracer_op_dir.replace(self.rcf.get('my.branch'), base)
        glob_patt = os.path.join(tracer_op_dir, 'sat-track_departures_[0-9]*.nc4')
        dep_files = glob.glob(glob_patt)
        ret_dict = {}
        for file_name in dep_files:
            ymd_string = os.path.splitext(file_name)[0].split('_')[-1]
            meas_list = []
            id_list = []
            print "Reading satellite observations from %s for tracer %s"%(file_name, tracer)
            with Dataset(file_name, 'r') as fid:
                for region in fid.groups.keys():
                    valid = fid.groups[region].variables['valid_departure'][:]
                    idx = where(valid == 1)[0]
                    ids = fid.groups[region].variables['input_index'][:][idx]
                    cols = fid.groups[region].variables['modeled_column'][:][idx]
                    meas_list.extend(cols)
                    id_list.extend(ids)
            ret_dict[ymd_string] = {'id': array(id_list), 'model': array(meas_list)}
        return ret_dict

    #def writePointObs(self, tracer, point_dict):
        ## write out the 'perfect' point observations to the 2-tracer point input file
        #print "Writing observations for tracer %s to %s"%(tracer, self.point_file)
        #fid = Dataset(self.point_file, 'a')
        #obs = fid.groups[tracer].variables['mixing_ratio'][:]
        #obs[point_dict['id']] = point_dict['model']
        #fid.groups[tracer].variables['mixing_ratio'][:] = obs
        #obs_err = fid.groups[tracer].variables['mixing_ratio_error'][:]
        #uniform_large_error = self.inflate_error_factor * abs(point_dict['model']).max()
        #obs_err[point_dict['id']] = uniform_large_error
        #fid.groups[tracer].variables['mixing_ratio_error'][:] = obs_err
        #fid.close()

    def writeSatObs(self, tracer, sat_dict):
        # write out the 'perfect' sat observtions to the 2-tracer sat input files
        input_dir = os.path.join(self.rcf.get('input.dir'), 'satellites')
        for ymd_tuple, sat_data in sat_dict.items():
            sat_file = os.path.join(input_dir, 'inputfile_'+ymd_tuple+'.nc4')
            print "Writing observations for tracer %s to %s"%(tracer, sat_file)
            fid = Dataset(sat_file, 'a')
            obs = fid.groups[tracer].variables['column_mixing'][:]
            obs[sat_data['id']] = sat_data['model']
            fid.groups[tracer].variables['column_mixing'][:] = obs
            obs_err = fid.groups[tracer].variables['sigma_column_mixing'][:]
            uniform_large_error = self.inflate_error_factor * abs(sat_data['model']).max()
            obs_err[sat_data['id']] = uniform_large_error
            fid.groups[tracer].variables['sigma_column_mixing'][:] = obs_err
            fid.close()

    def copyObs(self, tracer):
        data_dict = self.readPointObs(tracer)
        self.writePointObs(tracer, data_dict)
        data_dict = self.readSatObs(tracer)
        self.writeSatObs(tracer, data_dict)
        del data_dict
