#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, re, shutil, warnings
sys.dont_write_bytecode = True

import numpy as np
from scipy import interpolate
from glob import glob
from netCDF4 import Dataset, OrderedDict
from datetime import datetime, timedelta, time
from collections import defaultdict
from pyshell.base.helper.Utilities import checkDir
from calendar import monthrange
import cPickle as pickle
import h5py

warnings.simplefilter('error', UserWarning)

from pyshell.base.main.Observations import Observations
from pyshell.base.main.Observations import MakeSyntheticObs_point, del_time

# Import one-off modules, just in case
#from Obs_modules import Amazon_OSSE_create_point, Amazon_OSSE_tracks, Brazilian_future, OCO2_Odell_simul
#from Obs_validation import CONTRAIL_CO2, HIPPO_CO2, Validation

class Synthetic_CO2_point(MakeSyntheticObs_point):

    def __init__(self, *args, **kwargs):
        super(Synthetic_CO2_point, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.readRcKeys()

class Validate_points(Observations):

    def __init__(self, *args, **kwargs):
        super(Validate_points, self).__init__(*args, **kwargs)
        self.validate_methods = self.rcf.get('CO2.point.validate.methods').split()
        self.tracer = 'CO2'

    def writeStationFile(self):
        pass

    def consolidateObs(self):
        ret_dict = defaultdict(list)
        needed_vars = ['date', 'co2', 'lat', 'lon', 'alt', 'co2_err', 'samp_strat', 'site_id']

        for val_method in self.validate_methods:
            sys.stdout.write('Getting obs from method %s ... '%val_method) ; sys.stdout.flush()

            ValClass = self.get_class_from_name(val_method)
            val_inst = ValClass(self.StartDate.timetuple()[:6], self.EndDate.timetuple()[:6])
            data_dict = val_inst.createObservations()

            for var_name in needed_vars:
                ret_dict[var_name].extend(data_dict[var_name])

            sys.stdout.write('done!\n') ; sys.stdout.flush()

        # sort by time and convert into arrays
        sort_order = np.argsort(ret_dict['date'])
        for k,v in ret_dict.items():
            ret_dict[k] = np.array(v)[sort_order]

        return ret_dict

    def __call__(self):
        data_dict = self.consolidateObs()
        n_obs = len(data_dict['co2'])

        self.PointObservation[self.tracer]['dimensions'] = {'id': n_obs}
        self.PointObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.PointObservation[self.tracer]['attr_dict'] = OrderedDict()
        self.PointObservation[self.tracer]['variable_attrs'] = OrderedDict()

        self.PointObservation[self.tracer]['id'] = np.arange(1, n_obs+1, dtype=np.int32)
        self.PointObservation[self.tracer]['variable_shapes']['id'] = ('id',)

        self.PointObservation[self.tracer]['lat'] = np.float64(data_dict['lat'])
        self.PointObservation[self.tracer]['variable_shapes']['lat'] = ('id',)

        self.PointObservation[self.tracer]['lon'] = np.float64(data_dict['lon'])
        self.PointObservation[self.tracer]['variable_shapes']['lon'] = ('id',)

        self.PointObservation[self.tracer]['alt'] = np.float64(data_dict['alt'])
        self.PointObservation[self.tracer]['variable_shapes']['alt'] = ('id',)

        self.PointObservation[self.tracer]['mixing_ratio'] = np.float64(data_dict['co2'])
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio'] = ('id',)

        self.PointObservation[self.tracer]['mixing_ratio_error'] = np.float64(data_dict['co2_err'])
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio_error'] = ('id',)

        self.PointObservation[self.tracer]['date_components'] = np.array([d.timetuple()[:6] for d in data_dict['date']], np.int16)
        self.PointObservation[self.tracer]['variable_shapes']['date_components'] = ('id', 'idate')

        self.PointObservation[self.tracer]['sampling_strategy'] = np.int16(data_dict['samp_strat'])
        self.PointObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('id',)

        self.PointObservation[self.tracer]['station_id'] = np.int32(data_dict['site_id'])
        self.PointObservation[self.tracer]['variable_shapes']['station_id'] = ('id',)

class SEAC4RS_CO2(Observations):
    """
    This class parses CO2 observations from the SEAC4RS campaign
    """
    def __init__(self, *args, **kwargs):
        super(SEAC4RS_CO2, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.input_dir = self.rcf.get('SEAC4RS.CO2.data.folder')
        self.flights = self.rcf.get('SEAC4RS.CO2.flights').split() # 'DC8' and 'ER2', at most
        self.samp_strat = self.rcf.get('SEAC4RS.CO2.sampling.strategy', 'int')
        self.samp_err = self.rcf.get('SEAC4RS.CO2.error', 'float')
        self.assimilate = self.rcf.get('SEAC4RS.CO2.assimilate', 'bool', default=False)
        self.site_id = {}
        for aircraft in self.flights:
            self.site_id[aircraft] = self.rcf.get('SEAC4RS.CO2.%s.site_id'%aircraft, 'int')
        if not self.assimilate:
            self.unassim_mdm = self.rcf.get('CO2.point.unassimilate.mdm', 'float')

    def parseFile(self, file_name):
        ret_dict = {}

        # how many lines to skip?
        with open(file_name, 'r') as fid:
            comm_line = fid.readline()
            skiprows = int(comm_line.split(',')[0]) - 1 # need the header row with names of the columns

        temp_arr = np.genfromtxt(file_name, autostrip=True, names=True, skip_header=skiprows, delimiter=',', \
            usecols=('Fractional_Day', 'LATITUDE', 'LONGITUDE', 'GPS_ALT_NASDAT', 'PRESSURE', 'TEMPERATURE', 'CO2_HUPCRS', \
            'CO_HUPCRS', 'CH4_HUPCRS'), missing_values='-999999', usemask=True)

        # filter out lines of invalid CO2 data and coordinates
        valid_coords = ['Fractional_Day', 'LATITUDE', 'LONGITUDE', 'GPS_ALT_NASDAT', 'CO2_HUPCRS']
        valid_idx = np.zeros(len(temp_arr), bool)
        valid_idx[:] = True
        for col_name in valid_coords:
            if isinstance(temp_arr[col_name], np.ma.MaskedArray):
                valid_idx = np.logical_and(valid_idx, ~temp_arr[col_name].mask)
        temp_arr = temp_arr[valid_idx]

        # now filter out rows outside the time range of the run
        year = int(os.path.basename(file_name).split('_')[2][:4]) # f**king text format weirdos
        utc_times = np.array([datetime(year,1,1) + timedelta(days=d) for d in temp_arr['Fractional_Day']])
        valid_idx = np.array([i for i,d in enumerate(utc_times) if self.StartDate <= d <= self.EndDate])

        if len(valid_idx) > 0:
            temp_arr            = temp_arr[valid_idx]
            ret_dict['date']    = utc_times[valid_idx]
            ret_dict['lat']     = temp_arr['LATITUDE']
            # longitudes are given from 0 to 360
            ret_dict['lon']     = temp_arr['LONGITUDE'] - np.int32(temp_arr['LONGITUDE']/180.0) * 360.0
            # altitudes are given in km
            ret_dict['alt']     = temp_arr['GPS_ALT_NASDAT'] * 1000.0
            ret_dict['co2']     = temp_arr['CO2_HUPCRS']
            ret_dict['ch4']     = temp_arr['CH4_HUPCRS']
            ret_dict['co']      = temp_arr['CO_HUPCRS']
            ret_dict['temp']    = temp_arr['TEMPERATURE']
            # pressures are given in hPa
            ret_dict['pres']    = temp_arr['PRESSURE'] * 100.0
            if self.assimilate:
                ret_dict['co2_err'] = self.samp_err * np.ones(len(valid_idx), np.float64)
            else:
                ret_dict['co2_err'] = self.unassim_mdm * np.ones(len(valid_idx), np.float64)

        return ret_dict

    def createObservations(self):
        ret_dict = defaultdict(list)
        # What are the files to consider?
        for aircraft in self.flights:
            seac4rs_files = glob(os.path.join(self.input_dir, aircraft, 'SEAC4RS-mrg60-*_merge_*_R1_*.ict'))
            for file_name in seac4rs_files:
                data_dict = self.parseFile(file_name)

                n_obs = len(data_dict['co2']) if 'co2' in data_dict else 0
                if n_obs > 0:
                    data_dict['site_id'] = self.site_id[aircraft] * np.ones(n_obs, np.int32)
                    data_dict['samp_strat'] = self.samp_strat * np.ones(n_obs, np.int16)

                for k,v in data_dict.items():
                    ret_dict[k].extend(v)

        return ret_dict

class ZOTTO_CO2(Observations):

    def __init__(self, *args, **kwargs):
        super(ZOTTO_CO2, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.input_dir = self.rcf.get('ZOTTO.CO2.data.folder')
        self.samp_strat = self.rcf.get('ZOTTO.CO2.sampling.strategy', 'int')
        self.assimilate = self.rcf.get('ZOTTO.CO2.assimilate', 'bool', default=False)
        self.site_id = self.rcf.get('ZOTTO.CO2.site_id', 'int')
        if not self.assimilate:
            self.unassim_mdm = self.rcf.get('CO2.point.unassimilate.mdm', 'float')
        self.ground_level = self.rcf.get('ZOTTO.tower.ground_level', 'float')
        self.zotto_lat = self.rcf.get('ZOTTO.tower.latitude', 'float')
        self.zotto_lon = self.rcf.get('ZOTTO.tower.longitude', 'float')

    def parseRawFile(self, file_name):
        # Because the ZOTTO data are split by local time instead of UTC time, we'll convert everything into a single netcdf
        # file and then read from the netcdf file when we need to. This routine is for doing the first step.
        ret_dict = {}

        temp_arr = np.genfromtxt(file_name, autostrip=True, names=True, delimiter=',', \
            usecols=('UTCdate', 'UTCclock', 'AirLineNo', 'AirLine', 'CO2_ppm', 'CO2_stdev', 'CH4_ppb', 'CH4_stdev', \
            'H2Oraw_ppm', 'H2Oraw_stdev'), dtype=('|S12', '|S12', np.int8, '|S8', np.float32, np.float32, np.float32, \
            np.float32, np.float32, np.float32))

        # Filter out invalid entries
        valid_coords = ['CO2_ppm', 'CH4_ppb', 'CO2_stdev', 'CH4_stdev']
        valid_idx = np.zeros(len(temp_arr), bool)
        valid_idx[:] = True
        for col_name in valid_coords:
            valid_idx = np.logical_and(valid_idx, ~np.isnan(temp_arr[col_name]))
        temp_arr = temp_arr[valid_idx]

        # The UTC time can be in two different formats
        # The construct below assumes that within one file all dates are in the same format. If that's not the case,
        # I'm f**ked. Not really, I'll just have to write a loop :-(
        time_offsets = [[int(x) for x in d.split(':')] for d in temp_arr['UTCclock']]
        time_offsets = [timedelta(hours=h,minutes=m,seconds=s) for h,m,s in time_offsets]
        try:
            utc_times = [datetime.strptime(a, '%Y-%m-%d') + d for a,d in zip(temp_arr['UTCdate'], time_offsets)]
        except ValueError:
            utc_times = [datetime.strptime(a, '%d/%m/%Y') + d for a,d in zip(temp_arr['UTCdate'], time_offsets)]

        ret_dict['times'] = utc_times
        ret_dict['intake_height'] = np.array([int(s[:-1]) for s in temp_arr['AirLine']], np.int16)
        ret_dict['intake_line'] = temp_arr['AirLineNo']
        ret_dict['CO2'] = temp_arr['CO2_ppm']
        ret_dict['CO2_stdev'] = temp_arr['CO2_stdev']
        ret_dict['CH4'] = temp_arr['CH4_ppb']
        ret_dict['CH4_stdev'] = temp_arr['CH4_stdev']
        ret_dict['H2O'] = temp_arr['H2Oraw_ppm'] / 1.0E6
        ret_dict['H2O_stdev'] = temp_arr['H2Oraw_stdev'] / 1.0E6

        return ret_dict

    def convertToNetcdf(self):
        # What are the files that need conversion?
        all_files = glob(os.path.join(self.input_dir, 'ZOT*air.csv'))
        # Sort them by date, which is the same as sorting them by base name
        all_files.sort(key=lambda x: os.path.basename(x))

        # Create an output file
        op_file = self.rcf.get('ZOTTO.consolidated.obs.file')
        with Dataset(op_file, 'w') as fid:
            fid.createDimension('n_obs', 0)
            fid.createDimension('idate', 6)
            setattr(fid, 'ground_elevation', self.ground_level)
            setattr(fid, 'latitude', self.zotto_lat)
            setattr(fid, 'longitude', self.zotto_lon)

            v = fid.createVariable('time_components', np.int16, ('n_obs', 'idate'))

            v = fid.createVariable('intake_height', np.int16, ('n_obs',))
            v.unit = 'Meters above ground level'

            v = fid.createVariable('intake_line', np.int8, ('n_obs',))

            v = fid.createVariable('CO2', np.float32, ('n_obs',))
            v.unit = 'parts per million'

            v = fid.createVariable('CO2_stdev', np.float32, ('n_obs',))
            v.unit = 'parts per million'

            v = fid.createVariable('CH4', np.float32, ('n_obs',))
            v.unit = 'parts per billion'

            v = fid.createVariable('CH4_stdev', np.float32, ('n_obs',))
            v.unit = 'parts per billion'

            v = fid.createVariable('H2O', np.float32, ('n_obs',))
            v.unit = 'mole fraction'

            v = fid.createVariable('H2O_stdev', np.float32, ('n_obs',))
            v.unit = 'mole fraction'

        # Now fill it up
        for file_name in all_files:
            data_dict = self.parseRawFile(file_name)

            with Dataset(op_file, 'a') as fid:
                i_obs = len(fid.dimensions['n_obs'])
                n_obs = len(data_dict['CO2'])

                # The times have to be written separately, everything else can be in a for loop
                fid.variables['time_components'][i_obs:i_obs+n_obs,:] = np.array([d.timetuple()[:6] for d in data_dict['times']], np.int16)

                del data_dict['times']
                for var_name, var_value in data_dict.items():
                    fid.variables[var_name][i_obs:i_obs+n_obs] = var_value

            print '%s added'%os.path.basename(file_name)

    def createObservations(self):
        ret_dict = defaultdict(list)
        cons_file = self.rcf.get('ZOTTO.consolidated.obs.file')
        intakes = self.rcf.get('ZOTTO.CO2.accept.intakes').split()
        intakes = [int(i) for i in intakes]

        with Dataset(cons_file, 'r') as fid:
            all_intakes = fid.variables['intake_line'][:]
            all_times = np.array([datetime(*d) for d in fid.variables['time_components'][:]])

            valid_idx_1 = np.array([i in intakes for i in all_intakes], bool)
            valid_idx_2 = np.logical_and(all_times >= self.StartDate, all_times <= self.EndDate)
            valid_idx = np.logical_and(valid_idx_1, valid_idx_2)
            n_obs = valid_idx.sum()

            ret_dict['date'] = all_times[valid_idx]
            ret_dict['co2'] = fid.variables['CO2'][:][valid_idx]
            if self.assimilate:
                ret_dict['co2_err'] = fid.variables['CO2_stdev'][:][valid_idx]
            else:
                ret_dict['co2_err'] = self.unassim_mdm * np.ones(n_obs, np.float32)
            ret_dict['alt'] = fid.ground_elevation + fid.variables['intake_height'][:][valid_idx]
            ret_dict['ch4'] = fid.variables['CH4'][:][valid_idx]

            ret_dict['lat'] = fid.latitude * np.ones(n_obs, np.float32)
            ret_dict['lon'] = fid.longitude * np.ones(n_obs, np.float32)
            ret_dict['site_id'] = self.site_id * np.ones(n_obs, np.int32)
            ret_dict['samp_strat'] = self.samp_strat * np.ones(n_obs, np.int16)

        return ret_dict

class ATTREX_CO2(Observations):
    """
    This class parses CO2 observations from the ATTREX campaign.
    """
    def __init__(self, *args, **kwargs):
        super(ATTREX_CO2, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.input_dir = self.rcf.get('ATTREX.CO2.data.folder')
        self.samp_strat = self.rcf.get('ATTREX.CO2.sampling.strategy', 'int')
        self.samp_err = self.rcf.get('ATTREX.CO2.error', 'float')
        self.assimilate = self.rcf.get('ATTREX.CO2.assimilate', 'bool', default=False)
        self.site_id = self.rcf.get('ATTREX.CO2.site_id', 'int')
        if not self.assimilate:
            self.unassim_mdm = self.rcf.get('CO2.point.unassimilate.mdm', 'float')

    def parseFile(self, file_name):
        ret_dict = {}

        # ATTREX has different column names for different campaigns/files, joy oh joy!
        if os.path.basename(file_name).split('_')[2] == '2':
            lat_var = 'LN_LAT'
            lon_var = 'LN_LONG'
            alt_var = 'LN_ALT'
        elif os.path.basename(file_name).split('_')[2] == '3':
            lat_var = 'G_LAT'
            lon_var = 'G_LONG'
            alt_var = 'G_ALT'
        else:
            raise ValueError('You wanted to use ATTREX? Fine, code your header names!')

        temp_arr = np.genfromtxt(file_name, autostrip=True, comments='!', names=True, \
            usecols=("YYYYMMDD", "TIME_UTC", "P", "T", lat_var, lon_var, alt_var, "co2_HUPCRS", "co_HUPCRS", "ch4_HUPCRS"), \
            dtype=("|S8", float, float, float, float, float, float, float, float, float), missing_values='NaN')

        # filter out the lines of invalid CO2 data and dates, times, coordinates
        valid_coords = ['TIME_UTC', lat_var, lon_var, alt_var, 'co2_HUPCRS']
        valid_idx = np.zeros(len(temp_arr), bool)
        valid_idx[:] = True
        for col_name in valid_coords:
            valid_idx = np.logical_and(valid_idx, ~np.isnan(temp_arr[col_name]))
        temp_arr = temp_arr[valid_idx]

        # now filter out rows outside the time range of the run
        utc_times = np.array([datetime.strptime(d, '%Y%m%d') + timedelta(seconds=t) for d,t in zip(temp_arr['YYYYMMDD'], temp_arr['TIME_UTC'])])
        valid_idx = np.array([i for i,d in enumerate(utc_times) if self.StartDate <= d <= self.EndDate])

        if len(valid_idx) > 0:
            temp_arr            = temp_arr[valid_idx]
            ret_dict['date']    = utc_times[valid_idx]
            ret_dict['lat']     = temp_arr[lat_var]
            ret_dict['lon']     = temp_arr[lon_var]
            ret_dict['alt']     = temp_arr[alt_var]
            ret_dict['co2']     = temp_arr['co2_HUPCRS']
            ret_dict['ch4']     = temp_arr['ch4_HUPCRS']
            ret_dict['co']      = temp_arr['co_HUPCRS']
            ret_dict['temp']    = temp_arr['T']
            ret_dict['pres']    = temp_arr['P']
            if self.assimilate:
                ret_dict['co2_err'] = self.samp_err * np.ones(len(valid_idx), np.float64)
            else:
                ret_dict['co2_err'] = self.unassim_mdm * np.ones(len(valid_idx), np.float64)

        return ret_dict

    def createObservations(self):
        ret_dict = defaultdict(list)
        # What are the files to consider?
        attrex_files = glob(os.path.join(self.input_dir, '20151008_ATTREX_?_HUPCRS_merge_10s.tbl'))
        for file_name in attrex_files:
            data_dict = self.parseFile(file_name)
            for k,v in data_dict.items():
                ret_dict[k].extend(v)

        n_obs = len(ret_dict['co2'])
        if n_obs > 0:
            ret_dict['site_id'] = self.site_id * np.ones(n_obs, np.int32)
            ret_dict['samp_strat'] = self.samp_strat * np.ones(n_obs, np.int16)

        return ret_dict

class Gatti_CO2(Observations):
    """
    This class parses CO observations from Luciana Gatti's Brazilian profiles
    """
    def __init__(self, *args, **kwargs):
        super(Gatti_CO2, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.inputDir = self.rcf.get('gatti.co2.data.folder')
        self.site_list = self.rcf.get('gatti.co2.profile.sites').split()
        self.samp_strat = self.rcf.get('gatti.co2.sampling.strategy', 'int')
        sampling_err_line = self.rcf.get('gatti.co2.sampling.error').split(';')
        # The sampling error line has the form 'alt1,err1;alt2,err2,...;altN,errN', meaning that samples below alt1 get
        # err1, those between alt1 and alt2 get err2, etc.
        self.sampling_error_height = []
        self.sampling_error_mag = []
        for element in sampling_err_line:
            atoms = element.split(',')
            self.sampling_error_height.append(float(atoms[0]))
            self.sampling_error_mag.append(float(atoms[1]))
        self.sampling_error_height = np.array(self.sampling_error_height)
        self.sampling_error_mag = np.array(self.sampling_error_mag)
        # Should we assimilate Gatti's CO data?
        self.assimilate = self.rcf.get('gatti.co2.assimilate', 'bool', default=True)
        self.site_id_offset = self.rcf.get('gatti.co2.site_id.offset', 'int')
        if not self.assimilate:
            self.unassim_mdm = self.rcf.get('CO2.point.unassimilate.mdm', 'float')
        self.site_id_dict = {}
        self.site_names = {}
        with open(self.rcf.get('gatti.co2.site_names'), 'r') as fid:
            lines = fid.readlines()
        for line in lines:
            site_code, site_name = line.split(':')
            self.site_names[site_code.strip()] = site_name.strip()

    def parseFile(self, file_name):
        # Parse a file such as 'san.mrg'
        spec_col = 15 # This is for CO2, for CO or any other species, this needs to be changed
        with open(file_name, 'r') as fid:
            all_lines = fid.readlines()
        # Filter by quality flag
        all_lines = [l for l in all_lines if l.split()[spec_col+1].startswith('..')]
        # Get the dates
        date_times = [datetime.strptime(''.join(l.split()[1:6]), '%Y%m%d%H%M') for l in all_lines]
        # Get the indices corresponding to samples within the right time window
        idx = np.array([i for i,d in enumerate(date_times) if self.StartDate <= d < self.EndDate])
        ret_dict = defaultdict(list)
        if len(idx) > 0:
            ret_dict['date'] = np.array(date_times)[idx]
            # The coordinates
            ret_dict['lat'] = np.array([float(l.split()[8]) for l in all_lines])[idx]
            ret_dict['lon'] = np.array([float(l.split()[9]) for l in all_lines])[idx]
            ret_dict['alt'] = np.array([float(l.split()[10]) for l in all_lines])[idx]
            # Now the concentrations
            ret_dict['co2'] = np.array([float(l.split()[spec_col]) for l in all_lines])[idx]
            # Construct the MDM
            if self.assimilate:
                ret_dict['co2_error'] = np.zeros_like(ret_dict['co2'])
                for i, a in enumerate(ret_dict['alt']):
                    alt_idx = self.sampling_error_height.searchsorted(a)
                    ret_dict['co2_error'][i] = self.sampling_error_mag[alt_idx]
            else:
                ret_dict['co2_error'] = self.unassim_mdm * np.ones_like(ret_dict['co2'])
            # Also the event ID
            ret_dict['event_id'] = np.array([int(l.split()[14]) for l in all_lines])[idx]
        return ret_dict

    def createObservations(self):
        # The site ID has to be -100, for the first site, -101 for the second site, etc.
        write_dict = defaultdict(list)
        for i_site, site_code in enumerate(self.site_list):
            file_name = os.path.join(self.inputDir, '%s.mrg'%site_code)
            data_dict = self.parseFile(file_name)
            for k,v in data_dict.items():
                write_dict[k].extend(v)
            n_obs = len(data_dict['co2'])
            site_id = (self.site_id_offset - i_site) * np.ones(n_obs, np.int32)
            write_dict['site_id'].extend(list(site_id))

        # Now check if there are duplicate event IDs
        write_dict['event_id'] = np.array(write_dict['event_id'])
        u, i = np.unique(write_dict['event_id'], return_inverse=True)
        duplicates = u[np.bincount(i) > 1]
        if len(duplicates) > 0:
            sys.stderr.write("The following event IDs were found multiple times\n")
            for d in duplicates:
                sys.stderr.write("%i : %i times\n"%(d, np.count_nonzero(write_dict['event_id'] == d)))
            raise

        # Delete the event ID
        del write_dict['event_id']

        return write_dict

class RemoTeC(Observations):

    def GetFileList(self, ymd_tuple):
        if self.split_period == 'm':
            year, month = ymd_tuple
            datelist = [datetime(year,month,day) for day in range(1, monthrange(year,month)[1]+1)]
        elif self.split_period == 'd':
            datelist = [datetime(*ymd_tuple)]
        file_list = []
        for date in datelist:
            fileName = os.path.join(self.input_folder, date.strftime(self.filename_pattern))
            if os.path.isfile(fileName) and fileName not in file_list:
                file_list.append(fileName)
        return file_list

    def FilterData(self, fid):
        """
        Instead of just excluding all medium gain and glint soundings when they are not desired, we will keep them, but set
        their model data mismatches really high. That way, we can later compare how those soundings performed, compared to
        the optimised CO2 field. The only indices to be thrown out are going to be the np.ones QC-flagged by RemoTeC.
        """
        # /REMOTEC-OUTPUT/FLAG_QUAL should not be 'BAD'
        logical_mask = fid.groups['REMOTEC-OUTPUT'].variables['FLAG_QUAL'][:] == 'BAD'
        return np.where(logical_mask == False)[0]

    def ConsolidateData(self, ymd_tuple):
        self.data_dict = defaultdict(list)
        data_files = self.GetFileList(ymd_tuple)
        for file in data_files:
            self.AppendData(file)
        # If there are no valid data points, skip the rest
        if len(self.data_dict['times']) == 0:
            return

        # sort according to times
        sort_order = np.argsort(np.array(self.data_dict['times']))
        for var_name, var_value in self.data_dict.items():
            # sometimes variables like the aerosol size are not needed, and therefore not filled
            if len(var_value) > 0:
                self.data_dict[var_name] = np.array(var_value)[sort_order]
        if self.error_inflate:
            self.InflateErrors()

        # Now depending on whether we want glint and medium gain soundings or not, set the MDMs really high
        if self.exclude_medgain:
            mg_indices = np.where(self.data_dict['gain'] == 1)[0]
            self.data_dict['sigma XCO2'][mg_indices] = self.unassim_mdm

        if self.exclude_sunglint:
            sg_indices = np.where(self.data_dict['sunglint'] == 1)[0]
            self.data_dict['sigma XCO2'][sg_indices] = self.unassim_mdm

    def AppendData(self, fileName):
        data_dict = self.ReadData(fileName)
        for key, value in data_dict.items():
            self.data_dict[key].extend(value)

    def InflateErrors(self):
        """
        For inflating errors in total column measurements, based on some binning length and binning time
        """
        from tm5_utils import sat_utils

        CO2_sigma_column_mixing = self.data_dict['sigma XCO2']
        sample_times = self.data_dict['times']
        Latitudes = self.data_dict['latitude']
        Longitudes = self.data_dict['longitude']
        if len(CO2_sigma_column_mixing) > 0:
            # convert sample times to floating point seconds
            sample_seconds = [t-sample_times[0] for t in sample_times]
            sample_seconds = [dt.days * 86400 + dt.seconds for dt in sample_seconds]
            locations = np.zeros((len(Latitudes),2), Latitudes.dtype)
            locations[:,0] = Latitudes
            locations[:,1] = Longitudes
            self.data_dict['sigma XCO2'], num_samples = sat_utils.errorInflation(float64(sample_seconds), CO2_sigma_column_mixing, locations, \
                self.binning_time, self.binning_length, 'b', 0.0)
            del sample_seconds, locations, num_samples
        else:
            raise ValueError('No samples in this day/month, please investigate!')

class RemoTeC_2_1fluo(RemoTeC, Observations):
    """
    Class for Philippe's experiment with chlorophyll fluorescence. No bias correction is optimized.
    """
    def __init__(self, *args, **kwargs):
        # Since satellite observations are split by day or month, only the beginning time,
        # which specifies the day/month, matters
        super(RemoTeC_2_1fluo, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.instrument = 'satellite.GOSAT.SRON.CO2'
        # self.input_folder should contain a list of folders with input files
        # land and ocean pixels from the same file, for example input/GOSAT/SRON/V1.99/2010/09/2010_09_16_RemoTeC-GOSAT.nc
        self.input_folder = self.rcf.get(self.instrument + '.data.folder')
        self.filename_pattern = os.path.join(self.input_folder, self.rcf.get(self.instrument + '.inputfile.pattern'))
        self.exclude_medgain = self.rcf.get(self.instrument + '.exclude.mediumgain', 'bool')
        self.exclude_sunglint = self.rcf.get(self.instrument + '.exclude.sunglint', 'bool')
        self.error_inflate = self.rcf.get(self.instrument + '.errors.inflate', 'bool')
        if self.error_inflate:
            self.binning_time = self.rcf.get(self.instrument + '.binning.time', 'float')
            self.binning_length = self.rcf.get(self.instrument + '.binning.length', 'float')
        self.posterior_filters = self.rcf.get(self.instrument + '.posterior.filters', 'bool')
        self.xco2_var_name = self.rcf.get(self.instrument + '.varname')
        self.split_period = self.rcf.get('output.satellite.split.period')
        # Should we assimilate CO2 satellite data?
        self.assimilate = self.rcf.get('CO2.satellite.assimilate', 'bool', default=True)
        self.unassim_mdm = self.rcf.get('CO2.satellite.unassimilate.mdm', 'float')
        # Has fluorescence been retrieved?
        self.fluo_retrieved = self.rcf.get('remotec.CO2.include.fluorescence', 'bool')

    def ReadData(self, fileName):
        fid = Dataset(fileName, 'r')
        # filter first to get the necessary indices
        if self.posterior_filters:
            valid_indices = self.FilterData(fid)
        else:
            valid_indices = np.arange(len(fid.dimensions['time']))
        sample_time = []
        gdict = fid.groups['GOSAT'].variables
        sample_time = [datetime(int(y),int(m),int(d),int(H),int(M))+timedelta(seconds=S) for y,m,d,H,M,S in zip(gdict['TIME_YEAR'][:][valid_indices], gdict['TIME_MONTH'][:][valid_indices], gdict['TIME_DAY'][:][valid_indices], gdict['TIME_HOUR'][:][valid_indices], gdict['TIME_MIN'][:][valid_indices], gdict['TIME_SEC'][:][valid_indices])]
        nlev = len(fid.dimensions['nlev'])
        surface_pressure = 100.0 * fid.groups['METEO'].variables['PRESSURE'][:,nlev-1][valid_indices]
        # pressure_levels are mid-level pressures
        pressure_levels = 100.0 * 0.5 * (fid.groups['METEO'].variables['PRESSURE'][:,1:][valid_indices] + fid.groups['METEO'].variables['PRESSURE'][:,:-1][valid_indices])
        D_Xco2 = fid.groups['REMOTEC-OUTPUT'].variables['X_COLUMN_ERR'][:,0][valid_indices]
        AK = fid.groups['REMOTEC-OUTPUT'].variables['X_AK_COLUMN'][:,0,:][valid_indices]
        Prior = fid.groups['REMOTEC-OUTPUT'].variables['X_APR_PROF'][:,0,:][valid_indices]
        lat = fid.groups['GOSAT'].variables['LATITUDE'][:][valid_indices]
        lon = fid.groups['GOSAT'].variables['LONGITUDE'][:][valid_indices]
        glint = fid.groups['GOSAT'].variables['FLAG_SUNGLINT_CUSTOM'][:][valid_indices]
        gain = np.where(fid.groups['GOSAT'].variables['GAIN'][:][valid_indices] == 'HH', 0, 1) # 0 for high gain, 1 for medium gain
        Xco2 = fid.groups['REMOTEC-OUTPUT'].variables[self.xco2_var_name][:,0][valid_indices]
        unique_id_l1b = fid.groups['GOSAT'].variables['L1B_NAME'][:][valid_indices]
        unique_id_acos = fid.groups['GOSAT'].variables['IDACOS'][:][valid_indices]
        unique_id = np.array(zip(unique_id_l1b, unique_id_acos))
        fid.close()
        return_dict = {'times': sample_time, 'latitude': lat, 'longitude': lon, 'P surface': surface_pressure, 'P levels': pressure_levels, 'Avg kernel': AK,\
            'XCO2': Xco2, 'sigma XCO2': D_Xco2, 'prior profile': Prior, 'sunglint': glint, 'ID': unique_id, 'gain': gain}
        return return_dict

    def __call__(self):
        self.rejection_dict = defaultdict(int)
        if self.split_period == 'm':
            ymd_tuple = self.StartDate.timetuple()[:2]
        elif self.split_period == 'd':
            ymd_tuple = self.StartDate.timetuple()[:3]
        else:
            sys.stderr.write("Invalid split period specified for satellite data\n")
            raise
        self.ConsolidateData(ymd_tuple)
        n_obs = len(self.data_dict['times'])
        # maintain a dictionary of dimensions to write and variable shapes
        self.SatObservation[self.tracer]['dimensions'] = {'n_obs': n_obs}
        if n_obs == 0:
            return

        self.SatObservation[self.tracer]['dimensions']['n_levels'] = len(self.data_dict['prior profile'][0])
        self.SatObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.SatObservation[self.tracer]['variable_attrs'] = OrderedDict()
        self.SatObservation[self.tracer]['file_attrs'] = {}

        self.SatObservation[self.tracer]['cdate'] = np.array([d.timetuple()[:6] for d in self.data_dict['times']], np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['cdate'] = ('n_obs', 'idate')
        self.SatObservation[self.tracer]['variable_attrs']['cdate'] = [('description', 'time to be sampled')]

        samp_strat = self.rcf.get(self.instrument + '.sampling.strategy', 'int')
        self.SatObservation[self.tracer]['sampling_strategy'] = np.int16(samp_strat * np.ones(n_obs, np.int16))
        self.SatObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sampling_strategy'] = \
            [('description', 'sampling strategy, 3 => symmetric sampling, 2=> ndyn sampling, 1 => instantaneous sampling')]

        self.SatObservation[self.tracer]['p_surf'] = np.float32(self.data_dict['P surface'])
        self.SatObservation[self.tracer]['variable_shapes']['p_surf'] = ('n_obs',)

        self.SatObservation[self.tracer]['p_levels'] = np.float32(self.data_dict['P levels'])
        self.SatObservation[self.tracer]['variable_shapes']['p_levels'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['column_mixing'] = np.float64(self.data_dict['XCO2'])
        self.SatObservation[self.tracer]['variable_shapes']['column_mixing'] = ('n_obs',)

        if self.assimilate:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = np.float64(self.data_dict['sigma XCO2'])
        else:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('comment', 'These data will not be assimilated')]
        self.SatObservation[self.tracer]['variable_shapes']['sigma_column_mixing'] = ('n_obs',)

        self.SatObservation[self.tracer]['avg_kernel'] = np.float32(self.data_dict['Avg kernel'])
        self.SatObservation[self.tracer]['variable_shapes']['avg_kernel'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['prior_mixing'] = np.float32(self.data_dict['prior profile'])
        self.SatObservation[self.tracer]['variable_shapes']['prior_mixing'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['latitude'] = np.float32(self.data_dict['latitude'])
        self.SatObservation[self.tracer]['variable_shapes']['latitude'] = ('n_obs',)

        self.SatObservation[self.tracer]['longitude'] = np.float32(self.data_dict['longitude'])
        self.SatObservation[self.tracer]['variable_shapes']['longitude'] = ('n_obs',)

        self.SatObservation[self.tracer]['glint'] = np.int8(self.data_dict['sunglint'])
        self.SatObservation[self.tracer]['variable_shapes']['glint'] = ('n_obs',)

        self.SatObservation[self.tracer]['gain'] = np.int8(self.data_dict['gain'])
        self.SatObservation[self.tracer]['variable_shapes']['gain'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['gain'] = \
            [('legend', '0 for high gain, 1 for medium gain')]

        self.SatObservation[self.tracer]['acos_id'] = np.array(self.data_dict['ID'][:,1], int64)
        self.SatObservation[self.tracer]['variable_shapes']['acos_id'] = ('n_obs',)

class RemoTeC_2_1fluo_BC(RemoTeC_2_1fluo):
    """
    Class for Philippe's experiment with chlorophyll fluorescence. Optimize bias correction. For the retrieval without fluorescence,

    XCO2(bc) = XCO2(raw) * (1.0057425 - albedo_O2 * 0.01 + SOT * 0.015 + 0.015/alpha)

    or, in terms of variables stored in the netcdf files,

    X_COLUMN_BIASCOR(:,TAR0002) = X_COLUMN(:,TAR0002) * (1.0057425 - X_STATE(:,ALB_WIN01_ORDER00) * 0.01 + \
        SOT_WIN(:,WIN01) * 0.015 + 0.015/X_STATE(:,AER01SIZE1))

    For the retrieval with fluorescence,

    XCO2(bc) = XCO2(raw) * (1.005185 + albedo_O2 * 0.015 + SOT * 0.015 - height * 2.0E-7)

    or, in terms of variables stored in the netcdf files,

    X_COLUMN_BIASCOR(:,TAR0002) = X_COLUMN(:,TAR0002) * (1.005185 + X_STATE(:,ALB_WIN01_ORDER00) * 0.015 + \
        SOT_WIN(:,WIN01) * 0.015 - X_STATE(:,AER01HEIGHT1) * 2.0E-7)

    For the ocean glint,

    XCO2(bc) = XCO2(raw) * 1.00857

    or, in terms of variables stored in the netcdf files,

    X_COLUMN_BIASCOR(:,TAR0002) = X_COLUMN(:,TAR0002) * 1.00857
    """
    def ReadData(self, fileName):
        return_dict = {}
        with Dataset(fileName, 'r') as fid:
            # filter first to get the necessary indices
            if self.posterior_filters:
                valid_indices = self.FilterData(fid)
            else:
                valid_indices = np.arange(len(fid.dimensions['time']))

            # For this version of RemoTeC, all soundings should either be land H gain or sunglint. Check if there are any medium
            # gain soundings, and if there are, take them out.
            gdict = fid.groups['GOSAT'].variables
            glint = gdict['FLAG_SUNGLINT_CUSTOM'][:][valid_indices]

            gain  = np.where(gdict['GAIN'][:][valid_indices] == 'HH', 0, 1) # 0 for high gain, 1 for medium gain
            # Soundings that are not sunglint should have gain = 0
            land_soundings = np.where(glint == 0)[0]
            if any(gain[land_soundings] == 1):
                med_gain_indices = np.array([i for i,(ga,gl) in enumerate(zip(gain,glint)) if ga == 1 and gl == 0])
                gain = delete(gain, med_gain_indices)
                glint = delete(glint, med_gain_indices)
                valid_indices = delete(valid_indices, med_gain_indices)
                sys.stderr.write("WARNING: %i medium gain land pixels found in %s\n"%(len(med_gain_indices), fileName))

            return_dict['times'] = [datetime(int(y),int(m),int(d),int(H),int(M))+timedelta(seconds=S) for y,m,d,H,M,S in \
                zip(gdict['TIME_YEAR'][:][valid_indices], gdict['TIME_MONTH'][:][valid_indices], \
                gdict['TIME_DAY'][:][valid_indices], gdict['TIME_HOUR'][:][valid_indices], gdict['TIME_MIN'][:][valid_indices], \
                gdict['TIME_SEC'][:][valid_indices])]
            return_dict['latitude'] = gdict['LATITUDE'][:][valid_indices]
            return_dict['longitude'] = gdict['LONGITUDE'][:][valid_indices]
            return_dict['sunglint'] = gdict['FLAG_SUNGLINT_CUSTOM'][:][valid_indices]
            return_dict['gain'] = np.where(gdict['GAIN'][:][valid_indices] == 'HH', 0, 1) # 0 for high gain, 1 for medium gain
            return_dict['ID'] = gdict['IDACOS'][:][valid_indices]

            nlev = len(fid.dimensions['nlev'])
            gdict = fid.groups['METEO'].variables
            return_dict['P surface'] = 100.0 * gdict['PRESSURE'][:,nlev-1][valid_indices]
            # pressure_levels are mid-level pressures
            return_dict['P levels'] = 100.0 * 0.5 * (gdict['PRESSURE'][:,1:][valid_indices] + gdict['PRESSURE'][:,:-1][valid_indices])

            gdict = fid.groups['REMOTEC-OUTPUT'].variables
            ico2 = 0
            return_dict['sigma XCO2'] = gdict['X_COLUMN_ERR'][:,ico2][valid_indices]
            return_dict['Avg kernel'] = gdict['X_AK_COLUMN'][:,ico2,:][valid_indices]
            return_dict['prior profile'] = gdict['X_APR_PROF'][:,ico2,:][valid_indices]
            return_dict['XCO2'] = gdict[self.xco2_var_name][:,ico2][valid_indices] # this should be X_COLUMN, not X_COLUMN_BIASCOR

            # read some more variables for the bias correction
            dim_names = [s.strip() for s in gdict['X_STATE'].dim1.split(',')]
            var_idx = dim_names.index('ALB_WIN01_ORDER00')
            return_dict['O2 albedo'] = gdict['X_STATE'][:,var_idx][valid_indices]
            var_idx = dim_names.index('AER01SIZE1')
            return_dict['aerosol size'] = gdict['X_STATE'][:,var_idx][valid_indices]
            var_idx = dim_names.index('AER01HEIGHT1')
            return_dict['aerosol height'] = gdict['X_STATE'][:,var_idx][valid_indices]
            dim_names = [s.strip() for s in gdict['SOT_WIN'].dim1.split(',')]
            var_idx = dim_names.index('WIN01')
            return_dict['SOT'] = gdict['SOT_WIN'][:,var_idx][valid_indices]
            return_dict['XCO2 bias corrected'] = gdict['X_COLUMN_BIASCOR'][:,ico2][valid_indices]

        # Also store the difference of XCO2 bias corrected and XCO2 raw with bias correction, as a debug variable which should be zero
        Del_Xco2_bc = empty_like(return_dict['XCO2'])
        # What are the land H gain indices?
        idx = np.where(return_dict['sunglint'] == 0)[0]
        if self.fluo_retrieved:
            # XCO2(bc) = XCO2(raw) * (1.005185 + albedo_O2 * 0.015 + SOT * 0.015 - aer_height * 2.0E-7)
            Del_Xco2_bc[idx] = return_dict['XCO2 bias corrected'][idx] - return_dict['XCO2'][idx] * (1.005185 + \
                return_dict['O2 albedo'][idx] * 0.015 + return_dict['SOT'][idx] * 0.015 - return_dict['aerosol height'][idx] * 2.0E-7)
        else:
            # XCO2(bc) = XCO2(raw) * (1.0057425 - albedo_O2 * 0.01 + SOT * 0.015 + 0.015/aer_size)
            Del_Xco2_bc[idx] = return_dict['XCO2 bias corrected'][idx] - return_dict['XCO2'][idx] * (1.0057425 - \
                return_dict['O2 albedo'][idx] * 0.01 + return_dict['SOT'][idx] * 0.015 + 0.015/return_dict['aerosol size'][idx])
        # Now the ocean indices
        idx = np.where(return_dict['sunglint'] == 1)[0]
        # XCO2(bc) = XCO2(raw) * 1.00857
        Del_Xco2_bc[idx] = return_dict['XCO2 bias corrected'][idx] - return_dict['XCO2'][idx] * 1.00857

        return_dict['Del XCO2'] = Del_Xco2_bc

        return return_dict

    def __call__(self):
        if self.split_period == 'm':
            ymd_tuple = self.StartDate.timetuple()[:2]
        elif self.split_period == 'd':
            ymd_tuple = self.StartDate.timetuple()[:3]
        else:
            sys.stderr.write("Invalid split period specified for satellite data\n")
            raise
        self.ConsolidateData(ymd_tuple)
        n_obs = len(self.data_dict['times'])
        # maintain a dictionary of dimensions to write and variable shapes
        self.SatObservation[self.tracer]['dimensions'] = {'n_obs': n_obs}
        if n_obs == 0:
            return

        self.SatObservation[self.tracer]['dimensions']['n_levels'] = len(self.data_dict['prior profile'][0])
        self.SatObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.SatObservation[self.tracer]['variable_attrs'] = OrderedDict()
        self.SatObservation[self.tracer]['file_attrs'] = {}

        # create the bias correction formula and store it in the file
        if self.fluo_retrieved:
            biascorr_formula = 'XCO2(bc) = XCO2(raw) * (1.005185 + albedo_O2 * 0.015 + SOT * 0.015 - aer_height * 2.0E-7)'
        else:
            biascorr_formula = 'XCO2(bc) = XCO2(raw) * (1.0057425 - albedo_O2 * 0.01 + SOT * 0.015 + 0.015/aer_size)'
        self.SatObservation[self.tracer]['file_attrs']['land_H_biascorr'] = biascorr_formula
        biascorr_formula = 'XCO2(bc) = XCO2(raw) * 1.00857'
        self.SatObservation[self.tracer]['file_attrs']['glint_biascorr'] = biascorr_formula

        self.SatObservation[self.tracer]['cdate'] = np.array([d.timetuple()[:6] for d in self.data_dict['times']], np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['cdate'] = ('n_obs', 'idate')
        self.SatObservation[self.tracer]['variable_attrs']['cdate'] = [('description', 'time to be sampled')]

        samp_strat = self.rcf.get(self.instrument + '.sampling.strategy', 'int')
        self.SatObservation[self.tracer]['sampling_strategy'] = np.int16(samp_strat * np.ones(n_obs, np.int16))
        self.SatObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sampling_strategy'] = \
            [('description', 'sampling strategy, 3 => symmetric sampling, 2=> ndyn sampling, 1 => instantaneous sampling')]

        self.SatObservation[self.tracer]['p_surf'] = np.float32(self.data_dict['P surface'])
        self.SatObservation[self.tracer]['variable_shapes']['p_surf'] = ('n_obs',)

        self.SatObservation[self.tracer]['p_levels'] = np.float32(self.data_dict['P levels'])
        self.SatObservation[self.tracer]['variable_shapes']['p_levels'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['column_mixing'] = np.float64(self.data_dict['XCO2'])
        self.SatObservation[self.tracer]['variable_shapes']['column_mixing'] = ('n_obs',)

        if self.assimilate:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = np.float64(self.data_dict['sigma XCO2'])
        else:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('comment', 'These data will not be assimilated')]
        self.SatObservation[self.tracer]['variable_shapes']['sigma_column_mixing'] = ('n_obs',)

        self.SatObservation[self.tracer]['avg_kernel'] = np.float32(self.data_dict['Avg kernel'])
        self.SatObservation[self.tracer]['variable_shapes']['avg_kernel'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['prior_mixing'] = np.float32(self.data_dict['prior profile'])
        self.SatObservation[self.tracer]['variable_shapes']['prior_mixing'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['latitude'] = np.float32(self.data_dict['latitude'])
        self.SatObservation[self.tracer]['variable_shapes']['latitude'] = ('n_obs',)

        self.SatObservation[self.tracer]['longitude'] = np.float32(self.data_dict['longitude'])
        self.SatObservation[self.tracer]['variable_shapes']['longitude'] = ('n_obs',)

        self.SatObservation[self.tracer]['glint'] = np.int8(self.data_dict['sunglint'])
        self.SatObservation[self.tracer]['variable_shapes']['glint'] = ('n_obs',)

        self.SatObservation[self.tracer]['gain'] = np.int8(self.data_dict['gain'])
        self.SatObservation[self.tracer]['variable_shapes']['gain'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['gain'] = \
            [('legend', '0 for high gain, 1 for medium gain')]

        self.SatObservation[self.tracer]['acos_id'] = int64(self.data_dict['ID'])
        self.SatObservation[self.tracer]['variable_shapes']['acos_id'] = ('n_obs',)

        # now store some variables relating to bias correction
        self.SatObservation[self.tracer]['XCO2_bc'] = np.float64(self.data_dict['XCO2 bias corrected'])
        self.SatObservation[self.tracer]['variable_shapes']['XCO2_bc'] = ('n_obs',)

        self.SatObservation[self.tracer]['albedo_O2'] = np.float64(self.data_dict['O2 albedo'])
        self.SatObservation[self.tracer]['variable_shapes']['albedo_O2'] = ('n_obs',)

        self.SatObservation[self.tracer]['SOT'] = np.float64(self.data_dict['SOT'])
        self.SatObservation[self.tracer]['variable_shapes']['SOT'] = ('n_obs',)

        self.SatObservation[self.tracer]['aer_height'] = np.float64(self.data_dict['aerosol height'])
        self.SatObservation[self.tracer]['variable_shapes']['aer_height'] = ('n_obs',)

        self.SatObservation[self.tracer]['aer_size'] = np.float64(self.data_dict['aerosol size'])
        self.SatObservation[self.tracer]['variable_shapes']['aer_size'] = ('n_obs',)

        self.SatObservation[self.tracer]['del_XCO2_BC_debug'] = np.float64(self.data_dict['Del XCO2'])
        self.SatObservation[self.tracer]['variable_shapes']['del_XCO2_BC_debug'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['del_XCO2_BC_debug'] = [('description', \
            'difference between bias corrected XCO2 and raw XCO2 with RemoTeC bias correction applied, should be np.zeros'), \
            ('range', '%.5g, %.5g'%(self.data_dict['Del XCO2'].min(), self.data_dict['Del XCO2'].max()))]

        # Write out the initial bias correction parameters. This construction will only work for the case when the only
        # parameters we are optimizing are CO2/GOSAT related. If there are other parameters to be optimized as well, we
        # need to move this parameter apri and correlation writing business to somewhere more general. The CO2/GOSAT
        # part of that will be filled by this routine, and other parts will be filled by other routines.
        file_name = self.rcf.get('parameter.apri.filename')
        param_apri = np.array([float(s) for s in self.rcf.get('CO2.sat.bias.init_params').split()])
        if not os.path.isdir(os.path.dirname(file_name)):
            os.makedirs(os.path.dirname(file_name))
        savetxt(file_name, param_apri)

        file_name = self.rcf.get('parameter.prior.error.filename')
        param_dapri = np.array([float(s) for s in self.rcf.get('CO2.sat.bias.init_param.errs').split()])
        if not os.path.isdir(os.path.dirname(file_name)):
            os.makedirs(os.path.dirname(file_name))
        savetxt(file_name, param_dapri)

        file_name = self.rcf.get('parameter.correlation.filename')
        num_params = self.rcf.get('CO2.sat.bias.num_params', 'int')
        corr_matrix = diag(ones(num_params, dtype=float64))
        if not os.path.isdir(os.path.dirname(file_name)):
            os.makedirs(os.path.dirname(file_name))
        savetxt(file_name, corr_matrix, fmt='%10.5f')

class RemoTeC_2_11(Observations):
    """
    It seems that we will need to change the variables written to the TM5 satellite input
    files very frequently. Therefore, I'm rewriting the satellite observation class to be
    more flexible.

    We might decide to apply our inversion-optimized bias correction to the XCO2 before
    writing it in files. In that case, the form of the bias correction for V2.11 is:

    land M : X' = X*(b[0] + b[1]/size_parameter) # prior = 0.98298, 0.07616
    land H : X' = X*(b[2] + b[3]/size_parameter) # prior = 0.98747, 0.07586
    ocean  : X' = X*b[4] # prior = 1.00857
    """
    def __init__(self, *args, **kwargs):
        # Since satellite observations are split by day or month, only the beginning time,
        # which specifies the day/month, matters
        super(RemoTeC_2_11, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.instrument = 'satellite.GOSAT.SRON.CO2'
        # self.input_folder should contain a list of folders with input files
        # land and ocean pixels from the same file, for example input/GOSAT/SRON/V1.99/2010/09/2010_09_16_RemoTeC-GOSAT.nc
        self.input_folder = self.rcf.get(self.instrument + '.data.folder')
        self.filename_pattern = os.path.join(self.input_folder, self.rcf.get(self.instrument + '.inputfile.pattern'))
        self.exclude_medgain = self.rcf.get(self.instrument + '.exclude.mediumgain', 'bool')
        self.exclude_sunglint = self.rcf.get(self.instrument + '.exclude.sunglint', 'bool')
        self.error_inflate = self.rcf.get(self.instrument + '.errors.inflate', 'bool')
        if self.error_inflate:
            self.binning_time = self.rcf.get(self.instrument + '.binning.time', 'float')
            self.binning_length = self.rcf.get(self.instrument + '.binning.length', 'float')
        self.posterior_filters = self.rcf.get(self.instrument + '.posterior.filters', 'bool')
        self.apply_biascorr = self.rcf.get(self.instrument + '.bias.correct', 'bool', default=False)
        self.apply_optim_bias = self.rcf.get(self.instrument + '.apply.optim.bias', 'bool', default=False)
        if self.apply_optim_bias:
            # look at the comment on top of __init__ for the form of this bias correction.
            self.land_M_bias_1 = self.rcf.get('satellite.GOSAT.SRON.CO2.bias.land_M.1', 'float') # b[0]
            self.land_M_bias_2 = self.rcf.get('satellite.GOSAT.SRON.CO2.bias.land_M.2', 'float') # b[1]
            self.land_H_bias_1 = self.rcf.get('satellite.GOSAT.SRON.CO2.bias.land_H.1', 'float') # b[2]
            self.land_H_bias_2 = self.rcf.get('satellite.GOSAT.SRON.CO2.bias.land_H.2', 'float') # b[3]
            self.ocean_bias = self.rcf.get('satellite.GOSAT.SRON.CO2.bias.ocean', 'float') # b[4]
        self.split_period = self.rcf.get('output.satellite.split.period')
        # Should we assimilate CO2 satellite data?
        self.assimilate = self.rcf.get('CO2.satellite.assimilate', 'bool', default=True)
        self.unassim_mdm = self.rcf.get('CO2.satellite.unassimilate.mdm', 'float')

    def ReadData(self, fileName):
        fid = Dataset(fileName, 'r')
        # filter first to get the necessary indices
        if self.posterior_filters:
            valid_indices = self.FilterData(fid)
        else:
            valid_indices = np.arange(len(fid.dimensions['time']))
        sample_time = []
        gdict = fid.groups['GOSAT'].variables
        sample_time = [datetime(int(y),int(m),int(d),int(H),int(M))+timedelta(seconds=S) for y,m,d,H,M,S in zip(gdict['TIME_YEAR'][:][valid_indices], gdict['TIME_MONTH'][:][valid_indices], gdict['TIME_DAY'][:][valid_indices], gdict['TIME_HOUR'][:][valid_indices], gdict['TIME_MIN'][:][valid_indices], gdict['TIME_SEC'][:][valid_indices])]
        nlev = len(fid.dimensions['nlev'])
        surface_pressure = 100.0 * fid.groups['METEO'].variables['PRESSURE'][:,nlev-1][valid_indices]
        # pressure_levels are mid-level pressures
        pressure_levels = 100.0 * 0.5 * (fid.groups['METEO'].variables['PRESSURE'][:,1:][valid_indices] + fid.groups['METEO'].variables['PRESSURE'][:,:-1][valid_indices])
        D_Xco2 = fid.groups['REMOTEC-OUTPUT'].variables['X_COLUMN_ERR'][:,0][valid_indices]
        AK = fid.groups['REMOTEC-OUTPUT'].variables['X_AK_COLUMN'][:,0,:][valid_indices]
        Prior = fid.groups['REMOTEC-OUTPUT'].variables['X_APR_PROF'][:,0,:][valid_indices]
        lat = fid.groups['GOSAT'].variables['LATITUDE'][:][valid_indices]
        lon = fid.groups['GOSAT'].variables['LONGITUDE'][:][valid_indices]
        glint = fid.groups['GOSAT'].variables['FLAG_SUNGLINT_CUSTOM'][:][valid_indices]
        gain = np.where(fid.groups['GOSAT'].variables['GAIN'][:][valid_indices] == 'HH', 0, 1) # 0 for high gain, 1 for medium gain
        size_param = fid.groups['REMOTEC-OUTPUT'].variables['X_STATE'][:,49][valid_indices] # For optimizing the bias correction
        # From version 2.0 onward, the bias corrected product is already in the netcdf files
        if self.apply_biascorr:
            Xco2 = fid.groups['REMOTEC-OUTPUT'].variables['X_COLUMN_BIASCOR'][:,0][valid_indices]
        else:
            Xco2 = fid.groups['REMOTEC-OUTPUT'].variables['X_COLUMN'][:,0][valid_indices]
        if self.apply_optim_bias:
            # need to separate land H, land M and ocean pixels
            ocean_indices = np.where(glint == 1)[0]
            land_indices = np.where(glint == 0)[0]
            land_H_indices = np.array([i for i in land_indices if gain[i] == 0])
            land_M_indices = np.array([i for i in land_indices if gain[i] == 1])
            if len(ocean_indices) > 0:
                Xco2[ocean_indices] = Xco2[ocean_indices] * self.ocean_bias
            if len(land_M_indices) > 0:
                Xco2[land_M_indices] = Xco2[land_M_indices] * (self.land_M_bias_1 + self.land_M_bias_2/size_param[land_M_indices])
            if len(land_H_indices) > 0:
                Xco2[land_H_indices] = Xco2[land_H_indices] * (self.land_H_bias_1 + self.land_H_bias_2/size_param[land_H_indices])
        unique_id_l1b = fid.groups['GOSAT'].variables['L1B_NAME'][:][valid_indices]
        unique_id_acos = fid.groups['GOSAT'].variables['IDACOS'][:][valid_indices]
        unique_id = np.array(zip(unique_id_l1b, unique_id_acos))
        fid.close()
        return_dict = {'times': sample_time, 'latitude': lat, 'longitude': lon, 'P surface': surface_pressure, 'P levels': pressure_levels, 'Avg kernel': AK,\
            'XCO2': Xco2, 'sigma XCO2': D_Xco2, 'prior profile': Prior, 'sunglint': glint, 'ID': unique_id, 'aerosol size': size_param, 'gain': gain}
        return return_dict

    def __call__(self):
        self.rejection_dict = defaultdict(int)
        if self.split_period == 'm':
            ymd_tuple = self.StartDate.timetuple()[:2]
        elif self.split_period == 'd':
            ymd_tuple = self.StartDate.timetuple()[:3]
        else:
            sys.stderr.write("Invalid split period specified for satellite data\n")
            raise
        self.ConsolidateData(ymd_tuple)
        n_obs = len(self.data_dict['times'])
        # maintain a dictionary of dimensions to write and variable shapes
        self.SatObservation[self.tracer]['dimensions'] = {'n_obs': n_obs, 'n_levels': len(self.data_dict['prior profile'][0])}
        self.SatObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.SatObservation[self.tracer]['variable_attrs'] = OrderedDict()
        self.SatObservation[self.tracer]['file_attrs'] = {}

        self.SatObservation[self.tracer]['cdate'] = np.array([d.timetuple()[:6] for d in self.data_dict['times']], np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['cdate'] = ('n_obs', 'idate')
        self.SatObservation[self.tracer]['variable_attrs']['cdate'] = [('description', 'time to be sampled')]

        samp_strat = self.rcf.get(self.instrument + '.sampling.strategy', 'int')
        self.SatObservation[self.tracer]['sampling_strategy'] = np.int16(samp_strat * np.ones(n_obs, np.int16))
        self.SatObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sampling_strategy'] = \
            [('description', 'sampling strategy, 3 => symmetric sampling, 2=> ndyn sampling, 1 => instantaneous sampling')]

        self.SatObservation[self.tracer]['p_surf'] = np.float32(self.data_dict['P surface'])
        self.SatObservation[self.tracer]['variable_shapes']['p_surf'] = ('n_obs',)

        self.SatObservation[self.tracer]['p_levels'] = np.float32(self.data_dict['P levels'])
        self.SatObservation[self.tracer]['variable_shapes']['p_levels'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['column_mixing'] = np.float64(self.data_dict['XCO2'])
        self.SatObservation[self.tracer]['variable_shapes']['column_mixing'] = ('n_obs',)
        if self.apply_biascorr:
            self.SatObservation[self.tracer]['variable_attrs']['column_mixing'] = [('comment', 'bias corrected')]
        else:
            self.SatObservation[self.tracer]['variable_attrs']['column_mixing'] = [('comment', 'not bias corrected')]

        if self.assimilate:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = np.float64(self.data_dict['sigma XCO2'])
        else:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('comment', 'These data will not be assimilated')]
        self.SatObservation[self.tracer]['variable_shapes']['sigma_column_mixing'] = ('n_obs',)

        self.SatObservation[self.tracer]['avg_kernel'] = np.float32(self.data_dict['Avg kernel'])
        self.SatObservation[self.tracer]['variable_shapes']['avg_kernel'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['prior_mixing'] = np.float32(self.data_dict['prior profile'])
        self.SatObservation[self.tracer]['variable_shapes']['prior_mixing'] = ('n_obs', 'n_levels')

        self.SatObservation[self.tracer]['latitude'] = np.float32(self.data_dict['latitude'])
        self.SatObservation[self.tracer]['variable_shapes']['latitude'] = ('n_obs',)

        self.SatObservation[self.tracer]['longitude'] = np.float32(self.data_dict['longitude'])
        self.SatObservation[self.tracer]['variable_shapes']['longitude'] = ('n_obs',)

        self.SatObservation[self.tracer]['glint'] = np.int8(self.data_dict['sunglint'])
        self.SatObservation[self.tracer]['variable_shapes']['glint'] = ('n_obs',)

        self.SatObservation[self.tracer]['aerosol_size'] = np.float32(self.data_dict['aerosol size'])
        self.SatObservation[self.tracer]['variable_shapes']['aerosol_size'] = ('n_obs',)

        self.SatObservation[self.tracer]['gain'] = np.int8(self.data_dict['gain'])
        self.SatObservation[self.tracer]['variable_shapes']['gain'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['gain'] = \
            [('legend', '0 for high gain, 1 for medium gain')]
        # 0 for high gain, 1 for medium gain

        self.SatObservation[self.tracer]['acos_id'] = np.array(self.data_dict['ID'][:,1], int64)
        self.SatObservation[self.tracer]['variable_shapes']['acos_id'] = ('n_obs',)

        # Now also write the bias correction parameter files in the input folder, but use different file names,
        # in case there are more parameters to optimize, such as the initial concentration, etc.
        if not self.apply_biascorr:
            num_params = 5
            corr_file_name = self.rcf.get('parameter.correlation.filename')
            corr_file_name_base = 'remotec-2.11-' + os.path.basename(corr_file_name)
            corr_file_name = os.path.join(os.path.dirname(corr_file_name), corr_file_name_base)
            checkDir(corr_file_name)
            savetxt(corr_file_name, eye(num_params), fmt='%5.1f')

            apri_param_std = self.rcf.get('parameter.prior.error.filename')
            apri_param_std_base = 'remotec-2.11-' + os.path.basename(apri_param_std)
            apri_param_std = os.path.join(os.path.dirname(apri_param_std), apri_param_std_base)
            checkDir(apri_param_std)
            savetxt(apri_param_std, 0.03*ones(5, np.float32), fmt='%10.4f')

            apri_param = self.rcf.get('parameter.apri.filename')
            apri_param_base = 'remotec-2.11-' + os.path.basename(apri_param)
            apri_param = os.path.join(os.path.dirname(apri_param), apri_param_base)
            param_values = np.array([0.98298, 0.07616, 0.98747, 0.07586, 1.00857], np.float32)
            checkDir(apri_param)
            savetxt(apri_param, param_values, fmt='%10.5f')

class NIES_aircraft(Observations):
    """
    A class for preparing NIES and JMA aircrafts over Siberia for TM5. There are two sources of data,
    (a) Motoki Sasakawa, and (b) Toshinobu Machida.

    Machida has CO and CO2, at three locations, from aircraft flights. The files are in

    input/Siberia/Machida/CO2 and input/Siberia/Machida/CO

    one file per tower, although the CO2 files seem to contain CO as well.

    Motoki has tower and aircraft data. The aircraft data are in

    input/Siberia/Sasakawa/aircraft/CO2/in-situ/AMES

    and the tower data are in

    input/Siberia/Sasakawa/towers/CO2/AMES

    In fact, it seems that Machida's CO2 data are a subset of Motoki's, although Motoki doesn't have CO.

    There are two kinds of towers, those with two inputs (AZV, DEM, IGR, KRS, NOY, SVV, VGN, YAK)
    and those with more than two (BRZ)
    """
    def __init__(self, *args, **kwargs):
        super(NIES_aircraft, self).__init__(*args, **kwargs)
        self.sasakawa_aircraft_dir = self.rcf.get('nies.aircraft.data.folder')
        self.sasakawa_tower_dir = self.rcf.get('nies.tower.data.folder')
        self.sampling_strategy = self.rcf.get('nies.insitu.sampling.strategy', 'int')
        self.sampling_error = self.rcf.get('nies.insitu.sampling.error', 'float')
        self.data_dict = {}
        self.two_inlet_towers = ['AZV', 'DEM', 'IGR', 'KRS', 'NOY', 'SVV', 'VGN', 'YAK']
        self.multi_inlet_towers = ['BRZ']
        self.read_towers = self.rcf.get('nies.tower.data', 'bool', default=True)
        self.read_aircrafts = self.rcf.get('nies.aircraft.data', 'bool', default=False)
        self.site_id_dict = {}
        # This class will be called from, among others, ObsPack_umbrella. Therefore, we want to put stuff in a format
        # compliant with the site info structure of ObsPack_umbrella. If it's not used by the caller, no harm done.
        self.obspack_cats = {'nies-tower': None}
        site_id_offset = self.rcf.get('nies.co2.site_id.offset', 'int')
        all_tower_codes = self.two_inlet_towers + self.multi_inlet_towers
        self.obspack_cats['nies-tower'] = {}
        self.site_ids = {}
        for code in all_tower_codes:
            self.obspack_cats['nies-tower'][code.lower()] = {'n_obs': 0}
            self.site_ids[code] = site_id_offset
            site_id_offset = site_id_offset - 1

    def mid_afternoon(self, time_obj, lon):
        # return true if sample is mid-afternoon, false otherwise
        # mid-afternoon is defined as between 11 and 4, LT
        local_time = time_obj + timedelta(hours=lon/15.0)
        if time(11,0,0) <= local_time.time() <= time(16,0,0):
            return True
        else:
            return False

    def readAMEStower_two(self, file_name):
        """
        All AMES files have the same first seven lines.

        NLHEAD FFI      Line 1
        ONAME       Line 2
        ORG         Line 3
        SNAME       Line 4
        MNAME       Line 5
        IVOL NVOL       Line 6
        DATE RDATE      Line 7

        NLHEAD      Number of lines in file header
        FFI         File format index
        ONAME       List of author(s) in the format Lastname, Firstname separated by an arbitrary character (for example, a hyphen or a semi-colon).
        ORG         Organisation name (university, institute, etc). May include address and phone numbers.
        SNAME       Source of data, i.e. instrument, platform, model name, etc.
        MNAME       Name of mission, campaign, programme and/or project.
        NVOL        Total number of files belonging to the considered dataset (i.e. with same ONAME, ORG, SNAME, MNAME).
        IVOL        Number of the file in the above dataset (between 1 and NVOL).
        DATE        Universal Time date at which the data within the file begin.
        RDATE       Universal Time date at which the data within the file have been reduced or revised.
        """
        with open(file_name, 'r') as fid:
            all_lines = fid.readlines()
        num_header = int(all_lines[0].split()[0])
        header = all_lines[:num_header]
        data = all_lines[num_header:]
        r = re.compile('^SITE_NAME')
        # get the site name
        site_line = [l for l in header if r.search(l)][0]
        site_name = site_line.split(':')[1].strip()
        # get latitude and longitude
        r = re.compile('^LATITUDE')
        lat_line = [l for l in header if r.search(l)][0]
        lat = float(lat_line.split(':')[1].split()[0])
        r = re.compile('^LONGITUDE')
        lon_line = [l for l in header if r.search(l)][0]
        lon = float(lon_line.split(':')[1].split()[0])
        r = re.compile('^ELEVATION')
        elev_line = [l for l in header if r.search(l)][0]
        elev = float(elev_line.split(':')[1].split()[0])
        # What is the time zone of the measurements?
        r = re.compile('^TIME_ZONE')
        tz_line = [l for l in header if r.search(l)][0]
        tz_offset = int(tz_line.split(':')[1])
        tz_offset = timedelta(hours=tz_offset)
        # Each data line has two measurements, at the high and low inlets
        r = re.compile('Low inlet')
        inlet_line = [l for l in header if r.search(l)][0]
        low_inlet = float(inlet_line.split()[3]) + elev
        r = re.compile('High inlet')
        inlet_line = [l for l in header if r.search(l)][0]
        high_inlet = float(inlet_line.split()[3]) + elev
        # get the site ID in self.station_name_list
        stat_code = os.path.basename(file_name).split('_')[0].upper()
        stat_type = 'CM'
        stat_id = self.site_ids[stat_code]

        #if (stat_code, stat_type) in self.station_name_list:
            #if stat_id not in self.site_id_dict:
                #stat_lat = self.station_coords[(stat_code, stat_type)]['lat']
                #stat_alt = self.station_coords[(stat_code, stat_type)]['alt']
                #stat_lon = self.station_coords[(stat_code, stat_type)]['lon']
                #self.site_id_dict[stat_id] = {'type': stat_type, 'code': stat_code, \
                    #'name': self.station_name_dict[(stat_code, stat_type)], 'static': True}
                #self.site_id_dict[stat_id].update(self.station_coords[(stat_code, stat_type)])
        #else:
            #sys.stderr.write('Warning :: station code %s not found\n'%stat_code)

        # Add to self.obspack_cats as well
        temp_dict = {'code': stat_code.upper(), 'type': stat_type, 'lab': 'National Institute for Environmental Research', \
            'lab_abbr': 'NIES', 'altitude': elev, 'latitude': lat, 'longitude': lon, 'elev': elev, 'name': site_name, \
            'site_id': stat_id, 'site_num': stat_id, 'static': True, 'country': 'Russia'}
        self.obspack_cats['nies-tower'][stat_code.lower()].update(temp_dict)

        if site_name not in self.data_dict:
            self.data_dict[site_name] = defaultdict(list)

        n_obs = 0

        for line in data:
            fields = line.split()
            year = int(fields[1])
            month = int(fields[2])
            day = int(fields[3])
            hour = int(fields[4])
            minute = int(fields[5])
            co2_hi = float(fields[6])
            co2_lo = float(fields[7])
            obs_time = datetime(year, month, day, hour, minute) - tz_offset
            if self.StartDate <= obs_time <= self.EndDate and self.mid_afternoon(obs_time, lon):
                if 0.0 < co2_hi < 800.0:
                    self.data_dict[site_name]['time'].append(obs_time)
                    self.data_dict[site_name]['mix'].append(co2_hi)
                    self.data_dict[site_name]['latitude'].append(lat)
                    self.data_dict[site_name]['longitude'].append(lon)
                    self.data_dict[site_name]['altitude'].append(high_inlet)
                    self.data_dict[site_name]['site_id'].append(stat_id)
                    n_obs += 1
                # do not look at CO2 from the lower level
                #if 0.0 < co2_lo < 800.0:
                    #self.data_dict[site_name]['time'].append(obs_time)
                    #self.data_dict[site_name]['mix'].append(co2_lo)
                    #self.data_dict[site_name]['latitude'].append(lat)
                    #self.data_dict[site_name]['longitude'].append(lon)
                    #self.data_dict[site_name]['altitude'].append(low_inlet)

        self.obspack_cats['nies-tower'][stat_code.lower()]['n_obs'] += n_obs

    def readAMEStower_multi(self, file_name):
        with open(file_name, 'r') as fid:
            all_lines = fid.readlines()
        num_header = int(all_lines[0].split()[0])
        header = all_lines[:num_header]
        data = all_lines[num_header:]
        r = re.compile('^SITE_NAME')
        # get the site name
        site_line = [l for l in header if r.search(l)][0]
        site_name = site_line.split(':')[1].strip()
        # get latitude and longitude
        r = re.compile('^LATITUDE')
        lat_line = [l for l in header if r.search(l)][0]
        lat = float(lat_line.split(':')[1].split()[0])
        r = re.compile('^LONGITUDE')
        lon_line = [l for l in header if r.search(l)][0]
        lon = float(lon_line.split(':')[1].split()[0])
        r = re.compile('^ELEVATION')
        elev_line = [l for l in header if r.search(l)][0]
        elev = float(elev_line.split(':')[1].split()[0])
        # What is the time zone of the measurements?
        r = re.compile('^TIME_ZONE')
        tz_line = [l for l in header if r.search(l)][0]
        tz_offset = int(tz_line.split(':')[1])
        tz_offset = timedelta(hours=tz_offset)
        # there are multiple inlets
        header_line = header[-1]
        header_cols = header_line.split()
        co2_cols = []
        co2_alts = []
        for i_col, col_name in enumerate(header_cols):
            if col_name.startswith('CO2'):
                co2_cols.append(i_col)
                alt = float(col_name.split('_')[1].split('m')[0])
                co2_alts.append(alt)

        # which column is the highest level?
        i_max_alt = np.argmax(co2_alts)
        i_col = co2_cols[i_max_alt]
        alt = co2_alts[i_max_alt]

        # get the site ID in self.station_name_list
        stat_code = os.path.basename(file_name).split('_')[0].upper()
        stat_type = 'CM'
        stat_id = self.site_ids[stat_code]

        #if (stat_code, stat_type) in self.station_name_list:
            #if stat_id not in self.site_id_dict:
                #stat_lat = self.station_coords[(stat_code, stat_type)]['lat']
                #stat_alt = self.station_coords[(stat_code, stat_type)]['alt']
                #stat_lon = self.station_coords[(stat_code, stat_type)]['lon']
                #self.site_id_dict[stat_id] = {'type': stat_type, 'code': stat_code, \
                    #'name': self.station_name_dict[(stat_code, stat_type)], 'static': True}
                #self.site_id_dict[stat_id].update(self.station_coords[(stat_code, stat_type)])
        #else:
            #sys.stderr.write('Warning :: station code %s not found\n'%stat_code)

        # Add to self.obspack_cats as well
        temp_dict = {'code': stat_code.upper(), 'type': stat_type, 'lab': 'National Institute for Environmental Research', \
            'lab_abbr': 'NIES', 'altitude': elev, 'latitude': lat, 'longitude': lon, 'elev': elev, 'name': site_name, \
            'site_id': stat_id, 'site_num': stat_id, 'static': True, 'country': 'Russia'}
        self.obspack_cats['nies-tower'][stat_code.lower()].update(temp_dict)

        if site_name not in self.data_dict:
            self.data_dict[site_name] = defaultdict(list)

        n_obs = 0

        for line in data:
            fields = line.split()
            year = int(fields[1])
            month = int(fields[2])
            day = int(fields[3])
            hour = int(fields[4])
            minute = int(fields[5])
            obs_time = datetime(year, month, day, hour, minute) - tz_offset
            if self.StartDate <= obs_time <= self.EndDate and self.mid_afternoon(obs_time, lon):
                #for i_col, alt in zip(co2_cols, co2_alts):
                co2 = float(fields[i_col])
                if 0.0 < co2 < 800.0:
                    self.data_dict[site_name]['time'].append(obs_time)
                    self.data_dict[site_name]['mix'].append(co2)
                    self.data_dict[site_name]['latitude'].append(lat)
                    self.data_dict[site_name]['longitude'].append(lon)
                    self.data_dict[site_name]['altitude'].append(alt+elev)
                    self.data_dict[site_name]['site_id'].append(stat_id)
                    n_obs += 1

        self.obspack_cats['nies-tower'][stat_code.lower()]['n_obs'] += n_obs

    def readAMESaircraft(self, file_name):
        with open(file_name, 'r') as fid:
            all_lines = fid.readlines()
        num_header = int(all_lines[0].split()[0])
        header = all_lines[:num_header]
        data = all_lines[num_header:]
        r = re.compile('^SITE_NAME')
        # get the site name
        site_line = [l for l in header if r.search(l)][0]
        site_name = site_line.split(':')[1].strip()
        # get latitude and longitude
        r = re.compile('^LATITUDE')
        lat_line = [l for l in header if r.search(l)][0]
        lat = float(lat_line.split(':')[1].split()[0])
        r = re.compile('^LONGITUDE')
        lon_line = [l for l in header if r.search(l)][0]
        lon = float(lon_line.split(':')[1].split()[0])
        # What is the time zone of the measurements?
        r = re.compile('^TIME_ZONE')
        tz_line = [l for l in header if r.search(l)][0]
        tz_offset = int(tz_line.split(':')[1])
        tz_offset = timedelta(hours=tz_offset)
        if site_name not in self.data_dict:
            self.data_dict[site_name] = defaultdict(list)
        for line in data:
            fields = line.split()
            year = int(fields[1])
            month = int(fields[2])
            day = int(fields[3])
            hour = int(fields[4])
            minute = int(fields[5])
            second = int(fields[6])
            co2 = float(fields[7])
            alt = 1000.0 * float(fields[8])
            obs_time = datetime(year, month, day, hour, minute, second) - tz_offset
            if self.StartDate <= obs_time <= self.EndDate:
                if 0.0 < co2 < 800.0:
                    self.data_dict[site_name]['time'].append(obs_time)
                    self.data_dict[site_name]['mix'].append(co2)
                    self.data_dict[site_name]['latitude'].append(lat)
                    self.data_dict[site_name]['longitude'].append(lon)
                    self.data_dict[site_name]['altitude'].append(alt)
                    self.data_dict[site_name]['site_id'].append(-999) # aircraft does not have site id

    def consolidateObs(self):
        if self.read_towers:
            # get a list of tower files
            my_years = set()
            cur_date = self.StartDate
            while cur_date < self.EndDate:
                my_years.add(cur_date.year)
                cur_date += timedelta(days=1)
            for tower_code in self.two_inlet_towers: # most of them
                all_files = glob(os.path.join(self.sasakawa_tower_dir, tower_code+'_co2_*.txt'))
                relevant_files = [f for f in all_files if int(os.path.basename(f).split('_')[2]) in my_years]
                for file_name in relevant_files:
                    self.readAMEStower_two(file_name)
            for tower_code in self.multi_inlet_towers: # only Berezorechka
                all_files = glob(os.path.join(self.sasakawa_tower_dir, tower_code+'_co2_*.txt'))
                relevant_files = [f for f in all_files if int(os.path.basename(f).split('_')[2]) in my_years]
                for file_name in relevant_files:
                    self.readAMEStower_multi(file_name)
        if self.read_aircrafts:
            # get a list of aircraft files
            all_files = glob(os.path.join(self.sasakawa_aircraft_dir, '*.txt'))
            relevant_files = []
            for file_name in all_files:
                cur_date = datetime.strptime(os.path.basename(file_name), "%Y%m%d.txt")
                if self.StartDate <= cur_date < self.EndDate:
                    relevant_files.append(file_name)
            for file_name in relevant_files:
                self.readAMESaircraft(file_name)
        # now consolidate everything in self.data_dict into a bigger structure
        self.station_names = []
        self.station_codes = []
        self.write_dict = defaultdict(list)
        for stat_num, (stat_code, stat_data) in enumerate(self.data_dict.items()):
            n_obs = len(stat_data['mix'])
            for key, value in stat_data.items():
                self.write_dict[key].extend(value)
            self.station_names.append(stat_code)
            self.station_codes.append(stat_code)
        # Now we can sort according to time
        sort_order = np.argsort(np.array(self.write_dict['time']))
        for key, value in self.write_dict.items():
            self.write_dict[key] = np.array(value)[sort_order]

        # Add the 'n_obs' key to self.site_id_dict
        for id in self.site_id_dict.keys():
            n_obs = np.count_nonzero(self.write_dict['site_id'] == id)
            self.site_id_dict[id]['n_obs'] = n_obs

class ObsPack(Observations):
    """
    A class for parsing ObsPack observations. The files to be parsed are given as follows.

    (1) A list of categories, such as 'aircraft-pfp', 'aircraft-flask', 'surface-insitu', etc.
    (2) For each category, a list of three letter codes for which observations should be parsed
    """
    def __init__(self, *args, **kwargs):
        super(ObsPack, self).__init__(*args, **kwargs)
        self.obspack_dir = self.rcf.get('CO2.obspack.data.dir')
        self.site_id_offset = self.rcf.get('CO2.obspack.site_id.offset', 'int')
        site_id = self.site_id_offset
        self.obspack_cats = OrderedDict.fromkeys(self.rcf.get('CO2.obspack.categories').split())
        for obspack_cat in self.obspack_cats.keys():
            self.obspack_cats[obspack_cat] = OrderedDict.fromkeys(self.rcf.get('CO2.obspack.%s.sites'%obspack_cat).split())
            for site in self.obspack_cats[obspack_cat].keys():
                self.obspack_cats[obspack_cat][site] = {'site_id': site_id}
                site_id = site_id - 1
        self.sampling_error = OrderedDict.fromkeys(self.rcf.get('CO2.obspack.categories').split())
        for obspack_cat in self.sampling_error.keys():
            self.sampling_error[obspack_cat] = self.rcf.get('CO2.obspack.%s.sampling.error'%obspack_cat, 'float')
        self.sampling_strategy = self.rcf.get('CO2.obspack.sampling.strategy', 'int')
        # Right now this class has only been coded to handle all observations in a single file
        # So raise an error if that is not satisfied
        if self.point_split_period != 'a':
            sys.stderr.write('ObsPack class cannot handle point observations split over multiple files\n')
            raise

    def parseFile(self, file_name):
        site_code = os.path.basename(file_name).split('_')[1]
        site_cat = os.path.basename(file_name).split('_')[2]
        ret_dict = {}
        with Dataset(file_name, 'r') as fid:
            idx = np.where(fid.variables['obs_flag'][:] == 1)[0]
            # now filter by time
            all_times = [datetime(*d) for d in fid.variables['time_components'][idx]]
            t_idx = np.array([i for i,d in enumerate(all_times) if self.StartDate <= d <= self.EndDate])
            if len(t_idx) > 0:
                idx = idx[t_idx]
                ret_dict['times'] = [datetime(*d) for d in fid.variables['time_components'][idx]]
                ret_dict['co2'] = 1.0E6 * fid.variables['value'][idx]
                n_obs = len(ret_dict['co2'])
                ret_dict['co2_err'] = self.sampling_error[site_cat] * np.ones_like(ret_dict['co2'])
                ret_dict['samp_strat'] = self.sampling_strategy * np.ones(n_obs, np.int16)
                ret_dict['lat'] = fid.variables['latitude'][idx]
                ret_dict['alt'] = fid.variables['altitude'][idx]
                ret_dict['lon'] = fid.variables['longitude'][idx]
                ret_dict['obspack_id'] = fid.variables['obspack_num'][idx]
                ret_dict['site_id'] = self.obspack_cats[site_cat][site_code]['site_id'] * np.ones(n_obs, np.int32)
                loc_dict = {'site_name': fid.site_name, 'lab': fid.lab_name, 'lab_abbr': fid.lab_abbr}
                self.obspack_cats[site_cat][site_code].update(loc_dict)
        print '%s parsed'%os.path.basename(file_name)
        return ret_dict

    def reverseLookup(self):
        ret_dict = OrderedDict()
        for site_cat, site_cat_data in self.obspack_cats.items():
            for site_code, site_data in site_cat_data.items():
                site_id = site_data['site_id']
                ret_dict[site_id] = {'site_code': site_code, 'site_category': site_cat}
                for k,v in site_data.items():
                    if k != 'site_id':
                        ret_dict[site_id][k] = v
        return ret_dict

    def consolidateObs(self):
        all_files = glob(os.path.join(self.obspack_dir, 'co2*.nc'))
        self.write_dict = defaultdict(list)
        for file_name in all_files:
            site_code = os.path.basename(file_name).split('_')[1]
            site_cat = os.path.basename(file_name).split('_')[2]
            if site_cat in self.obspack_cats and site_code in self.obspack_cats[site_cat]:
                ret_dict = self.parseFile(file_name)
                for key, value in ret_dict.items():
                    self.write_dict[key].extend(value)
        # now write out the dictionary of site codes, etc.
        site_dict = self.reverseLookup()
        out_file = os.path.join(self.rcf.get('output.point.input.dir'), 'obspack_site_codes.dict')
        checkDir(out_file)
        with open(out_file, 'wb') as fid:
            pickle.dump(site_dict, fid, pickle.HIGHEST_PROTOCOL)

class ObsPack_umbrella(Observations):
    """
    A class for parsing ObsPack obs, for the case when we're getting obs solely from ObsPack, and not from CT2013.
    """
    def __init__(self, *args, **kwargs):
        super(ObsPack_umbrella, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.obspack_dir = self.rcf.get('CO2.obspack.data.dir')
        # Since ObsPack is the primary data source for us, no need to have a negative site ID, just start from 1 and go up
        site_id = 1
        # self.obspack_cats is a list like ['aircraft-pfp', 'aircraft-flask', 'shipboard-flask', ...]
        self.obspack_cats = OrderedDict.fromkeys(self.rcf.get('CO2.obspack.categories').split())
        for obspack_cat in self.obspack_cats.keys():
            # Read a list of (three-letter) site codes such as 'bgi bne car cma dnd esp etl ftl ...'
            self.obspack_cats[obspack_cat] = OrderedDict.fromkeys(self.rcf.get('CO2.obspack.%s.sites'%obspack_cat).split())
            for site in self.obspack_cats[obspack_cat].keys():
                self.obspack_cats[obspack_cat][site] = {'site_id': site_id, 'n_obs': 0}
                site_id = site_id + 1
        self.sampling_error = OrderedDict.fromkeys(self.rcf.get('CO2.obspack.categories').split())
        for obspack_cat in self.sampling_error.keys():
            self.sampling_error[obspack_cat] = self.rcf.get('CO2.obspack.%s.sampling.error'%obspack_cat, 'float')
        self.sampling_strategy = self.rcf.get('CO2.obspack.sampling.strategy', 'int')
        # Apart from the rc keys read in by ObsPack.__init__, we also need the following
        self.accept_types = {}
        for obspack_cat in self.obspack_cats:
            key = 'CO2.obspack.%s.accept'%obspack_cat
            # For each category, such as 'aircraft-pfp', there will be a list of event types we will accept,
            # such as ['allvalid', 'representative', ...]
            self.accept_types[obspack_cat] = self.rcf.get(key).split()
        self.mountaintop_sites = self.rcf.get('CO2.obspack.mountaintop.sites').split()
        self.allday_sites = self.rcf.get('CO2.obspack.allday.sites').split()
        # Some sites such as POC cruises and aircraft profiles need not be filtered by time of day, so get the site types that do
        self.tod_filter_site_types = self.rcf.get('CO2.obspack.time_of_day.subsample').split()
        # Some measurements need to be rejected. For each type of measurement (surface-flask, aircraft-pfp, etc.), they should
        # be provided as a list of site_code:lab_code in the rc file.
        self.reject_sites = defaultdict(list)
        for site_type in self.obspack_cats.keys():
            key = 'CO2.obspack.%s.reject_sites'%site_type
            reject_sites = self.rcf.get(key, default='').split()
            for site in reject_sites:
                atoms = site.split(':')
                site_code = atoms[0]
                lab_code = int(atoms[1])
                self.reject_sites[site_type].append((site_code, lab_code))

        # For sites classified as static platforms, we will output a continuous modeled time series
        self.static_platforms = self.rcf.get('CO2.obspack.static_platforms').split()
        # For the different static platforms, what are the types in the station list file?
        self.static_platform_types = {}
        for platform in self.static_platforms:
            key = 'CO2.obspack.%s.obs_type'%platform
            self.static_platform_types[platform] = self.rcf.get(key)

        # For ObsPack, keep track of which PI contributed data for which station
        self.stations_by_PI = defaultdict(list)

        # Also, define a list of keys to read in for each station
        self.site_info_mapping = {}
        map_key = self.rcf.get('CO2.obspack.site.metadata').split()
        for key_tuple in map_key:
            short,long = key_tuple.split(',')
            self.site_info_mapping[short] = long # e.g., self.site_info_mapping['lon'] = 'site_longitude'

        # Should we add NIES data?
        self.add_nies_data = self.rcf.get('CO2.add.nies.data', 'bool')
        # Should we add Gatti's data over the Amazon?
        self.add_gatti = self.rcf.get('CO2.add.gatti.data', 'bool', default=False)
        # We add flags to (optionally) disable the assimilation of these data
        self.assimilate = self.rcf.get('CO2.point.assimilate', 'bool', default=True)
        self.unassim_mdm = self.rcf.get('CO2.point.unassimilate.mdm', 'float')

        site_id_save_file = self.rcf.get('output.point.site_id.dict')
        bname = os.path.basename(site_id_save_file)
        bname = '%s_%s'%(self.tracer, bname)
        self.site_id_save_file = os.path.join(os.path.dirname(site_id_save_file), bname)

        # We need to keep track of the samples lost to the time of day filter
        self.lost_samples = {}

    def getStationList(self, file_name):
        # The getStationList from Observations needs to be overwritten with a dummy routine because we want to write it out
        self.station_name_list = []

    def writeStationFile(self):
        # Create the lines first
        write_lines = [' ID     LAT     LON     ALT TP STATION NAME\n']
        for site_cat, site_cat_data in self.obspack_cats.items():
            for site_code, site_data in site_cat_data.items():
                if site_data['n_obs'] > 0 and site_data['static']:
                    line = '%3s %7.2f %7.2f %7.1f %2s %s (%s)\n'%\
                        (site_code.upper(), site_data['latitude'], site_data['longitude'], site_data['altitude'], \
                        site_data['type'], site_data['name'], site_data['country'])
                    write_lines.append(line)

        write_lines.sort()
        file_name = self.rcf.get('output.station.timeseries.filename')
        with open(file_name, 'w') as fid:
            fid.writelines(write_lines)

    def consolidateObs(self):
        all_files = glob(os.path.join(self.obspack_dir, 'co2*.nc'))
        self.write_dict = defaultdict(list)

        for file_name in all_files:
            file_wo_ext = os.path.splitext(os.path.basename(file_name))[0]
            site_code = os.path.basename(file_wo_ext).split('_')[1] # oxk, mlo, etc.
            site_cat  = os.path.basename(file_wo_ext).split('_')[2] # tower-insitu, surface-flask, etc.
            lab_code = int(os.path.basename(file_wo_ext).split('_')[3]) # 1 for NOAA, 6 for EC, etc.
            meas_type = os.path.basename(file_wo_ext).split('_')[4] # representative, afternoon-379magl, etc.
            # self.reject_sites[site_type] is a list of tuples (site_code, lab_code)
            if (site_code, lab_code) in self.reject_sites[site_cat]:
                continue

            if site_cat in self.obspack_cats and site_code in self.obspack_cats[site_cat]:
                if meas_type in self.accept_types[site_cat]:
                    ret_dict = self.parseFile(file_name)
                    for key, value in ret_dict.items():
                        self.write_dict[key].extend(value)

        # Debug: print the sites from which samples have been lost
        lost_list = [(k,) + v for k,v in self.lost_samples.items()]
        lost_list.sort(key=lambda x: x[3])
        lost_list = [l for l in lost_list if l[1] != l[2]]
        for l in lost_list:
            print '%s :: kept %i of %i samples (%i %% loss)'%(l[0], l[2], l[1], l[3])
        # End debug

        # At this point, self.write_dict has the keys 'times', 'co2', 'co2_err', 'samp_strat', 'lat', 'alt', 'lon', 'obspack_id', 'site_id'

        # should we add NIES data? now is the time to do it.
        if self.add_nies_data:
            nies = NIES_aircraft(self.StartDate.timetuple()[:6], self.EndDate.timetuple()[:6])
            nies.consolidateObs()
            # now nies.write_dict has the keys 'time', 'mix', 'latitude', 'longitude', 'altitude', 'site_id'
            for k,v in nies.write_dict.items():
                nies.write_dict[k] = list(v)

            nies_obs = len(nies.write_dict['mix'])
            self.write_dict['times'].extend(nies.write_dict['time'])
            self.write_dict['lat'].extend(nies.write_dict['latitude'])
            self.write_dict['lon'].extend(nies.write_dict['longitude'])
            self.write_dict['alt'].extend(nies.write_dict['altitude'])
            self.write_dict['co2'].extend(nies.write_dict['mix'])
            self.write_dict['site_id'].extend(nies.write_dict['site_id'])
            # The following are not in nies.write_dict, so make up
            self.write_dict['co2_err'].extend(nies.sampling_error * np.ones(nies_obs, np.float64))
            self.write_dict['obspack_id'].extend(-999 * np.ones(nies_obs, np.int32))
            self.write_dict['samp_strat'].extend(nies.sampling_strategy * np.ones(nies_obs, np.int16))
            # save some memory
            self.obspack_cats.update(nies.obspack_cats)
            del nies

        # convert self.obspack_cats to a dictionary with site IDs as keys
        write_dict = self.reverseLookup()

        # Should we add Luciana Gatti's data over the Amazon?
        if self.add_gatti:
            gatti = Gatti_CO2(self.StartDate.timetuple()[:6], self.EndDate.timetuple()[:6])
            data_dict = gatti.createObservations()
            for k,v in data_dict.items():
                data_dict[k] = list(v)

            gatti_obs = len(data_dict['co2'])
            # Add stuff from data_dict to the main data structure
            self.write_dict['times'].extend(data_dict['date'])
            self.write_dict['lat'].extend(data_dict['lat'])
            self.write_dict['lon'].extend(data_dict['lon'])
            self.write_dict['alt'].extend(data_dict['alt'])
            self.write_dict['co2_err'].extend(data_dict['co2_error'])
            self.write_dict['co2'].extend(data_dict['co2'])
            self.write_dict['site_id'].extend(data_dict['site_id'])
            # The following are not in data_dict, so create
            self.write_dict['obspack_id'].extend(-999 * np.ones(gatti_obs, np.int32))
            self.write_dict['samp_strat'].extend(gatti.samp_strat * np.ones(gatti_obs, np.int16))
            # save some memory
            del gatti

        # now write out the dictionary of site codes, etc.

        checkDir(self.site_id_save_file)
        with open(self.site_id_save_file, 'wb') as fid:
            pickle.dump(write_dict, fid, pickle.HIGHEST_PROTOCOL)

        self.writeStationFile()

        # sort by time and put into dictionary
        sort_order = np.argsort(self.write_dict['times'])
        for k, v in self.write_dict.items():
            self.write_dict[k] = np.array(v)[sort_order]

    def parseFile(self, file_name):
        site_code = os.path.basename(file_name).split('_')[1] # such as 'drp'
        site_cat = os.path.basename(file_name).split('_')[2]  # such as 'shipboard-flask'

        ret_dict = {}
        with Dataset(file_name, 'r') as fid:
            # first filter by the obvious QC flag
            idx = np.where(fid.variables['obs_flag'][:] == 1)[0]
            # now filter by time
            all_times = np.array([datetime(*d) for d in fid.variables['time_components'][idx]])
            t_idx = np.array([i for i,d in enumerate(all_times) if self.StartDate <= d <= self.EndDate])
            if len(t_idx) == 0:
                return ret_dict

            # filter by mid-afternoon or nighttime samples
            utc_times = all_times[t_idx]
            longitudes = fid.variables['longitude'][idx][t_idx]

            if (site_cat not in self.tod_filter_site_types) or (site_code in self.allday_sites):
                # Accept measurements at any time of day
                lt_idx = np.arange(len(t_idx))
            elif site_code in self.mountaintop_sites:
                # For mountaintop sites, accept measurements from 0:00 to 07:00 local solar time
                lt_idx = self.filter_by_local_time(utc_times, longitudes, 'night')
            else:
                # For all other sites, accept measurements from 11:00 to 17:00 local solar time
                lt_idx = self.filter_by_local_time(utc_times, longitudes, 'afternoon')

            if len(lt_idx) > 0:
                n_old = len(t_idx)
                n_new = len(lt_idx)
                self.lost_samples[os.path.basename(file_name)] = (n_old, n_new, 100.0*(n_old-n_new)/n_old)
                t_idx = t_idx[lt_idx]
            else:
                print 'Check afternoon/night selection for %s (0 samples during requisite LST interval)'%os.path.basename(file_name)
                return ret_dict

            site_id = self.obspack_cats[site_cat][site_code]['site_id']
            idx = idx[t_idx]
            ret_dict['times'] = [datetime(*d) for d in fid.variables['time_components'][idx]]
            ret_dict['co2'] = 1.0E6 * fid.variables['value'][idx]
            n_obs = len(ret_dict['co2'])
            ret_dict['co2_err'] = self.sampling_error[site_cat] * np.ones_like(ret_dict['co2'])
            ret_dict['samp_strat'] = self.sampling_strategy * np.ones(n_obs, np.int16)
            ret_dict['lat'] = fid.variables['latitude'][idx]
            ret_dict['alt'] = fid.variables['altitude'][idx]
            ret_dict['lon'] = fid.variables['longitude'][idx]
            ret_dict['obspack_id'] = fid.variables['obspack_num'][idx]
            ret_dict['site_id'] = site_id * np.ones(n_obs, np.int32)

            loc_dict = {}
            for short,long in self.site_info_mapping.items():
                try:
                    loc_dict[short] = getattr(fid, long)
                except AttributeError:
                    raise AttributeError('File %s does not have attribute %s'%(file_name, long))
            loc_dict['n_obs'] = n_obs

            if site_cat in self.static_platforms:
                loc_dict['altitude'] = np.average(ret_dict['alt'])
                loc_dict['static'] = True
                stat_type = self.static_platform_types[site_cat]
                # Add this site to the list of sites, if it's not there already
                if not (site_code.upper(), stat_type) in self.station_name_list:
                    loc_dict['code'] = site_code.upper()
                    loc_dict['type'] = stat_type
                    loc_dict['site_num'] = site_id
                    self.station_name_list.append((site_code.upper(), stat_type))
            else:
                loc_dict['static'] = False

            self.obspack_cats[site_cat][site_code].update(loc_dict)

            # Store the PI's name and email
            pi_name = getattr(fid, 'provider_1_name')
            pi_email = getattr(fid, 'provider_1_email')
            pi_site_name = getattr(fid, 'site_name')
            pi_site_code = getattr(fid, 'site_code')
            pi_key = '%s (%s)'%(pi_name, pi_email)
            pi_value = '%s (%s)'%(pi_site_name, pi_site_code)
            self.stations_by_PI[pi_key].append(pi_value)

        return ret_dict

    def __call__(self):
        self.consolidateObs()
        n_obs = len(self.write_dict['co2'])

        self.PointObservation[self.tracer]['dimensions'] = {'id': n_obs}
        self.PointObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.PointObservation[self.tracer]['attr_dict'] = OrderedDict()
        self.PointObservation[self.tracer]['variable_attrs'] = OrderedDict()

        self.PointObservation[self.tracer]['id'] = np.arange(1, n_obs+1, dtype=np.int32)
        self.PointObservation[self.tracer]['variable_shapes']['id'] = ('id',)

        self.PointObservation[self.tracer]['lat'] = np.float64(self.write_dict['lat'])
        self.PointObservation[self.tracer]['variable_shapes']['lat'] = ('id',)

        self.PointObservation[self.tracer]['lon'] = np.float64(self.write_dict['lon'])
        self.PointObservation[self.tracer]['variable_shapes']['lon'] = ('id',)

        self.PointObservation[self.tracer]['alt'] = np.float64(self.write_dict['alt'])
        self.PointObservation[self.tracer]['variable_shapes']['alt'] = ('id',)

        self.PointObservation[self.tracer]['mixing_ratio'] = np.float64(self.write_dict['co2'])
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio'] = ('id',)

        if self.assimilate:
            self.PointObservation[self.tracer]['mixing_ratio_error'] = np.float64(self.write_dict['co2_err'])
        else:
            self.PointObservation[self.tracer]['mixing_ratio_error'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.PointObservation[self.tracer]['variable_attrs']['mixing_ratio_error'] = [('comment', 'These data will not be assimilated')]
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio_error'] = ('id',)

        self.PointObservation[self.tracer]['date_components'] = np.array([d.timetuple()[:6] for d in self.write_dict['times']], np.int16)
        self.PointObservation[self.tracer]['variable_shapes']['date_components'] = ('id', 'idate')

        self.PointObservation[self.tracer]['sampling_strategy'] = np.int16(self.write_dict['samp_strat'])
        self.PointObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('id',)

        self.PointObservation[self.tracer]['station_id'] = np.int32(self.write_dict['site_id'])
        self.PointObservation[self.tracer]['variable_shapes']['station_id'] = ('id',)

        self.PointObservation[self.tracer]['obspack_num'] = np.int32(self.write_dict['obspack_id'])
        self.PointObservation[self.tracer]['variable_shapes']['obspack_num'] = ('id',)

        if self.add_gatti:
            # Add the site IDs for Gatti's sites to the comment
            gatti_co_strings = []
            site_list = self.rcf.get('gatti.co2.profile.sites').split()
            site_offset = self.rcf.get('gatti.co2.site_id.offset', 'int')
            for i_site, site_code in enumerate(site_list):
                gatti_co_strings.append('%i => %s'%(site_offset-i_site, site_code.upper()))
            gatti_co_string = ', '.join(gatti_co_strings)
            self.PointObservation[self.tracer]['variable_attrs']['station_id'] = [('Gatti_site_indices', gatti_co_string)]

        # Add the list of PIs as attributes
        for i_pi, (pi_name, pi_stations) in enumerate(self.stations_by_PI.items()):
            attr_name = 'ObsPack_PI_%03i'%(i_pi+1)
            attr_value = pi_name
            self.PointObservation[self.tracer]['attr_dict'][attr_name] = attr_value
            attr_name = 'ObsPack_PI_%03i_stations'%(i_pi+1)
            attr_value = ', '.join(pi_stations)
            self.PointObservation[self.tracer]['attr_dict'][attr_name] = attr_value

class CT2013_CO2(Observations):
    """
    This is a class that reads preprocessed obs from CT2013.
    """
    def __init__(self, *args, **kwargs):
        super(CT2013_CO2, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.input_dir = self.rcf.get('point.ct2013.flask.data.folder')
        self.input_file = self.rcf.get('point.ct2013.flask.input.file')
        self.sampling_strat = self.rcf.get('point.ct2013.sampling.strategy', 'int')
        # Tolerances for matching observations to stations in the list
        self.lat_tol = 0.5
        self.lon_tol = 0.5
        self.alt_tol = 20.0
        # Should we add NIES data?
        self.add_nies_data = self.rcf.get('CO2.add.nies.data', 'bool')
        # Should we add ObsPack data that are not part of CT2013?
        self.add_obspack_data = self.rcf.get('CO2.add.obspack.data', 'bool')
        # Should we add Gatti's data over the Amazon?
        self.add_gatti = self.rcf.get('CO2.add.gatti.data', 'bool', default=False)
        # We add flags to (optionally) disable the assimilation of these data
        self.assimilate = self.rcf.get('CO2.point.assimilate', 'bool', default=True)
        if not self.assimilate:
            self.unassim_mdm = self.rcf.get('CO2.point.unassimilate.mdm', 'float')
        self.site_id_dict = {} # Translation from site_id to station type/code/name
        site_id_save_file = self.rcf.get('output.point.site_id.dict')
        bname = os.path.basename(site_id_save_file)
        bname = '%s_%s'%(self.tracer, bname)
        self.site_id_save_file = os.path.join(os.path.dirname(site_id_save_file), bname)

    def idxPeriod(self, time_start, time_end):
        # select all measurements from self.input_file between the two times, only for assimilated data
        with Dataset(self.input_file, 'r') as fid:
            mdm = fid.variables['mdm'][:]
            valid_idx = np.where(-mdm.mask)[0]
            all_times = [datetime(*d) for d in fid.variables['time_components'][:][valid_idx]]
        t_idx = np.array([i for i,d in enumerate(all_times) if time_start <= d < time_end])
        valid_idx = valid_idx[t_idx]
        return valid_idx

    def createObservations(self):
        valid_idx = self.idxPeriod(self.StartDate, self.EndDate)
        with Dataset(self.input_file, 'r') as fid:
            all_lats        = list(fid.variables['latitude'][:][valid_idx])
            all_lons        = list(fid.variables['longitude'][:][valid_idx])
            all_alts        = list(fid.variables['altitude'][:][valid_idx])
            all_mdm         = list(fid.variables['mdm'][:][valid_idx]/2.0) # divide by two because CT always errs on the side of caution
            all_time        = [datetime(*d) for d in fid.variables['time_components'][:][valid_idx]]
            all_co2         = list(1.0E6 * fid.variables['value'][:][valid_idx])
            all_obspack_id  = fid.variables['obspack_id'][:][valid_idx]
        all_obspack_id = [''.join(id.compressed()) for id in all_obspack_id]
        # now get the numerical obspack IDs
        all_obspack_num = [int(o.split('~')[-1]) for o in all_obspack_id]
        # create an array of sampling strategies
        all_sampling_strat = list(self.sampling_strat * np.ones(len(all_co2), np.int16))
        # figure out the station IDs
        all_stat_id = []
        for i, (alt, lat, lon, obs_id) in enumerate(zip(all_alts, all_lats, all_lons, all_obspack_id)):
            rel_part = obs_id.split('co2_')[-1]
            comps = rel_part.split('_')
            stat_code = comps[0].upper()
            seek_stat_id = True
            if comps[1] in ['aircraft-flask', 'aircraft-pfp', 'shipboard-flask', 'aircraft-insitu']:
                seek_stat_id = False
                stat_id = -999
            elif comps[1] in ['surface', 'surface-flask', 'surface-pfp']:
                stat_type = 'FM'
            elif comps[1] in ['surface-insitu']:
                stat_type = 'CM'
            elif comps[1] in ['tower-insitu']:
                stat_type = 'TM'
            else:
                sys.stderr.write('Unknown obs type %s\n'%comps[1])
                seek_stat_id = False
                stat_id = -999

            # Sometimes, it can happen at a site with both flask and continuous measurements, that a few flask
            # samples will be filled from the in situ intake, for the purpose of calibration. In that case, even
            # though the stat_type is 'FM', the measurement should be compared with the time series corresponding
            # to the 'CM' coordinate. An example is MLO, where FM is at 3397 masl, CM is at 3437 masl, but there
            # are some flask samples taken at 3437 masl.

            if seek_stat_id:
                if (stat_code, stat_type) not in self.station_name_list:
                    stat_id = -999
                    sys.stderr.write("WARNING: Station %s:%s not in %s\n"%(stat_type, stat_code, self.stationlist_file))
                else:
                    stat_lat = self.station_coords[(stat_code, stat_type)]['lat']
                    stat_alt = self.station_coords[(stat_code, stat_type)]['alt']
                    stat_lon = self.station_coords[(stat_code, stat_type)]['lon']

                    # Some stations had incorrect coordinates for their sampling locations. For example,
                    # the CSIRO sampling at Alert was at 210 masl, (200 m elevation + 10 magl), but in
                    # ObsPack they are recorded as 40 masl. Those observations need to be corrected.

                    if (stat_code, stat_type) == ('ALT', 'FM') and 39. < alt < 41.:
                        alt = 210.0
                        all_alts[i] = alt

                    if (stat_code, stat_type) == ('ASK', 'FM'):
                        alt = 2710.0
                        lat = 23.26
                        lon = 5.63
                        all_alts[i] = alt
                        all_lats[i] = lat
                        all_lons[i] = lon

                    # End manual correction block

                    if abs(stat_lat-lat) < self.lat_tol and abs(stat_lon-lon) < self.lon_tol and abs(stat_alt-alt) < self.alt_tol:
                        stat_id = self.station_name_list.index((stat_code, stat_type)) + 1
                        # Store the correpondence between the station ID and the station type/code in a dictionary
                        if stat_id not in self.site_id_dict:
                            self.site_id_dict[stat_id] = {'type': stat_type, 'code': stat_code, 'name': self.station_name_dict[(stat_code, stat_type)]}
                            self.site_id_dict[stat_id].update(self.station_coords[(stat_code, stat_type)])
                    else:
                        # Is this a case discussed above, with a station that shifted? To make sure that both the old and
                        # the new locations are sampled, I've made a 'FM' and a 'CM' instance of the same station in some
                        # cases. So the old location is a FM while the new location is a CM. It's a hack, I know.
                        special_case = False
                        if stat_type == 'FM' and (stat_code, 'CM') in self.station_coords:
                            new_stat_lat = self.station_coords[(stat_code, 'CM')]['lat']
                            new_stat_alt = self.station_coords[(stat_code, 'CM')]['alt']
                            new_stat_lon = self.station_coords[(stat_code, 'CM')]['lon']
                            if abs(new_stat_lat-lat) < self.lat_tol and abs(new_stat_lon-lon) < self.lon_tol \
                                and abs(new_stat_alt-alt) < self.alt_tol:
                                    stat_id = self.station_name_list.index((stat_code, 'CM')) + 1
                                    special_case = True
                                    # Store the correpondence between the station ID and the station type/code in a dictionary
                                    if stat_id not in self.site_id_dict:
                                        self.site_id_dict[stat_id] = {'type': 'CM', 'code': stat_code, 'name': self.station_name_dict[(stat_code, 'CM')], \
                                            'lat': new_stat_lat, 'lon': new_stat_lon, 'alt': new_stat_alt}
                        if not special_case:
                            stat_id = -999
                            sys.stderr.write("WARNING: Coordinate mismatch at %s:%s, %s\n"\
                                %(stat_type,stat_code,obs_id))
                            if abs(stat_lat-lat) > self.lat_tol:
                                sys.stderr.write("latitude %.2f in list, %.2f in obs file\n"%(stat_lat, lat))
                            if abs(stat_lon-lon) > self.lon_tol:
                                sys.stderr.write("longitude %.2f in list, %.2f in obs file\n"%(stat_lon, lon))
                            if abs(stat_alt-alt) > self.alt_tol:
                                sys.stderr.write("altitude %.2f in list, %.2f in obs file\n"%(stat_alt, alt))
                            sys.stderr.write("\n")

            all_stat_id.append(stat_id)

        # Write out the station ID <--> Site mapping to a file
        out_file = self.site_id_save_file
        checkDir(out_file)
        with open(out_file, 'wb') as fid:
            pickle.dump(self.site_id_dict, fid, pickle.HIGHEST_PROTOCOL)

        # should we add NIES data? now is the time to do it.
        if self.add_nies_data:
            nies = NIES_aircraft(self.StartDate.timetuple()[:6], self.EndDate.timetuple()[:6])
            nies.consolidateObs()
            # now nies.write_dict has all the data, add it to the main data structures
            for k,v in nies.write_dict.items():
                nies.write_dict[k] = list(v)
            nies_obs = len(nies.write_dict['mix'])
            all_time.extend(nies.write_dict['time'])
            all_lats.extend(nies.write_dict['latitude'])
            all_lons.extend(nies.write_dict['longitude'])
            all_alts.extend(nies.write_dict['altitude'])
            all_mdm.extend(nies.sampling_error * np.ones(nies_obs, np.float64))
            all_co2.extend(nies.write_dict['mix'])
            all_stat_id.extend(nies.write_dict['site_id'])
            all_obspack_num.extend(-999 * np.ones(nies_obs, np.int32))
            all_sampling_strat.extend(nies.sampling_strategy * np.ones(nies_obs, np.int16))
            # save some memory
            del nies

        # should we add some NOAA and partner lab aircraft data?
        if self.add_obspack_data:
            obspack = ObsPack(self.StartDate.timetuple()[:6], self.EndDate.timetuple()[:6])
            obspack.consolidateObs()
            # add data from obspack.write_dict to the main data structure
            # check for duplicate entries first
            unique_idx = nonzero(-np.in1d(np.array(obspack.write_dict['obspack_id']), np.array(all_obspack_num)))[0]
            for k,v in obspack.write_dict.items():
                obspack.write_dict[k] = list(np.array(v)[unique_idx])
            all_time.extend(obspack.write_dict['times'])
            all_lats.extend(obspack.write_dict['lat'])
            all_lons.extend(obspack.write_dict['lon'])
            all_alts.extend(obspack.write_dict['alt'])
            all_mdm.extend(obspack.write_dict['co2_err'])
            all_co2.extend(obspack.write_dict['co2'])
            all_stat_id.extend(obspack.write_dict['site_id'])
            all_obspack_num.extend(obspack.write_dict['obspack_id'])
            all_sampling_strat.extend(obspack.write_dict['samp_strat'])
            # save some memory
            del obspack

        # Should we add Luciana Gatti's data over the Amazon?
        if self.add_gatti:
            gatti = Gatti_CO2(self.StartDate.timetuple()[:6], self.EndDate.timetuple()[:6])
            data_dict = gatti.createObservations()
            gatti_obs = len(data_dict['co2'])
            # Add stuff from data_dict to the main data structure
            all_time.extend(data_dict['date'])
            all_lats.extend(data_dict['lat'])
            all_lons.extend(data_dict['lon'])
            all_alts.extend(data_dict['alt'])
            all_mdm.extend(data_dict['co2_error'])
            all_co2.extend(data_dict['co2'])
            all_stat_id.extend(data_dict['site_id'])
            all_obspack_num.extend(-999 * np.ones(gatti_obs, np.int32))
            all_sampling_strat.extend(gatti.samp_strat * np.ones(gatti_obs, np.int16))
            # save some memory
            del gatti

        # sort by time and put into dictionary
        ret_dict = {}
        sort_order = np.argsort(all_time)
        ret_dict['times'] = np.array([d.timetuple()[:6] for d in all_time], np.int16)[sort_order]
        ret_dict['lats'] = np.array(all_lats, np.float64)[sort_order]
        ret_dict['lons'] = np.array(all_lons, np.float64)[sort_order]
        ret_dict['alts'] = np.array(all_alts, np.float64)[sort_order]
        ret_dict['co2_err'] = np.array(all_mdm, np.float64)[sort_order]
        ret_dict['co2'] = np.array(all_co2, np.float64)[sort_order]
        ret_dict['stat_id'] = np.array(all_stat_id, np.int32)[sort_order]
        ret_dict['id'] = np.arange(1, len(all_co2)+1, dtype=np.int32)
        ret_dict['obspack_id'] = np.array(all_obspack_num, np.int32)[sort_order]
        ret_dict['sampling_strategy'] = np.array(all_sampling_strat, np.int16)[sort_order]

        return ret_dict

    def __call__(self):
        data_dict = self.createObservations()
        n_obs = len(data_dict['co2'])

        self.PointObservation[self.tracer]['dimensions'] = {'id': n_obs}
        self.PointObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.PointObservation[self.tracer]['attr_dict'] = OrderedDict()
        self.PointObservation[self.tracer]['variable_attrs'] = OrderedDict()

        self.PointObservation[self.tracer]['id'] = np.arange(1, n_obs+1, dtype=np.int32)
        self.PointObservation[self.tracer]['variable_shapes']['id'] = ('id',)

        self.PointObservation[self.tracer]['lat'] = data_dict['lats']
        self.PointObservation[self.tracer]['variable_shapes']['lat'] = ('id',)

        self.PointObservation[self.tracer]['lon'] = data_dict['lons']
        self.PointObservation[self.tracer]['variable_shapes']['lon'] = ('id',)

        self.PointObservation[self.tracer]['alt'] = data_dict['alts']
        self.PointObservation[self.tracer]['variable_shapes']['alt'] = ('id',)

        self.PointObservation[self.tracer]['mixing_ratio'] = data_dict['co2']
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio'] = ('id',)

        if self.assimilate:
            self.PointObservation[self.tracer]['mixing_ratio_error'] = data_dict['co2_err']
        else:
            self.PointObservation[self.tracer]['mixing_ratio_error'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.PointObservation[self.tracer]['variable_attrs']['mixing_ratio_error'] = [('comment', 'These data will not be assimilated')]
        self.PointObservation[self.tracer]['variable_shapes']['mixing_ratio_error'] = ('id',)

        self.PointObservation[self.tracer]['date_components'] = data_dict['times']
        self.PointObservation[self.tracer]['variable_shapes']['date_components'] = ('id', 'idate')

        self.PointObservation[self.tracer]['sampling_strategy'] = data_dict['sampling_strategy']
        self.PointObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('id',)

        self.PointObservation[self.tracer]['station_id'] = data_dict['stat_id']
        self.PointObservation[self.tracer]['variable_shapes']['station_id'] = ('id',)

        self.PointObservation[self.tracer]['obspack_num'] = data_dict['obspack_id']
        self.PointObservation[self.tracer]['variable_shapes']['obspack_num'] = ('id',)

        if self.add_gatti:
            # Add the site IDs for Gatti's sites to the comment
            gatti_co_strings = []
            site_list = self.rcf.get('gatti.co2.profile.sites').split()
            site_offset = self.rcf.get('gatti.co2.site_id.offset', 'int')
            for i_site, site_code in enumerate(site_list):
                gatti_co_strings.append('%i => %s'%(site_offset-i_site, site_code.upper()))
            gatti_co_string = ', '.join(gatti_co_strings)
            self.PointObservation[self.tracer]['variable_attrs']['station_id'] = [('Gatti_site_indices', gatti_co_string)]

class createACOSInputFile_b34(Observations):

    def __init__(self, *args, **kwargs):
        super(createACOSInputFile_b34, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.input_folder = self.rcf.get('satellite.GOSAT.ACOS.CO2.data.folder')
        # datetime_object.strtime(self.filename_pattern) should give a pattern that can be fed to glob.glob to get all files for a certain day
        self.filename_pattern = self.rcf.get('satellite.GOSAT.ACOS.CO2.filename.pattern')
        self.output_folder = self.rcf.get('output.satellite.output.directory')
        self.exclude_medgain = self.rcf.get('satellite.GOSAT.ACOS.CO2.exclude.mediumgain', 'bool')
        self.error_inflate = self.rcf.get('satellite.GOSAT.ACOS.CO2.errors.inflate', 'bool')
        # The data in the files are already filtered, yay!
        self.bias_correct = self.rcf.get('satellite.GOSAT.ACOS.CO2.bias.corrected', 'bool')
        self.exclude_sunglint = self.rcf.get('satellite.GOSAT.ACOS.CO2.exclude.sunglint', 'bool')
        if self.error_inflate:
            self.binning_time = self.rcf.get('satellite.GOSAT.ACOS.CO2.binning.time', 'float')
            self.binning_length = self.rcf.get('satellite.GOSAT.ACOS.CO2.binning.length', 'float')
        self.split_period = self.rcf.get('output.satellite.split.period')
        self.sampling_strategy = self.rcf.get('satellite.GOSAT.ACOS.CO2.sampling.strategy', 'int')
        self.filename_dict = defaultdict(list)
        if self.split_period == 'm':
            ym_list = self.get_YM_tuples()
            for year, month in ym_list:
                month_end = calendar.monthrange(year, month)[1]
                for i in range(1, month_end+1):
                    file_name = os.path.join(self.input_folder, datetime(year,month,i).strftime(self.filename_pattern))
                    if os.path.exists(file_name):
                        self.filename_dict[(year,month)].append(file_name)
        elif self.split_period == 'd':
            ymd_list = self.get_YMD_tuples()
            for year, month, day in ymd_list:
                file_name = os.path.join(self.input_folder, datetime(year,month,day).strftime(self.filename_pattern))
                if os.path.exists(file_name):
                    self.filename_dict[(year,month,day)].append(file_name)
        self.assimilate = self.rcf.get('CO2.satellite.assimilate', 'bool', default=True)
        if not self.assimilate:
            self.unassim_mdm = self.rcf.get('CO2.satellite.unassimilate.mdm', 'float')

    def readPeriod(self, ymd_tuple):
        # ymd_tuple can either be (year, month) or (year, month, day)
        return_dict = defaultdict(list)
        if len(ymd_tuple) == 2:
            print 'Reading files for %04i-%02i... '%(ymd_tuple),
        elif len(ymd_tuple) == 3:
            print 'Reading files for %04i-%02i-%02i... '%(ymd_tuple),
        for file_name in self.filename_dict[ymd_tuple]:
            new_dict = self.readData(file_name)
            for key, value in new_dict.items():
                return_dict[key].extend(value)
        # Sort according to time and convert to arrays
        sort_order = np.argsort(return_dict['times'])
        for key, value in return_dict.items():
            return_dict[key] = np.array(value)[sort_order]
        n_obs = len(return_dict['acos_id'])
        if self.error_inflate and n_obs>0:
            from tm5_utils import sat_utils

            # convert sample times to floating point seconds
            sample_seconds = [t-return_dict['times'][0] for t in return_dict['times']]
            sample_seconds = [dt.days * 86400 + dt.seconds for dt in sample_seconds]
            # convert lats and lons to a nx2 array
            locations = np.zeros((n_obs,2), return_dict['latitude'].dtype)
            locations[:,0] = return_dict['latitude']
            locations[:,1] = return_dict['longitude']
            return_dict['Xco2_err'], num_samples = sat_utils.errorInflation(float64(sample_seconds), return_dict['Xco2_err'], \
                locations, self.binning_time, self.binning_length, 'b', 0.0)
            del sample_seconds, locations, num_samples
        print 'done!'
        return return_dict

    def __call__(self):
        for ymd_tuple in sorted(self.filename_dict.keys()):
            self.writePeriod(ymd_tuple)
        # Now also write the bias correction parameter files in the input folder, but use different file names,
        # in case there are more parameters to optimize, such as the initial concentration, etc.
        if not self.bias_correct:
            self.writeBCparams()

class createACOSInputFile_b34_r03(createACOSInputFile_b34):
    """
    Parse ACOS b3.4/r03 'lite' files to create input files to be used by TM5
    """

    def readData(self, file_name):
        """
        Given a file_name, read the following fields and return them in a dictionary:

        /Retrieval/surface_type (0 for ocean, 1 for land)
        /Retrieval/psurf (retrieved surface pressure in hPa)
        /Retrieval/xco2_raw (XCO2 before bias correction)
        /Retrieval/SigmaB_Coeffiecient (psurf * these = pressure levels) (yes, the spelling mistake is there) (first element is TOA)
        /Retrieval/xco2_final (bias-corrected XCO2)
        /Retrieval/xco2_uncert (XCO2 error)
        /Retrieval/xco2_apriori (prior XCO2)
        /Retrieval/pwf (pressure weighting function)
        /Retrieval/xco2_ak (CO2 averaging kernel)
        /Sounding/latitude
        /Sounding/longitude
        /Sounding/sounding_id
        /Sounding/time (upto milliseconds)
        /Sounding/gain

        The following are needed for bias correction of land H/M:
        /Retrieval/albedo_3
        /Retrieval/delta_grad_co2
        /Retrieval/fs

        And for ocean:
        /Retrieval/s32
        /Retrieval/b1offset

        Chris O'Dell says we don't need to filter data, since that has been done by him already.
        """
        return_dict = {}
        with Dataset(file_name, 'r') as fid:
            # First check if some soundings should be thrown out
            n_obs = len(fid.dimensions['soundings_dimension'])
            n_lev = len(fid.dimensions['levels_dimension'])
            indices = np.ones((n_obs,2), bool)
            if self.exclude_medgain:
                indices[:,0] = fid.groups['Sounding'].variables['gain'][:] == 'H'
            if self.exclude_sunglint:
                indices[:,1] = fid.groups['Retrieval'].variables['surface_type'][:] == 1
            indices = all(indices, axis=1)
            n_obs = sum(indices)
            return_dict['latitude'] = fid.groups['Sounding'].variables['latitude'][:][indices]
            return_dict['longitude'] = fid.groups['Sounding'].variables['longitude'][:][indices]
            return_dict['times'] = np.array([datetime(*d) for d in fid.groups['Sounding'].variables['time'][:][indices]])
            return_dict['acos_id'] = fid.groups['Sounding'].variables['sounding_id'][:][indices]

            # For all soundings, the pressure boundaries are psurf * sigma
            return_dict['surf_pres'] = 100.0 * fid.groups['Retrieval'].variables['psurf'][:][indices] # in Pa
            sigma_pres = fid.groups['Retrieval'].variables['SigmaB_Coeffiecient'][:]
            return_dict['pres_levels'] = outer(return_dict['surf_pres'], sigma_pres)

            # Classify
            return_dict['gain'] = np.where(fid.groups['Sounding'].variables['gain'][:][indices] == 'H', 1, 0) # 1 is high gain
            return_dict['surface'] = fid.groups['Retrieval'].variables['surface_type'][:][indices] # 1 is land

            if self.bias_correct:
                return_dict['Xco2'] = fid.groups['Retrieval'].variables['xco2_final'][:][indices]
            else:
                return_dict['Xco2'] = fid.groups['Retrieval'].variables['xco2_raw'][:][indices]
            return_dict['Xco2_err'] = fid.groups['Retrieval'].variables['xco2_uncert'][:][indices]

            # The averaging kernel, prior profile and pressure weighting function have to be re-gridded to mid-levels
            raw_co2_apri = fid.groups['Retrieval'].variables['co2_profile_apriori'][:][indices]
            raw_ak = fid.groups['Retrieval'].variables['xco2_ak'][:][indices]
            raw_pwf = fid.groups['Retrieval'].variables['pwf'][:][indices]
            mid_levels = 0.5*(return_dict['pres_levels'][:,1:] + return_dict['pres_levels'][:,:-1])
            layer_thickness = np.diff(return_dict['pres_levels'], axis=1)
            return_dict['avg_kernel'] = np.zeros((n_obs, n_lev-1), np.float64)
            return_dict['co2_apri'] = np.zeros((n_obs, n_lev-1), np.float64)
            for i_obs in range(n_obs):
                ak_norm = raw_ak[i_obs]/raw_pwf[i_obs] # still 20 values
                spl = interpolate.InterpolatedUnivariateSpline(return_dict['pres_levels'][i_obs], ak_norm, k=1) # same as UnivariateSpline with smoothing forced to 0
                ak_new = spl(mid_levels[i_obs])
                pwf_new = layer_thickness[i_obs]/return_dict['surf_pres'][i_obs]
                ak_new = ak_new * pwf_new
                return_dict['avg_kernel'][i_obs] = ak_new
                spl = interpolate.InterpolatedUnivariateSpline(return_dict['pres_levels'][i_obs], raw_co2_apri[i_obs], k=1)
                return_dict['co2_apri'][i_obs] = spl(mid_levels[i_obs])

            # Now the variables for bias correction
            return_dict['albedo_3'] = fid.groups['Retrieval'].variables['albedo_3'][:][indices]
            return_dict['delta_grad_co2'] = fid.groups['Retrieval'].variables['delta_grad_co2'][:][indices]
            return_dict['Fs'] = fid.groups['Retrieval'].variables['fs'][:][indices]
            return_dict['s32'] = fid.groups['Retrieval'].variables['s32'][:][indices]
            return_dict['b1offset'] = fid.groups['Retrieval'].variables['b1offset'][:][indices]

        return return_dict

    def writePeriod(self, ymd_tuple):
        data_dict = self.readPeriod(ymd_tuple)

        n_obs = data_dict['acos_id'].shape[0]
        n_lev = data_dict['pres_levels'].shape[1]
        n_lay = data_dict['co2_apri'].shape[1]

        self.SatObservation[self.tracer]['dimensions'] = {'n_obs': n_obs, 'n_lev': n_lev, 'n_lay': n_lay} # idate already defined
        self.SatObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.SatObservation[self.tracer]['variable_attrs'] = OrderedDict()
        self.SatObservation[self.tracer]['file_attrs'] = {}

        self.SatObservation[self.tracer]['cdate'] = np.array([d.timetuple()[:6] for d in data_dict['times']], np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['cdate'] = ('n_obs', 'idate')
        self.SatObservation[self.tracer]['variable_attrs']['cdate'] = [('description', 'Observation date and time')]

        self.SatObservation[self.tracer]['latitude'] = np.float32(data_dict['latitude'])
        self.SatObservation[self.tracer]['variable_shapes']['latitude'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['latitude'] = [('unit', 'degrees north')]

        self.SatObservation[self.tracer]['longitude'] = np.float32(data_dict['longitude'])
        self.SatObservation[self.tracer]['variable_shapes']['longitude'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['lonitude'] = [('unit', 'degrees east')]

        self.SatObservation[self.tracer]['acos_id'] = int64(data_dict['acos_id'])
        self.SatObservation[self.tracer]['variable_shapes']['acos_id'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['acos_id'] = [('unit', 'YYYYMMDDhhmmss'), \
            ('description', 'ACOS sounding ID, constructed from scan start time in UTC')]

        self.SatObservation[self.tracer]['sampling_strategy'] = self.sampling_strategy * np.ones(n_obs, np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sampling_strategy'] = [('comment', '3 => symmetric sampling')]

        self.SatObservation[self.tracer]['gain'] = np.int8(data_dict['gain'])
        self.SatObservation[self.tracer]['variable_shapes']['gain'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['gain'] = [('comment', '1 => high gain, 0 => medium gain')]

        self.SatObservation[self.tracer]['surface'] = np.int8(data_dict['surface'])
        self.SatObservation[self.tracer]['variable_shapes']['surface'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['surface'] = [('comment', '1 => land, 0 => ocean')]

        self.SatObservation[self.tracer]['p_surf'] = np.float32(data_dict['surf_pres'])
        self.SatObservation[self.tracer]['variable_shapes']['p_surf'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['p_surf'] = [('unit', 'Pascal'), \
            ('description', 'Retrieved surface pressure')]

        self.SatObservation[self.tracer]['p_levels'] = np.float32(data_dict['pres_levels'])
        self.SatObservation[self.tracer]['variable_shapes']['p_levels'] = ('n_obs', 'n_lev')
        self.SatObservation[self.tracer]['variable_attrs']['p_levels'] = [('unit', 'Pascal'), \
            ('comment', 'Pressure boundaries for L2 retrieval, TOA first, ground last')]

        self.SatObservation[self.tracer]['column_mixing'] = np.float64(data_dict['Xco2'])
        self.SatObservation[self.tracer]['variable_shapes']['column_mixing'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['column_mixing'] = [('unit', 'Dry air mole fraction in parts per million')]
        if self.bias_correct:
            self.SatObservation[self.tracer]['variable_attrs']['column_mixing'].append(('comment', 'bias corrected'))
        else:
            self.SatObservation[self.tracer]['variable_attrs']['column_mixing'].append(('comment', 'not bias corrected'))

        if self.assimilate:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = np.float64(data_dict['Xco2_err'])
        else:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('comment', 'These data will not be assimilated')]
        self.SatObservation[self.tracer]['variable_shapes']['sigma_column_mixing'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('unit', 'Dry air mole fraction in parts per million')]

        self.SatObservation[self.tracer]['prior_mixing'] = np.float32(data_dict['co2_apri'])
        self.SatObservation[self.tracer]['variable_shapes']['prior_mixing'] = ('n_obs', 'n_lay')
        self.SatObservation[self.tracer]['variable_attrs']['prior_mixing'] = [('unit', 'Dry air mole fraction in parts per million'), \
            ('comment', 'First index is the TOA, last index is the surface')]

        self.SatObservation[self.tracer]['avg_kernel'] = np.float32(data_dict['avg_kernel'])
        self.SatObservation[self.tracer]['variable_shapes']['avg_kernel'] = ('n_obs', 'n_lay')
        self.SatObservation[self.tracer]['variable_attrs']['avg_kernel'] = [('description', 'Column averaging kernel'), \
            ('comment', 'First index is the TOA, last index is the surface')]

        self.SatObservation[self.tracer]['albedo_3'] = np.float32(data_dict['albedo_3'])
        self.SatObservation[self.tracer]['variable_shapes']['albedo_3'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['albedo_3'] = [('description', 'Retrieved Band 1 (2.04 micron) surface albedo')]

        self.SatObservation[self.tracer]['delta_grad_co2'] = np.float32(data_dict['delta_grad_co2'])
        self.SatObservation[self.tracer]['variable_shapes']['delta_grad_co2'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['delta_grad_co2'] = [('description', 'Change in CO2 vertical gradient (surface minus level 13), retrieved-apriori'), \
            ('comment', 'level 13 is at P/Psurf=0.631579')]

        self.SatObservation[self.tracer]['s32'] = np.float32(data_dict['s32'])
        self.SatObservation[self.tracer]['variable_shapes']['s32'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['s32'] = [('description', 'Ratio of Band 3 to Band 2 signal level')]

        self.SatObservation[self.tracer]['Fs'] = np.float32(data_dict['Fs'])
        self.SatObservation[self.tracer]['variable_shapes']['Fs'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['Fs'] = [('unit', 'W/m^2/micron/sr'), \
            ('description', 'Retrieved Fluorescence at 757 nm'), \
            ('comment', 'Simultaneous Fluorescence retrieval (at 757 nm) by L2 code; note this is different than the dedicated retrieval using only solar lines')]

        self.SatObservation[self.tracer]['b1offset'] = np.float32(data_dict['b1offset'])
        self.SatObservation[self.tracer]['variable_shapes']['b1offset'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['b1offset'] = [('unit', 'GOSAT radiance units * 1e7'), \
            ('description', 'Retrieved Band 1 (0.76 micron) radiance offset')]

        self.SatObservation[self.tracer]['file_attrs']['creation_time'] = datetime.now().strftime("%c")
        self.SatObservation[self.tracer]['file_attrs']['bias_correction_land_H'] = "XCO2_Bias_Corrected = XCO2_Raw + 9.5*(albedo_3-0.17) + 0.014*(delta_grad_co2-40) + 0.88*(Fs+0.04) - 0.10 ppm"
        self.SatObservation[self.tracer]['file_attrs']['bias_correction_land_M'] = "XCO2_Bias_Corrected = XCO2_Raw + (albedo_3-0.50)*4.4 + (delta_grad_co2-40)*0.0164 + 0.90 ppm"
        self.SatObservation[self.tracer]['file_attrs']['bias_correction_ocean'] = "XCO2_Bias_Corrected = XCO2_Raw + 0.62*(b1offset+1.0) - 42*(s32-0.61) - 1.0 ppm"
        self.SatObservation[self.tracer]['file_attrs']['param_errors_land_H'] = "albedo_3: 1.0, delta_grad_co2: 0.002, Fs: 0.1, offset: 0.2"
        self.SatObservation[self.tracer]['file_attrs']['param_errors_land_M'] = "albedo_3: 2.0, delta_grad_co2: 0.0015, offset: 0.4"
        self.SatObservation[self.tracer]['file_attrs']['param_errors_ocean'] = "s32: 8, b1offset: 0.1, offset: 0.25"

    def writeBCparams(self):
        num_params = 10
        corr_file_name = self.rcf.get('parameter.correlation.filename')
        corr_file_name_base = 'acos-3.4.r03-' + os.path.basename(corr_file_name)
        corr_file_name = os.path.join(os.path.dirname(corr_file_name), corr_file_name_base)
        checkDir(corr_file_name)
        savetxt(corr_file_name, eye(num_params), fmt='%5.1f')

        apri_param_std = self.rcf.get('parameter.prior.error.filename')
        apri_param_std_base = 'acos-3.4.r03-' + os.path.basename(apri_param_std)
        apri_param_std = os.path.join(os.path.dirname(apri_param_std), apri_param_std_base)
        param_std_values = np.array([1, 0.002, 0.1, 0.2, 2, 0.0015, 0.4, 0.1, 8, 0.25], np.float32)
        checkDir(apri_param_std)
        savetxt(apri_param_std, param_std_values, fmt='%10.5f')

        apri_param = self.rcf.get('parameter.apri.filename')
        apri_param_base = 'acos-3.4.r03-' + os.path.basename(apri_param)
        apri_param = os.path.join(os.path.dirname(apri_param), apri_param_base)
        param_values = np.array([9.5, 0.014, 0.88, -0.1, 4.4, 0.0164, 0.9, 0.62, -42, -1], np.float32)
        checkDir(apri_param)
        savetxt(apri_param, param_values, fmt='%10.5f')

class TCCON(Observations):

    def __init__(self, *args, **kwargs):
        super(TCCON, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.split_period = self.rcf.get('output.satellite.split.period')
        self.sampling_strategy = self.rcf.get('TCCON.CO2.sampling.strategy', 'int')
        self.assimilate = self.rcf.get('TCCON.CO2.assimilate', 'bool', default=False)
        self.unassim_mdm = self.rcf.get('CO2.satellite.unassimilate.mdm', 'float')
        self.tccon_files_dir = self.rcf.get('TCCON.CO2.folder')
        # In case we're assimilating TCCON, we should probably only assimilate close-to-zenith soundings
        self.tccon_max_sza = self.rcf.get('TCCON.CO2.max.sza', 'float', default=90.0)

    def __call__(self):
        if self.split_period == 'm':
            ymd_tuple = self.StartDate.timetuple()[:2]
        elif self.split_period == 'd':
            ymd_tuple = self.StartDate.timetuple()[:3]
        else:
            raise ValueError("Invalid split period specified for satellite data")

        data_dict, stat_dict = self.readPeriod(ymd_tuple)

        n_obs = len(data_dict['Xco2'])
        # maintain a dictionary of dimensions to write and variable shapes
        self.SatObservation[self.tracer]['dimensions'] = {'n_obs': n_obs}
        if n_obs == 0:
            return

        n_lev = data_dict['prior_pres_levels'].shape[1]
        self.SatObservation[self.tracer]['dimensions']['n_lev'] = n_lev

        self.SatObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.SatObservation[self.tracer]['variable_attrs'] = OrderedDict()
        self.SatObservation[self.tracer]['file_attrs'] = {}

        self.SatObservation[self.tracer]['file_attrs']['creation_time'] = datetime.now().strftime("%c")

        self.SatObservation[self.tracer]['cdate'] = np.array([d.timetuple()[:6] for d in data_dict['times']], np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['cdate'] = ('n_obs', 'idate')
        self.SatObservation[self.tracer]['variable_attrs']['cdate'] = [('description', 'Sounding date and time')]

        self.SatObservation[self.tracer]['latitude'] = np.float32(data_dict['latitude'])
        self.SatObservation[self.tracer]['variable_shapes']['latitude'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['latitude'] = [('unit', 'degrees north')]

        self.SatObservation[self.tracer]['longitude'] = np.float32(data_dict['longitude'])
        self.SatObservation[self.tracer]['variable_shapes']['longitude'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['lonitude'] = [('unit', 'degrees east')]

        self.SatObservation[self.tracer]['sampling_strategy'] = self.sampling_strategy * np.ones(n_obs, np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sampling_strategy'] = [('comment', '3 => symmetric sampling')]

        self.SatObservation[self.tracer]['p_surf'] = 100.0 * np.float32(data_dict['p_surf'])
        self.SatObservation[self.tracer]['variable_shapes']['p_surf'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['p_surf'] = [('unit', 'Pascals'), \
            ('description', 'Surface pressure from TCCON (pout_hpa)')]

        self.SatObservation[self.tracer]['p_levels_prior'] = 100.0 * np.float32(data_dict['prior_pres_levels'])
        self.SatObservation[self.tracer]['variable_shapes']['p_levels_prior'] = ('n_obs', 'n_lev')
        self.SatObservation[self.tracer]['variable_attrs']['p_levels_prior'] = [('unit', 'Pascals'), \
            ('comment', 'Pressure boundaries for TCCON prior, ground first, TOA last')]

        self.SatObservation[self.tracer]['column_mixing'] = np.float64(data_dict['Xco2'])
        self.SatObservation[self.tracer]['variable_shapes']['column_mixing'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['column_mixing'] = [('unit', 'Dry air mole fraction in parts per million')]

        if self.assimilate:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = np.float64(data_dict['Xco2_err'])
        else:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('comment', 'These data will not be assimilated')]
        self.SatObservation[self.tracer]['variable_shapes']['sigma_column_mixing'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('unit', 'Dry air mole fraction in parts per million')]

        self.SatObservation[self.tracer]['prior_mixing'] = np.float32(data_dict['co2_apri'])
        self.SatObservation[self.tracer]['variable_shapes']['prior_mixing'] = ('n_obs', 'n_lev')
        self.SatObservation[self.tracer]['variable_attrs']['prior_mixing'] = [('unit', 'Dry air mole fraction in parts per million'), \
            ('comment', 'First index is the surface layer, last index is TOA')]

        self.SatObservation[self.tracer]['avg_kernel'] = np.float32(data_dict['avg_kernel'])
        self.SatObservation[self.tracer]['variable_shapes']['avg_kernel'] = ('n_obs', 'n_lev')
        self.SatObservation[self.tracer]['variable_attrs']['avg_kernel'] = [('description', 'Column averaging kernel'), \
            ('comment', 'First index is the surface layer, last index is TOA')]

        self.SatObservation[self.tracer]['p_levels_ak'] = 100.0 * np.float32(data_dict['ak_pres_levels'])
        if self.SatObservation[self.tracer]['p_levels_ak'].ndim == 1:
            self.SatObservation[self.tracer]['variable_shapes']['p_levels_ak'] = ('n_lev',)
        elif self.SatObservation[self.tracer]['p_levels_ak'].ndim == 2:
            self.SatObservation[self.tracer]['variable_shapes']['p_levels_ak'] = ('n_obs', 'n_lev')
        self.SatObservation[self.tracer]['variable_attrs']['p_levels_ak'] = [('unit', 'Pascals'), \
            ('comment', 'Pressure boundaries for TCCON averaging kernel, ground first, TOA last')]

        self.SatObservation[self.tracer]['sza'] = np.float32(data_dict['sza'])
        self.SatObservation[self.tracer]['variable_shapes']['sza'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sza'] = [('description', 'Solar zenith angle')]

        self.SatObservation[self.tracer]['station_id'] = np.int8(data_dict['station_id'])
        self.SatObservation[self.tracer]['variable_shapes']['station_id'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['station_id'] = [('description', 'Mapping from station ID to name'), \
            ('num_stations', np.int32(len(stat_dict.keys())))]
        for k,v in stat_dict.items():
            self.SatObservation[self.tracer]['variable_attrs']['station_id'].append(('station_%03i'%k, v))

        if 'obs_num' in data_dict:
            self.SatObservation[self.tracer]['obs_num'] = np.int32(data_dict['obs_num'])
            self.SatObservation[self.tracer]['variable_shapes']['obs_num'] = ('n_obs',)
            self.SatObservation[self.tracer]['variable_attrs']['obs_num'] = [('description', 'Generic identifier for each sounding')]

class TCCON_GGG_2014(TCCON):

    def readPeriod(self, ymd_tuple):
        return_data_dict = defaultdict(list)
        station_indices = OrderedDict()

        # get a list of files
        tccon_files = glob(os.path.join(self.tccon_files_dir, '*.public.nc'))

        # make a mapping from file/station name to an index
        stat_idx = 0

        for file_name in tccon_files:
            new_dict = self.readData(file_name, ymd_tuple)

            var_names = new_dict['var names']
            for var_name in var_names:
                return_data_dict[var_name].extend(new_dict[var_name])

            # add station index
            if len(var_names) > 0:
                with Dataset(file_name, 'r') as fid:
                    stat_name = fid.longName
                station_indices[stat_idx] = stat_name
                n_obs = len(new_dict['Xco2'])
                stat_ids = stat_idx * np.ones(n_obs, np.int8)
                return_data_dict['station_id'].extend(stat_ids)
                stat_idx += 1

        sort_order = np.argsort(return_data_dict['times'])

        for var_name, var_value in return_data_dict.items():
            return_data_dict[var_name] = np.array(var_value)[sort_order]

        return return_data_dict, station_indices

    def readData(self, file_name, ymd_tuple):
        return_dict = {}
        return_dict['var names'] = []

        with Dataset(file_name, 'r') as fid:
            # filter by the correct time indices
            utc_times = np.array([datetime(1970,1,1) + timedelta(days=d) for d in fid.variables['time'][:]])
            year = np.array([u.year for u in utc_times])
            month = np.array([u.month for u in utc_times])
            valid_idx = np.logical_and(year == ymd_tuple[0], month == ymd_tuple[1])
            if len(ymd_tuple) == 3:
                day = np.array([u.day for u in utc_times])
                valid_idx = np.logical_and(valid_idx, day == ymd_tuple[2])

            n_obs = valid_idx.sum()
            if n_obs > 0:

                return_dict['times'] = utc_times[valid_idx]
                return_dict['latitude'] = fid.variables['lat_deg'][:][valid_idx]
                return_dict['longitude'] = fid.variables['long_deg'][:][valid_idx]
                return_dict['var names'].extend(['times', 'latitude', 'longitude'])

                return_dict['Xco2'] = fid.variables['xco2_ppm'][:][valid_idx]
                # We may or may not assimilate a TCCON sounding, depending on the SZA
                sza = fid.variables['asza_deg'][:][valid_idx] # the actual SZAs for the soundings
                return_dict['Xco2_err'] = np.where(sza <= self.tccon_max_sza, fid.variables['xco2_ppm_error'][:][valid_idx], self.unassim_mdm)
                return_dict['var names'].extend(['Xco2', 'Xco2_err'])

                # Read the surface pressure at the TCCON site (from a met station at the site)
                return_dict['p_surf'] = fid.variables['pout_hPa'][:][valid_idx]
                return_dict['var names'].append('p_surf')

                # The TCCON prior profiles and averaging kernels are specified on 71 pressure edges. They need to be cut off
                # near the ground depending on Pout. This will be done by Satellite.py/ApplyAveragingKernel, so no need to
                # do the interpolation to mid-levels here.

                # We need to read pressure levels, averaging kernels and prior profiles. THe prior profiles are in the form of a
                # lookup table, which the actual prior profiles being fewer than the number of soundings.
                prior_day_idx = fid.variables['prior_date_index'][:][valid_idx]
                return_dict['co2_apri'] = fid.variables['prior_co2'][:][prior_day_idx] # n_obs x 71
                return_dict['prior_pres_levels'] = fid.variables['prior_Pressure'][:][prior_day_idx] # n_obs x 71, boundary values
                return_dict['var names'].extend(['prior_pres_levels', 'co2_apri'])

                # The averaging kernels are only sensitive to the solar zenith angle. So there's a lookup table mapping each sounding
                # to the appropriate averaging kernel. Also, we need to interpolate to the pressure levels of the prior.
                sza_array = fid.variables['ak_zenith'][:] # the bin end-points for AK definition
                sza_idx = sza_array.searchsorted(sza) # all values will be between 0 and 15 (inclusive)
                return_dict['ak_pres_levels'] = fid.variables['ak_P_hPa'][:] # 71
                return_dict['avg_kernel'] = transpose(fid.variables['ak_co2'][:][:,sza_idx]) # n_obs x 71
                return_dict['var names'].extend(['avg_kernel', 'ak_pres_levels'])

                return_dict['sza'] = sza
                return_dict['var names'].append('sza')

        print 'Read ', n_obs, ' soundings for ', ymd_tuple, ' from ', os.path.basename(file_name)

        return return_dict

class TCCON_from_file(TCCON):
    """
    Class to create observations from TCCON.skeleton.nc and TCCON.auxiliary.nc, which has been distributed to all modellers.
    """
    def readPeriod(self, ymd_tuple):
        skel_file_name = os.path.join(self.tccon_files_dir, 'TCCON.skeleton.nc')
        aux_file_name = os.path.join(self.tccon_files_dir, 'TCCON.auxiliary.nc')

        # read the station ID to station name mapping
        station_indices = OrderedDict()
        with Dataset(skel_file_name, 'r') as fid:
            max_site_id = fid.variables['site_id'][:].max()
            for site_id in range(1, max_site_id+1):
                site_name = getattr(fid.variables['site_id'], 'site_%03i'%site_id)
                station_indices[site_id] = site_name

        return_data_dict = self.readData(skel_file_name, aux_file_name, ymd_tuple)

        return return_data_dict, station_indices

    def readData(self, skel_file_name, aux_file_name, ymd_tuple):
        return_dict = defaultdict(list)

        with Dataset(skel_file_name, 'r') as sfid, Dataset(aux_file_name, 'r') as afid:
            # filter by the correct time indices
            utc_times = sfid.variables['itime'][:]
            valid_idx = np.logical_and(utc_times[:,0] == ymd_tuple[0], utc_times[:,1] == ymd_tuple[1])
            if len(ymd_tuple) == 3:
                valid_idx = np.logical_and(valid_idx, utc_times[:,2] == ymd_tuple[2])

            n_obs = valid_idx.sum()
            if n_obs > 0:

                return_dict['times'] = np.array([datetime(*d) for d in utc_times[valid_idx]])
                return_dict['latitude'] = sfid.variables['latitude'][:][valid_idx]
                return_dict['longitude'] = sfid.variables['longitude'][:][valid_idx]

                return_dict['Xco2'] = afid.variables['XCO2'][:][valid_idx]
                sza = afid.variables['sza'][:][valid_idx]
                return_dict['Xco2_err'] = np.where(sza <= self.tccon_max_sza, afid.variables['XCO2_err'][:][valid_idx], self.unassim_mdm)

                return_dict['p_surf'] = sfid.variables['p_surf'][:][valid_idx]

                return_dict['prior_pres_levels'] = afid.variables['prior_pres_levels'][:][valid_idx]
                return_dict['co2_apri'] = afid.variables['co2_apri'][:][valid_idx]

                return_dict['ak_pres_levels'] = afid.variables['ak_pres_levels'][:][valid_idx]
                return_dict['avg_kernel'] = afid.variables['avg_kernel'][:][valid_idx]

                return_dict['sza'] = sza
                return_dict['station_id'] = sfid.variables['site_id'][:][valid_idx]
                return_dict['obs_num'] = np.arange(len(afid.dimensions['n_obs']))[valid_idx]

        print 'Read ', n_obs, ' soundings for ', ymd_tuple

        return return_dict

class TCCON_Pre_release(TCCON):

    def readPeriod(self, ymd_tuple):
        return_data_dict = defaultdict(list)
        station_indices = OrderedDict()

        # get a list of files
        tccon_files = glob(os.path.join(self.tccon_files_dir, 'tccon_pub.*.h5'))

        # make a mapping from file/station name to an index
        stat_idx = 0

        for file_name in tccon_files:
            new_dict = self.readData(file_name, ymd_tuple)

            var_names = new_dict['var names']
            for var_name in var_names:
                return_data_dict[var_name].extend(new_dict[var_name])

            # add station index
            if len(var_names) > 0:
                stat_name = os.path.basename(file_name).split('.')[1].split('_')[0]
                station_indices[stat_idx] = stat_name
                n_obs = len(new_dict['Xco2'])
                stat_ids = stat_idx * np.ones(n_obs, np.int8)
                return_data_dict['station_id'].extend(stat_ids)
                stat_idx += 1

        sort_order = np.argsort(return_data_dict['times'])

        for var_name, var_value in return_data_dict.items():
            return_data_dict[var_name] = np.array(var_value)[sort_order]

        return return_data_dict, station_indices

    def readData(self, file_name, ymd_tuple):
        return_dict = {}
        return_dict['var names'] = []

        with h5py.File(file_name, 'r') as fid:
            # filter by the correct time indices
            frac_days_epoch = fid['/TCCON/time'][:]
            utc_times = np.array([datetime(1970,1,1) + timedelta(days=d) for d in frac_days_epoch])

            year = np.array([u.year for u in utc_times])
            month = np.array([u.month for u in utc_times])
            valid_idx = np.logical_and(year == ymd_tuple[0], month == ymd_tuple[1])
            if len(ymd_tuple) == 3:
                day = np.array([u.day for u in utc_times])
                valid_idx = np.logical_and(valid_idx, day == ymd_tuple[2])

            n_obs = valid_idx.sum()
            if n_obs > 0:

                return_dict['times'] = utc_times[valid_idx]
                return_dict['latitude'] = fid['/TCCON/lat'][:][valid_idx]
                return_dict['longitude'] = fid['/TCCON/long'][:][valid_idx]
                return_dict['var names'].extend(['times', 'latitude', 'longitude'])

                return_dict['Xco2'] = fid['/TCCON/xco2'][:][valid_idx]
                sza = fid['/TCCON/asza'][:][valid_idx] # the actual SZAs for the soundings
                return_dict['Xco2_err'] = np.where(sza <= self.tccon_max_sza, fid['/TCCON/xco2_error'][:][valid_idx], self.unassim_mdm)
                return_dict['var names'].extend(['Xco2', 'Xco2_err'])

                # Read the surface pressure at the TCCON site (is this retrieved or prior?)
                return_dict['p_surf'] = fid['/TCCON/pout'][:][valid_idx]
                return_dict['var names'].append('p_surf')

                # The TCCON prior profiles and averaging kernels are specified on 71 pressure edges. They need to be cut off
                # near the ground depending on Pout. This will be done by Satellite.py/ApplyAveragingKernel, so no need to
                # do the interpolation to mid-levels here.

                # We need to read pressure levels, averaging kernels and prior profiles. THe prior profiles are in the form of a
                # lookup table, which the actual prior profiles being fewer than the number of soundings.
                prior_day_idx = fid['/TCCON/prior_date_index'][:][valid_idx]
                return_dict['co2_apri'] = fid['/PRIOR/co2'][:][prior_day_idx] # n_obs x 71
                return_dict['prior_pres_levels'] = fid['/PRIOR/pressure'][:][prior_day_idx] # n_obs x 71, boundary values
                return_dict['var names'].extend(['prior_pres_levels', 'co2_apri'])

                # The averaging kernels are only sensitive to the solar zenith angle. So there's a lookup table mapping each sounding
                # to the appropriate averaging kernel. Also, we need to interpolate to the pressure levels of the prior.
                sza_array = fid['/AK/zenith'][:] # the bin end-points for AK definition
                sza_idx = sza_array.searchsorted(sza) # all values will be between 0 and 15 (inclusive)
                return_dict['ak_pres_levels'] = fid['/AK/pressure'][:][sza_idx] # n_obs x 71
                return_dict['avg_kernel'] = fid['/AK/co2'][:][sza_idx] # n_obs x 71
                return_dict['var names'].extend(['avg_kernel', 'ak_pres_levels'])

                return_dict['sza'] = sza
                return_dict['var names'].append('sza')

        print 'Read ', n_obs, ' soundings for ', ymd_tuple, ' from ', os.path.basename(file_name)

        return return_dict

class OCO2(Observations):

    def __init__(self, *args, **kwargs):
        super(OCO2, self).__init__(*args, **kwargs)
        self.tracer = 'CO2'
        self.split_period = self.rcf.get('output.satellite.split.period')
        self.sampling_strategy = self.rcf.get('satellite.OCO2.ACOS.CO2.sampling.strategy', 'int')
        self.assimilate = self.rcf.get('CO2.satellite.assimilate', 'bool', default=True)
        if not self.assimilate:
            self.unassim_mdm = self.rcf.get('CO2.satellite.unassimilate.mdm', 'float')
        self.inflate_errors = self.rcf.get('satellite.OCO2.ACOS.CO2.error.inflate', 'bool', default=False)
        if self.inflate_errors:
            self.err_bin_time = self.rcf.get('satellite.OCO2.ACOS.CO2.error.binning.time', 'float')
            self.err_bin_length = self.rcf.get('satellite.OCO2.ACOS.CO2.binning.length', 'float')
        # Get a list of all L2 lite files
        self.lite_files_by_date = self.makeFileList()
        # Now classify them by day or month, depending on self.split_period
        self.filename_dict = defaultdict(list)
        if self.split_period == 'm':
            ym_list = self.get_YM_tuples()
            for year, month in ym_list:
                month_end = calendar.monthrange(year, month)[1]
                for i in range(1, month_end+1):
                    this_day = datetime(year, month, i)
                    if this_day in self.lite_files_by_date:
                        self.filename_dict[(year,month)].extend(self.lite_files_by_date[this_day])
        elif self.split_period == 'd':
            ymd_list = self.get_YMD_tuples()
            for year, month, day in ymd_list:
                this_day = datetime(year, month, day)
                if this_day in self.lite_files_by_date:
                    self.filename_dict[(year,month,day)].extend(self.lite_files_by_date[this_day])

    def readPeriod(self, ymd_tuple):
        # ymd_tuple can either be (year, month) or (year, month, day)
        return_data_dict = defaultdict(list)
        return_num_dict = defaultdict(int)
        aux_var_attrs = {}

        if len(self.filename_dict[ymd_tuple]) > 0:

            for file_name in self.filename_dict[ymd_tuple]:
                new_dict = self.readData(file_name)

                var_names = new_dict['var names']
                for var_name in var_names:
                    return_data_dict[var_name].extend(new_dict[var_name])

                for obs_type, n_obs in new_dict['num obs'].items():
                    return_num_dict[obs_type] += n_obs

                aux_var_attrs = new_dict['aux var attrs']

            # Sort according to time and convert to arrays
            if 'sounding_id' in return_data_dict:
                sort_order = np.argsort(return_data_dict['sounding_id'])
            else:
                sort_order = np.argsort(return_data_dict['times'])
            for var_name in var_names:
                return_data_dict[var_name] = np.array(return_data_dict[var_name])[sort_order]

            # Get rid of NaNs and Infs from some key variables
            invalid_idx = np.zeros(len(return_data_dict['times']), bool)
            invalid_idx[:] = False
            key_variables = self.rcf.get('satellite.OCO2.ACOS.check.valid_vars').split()
            for key_var in key_variables:
                v = return_data_dict[key_var]
                if v.ndim == 2:
                    invalid_idx = np.logical_or.reduce((invalid_idx, np.any(np.isnan(v), axis=1), np.any(np.isinf(v), axis=1)))
                elif v.ndim == 1:
                    invalid_idx = np.logical_or.reduce((invalid_idx, np.isnan(v), np.isinf(v)))
                else:
                    raise RuntimeError('Do not know how to check validity of rank %i variable %s'%(v.ndim, key_var))

            if np.any(invalid_idx):
                print 'WARNING :: %i invalid value(s) in '%invalid_idx.sum(), ymd_tuple
                for k, v in return_data_dict.items():
                    return_data_dict[k] = v[-invalid_idx]

            n_obs = len(return_data_dict['times'])

            # Do we want to inflate the errors?
            if self.inflate_errors and n_obs > 0:
                from tm5_utils import sat_utils

                # which algorithm to use to inflate errors?
                err_algo = self.rcf.get('satellite.OCO2.ACOS.CO2.err_inflate.choice', default='b')
                err_floor = self.rcf.get('satellite.OCO2.ACOS.CO2.err_inflate.floor', 'float', default=0.0)

                # convert sample times to floating point seconds
                sample_seconds = np.array([(t-return_data_dict['times'][0]).total_seconds() for t in return_data_dict['times']], np.float64)
                # convert lats and lons to a nx2 array
                locations = np.zeros((n_obs,2), return_data_dict['latitude'].dtype)
                locations[:,0] = return_data_dict['latitude']
                locations[:,1] = return_data_dict['longitude']

                num_err_samples = np.zeros(n_obs, np.int32)
                xco2_err = np.zeros(n_obs, np.float64)
                xco2_err[:] = return_data_dict['Xco2_err'] # initialise with the original errors
                orig_err = np.zeros(n_obs, np.float64)
                orig_err[:] = return_data_dict['Xco2_err'] # store the original errors

                # have to do this by sounding mode and surface type
                oper_mode = return_data_dict['operation_mode'] # 0=Nadir, 1=Glint, 2=Target, 3=Transition to/from Target
                surf_type = return_data_dict['surface_type'] # 1 => land, 0 => ocean
                for om in [0,1,2,3]:
                    for st in [0,1]:
                        idx = np.logical_and(oper_mode == om, surf_type == st)
                        if idx.sum() > 0:
                            this_locs = locations[idx]
                            this_errs = orig_err[idx]
                            this_secs = sample_seconds[idx]
                            new_errs, num_samples = sat_utils.errorInflation(this_secs, this_errs, this_locs, self.err_bin_time, \
                                self.err_bin_length, err_algo, err_floor)

                            xco2_err[idx] = new_errs
                            num_err_samples[idx] = num_samples

                return_data_dict['Xco2_err'] = xco2_err
                return_data_dict['Xco2_err_orig'] = orig_err
                return_data_dict['err_inflate_nsamples'] = num_err_samples

            if len(ymd_tuple) == 2:
                print_string = 'Read %i obs for %04i-%02i'%((n_obs,) + ymd_tuple)
            elif len(ymd_tuple) == 3:
                print_string = 'Read %i obs for %04i-%02i-%02i'%((n_obs,) + ymd_tuple)
            print print_string

        return return_data_dict, return_num_dict, aux_var_attrs

    def __call__(self):
        if self.split_period == 'm':
            ymd_tuple = self.StartDate.timetuple()[:2]
        elif self.split_period == 'd':
            ymd_tuple = self.StartDate.timetuple()[:3]
        else:
            raise ValueError("Invalid split period specified for satellite data")

        data_dict, num_dict, aux_var_attrs = self.readPeriod(ymd_tuple)

        n_obs = len(data_dict['Xco2'])
        # maintain a dictionary of dimensions to write and variable shapes
        self.SatObservation[self.tracer]['dimensions'] = {'n_obs': n_obs}
        if n_obs == 0:
            return

        n_lev = data_dict['pres_levels'].shape[1]
        n_lay = data_dict['co2_apri'].shape[1]
        self.SatObservation[self.tracer]['dimensions']['n_lev'] = n_lev
        self.SatObservation[self.tracer]['dimensions']['n_lay'] = n_lay

        written_vars = []

        self.SatObservation[self.tracer]['variable_shapes'] = OrderedDict()
        self.SatObservation[self.tracer]['variable_attrs'] = OrderedDict()
        self.SatObservation[self.tracer]['file_attrs'] = {}

        self.SatObservation[self.tracer]['file_attrs']['creation_time'] = datetime.now().strftime("%c")
        for obs_type, obs_num in num_dict.items():
            self.SatObservation[self.tracer]['file_attrs']['%s observations'%obs_type] = np.int32(obs_num)

        self.SatObservation[self.tracer]['cdate'] = np.array([d.timetuple()[:6] for d in data_dict['times']], np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['cdate'] = ('n_obs', 'idate')
        self.SatObservation[self.tracer]['variable_attrs']['cdate'] = [('description', 'Observation date and time')]
        written_vars.append('times')

        self.SatObservation[self.tracer]['latitude'] = np.float32(data_dict['latitude'])
        self.SatObservation[self.tracer]['variable_shapes']['latitude'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['latitude'] = [('unit', 'degrees north')]
        written_vars.append('latitude')

        self.SatObservation[self.tracer]['longitude'] = np.float32(data_dict['longitude'])
        self.SatObservation[self.tracer]['variable_shapes']['longitude'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['lonitude'] = [('unit', 'degrees east')]
        written_vars.append('longitude')

        self.SatObservation[self.tracer]['sampling_strategy'] = self.sampling_strategy * np.ones(n_obs, np.int16)
        self.SatObservation[self.tracer]['variable_shapes']['sampling_strategy'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sampling_strategy'] = [('comment', '3 => symmetric sampling')]
        written_vars.append('sampling_strategy')

        self.SatObservation[self.tracer]['surface_type'] = np.int8(data_dict['surface_type'])
        self.SatObservation[self.tracer]['variable_shapes']['surface_type'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['surface_type'] = [('comment', '1 => land, 0 => ocean')]
        written_vars.append('surface_type')

        self.SatObservation[self.tracer]['operation_mode'] = np.int8(data_dict['operation_mode'])
        self.SatObservation[self.tracer]['variable_shapes']['operation_mode'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['operation_mode'] = \
            [('comment', '0=Nadir, 1=Glint, 2=Target, 3=Transition to/from Target')]
        written_vars.append('operation_mode')

        self.SatObservation[self.tracer]['p_surf'] = np.float32(data_dict['surf_pres'])
        self.SatObservation[self.tracer]['variable_shapes']['p_surf'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['p_surf'] = [('unit', 'Pascal'), \
            ('description', 'Retrieved surface pressure')]
        written_vars.append('surf_pres')

        self.SatObservation[self.tracer]['p_levels'] = np.float32(data_dict['pres_levels'])
        self.SatObservation[self.tracer]['variable_shapes']['p_levels'] = ('n_obs', 'n_lev')
        self.SatObservation[self.tracer]['variable_attrs']['p_levels'] = [('unit', 'Pascal'), \
            ('comment', 'Pressure boundaries for L2 retrieval, TOA first, ground last')]
        written_vars.append('pres_levels')

        self.SatObservation[self.tracer]['column_mixing'] = np.float64(data_dict['Xco2'])
        self.SatObservation[self.tracer]['variable_shapes']['column_mixing'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['column_mixing'] = [('unit', 'Dry air mole fraction in parts per million')]
        written_vars.append('Xco2')

        if self.assimilate:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = np.float64(data_dict['Xco2_err'])
        else:
            self.SatObservation[self.tracer]['sigma_column_mixing'] = self.unassim_mdm * np.ones(n_obs, np.float64)
            self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('comment', 'These data will not be assimilated')]
        self.SatObservation[self.tracer]['variable_shapes']['sigma_column_mixing'] = ('n_obs',)
        self.SatObservation[self.tracer]['variable_attrs']['sigma_column_mixing'] = [('unit', 'Dry air mole fraction in parts per million')]
        written_vars.append('Xco2_err')

        if self.inflate_errors:
            self.SatObservation[self.tracer]['XCO2_error_original'] = data_dict['Xco2_err_orig']
            self.SatObservation[self.tracer]['variable_shapes']['XCO2_error_original'] = ('n_obs',)
            self.SatObservation[self.tracer]['variable_attrs']['XCO2_error_original'] = \
                [('unit', 'Dry air mole fraction in parts per million'), ('comment', 'Pre-inflation XCO2 errors from ACOS')]
            written_vars.append('Xco2_err_orig')

            self.SatObservation[self.tracer]['num_samples_in_err_bin'] = data_dict['err_inflate_nsamples']
            self.SatObservation[self.tracer]['variable_shapes']['num_samples_in_err_bin'] = ('n_obs',)
            self.SatObservation[self.tracer]['variable_attrs']['num_samples_in_err_bin'] = [('comment', 'Number of samples averaged for error inflation')]
            written_vars.append('err_inflate_nsamples')

        self.SatObservation[self.tracer]['prior_mixing'] = np.float32(data_dict['co2_apri'])
        self.SatObservation[self.tracer]['variable_shapes']['prior_mixing'] = ('n_obs', 'n_lay')
        self.SatObservation[self.tracer]['variable_attrs']['prior_mixing'] = [('unit', 'Dry air mole fraction in parts per million'), \
            ('comment', 'First index is the TOA, last index is the surface')]
        written_vars.append('co2_apri')

        self.SatObservation[self.tracer]['avg_kernel'] = np.float32(data_dict['avg_kernel'])
        self.SatObservation[self.tracer]['variable_shapes']['avg_kernel'] = ('n_obs', 'n_lay')
        self.SatObservation[self.tracer]['variable_attrs']['avg_kernel'] = [('description', 'Column averaging kernel'), \
            ('comment', 'First index is the TOA, last index is the surface')]
        written_vars.append('avg_kernel')

        # Now the auxiliary variables
        for var_name, var_value in data_dict.items():
            if not var_name in written_vars:
                self.SatObservation[self.tracer][var_name] = var_value
                self.SatObservation[self.tracer]['variable_shapes'][var_name] = ('n_obs',)
                if var_name in aux_var_attrs:
                    self.SatObservation[self.tracer]['variable_attrs'][var_name] = []
                    for attr_name, attr_value in aux_var_attrs[var_name].items():
                        self.SatObservation[self.tracer]['variable_attrs'][var_name].append((attr_name, attr_value))

    def makeFileList(self):
        """
        OCO2 lite files are named, e.g., oco2_LtCO2_150909_B7101A_150914134431s.nc4, where the last component is not predictable.
        Therefore, we should read in all the file names then associate them with their dates. This way, we don't have to look for
        oco2_LtCO2_yymmdd_*.nc4 each time we want the soundings for a particular day.
        """
        input_folder = self.rcf.get('satellite.OCO2.ACOS.CO2.data.folder')
        file_pattern = self.rcf.get('satellite.OCO2.ACOS.CO2.inputfile.pattern')
        next_day_hack = self.rcf.get('satellite.OCO2.ACOS.next_day_hack', 'bool', default=False)
        all_files = glob(os.path.join(input_folder, file_pattern))
        ret_dict = defaultdict(list)
        for file_name in all_files:
            base_name = os.path.basename(file_name)
            date_str  = base_name.split('_')[2]
            date_time = datetime.strptime(date_str, '%y%m%d') # this works because python's default pivot year is 1969
            ret_dict[date_time].append(file_name)
            if next_day_hack:
                # OCO-2 classifies soundings by orbit. Once an orbit is considered to be on a certain day, then all soundings
                # in that orbit are put in the lite file for that day. However, the 'day' for an orbit is calculated as its
                # day of equator crossing. So it's possible for the UTC day to roll over once the platform has crossed the
                # equator. All this is to say that the lite file for, say, June 22 2015, can contain a few soundings from June
                # 23, 2015. In other words, to get all soundings for June 23, 2015, we need to look at the lite files from June 23
                # as well as June 22. So each file has to be added to the next date as well.
                ret_dict[date_time+timedelta(days=1)].append(file_name)
        return ret_dict

class OCO2_B7101(OCO2):

    def correctResidualFPBias(self, xco2, oper_mode, surf_type, footprints):
        """
        In the v7101 release of the lite files, there are some residual biases in the bias corrected XCO2 that are depdenent
        on the footprint.(Table 5 of the users' guide). For coding this up, remember that glint is oper_mode 1, nadir is
        oper_mode 0, and target is oper_mode 2; land is surf_type 1, ocean is surf_type 0. The dictionary fp_bias has keys
        (surf_type, oper_mode).
        """
        fp_bias = { \
            (1, 1): np.array([-0.13, -0.09, -0.01, -0.14, 0.00,  0.06, 0.09, 0.21]), \
            (1, 0): np.array([-0.17, -0.10, -0.02, -0.09, 0.04,  0.03, 0.13, 0.18]), \
            (1, 2): np.array([-0.17, -0.04,  0.02, -0.10, 0.05,  0.01, 0.07, 0.16]), \
            (0, 1): np.array([ 0.03,  0.00,  0.02,  0.00, 0.00, -0.08, 0.02, 0.01]) \
            }

        xco2_ret = np.zeros_like(xco2)
        xco2_ret[:] = xco2

        for (st, om) in fp_bias.keys():
            sounding_idx = np.logical_and(np.ma.masked_equal(surf_type, st).mask, np.ma.masked_equal(oper_mode, om).mask)
            # these are the soundings for this surface type and operation mode
            if sounding_idx.sum() > 0:
                fp = footprints[sounding_idx] - 1 # 1 subtracted to serve as array indices
                residual_bias = fp_bias[(st,om)][fp]
                xco2_ret[sounding_idx] = xco2[sounding_idx] - residual_bias

        return xco2_ret

    def readData(self, file_name):
        use_warn = not self.rcf.get('satellite.OCO2.ACOS.quality_flag', 'bool')
        xco2_varname = self.rcf.get('satellite.OCO2.ACOS.CO2.varname') # whether we want to assimilate XCO2 with or without bias correction

        return_dict = {}
        with Dataset(file_name, 'r') as fid:

            n_sound = len(fid.dimensions['sounding_id'])

            # First, we need to classify the soundings
            # From the variable comment: "0=ocean;1=land"
            surface_type = fid.groups['Retrieval'].variables['surface_type'][:]
            # From the variable comment: "OCO-2 Operation Mode: 0=Nadir, 1=Glint, 2=Target, 3=Transition to/from Target"
            oper_mode = fid.groups['Sounding'].variables['operation_mode'][:]
            # We will reject the transition data, because it has not been looked at in detail by the L2 team. Target mode
            # data will have its own warn level threshold. Nadir ocean data should be rejected. Both glint land and ocean
            # data will be used, but with their individual warn levels.

            # We need to distinguish six types of soundings
            nadir_land  = np.logical_and(surface_type == 1, oper_mode == 0) ; N_nadir_land  = nadir_land.sum()
            nadir_ocean = np.logical_and(surface_type == 0, oper_mode == 0) ; N_nadir_ocean = nadir_ocean.sum() # Annmarie says these are horrible
            glint_land  = np.logical_and(surface_type == 1, oper_mode == 1) ; N_glint_land  = glint_land.sum()
            glint_ocean = np.logical_and(surface_type == 0, oper_mode == 1) ; N_glint_ocean = glint_ocean.sum()
            target_mode = (oper_mode == 2) ; N_target     = target_mode.sum()
            trans_mode  = (oper_mode == 3) ; N_transition = trans_mode.sum()

            # each sounding must be one of the above six types, so their sum should be the total number of soundings
            if n_sound != N_nadir_land + N_nadir_ocean + N_glint_land + N_glint_ocean + N_target + N_transition:
                sys.stderr.write('There are %i land nadir, %i ocean nadir, %i land glint, %i ocean glint, %i target and %i transition soundings, which do not add up to %i'%\
                    (N_nadir_land, N_nadir_ocean, N_glint_land, N_glint_ocean, N_target, N_transition, n_sound))
                raise RuntimeError

            if not use_warn:
                valid_idx = (fid.variables['xco2_quality_flag'][:] == 0) # boolean array
            else:
                warn_levels = fid.variables['warn_level'][:] # integers from 0 to 19
                # 19 means use all data, -1 means reject all data. By default, we'll reject all target and transition mode
                # soundings as well as all ocean nadir soundings, and accept all others (land glint, land nadir, ocean glint).
                gl_thresh = self.rcf.get('satellite.OCO2.ACOS.glint.land.warn_thres', 'int', default=19)
                go_thresh = self.rcf.get('satellite.OCO2.ACOS.glint.ocean.warn_thres', 'int', default=19)
                nl_thresh = self.rcf.get('satellite.OCO2.ACOS.nadir.land.warn_thres', 'int', default=19)
                no_thresh = self.rcf.get('satellite.OCO2.ACOS.nadir.ocean.warn_thres', 'int', default=-1)
                tg_thresh = self.rcf.get('satellite.OCO2.ACOS.target.warn_thres', 'int', default=-1)
                tr_thresh = self.rcf.get('satellite.OCO2.ACOS.transition.warn_thres', 'int', default=-1)

                valid_idx = np.logical_or.reduce( ( \
                    np.logical_and(nadir_land,  warn_levels <= nl_thresh), \
                    np.logical_and(nadir_ocean, warn_levels <= no_thresh), \
                    np.logical_and(glint_land,  warn_levels <= gl_thresh), \
                    np.logical_and(glint_ocean, warn_levels <= go_thresh), \
                    np.logical_and(target_mode, warn_levels <= tg_thresh), \
                    np.logical_and(trans_mode,  warn_levels <= tr_thresh), \
                    ) )

            # We are looking at lite files from this day as well as the previous day. So we need to filter for time.
            all_times = np.array([datetime(1970,1,1)+timedelta(seconds=d) for d in fid.variables['time'][:]])
            time_idx = np.array([self.StartDate <= d and d < self.EndDate for d in all_times])
            valid_idx = np.logical_and(time_idx, valid_idx)

            # Recompute the number of different obs types from the valid ones
            N_nadir_land  = nadir_land[valid_idx].sum()
            N_nadir_ocean = nadir_ocean[valid_idx].sum()
            N_glint_land  = glint_land[valid_idx].sum()
            N_glint_ocean = glint_ocean[valid_idx].sum()
            N_target      = target_mode[valid_idx].sum()
            N_transition  = trans_mode[valid_idx].sum()
            n_sound       = valid_idx.sum()

            return_dict['num obs'] = {}
            return_dict['num obs']['nadir land'] = N_nadir_land ; return_dict['num obs']['nadir ocean'] = N_nadir_ocean
            return_dict['num obs']['glint land'] = N_glint_land ; return_dict['num obs']['glint ocean'] = N_glint_ocean
            return_dict['num obs']['target'] = N_target ; return_dict['num obs']['transition'] = N_transition
            return_dict['num obs']['total'] = n_sound

            return_dict['var names'] = []

            # now read the relevant variables
            return_dict['latitude'] = fid.variables['latitude'][:][valid_idx]
            return_dict['longitude'] = fid.variables['longitude'][:][valid_idx]
            return_dict['times'] = np.array([datetime(1970,1,1)+timedelta(seconds=d) for d in fid.variables['time'][:][valid_idx]])

            # We need to correct for residual FP-dependent bias
            xco2_prebias = fid.variables[xco2_varname][:][valid_idx]
            footprints = fid.groups['Sounding'].variables['footprint'][:][valid_idx]
            return_dict['Xco2'] = self.correctResidualFPBias(xco2_prebias, oper_mode[valid_idx], surface_type[valid_idx], footprints)

            return_dict['Xco2_err'] = fid.variables['xco2_uncertainty'][:][valid_idx]

            return_dict['var names'].extend(['latitude', 'longitude', 'times', 'Xco2', 'Xco2_err'])

            # pressure levels (in the file they are oriented TOA to surface)
            return_dict['pres_levels'] = 100.0 * fid.variables['pressure_levels'][:][valid_idx] # in Pa, from TOA to surface
            return_dict['surf_pres'] = 100.0 * fid.groups['Retrieval'].variables['psurf'][:][valid_idx]

            return_dict['var names'].extend(['pres_levels', 'surf_pres'])

            # The averaging kernel, prior profile and pressure weighting function have to be re-gridded to mid-levels
            raw_co2_apri = fid.variables['co2_profile_apriori'][:][valid_idx]
            raw_ak = fid.variables['xco2_averaging_kernel'][:][valid_idx]
            raw_pwf = fid.variables['pressure_weight'][:][valid_idx]
            mid_levels = 0.5*(return_dict['pres_levels'][:,1:] + return_dict['pres_levels'][:,:-1])
            layer_thickness = np.diff(return_dict['pres_levels'], axis=1)
            n_obs, n_lev = raw_ak.shape

            return_dict['avg_kernel'] = np.zeros((n_obs, n_lev-1), np.float64)
            return_dict['co2_apri'] = np.zeros((n_obs, n_lev-1), np.float64)
            # What Chris calls the normalized averaging kernel -- standard since ACOS B3.5 -- is different from ACOS B3.4
            # and below. If A is the formal 20x20 averaging kernel in mixing ratio space, then a_i = (A^T h)_i/h_i, where
            # h = \delta_P/Psurf is the pressure weighting vector. Therefore, to interpolate to mid-levels, we do not need
            # to divide by the pressure weighting function.
            for i_obs in range(n_obs):
                ak_norm = raw_ak[i_obs] # still 20 values
                spl = interpolate.InterpolatedUnivariateSpline(return_dict['pres_levels'][i_obs], ak_norm, k=1) # same as UnivariateSpline with smoothing forced to 0
                ak_new = spl(mid_levels[i_obs])
                pwf_new = layer_thickness[i_obs]/return_dict['surf_pres'][i_obs]
                return_dict['avg_kernel'][i_obs] = ak_new
                spl = interpolate.InterpolatedUnivariateSpline(return_dict['pres_levels'][i_obs], raw_co2_apri[i_obs], k=1)
                return_dict['co2_apri'][i_obs] = spl(mid_levels[i_obs])

            return_dict['var names'].extend(['avg_kernel', 'co2_apri'])

            # Now read some auxiliary variables for diagnosing problems
            aux_vars = self.rcf.get('satellite.OCO2.ACOS.auxiliary.variables').split()
            aux_var_attrs = {}
            for aux_var in aux_vars:
                try:
                    grp_name, var_name = aux_var.split('/')
                    var_obj = fid.groups[grp_name].variables[var_name]
                except ValueError:
                    var_name = aux_var
                    var_obj = fid.variables[aux_var]

                # If there are missing values, raise a warning
                var_values = var_obj[:][valid_idx]
                if isinstance(var_values, np.ma.masked_array):
                    return_dict[var_name] = var_values.data
                    if any(var_values.mask):
                        missing_value = getattr(var_obj, 'missing_value')
                        sys.stderr.write('WARNING :: variable %s in file %s has %i masked/missing values\n'%\
                            (var_name, file_name, var_values.mask.sum()))
                else:
                    return_dict[var_name] = var_values

                return_dict['var names'].append(var_name)

                aux_var_attrs[var_name] = {}
                for attr_name in var_obj.ncattrs():
                    aux_var_attrs[var_name][attr_name] = getattr(var_obj, attr_name)

            return_dict['aux var attrs'] = aux_var_attrs

        return return_dict
