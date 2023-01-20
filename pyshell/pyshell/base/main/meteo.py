#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True

import re, os, shutil, glob, hashlib, pwd, progressbar, calendar
from pyshell.tmflex import rc
from pyshell.base.helper.Utilities import *
from numpy import *
from datetime import datetime, timedelta
from tarfile import open as TarFile
from dateutil.relativedelta import relativedelta
from collections import defaultdict
from stat import ST_UID
import tempfile

class LevelDefinitions(object):
    """
    Both the InitialConcentration and ApplySatelliteObs classes use level definitions from
    ECMWF for their pressure boundaries. We define them here, to be derived by any class that
    wants to use them.
    """
    def __init__(self):
        self.rcf = rc.RcFile(os.environ['pyshell.rc'])
        self.AT_coeffs = {} # in Pascals
        self.BT_coeffs = {}
        self.AT_coeffs['ml60/tropo25'] = array([0., 7.367743, 210.39389, 855.361755, 2063.779785, 3850.91333, 6144.314941, 8802.356444999999, 11632.758789, 14411.124023, 16899.46875, 18864.75, 20097.402344, 20429.863281, 19755.109375, 18045.183594, 15379.805664, 12077.446289, 8765.053711, 6018.019531, 3960.291504, 1680.640259, 713.218079, 298.495789, 95.63696299999999, 0.])
        self.BT_coeffs['ml60/tropo25'] = array([1., 0.99401945, 0.97966272, 0.95182151, 0.90788388, 0.84737492, 0.77159661, 0.6832686100000001, 0.58616841, 0.48477158, 0.38389215, 0.28832296, 0.2024759, 0.13002251, 0.07353382999999999, 0.03412116, 0.01114291, 0.00181516, 0.00007582, 0., 0., 0., 0., 0., 0., 0.])
        self.AT_coeffs['ml60/tropo36'] = array([0.0, 0.0, 7.367743, 210.393890, 467.333588, 855.361755, 2063.779785, 2887.696533, 3850.913330, 6144.314941, 7438.803223, 8802.356445, 11632.758789, 13043.218750, 14411.124023, 16899.468750, 17961.357422, 18864.75, 20097.402344, 20384.480469, 20429.863281, 19755.109375, 19027.695313, 18045.183594, 15379.805664, 13775.325195, 12077.446289, 8765.053711, 6018.019531, 3960.291504, 2579.888672, 1680.640259, 713.218079, 464.618134, 298.495789, 95.636963, 0.0])
        self.BT_coeffs['ml60/tropo36'] = array([1.0, 0.997630, 0.994019, 0.979663, 0.967645, 0.951822, 0.907884, 0.879657, 0.847375, 0.771597, 0.728786, 0.683269, 0.586168, 0.535710, 0.484772, 0.383892, 0.335155, 0.288323, 0.202476, 0.164384, 0.130023, 0.073534, 0.051690, 0.034121, 0.011143, 0.005081, 0.001815, 0.000076, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.AT_coeffs['ml60/tropo60'] = array([0., 20., 38.425343, 63.647804, 95.636963, 134.483307, 180.584351, 234.779053, 298.495789, 373.971924, 464.618134, 575.651001, 713.218079, 883.660522, 1094.834717, 1356.474609, 1680.640259, 2082.273926, 2579.888672, 3196.421631, 3960.291504, 4906.708496, 6018.019531, 7306.631348, 8765.053711, 10376.126953, 12077.446289, 13775.325195, 15379.805664, 16819.474609, 18045.183594, 19027.695313, 19755.109375, 20222.205078, 20429.863281, 20384.480469, 20097.402344, 19584.330078, 18864.75, 17961.357422, 16899.468750, 15706.447266, 14411.124023, 13043.218750, 11632.758789, 10209.500977, 8802.356445, 7438.803223, 6144.314941, 4941.778320, 3850.913330, 2887.696533, 2063.779785, 1385.912598, 855.361755, 467.333588, 210.393890, 65.889244, 7.367743, 0., 0.])[::-1]
        self.BT_coeffs['ml60/tropo60'] = array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.00007582, 0.00046139, 0.00181516, 0.00508112, 0.01114291, 0.02067788, 0.03412116, 0.05169041, 0.07353383, 0.09967469, 0.13002251, 0.16438432, 0.20247590, 0.24393314, 0.28832296, 0.33515489, 0.38389215, 0.43396294, 0.48477158, 0.53570992, 0.58616841, 0.63554746, 0.68326861, 0.72878581, 0.77159661, 0.81125343, 0.84737492, 0.87965691, 0.90788388, 0.93194032, 0.95182151, 0.96764523, 0.97966272, 0.98827010, 0.99401945, 0.99763012, 1.])[::-1]
        self.AT_coeffs['ml60/ml60']    = self.AT_coeffs['ml60/tropo60']
        self.BT_coeffs['ml60/ml60']    = self.BT_coeffs['ml60/tropo60']
        self.AT_coeffs['ml91/tropo25'] = array([0., 6.575628, 162.043427, 895.193542, 2356.202637, 3743.464355, 6353.920898, 8356.252930000001, 11543.166992, 14665.645508, 16544.585938, 18798.822266, 20107.03125, 20434.158203, 19587.513672, 18191.029297, 15638.053711, 11982.662109, 8564.624023, 6199.839355, 3767.196045, 1713.709595, 701.813354, 271.356506, 108.715561, 0.])
        self.BT_coeffs['ml91/tropo25'] = array([1., 0.994204, 0.9822379999999999, 0.950274, 0.897767, 0.85095, 0.764679, 0.698224, 0.589317, 0.475016, 0.399205, 0.291993, 0.20152, 0.131935, 0.067316, 0.036227, 0.012508, 0.001701, 5.5e-05, 0., 0., 0., 0., 0., 0., 0.])
        self.AT_coeffs['ml91/tropo34'] = array([0., 6.575628, 336.772369, 1297.656128, 3010.146973, 5422.802734, 8356.252930000001, 11543.166992, 14665.645508, 17385.595703, 19348.775391, 20319.011719, 20348.916016, 19919.796875, 19184.544922, 18191.029297, 16990.623047, 15638.053711, 14192.009766, 12713.897461, 11262.484375, 9873.560546999999, 8564.624023, 7341.469727, 6199.839355, 4663.776367, 3358.425781, 2292.155518, 1463.16394, 857.945801, 450.685791, 204.637451, 76.16765599999999, 21.413612, 0.])
        self.BT_coeffs['ml91/tropo34'] = array([1., 0.994204, 0.9734660000000001, 0.935157, 0.875518, 0.795385, 0.698224, 0.589317, 0.475016, 0.362203, 0.259554, 0.176091, 0.112979, 0.080777, 0.055474, 0.036227, 0.022189, 0.012508, 0.006322, 0.002765, 0.001, 0.000279, 0.000055, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
        self.AT_coeffs['ml91/ml91']    = array([0.000000,      2.000040,      3.980832,      7.387186,     12.908319,     21.413612,     33.952858,     51.746601,     76.167656,    108.715561,\
                                              150.986023,    204.637451,    271.356506,    352.824493,    450.685791,    566.519226,    701.813354,    857.945801,   1036.166504,   1237.585449,\
                                             1463.163940,   1713.709595,   1989.874390,   2292.155518,   2620.898438,   2976.302246,   3358.425781,   3767.196045,   4202.416504,   4663.776367,\
                                             5150.859863,   5663.156250,   6199.839355,   6759.727051,   7341.469727,   7942.926270,   8564.624023,   9208.305664,   9873.560547,  10558.881836,\
                                            11262.484375,  11982.662109,  12713.897461,  13453.225586,  14192.009766,  14922.685547,  15638.053711,  16329.560547,  16990.623047,  17613.281250,\
                                            18191.029297,  18716.968750,  19184.544922,  19587.513672,  19919.796875,  20175.394531,  20348.916016,  20434.158203,  20426.218750,  20319.011719,\
                                            20107.031250,  19785.357422,  19348.775391,  18798.822266,  18141.296875,  17385.595703,  16544.585938,  15633.566406,  14665.645508,  13653.219727,\
                                            12608.383789,  11543.166992,  10471.310547,   9405.222656,   8356.252930,   7335.164551,   6353.920898,   5422.802734,   4550.215820,   3743.464355,\
                                             3010.146973,   2356.202637,   1784.854614,   1297.656128,    895.193542,    576.314148,    336.772369,    162.043427,     54.208336,      6.575628,\
                                                0.003160,      0.000000])[::-1]
        self.BT_coeffs['ml91/ml91']    = array([0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,\
                                                0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,\
                                                0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000000,\
                                                0.000000,      0.000000,      0.000000,      0.000000,      0.000000,      0.000014,      0.000055,      0.000131,      0.000279,      0.000548,\
                                                0.001000,      0.001701,      0.002765,      0.004267,      0.006322,      0.009035,      0.012508,      0.016860,      0.022189,      0.028610,\
                                                0.036227,      0.045146,      0.055474,      0.067316,      0.080777,      0.095964,      0.112979,      0.131935,      0.152934,      0.176091,\
                                                0.201520,      0.229315,      0.259554,      0.291993,      0.326329,      0.362203,      0.399205,      0.436906,      0.475016,      0.513280,\
                                                0.551458,      0.589317,      0.626559,      0.662934,      0.698224,      0.732224,      0.764679,      0.795385,      0.824185,      0.850950,\
                                                0.875518,      0.897767,      0.917651,      0.935157,      0.950274,      0.963007,      0.973466,      0.982238,      0.989153,      0.994204,\
                                                0.997630,      1.000000])[::-1]
        level_key = '/'.join([self.rcf.get('ECLEVS'), self.rcf.get('LEVS')])
        self.AT = self.AT_coeffs[level_key]
        self.BT = self.BT_coeffs[level_key]

class Meteo(object):

    def __init__(self, tm5Obj):
        self.startDate = tm5Obj.StartTime.date()
        self.endDate = tm5Obj.EndTime.date()
        self.rcf = rc.RcFile(os.environ['pyshell.rc'])
        self.MissingFiles = []
        self.permanent_scratch = self.rcf.get('my.meteo.source.dir')

    def getTarFiles(self):
        meteo_archive_dir = self.rcf.get('my.meteo.archive')
        prelim_file_list = dirWalk(['*.tar'], meteo_archive_dir)
        prelim_file_list = filter(lambda x: len(os.path.basename(x).split('_')) > 1, prelim_file_list) # filter out oro.tar and lsm.tar
        # filter by month
        # first, construct all the months as strings
        self.month_string_list = []
        curDate = self.startDate.replace(day=1)
        while curDate <= self.endDate.replace(day=1):
            self.month_string_list.append(curDate.strftime("%Y%m"))
            curDate = curDate + relativedelta(months=1)
        final_file_list = []
        for filename in prelim_file_list:
            if os.path.basename(filename).split('_')[1][:6] in self.month_string_list:
                final_file_list.append(filename)
        return final_file_list

    def decompressTarFile(self, tar_file):
        arc_file = TarFile(tar_file, 'r')
        members = arc_file.getnames()
        if self.overwrite:
            arc_file.extractall(path=self.permanent_scratch)
            print "Success: all files in %s extracted to %s"%(tar_file, self.permanent_scratch)
        else:
            for member in members:
                if (not os.path.exists(os.path.join(self.permanent_scratch, member))):
                    try:
                        arc_file.extract(member, path=self.permanent_scratch)
                        print "Success: extracted %s"%(os.path.join(self.permanent_scratch, member))
                    except:
                        print "Failure: error extracting %s"%(os.path.join(self.permanent_scratch, member))
                else:
                    print "Warning: %s exists and I wasn't asked to overwrite"%(os.path.join(self.permanent_scratch, member))
        arc_file.close()

    def necessaryFiles(self):
        # Get a list of glb100x100 files that are needed for this run
        basename_patterns = []
        resol_names = self.rcf.get('my.meteo.field.keys').split()
        for resol_name in resol_names:
            field_names = self.rcf.get('my.meteo.'+resol_name+'.names').split()
            temp_resol = self.rcf.get('my.meteo.'+resol_name+'.resol')
            for field_name in field_names:
                basename_patterns.append('ec-ei-' + resol_name + '-glb100x100-' + field_name + '_%Y%m%d_' + temp_resol + '.hdf')
        cur_date = self.startDate - timedelta(days=1)
        end_date = self.endDate + timedelta(days=1)
        ret_list = []
        while cur_date <= end_date:
            for base_patt in basename_patterns:
                ret_list.append(os.path.join(self.permanent_scratch, cur_date.strftime(base_patt)))
            cur_date = cur_date + timedelta(days=1)
        return ret_list

    def checkCoarsenedFiles(self):
        """
        For this run, if we're going to run with coarsened meteo, check that all the coarsened meteo are there.
        If not, print a summary of the meteo files that are missing.
        """
        # What are the 3D fields needed? Don't add 'diffus' just yet.
        # Also remember that 'sp' needs to exist one day before and after.
        check_fields = ['cld', 'convec', 'mfuv', 'mfw', 'q', 'sp', 'sub', 't']
        # File names are , e.g., sp_19981231_00p03.hdf
        region_names = self.rcf.get('regions').split()
        meteo_dir = self.rcf.get('my.meteo.dir') # /scratch2/portfolios/BMC/co2/Sourish.Basu/var4d-mt/tm5_meteo_out/ct_nam
        ecclass_ecl = self.rcf.get('my.ecclass_ecl') # ei
        mlevs = self.rcf.get('my.mlevs') # tropo25
        big_missing_dict = dict.fromkeys(region_names)
        for region in region_names:
            src_dir = os.path.join(meteo_dir, 'ec', ecclass_ecl, 'fc012up2tr3', mlevs, region)
            missing_dict = defaultdict(list)
            cur_date = self.startDate
            while cur_date < self.endDate:
                for field in check_fields:
                    file_name = cur_date.strftime("%Y/%m/") + field + cur_date.strftime("_%Y%m%d_00p03.hdf")
                    file_name = os.path.join(src_dir, file_name)
                    if not os.path.exists(file_name):
                        missing_dict[(cur_date.year, cur_date.month)].append(file_name)
                cur_date += timedelta(days=1)
            # now print out the missing files
            print "For region %s, the following files are missing"%region
            print "========================================================"
            for year, month in sorted(missing_dict.keys()):
                print "%04i/%02i"%(year,month)
                print "-------"
                for file_name in sorted(missing_dict[(year,month)]):
                    print file_name
                print
            big_missing_dict[region] = missing_dict
        return big_missing_dict

    def checkGlobalFiles(self):
        # Check whether all glb100x100 files are present, and return a list of missing files
        required_files = self.necessaryFiles()
        ret_list = []
        for fname in required_files:
            if not os.path.exists(fname):
                ret_list.append(fname)
        return ret_list

    def __call__(self, overwrite=False):
        self.overwrite = overwrite
        file_list = self.getTarFiles()
        for file_name in file_list:
            self.decompressTarFile(file_name)

class MoveMeteo(object):
    """
    With the move from huygens to cartesius, the permanent scratch with all the meteo files have a new
    folder. This is a good time to reorganize the files a bit, instead of having ~100,000 files in one
    folder. The idea is that, for example, a file called

    /projects/huygens_projects/tm5meteo/mcn1/FILESET/fileset_mt5meteo/ec-ei-fc012up2tr3-ml60-glb100x100-t_20080101_00p03.hdf

    will become

    /project/tm5meteo/tmm/ec/ei/fc012up2tr3/ml60/glb100x100/2008/01/t_20080101_00p03.hdf
    """
    def __init__(self):
        self.rcf = rc.RcFile(os.environ['pyshell.rc'])
        self.source_folder = '/projects/huygens_projects/tm5meteo/mcn1/FILESET/fileset_mt5meteo'
        self.target_folder = self.rcf.get('my.meteo.source.dir') # permanent scratch for 1x1 meteo
        self.meteo_archive = self.rcf.get('my.meteo.archive')
        # construct the dictionary of variables for glb100x100 files
        self.var_dict = dict.fromkeys(self.rcf.get('my.meteo.time.resolution').split(','))
        for tkey in self.var_dict:
            self.var_dict[tkey] = dict.fromkeys(self.rcf.get('my.meteo.grids.'+tkey).split(','))
            for grid in self.var_dict[tkey]:
                self.var_dict[tkey][grid] = self.rcf.get('my.meteo.vars.'+tkey+'.'+grid).split(',')
        self.coarsened_archive = '/archive/sbasu/meteo/tm5_meteo_out'
        self.coarsened_scratch = self.rcf.get('tmm.output.dir')
        self.coarsened_levels = self.rcf.get('LEVS')
        self.scratch_folder = self.rcf.get('my.scratch_dir')
        self.meteo_class = self.rcf.get('my.meteo.class')

    def list_glb100x100_hdf_files(self, start_date, end_date, mode='all'):
        # List all the glb100x100 hdf files that should be present in the permanent scratch, between two dates
        startdate = datetime(*start_date)
        finaldate = datetime(*end_date)
        this_date = startdate
        file_list = []
        missing_dict = defaultdict(list)
        while this_date < finaldate:
            for base_name in self.build_hdf_file_names(this_date):
                fname = os.path.join(self.target_folder, 'ec', self.meteo_class, base_name)
                if mode == 'all':
                    file_list.append(fname)
                elif mode == 'missing':
                    if not os.path.isfile(fname):
                        missing_dict[(this_date.year, this_date.month)].append(fname)
            this_date += timedelta(days=1)
        if mode == 'all':
            return file_list
        elif mode == 'missing':
            return missing_dict

    def list_glb100x100_tar_files(self, start_date, end_date, mode='all', select_vars=None):
        # List all the glb100x100 tar files that should be present in the meteo archive, between two dates
        startdate = datetime(*start_date)
        finaldate = datetime(*end_date)
        date = startdate
        file_list = []
        missing_dict = defaultdict(list)
        while date < finaldate:
            yy=date.year
            mm=date.month
            dd=date.day
            for tres in self.var_dict:
                for grid in self.var_dict[tres]:
                    if select_vars == None:
                        var_list = self.var_dict[tres][grid]
                    else:
                        var_list = set(select_vars).intersection(self.var_dict[tres][grid])
                    for var_name in var_list:
                        filenames = [os.path.join(self.meteo_archive, f) for f in self.build_var_string(yy,mm,dd,var_name,tres,grid)]
                        for fname in filenames:
                            if mode == 'all':
                                file_list.append(fname)
                            elif mode == 'missing':
                                if not os.path.isfile(fname):
                                    missing_dict[(yy,mm)].append(fname)
            date = self.nextmonth(date)
        if mode == 'all':
            return file_list
        elif mode == 'missing':
            return missing_dict

    def nextmonth(self,date):
        # add one month to date, reset day of month to 1
        if date.month == 12:
            newdate=datetime(date.year+1,1,1)
        else:
            newdate=datetime(date.year,date.month+1,1)
        return newdate

    def build_var_string(self,year,month,day,var,tres,grid):
        """
        File names are five types, for either EI or OD meteo:

        OD:
        ec:/nlh/TM/meteo/ec/od_L91/fc012up2tr3/ml91/glb100x100/q_20090805_00p03.tar
        ec:/nlh/TM/meteo/ec/od_L91/fc012up2tr3/ml91/glb100x100/tsp_200908_00p03.tar
        ec:/nlh/TM/meteo/ec/od_L91/fc012up2tr3/sfc/glb100x100/sd_200908_00p03.tar
        ec:/nlh/TM/meteo/ec/od_L91/an0tr6/sfc/glb100x100/veg_200908_00p06.tar
        ec:/nlh/TM/meteo/ec/od_L91/an0tr6/sfc/glb100x100/srols_201004.tar

        EI:
        ec:/nlh/TM/meteo/ec/ei/fc012up2tr3/ml60/glb100x100/cld_20120100_00p03.tar (split each five days --> cld,mfuv,mfw,q,sub,t)
        ec:/nlh/TM/meteo/ec/ei/fc012up2tr3/ml60/glb100x100/sp_201201_00p03.tar (one for each month --> sp,tsp)
        ec:/nlh/TM/meteo/ec/ei/fc012up2tr3/sfc/glb100x100/sf_201201_00p03.tar (--> blh,ci,cp,d2m,ewss,g10m,lsp,nsss,sd,sf,skt,slhf,src,sshf,ssr,ssrd,sst,str,strd,swvl1,t2m,u10m,v10m)
        ec:/nlh/TM/meteo/ec/ei/an0tr6/sfc/glb100x100/sr_201201_00p06.tar (with the '00p06' --> albedo,sr,veg)
        ec:/nlh/TM/meteo/ec/ei/an0tr6/sfc/glb100x100/srols_201201.tar (without any time resolution identifier --> srols)

        We only need to return the part fc012up2tr3/ml91/glb100x100/q_20090805_00p03.tar, for example.
        """
        if tres == '00p03':
            tkey = 'fc012up2tr3'
        elif tres == '00p06':
            tkey = 'an0tr6'

        file_names = []
        # For some variables, the monthly data are split over five days
        if var in ['cld', 'mfuv', 'mfw', 'q', 'sub', 't']:
            for i in range(0,calendar.monthrange(year,month)[1]+1,5):
                # sometimes we want to start in the middle of the month, in which case there is no need to build a list of all files
                if day-i < 5:
                    filename = os.path.join(tkey, grid, 'glb100x100', "%s_%04i%02i%02i_%s.tar"%(var,year,month,i,tres))
                    file_names.append(filename)
        else:
            if var in ['srols']:
                filename = os.path.join(tkey, grid, 'glb100x100', "%s_%04i%02i.tar"%(var,year,month))
                file_names.append(filename)
            else:
                filename = os.path.join(tkey, grid, 'glb100x100', "%s_%04i%02i_%s.tar"%(var,year,month,tres))
                file_names.append(filename)

        return file_names

    def build_hdf_file_names(self, date_obj):
        """
        File names are of the following four types, for both EI and OD meteo:

        OD:
        an0tr6/sfc/glb100x100/2006/11/srols_200611.hdf --> one file per month
        an0tr6/sfc/glb100x100/2006/11/albedo_20061104_00p06.hdf --> one file per day, variables are self.var_dict['00p06']['sfc'] except for srols
        fc012up2tr3/ml91/glb100x100/2006/11/t_20061118_00p03.hdf --> one file per day, variables are self.var_dict['00p03']['ml91']
        fc012up2tr3/sfc/glb100x100/2010/10/v10m_20101003_00p03.hdf --> one file per day, variables are self.var_dict['00p03']['sfc']

        EI:
        an0tr6/sfc/glb100x100/2007/02/srols_200702.hdf --> one file per month
        an0tr6/sfc/glb100x100/2007/02/albedo_20070220_00p06.hdf --> one file per day
        fc012up2tr3/sfc/glb100x100/2001/12/sf_20011203_00p03.hdf --> one file per day
        fc012up2tr3/ml60/glb100x100/2007/08/q_20070802_00p03.hdf
        """

        file_names = []
        for tres in self.var_dict:
            if tres == '00p03':
                tkey = 'fc012up2tr3'
            elif tres == '00p06':
                tkey = 'an0tr6'
            for grid in self.var_dict[tres]:
                for var_name in self.var_dict[tres][grid]:
                    if var_name == 'srols' and date_obj.day != 1:
                        continue
                    if var_name == 'srols':
                        file_name = os.path.join(tkey, grid, 'glb100x100', date_obj.strftime("%Y/%m/srols_%Y%m.hdf"))
                    else:
                        file_name = os.path.join(tkey, grid, 'glb100x100', date_obj.strftime("%Y/%m"), "%s_%s_%s.hdf"%(var_name, date_obj.strftime("%Y%m%d"), tres))
                    file_names.append(file_name)

        return file_names

    def translateFileToPath(self, file_base_name):
        """
        Given the file name ec-ei-fc012up2tr3-ml60-glb100x100-t_20080101_00p03.hdf (only the base name),
        return the full path ec/ei/fc012up2tr3/ml60/glb100x100/2008/01/t_20080101_00p03.hdf. For files
        with names like ec-ei-an0tr6-sfc-glb100x100-oro.hdf, put them in ec/ei/an0tr6/sfc/glb100x100/oro.hdf.
        """
        file_parts = file_base_name.split('-')
        file_particles = file_parts[-1].split('_')
        if len(file_particles) > 1:
            # most files, such as t_20080101_00p03.hdf
            datetime_string = file_particles[1]
            file_parts.insert(-1, datetime_string[:4])
            file_parts.insert(-1, datetime_string[4:6])
        new_file_name = os.path.sep.join(file_parts)
        return new_file_name

    def translateCoarsenedFileToPath(self, file_base_name):
        """
        Given a file name like ec-ei-fc012up2tr3-tropo25-sea300x200-t_20110923_00p03.hdf, construct
        the name ec/ei/fc012up2tr3/tropo25/sea300x200/2011/09/t_20110923_00p03.hdf
        """
        file_parts = file_base_name.split('-')
        # file_parts is now ['ec', 'ei', 'fc012up2tr3', 'tropo25', 'sea300x200', 't_20110923_00p03.hdf']
        file_particles = file_parts[-1].split('_')
        # file_particles is now ['t', '20110923', '00p03.hdf']
        datetime_string = file_particles[1]
        file_parts.insert(-1, datetime_string[:4])
        # file_parts is now ['ec', 'ei', 'fc012up2tr3', 'tropo25', 'sea300x200', '2011', 't_20110923_00p03.hdf']
        file_parts.insert(-1, datetime_string[4:6])
        # file_parts is now ['ec', 'ei', 'fc012up2tr3', 'tropo25', 'sea300x200', '2011', '09', 't_20110923_00p03.hdf']
        new_file_name = os.path.sep.join(file_parts)
        # new_file_name is now ec/ei/fc012up2tr3/tropo25/sea300x200/2011/09/t_20110923_00p03.hdf
        return new_file_name

    def getTMPPfiles(self, pattern='*-glb100x100-*.hdf'):
        all_files = glob.glob(os.path.join(self.source_folder, pattern))
        return all_files

    def getCoarsenedFiles(self):
        source_folder = os.path.join(self.coarsened_archive, self.rcf.get('my.zoom'), self.rcf.get('ECLEVS'), self.rcf.get('LEVS'))
        all_files = glob.glob(os.path.join(source_folder, '*.hdf'))
        return all_files

    def copyFiles(self):
        print 'Getting list of source files'
        all_source_files = self.getTMPPfiles()
        print 'Finished getting list of source files'
        WidgetSet = ['Copying meteo : ', progressbar.Percentage(), ' ', progressbar.Bar(marker='#', left='|', right='|')]
        pbar = progressbar.ProgressBar(widgets=WidgetSet, maxval=len(all_source_files)).start()
        for i, source_file in enumerate(all_source_files):
            dirname, basename = os.path.split(source_file)
            new_file_name = self.translateFileToPath(basename)
            target_file = os.path.join(self.target_folder, new_file_name)
            dirname = os.path.dirname(target_file)
            if not os.path.isdir(dirname):
                os.makedirs(dirname)
            if not os.path.exists(target_file):
                shutil.copy(source_file, target_file)
            else:
                # Check if the source has a different size compared to the destination
                source_size = os.path.getsize(source_file)
                target_size = os.path.getsize(target_file)
                if source_size != target_size:
                    shutil.copy(source_file, target_file)
            pbar.update(i+1)
        pbar.finish()

    def untar_files(self, start_date, end_date, select_vars=None):
        """
        Untar files from /archive/sbasu into /projects/tm5meteo/tmm/hdf
        Sample files and the locations after untarring:

        /archive/sbasu/meteo/ec/od_L91/an0tr6/sfc/glb100x100/albedo_201010_00p06.tar    -->  /projects/tm5meteo/tmm/hdf/ec/od_L91/an0tr6/sfc/glb100x100/2010/10/albedo_201010??_00p06.hdf
        /archive/sbasu/meteo/ec/od_L91/fc012up2tr3/ml91/glb100x100/q_20101000_00p03.tar -->  /projects/tm5meteo/tmm/hdf/ec/od_L91/fc012up2tr3/ml91/glb100x100/2010/10/q_201010??_00p03.hdf
        /archive/sbasu/meteo/ec/od_L91/fc012up2tr3/sfc/glb100x100/t2m_201010_00p03.tar  -->  /projects/tm5meteo/tmm/hdf/ec/od_L91/fc012up2tr3/sfc/glb100x100/2010/10/t2m_201010??_00p03.hdf

        Within each tar file, the HDF file names are, for example, ec-od_L91-fc012up2tr3-ml91-glb100x100-q_20101001_00p03.hdf
        So we need to untar each file into a temporary folder, convert the base name (ec-od_L91-fc012up2tr3-ml91-glb100x100-q_20101001_00p03.hdf)
        to a full path name (ec/od_L91/fc012up2tr3/ml91/glb100x100/2010/10/q_20101001_00p03.hdf) and join that to self.target_folder (/projects/tm5meteo/tmm/hdf)
        to make up the destination file name. Then move the untarred files to their destinations. Finally, delete the temporary folder.
        """
        tar_files = self.list_glb100x100_tar_files(start_date, end_date, select_vars=select_vars)
        # first, check if all the files exist
        file_exist = True
        for file_name in tar_files:
            fe = os.path.isfile(file_name)
            file_exist = file_exist and fe
            if not fe:
                sys.stderr.write("File does not exist: %s\n"%file_name)
        if not file_exist:
            sys.stderr.write("Not all tar files exist, quitting\n")
            return
        # create a temporary folder
        temp_folder = tempfile.mkdtemp(prefix='tmp_meteo_dir_', dir=self.scratch_folder)
        WidgetSet = ['Extracting meteo : ', progressbar.Percentage(), ' ', progressbar.Bar(marker='#', left='|', right='|')]
        pbar = progressbar.ProgressBar(widgets=WidgetSet, maxval=len(tar_files)).start()
        for i,tar_file in enumerate(tar_files):
            arc_file = TarFile(tar_file, 'r')
            members = arc_file.getnames()
            # extract all the files
            arc_file.extractall(path=temp_folder)
            # now move the extracted files
            for base_name in members:
                source_file = os.path.join(temp_folder, base_name)
                target_file = os.path.join(self.target_folder, self.translateFileToPath(base_name))
                dirname = os.path.dirname(target_file)
                if not os.path.isdir(dirname):
                    os.makedirs(dirname)
                shutil.copy(source_file, target_file)
                os.remove(source_file)
            pbar.update(i+1)
        try:
            os.rmdir(temp_folder)
        except:
            sys.stderr.write("Error removing %s, please delete manually\n"%temp_folder)
        pbar.finish()

    def copyCoarsenedFiles(self):
        # like self.copyFiles, except that coarsened files are copied, not the original 1x1 files
        print 'Getting list of source files'
        all_source_files = self.getCoarsenedFiles()
        print 'Finished getting list of source files'
        WidgetSet = ['Copying meteo : ', progressbar.Percentage(), ' ', progressbar.Bar(marker='#', left='|', right='|')]
        pbar = progressbar.ProgressBar(widgets=WidgetSet, maxval=len(all_source_files)).start()
        for i, source_file in enumerate(all_source_files):
            # a source file will be called, for example, /archive/sbasu/meteo/tm5_meteo_out/se_asia/ml60/tropo25/ec-ei-fc012up2tr3-tropo25-sea300x200-t_20110923_00p03.hdf
            # it will have to be moved to /scratch/shared/sbasu/var4d/tm5_meteo_out/se_asia/ec/ei/fc012up2tr3/tropo25/sea300x200/2011/09/t_20110923_00p03.hdf
            # of which, self.coarsened_scratch contains /scratch/shared/sbasu/var4d/tm5_meteo_out/se_asia, and we need to make the
            # ec/ei/fc012up2tr3/tropo25/sea300x200/2011/09/t_20110923_00p03.hdf bit.
            dirname, basename = os.path.split(source_file)
            new_file_name = self.translateCoarsenedFileToPath(basename)
            target_file = os.path.join(self.coarsened_scratch, new_file_name)
            dirname = os.path.dirname(target_file)
            if not os.path.isdir(dirname):
                os.makedirs(dirname)
            if not os.path.exists(target_file):
                shutil.copy(source_file, target_file)
            else:
                # Check if the source has a different size compared to the destination
                source_size = os.path.getsize(source_file)
                target_size = os.path.getsize(target_file)
                if source_size != target_size:
                    shutil.copy(source_file, target_file)
            pbar.update(i+1)
        pbar.finish()

def hashfile(file_name, hasher, blocksize=134217728):
    """
    This is for taking the checksum of a file, usually used to check if two files are identical or not. The 'hasher'
    should be something like hashlib.md5() or hashlib.sha1().
    """
    afile = file(file_name, 'r')
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    afile.close()
    return hasher.digest()
