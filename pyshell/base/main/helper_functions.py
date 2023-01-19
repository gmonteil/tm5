from datetime import datetime, timedelta
import calendar, time, sys, os
import numpy as np
from netCDF4 import Dataset

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def SurfaceAreaGrid((lat_specs), (lon_specs)):
    """
    Returns the area of a rectangular lat/lon grid in square meters. 'lat_specs' is a tuple
    (lat_beg, lat_end, lat_divs) and 'lon_end' is a similar tuple for longitudes.
    """
    lat_beg, lat_end, lats = lat_specs
    lon_beg, lon_end, lons = lon_specs
    EarthRad = 6.371e6 # meters
    dLon = (np.pi/180.) * (lon_end-lon_beg)/lons
    dS = np.zeros((lats+1, lons), np.float64)
    Lats = (np.pi/180.) * np.linspace(lat_beg, lat_end, lats+1)
    for i, lat in enumerate(Lats):
        dS[i] = EarthRad * EarthRad * dLon * np.sin(lat)
    dS = np.diff(dS, axis=0)
    return dS

def PCDC_surface_area(nlat, nlon):
    # The surface area grid of a pole-centered, dateline-centered GEOS lat/lon grid
    EarthRad = 6.371e6 # meters
    # Construct the lat boundaries
    dlat = 180.0/(nlat-1) # in the bulk of the grid
    lats = np.zeros(nlat+1, dtype=np.float64)
    lats[1:nlat] = np.linspace(-90.0+dlat/2., 90.0-dlat/2., nlat-1)
    lats[0] = -90.0
    lats[nlat] = 90.0
    # convert to radians
    lats = (np.pi/180.) * lats
    # Longitudes are equi-spaced
    dlon = 2.0*np.pi/nlon
    # Now calculate the area
    # dS = R*R * (del_lon) * (sin(lat2)-sin(lat1)) # start indexing from South Pole
    dS = np.zeros((nlat+1, nlon), np.float64)
    for i, lat in enumerate(lats):
        dS[i] = EarthRad * EarthRad * dlon * np.sin(lat)
    dS = np.diff(dS, axis=0)
    return dS


def decimal_date_to_datetime(decimal_date):
    """
    Given a single decimal date such as 2003.4932, convert it to a datetime object datetime(2003, 6, 30, 0, 25, 55, 200000)
    """
    year = int(decimal_date)
    days = (365.0+calendar.isleap(year)) * (decimal_date - year)
    return datetime(year,1,1,0) + timedelta(days=days)

def decimal_date(date_array):
    """
    Given an array or list of datetime objects, convert them to decimal dates. For example, datetime(2009,3,15,9,12,30)
    would be converted to 2009.2010512.
    """
    # check if date_array is a single object
    if isinstance(date_array, datetime):
        return decimal_date_atomic(date_array)
    else:
        ret_array = np.array([decimal_date_atomic(d) for d in date_array], np.float64)
    return ret_array

def decimal_date_atomic(date_obj):
    # This is called by decimal_date to convert a single datetime object into a decimal date
    year = date_obj.year
    dt = date_obj - datetime(year,1,1)
    secs_in_year = dt.seconds + 86400.0 * dt.days
    total_secs = 86400.0 * (int(calendar.isleap(year)) + 365)
    return float(year) + secs_in_year/total_secs

class Timer(object):

    def __init__(self, timer_name=None):
        super(Timer, self).__init__()
        self.timer_name = timer_name

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start
        if self.timer_name is not None:
            print "Time taken for %s: %s"%(self.timer_name, self.repr_timedelta(self.interval))
            sys.stdout.flush()

    def repr_timedelta(self, s):
        hours, remainder = divmod(s, 3600)
        minutes, seconds = divmod(remainder, 60)
        cat_str = []
        if hours > 0:
            cat_str.append('%i hours'%hours)
        if minutes > 0:
            cat_str.append('%i minutes'%minutes)
        if seconds > 0:
            cat_str.append('%i seconds'%seconds)
        return ' '.join(cat_str)

class ExecEnvironment(object):
    # A class to make sure that foreground runs come back to the current folder
    def __init__(self, cur_dir, exec_dir):
        self.cur_dir = cur_dir
        self.exec_dir = exec_dir

    def __enter__(self):
        os.chdir(self.exec_dir)

    def __exit__(self, type, value, tb):
        os.chdir(self.cur_dir)

class del_time(object):

    def __init__(self, del_obj):
        super(del_time, self).__init__()
        # del_obj should be a timedelta object
        self.days = del_obj.days
        self.seconds = del_obj.seconds
        self.microseconds = del_obj.microseconds
        self.resolution = del_obj.resolution

    def to_seconds(self):
        # I'm sick of timedelta not implementing a to_seconds() method
        return (self.seconds + 1.0E-6 * self.microseconds + 86400.0 * self.days)

    def __div__(self, denom):
        self_sec = self.to_seconds()
        denom_sec = denom.to_seconds()
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

        # Also check if the folder in which the file resides exists
        dir_name = os.path.dirname(file_name)
        if not os.path.isdir(dir_name):
            raise RuntimeError('Folder %s does not exist or is not a folder'%dir_name)

        super(my_Dataset, self).__init__(file_name, *args, **kwargs)
        try:
            self.set_always_mask(False) # https://github.com/Unidata/netcdf4-python/issues/785#issuecomment-526711138 (in more recent versions)
        except AttributeError:
            pass

class Dummy_TQDM(object):

    def __init__(self, iterable, *args, **kwargs):
        super(Dummy_TQDM, self).__init__()
        self.iterator = iter(iterable)

    def __iter__(self):
        return self.iterator

    def __next__(self):
        return next(self.iterator)

    def update(self, *args, **kwargs):
        pass

class DummyProgress(object):

    def __init__(self, *args, **kwargs):
        super(DummyProgress, self).__init__()

    def __enter__(self, *args, **kwargs):
        pass

    def __exit__(self, *args, **kwargs):
        pass

    def next(self, *args, **kwargs):
        pass

    def finish(self, *args, **kwargs):
        pass
