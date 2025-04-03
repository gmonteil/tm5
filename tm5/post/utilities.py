#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Miscallaneous utility functions."""

#-- system packages
import os
import sys
from pathlib import Path
import datetime as dtm
import time
import calendar
import argparse
import pandas as pd
import numpy as np
from typing import Union


# WGS84 ellipsoid:
# https://en.wikipedia.org/wiki/Earth_radius
# Equatorial radius: a = (6378.1370 km)
# Polar radius:      b = (6356.7523 km)
_A = 6378.1370 * 1000
_B = 6356.7523 * 1000
# In geophysics, the International Union of Geodesy and Geophysics (IUGG) defines the Earth's mean radius (denoted R1) to be (2a + b)/3
#
_ERmeter  = (2*_A + _B)/3


def lonstr(lon, fmt="{:.3f}"):
    if lon<0:
        return fmt.format(abs(lon))+'W'
    else:
        return fmt.format(lon)+'E'
# ---end-of-lonstr


def latstr(lat, fmt="{:.3f}"):
    if lat<0:
        return fmt.format(abs(lat))+'S'
    else:
        return fmt.format(lat)+'N'
# ---end-of-latstr


def to_posix_ts(d: dtm.datetime, utc:bool=True) -> float:
    tt=d.timetuple()
    return (calendar.timegm(tt) if utc else time.mktime(tt)) + round(d.microsecond/1000000, 0)


def pd_timestamp_from_datetime(d: dtm.datetime) -> pd.Timestamp:
    return pd.to_datetime(to_posix_ts(d), unit='s')


def lst_intersect(l1,l2):
    """Returns ordered list of common entries in both lists/iterables.

    Parameters
    ----------
    l1 : list/tuple
        first list or tuple of elements
    l2 : list/tuple
        second list or tuple of elements

    Returns
    -------
    list
    """
    s1 = set(l1)
    s2 = set(l2)
    s = s1.intersection(s2)
    return sorted(list(s))
# ---end-of-lst_intersect


def X_is_running():
    """
    Function that checks whether the DISPLAY/X-server is running/accessible
    (will likely work only on Linux-like systems)
    """
    import os
    try:
        display = os.environ['DISPLAY']
        return True
    except KeyError:
        return False

    # from subprocess import Popen, PIPE
    # p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
    # p.communicate()
    # return p.returncode == 0
# ---end-of-X_is_running


def create_sha512(fname, blocksize=None):
    """Function to create sha512 checksum from file-content of filename given,
    in case of I/0 error(s) a value of -1 will be returned.

    Parameters
    ----------
    fname : str
       file name
    blocksize : int, optional
       reading the file in portions of kilobytes (by default file is read completly)

    Returns
    -------
    str
       SHA512 hash of the file if successful, None otherwise
    
    """
    import hashlib
    #-- sha512 instance
    hasher = hashlib.sha512()

    with open(fname,'rb') as afile:
        if blocksize==None:
            buf = afile.read()
            hasher.update(buf)
        else:
            while True:
                buf = afile.read(blocksize)
                if not buf:
                    break
                else:
                    hasher.update(buf)
        return hasher.hexdigest()

    #-- erroneous exit
    return None


def mkdirp_smart(newdir):
    """
    Function that creates new directory and tries to catch most of possible collisions,
    should act somehow similar to linux command 'mkdir -p',
    raises RuntimeError in case of failure.

    Parameters
    ----------
    newdir : str
        directory path which will be created.

    Raises
    ------
    RuntimeError
        In case the directory could not be created.
    """
    import os

    if os.path.isdir(newdir):
        pass
    elif os.path.exists(newdir):
        msg = "a non-directory path with the same name as the desired " \
              f"dir, '{newdir}', already exists."
        raise RuntimeError(msg)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            mkdirp_smart(head)
        if tail:
            try:
                os.mkdir(newdir)
            except OSError as exc:
                msg = f"newdir={newdir} could not be created on system (exc={exc})"
                raise RuntimeError(msg)


def area_on_earth( dlon, lat1, lat2,
                   angle_unit='deg', area_unit='km^2', ERmeter=_ERmeter ):
    """
    Calculate area of rectangular bounded region on a sphere.
    Default radius is the (approximated) radius of the Earth.

    Parameters
    ----------
    dlon : west-to-east extent
    lat1 : latitude of northern bound
    lat2 : latitude of southern bound
    angle_unit : unit of angular measure ('deg' or 'rad')
    ERmeter    : radius of Earth [m]
    area_unit  : unit of computed area ([km^2],[ha],[m^2])

    Returns
    -------
    area of rectangle on sphere (in [m^2], [ha] or [km^2])

    """
    #-- convert angular unit to 'rad':
    if angle_unit=='deg':
        dlon = np.radians(dlon)
        lat1 = np.radians(lat1)
        lat2 = np.radians(lat2)
    elif angle_unit=='rad':
        pass
    else:
        return RuntimeError(f"unexpected angle_unit={angle_unit}")

    #--
    dlon = np.abs(dlon)
    aweight = dlon*np.abs(np.sin(lat2)-np.sin(lat1))

    if area_unit in ['m^2','m2']:
        area = ERmeter*ERmeter*aweight
    elif area_unit=='ha':
        area = (ERmeter/100.)*(ERmeter/100.)*aweight
    elif area_unit=='km^2':
        area = (ERmeter/1000.)*(ERmeter/1000.)*aweight

    return area


def haversine( lon1, lat1, lon2, lat2,
               angle_unit='deg', rm=_ERmeter, dist_unit='m' ):
    """
    Calculate the great circle distance between two points 
    on a sphere according to the haversine formula.
    Default radius is the (approximated) of the Earth in meter.

    Parameters
    ----------
    lon1 : longitude of first point (on Earth)
    lat1 : latitude of first point
    lon2 : longitude of second point
    lat2 : latitude of second point
    angle_unit : unit of angular measure ('deg' or 'rad')
    rm   : radius of shpere [m]
    dist_unit  : unit used for the distance computation [m]/[km]

    Returns
    -------
    distance of two points on Earth (in [m] or [km])

    """
    # convert decimal degrees to radians
    if angle_unit=='deg':
        lon1 = np.deg2rad(lon1)
        lat1 = np.deg2rad(lat1)
        lon2 = np.deg2rad(lon2)
        lat2 = np.deg2rad(lat2)
    elif angle_unit=='rad':
        pass
    else:
        msg = f"unexpected angle_unit={angle_unit}"
        raise ValueError(msg)

    dlon = np.abs(lon2 - lon1)
    dlat = np.abs(lat2 - lat1)

    #-- haversine formula 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 

    #-- distance [m]
    dm = rm * c

    if dist_unit=='m':
        return dm
    elif dist_unit=='km':
        return dm/1000.
    else:
        msg = f"unexpected dist_unit +++{dist_unit}+++"
        raise RuntimeError(msg)


def sphere_grid_find_close(lonq, latq, longrd, latgrd, distm_max=None, radius_m=_ERmeter):
    """
    Determine indices of grid points with maximal distance 'distm_max'
    to query point based on spherical geometry (.

    Parameters
    ----------
    lonq : float
        longitude of query point (degrees)
    latq : float
        latitude of query point (degrees)
    longrd : numpy array of float (1 or 2 dimensional)
        longitude of grid points (degrees)
    latgrd : numpy array of float (1 or 2 dimensional)
        latitude of grid points (degrees)
    distm_max : float, optional
        maximal distance from query point (in meter)
        if none is given the grid-points closest to the pixel will be returned.
    radius_m : float, optional
        radius of sphere in meter (by default average radius of the Earth is taken)

    Returns:
    --------
    tuple
        first element are the indices of the close grid-points (result from np.where),
        second element are the distances of these grid-points to the query point.
    """
    #-- some consistency checks
    assert longrd.shape==latgrd.shape
    assert 1<=len(longrd.shape)<=2 #-- 1d/2d grid arrays expected

    #-- convert to radians
    longrd_rad = np.radians(longrd)
    latgrd_rad = np.radians(latgrd)
    lonq_rad = np.radians(lonq)
    latq_rad = np.radians(latq)
    dlat = np.abs( latgrd_rad - latq_rad )
    dlon = np.abs( longrd_rad - lonq_rad )
    a = np.sin(dlat/2)**2 + np.cos(latq_rad) * np.cos(latgrd_rad) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    dm = radius_m * c

    if distm_max!=None:
        idxs_min = np.where(dm<=distm_max)
    else:
        idxs_min = np.where(dm==dm.min())
    dist_min = dm[idxs_min]

    return (idxs_min, dist_min)


def set_outname(optionsORdir : Union[str, argparse.Namespace], aname : str, only_dir=False):
    """Function to assemble name of an output file according to a suggested name
    and settings made on command line.
    """
    from pathlib import Path

    outname = aname

    if type(optionsORdir)==argparse.Namespace:
        outdir = optionsORdir.outdir
        if only_dir:
            pass
        elif optionsORdir.outname!=None:
            outname = optionsORdir.outname
    else:
        outdir = optionsORdir

    #-- 
    if outdir!=None:
        if not os.path.isabs(outname):
            outname = os.path.join(outdir,outname)
        else:
            msg = f"outdir ***{outdir}*** ignored since outname is an absolute path " \
                f"+++{outname}+++"
            PkgLogger.warn(msg)
    #-- ensure that intermediate directories exist
    odir = Path(os.path.dirname(outname))
    odir.mkdir(exist_ok=True, parents=True)

    return outname
