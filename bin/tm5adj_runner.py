#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
======================================================================
DESCRIPTION



EXAMPLES


EXIT STATUS

AUTHOR

Michael Vossbeck <Michael.Vossbeck(at)Inversion-Lab.com>

======================================================================
"""

#- system packages
import os
from pathlib import Path
import sys
import os
import pathlib
import datetime as dtm
import calendar
import numpy as np
import numpy.ma as ma
import netCDF4 as nc4
import pandas as pd
import subprocess
from collections import OrderedDict
from loguru import logger
from argparse import ArgumentParser
from typing import List, Dict, Union

#
#-- iLab platform specific settings
#
_TM5DIR = '/srv/tm5' #-- TM5 top-level directory
_stafile_lst = [f"{_TM5DIR}/stationlist.verify2020.txt", f"{_TM5DIR}/co2m-stationlist_site_20231116_for-tm5.txt"]

#-- create parser
parser = ArgumentParser(description='build and run TM5 atmosphere tracer model')
parser.add_argument('--site',
                    default='so',
                    help="""observational site which is listed in station file (default: %(default)s).""")
parser.add_argument('--station_file',
                    choices=_stafile_lst,
                    default=_stafile_lst[1],
                    help="""csv file providing station coordinates (default: %(default)s).""")
parser.add_argument('--tpulse',
                    required=True,
                    help="""time-point of pulse UTC (any pandas Timestamp compliant format accepted, default: %(default)s).""")
parser.add_argument('--adj_trange',
                    nargs=2,
                    metavar=('yyyy-mm-dd_start', 'yyyy-mm-dd_end'),
                    help="""temporal range of adjoint run. These must be specified such that adjoint start is after *after* the adjoint end.""")
parser.add_argument('--mode',
                    choices=['point','satellite'],
                    default='point',
                    help="""which kind of departure file to generate (default: %(default)s).""")
parser.add_argument('--sat_interpolation',
                    choices=['gridbox','slopes','linear'],
                    default='linear',
                    help="""spatial interpolation applied in satellite mode (default: %(default)s).""")
parser.add_argument('--level_weights',
                    choices=['uniform','pressure'],
                    default='uniform',
                    help="""departure weights applied to the levels in satellite mode. Note, that for 'pressure' the TM5 forward simulation must have generated 'mix' output.""")
parser.add_argument('--rcfile',
                    help="""rc file for adjoint TM5 run (will be generated from forward.rc file, if not provided.""")
parser.add_argument('--forward_outdir',
                    required=True,
                    help="""top-level directory for outputs of forward run (required to happen before), adjoint outputs will also be placed in there (!!MANDATORY!!).""")
parser.add_argument('--adj_resdir',
                    help="""adjoint NetCDF files will be placed here.""")
parser.add_argument('--norun',
                    dest='run',
                    action='store_false',
                    help="""only prepare adjoin run, but do not start it (to inspect departure file etc...).""")


def move_file(source, destination):
    assert source.exists()
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open(mode="xb") as fp:
            fp.write(source.read_bytes())
    except FileExistsError:
        print(f"File {destination} exists already.")
    else:
        source.unlink()

def copy_file(source, destination):
    assert source.exists()
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open(mode="xb") as fp:
            fp.write(source.read_bytes())
    except FileExistsError:
        print(f"File {destination} exists already.")


def load_site( filepath : str, site_id : str ):
    """
    """
    header_c1 = ['ID', 'LAT', 'LON', 'ALT', 'TP', 'STATIONNAME']
    header_c2 = ['NUM','ID', 'LAT', 'LON', 'ALT', 'TP', 'STATIONNAME']
    site_dct = {}
    #-- the format is not well suitable for reading it with pandas...
    with open(filepath) as fp:
        for iline,line in enumerate(fp):
            tokens = line.split()
            if iline==0:
                assert tokens==header_c1 or tokens==header_c2, \
                    f"unexpected headers -->{tokens}<--"
            else:
                #-- running index at front of further lines
                inum = int(tokens[0])
                cur_id = tokens[1]
                if cur_id.lower() != site_id.lower():
                    continue
                else:
                    logger.info(f"selected site -->{site_id}<-- detected at inum={inum} (=={tokens}==)")
                    #-- getting coordinates and altitude plus station name
                    site_dct['lat'] = float(tokens[2])
                    site_dct['lon'] = float(tokens[3])
                    site_dct['alt'] = float(tokens[4])
                    site_dct['name'] = tokens[6].lower()
    #
    if len(site_dct)==0:
        raise RuntimeError(f"selected site {site_id} not found in station file -->{filepath}<--")
    return site_dct
# ---end-of-load_site


#
#-- Email Guillaume, 2024-05-17
#
# > The TM5 levels are probably not evenly distributed in terms of air mass. We would then need to adjust the factors 1/25 for this. You probably have used such adjusted weighting when you computed the XCO2/XCH4 in forward runs, so it is probably somewhere in the model ...
# >
# It's will depend on the location and time. The relative weight of each layer should be the difference between it's lower and upper pressure, divided by the surface pressure. And the pressure at the interface is obtained using the "a" and "b" coefficients of the sigma-pressure coordinates sytem (they are defined in the "dims_level.F90" file, which just picks the right levels from the "const_ec_v.F90" file, but you should probably find them also in some of the gridded output files). The equation is P(l) = a + b * Psurf
def load_pressurelevel_weights(mixfile : Path, lat : float, lon : float,
                               timepoint : pd.Timestamp, region_id : str = 'glb100x100'):
    """
    """
    fp = nc4.Dataset(str(mixfile))
    #
    #-- get sigma-pressure coordinates
    #
    at = fp.variables['at'][:]
    bt = fp.variables['bt'][:]

    #
    #-- spatial coordinates
    #
    lon_edges = fp[f'/{region_id}/lon_edges'][...]
    lat_edges = fp[f'/{region_id}/lat_edges'][...]
    cndlon = (lon>=lon_edges[:-1])&(lon<=lon_edges[1:])
    cndlat = (lat>=lat_edges[:-1])&(lat<=lat_edges[1:])
    ilon = np.where(cndlon)[0][0]
    ilat = np.where(cndlat)[0][0]
    grid_lon = 0.5*(lon_edges[ilon]+lon_edges[ilon+1])
    grid_lat = 0.5*(lat_edges[ilat]+lat_edges[ilat+1])
    msg = f"lon/lat = {lon}/{lat} detected at grid indices ilon/ilat = {ilon}/{ilat}"
    msg += f" (grid-cell centre {grid_lon}/{grid_lat})"
    logger.info(msg)
    #
    #-- temporal dimension
    #
    sample_times = fp[f'/{region_id}/sample_times'][...]
    timedelta_list = []
    time_list      = []
    for it,t in enumerate(sample_times):
        yr,mon,dy,hr,mn,sec = t
        tstr = f"{yr:04d}-{mon:02d}-{dy:02d} {hr:02d}:{mn:02d}:{sec:02d}"
        cur_time = pd.Timestamp(tstr)
        time_list.append(cur_time)
        if cur_time>=timepoint:
            timedelta_list.append( cur_time - timepoint )
        else:
            timedelta_list.append( timepoint - cur_time )
    timedelta_list = np.array(timedelta_list)
    it = np.argmin(timedelta_list)
    msg = f"time closest to pulse detected at it={it} -->{time_list[it]}<-- (diff: {timedelta_list[it]})"
    logger.info(msg)
    #
    #-- load surface pressure
    #
    ncvar = fp[f'/{region_id}/pressure']
    assert ncvar.dimensions==('times','latitude','longitude')
    assert ncvar.long_name=="surface pressure (in Pa)"
    psurf = ncvar[it,ilat,ilon]
    msg = f"detected surface pressure at coordinates and time {psurf} [Pa]"
    logger.info(msg)

    #
    #-- pressure levels: P(l) = a + b * Psurf
    pressure_levels = at + bt*psurf

    departures = (pressure_levels[:-1] - pressure_levels[1:])/psurf
    
    fp.close()

    return departures
# ---end-of-load_pressurelevel_weights


def create_point_departure(site_tag : str,
                           lat : float, lon : float, alt : float,
                           timepoint : pd.Timestamp,
                           outdir : Union[str, Path]):
    #
    #-- currently restricted to 1x1 degree global grid and CO2 tracer
    ntracer = 1
    tracer  = 'CO2'
    region  = 'glb100x100' #-- 1 by 1 degree global

    sampling_tag = 'instantaneous-sampling' #-- see below
    # site_tag = site_dict['id']
    time_tag = timepoint.strftime('%Y-%m-%dT%H%M%S')

    basename = f"point_departures_{sampling_tag}_{site_tag}_{time_tag}.nc4"
    outname = outdir / basename
    outname.parent.mkdir(parents=True, exist_ok=True)


    ##
    logger.info(f"start generation of point departure file...")
    fp = nc4.Dataset(outname, 'w')

    gid = fp.createGroup(region)

    ggid = gid.createGroup(tracer)

    #-- single sample
    nflasks = 1
    ggid.createDimension('samples', nflasks)
    ggid.createDimension('idate', 6)#-- year/mon/day/hour/minute/second

    #-- create variables (see adj_user_output_flask.F90, lines 202ff)

    #-- 6-component date and time of sampling
    ncvar = ggid.createVariable('date_components', 'i4', ('samples','idate'))
    t = timepoint
    ncvar[0,:] = [t.year, t.month, t.day, t.hour, t.minute, t.second]

    #-- number of samples accumulated
    ncvar = ggid.createVariable('nsamples', 'i4', ('samples',))
    ncvar[:] = 1

    #-- total weight accumulated during the forward run
    ncvar = ggid.createVariable('total_weight', 'f8', ('samples',))
    ncvar[:] = 1

    #-- sampling strategy (1 for 4-hour average, 2 for instantaneous samples, 3 for dT sampling, 4 for custom)
    ncvar = ggid.createVariable('sampling_strategy', 'i4', ('samples',))
    ncvar[:] = 2

    #-- custom time window length for averaging
    ncvar = ggid.createVariable('time_window_length', 'i4', ('samples',))
    ncvar[:] = 3600 #-- looked up in output/point_input.nc4

    #-- mixing ratio forcing of interest
    ncvar = ggid.createVariable('forcing', 'f8', ('samples',))
    ncvar[:] = 1.

    ncvar = ggid.createVariable('lat', 'f8', ('samples',))
    ncvar[:] = lat

    #-- longitude of sampling
    ncvar = ggid.createVariable('lon', 'f8', ('samples',))
    ncvar[:] = lon

    #-- altitude of sampling
    ncvar = ggid.createVariable('alt', 'f8', ('samples',))
    ncvar[:] = alt

    fp.close()

    logger.info(f"...generated departure file {outname}")

    #-- filemame in TM5
    default_outname = outdir / 'point_departures.nc4'
    try:
        default_outname.unlink()
    except FileNotFoundError as e:
        logger.info(f"no need to unlink -->{default_outname}<--")
        pass
    default_outname.symlink_to(basename)
    logger.info(f"...set system link for file -->{default_outname}<--")
# ---end-of-create_point_departure


def create_satellite_departure(site_tag : str,
                               lat : float, lon : float, alt : float,
                               timepoint : pd.Timestamp,
                               departures : np.ndarray,
                               departure_tag : str,
                               outdir : Union[str, Path]):

    """function to create satellite track departure file for adjoint run
    with TM5.
    TM5 will read from
    outdir_satellite/satellite/sat-track_departures_yyyymm.nc4   (if split_period=='m')
    outdir_satellite/satellite/sat-track_departures_yyyymmdd.nc4 (if split_period=='d')

    This routine is currently limited to
    - 1x1 degree global grid ('glb100x100')
    - tracer name 'CO2'
    - 25 levels ('tropo25')
    """

    #
    #-- tracer/region settings
    #
    ntracer = 1
    tracer  = 'CO2'
    region  = 'glb100x100' #-- 1 by 1 degree global
    nlevel = len(departures)
    assert nlevel==25, \
        f"currently expecting 25 vertical levels only ('tropo25')"

    time_tag = timepoint.strftime('%Y-%m-%dT%H%M%S')

    #-- sampling strategy
    #   2: instantaneous
    #   3: wihin ndyn/tref
    #   -> see more in user_output_satellite.F90, line 39ff 
    sampling_strategy = ('instantaneous-sampling', 2)
    
    basename = f"sat-track_departures_{departure_tag}_{sampling_strategy[0]}_{site_tag}_{time_tag}.nc4"
    outname = outdir / basename
    outname.parent.mkdir(parents=True, exist_ok=True)


    ##
    logger.info(f"start generation of satellite-track departure file...")
    fp = nc4.Dataset(outname, 'w')

    gid = fp.createGroup(region)

    ggid = gid.createGroup(tracer)

    #
    #-- single pulse
    #
    n_obs = 1
    ggid.createDimension('n_obs', n_obs)
    ggid.createDimension('idate', 6) #-- year/mon/day/hour/minute/second
    ggid.createDimension('nlevel', nlevel)

    #
    #-- create variables (see adj_user_output_satellite.F90, lines 248ff)
    #
    #
    #-- Note: look into 'type sat_tracer_forcing' (user_output_satellite_data.F90, lines 62ff)
    #         nsamples, sampling_strategy, idates are all integer(kind=2)

    
    #
    #-- 6-component date and time of sampling
    #   MVO NOTE:naming 'idate' differs from adj_user_output_flask.F90 ('date_components')
    #
    ncvar = ggid.createVariable('idate', 'i2', ('n_obs','idate'))
    t = timepoint
    ncvar[0,:] = [t.year, t.month, t.day, t.hour, t.minute, t.second]

    #-- number of samples accumulated, MVO::unclear if really needed, but it is read by TM5
    ncvar = ggid.createVariable('nsamples', 'i2', ('n_obs',))
    ncvar[:] = 1

    #-- total weight accumulated during the forward run
    ncvar = ggid.createVariable('total_weight', 'f8', ('n_obs',))
    ncvar[:] = 1

    #-- sampling strategy
    #   2: instantaneous
    #   3: wihin ndyn/tref
    #   -> see more in user_output_satellite.F90, line 39ff
    #  NOTE: 
    ncvar = ggid.createVariable('sampling_strategy', 'i2', ('n_obs',))
    ncvar[:] = sampling_strategy[1]

    #
    #-- 
    #
    ncvar = ggid.createVariable('departures', 'f8', ('n_obs','nlevel'))
    ncvar[:] = departures
   
    ncvar = ggid.createVariable('lat', 'f8', ('n_obs',))
    ncvar[:] = lat

    #-- longitude of sampling
    ncvar = ggid.createVariable('lon', 'f8', ('n_obs',))
    ncvar[:] = lon

    #-- altitude of sampling !-- *NOT* read by adjoint code
    ncvar = ggid.createVariable('alt', 'f8', ('n_obs',))
    ncvar[:] = alt

    fp.close()

    logger.info(f"...generated departure file {outname}")

    #
    #-- currently we fix to split_period='m'
    #-- filemame in TM5
    #
    default_outname = outdir / f"sat-track_departures_{timepoint.strftime('%Y%m')}.nc4"
    try:
        default_outname.unlink()
    except FileNotFoundError as e:
        logger.info(f"no need to unlink -->{default_outname}<--")
        pass
    default_outname.symlink_to(basename)
    logger.info(f"...set system link for file -->{default_outname}<--")
# ---end-of-create_satellite_departure


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#                    M A I N
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-- parse command line
args = parser.parse_args(sys.argv[1:])


msg = f"running command ***{' '.join(sys.argv)}***"
logger.info(msg)

forward_outdir = args.forward_outdir

#
#-- site information
#
site_info = load_site(args.station_file, args.site)
site_tag = f"{args.site}-{site_info['name']}"
lat, lon, alt = site_info['lat'], site_info['lon'], site_info['alt']

#
#-- time point of backward pulse LST
#   CHANGED:2024-06-03 time has to be passed in UTC now,
#           conversion code (below) for conversion LST to UTC no longer required
#
tpulse_utc     = pd.Timestamp(args.tpulse)
# #
# #-- time point of backward pulse UTC,
# #   rounded to full hour
# #
# tpulse_utc = tpulse - pd.Timedelta(hours=lon/15.) #-- time of pulse [UTC]
# msg = f"time-of-pulse {tpulse}[LST] yields {tpulse_utc} [UTC]"
# logger.info(msg)
# tpulse_utc = tpulse_utc.round(freq='H')
# msg = f"time-of-pulse, rounded to full hour {tpulse_utc} [UTC]"
# logger.info(msg)

if args.rcfile==None:
    #
    #-- no adjoint rcfile is provided,
    #   try using the rcfile from forward run
    #   and adjust/add keys
    #
    rcfile = Path(f'{forward_outdir}/forward.rc')
    assert rcfile.exists()

    rc_dct = OrderedDict()
    with open(rcfile, 'r') as fp:
        for line in fp:
            tokens = line.replace('\n','').split(':', 1) #-- split at first ':' only
            tokens = [t.strip() for t in tokens]
            k,v = tokens
            rc_dct[k] = v
    #
    #-- temporal range of adjoint run
    #
    fwd_tstart = pd.Timestamp(rc_dct['jobstep.timerange.start'])
    fwd_tend   = pd.Timestamp(rc_dct['jobstep.timerange.end'])
    if args.adj_trange!=None:
        adj_tstart = pd.Timestamp(args.adj_trange[0])
        adj_tend   = pd.Timestamp(args.adj_trange[1])
        assert adj_tstart>adj_tend
        assert adj_tstart<=fwd_tend
        assert adj_tend>=fwd_tstart
        rc_dct['jobstep.timerange.start'] = adj_tend.strftime('%Y-%m-%d %H:%M:%S')
        rc_dct['jobstep.timerange.end']   = adj_tstart.strftime('%Y-%m-%d %H:%M:%S')
    else:
        adj_tstart = fwd_tend
        adj_tend   = fwd_tstart
        #-- now using 1 hour
        if tpulse_utc - pd.Timedelta(days=1) >= fwd_tstart:
            adj_tend = tpulse_utc - pd.Timedelta(days=1)
            rc_dct['jobstep.timerange.start'] = adj_tend.strftime('%Y-%m-%d %H:%M:%S')
    #
    #-- consistencty for time of pulse
    #
    assert adj_tend<=tpulse_utc<=adj_tstart, \
        f"selected timepoint of backward pulse {args.timepoint} not in temporal range of rcfile {adj_rcfile}"

    assert adj_tstart==pd.Timestamp(adj_tstart.strftime('%Y-%m-%d')), \
        f"adjoint TM5 run must start at H:M:S 00:00:00 but have adj_tend={adj_tend}"
    #
    #-- swith to run mode 2 (for adjoint)
    #
    rc_dct['my.runmode'] = "2"

    #MVO::maybe not required for adjoint run
    #-- output.satellite.output.directory

    #
    #--  output.satellite.interpolation
    # user_output_satellite_data.F90:
    # integer, parameter    :: SAT_INTERPOLATION_GRIDBOX = 1
    # integer, parameter    :: SAT_INTERPOLATION_SLOPES = 2
    # integer, parameter    :: SAT_INTERPOLATION_LINEAR = 3
    satinterpolation_dct = {'gridbox':1, 'slopes':2, 'linear':3}

    sat_ipmode = str(satinterpolation_dct[args.sat_interpolation])
    if args.mode=='point':
        adj_tag = args.mode
    else:
        adj_tag = f"{args.mode}-{args.sat_interpolation}"
    
    adj_rcfile = f"{forward_outdir}/adjoint_{site_tag}_{adj_tag}_{tpulse.strftime('%Y%m%dT%H')}-LST.rc"
    if args.mode=='point':
        rc_dct['adjoint.input.point']     = 'T'
        rc_dct['adjoint.input.satellite'] = 'F'
        if not 'output.point.verbose' in rc_dct:
            rc_dct['output.point.verbose'] = 'T'
    elif args.mode=='satellite':
        rc_dct['output.point']                   = 'F'
        rc_dct['adjoint.input.point']            = 'F'
        rc_dct['adjoint.input.satellite']        = 'T'      
        rc_dct['output.satellite.verbose']       = 'T'
        rc_dct['output.satellite.split.period']  = 'm'  #-- sat-track_yyyymm.nc4
        rc_dct['output.satellite.interpolation'] = sat_ipmode
    with open(adj_rcfile, 'w') as fpout:
        for k,v in rc_dct.items():
            line = f"{k} : {v}"
            fpout.write(line + '\n')
        logger.info(f"generated rcfile for adjoint run ***{adj_rcfile}***")
else:
    #
    #-- an rcfile for adjoint is provided,
    #   only copy to default name
    #
    copy_file(Path(args.rcfile), Path(adj_rcfile))




#
#-- create departure file
#
if args.mode=='point':
    ptoutdir = Path(f'{args.forward_outdir}/point')
    create_point_departure(site_tag, lat, lon, alt, tpulse_utc, ptoutdir)
elif args.mode=='satellite':
    if args.level_weights=='uniform':
        #
        #-- initial approach: uniform over levels
        #
        #
        #-- fixed to 25 level simulations ('tropo25')
        #
        nlevel = 25
        departures = np.full(nlevel, 1/nlevel)
        departure_tag = 'uniform-level-weights'
        logger.info(f"created uniform departures over levels")
    else:
        #
        #-- naming of example output for mixing ratio as generated by forward run:
        #   mix/2018/06/mix20180601.nc4
        #
        mixoutdir = Path(f"{args.forward_outdir}/mix/{tpulse_utc.year}/{tpulse_utc.month:02d}")
        mixfile = mixoutdir / f"mix{tpulse_utc.strftime('%Y%m%d')}.nc4"
        if mixfile.exists():
            msg = f"derive level-dependent departures from file ***{mixfile}***"
            logger.info(msg)
            departures = load_pressurelevel_weights( mixfile, lat, lon, tpulse_utc)
            departure_tag = 'pressure-level-weights'
        else:
            msg = f"required input for computing pressure level weights is missing ==>{mixfile}<=="
            raise RuntimeError(msg)
    satoutdir = Path(f"{args.forward_outdir}/satellite")
    create_satellite_departure(site_tag, lat, lon, alt, tpulse_utc, departures, departure_tag, satoutdir)
#
#-- TM5 command line
#
exe = f'./{args.forward_outdir}/tm5.x'
cmd = [exe, str(adj_rcfile)]
logger.info(f"starting running  TM5 with command {'_'.join(cmd)}...")
if args.run:
    try:
        _ = subprocess.run(cmd, stdout=None)
        if _.returncode != 0 :
            raise RuntimeError(_)
    except Exception as e:
        print(f"exception -->{e}<--")
        raise

#
#-- expected NetCDF files generated by TM5 adjoint run
#
adj_emis_file  = Path(args.forward_outdir) / "adj_emissions.nc4"
adj_sav_file   = Path(args.forward_outdir) / f"adj_save_{adj_tend.strftime('%Y%m%d%H%M')}.nc4"
adj_tprof_file = Path(args.forward_outdir) / 'timing' / f"{rcfile.stem}.prf"

if args.run:
    assert adj_emis_file.exists(), \
        f"file ***{adj_emis_file}*** was not generated"
    assert adj_sav_file.exists(), \
        f"file ***{adj_sav_file}*** was not generated"
else:
    adj_emis_file.touch(exist_ok=True)
    adj_sav_file.touch(exist_ok=True)

#-- move files to target directory
if args.adj_resdir!=None:
    adj_resdir = Path(args.adj_resdir)
else:
    adj_resdir = Path(args.forward_outdir) / f"adjrun_{args.site}_{tpulse_utc.strftime('%Y%m%dT%H%M%S')}"
    adj_resdir.mkdir(parents=True, exist_ok=True)
for src in [adj_emis_file, adj_sav_file, adj_tprof_file]:
    if src.exists():
        dst = adj_resdir / src.name
        print(f"moving file -->{src}<-- to -->{dst}<--")
        move_file(src, dst)
