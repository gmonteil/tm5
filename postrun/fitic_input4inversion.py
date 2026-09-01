#!/usr/bin/env python

from argparse import ArgumentParser, Namespace as ArgumentNamespace
import sys
import os
from omegaconf import OmegaConf, DictConfig
from pathlib import Path
from collections import OrderedDict
import datetime as dtm
from loguru import logger
import pandas as pd
from pandas import date_range, DatetimeIndex,DataFrame
from pandas import Timestamp, Timedelta, concat
import xarray as xr
import numpy as np
from numpy import zeros, tile
from netCDF4 import Dataset, stringtochar
import xesmf
from types import SimpleNamespace
import pickle
import lzma

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots,colorbar
import matplotlib.patches as mpatches
from cartopy import crs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
from tm5.fitic import read_obs_table
from tm5.gridtools import TM5Grids
from tm5.observations import read_obspack_file
from tm5.post.footprint_io import load_adjoint_fwd #-- this was for earlier diagnostics
from tm5.fitic import get_fitic_region_table
# from tm5.post.footprint_io import tm5rundir_obstable, tm5rundir_iniconc_1obs
# from tm5.post.footprint_io import regions1D_info
# from tm5.post.footprint_io import tm5rundir_jacobian3D, tm5rundir_simustart
# from tm5.post.footprint_io import tm5rundir_jacobian3D_old, tm5_fitic_adjoint_corrected_halos
# from tm5.post.footprint_io import tm5rundir_emissions2D, tm5emisdir_load_emissions2D
# from tm5.post.footprint_io import jacobian_redistribute_glb6x4_to_avengers_zoom
from tm5.post.plot_util import cnorm_set
from tm5.post.utilities import lonstr,latstr,set_outname,create_sha512


# #-- glb6x4
# dom_glb6x4 = SimpleNamespace(name='glb6x4', west=-180, east=180, south=-90, north=90, dlon=6, dlat=4)
# #-- eur3x2
# dom_eur3x2 = SimpleNamespace(name='eur3x2', west=-36, east=54, south=22, north=74, dlon=3, dlat=2)



# #
# #-- drop 30W-48E x 26N-70N (glb6x4 indices 29-39, 25-38, end included)
# #
# glb6x4_dropmask = np.full((45,60), False)
# glb6x4_dropmask[29:40,25:38] = True #-- 143 grid-cells
# glb6x4_ndrop = np.count_nonzero(glb6x4_dropmask)
# glb6x4_ng = (glb6x4_dropmask.size-glb6x4_ndrop)
# msg = f"glb6x4: ndrop={glb6x4_ndrop} (of {glb6x4_dropmask.size}) ng={glb6x4_ng}"
# print(msg)
# #
# #-- drop 3E-15E x 44N-56N  (eur3x2 indices 11-16,13-17)
# #
# eur3x2_dropmask = np.full((26,30), False)
# eur3x2_dropmask[11:17,13:17] = True #-- 24 grid-cells
# eur3x2_halomask = np.full((26,30), True)
# eur3x2_halomask[2:26-2,2:30-2] = False
# eur3x2_dropmask |= eur3x2_halomask
# eur3x2_ndrop = np.count_nonzero(eur3x2_dropmask)
# eur3x2_ng = (eur3x2_dropmask.size-eur3x2_ndrop)
# msg = f"eur3x2: ndrop={eur3x2_ndrop} (of {eur3x2_dropmask.size}) ng={eur3x2_ng}"
# print(msg)
# #
# #-- drop HALO part of gns1x1 domain
# #
# gns1x1_halomask = np.full((16,18), True)
# gns1x1_halomask[2:16-2,3:18-3] = False
# gns1x1_dropmask = gns1x1_halomask
# gns1x1_ndrop = np.count_nonzero(gns1x1_dropmask)
# gns1x1_ng = (gns1x1_dropmask.size-gns1x1_ndrop)
# msg = f"gns1x1: ndrop={gns1x1_ndrop} (of {gns1x1_dropmask.size}) ng={gns1x1_ng}"
# print(msg)
# ng = glb6x4_ng + eur3x2_ng + gns1x1_ng
# msg = f"number of grid-cells contributing to emissions and Jacobian, ng={ng}"
# logger.info(msg)

# drop_table = {
#     'glb600x400': glb6x4_dropmask,
#     'eur300x200': eur3x2_dropmask,
#     'gns100x100': gns1x1_dropmask
#     }


def subcmd_prepare_daily_obsjacobian(args : ArgumentNamespace) -> None: 
    """
    """
    complevel = args.__dict__.get('complevel',4)
    #
    #-- load obstable
    #2021-01-01 01:00:00 in second column
    msg = f"reading obstabe from file ***{str(args.obsfile_filepath)}***"
    logger.info(msg)
    obstable = pd.read_csv(args.obsfile_filepath, index_col=0, parse_dates=[1], date_format='%Y-%m-%d %H:%M:%S')
    region_list = obstable['region'].unique()
    if 'gns100x100' in region_list:
        regions = ['glb600x400','eur300x200','gns100x100',]
        domain_tag = 'gns1x1'
    else:
        regions = ['glb600x400',]
        domain_tag = 'glb6x4'
    msg = f"...yields nobs={len(obstable)} (regions ==>{regions}<==)"
    logger.info(msg)
    days = sorted(obstable.time.dt.date.drop_duplicates())
    obs_firstday = Timestamp(days[0])
    obs_lastday  = Timestamp(days[-1])
    msg = f"last observational day in obstable -->{obs_lastday}<--"
    logger.info(msg)
    #
    #-- drop entries where actual observation is missing
    #
    msg = f"dropping entries where actual observation is missing..."
    logger.info(msg)
    obs_missval = -99999.
    cnd_keep = obstable.loc[:,'mixing_ratio']!=obs_missval
    obstable = obstable.loc[cnd_keep,:]
    msg = f"...yields nobs={len(obstable)}"
    logger.info(msg)
    #
    #-- restrict temporal range for obs (potentially)
    #
    if args.obs_lastday!=None:
        if args.obs_lastday<obs_lastday:
            obs_lastday = args.obs_lastday
        if obs_lastday!=days[-1]:
            msg = f"restricting obstable until {obs_lastday}"
            logger.info(msg)
            cnd_date = obstable.time.dt.date<=obs_lastday.date()
            obstable = obstable.loc[cnd_date,:]
            msg = f"...yields nobs={len(obstable)}"
            logger.info(msg)
    else:
        obs_lastday = days[-1]
    #
    #-- restrict to selected stations (potentially)
    #
    if args.stations!=None:
        msg = f"restricting obstable to stations ==>{args.stations}<=="
        logger.info(msg)
        cnd_sta = obstable.loc[:,'obsid'].isin(args.stations)
        obstable = obstable.loc[cnd_sta,:]
        msg = f"...yields nobs={len(obstable)}"
        logger.info(msg)
    #
    #--
    #
    obsid_values = obstable.loc[:,'obsid'].values
    if len(obsid_values[0].split('_'))==2:
        obstable.loc[:,'obs_stationid'] = obsid_values
    elif len(obsid_values[0].split('_'))==3:
        obstable.loc[:,'obs_stationid'] = ['_'.join(_.split('_')[:2]) for _ in obsid_values]
    #
    #-- assume coordinates and alitude do not depend on time
    #
    station_table = obstable.sort_values('obs_stationid')[['obs_stationid','time','lon','lat','alt',]].groupby('obs_stationid').first()
    staname_list = list(station_table.index)
    nsta = len(station_table)
    msg = f"preparing for nsta={nsta} (==>{staname_list}<==)"
    logger.info(msg)
    if nsta<=4:
        obsid_tag = '--'.join(staname_list)
    else:
        obsid_tag = f"{nsta}-obslocations"
    #
    #--
    #
    nobs = len(obstable)
    obs_dates = sorted(obstable.time.dt.date.drop_duplicates())
    obsday_tag = f"obs-{obs_dates[0]}--{obs_dates[-1]}"
    #
    #-- emissions are always from Oct 1, 2020
    #
    emis_start = Timestamp(2020, 10, 1) #--SHOULD NOT BE HARD-CODED
    emisday_range = date_range(emis_start, obs_lastday, freq='1D')
    nemisday = len(emisday_range)
    msg = f"detected maximal nemisday={nemisday}"
    logger.info(msg)
    emis_tag = emis_start.strftime(f"emis-start-%Y%m%d")
    #-- prepare for monthly emissions (if possible)
    emismon_range = date_range(emis_start, obs_lastday, freq='MS')
    nemismon = len(emismon_range)
    #
    #-- load region table
    #
    region_table = get_fitic_region_table()
    ng = 0
    regionid_1D = []
    lon_1D = None
    lat_1D = None
    for region,region_info in region_table.items():
        keep_mask = ~region_info.drop_mask
        ng_reg = np.count_nonzero(keep_mask)
        ng += ng_reg
        regionid_1D += [region,]*ng_reg
        lon_reg = region_info.lonmesh[keep_mask]
        lat_reg = region_info.latmesh[keep_mask]
        if lon_1D is None:
            lon_1D = lon_reg
        else:
            lon_1D= np.hstack((lon_1D,lon_reg))
        if lat_1D is None:
            lat_1D = lat_reg
        else:
            lat_1D= np.hstack((lat_1D,lat_reg))
    #
    #-- load footprints
    #
    if not args.pickle_filepath.exists():
        msg = f"serialised footprints ***{str(args.pickle_filepath)}*** not found."
        raise RuntimeError(msg)
    elif args.pickle_filepath.suffix=='.pickle':
        with open(args.pickle_filepath, 'rb') as fid:
            msg = f"start reading ***{str(args.pickle_filepath)}***"
            logger.info(msg)
            footprints = pickle.load(fid)
            msg = f"...pickle file loaded"
            logger.info(msg)
    elif args.pickle_filepath.suffix=='.xz':
        with lzma.open(args.pickle_filepath, 'rb') as fid:
            msg = f"start reading ***{str(args.pickle_filepath)}***"
            logger.info(msg)
            footprints = pickle.load(fid)
            msg = f"...pickle file loaded"
            logger.info(msg)
    else:
        msg = f"serialised footprints ***{str(args.pickle_filepath)}*** with unexpected suffix!"
        raise RuntimeError(msg)
    ##################################################
    #
    #--       c o l l e c t   f o o t p r i n t s   f o r   J a c o b i a n
    #
    msg = f"start preparation of obs-Jacobian for nobs/nemisday/ng={nobs}/{nemisday}/{ng}..."
    logger.info(msg)
    #
    missval = -99999.0
    obs_jacobian_units = 'ppb/(kgCH4/cell/s)' #-- by TM5
    obs_jacobian = np.zeros((nobs,nemisday,ng))
    iniconc_1D = np.full((nobs,), missval)
    obsconc_1D = np.full((nobs,), missval)
    tm5fwd_1D  = np.full((nobs,), missval)
    obsid_1D = []
    stationid_1D = []
    obstime_1D = []
    iobs = -1
    for obs_day in obs_dates:
        cnd_day = obstable.time.dt.date==obs_day
        obstable_day = obstable.loc[cnd_day,:].sort_values('obsid')
        nobs_day = len(obstable_day)
        # print(f"@{obs_day}, nobs={nobs_day}")
        for obs in obstable_day.itertuples():
            iobs += 1
            iniconc_1D[iobs] = obs.iniconc
            obsconc_1D[iobs] = obs.mixing_ratio
            tm5fwd_1D[iobs]  = obs.tm5_fwd
            obstime_1D.append(obs.time)
            obsid_1D.append(obs.obsid)
            stationid_1D.append(obs.obs_stationid)
            # print(f"iobs={iobs} for {obs.obsid}/{obs.time}")
            for idate, date in enumerate(date_range(emis_start, obs.time, freq='D')):
                cur_footplist = []
                for region in regions:
                    cur_footp = footprints[obs.Index][region][idate,:]
                    #-- restrict to relevant grid-cells
                    drop_mask = region_table[region].drop_mask
                    cur_footp = cur_footp[~drop_mask]
                    cur_footplist.append(cur_footp)
                obs_jacobian[iobs,idate,:] = np.hstack(cur_footplist)
    #--
    msg = f"...reading footprint data done."
    logger.info(msg)
     
    ##################################################
    #
    #--       J a c o b i a n   w.r.t.   m o n t h l y   e m i s s i o n s
    #
    obs_jacobian_mm = np.zeros((nobs,nemismon,ng))
    obs_jacobian_mm_units = "ppb/(kgCH4/cell/month)"
    for imon,emismondayf in enumerate(emismon_range):
        emismondayl = (emismondayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
        monday_range = date_range(emismondayf, emismondayl)
        #-- unit conversion [ppb/kgCH4/cell/s] --> [ppb/kgCH4/cell/month]
        nsecmon = len(monday_range)*86400
        cnd_emismon = (emisday_range>=emismondayf)&(emisday_range<=emismondayl)
        ndpmon = np.count_nonzero(cnd_emismon)
        nsecmon = ndpmon*86400 #-- 
        idxs_emismon = np.where(cnd_emismon)[0]
        jac_dd = obs_jacobian[:,idxs_emismon,:]
        jac_mm =  jac_dd.sum(axis=1)/nsecmon
        obs_jacobian_mm[:,imon,:] = jac_mm[:]
        
    ##################################################
    #
    #--       V E R I F I C A T I O N
    #
    if args.emission_dir!=None:
        msg = f"start reading emissions from ***{str(args.emission_dir)}***..."
        logger.info(msg)
        emis2D = np.full((nemisday,ng), missval)
        for idate, date in enumerate(emisday_range):
            emis_list = []
            for region in regions:
                drop_mask = region_table[region].drop_mask
                filepath = args.emission_dir / f"ch4emis.CH4.{region}.{date.strftime('%Y%m%d')}.nc"
                ds = xr.open_dataset(filepath)
                cur_emis = ds.to_array().sum("variable").values # [kgCH4/cell/s]
                cur_emis = cur_emis[~drop_mask]
                emis_list.append(cur_emis)
            emis2D[idate,:] = np.hstack(emis_list)
        nmiss = np.count_nonzero(emis2D==missval)
        nnan = np.count_nonzero(np.isnan(emis2D))
        msg = f"...reading emissions done (nmiss={nmiss}, nnan={nnan})"
        logger.info(msg)
        emis2D_mm = np.full((nemismon,ng), missval)
        for imon,emismondayf in enumerate(emismon_range):
            emismondayl = (emismondayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
            monday_range = date_range(emismondayf, emismondayl)
            nsecday = 86400
            cnd_emismon = (emisday_range>=emismondayf)&(emisday_range<=emismondayl)
            idxs_emismon = np.where(cnd_emismon)[0]
            emis2D_mm[imon,:] = np.sum(emis2D[idxs_emismon,:]*nsecday, axis=0)
        #
        #-- propagate emissions forward with Jacobian
        #
        obs_jac2D = obs_jacobian.reshape((nobs,nemisday*ng))
        emis1D    = emis2D.reshape(nemisday*ng)
        linfwd_1D = np.dot(obs_jac2D, emis1D) + iniconc_1D
        #
        #-- propagate monthly emissions forward with Jacobian
        #
        obs_jac2D_mm = obs_jacobian_mm.reshape((nobs,nemismon*ng))
        emis1D_mm    = emis2D_mm.reshape(nemismon*ng)
        linfwd_1D_mm = np.dot(obs_jac2D_mm, emis1D_mm) + iniconc_1D
        #
        #--
        #
        cmp_dict = {'time':obstime_1D,
                    'station':obsid_1D,
                    'tm5fwd': tm5fwd_1D,
                    'linfwd':linfwd_1D,
                    'linfwd_mm':linfwd_1D_mm,
                    'diff_lin-full':linfwd_1D-tm5fwd_1D,
                    'diff_linmm-lin':linfwd_1D_mm-linfwd_1D
                    }
        dfcmp = DataFrame.from_dict(cmp_dict)
        for col in ['diff_lin-full','diff_linmm-lin',]:
            msg = f"{col}: min/mean/max = " \
                f"{dfcmp[col].min()}/{dfcmp[col].mean()}/{dfcmp[col].max()}"
            logger.info(msg)
        outname_tokens = [f"obsjac-forward-comparison", obsday_tag, emis_tag]
        outname = '_'.join(outname_tokens) + '.csv'
        if args.outdir!=None:
            outname = args.outdir / outname
            outname.parent.mkdir(parents=True, exist_ok=True)
        dfcmp.to_csv(outname, index=False)
    ##################################################
    #
    #--       o u t p u t   g e n e r a t i o n
    #
    outname_tokens = ["fitic-inversion-input", obsid_tag, domain_tag, obsday_tag, emis_tag]
    outname = '_'.join(outname_tokens) + '.nc'
    if args.outdir!=None:
        outname = args.outdir / outname
        outname.parent.mkdir(parents=True, exist_ok=True)
    #
    #-- create dimensions
    #
    n_strlen = 32
    fp = Dataset(outname, 'w')
    fp.createDimension('ng', ng)
    if args.add_daily_obsjac:
        fp.createDimension('nemisday', nemisday)
    fp.createDimension('nemismon', nemismon)
    fp.createDimension('nobs', nobs)
    fp.createDimension('nsta', nsta)
    fp.createDimension('ntc', 6) #-- year/mon/day/hour/minute/second for calendar type variable(s)
    fp.createDimension('nstrlen', n_strlen)
    #
    #-- longitude
    #
    ncvar = fp.createVariable('lon', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'longitude'
    ncvar.units = 'degrees_east'
    ncvar.comment = 'references center of grid-cell in related zoom domain'
    ncvar[:] = lon_1D[:]
    #
    #-- latitude
    #
    ncvar = fp.createVariable('lat', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'latitude'
    ncvar.units = 'degrees_north'
    ncvar.comment = 'references center of grid-cell in related zoom domain'
    ncvar[:] = lat_1D[:]
    #
    #-- region identifier
    #
    ncvar = fp.createVariable('region', str, ('ng',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.units = ''
    ncvar[:] = np.array(regionid_1D[:])
    #
    #-- region identifier
    #
    ncvar = fp.createVariable('region_ftn', 'S1', ('ng','nstrlen',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.comment = f"region identifer in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    ncvar[:] = stringtochar(np.array(regionid_1D), n_strlen=n_strlen)
    #
    #-- observed concentration
    #
    ncvar = fp.createVariable('obs', 'f8', ('nobs',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = obsconc_1D[:]
    ncvar.long_name = f"observed CH4 concentration"
    ncvar.units = 'ppb'
    #
    #-- initial concentration
    #
    ncvar = fp.createVariable('iniconc', 'f8', ('nobs',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = iniconc_1D[:]
    ncvar.long_name = f"initial CH4 concentration"
    ncvar.units = 'ppb'
    #
    #-- station identifier (per observation)
    #
    stationid_1D = np.array(stationid_1D[:])
    ncvar = fp.createVariable('station', str, ('nobs',) )
    ncvar[:] = stationid_1D[:]
    ncvar.long_name = 'station_identifier'
    ncvar.units = ''
    ncvar = fp.createVariable('station_ftn', 'S1', ('nobs','nstrlen',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = stringtochar(stationid_1D[:], n_strlen=n_strlen)
    ncvar.long_name = 'station_identifier'
    ncvar.comment = f"station identifier in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    #
    #-- observational time points
    #
    ncvar = fp.createVariable('obstime', str, ('nobs',) )
    ncvar[:] = np.array([ _.strftime('%Y%m%dT%H%M%S') for _ in obstime_1D ])
    ncvar.long_name = 'time_of_observation'
    ncvar.units = ''
    #
    #-- observational calendar (to ease integration in Fortran inversion environment)
    #
    ncvar = fp.createVariable('obs_calendar', 'i4', ('nobs','ntc'),
                              compression='zlib', complevel=complevel)
    for iobs,_obst in enumerate(obstime_1D):
        ncvar[iobs,:] = [_obst.year,_obst.month,_obst.day,_obst.hour,_obst.minute,_obst.second]
    ncvar.long_name = 'time_of_observation'
    ncvar.units = ''
    ncvar.comment = f"observational time points in a format which is suitable " \
        f"for Fortran based I/O"
    #
    #-- unique list of stations
    #
    ncvar = fp.createVariable('station_id', str, ('nsta',))
    ncvar[:] = np.array(staname_list[:])
    ncvar.long_name = f"station_identifier_list"
    ncvar.units = ''
    ncvar.comment = f"Comprises the overall list of stations. Note, that there may be no observations for a station on certain day(s)."
    ncvar = fp.createVariable('station_id_ftn', 'S1', ('nsta','nstrlen',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = stringtochar(np.array(staname_list[:]), n_strlen=n_strlen)
    ncvar.long_name = 'station_identifier'
    ncvar.comment = f"station identifier in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    #-- longitude
    ncvar = fp.createVariable('station_lon', 'f8', ('nsta',))
    ncvar[:] = station_table.loc[:,'lon']
    ncvar.long_name = 'station_longitude'
    ncvar.units = 'degrees_east'
    #-- latitude
    ncvar = fp.createVariable('station_lat', 'f8', ('nsta',))
    ncvar[:] = station_table.loc[:,'lat']
    ncvar.long_name = 'station_longitude'
    ncvar.units = 'degrees_north'
    #-- altitude
    ncvar = fp.createVariable('station_alt', 'f8', ('nsta',))
    ncvar[:] = station_table.loc[:,'alt']
    ncvar.long_name = 'station_altitude'
    ncvar.units = 'm'
    if args.add_daily_obsjac:
        #
        #-- emission day
        #
        ncvar = fp.createVariable('emisday', str, ('nemisday',))
        ncvar.long_name = 'emission_month'
        ncvar.units = ''
        ncvar[:] = np.array([ _.strftime('%Y%m%d') for _ in emisday_range ])
        #
        #-- emission day (as calendar variable)
        #
        ncvar = fp.createVariable('emisday_calendar', 'i4', ('nemisday','ntc',),
                                  compression='zlib', complevel=complevel)
        for imon,_mon in enumerate(emisday_range):
            ncvar[imon,:] = [_mon.year,_mon.month,_mon.day,0,0,0]
        ncvar.long_name = 'emission_month_calendar'
        ncvar.comment = f"emission month information in a format which is suitable " \
            f"for Fortran based I/O"
        ncvar.units = ''
    #
    #-- emission month
    #
    ncvar = fp.createVariable('emismon', str, ('nemismon',))
    ncvar.long_name = 'emission_month'
    ncvar.units = ''
    ncvar[:] = np.array([ _.strftime('%Y%m%d') for _ in emismon_range ])
    #
    #-- emission month (as calendar variable)
    #
    ncvar = fp.createVariable('emismon_calendar', 'i4', ('nemismon','ntc',),
                              compression='zlib', complevel=complevel)
    for imon,_mon in enumerate(emismon_range):
        ncvar[imon,:] = [_mon.year,_mon.month,_mon.day,0,0,0]
    ncvar.long_name = 'emission_month_calendar'
    ncvar.comment = f"emission month information in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    #
    #-- (monthly) observational Jacobian
    #
    ncvar = fp.createVariable('obs_jacobian', 'f8', ('nobs','nemismon','ng',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = obs_jacobian_mm[:,:,:]
    ncvar.units = obs_jacobian_units
    ncvar.comment = f"Jacobian quantifies the sensitivity of concentration at " \
        f"observed times and locations w.r.t. to monthly total emissions."
    if args.add_daily_obsjac:
        #
        #-- daily observational Jacobian
        #
        ncvar = fp.createVariable('obs_jacobian_daily', 'f8', ('nobs','nemisday','ng',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = obs_jacobian[:,:,:]
        ncvar.units = obs_jacobian_units
        ncvar.comment = f"Jacobian quantifies the sensitivity of concentration at " \
            f"observed times and locations w.r.t. to daily emission rates."
    #
    #-- global attributes
    #
    fp.obstable_filepath = str(args.obsfile_filepath)
    fp.footprint_pickle_filepath = str(args.pickle_filepath)
    # fp.time_coverage_start = day_first.strftime('%Y-%m-%d')
    # fp.time_coverage_end   = day_last.strftime('%Y-%m-%d')
    # fp.time_coverage_resolution = "P1M"
    try:
        fp.processing_platform = f"{os.environ['USER']}@{os.environ['HOSTNAME']}"
    except KeyError:
        pass
    fp.history = f"{' '.join(sys.argv)}"
    fp.date_created = Timestamp.now('UTC').isoformat()
    #
    #-- close
    #
    fp.close()
    msg = f"generated file ***{outname}***"
    logger.info(msg)

    
################################################################################
#
#                   p a r s e r
#
parser = ArgumentParser(usage=globals()['__doc__'])
parser = ArgumentParser()
#----------------------------
#     s u b c o m m a n d s
#
subparsers = parser.add_subparsers( title='Available Subcommands',
                                    metavar='CMDS',
                                    description='',
                                    dest='subcmds',
                                    help='')

#
#--       prepare_daily_obsjacobian
#
sparser = subparsers.add_parser('prepare_daily_obsjacobian')
sparser.add_argument('pickle_filepath',
                    type=Path,
                    help="""use serialized footprint information from previous run (and don't reparse footprints.""")
sparser.add_argument('obsfile_filepath',
                    type=Path,
                    help="""csv file providing observational information (time and station identifier) as well as TM5 initial concentration at observational sites.""")
sparser.add_argument('--obs_lastday',
                     type=Timestamp,
                     default=Timestamp(2021,1,1),
                     help="""last observational day (default:%(default)s).""")
sparser.add_argument('--stations',
                     nargs='+',
                     help="""restrict to selected stations.""")
sparser.add_argument('--emission_dir',
                     type=Path,
                     help="""propagate emissions forward (and compare against reference forward results).""")
sparser.add_argument('--add_daily_obsjac',
                     action='store_true',
                     help="""whether to add the daily observational Jacobian to NetCDF ouput (which is currently not used in the inversion environment).""")
sparser.add_argument('--outdir',
                     type=Path,
                     help="""destination directory for all generated outputs.""")



################################################################################
#
#                   p r o g r a m   s t a r t
#
def main(args):

    ts = Timestamp.now('UTC')

    if args.subcmds=='prepare_daily_obsjacobian':
        subcmd_prepare_daily_obsjacobian(args)

    #
    te = Timestamp.now('UTC')
    msg = f"...subcommand +++{args.subcmds}+++ DONE (time_elapsed={te-ts})"
    logger.info(msg)
#
if __name__ == '__main__':
    import datetime as dtm

    progname = os.path.basename(__file__)

    #-----------------------------
    #          P R O G R A M   S T A R T
    #
    ttstart = Timestamp.now('UTC').isoformat()
    logger.info(f"{progname}::PROGRAM START::{ttstart}")
    argv = ' '.join(sys.argv)
    logger.info(f"  command-line -->{argv}<--")

    #
    #          p a r s e   c o m m a n d   l i n e
    #
    args = parser.parse_args()

    #
    #--        s t a r t   e x e c u t i o n
    #
    main(args)
