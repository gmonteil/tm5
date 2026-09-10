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
import hvplot.pandas
import hvplot.xarray
from holoviews import opts, Overlay
import holoviews as hv
import geoviews.feature as gf

from tm5.fitic import read_obs_table
from tm5.gridtools import TM5Grids
from tm5.observations import read_obspack_file
from tm5.post.footprint_io import load_adjoint_fwd #-- this was for earlier diagnostics
from tm5.fitic import tm5emisdir_load_emissions2D
from tm5.fitic import get_fitic_region_table
from tm5.fitic import ojac_glb6x4_redistribute_to_fitic
from tm5.fitic import ojac_glb6x4_redistribute_to_fitic_sqm
from tm5.post.plot_util import cnorm_set
from tm5.post.utilities import lonstr,latstr,set_outname,create_sha512


def subcmd_prepare_obsjacobian(args : ArgumentNamespace) -> None: 
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
    #-- assume coordinates and altitude do not depend on time
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
    #-- load region table ***FIT-IC compliant***
    #
    fitic_region_table = get_fitic_region_table()
    fitic_regions = list(fitic_region_table.keys())
    ng = 0
    regionid_1D = []
    lon_1D = None
    lat_1D = None
    for region,region_info in fitic_region_table.items():
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
    if domain_tag=='gns1x1':
        obs_jacobian = np.zeros((nobs,nemisday,ng))
    elif domain_tag=='glb6x4':
        _nlat = fitic_region_table['glb600x400'].grid.nlat
        _nlon = fitic_region_table['glb600x400'].grid.nlon
        ng_6x4 = _nlat*_nlon
        obs_jacobian_6x4 = np.zeros((nobs,nemisday,_nlat,_nlon))
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
                if domain_tag=='gns1x1':
                    cur_footplist = []
                    for region in regions:
                        cur_footp = footprints[obs.Index][region][idate,:]
                        #-- restrict to relevant grid-cells
                        drop_mask = fitic_region_table[region].drop_mask
                        keep_mask = ~drop_mask
                        cur_footp = cur_footp[keep_mask]
                        cur_footplist.append(cur_footp)
                    obs_jacobian[iobs,idate,:] = np.hstack(cur_footplist)
                elif domain_tag=='glb6x4':
                    cur_footp = footprints[obs.Index]['glb600x400'][idate,:]
                    obs_jacobian_6x4[iobs,idate,:] = cur_footp[:]
    #--
    msg = f"...reading footprint data done."
    logger.info(msg)
     
    ##################################################
    #
    #--       r e d i s t r i b u t i o n   o f   f l a s k   J a c o b i a n
    #
    if domain_tag=='glb6x4':
        msg = f"flask footprints computed for {domain_tag} require " \
            f"spatial re-distribution to FIT-IC grid-cells."
        logger.info(msg)
        #
        #--
        #
        ojac6x4_da = xr.DataArray(
            obs_jacobian_6x4,
            dims=('obs','emisday','lat','lon'),
            coords={'obs':obsid_1D,
                    'emisday':emisday_range,
                    'lat':fitic_region_table['glb600x400'].grid.latc,
                    'lon':fitic_region_table['glb600x400'].grid.lonc
                    },
            attrs = {'units': 'ppb/(kgCH4/cell/s)'}
            )
        # obs_jacobian = ojac_glb6x4_redistribute_to_fitic_sqm(ojac6x4_da, fitic_region_table)
        obs_jacobian = ojac_glb6x4_redistribute_to_fitic(ojac6x4_da, fitic_region_table)
        obs_jacobian_6x4 = obs_jacobian_6x4.reshape((nobs,nemisday,ng_6x4))

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
    if domain_tag=='glb6x4':
        obs_jacobian_6x4_mm = np.zeros((nobs,nemismon,ng_6x4))
        obs_jacobian_6x4_mm_units = "ppb/(kgCH4/cell/month)"
        for imon,emismondayf in enumerate(emismon_range):
            emismondayl = (emismondayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
            monday_range = date_range(emismondayf, emismondayl)
            #-- unit conversion [ppb/kgCH4/cell/s] --> [ppb/kgCH4/cell/month]
            nsecmon = len(monday_range)*86400
            cnd_emismon = (emisday_range>=emismondayf)&(emisday_range<=emismondayl)
            ndpmon = np.count_nonzero(cnd_emismon)
            nsecmon = ndpmon*86400 #-- 
            idxs_emismon = np.where(cnd_emismon)[0]
            jac_dd = obs_jacobian_6x4[:,idxs_emismon,:]
            jac_mm =  jac_dd.sum(axis=1)/nsecmon
            obs_jacobian_6x4_mm[:,imon,:] = jac_mm[:]
        
        
    ##################################################
    #
    #--       V E R I F I C A T I O N
    #
    if args.emission_dir!=None:
        msg = f"start reading emissions from ***{str(args.emission_dir)}***..."
        logger.info(msg)
        #
        #-- emissions on FIT-IC grid-cells
        #
        msg = f"fitic_regions -->{fitic_regions}<--"
        logger.debug(msg)
        emis_info = tm5emisdir_load_emissions2D(args.emission_dir, 'ch4emis', emisday_range, fitic_regions, drop=True)
        emis2D = emis_info.emis2D
        lonc1D_fitic = emis_info.lonc1D
        latc1D_fitic = emis_info.latc1D
        nnan = np.count_nonzero(np.isnan(emis2D))
        msg = f"...reading emissions done nnan={nnan})"
        logger.info(msg)
        emis2D_mm = np.full((nemismon,ng), missval)
        for imon,emismondayf in enumerate(emismon_range):
            emismondayl = (emismondayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
            monday_range = date_range(emismondayf, emismondayl)
            nsecday = 86400
            cnd_emismon = (emisday_range>=emismondayf)&(emisday_range<=emismondayl)
            idxs_emismon = np.where(cnd_emismon)[0]
            emis2D_mm[imon,:] = np.sum(emis2D[idxs_emismon,:]*nsecday, axis=0)

        if domain_tag=='glb6x4':
            emis_info = tm5emisdir_load_emissions2D(args.emission_dir, 'ch4emis', emisday_range, regions, drop=False)
            emis2D_6x4 = emis_info.emis2D
            lonc1D_6x4 = emis_info.lonc1D
            latc1D_6x4 = emis_info.latc1D
            nnan = np.count_nonzero(np.isnan(emis2D_6x4))
            msg = f"...reading emissions done nnan={nnan})"
            logger.info(msg)
            emis2D_6x4_mm = np.full((nemismon,ng_6x4), missval)
            for imon,emismondayf in enumerate(emismon_range):
                emismondayl = (emismondayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
                monday_range = date_range(emismondayf, emismondayl)
                nsecday = 86400
                cnd_emismon = (emisday_range>=emismondayf)&(emisday_range<=emismondayl)
                idxs_emismon = np.where(cnd_emismon)[0]
                emis2D_6x4_mm[imon,:] = np.sum(emis2D_6x4[idxs_emismon,:]*nsecday, axis=0)
                #--- DEBUG
                msg = f"emissions@imon={imon}: fitic/glb6x4 = " \
                    f"{emis2D_mm[imon,:].sum()}/{emis2D_6x4_mm[imon,:].sum()}" \
                    f"[kgCH4/month]"
                logger.debug(msg)
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
        #-- propagation with raw 6x4 emissions
        #
        if domain_tag=='glb6x4':
            obs_jac2D_6x4 = obs_jacobian_6x4.reshape((nobs,nemisday*ng_6x4))
            emis1D_6x4    = emis2D_6x4.reshape(nemisday*ng_6x4)
            linfwd_6x4_1D = np.dot(obs_jac2D_6x4, emis1D_6x4) + iniconc_1D
            # monthly
            obs_jac2D_6x4_mm = obs_jacobian_6x4_mm.reshape((nobs,nemismon*ng_6x4))
            emis1D_6x4_mm    = emis2D_6x4_mm.reshape(nemismon*ng_6x4)
            linfwd_6x4_1D_mm = np.dot(obs_jac2D_6x4_mm, emis1D_6x4_mm) + iniconc_1D
            # ##################################################
            # ### MVO-DEBUG tracing delta-concentration differences...
            # ###
            # iemisday = 0
            # iobs = 0
            # cur_emis2D_fitic = emis2D[iemisday,:]
            # cur_ojac_fitic = obs_jacobian[iobs,iemisday,:]
            # cur_emis2D_6x4 = emis2D_6x4[iemisday,:]
            # cur_ojac_6x4 = obs_jacobian_6x4[iobs,iemisday,:]
            # ### MVO-DEBUG
            # dconc_6x4 = []
            # dconc_fitic = []
            # ngc_fitic = []
            # idxs_fitic = np.array([],dtype='i4')
            # for ig_6x4 in range(ng_6x4):
            #     _lonc = lonc1D_6x4[ig_6x4]
            #     _latc = latc1D_6x4[ig_6x4]
            #     cnd_lon = (lonc1D_fitic>=_lonc-3)&(lonc1D_fitic<=_lonc+3)
            #     cnd_lat = (latc1D_fitic>=_latc-2)&(latc1D_fitic<=_latc+2)
            #     cur_idxs_fitic = np.where(cnd_lon&cnd_lat)[0]
            #     idxs_fitic = np.hstack((idxs_fitic, cur_idxs_fitic))
            #     _ng_fitic = len(cur_idxs_fitic)
            #     ngc_fitic.append(_ng_fitic)
            #     msg = f"@ig_6x4={ig_6x4} lon/lat = {_lonc}/{_latc}, iemisday={iemisday} " \
            #         f" _ng_fitic={_ng_fitic}"
            #     _emis_6x4 = cur_emis2D_6x4[ig_6x4]
            #     _emis_fitic = cur_emis2D_fitic[cur_idxs_fitic]
            #     _emis_fitic_sum = np.sum(_emis_fitic)
            #     if max(_emis_6x4,_emis_fitic_sum) > 0:
            #         _rdiff = abs(_emis_fitic_sum-_emis_6x4)/_emis_fitic_sum
            #         if _rdiff>1e-10:
            #             msg = f"@ig_6x4={ig_6x4} lon/lat = {_lonc}/{_latc}" \
            #                 f"emission glb6x4/fitic/rdiff = {_emis_6x4}/{_emis_fitic_sum}/{_rdiff}"
            #             print(msg)
            #     else:
            #         pass
            #     _ojac_6x4 = cur_ojac_6x4[ig_6x4]
            #     _ojac_fitic = cur_ojac_fitic[cur_idxs_fitic]
            #     _ojac_fitic_sum = np.sum(_ojac_fitic)
            #     _rdiff = abs(_ojac_fitic_sum-_ojac_6x4)
            #     # msg = f"@ig_6x4={ig_6x4} lon/lat = {_lonc}/{_latc}" \
            #     #     f"iobs/iemisday={iobs}/{iemisday}: " \
            #     #     f"Jacobian values glb6x4/fitic = {_ojac_6x4}/{_ojac_fitic_sum}"
            #     # print(msg)
            #     _dconc_6x4 = _ojac_6x4*_emis_6x4
            #     _dconc_fitic = np.sum(_ojac_fitic*_emis_fitic)
            #     if _dconc_6x4!=_dconc_fitic:
            #         msg = f"@ig_6x4={ig_6x4} lon/lat = {_lonc}/{_latc}" \
            #             f"iobs/iemisday={iobs}/{iemisday}: " \
            #             f"Jacobian values glb6x4/fitic = {_ojac_6x4}/{_ojac_fitic_sum}"
            #         print(msg)
            #         msg = f"@ig_6x4={ig_6x4} lon/lat = {_lonc}/{_latc}" \
            #             f"iobs/iemisday={iobs}/{iemisday}: " \
            #             f"_dconc_6x4/_dconc_fitic = {_dconc_6x4}/{_dconc_fitic}"
            #         print(msg)
            #         msg = f"...cur_idxs_fitic={cur_idxs_fitic}"
            #         print(msg)
            #     dconc_6x4.append(_dconc_6x4)
            #     dconc_fitic.append(_dconc_fitic)
            # #
            # _data_dict = {'lonc_6x4': lonc1D_6x4,
            #               'latc_6x4': latc1D_6x4,
            #               'ngc_fitic': ngc_fitic,
            #               'dconc_6x4':dconc_6x4,
            #               'dconc_fitic':dconc_fitic,
            #               }
            # _df = pd.DataFrame.from_dict(_data_dict)
            # _df.loc[:,'dconc_diff'] = _df.loc[:,'dconc_fitic'] - _df.loc[:,'dconc_6x4']
            # print(_df[['dconc_6x4','dconc_fitic','dconc_diff']].describe())
            # msg = f"iobs/iemisday={iobs}/{iemisday}: " \
            #     f"ngc_fitic/dconc_fitic/dconc_6x4 = " \
            #     f"{_df['ngc_fitic'].sum()}/{_df['dconc_fitic'].sum()}/{_df['dconc_6x4'].sum()}"
            # print(msg)
            # outname = f"dconc_debug-comparison_iobs{iobs}_iemisday{iemisday}.csv"
            # _df.to_csv(outname, index=False)
            # print(f"-"*50)
            # print(f"-"*50)
            # dconc_fitic = 0
            # msg = f"len(idxs_fitic)={len(idxs_fitic)}"
            # print(msg)
            # for ig in idxs_fitic:
            #     dconc_fitic += cur_ojac_fitic[ig]*cur_emis2D_fitic[ig]
            # msg = f"loop-computed dconc_fitic = {dconc_fitic}"
            # print(msg)
            # print(f"-"*50)
            # print(f"-"*50)
            # msg = f"iemisday={iemisday}: " \
            #     f"np.sum(cur_emis2D_fitic)={np.sum(cur_emis2D_fitic[:])}"
            # print(msg)
            # msg = f"iemisday={iemisday}: " \
            #     f"np.sum(cur_emis2D_6x4[:])={np.sum(cur_emis2D_6x4[:])}"
            # print(msg)
            # msg = f"iobs/iemisday={iobs}/{iemisday}: " \
            #     f"np.sum(cur_ojac_fitic[:])={np.sum(cur_ojac_fitic[:])}"
            # print(msg)
            # msg = f"iobs/iemisday={iobs}/{iemisday}: " \
            #     f"np.sum(cur_ojac_6x4[:])={np.sum(cur_ojac_6x4[:])}"
            # print(msg)
            # dconc = np.sum(cur_ojac_fitic[:]*cur_emis2D_fitic[:])
            # dconc_6x4 = np.sum(cur_ojac_6x4[:]*cur_emis2D_6x4[:])
            # msg = f"iobs/iemisday={iobs}/{iemisday}: " \
            #     f"dconc/dconc_6x4 = {dconc}/{dconc_6x4}"
            # print(msg)
            # ### MVO-DEBUG-END tracing delta-concentration differences...
            # ##################################################
        #
        #--
        #
        if domain_tag=='gns1x1':
            cmp_dict = {'time':obstime_1D,
                        'station':obsid_1D,
                        'iniconc':iniconc_1D,
                        'tm5fwd': tm5fwd_1D,
                        'linfwd':linfwd_1D,
                        'linfwd_mm':linfwd_1D_mm,
                        'diff_lin-full':linfwd_1D-tm5fwd_1D,
                        'diff_linmm-lin':linfwd_1D_mm-linfwd_1D
                        }
        elif domain_tag=='glb6x4':
             cmp_dict = {'time':obstime_1D,
                         'station':stationid_1D,
                         'iniconc':iniconc_1D,
                         'tm5fwd': tm5fwd_1D,
                         'linfwd':linfwd_1D,
                         'linfwd_mm':linfwd_1D_mm,
                         'diff_lin-full':linfwd_1D-tm5fwd_1D,
                         'diff_linmm-lin':linfwd_1D_mm-linfwd_1D,
                         'linfwd_6x4': linfwd_6x4_1D,
                         'linfwd_mm_6x4': linfwd_6x4_1D_mm,
                         'diff_lin-full_6x4': linfwd_6x4_1D-tm5fwd_1D,
                         'diff_linmm-lin_6x4': linfwd_6x4_1D_mm-linfwd_6x4_1D
                        }
        dfcmp = DataFrame.from_dict(cmp_dict)
        for col in ['diff_lin-full','diff_linmm-lin',]:
            msg = f"{col}: min/mean/max = " \
                f"{dfcmp[col].min()}/{dfcmp[col].mean()}/{dfcmp[col].max()}"
            logger.info(msg)
        outname_tokens = [f"obsjac-forward-comparison", obsid_tag, domain_tag, obsday_tag, emis_tag]
        outname = '_'.join(outname_tokens) + '.csv'
        if args.outdir!=None:
            outname = args.outdir / outname
            outname.parent.mkdir(parents=True, exist_ok=True)
        dfcmp.to_csv(outname, index=False)
        msg = f"generated ***{str(outname)}***"
        logger.info(msg)
        #
        #-- and now plotting
        #
        rename_dict = {'linfwd':'linfwd_daily', 'linfwd_mm':'linfwd_monthly',
                       'diff_lin-full':'linfwd-tm5fwd',
                       'diff_linmm-lin':'linfwd_monthly-daily'}
        dfcmp = dfcmp.rename(rename_dict,axis=1)
        cols_fwd = ['tm5fwd', 'iniconc','linfwd_daily','linfwd_monthly']
        cols_diff = ['linfwd-tm5fwd','linfwd_monthly-daily']
        min_fwd = dfcmp.loc[:,cols_fwd].min().min()
        max_fwd = dfcmp.loc[:,cols_fwd].max().max()
        min_diff = dfcmp.loc[:,cols_diff].min().min()
        max_diff = dfcmp.loc[:,cols_diff].max().max()
        #
        outname = '_'.join(outname_tokens) + '.html'
        if args.outdir!=None:
            outname = args.outdir / outname
            outname.parent.mkdir(parents=True, exist_ok=True)
        p = (
            dfcmp.hvplot(x='time', y=cols_fwd, groupby='station', width=1500, height=600, grid=True, ylim=(min_fwd,max_fwd)) +
            dfcmp.hvplot(x='time', y=cols_diff, groupby='station', width=1500, height=600, ylim=(min_diff,max_diff))
        ).cols(1)
        hv.save(p, outname)
        msg = f"generated ***{str(outname)}***"
        logger.info(msg)
        #>> only forward
        xoutname_tokens = [f"obsjac-forward", obsid_tag, domain_tag, obsday_tag, emis_tag]
        outname = '_'.join(xoutname_tokens) + '.html'
        if args.outdir!=None:
            outname = args.outdir / outname
            outname.parent.mkdir(parents=True, exist_ok=True)
        p = (
            dfcmp.hvplot(x='time', y=cols_fwd, groupby='station', width=1500, height=600, grid=True)
        )
        hv.save(p, outname)
        msg = f"generated ***{str(outname)}***"
        logger.info(msg)
        #>> only differences
        xoutname_tokens = [f"obsjac-differences", obsid_tag, domain_tag, obsday_tag, emis_tag]
        outname = '_'.join(xoutname_tokens) + '.html'
        if args.outdir!=None:
            outname = args.outdir / outname
            outname.parent.mkdir(parents=True, exist_ok=True)
        p = (
            dfcmp.hvplot(x='time', y=cols_diff, groupby='station', width=1500, height=600, grid=True)
        )
        hv.save(p, outname)
        msg = f"generated ***{str(outname)}***"
        logger.info(msg)

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
        ncvar.long_name = 'day_of_emission'
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
    ncvar.units = obs_jacobian_mm_units
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

    
def subcmd_monthly_emissions_for_inversion(args : ArgumentNamespace) -> None:
    """
    Preparation of monthly averaged emissions suitable as input for (Fortran based)
    inversion system.
    """
    tm5emisdir = args.tm5emisdir
    time_start, time_end = args.time_range
    regions = args.regions
    complevel = args.__dict__.get('complevel',4)

    #
    #--
    #
    mon_range = date_range(time_start, time_end, freq='MS')
    monend_range = date_range(time_start, time_end, freq='ME')
    nmon = len(mon_range)
    ntc = 3 #-- recording year/month/day
    time_data = np.full((nmon,ntc), -1)
    for imon,mon in enumerate(mon_range):
        time_data[imon,:] = [mon.year, mon.month, mon.day] #-- deliberately take first day
    day_first = mon_range[0]
    day_last  = monend_range[-1]
    month_tag = f"{day_first.strftime('%Y%m%d')}--{day_last.strftime('%Y%m%d')}"
    day_range = date_range(day_first, day_last, freq='1D')

    #
    #-- initialise array for emissions
    #
    nsecday = 86400
    emis_miss = -99999.
    emis_data = None
    reginfo = None
    ng = None
    #
    #-- monthly total emissions [kgCH4/cell/month
    #
    for imon,dayf in enumerate(mon_range):
        dayl = (dayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
        day_range = date_range(dayf,dayl)
        msg = f"...loading emissions for {dayf.strftime('%Y%m%d')} to {dayl.strftime('%Y%m%d')}"
        logger.info(msg)
        #-- load daily emissions for every day in month
        drop = len(regions)>1
        emis_info = tm5emisdir_load_emissions2D(tm5emisdir, 'ch4emis', day_range, regions, drop=drop)
        #-- collect ancillary infos, allocate array
        if imon==0:
            _,ng = emis_info.emis2D.shape
            emis_data = np.full((nmon,ng), emis_miss)
            reginfo = emis_info
        #-- convert daily emission rates [kgCH4/cell/s] to [kgCH4/cell/month]
        emis_mm =  np.sum(emis_info.emis2D*nsecday, axis=0)
        #-- insert current month into buffer[mon,grid]
        emis_data[imon,:] = emis_mm[:]
        
    msg = f"...monthly emission data ready."
    logger.info(msg)

    #
    #-- output preparation
    #
    region_tag = "-".join(regions)
    outname_tokens = [f"fitic-monthly-emissions", month_tag, region_tag,]
    outname = '_'.join(outname_tokens) + '.nc'
    outname = set_outname(args, outname)
    msg = f"writing emission inputs for inversion inputs to file ***{outname}***..."
    logger.info(msg)
    #
    #-- spatial dimensions
    #
    fp = Dataset(outname, 'w')
    n_strlen = 32
    fp.createDimension('ntc', ntc)
    fp.createDimension('nstrlen', n_strlen)
    fp.createDimension('ng', ng)
    fp.createDimension('nmon', nmon)
    #-- time variable
    ncvar = fp.createVariable('time', 'i4', ('nmon','ntc',))
    ncvar.long_name = "date_of_first_day_in_month"
    ncvar.units = ''
    ncvar[:] = time_data[:]
    #
    #-- longitude
    #
    ncvar = fp.createVariable('lon', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'longitude'
    ncvar.units = 'degrees_east'
    ncvar.comment = 'references center of grid-cell in underlying domain'
    ncvar[:] = reginfo.lonc1D
    #
    #-- latitude
    #
    ncvar = fp.createVariable('lat', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'latitude'
    ncvar.units = 'degrees_north'
    ncvar.comment = 'references center of grid-cell in underlying domain'
    ncvar[:] = reginfo.latc1D
    #
    #-- area
    #
    ncvar = fp.createVariable('area', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'gridcell_area'
    ncvar.units = 'm2'
    ncvar[:] = reginfo.area1D
    #
    #-- region identifier
    #
    ncvar = fp.createVariable('region', reginfo.reg1D.dtype, ('ng',))
    ncvar.long_name = f"gridcell_region_identifier"
    ncvar.units = ''
    ncvar[:] = reginfo.reg1D[:]
    #
    #-- region identifier (Fortran compliant)
    #
    ncvar = fp.createVariable('region_ftn', 'S1', ('ng','nstrlen'))
    ncvar.long_name = f"gridcell_region_identifier"
    ncvar.comment = f"region identifer in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    ncvar[:] = stringtochar(reginfo.reg1D[:], n_strlen=n_strlen)
    #
    #-- emission variable
    #
    if nmon>1:
        ncvar = fp.createVariable('emission', 'f8', ('nmon','ng'),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = emis_data[:]
    else:
        ncvar = fp.createVariable('emission', 'f8', ('ng',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = emis_data[month-1,:]
    ncvar.long_name = "CH4 emissions"
    ncvar.units = 'kgCH4/cell/month'

    #
    #-- global attributes
    #
    fp.emission_directory = str(tm5emisdir)
    fp.time_coverage_start = day_first.strftime('%Y-%m-%d')
    fp.time_coverage_end   = day_last.strftime('%Y-%m-%d')
    fp.time_coverage_resolution = "P1M"
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


def subcmd_prepare_tgtjacobian(args : ArgumentNamespace) -> None:
    """
    """
    complevel = args.__dict__.get('complevel',4)
    #
    #--
    #
    #
    #-- load region table
    #
    region_table = get_fitic_region_table()
    regions = list(region_table.keys())
    ng = 0
    regionid_1D = []
    area_1D = None
    lon_1D = None
    lat_1D = None
    for region,region_info in region_table.items():
        keep_mask = ~region_info.drop_mask
        ng_reg = np.count_nonzero(keep_mask)
        ng += ng_reg
        regionid_1D += [region,]*ng_reg
        lon_reg = region_info.lonmesh[keep_mask]
        lat_reg = region_info.latmesh[keep_mask]
        area_reg = region_info.grid.area[keep_mask]
        if lon_1D is None:
            lon_1D = lon_reg
        else:
            lon_1D= np.hstack((lon_1D,lon_reg))
        if lat_1D is None:
            lat_1D = lat_reg
        else:
            lat_1D= np.hstack((lat_1D,lat_reg))
        if area_1D is None:
            area_1D = area_reg
        else:
            area_1D= np.hstack((area_1D,area_reg))
    regionid_1D = np.array(regionid_1D)
    msg = f"-->{regions}<-- yield overall {ng} grid-cells"
    logger.info(msg)

    #
    #--
    #
    cnd_gns_nohalo = (regionid_1D=='gns100x100')
    gns_reso = 1.
    gns_nohalo_w = np.min(lon_1D[cnd_gns_nohalo]) - gns_reso/2
    gns_nohalo_e = np.max(lon_1D[cnd_gns_nohalo]) + gns_reso/2
    gns_nohalo_s = np.min(lat_1D[cnd_gns_nohalo]) - gns_reso/2
    gns_nohalo_n = np.max(lat_1D[cnd_gns_nohalo]) + gns_reso/2
    gns_nohalo_nlon = int((gns_nohalo_e-gns_nohalo_w)//gns_reso)
    gns_nohalo_nlat = int((gns_nohalo_n-gns_nohalo_s)//gns_reso)
    gns_nohalo_area2D = area_1D[cnd_gns_nohalo].reshape(gns_nohalo_nlat,gns_nohalo_nlat)
    #
    #-- dedicated country targets *only* in the 1x1 innermost zoom domain
    #   with HALO parts removed
    #
    tgt_country_table = OrderedDict()
    if args.countryfrct_filepath!=None:
        cfrctfile = args.countryfrct_filepath
        # print(f"==>{cfrctfile}<== ({type(cfrctfile)}")
        frctds = xr.open_dataset(cfrctfile)
        reso = 1
        w = int(frctds.lon.values.min()-reso/2)
        e = int(frctds.lon.values.max()+reso/2)
        s = int(frctds.lat.values.min()-reso/2)
        n = int(frctds.lat.values.max()+reso/2)
        cgrid = TM5Grids.from_corners(west=w,east=e,south=s,north=n,dlon=reso,dlat=reso)
        country_id = frctds.country_ID.values
        msg = f"country identifiers in file ==>{country_id}<=="
        logger.info(msg)
        country_frct = frctds.country_fraction.values
        country_area = country_frct*cgrid.area[np.newaxis,:,:]
        # for ic,cid in enumerate(country_id):
        #     _area = country_area[ic,:].sum() # [m2]
        #     print(f"{cid}: {_area/1e6:.2f}[km2]")
        #
        #-- restrict to innermost zoom domain
        #
        frctds_inner = frctds.sel(lon=(frctds.lon>=gns_nohalo_w)&(frctds.lon<=gns_nohalo_e),lat=(frctds.lat>=gns_nohalo_s)&(frctds.lat<=gns_nohalo_n))
        country_frct_inner = frctds_inner.country_fraction.values
        
        country_area_inner = country_frct_inner*gns_nohalo_area2D[np.newaxis,:,:]
        for ic,cid in enumerate(country_id):
            _area    = country_area_inner[ic,:].sum()/1e6 #[km2]
            _areatot = country_area[ic,:].sum()/1e6     #[km2]
            if _area>0:
                msg = f"{cid}: area_gns-nohalo={_area:.2f}[km2] (area={_areatot:.2f}[km2])"
                print(msg)
        #
        #--
        #
        if args.countries!=None:
            for _cid in args.countries:
                #-- index
                cid = _cid.upper()
                ic = np.where(country_id==cid)[0]
                if len(ic)==0:
                    msg = f"country identifier -->{cid}<-- not found"
                    raise RuntimeError(msg)
                else:
                    if country_frct_inner[ic,:].sum()!=country_frct[ic,:].sum():
                        msg = f"country -->{cid}<-- not fully covered by innermost zoom domain, " \
                            f"needs to be ignored!"
                        logger.warning(msg)
                        continue
                    else:
                        cur_fraction = country_frct_inner[ic,:]
                        tgt_country_table[cid] = country_frct_inner[ic,:].ravel()
                        msg = f"@{cid}, #non-zero={np.count_nonzero(cur_fraction>0)}, fraction-sum={np.sum(cur_fraction)}"
                        logger.info(msg)
    #
    #-- target jacobian
    #
    target_list = ['global', 'target domain',] + list(tgt_country_table.keys())
    ntgt = len(target_list)
    tjac2D = zeros((ntgt,ng), dtype='f8')
    for itgt,tgt in enumerate(target_list):
        if tgt=='global':
            #
            #-- updated emission generation does no longer
            #   have double counted emissions (!), i.e.
            #   *every* grid-cell must now be accounted.
            #
            tjac2D[itgt,:] = 1.
        elif tgt=='target domain':
            #
            #-- contribution of all grid-cells that are simulated at 1x1 degree
            #
            tjac2D[itgt,cnd_gns_nohalo] = 1.
        elif tgt in tgt_country_table:
            tgt_frct = tgt_country_table[tgt]
            #-- just one more consistency check
            tgt_area = np.sum(tgt_frct*area_1D[cnd_gns_nohalo])/1e6
            msg = f"inserting grid cell fractions for target -->{tgt}<-- (area: {tgt_area:.2f}[km2])"
            logger.info(msg)
            tjac2D[itgt,cnd_gns_nohalo] = tgt_frct
    #
    #--
    #
    if args.emission_filepath!=None:
        emis_ds = xr.open_dataset(args.emission_filepath)
        units = "kgCH4/cell/month"
        assert emis_ds.emission.units==units
        target_unit = "kgCH4"
        emisvec = emis_ds.emission.sel(nmon=0).values
        if len(emisvec)!=ng:
            msg = f"emissions from file ***{str(args.emission_filepath)}*** are not compliant " \
                f"with target Jacobian ng={len(emisvec)} instead of expected ng={ng}"
            raise RuntimeError(msg)
        emisvec_gns = emisvec[cnd_gns_nohalo]
        msg = f"sum(emisvec_gns)={sum(emisvec_gns)}"
        logger.info(msg)
        for itgt,tgt in enumerate(target_list):
            tjac_vec = tjac2D[itgt,:]
            msg = f"@itgt={itgt}/{tgt}, sum(tjac)={np.sum(tjac_vec)}"
            logger.info(msg)
            tgt_emis = np.dot(tjac_vec, emisvec)
            msg = f"{tgt}, target_emissions={tgt_emis:.0f}[{target_unit}]"
            logger.info(msg)
        sys.exit(0)
        
    #
    #-- prepare output
    #
    outname_tokens = [f"fitic-tarjac_ntgt{ntgt}",]
    if len(tgt_country_table)>0:
        country_tag = 'with-' + '-'.join(list(tgt_country_table.keys()))
        outname_tokens.append(country_tag)
    outname = '_'.join(outname_tokens) + '.nc'
    outname = set_outname(args, outname)
    msg = f"writing inversion inputs to file ***{outname}***..."
    logger.info(msg)
    #
    #-- spatial dimensions
    #
    n_strlen = 32
    fp = Dataset(outname, 'w')
    fp.createDimension('ng', ng)
    fp.createDimension('ntgt', ntgt)
    fp.createDimension('nstrlen', n_strlen)
    #
    ncvar = fp.createVariable('lon', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'longitude'
    ncvar.units = 'degrees_east'
    ncvar.comment = 'references center of grid-cell in related target domain'
    ncvar[:] = lon_1D
    #
    ncvar = fp.createVariable('lat', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'latitude'
    ncvar.units = 'degrees_north'
    ncvar.comment = 'references center of grid-cell in related target domain'
    ncvar[:] = lat_1D
    #
    ncvar = fp.createVariable('region', regionid_1D.dtype, ('ng',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.units = ''
    ncvar[:] = regionid_1D[:]
    ncvar = fp.createVariable('region_ftn', 'S1', ('ng','nstrlen',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.comment = f"region identifer in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    ncvar[:] = stringtochar(regionid_1D[:], n_strlen=n_strlen)
    #
    ncvar = fp.createVariable('area', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.units = 'm2'
    ncvar[:] = area_1D
    #
    ncvar = fp.createVariable('targets', str, ('ntgt',))
    ncvar.long_name = f"target_identifier"
    ncvar.units = ''
    ncvar[:] =  np.array(target_list)
    #
    ncvar = fp.createVariable('targets_ftn', 'S1', ('ntgt','nstrlen',))
    ncvar.long_name = f"target_identifier"
    ncvar.comment = f"region identifer in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    ncvar[:] = stringtochar(np.array(target_list), n_strlen=n_strlen)
    #
    ncvar = fp.createVariable('tgt_jacobian', 'f8', ('ntgt','ng',),
                              compression='zlib', complevel=complevel)
    ncvar.units = ''
    ncvar[:] = tjac2D[:]

    #
    #-- global attributes
    #
    fp.description = f"Target Jacobian for Fortran inversion environment within FIT-IC"
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


def subcmd_fitic_monthly_emissions_visu(args):
    filepath = args.filepath
    varname = args.variable
    cmap = 'RdBu_r'
    clim = (-0.005, 0.005)

    emis_tag = filepath.stem.split('_')[0]
    
    ds = xr.open_dataset(filepath)

    imon = ds.nmon.values[-1]
    ds = ds.isel(nmon=imon)
    emis_mon = Timestamp(*ds.time.values)

    #
    #-- loop over regions
    #
    region_list = np.unique(ds.region.values)
    for region in region_list:
        ds_cur = ds.sel( ng=(ds.region==region) )
        emis_in = ds_cur[varname]
        assert emis_in.attrs['units'] in ['kgCH4/cell',"kgCH4/cell/month"], \
            f"unexpected emission units -->{emis_in.attrs['units']}<--"
        emis_out = (emis_in / ds_cur.area).values
        ilat = ((ds_cur.lat - ds_cur.lat.min()) / int(region[7])).values.astype(int)
        ilon = ((ds_cur.lon - ds_cur.lon.min()) / int(region[3])).values.astype(int)
        
        emis_visu = zeros((ilat.max() + 1, ilon.max() + 1))
        emis_visu[ilat, ilon] = emis_out
        emis_visu = xr.DataArray(data = emis_visu,
                                 dims=('lat', 'lon'),
                                 attrs={'units':'kgCH4/m2'})
        emis_visu['lat'] = ('lat', sorted(set(ds_cur.lat.values)))
        emis_visu['lon'] = ('lon', sorted(set(ds_cur.lon.values)))
        projection = crs.PlateCarree()
        p = emis_visu.hvplot.quadmesh(cmap=cmap, clim=clim,
                                      geo=True,
                                      crs=crs.PlateCarree(),
                                      features=['borders', 'coastline'],
                                      projection=projection,
                                      width=1500, height=600)
        title = emis_mon.strftime(f"{emis_tag} emissions ({region}, %b %Y)")
        plotcfg = opts.Overlay(title=title, ylabel=f"{emis_visu.attrs['units']}")
        p.opts(plotcfg)
        outname = f"{emis_tag}_{region}-emissions_{emis_mon.strftime('%Y-%m')}.html"
        if args.outdir!=None:
            outname = args.outdir / outname
            outname.parent.mkdir(parents=True, exist_ok=True)
        hv.save(p, outname)
        msg = f"generated ***{str(outname)}***"
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
#--       prepare_obsjacobian
#
sparser = subparsers.add_parser('prepare_obsjacobian',
                                help="""preparation of NetCDF file providing an observational Jacobian and further inputs required for FIT-IC inversion system.""")
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

#
#--       monthly_emissions_for_inversion
#
sparser = subparsers.add_parser('monthly_emissions_for_inversion',
                                help="""test preparation of inputs for Fortran inversion system.""")
sparser.add_argument('tm5emisdir',
                     help="""name of directory containing daily emissions files as prepared for TM5 simulations plus the initial part of the file name pattern.""")
sparser.add_argument('--time_range',
                     type=Timestamp,
                     nargs=2,
                     default=[Timestamp(2021,1,1), Timestamp(2021,12,31)],
                     help="""temporal range, only year and month are significant (default: %(default)s).""")
sparser.add_argument('--regions',
                     nargs='+',
                     choices=['glb600x400','eur300x200','gns100x100',],
                     default=['glb600x400','eur300x200','gns100x100',],
                     help="""selected regions (default: %(default)s), better only change for test purposes.""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       prepare_tgtjacobian
#
sparser = subparsers.add_parser('prepare_tgtjacobian',
                                help="""preparation of dedicated target Jacobian for Fortran inversion system.""")
sparser.add_argument('--countryfrct_filepath',
                     type=Path,
                     help="""provide 1degree gridded NetCDF file with country fractions.""")
sparser.add_argument('--countries',
                     nargs='+',
                     help="""select list of countries (ISO3 identifier) which must be fully covered within the FIT-IC innermost target domain.""")
sparser.add_argument('--emission_filepath',
                     type=Path,
                     help="""for debugging: provide a compliant emissions NetCDF file to check the target values.""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--
#
sparser = subparsers.add_parser('fitic_monthly_emissions_visu',
                                help="""visualisation of emissions prepared for or generated by the Fortran inversion environment.""")
sparser.add_argument('filepath',
                     type=Path,
                     help="""NetCDF file providing the emissions.""")
sparser.add_argument('--variable',
                     default='emission',
                     help="""name of NetCDF variable for the emissions (default: %(default)s).""")
sparser.add_argument('--outdir',
                     type=Path,
                     help="""top-level directory for any generated outputs..""")



################################################################################
#
#                   p r o g r a m   s t a r t
#
def main(args):

    ts = Timestamp.now('UTC')

    if args.subcmds=='prepare_obsjacobian':
        subcmd_prepare_obsjacobian(args)

    if args.subcmds=='monthly_emissions_for_inversion':
        subcmd_monthly_emissions_for_inversion(args)

    if args.subcmds=='prepare_tgtjacobian':
        subcmd_prepare_tgtjacobian(args)

    if args.subcmds=='fitic_monthly_emissions_visu':
        subcmd_fitic_monthly_emissions_visu(args)

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
