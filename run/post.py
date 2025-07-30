#!/usr/bin/env python

"""
M. Vossbeck, March-July 2025:
  initial script used for simple post-processing (visualisation, analysis,...)
  outputs of a TM5 simulation,

July 30, 2025: added under version control 
"""
#
#-- system packages
#
import sys
import os
import shutil
from omegaconf import OmegaConf
from argparse import ArgumentParser
from pathlib import Path
import datetime as dtm
from collections import OrderedDict
import netCDF4 as nc4
import xarray as xr
import pandas as pd
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from loguru import logger
from typing import Union, List, Tuple

#
#-- local packages
#
import tm5
from tm5.gridtools import TM5Grids
from tm5.post.utilities import get_hostname,  set_outname
from tm5.post.utilities import sphere_grid_find_close
from tm5.post.cams_io   import cams_at_obspack_load_conctseries


def obspack_load_conctseries( obspackdir :  Union[str, Path],
                              tcover_start : pd.Timestamp,
                              tcover_end : pd.Timestamp,
                              sta_tag : str,
                              lonq : float,
                              latq : float,
                              altq : float,
                              altdif_threshold : float = 2 ) -> Union[None,pd.DataFrame]:
    #
    #-- station tag is like: [id]_[alt]
    #   (with alt likely level above ground
    #
    obs_dict = None
    id_tag,alt_tag = sta_tag.split('_')
    fptn = f'ch4_{id_tag}*.nc'
    file_lst = list(Path(obspackdir).glob(fptn))
    if len(file_lst)==0:
        msg = f"...@{sta_tag} no matching obspack files found."
        logger.warning(msg)
        return obs_dict
    else:
        nfiles = len(file_lst)
        msg = f"...@{sta_tag}, {len(file_lst)} candidate files " \
            f"(***{[_.name for _ in file_lst]}***)"
        logger.info(msg)
        obspack_dict = None
        obspack_fpath = None
        obspack_ch4_lst = [None]*nfiles
        altdif_lst = np.full(nfiles, np.inf)
        obsalt_lst = np.full(nfiles, -1)
        #
        #-- loop over obspack files with matching station identifier
        #   -> find file where altitude is closest to request
        #
        for ipath,fpath in enumerate(file_lst):
            ch4_dict = load_obspack_ch4(fpath,
                                        date_first=tcover_start,
                                        date_last=tcover_end)
            # print(f"fpath ***{fpath}*** yields\n{ch4_dict}")
            obspack_ch4_lst[ipath] = ch4_dict
            if ch4_dict==None:
                continue
            elev = ch4_dict.get('elevation',None)
            alt  = ch4_dict.get('altitude',None) #-- that should be the measurement height
            # print(f"MVODEBUG::@{sta_tag} (height={altq}), ***{fpath.name}*** elev ==>{elev}<== alt ==>{alt}<==")
            if alt is None:
                continue #-- Hmm, obspack data file without measurement height...
            elif len(alt)==0:
                continue #-- Hmm, obspack data file without measurement height...
            #
            #-- maximal difference in altitude
            #
            if np.all(np.diff(alt)==0):
                obsalt = alt[0]
            else:
                obsalt = alt.mean()
            obsalt_lst[ipath] = obsalt
            alt_dif = np.abs(altq-obsalt)
            altdif_lst[ipath] = alt_dif
            # msg = f"...obsalt={obsalt} @{fpath}"
            # logger.debug(msg)
        #
        #-- finding best altitude match
        #
        ialtmin = np.argmin(altdif_lst)
        altdif_min = altdif_lst[ialtmin]
        msg = f"...@{sta_tag} with altq={altq}, " \
            f"yields best altitude match {obsalt_lst[ialtmin]} " \
            f"(difference={altdif_min}) with file ***{file_lst[ialtmin]}***"
        logger.info(msg)
        #
        #--
        #
        if altdif_min<altdif_threshold: #
            obspack_fpath = file_lst[ialtmin]
            obspack_dict  = obspack_ch4_lst[ialtmin]
        if obspack_fpath==None:
            msg = f"...@{sta_tag} (height={altq}), no matching obspack observation found."
            logger.warning(msg)
        else:
            msg = f"...using obspack data from file ***{obspack_fpath}*** for comparison."
            logger.info(msg)
            dfobs = pd.DataFrame.from_dict({'time':obspack_dict['time'],
                                            'obspack_ch4':obspack_dict['ch4']})
            obs_dict = {}
            obs_dict['df'] = dfobs
            obs_dict['filepath'] = obspack_fpath
            obs_dict['altitude'] = obsalt_lst[ialtmin]
    return obs_dict


def load_obspack_ch4(filepath : Union[str,Path],
                     date_first : Union[pd.Timestamp,None] = None,
                     date_last  : Union[pd.Timestamp,None] = None):
    """Reading CH4 concentrations from NetCDF file provided by obspack.
    """
    if date_first!=None and date_last!=None and not date_first<date_last:
        msg = f"inconsistent temporal selection {date_first} --- {date_last}"
        raise RuntimeError(msg)

    # msg = f"loading obspack data from ***{filepath}***..."
    # logger.debug(msg)
    fp = nc4.Dataset(str(filepath))
    nctime = fp.variables['time']
    date_lst = nc4.num2date(nctime[:], nctime.units,
                            only_use_cftime_datetimes=False,
                            only_use_python_datetimes=True)
    #-- convert to pd.Timestamp
    date_lst = np.array([pd.Timestamp(_) for _ in date_lst])
    assert np.all(date_lst==np.sort(date_lst))
    nobs = len(date_lst)

    #
    #-- temporal selection
    #
    tcnd = None
    if date_first!=None and date_lst[-1]<date_first:
        return None
    elif date_first!=None:
        if tcnd is None:
            tcnd = date_lst>=date_first
        else:
            tcnd &= (date_lst>=date_first)
    if date_last!=None and date_lst[0]>date_last:
        return None
    elif date_last!=None:
        if tcnd is None:
            tcnd = date_lst<=date_last
        else:
            tcnd &= (date_lst<=date_last)
    if tcnd is None:
        i0 = 0
        i1 = len(date_lst)+1
    else:
        tidxs = np.where(tcnd)[0]
        if len(tidxs)==0:
            return None
        else:
            i0 = tidxs[0]
            i1 = tidxs[-1]+1
    date_lst = date_lst[i0:i1]

    nnobs = len(date_lst)
    # print(f"MVODEBUG:: nnobs={nnobs}, detected date range {date_lst[0]} --- {date_lst[-1]})")
    #
    #-- prepare output dictionary
    #
    ch4_dict = {}
    ch4_dict['time'] = date_lst
    #-- CH4 concentration
    ncvar = fp.variables['value']
    assert ncvar.units=='mol mol-1'
    ch4_data = ncvar[i0:i1]
    ch4_dict['units'] = 'ppb'
    ch4_dict['ch4'] = ch4_data*1e9
    #-- CH4 uncertainty
    if 'value_unc' in fp.variables:
        ncvar = fp.variables['value_unc']
        assert ncvar.units=='mol mol-1'
        ch4_dict['dataunc'] = ncvar[i0:i1]*1e9
    else:
        ch4_dict['dataunc'] = None
    if 'altitude' in fp.variables:
        ncvar = fp.variables['altitude']
        assert ncvar.units=='m'
        ch4_dict['altitude'] = ncvar[i0:i1]
    else:
        ch4_dict['altitude'] = None
    if 'elevation' in fp.variables:
        ncvar = fp.variables['elevation']
        assert ncvar.units=='m'
        ch4_dict['elevation'] = ncvar[i0:i1]
    else:
        ch4_dict['elevation'] = None
    if 'elevation' in fp.variables:
        ncvar = fp.variables['elevation']
        assert ncvar.units=='m'
        ch4_dict['elevation'] = ncvar[i0:i1]
    else:
        ch4_dict['elevation'] = None
    if 'intake_height' in fp.variables:
        ncvar = fp.variables['intake_height']
        assert ncvar.units=='m'
        ch4_dict['intake_height'] = ncvar[i0:i1]
    else:
        ch4_dict['intake_height'] = None
    fp.close()

    return ch4_dict


def load_obspack_ch4(filepath : Union[str,Path],
                     date_first : Union[pd.Timestamp,None] = None,
                     date_last  : Union[pd.Timestamp,None] = None):
    """Reading CH4 concentrations from NetCDF file provided by obspack.
    """
    if date_first!=None and date_last!=None and not date_first<date_last:
        msg = f"inconsistent temporal selection {date_first} --- {date_last}"
        raise RuntimeError(msg)

    # msg = f"loading obspack data from ***{filepath}***..."
    # logger.debug(msg)
    fp = nc4.Dataset(str(filepath))
    nctime = fp.variables['time']
    date_lst = nc4.num2date(nctime[:], nctime.units,
                            only_use_cftime_datetimes=False,
                            only_use_python_datetimes=True)
    #-- convert to pd.Timestamp
    date_lst = np.array([pd.Timestamp(_) for _ in date_lst])
    assert np.all(date_lst==np.sort(date_lst))
    nobs = len(date_lst)

    #
    #-- temporal selection
    #
    tcnd = None
    if date_first!=None and date_lst[-1]<date_first:
        return None
    elif date_first!=None:
        if tcnd is None:
            tcnd = date_lst>=date_first
        else:
            tcnd &= (date_lst>=date_first)
    if date_last!=None and date_lst[0]>date_last:
        return None
    elif date_last!=None:
        if tcnd is None:
            tcnd = date_lst<=date_last
        else:
            tcnd &= (date_lst<=date_last)
    if tcnd is None:
        i0 = 0
        i1 = len(date_lst)+1
    else:
        tidxs = np.where(tcnd)[0]
        if len(tidxs)==0:
            return None
        else:
            i0 = tidxs[0]
            i1 = tidxs[-1]+1
    date_lst = date_lst[i0:i1]

    nnobs = len(date_lst)
    # print(f"MVODEBUG:: nnobs={nnobs}, detected date range {date_lst[0]} --- {date_lst[-1]})")
    #
    #-- prepare output dictionary
    #
    ch4_dict = {}
    ch4_dict['time'] = date_lst
    #-- CH4 concentration
    ncvar = fp.variables['value']
    assert ncvar.units=='mol mol-1'
    ch4_data = ncvar[i0:i1]
    ch4_dict['units'] = 'ppb'
    ch4_dict['ch4'] = ch4_data*1e9
    #-- CH4 uncertainty
    if 'value_unc' in fp.variables:
        ncvar = fp.variables['value_unc']
        assert ncvar.units=='mol mol-1'
        ch4_dict['dataunc'] = ncvar[i0:i1]*1e9
    else:
        ch4_dict['dataunc'] = None
    if 'altitude' in fp.variables:
        ncvar = fp.variables['altitude']
        assert ncvar.units=='m'
        ch4_dict['altitude'] = ncvar[i0:i1]
    else:
        ch4_dict['altitude'] = None
    if 'elevation' in fp.variables:
        ncvar = fp.variables['elevation']
        assert ncvar.units=='m'
        ch4_dict['elevation'] = ncvar[i0:i1]
    else:
        ch4_dict['elevation'] = None
    if 'elevation' in fp.variables:
        ncvar = fp.variables['elevation']
        assert ncvar.units=='m'
        ch4_dict['elevation'] = ncvar[i0:i1]
    else:
        ch4_dict['elevation'] = None
    if 'intake_height' in fp.variables:
        ncvar = fp.variables['intake_height']
        assert ncvar.units=='m'
        ch4_dict['intake_height'] = ncvar[i0:i1]
    else:
        ch4_dict['intake_height'] = None
    fp.close()

    return ch4_dict


def subcmd_stations_visu(args):
    """Visualisation of TM5 simulated concentrations at stations.
    """
    expdir = Path(args.expdir)
    station = args.station
    obspackdir = args.obspackdir
    altdif_threshold = args.altdif_threshold
    #-
    figsize = args.figsize
    dpi     = args.dpi
    markersize = args.markersize

    if not expdir.is_dir():
        msg = f"experiment directory ==>{expdir}<== not existing."
        raise RuntimeError(msg)

    stafile = expdir / 'stations' / 'stations.nc4'
    if not stafile.exists():
        msg = f"station concentration file ==>{stafile}<== not found."
        raise RuntimeError(msg)

    #
    #-- comparison against CAMS concentrations (?)
    #
    with_cams = (args.camsfile!=None)
    #
    #-- comparison against obspack (?)
    #
    if obspackdir!=None:
        obspackdir = Path(obspackdir)
        if not obspackdir.is_dir():
            msg = f"provided obspack directory ***{args.obspackdir}*** not accessible"
            raise RuntimeError(msg)
    with_obspack = (obspackdir!=None)
    #
    #--
    #
    ds = nc4.Dataset(str(stafile))
    #
    #-- access global attributes
    #
    ntrac = ds.dimensions['tracers'].size
    date_lst = [pd.Timestamp(*_) for _ in ds['date_midpoints'][:]]
    tracers = [getattr(ds, f'tracer_{itrac+1:03d}') for itrac in range(ntrac)]
    tcover_start = pd.Timestamp(ds.getncattr('starting time'))
    tcover_end   = pd.Timestamp(ds.getncattr('ending time'))
    msg = f"detected temporal coverage {tcover_start} -- {tcover_end}"
    logger.info(msg)
    tcover_tag = f"{tcover_start.strftime('%Y%m%d')}-{tcover_end.strftime('%Y%m%d')}"
    #
    #-- each station time-serious is in separate NetCDF group
    #
    ngrp = len(ds.groups)
    stations = sorted(set([ds[_].getncattr('name') for _ in ds.groups]))
    nsta = len(stations)
    msg = f"detected {ngrp} NetCDF groups with {nsta} different station names"
    logger.info(msg)
    if station!=None:
        sta_lst = station
    else:
        sta_lst = stations
    # #
    # #--
    # #
    # gns_stations = set()
    # for _ in ds.groups:
    #     stagrp = ds[_]
    #     abbr_tag = stagrp.abbr.replace('FM/','')
    #     stareg = stagrp.region
    #     if stareg=='gns100x100':
    #         gns_stations.add(abbr_tag)
    # gns_stations = sorted(list(gns_stations))
    # with open('gns_stations.txt', 'w') as fp:
    #     for _ in gns_stations:
    #         fp.write(f"{_}" + '\n')
    # sys.exit(0)
    #
    #-- loop over stations (station names)
    #
    for sta in sta_lst:
        #
        #-- detect entries starting with station name
        #
        station_ids = [_ for _ in ds.groups if ds[_].getncattr('name').lower().startswith(sta.lower())]
        nstaid = len(station_ids)
        if nstaid==0:
            msg = f"no entries detected for station -->{sta}<--"
            logger.warning(msg)
            continue
        msg = f"detected {nstaid} entries for @station -->{sta}<--"
        logger.debug(msg)
        for staid in station_ids:
            stagrp = ds[staid]
            abbr_tag = stagrp.abbr.replace('FM/','')
            sta_alt = stagrp.altitude
            stalon = stagrp.longitude
            stalat = stagrp.latitude
            stareg = stagrp.region
            ncmix = stagrp['mixing_ratio']
            assert ncmix.dimensions==('tracers','samples')
            mix_unit = ncmix.unit
            assert mix_unit.startswith('mole fraction (')
            assert mix_unit.endswith(')')
            mix_unit = mix_unit[:-1].replace('mole fraction (','')
            #-- each tracer has (potentially) it' own unit
            mix_unit = [_.strip() for _ in mix_unit.split(',')]
            #-- but here make sure all have same unit
            assert np.all(np.array(mix_unit)==mix_unit[0])
            mix_unit = mix_unit[0]
            mix_ratio = ncmix[:]
            assert ma.count_masked(mix_ratio)==0
            mix_ratio = mix_ratio.data
            data_dict = {'time': date_lst}
            #
            #-- insert mixing ratio(s)
            #
            if ntrac==1:
                tracer_tag = f"fit-ic {tracers[0].lower()}"
                data_dict[tracer_tag] = mix_ratio[0,:]
            elif args.tracersum:
                tracer_tag = f"fit-ic tracer-sum"
                data_dict[tracer_tag] = mix_ratio.sum(axis=0)
            else:
                for itrac,tracer in enumerate(tracers):
                    tracer_tag = f"fit-ic {tracer.lower()}"
                    data_dict[tracer_tag] = mix_ratio[itrac,:]

            #
            #-- build data frame(s) (per current station and station level)
            #
            df = pd.DataFrame.from_dict(data_dict)
            if args.hour!=None:
                msg = f"...restrict to simulations in hour={args.hour}"
                logger.info(msg)
                df = df[df.time.dt.hour==args.hour]
            # print(f"FIT-IC simulated")
            # print(df.head())
            # print(f"-"*30)
            # dff = df[df.time.dt.hour==15]
            # print(f"FIT-IC@hour=15")
            # print(dff.head())
            #
            #-- comparison against CAMS concentrations
            #
            cams_df = None
            if with_cams:
                cams_df = \
                    cams_at_obspack_load_conctseries(args.camsfile,
                                                     stalon,
                                                     stalat,
                                                     sta_alt)
                # print(cams_df.head())
                if args.hour!=None:
                    cams_df = cams_df[cams_df.time.dt.hour==args.hour]
                # print(cams_df.head())
                # sys.exit(0)
            #
            #-- comparison versus obspack (?)
            #
            dfobspack = None
            obspack_info = None
            if obspackdir!=None:
                obspack_info = obspack_load_conctseries(obspackdir,
                                                        tcover_start,
                                                        tcover_end,
                                                        abbr_tag,
                                                        stalon,
                                                        stalat,
                                                        sta_alt,
                                                        altdif_threshold)
                if obspack_info==None:
                    msg = f"...no matching obspack data found"
                    logger.warning(msg)
                else:
                    msg = f"obspack_info loaded, obsalt={obspack_info['altitude']}"
                    logger.debug(msg)
                    dfobspack = obspack_info['df']
                    if args.hour!=None:
                        msg = f"...restrict obspack to/in hour={args.hour}"
                        logger.info(msg)
                        dfobspack = dfobspack[dfobspack.time.dt.hour==args.hour]
                    # print(f"obspack")
                    # print(dfobspack.head())
            #--
            #
            #-- export to csv
            #
            if args.csv_output:
                dfcsv = df.copy()
                dfcsv.index = dfcsv['time']
                dfcsv = dfcsv.drop(['time',], axis=1)
                if len(dfcsv.columns)>0:
                    tracer_tag = 'multiple-tracers'
                else:
                    tracer_tag = dfcsv.columns[0].lower()
                outname_tokens = [abbr_tag, 'tm5-simu', tracer_tag, tcover_tag,]
                if args.hour!=None:
                    outname_tokens.append(f'hour{args.hour:02d}')
                outname = '_'.join(outname_tokens) + '.csv'
                outname = set_outname(args, outname)
                with open(outname, 'w') as fp:
                    fp.write(f"## input_file: {str(stafile.absolute())}" + '\n')
                    fp.write(f"## station: {sta}" + '\n')
                    fp.write(f"## longitude: {stalon}" + '\n')
                    fp.write(f"## latitude:  {stalat}" + '\n')
                    fp.write(f"## altitude:  {sta_alt}" + '\n')
                    fp.write(f"## tm5_region: {stareg}" + '\n')
                    df.to_csv(fp, index=False)
                    msg = f"generated file ***{outname}***"
                    logger.info(msg)
                #
                #--
                #
                if obspack_info!=None:
                    obspack_file = obspack_info['filepath']
                    dfobs = obspack_info['df']
                    obspack_alt = obspack_info['altitude']
                    outname_tokens = [abbr_tag, 'obspack', tcover_tag,]
                    if args.hour!=None:
                        outname_tokens.append(f'hour{args.hour:02d}')
                    outname = '_'.join(outname_tokens) + '.csv'
                    outname = set_outname(args, outname)
                    with open(outname, 'w') as fp:
                        fp.write(f"## input_file: {str(stafile.absolute())}" + '\n')
                        fp.write(f"## station: {sta}" + '\n')
                        fp.write(f"## longitude: {stalon}" + '\n')
                        fp.write(f"## latitude:  {stalat}" + '\n')
                        fp.write(f"## altitude_mean:  {obspack_alt}" + '\n')
                        fp.write(f"## obspack_filepath: {obspack_file}" + '\n')
                        dfobs.to_csv(fp, index=False)
                        msg = f"generated file ***{outname}***"
                        logger.info(msg)
            #
            #-- currently one plot per tracer
            #
            for curtrac in df.columns:
                if curtrac=='time':
                    continue
                dfplot = df[['time',curtrac]]
                tracer_tag = f"tm5-simu-{curtrac.lower()}"
                if with_cams and len(cams_df)>0:
                    tracer_tag += '-vs-cams'
                if obspackdir!=None:
                    tracer_tag += '-vs-obspack'
                outname_tokens = [abbr_tag, tracer_tag, tcover_tag,]
                if args.hour!=None:
                    outname_tokens.append(f'hour{args.hour:02d}')
                outname = '_'.join(outname_tokens) + '.csv'
                outname = set_outname(args, outname)
                fig = plt.figure(figsize=figsize, dpi=dpi)#, tight_layout=True)
                ax = fig.add_subplot(111)
                #
                #-- plot FIT-IC simulation
                #
                #-- 2025-07-30: switchec back to only line for simulation
                dfplot.plot(ax=ax, x='time', grid=True,
                            ylabel=','.join(tracers), xlabel='time',
                            color='blue')
                # if args.style=='line':
                #     dfplot.plot(ax=ax, x='time', grid=True,
                #                 ylabel=','.join(tracers), xlabel='time')
                # elif args.style=='marker':
                #     dfplot.plot(ax=ax, x='time', grid=True,
                #                 ls='', marker='.', markersize=markersize,
                #                 color='blue',
                #                 ylabel=','.join(tracers), xlabel='time')
                #
                #-- plot obspack concentrations
                #
                if not dfobspack is None:
                    if  args.style=='line':
                        dfobspack.plot(ax=ax, x='time', grid='true', alpha=0.5)
                    elif args.style=='marker':
                        dfobspack.plot(ax=ax, x='time', grid='true',
                                       ls='', marker='x', markersize=markersize,
                                       color='orange')
                #
                #-- plot CAMS concentrations
                #
                if not cams_df is None and len(cams_df)>0:
                    #-- 2025-07-30: switched back to only line plot
                    cams_df.plot(ax=ax, x='time', grid=True,
                                 color='green', alpha=0.5, )
                    # if args.style=='line':
                    #     cams_df.plot(ax=ax, x='time', grid=True, alpha=0.5)
                    # elif args.style=='marker':
                    #     cams_df.plot(ax=ax, x='time', grid=True,
                    #                  ls='', marker='+', markersize=markersize,
                    #                  color='green')
                #
                #-- ylimits
                #
                if args.ylimits!=None:
                    ax.set_ylim(args.ylimits)
                #
                #-- title
                #
                title = f"CH4@{sta} (alt: {sta_alt}[m], {stareg})"
                if obspack_info!=None:
                    obspack_alt = obspack_info['altitude']
                    if sta_alt!=obspack_alt:
                        title += f" obspack_alt={obspack_alt}[m]"
                ax.set_title(title)
                ax.set_xlabel('time')
                ax.set_ylabel(f"conc [{mix_unit}]")
                outname = '_'.join(outname_tokens) + '.png'
                outname = set_outname(args, outname)
                plt.savefig(outname, dpi=dpi)
                plt.close()
                msg = f"generated file ***{outname}***"
                logger.info(msg)


def subcmd_stations_cmpvisu(args):
    """Visualisation of TM5 simulated concentrations at stations
    from two simulations.
    """
    expdir_lst = [Path(_) for _ in args.expdir]
    exptag_lst = args.exptag
    station = args.station
    obspackdir = args.obspackdir
    figsize = args.figsize
    dpi     = args.dpi

    tag1,tag2 = exptag_lst
    exp_map = OrderedDict()
    for exp in zip(exptag_lst, expdir_lst):
        _tag, expdir = exp

        if not expdir.is_dir():
            msg = f"experiment directory ==>{expdir}<== not existing."
            raise RuntimeError(msg)

        stafile = expdir / 'stations' / 'stations.nc4'
        if not stafile.exists():
            msg = f"station concentration file ==>{stafile}<== not found."
            raise RuntimeError(msg)
        exp_map[_tag] = {'stafile':stafile}
        #
        #--
        #
        ds = nc4.Dataset(str(stafile))
        exp_map[_tag]['ds'] = ds
        #
        #-- access global attributes
        #
        ntrac = ds.dimensions['tracers'].size
        date_lst = [pd.Timestamp(*_) for _ in ds['date_midpoints'][:]]
        tracers = [getattr(ds, f'tracer_{itrac+1:03d}') for itrac in range(ntrac)]
        tcover_start = pd.Timestamp(ds.getncattr('starting time'))
        tcover_end   = pd.Timestamp(ds.getncattr('ending time'))
        msg = f"detected temporal coverage {tcover_start} -- {tcover_end}"
        logger.info(msg)
        tcover_tag = f"{tcover_start.strftime('%Y%m%d')}-{tcover_end.strftime('%Y%m%d')}"
        #
        #-- each station time-serious is in separate NetCDF group
        #
        ngrp = len(ds.groups)
        stations = sorted(set([ds[_].getncattr('name') for _ in ds.groups]))
        #
        #--
        #
        exp_map[_tag]['ngrp'] = ngrp
        exp_map[_tag]['date_lst'] = date_lst
        exp_map[_tag]['tracers'] = tracers
        exp_map[_tag]['tcover_tag'] = tcover_tag
        exp_map[_tag]['stations'] = stations
    #
    assert exp_map[tag1]['ngrp']==exp_map[tag2]['ngrp']
    assert exp_map[tag1]['stations']==exp_map[tag2]['stations']
    assert exp_map[tag1]['tracers']==exp_map[tag2]['tracers']
    stations = exp_map[tag1]['stations']
    nsta = len(stations)
    msg = f"detected {ngrp} NetCDF groups with {nsta} different station names"
    logger.info(msg)
    tracers  = exp_map[tag1]['tracers']
    ntrac = len(tracers)
    msg = f"detected {ntrac} tracers in file -->{tracers}<--"
    logger.info(msg)
    #
    #-- comparison against CAMS concentrations (?)
    #
    with_cams = (args.camsfile!=None)
    #
    #-- comparison against obspack (?)
    #
    if obspackdir!=None:
        obspackdir = Path(obspackdir)
        if not obspackdir.is_dir():
            msg = f"provided obspack directory ***{args.obspackdir}*** not accessible"
            raise RuntimeError(msg)
    with_obspack = (obspackdir!=None)

    if station!=None:
        sta_lst = station
    else:
        sta_lst = stations
    # #
    # #--
    # #
    # gns_stations = set()
    # for _ in ds.groups:
    #     stagrp = ds[_]
    #     abbr_tag = stagrp.abbr.replace('FM/','')
    #     stareg = stagrp.region
    #     if stareg=='gns100x100':
    #         gns_stations.add(abbr_tag)
    # gns_stations = sorted(list(gns_stations))
    # with open('gns_stations.txt', 'w') as fp:
    #     for _ in gns_stations:
    #         fp.write(f"{_}" + '\n')
    # sys.exit(0)
    #
    #-- loop over stations (station names)
    #
    ds1 = exp_map[tag1]['ds']
    ds2 = exp_map[tag2]['ds']
    for sta in sta_lst:
        #
        #-- detect entries starting with station name
        #   -> can take from first experiment directory
        station_ids = [_ for _ in ds1.groups if ds1[_].getncattr('name').lower().startswith(sta.lower())]
        nstaid = len(station_ids)
        if nstaid==0:
            msg = f"no entries detected for station -->{sta}<--"
            logger.warning(msg)
            continue
        msg = f"detected {nstaid} entries for @station -->{sta}<--"
        logger.debug(msg)
        for staid in station_ids:
            stagrp1 = ds1[staid]
            stagrp2 = ds2[staid]
            abbr_tag = stagrp1.abbr.replace('FM/','')
            sta_alt = stagrp1.altitude
            stalon = stagrp1.longitude
            stalat = stagrp1.latitude
            stareg = stagrp1.region
            ncmix1 = stagrp1['mixing_ratio']
            ncmix2 = stagrp2['mixing_ratio']
            #-- NOTE: no checking of second experiment outputs yet
            assert ncmix1.dimensions==('tracers','samples')
            mix_unit = ncmix1.unit
            assert mix_unit.startswith('mole fraction (')
            assert mix_unit.endswith(')')
            mix_unit = mix_unit[:-1].replace('mole fraction (','')
            #-- each tracer has (potentially) it' own unit
            mix_unit = [_.strip() for _ in mix_unit.split(',')]
            #-- but here make sure all have same unit
            assert np.all(np.array(mix_unit)==mix_unit[0])
            mix_unit = mix_unit[0]
            mix_ratio1 = ncmix1[:]
            assert ma.count_masked(mix_ratio1)==0
            mix_ratio1 = mix_ratio1.data
            mix_ratio2 = ncmix2[:]
            assert ma.count_masked(mix_ratio2)==0
            mix_ratio2 = mix_ratio2.data
            #
            #-- build data dictionary
            #
            data_dict = {'time': date_lst}
            #
            #-- insert mixing ratio(s)
            #
            if ntrac==1:
                tracer_tag = f"{tag1} {tracers[0].lower()}"
                data_dict[tracer_tag] = mix_ratio1[0,:]
                tracer_tag = f"{tag2} {tracers[0].lower()}"
                data_dict[tracer_tag] = mix_ratio2[0,:]
                tracer_ftag = f"{tracers[0].lower()}_{tag1}-vs{tag2}"
            elif args.tracersum:
                tracer_tag = f"{tag1} tracer-sum"
                data_dict[tracer_tag] = mix_ratio1.sum(axis=0)
                tracer_tag = f"{tag2} tracer-sum"
                data_dict[tracer_tag] = mix_ratio2.sum(axis=0)
                tracer_ftag = f"sum-of-tracers_{tag1}-vs{tag2}"
            else:
                tracer_ftag = f"multiple-tracers_{tag1}-vs{tag2}"
                for itrac,tracer in enumerate(tracers):
                    tracer_tag = f"{tag1} {tracer.lower()}"
                    data_dict[tracer_tag] = mix_ratio1[itrac,:]
                    tracer_tag = f"{tag2} {tracer.lower()}"
                    data_dict[tracer_tag] = mix_ratio2[itrac,:]

            #
            #-- build data frame(s) (per current station and station level)
            #
            df = pd.DataFrame.from_dict(data_dict)
            cams_df = None
            if with_cams:
                cams_df = \
                    cams_at_obspack_load_conctseries(args.camsfile,
                                                     stalon,
                                                     stalat,
                                                     sta_alt)
            #-- obspack(?)
            dfobspack = None
            if obspackdir!=None:
                dfobspack = obspack_load_conctseries(obspackdir,
                                                     tcover_start,
                                                     tcover_end,
                                                     abbr_tag,
                                                     stalon,
                                                     stalat,
                                                     sta_alt)
            #--
            #
            #-- export to csv
            #
            if args.csv_output:
                dfcsv = df.copy()
                dfcsv.index = dfcsv['time']
                dfcsv = dfcsv.drop(['time',], axis=1)
                outname_tokens = [abbr_tag, '', tracer_ftag, tcover_tag,]
                outname = '_'.join(outname_tokens) + '.csv'
                outname = set_outname(args, outname)
                with open(outname, 'w') as fp:
                    fp.write(f"## input_file: {str(stafile.absolute())}" + '\n')
                    fp.write(f"## station: {sta}" + '\n')
                    fp.write(f"## longitude: {stalon}" + '\n')
                    fp.write(f"## latitude:  {stalat}" + '\n')
                    fp.write(f"## altitude:  {sta_alt}" + '\n')
                    fp.write(f"## tm5_region: {stareg}" + '\n')
                    df.to_csv(fp, index=False)
                    msg = f"generated file ***{outname}***"
                    logger.info(msg)
            #
            #-- currently one plot per tracer
            #
            tracer_plot_lst = [ _.replace(f'{tag1} ','').replace(f'{tag2} ','') for _ in df.columns if _!='time' ]
            for curtrac in tracer_plot_lst:
                _c1 = f"{tag1} {curtrac}"
                _c2 = f"{tag2} {curtrac}"
                dfplot = df[['time',_c1, _c2]]
                tracer_tag = f"{curtrac.lower()}-{tag1}-vs-{tag2}"
                if with_cams:
                    tracer_tag += '-vs-cams'
                if obspackdir!=None:
                    tracer_tag += '-vs-obspack'
                outname_tokens = [abbr_tag, tracer_tag, tcover_tag,]
                outname = '_'.join(outname_tokens) + '.csv'
                outname = set_outname(args, outname)
                fig = plt.figure(figsize=figsize, dpi=dpi)#, tight_layout=True)
                ax = fig.add_subplot(111)
                dfplot.plot(ax=ax, x='time', grid=True)# ,
                            # ylabel=','.join(tracers), xlabel='time')
                if not cams_df is None:
                    cams_df.plot(ax=ax, x='time', grid=True, alpha=0.5)
                if not dfobspack is None:
                    dfobspack.plot(ax=ax, x='time', grid='true', alpha=0.5)
                ax.set_title(f"CH4@{sta} (alt: {sta_alt}[m], {stareg})")
                ax.set_xlabel('time')
                ax.set_ylabel(f"conc [{mix_unit}]")
                outname = '_'.join(outname_tokens) + '.png'
                outname = set_outname(args, outname)
                plt.savefig(outname, dpi=dpi)
                plt.close()
                msg = f"generated file ***{outname}***"
                logger.info(msg)


def subcmd_mixdir_inspect(args):
    """
    """
    expdir = Path(args.expdir)

    mixdir = expdir / 'mix'

    if not mixdir.is_dir():
        msg = f"mixing files top-level directory ==>{mixdir}<== not existing."
        raise RuntimeError(msg)

    year_dirlst = sorted(mixdir.glob('????'))
    print(f"year_dirlst ==>{year_dirlst}<==")
    assert len(year_dirlst)==1
    yeardir = year_dirlst[0]
    mon_dirlst = sorted(yeardir.glob('??'))
    ifile = 0
    for mdir in mon_dirlst:
        print(f"@{mdir} now...")
        mfile_lst = sorted(mdir.glob('mix_????????.nc4'))
        for _mfile in mfile_lst:
            ds = nc4.Dataset(str(_mfile))
            for grp in ds.groups:
                ncgrp = ds[grp]
                cur_mix = ncgrp['mix'][:]
                msg = f"...@{_mfile.name},{grp}: mix min/mean/max = " \
                    f"{cur_mix.min()}/{cur_mix.mean()}/{cur_mix.max()}"
                logger.info(msg)
            ds.close()
            ifile += 1

def subcmd_mix_inspect(args):
    """
    """
    filepath = Path(args.filepath)
    if not filepath.exists():
        msg = f"mixing file ***{str(filepath)}*** not found."
        raise RuntimeError(msg)
    ds = nc4.Dataset(filepath)
    #
    #-- global dimensions
    #
    ntime = ds.dimensions['times'].size
    ntrac = ds.dimensions['tracers'].size
    nlev  = ds.dimensions['levels'].size
    #-- at least for now: assume daily output
    assert ntime==1
    #
    #-- regions are stored in groups
    for grp in ds.groups:
        ncgrp = ds.groups[grp]
        nlon = ncgrp.dimensions['longitude'].size
        nlat = ncgrp.dimensions['latitude'].size
        #
        #--
        #
        date_lst = [pd.Timestamp(*_) for _ in ncgrp['sample_times'][:]]
        print(date_lst) 
        ncmix = ncgrp['mix']
        assert ncmix.dimensions==('tracers','times','levels','latitude','longitude')
        assert ncmix.shape==(ntrac,ntime,nlev,nlat,nlon)
        ncgph = ncgrp['gph']
        assert ncgph.dimensions==('times','boundaries','latitude','longitude')
        for ilev in range(nlev):
            mix = ncmix[:,0,ilev,:,:]
            msg = f"@{grp}, ilev={ilev}, mixing ratio min/mean/max = " \
                f"{mix.min()}/{mix.mean()}/{mix.max()}"
            print(msg)
    ds.close()


def subcmd_emis_visu(args):
    """
    """
    filepath = args.filepath
    figsize = args.figsize
    dpi     = args.dpi
    ds = xr.open_dataset(filepath)

    dschk = ds.sel(lat=53.5,lon=10, method='nearest')
    for _v in ds.data_vars:
        print(dschk[_v].data)
    sys.exit(0)
    fig, axis = plt.subplots(1, 1,
                             subplot_kw=dict(projection=ccrs.PlateCarree()))

    ds.cams_total.plot(
        ax=axis,
        transform=ccrs.PlateCarree(),  # this is important!
        # usual xarray stuff
        cbar_kwargs={"orientation": "horizontal", "shrink": 0.7},
        robust=True,
    )
    axis.coastlines()  # cartopy function

    outname_tokens = [Path(filepath).stem, 'visu']
    outname = '_'.join(outname_tokens) + '.png'
    outname = set_outname(args, outname)
    plt.savefig(outname, dpi=dpi)
    plt.close()
    msg = f"generated file ***{outname}***"
    print(msg)

def parser():

    def _add_io_options(aparser):
        aparser.add_argument( '--outname',
                              help="""write output to this file""" )
        aparser.add_argument( '--outdir',
                              help="""where to place generated file(s).""" )
    def _add_plot_options(aparser):
        aparser.add_argument('--figsize',
                             type=float,
                             nargs=2,
                             default=(10,6),
                             help="""figure size, width/height[inch] (default: %(default)s).""")
        aparser.add_argument('--dpi',
                             type=int,
                             default=150,
                             help="""dots-per-inch (default: %(default)s).""")

    #
    #--
    #
    parser = ArgumentParser(usage=globals()['__doc__'])

    #----------------------------
    #     s u b c o m m a n d s
    #
    subparsers = parser.add_subparsers( title='Available Subcommands',
                                        metavar='CMDS',
                                        description='',
                                        dest='subcmds',
                                        help='')
    #
    #--       stations_visu
    #
    xparser = subparsers.add_parser('stations_visu',
                                    help="""showing concentration time-series at stations.""")
    xparser.add_argument('expdir', help="""TM5 simulation directory.""")
    xparser.add_argument('--station',
                         nargs='+',
                         help="""select stations by name.""")
    xparser.add_argument('--tracersum',
                         action='store_true',
                         help="""whether to show sum of simulated tracers (in case of multiple-tracer run.""")
    xparser.add_argument('--camsfile',
                         help="""provision of NetCDF file with CAMS concentrations at obspack sites triggers comparison plot.""")
    xparser.add_argument('--obspackdir',
                         help="""comparison against obspack CH4 measurements will be added (if matching obspack file for location is found).""")
    xparser.add_argument('--altdif_threshold',
                         type=float,
                         default=2,
                         help="""maximal allowed difference in altitude when picking obspack concentrations at station altitude.""")
    xparser.add_argument('--csv_output',
                         action='store_true',
                         help="""whether to write station time-series also to csv file.""")
    xparser.add_argument('--hour',
                         type=int,
                         help="""restrict to simulations/observations at/in this hour.""")
    xparser.add_argument('--ylimits',
                         type=float,
                         nargs=2,
                         help="""explicitly set limits for the y-axis.""")
    xparser.add_argument('--style',
                         choices=['line','marker'],
                         default='line',
                         help="""style used for plotting the time-series (default: %(default)s).""")
    xparser.add_argument('--markersize',
                         type=float,
                         default=2,
                         help="""marker size passed to matplotlib (default: %(default)s).""")
    _add_plot_options(xparser)
    _add_io_options(xparser)

    #
    #--       stations_cmpvisu
    #
    xparser = subparsers.add_parser('stations_cmpvisu',
                                    help="""showing concentration time-series at stations.""")
    xparser.add_argument('expdir', nargs=2,
                         help="""directory from two similar TM5 simulations""")
    xparser.add_argument('--exptag', nargs=2,
                         default=['exp1','exp2'],
                         help="""descriptive tags for experiments (default: %(default)s).""")
    xparser.add_argument('--station',
                         nargs='+',
                         help="""select stations by name.""")
    xparser.add_argument('--tracersum',
                         action='store_true',
                         help="""whether to show sum of simulated tracers (in case of multiple-tracer run.""")
    xparser.add_argument('--camsfile',
                         help="""provision of NetCDF file with CAMS concentrations at obspack sites triggers comparison plot.""")
    xparser.add_argument('--obspackdir',
                         help="""comparison against obspack CH4 measurements will be added (if matching obspack file for location is found).""")
    xparser.add_argument('--csv_output',
                         action='store_true',
                         help="""whether to write station time-series also to csv file.""")
    _add_plot_options(xparser)
    _add_io_options(xparser)

    #
    #--
    #
    xparser = subparsers.add_parser('mixdir_inspect',
                                    help="""check TM5 generated mixing ratio files.""")
    xparser.add_argument('expdir', help="""TM5 simulation directory.""")
    _add_io_options(xparser)

    #
    #--
    #
    xparser = subparsers.add_parser('mix_inspect',
                                    help="""check TM5 generated mixing ratio single file.""")
    xparser.add_argument('filepath', help="""TM5 simulation generated mix file.""")
    _add_io_options(xparser)

    #
    #--
    #
    xparser = subparsers.add_parser('emis_visu',
                                    help="""visualisation of single emissions file.""")
    xparser.add_argument('filepath',
                         help="""NetCDF file generated by TM5 pre-processor.""")
    _add_plot_options(xparser)
    _add_io_options(xparser)
    
    return parser


def main(args):

    if args.subcmds=='stations_visu':
        subcmd_stations_visu(args)

    if args.subcmds=='stations_cmpvisu':
        subcmd_stations_cmpvisu(args)

    if args.subcmds=='mixdir_inspect':
        subcmd_mixdir_inspect(args)

    if args.subcmds=='mix_inspect':
        subcmd_mix_inspect(args)

    if args.subcmds=='emis_visu':
        subcmd_emis_visu(args)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#                    M A I N
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':
    import datetime as dtm

    progname = os.path.basename(__file__)

    #-----------------------------
    #          P R O G R A M   S T A R T
    #
    fmt = "%Y-%m-%dT%H:%M:%S.%f"
    ttstart = dtm.datetime.now()
    logger.info(f"{progname}::PROGRAM START::{ttstart.strftime(fmt)}")
    argv = ' '.join(sys.argv)
    logger.info(f"  command-line -->{argv}<--")

    #
    #          p a r s e   c o m m a n d   l i n e
    #
    parser = parser()
    args = parser.parse_args()

    #
    #--        s t a r t   e x e c u t i o n
    #
    main(args)
