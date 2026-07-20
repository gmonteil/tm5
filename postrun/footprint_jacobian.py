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
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots,colorbar
from cartopy import crs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from tm5.fitic import read_obs_table
from tm5.gridtools import TM5Grids
from tm5.observations import read_obspack_file
from tm5.post.footprint_io import load_adjoint_fwd #-- this was for earlier diagnostics
from tm5.post.footprint_io import region_table, _init_region_table
from tm5.post.footprint_io import tm5rundir_obstable, tm5rundir_iniconc_1obs
from tm5.post.footprint_io import regions1D_info
from tm5.post.footprint_io import tm5rundir_jacobian3D, tm5rundir_simustart
from tm5.post.footprint_io import tm5rundir_jacobian3D_old, tm5_fitic_adjoint_corrected_halos
from tm5.post.footprint_io import tm5rundir_emissions2D, tm5emisdir_load_emissions2D
from tm5.post.footprint_io import jacobian_redistribute_glb6x4_to_avengers_zoom
from tm5.post.plot_util import cnorm_set
from tm5.post.utilities import lonstr,latstr,set_outname,create_sha512

#
#-- mapping of station identifiers to obspack NetCDF file names
#   (for identifiers where the mapping is not "generic")
#
obspack_lookup_table = OrderedDict()
obspack_lookup_table['brm_45']   = 'ch4_brm_tower-insitu_49_allvalid-44magl.nc'
obspack_lookup_table['cbw_35']   = 'ch4_cbw_surface-insitu_118_allvalid.nc'
obspack_lookup_table['cmn_7']    = 'ch4_cmn_surface-insitu_106_allvalid.nc'
obspack_lookup_table['cmn_8']    = 'ch4_cmn_surface-insitu_443_allvalid.nc'
obspack_lookup_table['ers_40']   = 'ch4_ers_surface-insitu_11_allvalid.nc'
obspack_lookup_table['hei_30']   = 'ch4_hei_surface-insitu_22_allvalid.nc'
obspack_lookup_table['hel_110']  = 'ch4_hel_surface-insitu_147_allvalid.nc'
obspack_lookup_table['hpb_5']    = 'ch4_hpb_surface-flask_1_representative.nc'
obspack_lookup_table['hun_96']   = 'ch4_hun_surface-flask_1_representative.nc'
obspack_lookup_table['jfj_0']    = 'ch4_jfj_surface-flask_45_representative.nc'
obspack_lookup_table['jfj_10']   = 'ch4_jfj_surface-insitu_5_allvalid.nc'
obspack_lookup_table['jfj_14']   = 'ch4_jfj_tower-insitu_49_allvalid-13magl.nc'
obspack_lookup_table['lhw_32']   = 'ch4_lhw_surface-insitu_5_allvalid.nc'
obspack_lookup_table['lin_2']    = 'ch4_lin_tower-insitu_147_allvalid-2.5magl.nc'
obspack_lookup_table['lut_60']   = 'ch4_lut_surface-insitu_44_allvalid.nc'
obspack_lookup_table['oxk_150']  = 'ch4_oxk_surface-flask_1_representative.nc'
obspack_lookup_table['pdm_28']   = 'ch4_pdm_surface-insitu_11_allvalid.nc'
obspack_lookup_table['prs_10']   = 'ch4_prs_surface-insitu_21_allvalid.nc'
obspack_lookup_table['puy_10']   = 'ch4_puy_surface-insitu_11_allvalid.nc'
obspack_lookup_table['tac_180']  = 'ch4_tac_surface-flask_1_representative.nc'
obspack_lookup_table['wao_10']   = 'ch4_wao_surface-insitu_13_allvalid.nc'
obspack_lookup_table['wes_14']   = 'ch4_wes_surface-insitu_25_allvalid.nc'
obspack_lookup_table['zsf_3']    = 'ch4_zsf_surface-insitu_25_allvalid.nc'


nsecday = 86400

def tm5refdir_load_stationconc( refdir : str | Path, obsid : str ) -> xr.DataArray:
    station_file = Path(refdir) / 'stations' / 'stations.nc4'
    if not station_file.exists():
        msg = f"stations result file ***{str(station_file)}*** not found on system."
        raise FileNotFoundError(msg)
    #
    #--
    #
    msg = f"...try loading {obsid} concentration from file {str(station_file)}"
    logger.debug(msg)
    ds = Dataset(station_file)
    assert ds.dimensions['tracers'].size==1, \
        f"multiple tracers in station file not yet supported " \
        f"({str(station_file)})"
    grp_fnd = None
    for g,gg in ds.groups.items():
        if gg.abbr.endswith(obsid):
            grp_fnd = gg
    #
    if grp_fnd is None:
        raise RuntimeError(f"-->{obsid}<-- not found in station file ***{str(station_file)}***")
    else:
        msg = f"obsid={obsid} found at group {grp_fnd.abbr}"
        logger.info(msg)
    #
    time_list = np.array([Timestamp(*_) for _ in ds['date_midpoints'][:]])
    ncvar = grp_fnd['/mixing_ratio']
    if not ncvar.dimensions!=('tracer','samples'):
        msg = f"mixing_ratio with unexpected dimensions -->{ncvar.dimensions}<--"
        raise RuntimeError(msg)
    conc = ncvar[0,:]
    df_conc = DataFrame.from_dict( {'time':time_list, 'conc': conc} )
    # conc_da = xr.DataArray(
    #     conc,
    #     dims = ('time',),
    #     coords = {
    #         'time': time_list
    #         }
    #     )
    # return conc_da
    return df_conc
  

def collect_input4inversion_obs1D( topdir : Path, domain_tag : str,
                                   obsday_range : DatetimeIndex,
                                   emisday_range : DatetimeIndex,
                                   remove_halo : bool = True,
                                   obsid : list|None = None,
                                   obsdir : Path|None = None,
                                   refdir : Path|None = None,
                                   outdir : Path|None = None) -> SimpleNamespace:

    """Routine to collect results from (a series of) TM5 adjoint footprint simulations
    that were run for observations within a dedicated domain.

    Expected directory naming in the toplevel directory is
    - footprints_gns100x100_yyyymmdd (for continuous measurements at stations in the innermost zoom)
    - footprints_glb600x400_yyyymmdd (for flask measurements)

    In the observations space the data will be concatenated to a 1D vector
    - day1
      - sta1, sta2, sta3,...
    - day2
      - sta1, sta2, ...

    Note, in particular for the flask measurements, the list of stations varies per day.
    """
    #-- arguments
    nobsday = len(obsday_range)
    obsdayf = obsday_range[0]
    obsdayl = obsday_range[-1]
    msg = f"...collecting footprint information for (daily) observations " \
        f"in temporal range {obsdayf.strftime('%Y%m%d')} - {obsdayl.strftime('%Y%m%d')}"
    logger.info(msg)

    if domain_tag=='gns100x100':
        region_list = ['glb600x400','eur300x200','gns100x100',]
    elif domain_tag=='glb600x400':
        region_list = ['glb600x400',]
    else:
        msg = f"unexpected domain selection -->{domain_tag}<--"
        raise RuntimeError(msg)
    
    #
    #-- loop over (obs) days
    #
    obsfile_info_cache = OrderedDict()
    stationid_1D = []
    obstime_1D = []
    obslon_1D = np.empty(0)
    obslat_1D = np.empty(0)
    obsalt_1D = np.empty(0)
    obsmix_1D = np.empty(0)
    inic_array1D = np.empty(0)
    jac_da    = None
    emis_lonc1D = None
    emis_latc1D = None
    emis_reg1D  = None
    rundir_list = []
    obs_extrainfo_dict = {} #-- to collect hour of observation and time-window length
    simu_start = None
    for iday,obsday in enumerate(obsday_range):
        #
        #-- naming pattern used by Guillaume: footprints_<domain_tag>_%Y%m%d
        #   -  where '%Y%m%d' refers to last day of simulation at 0:00 !
        #      E.g. to assemble inputs for February (Feb 1 to Feb 28)
        #      the footprint directory names are ranging from
        #      footprints_<domain_tag>_20210202, ..., footprints_<domain_tag>_20210301
        #
        dirday = obsday + Timedelta(days=1)
        dirpattern = dirday.strftime(f"footprints_{domain_tag}_*%Y%m%d")
        cur_rundir_list = sorted(topdir.glob(dirpattern))
        if len(cur_rundir_list)==0:
            msg = f"@obsday={obsday}, no TM5 run directory found " \
                f"(which may happen if there were not observations on that day)."
            logger.warning(msg)
            continue
        elif len(cur_rundir_list)>1:
            msg = f"@obsday={obsday}, expected single TM5 run directory only, but " \
                f"found ==>{cur_rundir_list}<=="
            raise RuntimeError(msg)
        rundir = cur_rundir_list[0]
        #
        #--
        #
        cur_simu_start = tm5rundir_simustart(rundir)
        if iday==0:
            simu_start = cur_simu_start
        else:
            if cur_simu_start!=simu_start:
                msg = f"expected simulation start {simu_start} but found " \
                    f"{cur_simu_start} (@***{rundir}***)"
                raise RuntimeError(msg)
        #
        #-- all footprint run directories should contain (a copy of) the
        #   yaml configuration file involved
        #
        yamlfile = rundir / 'tm5.yaml'
        if not yamlfile.exists():
            msg = f"...@{rundir}, ATTENTION yaml configuration file not present! " \
                f"Should usually not happen, here we will just skip!"
            logger.warning(msg)
            continue
        msg = f"@obsday={obsday.strftime('%Y-%m-%d')}, reading from directory -->{rundir}<--"
        logger.debug(msg)
        rundir_list.append(rundir)
        #
        #-- load observations
        #   MVO-NOTE: the footprint observation files for the zoom domain
        #             were prepared with all stations on every day
        #             (also if an actual observation was missing).
        #             But here we extract only those stations with
        #             observations on the current observational day.
        #
        obs_table = tm5rundir_obstable(rundir, drop_missing_value=True)
        if obsid==None:
            obsinfo_curday = obs_table.copy()
        else:
            #-- ATTENTION:problematic for flask sites (as these have the date appended)
            # cnd_sta = obs_table.index.isin(obsid)
            # obsinfo_curday = obs_table.loc[cnd_sta,:]
            if domain_tag=='gns100x100':
                cnd_sta = obs_table.index.isin(obsid)
                if np.count_nonzero(cnd_sta)==0:
                    msg = f"..@obsday={obsday.strftime('%Y-%m-%d')} no available observations for -->{obsid}<--"
                    logger.warning(msg)
                    continue
                #
                obsinfo_curday = obs_table.loc[obsid,:] #-- keep station ordering from command line
            elif domain_tag=='glb600x400':
                #-- CBW_60_20210131
                obsids_shortened = [ '_'.join(_.lower().split('_')[:2]) for _ in obs_table.index ]
                xobsids = []
                for _obsid in obsid:
                    if _obsid in obsids_shortened:
                        try:
                            idx = obsids_shortened.index(_obsid)
                            xobsids.append(obs_table.index[idx])
                        except KeyError:
                            continue
                obsinfo_curday = obs_table.loc[xobsids,:]
        #
        #-- drop duplicated 'obsid' in obstable index
        #   NOTE: e.g. on January 16, 2021 there were two entries
        #         for 'SPO_2847_2021011520' (with measurement times 20:43 and 20:56)
        #
        _n1 = len(obsinfo_curday.index)
        _n2 = len(list(set(obsinfo_curday.index)))
        if _n2<_n1:
            msg = f"@{rundir}, will drop duplicated observation identifiers in index!"
            logger.error(msg)
            logger.info(f"initially: -->{obsinfo_curday.index}<--")
            obsinfo_curday = obsinfo_curday.loc[~obsinfo_curday.index.duplicated(keep=False)]
            logger.info(f"after droping duplicates: -->{obsinfo_curday.index}<--")
        #
        #-- read observations
        #
        obslist_curday = []
        stalist_curday = []
        obslonlist_curday = []
        obslatlist_curday = []
        obsaltlist_curday = []
        for ista,_obsid in enumerate(obsinfo_curday.index):
            #
            #-- ATTENTION:
            #   - Guillaume has prepared flask obs table files
            #     such that the actual measurement is included and we can take if directly
            #     from this file.
            #   - !!!Contrary!!! the continuous obs table files that were used
            #     -->SO FAR<-- for the footprint simulations only contain dummy
            #     observations.
            #     In this case, getting the measurement value is much more cumbersome:
            #     - we have to go through the original observation files
            #       (the user must have specified the observation directory on invocation!),
            #       and extract the observation averaged over the station-specific time-window.
            #     
            #
            #
            if args.obsdir==None:
                _obsid_data = obsinfo_curday.loc[_obsid,:]
                stalist_curday.append(_obsid)
                obslist_curday.append(_obsid_data.mixing_ratio)
                obslonlist_curday.append(_obsid_data.lon)
                obslatlist_curday.append(_obsid_data.lat)
                obsaltlist_curday.append(_obsid_data.alt)
                obstime_1D.append(_obsid_data.time)
                #
                #-- record hour of observation and time-window
                #   (make sure these do not change over time)
                #
                if not _obsid in obs_extrainfo_dict:
                    obs_extrainfo_dict[_obsid] = {
                        'obshour': _obsid_data.time.hour,
                        'time_window_length': _obsid_data.time_window_length
                    }
                else:
                    assert obs_extrainfo_dict[_obsid]['obshour']==_obsid_data.time.hour
                    assert obs_extrainfo_dict[_obsid]['time_window_length']==_obsid_data.time_window_length
            else: #-- generally this should not be used but it is/was required for the early continuous footprint simulations where the point_input.nc4 files contained dummy and not real obspack observations.
                #
                #-- cumbersome case: extract averaged observation from original observation file
                #
                _obstime = obsinfo_curday.loc[_obsid,'time']
                obs_hr = _obstime.hour
                obs_tw = obsinfo_curday.loc[_obsid,'time_window_length'] #-- time-window length [s]
                if not _obsid in obsfile_info_cache:
                    #
                    #-- extract station identifier (e.g. from 'cbw_207' or xxx)
                    #
                    _obsid_tokens = _obsid.split('_')
                    if len(_obsid_tokens)==2:
                        staid,sta_alt = _obsid_tokens
                    elif len(_obsid_tokens)==3:
                        staid,sta_alt,obs_time = _obsid_tokens
                    else:
                        msg = f"unexpected obsid -->{_obsid}<--"
                        raise RuntimeError(msg)
                    if _obsid in obspack_lookup_table:
                        obsfile = Path(args.obsdir) / obspack_lookup_table[_obsid]
                    else:
                        #
                        #-- gns100x100: should be the continuous measurements
                        #
                        if domain_tag=='gns100x100':
                            ptn = f"ch4_{staid}*-{sta_alt}magl*.nc"
                        elif domain_tag=='glb600x400':
                            ptn = f"ch4_{staid.lower()}_surface-flask_*_representative.nc"
                        obspack_list = list(obsdir.glob(ptn))
                        if len(obspack_list)==0:
                            msg = f"no matching observation file found for staid -->{staid}<--"
                            raise RuntimeError(msg)
                        elif len(obspack_list)>1:
                            msg = f"staid -->{staid}<-- yields multiple observation files " \
                                f"==>{obspack_list}<==. This cannot be handled yet and station is " \
                                f"ignored."
                            logger.error(msg)
                            continue
                        else:
                            obsfile = obspack_list[0]
                    msg = f"...reading observed concentrations from file ***{str(obsfile)}***..."
                    logger.info(msg)
                    tstart = obsdayf 
                    tend   = obsdayl+Timedelta(seconds=86399)
                    #-- ATTENTION:: restricting to the day is unsufficient,
                    #               i.e. we have sites (at high altitude) where
                    #               1am is in the mid of the time-window which
                    #               thus exceeds to the day before.
                    #               we thus extend the period that is extracted
                    tstart -= Timedelta(days=1)
                    tend   += Timedelta(days=1)
                    obspack_info = read_obspack_file(obsfile, start=tstart, end=tend)
                    obs_df = obspack_info.data
                    obsfile_info_cache[_obsid] = obspack_info
                else:
                    obs_df = obsfile_info_cache[_obsid].data
                #
                #--
                #
                _ostart = obsday + Timedelta(hours=obs_hr) - Timedelta(seconds=obs_tw)
                _oend   = obsday + Timedelta(hours=obs_hr) + Timedelta(seconds=obs_tw)
                cnd_time = (obs_df['time']>=_ostart)&(obs_df['time']<=_oend)
                if np.count_nonzero(cnd_time)>0:
                    mix_day = obs_df.loc[cnd_time,'value'].mean() * 1.e9 #-- [mol/mol] to [ppb]
                    ## mix_day *= 1.e9 #-- [mol/mol] to [ppb]
                    # DEBUGGING (to check tiny differences in observed concentrations
                    #            when prepared with the prepare_obs4footp.py utility)
                    # msg = f"@{_obsid}/{obsday}: conc={mix_day}[ppb] (==>{obs_df.loc[cnd_time,'value'].values}<==)"
                    # logger.debug(msg)
                    stalist_curday.append(_obsid)
                    obslist_curday.append(mix_day)
                    obslonlist_curday.append(obs_df.loc[cnd_time,'longitude'].mean())
                    obslatlist_curday.append(obs_df.loc[cnd_time,'latitude'].mean())
                    obsaltlist_curday.append(obs_df.loc[cnd_time,'altitude'].mean())
                    obstime_1D.append(_obstime)
                else:
                    msg = f"...@{staid}, no observations found in time window {_ostart}==>{_oend}"
                    logger.debug(msg)
        #
        #-- restrict obsdata frame for current day to those stations without missing obs.
        #
        if len(stalist_curday)!=len(obsinfo_curday):
            cnd_sta = obsinfo_curday.index.isin(stalist_curday)
            obsinfo_curday = obsinfo_curday.loc[cnd_sta,:]
        #
        #-- observations done
        #
        nobs = len(obslist_curday)
        msg = f"...@obsday={obsday.strftime('%Y-%m-%d')}, nobs={nobs}"
        logger.info(msg)
        if nobs==0:
            continue
        #
        #-- extend 1D arrays
        #
        if domain_tag=='glb600x400':
            stationid_1D += [ '_'.join(_.lower().split('_')[:2]) for _ in obsinfo_curday.index ]
        elif domain_tag=='gns100x100':
            stationid_1D += list(obsinfo_curday.index)
        else:
            raise RuntimeError(f"unexpected domain -->{domain_tag}<--")
        obsmix_1D = np.concat((obsmix_1D,np.array(obslist_curday)))
        obslon_1D = np.concat((obslon_1D,np.array(obslonlist_curday)))
        obslat_1D = np.concat((obslat_1D,np.array(obslatlist_curday)))
        obsalt_1D = np.concat((obsalt_1D,np.array(obsaltlist_curday)))
        #
        #-- initial concentration
        #
        iniclist_curday = []
        for ista,staid in enumerate(obsinfo_curday.index):
            sta_obs_info = obsinfo_curday.loc[staid,:]
            # print(f"...@{staid}, sta_obs_info\n{sta_obs_info}\n")
            msg = f"...@{staid}, reading initial concentration."
            logger.info(msg)
            inic_info = tm5rundir_iniconc_1obs(rundir, sta_obs_info)
            point_out = inic_info.point_output
            # logger.info(f"-->{point_out}<--")
            if 'mixing_ratio' in point_out.index:
                inic = point_out.mixing_ratio
            else:
                #
                #-- different station identifiers being mapped to same numerical ID
                #   e.g. on January 29, 2021 the obs identifiers
                #   "BRW_27_2021012919", "BRW_16_2021012920" were both mapped to
                #   station_id=75.
                #
                _station_id  = np.unique(np.array(point_out.index))
                assert len(_station_id)==1, \
                    f"expected one single station_id but found ==>{_station_id}<=="
                _station_id = _station_id[0]
                #
                #-- extract entries with this station_id
                #
                _cnd = obsinfo_curday['station_id']==_station_id
                _obsinfo_cnd = obsinfo_curday.loc[_cnd,:]
                # logger.info(f"\n{_obsinfo_cnd[['time','lat','lon','mixing_ratio']]}")
                #
                #-- select index where observation identifier matches current station
                #
                inic = _obsinfo_cnd.loc[staid,'mixing_ratio']
            iniclist_curday.append(inic)
        cur_inic_array = np.array(iniclist_curday)
        inic_array1D = np.concat((inic_array1D,cur_inic_array))
        #
        #-- read Jacobian for
        #   - current obsday, but possibly restricted to selected stations
        #   - w.r.t. to emissions in the selected range
        #
        # jac_info = tm5rundir_jacobian3D_old(rundir, emisday_range=emisday_range,
        #                                 obsid=list(obsinfo_curday.index),
        #                                 remove_halo=remove_halo, clip_child=False)
        jac_info = tm5rundir_jacobian3D(rundir, emisday_range=emisday_range,
                                        obsid=list(obsinfo_curday.index),
                                        remove_halo=remove_halo, clip_child=False)
        #-- Jacobian shape: [nobs,nemisday,ng]
        if jac_da is None:
            jac_da = jac_info.jac3D
        else:
            jac_da = xr.concat([jac_da,jac_info.jac3D], dim='obs')
        #
        #-- emission location (required only once)
        #
        if emis_lonc1D is None:
            emis_lonc1D = jac_info.lonc1D
            emis_latc1D = jac_info.latc1D
            emis_reg1D  = jac_info.reg1D
    #
    #-- loop over obs days terminated
    #
    nobs = len(obsmix_1D)
    msg = f"...collected {nobs} observations overall."
    logger.info(msg)
    #
    #-- turn lists into arrays
    #
    stationid_1D = np.array(stationid_1D)
    stationid_uniq = np.unique(stationid_1D)
    msg = f"...collected inputs for {len(stationid_uniq)} stations " \
        f"(-->{stationid_uniq}<--"
    logger.debug(msg)
    obstime_1D = np.array(obstime_1D)
    # print(f"stationid_1D -->{stationid_1D}<--")
    # print(f"obstime_1D -->{obstime_1D}<--")
    # print(f"obsmix_1D -->{obsmix_1D}<--")
    # print(f"inic_array1D -->{inic_array1D}<--")
    # print(f"jac_da.shape={jac_da.shape}")

    #->
    #
    #-- instantiate namespace to be returned
    #
    input4inv = SimpleNamespace(
        simu_start=simu_start,
        rundir_list=rundir_list,
        obs_extrainfo=obs_extrainfo_dict,
        stationid_1D=stationid_1D,
        obstime_1D=obstime_1D,
        obslon_1D=obslon_1D,
        obslat_1D=obslat_1D,
        obsalt_1D=obsalt_1D,
        obsmix_1D=obsmix_1D,
        inic_array1D=inic_array1D,
        emis_lonc1D=emis_lonc1D,
        emis_latc1D=emis_latc1D,
        emis_reg1D=emis_reg1D,
        jac_da=jac_da)

    #
    #-- verification
    #
    if refdir!=None:
        # if obsid==None or len(obsid)!=1:
        #     msg = f"...verification currently only supported in case one single " \
        #         f"station is being processed."
        #     raise RuntimeError(msg)
        # staid = obsid[0]
        # assert staid in obs_extrainfo_dict
        for staid in  stationid_uniq:
            idxs_staid = np.where(stationid_1D==staid)[0]
            noday = len(idxs_staid)
            # msg = f"...verification@{staid}, nobs={noday}"
            # logger.debug(msg)
            #
            #--
            #
            cur_inic    = inic_array1D[idxs_staid]
            cur_obstime = obstime_1D[idxs_staid]
            cur_obsc    = obsmix_1D[idxs_staid]
            _nobs,_nemis,_ng = jac_da.shape
            cur_ojac2D = jac_da.values.reshape(_nobs,_nemis*_ng)
            cur_ojac2D = cur_ojac2D[idxs_staid,:]
            #
            #--
            #
            obs_hr  = obs_extrainfo_dict[staid]['obshour']
            obs_twl = obs_extrainfo_dict[staid]['time_window_length']
            msg = f"...@{staid} (obs_hr={obs_hr}, obs_twl={obs_twl}, " \
                f"starting verification against reference run."
            logger.debug(msg)
            #
            #-- emissions used in reference run
            #
            refemisdir = refdir / 'emissions'
            msg = f"reading emissions from reference run directory ***{refemisdir}***..."
            logger.info(msg)
            refemis_info = tm5emisdir_load_emissions2D(refemisdir, 'ch4emis', emisday_range, region_list, remove_halo=remove_halo)
            refemis2D = refemis_info.emis2D
            msg = f"...refemis2D read, (shape={refemis2D.shape})"
            logger.debug(msg)
            msg = f"jac_da.shape = {jac_da.shape}"
            logger.debug(msg)
            if refemis2D.shape!=jac_da.shape[1:]:
                msg = f"...inconsistent shape of observational Jacobian!"
                raise RuntimeError(msg)
            refemis1D = refemis2D.reshape(_nemis*_ng)
            #
            #-- load reference concentrations
            #
            refconc = tm5refdir_load_stationconc(refdir, staid)           
            #
            #-- propagate emissions with linear model
            #
            lindconc = np.dot(cur_ojac2D, refemis1D)
            linsimu_conc = lindconc + cur_inic #-- yields nobs values
            #
            #--
            #
            obstime_tag = cur_obstime[0].strftime('%Y%m%d') + '-' + cur_obstime[-1].strftime('%Y%m%d')
            allobstime_tag = min(obstime_1D).strftime('%Y%m%d') + '-' + max(obstime_1D).strftime('%Y%m%d')
            msg = f"...loading concentration at {obsid} from reference run"
            logger.info(msg)
            logger.debug(f"noday={noday}, cur_obstime={cur_obstime}")
            cmp_df = DataFrame.from_dict(
                {
                    'time': cur_obstime,
                    'refconc': [np.nan,]*noday,
                    'lindeltaconc': [np.nan,]*noday,
                    'linconc': [np.nan,]*noday,
                    'iniconc': [np.nan,]*noday,
                    'obsconc': [np.nan,]*noday
                    }
                )
            cmp_df = cmp_df.set_index('time')
            for iday,obstime in enumerate(cur_obstime):
                _ostart = obstime - Timedelta(seconds=obs_twl)
                _oend   = obstime + Timedelta(seconds=obs_twl)
                cnd_day = (refconc['time']>=_ostart)&(refconc['time']<=_oend)
                refconc_tw = refconc.loc[cnd_day,'conc']
                # msg = f"...refconc --> {refconc_tw.values} ({_ostart} -- {_oend})"
                # logger.debug(msg)
                refconc_day = refconc_tw.mean()
                lindconc_day = lindconc[iday]
                linconc_day = linsimu_conc[iday]
                iniconc_day = cur_inic[iday]
                obsconc_day = cur_obsc[iday]
                cmp_df.loc[obstime,:] = [refconc_day, lindconc_day, linconc_day, iniconc_day, obsconc_day]
                msg = f"...@{staid}/{obstime.strftime('%Y%m%dT%H')} " \
                    f"refconc/linconc/iniconc = {refconc_day}/{linconc_day}/{iniconc_day}"
                logger.info(msg)
            #
            #-- write comparions to csv (done per station)
            #
            outname =  f"obsjac-verification_{staid.replace('_','-')}_simustart-{simu_start.strftime('%Y%m%d')}_{allobstime_tag}.csv"
            if outdir!=None:
                outname = outdir / outname
                outname.parent.mkdir(parents=True, exist_ok=True)
            cmp_df.to_csv(outname, index=True)
            msg = f"...generated comparison file ***{outname}***"
            
    #
    #-- return collected results
    #
    return input4inv
    

def subcmd_test_jacobianfwd_1day(args : ArgumentNamespace) -> None:
    tm5rundir = args.tm5rundir
    obsid     = args.obsid
    if args.trange!=None:
        trange = date_range(args.trange[0], args.trange[1], freq='1D')
        print(trange)
    else:
        trange = args.trange

    df_fwd = load_adjoint_fwd(tm5rundir, trange=trange)

    cndobs = df_fwd.loc[:,'obs']==obsid
    df_obs = df_fwd.loc[cndobs,:]
    print(df_obs)
    #-- sum over regions
    df_out = df_obs.groupby(['obs', 'region']).sum(numeric_only=True).reset_index()
    print(df_out)
    csim = df_obs.loc[:,'mix'].sum()
    print(f"csim={csim}")


def subcmd_debug_footprint_io(args):
    tm5rundir = args.tm5rundir
    if args.emisday_firstlast!=None:
        emisdayf,emisdayl = args.emisday_firstlast
        emisday_range = date_range(emisdayf, emisdayl, freq='1d')
    else:
        emisday_range = None
    #
    #--
    #
    obs_table = tm5rundir_obstable(tm5rundir, drop_missing_value=False)
    obsid = ['cbw_207','cra_60',]
    print(obs_table.head())
    x = obs_table.loc[obsid,:]
    print(x)
    footprint_info = tm5rundir_jacobian3D_update(tm5rundir, emisday_range, obsid=obsid,
                                                 remove_halo=True)
    sys.exit(0)
    #
    #--
    #
    footprint_info = tm5_fitic_adjoint_corrected_halos(tm5rundir, emisday_range=emisday_range)
    

def subcmd_build_jacobian_period_obs1D(args : ArgumentNamespace) -> None:
    """Updated approach for preparation of inputs for Fortran-based inversion environment
    based on TM5 adjoint runs for a selected period of observations.
    Observations (i.e. the combination of station and observational day) are being
    stored as 1D array, and similar for the initial concentration and the first
    dimension (row index) of the Jacobian.
    Sensitivities are being computed w.r.t. to emissions on each day of the temporal
    period.

    On invocation the user has the option to condense the Jacobian to
    reflect sensitivity of the total emission within the selected temporal range.
    """
    topdir = args.outpath_tm5
    domain = args.domain
    obsday_firstlast = args.obsday_firstlast
    emisday_firstlast = args.emisday_firstlast
    obsid = args.obsid
    difdeg_max = args.difdeg_max
    difalt_max = args.difalt_max
    complevel = args.__dict__.get('complevel',4)

    #
    #-- collect inputs
    #
    msg = f"START collecting sensitiy inputs..."
    logger.info(msg)
    obsdayf, obsdayl = obsday_firstlast
    obsday_range = date_range(obsdayf, obsdayl, freq='1d')
    nobsday = len(obsday_range)
    emisdayf, emisdayl = emisday_firstlast
    emisday_range = date_range(emisdayf, emisdayl, freq='1d')
    nemisday = len(emisday_range)
    input4inv = collect_input4inversion_obs1D(topdir, domain, obsday_range, emisday_range,
                                              remove_halo=args.remove_halo,
                                              obsid=args.obsid,
                                              obsdir=args.obsdir)
    msg = f"...input collection FINISHED."
    logger.info(msg)
    #
    #--
    #
    stationid_1D = input4inv.stationid_1D
    obstime_1D = input4inv.obstime_1D
    #--
    obsmix_1D    = input4inv.obsmix_1D
    inic_array1D = input4inv.inic_array1D
    jac_da       = input4inv.jac_da
    nobs, _nemisday, ng = jac_da.shape
    if nemisday!=_nemisday:
        msg = f"inversion input collection supposed for {nemisday} emission days, " \
            f"but shape of resulting Jacobian {jac_array.shape} is unexpected"
        raise RuntimeError(msg)
    #
    #-- Jacobian quantifies deltac [ppb] w.r.t. daily emission rates [kgCH4/cell/s]
    #
    sensitivity_units = jac_da.units
    msg = f"...detected sensitivity units -->{sensitivity_units}<--"
    logger.debug(msg)
    #
    #-- unit convert to
    #   - total emissions (per grid-cell) for the complete period
    #   - need to scale and average entries in Jacobian accordingly
    #
    if args.jac4totemis:
        ojac_daily = jac_da.values
        sensitivity_units = 'ppb/(kgCH4/cell)'
        ojac_da = jac_da.sum(dim='emisday') / (nemisday*nsecday)
        #
        #-- verification of monthly total aggregation
        #
        #-- select emissions from from first rundir
        #
        rundir_first = input4inv.rundir_list[0]
        emis_info = tm5rundir_emissions2D(rundir_first, trange=obsday_range, remove_halo=args.remove_halo)
        emis2D = emis_info.emis2D
        emis_tot = np.sum(emis2D*nsecday, axis=0) #-- overall emissions in temporal range
        msg = f"emis_tot min/mean/max = {emis_tot.min()}/{emis_tot.mean()}/{emis_tot.max()}"
        logger.debug(msg)
        ojac_tot = ojac_da.values
        #
        for iobs in range(nobs):
            _staid  = stationid_1D[iobs]
            _obsday = obstime_1D[iobs].strftime('%Y%m%d')
            dc_tot = np.dot(ojac_tot[iobs,:], emis_tot)
            dc     = np.dot(ojac_daily[iobs,:].ravel(), emis2D.ravel())
            msg = f"@{_obsday},{_staid}: deltacconc derived by daily-rate/temporal-total = " \
                f"{dc}/{dc_tot}"
            logger.info(msg)
        ojac_out_da = ojac_da
        ojac_out = ojac_da.values
    else:
        ojac_out_da = jac_da
        ojac_out = jac_da.values

    #
    #-- re-distribute global flask Jacobian (which has been computed *only* on global 6x4 grid)
    #
    if args.domain=='glb600x400' and args.glb6x4_to_avengers_zoom:
        #
        #-- spatially re-distribute sensitivities to the 1D avengers zoom vector
        #
        jacobian_redistributed = jacobian_redistribute_glb6x4_to_avengers_zoom(ojac_out_da)
        #
        #-- extract the relevant bits
        #
        ojac_out = jacobian_redistributed.jacobian.values
        _, ng = ojac_out.shape #-- update number of grid-cells (!)
        emis_lonc1D = jacobian_redistributed.lonc1D
        emis_latc1D = jacobian_redistributed.latc1D
        emis_reg1D = jacobian_redistributed.reg1D
    else: #-- spatial information (for emissions) can remain as they have been read
        emis_lonc1D = input4inv.emis_lonc1D
        emis_latc1D = input4inv.emis_latc1D
        emis_reg1D  = input4inv.emis_reg1D
    #
    #-- prepare output
    #
    if args.obsid==None:
        station_list = np.unique(stationid_1D)
    else:
        station_list = np.array(args.obsid)
    nsta = len(station_list)
    #
    #-- station coordinates were collected per-observation
    #   (but we expect them not to change from obs to obs at the same station)
    #
    coords_fill = -9999.
    station_coords = np.full((nsta,3), coords_fill) #-- lon/lat/alt
    for ista,sta in enumerate(station_list):
        idxs_sta = np.where(stationid_1D==sta)
        nidxs = len(idxs_sta[0])
        _lon_sta = input4inv.obslon_1D[idxs_sta]
        _lat_sta = input4inv.obslat_1D[idxs_sta]
        _alt_sta = input4inv.obsalt_1D[idxs_sta]
        if nidxs==1:
            station_coords[ista,:] = (_lon_sta[0],_lat_sta[0],_alt_sta[0])
        else:
            #-- longitude
            diflon_max = np.max(np.diff(np.abs(_lon_sta-_lon_sta[0])))
            if diflon_max<=difdeg_max:
                station_coords[ista,0] = _lon_sta[0]
            else:
                msg = f"@{sta}, varying longitudes diflon_max={diflon_max} exceeds threshold " \
                    f"{difdeg_max}"
                logger.debug(msg)
            #-- latitude
            diflat_max = np.max(np.diff(np.abs(_lat_sta-_lat_sta[0])))
            if diflat_max<=difdeg_max:
                station_coords[ista,1] = _lat_sta[0]
            else:
                msg = f"@{sta}, varying latitudes diflat_max={diflat_max} exceeds threshold " \
                    f"{difdeg_max}" 
                logger.debug(msg)
            #-- altitude
            sta_difalt_max = np.max(np.diff(np.abs(_alt_sta-_alt_sta[0])))
            if sta_difalt_max<=difalt_max:
                station_coords[ista,2] = _alt_sta[0]
            else:
                msg = f"@{sta}, varying altitudes sta_difalt_max={sta_difalt_max} exceeds threshold " \
                    f"{difalt_max}"
                logger.debug(msg)
    #--
    stacoords_per_sta = not (coords_fill in station_coords)
    if stacoords_per_sta:
        msg = f"observation coordinates will be written --per-station--"
        logger.info(msg)
    else:
        msg = f"observation coordinates will be written --per-observation--"
        logger.info(msg)
    if nobsday==1:
        obsday_tag = f"obsday-{obsdayf.strftime('%Y%m%d')}"
    else:
        obsday_tag = f"obsdays-{obsdayf.strftime('%Y%m%d')}--{obsdayl.strftime('%Y%m%d')}"
    if nemisday==1:
        emisday_tag = f"{emisdayf.strftime('%Y%m%d')}"
    else:
        emisday_tag = f"{emisdayf.strftime('%Y%m%d')}--{emisdayl.strftime('%Y%m%d')}"
    if args.jac4totemis:
        emisday_tag = f"wrt-totalemis-{emisday_tag}"
    else:
        emisday_tag = f"wrt-dailyemis-{emisday_tag}"
    if nsta==1:
        obsid_tag = station_list[0]
    elif nsta<=5:
        obsid_tag = '--' + '--'.join([_ for _ in station_list]) + '--'
    else:
        obsid_tag = f"{nsta}-obslocations"
    if args.domain=='glb600x400' and args.glb6x4_to_avengers_zoom:
        domain_tag = f"{args.domain}-to-gns100x100"
    else:
        domain_tag = args.domain
    outname_tokens = ["fitic-inversion-input-obs1D", obsid_tag, domain_tag, obsday_tag, emisday_tag]
    if args.remove_halo:
        outname_tokens.append('removed-halos')
    outname = '_'.join(outname_tokens) + '.nc'
    outname = set_outname(args, outname)
    msg = f"writing inversion inputs to file ***{outname}***..."
    logger.info(msg)
    #
    #-- spatial dimensions
    #
    fp = Dataset(outname, 'w')
    fp.createDimension('ng', ng)
    fp.createDimension('nobs', nobs)
    if not args.jac4totemis:
        fp.createDimension('nemisday', nemisday)
    fp.createDimension('nsta', nsta)
    #
    #-- unique list of stations
    #
    ncvar = fp.createVariable('station_id', str, ('nsta',))
    ncvar[:] = station_list[:]
    ncvar.long_name = f"station_identifier_list"
    ncvar.units = ''
    ncvar.comment = f"Comprises the overall list of stations. Note, that there may be no observations for a station on certain day(s)."
    if stacoords_per_sta:
        #-- longitude
        ncvar = fp.createVariable('station_lon', 'f8', ('nsta',))
        ncvar[:] = station_coords[:,0]
        ncvar.long_name = 'station_longitude'
        ncvar.units = 'degrees_east'
        #-- latitude
        ncvar = fp.createVariable('station_lat', 'f8', ('nsta',))
        ncvar[:] = station_coords[:,1]
        ncvar.long_name = 'station_longitude'
        ncvar.units = 'degrees_north'
        #-- altitude
        ncvar = fp.createVariable('station_alt', 'f8', ('nsta',))
        ncvar[:] = station_coords[:,2]
        ncvar.long_name = 'station_altitude'
        ncvar.units = 'm'

    #
    ncvar = fp.createVariable('lon', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'longitude'
    ncvar.units = 'degrees_east'
    ncvar.comment = 'references center of grid-cell in related zoom domain'
    ncvar[:] = emis_lonc1D[:]
    #
    ncvar = fp.createVariable('lat', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'latitude'
    ncvar.units = 'degrees_north'
    ncvar.comment = 'references center of grid-cell in related zoom domain'
    ncvar[:] = emis_latc1D[:]
    #
    ncvar = fp.createVariable('region', input4inv.emis_reg1D.dtype, ('ng',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.units = ''
    ncvar[:] = emis_reg1D[:]
    #
    ncvar = fp.createVariable('obs', 'f8', ('nobs',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = obsmix_1D[:]
    ncvar.long_name = f"observed CH4 concentration"
    ncvar.units = 'ppb'
    #
    ncvar = fp.createVariable('iniconc', 'f8', ('nobs',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = inic_array1D[:]
    ncvar.long_name = f"initial_concentration"
    ncvar.units = 'ppb'
    #
    ncvar = fp.createVariable('station', str, ('nobs',) )
    ncvar[:] = stationid_1D[:]
    ncvar.long_name = 'station_identifier'
    ncvar.units = ''
    #
    ncvar = fp.createVariable('obstime', str, ('nobs',) )
    ncvar[:] = np.array([ _.strftime('%Y%m%dT%H%M%S') for _ in obstime_1D ])
    # for iobs in range(nobs):
    #     ncvar[iobs] = obstime_1D[iobs].strftime('%Y%m%dT%H')
    ncvar.long_name = 'time_of_observation'
    ncvar.units = ''
    if not stacoords_per_sta:
        #
        ncvar = fp.createVariable('obslon', 'f8', ('nobs',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = input4inv.obslon_1D[:]
        ncvar.long_name = 'longitude_of_observation'
        ncvar.units = 'degrees_east'
        #
        ncvar = fp.createVariable('obslat', 'f8', ('nobs',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = input4inv.obslat_1D[:]
        ncvar.long_name = 'latitude_of_observation'
        ncvar.units = 'degrees_north'
        #
        ncvar = fp.createVariable('obsalt', 'f8', ('nobs',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = input4inv.obsalt_1D[:]
        ncvar.long_name = 'altitude_of_observation'
        ncvar.units = 'm'

    #
    if args.jac4totemis:
        ncvar = fp.createVariable('obs_jacobian', 'f8', ('nobs','ng',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = ojac_out[:]
        ncvar.units = sensitivity_units
        ncvar.comment = f"Jacobian quantifies the sensitivity of concentration at " \
            f"observed times and locations w.r.t. to the total emission field in the " \
            f"temporal range from {emisdayf.strftime('%Y%m%d')} to {emisdayl.strftime('%Y%m%d')}"
    else:
        ncvar = fp.createVariable('obs_jacobian', 'f8', ('nobs','nemisday','ng',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = ojac_out[:]
        ncvar.units = sensitivity_units
        ncvar.comment = f"Jacobian quantifies the sensitivity of concentration at " \
            f"observed times and locations w.r.t. to daily emission rates."

    #
    #-- global attributes
    #
    fp.description = f"Inputs for FIT-IC inversion environment comprising the observational" \
        f"Jacobian, the initial concentrations, and the observed concentrations " \
        f"at the selected stations and for the selected days."
    fp.footprint_directory = str(topdir.absolute())
    fp.removed_halos = np.int32(args.remove_halo)
    try:
        fp.processing_platform = f"{os.environ['USER']}@{os.environ['HOSTNAME']}"
    except KeyError:
        pass
    fp.history = f"{' '.join(sys.argv)}"
    fp.date_created = Timestamp.utcnow().isoformat()
    #
    #-- close
    #
    fp.close()
    msg = f"generated file ***{outname}***"
    logger.info(msg)


def subcmd_monthly_obsjacobian_for_inversion(args  : ArgumentNamespace) -> None:
    """
    """
    topdir = args.outpath_tm5
    domain = args.domain
    obsid = args.obsid
    obsmon_firstlast = args.obsmon_firstlast
    emisday_first = args.emisday_first
    difdeg_max = args.difdeg_max
    difalt_max = args.difalt_max
    complevel = args.__dict__.get('complevel',4)
    #
    #--
    #
    if domain=='gns100x100':
        region_list = ['glb600x400','eur300x200','gns100x100']
    elif domain=='glb600x400':
        region_list = ['glb600x400']
    else:
        msg = f"unexpected domain ==>{domain}<--"
        raise RuntimeError(msg)
    reginfo = regions1D_info(region_list, remove_halo=args.remove_halo)
    ng = reginfo.ng
    #
    #-- range of observation months
    #
    monf,monl = obsmon_firstlast
    #- normalise to day1 in month
    obsdayf = monf.replace(day=1)
    obsdayl = (monl.replace(day=1) + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
    obsdayf, obsdayl = obsmon_firstlast
    obsmon_range = date_range(obsdayf, obsdayl, freq='MS')
    nobsmon = len(obsmon_range)
    msg = f"obsmon_range -->{obsmon_range}<--"
    logger.debug(msg)
    obsdayf, obsdayl = obsmon_firstlast
    obsday_range = date_range(obsdayf, obsdayl, freq='1D')
    nobsday = len(obsday_range)
    #
    #-- 
    #
    if emisday_first.day!=1:
        msg = f"...first day of emissions must be on the first day of month, " \
            f"but got -->{emisday_first}<--"
        raise RuntimeError(msg)
    emismon_range = date_range(emisday_first, obsdayl, freq='MS')
    nemismon = len(emismon_range)
    msg = f"emismon_range -->{emismon_range}<--"
    logger.debug(msg)
    njaccol = nemismon*ng
    #
    #--
    #
    sensitivity_units = 'ppb/(kgCH4/cell/month)'
    stationid_1D   = None
    obslon_1D      = None
    obslat_1D      = None
    obsalt_1D      = None
    obstime_1D     = None
    obsmix_1D      = None
    inic_array1D   = None
    jac_da         = None
    jac3D_da       = None
    jaccol_emismon = np.array([_.strftime('%Y-%m-%d') for _ in emismon_range]*ng)
    jaccol_emisgc  = np.tile(np.arange(ng,dtype='i4'), nemismon)
    emis_lonc1D    = None
    emis_latc1D    = None
    emis_reg1D     = None
    #
    #-- loop in monthly bins for the observations
    #
    for imon,curobsdayf in enumerate(obsmon_range):
        curobsdayl = (curobsdayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
        curobsday_range = date_range(curobsdayf, curobsdayl, freq='1D')
        # #-- define shorter time-span when debugging
        # dayf = max(curobsdayf, obsdayf)
        # dayl = min(curobsdayl, obsdayl)
        # curobsday_range = date_range(dayf, dayl, freq='1D')
        curemisdayf = emisday_first
        curemisdayl = (curobsdayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
        curemisday_range = date_range(curemisdayf, curemisdayl, freq='1D')
        nemisday = len(curemisday_range)
        msg = f"start collection for month starting on {curobsdayf.strftime('%Y-%m-%d')})"
        logger.info(msg)
        input4inv = collect_input4inversion_obs1D(topdir, domain, curobsday_range, curemisday_range,
                                                  remove_halo=args.remove_halo,
                                                  obsid=args.obsid,
                                                  obsdir=args.obsdir,
                                                  refdir=args.refdir,
                                                  outdir=args.outdir)
        # #-- 
        curjacdaily_da = input4inv.jac_da
        curnobs,curnemisday,_ = curjacdaily_da.shape
        if nemisday!=curnemisday:
            msg = f"inversion input collection supposed for {nemisday} emission days, " \
                f"but shape of resulting Jacobian {jac_array.shape} is unexpected"
            raise RuntimeError(msg)
        #
        #-- aggregate to sensitivities per month
        #
        curjac_array = zeros((curnobs,njaccol))
        #
        #-- allocate buffer for current observation month
        #
        #   TODO::currently two type of buffers are created,
        #         buf[obs,njaccol] and buf[obs,emismon,ng],
        #         Eventually we need to select which variant is best for inversion environment
        #
        # curjac_da = xr.DataArray(
        #     zeros((curnobs,njaccol)),
        #     dims=('obs','njaccol'),
        #     coords={'obs':curjacdaily_da.obs, 'njaccol': np.arange(njaccol)},
        #     attrs = {'units': sensitivity_units}
        #     )
        curjac3D_da = xr.DataArray(
            zeros((curnobs,nemismon,ng)),
            dims=('obs','emismon','ng'),
            coords={'obs':curjacdaily_da.obs, 'emismon':emismon_range, 'ng': np.arange(ng)},
            attrs = {'units': sensitivity_units}
            )
        #
        #-- loop over (daily) emissions per month
        #
        for iemismon,emismondayf in enumerate(emismon_range):
            #-- no sensitivity for emission months that are past the last observation
            if emismondayf>max(curemisday_range):
                break
            #
            msg = f"...getting sensitivities for emismondayf={emismondayf.strftime('%Y-%m-%d')}"
            logger.info(msg)
            #
            #-- extraction sensitivities for current emission month
            #
            cnd_emismon = (curjacdaily_da.emisday>=curemisdayf)&(curjacdaily_da.emisday<=curemisdayl)
            _jac_da = curjacdaily_da.sel(emisday=cnd_emismon)
            msg = f"......_jac_da.shape={_jac_da.shape}"
            logger.debug(msg)
            #
            #-- unit conversion [ppb/kgCH4/cell/s] --> [ppb/kgCH4/cell/month]
            #
            nsec = _jac_da.shape[1]*nsecday #-- number of emission seconds [in month]
            _jac_da = _jac_da.sum(dim='emisday') / nsec
            msg = f"......after summing-up for month, _jac_da.shape={_jac_da.shape}"
            logger.debug(msg)
            #
            #-- insert into overall array
            #
            icolf = iemismon*ng
            icoll = icolf+ng-1
            # curjac_da.loc[dict(njaccol=slice(icolf,icoll))] = _jac_da.values
            curjac3D_da.loc[dict(emismon=slice(curemisdayf,curemisdayl))] = _jac_da.values.reshape((curnobs,1,ng))
        #
        #--
        #
        if imon==0:
            stationid_1D  = input4inv.stationid_1D
            obslon_1D     = input4inv.obslon_1D
            obslat_1D     = input4inv.obslat_1D
            obsalt_1D     = input4inv.obsalt_1D
            obstime_1D    = input4inv.obstime_1D
            inic_array1D  = input4inv.inic_array1D
            obsmix_1D     = input4inv.obsmix_1D
            # jac_da        = curjac_da
            jac3D_da      = curjac3D_da
            emis_lonc1D   = input4inv.emis_lonc1D
            emis_latc1D   = input4inv.emis_latc1D
            emis_reg1D    = input4inv.emis_reg1D
        else:
            stationid_1D   = np.hstack((stationid_1D, input4inv.stationid_1D))
            obslon_1D      = np.hstack((obslon_1D, input4inv.obslon_1D))
            obslat_1D      = np.hstack((obslat_1D, input4inv.obslat_1D))
            obsalt_1D      = np.hstack((obsalt_1D, input4inv.obsalt_1D))
            obstime_1D     = np.hstack((obstime_1D, input4inv.obstime_1D))
            inic_array1D   = np.hstack((inic_array1D, input4inv.inic_array1D))
            obsmix_1D      = np.hstack((obsmix_1D, input4inv.obsmix_1D))
            # jac_da         = xr.concat([jac_da, curjac_da], dim='obs')
            jac3D_da       = xr.concat([jac3D_da, curjac3D_da], dim='obs')
            if np.any(emis_lonc1D!=input4inv.emis_lonc1D):
                msg = f"@{curobsdayf},imon={imon}, inconsistency in emission longitudes"
                raise RuntimeError(msg)
            if np.any(emis_latc1D!=input4inv.emis_latc1D):
                msg = f"@{curobsdayf},imon={imon}, inconsistency in emission latitudes"
                raise RuntimeError(msg)
            if np.any(emis_reg1D!=input4inv.emis_reg1D):
                msg = f"@{curobsdayf},imon={imon}, inconsistency in emission region identifiers"
                raise RuntimeError(msg)
    #
    #--
    #
    nobs, _nemisdays,_ng = jac3D_da.shape
    msg = f"...Jacobian array created (dimensions={jac3D_da.dims}, shape={jac3D_da.shape})."
    logger.debug(msg)
    #
    #-- re-distribute global flask Jacobian (which has been computed *only* on global 6x4 grid)
    #
    if args.domain=='glb600x400' and args.glb6x4_to_avengers_zoom:
        #
        #-- perfrom the re-distribution to the AVENGERS zoom domain
        #
        jac3D_redistributed = jacobian_redistribute_glb6x4_to_avengers_zoom(jac3D_da)
        #
        #-- update local variables for output accordingly
        #
        jac3D_da    = jac3D_redistributed.jacobian
        emis_lonc1D = jac3D_redistributed.lonc1D
        emis_latc1D = jac3D_redistributed.latc1D
        emis_reg1D  = jac3D_redistributed.reg1D
        ngout = len(emis_lonc1D)       
    else:
        ngout = ng
        
    #
    #-- determine unique stations
    #
    if args.obsid==None:
        station_list = np.unique(stationid_1D)
    else:
        station_list = np.array(args.obsid)
    nsta = len(station_list)
    msg = f"...Jacobian prepared for {nsta} different stations (-->{station_list}<--)"
    logger.info(msg)
    #
    #-- station coordinates were collected per-observation
    #   - we expect them not to change from obs to obs at the same station,
    #     but there occurred small differences (at least in case the observations
    #     were extracted from the origin obspack NetCDF files)
    #
    coords_fill = -9999.
    station_coords = np.full((nsta,3), coords_fill) #-- lon/lat/alt
    for ista,sta in enumerate(station_list):
        idxs_sta = np.where(stationid_1D==sta)
        nidxs = len(idxs_sta[0])
        _lon_sta = obslon_1D[idxs_sta]
        _lat_sta = obslat_1D[idxs_sta]
        _alt_sta = obsalt_1D[idxs_sta]
        if nidxs==1:
            station_coords[ista,:] = (_lon_sta[0],_lat_sta[0],_alt_sta[0])
        else:
            #-- longitude
            diflon_max = np.max(np.diff(np.abs(_lon_sta-_lon_sta[0])))
            if diflon_max<=difdeg_max:
                station_coords[ista,0] = _lon_sta[0]
            else:
                msg = f"@{sta}, varying longitudes diflon_max={diflon_max} exceeds threshold " \
                    f"{difdeg_max}"
                logger.debug(msg)
            #-- latitude
            diflat_max = np.max(np.diff(np.abs(_lat_sta-_lat_sta[0])))
            if diflat_max<=difdeg_max:
                station_coords[ista,1] = _lat_sta[0]
            else:
                msg = f"@{sta}, varying latitudes diflat_max={diflat_max} exceeds threshold " \
                    f"{difdeg_max}" 
                logger.debug(msg)
            #-- altitude
            sta_difalt_max = np.max(np.diff(np.abs(_alt_sta-_alt_sta[0])))
            if sta_difalt_max<=difalt_max:
                station_coords[ista,2] = _alt_sta[0]
            else:
                msg = f"@{sta}, varying altitudes sta_difalt_max={sta_difalt_max} exceeds threshold " \
                    f"{difalt_max}"
                logger.debug(msg)
    #--
    stacoords_per_sta = not (coords_fill in station_coords)
    if stacoords_per_sta:
        msg = f"observation coordinates will be written --per-station--"
        logger.info(msg)
    else:
        msg = f"observation coordinates will be written --per-observation--"
        logger.info(msg)
    #
    #-- prepare output
    #
    if nobsday==1:
        obsday_tag = f"obsday-{obsdayf.strftime('%Y%m%d')}"
    else:
        obsday_tag = f"obsdays-{obsdayf.strftime('%Y%m%d')}--{obsdayl.strftime('%Y%m%d')}"
    emismon_tag = f"wrt-monthlyemis-{emismon_range[0].strftime('%Y%m')}--{emismon_range[-1].strftime('%Y%m')}"
    if nsta==1:
        obsid_tag = station_list[0]
    elif nsta<=5:
        obsid_tag = '--' + '--'.join([_ for _ in station_list]) + '--'
    else:
        obsid_tag = f"{nsta}-obslocations"
    if args.domain=='glb600x400' and args.glb6x4_to_avengers_zoom:
        domain_tag = f"{args.domain}-to-gns100x100"
    else:
        domain_tag = args.domain
    outname_tokens = ["fitic-inversion-input", obsid_tag, domain_tag, obsday_tag, emismon_tag]
    if args.remove_halo:
        outname_tokens.append('removed-halos')
    outname = '_'.join(outname_tokens) + '.nc'
    outname = set_outname(args, outname)
    msg = f"writing inversion inputs to file ***{outname}***..."
    logger.info(msg)
    #
    #-- create dimensions
    #
    n_strlen = 32
    fp = Dataset(outname, 'w')
    fp.createDimension('ng', ngout)
    fp.createDimension('nemismon', nemismon)
    fp.createDimension('nobs', nobs)
    # fp.createDimension('njaccol', njaccol)
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
    ncvar[:] = emis_lonc1D[:]
    #
    #-- latitude
    #
    ncvar = fp.createVariable('lat', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'latitude'
    ncvar.units = 'degrees_north'
    ncvar.comment = 'references center of grid-cell in related zoom domain'
    ncvar[:] = emis_latc1D[:]
    #
    #-- region identifier
    #
    ncvar = fp.createVariable('region', emis_reg1D.dtype, ('ng',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.units = ''
    ncvar[:] = emis_reg1D[:]
    #
    #-- region identifier
    #
    ncvar = fp.createVariable('region_ftn', 'S1', ('ng','nstrlen',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.comment = f"region identifer in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    ncvar[:] = stringtochar(emis_reg1D[:], n_strlen=n_strlen)
    #
    #-- observed concentration
    #
    ncvar = fp.createVariable('obs', 'f8', ('nobs',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = obsmix_1D[:]
    ncvar.long_name = f"observed CH4 concentration"
    ncvar.units = 'ppb'
    #
    #-- initial concentration
    #
    ncvar = fp.createVariable('iniconc', 'f8', ('nobs',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = inic_array1D[:]
    ncvar.long_name = f"initial_concentration"
    ncvar.units = 'ppb'
    #
    #-- station identifier (per observation)
    #
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
    ncvar[:] = station_list[:]
    ncvar.long_name = f"station_identifier_list"
    ncvar.units = ''
    ncvar.comment = f"Comprises the overall list of stations. Note, that there may be no observations for a station on certain day(s)."
    ncvar = fp.createVariable('station_id_ftn', 'S1', ('nsta','nstrlen',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = stringtochar(station_list[:], n_strlen=n_strlen)
    ncvar.long_name = 'station_identifier'
    ncvar.comment = f"station identifier in a format which is suitable " \
        f"for Fortran based I/O"
    ncvar.units = ''
    
    if stacoords_per_sta:
        #-- longitude
        ncvar = fp.createVariable('station_lon', 'f8', ('nsta',))
        ncvar[:] = station_coords[:,0]
        ncvar.long_name = 'station_longitude'
        ncvar.units = 'degrees_east'
        #-- latitude
        ncvar = fp.createVariable('station_lat', 'f8', ('nsta',))
        ncvar[:] = station_coords[:,1]
        ncvar.long_name = 'station_longitude'
        ncvar.units = 'degrees_north'
        #-- altitude
        ncvar = fp.createVariable('station_alt', 'f8', ('nsta',))
        ncvar[:] = station_coords[:,2]
        ncvar.long_name = 'station_altitude'
        ncvar.units = 'm'

    else:#if not stacoords_per_sta:
        #
        ncvar = fp.createVariable('obslon', 'f8', ('nobs',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = obslon_1D[:]
        ncvar.long_name = 'longitude_of_observation'
        ncvar.units = 'degrees_east'
        #
        ncvar = fp.createVariable('obslat', 'f8', ('nobs',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = obslat_1D[:]
        ncvar.long_name = 'latitude_of_observation'
        ncvar.units = 'degrees_north'
        #
        ncvar = fp.createVariable('obsalt', 'f8', ('nobs',),
                                  compression='zlib', complevel=complevel)
        ncvar[:] = obsalt_1D[:]
        ncvar.long_name = 'altitude_of_observation'
        ncvar.units = 'm'
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
    # #
    # #-- 2D Jacobian dataset
    # #
    # ncvar = fp.createVariable('obs_jacobian2D', 'f8', ('nobs','njaccol',),
    #                           compression='zlib', complevel=complevel)
    # ncvar[:] = jac_da[:]
    # ncvar.units = jac_da.attrs['units']
    # ncvar.comment = f"Jacobian quantifies the sensitivity of concentration at " \
    #     f"observed times and locations w.r.t. to monthly emissions."
    # #-- 2D Jacobian ancillary information for columns
    # #> NOTE: compression seems  problematic for 'str' typed variable
    # ncvar = fp.createVariable('jaccol_emismon', jaccol_emismon.dtype, ('njaccol',))
    # ncvar[:] = jaccol_emismon[:]
    # ncvar.long_name = "emission_month_of_column_in_jacobian"
    # ncvar.units = ''
    # ncvar = fp.createVariable('jaccol_gridcell_index', jaccol_emisgc.dtype, ('njaccol',),
    #                           compression='zlib', complevel=complevel)
    # ncvar[:] = jaccol_emisgc[:]
    # ncvar.long_name = "0-based index of emission grid-cell of column in Jacobian"
    # ncvar.units = ''
    #
    #-- 3D Jabobian dataset
    #
    ncvar = fp.createVariable('obs_jacobian', 'f8', ('nobs','nemismon','ng',),
                              compression='zlib', complevel=complevel)
    ncvar[:] = jac3D_da[:]
    ncvar.units = jac3D_da.attrs['units']
    ncvar.comment = f"Jacobian quantifies the sensitivity of concentration at " \
        f"observed times and locations w.r.t. to monthly emissions."
    #
    #-- global attributes
    #
    fp.description = f"Inputs for FIT-IC inversion environment comprising the observational" \
        f"Jacobian, as well as initial and observed concentrations " \
        f"at selected stations and for the selected days."
    fp.footprint_directory = str(topdir.absolute())
    fp.removed_halos = np.int32(args.remove_halo)
    try:
        fp.processing_platform = f"{os.environ['USER']}@{os.environ['HOSTNAME']}"
    except KeyError:
        pass
    fp.history = f"{' '.join(sys.argv)}"
    fp.date_created = Timestamp.utcnow().isoformat()
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
    nmon = len(mon_range)
    ntc = 3 #-- recording year/month/day
    time_data = np.full((nmon,ntc), -1)
    for imon,mon in enumerate(mon_range):
        time_data[imon,:] = [mon.year, mon.month, mon.day] #-- deliberately take first day
    day_first = mon_range[0]
    day_last  = (mon_range[-1] + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
    month_tag = f"{day_first.strftime('%Y%m%d')}--{day_last.strftime('%Y%m%d')}"
    #
    #-- spatial configuration
    #
    reginfo = regions1D_info(regions, remove_halo=args.remove_halo)
    ng = reginfo.ng
    msg = f"-->{regions}<-- yield overall {ng} grid-cells"
    logger.info(msg)

    #
    #-- debugging / control
    #
    reg_table = reginfo.table
    for reg,reg_data in reg_table.items():
        _ng = reg_data.ng1D
        _ng_nohalo = reg_data.ng1D_nohalo
        msg = f"@{reg}, ng/ng_nohalo/difference = {_ng}/{_ng_nohalo}/{_ng-_ng_nohalo}"
        logger.info(msg)
    #
    #-- determine list of months
    #
    region_file_table = OrderedDict()
    nday = None
    for ireg,reg in enumerate(regions):
        region_file_table[reg] = []
        for mondate in mon_range:
            file_list = sorted(Path(tm5emisdir).glob(f"ch4emis.CH4.{reg}.{mondate.year}{mondate.month}??.nc"))
            region_file_table[reg] += file_list
        if ireg==0:
            nday = len(region_file_table[reg])
        else:
            assert nday==len(region_file_table[reg])
    #
    #-- initialise array for emissions
    #
    nsecday = 86400
    emis_miss = -99999.
    emis_data = np.full((nmon,ng), emis_miss)
    #-- fill array
    for imon,dayf in enumerate(mon_range):
        dayl = (dayf + Timedelta(days=32)).replace(day=1) - Timedelta(days=1)
        day_range = date_range(dayf,dayl)
        msg = f"...loading emissions for {dayf.strftime('%Y%m%d')} to {dayl.strftime('%Y%m%d')}"
        logger.info(msg)
        #-- load daily emissions for every day in month
        emis_info = tm5emisdir_load_emissions2D(tm5emisdir, tm5emisdir+'/ch4emis', day_range, regions, remove_halo=args.remove_halo)
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
    fp.removed_halos = np.int32(args.remove_halo)
    fp.time_coverage_start = day_first.strftime('%Y-%m-%d')
    fp.time_coverage_end   = day_last.strftime('%Y-%m-%d')
    fp.time_coverage_resolution = "P1M"
    try:
        fp.processing_platform = f"{os.environ['USER']}@{os.environ['HOSTNAME']}"
    except KeyError:
        pass
    fp.history = f"{' '.join(sys.argv)}"
    fp.date_created = Timestamp.utcnow().isoformat()
    #
    #-- close
    #
    fp.close()
    msg = f"generated file ***{outname}***"
    logger.info(msg)

    
def subcmd_create_target_jacobian(args : ArgumentNamespace) -> None:
    """
    """
    complevel = args.__dict__.get('complevel',4)
    #
    #--
    #
    regions = ['glb600x400', 'eur300x200', 'gns100x100',]
    reginfo = regions1D_info(regions, remove_halo=args.remove_halo)
    ng = reginfo.ng
    region_table = reginfo.table
    area1D = reginfo.area1D
    msg = f"-->{regions}<-- yield overall {ng} grid-cells"
    logger.info(msg)

    eur_nohalo_w = region_table['eur300x200'].grid.west + region_table['glb600x400'].grid.dlon
    eur_nohalo_e = region_table['eur300x200'].grid.east - region_table['glb600x400'].grid.dlon
    eur_nohalo_s = region_table['eur300x200'].grid.south + region_table['glb600x400'].grid.dlat
    eur_nohalo_n = region_table['eur300x200'].grid.north - region_table['glb600x400'].grid.dlat

    #
    #-- dedicated country targets processed
    #   *only* in the 1x1 innermost zoom domain with HALO parts removed
    #
    gns_grid = region_table['gns100x100'].grid
    gns_w,gns_e,gns_s,gns_n = gns_grid.domain
    gns_nohalo_w = gns_w + region_table['eur300x200'].grid.dlon
    gns_nohalo_e = gns_e - region_table['eur300x200'].grid.dlon
    gns_nohalo_s = gns_s + region_table['eur300x200'].grid.dlat
    gns_nohalo_n = gns_n - region_table['eur300x200'].grid.dlat
    gns_nohalo_grid = TM5Grids.from_corners(west=gns_nohalo_w,east=gns_nohalo_e,
                                            south=gns_nohalo_s,north=gns_nohalo_n,
                                            dlon=gns_grid.dlon,dlat=gns_grid.dlat)
    #
    cnd_gns_nohalo = (reginfo.reg1D=='gns100x100') & \
        (reginfo.lonc1D>=gns_nohalo_grid.west)&(reginfo.lonc1D<=gns_nohalo_grid.east) & \
        (reginfo.latc1D>=gns_nohalo_grid.south)&(reginfo.latc1D<=gns_nohalo_grid.north)
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
        country_area_inner = country_frct_inner*gns_nohalo_grid.area[np.newaxis,:,:]
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
                        tgt_country_table[cid] = country_frct_inner[ic,:].ravel()
    #
    #--
    #
    
    #
    #-- target jacobian
    #
    target_list = ['global', 'zoom_domain',] + list(tgt_country_table.keys())
    ntgt = len(target_list)
    tjac2D = zeros((ntgt,ng), dtype='f8')
    for itgt,tgt in enumerate(target_list):
        if tgt=='global':
            #-- ATTENTION:: avoid double counting of grid-cells,
            #               setting all entries to 1. is !WRONG!
            # tjac2D[itgt,:] = 1.
            #
            #-- from coarse to finer
            #
            # (1) global domain  --> remove inner part (except for HALO)
            cnd = (reginfo.reg1D=='glb600x400') & \
                ( (reginfo.lonc1D<=eur_nohalo_w) | (reginfo.lonc1D>=eur_nohalo_e) | \
                  (reginfo.latc1D<=eur_nohalo_s) | (reginfo.latc1D>=eur_nohalo_n) )
            ncnd = np.count_nonzero(cnd)
            msg = f"...glb600x400 part contributing with {ncnd} grid-cells to global target Jacobian"
            logger.info(msg)
            tjac2D[itgt,cnd] = 1.
            #
            # (2) european domain --> remove inner part covered by zoom (except for HALO)
            cnd = (reginfo.reg1D=='eur300x200') & \
                ( (reginfo.lonc1D<=gns_nohalo_w) | (reginfo.lonc1D>=gns_nohalo_e) | \
                  (reginfo.latc1D<=gns_nohalo_s) | (reginfo.latc1D>=gns_nohalo_n) )
            #-- but also remove the HALO part in the domain itself
            cnd &= (reginfo.reg1D=='eur300x200') & \
                ( (reginfo.lonc1D>=eur_nohalo_w)&(reginfo.lonc1D<=eur_nohalo_e) & \
                  (reginfo.latc1D>=eur_nohalo_s)&(reginfo.latc1D<=eur_nohalo_n) )
            ncnd = np.count_nonzero(cnd)
            msg = f"...eur300x200 part contributing with {ncnd} grid-cells to global target Jacobian"
            logger.info(msg)
            tjac2D[itgt,cnd] = 1.
            # (3) inner-most domain
            ncnd = np.count_nonzero(cnd_gns_nohalo)
            msg = f"...gns100x100 part contributing with {ncnd} grid-cells to global target Jacobian"
            logger.info(msg)
            tjac2D[itgt, cnd_gns_nohalo] = 1.
        elif tgt=='zoom_domain':
            #
            #-- part1: grid-cells in inner-most zoom domain with HALO removed
            #
            tjac2D[itgt,cnd_gns_nohalo] = 1.
            #
            #-- part2: grid-cells from eur300x200 for the HALO part of inner zoom domain
            #          -> this is done in two steps
            cnd_1 = (reginfo.reg1D=='eur300x200') & \
                ( (reginfo.lonc1D>=gns_w) & (reginfo.lonc1D<=gns_e) ) & \
                ( (reginfo.latc1D>=gns_s) & (reginfo.latc1D<=gns_n) )
            tjac2D[itgt, cnd_1] = 1.
            cnd_2 = (reginfo.reg1D=='eur300x200') & \
                ( (reginfo.lonc1D>=gns_nohalo_w) & (reginfo.lonc1D<=gns_nohalo_e) ) & \
                ( (reginfo.latc1D>=gns_nohalo_s) & (reginfo.latc1D<=gns_nohalo_n) )
            tjac2D[itgt, cnd_2] = 0.
        elif tgt in tgt_country_table:
            tgt_frct = tgt_country_table[tgt]
            #-- just one more consistency check
            tgt_area = np.sum(tgt_frct*area1D[cnd_gns_nohalo])/1e6
            msg = f"inserting grid cell fractions for target -->{tgt}<-- (area: {tgt_area:.2f}[km2])"
            logger.info(msg)
            tjac2D[itgt,cnd_gns_nohalo] = tgt_frct
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
    fp = Dataset(outname, 'w')
    fp.createDimension('ng', ng)
    fp.createDimension('ntgt', ntgt)
    #
    ncvar = fp.createVariable('lon', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'longitude'
    ncvar.units = 'degrees_east'
    ncvar.comment = 'references center of grid-cell in related zoom domain'
    ncvar[:] = reginfo.lonc1D
    #
    ncvar = fp.createVariable('lat', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'latitude'
    ncvar.units = 'degrees_north'
    ncvar.comment = 'references center of grid-cell in related zoom domain'
    ncvar[:] = reginfo.latc1D
    #
    ncvar = fp.createVariable('region', reginfo.reg1D.dtype, ('ng',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.units = ''
    ncvar[:] = reginfo.reg1D[:]
    #
    ncvar = fp.createVariable('area', 'f8', ('ng',),
                              compression='zlib', complevel=complevel)
    ncvar.units = 'm2'
    ncvar[:] = area1D
    #
    ncvar = fp.createVariable('targets', str, ('ntgt',))
    ncvar.long_name = f"target_identifier"
    ncvar.units = ''
    for itgt,tgt in enumerate(target_list):
        ncvar[itgt] = tgt
    #
    ncvar = fp.createVariable('tgt_jacobian', 'f8', ('ntgt','ng',),
                              compression='zlib', complevel=complevel)
    ncvar.units = ''
    ncvar[:] = tjac2D[:]

    #
    #-- global attributes
    #
    fp.description = f"Target Jacobian for Fortran inversion environment within FIT-IC"
    fp.removed_halos = np.int32(args.remove_halo)
    try:
        fp.processing_platform = f"{os.environ['USER']}@{os.environ['HOSTNAME']}"
    except KeyError:
        pass
    fp.history = f"{' '.join(sys.argv)}"
    fp.date_created = Timestamp.utcnow().isoformat()
    #
    #-- close
    #
    fp.close()
    msg = f"generated file ***{outname}***"
    logger.info(msg)


def subcmd_concat_ojac_obs1D(args):
    """
    """
    filepath_ojac_cont = args.filepath_ojac_cont
    filepath_ojac_flask = args.filepath_ojac_flask
    stations_flask = args.stations_flask
    complevel = args.__dict__.get('complevel',4)

    #
    #-- continuous obs Jacobian
    #
    ds_ojac_cont = xr.open_dataset(filepath_ojac_cont)

    #
    #-- flask obs Jacobian
    #
    ds_ojac_flask = xr.open_dataset(filepath_ojac_flask)
#    print(ds_ojac_flask)
    uniq_station_ids = ds_ojac_flask['station_id'].values
    msg = f"flask station  identifiers -->{uniq_station_ids}<-- " \
        f"obs_jacobian shape={ds_ojac_flask.obs_jacobian.shape}"
    logger.debug(msg)
    #
    #-- filter for only selected stations
    #
    idxs_uniq_station_ids = np.where(np.isin(uniq_station_ids, stations_flask))
    idxs_uniq_station_ids = idxs_uniq_station_ids[0]
    obs_station_ids = ds_ojac_flask['station'].values
    idxs_stations_flask = np.where(np.isin(obs_station_ids, stations_flask))
    idxs_stations_flask = idxs_stations_flask[0]
    ds_ojac_flask = ds_ojac_flask.sel(nsta=idxs_uniq_station_ids, nobs=idxs_stations_flask)
    # print(ds_ojac_flask)
    msg = f"station_id values after filtering selected stations " \
          f"-->{ds_ojac_flask['station_id'].values}<--, " \
          f"yields shape {ds_ojac_flask.obs_jacobian.shape}"
    logger.debug(msg)

    #
    #-- some consistency checks
    #
    ojac_cont  = ds_ojac_cont.obs_jacobian
    ojac_flask = ds_ojac_flask.obs_jacobian
    assert ojac_cont.dims==('nobs','ng')
    assert ojac_flask.dims==ojac_cont.dims
    assert ojac_cont.shape[1]==ojac_flask.shape[1]
    assert ojac_cont.units=="ppb/(kgCH4/cell)"
    assert ojac_flask.units==ojac_cont.units
    nobs_cont, ng = ojac_cont.shape
    nobs_flask, _ng = ojac_flask.shape
    
    dsmerged_table = OrderedDict()
    for v in ds_ojac_cont.variables:
        dims = ds_ojac_cont[v].dims
        ds_cont = ds_ojac_cont[v].values
        ds_flask = ds_ojac_flask[v].values
        if dims==('nsta',):
            print(f"-->{v}<-- depdends on nsta only")
            dsmerged_table[v] = (np.hstack((ds_cont,ds_flask)), dims)
        elif dims==('nobs',):
            print(f"-->{v}<-- depdends on nobs only ({type(v)})")
            print(type(ds_cont), type(ds_flask))
            print(ds_cont.shape, ds_flask.shape)
            dsmerged_table[v] = (np.hstack((ds_cont,ds_flask)), dims)
        elif dims==('ng',):
            print(f"-->{v}<-- depdends on ng only")
            dsmerged_table[v] = (ds_cont,dims) #-- 1D grid vector taken from
        elif v=='obs_jacobian':
            #-- stack Jacobians along first ('nobs') axis
            ojac_out = np.concat((ds_cont,ojac_flask), axis=0)
            dsmerged_table[v] = (ojac_out, dims)
        else:
            raise RuntimeError(f"-->{v}<-- unexpected dimensions ==>{dims}<== or variable ==>{v}<==")
    #
    #--
    #
    nobsout = nobs_cont + nobs_flask
    ngout   = ng
    nstaout = dsmerged_table['station_id'][0].size
    #
    #--
    #
    if len(stations_flask)<=4:
        flask_station_tag = 'flask-stations-' + '--'.join(stations_flask)
    else:
        flask_station_tag = f'{len(stations_flask)}-flask-stations'
    outname_tokens = [filepath_ojac_cont.stem, '---merged---', filepath_ojac_flask.stem, flask_station_tag]
    outname = '_'.join(outname_tokens) + '.nc'
    outname = set_outname(args, outname)
    msg = f"writing inversion inputs to file ***{outname}***..."
    logger.info(msg)
    fp = Dataset(outname, 'w')
    fp.createDimension('ng', ngout)
    fp.createDimension('nobs', nobsout)
    fp.createDimension('nsta', nstaout)
    for v,v_info in dsmerged_table.items():
        msg = f"@{v}, start with output..."
        logger.debug(msg)
        v_data, v_dims = v_info
        if v in ['station_id','station','region','obstime']:
            v_dtype = str
            ncvar = fp.createVariable(v, v_dtype, v_dims)
        else:
            v_dtype = 'f8'
            ncvar = fp.createVariable(v, v_dtype, v_dims, compression='zlib', complevel=complevel)
        ncvar[:] = v_data[:]
        attr_table = ds_ojac_cont[v].attrs
        for k,v in ds_ojac_cont[v].attrs.items():
            ncvar.setncattr(k, v)
    #
    #-- global attributes
    #
    fp.description = f"Observational Jacobian merged from continuous and flask observation contributions prepared for use within FIT-IC inversion environment"
    fp.filepath_ojac_cont = str(filepath_ojac_cont.absolute())
    fp.filepath_ojac_cont_sha512 = create_sha512(str(filepath_ojac_cont))
    fp.filepath_ojac_flask = str(filepath_ojac_flask.absolute())
    fp.filepath_ojac_flask_sha512 = create_sha512(str(filepath_ojac_flask))
    fp.history = f"{' '.join(sys.argv)}"
    fp.date_created = Timestamp.utcnow().isoformat()
    fp.close()
    msg = f"generated file ***{outname}***"
    logger.info(msg)


def subcmd_inspect_ojac_obs1D(args):
    """
    """
    filepath_ojac =  args.filepath_ojac

    #-- initialise region table
    _init_region_table()
    ds_ojac = xr.open_dataset(filepath_ojac)
    #-- determine regions
    region_uniq = np.unique(ds_ojac.region.values)
    msg = f"detected regions -->{region_uniq}<--"
    logger.info(msg)
    for region,region_info in region_table.items():
        if not region in region_uniq:
            msg = f"expected region -->{region}<-- not present in NetCDF file?"
            raise RuntimeError(msg)
        match region:
            case 'glb600x400':
                #-- check that contribution is zero within child taking into account halo band
                child_info = region_table[region_info.child]
                dlon = region_info.grid.dlon
                dlat = region_info.grid.dlat
                lonmin = child_info.grid.west + dlon
                lonmax = child_info.grid.east - dlon
                latmin = child_info.grid.south + dlat
                latmax = child_info.grid.north - dlat
                xtag = f"{lonstr(lonmin)}-{lonstr(lonmax)}x{latstr(latmin)}-{latstr(latmax)}"
                xcnd = (ds_ojac.region==region) & \
                    (ds_ojac.lon>=lonmin)&(ds_ojac.lon<=lonmax) & \
                    (ds_ojac.lat>=latmin)&(ds_ojac.lat<=latmax)
                xds = ds_ojac.sel(ng=xcnd)
                xojac = xds.obs_jacobian.values
                msg = f"ojac@{region}, restricted to {xtag} yields min/mean/max = " \
                    f"{xojac.min()}/{xojac.mean()}/{xojac.max()}"
                logger.info(msg)
            case 'eur300x200':
                #
                child_info = region_table[region_info.child]
                dlon = region_info.grid.dlon
                dlat = region_info.grid.dlat
                lonmin = child_info.grid.west + dlon
                lonmax = child_info.grid.east - dlon
                latmin = child_info.grid.south + dlat
                latmax = child_info.grid.north - dlat
                xtag = f"{lonstr(lonmin)}-{lonstr(lonmax)}x{latstr(latmin)}-{latstr(latmax)}"
                xcnd = (ds_ojac.region==region) & \
                    (ds_ojac.lon>=lonmin)&(ds_ojac.lon<=lonmax) & \
                    (ds_ojac.lat>=latmin)&(ds_ojac.lat<=latmax)
                xds = ds_ojac.sel(ng=xcnd)
                xojac = xds.obs_jacobian.values
                msg = f"ojac@{region}, restricted to {xtag} yields min/mean/max = " \
                    f"{xojac.min()}/{xojac.mean()}/{xojac.max()}"
                logger.info(msg)
                #-- check that senstivities in halo band are zero
                parent_info = region_table[region_info.parent]
                dlon = parent_info.grid.dlon
                dlat = parent_info.grid.dlat
                lonmin = region_info.grid.west + dlon
                lonmax = region_info.grid.east - dlon
                latmin = region_info.grid.south + dlat
                latmax = region_info.grid.north - dlat
                halo_tag = f"lon<={lonstr(lonmin)} or lon>={lonstr(lonmax)} " \
                    f"or lat<={latstr(latmin)} or lat>={latstr(latmax)}"
                xcnd = (ds_ojac.region==region) & \
                    ( \
                      (ds_ojac.lon<=lonmin) | (ds_ojac.lon>=lonmax) \
                      | (ds_ojac.lat<=latmin) | (ds_ojac.lat>=latmax) \
                      )
                xds = ds_ojac.sel(ng=xcnd)
                xojac = xds.obs_jacobian.values
                msg = f"ojac@{region}, restricted to HALO band ({halo_tag}) yields min/mean/max = " \
                    f"{xojac.min()}/{xojac.mean()}/{xojac.max()}"
                logger.info(msg)
            case 'gns100x100':
                #-- check that senstivities in halo band are zero
                parent_info = region_table[region_info.parent]
                dlon = parent_info.grid.dlon
                dlat = parent_info.grid.dlat
                lonmin = region_info.grid.west + dlon
                lonmax = region_info.grid.east - dlon
                latmin = region_info.grid.south + dlat
                latmax = region_info.grid.north - dlat
                halo_tag = f"lon<={lonstr(lonmin)} or lon>={lonstr(lonmax)} " \
                    f"or lat<={latstr(latmin)} or lat>={latstr(latmax)}"
                xcnd = (ds_ojac.region==region) & \
                    ( \
                      (ds_ojac.lon<=lonmin) | (ds_ojac.lon>=lonmax) \
                      | (ds_ojac.lat<=latmin) | (ds_ojac.lat>=latmax) \
                      )
                xds = ds_ojac.sel(ng=xcnd)
                xojac = xds.obs_jacobian.values
                msg = f"ojac@{region}, restricted to HALO band ({halo_tag}) yields min/mean/max = " \
                    f"{xojac.min()}/{xojac.mean()}/{xojac.max()}"
                logger.info(msg)


def subcmd_compare_ojac_obs1D(args):
    filepath1, filepath2 = args.filepath_ojac
    ds1 = xr.open_dataset(filepath1)
    ds2 = xr.open_dataset(filepath2)
    #
    #-- compare global part
    #
    ds1_glb = ds1.sel(ng=(ds1.region=='glb600x400'))
    ds2_glb = ds2.sel(ng=(ds2.region=='glb600x400'))
    ojac1_glb = ds1_glb.obs_jacobian
    ojac2_glb = ds2_glb.obs_jacobian
    msg = f"shapes of global part of Jacobian ojac1/ojac2 = {ojac1_glb.shape}/{ojac2_glb.shape}"
    logger.info(msg)
    print(np.all(ojac1_glb.values==ojac2_glb.values))
    ojacdif_glb = ojac1_glb.values - ojac2_glb.values
    msg = f"ojac differences (1 minus 2), min/mean/max = " \
        f"{ojacdif_glb.min()}/{ojacdif_glb.mean()}/{ojacdif_glb.max()}"
    logger.info(msg)


def subcmd_mmojac_propagate(args):
    """
    """
    obsjac_filepath = args.obsjac_filepath
    emis_filepath = args.emis_filepath
    figsize = args.__dict__.get('figsize',(10,6))
    dpi     = args.__dict__.get('dpi', 150)

    #
    #--
    #
    fpobs = Dataset('/srv/data/AVENGERS/fit_ic/footprint_obsdir/gns1x1_with-real-obs/tm5-obstable_gns1x1_2021-01-02_obs-fitic.nc')
    obshour = fpobs['/time'][:]
    obsids  = fpobs['/obsid'][:]
    twlength = fpobs['/time_window_length'][:] #-- secnod
    fpobs.close()
    obs_dct = { 'obsid': obsids, 'obshour': obshour, 'twlength_s': twlength}
    obshour_df = DataFrame.from_dict(obs_dct)
    obshour_df = obshour_df.set_index('obsid')
    #
    #--
    #
    emis_tokens = emis_filepath.stem.split('_')
    emis_tag = emis_tokens[0]
    ds_ojac = xr.open_dataset(obsjac_filepath)
    ds_emis = xr.open_dataset(emis_filepath)
    # ds_emis['emis_time'] = [Timestamp(*ds_emis.time.loc[imon,:].values) for imon in ds_emis.nmon]
    # ds_ojac['emis_time'] = [Timestamp(*ds_ojac.emismon_calendar.loc[imon,:].values) for imon in ds_ojac.nemismon]
    # ds_emis = ds_emis.assign(emis_time=[Timestamp(*ds_emis.time.loc[imon,:].values) for imon in ds_emis.nmon])
    # print(ds_emis)
    # sys.exit(0)
    #
    #-- ensure emissions start from same month as emission senstivities
    #
    ojac_emismon = ds_ojac.emismon_calendar.values
    emis_emismon = ds_emis.time.values
    ntc = min(ojac_emismon.shape[1],emis_emismon.shape[1])
    assert np.all(ojac_emismon[0,:ntc]==emis_emismon[0,:ntc])
    #
    #-- ensure consistent dimensions
    #
    nobs,nmon,ng = ds_ojac.obs_jacobian.shape
    _nmon,_ng   = ds_emis.emission.shape
    assert ng==_ng
    assert nmon<=_nmon
    emis = ds_emis.emission.values[:nmon,:]
    emis1D = emis.reshape(ng*nmon)
    #
    #--
    #
    if args.refsimu_filepath!=None:
        refsimu_df = pd.read_csv(args.refsimu_filepath)
        refsimu_df = refsimu_df.set_index('time')
        refsimu_df.index = [Timestamp(_) for _ in refsimu_df.index]
    else:
        refsimu_df = None
    #
    #--
    #
    for sta in ds_ojac.station_id.values:
        cur_ojac = ds_ojac.where(ds_ojac.station==sta, drop=True)
        cur_nobs = len(cur_ojac.nobs.values)
        msg = f"@{sta}, compute simulated concentrations, cur_nobs={cur_nobs}..."
        logger.info(msg)
        iniconc = cur_ojac.iniconc.values
        obs     = cur_ojac.obs.values
        obstime = cur_ojac.obstime.values
        ojac2D = cur_ojac.obs_jacobian.values.reshape(cur_nobs,ng*nmon)
        csimu = np.matmul(ojac2D, emis1D) + iniconc
        # for iobs in range(cur_nobs):
        #     print(f"...{sta}/{obstime[iobs]}: simulated/observed {csimu[iobs]}/{obs[iobs]}")
        obstime = [Timestamp(_) for _ in obstime]
        conc_dict = {'time': obstime,
                     'csimu': csimu,
                     'obs': obs }
        if not refsimu_df is None:
            cnd_reftime = (refsimu_df.index>=obstime[0])&(refsimu_df.index<=obstime[-1])
            refsimu = refsimu_df.loc[cnd_reftime,sta]
            #-- adjust hour of day
            obshour = obshour_df.loc[sta,'obshour']
            refsimu.index = [ (_ + pd.Timedelta(hours=obshour)) for _ in refsimu.index]
            # conc_dict['refsimu'] = refsimu.values
        cur_df = DataFrame.from_dict(conc_dict).set_index('time')
        fig, ax = subplots(1, 1, figsize=figsize)
        cur_df.plot(y='obs', ax=ax, kind='line',
                    color='black', ls='', marker='D', markersize=4, label='observed')
        cur_df.plot(y='csimu', ax=ax, kind='line',
                    color='red', ls='-', marker='+', markersize=4, label='simulated')
        if not refsimu_df is None:
            refsimu.plot(ax=ax, kind='line',
                    color='green', ls='-', marker='x', markersize=4, label='reference-simulation')
        title = f"{sta.upper()}, observed/simulated concentration ({emis_tag} emissions)"
        ax.set_title(title)
        ax.set_ylabel(f"[ppb]")
        ax.grid(which='major', axis='both',color='grey', alpha=0.5)
        trange_tag = f"{obstime[0].strftime('%Y%m%d')}-{obstime[-1].strftime('%Y%m%d')}"
        emis_ftag = f"{emis_tag}-mm-emissions"
        outname_tokens = [sta, 'ch4-concentration', trange_tag, emis_ftag]
        outname = Path('_'.join(outname_tokens) + '.png')
        if args.outdir!=None:
            outname = args.outdir / outname
        outname.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outname, dpi=dpi, bbox_inches='tight')
        plt.close()
        msg = f"...generated ***{str(outname)}***"
        logger.info(msg)
        # print(cur_ojac)
        # print(f"="*30)
    sys.exit(0)
    

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
#--       test_jacobianfwd_1day
#
sparser = subparsers.add_parser('test_jacobianfwd_1day',
                                help="""compute contribution of emissions to selected observation""")
sparser.add_argument('tm5rundir',
                     help="""top-level directory of TM5 adjoint run for footprint creation.""")
sparser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
sparser.add_argument('--obsid',
                     default='cbw_207',
                     help="""restrict Jacobian to one single observational location (default: %(default)s).""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--
#
sparser = subparsers.add_parser('debug_footprint_io',
                                help="""some testing and inspection of the lower level footprint I/O""")
sparser.add_argument('tm5rundir',
                     help="""top-level directory of TM5 adjoint run for footprint creation.""")
sparser.add_argument('--emisday_firstlast',
                     metavar=('tstart','tend'),
                     nargs=2,
                     type=Timestamp,
                     help="""potentially restrict to footprint on certain emission days""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       build_jacobian_period_obs1D
#
sparser = subparsers.add_parser('build_jacobian_period_obs1D',
                                help="""test preparation of inputs for Fortran inversion system.""")
sparser.add_argument('outpath_tm5',
                     type=Path,
                     help="""top-level directory of series of TM5 adjoint runs for footprint creation each of those for one observation day.""")
sparser.add_argument('domain',
                     choices=['gns100x100','glb600x400',],
                     help="""observational domain. As of now the TM5 footprint simulations have been run with two different setups. gns100x100 indicates footprints calculated with the zoom setup (for continuous measurements), while glb600x400 has been used to calculate footprints for flask measurements at only the coarse 6x4 degree resolution globally.""")
sparser.add_argument('--obsid',
                     nargs='+',
                     help="""select one single observational location (default: %(default)s).""")
sparser.add_argument('--obsday_firstlast',
                     nargs=2,
                     type=Timestamp,
                     default=[Timestamp("20210101"),Timestamp("20210102")],
                     help="""first/last day of observational period (default: %(default)s).""")
sparser.add_argument('--emisday_firstlast',
                     nargs=2,
                     type=Timestamp,
                     default=[Timestamp("20210101"),Timestamp("20210102")],
                     help="""first/last day of emissions included for sensitivities (default: %(default)s).""")
sparser.add_argument('--jac4totemis',
                     action='store_true',
                     help="""whether to condense the Jacobian to the sensitivity w.r.t. to the total emissions within the emissions temporal domain.""")
sparser.add_argument('--glb6x4_to_avengers-zoom',
                     action='store_true',
                     help="""Option to be used only in conjunction with domain=='glb600x400'! This will re-distribute the globally computed sensitivities at 6x4 degrees to sensitivities w.r.t. to the grid-cells used for the (gns1x1) zoom  domain.""")
sparser.add_argument('--no-remove-halo','--no-halo_correction',
                     action='store_false',
                     dest='remove_halo',
                     help="""meanwhile, by default the HALO part in zoom domains with parent are being removed. Use this option to activate the old behaviour.""")
sparser.add_argument('--obsdir',
                     type=Path,
                     help="""explicitly provided directory that contains obspack NetCDF data files with CH4 observations for selected station/site (!!!NOTE: this is only required in case the obs table files used for the underlying footprint computations still contain dummy measurement values instead of the real ones.!!!).""")
sparser.add_argument('--difdeg_max',
                     type=float,
                     default=0.00001,
                     help="""geographical coordinates may be varying per observation, in case their differences are below this threshold a single value is associated to only the station in the NetCDF file (default: %(default)s).""")
sparser.add_argument('--difalt_max',
                     type=float,
                     default=0.1,
                     help="""measurement altitude may be varying per observation, in case their differences are below this threshold a single value is associated to only the station in the NetCDF file (default: %(default)s[m]).""")
sparser.add_argument('--refdir',
                     type=Path,
                     help="""TM5 forward simulation, can be used to verify the Jacobian approach.""") 
sparser.add_argument('--outdir',
                     type=Path,
                     help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                     help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       monthly_obsjacobian_for_inversion
#
sparser = subparsers.add_parser('monthly_obsjacobian_for_inversion',
                                help="""generation of an observational Jacobian (and provision of daily inicial concentrations and observations) quantifying the sensitivity of (daily) observations at selected sites with respect to monthly emisssions.""")
sparser.add_argument('outpath_tm5',
                     type=Path,
                     help="""top-level directory of series of TM5 adjoint runs for footprint creation each of those for one observation day.""")
sparser.add_argument('domain',
                     choices=['gns100x100','glb600x400',],
                     help="""observational domain. As of now the TM5 footprint simulations have been run with two different setups. gns100x100 indicates footprints calculated with the zoom setup (for continuous measurements), while glb600x400 has been used to calculate footprints for flask measurements at only the coarse 6x4 degree resolution globally.""")
sparser.add_argument('--obsid',
                     nargs='+',
                     help="""select one single observational location (default: %(default)s).""")
sparser.add_argument('--obsmon_firstlast',
                     nargs=2,
                     type=Timestamp,
                     default=[Timestamp('2021-01-01'),Timestamp('2021-02-28')],
                     help="""selected temporal range of observation months (default: %(default)s), from the specied timestamps only the year/month will actually been used.""")
sparser.add_argument('--emisday_first',
                     type=Timestamp,
                     default=[Timestamp("20210101")],
                     help="""first day of emissions included for sensitivities, must start on day 1 of a month (default: %(default)s).""")
sparser.add_argument('--no-glb6x4_to_avengers-zoom',
                     dest='glb6x4_to_avengers_zoom',
                     action='store_false',
                     help="""Option to be used only in conjunction with domain=='glb600x400'! This will re-distribute the globally computed sensitivities at 6x4 degrees to sensitivities w.r.t. to the grid-cells used for the (gns1x1) zoom  domain.""")
sparser.add_argument('--no-remove-halo','--no-halo_correction',
                     action='store_false',
                     dest='remove_halo',
                     help="""meanwhile, by default the HALO part in zoom domains with parent are being removed. Use this option to activate the old behaviour.""")
sparser.add_argument('--obsdir',
                     type=Path,
                     help="""explicitly provided directory that contains obspack NetCDF data files with CH4 observations for selected station/site (!!!NOTE: this is only required in case the obs table files used for the underlying footprint computations still contain dummy measurement values instead of the real ones.!!!).""")
sparser.add_argument('--difdeg_max',
                     type=float,
                     default=0.00001,
                     help="""geographical coordinates may be varying per observation, in case their differences are below this threshold a single value is associated to only the station in the NetCDF file (default: %(default)s).""")
sparser.add_argument('--difalt_max',
                     type=float,
                     default=0.1,
                     help="""measurement altitude may be varying per observation, in case their differences are below this threshold a single value is associated to only the station in the NetCDF file (default: %(default)s[m]).""")
sparser.add_argument('--refdir',
                     type=Path,
                     help="""TM5 forward simulation, can be used to verify the Jacobian approach.""") 
sparser.add_argument('--outdir',
                     type=Path,
                     help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                     help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       concat_ojac_obs1D
#
sparser = subparsers.add_parser('concat_ojac_obs1D',
                                help="""combine observational Jacobians prepared for continuous measurements and flask measurements.""")
sparser.add_argument('filepath_ojac_cont',
                     type=Path,
                     help="""NetCDF file prepared for continous observations (it is assumed that it was built for grid-cells based on the AVENGERS 3-level zoom).""")
sparser.add_argument('filepath_ojac_flask',
                     type=Path,
                     help="""NetCDF file prepared for flask observations (it is assumed that it was built for only the grid-cells in the global glb600x400 domain).""")
sparser.add_argument('--stations_flask',
                     nargs='+',
                     default=['asc_90','brw_16','cgo_164','izo_2377','mhd_26','mlo_3437','spo_2815',],
                     help="""restrict to selected flask stations (default: %(default)s).""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

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
# sparser.add_argument('--month',
#                      type=int,
#                      choices=list(np.arange(1,13)),
#                      help="""restrict to single month.""")
sparser.add_argument('--regions',
                     nargs='+',
                     choices=['glb600x400','eur300x200','gns100x100',],
                     default=['glb600x400','eur300x200','gns100x100',],
                     help="""selected regions (default: %(default)s), better only change for test purposes.""")
sparser.add_argument('--no-remove-halo','--no-halo_correction',
                     action='store_false',
                     dest='remove_halo',
                     help="""meanwhile, by default the HALO part in zoom domains with parent are being removed. Use this option to activate the old behaviour.""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       create_target_jacobian
#
sparser = subparsers.add_parser('create_target_jacobian',
                                help="""preparation of dedicated target Jacobian for Fortran inversion system.""")
sparser.add_argument('--no-remove-halo','--no-halo_correction',
                     action='store_false',
                     dest='remove_halo',
                     help="""meanwhile, by default the HALO part in zoom domains with parent are being removed. Use this option to activate the old behaviour.""")
sparser.add_argument('--countryfrct_filepath',
                     type=Path,
                     help="""provide 1degree gridded NetCDF file with country fractions.""")
sparser.add_argument('--countries',
                     nargs='+',
                     help="""select list of countries (ISO3 identifier) which must be fully covered within the FIT-IC innermost zoom domain.""")
sparser.add_argument('--emission_filepath',
                     type=Path,
                     help="""for debugging: provide a compliant emissions NetCDF file to check the target values.""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       inspect_ojac_obs1D
#
sparser = subparsers.add_parser('inspect_ojac_obs1D',
                                help="""some inspection and consistency check on generated observational Jacobian (should be applicable to both, "global/flask" and "zoomed/continuous)""")
sparser.add_argument('filepath_ojac',
                     type=Path,
                     help="""NetCDF file prepared for observational Jacobian""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       compare_ojac_obs1D
#
sparser = subparsers.add_parser('compare_ojac_obs1D',
                                help="""some comparison of generated observational Jacobian (primarily to trace potential issues with Jacobians when HALOs are removed).""")
sparser.add_argument('filepath_ojac',
                     nargs=2,
                     type=Path,
                     help="""NetCDF file prepared for observational Jacobian""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       mmojac_propagate
#
sparser = subparsers.add_parser('mmojac_propagate',
                                help="""compute concentrations by matrix multiplication of observational Jacobian with a monthly emission field.""")
sparser.add_argument('obsjac_filepath',
                     type=Path,
                     help="""NetCDF file providing observational Jacobian as well as observed and inicial concentrations.""")
sparser.add_argument('emis_filepath',
                     type=Path,
                     help="""NetCDF file providing monthly emissions (must cover the range within the observational Jacobian and be on the same grid.""")
sparser.add_argument('--refsimu_filepath',
                     type=Path,
                     help="""csv file providing fit-ic averaged daily concentrations simulated with TM5 in forward mode.""")
sparser.add_argument('--outdir',
                     type=Path,
                     help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")




################################################################################
#
#                   p r o g r a m   s t a r t
#
def main(args):

    ts = Timestamp.utcnow()
    
    if args.subcmds=='test_jacobianfwd_1day':
        subcmd_test_jacobianfwd_1day(args)

    if args.subcmds=='debug_footprint_io':
        subcmd_debug_footprint_io(args)

    if args.subcmds=='build_jacobian_period_obs1D':
        subcmd_build_jacobian_period_obs1D(args)

    if args.subcmds=='monthly_obsjacobian_for_inversion':
        subcmd_monthly_obsjacobian_for_inversion(args)

    if args.subcmds=='monthly_emissions_for_inversion':
        subcmd_monthly_emissions_for_inversion(args)

    if args.subcmds=='create_target_jacobian':
        subcmd_create_target_jacobian(args)

    if args.subcmds=='concat_ojac_obs1D':
        subcmd_concat_ojac_obs1D(args)

    if args.subcmds=='inspect_ojac_obs1D':
        subcmd_inspect_ojac_obs1D(args)

    if args.subcmds=='compare_ojac_obs1D':
        subcmd_compare_ojac_obs1D(args)

    if args.subcmds=='mmojac_propagate':
        subcmd_mmojac_propagate(args)

    #
    te = Timestamp.utcnow()
    msg = f"...subcommand +++{args.subcmds}+++ DONE (time_elapsed={te-ts})"
    logger.info(msg)
#
if __name__ == '__main__':
    import datetime as dtm

    progname = os.path.basename(__file__)

    #-----------------------------
    #          P R O G R A M   S T A R T
    #
    ttstart = Timestamp.utcnow().isoformat()
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
