#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface to read obspack CH4 concentrations from NetCDF file(s)."""

#-- system packages
import os
import sys
from pathlib import Path
import numpy as np
import netCDF4 as nc4
import xarray as xr
import pandas as pd
from loguru import logger
from typing import Union, List, Tuple

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
