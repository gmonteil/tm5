#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""simple interface(s) to read CAMS posterior concentrations at obspack locations. The hat were extracted at 1x1 degree ."""

#-- system packages
import os
import sys
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
from loguru import logger
from typing import Union, List, Tuple
#-- local packages
from tm5.post.utilities import sphere_grid_find_close


def camsdir_load_concentrations( cams_dir : Union[str, Path],
                                 year     : int = 2021 ) -> xr.Dataset:
    """Function that returns xarray dataset derived
    from concatenating monthly CAMS CH4 concentration files
    along the time-dimensions.
    """
    pattern = f"cams73_latest_ch4_conc_surface_inst_{year}??.nc"
    file_lst = sorted(Path(cams_dir).glob(pattern))
    # print(file_lst)
    msg = f"loading CAMS concentrations from {len(file_lst)} input files..."
    logger.info(msg)
    data = xr.open_mfdataset(file_lst)
    return data


def cams_at_obspack_load_conctseries( camsfile : str,
                                      lonq : float,
                                      latq : float,
                                      altq : float ) -> pd.DataFrame:
    """Read concentration time-series from a dedicated NetCDF file
    providing CAMS (posterior) concentrations extracted from a global
    CAMS inversion (at 1x1 degree) from where the concentrations at a list
    of obspack station locations were picked from the grid-cells of the
    station location.
    """
    cams_path = Path(camsfile)
    if not cams_path.exists():
        msg = f"...provided CAMS file ***{camsfile}*** not accessible."
        raise IOError(msg)
    cams_xds = xr.open_dataset(camsfile)
    cams_lon = cams_xds.longitude.values
    cams_lat = cams_xds.latitude.values
    cams_date_lst = cams_xds.time.values

    #
    #-- detect index in file closest to lonq/latq
    #
    icams, dstcams = sphere_grid_find_close(lonq, latq, cams_lon, cams_lat)
    #-- get rid of tuple (cams_lon/cams_lat have single dimension only),
    #   indices are single tuple element here
    icams = icams[0]
    msg = f"...lon/lat={lonq}/{latq} detected at index={icams} in data file"
    logger.info(msg)
    #
    #-- unfortunately file is prepared such that
    #   there can be multiple entries at the same location
    #   (for obspack stations with different vertical levels)
    #   - horizontally all are at the same location and we can
    #     pick the first index.
    #
    if len(icams)==0:
        msg = f"...no proper location found for lonq/latq = {lonq}/{latq} " \
            "in file ***{camsfile}***"
        logger.warning(msg)
        cams_xds.close()
        return None
    else:
        icams = icams[0] #-- first index!
    #
    #-- CAMS concentration closest to requested coordinates
    #
    cams_loc = cams_xds.sel(obspack_location=icams)
    #
    #-- need altitude and concentration
    #
    cams_alt = cams_loc.altitude
    cams_conc = cams_loc.CH4
    #
    #-- altitude provides the layer-to-layer values
    #   (i.e. should be on the vertical box boundaries)
    #   and also depends on time
    #
    assert cams_alt.dims==('time','hlevel')
    assert cams_conc.dims==('time','level')
    assert cams_alt.shape[1]==cams_conc.shape[1]+1 #-- one more
    nt,nlev = cams_conc.shape
    #
    #-- altitude/concentration as numpy arrays
    #
    cams_alt  = cams_alt.data
    cams_conc = cams_conc.data
    msg = f"...requested altq={altq}, average cams_alt[0] = {cams_alt[0,:].mean()}"
    logger.debug(msg)
    #
    #-- determine level (and time) indices with
    #   cams_lower<=altq<=cams_upper
    #
    cnd_alt = (altq>=cams_alt[:,:nlev])&(altq<=cams_alt[:,1:nlev+1])
    iit,iilev = np.where(cnd_alt)
    #
    #-- ATTENTION: very likely the level index will not change over time...
    #
    # assert np.all(np.diff(iilev)==0)
    if len(iit)!=nt:
        msg = f"...requested altitude altq={altq} did not match with properly " \
            f"with altitude of CAMS concentrations. CAMS will be discarded."
        logger.warning(msg)
        cams_xds.close()
        return None
    #
    #--
    #
    cams_conc = cams_conc[iit,iilev] #-- yields 1D array
    cams_df = pd.DataFrame.from_dict({'time':cams_date_lst,
                                      'cams':cams_conc})
    #-- dispose resources
    cams_xds.close()
    return cams_df
