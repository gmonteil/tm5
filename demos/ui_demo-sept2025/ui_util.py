import os
import sys
import xarray as xr
import pandas as pd
import datetime as dtm
import netCDF4 as nc4
import numpy as np
from pathlib import Path
from collections import OrderedDict
from loguru import logger
from typing import Union, List, Tuple

#
# WGS84 ellipsoid:
# https://en.wikipedia.org/wiki/Earth_radius
# Equatorial radius: a = (6378.1370 km)
# Polar radius:      b = (6356.7523 km)
_A = 6378.1370 * 1000
_B = 6356.7523 * 1000
# In geophysics, the International Union of Geodesy and Geophysics (IUGG) defines the Earth's mean radius (denoted R1) to be (2a + b)/3
#
_ERmeter  = (2*_A + _B)/3


def get_hostname():
    import socket
    hostname = socket.gethostname()
    return hostname

def is_jupyterhub():
    import os
    return 'JUPYTERHUB_API_TOKEN' in os.environ

# Fix env variables for cartopy:
def fix_env() -> None:
    import sys, os
    from pathlib import Path
    env_base_path = Path(sys.executable).parents[1]
    os.environ["SSL_CERT_FILE"] = str(env_base_path / 'ssl' / 'cert.pem')
    os.environ["SSL_CERT_DIR"] = str(env_base_path / 'ssl' / 'certs')
    os.environ["REQUESTS_CA_BUNDLE"] = str(env_base_path / 'ssl' / 'cert.pem')
    os.environ["PROJ_LIB"] = str(env_base_path / 'share' / 'proj')
    if get_hostname().find('cosmos')>=0:
        pass
    if is_jupyterhub():
        import cartopy
        cartopy.config['data_dir'] = str(env_base_path/'share'/'cartopy')

#---------------------------------------
#
#         h a r d - c o d e d   p a t h s
#
#---------------------------------------
#
#-- file/directory selection depends on host
#   - COSMOS
#   - ICOS Jupyter Hub
#
#
#
#       COSMOS
#
#
if get_hostname().find('cosmos')>=0:
    #
    #-- top-level directory with pre-computed TM5 simulations
    #
    demo_precompdir = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/runs_with-convdiff_repeated')
    if not demo_precompdir.exists():
        msg = f"directory with pre-computed simulations not found -->{demo_precompdir}<--"
        raise RuntimeError(msg)
    #
    #--  CAMS/obspack concentrations
    #
    camsfile = '/lunarc/nobackup/projects/ghg_inv/michael/CAMS/ch4/cams_ch4conc_at-obspack-locations_2021.nc'
    obspackdir = '/lunarc/nobackup/projects/ghg_inv/michael/FIT-IC/obspack_ch4_1_GLOBALVIEWplus_v6.0_2023-12-01/data/nc'
    logger.remove()
    logger.add(sys.stdout, level="DEBUG")
#
#
#       ICOS jupyter lab
#
#
elif is_jupyterhub(): #-- ICOS jupyter lab
    #
    #-- top-level directory with pre-computed TM5 simulations
    #
    demo_precompdir = Path('/data/avengers/precomputed_output_sept2025')
    if not demo_precompdir.exists():
        msg = f"directory with pre-computed simulations not found -->{demo_precompdir}<--"
        raise RuntimeError(msg)
    #
    #--  CAMS/obspack concentrations
    #
    camsfile   = '/data/avengers/fit_ic/validation/cams_ch4conc_at-obspack-locations_2021.nc'
    obspackdir = '/data/avengers/fit_ic/validation/obspack_ch4_1_GLOBALVIEWplus_v6.0_2023-12-01/data/nc'
    #-- no loguru logging on ICOS
    logger.remove()
    logger.add(sys.stdout, level="WARNING")
#
#
#       Laptop
#
#
elif get_hostname().find('mvobook2')>=0:
    #
    #-- top-level directory with pre-computed TM5 simulations
    #
    demo_precompdir = Path('/srv/data/avengers/geneva2025/fit-ic_precomputed-output')
    if not demo_precompdir.exists():
        msg = f"directory with pre-computed simulations not found -->{demo_precompdir}<--"
        raise RuntimeError(msg)
    camsfile = '/data/avengers/fit_ic/validation/cams_ch4conc_at-obspack-locations_2021.nc'
    obspackdir = '/data/avengers/fit_ic/validation/obspack_ch4_1_GLOBALVIEWplus_v6.0_2023-12-01/data/nc'
else:
    msg = f"detected platform not yet supported"
    raise RuntimeError(msg)


def get_precomputed_expdirs() -> OrderedDict:
    """preliminary ad-hoc routine that builds a dictionary
    of tags (identifying the experiment) to the output directory
    of the associated TM5 simulation.

    The current listing looks like this (Sep 8, 2025):
    ------------------------------------
    fitic-simu-default_oh-cams
    fitic-simu-edgarflat_oh-cams
    fitic-simu-regional_anthro-no-agri_oh-cams
    fitic-simu-regional_anthro-no-fossil_oh-cams
    fitic-simu-regional_anthro-no-france_oh-cams
    fitic-simu-regional_anthro-no-netherlands_oh-cams
    fitic-simu-regional_anthro-no-waste_oh-cams
    fitic-simu-regional_oh-cams

    Extended by further experiments (Sept 18, 2025)
    fitic-simu-half-oh_oh-cams
    fitic-simu-no-germany_oh-cams
    fitic-simu-no-northamerica_oh-cams
    """
    exp_dict = OrderedDict()
    #-- first level:
    exp_dirs = sorted(demo_precompdir.glob('fitic-simu-*'))
    for _d in exp_dirs:
        # print(f"-->{_d}<--")
        exp_outdir = _d / 'output_2021-01-01--2022-01-01'
        if not exp_outdir.is_dir():
            #-- expected simulation does not exist
            msg = f"expected TM5 simulation output directory " \
                f"***{str(exp_outdir)}*** does not exist, skipping..."
            logger.warn(msg)
            continue
        _basename = _d.stem.replace('fitic-simu-','').replace('_oh-cams','')
        # print(f"-->{_basename}<--")
        if _basename.startswith('regional_'):
            if _basename.find('anthro-no'):
                exp_tag = 'regional,' + \
                    _basename.replace('regional_','').replace('anthro-','')
                exp_tag = 'regional_' + \
                    _basename.replace('regional_','').replace('anthro-','')
                is_country = _basename.find('france')>=0 or _basename.find('netherlands')>=0
                if is_country:
                    exp_tag = _basename
                else:
                    exp_tag = _basename.replace('anthro-','')
            else:
                exp_tag = 'regional'
        else:
            exp_tag = _basename
        exp_dict[exp_tag] = exp_outdir
    # #-- DEBUG only
    # for k,v in exp_dict.items():
    #     print(f"{k} --> {v}")
    return exp_dict


def fitic_inputemisdir() -> str:
    if get_hostname().find('cosmos')>=0:
        emisdir = '/lunarc/nobackup/projects/ghg_inv/michael/TM5/input/ch4/emissions_fitic'
        emisdir = '/lunarc/nobackup/projects/ghg_inv/michael/TM5/input/ch4/emissions'
    else:
        #-- on ICOS jupyter lab
        emisdir = '/project/fit_ic/data/input/emissions/CH4'
        #-- 2025-04-14: switchted to directory path as used on VM
        #               (and mounted as such on the Jupyter-Hub)
        emisdir = '/data/avengers/fit_ic/input/emissions'

    emisdir = Path(emisdir)
    if not emisdir.is_dir():
        msg = f"fit-ic input emissions directory not found on system " \
            f"***{emisdir}***"
        raise RuntimeError(msg)

    return emisdir

#===========================================================
#
#                   s h o u l d   b e   m o v e d   t o   m o d u l e ....
#
#===========================================================
def sphere_grid_find_close(lonq : float, latq : float,
                           longrd : np.ndarray, latgrd : np.ndarray,
                           distm_max : float = None, radius_m : float = _ERmeter) -> tuple:
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
        if len(longrd.shape)==1:
            idxs_min = idxs_min[0]
    dist_min = dm[idxs_min]

    return (idxs_min, dist_min)


def cams_at_obspack_load_conctseries( camsfile : str,
                                      lonq : float,
                                      latq : float,
                                      altq : float ) -> pd.DataFrame:
    """Read concentration time-series from a dedicated NetCDF file
    providing CAMS concentrations extracted at a list of obspack station
    locations.
    Concentrations are extracted for specific spatial coordinates,
    including altitude (to be provided in [m]).
    """
    cams_path = Path(camsfile)
    if not cams_path.exists():
        msg = f"...provided CAMS file ***{camsfile}*** not accessible."
        raise IOError(msg)
    cams_xds = xr.open_dataset(camsfile)
    #
    #--
    #
    try:
        cams_release = cams_xds.attrs['release']
    except KeyError:
        cams_release = 'rev. unknown'
    #
    #-- spatial coordinate values
    #
    cams_lon = cams_xds.longitude.values
    cams_lat = cams_xds.latitude.values
    cams_date_lst = cams_xds.time.values

    #
    #-- detect index in file closest to lonq/latq
    #
    icams, dstcams = sphere_grid_find_close(lonq, latq, cams_lon, cams_lat)
    #
    #-- unfortunately file is prepared such that
    #   there can be multiple entries at the same location
    #   (for obspack stations with different vertical levels)
    #   - we can select the first index here
    #
    if len(icams)==0:
        msg = f"...no proper location found for lonq/latq = {lonq}/{latq} " \
            "in file ***{camsfile}***"
        logger.warning(msg)
        cams_xds.close()
        return None
    #
    #-- select this location
    #
    icams = icams[0]
    cams_loc = cams_xds.sel(obspack_location=icams)
    #
    #-- altitude is time-dependent!
    #
    # print(f"-->{cams_loc.altitude.coords}<--")
    cams_alt = cams_loc.altitude.data
    nt,nlevalt = cams_alt.shape
    req_alt = np.array([altq,])
    req_alt = req_alt[np.newaxis,:]
    #
    #-- altitude difference  query-cams
    #
    dif_alt = req_alt - cams_alt
    dif_alt[dif_alt<0] = 1000000 #-- large value
    ilev = np.argmin(dif_alt, axis=1)
    itime = np.arange(nt)
    #
    #--
    #
    cams_conc = cams_loc.CH4.data
    _,nlevconc = cams_conc.shape
    assert nlevconc==nlevalt-1
    cams_conc = cams_conc[itime,ilev]
    cams_tag = f"cams-inversion ({cams_release})"
    cams_df = pd.DataFrame.from_dict({'time':cams_date_lst,
                                      cams_tag:cams_conc})
    #-- dispose resources
    cams_xds.close()
    return cams_df


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


if __name__ == '__main__':

    obspackdir = '/srv/data/avengers/fit_ic/validation/obspack_ch4_1_GLOBALVIEWplus_v6.0_2023-12-01/data/nc'
    # // group attributes:
    #               :name = "South Pole, Antarctica" ;
    #               :latitude = -89.98 ;
    #               :longitude = -24.8 ;
    #               :altitude = 2821. ;
    #               :region = "glb600x400" ;
    #               :abbr = "FM/spo_11" ;
    sta_tag = 'spo_11'
    lonq = -24.8
    latq = -89.98
    altq = 2821.
    tstart = pd.Timestamp('20210101')
    tend   = pd.Timestamp('20211231')
    ch4_dict = obspack_load_conctseries(obspackdir, tstart, tend,
                                        sta_tag, lonq, latq, altq)
