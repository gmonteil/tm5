import os
import sys
import param
import xarray as xr
import pandas as pd
import datetime as dtm
from pandas import date_range, DatetimeIndex, Timestamp, DataFrame
from functools import lru_cache
import  netCDF4 as nc4
import numpy as np
from pathlib import Path
from collections import OrderedDict
import holoviews as hv
import hvplot.xarray
import panel as pn
import hvplot.pandas
from loguru import logger
from typing import Union, List, Tuple
#-- to avoid PROJ path warning message on JupyterHub,
#   fix_env would need to be called prior import of ccrs.
#   For now, this import of crs is made only in
#   routines where needed!
# import cartopy.crs as ccrs

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

# #
# #-- fix to single output directory
# #
# outdir = '/project/fit_ic/data/output/guillaume_30jan' #-- initial output directory (ICOS Jupyter Hub only)

#
#-- FIT-IC simulations limited to CH4 species, only
#
species = 'CH4'

#
#-- file/directory selection depends on host (COSMOS or ICOS Jupyter Hub)
#
if get_hostname().find('cosmos')>=0:
    outdir = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/testruns_output-cosmos-home_until-2025-03-14/testruns.v5/tm5simu-tropo34-avengers-1_meteo-coarsened-True_Makefile.singularity.ifort_platform-cx03/output_2021-01-01--2022-01-01')
    
    camsfile = '/lunarc/nobackup/projects/ghg_inv/michael/CAMS/ch4/cams_ch4conc_at-obspack-locations_2021.nc'
    obspackdir = '/lunarc/nobackup/projects/ghg_inv/michael/FIT-IC/obspack_ch4_1_GLOBALVIEWplus_v6.0_2023-12-01/data/nc'
    logger.remove()
    logger.add(sys.stdout, level="DEBUG")
else: #-- ICOS jupyter lab
    outdir = Path('/project/fit_ic/data/output_misc/output_2021-01-01--2022-01-01_avengers-1_singletracer_all-emis-default')
    outdir = Path('/data/avengers/ga2025/fit-ic_precomputed-output/fitic-simu-default_20210101--20211231')
    camsfile = '/project/fit_ic/data/validation/cams_ch4conc_at-obspack-locations_2021.nc'
    obspackdir = '/project/fit_ic/data/validation/obspack_ch4_1_GLOBALVIEWplus_v6.0_2023-12-01/data/nc'
    #-- no loguru logging on ICOS
    logger.remove()
    logger.add(sys.stdout, level="WARNING")

emisdir       = outdir / 'emissions'
stations_file = outdir / 'stations/stations.nc4'

if not stations_file.exists():
    msg = f"file ***{stations_file}*** not accessbible"
    raise RuntimeError(msg)


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
                              altq : float) -> Union[None,pd.DataFrame]:
    #
    #-- station tag is like: [id]_[alt]
    #   (with alt likely level above ground
    #
    dfobs = None
    id_tag,alt_tag = sta_tag.split('_')
    fptn = f'ch4_{id_tag}*.nc'
    file_lst = list(Path(obspackdir).glob(fptn))
    if len(file_lst)==0:
        msg = f"...@{sta_tag} no matching obspack files found."
        logger.warning(msg)
        return dfobs
    else:
        # msg = f"...@{sta_tag}, {len(file_lst)} candidate files " \
        #     f"(***{[_.name for _ in file_lst]}***)"
        # logger.info(msg)
        obspack_dict = None
        obspack_fpath = None
        for ipath,fpath in enumerate(file_lst):
            ch4_dict = load_obspack_ch4(fpath,
                                        date_first=tcover_start,
                                        date_last=tcover_end)
            if ch4_dict==None:
                continue
            elev = ch4_dict.get('elevation',None)
            alt  = ch4_dict.get('altitude',None) #-- that should be the measurement height
            # print(f"MVODEBUG::@{sta_tag} (height={altq}), ***{fpath.name}*** elev ==>{elev}<== alt ==>{alt}<==")
            if alt is None:
                continue #-- Hmm, obspack data file without measurement height...
            elif len(alt)==0:
                continue #-- Hmm, obspack data file without measurement height...
            alt_dif = np.max(np.abs(altq-alt))
            if alt_dif<=1: #-- less than 1meter difference
                obspack_fpath = fpath
                obspack_dict = ch4_dict
                msg = f"...@{sta_tag} (height={altq}), selected file ***{obspack_fpath}*** (alt_dif={alt_dif}[m])"
                logger.info(msg)
                break
        if obspack_fpath==None:
            msg = f"...@{sta_tag} (height={altq}), no matching obspack observation found."
            logger.warning(msg)
        else:
            dfobs = pd.DataFrame.from_dict({'time':obspack_dict['time'],
                                            'obspack_ch4':obspack_dict['ch4']})
    return dfobs


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



# Fix env variables for cartopy:
def fix_env() -> None:
    import sys, os
    from pathlib import Path
    env_base_path = Path(sys.executable).parents[1]
    os.environ["SSL_CERT_FILE"] = str(env_base_path / 'ssl' / 'cert.pem')
    os.environ["SSL_CERT_DIR"] = str(env_base_path / 'ssl' / 'certs')
    os.environ["REQUESTS_CA_BUNDLE"] = str(env_base_path / 'ssl' / 'cert.pem')
    os.environ["PROJ_LIB"] = str(env_base_path / 'share' / 'proj')

@lru_cache
def get_file_list(path: Path, pattern: str) -> List[Path]:
    return list(Path(path).glob(pattern))
        

@lru_cache
def load_emis(pattern: str) -> xr.Dataset:
    return xr.open_mfdataset(pattern, concat_dim='time', combine='nested')


class EmissionExplorer(pn.viewable.Viewer):
    time_index = param.Integer()
    data = param.ObjectSelector()
    category = param.Selector()
    region   = param.Selector()
    
    def __init__(self, settings):
        super().__init__()
        fix_env()
        self.settings = settings
#        emis_ptn = f"{str(emisdir)}/ch4emis.CH4.{region}.*.nc"
        # emis_ptn = f"{str(emisdir)}/ch4emis.CH4.glb600x400.*.nc"
        # self.data = load_emis(emis_ptn)
        # self.param.category.objects = list(self.data.data_vars)
        # self.category = self.param.category.objects[0]
        self.setup_emis()
        #-- self.setup_emis_mp() #-- multiprocessor version to speed-up
        #
        #--
        #
        msg = f"{dtm.datetime.utcnow()}, @__init__, setting self.param.region.objects"
        # print(f"DEBUG::{msg}")
        self.param.region.objects = list(self.data.keys())
        self.region = self.param.region.objects[0]
        msg = f"{dtm.datetime.utcnow()}, @__init__, setting self.param.category.objects"
        # print(f"DEBUG::{msg}")
        self.param.category.objects = list(self.data[self.region].data_vars)
        self.category = self.param.category.objects[0]

        # dates = {Timestamp(v).strftime('%B %Y'): iv for (iv, v) in enumerate(self.data.time.values)}
        self.widgets = {
            'date_selector': pn.widgets.DiscreteSlider.from_param(self.param.time_index, options=self.dates, name=''),
            'field_selector': pn.widgets.Select.from_param(self.param.category),
            'region_selector': pn.widgets.Select.from_param(self.param.region)
        }
        
    def __panel__(self):
        msg = f"{dtm.datetime.utcnow()}, start@__panel__ with self.region={self.region}"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        widget_regsel = self.widgets['region_selector']
        widget_fldsel = self.widgets['field_selector']
        widget_datesel = self.widgets['date_selector']
        # widget_map1 = hv.DynamicMap(pn.bind(self.map_emis,widget_datesel,widget_fldsel,widget_regsel))
        # widget_map2 = hv.DynamicMap(pn.bind(self.plot_domain_emis,widget_datesel,widget_fldsel,widget_regsel))

        widget_map1 = hv.DynamicMap(self.map_emis)
        widget_map2 = hv.DynamicMap(self.plot_domain_emis)
        widget = pn.Column(widget_regsel,widget_fldsel,widget_datesel,
                           pn.Row(widget_map1,widget_map2))
        return widget
        # return pn.Column(
        #     self.widgets['region_selector'],
        #     self.widgets['field_selector'],
        #     self.widgets['date_selector'],
        #     pn.Row(
        #         hv.DynamicMap(self.map_emis),
        #         hv.DynamicMap(self.plot_domain_emis)
        #     )
        # )

    #-- probably better drop 'watch=True' (?)
    #   (see https://panel.holoviz.org/how_to/param/dependencies.html)
    @property
    # @pn.depends('time_index', watch=True)
    @pn.depends('time_index')
    def current_date(self):
        return self.data[self.region].time.values[self.time_index]

    @property
    # @pn.depends('region', watch=True)
    @pn.depends('region')
    def current_extent(self):
        if self.region=='glb600x400':
            return [-180, 180, -90, 90]
        elif self.region=='eur300x200':
            return [-36, 54, 22, 74]
        elif self.region=='gns100x100':
            return [0, 18, 42, 58]
        else:
            return RuntimeError(f"unexpected region -->{self.region}<--")
    
    # @pn.depends('time_index', 'region', 'category', watch=True)
    @pn.depends('time_index', 'region', 'category')
    def map_emis(self):
        import cartopy.crs as ccrs
        #-- MVO::can end-up here with self.category==None,
        #        so need to catch this case
        msg = f"{dtm.datetime.utcnow()}, start@map_emis -->{self.region}<-- -->{self.category}<-- itime={self.time_index}"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        cur_cat = self.category if self.category!=None else 'wetland'
        cur_date = pd.to_datetime(str(self.current_date)).strftime('%Y-%m-%d')
        lonmin,lonmax,latmin,latmax = self.current_extent
        title = f"{cur_cat}@{self.region} ([kgCH4/cell/s], {cur_date})"
        clabel = f"[kgCH4/cell/s]"
        # print(f"DEBUG::@map_emis title -->{title}<--")
        # features features (default=None): A list of features or a dictionary of features and the scale at which to render it. Available features include ‘borders’, ‘coastline’, ‘lakes’, ‘land’, ‘ocean’, ‘rivers’ and ‘states’. Available scales include ‘10m’/’50m’/’110m’.
        if self.region=='glb600x400':
            cproj = ccrs.PlateCarree()
            cfeatures = {'borders':'110m', 'coastline':'110m'}
            coastline = '110m'
        elif self.region=='eur300x200':
            cproj = ccrs.GOOGLE_MERCATOR
            cproj = ccrs.PlateCarree()
            # cfeatures = ['borders',]
            cfeatures = {'borders':'50m', 'coastline':'50m'}
            coastline = '50m'
        elif self.region=='gns100x100':
            cproj = ccrs.GOOGLE_MERCATOR
            cproj = ccrs.PlateCarree()
            cfeatures = {'borders':'10m', 'coastline':'10m'}
            coastline = '10m'
        cur_emis = self.data[self.region][cur_cat].isel(time=self.time_index)
        msg = f"{dtm.datetime.utcnow()}, @map_emis, current data ready"
        logger.debug(msg)
        # emis_hvplot = cur_emis.hvplot.quadmesh(geo=True,
        #                                        coastline=coastline,
        #                                        features=cfeatures,
        #                                        xlim=(lonmin,lonmax),
        #                                        ylim=(latmin,latmax),
        #                                        clabel=clabel, title=title)
        emis_hvplot = cur_emis.hvplot.quadmesh(
            xlim=(lonmin,lonmax), ylim=(latmin,latmax),
            coastline=coastline, title=title,
            crs=ccrs.PlateCarree(),
            projection=cproj, features=cfeatures, project=True)
        msg = f"{dtm.datetime.utcnow()}, @map_emis, hvplot ready"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        return emis_hvplot.opts(backend_opts={"plot.toolbar.autohide": True})

    # @pn.depends('time_index', 'region', 'category', watch=True)
    @pn.depends('time_index', 'region', 'category')
    def plot_domain_emis(self):
        msg = f"{dtm.datetime.utcnow()}, start@plot_domain_emis with time_index={self.time_index}, " \
            f"region -->{self.region}<--  category -->{self.category}<--"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        if self.region==None:
            return
        # cat_data = self.data[self.region].sum(('lat', 'lon'))
        # cat_df = cat_data.to_dataframe()
        # cat_df = cat_df.reset_index()
        cat_df = self.glob_timeseries[self.region]
        ylabel = f"[kgCH4/s]"
        title = f"sectorial CH4 emissions (@{self.region})"
        # print(cat_df.head())
        cat_hvplot = cat_df.hvplot(grid=True, x='time', xlabel='time', ylabel=ylabel, title=title)
        msg = f"{dtm.datetime.utcnow()}, @plot_domain_emis, hvplot prepared"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        # print(f"DEBUG::type(cat_hvplot)={type(cat_hvplot)}")
        return cat_hvplot.opts(backend_opts={"plot.toolbar.autohide": True})


    @lru_cache
    def setup_emis(self):
        msg = f"{dtm.datetime.utcnow()}, @setup_emis, start"
        logger.debug(msg)
        self.data = OrderedDict()
        self.glob_timeseries = OrderedDict()
        self.dates = None
        for reg in ['glb600x400','eur300x200','gns100x100']:
            msg = f"...@setup_emis@{reg} reading input"
            logger.debug(msg)
            #-- check whether emissions for region were generated
            _file_lst = sorted(list(emisdir.glob(f"ch4emis.CH4.{reg}.*.nc")))
            if len(_file_lst)>0:
                emis_ptn = f"{str(emisdir)}/ch4emis.CH4.{reg}.*.nc"
                ### shorten list of files for speeding up when devloping
                emis_ptn = f"{str(emisdir)}/ch4emis.CH4.{reg}.202101??.nc"
                ###
                cur_emis = load_emis(emis_ptn)
                cur_dates = {Timestamp(v).strftime('%B %Y'): iv for (iv, v) in enumerate(cur_emis.time.values)}
                if self.dates is None:
                    self.dates = cur_dates
                else:
                    #-- make sure all regions are for the same dates
                    assert np.all(self.dates==cur_dates)
                #
                #-- add cat emissions field
                #
                self.data[reg] = cur_emis
                #
                #-- caching regional-sum time-series (per category),
                #   otherwise we have long delay when visualising...
                #
                msg = f"...@setup_emis@{reg}, prepare overall time-series"
                logger.debug(msg)
                cur_glob = cur_emis.sum(('lat','lon')).to_dataframe()
                cur_glob = cur_glob.reset_index()
                self.glob_timeseries[reg] = cur_glob
                msg = f"...@setup_emis@{reg} done"
                logger.debug(msg)
        # print(f"MVODEBUG:setup_emis: list(self.data.keys()) ***{list(self.data.keys())}***")
        msg = f"@setup_emis, finished"
        logger.debug(msg)

    @lru_cache
    def setup_emis_mp(self):
        from multiprocessing import Process, Manager

        def do_region(reg):
            msg = f"...@setup_emis@{reg} reading input"
            logger.debug(msg)
            #-- check whether emissions for region were generated
            _file_lst = sorted(list(emisdir.glob(f"ch4emis.CH4.{reg}.*.nc")))
            if len(_file_lst)>0:
                emis_ptn = f"{str(emisdir)}/ch4emis.CH4.{reg}.*.nc"
                cur_emis = load_emis(emis_ptn)
                cur_dates = {Timestamp(v).strftime('%B %Y'): iv for (iv, v) in enumerate(cur_emis.time.values)}
                mp_dates[reg] = cur_dates
                # print(f"@do_region for -->{reg}<--, mp_dates -->{mp_dates.keys()}")
                # if ireg==1:
                #     self.dates = cur_dates
                # else:
                #     #-- make sure all regions are for the same dates
                #     assert np.all(self.dates==cur_dates)
                #
                #-- add cat emissions field
                #
                mp_data[reg] = cur_emis
                # self.data[reg] = cur_emis
                #
                #-- caching regional-sum time-series (per category),
                #   otherwise we have long delay when visualising...
                #
                msg = f"...@setup_emis@{reg}, prepare overall time-series"
                logger.debug(msg)
                cur_glob = cur_emis.sum(('lat','lon')).to_dataframe()
                cur_glob = cur_glob.reset_index()
                # self.glob_timeseries[reg] = cur_glob
                mp_domtseries[reg] = cur_glob
                msg = f"...@setup_emis@{reg} done"
                logger.debug(msg)

        msg = f"{dtm.datetime.utcnow()}, @setup_emis, start"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        dom_list = ['glb600x400','eur300x200','gns100x100']
        manager = Manager()
        mp_dates = manager.dict()
        mp_data  = manager.dict()
        mp_domtseries = manager.dict()
        processes = [Process(target=do_region, args=(_reg,)) for _reg in dom_list ]
        for process in processes:
            process.start()
        #-- wait for all process to complete
        for process in processes:
            process.join()
        # print(f"...processes joined.")
        self.data = OrderedDict()
        self.glob_timeseries = OrderedDict()
        for _dom in dom_list:
            self.data[_dom] = mp_data[_dom]
            self.glob_timeseries[_dom] = mp_domtseries[_dom]
        #-- TODO: should still make sure dates are equal for all regions
        self.dates = mp_dates['glb600x400']

class StationExplorer(pn.viewable.Viewer):
    tracer = param.Selector()
    station = param.Selector()
    data = param.ObjectSelector()
    
    def __init__(self, settings):
        super().__init__()
        self.settings = settings
        self.data = nc4.Dataset(str(stations_file))
        self._dates = [Timestamp(*_) for _ in self.data['date_midpoints'][:]]
        
        ntrac = self.data.dimensions['tracers']
        tracers = [getattr(self.data, f'tracer_{itrac+1:03.0f}') for itrac in range(self.data.dimensions['tracers'].size)]
        self.param.tracer.objects = tracers
        self.tracer = tracers[0]
        
        stations = sorted(set([self.data[_].getncattr('name') for _ in self.data.groups]))
        self.param.station.objects = sorted(stations)
        self.station = stations[-1]
        
    def __panel__(self):
        return pn.Column(
            pn.widgets.Select.from_param(self.param.tracer),
            pn.widgets.Select.from_param(self.param.station),
            hv.DynamicMap(self.plot_timeseries),
        )

    def _get_unit(self, station_id):
        stagrp = self.data[station_id]
        ncmix = stagrp['mixing_ratio']
        assert ncmix.dimensions==('tracers','samples')
        mix_unit = ncmix.unit
        assert mix_unit.startswith('mole fraction (')
        assert mix_unit.endswith(')')
        mix_unit = mix_unit[:-1].replace('mole fraction (','')
        return mix_unit
    
    # @pn.depends('tracer', 'station', watch=True)
    @pn.depends('tracer', 'station')
    def plot_timeseries(self):
        if self.tracer is None or self.station is None:
            return
        #-- temporal coverage
        tcover_start = pd.Timestamp(self.data.getncattr('starting time'))
        tcover_end   = pd.Timestamp(self.data.getncattr('ending time'))
        # print(tcover_start)
        # print(tcover_end)
        #-- pick stations with matching name
        station_ids = [_ for _ in self.data.groups if self.data[_].getncattr('name') == self.station]
        #
        #-- currently restrict to upper-most vertical level
        #
        station_level = [self.data[_].getncattr('altitude') for _ in station_ids]
        station_ids = station_ids[station_level.index(max(station_level))]
        #
        #--
        #
        # print(f"DEBUG ***{station_ids}***")
        stagrp = self.data[station_ids]
        abbr_tag = stagrp.abbr.replace('FM/','')
        sta_alt = stagrp.altitude
        stalon = stagrp.longitude
        stalat = stagrp.latitude
        sta_region   = stagrp.region
        ncmix = stagrp['mixing_ratio']
        itrac = self.param.tracer.objects.index(self.tracer)
        staconc = ncmix[itrac,:].data
        data_dict = {'time':self._dates, self.tracer: staconc}
        df = DataFrame.from_dict(data_dict)
        #
        #-- compare against cams
        #
        dfcams = cams_at_obspack_load_conctseries(camsfile,
                                                  stalon,
                                                  stalat,
                                                  sta_alt)
        if not dfcams is None:
            #
            #-- merging daily averages from TM5 simulation and CAMS
            #   (
            dfd = df.copy()
            dfd.index = dfd['time']
            dfd = dfd.drop(['time',], axis=1)
            dfd = dfd.resample('D').mean()
            dfcams.index = dfcams['time']
            dfcams = dfcams.drop(['time',], axis=1)
            dfd_cams = dfcams.resample('D').mean()
            dfplot = pd.merge(dfd, dfd_cams,
                              left_index=True, right_index=True, how='inner')
            dfplot = dfplot.reset_index()
            ### TESTING ONLY
            # dfplot = dfd.reset_index()
        else:
            dfplot = df
        #
        #-- add comparison against obspack
        #
        if obspackdir!=None:
            dfobspack = obspack_load_conctseries(obspackdir,
                                                 tcover_start,
                                                 tcover_end,
                                                 abbr_tag,
                                                 stalon,
                                                 stalat,
                                                 sta_alt)
            #-- there seem to be simulated stations
            #   where no appropriate obspack counterpart is available...
            if not dfobspack is None:
                # msg = f"DEBUG @{abbr_tag}, obspack preparation ({dfobspack.shape}), " \
                #     f"{dfobspack['time'].min()} to {dfobspack['time'].max()}"
                # print(msg)
                #
                #-- for now: convert to daily means
                #
                dfobspack.index = dfobspack['time']
                dfobspack = dfobspack.drop(['time',], axis=1)
                dfd_obspack = dfobspack.resample('D').mean()
                dfd_obspack = dfd_obspack[dfd_obspack.notnull()]
                #-- dfplot: make 'time' become index again
                # print(dfd_obspack.head())
                dfplot.index = dfplot['time']
                dfplot = dfplot.drop(['time',], axis=1)
                #-- insert emtpy 'obspack' column
                dfplot['obspack'] = np.nan
                #-- insert obspack data at dates where available
                for _date,_row in dfd_obspack.iterrows():
                    try:
                        _ch4 = float(_row['obspack_ch4'])
                        dfplot.loc[_date,'obspack'] = _ch4
                    except TypeError:
                        pass
                    # if not type(_ch4)==pd.NaT:
                    #     print(f"xxx{_date}xxx")
                    #     print(f"zzz{type(_ch4)}zzz {type(_ch4)==pd.NaT}")
                    #     print(f"yyy{_ch4}yyy")
                    #     dfplot.loc[_date,'obspack'] = _ch4
                # print(dfplot.head())
                #-- 'time' becomes column again
                dfplot = dfplot.reset_index()
                # print(dfplot.head(n=3))
        #--
        title = f"{species}@{self.station} ({sta_region}, {sta_alt}[m])"
        mixunit = self._get_unit(station_ids)
        ylabel = f"concentration [{mixunit}]"

        #
        #-- make more descriptive name for plotting
        #
        if self.tracer==species:
            rename_map = { self.tracer : 'FIT-IC' }
        else:
            rename_map = { self.tracer : f"FIT-IC ({self.tracer})" }
        dfplot = dfplot.rename(columns=rename_map)

        #--
        plot_columns = ['time', rename_map[self.tracer],]
        for _c in dfplot.columns:
            if _c.startswith('cams') or _c.startswith('obspack'):
                plot_columns.append(_c)
        # for _c in ['cams','obspack',]:
        #     if _c in dfplot.columns:
        #         plot_columns.append(_c)
        dfplot = dfplot[plot_columns]
        # print(dfplot.head(n=3))
        # print(f"TM5 entries: {dfplot[self.tracer].notnull().sum()}")
        # print(f"obspack:     {dfplot['obspack'].notnull().sum()}")
        hvret = dfplot.hvplot(x='time', grid=True,
                              ylabel=ylabel, xlabel='time', title=title)
        # print(f"@plot_timeseries, returning type -->{type(hvret)}<--")
        return hvret
