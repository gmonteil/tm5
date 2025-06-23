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
#
#-- to avoid PROJ path warning message on JupyterHub,
#   fix_env would need to be called prior import of ccrs.
#   For now, this import of crs is made only in
#   routines where needed!
# import cartopy.crs as ccrs
hv.extension('bokeh')


#-- local packages
#   NOTE: preliminary imports currently,
#         will/need to be accessible eventually via the FIT-IC environment.
from ui_util import fix_env
from ui_util import cams_at_obspack_load_conctseries
from ui_util import obspack_load_conctseries
from ui_util import outdir_default, outdir_overwrite, outdir_edgarflat
from uu_util import camsfile, obspackdir #-- Note: camsfile not required in Geneva25 demo

#
#
#-- FIT-IC simulations currently limited to CH4 as single species
#
species = 'CH4'
supported_domain_list = ['glb600x400','eur300x200','gns100x100']

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
    
    def __init__(self, settings, pattern : str = None, mode : str = 'precomputed_default', load_parallel : bool = True):
        super().__init__()
        fix_env()
        self.settings = settings
        #-- emission file naming convention is
        #   ch4emis.{species}.{reg}.yyyymmdd.nc
        self.yyyymmdd_ptn = pattern if pattern!=None else '????????'
        #
        #--
        #
        if mode=='precomputed_default':
            self.outdir_precomputed = outdir_default
        elif mode=='precomputed_regional':
            self.outdir_precomputed = outdir_overwrite
        else:
            msg = f"selected mode -->{mode}<-- not supported."
            raise RuntimeError(msg)
        #
        self.emisdir       = self.outdir_precomputed / 'emissions'
        #-- basic consistency tests
        if not self.emisdir.is_dir():
            msg = f"directory for emissions ==>{self.emisdir}<== is not accessible"
            raise RuntimeError(msg)
        else:
            msg = f"emissions directory ***{self.emisdir}***"
            logger.debug(msg)

        #
        #-- read and prepare emissions to be ready for plotting
        #
        if load_parallel:
            self.setup_emis_mp() #-- multiprocessor version to speed-up
        else:
            self.setup_emis()
        #
        #--
        #
        msg = f"{dtm.datetime.utcnow()}, @__init__, setting self.param.region.objects"
        # print(f"DEBUG::{msg}")
        self.param.region.objects = list(self.data.keys())
        self.region = self.param.region.objects[-1]
        self.region = self.param.region.objects[1]
        msg = f"{dtm.datetime.utcnow()}, @__init__, setting self.param.category.objects"
        # print(f"DEBUG::{msg}")
        datavar_list = sorted(list(self.data[self.region].data_vars))
        self.param.category.objects = datavar_list
        #
        #-- start with 'fossil' as initial category (if present),
        #   otherwise start with first in list
        #
        istart = 0
        for i,varname in enumerate(datavar_list):
            if varname=='fossil':
                istart = i
                break
        self.category = self.param.category.objects[istart]
        

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

        # dyn_params = {'cache_size':1}
        # widget_map1 = hv.DynamicMap(self.map_emis, **dyn_params)
        # widget_map1 = self.map_emis
        # widget_map1 = hv.DynamicMap(self.map_emis, **dyn_params)
        widget_map1 = hv.DynamicMap(self.map_emis)
        widget_map2 = hv.DynamicMap(self.plot_domain_emis)
        # widget_map2 = self.plot_domain_emis
        widget = pn.Column(widget_regsel,widget_fldsel,widget_datesel,
                           pn.Row(widget_map1,widget_map2))
        return widget

    #-- probably better drop 'watch=True' (?)
    #   (see https://panel.holoviz.org/how_to/param/dependencies.html)
    @property
    # @pn.depends('time_index', watch=True)
    @pn.depends('time_index')
    def current_date(self):
        return self.data[self.region].time.values[self.time_index]

    # @property
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
        lonmin,lonmax,latmin,latmax = self.current_extent()
        # print(f"DEBUG::@map_emis title -->{title}<--")
        # features features (default=None): A list of features or a dictionary of features and the scale at which to render it. Available features include ‘borders’, ‘coastline’, ‘lakes’, ‘land’, ‘ocean’, ‘rivers’ and ‘states’. Available scales include ‘10m’/’50m’/’110m’.
        if self.region=='glb600x400':
            cproj = ccrs.PlateCarree()
            cfeatures = {'borders':'110m', 'coastline':'110m'}
            coastline = '110m'
        elif self.region=='eur300x200':
            cproj = ccrs.GOOGLE_MERCATOR
            cproj = ccrs.PlateCarree()
            cfeatures = {'borders':'50m', 'coastline':'50m'}
            coastline = '50m'
        elif self.region=='gns100x100':
            cproj = ccrs.GOOGLE_MERCATOR
            cproj = ccrs.PlateCarree()
            cfeatures = {'borders':'10m', 'coastline':'10m'}
            coastline = '10m'
        cur_emis = self.data[self.region][cur_cat].isel(time=self.time_index)
        unitlabel = f"[kg{species}/cell/s]"
        #
        #-- compute total
        #
        emis_tot = cur_emis.sum(('lat','lon')).values #- kgTRACER/s
        emis_tot = emis_tot / 1e6 * (365*86400)       #- ktTRACER/year
        emis_tot_unit = f"[kt{species}/year]"
        if emis_tot>=1.e3:
            emis_tot = emis_tot / 1e3
            emis_tot_unit = f"[Mt{species}/year]"
        msg = f"{dtm.datetime.utcnow()}, @map_emis, current data ready"
        logger.debug(msg)
        title = f"{cur_cat}@{self.region} " \
            f"({cur_date}, total={emis_tot:.3f}{emis_tot_unit})"
        title = f"{cur_cat}@{self.region} {unitlabel} ({cur_date})" + '\n' \
            f"domain total: {emis_tot:.3f}{emis_tot_unit}"
        # with open(f"map_emis_{dtm.datetime.utcnow().isoformat()}.log", 'w') as fp:
        #     msg = f"***{title}*** coastline={coastline} cfeatures={cfeatures}"
        #     fp.write(f"{msg}" +'\n')
        # emis_hvplot = cur_emis.hvplot.quadmesh(geo=True,
        #                                        coastline=coastline,
        #                                        features=cfeatures,
        #                                        xlim=(lonmin,lonmax),
        #                                        ylim=(latmin,latmax),
        #                                        clabel=unitlabel, title=title)
        emis_hvplot = cur_emis.hvplot.quadmesh(
            xlim=(lonmin,lonmax), ylim=(latmin,latmax),
            coastline=coastline, title=title,
            crs=ccrs.PlateCarree(),
            projection=cproj, features=cfeatures, project=True)#, rasterize=True)
        msg = f"{dtm.datetime.utcnow()}, @map_emis, hvplot ready"
        logger.debug(msg)
        #
        #-- autohide toolbar
        #
        emis_hvplot = emis_hvplot.opts(
            # colorbar_options={'label':unitlabel},
            backend_opts={"plot.toolbar.autohide": True})
        # #-- 'framewise=True' did not help to update maps
        # #   when changing region...
        # emis_hvplot = emis_hvplot.opts(framewise=True, backend_opts={"plot.toolbar.autohide": True})
        return emis_hvplot

    # @pn.depends('time_index', 'region', 'category', watch=True)
    @pn.depends('time_index', 'region', 'category')
    def plot_domain_emis(self):
        msg = f"{dtm.datetime.utcnow()}, start@plot_domain_emis with time_index={self.time_index}, " \
            f"region -->{self.region}<--  category -->{self.category}<--"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        if self.region==None:
            return
        cat_df = self.glob_timeseries[self.region]
        ylabel = f"[kg{species}/s]"
        if self.region=='glb600x400':
            title = f"global sectoral {species} emission totals (@{self.region})"
        else:
            title = f"domain sectoral {species} emission totals (@{self.region})"
        # print(cat_df.head())
        # cat_hvplot = cat_df.hvplot(grid=True, x='time',
        #                            xlabel='time', ylabel=ylabel, title=title,
        #                            fontsize={'legend':6})
        cat_hvplot = cat_df.hvplot(grid=True, x='time',
                                   xlabel='time', ylabel=ylabel, title=title)
        # cat_hvplot.opts(backend_opts={"plot.toolbar.autohide": True}, legend_position='bottom_right', legend_offset=(0,0), fontsize={'legend':6,'legend_title':6})
        cat_hvplot.opts(backend_opts={"plot.toolbar.autohide": True}, fontsize={'legend':6,'legend_title':6})
        # cat_hvplot = cat_df.plot(grid=True, x='time',
        #                          xlabel='time', ylabel=ylabel, title=title)
        msg = f"{dtm.datetime.utcnow()}, @plot_domain_emis, hvplot prepared"
        logger.debug(msg)
        # print(f"DEBUG::{msg}")
        # print(f"DEBUG::type(cat_hvplot)={type(cat_hvplot)}")
        return cat_hvplot


    @lru_cache
    def setup_emis(self):
        msg = f"{dtm.datetime.utcnow()}, @setup_emis, start"
        logger.debug(msg)
        self.data = OrderedDict()
        self.glob_timeseries = OrderedDict()
        self.dates = None
        for reg in supported_domain_list:
            msg = f"...@setup_emis@{reg} reading input"
            logger.debug(msg)
            #-- check whether emissions for region were generated
            fptn = f"ch4emis.{species}.{reg}.{self.yyyymmdd_ptn}.nc"
            _file_lst = sorted(list(self.emisdir.glob(fptn)))
            if len(_file_lst)==0:
                msg = f"no emissions files detected with pattern " \
                    f"==>{fptn}<=="
                raise RuntimeError(msg)
            else:
                emis_ptn = f"{str(self.emisdir)}/{fptn}"
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
            _file_lst = sorted(list(self.emisdir.glob(f"ch4emis.{species}.{reg}.{self.yyyymmdd_ptn}.nc")))
            if len(_file_lst)>0:
                emis_ptn = f"{str(self.emisdir)}/ch4emis.{species}.{reg}.{self.yyyymmdd_ptn}.nc"
                cur_emis = load_emis(emis_ptn)
                cur_dates = {Timestamp(v).strftime('%B %Y'): iv for (iv, v) in enumerate(cur_emis.time.values)}
                mp_dates[reg] = cur_dates
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
        manager = Manager()
        mp_dates = manager.dict()
        mp_data  = manager.dict()
        mp_domtseries = manager.dict()
        processes = [Process(target=do_region, args=(_reg,)) for _reg in supported_domain_list ]
        for process in processes:
            process.start()
        #-- wait for all process to complete
        for process in processes:
            process.join()
        # print(f"...processes joined.")
        self.data = OrderedDict()
        self.glob_timeseries = OrderedDict()
        for _dom in supported_domain_list:
            self.data[_dom] = mp_data[_dom]
            self.glob_timeseries[_dom] = mp_domtseries[_dom]
        #-- TODO: should still make sure dates are equal for all regions
        self.dates = mp_dates['glb600x400']

class StationExplorer(pn.viewable.Viewer):
    tracer  = param.Selector()
    station = param.Selector()
    data    = param.ObjectSelector()
    
    def __init__(self, settings, mode : str = 'precomputed_default'):
        super().__init__()
        self.settings = settings
        #
        #--
        #
        self.precompoutdir_default = outdir_default
        self.precompoutdir_flat    = outdir_edgarflat
        self.precompoutdir_ovr     = outdir_overwrite
        #
        #-- up to 3 files with simulated concentrations at stations
        #
        self.stations_file_default = self.precompoutdir_default / 'stations' / 'stations.nc4'
        self.stations_file_flat    = self.precompoutdir_flat / 'stations' / 'stations.nc4'
        #-- basic consistency tests
        if not self.stations_file_default.exists():
            msg =f"file with simulated concentrations at stations " \
                f"==>{self.stations_file_default}<== is not accessible."
            raise RuntimeError(msg)
        elif not self.stations_file_flat.exists():
            msg =f"file with simulated concentrations at stations " \
                f"==>{self.stations_file_flat}<== is not accessible."
            raise RuntimeError(msg)
        #
        #-- open station NetCDF files
        #
        self.data_default = nc4.Dataset(str(self.stations_file_default))
        self.data_flat    = nc4.Dataset(str(self.stations_file_flat))
        #
        #--
        #
        assert np.all(self.data_default['date_midpoints'][:]==self.data_flat['date_midpoints'][:])
        #
        #-- store concentration dates
        #
        self._dates = [Timestamp(*_) for _ in self.data_default['date_midpoints'][:]]
        #
        #-- store tracer names (expected to be equal among simulations)
        #
        ntrac = self.data_default.dimensions['tracers']
        tracers = [getattr(self.data_default, f'tracer_{itrac+1:03.0f}') for itrac in range(self.data_default.dimensions['tracers'].size)]
        self.param.tracer.objects = tracers
        self.tracer = tracers[0]
        #

        #
        #-- for the panel
        #
        stations = sorted(set([self.data_default[_].getncattr('name') for _ in self.data_default.groups]))
        self.param.station.objects = sorted(stations)
        self.station = stations[-1]
        
    def __panel__(self):
        return pn.Column(
            pn.widgets.Select.from_param(self.param.tracer),
            pn.widgets.Select.from_param(self.param.station),
            hv.DynamicMap(self.plot_timeseries),
        )

    def _get_unit(self, station_id):
        #
        #-- mixing unit is equal among default/overwrite/edgarflat
        #
        stagrp = self.data_default[station_id]
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
        tcover_start = pd.Timestamp(self.data_default.getncattr('starting time'))
        tcover_end   = pd.Timestamp(self.data_default.getncattr('ending time'))
        # print(tcover_start)
        # print(tcover_end)
        #
        #-- pick stations with matching name
        #
        station_ids = [_ for _ in self.data_default.groups if self.data_default[_].getncattr('name') == self.station]
        #
        #-- currently restrict to upper-most vertical level
        #
        station_level = [self.data_default[_].getncattr('altitude') for _ in station_ids]
        station_ids = station_ids[station_level.index(max(station_level))]
        #
        #--
        #
        # print(f"DEBUG ***{station_ids}***")
        stagrp = self.data_default[station_ids]
        abbr_tag = stagrp.abbr.replace('FM/','')
        sta_alt = stagrp.altitude
        stalon = stagrp.longitude
        stalat = stagrp.latitude
        sta_region   = stagrp.region
        #
        #-- load (default) concentration at station (and level),
        #   insert into data frame
        #
        ncmix = stagrp['mixing_ratio']
        itrac = self.param.tracer.objects.index(self.tracer)
        staconc = ncmix[itrac,:].data
        data_dict = {'time':self._dates, self.tracer: staconc}
        df = DataFrame.from_dict(data_dict)
        #
        #-- add simulation with flat profile
        #
        stagrp_flat = self.data_flat[station_ids]
        ncmix_flat = stagrp_flat['mixing_ratio']
        itrac = self.param.tracer.objects.index(self.tracer)
        staconc_flat = ncmix_flat[itrac,:].data
        df['FIT-IC (flat anthropogenic)'] = staconc_flat[:]
        #-- MVO-20250622:don't show CAMS in Geneva demo
        # #
        # #-- compare against cams
        # #
        # dfcams = cams_at_obspack_load_conctseries(camsfile,
        #                                           stalon,
        #                                           stalat,
        #                                           sta_alt)
        dfcams = None #-- i.e. not showing CAMS concentrations
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
        hvplot = dfplot.hvplot(x='time', grid=True,
                              ylabel=ylabel, xlabel='time', title=title)
        # print(f"@plot_timeseries, returning type -->{type(hvplot)}<--")
        #
        #
        #-- autohide toolbar
        #
        hvplot = hvplot.opts(backend_opts={"plot.toolbar.autohide": True})

        return hvplot
