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
from ui_util import get_precomputed_expdirs
from ui_util import obspackdir

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
    data       = param.ObjectSelector()
    category   = param.Selector()
    region     = param.Selector()
    
    def __init__(self, settings, pattern : str = None,
                 mode : str = 'precomputed_default',
                 load_parallel : bool = True):
        super().__init__()
        fix_env()
        self.settings = settings
        #-- emission file naming convention is
        #   ch4emis.{species}.{reg}.yyyymmdd.nc
        self.yyyymmdd_ptn = pattern if pattern!=None else '????????'
        #
        #-- directories with pre-computed output
        #
        precomp_dict = get_precomputed_expdirs()
        if mode=='precomputed_default':
            self.outdir_precomputed = precomp_dict['default']
        elif mode=='precomputed_regional':
            self.outdir_precomputed = precomp_dict['regional']
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
    #-- 2025-09-09: vertical level not selectable anymore
    # vlevel  = param.Selector(label='vertical level at station')
    exp1 = param.Selector(label='first TM5 simulation')
    exp2 = param.Selector(label='second TM5 simulation')
    #-- swichted back to daily-mean concentration
    #hour    = param.Selector()
    data    = param.ObjectSelector()

    def __init__(self, settings):
        super().__init__()
        self.settings = settings
        #
        #--
        #
        self.plot_width = 1000
        self.plot_height = 500
        #
        #-- station plot comparison against obspack
        #   currently limited to single level
        #   (no longer switch to select highest/lowest level)
        #
        self.vlevel = 'highest'
        self.default_site = 'Cabauw' #None
        #
        #-- directories with pre-computed output
        #
        precomp_dict = get_precomputed_expdirs()
        #
        #-- build table of station-simulation file handles
        #
        self.precomp_table = OrderedDict()
        simu_dates = None
        tcover_start = None
        tcover_end   = None
        tracers    = None
        stations   = None
        for tag,expdir in precomp_dict.items():
            stations_file = precomp_dict[tag] / 'stations' / 'stations.nc4'
            if not stations_file.exists():
                msg = f"@{tag}, did not find 'stations.nc4' in simulation directory"
                logger.warning(msg)
                continue
            else:
                #-- NetCDF file handle
                _cur_handle = nc4.Dataset(str(stations_file))
                self.precomp_table[tag] = _cur_handle
            #
            #-- save simulation dates, assure consinstency among simulations
            #
            _cur_tstart = _cur_handle.getncattr('starting time')
            _cur_tend   = _cur_handle.getncattr('ending time')
            _cur_dates = _cur_handle['date_midpoints'][:]
            if len(self.precomp_table)==1:
                simu_dates = np.copy(_cur_dates)
                tcover_start = _cur_tstart
                tcover_end   = _cur_tend
            else:
                assert np.all(simu_dates[:]==_cur_dates[:])
            #
            #-- make sure all simulations are for the same (list of) tracers
            #
            _cur_tracers = [getattr(_cur_handle, f'tracer_{itrac+1:03.0f}') for itrac in range(_cur_handle.dimensions['tracers'].size)]
            if len(self.precomp_table)==1:
                tracers = _cur_tracers
            else:
                assert _cur_tracers==tracers
            #
            #-- build list of stations
            #
            _cur_stations = sorted(set([_cur_handle[_].getncattr('name') for _ in _cur_handle.groups]))
            if len(self.precomp_table)==1:
                stations = _cur_stations
            else:
                assert _cur_stations==stations
        #
        #-- store temporal settings
        #
        self._dates = [Timestamp(*_) for _ in simu_dates]
        self.tcover_start = pd.Timestamp(tcover_start)
        self.tcover_end   = pd.Timestamp(tcover_end)
        #
        #-- set widget selectors
        #
        self.param.tracer.objects = tracers
        self.tracer = tracers[0]
        self.param.station.objects = sorted(stations)
        #-- vertical level not selectable anymore (stick with highest level)
        # self.param.vlevel.objects  = ['highest','lowest']
        # self.vlevel                = 'highest'
        #
        #-- defaults to last site in list
        #
        if self.default_site==None:
            try:
                idx = stations.index("Mauna Loa, Hawaii")
                idx = stations.index("Cabauw")
                # idx = stations.index("Palmer Station, Antarctica")
            except ValueError:
                idx = -1
            self.station = stations[idx]
        else:
            try:
                idx = stations.index(self.default_site)
            except ValueError:
                idx = -1
            self.station = stations[idx]
        #
        # self.widgets = dict(
        #     info1 = pn.pane.Markdown(width=250),
        #     info2 = pn.pane.Markdown(width=250)
        #     )
        self.widgets = dict(
            info = pn.pane.Markdown(width=self.plot_width,
                                    styles={'font-size': '14px'}),
            )
        #
        #-- re-order experiments for the selection:
        #   - ad-hoc solution here, but for the demonstration
        #     it is favoutable to start with default/edgarflat/regional
        #
        exp_list = list( self.precomp_table.keys())
        self.exp_list = []
        #
        #-- standard simulations first
        #
        exp_list_first = ['default','edgarflat','regional',
                          'half-oh','no-germany','no-northamerica',]
        for exp in exp_list_first:
            if exp in exp_list:
                self.exp_list.append(exp)
        #
        #-- sensitivity experiments based on regional changes
        #
        exp_list_regional = ['no-agri','no-fossil','no-waste',
                             'anthro-no-france','anthro-no-netherlands',]
        exp_list_regional = [ f'regional_{_}' for _ in exp_list_regional]
        for exp in exp_list_regional:
            if exp in self.exp_list:
                continue
            elif exp in exp_list:
                self.exp_list.append(exp)
            else:
                pass
        #
        #-- any remaining sensitivity experiments
        #
        for exp in exp_list:
            if exp in self.exp_list:
                continue
            else:
                self.exp_list.append(exp)

        #
        self.param.exp1.objects = self.exp_list
        self.param.exp2.objects = self.exp_list
        iexp1 = self.exp_list.index('default')
        iexp2 = self.exp_list.index('regional')
        self.exp1 = self.exp_list[iexp1]
        self.exp2 = self.exp_list[iexp2]
        # self.update_desc_exp()
        self.set_experiments_desc_table()
        
    def __panel__(self):
        return pn.Column(
            self.widgets['info'], #-- list all available experiments
            pn.widgets.Select.from_param(self.param.tracer),
            pn.widgets.Select.from_param(self.param.station),
            #-- 2025-09-09: vertical level not selectable anymore
            # pn.Row(pn.widgets.Select.from_param(self.param.station),
            #        pn.widgets.Select.from_param(self.param.vlevel)),
            pn.Row(pn.widgets.Select.from_param(self.param.exp1),
                   pn.widgets.Select.from_param(self.param.exp2)),
            # pn.Row(pn.Column(pn.widgets.Select.from_param(self.param.exp1),
            #                  self.widgets['info1']),
            #        pn.Column(pn.widgets.Select.from_param(self.param.exp2),
            #                  self.widgets['info2'])
            #        ),
            #-- was here when 'info' widget only showed the two
            #   active experiments
            # self.widgets['info'],
            # pn.widgets.Select.from_param(self.param.hour),
            hv.DynamicMap(self.plot_timeseries),
        )

    def _get_desc(self, exp):
        match exp:
            case 'default':
                desc = f"TM5 run using global prior default emissions"
            case 'edgarflat':
                desc = f"Similar to the default case, but using a flat " \
                    "temporal profile for EDGAR anthropogenic emissions."
            case 'regional':
                desc = f"TM5 run using wetland, mineral-soils, " \
                    f"and anthropogenic emissions provided by " \
                    f"AVENGERS WP2 over the European domain, " \
                    f"and with the global default emissions elsewhere."
            case 'regional_no-agri' | 'regional_anthro-no-agri':
                desc = f"Prior emissions similar to the regional case " \
                    f"but without emissions from agriculture sector " \
                    f"in the European domain."
            case 'regional_no-waste' | 'regional_anthro-no-waste':
                desc = f"Prior emissions similar to the regional case " \
                    f"but without emissions from waste sector " \
                    f"in the European domain."
            case 'regional_no-fossil' | 'regional_anthro-no-fossil':
                desc = f"Prior emissions similar to the regional case " \
                    f"but without emissions from fossil sector " \
                    f"in the European domain."
            case 'regional_no-france' | 'regional_anthro-no-france':
                desc = f"Prior emissions similar to the regional case " \
                    f"but without anthropogenic emissions over France."
            case 'regional_no-netherlands' | 'regional_anthro-no-netherlands':
                desc = f"Prior emissions similar to the regional case " \
                    f"but without anthropogenic emissions over the Netherlands. "
            case 'no-northamerica':
                desc = f"Similar to default case " \
                    f"but without emissions over Northern America."
            case 'no-germany':
                desc = f"Similar to default case " \
                    f"but without emissions over domain around Germany (6E-15E,47N-55N)."
            case 'half-oh':
                desc = f"Similar to default case " \
                    f"but TM5 only 50% of the CAMS OH field is used chemistry."
            case _:
                desc = f"no description available yet."

        return desc


    def set_experiments_desc(self):
        def _get_desc(exp):
            desc_html = "<dt>"
            desc_html += '\n' + f'<b><span style="text-decoration:underline">{exp}</span></b> '
            desc_html += '\n' +"</dt>"
            desc = self._get_desc(exp)
            desc_html += f"<dd>{desc}</dd>"
            return desc_html

    def set_experiments_desc_list(self):
        def _get_desc(exp):
            desc_html = "<li>"
            desc_html += f'<b>{exp}: </b>'
            desc_html += '<span style="margin-left:30px">'
            desc = self._get_desc(exp)
            desc_html += f"{desc}</span></li>"
            return desc_html
        #
        #--
        #
        desc = f'<b><span style="text-decoration:underline">List of precomputed TM5 simulations:</span></b>'
        desc +="<br>"
        desc += f"<ul>"
        for iexp,exp in enumerate(self.exp_list):
            desc += _get_desc(exp)
            if iexp==len(self.exp_list)-1:
                pass
            else:
                desc += '<br>'
        desc += f"</ul>"
        self.widgets['info'].object = desc
            
    def set_experiments_desc_table(self):
        def _get_desc(exp):
            desc = self._get_desc(exp)
            desc_html = "<tr>"
            # desc_html += f'<td><b>{exp}</b></td>'
            desc_html += f'<td>{exp}</td>'
            desc_html += f'<td>{desc}</td>'
            desc_html += f'</tr>'

            return desc_html
        #
        #--
        #
        desc = f'<table>'
        desc += f'<caption><b>List of precomputed TM5 simulations<b></caption>'
        desc += '<thead><tr><th>Identifier</th><th>Description</th></tr></thead>'
        for iexp,exp in enumerate(self.exp_list):
            desc += _get_desc(exp)
            if iexp==len(self.exp_list)-1:
                pass
            else:
                desc += '<br>'
        desc += f"</table>"
        self.widgets['info'].object = desc

    # @pn.depends('exp1','exp2',watch=True)
    # def update_desc_exp(self):
    #     def _get_desc(exp):
    #         desc_html = "<dt>"
    #         # desc_html += '\n' + f'<b><span style="text-decoration:underline">{exp}</span></b> '
    #         desc_html += '\n' + f'<b>{exp}</b>'
    #         desc_html += '\n' +"</dt>"
    #         match exp:
    #             case 'default':
    #                 desc = f"TM5 run using global prior default emissions"
    #             case 'edgarflat':
    #                 desc = f"similar to the default case, but using a flat " \
    #                     "temporal profile for EDGAR anthropogenic emissions."
    #             case 'regional':
    #                 desc = f"TM5 run using wetland, mineral-soils, " \
    #                     f"and anthropogenic emissions provided by " \
    #                     f"AVENGERS WP2 over the European domain, " \
    #                     f"and with the global default emissions elsewhere."
    #             case 'regional_no-agri':
    #                 desc = f"Prior emissions similar to the regional case " \
    #                     f"but without emissions from agriculture sector " \
    #                     f"in the European domain."
    #             case 'regional_no-waste':
    #                 desc = f"Prior emissions similar to the regional case " \
    #                     f"but without emissions from waste sector " \
    #                     f"in the European domain."
    #             case 'regional_no-fossil':
    #                 desc = f"Prior emissions similar to the regional case " \
    #                     f"but without emissions from fossil sector " \
    #                     f"in the European domain."
    #             case 'regional_no-france':
    #                 desc = f"Prior emissions similar to the regional case " \
    #                     f"but without anthropogenic emissions over France."
    #             case 'regional_no-netherlands':
    #                 desc = f"Prior emissions similar to the regional case " \
    #                     f"but without anthropogenic emissions over the Netherlands. "
    #             case _:
    #                 desc = f"no description available yet."
    #         desc_html += f"<dd>{desc}</dd>"
    #         return desc_html
    #     #--
    #     desc = f"<dl>"
    #     desc += _get_desc(self.exp1)
    #     desc += '<br>'
    #     desc += _get_desc(self.exp2)
    #     desc += f"</dl>"
    #     self.widgets['info'].object = desc

    def _get_unit(self, station_id):
        def first(s):
            '''Return the first element from an ordered collection
            or an arbitrary element from an unordered collection.
            Raise StopIteration if the collection is empty.
            '''
            return next(iter(s))

        #
        #-- mixing unit is equal among TM5 simulations
        #
        #-> Note: expected to be equal for all simulations,
        #         so can be extracted from any of the table entries
        exp1 = first(self.precomp_table.keys())
        stagrp = self.precomp_table[exp1][station_id]
        ncmix = stagrp['mixing_ratio']
        assert ncmix.dimensions==('tracers','samples')
        mix_unit = ncmix.unit
        assert mix_unit.startswith('mole fraction (')
        assert mix_unit.endswith(')')
        mix_unit = mix_unit[:-1].replace('mole fraction (','')
        return mix_unit
    
    # @pn.depends('tracer', 'station', watch=True)
    #-- switched back to daily mean concentration
    # @pn.depends('tracer', 'station', 'hour')
    #-- added selection of experiments
    # @pn.depends('tracer', 'station','vlevel','exp1','exp2')
    #-- 2025-09-09: disabled vertical level
    @pn.depends('tracer', 'station','exp1','exp2', watch=True)
    def plot_timeseries(self):
        if self.tracer is None or self.station is None:# or self.hour is None:
            return
        elif self.exp1 is None or self.exp2 is None:
            return
        # show_hour = self.hour
        simu_colors = ['red','blue',]
        obs_color   = ['orange',]
        obs_color   = 'k'
        obs_marker  = 'o'#'+'
        obs_ms      = 5 #-- marker size
        legend_loc  = 'top_right'
        #
        #-- temporal coverage
        #
        tcover_start = self.tcover_start
        tcover_end   = self.tcover_end
        #
        #--
        #
        simu_tag1 = self.exp1
        simu_tag2 = self.exp2
        simu1 = self.precomp_table[simu_tag1]
        simu2 = self.precomp_table[simu_tag2]
        #
        #-- pick stations with matching name
        #
        station_ids = [_ for _ in simu1.groups if simu1[_].getncattr('name') == self.station]
        #
        #-- vertical level of simulation
        #
        station_level = [simu1[_].getncattr('altitude') for _ in station_ids]
        #
        #-- currently restricted to lowest/highest level
        #
        if self.vlevel=='lowest':
            station_id = station_ids[station_level.index(min(station_level))]
        elif self.vlevel=='highest':
            station_id = station_ids[station_level.index(max(station_level))]
        else:
            msg = f"unexpected vertical level -->{self.vlevel}<--"
            raise RuntimeError(msg)
        #
        #--
        #
        # print(f"DEBUG station_id ***{station_id}***")
        stagrp = simu1[station_id]
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
        coltag = f"{self.tracer}_simu1"
        data_dict = {'time':self._dates, coltag: staconc}
        df = DataFrame.from_dict(data_dict)
        #
        #-- add second simulation, currently the one with regional emissions
        #
        stagrp_cmp = simu2[station_id]
        ncmix_cmp = stagrp_cmp['mixing_ratio']
        itrac = self.param.tracer.objects.index(self.tracer)
        staconc_cmp = ncmix_cmp[itrac,:].data
        coltag = f"{self.tracer}_simu2"
        df[coltag] = staconc_cmp[:]
        # #
        # #-- restrict to selected hour
        # #
        # dfplot = df[df['time'].dt.hour==show_hour]
        #
        #-- back again to daily concentrations
        #
        dfd = df.copy()
        dfd.index = dfd['time']
        dfd = dfd.drop(['time',], axis=1)
        dfd = dfd.resample('D').mean()
        dfd = dfd.reset_index()
        dfplot = dfd

        #
        #-- make more descriptive name for plotting
        #
        if self.tracer==species:
            if simu_tag1!=simu_tag2:
                simu_columns = [f'FIT-IC-{simu_tag1}', f'FIT-IC-{simu_tag2}',]
            else:
                simu_columns = [f'FIT-IC-{simu_tag1}', f'FIT-IC-also-{simu_tag2}',]
        else:
            if simu_tag1!=simu_tag2:
                simu_columns = [f'FIT-IC-{simu_tag1} (f{self.tracer}',
                                f'FIT-IC-{simu_tag2} (f{self.tracer}',]
            else:
                simu_columns = [f'FIT-IC-{simu_tag1} (f{self.tracer}',
                                f'FIT-IC-also-{simu_tag2} (f{self.tracer}',]
        #--
        rename_map = { f"{self.tracer}_simu1" : simu_columns[0],
                       f"{self.tracer}_simu2" : simu_columns[1]  }
        dfplot = dfplot.rename(columns=rename_map)
        # print(f"DEBUG: dfplot columns (after renaming) ***{dfplot.columns}***")
        #
        #--     title and unit
        #
        title = f"{species}@{self.station} (daily-mean, {sta_region}, {sta_alt}[m])"
        mixunit = self._get_unit(station_id)
        xlabel = 'time'
        ylabel = f"concentration [{mixunit}]"
        #
        #--
        #
        #-- add comparison against obspack (if available)
        #
        if obspackdir!=None:
            obspack_info = obspack_load_conctseries(obspackdir,
                                                    tcover_start,
                                                    tcover_end,
                                                    abbr_tag,
                                                    stalon,
                                                    stalat,
                                                    sta_alt)
            if obspack_info!=None:
                dfobspack = obspack_info['df']
            else:
                dfobspack = None
            #
            #-- there seem to be simulated stations
            #   where no appropriate obspack counterpart is available...
            #
            if not dfobspack is None:
                # msg = f"DEBUG @{abbr_tag}, obspack preparation ({dfobspack.shape}), " \
                #     f"{dfobspack['time'].min()} to {dfobspack['time'].max()}"
                # print(msg)
                # #
                # #-- when selecting hour-of-day
                # #
                # dfobspack = dfobspack[dfobspack['time'].dt.hour==show_hour]
                # dfobspack.index = dfobspack['time']
                # dfobspack = dfobspack.drop(['time',], axis=1)
                #
                #--
                #
                dfd_obspack = dfobspack.copy()
                dfd_obspack.index = dfd_obspack['time']
                dfd_obspack = dfd_obspack.drop(['time',], axis=1)
                dfd_obspack = dfd_obspack.resample('D').mean()
                #
                #-- insert obspack concentrations as column 'obspack'
                #   (where available and at correct time-points)
                #
                dfplot.index = dfplot['time']
                dfplot = dfplot.drop(['time',], axis=1)
                # dfplot.loc[dfobspack.index,'obspack'] = dfobspack.loc[:,'obspack_ch4']
                dfplot.loc[dfd_obspack.index,'obspack'] = dfd_obspack.loc[:,'obspack_ch4']
                dfplot = dfplot.reset_index()
                #
                #-- add RMSE to title
                #
                rmse_title = "RMSE: "
                for _sim in simu_columns:
                    rmse = ((dfplot[_sim]-dfplot['obspack'])**2).mean() ** 0.5
                    exp_tag = _sim.replace('FIT-IC-','')
                    rmse_title += f"{exp_tag}={rmse:.3f} "
                if self.plot_width>=1000:
                    title += f",   {rmse_title}"
                else:
                    title += '\n' + rmse_title

        #
        #--
        #
        # print(f"prior plotting, columns -->{simu_columns}<--, colors -->{simu_colors}<--")
        # simplot = dfplot.hvplot(x='time', y=simu_columns,
        #                         color=simu_colors,
        #                         grid=True,
        #                         ylabel=ylabel, xlabel=xlabel, title=title)
        #
        #-- first simulation
        #
        p1 = dfplot.hvplot(x='time', y=simu_columns[0],
                           color=simu_colors[0],
                           label=simu_columns[0], legend=legend_loc,
                           width=self.plot_width, height=self.plot_height,
                           ylabel=ylabel, xlabel=xlabel, title=title)
        #
        #-- second simulation
        #
        p1 *= dfplot.hvplot(x='time', y=simu_columns[1],
                            color=simu_colors[1],
                            label=simu_columns[1], legend=legend_loc)
        #
        #-- observations
        #
        if 'obspack' in dfplot.columns:
            obsplot = dfplot.hvplot.points(x='time', y='obspack',
                                           marker=obs_marker,
                                           label='obspack',
                                           color=obs_color)
            # p1 *= obsplot
            #
            #-- update by Guillaume, 2025-09-09
            #
            p1 *= dfplot.hvplot.scatter(x='time', y='obspack',
                                        label='Observations (obspack)',
                                        color=obs_color,
                                        marker=obs_marker, s=obs_ms,
                                        legend=legend_loc)
        #
        #-- autohide toolbar
        #
        p1 = p1.opts(backend_opts={"plot.toolbar.autohide": True})

        return p1
