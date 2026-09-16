#!/usr/bin/env python

from omegaconf import OmegaConf, DictConfig
import requests
from pathlib import Path
import sys
from loguru import logger
import xarray as xr
from pandas import read_csv, DataFrame, concat, Timestamp
import numpy as np
from collections import OrderedDict
from glob import glob
import panel as pn
import param
import hvplot.xarray
from holoviews import opts, Overlay
from typing import Tuple, Dict, List
import io
from numpy import zeros
from itertools import cycle
from bokeh.palettes import Category10
import itertools
import geoviews.feature as gf
from cartopy import crs

from tm5 import debug
from tm5.gui.css import *
from tm5.gui.widgets.emissions import EmissionSettings
from tm5.gui.widgets.stations import calc_statistics
from tm5.gui.widgets.widget_utils import experiment_desc, plot_site_info, load_observations_metadata




def get_exp_label(expfile: str | Path) -> str:
    return Path(expfile).name.replace('fitic-', '').rsplit('_monthly')[0]


@debug.timer
def simulation_read_targets(output_path: Path) -> Tuple[np.ndarray, DataFrame]:
    """
    """
    #
    #-- read target identifier
    #
    tgt = xr.open_dataset(output_path / 'ftj.nc')
    tgt_dict = { 'target': tgt.targets.values,
                 'prior': [],
                 'posterior': [],
                 'posterior_uncertainty': []
                }
    #
    #-- output/t0.dat --> prior target
    #   output/t.dat  --> posterior target
    #   output/ct.dat --> posterior correlation of targets
    #                     (uncertainties are square-root of diagonal elements)
    #
    #-- prior
    #
    tfile = output_path / 'output' / 't0.dat'
    if not tfile.exists():
        msg = f"expected prior target file ***{tfile.name}*** not present"
        raise RuntimeError(msg)
    else:
        msg = f"reading from ***{tfile.name}***"
        logger.debug(msg)
        dft = read_csv(tfile, sep=r"\s+|\t+", engine='python')
        # dft.to_csv('tgt-apri.csv')
        # msg = f"prior targets -->{dft}<--"
        # logger.debug(msg)
        tgt_dict['prior'] = dft.loc[:,'col_1'].values
    #
    #-- posterior
    #
    tfile = output_path / 'output' / 't.dat'
    if not tfile.exists():
        msg = f"expected prior target file ***{tfile.name}*** not present"
        raise RuntimeError(msg)
    else:
        msg = f"reading from ***{tfile.name}***"
        logger.debug(msg)
        dft = read_csv(tfile, sep=r"\s+|\t+", engine='python')
        # dft.to_csv('tgt-apos.csv')
        # msg = f"prior targets -->{dft}<--"
        # logger.debug(msg)
        tgt_dict['posterior'] = dft.loc[:,'col_1'].values
    #
    #-- posterior uncertainty
    #
    tfile = output_path / 'output' / 'ct.dat'
    if not tfile.exists():
        msg = f"expected prior target file ***{tfile.name}*** not present"
        raise RuntimeError(msg)
    else:
        msg = f"reading from ***{tfile.name}***"
        logger.debug(msg)
        dft = read_csv(tfile, sep=r"\s+|\t+", engine='python')
        # dft.to_csv('tgtunc-apos.csv')
        for itgt, tgt in enumerate(tgt.targets):
            #-- correlation matrix (!), need to take square root of diagonal
            #-- MVO-ATTENTION:column indexing is Fortran based (col_1,col_2,col_3,...)
            _tgtunc = np.sqrt(dft.loc[itgt, f'col_{itgt+1}'])
            tgt_dict['posterior_uncertainty'].append(_tgtunc)
        # msg = f"prior targets -->{dft}<--"
        # logger.debug(msg)
    return DataFrame.from_dict(tgt_dict).set_index('target')


@debug.timer
def conc_statistics(conc: xr.Dataset, label: str) -> DataFrame:
    """
    """
    dfc = conc.to_dataframe()
    stations = sorted(list(set(dfc.station.values)))

    #-- mapping of cases onto identifiers in data frame
    case_table = OrderedDict()
    case_table['prior'] = 'apri'
    case_table['post'] = 'apos'
    # msg = f"stations -->{stations}<--"
    # logger.debug(stations)
    stats = OrderedDict()
    stats['station'] = []
    for case_name in case_table.keys():
        stats[f"Mean bias ({case_name})"] = []
    for case_name in case_table.keys():
        stats[f"RMSE ({case_name})"] = []
    for case_name in case_table.keys():
        stats[f"Correlation coefficient ({case_name})"] = []
        
    for sta in stations:
        # msg = f"now @station={sta}"
        # logger.debug(msg)
        #
        #-- select current station (for all days)
        #
        cnd = dfc['station']==sta
        _df = dfc.loc[cnd,:]
        _cobs  = _df.loc[:,'obs']
        stats['station'].append(sta)
        #-- stats for prior and posteror
        for case_name, case_tag in case_table.items():
            _csim = _df.loc[:,f"{case_tag}_{label}"]
            #
            #-- compute statistics
            #
            _bias  = _csim - _cobs
            _meanbias = _bias.mean()
            _rmse = (_bias ** 2).mean() ** .5
            _corrcoef = np.corrcoef(_csim.values,_cobs.values)[0,1]
            # msg = f"@{sta}/{case_name}/{label}: meanbias/rmse/corrcoef = {_meanbias}/{_rmse}/{_corrcoef}"
            # logger.debug(msg)
            #
            #-- insert statistics
            #
            stats[f'Mean bias ({case_name})'].append(_meanbias)
            stats[f'RMSE ({case_name})'].append(_rmse)
            stats[f'Correlation coefficient ({case_name})'].append(_corrcoef)
    # msg = f"...loop terminated, stats -->{stats}<--"
    # logger.info(msg)
    #
    #-- turn into dataframe
    #
    stats = DataFrame.from_dict(stats).set_index('station')
    #-- reorder columns
    
    # #-- DEBUGGING
    # stats.to_csv('stats.csv', index=True)

    return stats


def load_inversion_concentrations(path: Path, label: str) -> xr.Dataset:
    fc = xr.open_dataset(path / 'fcpost.nc')
    # Reformat the data as a dataframe, consistent with _conc_plot
    conc = fc[['obs', 'cprior', 'cpost', 'station', 'obstime']].to_dataframe()
    for station_id in fc.station_id.values:
        stat = fc.sel(nsta = fc.station_id == station_id)
        conc.loc[conc.station == station_id, 'station_lon'] = float(stat.station_lon.values[0])
        conc.loc[conc.station == station_id, 'station_lat'] = float(stat.station_lat.values[0])
    conc['time'] = [Timestamp(_) for _ in conc.loc[:,'obstime']]
    conc = conc.rename(columns={'cprior':f'apri_{label}', 'cpost':f'apos_{label}'})
    conc = conc.to_xarray()
    conc.attrs['units'] = fc.obs.units
    return conc


def load_forward_concentrations(path: Path, label: str) -> xr.Dataset:
    fc = xr.open_dataset(path / 'fc.nc')
    conc = fc[['obs', 'conc', 'station', 'obstime']].to_dataframe()
    conc.loc[:,'station'] = conc.loc[:,'station'].astype('string')
    for station_id in fc.station_id.values:
        stat = fc.sel(nsta = fc.station_id == station_id)
        conc.loc[conc.station == station_id, 'station_lon'] = float(stat.station_lon.values[0])
        conc.loc[conc.station == station_id, 'station_lat'] = float(stat.station_lat.values[0])
        conc.loc[conc.station == station_id, 'station_lon'] = float(stat.station_lon.values[0])
        conc.loc[conc.station == station_id, 'station_alt'] = float(stat.station_alt.values[0])
    conc['time'] = [Timestamp(_) for _ in conc.loc[:,'obstime']]
    conc = conc.rename(columns={'conc':f'forward_{label}'})
    conc = conc.to_xarray()
    conc.attrs['units'] = fc.obs.units
    return conc


def extract_emis_map(ds, region: str, varname : str) -> xr.DataArray:
    ds = ds.sel(ng = (ds.region == region))
    # #-- MVO-20260728::
    # #   fe.nc:      data variable is 'emissions'
    # #   fepost.nc:  data variables are 'emission_prior' *and* 'emission_post'
    # if 'emission' in ds:
    #     emis_in = ds['emission']
    # elif 'emission_post' in ds:
    #     emis_in = ds['emission_post']
    emis_in = ds[varname]
    assert emis_in.attrs['units'] in ['kgCH4/cell',"kgCH4/cell/month"], \
        f"unexpected emission units -->{emis_in.attrs['units']}<--"
    emis_out = (emis_in / ds.area).values
    ilat = ((ds.lat - ds.lat.min()) / int(region[7])).values.astype(int)
    ilon = ((ds.lon - ds.lon.min()) / int(region[3])).values.astype(int)
    emis = zeros((ilat.max() + 1, ilon.max() + 1))
    emis[ilat, ilon] = emis_out
    emis = xr.DataArray(data = emis,
                        dims=('lat', 'lon'),
                        attrs={'units':'kgCH4/m2'})
    emis['lat'] = ('lat', sorted(set(ds.lat.values)))
    emis['lon'] = ('lon', sorted(set(ds.lon.values)))
    # logger.debug(f"@{region}, emis={emis}")
    return emis


def load_emissions(path: Path) -> Dict[str, xr.Dataset]:
    apri = xr.open_dataset(path / 'fe.nc')
    apos = xr.open_dataset(path / 'fepost.nc')
    #-- deliberately selecting last month
    ilastmon = apri.nmon.values[-1]

    apri = apri.isel(nmon=ilastmon)
    apos = apos.isel(nmon=ilastmon)
    emis_mon = Timestamp(*apri.time.values)
    # logger.debug(f"emis_mon ->{emis_mon}<-")
    #-- meanwhile 'area' available in fepost.nc:
    if not 'area' in apos:
        # "area" seems missing in the posterior file ... copy it from the prior
        apos['area'] = apri['area']
    
    emis = {}
    for region in ['glb600x400', 'eur300x200', 'gns100x100']:
        emis[region] = xr.Dataset(
            dict(
                apri = extract_emis_map(apri, region, 'emission'), 
                apos = extract_emis_map(apos, region, 'emission_post')
            ),
            attrs={'emis_month':emis_mon.strftime("%b %Y")}
        )
    return emis


def plot_conc_timeseries(df: DataFrame, simul_type: str, cur_exp: str):
    p = df.hvplot.points(x='time', y='obs', grid=True, c='k', label='obs', width=1200, height=400)
    color_palette = itertools.cycle(Category10[10])
    # print(cur_exp, dfc.columns)

    # Find all "forward" experiments
    if simul_type == 'fwd':
        experiments = [c.split('_', maxsplit=1)[1] for c in df.columns if c.startswith('forward_')]
        # print(experiments)
        for iexp, exp in enumerate(experiments):
            col = next(color_palette) # Category10[10][iexp]
            # print(exp, iexp, col)
            if exp == cur_exp:
                p *= df.hvplot.line(x='time', y=f'forward_{exp}', c=col, label=exp, muted_alpha=0, line_width=4)
            else:
                p *= df.hvplot.line(x='time', y=f'forward_{exp}', c=col, label=exp, muted_alpha=0, line_width=2)

    # Find all "inversion" experiments
    elif simul_type == 'inv':
        experiments = [c[5:] for c in df.columns if c.startswith('apri_')]
        # Plot the experiments:
        for iexp, exp in enumerate(experiments):
            col = next(color_palette) # Category10[10][iexp]
            if exp == cur_exp:
                p *= df.hvplot.line(x='time', y=f'apri_{exp}', c=col, line_dash='dashed', label=f'prior_{exp}', line_width=4, muted_alpha=0)
                p *= df.hvplot.line(x='time', y=f'apos_{exp}', c=col, label=f'posterior_{exp}', line_width=4, muted_alpha=0)
            else:
                p *= df.hvplot.line(x='time', y=f'apri_{exp}', c=col, line_dash='dashed', label=f'prior_{exp}', line_width=2, muted_alpha=0)
                p *= df.hvplot.line(x='time', y=f'apos_{exp}', c=col, label=f'posterior_{exp}', line_width=2, muted_alpha=0)
    return p


def plot_stats_table(df: DataFrame, emis_dataset: str):
    #-- reduce decimals for visualisation
    df = df.round(decimals=2)
    # msg = f"df.columns ==>{df.columns}<=="
    # logger.debug(msg)
    # columns = [
    #     {
    #         "title": col,
    #         "field": col,
    #         "headerSort": col != "station",
    #     }
    #     for col in df.columns
    # ]
    # msg = f"columns ***{columns}***"
    # logger.debug(msg)
    formatters = {
        col: {"type": "number", "precision": 2}
        for col in df.columns
        if df[col].dtype.kind in "fiu"
    }
    # msg = f"formatters ***{formatters}***"
    # logger.debug(msg)
    
    table = pn.widgets.Tabulator(
        df,
        text_align={'station':"left"},
        # text_align="center",
        # columns=columns,
        # formatters=formatters,
    )
    # msg = f"tabulator widget generated table ***{table}***"
    # logger.debug(msg)
    title = f"# Fit statistics for all stations (emissions scenario: {get_exp_label(emis_dataset)})"
    return pn.Column(pn.pane.Markdown(title), table)

# def plot_stats_table(df: DataFrame, emis_dataset: str):
#     nc = len(df.columns)
#     formatters = [lambda x: f'{x:.2f}'] * nc
#     table = pn.pane.DataFrame(df, text_align='center', formatters=formatters)
#     title = f'# Fit statistics for all stations ({get_exp_label(emis_dataset)})'
#     return pn.Column(pn.pane.Markdown(title), table)


def plot_emis_table_md(emis_datasets: List[str]):
    lines = ['| **Emissions setup** | **Description** |']
    lines.append('| --- | --- |')
    for exp in emis_datasets:
        exp = get_exp_label(exp)
        if desc := experiment_desc(exp):
            lines.append(f'| {exp} | {desc} |')
    return '\n'.join(lines)


def plot_emission_map(emissions: xr.Dataset, emis_dataset: str):
    # print("computing emission map")
    logger.debug("computing emission map")
            
    cmap = 'RdBu_r'
    clim = (-0.005, 0.005)

    projection = crs.RotatedPole(pole_longitude=185, pole_latitude=50)
    xlim = (-15, 35)
    ylim = (33, 73)

    projections = {
        'glb600x400': crs.PlateCarree(),
        'eur300x200': crs.PlateCarree(),
        'gns100x100': crs.PlateCarree()
    }

    xlims = {
        'glb600x400': (-180, 180),
        'eur300x200': (-36, 54),
        'gns100x100': (0, 18)
    }

    ylims = {
        'glb600x400': (-90, 90),
        'eur300x200': (22, 74),
        'gns100x100': (42, 58)
    }
        
    logger.debug("...returning emissions map now!!!")
    
    emis_units = emissions['glb600x400'].apos.attrs['units']
    emis_month = emissions['glb600x400'].attrs['emis_month']

    logger.debug(f"emis_units ==>{emis_units}<==")

    mode = 'guillaume'
    if mode in ['glb100x100', 'eur300x200', 'gns100x100']:
        prow = pn.Row(
            emissions[mode].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projections[mode], xlim=xlims[mode], ylim=ylims[mode]),
            emissions[mode].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projections[mode], xlim=xlims[mode], ylim=ylims[mode]),
            (emissions[mode].apos - emissions[mode].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection[mode], xlim=xlims[mode], ylim=ylims[mode])
        )
        title = f"# emission maps (posterior, prior, posterior-prior) ({get_exp_label(emis_dataset)})"
        p = pn.Column(pn.pane.Markdown(title), prow)
    elif mode == 'row_merged':
        prow = pn.Row(
            emissions['glb600x400'].apos.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=projections['glb600x400'], xlim=xlims['glb600x400'], ylim=ylims['glb600x400']) *
            emissions['eur300x200'].apos.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim) *
            emissions['gns100x100'].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim),
            emissions['glb600x400'].apri.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=projections['glb600x400'], xlim=xlim['glb600x400'], ylim=ylim['glb600x400']) *
            emissions['eur300x200'].apri.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim) *
            emissions['gns100x100'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim),
            (emissions['glb600x400'].apos - emissions['glb600x400'].apri).hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=projections['glb600x400'], xlim=xlims['glb600x400'], ylim=ylims['glb600x400']) *
            (emissions['eur300x200'].apos - emissions['eur300x200'].apri).hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim) *
            (emissions['gns100x100'].apos - emissions['gns100x100'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim)
        )
        title = f"# emssions maps (posterior, prior, posterior-prior) ({get_exp_label(emis_dataset)})"
        p = pn.Column(pn.pane.Markdown(title), prow)
    elif mode == 'guillaume':
        ppost = (
            emissions['glb600x400'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['eur300x200'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['gns100x100'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
        )
        plotcfg = opts.Overlay(title=f'posterior emissions ({emis_month}, {get_exp_label(emis_dataset)})', ylabel=f"[{emis_units}]" )
        ppost.opts(plotcfg)
        pprior = (
            emissions['glb600x400'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['eur300x200'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['gns100x100'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
            )
        plotcfg = opts.Overlay(title=f'prior emissions ({emis_month}, {get_exp_label(emis_dataset)})' )
        pprior.opts(plotcfg)
        pdiff = (
            (emissions['glb600x400'].apos - emissions['glb600x400'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            (emissions['eur300x200'].apos - emissions['eur300x200'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            (emissions['gns100x100'].apos - emissions['gns100x100'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
            )
        plotcfg = opts.Overlay(title=f'posterior-prior emissions ({emis_month}, {get_exp_label(emis_dataset)})')
        pdiff.opts(plotcfg)
        
        prow = (
            ppost + pprior + pdiff
        )
        p = prow
    return p
    

def plot_target_table(df: DataFrame, emis_dataset: str):
    def get_csv_file():
        # You can re-compute data here if it's dynamic
        fid = io.BytesIO()
        fid.write(df.to_csv(index=False).encode('utf-8'))
        fid.seek(0)
        return fid

    nc = len(df.columns)
    formatters = [lambda x: f'{x:.3f}'] * nc
    # outname = f"targets_{get_exp_label(self.emis_dataset)}.csv"
    # outname = f"targets.csv"
    # button = pn.widgets.FileDownload(callback=get_csv_file, filename=outname, label="Download Data (CSV)", button_type="primary")
    button = pn.widgets.FileDownload(callback=get_csv_file, filename="targets.csv", label="Download Data (CSV)", button_type="primary")

    p = pn.Column(pn.pane.DataFrame(df, text_align='center', formatters=formatters), button)
    # return self.tgt_table
    #-- MVO-TODO::units [MtCH4] should not be hard-coded here
    logger.debug("...returning target table now")
    title = f'# Target emission quantities [MtCH4] ({get_exp_label(emis_dataset)})'
    return pn.Column(pn.pane.Markdown(title), p)


def plot_map_sites(df: DataFrame, current_site: str):
    # msg = f"...@{current_site}"
    # logger.debug(msg)
    logger.debug(df.head())
    df.loc[:, 'cur_site'] = 0
    df.loc[df.station == current_site, 'cur_site'] = 1
    # msg = f"df.columns -->{df.columns}<--"
    # logger.debug(msg)
    lonmin = df.loc[:,'station_lon'].values.min()
    lonmax = df.loc[:,'station_lon'].values.max()
    latmin = df.loc[:,'station_lat'].values.min()
    latmax = df.loc[:,'station_lat'].values.max()
    # msg = f"lonmin/lonmax = {lonmin}/{lonmax}, latmin/latmax = {latmin}/{latmax}"
    # logger.debug(msg)
    #
    cnd_gns = (0<=lonmin<=lonmax<=18 ) and (42<=latmin<=latmax<=58)
    if cnd_gns:
        xlim = (-15, 35)
        ylim = (33, 73)
        tiles = 'EsriTerrain'
    else:
        xlim = (-180,180)
        ylim = (-90,90)
        #-- MVO-TODO:global map does not show up when using 'EsriTerrain'
        tiles = None
    #
    if 'station_alt' in df.columns:
        hover_cols = ['station','station_alt',]
    else:
        hover_cols = ['station',]
    # msg = f"...@{current_site}, xlim/ylim = {xlim}/{ylim}"
    # logger.debug(msg)
    emis_map = df.hvplot.points(
        x='station_lon', y='station_lat', color='cur_site', cmap=['LightSlateGray', 'red'],
        geo=True, coastline=True,
        hover_cols=hover_cols,
        xlim=xlim, ylim=ylim, colorbar=False, tiles=tiles
    )
    return emis_map
    #-- widgets-boarders currently make problems
    #   on exploredata.icos-cp.eu,
    #   disabled for NCGG10
    # return df.hvplot.points(
    #     x='station_lon', y='station_lat', color='cur_site', cmap=['LightSlateGray', 'red'],
    #     geo=True, coastline=True, xlim=xlim, ylim=ylim, colorbar=False, tiles=tiles
    # ) * self.widgets['borders']


class PreconfExperimentGUI(pn.viewable.Viewer):
    emis_dataset = param.FileSelector(doc='Prior emission dataset')
    run_forward = param.Event(doc='Do a forward run', label='Perform a forward simulation')
    run_inv = param.Event(doc='Do an inversion', label='Perform an inversion')
    alert = param.String(doc='Generic object for error messages or others ...', default='')
    current_site = param.Selector(doc='Current site to be displayed', default=None)
    sites_list = param.List(default=[], doc='List of observation sites available (for internal use ...)')
    simul_type = param.Selector(objects=['fwd', 'inv'], allow_None=True, default=None)
    correlation_switch = param.Boolean(doc='Switch to enable/disable correlated emission adjustments', default=False, label='Spatially correlated prior emissions uncertainty')
    add_category_event = param.Event(doc='Add a new emission category', label='Add category')
    preconf_scenario = param.Selector(default=None, allow_None=True, doc='Preconfigured emission scenario')

    # Data containers:
    conc        = param.ClassSelector(class_=xr.Dataset)
    stats4conc  = param.DataFrame()
    tgt_table   = param.DataFrame()
    emissions   = param.Dict()

    def __init__(self, gui_settings: DictConfig):

        super().__init__()

        self.cache_fwd = OrderedDict()
        self.cache_inv = OrderedDict()

        self._message = ''
        self.gui_settings = gui_settings

        # Load the file list
        self.param.emis_dataset.path = self.gui_settings.emissions.glob_pattern
        self.emis_dataset = self.param.emis_dataset.objects[0]

        self.emission_scenario = []
        self.emission_scenario_widgets = pn.Column()

        scenarios = self.gui_settings.emissions.get('scenarios', {})
        self.param.preconf_scenario.objects = {v['title']: k for k, v in scenarios.items()}
        self.preconf_scenario = 'default'

        # Globally accessible widgets
        self.widgets = {
            'station_selector': pn.widgets.Select.from_param(self.param.current_site),
            'borders': gf.borders()
        }
        self.widgets['station_selector'].visible = False

    def __panel__(self):
        header_pane = pn.pane.Markdown('# Preconfigured prior emission scenarios')
        expdesc_pane = pn.pane.Markdown(
            plot_emis_table_md(self.param.emis_dataset.objects),
            stylesheets=[preconfsim_stylesheet], 
            css_classes=['precomp-right']
        )

        widgets = [
            header_pane,
            pn.Row(pn.widgets.Select.from_param(self.param.emis_dataset), expdesc_pane),
            pn.Column(
                pn.Row(
                    pn.widgets.Select.from_param(self.param.preconf_scenario, name='Preconfigured scenario'), 
                    self.scenario_description
                ),
                self.emission_scenario_widgets,
                pn.widgets.Button.from_param(self.param.add_category_event)
            ),
            pn.Row(
                    pn.widgets.Button.from_param(self.param.run_forward),
                    pn.Column(
                        pn.widgets.Button.from_param(self.param.run_inv),
                        pn.widgets.Switch.from_param(self.param.correlation_switch)
                    ),
            ),
            self._alert,
            self.widgets['station_selector'],
            pn.Row(self.conc_plot, self.map_sites),
            pn.Row(self.conc_stats_table, self.target_table)
        ]
        if self.gui_settings.get('show_emismap', False):
            widgets.append(self.map_emissions)
        return pn.Column(*widgets)

    # ------------------------------------------------
    # Interactive panels/widgets

    @param.depends('run_forward', watch=True)
    def _run_forward(self):
        output_path = self._call_backend('forward')
        self.simul_type = 'fwd'
        if output_path is not None:
            self._read_concentrations(output_path, 'forward')
            self.stats4conc = None
            self.tgt_table = None
            self.emissions = None

    @param.depends('run_inv', watch=True)
    def _run_inv(self):
        output_path = self._call_backend('inversion')
        # msg = f"...output_path ***{output_path}*** ('output_path is not None': {output_path is not None})"
        # logger.debug(msg)
        self.simul_type = 'inv'
        if output_path is not None:
            # msg = f"...start reading concentrations"
            # logger.debug(msg)
            self._read_concentrations(output_path, 'inversion')
            # msg = f"...start reading emissions from -->{output_path}<--"
            # logger.debug(msg)
            self.emissions = load_emissions(output_path)
            # msg = f"...computing conc statistics"
            # logger.debug(msg)
            self.stats4conc = conc_statistics(self.conc, get_exp_label(self.emis_dataset))
            # msg = f"...reading simulation targets from directory ***{str(output_path)}***"
            # logger.debug(msg)
            self.tgt_table = simulation_read_targets(output_path)

    @param.depends('add_category_event', watch=True)
    def _add_emission_category(self, catname: str = None):
        if catname is None:
            catname = f'category_{len(self.emission_scenario) + 1}'
        es = EmissionSettings(
            catname=catname,
            regions=['global', 'regional'],
            path=self.gui_settings.emissions.path,
            remove_callback=self._remove_emission_category,
            locked=(self.preconf_scenario != 'custom'),
        )
        self.emission_scenario.append(es)
        self.emission_scenario_widgets.append(es.__panel__())

    def _remove_emission_category(self, es: EmissionSettings):
        self.emission_scenario.remove(es)
        self.emission_scenario_widgets.objects = [e.__panel__() for e in self.emission_scenario]

    @param.depends('preconf_scenario', watch=True)
    def _load_scenario(self):
        if self.preconf_scenario is None:
            return
        scenario = self.gui_settings.emissions.scenarios[self.preconf_scenario]
        self.emission_scenario = []
        for catname, spec in scenario.get('categories', {}).items():
            self._add_emission_category(catname)
            self.emission_scenario[-1].set_category(spec)
        self.emission_scenario_widgets.objects = [e.__panel__() for e in self.emission_scenario]

    @param.depends('preconf_scenario')
    def scenario_description(self):
        if self.preconf_scenario is None:
            return ''
        return pn.pane.Markdown(
            self.gui_settings.emissions.scenarios[self.preconf_scenario].get('description', ''),
            stylesheets=[preconfsim_stylesheet], css_classes=['precomp-right']
        )

    @param.depends('alert')
    def _alert(self):
        if self.alert == '':
            return ''
        return pn.pane.Alert(self.alert, alert_type='danger')

    @param.depends('conc', 'current_site')
    def conc_plot(self):
        # msg = f"...current_site={self.current_site}"
        # logger.debug(msg)
        if self.conc is None:
            return ''
        if self.current_site is None:
            return ''
        cur_exp = get_exp_label(self.emis_dataset)
        dfc = self.conc.to_dataframe()
        dfc = dfc[dfc.station == self.current_site]
        # msg = f"...@{self.simul_type},cur_exp={cur_exp}: calling plot_conc_timeseries..."
        # logger.debug(msg)
        return plot_conc_timeseries(dfc, self.simul_type, cur_exp)

    @param.depends('stats4conc')
    def conc_stats_table(self):
        if self.stats4conc is None:
            return ''
        return plot_stats_table(self.stats4conc, self.emis_dataset)

    @param.depends('emissions')
    def map_emissions(self):
        if self.emissions is None:
            # print("resetting emission map")
            logger.debug("resetting emission map")
            return 
        return plot_emission_map(self.emissions, self.emis_dataset)

    @param.depends('tgt_table')
    def target_table(self):
        if self.tgt_table is None:
            logger.debug("returning None")
            return ''
        return plot_target_table(self.tgt_table, self.emis_dataset)

    @param.depends('current_site', 'sites_list')
    def map_sites(self):
        # msg = f"...current_site={self.current_site}"
        # logger.debug(msg)
        if self.conc is None or self.current_site is None:
            return ''
        # msg = f"...calling plot_map_sites"
        # logger.debug(msg)
        site_columns = ['station', 'station_lon', 'station_lat','station_alt']
        site_columns = ['station', 'station_lon', 'station_lat',]
        site_map = plot_map_sites(
            self.conc.to_dataframe().loc[:, site_columns].drop_duplicates(),
            self.current_site
        )
        # msg = f"...returning site_map  ==>{site_map}<=="
        # logger.debug(msg)
        return site_map 

    # ------------------------------------------------
    # Internal methods (communication with the backend)

    def _call_backend(self, task: str) -> Path | None :
        """
        Triggers an inversion on the VM, and return either the output path, or an error code:
            - 100: something went wrong ...
            - 101: result is not valid json
        """
        # if self.emis_dataset in self.cache_inv and task == 'inversion':
        #     return self.cache_inv[self.emis_dataset]
        # elif self.emis_dataset in self.cache_fwd and task == 'forward':
        #     return self.cache_fwd[self.emis_dataset]

        url = f"{self.gui_settings.backend_url}/forward"

        settings = {
            'emis': self.emis_dataset,
            'emissions': {
                es.catname: {
                    'global': {
                        'file': es.emis_glo.filename, 
                        'field': es.emis_glo.fieldname},
                    **({'regional': {
                        'file': es.emis_reg.filename, 
                        'field': es.emis_reg.fieldname}} if es.switch_reg else {}
                    )
                }
                for es in self.emission_scenario
            },
            'task': task,
            'namelist': {
                'fix': self.correlation_switch
            }
        }
        r = requests.post(url, data={'conf': OmegaConf.to_yaml(settings)})
        if not r.ok:
            self.alert = f"{task} run failed: backend returned {r.status_code} for emis={self.emis_dataset} at {url}. Body: {r.text[:500]}"
            return
        try:
            payload = r.json()
        except requests.exceptions.JSONDecodeError:
            self.alert = f"{task} run failed: backend returned non-JSON for emis={self.emis_dataset} at {url}. Body: {r.text[:500]}"
            return
        self.alert = ''

        output_path = Path(payload['output'])
        # msg = f"@task={task} for {self.emis_dataset} yields output_path ***{str(output_path)}***"
        # logger.debug(msg)
        if task == 'inversion':
            self.cache_inv[self.emis_dataset] = output_path
        else:
            self.cache_fwd[self.emis_dataset] = output_path
        return output_path
        
    def _read_concentrations(self, path: Path, task: str):
        label = get_exp_label(self.emis_dataset)
        # msg = f"@task={task}, emissions_label -->{label}<--"
        # logger.debug(msg)
        if task == 'inversion':
            conc = load_inversion_concentrations(path, label)
        else:
            conc = load_forward_concentrations(path, label)
        # msg = f"...@task={task}, loading concentrations done."
        # logger.debug(msg)
        if self.conc is None:
            conc_update = conc
        else:
            conc_update = xr.merge([self.conc, conc], compat='override')

        # Now update the "sites_list", if needed:
        sites_available = set(conc_update.station.values.reshape(-1))
        update_sites = (sites_available != set(self.sites_list))
        # msg = f"@task={task}, emissions_label -->{label}<-- " \
        #     f"sites_available -->{sites_available}<--, update_sites={update_sites}"
        # logger.info(msg)
        if update_sites:
            updated_site_list = sorted(list(sites_available))
            # msg = f"...self.current_site before change ***{self.current_site}***"
            # logger.debug(msg)
            self.param.current_site.objects = updated_site_list
            # msg = f"...self.param.current_site.objects now set"
            # logger.debug(msg)
            self.current_site = self.param.current_site.objects[0]
            # msg = f"...self.current_site now set to -->{self.current_site}<--"
            # logger.debug(msg)
            self.widgets['station_selector'].visible = True
            # msg = f"...self.widgets['station_selector'].visible now set " \
            #     f"value={self.widgets['station_selector'].visible}"
            # logger.debug(msg)
            #-- MVO-NOTE::setting self.conc seems to slightly improve GUI performance
            #             (Note that self.conc triggers member routines conc_plot/map_sites to
            #              be executed...)
            self.conc = conc_update
            # msg = f"...updated_site_list ==>{updated_site_list}<=="
            # logger.debug(msg)
            self.sites_list = updated_site_list
            # msg = f"...self.sites_list now set!"
            # logger.debug(msg)
        else:
            self.conc = conc_update
