#!/usr/bin/env python

import asyncio
import copy
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
from tm5.util import get_dict_checksum as dict_checksum
from tm5.gui.css import *
from tm5.gui.widgets.emissions import EmissionSettings
from tm5.gui.widgets.stations import calc_statistics
from tm5.gui.widgets.widget_utils import experiment_desc, plot_site_info, load_observations_metadata


simu_stylesheets = [
    """
    :host(.forward-button) .bk-btn {
    background-color: #cfd9e2 !important;
    border-color: #cfd9e2 !important;
    color: white;
    /*border: none;*/
    font-size: 1.25em;
    border: 1px solid #388E3C;
    border-radius: 6px;
    }
    :host(.forward-button) .bk-btn:hover {
    background-color: #00cc66 !important;
    border-color: #00cc66 !important;
    }
    :host(.inversion-button) .bk-btn {
    background-color: #cfd9e2 !important;
    border-color: #cfd9e2 !important;
    color: white;
    border: none;
    font-size: 1.25em;
    border: 1px solid #388E3C;
    border-radius: 6px;
    }
    :host(.inversion-button) .bk-btn:hover {
    background-color: #00cc66 !important;
    border-color: #00cc66 !important;
    }
    :host(.duplicate-button) .bk-btn {
    background-color: #c8d9ea !important;
    border-color: #c8d9ea !important;
    color: white;
    border: none;
    font-size: 1.25em;
    border: 1px solid #388E3C;
    border-radius: 6px;
    }
    :host(.duplicate-button) .bk-btn:hover {
    background-color: #00cc66 !important;
    border-color: #00cc66 !important;
    }
    :host(.category-button) .bk-btn {
    background-color: #c8d9ea !important;
    border-color: #c8d9ea !important;
    color: white;
    border: none;
    font-size: 1.25em;
    border: 1px solid #388E3C;
    border-radius: 6px;
    }
    :host(.category-button) .bk-btn:hover {
    background-color: #00cc66 !important;
    border-color: #00cc66 !important;
    }
    :host(.simu-select) .bk-input {
/*    background-color: #f0f8ff;*/
    color: #222;
    border: 1px solid #1976d2;
    border: 1px solid #000000;
    border-radius: 6px;
    font-weight: 500;
    font-size: 1.25em;
    }
    :host(.simu-select) .bk-input:focus {
    border-color: #0d47a1;
    box-shadow: 0 0 0 2px rgba(25, 118, 210, 0.2);
    }
    :host(.simu-select) .bk-input-group > label {
    font-size: 1.25em;
    font-weight: 600;
    color: #444;
    }
    """]

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
    # for case_name in case_table.keys():
    #     stats[f"Correlation coefficient ({case_name})"] = []
        
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
            # stats[f'Correlation coefficient ({case_name})'].append(_corrcoef)
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


def plot_stats_table(df: DataFrame, scenario_name: str):
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
    title = f"# Fit statistics for all stations (emissions scenario: {scenario_name})"
    return pn.Column(pn.pane.Markdown(title), table)

# def plot_stats_table(df: DataFrame, scenario_name: str):
#     nc = len(df.columns)
#     formatters = [lambda x: f'{x:.2f}'] * nc
#     table = pn.pane.DataFrame(df, text_align='center', formatters=formatters)
#     title = f'# Fit statistics for all stations ({scenario_name})'
#     return pn.Column(pn.pane.Markdown(title), table)


def plot_emis_table_md(emis_datasets: List[str]):
    lines = ['| **Emissions setup** | **Description** |']
    lines.append('| --- | --- |')
    for exp in emis_datasets:
        if desc := experiment_desc(exp):
            lines.append(f'| {exp} | {desc} |')
    return '\n'.join(lines)


def plot_scenario_table_md(scenarios: DictConfig):
    lines = ['| **Scenario** | **Description** |']
    lines.append('| --- | --- |')
    for spec in scenarios.values():
        title = spec.get('title', '')
        desc = spec.get('description', '')
        editable = spec.get('editable',False)
        if editable:
            title = f'<div style="background-color: #e8f5e9;">{title}</div>'
        else:
            title = f'<div style="background-color: #F2D2B8;">{title}</div>'
        cur_line = f"| {title} | {desc} |"
        lines.append(cur_line)
    return '\n'.join(lines)


def plot_emission_map(emissions: xr.Dataset, emis_dataset: str):
    # print("computing emission map")
    msg = f"computing emission map for emis_dataset -->{emis_dataset}<--"
    logger.debug(msg)
            
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
        
    emis_units = emissions['glb600x400'].apos.attrs['units']
    emis_month = emissions['glb600x400'].attrs['emis_month']

    msg = f"detected emis_units={emis_units}, emis_month={emis_month}"
    logger.debug(msg)
    
    mode = 'guillaume'
    if mode in ['glb100x100', 'eur300x200', 'gns100x100']:
        prow = pn.Row(
            emissions[mode].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projections[mode], xlim=xlims[mode], ylim=ylims[mode]),
            emissions[mode].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projections[mode], xlim=xlims[mode], ylim=ylims[mode]),
            (emissions[mode].apos - emissions[mode].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection[mode], xlim=xlims[mode], ylim=ylims[mode])
        )
        title = f"# emission maps (posterior, prior, posterior-prior) ({emis_dataset})"
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
        title = f"# emssions maps (posterior, prior, posterior-prior) ({emis_dataset})"
        p = pn.Column(pn.pane.Markdown(title), prow)
    elif mode == 'guillaume':
        #
        #-- posterior
        #
        ppost = (
            emissions['glb600x400'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['eur300x200'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['gns100x100'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
        )
        title = f'posterior emissions ({emis_month}, {emis_dataset})'
        plotcfg = opts.Overlay(title=title, ylabel=f"[{emis_units}]" )
        ppost.opts(plotcfg)
        msg = f"...hvplot for posterior done"
        logger.debug(msg)
        #
        #-- prior
        #
        pprior = (
            emissions['glb600x400'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['eur300x200'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            emissions['gns100x100'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
            )
        title = f'prior emissions ({emis_month}, {emis_dataset})'
        plotcfg = opts.Overlay(title=title)
        pprior.opts(plotcfg)
        msg = f"...hvplot for prior done"
        logger.debug(msg)
        #
        #-- posterior - prior
        #
        pdiff = (
            (emissions['glb600x400'].apos - emissions['glb600x400'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            (emissions['eur300x200'].apos - emissions['eur300x200'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
            (emissions['gns100x100'].apos - emissions['gns100x100'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
            )
        title = title=f'posterior-prior emissions ({emis_month}, {emis_dataset})'
        plotcfg = opts.Overlay(title)
        pdiff.opts(plotcfg)
        msg = f"...hvplot for differences done"
        logger.debug(msg)
        #
        #-- combining hvplots
        #
        prow = (
            ppost + pprior + pdiff
        )
        p = prow
    msg = f"...returning emissions map for mode -->{mode}<-- (title ==>{title}<==)"
    logger.debug(msg)
    
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
    # outname = f"targets_{emis_dataset}.csv"
    # outname = f"targets.csv"
    # button = pn.widgets.FileDownload(callback=get_csv_file, filename=outname, label="Download Data (CSV)", button_type="primary")
    button = pn.widgets.FileDownload(callback=get_csv_file, filename="targets.csv", label="Download Data (CSV)", button_type="primary")

    p = pn.Column(pn.pane.DataFrame(df, text_align='center', formatters=formatters), button)
    # return self.tgt_table
    #-- MVO-TODO::units [MtCH4] should not be hard-coded here
    logger.debug("...returning target table now")
    title = f'# Target emission quantities [MtCH4] ({emis_dataset})'
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
    # emis_dataset = param.FileSelector(doc='Prior emission dataset')
    # run_forward = param.Event(doc='Do a forward run', label='Perform a forward simulation')
    # run_inv = param.Event(doc='Do an inversion', label='Perform an inversion')
    run_forward = param.Event(doc='', label='Perform a forward simulation')
    run_inv = param.Event(doc='', label='Perform an inversion')
    duplicate_scenario_event = param.Event(doc='Create a new scenario based on an existing one.', label='Duplicate scenario')
    confirm_duplicate_event = param.Event(doc='Click to create the new scenario under that name', label='Create')
    new_scenario_name = param.String(default='', doc='Name for the new scenario')
    alert = param.String(doc='Generic object for error messages or others ...', default='')
    # current_site = param.Selector(doc='Current site to be displayed', default=None)
    current_site = param.Selector(doc='d', default=None)
    sites_list = param.List(default=[], doc='List of observation sites available (for internal use ...)')
    simul_type = param.Selector(objects=['fwd', 'inv'], allow_None=True, default=None)
    correlation_switch = param.Selector(
        default="full grid",
        objects=["fixed patterns", "full grid"],
        label="Resolution of Emission space (Please note that option 'fixed patterns' is not implemented yet)",
    )
    add_category_event = param.Event(doc='Add a new emission category', label='Add category')
    hide_categories_event = param.Event(doc='Hide emissione categories', label='Hide categories')
    show_categories_event = param.Event(doc='Show emissione categories', label='Show categories')
    select_scenario = param.Selector(default=None, allow_None=True, doc='Selection or configurationtion of emission scenario')
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
        # self.param.emis_dataset.path = self.gui_settings.emissions.glob_pattern
        # self.emis_dataset = self.param.emis_dataset.objects[0]

        self.emission_scenario_widgets = pn.Column()

        # Copy the settings since this may get edited (I think it's cleaner to keep 
        # settings read-only).
        self.emission_scenarios = copy.deepcopy(self.gui_settings.emissions.scenarios)

        # Globally accessible widgets
        self.widgets = {
            'station_selector': pn.widgets.Select.from_param(self.param.current_site,
                                                             css_classes=["simu-select"],
                                                             stylesheets=simu_stylesheets,),
            'borders': gf.borders(),
            'add_category': pn.widgets.Button.from_param(self.param.add_category_event),
            'hide_categories': pn.widgets.Button.from_param(self.param.hide_categories_event,
                                                            css_classes=["category-button"],
                                                            stylesheets=simu_stylesheets,),
            'show_categories': pn.widgets.Button.from_param(self.param.show_categories_event,
                                                            css_classes=["category-button"],
                                                            stylesheets=simu_stylesheets,),
            'select_scenario': pn.widgets.Select.from_param(self.param.select_scenario,
                                                            name='Select scenario',
                                                            css_classes=["simu-select"],
                                                            stylesheets=simu_stylesheets,),
        }
        self.widgets['station_selector'].visible = False
        self.widgets['duplicate_prompt'] = pn.Row(
            pn.widgets.TextInput.from_param(self.param.new_scenario_name, name='New scenario name',
                                            css_classes=["simu-select"],
                                            stylesheets=simu_stylesheets,),
            pn.widgets.Button.from_param(self.param.confirm_duplicate_event,
                                         css_classes=["duplicate-button"],
                                         stylesheets=simu_stylesheets,),
            visible=False,
        )

        self.param.select_scenario.objects = list(self.emission_scenarios.keys())
        self._fix_scenario_display_name()
        self.select_scenario = 'default'

    def _fix_scenario_display_name(self):
        """
        Ensure that the "preconfigured scenarios" menu items reflect the "title"
        key of the scenarios (from the yaml file), and not the section name.
        """
        self.widgets['select_scenario'].options = {v['title']: k for k, v in self.emission_scenarios.items()}

    def __panel__(self):
        #        #--
        #
        intro_text =  """
        ## Introduction
         You are running a fast demo configuration of the Flexible Inversion Tool for Inventory Compilers (FIT-IC) with a focus on central Europe and for January 2021.<br>
        This demo allows you<br>
        <ul>
        <li>to select preconfigured emission scenarios (marked in orange below)</li>
        <li>to modify a preconfigured emission scenario (using the duplicate button)</li>
        <li>to create own emission scenarios starting from scratch or from FIT-IC default (marked in green below)</li>
        <li>to upload a prepared user defined emission scenario (by switchting to the 'upload emissions tab')</li>
        </ul>
        and to perform
        <ul>
        <li>a forward simulation based on the selected scenario(s) and compare the simulated atmospheric signal(s) to observed methane concentrations</li>
        <li>perform an atmospheric transport inversion using the selected scenario as prior emission field.</li>
        </ul>
        <br>
         For further background on the tool see <a href="https://fit-ic.inversion-lab.com">FIT-IC website</a>.
         """
        intro_pane = pn.pane.Markdown(
            intro_text,
            stylesheets=[preconfsim_stylesheet], 
            css_classes=['precomp-intro']
        )
        scenario_table_md = plot_scenario_table_md(self.emission_scenarios)
        scenario_table_md = f"## Description of prior emission scenarios\n{scenario_table_md}"
        scenario_table_pane = pn.pane.Markdown(
            scenario_table_md,
            stylesheets=[preconfsim_stylesheet],
            css_classes=['precomp-right']
        )

        widgets = [
            intro_pane,
            # header_pane,
            scenario_table_pane,
            pn.Column(
                pn.Column(
                    pn.pane.Markdown("## Select or configure prior emission scenario"),
                    pn.Row(
                        
                    self.widgets['select_scenario'],
                    pn.widgets.Button.from_param(self.param.duplicate_scenario_event,
                                                 css_classes=["duplicate-button"],
                                                 stylesheets=simu_stylesheets,),
                    self.widgets['duplicate_prompt'],
                ),
                self.emission_scenario_widgets,
                pn.Row(
                    self.widgets['add_category'],
                    self.widgets['hide_categories'],
                    self.widgets['show_categories'],
                    )
                ),
                sizing_mode="stretch_width",
                styles={
                    "min-width": "0",
                    "background": "#f5f7fa",
                    "border": "2px solid #ddd",
                    "border-radius": "8px",
                    "max-width": "75%",
                    "padding": "15px",
                }
            ),
            pn.Column(
                pn.pane.Markdown("## Running a prior emission scenario"),
                pn.Row(
                    pn.widgets.Button.from_param(self.param.run_forward,
                                                 css_classes=["forward-button"],
                                                 stylesheets=simu_stylesheets,
                                                 ),
                    pn.Column(
                        pn.widgets.Button.from_param(self.param.run_inv,
                                                     css_classes=["inversion-button"],
                                                     stylesheets=simu_stylesheets),
                        pn.widgets.Select.from_param(self.param.correlation_switch,
                                                     css_classes=["simu-select"],
                                                     stylesheets=simu_stylesheets),
                    ),
                    # sizing_mode="stretch_width",
                    # styles={
                    #     "min-width": "0",
                    #     "background": "#f5f7fa",
                    #     "border": "2px solid #ddd",
                    #     "border-radius": "8px",
                    #     "max-width": "75%",
                    #     "padding": "15px",
                    # }
                ),
                sizing_mode="stretch_width",
                styles={
                    "min-width": "0",
                    "background": "#f5f7fa",
                    "border": "2px solid #ddd",
                    "border-radius": "8px",
                    "max-width": "75%",
                    "padding": "15px",
              }  
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
            self.stats4conc = conc_statistics(self.conc, self.select_scenario)
            # msg = f"...reading simulation targets from directory ***{str(output_path)}***"
            # logger.debug(msg)
            self.tgt_table = simulation_read_targets(output_path)

    def _build_emission_category(self, catname: str, visible : bool|None = None) -> EmissionSettings:
        if visible==None:
            visible = self.emission_scenarios[self.select_scenario].get('editable', False)
        return EmissionSettings(
            catname=catname,
            regions=['global', 'regional'],
            path=self.gui_settings.emissions.path,
            remove_callback=self._remove_emission_category,
            visible=visible,
        )

    @param.depends('duplicate_scenario_event', watch=True)
    def _show_duplicate_prompt(self):
        self.widgets['duplicate_prompt'].visible = True

    @param.depends('confirm_duplicate_event', watch=True)
    def _duplicate_scenario(self):
        self.alert = ''

        # Get the value from the "new_scenario_name" widget, once the user has clicked to confirm.
        new_key = self.new_scenario_name.strip()
        if self.select_scenario is None or new_key == '':
            return

        # Prevent editing an existing scenario
        if new_key in self.emission_scenarios:
            self.alert = f"a scenario named '{new_key}' already exists, pick another name"
            return

        # Copy the "source scenario"
        source = self.emission_scenarios[self.select_scenario]

        # Create the new one based on it
        self.emission_scenarios[new_key] = {
            'title': new_key,
            'description': f"Duplicated from \"{source.get('title', self.select_scenario)}\"",
            'categories': source.get('categories', {}),
            'editable': True,
        }

        # Update the "select_scenario" object (re-create it completely in fact)
        self.param.select_scenario.objects = list(self.param.select_scenario.objects) + [new_key]

        # Fix the drop-down menu:
        self._fix_scenario_display_name()

        # Reset "edit" widgets
        self.new_scenario_name = ''
        self.widgets['duplicate_prompt'].visible = False

        #
        self.widgets['show_categories'].visible = False
        # Set the current scenario to the scenario we just created
        self.select_scenario = new_key

    @param.depends('add_category_event', watch=True)
    def _add_emission_category(self, catname: str = None):
        if catname is None:
            catname = f'category_{len(self.emission_scenario) + 1}'
        es = self._build_emission_category(catname)
        self.emission_scenario.append(es)
        self.emission_scenario_widgets.append(es.__panel__())

    @param.depends('hide_categories_event', watch=True)
    def _hide_emission_categories(self):
        self._load_scenario(visible=False)
        self.widgets['show_categories'].visible = True
        self.widgets['hide_categories'].visible = False

    @param.depends('show_categories_event', watch=True)
    def _show_emission_categories(self):
        self._load_scenario(visible=True)
        self.widgets['show_categories'].visible = False
        self.widgets['hide_categories'].visible = True

    def _remove_emission_category(self, es: EmissionSettings):
        self.emission_scenario.remove(es)
        self.emission_scenario_widgets.objects = [e.__panel__() for e in self.emission_scenario]

    @param.depends('select_scenario', watch=True)
    def _load_scenario(self, visible : bool|None = None):
        if self.select_scenario is None:
            return
        scenario = self.emission_scenarios[self.select_scenario]
        editable = scenario.get('editable', False)
        if visible==None: 
            visible = editable
        self.widgets['add_category'].visible = editable
        self.widgets['hide_categories'].visible = visible
        self.widgets['show_categories'].visible = not visible
        emission_scenario = []
        for catname, spec in scenario.get('categories', {}).items():
            es = self._build_emission_category(catname, visible=visible)
            es.set_category(spec)
            emission_scenario.append(es)
        self.emission_scenario = emission_scenario
        self.emission_scenario_widgets.objects = [e.__panel__() for e in self.emission_scenario]

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
        dfc = self.conc.to_dataframe()
        dfc = dfc[dfc.station == self.current_site]
        # msg = f"...@{self.simul_type},cur_exp={cur_exp}: calling plot_conc_timeseries..."
        # logger.debug(msg)
        return plot_conc_timeseries(dfc, self.simul_type, self.select_scenario)

    @param.depends('stats4conc')
    def conc_stats_table(self):
        if self.stats4conc is None:
            return ''
        return plot_stats_table(self.stats4conc, self.select_scenario)

    @param.depends('emissions')
    def map_emissions(self):
        if self.emissions is None:
            # print("resetting emission map")
            logger.debug("resetting emission map")
            return 
        return plot_emission_map(self.emissions, self.select_scenario)

    @param.depends('tgt_table')
    def target_table(self):
        if self.tgt_table is None:
            logger.debug("returning None")
            return ''
        return plot_target_table(self.tgt_table, self.select_scenario)

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

        msg = f"self.emission_scenario ***{self.emission_scenario}***"
        logger.debug(msg)
        
        settings = {
            # 'emis': self.emis_dataset,
            'emissions': {
                'name': self.select_scenario,
                'start': self.gui_settings.start,
                'end': self.gui_settings.end,
                'regions': OmegaConf.to_container(self.gui_settings.regions),
                'categories': {
                    es.catname: {
                        'global': {
                            'file': f'{es.emis_glo.path}/{es.emis_glo.filename}' if f'{es.emis_glo.path}/{es.emis_glo.filename}'.endswith('.nc') else f'{es.emis_glo.path}/{es.emis_glo.filename}*.nc',
                            'field': es.emis_glo.fieldname},
                        **({'regional': {
                            'file': f'{es.emis_reg.path}/{es.emis_reg.filename}' if f'{es.emis_reg.path}/{es.emis_reg.filename}'.endswith('.nc') else f'{es.emis_reg.path}/{es.emis_reg.filename}*.nc',
                            'field': es.emis_reg.fieldname}} if es.switch_reg else {}
                        )
                    }
                    for es in self.emission_scenario
                },
            },
            'task': task,
            'namelist': {
                'fix': (self.correlation_switch=='fixed patterns')
            }
        }

        msg = f"@self.emission_scenario ***{self.emission_scenario}*** yields checksum " \
            f"-->{dict_checksum(settings['emissions'])}<--"
        logger.debug(msg)

        yaml_conf =  OmegaConf.to_yaml(settings)
        try:
            r = requests.post(url, data={'conf':yaml_conf})
        except requests.exceptions.ConnectionError:
            self.alert = f"{task} run failed: could not connect to server -->{url}<--"
            return
        if not r.ok:
            self.alert = f"{task} run failed: backend returned {r.status_code} " \
                f"for configuration *****{yaml_conf}***** at {url}. Body: {r.text[:500]}"
            return
        try:
            payload = r.json()
        except requests.exceptions.JSONDecodeError:
            self.alert = f"{task} run failed: backend returned non-JSON for " \
                f"configuration *****{yaml_conf}**** at {url}. Body: {r.text[:500]}"
            return
        self.alert = ''

        output_path = Path(payload['output'])
        # msg = f"@task={task} for {self.emis_dataset} yields output_path ***{str(output_path)}***"
        # logger.debug(msg)
        if task == 'inversion':
            self.cache_inv[self.select_scenario] = output_path
        else:
            self.cache_fwd[self.select_scenario] = output_path
        return output_path
        
    def _read_concentrations(self, path: Path, task: str):
        if task == 'inversion':
            conc = load_inversion_concentrations(path, self.select_scenario)
        else:
            conc = load_forward_concentrations(path, self.select_scenario)
        # msg = f"...@task={task}, loading concentrations done."
        # logger.debug(msg)
        if self.conc is None:
            conc_update = conc
        else:
            conc_update = xr.merge([self.conc, conc], compat='override')

        # Now update the "sites_list", if needed:
        sites_available = set(conc_update.station.values.reshape(-1))
        update_sites = (sites_available != set(self.sites_list))
        # msg = f"@task={task}, emissions_label -->{self.prefonf_scenario}<-- " \
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
