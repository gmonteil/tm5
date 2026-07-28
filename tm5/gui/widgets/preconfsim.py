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
from typing import Tuple, Dict
import io
from numpy import zeros

from tm5 import debug
from tm5.gui.css import *
from tm5.gui.widgets.stations import calc_statistics
from tm5.gui.widgets.widget_utils import experiment_desc, plot_site_info, load_observations_metadata

from itertools import cycle
from bokeh.palettes import Category10
import itertools
import geoviews.feature as gf
from cartopy import crs



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
    stations = set(dfc.station.values)

    # msg = f"stations -->{stations}<--"
    # logger.debug(stations)
    stats = {
        'station': [],
        'Mean bias (prior)': [],
        'Mean bias (post)' : [],
        'RMSE (prior)': [],
        'RMSE (post)': [],
        'Correlation coefficient (prior)': [],
        'Correlation coefficient (post)': [],
        }
    # msg = f"stats initial -->{stats}<-- (==>{stations}<==)"
    # logger.debug(msg)
    for sta in stations:
        # msg = f"now @station={sta}"
        # logger.debug(msg)
        #
        #-- select current station (for all days)
        #
        cnd = dfc['station']==sta
        _df = dfc.loc[cnd,:]
        stats['station'].append(sta)
        _capri = _df.loc[:,f'apri_{label}']
        _capos = _df.loc[:,f'apos_{label}']
        _cobs  = _df.loc[:,'obs']
        # msg = f"...@{sta} RMSE ...({_df.shape})"
        # logger.debug(msg)
        #
        #-- prior statistics
        #
        _bias_prior  = _capri - _cobs
        meanbias_prior = _bias_prior.mean()
        rmse_prior = (_bias_prior ** 2).mean() ** .5
        corrcoef_prior = np.corrcoef(_capri.values,_cobs.values)[0,1]
        # msg = f"@{sta}, rmse prior ==>{rmse_prior}<=="
        # logger.debug(msg)
        #
        #-- posterior statistics
        #
        _bias_post  = _capos - _cobs
        meanbias_post = _bias_post.mean()
        rmse_post  = (_bias_post **2).mean() ** .5
        corrcoef_post = np.corrcoef(_capos.values,_cobs.values)[0,1]
        # msg = f"@{sta}, rmse post ==>{rmse_post}<=="
        # logger.debug(msg)
        #
        #-- fill into dictionary
        #
        stats['Mean bias (prior)'].append(meanbias_prior)
        stats['RMSE (prior)'].append(rmse_prior)
        stats['Correlation coefficient (prior)'].append(corrcoef_prior)
        stats['Mean bias (post)'].append(meanbias_post)
        stats['RMSE (post)'].append(rmse_post)
        stats['Correlation coefficient (post)'].append(corrcoef_post)
    # msg = f"...loop terminated, stats -->{stats}<--"
    # logger.info(msg)
    #
    #-- turn into dataframe
    #
    stats = DataFrame.from_dict(stats).set_index('station')
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
    for station_id in fc.station_id.values:
        stat = fc.sel(nsta = fc.station_id == station_id)
        conc.loc[conc.station == station_id, 'station_lon'] = float(stat.station_lon.values[0])
        conc.loc[conc.station == station_id, 'station_lat'] = float(stat.station_lat.values[0])
    conc['time'] = [Timestamp(_) for _ in conc.loc[:,'obstime']]
    conc = conc.rename(columns={'conc':f'forward_{label}'})
    conc = conc.to_xarray()
    conc.attrs['units'] = fc.obs.units
    return conc


def extract_emis_map(ds, region: str) -> xr.DataArray:
    ds = ds.sel(ng = ds.region == region)
    #-- MVO-20260728::
    #   fe.nc:      data variable is 'emissions'
    #   fepost.nc:  data variables are 'emission_prior' *and* 'emission_post'
    if 'emission' in ds:
        emis_in = ds['emission']
    elif 'emission_post' in ds:
        emis_in = ds['emission_post']
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
    apri = xr.open_dataset(path / 'fe.nc').isel(nmon=0)
    apos = xr.open_dataset(path / 'fepost.nc')
    apos['area'] = apri['area']  # "area" seems missing in the posterior file ... copy it from the prior
    emis = {}
    for region in ['glb600x400', 'eur300x200', 'gns100x100']:
        emis[region] = xr.Dataset(dict(
            apri = extract_emis_map(apri, region), 
            apos = extract_emis_map(apos, region)
        ))
    return emis


class PreconfExperimentGUI(pn.viewable.Viewer):
    experiment = param.FileSelector(doc='Prior emission dataset')
    run_forward = param.Event(doc='Do a forward run', label='Perform a forward simulation')
    run_inv = param.Event(doc='Do an inversion', label='Perform an inversion')
    alert = param.String(doc='Generic object for error messages or others ...', default='')
    current_site = param.Selector(doc='Current site to be displayed', default=None)
    sites_list = param.List(default=[], doc='List of observation sites available (for internal use ...)')
    simul_type = param.Selector(objects=['fwd', 'inv'], allow_None=True, default=None)

    # Data containers:
    conc        = param.ClassSelector(class_=xr.Dataset)
    stats4conc  = param.DataFrame()
    tgt_table   = param.DataFrame()
    emissions   = param.Dict()

    def __init__(self, gui_settings: DictConfig):

        self._message = ''

        super().__init__()
        self.gui_settings = gui_settings

        # Load the file list
        self.param.experiment.path = self.gui_settings.emissions.glob_pattern
        self.experiment = self.param.experiment.objects[0]
        self.cache_fwd = OrderedDict()
        self.cache_inv = OrderedDict()

        # Globally accessible widgets
        self.widgets = {
            'station_selector': pn.widgets.Select.from_param(self.param.current_site)
        }
        self.widgets['station_selector'].visible = False
        self.widgets['borders'] = gf.borders()

        # TEMPORARY:
        self.namelist = {}

    def __panel__(self):
        header_pane = pn.pane.Markdown('# Preconfigured experiments')
        expdesc_pane = pn.pane.Markdown(
                    self._emistable_md(),
                    stylesheets=[preconfsim_stylesheet],
                    css_classes=['precomp-right'])
        logger.debug(f"==>{self.gui_settings.keys()}<==")
        
        show_emismap = self.gui_settings.get('show_emismap', False)
        # logger.debug(f"show_emismap = {show_emismap}")
        if show_emismap:
            return pn.Column(
                header_pane,
                pn.Row(pn.widgets.Select.from_param(self.param.experiment), expdesc_pane),
                pn.Row(
                    pn.widgets.Button.from_param(self.param.run_forward),
                    pn.widgets.Button.from_param(self.param.run_inv)
                ),
                self._alert,
                # self.conc_plot,
                self.widgets['station_selector'],
                pn.Row(self.conc_plot, self.map_sites),
                pn.Row(self.conc_stats_table, self.target_table),
                self.map_emissions
            )
        else:
             return pn.Column(
                header_pane,
                pn.Row(pn.widgets.Select.from_param(self.param.experiment), expdesc_pane),
                pn.Row(
                    pn.widgets.Button.from_param(self.param.run_forward),
                    pn.widgets.Button.from_param(self.param.run_inv)
                ),
                self._alert,
                # self.conc_plot,
                self.widgets['station_selector'],
                pn.Row(self.conc_plot, self.map_sites),
                pn.Row(self.conc_stats_table, self.target_table)
            )

    def _emistable_md(self):
        lines = ['| **Emissions setup** | **Description** |']
        lines.append('| --- | --- |')
        for exp in self.param.experiment.objects:
            exp = get_exp_label(exp)
            if desc := experiment_desc(exp):
                lines.append(f'| {exp} | {desc} |')
        return '\n'.join(lines)

    def _get_output_path(self, task: str) -> Path | None :
        """
        Triggers an inversion on the VM, and return either the output path, or an error code:
            - 100: something went wrong ...
            - 101: result is not valid json
        """
        # if self.experiment in self.cache_inv and task == 'inversion':
        #     return self.cache_inv[self.experiment]
        # elif self.experiment in self.cache_fwd and task == 'forward':
        #     return self.cache_fwd[self.experiment]

        url = f"{self.gui_settings.backend_url}/forward"

        settings = {
            'emis': self.experiment,
            'task': task,
            'namelist': self.namelist
        }

        r = requests.post(url, data={'conf': OmegaConf.to_yaml(settings)})
        if not r.ok:
            self.alert = f"{task} run failed: backend returned {r.status_code} for emis={self.experiment} at {url}. Body: {r.text[:500]}"
            return
        try:
            payload = r.json()
        except requests.exceptions.JSONDecodeError:
            self.alert = f"{task} run failed: backend returned non-JSON for emis={self.experiment} at {url}. Body: {r.text[:500]}"
            return
        self.alert = ''

        output_path = Path(payload['output'])
        # msg = f"@task={task} for {self.experiment} yields output_path ***{str(output_path)}***"
        # logger.debug(msg)
        if task == 'inversion':
            self.cache_inv[self.experiment] = output_path
        else:
            self.cache_fwd[self.experiment] = output_path
        return output_path

    @param.depends('run_forward', watch=True)
    def _run_forward(self):
        output_path = self._get_output_path('forward')
        self.simul_type = 'fwd'
        if output_path is not None:
            self._read_concentrations(output_path, 'forward')
            self.stats4conc = None
            self.tgt_table = None
            self.emissions = None

    @param.depends('run_inv', watch=True)
    def _run_inv(self):
        output_path = self._get_output_path('inversion')
        # msg = f"...output_path ***{output_path}*** ('output_path is not None': {output_path is not None})"
        # logger.debug(msg)
        self.simul_type = 'inv'
        if output_path is not None:
            # msg = f"...start reading concentrations"
            # logger.debug(msg)
            self._read_concentrations(output_path, 'inversion')
            # msg = f"...start reading emissions"
            # logger.debug(msg)
            self.emissions = load_emissions(output_path)
            # msg = f"...computing conc statistics"
            # logger.debug(msg)
            self.stats4conc = conc_statistics(self.conc, get_exp_label(self.experiment))
            # msg = f"...reading simulation targets from directory ***{str(output_path)}***"
            # logger.debug(msg)
            self.tgt_table = simulation_read_targets(output_path)

    def _read_concentrations(self, path: Path, task: str):
        label = get_exp_label(self.experiment)
        if task == 'inversion':
            conc = load_inversion_concentrations(path, label)
        else:
            conc = load_forward_concentrations(path, label)

        if self.conc is None:
            self.conc = conc
        else:
            self.conc = xr.merge([self.conc, conc], compat='override')

        # Now update the "sites_list", if needed:
        sites_available = set(self.conc.station.values.reshape(-1))
        if sites_available != set(self.sites_list):
            self.sites_list = list(sites_available)
            self.param.current_site.objects = set(self.conc.station.values.reshape(-1))
            self.current_site = self.param.current_site.objects[0]
            self.widgets['station_selector'].visible = True

    @param.depends('alert')
    def _alert(self):
        if self.alert == '':
            return ''
        return pn.pane.Alert(self.alert, alert_type='danger')

    @param.depends('conc', 'current_site')
    def conc_plot(self):
        if self.conc is None:
            return ''
        if self.current_site is None:
            return ''
        cur_exp = get_exp_label(self.experiment)
        dfc = self.conc.to_dataframe()
        dfc = dfc[dfc.station == self.current_site]
        p = dfc.hvplot.points(x='time', y='obs', grid=True, c='k', label='obs', width=1200, height=400)
        color_palette = itertools.cycle(Category10[10])
        # print(cur_exp, dfc.columns)

        # Find all "forward" experiments
        if self.simul_type == 'fwd':
            experiments = [c.split('_', maxsplit=1)[1] for c in dfc.columns if c.startswith('forward_')]
            # print(experiments)
            for iexp, exp in enumerate(experiments):
                col = next(color_palette) # Category10[10][iexp]
                # print(exp, iexp, col)
                if exp == cur_exp:
                    p *= dfc.hvplot.line(x='time', y=f'forward_{exp}', c=col, label=exp, muted_alpha=0, line_width=4)
                else:
                    p *= dfc.hvplot.line(x='time', y=f'forward_{exp}', c=col, label=exp, muted_alpha=0, line_width=2)

        # Find all "inversion" experiments
        elif self.simul_type == 'inv':
            experiments = [c[5:] for c in dfc.columns if c.startswith('apri_')]
            # Plot the experiments:
            for iexp, exp in enumerate(experiments):
                col = next(color_palette) # Category10[10][iexp]
                if exp == cur_exp:
                    p *= dfc.hvplot.line(x='time', y=f'apri_{exp}', c=col, line_dash='dashed', label=f'prior_{exp}', line_width=4, muted_alpha=0)
                    p *= dfc.hvplot.line(x='time', y=f'apos_{exp}', c=col, label=f'posterior_{exp}', line_width=4, muted_alpha=0)
                else:
                    p *= dfc.hvplot.line(x='time', y=f'apri_{exp}', c=col, line_dash='dashed', label=f'prior_{exp}', line_width=2, muted_alpha=0)
                    p *= dfc.hvplot.line(x='time', y=f'apos_{exp}', c=col, label=f'posterior_{exp}', line_width=2, muted_alpha=0)
        return p

    @param.depends('stats4conc')
    def conc_stats_table(self):
        if self.stats4conc is None:
            return ''
        else:
            df = self.stats4conc
            nc = len(df.columns)
            formatters = [lambda x: f'{x:.2f}'] * nc
            p = pn.pane.DataFrame(df, text_align='center', formatters=formatters)
            title = f'# Fit statistics for all stations ({get_exp_label(self.experiment)})'
            return pn.Column(pn.pane.Markdown(title), p)

    @param.depends('tgt_table')
    def target_table(self):

        if self.tgt_table is None:
            logger.debug("returning None")
            return ''
        else:
            df = self.tgt_table
            def get_csv_file():
                # You can re-compute data here if it's dynamic
                fid = io.BytesIO()
                fid.write(df.to_csv(index=False).encode('utf-8'))
                fid.seek(0)
                return fid

            nc = len(df.columns)
            formatters = [lambda x: f'{x:.3f}'] * nc
            # outname = f"targets_{get_exp_label(self.experiment)}.csv"
            # outname = f"targets.csv"
            # button = pn.widgets.FileDownload(callback=get_csv_file, filename=outname, label="Download Data (CSV)", button_type="primary")
            button = pn.widgets.FileDownload(callback=get_csv_file, filename="targets.csv", label="Download Data (CSV)", button_type="primary")

            p = pn.Column(pn.pane.DataFrame(df, text_align='center', formatters=formatters), button)
            # return self.tgt_table
            #-- MVO-TODO::units [MtCH4] should not be hard-coded here
            logger.debug("...returning target table now")
            title = f'# Target emission quantities [MtCH4] ({get_exp_label(self.experiment)})'
            return pn.Column(pn.pane.Markdown(title), p)

    @param.depends('current_site', 'sites_list')
    def map_sites(self):
        if self.conc is None:
            return ''
        if self.current_site is None:
            return ''
        df = self.conc.to_dataframe().loc[:, ['station', 'station_lon', 'station_lat']].drop_duplicates()
        df.loc[:, 'cur_site'] = 0
        df.loc[df.station == self.current_site, 'cur_site'] = 1
        return df.hvplot.points(
            x='station_lon', y='station_lat', color='cur_site', cmap=['LightSlateGray', 'red'],
            geo=True, coastline=True, xlim=(-15, 35), ylim=(33, 73), colorbar=False, tiles='EsriTerrain'
        )
        #-- widgets-boarders currently make problems
        #   on exploredata.icos-cp.eu,
        #   disabled for NCGG10
        return df.hvplot.points(
            x='station_lon', y='station_lat', color='cur_site', cmap=['LightSlateGray', 'red'],
            geo=True, coastline=True, xlim=(-15, 35), ylim=(33, 73), colorbar=False, tiles='EsriTerrain'
        ) * self.widgets['borders']

    @param.depends('emissions')
    def map_emissions(self):
        if self.emissions is None:
            # print("resetting emission map")
            logger.debug("resetting emission map")
            return
        # print("computing emission map")
        logger.debug("computing emission map")
            
        cmap = 'RdBu_r'
        clim = (-0.005, 0.005)

        projection = crs.RotatedPole(pole_longitude=185, pole_latitude=50)
        xlim = (-15, 35)
        ylim = (33, 73)
        
        # projection = crs.PlateCarree()
        # xlim = (-180,180)
        # ylim = (-90,90)
        # projection = crs.PlateCarree()
        # xlim = (-180,180)
        # ylim = (-90,90)
        proj_glb600x400 = crs.PlateCarree()
        xlim_glb600x400 = (-180,180)
        ylim_glb600x400 = (-90,90)
        # 36W-54E, 22N-74N
        proj_eur300x200 = crs.PlateCarree()
        xlim_eur300x200 = (-36,54)
        ylim_eur300x200 = (22,74)
        # 0E-18E, 42N-58N
        proj_gns100x100 = crs.RotatedPole(pole_longitude=185, pole_latitude=50)
        proj_gns100x100 = crs.PlateCarree()
        xlim_gns100x100 = (0,18)
        ylim_gns100x100 = (42,58)
        
        logger.debug("...returning emissions map now!!!")
        # apri_glb600x400 = self.emissions['glb600x400'].apri
        # apos_glb600x400 = self.emissions['glb600x400'].apos
        # apri_eur300x200 = self.emissions['eur300x200'].apri
        # apos_eur300x200 = self.emissions['eur300x200'].apos
        # apri_gns100x100 = self.emissions['gns100x100'].apri
        # apos_gns100x100 = self.emissions['gns100x100'].apos

        # p = pn.Row(
        #     apri_glb600x400.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_600x400, xlim=xlim_600x400, ylim=ylim_600x400),
        #     apos_glb600x400.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_600x400, xlim=xlim_600x400, ylim=ylim_600x400),
        #     (apos_glb600x400-apri_glb600x400).quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_600x400, xlim=xlim_600x400, ylim=ylim_600x400)     
        # )
        # p = pn.Row(
        #     apri_glb600x400.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_600x400, xlim=xlim_600x400, ylim=ylim_600x400),
        #     apos_glb600x400.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_600x400, xlim=xlim_600x400, ylim=ylim_600x400)
        # )
                   
        
        # p = (self.emissions['glb600x400'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) +
        #     self.emissions['glb600x400'].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim)
        #     )
        emis_units = self.emissions['glb600x400'].apos.attrs['units']
        logger.debug(f"emis_units ==>{emis_units}<==")
        mode = 'guillaume'
        if mode=='glb600x400':
            prow = pn.Row(
                self.emissions['glb600x400'].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400,),
                self.emissions['glb600x400'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400),
                (self.emissions['glb600x400'].apos-self.emissions['glb600x400'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400,)
            )
            title = f"# emssions maps (posterior, prior, posterior-prior) ({get_exp_label(self.experiment)})"
            p = pn.Column(pn.pane.Markdown(title), prow)
        elif mode=='eur300x200':
            prow_eur300x200 = pn.Row(
                self.emissions['eur300x200'].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_eur300x200, xlim=xlim_eur300x200, ylim=ylim_eur300x200,),
                self.emissions['eur300x200'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_eur300x200, xlim=xlim_eur300x200, ylim=ylim_eur300x200),
                (self.emissions['eur300x200'].apos-self.emissions['eur300x200'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_eur300x200, xlim=xlim_eur300x200, ylim=ylim_eur300x200,)
            )
            title = f"# emssions maps (posterior, prior, posterior-prior) ({get_exp_label(self.experiment)})"
            p = pn.Column(pn.pane.Markdown(title), prow)
        elif mode=='gns100x100':
            prow = pn.Row(
                self.emissions['gns100x100'].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_gns100x100, xlim=xlim_gns100x100, ylim=ylim_gns100x100,),
                self.emissions['gns100x100'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_gns100x100, xlim=xlim_gns100x100, ylim=ylim_gns100x100),
                (self.emissions['gns100x100'].apos-self.emissions['gns100x100'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_gns100x100, xlim=xlim_gns100x100, ylim=ylim_gns100x100,)
            )
            title = f"# emssions maps (posterior, prior, posterior-prior) ({get_exp_label(self.experiment)})"
            p = pn.Column(pn.pane.Markdown(title), prow)
        elif mode=='row_merged':
            prow = pn.Row(
                self.emissions['glb600x400'].apos.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400) *
                self.emissions['eur300x200'].apos.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400) *
                self.emissions['gns100x100'].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400,),
                self.emissions['glb600x400'].apri.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400) *
                self.emissions['eur300x200'].apri.hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400) *
                self.emissions['gns100x100'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400,),
                (self.emissions['glb600x400'].apos - self.emissions['glb600x400'].apri).hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400) *
                (self.emissions['eur300x200'].apos-self.emissions['eur300x200'].apri).hvplot.quadmesh(rasterize=True, geo=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400) *
                (self.emissions['gns100x100'].apos-self.emissions['gns100x100'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=proj_glb600x400, xlim=xlim_glb600x400, ylim=ylim_glb600x400,)
            )
            title = f"# emssions maps (posterior, prior, posterior-prior) ({get_exp_label(self.experiment)})"
            p = pn.Column(pn.pane.Markdown(title), prow)
        elif mode=='guillaume':
            ppost = (
                self.emissions['glb600x400'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
                self.emissions['eur300x200'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
                self.emissions['gns100x100'].apos.hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
            )
            plotcfg = opts.Overlay(title=f'posterior emissions ({get_exp_label(self.experiment)})', ylabel=f"[{emis_units}]" )
            ppost.opts(plotcfg)
            pprior = (
                self.emissions['glb600x400'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
                self.emissions['eur300x200'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
                self.emissions['gns100x100'].apri.hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
                )
            plotcfg = opts.Overlay(title=f'prior emissions ({get_exp_label(self.experiment)})' )
            pprior.opts(plotcfg)
            pdiff = (
                (self.emissions['glb600x400'].apos - self.emissions['glb600x400'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
                (self.emissions['eur300x200'].apos - self.emissions['eur300x200'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, xlim=xlim, ylim=ylim, projection=projection) *
                (self.emissions['gns100x100'].apos - self.emissions['gns100x100'].apri).hvplot.quadmesh(cmap=cmap, clim=clim, coastline=True, xlim=xlim, ylim=ylim, projection=projection)
                )
            plotcfg = opts.Overlay(title=f'posterior-prior emissions ({get_exp_label(self.experiment)})')
            pdiff.opts(plotcfg)
            
            prow = (
                ppost + pprior + pdiff
            )
            p = prow
        ##
 
        # logger.debug(f"...returning p ({type(p)})")
        return p
    
        # return (
        #     self.emissions['glb600x400'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) +
        #     self.emissions['glb600x400'].apos.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim)
        #     )
        # return (
        #     self.emissions['glb600x400'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     self.emissions['eur300x200'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     self.emissions['gns100x100'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim)#  *
        #     # self.widgets['borders']
        #     )
        # return (
        #     self.emissions['glb600x400'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     self.emissions['eur300x200'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     self.emissions['gns100x100'].apri.hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     self.widgets['borders'] +
        #     (self.emissions['glb600x400'].apos - self.emissions['glb600x400'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     (self.emissions['eur300x200'].apos - self.emissions['eur300x200'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     (self.emissions['gns100x100'].apos - self.emissions['gns100x100'].apri).hvplot.quadmesh(rasterize=True, geo=True, coastline=True, cmap=cmap, clim=clim, projection=projection, xlim=xlim, ylim=ylim) *
        #     self.widgets['borders']
        # )

# class PreconfExperimentGUI_(pn.viewable.Viewer):
#     experiment = param.FileSelector(doc='Prior emission dataset')
#     run_forward = param.Event(doc='Do a forward run', label='Perform a forward simulation')
#     run_inv = param.Event(doc='Do an inversion', label='Perform an inversion')

#     # Data containers:
#     conc        = param.ClassSelector(class_=DataFrame, precedence=-1)
#     stats4conc  = param.ClassSelector(class_=DataFrame, precedence=-1)
#     tgt_table   = param.ClassSelector(class_=DataFrame, precedence=-1)

#     def __init__(self, gui_settings: DictConfig):
#         super().__init__()
#         self.button_fwd = pn.widgets.Button.from_param(self.param.run_forward)
#         self.button_inv = pn.widgets.Button.from_param(self.param.run_inv)
#         self.gui_settings = gui_settings

#         # Load the file list
#         self.param.experiment.path = self.gui_settings.emissions.glob_pattern
#         self.experiment = self.param.experiment.objects[0]
#         self.cache_fwd = OrderedDict()
#         self.cache_inv = OrderedDict()
#         self.stations = None
#         # self.stations_widgets = pn.Column()
#         self.obs_table = None

#     def __panel__(self):
#         header_pane = pn.pane.Markdown('# Preconfigured experiments',
#                                         stylesheets=[preconfsim_stylesheet],
#                                         css_classes=['precomp-right'])
#         stats_pane = pn.pane.Markdown('Fit statistics for all stations',
#                                       stylesheets=[preconfsim_stylesheet],
#                                       css_classes=['precomp-right'])
#         header_pane = pn.pane.Markdown('# Preconfigured experiments')
#         stats_pane = pn.pane.Markdown('Fit statistics for all stations')
#         expdesc_pane = pn.pane.Markdown(
#             f"{self._emistable_html()}",
#             stylesheets=[preconfsim_stylesheet],
#             css_classes=['precomp-right'])
#         return pn.Column(
#             header_pane,
#             pn.Row(pn.widgets.Select.from_param(self.param.experiment),expdesc_pane),
#             # pn.pane.Markdown("== some description of the selected experiment =="),
#             pn.Row(self.button_fwd, self.button_inv),
#             # self.stations_widgets,
#             self._results_pane
#             # self._conc_plot,
#             # self._conc_stats_table,
#             # self._target_table
#             # pn.Column(
#             #     stats_pane,
#             #     self._conc_stats_table)
#         )

#     def _emistable_html(self):
#         def _get_desc(exp):
#             desc = experiment_desc(exp)
#             if desc is None:
#                 return desc
#             desc_html = '<tr>'
#             desc_html += f'<td>{exp}</td>'
#             desc_html += f'<td>{desc}</td>'
#             desc_html += f'</tr>'
#             return desc_html

#         exp_list = list(self.param.experiment.objects)
#         tbl = f'<table>'
#         #-- header
#         tbl += f'<thead><tr><th>Emissions setup</th><th>Description</th></tr></thead>'
#         for iexp,exp in enumerate(exp_list):
#             p = Path(exp).stem
#             ptokens = p.split('_')
#             _exp = ptokens[0].replace('fitic-','')
#             if ptokens[0]=='fitic-regional':
#                 if ptokens[1]=='monthly-emissions':
#                     _exp = 'regional'
#                 else:
#                     _exp = _exp + '_' + ptokens[1]
#             desc = _get_desc(_exp)
#             if not desc is None:
#                 tbl += desc
#                 # if iexp<len(exp_list)-1:
#                 #     tbl += '<br>'
#         tbl += f'</table>'

#         return tbl

#     # def _set_obstable(self):
#     #     #-- code below to visualise station info next to time-series
#     #     #   does *NOT* work yet.
#     #     #   Thus need to keep value None for self.obstable !
#     #     return
#     #     ####
#     #     logger.debug(f"--> start creation of self.obs_table...")
#     #     obsfile_list_all = glob(self.gui_settings.observations.files)
#     #     obsinfo_list = []
#     #     for ista,_staid in enumerate(self.stations):
#     #         staid,sta_alt = _staid.split('_')
#     #         for o in obsfile_list_all:
#     #             p = Path(o)
#     #             # logger.debug(f"@{_staid}: staid-->{staid}<-- {p.name}")
#     #             if p.name.startswith(f'ch4_{staid}'):
#     #                 # msg = f"station -->{_staid}<-- with obsfile ***{o}***"
#     #                 # logger.info(msg)
#     #                 df = load_observations_metadata(p)
#     #                 obsinfo_list.append(df)
#     #     self.obs_table = concat(obsinfo_list)
#     #     logger.debug(f"self.obs_table -->{self.obs_table.shape}<--")

#     @param.depends('run_forward', watch=True)
#     def _run_forward(self):
#         #
#         #-- clear previous results
#         #
#         self.conc = None
#         self.stats4conc = None
#         self.tgt_table = None
#         # # Here "emis" should point to the file from the "Experiment" selector
#         # r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis':self.experiment, 'task':'forward'})

#         # # Retrieve results (here just the concentrations):
#         # output_path = Path(r.json()['output'])
#         #-- modification provided by Zois (2026-05-25)
#         #>>MVO:: it must be the full path (otherwise symlink will fail on the backend)
#         # emis = Path(self.experiment).name
#         emis = str(self.experiment)
#         if emis in self.cache_fwd:
#             output_path = self.cache_fwd[emis]
#             msg = f"forward simulation already cached for emis -->{emis}<--"
#             logger.debug(msg)
#         else:
#             msg = f"running forward simulation for emissions -->{emis}<--"
#             logger.debug(msg)
#             url = f"{self.gui_settings.backend_url}/forward"
#             r = requests.get(url, params={'emis': emis, 'task': 'forward'})

#             if not r.ok:
#                 logger.error(
#                     f"forward run failed: backend returned {r.status_code} for "
#                     f"emis={emis!r} at {url}. Body: {r.text[:500]}"
#                 )
#                 return
#             try:
#                 payload = r.json()
#             except requests.exceptions.JSONDecodeError:
#                 logger.error(
#                     f"forward run failed: backend returned non-JSON for "
#                     f"emis={emis!r} at {url}. Body: {r.text[:500]}"
#                 )
#                 return
#             output_path = Path(payload['output'])
#             #-- store in cache
#             self.cache_fwd[emis] = output_path
#         #--
#         logger.debug(f"reading from ouput directory {str(output_path)}")

#         #
#         #-- processing result/output folder
#         #
#         #-- fc.nc: simulated concentrations using the ingoing emissions
#         #          c = iniconc + ojac*emis
#         #
#         fc = xr.open_dataset(output_path / 'fc.nc')
#         #-- store list of stations
#         self.stations = fc.station.values #-- get station identifiers
#         self.conc_units = fc.obs.units
#         conc = fc[['obs', 'conc', 'station','obstime',]]
#         conc_df = conc.to_dataframe()
#         conc_df['time'] = [Timestamp(_) for _ in conc_df.loc[:,'obstime']]
#         #-- rename (NOTE: must be consistent with _conc_plot!)
#         conc_df = conc_df.rename(columns={'conc':'forward'})
#         # tag = self.experiment
#         # self.conc[f'forward_{tag}'] = obs['forward']
#         # self.conc['forward'] = obs['forward']
#         self.conc = conc_df
#         msg = f"self.conc set with columns -->{self.conc.columns}<--"
#         logger.debug(msg)
#         #
#         #--
#         #
#         logger.debug(f"self.obs_table ==>{self.obs_table}<== ({self.obs_table is None})")
#         if self.obs_table is None:
#             logger.debug(f"--> setting self.obs_table...")
#             self._set_obstable()

#     def _get_inv_path(self) -> Path | int :
#         """
#         Triggers an inversion on the VM, and return either the output path, or an error code:
#             - 100: something went wrong ...
#             - 101: result is not valid json
#         """
#         if self.experiment in self.cache_inv:
#             return self.cache_inv[self.experiment]

#         r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis': self.experiment, 'task': 'inversion'})
#         if not r.ok:
#             return 100

#         try:
#             payload = r.json()
#         except requests.exceptions.JSONDecodeError:
#             return 101

#         output_path = Path(payload['output'])
#         self.cache_inv[self.experiment] = output_path
#         return output_path

#     def _load_inversion_results(self, path: Path):
#         fc = xr.open_dataset(path / 'fcpost.nc')

#         # Reformat the data as a dataframe, consistent with _conc_plot
#         conc = fc[['obs', 'cprior', 'cpost','station','obstime',]].to_dataframe()
#         conc['time'] = [Timestamp(_) for _ in conc.loc[:,'obstime']]
#         conc = conc.rename(columns={'cprior':'apri','cpost':'apos'})

#         # Save ...
#         self.stations = fc.station_id.values #-- get station identifiers
#         self.conc = conc
#         self.conc_units = fc.obs.units

#     @param.depends('run_inv', watch=True)
#     def _run_inv(self):
#         self.tgt_table = None

#         output_path = self._get_inv_path()
#         if isinstance(output_path, int):
#             return
#         self._load_inversion_results(output_path)
#         self.stats4conc = conc_statistics(self.conc, self.stations)
#         self._read_targets(output_path)
#         self.targets, self.tgt_table = simulation_read_targets(output_path)
#         # if self.obs_table is None:
#         #     self._set_obstable()

#     @param.depends('conc')
#     def _conc_plot(self):
#         # msg = f"self.conc -->{self.conc}<--"
#         # logger.debug(msg)
#         if self.conc is None:
#             return ''
#         #--
#         dfc = self.conc
#         p = dfc.hvplot.points(x='time', y='obs', grid=True, c='k', label='obs', groupby='station')
#         if 'forward' in dfc.columns:  #-- only forward simulation
#             p *= dfc.hvplot.line(x='time', y='forward', c='r', label='forward', groupby='station')
#             # Get all the other (older) forwards:
#             # for fwd in [_ for _ in dfc.columns if _.startswith('forward') and not _.endswith(self.experiment)]:
#             #     col = next(color_palette)
#             #     p *= dfc.hvplot.line(x='nobsday', y=fwd, label=fwd.rsplit('_', maxsplit=1)[0], groupby='station', line_width=1, color=col)
#         elif 'apos' in dfc.columns:   #-- result from inversion
#             p *= dfc.hvplot.line(x='time', y='apri', label='prior', groupby='station', c='r')
#             p *= dfc.hvplot.line(x='time', y='apos', label='posterior', groupby='station', c='c')
#             # for apri in [_ for _ in dfc.columns if _.startswith('apri') and not _.endswith(self.experiment)]:
#             #     apos = apri.replace('apri', 'apos')
#             #     col = next(color_palette)
#             #     p *= dfc.hvplot.line(x='time', y='apri', label=apri.rsplit('_', maxsplit=1)[0], groupby='station', color=col, line_dash='time')
#             #     p *= dfc.hvplot.line(x='nobsday', y='apos', label=apos.rsplit('_', maxsplit=1)[0], groupby='station', color=col)
#         title = "Time-series of observed and simulated concentration at selected station"
#         plotcfg = opts.Overlay(title=title, ylabel="[{self.conc_units}]")
#         p.opts(plotcfg)
#         msg = f"self.obs_table -->{self.obs_table}<--"
#         logger.debug(msg)
#         return p
#         # if self.obs_table is None:
#         #     return p
#         # else:
#         #     return pn.Row(p, plot_site_info(self.obs_table, 'cbw_207'))
#         # return p
#         # return pn.Column(
#         #     pn.pane.Markdown('# Concentration time-series at selected station'),
#         #     p
#         #     )
#         # return p

#     @param.depends('stats4conc')
#     def _conc_stats_table(self):
#         # msg = f"self.stats4conc -->{self.stats4conc}<--"
#         # logger.debug(msg)
#         if self.stats4conc is None:
#             return ''
#         else:
#             df = self.stats4conc
#             nc = len(df.columns)
#             formatters = [lambda x: f'{x:.2f}'] * nc
#             p = pn.pane.DataFrame(df, text_align='center', formatters=formatters)
#             # title = "Fit statistics for all stations"
#             # plotcfg = opts.Overlay(title=title)
#             # p.opts(plotcfg)
#             # return p
#             return pn.Column(pn.pane.Markdown('# Fit statistics for all stations'), p)
#             # return pn.pane.DataFrame(df, text_align='center', formatters=formatters)

#     @param.depends('tgt_table')
#     def _target_table(self):
#         # msg = f"self.tgt_table -->{self.tgt_table}<--"
#         # logger.debug(msg)
#         if self.tgt_table is None:
#             return ''
#         else:
#             df = self.tgt_table
#             nc = len(df.columns)
#             formatters = [lambda x: f'{x:.3f}'] * nc
#             p = pn.pane.DataFrame(df, text_align='center', formatters=formatters)
#             #-- MVO-TODO::units [MtCH4] should not be hard-coded here
#             return pn.Column(pn.pane.Markdown('# Target emission quantities [MtCH4]'), p)
#             # return pn.pane.DataFrame(df, text_align='center', formatters=formatters)
