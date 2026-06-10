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
from holoviews import opts

from tm5 import debug
from tm5.gui.css import *
from tm5.gui.widgets.widget_utils import experiment_desc, plot_site_info, load_observations_metadata

from itertools import cycle
from bokeh.palettes import Category10 as color_palette


@debug.timer
def simulation_read_targets( simu : pn.viewable.Viewer, output_path : str|Path ) -> None:
    """
    """
    #
    #-- read target identifier
    #
    tgt = xr.open_dataset(output_path / 'ftj.nc')
    simu.targets = tgt.targets.values
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
    tfile = Path(output_path) / 'output' / 't0.dat'
    msg = f"reading from ***{tfile.name}***"
    logger.debug(msg)
    if not tfile.exists():
        msg = f"expected prior target file ***{tfile.name}*** not present"
        raise RuntimeError(msg)
    else:
        dft = read_csv(tfile, sep=r"\s+|\t+", engine='python')
        # dft.to_csv('tgt-apri.csv')
        # msg = f"prior targets -->{dft}<--"
        # logger.debug(msg)
        tgt_dict['prior'] = dft.loc[:,'col_1'].values
    #
    #-- posterior
    #
    tfile = Path(output_path) / 'output' / 't.dat'
    msg = f"reading from ***{tfile.name}***"
    logger.debug(msg)
    if not tfile.exists():
        msg = f"expected prior target file ***{tfile.name}*** not present"
        raise RuntimeError(msg)
    else:
        dft = read_csv(tfile, sep=r"\s+|\t+", engine='python')
        # dft.to_csv('tgt-apos.csv')
        # msg = f"prior targets -->{dft}<--"
        # logger.debug(msg)
        tgt_dict['posterior'] = dft.loc[:,'col_1'].values
    #
    #-- posterior uncertainty
    #
    tfile = Path(output_path) / 'output' / 'ct.dat'
    msg = f"reading from ***{tfile.name}***"
    logger.debug(msg)
    if not tfile.exists():
        msg = f"expected prior target file ***{tfile.name}*** not present"
        raise RuntimeError(msg)
    else:
        dft = read_csv(tfile, sep=r"\s+|\t+", engine='python')
        # dft.to_csv('tgtunc-apos.csv')
        for itgt,tgt in enumerate(simu.targets):
            #-- correlation matrix (!), need to take square root of diagonal
            #-- MVO-ATTENTION:column indexing is Fortran based (col_1,col_2,col_3,...)
            _tgtunc = np.sqrt(dft.loc[itgt, f'col_{itgt+1}'])
            tgt_dict['posterior_uncertainty'].append(_tgtunc)
        # msg = f"prior targets -->{dft}<--"
        # logger.debug(msg)
    #
    #-- store target table
    #
    simu.tgt_table = DataFrame.from_dict(tgt_dict).set_index('target')


@debug.timer
def conc_statistics( conc : DataFrame, stations : list[str] ) -> DataFrame:
    """
    """
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
    dfc = conc#.to_dataframe().reset_index()
    for sta in stations:
        # msg = f"now @station={sta}"
        # logger.debug(msg)
        #
        #-- select current station (for all days)
        #
        cnd = dfc['station']==sta
        _df = dfc.loc[cnd,:]
        stats['station'].append(sta)
        _capri = _df.loc[:,'apri']
        _capos = _df.loc[:,'apos']
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


class PreconfExperimentGUI(pn.viewable.Viewer):
    experiment = param.FileSelector(doc='Prior emission dataset')
    run_forward = param.Event(doc='Do a forward run', label='Perform a forward simulation')
    run_inv = param.Event(doc='Do an inversion', label='Perform an inversion')

    # Data containers:
    conc        = param.ClassSelector(class_=DataFrame, precedence=-1)
    stats4conc  = param.ClassSelector(class_=DataFrame, precedence=-1)
    tgt_table   = param.ClassSelector(class_=DataFrame, precedence=-1)

    def __init__(self, gui_settings: DictConfig):
        super().__init__()
        self.button_fwd = pn.widgets.Button.from_param(self.param.run_forward)
        self.button_inv = pn.widgets.Button.from_param(self.param.run_inv)
        self.gui_settings = gui_settings

        # Load the file list
        self.param.experiment.path = self.gui_settings.emissions.glob_pattern
        self.experiment = self.param.experiment.objects[0]
        self.cache_fwd = OrderedDict()
        self.cache_inv = OrderedDict()
        self.stations = None
        # self.stations_widgets = pn.Column()
        self.obs_table = None

    def __panel__(self):
        header_pane = pn.pane.Markdown('# Preconfigured experiments',
                                        stylesheets=[preconfsim_stylesheet],
                                        css_classes=['precomp-right'])
        stats_pane = pn.pane.Markdown('Fit statistics for all stations',
                                      stylesheets=[preconfsim_stylesheet],
                                      css_classes=['precomp-right'])
        header_pane = pn.pane.Markdown('# Preconfigured experiments')
        stats_pane = pn.pane.Markdown('Fit statistics for all stations')
        expdesc_pane = pn.pane.Markdown(
            f"{self._emistable_html()}",
            stylesheets=[preconfsim_stylesheet],
            css_classes=['precomp-right'])
        return pn.Column(
            header_pane,
            pn.Row(pn.widgets.Select.from_param(self.param.experiment),expdesc_pane), 
            # pn.pane.Markdown("== some description of the selected experiment =="),
            pn.Row(self.button_fwd, self.button_inv),
            # self.stations_widgets,
            self._conc_plot,
            self._conc_stats_table,
            self._target_table
            # pn.Column(
            #     stats_pane,
            #     self._conc_stats_table)
        )

    def _emistable_html(self):
        def _get_desc(exp):
            desc = experiment_desc(exp)
            if desc is None:
                return desc
            desc_html = '<tr>'
            desc_html += f'<td>{exp}</td>'
            desc_html += f'<td>{desc}</td>'
            desc_html += f'</tr>'
            return desc_html
        
        exp_list = list(self.param.experiment.objects)
        tbl = f'<table>'
        #-- header
        tbl += f'<thead><tr><th>Emissions setup</th><th>Description</th></tr></thead>'
        for iexp,exp in enumerate(exp_list):
            p = Path(exp).stem
            ptokens = p.split('_')
            _exp = ptokens[0].replace('fitic-','')
            if ptokens[0]=='fitic-regional':
                if ptokens[1]=='monthly-emissions':
                    _exp = 'regional'
                else:
                    _exp = _exp + '_' + ptokens[1]
            desc = _get_desc(_exp)
            if not desc is None:
                tbl += desc
                # if iexp<len(exp_list)-1:
                #     tbl += '<br>'
        tbl += f'</table>'

        return tbl

    def _set_obstable(self):
        #-- code below to visualise station info next to time-series
        #   does *NOT* work yet.
        #   Thus need to keep value None for self.obstable !
        return
        ####
        logger.debug(f"--> start creation of self.obs_table...")
        obsfile_list_all = glob(self.gui_settings.observations.files)
        obsinfo_list = []
        for ista,_staid in enumerate(self.stations):
            staid,sta_alt = _staid.split('_')
            for o in obsfile_list_all:
                p = Path(o)
                # logger.debug(f"@{_staid}: staid-->{staid}<-- {p.name}")
                if p.name.startswith(f'ch4_{staid}'):
                    # msg = f"station -->{_staid}<-- with obsfile ***{o}***"
                    # logger.info(msg)
                    df = load_observations_metadata(p)
                    obsinfo_list.append(df)
        self.obs_table = concat(obsinfo_list)
        logger.debug(f"self.obs_table -->{self.obs_table.shape}<--")
        
    
    @param.depends('run_forward', watch=True)
    def _run_forward(self):
        #
        #-- clear previous results
        #
        self.conc = None
        self.stats4conc = None
        self.tgt_table = None
        # # Here "emis" should point to the file from the "Experiment" selector
        # r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis':self.experiment, 'task':'forward'})

        # # Retrieve results (here just the concentrations):
        # output_path = Path(r.json()['output'])
        #-- modification provided by Zois (2026-05-25)
        #>>MVO:: it must be the full path (otherwise symlink will fail on the backend)
        # emis = Path(self.experiment).name
        emis = str(self.experiment)
        if emis in self.cache_fwd:
            output_path = self.cache_fwd[emis]
            msg = f"forward simulation already cached for emis -->{emis}<--"
            logger.debug(msg)
        else:
            msg = f"running forward simulation for emissions -->{emis}<--"
            logger.debug(msg)
            url = f"{self.gui_settings.backend_url}/forward"
            r = requests.get(url, params={'emis': emis, 'task': 'forward'})

            if not r.ok:
                logger.error(
                    f"forward run failed: backend returned {r.status_code} for "
                    f"emis={emis!r} at {url}. Body: {r.text[:500]}"
                )
                return
            try:
                payload = r.json()
            except requests.exceptions.JSONDecodeError:
                logger.error(
                    f"forward run failed: backend returned non-JSON for "
                    f"emis={emis!r} at {url}. Body: {r.text[:500]}"
                )
                return
            output_path = Path(payload['output'])
            #-- store in cache
            self.cache_fwd[emis] = output_path
        #--
        logger.debug(f"reading from ouput directory {str(output_path)}")

        #
        #-- processing result/output folder
        #
        #-- fc.nc: simulated concentrations using the ingoing emissions
        #          c = iniconc + ojac*emis
        #
        fc = xr.open_dataset(output_path / 'fc.nc')
        #-- store list of stations
        self.stations = fc.station.values #-- get station identifiers
        self.conc_units = fc.obs.units
        conc = fc[['obs', 'conc', 'station','obstime',]]
        conc_df = conc.to_dataframe()
        conc_df['time'] = [Timestamp(_) for _ in conc_df.loc[:,'obstime']]
        #-- rename (NOTE: must be consistent with _conc_plot!)
        conc_df = conc_df.rename(columns={'conc':'forward'})
        # tag = self.experiment
        # self.conc[f'forward_{tag}'] = obs['forward']
        # self.conc['forward'] = obs['forward']
        self.conc = conc_df
        msg = f"self.conc set with columns -->{self.conc.columns}<--"
        logger.debug(msg)
        #
        #--
        #
        logger.debug(f"self.obs_table ==>{self.obs_table}<== ({self.obs_table is None})")
        if self.obs_table is None:
            logger.debug(f"--> setting self.obs_table...")
            self._set_obstable()
        
        
    @param.depends('run_inv', watch=True)
    def _run_inv(self):
        #
        #-- clear previous results
        #
        self.conc = None
        self.stats4conc = None
        self.tgt_table = None
        # r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis':self.experiment, 'task':'inversion'})

        # # Retrieve results (here just the concentrations):
        # # This file seems to not have the info on which station the files come from ...
        # # The is the only thing needed to make the plots site-specific.
        # output_path = Path(r.json()['output'])
        #
        #-- as modified by Zois (2026-05-25)
        #
        #>>MVO:: it must be the full path (otherwise symlink will fail on the backend)
        # emis = Path(self.experiment).name
        
        emis = str(self.experiment)
        if emis in self.cache_inv:
            output_path = self.cache_inv[emis]
            msg = f"forward simulation already cached for emis -->{emis}<--"
            logger.debug(msg)
        else:
            msg = f"running inversion for emissions -->{emis}<--"
            logger.debug(msg)
            url = f"{self.gui_settings.backend_url}/forward"
            r = requests.get(url, params={'emis': emis, 'task': 'inversion'})

            if not r.ok:
                logger.error(
                    f"inversion failed: backend returned {r.status_code} for "
                    f"emis={emis!r} at {url}. Body: {r.text[:500]}"
                )
                return
            try:
                payload = r.json()
            except requests.exceptions.JSONDecodeError:
                logger.error(
                    f"inversion failed: backend returned non-JSON for "
                    f"emis={emis!r} at {url}. Body: {r.text[:500]}"
                )
                return
            output_path = Path(payload['output'])
            #-- store in cache
            self.cache_inv[emis] = output_path
        #--
        logger.debug(f"reading from ouput directory {str(output_path)}")

        #
        #-- processing result/output folder
        #
        #-- 20260526: txk had changed code such that prior/posterior
        #             simulated concentrations (including the signal from
        #             the initial concentration) both are in file
        #             fcpost.nc
        #-- 20260609: mvo start adaption to 1d obs vector output
        fc = xr.open_dataset(output_path / 'fcpost.nc')
        self.stations = fc.station_id.values #-- get station identifiers
        conc = fc[['obs', 'cprior', 'cpost','station','obstime',]]
        #-- turn into data frame
        conc_df = conc.to_dataframe()
        conc_df['time'] = [Timestamp(_) for _ in conc_df.loc[:,'obstime']]
        #-- renaming columns to be consistent with _conc_plot
        conc_df = conc_df.rename(columns={'cprior':'apri','cpost':'apos'})
        
        self.conc = conc_df
        self.conc_units = fc.obs.units
        self.stats4conc = conc_statistics(self.conc, self.stations)
        # msg = f"...setting done."
        # logger.debug(msg)
        #
        #-- read target identifier
        #
        # msg = f"...start reading targets..."
        # logger.info(msg)
        simulation_read_targets(self, output_path)
        # msg = f"...detected targets -->{self.targets}<--"
        # logger.debug(msg)
        # self.tgt_table.to_csv('yy.csv', index=True)
        # msg = f"...generated target table -->\n{self.tgt_table}\n<--"
        # logger.info(msg)
        #
        #--
        #
        if self.obs_table is None:
            self._set_obstable()


    @param.depends('conc')
    def _conc_plot(self):
        # msg = f"self.conc -->{self.conc}<--"
        # logger.debug(msg)
        if self.conc is None:
            return ''
        #--
        dfc = self.conc
        p = dfc.hvplot.points(x='time', y='obs', grid=True, c='k', label='obs', groupby='station')
        if 'forward' in dfc.columns:  #-- only forward simulation
            p *= dfc.hvplot.line(x='time', y='forward', c='r', label='forward', groupby='station')
            # Get all the other (older) forwards:
            # for fwd in [_ for _ in dfc.columns if _.startswith('forward') and not _.endswith(self.experiment)]:
            #     col = next(color_palette)
            #     p *= dfc.hvplot.line(x='nobsday', y=fwd, label=fwd.rsplit('_', maxsplit=1)[0], groupby='station', line_width=1, color=col)
        elif 'apos' in dfc.columns:   #-- result from inversion
            p *= dfc.hvplot.line(x='time', y='apri', label='prior', groupby='station', c='r')
            p *= dfc.hvplot.line(x='time', y='apos', label='posterior', groupby='station', c='c')
            # for apri in [_ for _ in dfc.columns if _.startswith('apri') and not _.endswith(self.experiment)]:
            #     apos = apri.replace('apri', 'apos')
            #     col = next(color_palette)
            #     p *= dfc.hvplot.line(x='time', y='apri', label=apri.rsplit('_', maxsplit=1)[0], groupby='station', color=col, line_dash='time')
            #     p *= dfc.hvplot.line(x='nobsday', y='apos', label=apos.rsplit('_', maxsplit=1)[0], groupby='station', color=col)
        title = "Time-series of observed and simulated concentration at selected station"
        plotcfg = opts.Overlay(title=title, ylabel="[{self.conc_units}]")
        p.opts(plotcfg)
        msg = f"self.obs_table -->{self.obs_table}<--"
        logger.debug(msg)
        return p
        # if self.obs_table is None:
        #     return p
        # else:
        #     return pn.Row(p, plot_site_info(self.obs_table, 'cbw_207'))
        # return p
        # return pn.Column(
        #     pn.pane.Markdown('# Concentration time-series at selected station'),
        #     p
        #     )
        # return p

    @param.depends('stats4conc')
    def _conc_stats_table(self):
        # msg = f"self.stats4conc -->{self.stats4conc}<--"
        # logger.debug(msg)
        if self.stats4conc is None:
            return ''
        else:
            df = self.stats4conc
            nc = len(df.columns)
            formatters = [lambda x: f'{x:.2f}'] * nc
            p = pn.pane.DataFrame(df, text_align='center', formatters=formatters)
            # title = "Fit statistics for all stations"
            # plotcfg = opts.Overlay(title=title)
            # p.opts(plotcfg)
            # return p
            return pn.Column(pn.pane.Markdown('# Fit statistics for all stations'), p)
            # return pn.pane.DataFrame(df, text_align='center', formatters=formatters)

    @param.depends('tgt_table')
    def _target_table(self):
        # msg = f"self.tgt_table -->{self.tgt_table}<--"
        # logger.debug(msg)
        if self.tgt_table is None:
            return ''
        else:
            df = self.tgt_table
            nc = len(df.columns)
            formatters = [lambda x: f'{x:.3f}'] * nc
            p = pn.pane.DataFrame(df, text_align='center', formatters=formatters)
            #-- MVO-TODO::units [MtCH4] should not be hard-coded here
            return pn.Column(pn.pane.Markdown('# Target emission quantities [MtCH4]'), p)
            # return pn.pane.DataFrame(df, text_align='center', formatters=formatters)
