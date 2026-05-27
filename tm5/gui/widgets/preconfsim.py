#!/usr/bin/env python

import panel as pn
import param
from omegaconf import OmegaConf, DictConfig
import requests
from pathlib import Path
import sys
from loguru import logger
import numpy as np
import xarray as xr
import hvplot.xarray
from pandas import read_csv, DataFrame
from numpy import corrcoef

from tm5 import debug
from tm5.gui.css import *

@debug.timer
def conc_statistics( conc : xr.Dataset, stations : list[str] ) -> DataFrame:
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
    dfc = conc.to_dataframe().reset_index()
    #-- prefer station identifier (rather than station index) for output
    dfc.loc[:,'station'] = stations[dfc.loc[:,'nsta'].values]
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
        corrcoef_prior = corrcoef(_capri.values,_cobs.values)[0,1]
        # msg = f"@{sta}, rmse prior ==>{rmse_prior}<=="
        # logger.debug(msg)
        #
        #-- posterior statistics
        #
        _bias_post  = _capos - _cobs
        meanbias_post = _bias_post.mean()       
        rmse_post  = (_bias_post **2).mean() ** .5
        corrcoef_post = corrcoef(_capos.values,_cobs.values)[0,1]
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
    run_forward = param.Event(doc='Do a forward run', label='Submit a new forward run')
    run_inv = param.Event(doc='Do an inversion', label='Perform an inversion')

    # Data containers:
    conc       = param.ClassSelector(class_=xr.Dataset, precedence=-1)
    stats4conc = param.ClassSelector(class_=DataFrame, precedence=-1)

    def __init__(self, gui_settings: DictConfig):
        super().__init__()
        self.button_fwd = pn.widgets.Button.from_param(self.param.run_forward)
        self.button_inv = pn.widgets.Button.from_param(self.param.run_inv)
        self.gui_settings = gui_settings

        # Load the file list
        self.param.experiment.path = self.gui_settings.emissions.glob_pattern
        self.experiment = self.param.experiment.objects[0]
        self.stations = None
        # self.stations_widgets = pn.Column()

    def __panel__(self):
        return pn.Column(
            pn.pane.Markdown("# Preconfigured experiments"),
            pn.widgets.Select.from_param(self.param.experiment), 
            pn.pane.Markdown("== some description of the selected experiment =="),
            pn.Row(self.button_fwd, self.button_inv),
            # self.stations_widgets,
            self._conc_plot,
            self._conc_stats_table
        )

    @param.depends('run_forward', watch=True)
    def _run_forward(self):
        # # Here "emis" should point to the file from the "Experiment" selector
        # r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis':self.experiment, 'task':'forward'})

        # # Retrieve results (here just the concentrations):
        # output_path = Path(r.json()['output'])
        #-- modification provided by Zois (2026-05-25)
        #>>MVO:: it must be the full path (otherwise symlink will fail on the backend)
        # emis = Path(self.experiment).name
        emis = str(self.experiment)
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
      
        logger.debug(f"reading from ouput directory {str(output_path)}")

        #
        #-- processing result/output folder
        #
        #-- fc.nc: simulated concentrations using the ingoing emissions
        #          c = iniconc + ojac*emis
        fc = xr.open_dataset(output_path / 'fc.nc')
        #-- for now extracting the observations from file foj.nc
        #   NOTE: this will certainly change, but currently all
        #         input for the inversion (including obs.) are present in this file
        obs = xr.open_dataset(output_path / 'foj.nc')
        #-- store list of stations
        self.stations = obs.station.values #-- get station identifiers
        # logger.debug(f"stations -->{self.stations}<--")
        obs['forward'] = fc['conc']
        self.conc = obs[['obs', 'forward', 'iniconc']]

    @param.depends('run_inv', watch=True)
    def _run_inv(self):
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

        #
        #-- processing result/output folder
        #
        #-- 20260526: txk had changed code such that prior/posterior
        #             simulated concentrations (including the signal from
        #             the initial concentration) both are in file
        #             fcpost.nc
        fc = xr.open_dataset(output_path / 'fcpost.nc')
        
        obs = xr.open_dataset(output_path / 'foj.nc')
        self.stations = obs.station.values #-- get station identifiers
        obs['apri'] = fc['cprior']
        obs['apos'] = fc['cpost']
        # msg = f"setting self.conc/self.stats4conc"
        # logger.debug(msg)
        self.conc = obs[['obs', 'apri', 'apos']]
        self.stats4conc = conc_statistics(self.conc, self.stations)
        # msg = f"...setting done."
        # logger.debug(msg)

    @param.depends('conc')
    def _conc_plot(self):
        # msg = f"self.conc -->{self.conc}<--"
        # logger.debug(msg)
        if self.conc is None:
            return ''
        dfc = self.conc.to_dataframe().reset_index()
        # logger.debug(f"dfc.columns -->{dfc.columns}<--")
        #-- prefer station identifier (rather than station index) for output
        dfc.loc[:,'station'] = self.stations[dfc.loc[:,'nsta'].values]
        p = dfc.hvplot.points(x='nobsday', y='obs', grid=True, c='k', label='obs', groupby='station')
        if 'forward' in self.conc:  #-- only forward simulation
            p *= dfc.hvplot.line(x='nobsday', y='forward', c='r', label='forward', groupby='station')
        elif 'apos' in self.conc:   #-- result from inversion
            p *= dfc.hvplot.line(x='nobsday', y='apri', c='r', label='prior', groupby='station')
            p *= dfc.hvplot.line(x='nobsday', y='apos', c='c', label='posterior', groupby='station')
        return p

    @param.depends('stats4conc')
    def _conc_stats_table(self):
        # msg = f"self.stats4conc -->{self.stats4conc}<--"
        # logger.debug(msg)
        if self.stats4conc is None:
            return ''
        else:
            df = self.stats4conc
            return pn.pane.DataFrame(df, text_align='center')
            # return pn.pane.DataFrame(df, formatters=[lambda x: f'{x:.2f}'] * 3, text_align='center')
