#!/usr/bin/env python

import hashlib
from pandas import DataFrame, Timestamp, concat, Grouper, merge
from pathlib import Path
import xarray as xr
from netCDF4 import Dataset
from numpy import zeros, corrcoef, where
import panel as pn
import param
from omegaconf import DictConfig
from glob import glob
import holoviews as hv
from holoviews import opts
import hvplot.pandas
from functools import partial
from typing import List
from tqdm.contrib.concurrent import process_map
import xxhash
import numpy as np
import multiprocessing.pool as mp
from loguru import logger

from tm5 import debug
from tm5.gui.css import *
from tm5.gui.widgets.widget_utils import experiment_desc

CSS2 = """
:host {
    --font-size: 64px;
}
"""

class PrecomputedInfo(pn.viewable.Viewer):

    def __init__(self, settings: DictConfig):
        super().__init__()
        self.settings = settings

    def __panel__(self):
        info1 = """
        ## <ins>Precomputed experiments</ins>:

        A series of TM5 simulations for 2021 have been
        carried out with the hourly simulated Methane concentrations
        recorded at a list of sites where obspack observations
        are available for comparison.

        The standard setup (referred to as 'default' below) applies CAMS
        global atmospheric concentrations derived from inversion of surface observations (v23r1) as initial condition, and the following Methane fluxes are ingested:
        - wetland emissions at daily time-scale (LPJ-GUESS, AVENGERS WP2)
        - emissions from mineral soils at daily time-scale (LPJ-GUESS, AVENGERS WP2)
        - anthropogenic emissions (agriculture, fossil, and waste) at monthly time-scale taken from EDGAR (v8)
        - emissions from biofuels and biomass burning at daily time-scale taken from GFAS
        - emissions from inland water (annual climatology, Johnson et al. 2022)
        - emissions from termites (annual climatology, Castaldi 2013)
        - emissions from the oceans (annual climatology, Weber 2021)
        - geological emissions (annual climatology, Etiope 2015) 

        The further precomputed simulations differ from the default case in that different inputs (primarily in terms of the emission scenarios) are applied:
        """
        pane1 = pn.pane.Markdown(
            info1,
            stylesheets=[precomp_stylesheet,], css_classes=['precomp-right']
            )
        pane2 = pn.pane.Markdown(
            f"{self._exptable_html()}",
            stylesheets=[precomp_stylesheet,], css_classes=['precomp-right']
        )
        pane3 = pn.pane.Markdown(
            """
            ## <ins>Statistics and Comparison</ins>

            Two principal types of statistics and comparison of simulation experiments against observed concentrations are available, which can be selected via the corresponding tab on the left-hand side:
            - **Fit statistics**
              - visualisation of statistics or concentrations for one selected experiment 
              - temporally averaged simulated or observed concentration differentiated by station (as map and as table)
              - metrics based on complete time-series differentiated by station
              - available metrics include bias, RMSE (Root-mean squared error), and Pearson correlation coefficient
            - **Modelled timeseries**
              - comparison/visualsation of observed and simulated concentrations (potentially from multiple experiments) for one selected site
              - time-series of hourly concentrations
              - time-series of weekly-averaged biases
              - histogram of hourly biases
              - table with metrics (RMSE, bias, correlation coeffient) based on complete time-series
            """,
            stylesheets=[precomp_stylesheet,], css_classes=['precomp-right']
            )
        return pn.Column(
            pane1,
            pane2,
            pane3,
            )

    def _exptable_html(self):
        def _get_desc(exp):
            desc = experiment_desc(exp)
            if desc is None:
                return desc
            desc_html = '<tr>'
            desc_html += f'<td>{exp}</td>'
            desc_html += f'<td>{desc}</td>'
            desc_html += f'</tr>'
            return desc_html
        
        exp_list = list(self.settings.experiments.list)
        tbl = f'<table>'
        #-- header
        tbl += f'<thead><tr><th>Experiment</th><th>Description</th></tr></thead>'
        for iexp,exp in enumerate(exp_list):
            desc = _get_desc(exp)
            if not desc is None:
                tbl += desc
                # if iexp<len(exp_list)-1:
                #     tbl += '<br>'
        tbl += f'</table>'

        return tbl
        
