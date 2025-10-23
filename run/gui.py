#!/usr/bin/env python

from tm5.gui.main import FitIC_UI
from tm5.gui.widgets.stations import StationExplorer, StatisticsViewer
import panel as pn
from omegaconf import OmegaConf
pn.extension()
pn.extension('terminal')

# Setup TM5 tab
setup_gui = FitIC_UI()


# Results tab (stations explorer)

# Normally the settings should be passed from the gui, but for now, I just hardcoded them ...
conf = OmegaConf.create()
conf.observations = {}
conf.observations.files = '/home/gmonteil/data/iLab/TM5/observations/ch4_*.nc'
conf.experiments = {}
conf.experiments.list = {
     'default': 'fitic-simu-default_oh-cams/output_2021-01-01--2022-01-01',
     'edgarflat': 'fitic-simu-edgarflat_oh-cams/output_2021-01-01--2022-01-01',
     'half-oh': 'fitic-simu-half-oh_oh-cams/output_2021-01-01--2022-01-01',
     'no-germany': 'fitic-simu-no-germany_oh-cams/output_2021-01-01--2022-01-01',
     'no-gns': 'fitic-simu-no-gns_oh-cams/output_2021-01-01--2022-01-01',
     'no-northamerica': 'fitic-simu-no-northamerica_oh-cams/output_2021-01-01--2022-01-01',
     'regional_no-agri': 'fitic-simu-regional_anthro-no-agri_oh-cams/output_2021-01-01--2022-01-01',
     'regional_no-fossil': 'fitic-simu-regional_anthro-no-fossil_oh-cams/output_2021-01-01--2022-01-01',
     'regional_no-france': 'fitic-simu-regional_anthro-no-france_oh-cams/output_2021-01-01--2022-01-01',
     'regional_no-netherlands': 'fitic-simu-regional_anthro-no-netherlands_oh-cams/output_2021-01-01--2022-01-01',
     'regional_no-waste': 'fitic-simu-regional_anthro-no-waste_oh-cams/output_2021-01-01--2022-01-01',
     'regional': 'fitic-simu-regional_oh-cams/output_2021-01-01--2022-01-01'
}
conf.experiments.path = '/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/runs_with-convdiff_repeated'
conf.experiments.cache = '/lunarc/nobackup/projects/ghg_inv/guillaume/iLab/TM5/experiments'

stations = StationExplorer(conf)
station_statistics = StatisticsViewer(conf)

gui = pn.Tabs(
    ("Setup simulation", setup_gui),
    ("Results", pn.Tabs(
        ('Fit statistics', station_statistics),
        ('Modelled timeseries', stations), 
        tabs_location='left'),
    ),
    dynamic=True
    )
gui.servable()