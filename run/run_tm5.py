#!/usr/bin/env python

"""
G. Monteil, 15 May 2024: Basic script to do a forward run with TM5
"""

from omegaconf import OmegaConf
import tm5

# 1. Build the model
tm = tm5.TM5('forward_co2.yaml')
tm.build()

# 2. setup input files
tm.setup_emissions()
tm.setup_meteo()
tm.setup_observations()
tm.setup_tm5_optim()
tm.setup_run('forward')
tm.setup_output()
tm.setup_regions()
tm.setup_iniconc()
tm.setup_optim()
tm.setup_tracers()
tm.setup_system()