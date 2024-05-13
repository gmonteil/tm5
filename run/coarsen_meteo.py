#!/usr/bin/env python

"""
G. Monteil, 15 May 2024: Simplified run script to compute meteo on a coarsened grid 
"""

from omegaconf import OmegaConf
from pathlib import Path
import tm5
from tm5.system import runcmd

# 1. Build the model
tm = tm5.TM5('coarsen_meteo.yaml', host='donkey')
#tm.build()

# 2. setup input files
tm.setup_meteo()
#tm.setup_run('forward')
#tm.setup_iniconc()
#tm.setup_regions()
#tm.setup_system()
tm.settings['proces.source'] = 'F'

# Write the rc-file
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'forward.rc')

# Run TM5
runcmd(tm.dconf.run.run_cmd.split() + [str(rcf)])
#runcmd(tm.settings.run.cmd.split() + ['forward.rc'])