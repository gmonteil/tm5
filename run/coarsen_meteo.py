#!/usr/bin/env python

"""
G. Monteil, 15 May 2024: Simplified run script to compute meteo on a coarsened grid 
"""
import os

from omegaconf import OmegaConf
from pathlib import Path
import tm5
from tm5.system import runcmd
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument('-b', '--build', action='store_true', default=False, help='Use this option to compile the code')
parser.add_argument('-m', '--host', default=os.environ['TM5_HOST'])
args = parser.parse_args(sys.argv[1:])

# 1. Build the model
tm = tm5.TM5('coarsen_meteo.yaml', host=args.host)
if args.build :
    tm.build()

# 2. setup input files
tm.setup_meteo()
tm.setup_run('forward')
tm.setup_iniconc('zero')
tm.setup_regions()
tm.setup_system()
tm.setup_tracers()
tm.setup_output()
tm.settings['proces.source'] = 'F'

# Write the rc-file
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'forward.rc')

# Run TM5
runcmd(tm.dconf.run.run_cmd.split() + [str(rcf)])