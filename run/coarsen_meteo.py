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
parser.add_argument('--config_file', default='coarsen_meteo.yaml')
args = parser.parse_args(sys.argv[1:])

# keys = sorted(os.environ.keys())
# for k in keys:
#     print(f"{k:<20}: {os.environ[k]}")
try:
    platform = os.environ['HOSTNAME']
except KeyError:
    platform = 'unknown'

# 1. Build the model
tm = tm5.TM5(args.config_file, host=args.host, platform=platform)
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
#-- MVO-20250205: when running this coarsen_meteo.py I meanwhile got
#                 the error messages below at TM5 run-time about missing
#                 emissions keys,
#                 thus I inserted the setup of emissions:
tm.setup_emissions2(skip_emis_gen=True)
      # ERROR - key not found and no default specified ...
      # ERROR -   rcfile : /home/mivo/work/TM5/expdir-coarsen-meteo_tm5-iLab/coarsen-meteo-tropo34-fitic_coarsened-False_Makefile.singularity.ifort_platform-cosmos1.int.lunarc/output_2021-01-01--2021-01-03/forward.rc
      # ERROR -   key    : emissions.CH4.glb600x400.ncats
      # ERROR - in GO_Rc/ReadRc_i (go_rc.F90, line  424)
      # ERROR - in Emission_Data/Emission_Data_Init (emission_data.F90, line  185)
      # ERROR - in Emission/Emission_Init (emission.F90, line  103)
      # ERROR - in ModelIntegration/Proces_Init (modelIntegration.F90, line  182)
      # ERROR - in Tracer/Tracer_Model (tracer.F90, line  281)
      # ERROR - in TM5var4D/TM5var4D_Run (tm5var4d.F90, line  312)

tm.settings['proces.source'] = 'F'

# Write the rc-file
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'forward.rc')

# Run TM5
runcmd(tm.dconf.run.run_cmd.split() + [str(rcf)])
