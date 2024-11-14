#!/usr/bin/env python

"""
G. Monteil, 15 May 2024: Basic script to do a forward run with TM5
"""

import os
from omegaconf import OmegaConf
import tm5
from argparse import ArgumentParser
import sys
from pathlib import Path
from tm5.system import runcmd
from loguru import logger

parser = ArgumentParser()
parser.add_argument('-b', '--build', action='store_true', default=False, help='Use this option to compile the code')
parser.add_argument('--build-only', action='store_true', default=False)
parser.add_argument('-m', '--host', default=os.environ['TM5_HOST'])
parser.add_argument('config_file')
args = parser.parse_args(sys.argv[1:])

#=====================================================
# 1. Build the model
#=====================================================
tm = tm5.TM5(args.config_file, host=args.host)
if args.build :
    tm.build()

if args.build_only:
    sys.exit()

#=====================================================
# 2. Setup input files:
#=====================================================

# Fill in all the meteo-related rc-keys + download the meteo files (if needed)
tm.setup_meteo()

# Fill in the rc-keys related to run type and duration
tm.setup_run('forward')

# Set the keys related to regions
#TODO: check that the redgrid keys are handled ok!
tm.setup_regions()

# Set machine-dependent rc-keys?? ==> should be merged with something else ...
tm.setup_system()

# Set rc-keys for tracers (for the tracer-generic chem_params.F90 module)
tm.setup_tracers()

# Set the keys related to TM5 "output" modules (observations, ...)
tm.setup_output()

# Set the keys related to initial condition
tm.setup_iniconc()

# Set the emissions
#tm.setup_emissions()
tm.setup_emissions2()

#=====================================================
# Write the rc-file
#=====================================================
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'forward.rc')

#=====================================================
# Run TM5
#=====================================================
runcmd(tm.dconf.run.run_cmd.split() + [str(rcf)])

# 2. setup input files
#tm.setup_emissions()
#tm.setup_meteo()
#tm.setup_observations()
#tm.setup_tm5_optim()
#tm.setup_run('forward')
#tm.setup_output()
#tm.setup_regions()
#tm.setup_iniconc()
#tm.setup_optim()
#tm.setup_tracers()
#tm.setup_system()
