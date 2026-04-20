#!/usr/bin/env python

# Prepare daily emission files for the requested grid
import os
import sys
from omegaconf import OmegaConf
import xarray as xr
from pandas import Timestamp, date_range
from tm5.settings import load_config
from pathlib import Path
from loguru import logger
from argparse import ArgumentParser
from tm5.gridtools import TM5Grids
from tm5.emissions import prepare_emissions, coarsen_file

try:
    TM5_HOST = os.environ['TM5_HOST']
except KeyError:
    TM5_HOST = None

parser = ArgumentParser()
parser.add_argument('--config_file', default='forward.yaml')
parser.add_argument('-m', '--host', default=TM5_HOST)
parser.add_argument('--expdir',
                    help="""override target dir from host section in yaml file.""")
parser.add_argument('--outdir',
                    help="""override output directory of selected machine in yaml file.""")
parser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
args = parser.parse_args(sys.argv[1:])

if args.host==None:
    msg = f"host must be specified, either via TM5_HOST environment " \
        "variable or by using option --host. Cannot continue!"
    raise RuntimeError(msg)

yaml_file = Path(args.config_file)
if not yaml_file.exists():
    msg = f"provided yaml file ***{str(yaml_file)}*** not accessible!"
    raise RuntimeError(msg)

#=====================================================
# 1. read configuration, potentiall override part of the settings, set 'host' key
#=====================================================
dconf = OmegaConf.load(str(yaml_file))
#-- consistency
if not args.host in dconf.keys():
    msg = f"selected machine -->{args.host}<-- not found in configuration file."
    raise RuntimeError(msg)
#
#-- override settings
#
if args.trange!=None:
    tstart, tend = args.trange
    dconf.run['start'] = tstart
    dconf.run['end']   = tend
if args.outdir!=None:
    outdir_sav = dconf[args.host].paths.output
    msg = f"setting experiment directory -->{args.outdir}<-- " \
        f"(overriding configuration file ==>{outdir_sav}<=="
    logger.info(msg)
    dconf[args.host].paths.output = args.outdir
dconf['host'] = dconf[args.host]

#=====================================================
# copy originating yaml file
#=====================================================
dst = Path(dconf.run.paths.output) / yaml_file.name
dst.parent.mkdir(parents=True, exist_ok=True)
OmegaConf.save(config=dconf, f=str(dst))
msg = f"...stored yaml configuration file ***{str(dst)}***"
logger.info(msg)

#=====================================================
# 2. prepare (daily) emissions files 
#=====================================================
prepare_emissions(dconf)
