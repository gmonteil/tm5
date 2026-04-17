#!/usr/bin/env python

"""
G. Monteil, 15 May 2024: Basic script to do a forward run with TM5
"""

import sys
import os
import shutil
from omegaconf import OmegaConf
import tm5
from argparse import ArgumentParser
from pathlib import Path
from tm5.system import runcmd
from loguru import logger

def get_hostname():
    import socket
    hostname = socket.gethostname()
    return hostname

parser = ArgumentParser()
parser.add_argument('-b', '--build', action='store_true', default=False, help='Use this option to compile the code')
parser.add_argument('--build-only', action='store_true', default=False)
parser.add_argument('--build-dir',
                    help="""over-ride default build directory as it is defined in configuration file.""")
parser.add_argument('--rcfile-only', action='store_true', default=False, help="""only create the TM5 rcfile (nor compiling neither running TM5).""")
parser.add_argument('-m', '--host', default=os.environ['TM5_HOST'])
parser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
parser.add_argument('--expdir',
                    help="""over-ride top-level experiment (output) directory as defined in configuration file.""")
parser.add_argument('config_file')
args = parser.parse_args(sys.argv[1:])


yaml_file = Path(args.config_file)
if not yaml_file.exists():
    msg = f"provided yaml file ***{str(yaml_file)}*** not accessible!"
    raise RuntimeError(msg)

#=====================================================
# -1. read configuration, potentiall override part of the settings
#=====================================================
machine = args.host
platform = get_hostname()
dconf = OmegaConf.load(str(yaml_file))
#-- consistency
if not machine in dconf.keys():
    msg = f"selected machine -->{args.host}<-- not found in configuration file."
    raise RuntimeError(msg)
#-- override settings
if args.trange!=None:
    tstart, tend = args.trange
    dconf.run['start'] = tstart
    dconf.run['end']   = tend
if args.expdir!=None: #-- note, expdir is currently set in selected machine/host
    expdir_sav = dconf[args.host].paths.expdir
    msg = f"setting experiment directory -->{args.expdir}<-- " \
        f"(overriding configuration file ==>{expdir_sav}<=="
    logger.info(msg)
    dconf[args.host].paths.expdir = args.expdir
if args.build_dir!=None:
    dconf[machine].paths.build = args.build_dir
if platform!=None and 'platform' in dconf.run:
        if dconf.run.platform!=platform:
            msg = f"overriding dconf.run.platform, {dconf.run.platform} " \
                f"{platform}"
            logger.info(msg)
            dconf.run['platform'] = platform
# OmegaConf.save(config=dconf, f='xxx.yml')
# sys.exit(0)

#=====================================================
# 0. create pre-processor instance
#=====================================================
tm = tm5.TM5(dconf, host=args.host)

#=====================================================
# 1. Build the model
#=====================================================
if args.build or args.build_only and not args.rcfile_only:
    tm.build()

if args.build_only:
    logger.info(f"TM5 compilation done, exiting now")
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
logger.info(f"start emissions preparation...")
tm.setup_emissions2()
logger.info(f"...emissions done.")

#=====================================================
# Write the rc-file
#=====================================================
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'forward.rc')

#=====================================================
# copy originating yaml file
#=====================================================
dst = Path(tm.dconf.run.paths.output) / yaml_file.name
# #-- copy2 tries to preserve file attributes (as far as possible)
# shutil.copy2(args.config_file, str(dst))
OmegaConf.save(config=dconf, f=str(dst))

#=====================================================
# Run TM5
#=====================================================
run_cmd = tm.dconf.run.run_cmd.split() + [str(rcf)]
if args.rcfile_only:
    msg = f"TM5 rcfile ***{str(rcf)}*** has been created"
    logger.info(msg)
    msg = f"...skipping run command \n-->{' '.join(run_cmd)}<--\n"
    logger.info(msg)
else:
    runcmd(run_cmd)
logger.info(f"regular program termination")
