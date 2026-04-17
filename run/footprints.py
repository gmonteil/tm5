#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
from omegaconf import OmegaConf
from pathlib import Path
from loguru import logger
import tm5
from tm5.system import runcmd
from tm5.fitic import read_obs_table, create_departure_files

                
parser = ArgumentParser()
parser.add_argument('-b', '--build', action='store_true', default=False, help='Use this option to compile the code')
parser.add_argument('--build-only', action='store_true', default=False)
parser.add_argument('-m', '--host', default=os.environ['TM5_HOST'])
parser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
parser.add_argument('--obsfile',
                    help="""whether to override the observations file used to create the departures for the adjoint run.""")
parser.add_argument('--outdir',
                    help="""override top-level experiment (output) directory as defined in configuration file.""")
parser.add_argument('config_file')
args = parser.parse_args(sys.argv[1:])

yaml_file = Path(args.config_file)
if not yaml_file.exists():
    msg = f"provided yaml file ***{str(yaml_file)}*** not accessible!"
    raise RuntimeError(msg)

#=====================================================
# -1. read configuration, potentiall override part of the settings
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
if args.obsfile!=None:
    obsfile_sav = dconf.observations.file
    msg = f"setting observations file -->{args.obsfile}<-- " \
        f"(overriding entry in configuration file ==>{obsfile_sav}<=="
    logger.info(msg)
    dconf.observations.file = args.obsfile
if args.outdir!=None:
    outdir_sav = dconf[args.host].paths.output
    msg = f"setting experiment directory -->{args.outdir}<-- " \
        f"(overriding configuration file ==>{outdir_sav}<=="
    logger.info(msg)
    dconf[args.host].paths.output = args.outdir

# =====================================================
# 1. Build the model
# =====================================================
tm = tm5.TM5(dconf, host=args.host)
if args.build :
    tm.build()

if args.build_only:
    sys.exit()
    

# =====================================================
# 2. Setup the observations
# =====================================================
observations_table = read_obs_table(tm.dconf.observations.file)
tm.setup_observations(observations_table)
   
    
# =====================================================
# 3. Setup input files:
# =====================================================

# Fill in all the meteo-related rc-keys + download the meteo files (if needed)
tm.setup_meteo()
tm.setup_run('forward')
tm.setup_regions()
tm.setup_system()
tm.setup_tracers()
tm.setup_output()
tm.setup_iniconc()
logger.info("start emissions preparation...")
tm.setup_emissions2()
logger.info("...emissions done.")


# =====================================================
# Write the rc-file
# =====================================================
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'forward.rc')

#=====================================================
# copy originating yaml file
#=====================================================
dst = Path(tm.dconf.run.paths.output) / yaml_file.name
OmegaConf.save(config=dconf, f=str(dst))

# =====================================================
# Do a forward run
# =====================================================
runcmd(tm.dconf.run.run_cmd.split() + [str(rcf)])


# =====================================================
# Setup the adjoint run
# =====================================================

# Create a departure file
create_departure_files(tm.dconf)

# Setup one tracer per observation:
tm.setup_adjoint_tracers(observations_table)

# Set the main keys to adjoint mode 
tm.setup_run('adjoint')

# Setup the emissions: we still need some emission rc-keys for each new tracer
# This could probably be simplified since we write a simplified adjoint structure, but for now it works ...
tm.setup_emissions2(skip_emis_gen=True)


# =====================================================
# Do an adjoint run
# =====================================================

rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'adjoint.rc')
runcmd(tm.dconf.run.run_cmd.split() + [str(rcf)])
