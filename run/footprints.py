#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
from omegaconf import OmegaConf
from pathlib import Path
from loguru import logger
from pandas import Timestamp, Timedelta
import tm5
from tm5.system import runcmd
from tm5.fitic import read_obs_table, create_departure_files

try:
    default_host = os.environ['TM5_HOST']
except KeyError:
    default_host = None
                
parser = ArgumentParser()
parser.add_argument('-b', '--build', action='store_true', default=False, help='Use this option to compile the code')
parser.add_argument('--build-dir',
                    help="""over-ride default build directory as it is defined in configuration file.""")
parser.add_argument('--build-only', action='store_true', default=False)
parser.add_argument('-m', '--host', default=default_host)
parser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
parser.add_argument('--ndays',
                    type=int,
                    help="""number of days to simulate (ATTENTION: overrules 'tend' provided by --trange).""")
parser.add_argument('--obsfile',
                    help="""whether to override the observations file used to create the departures for the adjoint run.""")
parser.add_argument('--outdir',
                    help="""override top-level experiment (output) directory as defined in configuration file.""")
parser.add_argument('--tm5exec',
                    help="""specify an already compiled TM5 executablea as an *absolute* path, user is responsible to ensure it is consistent with the spatial settings of the provided yaml configuration file. NOTE, that the provided executable will *not* be used in case one of the options '--build' or --build-only' war provided, too!""")
parser.add_argument('config_file')
args = parser.parse_args(sys.argv[1:])

#-- check yaml file exists
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
if args.ndays!=None:
    dend = Timestamp(dconf.run.start) + Timedelta(days=args.ndays)
    dconf.run['end'] = dend.strftime('%Y%m%d')
if args.obsfile!=None:
    obsfile_sav = dconf.observations.file
    msg = f"setting observations file -->{args.obsfile}<-- " \
        f"(overriding entry in configuration file ==>{obsfile_sav}<=="
    logger.info(msg)
    dconf.observations.file = args.obsfile
if args.build_dir!=None:
    dconf[args.host].paths.build = args.build_dir
if args.outdir!=None:
    outdir_sav = dconf[args.host].paths.output
    msg = f"setting experiment directory -->{args.outdir}<-- " \
        f"(overriding configuration file ==>{outdir_sav}<=="
    logger.info(msg)
    dconf[args.host].paths.output = args.outdir
if args.tm5exec!=None:
    if args.build or args.build_only:
        msg = f"discarding provided executable ***{args.tm5exec}*** " \
            f"because request for building it is also triggered via command line!"
        logger.warn(msg)
    else:
        exe = Path(args.tm5exec)
        if not exe.is_absolute():
            msg = f"provided executable must be an absolute path (-->{str(exe)}<--)"
            raise RuntimeError(msg)
        elif not exe.exists():
            msg = f"provided executable not found on sysstem (-->{str(exe)}<--)"
            raise RuntimeError(msg)
        else:
            dconf.run.paths['tm5exec'] = str(exe)
#
#-- adapt obsfile
#
if dconf.observations.file.find('%Y-%m-%d')>0:
    date_end = Timestamp(dconf.run.end) #-- end-of-simulation (yyyymmddT00:00:00)
    date_obs = date_end - Timedelta(days=1)
    dconf.observations.file = date_obs.strftime(dconf.observations.file)

# =====================================================
# 0. create configuration instance
# =====================================================
tm = tm5.TM5(dconf, host=args.host)


# =====================================================
# 1. Build the model
# =====================================================
if args.build or args.build_only:
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
dst = Path(tm.dconf.run.paths.output) / 'tm5.yaml'
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
