#!/usr/bin/env python
import sys
import os
from pathlib import Path
import pandas as pd
from omegaconf import OmegaConf
from argparse import ArgumentParser
from loguru import logger

from tm5.main import TM5

#-- command line arguments
parser = ArgumentParser(description='build and run TM5 atmosphere tracer model')
parser.add_argument('--yaml',
                    required=True,
                    # default='tm5.yaml',
                    help="""yaml configuration file for TM5 (!mandatory!) (default: %(default)s).""")
parser.add_argument('--run_mode',
                    choices=['only_dump', 'only_build', 'only_prepforward',
                             'only_forward', 'coarsen_meteo',],
                    # default='only_dump',
                    help="""selected run mode (default: %(default)s).""")
parser.add_argument('--iniconc',
                    choices=['zero'],
                    help="""initial concentration type (to override settings from the yaml file).""")
parser.add_argument('--trange',
                    nargs=2,
                    metavar=('yyyy-mm-dd_start', 'yyyy-mm-dd_end'),
                    help="""temporal range of simulation (overrides settings in yaml file), start/end of simulation period to be provided in format yyyy-mm-dd.""")
parser.add_argument('--ndyn',
                    type=int,
                    help="""length of time-step [s] of TM5 ('ndyn'), note that cfl.outputstep will be set to the same value.""")
#-- parse command line
args = parser.parse_args(sys.argv[1:])


if args.trange!=None:
    tstart, tend = args.trange
    msg = f"overriding simulation period to {tstart} -- {tend}"
    logger.info(msg)
    tstart = pd.Timestamp(tstart)
    tend   = pd.Timestamp(tend)
    yamlfile_org = Path(args.yaml)
    #-- load original yaml file
    dconf = OmegaConf.load(yamlfile_org)
    #-- override settings
    if tstart.hour==0:
        dconf.run.start = tstart.strftime('%Y-%m-%d')
    else:
        dconf.run.start = tstart.strftime('%Y-%m-%dT%H')
    if tend.hour==0:
        dconf.run.end   = tend.strftime('%Y-%m-%d')
    else:
        dconf.run.end   = tend.strftime('%Y-%m-%dT%H')
    #--
    bname = os.path.splitext(yamlfile_org)[0]
    yamlfile = f"{bname}_{tstart.strftime('%Y%m%d')}-{tend.strftime('%Y%m%d')}.yaml"
    with open(yamlfile, 'w') as fpout:
        OmegaConf.save(config=dconf, f=fpout.name)
    #-- consistency check
    loaded = OmegaConf.load(yamlfile)
    assert dconf==loaded, \
        f"unexpected differences in generated yaml configuration from ***{yamlfile}***"
    msg = "...generated yaml file ***{yamlfile}***"
    logger.info(msg)
else:
    yamlfile = args.yaml


#
#-- create TM5 instance from yaml file
#
tm = TM5(yamlfile, machine='ilab')


if args.run_mode=='only_dump':
    logger.info(f"dumping TM5 configuration only")
    print(f"{tm5.dconf}")
    sys.exit(0)


if args.run_mode=='coarsen_meteo':
    tm.coarsen_meteo()
    sys.exit(0)

#--
if args.run_mode in [None,'only_build']:
    tm.build()
    if args.run_mode=='only_build':
        sys.exit(0)
#--
if args.run_mode in [None, 'only_prepforward', 'only_forward']:
    #-- must ensure
    #   - TM5 executable exists
    #     (but potentially difficult to ensure it is consistent with current yaml)
    #   - output directory must be present, tm5
    # #MVDEBUG::
    # print(f"MVMV::tm.dconf -->{tm.dconf}<--")
    # print(f"MVMV::tm.dconf.run.paths.output -->{tm.dconf.run.paths.output}<--")
    # print(f"MVMV::tm.dconf.build.directory -->{tm.dconf.build.directory}<--")
    output = Path(tm.dconf.run.paths.output)
    tm5exe = output / 'tm5.x'
    if not tm5exe.exists():
        if not output.exists():
            logger.info(f"output directory -->{output}<-- does not yet exist.")
            output.mkdir(parents=True, exist_ok=True)
        #--
        logger.info(f"need to create link to TM5 executable.")
        build_dir = Path(tm.dconf.build.directory)
        tm5exec_built = build_dir / 'tm5.x'
        assert tm5exec_built.exists(), \
            f"TM5 exectuable -->{tm5exec_built}<-- not present"
        # print(f"-->{tm5exe}<--")
        # print(f"-->{tm5exec_built}<--   ==>{tm5exec_built.absolute()}<==")
        os.symlink(tm5exec_built.absolute(), tm5exe)
    dry_run = args.run_mode=='only_prepforward'
    ndyn = None if args.ndyn==None else str(args.ndyn)
    #-- for now taking default values selected by Guillaume,
    #   but potentially only do a "dry run"
    #   (i.e. preparing all inputs for the run but eventually only dump the command that would
    #    be executed.)
    #-- MVO-ADDED-20240104: allow to set zero initial concentrations (and bypass any input files)
    #-- MVO-ADDED-20240529: allow overriding 'ndyn' (which is applied to cfl.outputstep, too)
    tm.forward(iniconc=args.iniconc, dry_run=dry_run, ndyn=ndyn)
