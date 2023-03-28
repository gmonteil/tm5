#!/usr/bin/env python

from argparse import ArgumentParser
import tm5.setup
import tm5.build
from pathlib import Path
from omegaconf import OmegaConf
import subprocess
from loguru import logger


# tm5 --dev pyshell forward --legacy --build --rc forward.yaml

p = ArgumentParser()
p.add_argument('--build', action='store_true', default=False, help='compile tm5')
p.add_argument('--pyshell', action='store_true', default=False)
p.add_argument('action', choices=['build', 'forward', 'background', 'optim'])
p.add_argument('config_file', nargs=1)
args = p.parse_args()


#TODO: replace all the calls to "tm5" by direct calls to singularity (and delete the "tm5" script).


conf = OmegaConf.load(args.config_file[0])

if args.build or args.action == 'build':
    exec = tm5.build.build_tm5(conf)


if args.pyshell and args.action == 'forward' and False:

    # Setup TM5 using pyshell
    output = Path(conf.run.paths.output).parent
    cmd = f'tm5 --output {output} --dev pyshell setup --rc {args.config_file[0]}'
    subprocess.run(cmd.split())

    # Run TM5
    tm5.run.run_tm5([exec, Path(conf.run.paths.output) / 'tm5.rc'], output='.')


if args.pyshell and args.action == 'background':
    #TODO: this one looks wrong, needs testing ...

    # Setup TM5 using pyshell
    output = Path(conf.run.paths.output).parent
    cmd = f'tm5 --output {output} --dev pyshell background --rc {args.config_file[0]}'
    subprocess.run(cmd.split())

    # Run TM5
    tm5.run.run_tm5([exec, Path(conf.run.paths.output) / 'tm5.rc'], output=output)


if args.pyshell and args.action == 'optim':

    # Setup TM5 using pyshell
    output = Path(conf.run.paths.output).parent
    # cmd = f'tm5 --output {output} --dev pyshell optim --rc {args.config_file[0]}'
    cmd = 'singularity run '\
        f'--bind {conf.machine.paths.meteo}:/meteo '\
        f'--bind {output}:/output '\
        f'--bind {conf.machine.paths.input}:/input '\
        f'--bind {conf.machine.paths.scratch}:/scratch '\
        f'--bind {conf.machine.paths.data}:/data '\
        f'--bind {conf.machine.paths.tm5}:/tm5 '\
        f'{conf.machine.paths.container} pyshell optim --rc {args.config_file[0]}'
    logger.info(cmd)
    
    tm5.setup.setup_tm5(conf)
    subprocess.run(cmd.split())


if args.pyshell and args.action == 'forward':

    # Setup TM5 using pyshell
    output = Path(conf.run.paths.output).parent
    # cmd = f'tm5 --output {output} --dev pyshell optim --rc {args.config_file[0]}'
    cmd = 'singularity run '\
        f'--bind {conf.machine.paths.meteo}:/meteo '\
        f'--bind {output}:/output '\
        f'--bind {conf.machine.paths.input}:/input '\
        f'--bind {conf.machine.paths.scratch}:/scratch '\
        f'--bind {conf.machine.paths.data}:/data '\
        f'--bind {conf.machine.paths.tm5}:/tm5 '\
        f'{conf.machine.paths.container} pyshell forward --rc {args.config_file[0]}'
    logger.info(cmd)
    
    tm5.setup.setup_tm5(conf)
    subprocess.run(cmd.split())
