#!/usr/bin/env python

from argparse import ArgumentParser
import tm5.setup
import tm5.run
import tm5.build
from pathlib import Path
from omegaconf import OmegaConf
import subprocess


p = ArgumentParser()
p.add_argument('--build', action='store_true', default=False, help='compile tm5')
p.add_argument('--pyshell', action='store_true', default=False)
p.add_argument('action', choices=['build', 'forward', 'background', 'optim'])
p.add_argument('config_file', nargs=1)
args = p.parse_args()


conf = OmegaConf.load(args.config_file[0])

if args.build or args.action == 'build':
    exec = tm5.build.build_tm5(conf)

tm5.setup.setup_tm5(conf)

if args.pyshell and args.action == 'background':
    # output = Path(conf.run.paths.output).parent
    tm5.run.run_tm5(f'pyshell background --rc {args.config_file[0]}', settings=conf.machine.host)


if args.pyshell and args.action == 'optim':
    # output = Path(conf.run.paths.output).parent
    tm5.run.run_tm5(f'pyshell optim --rc {args.config_file[0]}', settings=conf.machine.host)


if args.pyshell and args.action == 'forward':
    # output = Path(conf.run.paths.output).parent
    tm5.run.run_tm5(f'pyshell forward --rc {args.config_file[0]}', settings=conf.machine.host)
