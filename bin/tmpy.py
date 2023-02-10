#!/usr/bin/env python

from argparse import ArgumentParser
from tm5 import build as tm5
from pathlib import Path
from omegaconf import OmegaConf
import subprocess


# tm5 --dev pyshell forward --legacy --build --rc forward.yaml

p = ArgumentParser()
p.add_argument('--build', action='store_true', default=False, help='compile tm5')
p.add_argument('--pyshell', action='store_true', default=False)
p.add_argument('action', choices=['build', 'forward', 'background'])
p.add_argument('config_file', nargs=1)
args = p.parse_args()

conf = OmegaConf.load(args.config_file[0])

if args.build or args.action == 'build':
    exec = tm5.build_tm5(conf)

if args.pyshell and args.action == 'forward':

    # Setup TM5 using pyshell
    output = Path(conf.run.paths.output).parent
    cmd = f'tm5 --output {output} --dev pyshell setup --rc {args.config_file[0]}'
    subprocess.run(cmd.split())

    # Run TM5
    tm5.run_tm5([exec, Path(conf.run.paths.output) / 'tm5.rc'], output='.')

if args.pyshell and args.action == 'background':

    # Setup TM5 using pyshell
    output = Path(conf.run.paths.output).parent
    cmd = f'tm5 --output {output} --dev pyshell background --rc {args.config_file[0]}'
    subprocess.run(cmd.split())

    # Run TM5
    tm5.run_tm5([exec, Path(conf.run.paths.output) / 'tm5.rc'], output=output)

