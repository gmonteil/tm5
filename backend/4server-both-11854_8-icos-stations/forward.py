#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import subprocess
from pathlib import Path
import shutil


parser = ArgumentParser()
parser.add_argument('--output', help='Path where the code will be run and the output written', type=Path)
parser.add_argument('--data', help='Path where the input files are located. Should be visible by both the VM and the JHub', type=Path)
parser.add_argument('--emis', help='Name of the file to use for the emissions (should be inside the data path)')
parser.add_argument('--task', choices=['forward', 'inversion'])
args = parser.parse_args(sys.argv[1:])

outpath = args.output
datapath = args.data
emfile = args.emis

# Create the run directory and copy the files in it
if outpath.exists():
    shutil.rmtree(outpath)
outpath.mkdir(parents=True, exist_ok=True)

(outpath / 'forward').symlink_to(datapath / 'forward')
(outpath / 'inversion').symlink_to(datapath / 'inversion')
(outpath / 'foj.nc').symlink_to(datapath / 'foj.nc')
(outpath / 'ftj.nc').symlink_to(datapath / 'ftj.nc')
(outpath / 'fe.nc').symlink_to(datapath / emfile)
(outpath / 'fesigma.nc').symlink_to(datapath / 'fesigma.nc')
shutil.copy(datapath / 'foj.nc', outpath / 'fc.nc')
shutil.copy(datapath / 'foj.nc', outpath / 'fepost.nc')
shutil.copy(datapath / 'foj.nc', outpath / 'fcpost.nc')

(outpath / 'output').mkdir()

# Run the code
p = subprocess.run('./forward', check=True, capture_output=True, cwd=outpath)
with open(outpath / 'output/forward.txt', 'w') as fid:
    fid.writelines(p.stdout.decode())

if args.task == 'inversion':
    shutil.copy(outpath / 'fc.nc', outpath / 'fc_apri.nc')
    p = subprocess.run('./inversion', check=True, capture_output=True, cwd=outpath)
    with open(outpath / 'output/inversion.txt', 'w') as fid:
        fid.writelines(p.stdout.decode())
        shutil.move(outpath / 'fc.nc', outpath / 'fc_apos.nc')
