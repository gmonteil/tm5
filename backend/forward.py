#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
import subprocess
from pathlib import Path
import shutil
from loguru import logger

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
# MVO-CHANGED: this is now already done in server.py
# if outpath.exists():
#     shutil.rmtree(outpath)
# outpath.mkdir(parents=True, exist_ok=True)
if not outpath.exists():
    msg = f"...expected output directory {str(outpath)} not found on system!"
    raise RuntimeError(msg)

(outpath / 'forward').symlink_to(datapath / 'forward')
(outpath / 'inversion').symlink_to(datapath / 'inversion')
(outpath / 'foj.nc').symlink_to(datapath / 'foj.nc')
(outpath / 'ftj.nc').symlink_to(datapath / 'ftj.nc')
(outpath / 'fe.nc').symlink_to(datapath / emfile)
(outpath / 'fesigma.nc').symlink_to(datapath / 'fesigma.nc')

#
#-- Fortran inversion system will write to this directory
#
ftn_outdir = 'output'
outdir = outpath / f'{ftn_outdir}'
outdir.mkdir(parents=True, exist_ok=True)
#
#-- Fortran system operates on existing files for I/O,
#   we want them in the output folder, but Fortran expects those in the
#   current run direcgt
#
shutil.copy(outpath / 'foj.nc', outdir / 'fc.nc')
shutil.copy(outpath / 'fe.nc', outdir / 'fepost.nc')
shutil.copy(outpath / 'foj.nc', outdir / 'fcpost.nc')
#-- and these files must be available in current working directory!
# (outpath / 'fc.nc').symlink_to(outdir / 'fc.nc')
# (outpath / 'fcpost.nc').symlink_to(outdir / 'fcpost.nc')
# (outpath / 'fepost.nc').symlink_to(outdir / 'fepost.nc')
#
#-- prefer "local" links (e.g. fc.nc -> output/fc.nc)
#
dir_fd = os.open(str(outpath), os.O_RDONLY)
for p in ['fc.nc', 'fepost.nc', 'fcpost.nc',]:
    src = f"{ftn_outdir}/{p}"
    dst = p
    os.symlink(src=src, dst=dst, dir_fd=dir_fd)


# Run the code
p = subprocess.run('./forward', check=True, capture_output=True, cwd=outpath)
with open(outpath / 'output/forward.txt', 'w') as fid:
    fid.writelines(p.stdout.decode())

if args.task == 'inversion':
    #-- MVO,20260729:
    #   - Fortran code writes prior and posterior concentration to file 'fcpost.nc'
    #   - there is no longer need to copy (and the copy to fc_apos.nc even did not
    #     contain the posterior but still the prior concentrations)
    #
    # shutil.copy(outpath / 'fc.nc', outpath / 'fc_apri.nc')
    p = subprocess.run('./inversion', check=True, capture_output=True, cwd=outpath)
    with open(outpath / 'output/inversion.txt', 'w') as fid:
        fid.writelines(p.stdout.decode())
        # shutil.move(outpath / 'fc.nc', outpath / 'fc_apos.nc')
