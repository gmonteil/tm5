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

#
#-- required links from the datapath
#   - both exectuables
#   - observational and target Jacobian
#   - emissions and emission uncertainties
#
(outpath / 'forward').symlink_to(datapath / 'forward')
(outpath / 'inversion').symlink_to(datapath / 'inversion')
(outpath / 'foj.nc').symlink_to(datapath / 'foj.nc')
(outpath / 'ftj.nc').symlink_to(datapath / 'ftj.nc')
(outpath / 'fe.nc').symlink_to(datapath / emfile)
(outpath / 'fesigma.nc').symlink_to(datapath / 'fesigma.nc')
#
#-- Fortran inversion system outputs go to dedicated subdirectory
#
ftn_outdir = 'output'
outdir = outpath / f'{ftn_outdir}'
outdir.mkdir(parents=True, exist_ok=True)

#
#-- required NetCDF inputs depend on task
#
if args.task=='forward':
    mapping_dict = {
        'fc.nc':     'foj.nc',
    }
elif args.task=='inversion':
    mapping_dict = {
        'fepost.nc': 'fe.nc',
        'fcpost.nc': 'foj.nc',
    }


#
#-- Fortran system operates on existing files for I/O,
#   we want them in the output folder,
#   but the Fortran code expects those in the current run directory
#   and we set symbolic links
#
dir_fd = os.open(str(outpath), os.O_RDONLY)
for dstfile,srcfile in mapping_dict.items():
    shutil.copy(outpath / srcfile, outdir / dstfile)
    #
    #-- prefer "local" links in run directory (e.g. fc.nc -> output/fc.nc)
    #
    srcfile_loc = (Path(ftn_outdir) / dstfile)
    os.symlink(src=str(srcfile_loc), dst=dstfile, dir_fd=dir_fd)


# Run the code
if args.task=='forward':
    p = subprocess.run('./forward', check=True, capture_output=True, cwd=outpath)
    with open(outpath / 'output/forward.txt', 'w') as fid:
        fid.writelines(p.stdout.decode())

elif args.task == 'inversion':
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
