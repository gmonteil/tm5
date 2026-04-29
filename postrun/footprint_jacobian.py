#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
import datetime as dtm
from omegaconf import OmegaConf, DictConfig
from pathlib import Path
import datetime as dtm
from loguru import logger
from pandas import date_range, DataFrame
from pandas import Timestamp, Timedelta, concat
import xarray as xr
import numpy as np
from numpy import zeros, tile
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots,colorbar
from cartopy import crs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from tm5.fitic import read_obs_table
from tm5.post.footprint_io import tm5_fitic_adjoint_corrected_halos
from tm5.post.footprint_io import tm5_fitic_footprint4jacobian_v1
from tm5.post.footprint_io import tm5rundir_obsids_extra
from tm5.post.footprint_io import tm5rundir_emissions1D
from tm5.post.footprint_io import tm5rundir_jacobian2D
from tm5.gridtools import TM5Grids
from tm5.post.plot_util import cnorm_set
from tm5.post.utilities import lonstr,latstr



def subcmd_testbuild_jacobian_1day( args ):
    outpath_tm5 = args.outpath_tm5
    obsid       = args.obsid
    
    # footp_info = tm5_fitic_adjoint_corrected_halos(outpath_tm5)
    # print(footp.loc[:,'itime'].unique())
    if args.trange!=None:
        trange = date_range(args.trange[0], args.trange[1], freq='1d')
        print(trange)
    else:
        trange = args.trange
    #
    #-- load obs table used for this simulation
    #
    obs_table = tm5rundir_obsids_extra(outpath_tm5)
    if not obsid in obs_table:
        msg = f"selected observation identifier -->{obsid}<-- not present in simulation " \
            f"(==>{list(obs_table.keys())}<==)"
        raise RuntimeError(msg)
    obs_info = obs_table[obsid]
    print(obs_info)
    #
    #-- load emissions
    #
    msg = f"start preparing emissions..."
    logger.info(msg)
    emis1D = tm5rundir_emissions1D(outpath_tm5, trange=trange)
    msg = f"...emissions done."
    logger.info(msg)
    #
    #-- load Jacobian
    #
    # jacobian_info = tm5_fitic_footprint4jacobian_v1(outpath_tm5, trange=trange, obsid=obsid)
    # jac_table = jacobian_info.jac_table
    # days = jacobian_info.days #-- [emission days]
    # nday = len(days)
    # ngc = 0
    # for region,region_jac in jac_table.items():
    #     iday = len(region_jac.coords['iday'])
    #     nlon = len(region_jac.coords['lon'])
    #     nlat = len(region_jac.coords['lat'])
    #     if iday!=nday:
    #         msg = f"unexpected dimension size iday={iday} (expected {nday})"
    #         raise RuntimeError(msg)
    #     ngc += nlat*nlon #_ngc
    #     for iday,day in enumerate(days):
    #         _jac = region_jac.sel(iday=iday)
    #         msg = f"@{region},{day.strftime('%Y-%m-%d')}: jac min/mean/max = " \
    #             f"{_jac.min().values}/{_jac.mean().values}/{_jac.max().values}"
    #         print(msg)
    # msg = f"overall number of grid-cells ngc={ngc}"
    # print(msg)
    # #
    # #-- reformat Jacobian values to 1D array
    # #   ordering:
    # #   - emisday
    # #     - grid-cells glb600x400 (flattened)
    # #     - grid-cells eur300x200 (flattened)
    # #     - grid-cells gns100x100 (flattened)
    # #
    # jacobian2D = np.array([]) #empty()
    # for iday in range(nday):
    #     for reg,reg_jac in jac_table.items():
    #         dayreg_jac = reg_jac.values[iday,:,:]
    #         jacobian2D = np.concat( (jacobian2D, dayreg_jac.ravel()))
    #     # print(f"@iday={iday}, jacobian2D.shape={jacobian1D.shape}")
    #-- consistency check
    # assert emis1D.shape==jacobian1D.shape, \
    #     f"{emis1D.shape} vs {jacobian1D.shape}"
    # #
    # #--
    # #
    # jacobian2D = jacobian1D.reshape(1,len(jacobian1D))
    #
    jac_info = tm5rundir_jacobian2D(outpath_tm5, trange=trange, obsid=obsid)
    jacobian2D = jac_info.jac2D
    nobs,nemis = jacobian2D.shape
    # print(f"jacobian2D.shape={jacobian2D.shape}")
    for iobs in range(nobs):
        msg = f"@iobs={iobs}, jac min/mean/max = " \
            f"{jacobian2D[iobs,:].min()}/{jacobian2D[iobs,:].mean()}/{jacobian2D[iobs,:].max()}"
        print(msg)
    assert len(emis1D)==nemis, \
        f"expected {len(emis1D)} emissions, but Jacobian shape is {jacobian2D.shape}"
    cobs = np.dot(jacobian2D,emis1D)
    print(f"cobs -->{cobs}<--")


def subcmd_testbuild_jacobian_1month(args):
    outpath_tm5 = Path(args.outpath_tm5)
    month = args.month
    obsid = args.obsid
    #--
    dayf = Timestamp(f"2021-{month:02d}-01")          #-- 
    dayl = (dayf + Timedelta(days=31)).replace(day=1) - Timedelta(days=1) #-- last day in month
    #
    #-- naming pattern used by Guillaume: footprints_gns100x100_%Y%m%d,
    #   where '%Y%m%d' refers to last day of simulation at 0:00 !
    # e.g. for February we need to pick directories
    # footprints_gns100x100_20210202 to footprints_gns100x100_20210301
    odayf = dayf + Timedelta(days=1)
    odayl = dayf.replace(month=month+1)
    obs_trange = date_range(odayf, odayl, freq='1d')
    flx_trange = date_range(dayf, dayl, freq='1d')
    nobs = len(obs_trange)
    nday = len(flx_trange)  #-- number of emission days
    #
    #-- load Jacobians (differentiated by zoom region) for each observational time
    #
    jacobian_list = []
    for iobs,obsday in enumerate(obs_trange):
        msg = f"obs@{obsday.strftime('%Y-%m-%d')}, starting..."
        logger.info(msg)
        rundir = outpath_tm5 / f"footprints_gns100x100_{obsday.strftime('%Y%m%d')}"
        if not rundir.is_dir():
            msg = f"...expected directory -->{str(rundir)}<-- not found!"
            raise RuntimeError(msg)
        #
        #--
        #
        if iobs==1:
            #
            #-- load emissions
            #
            msg = f"start preparing emissions..."
            logger.info(msg)
            emis1D = tm5rundir_emissions1D(rundir, trange=flx_trange)
            msg = f"...emissions done."
            logger.info(msg)
        # #
        # #--
        # #
        # jac_info = tm5_fitic_footprint4jacobian_v1(rundir, trange=flx_trange, obsid=obsid)
        # jac_table = jac_info.jac_table
        # #
        # #-- consistency checks
        # #
        # assert list(jac_table.keys())==['glb600x400','eur300x200','gns100x100',], \
        #             f"unexpected regions in jacobian table -->{jac_table.keys()}<--"
        # for reg,reg_jac in jac_table.items():
        #     assert reg_jac.dims==('iday','lat','lon'), \
        #         f"unexpected dimensions in jacobian table -->{jac_table.dims}<--"
        # #
        # #-- reformat Jacobian values to 1D in spatial dimensions
        # #
        # jacobian_a = np.array([]) #empty()
        # for iday in range(nday):
        #     for reg,reg_jac in jac_table.items():
        #         dayreg_jac = reg_jac.values[iday,:,:]
        #         jacobian_a = np.concat( (jacobian_a, dayreg_jac.ravel()))
        #     print(f"@iday={iday}, jacobian_a.shape={jacobian_a.shape}")
        #
        # jacobian_list.append(jacobian_a)
        #--
        #
        jac_info = tm5rundir_jacobian2D(rundir, trange=flx_trange, obsid=obsid)
        jacobian_list.append(jac_info.jac2D)
    #
    #--
    #
    jac2D = np.vstack(jacobian_list)
    print(f"jac2D.shape={jac2D.shape}")

    csim = np.dot(jac2D, emis1D)
    print(f"csim -->{csim}<--")

################################################################################
#
#                   p a r s e r
#
parser = ArgumentParser(usage=globals()['__doc__'])
parser = ArgumentParser()
#----------------------------
#     s u b c o m m a n d s
#
subparsers = parser.add_subparsers( title='Available Subcommands',
                                    metavar='CMDS',
                                    description='',
                                    dest='subcmds',
                                    help='')
#
#--       testbuild_jacobian_1day
#
sparser = subparsers.add_parser('testbuild_jacobian_1day')
sparser.add_argument('outpath_tm5',
                     help="""top-level directory of TM5 adjoint run for footprint creation.""")
sparser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
sparser.add_argument('--obsid',
                     default='cbw_207',
                     help="""restrict Jacobian to one single observational location (default: %(default)s).""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")
#
#--       testbuild_jacobian_1month
#
sparser = subparsers.add_parser('testbuild_jacobian_1month')
sparser.add_argument('outpath_tm5',
                     help="""top-level directory of series of TM5 adjoint runs for footprint creation each of those for one observation day.""")
sparser.add_argument('--month',
                     type=int,
                     choices=list(np.arange(1,13)),
                     default=2,
                     help="""selected month in 2021 (default: %(default)s).""")
# sparser.add_argument('--trange',
#                     metavar=('tstart','tend'),
#                     nargs=2,
#                     help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
sparser.add_argument('--obsid',
                     default='cbw_207',
                     help="""restrict Jacobian to one single observational location (default: %(default)s).""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")



################################################################################
#
#                   p r o g r a m   s t a r t
#
def main(args):

    if args.subcmds=='testbuild_jacobian_1day':
        subcmd_testbuild_jacobian_1day(args)

    if args.subcmds=='testbuild_jacobian_1month':
        subcmd_testbuild_jacobian_1month(args)

if __name__ == '__main__':
    import datetime as dtm

    progname = os.path.basename(__file__)

    #-----------------------------
    #          P R O G R A M   S T A R T
    #
    fmt = "%Y-%m-%dT%H:%M:%S.%f"
    ttstart = dtm.datetime.now()
    logger.info(f"{progname}::PROGRAM START::{ttstart.strftime(fmt)}")
    argv = ' '.join(sys.argv)
    logger.info(f"  command-line -->{argv}<--")

    #
    #          p a r s e   c o m m a n d   l i n e
    #
    args = parser.parse_args()

    #
    #--        s t a r t   e x e c u t i o n
    #
    main(args)



# ################################################################################
# #
# #                   p r o g r a m   s t a r t
# #
# args = parser.parse_args(sys.argv[1:])

