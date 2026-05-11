#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
import datetime as dtm
from omegaconf import OmegaConf, DictConfig
from pathlib import Path
from collections import OrderedDict
import datetime as dtm
from loguru import logger
from pandas import date_range, DataFrame
from pandas import Timestamp, Timedelta, concat
import xarray as xr
import numpy as np
from numpy import zeros, tile
from netCDF4 import Dataset
from types import SimpleNamespace
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots,colorbar
from cartopy import crs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from tm5.fitic import read_obs_table
from tm5.gridtools import TM5Grids
from tm5.observations import read_obspack_file
from tm5.post.footprint_io import load_adjoint_fwd #-- this was for earlier diagnostics
from tm5.post.footprint_io import tm5rundir_obstable, tm5rundir_iniconc_1obs
from tm5.post.footprint_io import regions1D_info
from tm5.post.footprint_io import tm5rundir_jacobian3D, tm5rundir_emissions2D
from tm5.post.footprint_io import tm5emisdir_load_emissions2D, emisvector_to_global1x1
from tm5.post.plot_util import cnorm_set
from tm5.post.utilities import lonstr,latstr,set_outname


#
#-- restrict to current AVENGERS zoom setup
#
regions_expect = ['glb600x400','eur300x200','gns100x100',]


def subcmd_emis_visu_debug(args):
    filepath_emis = args.emisfile
    regions_select = args.region

    #
    #-- load meta information from prior emissions file
    #
    if not Path(filepath_emis).exists():
        msg = f"prior emissions file ***{filepath_emis}*** not found!"
        raise FileNotFoundError(msg)
    dsemis = Dataset(filepath_emis)
    ng = dsemis.dimensions['ng'].size
    reg1D = dsemis['/region'][:]
    reg1D_uniq = np.unique(reg1D)
    if len(reg1D_uniq)!=len(regions_expect):
        msg = f"...detected unexpected regions ***{reg1D_uniq}***"
        raise RuntimeError(msg)
    elif not set(reg1D_uniq)==set(regions_expect):
        msg = f"...detected unexpected regions ***{reg1D_uniq}***"
        raise RuntimeError(msg)
    msg = f"...emissions vector is for regions -->{regions_expect}<--"
    logger.info(msg)
    region_info = regions1D_info(regions_expect)
    region_table = region_info.table

    #
    #-- extent for plotting
    #
    lonw, lone, lats, latn = args.extent
    if lonw==-180 and lone==180 and lats==-90 and latn==90:
        domain_tag = 'global'
    else:
        domain_tag = f"{lonstr(lonw)}-{lonstr(lone)}x{latstr(lats)}-{latstr(latn)}"
    msg = f"...using domain tag -->{domain_tag}<--"
    logger.info(msg)

    #
    #-- load emissions
    #
    emission_table = {}
    ncemis = dsemis['/emission']
    try:
        emis_units = ncemis.units
    except AttributeError:
        msg = f"expected attribute -->{units}<-- missing for variable 'emissions'"
        raise AttributeError(msg)
    for reg in regions_select:
        grid = region_table[reg].grid
        cnd_reg = reg1D==reg
        assert np.count_nonzero(cnd_reg)==region_table[reg].ng1D
        if ncemis.dimensions==('nmon','ng'):
            emis_data = ncemis[:][:,cnd_reg]
            nt,_ng = emis_data.shape
            emis_data = emis_data.reshape(nt,grid.nlat,grid.nlon)
        elif ncemis.dimensions==('ng',):
            emis_data = ncemis[:][cnd_reg]
            emis_data = emis_data.reshape(grid.nlat,grid.nlon)
        else:
            raise RuntimeError(f"unexpected dimensions {ncemis.dimensions}")
        if 'nmon' in ncemis.dimensions:
            emis_plot = emis_data[0,:]
        else:
            emis_plot = emis_data[:]
        #-- turn into data array which eases plotting
        da_plot = xr.DataArray(
            emis_plot,
            dims=('lat','lon'),
            coords = {'lon': grid.lonc, 'lat': grid.latc }
        )
        #-- restrict to extent
        da_plot = da_plot.sel(lat=slice(lats,latn),lon=slice(lonw,lone))
        emis_tot = da_plot.sum().values
        pltmin = da_plot.min().values
        pltmean = da_plot.mean().values
        pltmax = da_plot.max().values
        pkw = { 'cbmin':args.cbmin, 'cbmax':args.cbmax }
        pkw, cnorm = cnorm_set(pkw, pltmin, pltmax)
        cmap = 'Reds'
        f, ax = subplots(1, 1, figsize=args.figsize, subplot_kw=dict(projection=crs.PlateCarree()))
        img = ax.imshow(da_plot, origin='lower', extent=args.extent, norm=cnorm, cmap=cmap)
        #-- add colorbar
        cbar = colorbar(img)
        cbar.set_label(f"[{emis_units}]")
        #--
        ax.set_extent(args.extent)
        ax.coastlines()
        #-- title
        title = f"emissions@{reg} ({str(filepath_emis)}), total: {emis_tot:.1f}"
        # title += f"min/mean/max = {pltmin}/{pltmean}/{pltmax}"
        ax.set_title(title)
        #-- add gridlines
        gl = ax.gridlines(crs=crs.PlateCarree(),
                          draw_labels=True, linewidth=0.5, color='gray')
        gl.top_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 8, 'color':'gray'}
        gl.ylabel_style = {'size': 8, 'color':'gray'}
        #
        #-- create file
        #
        outname_tokens = ['fitic-emissions', reg, domain_tag,]
        outname = '_'.join(outname_tokens) + '.png'
        outname = set_outname(args, outname)
        plt.tight_layout()
        plt.savefig(str(outname), dpi=args.dpi)
        plt.close()
        logger.info(f"generated ***{outname}***")
    #-- close handles
    dsemis.close()


def subcmd_emis_visu(args):
    filepath_emis = args.emisfile
    #
    #-- read in emissions onto 1x1 global grid
    #
    emis_info = emisvector_to_global1x1(filepath_emis)
    #-- use the combined field
    emis_select = emis_info.emis_glb1x1
    #-- this was for debugging only
    # emis_select = emis_info.table_native['glb600x400']
    # emis_select = emis_info.table_1x1['glb600x400']
    try:
        emis_units = emis_select.attrs['units']
    except AttributeError:
        emis_units = "unknown"
    if args.tosqm and emis_units=='kgCH4/cell':
        grid_glb1x1 = emis_info.grid_glb1x1
        emis_select = emis_select / grid_glb1x1.area
        emis_select.attrs['units'] = 'kgCH4/m2'
    #
    #-- extent for plotting
    #
    lonw, lone, lats, latn = args.extent
    if lonw==-180 and lone==180 and lats==-90 and latn==90:
        domain_tag = 'global'
    else:
        domain_tag = f"{lonstr(lonw)}-{lonstr(lone)}x{latstr(lats)}-{latstr(latn)}"
    msg = f"...using domain tag -->{domain_tag}<--"
    logger.info(msg)

    #
    #-- restrict to extent
    #
    emis_plot = emis_select.sel(lat=slice(lats,latn),lon=slice(lonw,lone))
    emis_tot = emis_plot.sum().values
    pltmin = emis_plot.min().values
    pltmean = emis_plot.mean().values
    pltmax = emis_plot.max().values
    pkw = { 'cbmin':args.cbmin, 'cbmax':args.cbmax }
    pkw, cnorm = cnorm_set(pkw, pltmin, pltmax)
    cmap = 'Reds'
    f, ax = subplots(1, 1, figsize=args.figsize, subplot_kw=dict(projection=crs.PlateCarree()))
    img = ax.imshow(emis_plot, origin='lower', extent=args.extent, norm=cnorm, cmap=cmap)
    #-- add colorbar
    cbar = colorbar(img)
    cbar.set_label(f"[{emis_units}]")
    #--
    ax.set_extent(args.extent)
    ax.coastlines()
    #-- title
    if args.title:
        title = args.title
    elif args.tosqm:
        title = f"emissions ({str(filepath_emis)})"
    else:
        title = f"emissions ({str(filepath_emis)}), total: {emis_tot:.1f}"
    # title += f"min/mean/max = {pltmin}/{pltmean}/{pltmax}"
    ax.set_title(title)
    #-- add gridlines
    gl = ax.gridlines(crs=crs.PlateCarree(),
                      draw_labels=True, linewidth=0.5, color='gray')
    gl.top_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 8, 'color':'gray'}
    gl.ylabel_style = {'size': 8, 'color':'gray'}
    #
    #-- create file
    #
    outname_tokens = ['fitic-emissions', domain_tag,]
    outname = '_'.join(outname_tokens) + '.png'
    outname = set_outname(args, outname)
    plt.tight_layout()
    plt.savefig(str(outname), dpi=args.dpi)
    plt.close()
    logger.info(f"generated ***{outname}***")


    
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
#--       emis_visu
#
sparser = subparsers.add_parser('emis_visu',
                                help="""visualisation of posterior(or prior) emissions.""")
sparser.add_argument('emisfile',
                     help="""NetCDF file with prior or posterior emissions from FIT-IC inversion experiment directory.""")
sparser.add_argument('--extent',
                     type=float,
                     nargs=4,
                     default=(-180,180,-90,90),
                    help="""selected domain for spatial footprint plots (default: %(default)s).""")
sparser.add_argument('--figsize',
                     nargs=2,
                     type=float,
                     default=(16,10),
                     help="""figure size (width,height) in inches (default: %(default)s).""")
sparser.add_argument('--dpi',
                     type=int,
                     default=150,
                     help="""dots-per-inch (default: %(default)s).""")
sparser.add_argument('--tosqm',
                     action='store_true',
                     help="""whether to convert from kgCH4/cell to kgCH4/m2.""")
sparser.add_argument('--title',
                     help="""explicitly set title of plot.""")
sparser.add_argument('--cbmin',
                     type=float,
                     help="""explicit minimum at colorbar.""")
sparser.add_argument('--cbmax',
                     type=float,
                     help="""explicit maximum at colorbar.""")
sparser.add_argument('--outdir',
                     help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                     help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")


#
#--       emis_visu_debug
#
sparser = subparsers.add_parser('emis_visu_debug',
                                help="""visualisation of posterior(or prior) emissions.""")
sparser.add_argument('emisfile',
                     help="""NetCDF file with prior or posterior emissions from FIT-IC inversion experiment directory.""")
sparser.add_argument('--region',
                     nargs='+',
                     default=['glb600x400',],
                     help="""restrict to emissions assigned to selected regions (default: %(default)s).""")
sparser.add_argument('--extent',
                     type=float,
                     nargs=4,
                     default=(-180,180,-90,90),
                    help="""selected domain for spatial footprint plots (default: %(default)s).""")
sparser.add_argument('--figsize',
                     nargs=2,
                     type=float,
                     default=(16,10),
                     help="""figure size (width,height) in inches (default: %(default)s).""")
sparser.add_argument('--dpi',
                     type=int,
                     default=150,
                     help="""dots-per-inch (default: %(default)s).""")
sparser.add_argument('--cbmin',
                     type=float,
                     help="""explicit minimum at colorbar.""")
sparser.add_argument('--cbmax',
                     type=float,
                     help="""explicit maximum at colorbar.""")
sparser.add_argument('--outdir',
                     help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                     help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")



################################################################################
#
#                   p r o g r a m   s t a r t
#
def main(args):

    if args.subcmds=='emis_visu':
        subcmd_emis_visu(args)

    if args.subcmds=='emis_visu_debug':
        subcmd_emis_visu_debug(args)

if __name__ == '__main__':
    import datetime as dtm

    progname = os.path.basename(__file__)

    #-----------------------------
    #          P R O G R A M   S T A R T
    #
    ttstart = Timestamp.utcnow().isoformat()
    logger.info(f"{progname}::PROGRAM START::{ttstart}")
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
