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
import cartopy.feature as cfeature
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
    emisname = args.emisvar
    regions_select = args.region
    boarders_lw = args.__dict__.get('boarders_lw', 1.5)
    if args.extent!=None:
        lonw, lone, lats, latn = args.extent
        domain_tag = f"{lonstr(lonw)}-{lonstr(lone)}x{latstr(lats)}-{latstr(latn)}"
    else:
        lonw, lone, lats, latn = -180, 180, -90, 90
        domain_tag = 'global'
    lat_slice = slice(lats,latn)
    lon_slice = slice(lonw,lone)
    #
    #--
    #
    emis_results = emisvector_to_global1x1(filepath_emis, varname=emisname, tosqm=False)
    emistable_native = emis_results.table_native
    for reg,reg_data in emistable_native.items():
        sel_data = reg_data.loc[dict(lat=lat_slice,lon=lon_slice)]
        # print(f"@{reg} sel_data lon/lat = {sel_data.lon} / {sel_data.lat}")
        emtot = sel_data.sum().values
        units = sel_data.units
        units = units.replace('/cell','')
        print(f"@{reg}, {domain_tag}: emission shape/total = {sel_data.shape} {emtot}[{units}]")

    # #
    # #-- extent for plotting
    # #
    # lonw, lone, lats, latn = args.extent
    # if lonw==-180 and lone==180 and lats==-90 and latn==90:
    #     domain_tag = 'global'
    # else:
    #     domain_tag = f"{lonstr(lonw)}-{lonstr(lone)}x{latstr(lats)}-{latstr(latn)}"
    # msg = f"...using domain tag -->{domain_tag}<--"
    # logger.info(msg)

    #
    #-- load emissions (at their native resolution and extent)
    #
    for reg in regions_select:
        if not reg in emistable_native:
            msg = f"no emissions at native resolution for {reg}"
            logger.error(msg)
            continue
        da_plot = emistable_native[reg]
        #
        #-- restrict to extent
        #
        lon = da_plot.lon.values
        lat = da_plot.lat.values
        dlon = np.unique(np.diff(lon))
        dlat = np.unique(np.diff(lat))
        assert len(dlon)==1 and len(dlat)==1
        dlon = dlon[0]
        dlat = dlat[0]
        lonw = lon.min() - dlon/2
        lone = lon.max() + dlon/2
        lats = lat.min() - dlat/2
        latn = lat.max() + dlat/2
        domain_tag = f"{lonstr(lonw)}-{lonstr(lone)}x{latstr(lats)}-{latstr(latn)}"
        domain_extent = [lonw,lone,lats,latn]
        # da_plot = da_plot.sel(lat=slice(lats,latn),lon=slice(lonw,lone))
        # print(da_plot)
        #--
        f, ax = subplots(1, 1, figsize=args.figsize, subplot_kw=dict(projection=crs.PlateCarree()))
        if args.extent!=None:
            lonw,lone,lats,latn = args.extent
            if list(args.extent)==[-180,180,-90,90]:
                domain_tag = 'global'
            else:
                domain_tag = f"{lonstr(lonw)}-{lonstr(lone)}x{latstr(lats)}-{latstr(latn)}"
            ax.set_extent(args.extent)
            #-- compute totals for only the visible domain
            da_visible = da_plot.sel(lat=slice(lats,latn),lon=slice(lonw,lone))
            emis_tot = da_visible.sum().values
            pltmin = da_visible.min().values
            pltmean = da_visible.mean().values
            pltmax = da_visible.max().values           
        else:
            #-- compute totals for these regional emissions
            da_visible = da_plot
            emis_tot = da_visible.sum().values
            pltmin = da_visible.min().values
            pltmean = da_visible.mean().values
            pltmax = da_visible.max().values
        pkw = { 'cbmin':args.cbmin, 'cbmax':args.cbmax }
        pkw, cnorm = cnorm_set(pkw, pltmin, pltmax)
        cmap = 'Reds'
        img = ax.imshow(da_plot, origin='lower', extent=domain_extent, norm=cnorm, cmap=cmap)
        #-- add colorbar
        cbar = colorbar(img)
        cbar.set_label(f"[{da_plot.units}]")

        #
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, lw=boarders_lw)
        #-- title
        title = f"emissions@{reg} ({filepath_emis.name}), total: {emis_tot:.1f}"
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


def subcmd_emis_visu(args):
    filepath_emis = args.emisfile
    boarders_lw = args.__dict__.get('boarders_lw', 1.5)
    #
    #-- read in emissions onto 1x1 global grid
    #
    emis_info = emisvector_to_global1x1(filepath_emis, tosqm=args.tosqm)
    #
    #-- use the merged emissions field on 1x1 degree
    #
    emis_select = emis_info.emis_glb1x1
    #-- this was for debugging only
    # emis_select = emis_info.table_native['glb600x400']
    # emis_select = emis_info.table_1x1['glb600x400']
    try:
        emis_units = emis_select.attrs['units']
    except AttributeError:
        emis_units = "unknown"
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
    #
    #-- coast-lines, country boarders
    #
    ax.add_feature(cfeature.BORDERS, lw=boarders_lw)
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


def subcmd_inversion_inspect(args):
    fcpost = args.fcpost
    obsfile = args.obsfile
    if not fcpost.exists():
        msg = f"fcpost file not present ***{fcpost.name}***"
        raise RuntimeError(msg)
    if not obsfile.exists():
        msg = f"obsfile not present ***{obsfile.name}***"
        raise RuntimeError(msg)
    dspost = xr.open_dataset(fcpost)
    dsobs  = xr.open_dataset(obsfile)
    stations = dspost['station'].values
    #
    #-- concentrations as numpy arrays
    #
    cprior = dspost['cprior'].values
    cpost  = dspost['cpost'].values
    cobs   = dsobs['obs'].values
    nobsday,nsta = cprior.shape
    msg = f"nobsday={nobsday} nsta={nsta}"
    logger.info(msg)
    #
    #-- over all stations
    #
    bias_prior = (cprior - cobs).ravel()
    bias_post  = (cpost - cobs).ravel()
    rmse_prior = (bias_prior**2).mean()**0.5
    rmse_post  = (bias_post**2).mean()**0.5
    print(f"RMSE prior/post = {rmse_prior}/{rmse_post}")
    for ista in range(nsta):
        sta = stations[ista]
        bias_prior = cprior[:,ista] - cobs[:,ista]
        bias_post = cpost[:,ista] - cobs[:,ista]
        rmse_prior = (bias_prior**2).mean()**0.5
        rmse_post  = (bias_post**2).mean()**0.5
        print(f"@{sta}: RMSE prior/post = {rmse_prior}/{rmse_post}")
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
                     type=Path,
                     help="""NetCDF file with prior or posterior emissions from FIT-IC inversion experiment directory.""")
sparser.add_argument('--emisvar',
                     default='emission',
                     help="""name of variable quantifiying the emissions (default: %(default)s).""")
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


sparser = subparsers.add_parser('inversion_inspect',
                                help="""debugging inversion results.""")
sparser.add_argument('fcpost',
                     type=Path,
                     help="""NetCDF file providing prior and posterior concentrations.""")
sparser.add_argument('obsfile',
                     type=Path,
                     help="""NetCDF file provding observed concentrations (currently those are still integrated in the foj.nc file.""")



################################################################################
#
#                   p r o g r a m   s t a r t
#
def main(args):

    if args.subcmds=='emis_visu':
        subcmd_emis_visu(args)

    if args.subcmds=='emis_visu_debug':
        subcmd_emis_visu_debug(args)

    if args.subcmds=='inversion_inspect':
        subcmd_inversion_inspect(args)

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
