#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
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
from tm5.post.footprint_io import load_adjoint_fitic_footprint
from tm5.gridtools import TM5Grids
from tm5.post.plot_util import cnorm_set
from tm5.post.utilities import lonstr,latstr

##plt.rcParams['text.usetex'] = True

def load_adjemis_region_total( conf : DictConfig, obs_table : DataFrame ) -> DataFrame:

    # Load the observations tabled
    tracers = set(obs_table.tracer)

    # Accumulators:
    date, region, obs, dmix = [], [], [], []

    for day in date_range(conf.run.start, conf.run.end, freq='1d'):
        for reg in conf.run.regions:
            adjfile = day.strftime(f'{conf.run.paths.output}/adjemis/adjemis.{reg}.%Y%m%d.nc')

            # Adjoint files are written only if needed, so it may not exist for a given day/region
            if not Path(adjfile).exists():
                continue

            # I don't know which key should code for that exact path ... for now I guessed ...
            adj = xr.open_dataset(day.strftime(f'{conf.run.paths.output}/adjemis/adjemis.{reg}.%Y%m%d.nc'))

            ilons = adj.ilon.values
            ilats = adj.ilat.values
            sensi = adj['values'].values
            itrac = adj.itrac.values

            # Read the emissions. For now, I'm hardcoding this to CH4 ...
            em = xr.open_dataset(day.strftime(f'{conf.emissions.CH4.prefix}.CH4.{reg}.%Y%m%d.nc'))
            emtot = em.to_array().sum('variable').values

            for itr in set(itrac):
                date.append(day)
                region.append(reg)
                obs.append(adj.tracer.values[itr].strip().decode())

                subset = itrac == itr
                dmix.append((emtot[ilats[subset], ilons[subset]] * sensi[subset]).sum())

    # Store everything in a common dataframe. Each row contains the contribution of one date and region to one specific obs
    fwd = DataFrame.from_dict({'obs': obs, 'region': region, 'date': date, 'mix': dmix})
    return fwd




################################################################################
#
#                   p a r s e r
#
parser = ArgumentParser()
parser.add_argument('config_file',
                    help="""configuration file used for TM5 forward/adjoint run.""")
parser.add_argument('-m', '--host', default=os.environ['TM5_HOST'],
                    help="""selected machine/host from configuration file (default: %(default)s).""")
parser.add_argument('--mode',
                    choices=['write_footp', 'visu_footp','visu_timeseries',],
                    default='visu_footp',
                    help="""requested analysis or visualisation (default: %(default)s).""")
# parser.add_argument('--trange',
#                     metavar=('tstart','tend'),
#                     nargs=2,
#                     help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
# parser.add_argument('--obsfile',
#                     help="""whether to override the observations file used to create the departures for the adjoint run.""")
parser.add_argument('--stationids',
                    nargs='+',
                    help="""restrict plotting to selected observations/stations (expected format similar to, e.g., cbw_207).""")
parser.add_argument('--footp_tmode',
                    choices=['total','daily',],
                    default='total',
                    help="""whether to plot footprints for the observation w.r.t. contributions from all days before *or* w.r.t. to the individual days.""")
parser.add_argument('--figsize',
                    nargs=2,
                    type=float,
                    default=(10,6),
                    help="""figure size, width/height [inches] (default: %(default)s).""")
parser.add_argument('--dpi',
                    type=int,
                    default=300,
                    help="""dots-per-inch (default: %(default)s).""")
parser.add_argument('--cbmin',
                    type=float,
                    help="""explicit minimum at colorbar.""")
parser.add_argument('--cbmax',
                    type=float,
                    help="""explicit maximum at colorbar.""")
parser.add_argument('--extent',
                    type=float,
                    nargs=4,
                    default=(-180,180,-90,90),
                    help="""selected domain for spatial footprint plots (default: %(default)s).""")
parser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
parser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

################################################################################
#
#                   p r o g r a m   s t a r t
#
args = parser.parse_args(sys.argv[1:])

yaml_file = args.config_file
host      = args.host

#=====================================================
# load configuration, set host
#=====================================================
conf = OmegaConf.load(str(yaml_file))
conf.host = conf[host]
tstart = Timestamp(conf.run.start)
tend   = Timestamp(conf.run.end)
#-- temporal range (simulation ends at last day 0:00, this day is (typically) not simulated!)
trange = list(date_range(conf.run.start, conf.run.end, freq='1d'))[:-1]
nday = len(trange)
# Load the observations tabled
msg = f"load observation table from file ***{conf.observations.file}***"
logger.info(msg)
obs_table = read_obs_table(conf.observations.file).set_index('obsid')
obsids = list(obs_table.index)
nobs = len(obsids)
msg = f"...observations read (#obs = {len(obsids)})."
logger.info(msg)
# print(obs_table)
# print("-"*30)
#
#-- add descriptor for observational sampling
#
obs_table['obstime_info'] = "---"
obs_table['sampling_tag'] = "---"
for row in obs_table.iterrows():
    _idx,_cols = row
    _tobs = _cols.time
    _dts  = _cols.time_window_length #-- [s], _tobs-_dts -- _tobs+_dts
    match _cols.sampling_strategy:
        case 2:
            _tinfo = _tobs.strftime('%Y-%m-%dT%H')
            obs_table.loc[_idx,'obstime_info'] = f"{_tinfo}, instantaneous"
            obs_table.loc[_idx,'sampling_tag'] = f"inst-{_tobs.strftime('%Y%m%dT%H')}"
        case 4:
            _t1 = _tobs-Timedelta(seconds=_dts)
            _t2 = _tobs+Timedelta(seconds=_dts)
            _tinfo = f"{_t1.strftime('%Y-%m-%dT%H')}--{_t2.strftime('%Y-%m-%dT%H')}"
            obs_table.loc[_idx,'obstime_info'] = f"{_tinfo}, average"
            obs_table.loc[_idx,'sampling_tag'] = f"avg-{_t1.strftime('%Y%m%dT%H')}--{_t2.strftime('%Y%m%dT%H')}"
        case other:
            msg = f"sampling strategy ***{_cols.sampling_strategy}*** not handled yet."
            logger.error(msg)
#
#--
#
if args.mode=='visu_timeseries':
    fwd = load_adjemis_region_total(conf, obs_table)
    print(fwd.head(n=3))
    msg = f"time-series plot not yet completed!"
    raise NotImplementedError(msg)

elif args.mode=='write_footp':
    #=====================================================
    # read footprints
    #=====================================================
    logger.info(f"start reading footprints...")
    fpresult = load_adjoint_fitic_footprint(conf)
    fp = fpresult.data
    logger.info(f"...footprints done.")
    if not obsids==fpresult.obsids:
        msg = f"...inconsistency between observation identifiers from obstable " \
            f"({obsids}, {conf.observations.file}) and those from adjoint emissions " \
            f"({fpresult.tracer})"
        raise RuntimeError(msg)
    adjemis_units = fp.attrs['units']
    #-- dump some statistics
    msg = f"...global adjoint emissions assembled, " \
        f"min/mean/max = {fp.values.min()}/{fp.values.mean()}/{fp.values.max()} [{adjemis_units}]"
    logger.info(msg)
    #
    #-- turn data array to dataset
    #
    ds = fp.to_dataset(name='footprint')
    #
    #-- add time variable (*this refers to the adj. emissions!*)
    #
    refday = fpresult.days[0]
    days_since = [ (_-refday).days for _ in fpresult.days ]
    da = xr.DataArray(
        days_since,
        dims=('iday',),
        # coords = {
        #     'iday': fp.iday
        #     },
        attrs = {
            'units': refday.strftime(f"days since %Y-%m-%d")
            }
        )
    ds['time'] = da
    #
    #-- add obsids (very likely corresponding to stations)
    #
    da = xr.DataArray(
        fpresult.obsids,
        dims=('iobs',),
        )
    ds['obsids'] = da
    #
    #-- add observational time-points
    #
    obs_secsince = [ (_-refday).total_seconds() for _ in obs_table.time ]
    da = xr.DataArray(
        obs_secsince,
        dims=('iobs',),
        attrs = {
            'units': refday.strftime(f"seconds since %Y-%m-%d")
            }
        )
    ds['obstime'] = da
    #
    #-- global attributes
    #
    ds.attrs['tm5_rundir']  = conf.run.paths.output
    ds.attrs['date_created'] = dtm.datetime.utcnow().isoformat()
    #
    #-- write to file
    #
    if len(obsids)>4:
        obs_tag = f"{len(obsids)}obs"
    else:
        obs_tag = '--'.join( [_.replace('_','-') for _ in fpresult.obsids] )
    time_tag = "{fpresult.days[0].strftime('%Y%m%d')}-{fpresult.days[-1].strftime('%Y%m%d')}"
    outname_tokens = ['tm5-footprint', obs_tag, time_tag,]
    outfile = args.outname if args.outname!=None else '_'.join(outname_tokens)+'.nc'
    if args.outdir!=None:
        if Path(outfile).is_absolute():
            pass
        else:
            outfile = Path(args.outdir) / outfile
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(str(outfile))
    logger.info(f"generated file ***{outfile}***")

elif args.mode=='visu_footp':
    lonw, lone, lats, latn = args.extent
    # lonw, lone = (-15, 35)
    # lats, latn = (33, 73)
    plotlim_lon = (lonw,lone)
    plotlim_lat = (lats,latn)
    if lonw==-180 and lone==180 and lats==-90 and latn==90:
        domain_tag = 'global'
    else:
        domain_tag = f"{lonstr(lonw)}-{lonstr(lone)}x{latstr(lats)}-{latstr(latn)}"
    msg = f"...using domain tag -->{domain_tag}<--"
    logger.info(msg)
    #=====================================================
    # read footprints
    #=====================================================
    logger.info(f"start reading footprints...")
    # df_tot = load_adjoint_footprint(conf, obs_table)
    fpresult = load_adjoint_fitic_footprint(conf)
    fp = fpresult.data
    logger.info(f"...footprints done.")
    if not obsids==fpresult.obsids:
        msg = f"...inconsistency between observation identifiers from obstable " \
            f"({obsids}, {conf.observations.file}) and those from adjoint emissions " \
            f"({fpresult.tracer})"
        raise RuntimeError(msg)
    adjemis_units = fp.attrs['units']
    #-- dump some statistics
    msg = f"...global adjoint emissions assembled, " \
        f"min/mean/max = {fp.values.min()}/{fp.values.mean()}/{fp.values.max()} [{adjemis_units}]"
    logger.info(msg)
    #
    #--
    #
    footp_units   = "ppb/kgCH4/m2/s"
    grid1x1 = TM5Grids.global1x1()
    msg = f"...apply unit conversion {adjemis_units} --> {footp_units}"
    logger.info(msg)
    fp = fp / grid1x1.area[np.newaxis,:,:,np.newaxis]
    msg = f"...global footprint statistics after unit conversion, " \
        f"min/mean/max = {fp.values.min()}/{fp.values.mean()}/{fp.values.max()} [{footp_units}]"
    logger.info(msg)
    glbmin = fp.min().values
    glbmean = fp.mean().values
    glbmax  = fp.max().values
    for iobs,obsid in enumerate(obsids):# in range(len(obsids)):
        obsid = obsids[iobs]
        obs_info = obs_table.loc[obsid,'obstime_info']
        if 'station_name' in obs_table.columns:
            staid = obs_table.loc[obsid,'station_name']
        else:
            staid = obsid
        #-- coordinates of current obs/station
        coords = obs_table.loc[obsid,['lat', 'lon']]##.drop_duplicates()
        if args.stationids!=None and not staid in args.stationids:
            msg = f"...current footprint for -->{staid}<-- not within selected stations."
            logger.info(msg)
            continue
        #
        #-- footprint restricted to plotted domain
        #
        datrac = fp.sel(iday=np.arange(nday),lat=slice(lats,latn),lon=slice(lonw,lone),iobs=iobs)
        dmin = datrac.min().values
        dmean = datrac.mean().values
        dmax  = datrac.max().values
        print(f"dmin/dmean/dmax =  {dmin:.15f}/{dmean:.15f}/{dmax:.15f}")
        #
        #-- footprint w.r.t. to contributions from all time(s)
        #
        if args.footp_tmode=='total':
            fpplot = fp[..., iobs].sum(0) #-- add contribution from all (daily) sensitivity
            f, ax = subplots(1, 1, figsize=args.figsize, subplot_kw=dict(projection=crs.PlateCarree()))
            ax.set_global()
            ax.coastlines()
            img = ax.imshow(fpplot, origin='lower', extent=(-180, 180, -90, 90), cmap='Reds')
            cbar = colorbar(img)
            cbar.set_label(f"[{footp_units}]")
            #
            # Add the obs coordinates
            ax.plot(coords.lon, coords.lat, 'c+', ms=30, alpha=1)
            # title = rf"$\frac{\delta conc}{\delta emis}$, {staid} ({obs_info})"
            # title = r'$\delta conc$'
            title = f"footprint@{staid} ({obs_info})"
            ax.set_title(title)
            # Restrict to a smaller domain:
            ax.set_xlim(*plotlim_lon)
            ax.set_ylim(*plotlim_lat)
            outdir = args.outdir if args.outdir!=None else 'tm5_footprint-plots_test'
            sampling_tag = obs_table.loc[obsid,'sampling_tag']
            location_tag = f'{staid.replace('_','-')}'
            outname_tokens = ['tm5_footprint', location_tag, sampling_tag,]
            outname = '_'.join(outname_tokens) + '.png'
            outname = Path(outdir) / outname
            outname.parent.mkdir(parents=True, exist_ok=True)
            plt.tight_layout()
            plt.savefig(str(outname), dpi=args.dpi)
            plt.close()
            logger.info(f"generated ***{outname}***")
        #
        #-- footprint w.r.t. to contributions from all time(s)
        #
        elif args.footp_tmode=='daily':
            for iday,day in enumerate(trange):
                #-- contribution from current day
                fpplot = fp[iday,..., iobs] #-- global field(!)
                #-- dataset restricted to domain for plot
                daplot = fp.sel(iday=iday,lat=slice(lats,latn),lon=slice(lonw,lone),iobs=iobs)
                vmin = fpplot.min()
                vmean = fpplot.mean()
                vmax = fpplot.max()
                pltmin  = daplot.min().values
                pltmean = daplot.mean().values
                pltmax  = daplot.max().values
                # msg = f"...@{day},iobs={iobs}: vmin/vmean/vmax = {vmin}/{vmean}/{vmax}"
                # logger.info(msg)
                msg = f"...@{day},iobs={iobs}: {pltmin}/{pltmean}/{pltmax}"
                logger.info(msg)
                f, ax = subplots(1, 1, figsize=args.figsize,
                                 subplot_kw=dict(projection=crs.PlateCarree()))
                ax.set_global()
                ax.coastlines()
                # pkw = {'lognorm':True}
                pkw = { 'cbmin':args.cbmin, 'cbmax':args.cbmax }
                pkw, cnorm = cnorm_set(pkw, pltmin, pltmax)
                # cnorm = None #-- don't yet impose user-defined norm
                # img = ax.imshow(fpplot, norm=cnorm,
                #                 origin='lower', extent=(-180, 180, -90, 90), cmap='Reds')
                img = ax.imshow(daplot, norm=cnorm,
                                origin='lower', extent=(lonw, lone, lats, latn), cmap='Reds')
                # Restrict to a smaller domain:
                ax.set_xlim(*plotlim_lon)
                ax.set_ylim(*plotlim_lat)
                #-- add colorbar
                cbar = colorbar(img)
                cbar.set_label(f"[{footp_units}]")
                #
                # Add the obs coordinates
                ax.plot(coords.lon, coords.lat, 'c+', ms=30, alpha=1)
                #-- add title
                title = f"footprint@{staid} w.r.t. emissions on " \
                    f"{day.strftime('%Y-%m-%d')} ({obs_info})"
                title += '\n' + \
                    f"min/mean/max = {pltmin:.4E}/{pltmean:.4E}/{pltmax:.4E}"
                ax.set_title(title)
                gl = ax.gridlines(crs=crs.PlateCarree(),
                                  draw_labels=True, linewidth=0.5, color='gray')
                gl.top_labels = False
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                gl.xlabel_style = {'size': 8, 'color':'gray'}
                gl.ylabel_style = {'size': 8, 'color':'gray'}
                #--
                outdir = args.outdir if args.outdir!=None else 'tm5_footprint-plots_test'
                sampling_tag = obs_table.loc[obsid,'sampling_tag']
                location_tag = f'{staid.replace('_','-')}'
                emis_tag = f"{day.strftime('emis-%Y%m%d')}"
                outname_tokens = ['tm5_footprint', location_tag, emis_tag, sampling_tag, domain_tag]
                outname = '_'.join(outname_tokens) + '.png'
                outname = Path(outdir) / outname
                outname.parent.mkdir(parents=True, exist_ok=True)
                plt.tight_layout()
                plt.savefig(str(outname), dpi=args.dpi)
                plt.close()
                logger.info(f"generated ***{outname}***")
else:
    msg = f"mode {args.mode} not supported."
    raise RuntimeError(msg)
