#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
from omegaconf import OmegaConf, DictConfig
from pathlib import Path
from loguru import logger
from omegaconf import OmegaConf
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


def load_adjoint_footprint(conf : DictConfig, obs_table : DataFrame) -> DataFrame:
    outpath = Path(conf.run.paths.output)/'adjemis'
    regions = conf.run.regions
    start = Timestamp(conf.run.start)
    end = Timestamp(conf.run.end)
    # Load the observations tabled
    tracers = set(obs_table.tracer)

    footprints = { 'ilat':[],
                   'ilon':[],
                   'itime':[],
                   'value':[],
                   'region':[],
                   'itrac':[]   }
    coords_table = {}
    for itrac, tracer in enumerate(tracers):
        for iday, day in enumerate(date_range(start, end, freq='d')):
            for ireg, region in enumerate(regions):
                fname = outpath / day.strftime(f'adjemis.{region}.%Y%m%d.nc')
                if not fname.exists():
                    continue
                with Dataset(fname, 'r') as ds:
                    ilat = ds['ilat'][:].data
                    ilon = ds['ilon'][:].data
                    itrac = ds['itrac'][:].data
                    values = ds['values'][:].data
                    lonc = ds['lon'][:].data
                    latc = ds['lat'][:].data
                ndat = len(ilat)

                footprints['ilat'].extend(ilat)
                footprints['ilon'].extend(ilon)
                footprints['itime'].extend([iday] * ndat)
                footprints['value'].extend(values)
                footprints['itrac'].extend(itrac)
                footprints['region'].extend([region] * ndat)
                #-- store spatial coordinates for region (only once)
                if not region in coords_table:
                    coords_table[region] = (lonc,latc)
    footprints = DataFrame.from_dict(footprints)

    # print(f"-"*30)
    # for k,v in coords_table.items():
    #     lonc,latc = v
    #     dlon = np.unique(np.diff(lonc))
    #     dlat = np.unique(np.diff(latc))
    #     print(f"{k}, dlon/dlat = {dlon}/{dlat}:")
    #     print(f"   lon_centers: {lonc}")
    #     print(f"   lat_centers: {latc}")
    # print(f"-"*30)
    # Correct the halos:
    r2 = footprints[footprints.region == 'eur300x200']
    # print("r2.ilon=",r2.ilon.values)
    # print("r2.ilat=",r2.ilat.values)
    footprints.loc[footprints.region == 'eur300x200', 'ilon'] = r2.ilon.values - 2
    footprints.loc[footprints.region == 'eur300x200', 'ilat'] = r2.ilat.values - 2
    r3 = footprints[footprints.region == 'gns100x100']
    # print("r3.ilon=",r3.ilon.values)
    # print("r3.ilat=",r3.ilat.values)
    footprints.loc[footprints.region == 'gns100x100', 'ilon'] = r3.ilon.values - 3
    footprints.loc[footprints.region == 'gns100x100', 'ilat'] = r3.ilat.values - 2

    # ------------------
    # Reproject the eur3x2 data on a 1x1 grid
    fpreg = footprints[footprints.region == 'eur300x200']
    grid0 = TM5Grids.from_corners(west=-36+6, east=54-6, south=22+4, north=74-4, dlon=3, dlat=2)
    grid1 = TM5Grids.from_corners(west=-36+6, east=54-6, south=22+4, north=74-4, dlon=1, dlat=1)

    # Duplicate the points
    rlon, rlat = 3, 2
    # print(f"#ilat={len(fpreg.ilat.values)}, #ilon={len(fpreg.ilon.values)}")
    # x = tile(fpreg.value.values, (rlat, rlon, 1))
    # print(f"x.shape={x.shape}")
    # x = x.reshape(-1)
    # print(f"x.shape={x.shape}")
    ilat = tile(fpreg.ilat.values * rlat, (rlat, rlon, 1))
    ilon = tile(fpreg.ilon.values * rlon, (rlat, rlon, 1))
    values = tile(fpreg.value.values, (rlat, rlon, 1)).reshape(-1)
    itime = tile(fpreg.itime.values, (rlat, rlon, 1)).reshape(-1)
    itrac = tile(fpreg.itrac.values, (rlat, rlon, 1)).reshape(-1)

    # Adjust the indices
    for jj in range(rlat):
        ilat[jj, :, :] += jj
    for ii in range(rlon):
        ilon[:, ii, :] += ii

    # Re-create a DataFrame:
    df_eur = DataFrame.from_dict(dict(ilat=ilat.reshape(-1), ilon=ilon.reshape(-1), value=values, itime=itime, itrac=itrac))
    df_eur.loc[:, 'region'] = 'eur300x200'

    # Adjust the indices to those of a glb1x1 grid:
    df_eur.loc[:, 'ilat'] += grid1.south + 90
    df_eur.loc[:, 'ilon'] += grid1.west + 180

    # ------------------
    # Reproject the glb6x4 grid on glo1x1

    fpreg = footprints[footprints.region == 'glb600x400']
    grid0 = TM5Grids.from_corners(west=-180, east=180, south=-90, north=90, dlon=6, dlat=4)
    grid1 = TM5Grids.from_corners(west=-180, east=180, south=-90, north=90, dlon=1, dlat=1)

    rlon, rlat = 6, 4
    ilat = tile(fpreg.ilat.values * rlat, (rlat, rlon, 1))
    ilon = tile(fpreg.ilon.values * rlon, (rlat, rlon, 1))
    values = tile(fpreg.value.values, (rlat, rlon, 1)).reshape(-1)
    itime = tile(fpreg.itime.values, (rlat, rlon, 1)).reshape(-1)
    itrac = tile(fpreg.itrac.values, (rlat, rlon, 1)).reshape(-1)

    for jj in range(rlat):
        ilat[jj, :, :] += jj
    for ii in range(rlon):
        ilon[:, ii, :] += ii

    df_glb = DataFrame.from_dict(dict(ilat=ilat.reshape(-1), ilon=ilon.reshape(-1), value=values, itime=itime, itrac=itrac))
    df_glb.loc[:, 'region'] = 'eur600x400'

    # ---------------------
    # Finally reproject gns1x1 on glo1x1
    df_gns = footprints[footprints.region == 'gns100x100']
    grid0 = TM5Grids.from_corners(west=0+3, east=18-3, south=42+2, north=58-2, dlon=1, dlat=1)
    df_gns.loc[:, 'ilat'] += grid0.south + 90
    df_gns.loc[:, 'ilon'] += grid0.west + 180

    # ----------------------
    # Concatenate all 3 dataframes

    df_tot = concat([df_eur, df_glb, df_gns])
    return df_tot


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
                    choices=['visu_footp','visu_timeseries',],
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
#-- temporal range (last day excluded)
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
# print(obs_table.columns)
# print(obs_table.dtypes)
# sys.exit(0)

#
#--
#
if args.mode=='visu_timeseries':
    fwd = load_adjemis_region_total(conf, obs_table)
    print(fwd.head(n=3))
    msg = f"time-series plot not yet completed!"
    raise NotImplementedError(msg)

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
    df_tot = load_adjoint_footprint(conf, obs_table)
    logger.info(f"...footprints done.")

    # Store all footprints in a single array
    # - array is global with dimensions time/lat/lon/obs
    #
    fp = zeros((nday, 180, 360, nobs))
    fp[df_tot.itime, df_tot.ilat, df_tot.ilon, df_tot.itrac] = df_tot.value
    adjemis_units = "ppb/kgCH4/cell/s"
    msg = f"...global adjoint emissions assembled, " \
        f"min/mean/max = {fp.min()}/{fp.mean()}/{fp.max()} [{adjemis_units}]"
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
        f"min/mean/max = {fp.min()}/{fp.mean()}/{fp.max()} [{footp_units}]"
    logger.info(msg)
    da = xr.DataArray(
        fp,
        dims=('iday','lat','lon','itrac'),
        coords = {
            'iday': np.arange(nday),
            'lat': grid1x1.latc,
            'lon': grid1x1.lonc,
            'itrac': np.arange(nobs)
            }
        )
    #
    glbmin = da.min().values
    glbmean = da.mean().values
    glbmax  = da.max().values
    print(f"glbmin/glbmean/glbmax =  {glbmin:.15f}/{glbmean:.15f}/{glbmax:.15f}")
    
    for itrac,obsid in enumerate(obsids):# in range(len(obsids)):
        obsid = obsids[itrac]
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
        datrac = da.sel(iday=np.arange(nday),lat=slice(lats,latn),lon=slice(lonw,lone),itrac=itrac)
        dmin = datrac.min().values
        dmean = datrac.mean().values
        dmax  = datrac.max().values
        print(f"dmin/dmean/dmax =  {dmin:.15f}/{dmean:.15f}/{dmax:.15f}")
        #
        #-- footprint w.r.t. to contributions from all time(s)
        #
        if args.footp_tmode=='total':
            fpplot = fp[..., itrac].sum(0) #-- add contribution from all (daily) sensitivity
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
                fpplot = fp[iday,..., itrac] #-- global field(!)
                #-- dataset restricted to domain for plot
                daplot = da.sel(iday=iday,lat=slice(lats,latn),lon=slice(lonw,lone),itrac=itrac)
                vmin = fpplot.min()
                vmean = fpplot.mean()
                vmax = fpplot.max()
                pltmin  = daplot.min().values
                pltmean = daplot.mean().values
                pltmax  = daplot.max().values
                # msg = f"...@{day},itrac={itrac}: vmin/vmean/vmax = {vmin}/{vmean}/{vmax}"
                # logger.info(msg)
                msg = f"...@{day},itrac={itrac}: {pltmin}/{pltmean}/{pltmax}"
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
