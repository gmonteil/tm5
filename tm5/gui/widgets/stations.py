#!/usr/bin/env python

import hashlib
from pandas import DataFrame, Timestamp, concat, Grouper, merge
from pathlib import Path
import xarray as xr
from netCDF4 import Dataset
import numpy as np
from numpy import zeros, corrcoef, where
import panel as pn
import param
from omegaconf import DictConfig
from glob import glob
import holoviews as hv
from holoviews import opts
import hvplot.pandas
from functools import partial
from typing import List
from tqdm.contrib.concurrent import process_map
import xxhash
import numpy as np
import multiprocessing.pool as mp
from loguru import logger

from tm5 import debug
from tm5.gui.css import *


# ----- utilities -----


# @debug.timer
# def md5_(fname: str) -> str:
#     """
#     Returns the md5sum of a file (to check if it has changed or not)
#     """
#     hash_md5 = hashlib.md5()
#     with open(fname, "rb") as f:
#         for chunk in iter(lambda: f.read(4096), b""):
#             hash_md5.update(chunk)
#     return hash_md5.hexdigest()


def md5(fname: str, chunk_size: int=1024 * 1024):
    h = xxhash.xxh32()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def utc_to_lst( time_utc : np.ndarray | np.datetime64, longitude : np.ndarray | float ):
    """conversion of UTC time to local solar time (LST),
    longitude coordinate may be provided as
    scalar (equally applied for all time points)
    array  (individual longitude for each time point)
    """
    #-- ensure same length (if both are arrays)
    if hasattr(time_utc, "__len__") and hasattr(longitude, "__len"):
        assert len(time_utc)==len(longitude)
        
    #-- UTC/LST time difference (15 degrees equals 1 hour)
    if hasattr(longitude, "__len__"):
        td = np.asarray(longitude/15. * 3600, dtype='i4') # [s]
        #-- as timedelta64
        td = np.array(td, dtype='timedelta64[s]')
    else:
        td = int(longitude/15. * 3600) # [s]
        td = np.timedelta64(td, 's')
    #--
    time_LST = time_utc + td

    return time_LST


#@debug.timer
def load_experiment(conf, expname, outmode : str = 'dataframe') -> DataFrame | xr.Dataset:
    """
    Load the TM5 stations output as a single DataFrame.
    Since this is time-consuming, this function will try to cache the results into a netCDF file.
    The "conf" argument is an "omegaconf.DictConfig" object, and should contain the following keys:
    - experiments.path : path where the TM5 "stations.
    - experiments.cache :
    - experiments.list:

    """
    stations_file = Path(conf.experiments.path) / conf.experiments.list[expname] / 'stations/stations.nc4'
    fitic_file = Path(conf.experiments.cache) / expname / 'stations.nc'
    checksum = md5(stations_file)

    if fitic_file.exists():
        out = xr.open_dataset(fitic_file)
        if out.attrs['sourcefile_md5sum'] == checksum:
            #
            #--
            #
            # logger.debug(f"...cached {fitic_file} can be used!")

            #-- conversion to DataFrame
            if outmode=='dataframe':
                logger.debug(f"...converting to dataframe")
                df = out.to_dataframe().reset_index()
                logger.debug(f"......conversion done.")
                return df
            elif outmode in ['xarray','dataset',]:#-- keep as xarray.dataset
                return out
            else:
                msg = f"unexpected outmode ==>{outmode}<=="
                raise RuntimeError(msg)

    with Dataset(stations_file) as ds:

        # Dimensions of the output file should be (time, tracer, site)
        times = [Timestamp(*_) for _ in ds['date_midpoints'][:]]
        tracers = [getattr(ds, f'tracer_{itr + 1:03.0f}') for itr in range(ds.dimensions['tracers'].size)]
        sites = [ds[_].getncattr('name') for _ in ds.groups]
        sitelat = [ds[_].getncattr('latitude') for _ in ds.groups]
        sitelon = [ds[_].getncattr('longitude') for _ in ds.groups]
        #-- MVO::this should be height of simulation (above sea-level)
        #        (had checked, e.g. Zugspitze 2669 which is 2666m + 3m)
        sitealt = [ds[_].getncattr('altitude') for _ in ds.groups]
        abbr = [ds[_].getncattr('abbr').split('/')[1] for _ in ds.groups]

        # Create output dataset:
        out = xr.Dataset(coords=dict(time=times, tracer=tracers, site=abbr))
        out['station'] = (('site',), sites)
        out['latitude'] = (('site',), sitelat)
        out['longitude'] = (('site',), sitelon)
        out['altitude']  = (('site',), sitealt)

        mix = zeros((len(sites), len(times), len(tracers)))
        #
        #-- add local solar time (differntiated by station)
        #
        times_LST = np.empty((len(sites), len(times)), dtype='datetime64[s]')
        times_utc = np.array(times).astype('datetime64[s]')
        for isite, site in enumerate(ds.groups):
            mix[isite, :, :] = ds[site]['mixing_ratio'][:].transpose()
            _lon = float(sitelon[isite])
            # logger.debug(f"@isite={isite} (site={site}, {sites[isite]} lon={_lon})...")
            times_LST[isite,:] = utc_to_lst(times_utc, _lon)
            # logger.debug(f"UTC first/last = {times_utc[0]}/{times_utc[-1]}, " \
            #              f"LST first/last = {times_LST[isite,0]}/{times_LST[isite,-1]}")
        
    #
    #-- additional data variables
    #
    out['mixing_ratio'] = (('site', 'time', 'tracer'), mix)
    #-- pass LST with conversion to datetime64[ns] (nano seconds)
    #   to avoid warning raised by xarray
    #   (which is not severe but should better not pop-up in the GUI)
    #
    out['LST']          = (('site', 'time'), np.asarray(times_LST, dtype='datetime64[ns]'))
    out['sampling_height'] = (('site'), [int(_.split('_')[1]) for _ in out.site.values])
    out['sitecode'] = (('site'), [_.split('_')[0] for _ in out.site.values])
    
    # Select the top vertical level for each site (#TODO: we need to find something better ...)
    df = out[['sampling_height', 'sitecode']].to_dataframe()
    max_height = df.groupby('sitecode')['sampling_height'].transform('max').values
    selected = where(df.sampling_height.values == max_height)[0]
    out = out.isel(site=selected)    
    #
    #-- save FIT-IC enhanced station file to cache
    #
    out.attrs['sourcefile'] = str(stations_file)
    out.attrs['sourcefile_md5sum'] = checksum
    fitic_file.parent.mkdir(exist_ok=True, parents=True)
    out.to_netcdf(fitic_file)
    #
    #-- need to return dataframe
    #
    if outmode=='dataframe':
        df = out.to_dataframe().reset_index()
        return df
    else:
        return out


def load_observations_data(fname: Path) -> DataFrame:
    """
    Read one observation file. Return the content as a dataframe, with the following columns:
    - time
    !dropped - value: observation, in the original units (mol/mol)
    - altitude: altitude (above sea level) of the measurement/observation (I think)
    - latitude
    - longitude
    !dropped - elevation: ground altitude (above sea level) of the observation site (I think)
    !dropped - intake_height: sampling height (above ground)
    - code: site code
    !dropped - filename: name of the observation file
    - site_name: full name of the site
    - obs: observation, in ppb

    The function is very ad-hoc for FitIC => CH4-specific, hardcoded for 2021. But this should be easy enough to fix ...
    I tried making it faster using h5py instead of xarray but it was the same
    """
    vars_select = ['time', 'value', 'altitude', 'latitude', 'longitude', 'elevation', 'intake_height']
    #-- MVO::restrict to only variables required in upstream use
    #        - altitude: for switch night-time/afternoon period
    #        - longitude: for (fast) UTC to LST conversion
    vars_select = ['time', 'value', 'altitude', 'latitude', 'longitude',]
    ds = xr.open_dataset(fname, decode_timedelta=False)[vars_select]
    df = ds.sel(obs=ds.time.dt.year == 2021).to_dataframe()
    if len(df) > 0:
        df['code'] = ds.site_code.lower()
        #-- MVO::(long) filename also present in meta dataframe
        # df['filename'] = ds.encoding['source']
        df.loc[:, 'site_name'] = ds.attrs['site_name']
    #-- unit conversion [mol mol-1] to [ppb]
    df.loc[:, 'value'] = df.loc[:,'value'] * 1.e9
    #-- rename
    df = df.rename(columns={'value':'obs'})
    #
    #-- add LOCAL solar time
    #
    df.loc[:,'LST'] = utc_to_lst(df.loc[:,'time'].values, df.loc[:,'longitude'].values)

    return df


def load_all_observations(observation_files: str, nproc : int | None = None) -> DataFrame:
    with mp.Pool(processes=nproc) as pp:
        return concat(pp.map(load_observations_data, list(glob(observation_files))))


def load_observations_metadata(fname: Path) -> DataFrame:
    """
    Complementary function to "load_observations_data": this one loads a bunch of metadata, for each site:
    - site_name
    - site_code
    - country
    - latitude
    - longitude
    - elevation
    - doi
    - filename
    """
    vars_select = ['time', 'value', 'altitude', 'latitude', 'longitude', 'elevation', 'intake_height']
    ds = xr.open_dataset(fname, decode_timedelta=False)[vars_select]
    return DataFrame({
        'site_name': ds.attrs['site_name'],
        'site_code': ds.attrs['site_code'],
        'country': ds.attrs['site_country'],
        'latitude': ds.attrs['site_latitude'],
        'longitude': ds.attrs['site_longitude'],
        'elevation': ds.attrs['site_elevation'],
        'doi': ds.attrs['obspack_identifier_link'],
        'filename': fname
    }, index=[ds.attrs['site_name']])

@debug.timer
def extract_timeperiod_during_day( obs_model : DataFrame,
                                   high_altitude : float = 1000.,
                                   hours_high : List = ['00:00','04:00'],
                                   hours_low  : List = ['12:00','16:00'],
                                   aggregate : bool = False               ) -> DataFrame:
    """Function to extract (and average) observed and simulated
    concentrations for selected times-of-day (local solar time),
    the times-of-day are differentiated per station according to
    their altitude.

    The input DataFrame is expected to provide columns
    'time': measurement times in UTC
    'longitude': zonal location of stations [degrees East]
    'altitude': altitude of stations [m]
    """
    #
    #-- initial testing to differentiate
    #   stations at high and lower altitude with differing hours-of-day
    #   for the comparison
    # MVO-NOTE::too expensive to do the work here every time!!
    # self.model.to_csv(f"~/obs-and-experiments.csv", sep=';', index=False)
    #
    #
    logger.debug(f"start temporal filtering differentiated by station altitude")
    dfx = obs_model.copy()
    #
    #-- MVO::LST should now have been set already
    #        (see load_observations)
    if not 'LST' in dfx.columns:
        utc_values = dfx.loc[:,'time'].values
        lon_values = dfx.loc[:,'longitude'].values
        dfx.loc[:,'LST'] = utc_to_lst(utc_values, lon_values)
    
    #- Step3:
    #  - make LST index
    dfx = dfx.set_index('LST')
    # dfx.index = dfx.loc[:,'LST']
    # if 'LST' in dfx:
    #     dfx = dfx.drop(['LST',], axis=1)
    #- Step4:
    #  - separate high/low altitude stations
    cnd_hi = dfx['altitude']>=high_altitude
    dfx_hi = dfx.loc[cnd_hi, :]
    dfx_lo = dfx.loc[~cnd_hi, :]
    nhi = len(dfx_hi.loc[:,'code'].unique())
    nlo = len(dfx_lo.loc[:,'code'].unique())
    msg = f"found {nhi}/{nlo} high/low altitude stations " \
        f"(high_altitude={high_altitude})"
    logger.debug(msg)
    #- Step5:
    #  - high/low differentiated selected hours
    hrs,hre = hours_high
    dfx_hi = dfx_hi.between_time(hrs,hre, inclusive='both')
    hrs,hre = hours_low
    dfx_lo = dfx_lo.between_time(hrs,hre, inclusive='both')
    #
    #-- MVO-TODO::for upstream use it would be more flexible
    #             not to average (to daily) already here.
    #             For the time-series plots the daily-averages are good,
    #             but for the weekly biases plots or the total comparison
    #             table it might be more consistent to compute those
    #             on the raw (hourly frequency data in the respective hours
    #             of day) rather than from the averages over these hours.
    #- Step6:
    #  - aggregate to daily averages per station (code)
    if aggregate:
        aggr_dict = {}
        for _ in dfx_hi:
            if _ in ['altitude', 'latitude', 'longitude', 'elevation', 'intake_height', 'filename', 'site_name']:
                aggr_dict[_] = 'first'
            elif _=='code':
                pass
            else:
                aggr_dict[_] = 'mean'
        dfx_hi = dfx_hi.groupby(['code',]).resample('D').agg(aggr_dict).reset_index()
        dfx_lo = dfx_lo.groupby(['code',]).resample('D').agg(aggr_dict).reset_index()
        logger.debug(f"...temporal aggregation done.")
    else:
        dfx_lo = dfx_lo.reset_index()
        dfx_hi = dfx_hi.reset_index()
    #- Step7:
    #  - recombine both data frames
    dfout = concat([dfx_lo, dfx_hi])
    logger.debug(f"...temporal filtering done.")
    return dfout

@debug.timer
def interp_model_1sta(observations: DataFrame, model: DataFrame, station: str, experiments: List[str]) -> DataFrame:
    # TODO: bad stuff will happen if we have more than one tracer ...
    obs = observations[observations.site_name == station].copy()
    mod = model[model.station == station]

    # Interpolate along the time dimension
    mod = mod.set_index(['site', 'time', 'tracer', 'station']).to_xarray().squeeze().interp(time=obs.time).to_dataframe()

    # Merge the two arrays:
    for exp in experiments:
        obs.loc[:, exp] = mod[exp].values

        # bias:
        obs.loc[:, f'bias_{exp}'] = obs[exp] - obs.obs

    return obs


# ----- Plots -----


@debug.timer
def plot_stations(model: DataFrame, observations: DataFrame, station: str | None, experiments: List[str] | None, title : str | None = None):
    """
    Returns a plot comparing observations to model data at a specific station.
    - if no station is provided, the function will return "None"
    - if no experiment is provided, the function will just return a plot of the observations
    """

    if station is None:
        return

    plot = observations[observations.site_name == station].hvplot.points(
        x='time', y='obs', color='k', s=1, label='observations', width=1200, grid=True
    )

    # If no experiment has been requested, return the plot with the obs only
    if experiments is None:
        plotcfg = opts.Overlay(title="observed conentrations",
                               ylabel="[ppb]")
        plot.opts(plotcfg)
        return plot

    for exp in experiments:
        plot *= model.loc[model.station == station].hvplot(x='time', y=exp, label=exp)

    if title==None:
        title = "hourly observed and simulated concentrations"
    plotcfg = opts.Overlay(title=title, ylabel="[ppb]")
    plot.opts(plotcfg)
    return plot

@debug.timer
def plot_stations_v3(obs_model: DataFrame, station: str | None, experiments: List[str] | None, title : str | None = None):
    """
    Returns a plot comparing observations to model data at a specific station.
    - if no station is provided, the function will return "None"
    - if no experiment is provided, the function will just return a plot of the observations
    """

    if station is None:
        return

    cnd_station = (obs_model.site_name == station)
    dfplot = obs_model.loc[cnd_station,:]
    if title==None:
        if experiments is None or len(experiments)==0:
            title = 'observed concentrations'
        else:
            title = "simulated and observed concentrations"
    #-- MVO::should not hard-code the concentration units!
    #        (though, for CH4 simulation we now it is in fact ppb)
    units = '[ppb]'
    plot = dfplot.hvplot.points(
        x='time', y='obs', color='k', s=1, label='observations', width=1200, grid=True, ylabel=units, title=title
    )
    # If no experiment has been requested, return the plot with the obs only
    if experiments is None or len(experiments)==0:
        return plot
    #-- otherwise: add simulated concentrations from selected experiments
    for exp in experiments:
        plot *= dfplot.hvplot(x='time', y=exp, label=exp)

    #
    #-- MVO::this did not work, and is now achieved
    #        by providing title and y-label already
    #        when creating the initial plot.
    #
    # if title==None:
    #     title = "hourly observed and simulated concentrations"
    # plotcfg = opts.Overlay(title=title, ylabel="[ppb]")
    # plot.opts(plotcfg)

    return plot

@debug.timer
def plot_weekly_bias(observations: DataFrame, model: DataFrame, station: str | None, experiments, title : str | None = None):
    if station is None or experiments is None:
        return
    wmean = calc_weekly_bias(observations, model, station, experiments)
    if len(wmean)==0:
        return
    if title==None:
        title = "Weekly averaged bias"
    plot = wmean.hvplot(
        x='time', y=[f'bias_{exp}' for exp in experiments], grid=True,
        width=1200, ylabel='[ppb]', title=title 
    )
    return plot

@debug.timer
def plot_weekly_bias_v3(obs_model: DataFrame, station: str | None, experiments, title : str | None = None):
    if station is None or experiments is None:
        return
    wmean = calc_weekly_bias_v3(obs_model, station, experiments)
    if len(wmean)==0:
        return
    if title is None:
        title = "Weekly averaged bias"
    plot = wmean.hvplot(
        x='time', y=[f'bias_{exp}' for exp in experiments], grid=True,
        width=1200, ylabel='[ppb]', title=title 
    )
    return plot

@debug.timer
def plot_site_info(sites: DataFrame, station: str | None):
    site = sites.loc[station]

    text = pn.pane.Markdown(f"""
    ### {site.site_name}

    - latitude: {site.latitude}
    - longitude: {site.longitude}
    - elevation: {site.elevation}
    - DOI: {site.doi}
    """)

    return pn.Column(
        text,
        sites.hvplot.points(
            x='longitude', y='latitude', geo=True, coastline=True, xlim=(-180, 180), ylim=(-90, 90),
            frame_width=300, hover_cols=['site_name']
        ) *
        sites.loc[[station]].hvplot.points(
            x='longitude', y='latitude', geo=True, coastline=True, xlim=(-180, 180), ylim=(-90, 90), frame_width=300, color='r'
        )
    )


@debug.timer
def plot_histogram_of_fit_residuals(
    observations: DataFrame,
    model: DataFrame,
    station: str,
    experiments: List[str]
) -> hv.Overlay:
    if station is None or experiments is None:
        return
    df = interp_model_1sta(observations, model, station, experiments)
    histplots = []
    for exp in experiments:
        histplots.append(df[f'bias_{exp}'].hvplot.hist(bins=100, alpha=.5, line_width=0, label=exp))
    hvplot = hv.Overlay(histplots).opts(
        title="Histogram of hourly biases",
        xlabel="[ppb]", ylabel="count")
    return hvplot

@debug.timer
def plot_histogram_of_fit_residuals_v3(
        obs_model: DataFrame,
        station: str,
        experiments: List[str],
        title : str = None
) -> hv.Overlay:
    if station is None or experiments is None:
        return
    if title is None:
        title = "Histogram of biases"
    cnd_station = (obs_model.site_name == station)
    dfhist = obs_model.loc[cnd_station,:]
    histplots = []
    for exp in experiments:
        histplots.append(dfhist[f'bias_{exp}'].hvplot.hist(bins=100, alpha=.5, line_width=0, label=exp))
    hvplot = hv.Overlay(histplots).opts(
        title=title,
        xlabel="[ppb]", ylabel="count")
    return hvplot


@debug.timer
def plot_table_statistics(observations: DataFrame, model: DataFrame, station: str, experiments: List[str]):
    df = calc_statistics(observations, model, station, experiments)
    return pn.pane.DataFrame(df, formatters=[lambda x: f'{x:.2f}'] * 3, text_align='center')

@debug.timer
def plot_table_statistics_v3(obs_model: DataFrame, station: str, experiments: List[str]):
    df = calc_statistics_v3(obs_model, station, experiments)
    return pn.pane.DataFrame(df, formatters=[lambda x: f'{x:.2f}'] * 3, text_align='center')


@debug.timer
def plot_stat_maps(model: DataFrame, experiment: str, statistics_type: str):
    #
    #-- make sure only statistics based on total time-series is shown
    #   (see routine calc_fit_statistics2)
    #
    df = model[model.Month.isna() & (model.experiment == experiment)]
    title = f"{statistics_type} ({experiment}, based on overall time-series)"
    return df.hvplot.points(
        x='longitude', y='latitude', c=statistics_type, s='count', geo=True, title=title,
        hover_cols=['site_name', 'Bias', 'RMSE', 'Correlation coefficient'],
        coastline=True, xlim=(-180, 180), ylim=(-90, 90), frame_width=1200)


@debug.timer
def plot_stats_table(model: DataFrame, statistics_type: str, experiment_list: List[str], highlighted_experiment: str):
    #
    #-- make sure only statistics based on total time-series is shown
    #   (see routine calc_fit_statistics2)
    #
    df = model[model.Month.isna()]  # & (self.model.experiment == self.experiment)]
    title = f"{statistics_type} (based on overall time-series)"
    fig = df[df.experiment == highlighted_experiment].sort_values(statistics_type).hvplot.scatter(
        x='site_name',
        s='count',
        y=statistics_type,
        title=title,
        rot=90,
        frame_width=1200,
        height=800,
        grid=True,
        label=highlighted_experiment,
        muted=False
    )

    for exp in experiment_list:
        if exp != highlighted_experiment:
            fig *= df[df.experiment == exp].hvplot.scatter(
                x='site_name',
                s='count',
                y=statistics_type,
                label=exp,
                muted_alpha=0,
                muted=True
            )
    plotcfg = opts.Overlay(ylabel="[ppb]")
    return fig.opts(plotcfg)


# ----- Statistics -----

# There are two functions doing more or less the same stuff (calc_statistics and calc_statistics2). The first one is used by "StationExplorer", the second one by "StatisticsViewer".
# The seconf one (calc_fit_statistics2) is much more efficient, and also computes monthly statistics, but needs to be adapted a bit before being usable in "StationExplorer". 
# So keeping both for now ...


@debug.timer
def calc_weekly_bias(observations: DataFrame, model: DataFrame, station: str | None, experiments: List[str] | None):
    df = interp_model_1sta(observations, model, station, experiments)
    if df is None:
        return

    # Calculate weekly mean:
    return df.groupby(Grouper(key='time', freq='7D')).mean(numeric_only=True)

@debug.timer
def calc_weekly_bias_v3(obs_model: DataFrame, station: str | None, experiments: List[str] | None):
    if obs_model is None:
        return
    # logger.debug(f"obs_model, columns -->{list(obs_model.columns)}<--, station={station}")
    df = obs_model.loc[obs_model.site_name==station,:]
    if len(df)==0:
        return
    #
    #-- lazy computation, calc weakly mean for all numeric columns
    #
    dfw = df.groupby(Grouper(key='time', freq='7D')).mean(numeric_only=True)
    # Calculate weekly mean:
    # logger.debug(f"dfw,columns -->{list(dfw.columns)}<-- shape={dfw.shape}")
    return dfw


@debug.timer
def calc_statistics(observations: DataFrame, model: DataFrame, station: str | None, experiments: List[str] | None):
    """
    - RMSE
    - corr
    - bias
    """
    if station is None or experiments is None:
        return

    df = interp_model_1sta(observations, model, station, experiments)
    stats = {
        'Mean bias': [],
        'Correlation coefficient': [],
        'RMSE': [],
        'experiment': []
    }
    for exp in experiments:
        stats['experiment'].append(exp)
        stats['Mean bias'].append((df[exp] - df.obs).mean())
        stats['RMSE'].append(((df[exp] - df.obs) ** 2).mean() ** .5)
        stats['Correlation coefficient'].append(corrcoef(df[exp].values, df.obs.values)[0, 1])

    return DataFrame.from_dict(stats).set_index('experiment')

@debug.timer
def calc_statistics_v3(obs_model: DataFrame, station: str | None, experiments: List[str] | None):
    """
    - RMSE
    - corr
    - bias
    """
    if station is None or experiments is None:
        return
    df = obs_model.loc[obs_model.site_name==station,:]
    stats = {
        'Mean bias': [],
        'Correlation coefficient': [],
        'RMSE': [],
        'experiment': []
    }
    for exp in experiments:
        bias_column = f"bias_{exp}"
        stats['experiment'].append(exp)
        stats['Mean bias'].append(df.loc[:,bias_column].mean())
        stats['RMSE'].append((df.loc[:,bias_column] ** 2).mean() ** .5)
        stats['Correlation coefficient'].append(corrcoef(df[exp].values, df.obs.values)[0, 1])

    return DataFrame.from_dict(stats).set_index('experiment')


#def interp_model(model: DataFrame, obs: DataFrame) -> DataFrame:
#@debug.timer
def interp_model(model: xr.Dataset, obs: DataFrame, column : str = 'model') -> DataFrame:
    obs = obs.sort_values(["site_name", "time"])
    out = []
    # site_list = []
    for site, g in obs.groupby("site_name", sort=False):
        # #-- if model is dataframe...
        # cnd_site = model['station']==site
        # m = model.loc[cnd_site,:]
        m = model.sel(site=model.station == site)
        mmix_interp = np.interp(
            g["time"].values.astype("datetime64[ns]").astype(float),
            m["time"].values.astype("datetime64[ns]").astype(float),
            m["mixing_ratio"].values.squeeze())
        out.append(mmix_interp)
        # cur_sites = np.full(len(mmix_interp), site)
        # if site_list is None:
        #     site_list = cur_sites
        # else:
        #     site_list = np.concatenate((site_list,cur_sites))
        # site_list.append([site,]*len(mmix_interp))
        # logger.info(f"...appending at site {site}")
        # site_list.append(cur_sites)
        # logger.info(f"...appending done.")
    # logger.trace(f"...concatenating site list")
    # site_list = np.concatenate(site_list)
    # logger.trace(f"...concatenation done.")
    #-- insert simulated mixing ratio
    obs.loc[:, column] = np.concatenate(out)
    #--
    # if 'station' in obs:
    #     assert np.all(obs.loc[:,'station'].values==site_list)
    # else:
    #     obs.loc[:,'station'] = site_list
    return obs

#@debug.timer
def calc_fit_statistics2(model: xr.Dataset, obs: DataFrame) -> DataFrame:
    obs = interp_model(model, obs)
    obs.loc[:, 'Bias'] = obs.model - obs.obs
    obs.loc[:, 'RMSE'] = obs.Bias ** 2
    # obs.loc[:, 'chi2'] = obs.rmse / (obs.err ** 2)

    total_statistics = obs.groupby('site_name').mean(numeric_only=True)
    total_statistics.loc[:, 'RMSE'] = total_statistics.RMSE ** .5
    total_statistics.loc[:, 'Correlation coefficient'] = obs.groupby('site_name')[['obs', 'model']].corr().obs.values[1::2]
    totvar_list = ['site_name', 'latitude', 'longitude', 'altitude', 'elevation', 'intake_height', 'obs', 'model', 'Bias', 'RMSE', 'Correlation coefficient']
    totvar_list = ['site_name', 'latitude', 'longitude', 'obs', 'model', 'Bias', 'RMSE', 'Correlation coefficient']
    total_statistics = total_statistics.reset_index()[totvar_list]
    total_statistics.loc[:, 'count'] = obs.groupby('site_name').obs.count().values ** .5

    monthly_data = obs.groupby([obs.time.to_numpy().astype('datetime64[M]'), 'site_name'])
    monthly_statistics = monthly_data.mean(numeric_only=True)
    monthly_statistics.loc[:, 'RMSE'] = monthly_statistics.RMSE ** .5
    monthly_statistics.loc[:, 'Correlation coefficient'] = monthly_data[['obs', 'model']].corr().obs.values[1::2]
    monthly_statistics.index.rename('Month', level=0, inplace=True)
    monthly_statistics.loc[:, 'count'] = obs.groupby([obs.time.to_numpy().astype('datetime64[M]'), 'site_name']).obs.count() ** .5
    monvar_list = ['Month', 'site_name', 'latitude', 'longitude', 'altitude', 'elevation', 'intake_height', 'obs', 'model', 'Bias', 'RMSE', 'Correlation coefficient', 'count']
    monvar_list = ['Month', 'site_name', 'latitude', 'longitude', 'obs', 'model', 'Bias', 'RMSE', 'Correlation coefficient', 'count']
    monthly_statistics = monthly_statistics.reset_index()[monvar_list]

    return concat([monthly_statistics, total_statistics]).rename(columns={'obs': 'Observations', 'model':'Model'})


#@debug.timer
def calc_statistics2(exp: str, obs: DataFrame, conf: DictConfig) -> DataFrame:
    model = load_experiment(conf, exp, outmode='xarray')
    stats = calc_fit_statistics2(model, obs)
    stats.loc[:, 'experiment'] = exp
    return stats

    
def load_experiments(observations: DataFrame, settings: DictConfig, experiments: List[str]) -> DataFrame:
    fn = partial(calc_statistics2, obs=observations, conf=settings)

    nproc = settings.maxnproc if 'maxnproc' in settings else None
    with mp.Pool(processes=nproc) as pp:
        stats = pp.map(fn, experiments)
    return concat(stats)


# ===== Main widgets =====


class StationExplorer(pn.viewable.Viewer):
    station = param.Selector()
    experiments = param.ListSelector()
    data = param.DataFrame()
    #-- MVO:model will/may represent combined dataframe
    #       containing observations *and* model simulations
    #       Thus naming 'model' is bit misleading...
    model = param.DataFrame()
    sites = param.DataFrame()

    @debug.timer
    def __init__(self, settings: DictConfig):
        super().__init__()
        self.settings = settings
        # Loading mode for experiments
        # - v3:              use new experiment handling
        # - any other value (e.g. v0): stick with previous implementation
        # NOTE:
        # - once thoroughly tested we would stick with implementation
        #   of 'v3' and eventually drop this switch!
        # - there is currently one drawback left with 'v3',
        #   which is that the simulated time-series are only shown
        #   for the time-points interpolated to the observation times...
        self.expmode = 'v0'
        self.expmode = 'v3'
        if 'expmode' in self.settings:
            self.expmode = self.settings.expmode
        #
        #-- dedicated comparison of simulated vs observed concentrations
        #   (1) for high altitude stations using averaged concentration
        #       of midnight data (midnight to 4am local solar time)
        #   (2) for lower altitude stations the averaged concentrations
        #       in the afternoon (12am to 4pm) are taken
        #
        self.obs_model_cmp = None
        self.high_altitude = 1000. #
        self.hours_of_day = {'high':['00:00','04:00'],
                             'low':['12:00','16:00'] }
        self.use_hours_of_day = self.settings.get('use_hours_of_day',False)
        # Initialize widgets
        self.param.experiments.objects = list(self.settings.experiments.list)

        # Preload observations
        self.load_observations()


    @debug.timer
    def __panel__(self):
        station_widget = pn.widgets.Select.from_param(self.param.station,
                                                      description="Only one single station can be selected at a time.")
        exp_widget =  pn.widgets.MultiSelect.from_param(self.param.experiments,
                                                        align='end', height=200, height_policy='max', width_policy='max', description="Multiple experiments can be selected at a time (keep <ctrl> key pressed while selecting further experiments).")
        widgets = pn.Column(
            pn.Row(
                station_widget,
                exp_widget,
                self.table_statistics
            ),
            pn.layout.Divider(),
            pn.Column(
                self.plot_timeseries,
                self.plot_weekly_bias
            ),
            pn.layout.Divider(),
            pn.Row(
                self.site_info,
                self.histogram_of_fit_residuals
            ),
            stylesheets=[precomp_stylesheet,], css_classes=['precomp-right']
        )
        return widgets

    def _maxnproc(self):
        nproc = self.settings.maxnproc if 'maxnproc' in self.settings else None
        return nproc

    def _startup_station(self):
        site_start = self.settings.station_first if 'station_first' in self.settings else 'Cabauw'
        return site_start

    def _station_is_high(self):
        if self.data is None:
            return False
        if self.station is None:
            return False
        cnd_station = (self.data.site_name==self.station)
        # logger.info(f"self.data columns ==>{self.data.columns}<==")
        df = self.data.loc[cnd_station,['altitude']]
        return np.all(df.loc[:,'altitude']>=self.high_altitude)

    def _station_hours_of_day(self):
        is_high = self._station_is_high()
        if is_high:
            return self.hours_of_day['high']
        else:
            return self.hours_of_day['low']

    @debug.timer
    def load_observations(self):
        """
        Load the observations in memory (in the "self.station" DataFrame.
        This is normally called just once, during the "__init__". All the observations in the obs folder will be read => The site selection is done by selecting which file(s) go in that folder!
        """
        # TODO: This is slow. We can make it ways faster by pre-computing and storing the concatenated dataframe
        
        with mp.Pool(processes=self._maxnproc()) as pp:
            self.data = concat(pp.map(load_observations_data, glob(self.settings.observations.files)))
            # self.data.to_csv('stationexplorer_obsdata.csv', index=True)
            self.sites = concat(pp.map(load_observations_metadata, glob(self.settings.observations.files)))
        #
        #-- list of stations
        #-- selected station shown on start-up
        #
        site_list = sorted(set(self.data.site_name))
        site_start = self._startup_station()
        try:
            isite0 = site_list.index(site_start)
        except ValueError:
            isite0 = 0 #-- take first in list if selected station is not present
        self.param.station.objects = site_list
        self.station = self.param.station.objects[isite0]

        #
        #-- self.model becomes DataFrame containing
        #   both, observations and model simulations from
        #   (potentially) multiple experiments
        #
        if self.expmode=='v3':
            self.model = self.data.copy()
        #
        #-- build DataFrame with averaged concentrations over
        #   dedicated periods of time-of-day.
        #
        if self.use_hours_of_day:
            hours_high = self.hours_of_day['high']
            hours_low = self.hours_of_day['low']
            self.obs_model_cmp = \
                extract_timeperiod_during_day(self.model,
                                              high_altitude=self.high_altitude,
                                              hours_high=hours_high,
                                              hours_low=hours_low,
                                              aggregate=True)

    @pn.depends('experiments', watch=True)
    def load_experiments(self):
        """
        Load the model timeseries. The data is loaded somewhat lazily:
        - if the model is already loaded, it stays in memory, even if it is de-selected by the user.
        - if not, it is loaded using "load_experiment", which computes the temporal interpolation based on the "stations.nc" files. Further speedups are embedded in "load_experiment" to avoid recomputing that interpolation multiple times.
        """
        for exp in self.experiments:
            if self.expmode=='v3':
                # If no modelled timeseries has been loaded yet:
                if self.model is None:
                    self.model = self.data.copy()
                    m = load_experiment(self.settings, exp, outmode='xarray')
                    self.model = interp_model(m, self.model, column=exp)
                # Load the remaining data:
                else:
                    if exp not in self.model.columns:
                        m = load_experiment(self.settings, exp, outmode='xarray')
                        self.model = interp_model(m, self.model, column=exp)
                #
                #-- pre-calculate observation-model biases
                #
                bias_column = f"bias_{exp}"
                if not bias_column in self.model:
                    _obs = self.model.loc[:,'obs']
                    _sim = self.model.loc[:,exp]
                    self.model.loc[:, bias_column] = _sim - _obs
            else:
                # If no modelled timeseries has been loaded yet:
                if self.model is None:
                    self.model = load_experiment(self.settings, exp).rename(columns={'mixing_ratio': exp})

                # Load the remaining data:
                else:
                    if exp not in self.model.columns:
                        df = load_experiment(self.settings, exp).rename(columns={'mixing_ratio': exp})
                        self.model = merge(self.model, df[[exp, 'site', 'time', 'tracer', 'station']], on=['site', 'time', 'tracer', 'station'])
        #
        #
        #
        # if not self.model is None:
        #     logger.debug(f"self.model, columns -->{list(self.model.columns)}<--")
        #
        #--
        #
        if self.use_hours_of_day:
            hours_high = self.hours_of_day['high']
            hours_low  = self.hours_of_day['low']
            self.obs_model_cmp = \
                extract_timeperiod_during_day(self.model,
                                              high_altitude=self.high_altitude,
                                              hours_high=hours_high,
                                              hours_low=hours_low,
                                              aggregate=True)

    @pn.depends('station', 'experiments')
    def plot_timeseries(self):
        if self.expmode=='v3':
            if self.obs_model_cmp is None:
                title = "hourly observed and simulated concentrations"
                all_plots = plot_stations_v3(self.model, self.station, self.experiments, title=title)
            else:
                hr1,hr2 = self._station_hours_of_day()
                avg_tag = f"(averaged from {hr1} - {hr2} LST)"
                title = f"observed and simulated concentrations {avg_tag}"
                all_plots = plot_stations_v3(self.obs_model_cmp, self.station, self.experiments, title=title)
        else:
            title = "hourly observed and simulated concentrations"
            all_plots = plot_stations(self.model, self.data, self.station, self.experiments, title=title)
        return all_plots

    @pn.depends('station', 'experiments')
    def plot_weekly_bias(self):
        if self.expmode=='v3':
            if self.obs_model_cmp is None:
               title = "Weekly averaged bias"
               dfvisu = plot_weekly_bias_v3(self.model, self.station, self.experiments, title=title)
            else:
                hr1,hr2 = self._station_hours_of_day()
                avg_tag = f"(concentrations taken from {hr1} - {hr2} LST)"
                title = f"Weekly averaged bias {avg_tag}"
                dfvisu = plot_weekly_bias_v3(self.obs_model_cmp, self.station, self.experiments, title=title)
            return dfvisu
        else:
            return plot_weekly_bias(self.data, self.model, self.station, self.experiments)

    @pn.depends('station')
    def site_info(self):
        return plot_site_info(self.sites, self.station)

    @pn.depends('station', 'experiments')
    def histogram_of_fit_residuals(self):
        if self.expmode=='v3':
            if self.obs_model_cmp is None:
                title = "Histogramm of biases (on hourly basis)"
                return plot_histogram_of_fit_residuals_v3(self.model, self.station, self.experiments, title=title)
            else:
                hr1,hr2 = self._station_hours_of_day()
                avg_tag = f"(concentration taken from {hr1} - {hr2} LST)"
                title = f"Histogramm of biases {avg_tag}"
                return plot_histogram_of_fit_residuals_v3(self.obs_model_cmp, self.station, self.experiments, title=title)
        else:
            return plot_histogram_of_fit_residuals(self.data, self.model, self.station, self.experiments)

    @pn.depends('station', 'experiments')
    def table_statistics(self):
        if self.expmode=='v3':
            if self.obs_model_cmp is None:
                return plot_table_statistics_v3(self.model, self.station, self.experiments)
            else:
                return plot_table_statistics_v3(self.obs_model_cmp, self.station, self.experiments)
        else:
            return plot_table_statistics(self.data, self.model, self.station, self.experiments)
        


class StatisticsViewer(pn.viewable.Viewer):

    experiment = param.Selector()
    data = param.DataFrame()
    model = param.DataFrame()
    statistics_type = param.Selector(objects=['RMSE', 'Bias', 'Correlation coefficient', 'Model', 'Observations'])

    @debug.timer
    def __init__(self, settings: DictConfig):
        super().__init__()
        self.settings = settings

        # Initialize widgets
        self.param.experiment.objects = list(self.settings.experiments.list)

        # Pre-load observations
        self.load_observations()

        # Pre-load experiments data
        self.load_experiments()

        # Set default values:
        self.experiment = self.param.experiment.objects[0]

        #-- MVO::copied from StationExplorer,
        #        but not used yet!
        #
        #-- dedicated comparison of simulated vs observed concentrations
        #   (1) for high altitude stations using averaged concentration
        #       of midnight data (midnight to 4am local solar time)
        #   (2) for lower altitude stations the averaged concentrations
        #       in the afternoon (12am to 4pm) are taken
        #
        self.obs_model_cmp = None
        self.high_altitude = 1000. #
        self.hours_of_day = {'high':['00:00','04:00'],
                             'low':['12:00','16:00'] }
        self.use_hours_of_day = self.settings.get('use_hours_of_day',False)


    @debug.timer
    def __panel__(self):
        exp_desc = f"Only one single experiment can be selected at a time."
        exp_widget =  pn.widgets.Select.from_param(self.param.experiment,
                                                   description=exp_desc)
        metric_desc = f""
        metric_widget = pn.widgets.RadioButtonGroup.from_param(self.param.statistics_type,
                                                               description=metric_desc)
        return pn.Column(
            pn.Row(
                exp_widget,
                metric_widget
            ),
            self.plot_stat_maps,
            self.plot_stats_table,
            stylesheets=[precomp_stylesheet,], css_classes=['precomp-right']
        )

    def _maxnproc(self):
        nproc = self.settings.maxnproc if 'maxnproc' in self.settings else None
        return nproc

    @debug.timer
    def load_observations(self):
        self.data = load_all_observations(self.settings.observations.files,
                                          nproc=self._maxnproc())
        #-- for debugging
        # self.data.to_csv('statisticsviewer_obsdata.csv', index=True)
    @debug.timer
    def load_experiments(self):
        self.model = load_experiments(self.data, self.settings, self.param.experiment.objects)
        #-- for debugging
        # self.model.to_csv(f"statistics-viewer_model.csv", sep=';')

    @pn.depends('statistics_type', 'experiment')
    def plot_stat_maps(self):
        return plot_stat_maps(self.model, self.experiment, self.statistics_type)

    @pn.depends('statistics_type', 'experiment')
    def plot_stats_table(self):
        return plot_stats_table(self.model, self.statistics_type, self.param.experiment.objects, self.experiment)
    
