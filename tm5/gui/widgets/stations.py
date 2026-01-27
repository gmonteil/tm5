#!/usr/bin/env python

import hashlib
from pandas import DataFrame, Timestamp, concat, Grouper, merge
from pathlib import Path
import xarray as xr
from netCDF4 import Dataset
from numpy import zeros, corrcoef, where
import panel as pn
import param
from omegaconf import DictConfig
from glob import glob
import holoviews as hv
from typing import List
import hvplot.pandas
from functools import partial
from tqdm.contrib.concurrent import process_map
from tm5 import debug
import xxhash
import numpy as np
import multiprocessing.pool as mp

#-- background color covering (most) elements of
#   - StationExplorer
#   - StatisticsViewer
_bgcolor = '#E4EDED'


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
            if outmode=='dataframe':
                df = out.to_dataframe().reset_index()
                return df
            else:#-- assuming xarray.dataset
                return out

    with Dataset(stations_file) as ds:

        # Dimensions of the output file should be (time, tracer, site)
        times = [Timestamp(*_) for _ in ds['date_midpoints'][:]]
        tracers = [getattr(ds, f'tracer_{itr + 1:03.0f}') for itr in range(ds.dimensions['tracers'].size)]
        sites = [ds[_].getncattr('name') for _ in ds.groups]
        abbr = [ds[_].getncattr('abbr').split('/')[1] for _ in ds.groups]

        # Create output dataset:
        out = xr.Dataset(coords=dict(time=times, tracer=tracers, site=abbr))
        out['station'] = (('site',), sites)

        mix = zeros((len(sites), len(times), len(tracers)))
        for isite, site in enumerate(ds.groups):
            mix[isite, :, :] = ds[site]['mixing_ratio'][:].transpose()
    #
    #-- additional data variables
    #
    out['mixing_ratio'] = (('site', 'time', 'tracer'), mix)
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
    - value: observation, in the original units (mol/mol)
    - altitude: ground altitude (above sea level) of the observation site
    - latitude
    - longitude
    - elevation: altitude (above sea level) of the observations (I think!)
    - intake_height: sampling height (above ground)
    - code: site code
    - filename: name of the observation file
    - site_name: full name of the site
    - obs: observation, in ppb

    The function is very ad-hoc for FitIC => CH4-specific, hardcoded for 2021. But this should be easy enough to fix ...
    I tried making it faster using h5py instead of xarray but it was the same
    """
    ds = xr.open_dataset(fname, decode_timedelta=False)[['time', 'value', 'altitude', 'latitude', 'longitude', 'elevation', 'intake_height']]
    df = ds.sel(obs=ds.time.dt.year == 2021).to_dataframe()
    if len(df) > 0:
        df['code'] = ds.site_code.lower()
        df['filename'] = ds.encoding['source']
        df.loc[:, 'site_name'] = ds.attrs['site_name']
    df.loc[:, 'obs'] = df['value'] * 1.e9
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
    ds = xr.open_dataset(fname, decode_timedelta=False)[['time', 'value', 'altitude', 'latitude', 'longitude', 'elevation', 'intake_height']]
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
def plot_stations(model: DataFrame, observations: DataFrame, station: str | None, experiments: List[str] | None):
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
        return plot

    for exp in experiments:
        plot *= model.loc[model.station == station].hvplot(x='time', y=exp, label=exp)

    return plot


@debug.timer
def plot_weekly_bias(observations: DataFrame, model: DataFrame, station: str | None, experiments):
    if station is None or experiments is None:
        return
    wmean = calc_weekly_bias(observations, model, station, experiments)
    plot = wmean.hvplot(
        x='time', y=[f'bias_{exp}' for exp in experiments], grid=True, width=1200, ylabel='bias (ppb)'
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
    return hv.Overlay(histplots)


@debug.timer
def plot_table_statistics(observations: DataFrame, model: DataFrame, station: str, experiments: List[str]):
    df = calc_statistics(observations, model, station, experiments)
    return pn.pane.DataFrame(df, formatters=[lambda x: f'{x:.2f}'] * 3, text_align='center')


@debug.timer
def plot_stat_maps(model: DataFrame, experiment: str, statistics_type: str):
    df = model[model.Month.isna() & (model.experiment == experiment)]
    return df.hvplot.points(
        x='longitude', y='latitude', c=statistics_type, s='count', geo=True, title=statistics_type,
        hover_cols=['site_name', 'Bias', 'RMSE', 'Correlation coefficient'],
        coastline=True, xlim=(-180, 180), ylim=(-90, 90), frame_width=1200)


@debug.timer
def plot_stats_table(model: DataFrame, statistics_type: str, experiment_list: List[str], highlighted_experiment: str):
    df = model[model.Month.isna()]  # & (self.model.experiment == self.experiment)]
    fig = df[df.experiment == highlighted_experiment].sort_values(statistics_type).hvplot.scatter(
        x='site_name',
        s='count',
        y=statistics_type,
        title=statistics_type,
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
    return fig


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


#@debug.timer
#def interp_model(model: DataFrame, obs: DataFrame) -> DataFrame:
def interp_model(model: xr.Dataset, obs: DataFrame) -> DataFrame:
    obs = obs.sort_values(["site_name", "time"])
    out = []
    for site, g in obs.groupby("site_name", sort=False):
        # #-- if model is dataframe...
        # cnd_site = model['station']==site
        # m = model.loc[cnd_site,:]
        m = model.sel(site=model.station == site)
        out.append(np.interp(
            g["time"].values.astype("datetime64[ns]").astype(float),
            m["time"].values.astype("datetime64[ns]").astype(float),
            m["mixing_ratio"].values.squeeze()
        ))

    obs.loc[:, "model"] = np.concatenate(out)
    return obs


#@debug.timer
##def calc_fit_statistics2(model: DataFrame, obs: DataFrame) -> DataFrame:
def calc_fit_statistics2(model: xr.Dataset, obs: DataFrame) -> DataFrame:
    obs = interp_model(model, obs)
    obs.loc[:, 'Bias'] = obs.model - obs.obs
    obs.loc[:, 'RMSE'] = obs.Bias ** 2
    # obs.loc[:, 'chi2'] = obs.rmse / (obs.err ** 2)

    total_statistics = obs.groupby('site_name').mean(numeric_only=True)
    total_statistics.loc[:, 'RMSE'] = total_statistics.RMSE ** .5
    total_statistics.loc[:, 'Correlation coefficient'] = obs.groupby('site_name')[['obs', 'model']].corr().obs.values[1::2]
    total_statistics = total_statistics.reset_index()[['site_name', 'latitude', 'longitude', 'altitude', 'elevation', 'intake_height', 'obs', 'model', 'Bias', 'RMSE', 'Correlation coefficient']]
    total_statistics.loc[:, 'count'] = obs.groupby('site_name').obs.count().values ** .5

    monthly_data = obs.groupby([obs.time.to_numpy().astype('datetime64[M]'), 'site_name'])
    monthly_statistics = monthly_data.mean(numeric_only=True)
    monthly_statistics.loc[:, 'RMSE'] = monthly_statistics.RMSE ** .5
    monthly_statistics.loc[:, 'Correlation coefficient'] = monthly_data[['obs', 'model']].corr().obs.values[1::2]
    monthly_statistics.index.rename('Month', level=0, inplace=True)
    monthly_statistics.loc[:, 'count'] = obs.groupby([obs.time.to_numpy().astype('datetime64[M]'), 'site_name']).obs.count() ** .5
    monthly_statistics = monthly_statistics.reset_index()[['Month', 'site_name', 'latitude', 'longitude', 'altitude', 'elevation', 'intake_height', 'obs', 'model', 'Bias', 'RMSE', 'Correlation coefficient', 'count']]

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
    model = param.DataFrame()
    sites = param.DataFrame()

    @debug.timer
    def __init__(self, settings: DictConfig):
        super().__init__()
        self.settings = settings
        # Initialize widgets
        self.param.experiments.objects = list(self.settings.experiments.list)

        # Preload observations
        self.load_observations()

    @debug.timer
    def __panel__(self):
        return pn.Column(
            pn.Row(
                pn.widgets.Select.from_param(self.param.station),
                pn.widgets.MultiSelect.from_param(self.param.experiments, align='end', height=200, height_policy='max', width_policy='max'),
                self.table_statistics
            ),
            pn.Column(
                self.plot_timeseries,
                self.plot_weekly_bias
            ),
            pn.Row(
                self.site_info,
                self.histogram_of_fit_residuals
            ),
            styles=dict(background=_bgcolor)
        )

    def _maxnproc(self):
        nproc = self.settings.maxnproc if 'maxnproc' in self.settings else None
        return nproc

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
        self.param.station.objects = sorted(set(self.data.site_name))
        self.station = self.param.station.objects[0]

    @pn.depends('experiments', watch=True)
    def load_experiments(self):
        """
        Load the model timeseries. The data is loaded somewhat lazily:
        - if the model is already loaded, it stays in memory, even if it is de-selected by the user.
        - if not, it is loaded using "load_experiment", which computes the temporal interpolation based on the "stations.nc" files. Further speedups are embedded in "load_experiment" to avoid recomputing that interpolation multiple times.
        """
        for exp in self.experiments:

            # If no modelled timeseries has been loaded yet:
            if self.model is None:
                self.model = load_experiment(self.settings, exp).rename(columns={'mixing_ratio': exp})

            # Load the remaining data:
            else:
                if exp not in self.model.columns:
                    df = load_experiment(self.settings, exp).rename(columns={'mixing_ratio': exp})
                    self.model = merge(self.model, df[[exp, 'site', 'time', 'tracer', 'station']], on=['site', 'time', 'tracer', 'station'])

    @pn.depends('station', 'experiments')
    def plot_timeseries(self):
        return plot_stations(self.model, self.data, self.station, self.experiments)

    @pn.depends('station', 'experiments')
    def plot_weekly_bias(self):
        return plot_weekly_bias(self.data, self.model, self.station, self.experiments)

    @pn.depends('station')
    def site_info(self):
        return plot_site_info(self.sites, self.station)

    @pn.depends('station', 'experiments')
    def histogram_of_fit_residuals(self):
        return plot_histogram_of_fit_residuals(self.data, self.model, self.station, self.experiments)

    @pn.depends('station', 'experiments')
    def table_statistics(self):
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

    @debug.timer
    def __panel__(self):
        return pn.Column(
            pn.Row(
                pn.widgets.Select.from_param(self.param.experiment),
                pn.widgets.RadioButtonGroup.from_param(self.param.statistics_type)
            ),
            self.plot_stat_maps,
            self.plot_stats_table,
            styles=dict(background=_bgcolor)
        )

    def _maxnproc(self):
        nproc = self.settings.maxnproc if 'maxnproc' in self.settings else None
        return nproc

    @debug.timer
    def load_observations(self):
        self.data = load_all_observations(self.settings.observations.files,
                                          nproc=self._maxnproc())
        # self.data.to_csv('statisticsviewer_obsdata.csv', index=True)
    @debug.timer
    def load_experiments(self):
        self.model = load_experiments(self.data, self.settings, self.param.experiment.objects)

    @pn.depends('statistics_type', 'experiment')
    def plot_stat_maps(self):
        return plot_stat_maps(self.model, self.experiment, self.statistics_type)

    @pn.depends('statistics_type', 'experiment')
    def plot_stats_table(self):
        return plot_stats_table(self.model, self.statistics_type, self.param.experiment.objects, self.experiment)
    
