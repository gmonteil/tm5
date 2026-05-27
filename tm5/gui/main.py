import panel as pn
import param
from omegaconf import OmegaConf, DictConfig
import requests
from pathlib import Path
import sys
from loguru import logger
import numpy as np
from tm5.gui.widgets import RunSettings
from tm5.gui.widgets.stations import StationExplorer, StatisticsViewer
from tm5.gui.widgets.precomputed import PrecomputedInfo
from tm5.gui.css import *
# from tm5.gui.widgets.emis import EmissionExplorer
from tm5 import debug
import xarray as xr
import hvplot.xarray
from pandas import read_csv, DataFrame
from numpy import corrcoef

pn.extension()
pn.extension('terminal')
pn.extension('floatpanel')
#-- MVO::would potentially be nicer to differentiate (widget) loading
#        indicator more individually...
#   Reference: https://panel.holoviz.org/how_to/param/examples/loading.html
#
pn.extension(loading_spinner='dots', loading_color='#00aa41', template='bootstrap')
pn.extension(loading_spinner='petal', loading_color='black', template='bootstrap')
pn.param.ParamMethod.loading_indicator = True


def fix_env() -> None:
    import sys, os
    from pathlib import Path
    env_base_path = Path(sys.executable).parents[1]
    # os.environ["SSL_CERT_FILE"] = str(env_base_path / 'ssl' / 'cert.pem')
    # os.environ["SSL_CERT_DIR"] = str(env_base_path / 'ssl' / 'certs')
    # os.environ["REQUESTS_CA_BUNDLE"] = str(env_base_path / 'ssl' / 'cert.pem')
    os.environ["PROJ_LIB"] = str(env_base_path / 'share' / 'proj')


@debug.timer
def conc_statistics( conc : xr.Dataset, stations : list[str] ) -> DataFrame:
    """
    """
    # msg = f"stations -->{stations}<--"
    # logger.debug(stations)
    stats = {
        'station': [],
        'Mean bias (prior)': [],
        'Mean bias (post)' : [],
        'RMSE (prior)': [],
        'RMSE (post)': [],
        'Correlation coefficient (prior)': [],
        'Correlation coefficient (post)': [],
        }
    # msg = f"stats initial -->{stats}<-- (==>{stations}<==)"
    # logger.debug(msg)
    dfc = conc.to_dataframe().reset_index()
    #-- prefer station identifier (rather than station index) for output
    dfc.loc[:,'station'] = stations[dfc.loc[:,'nsta'].values]
    for sta in stations:
        # msg = f"now @station={sta}"
        # logger.debug(msg)
        #
        #-- select current station (for all days)
        #
        cnd = dfc['station']==sta
        _df = dfc.loc[cnd,:]
        stats['station'].append(sta)
        _capri = _df.loc[:,'apri']
        _capos = _df.loc[:,'apos']
        _cobs  = _df.loc[:,'obs']
        # msg = f"...@{sta} RMSE ...({_df.shape})"
        # logger.debug(msg)
        #
        #-- prior statistics
        #
        _bias_prior  = _capri - _cobs
        meanbias_prior = _bias_prior.mean()
        rmse_prior = (_bias_prior ** 2).mean() ** .5
        corrcoef_prior = corrcoef(_capri.values,_cobs.values)[0,1]
        # msg = f"@{sta}, rmse prior ==>{rmse_prior}<=="
        # logger.debug(msg)
        #
        #-- posterior statistics
        #
        _bias_post  = _capos - _cobs
        meanbias_post = _bias_post.mean()       
        rmse_post  = (_bias_post **2).mean() ** .5
        corrcoef_post = corrcoef(_capos.values,_cobs.values)[0,1]
        # msg = f"@{sta}, rmse post ==>{rmse_post}<=="
        # logger.debug(msg)
        #
        #-- fill into dictionary
        #
        stats['Mean bias (prior)'].append(meanbias_prior)
        stats['RMSE (prior)'].append(rmse_prior)
        stats['Correlation coefficient (prior)'].append(corrcoef_prior)
        stats['Mean bias (post)'].append(meanbias_post)
        stats['RMSE (post)'].append(rmse_post)
        stats['Correlation coefficient (post)'].append(corrcoef_post)
    # msg = f"...loop terminated, stats -->{stats}<--"
    # logger.info(msg)
    #
    #-- turn into dataframe
    #
    stats = DataFrame.from_dict(stats).set_index('station')
    # stats.to_csv('stats.csv', index=True)

    return stats


class ExperimentSetupGUI(pn.viewable.Viewer):
    """
    Top-level container for the GUI. Contains the widgets that aren't related to the settings (i.e. terminal, buttons, etc.).
    """

    rcfile = param.FileSelector(path='*.yaml', doc='TM5 settings (yaml file) to be used as a template', label='TM5 config file')
    # run_tm5_button = param.Event(doc='run TM5', label='Run TM5')
    # build_tm5_button = param.Event(doc='compile TM5', label='Compile TM5')
    submit_event = param.Event(doc='submit tm5', label='Submit a new run')
    check_status_event = param.Event(doc='Check status', label='Check status')
    jobid = param.Integer(doc='TM5 job ID')

    @debug.timer
    def __init__(self, gui_settings: DictConfig, url_tm5: str = 'http://pancake.nebula:5000', **params):
        super().__init__(**params)
        self.gui_settings = gui_settings
        self.settings = RunSettings(self.gui_settings)
        self.textbox = pn.pane.Alert(visible=False, width=300)
        self.terminal = pn.widgets.Terminal(options={"cursorBlink": True}, height=300, sizing_mode='stretch_width', write_to_console=True, visible=False)
        self.submit_button = pn.widgets.Button.from_param(self.param.submit_event)
        self.check_status_button = pn.widgets.Button.from_param(self.param.check_status_event, visible=False)
        self.job_selector = pn.widgets.IntInput.from_param(self.param.jobid, visible=False, name='', width=60)
        self.alert = pn.pane.Alert(visible=False)

        self.url_tm5 = url_tm5

    def __panel__(self):
        return pn.Column(
            pn.pane.Markdown("# Setup experiment"),
            # pn.widgets.Select.from_param(self.param.rcfile,
            #                              stylesheets=[setup_stylesheet,]),
            pn.widgets.Select.from_param(self.param.rcfile, styles=runsettings_style),
            self.settings,
            pn.layout.Divider(),
            pn.pane.Markdown("# Submit and monitor experiment"),
            pn.Row(
                self.submit_button, self.check_status_button, self.job_selector
            ),
            self.alert,
            pn.Row(
                self.textbox,
                self.terminal
            ),
            stylesheets=[setup_stylesheet,], css_classes=['setup-base']
        )

    @param.depends('submit_event', watch=True)
    def submit(self):

        # Just in case, clear and hide the terminal and "textbox"
        self.terminal.clear()
        self.terminal.visible = False
        self.textbox.object = None
        self.textbox.visible = False

        try:
            url = f'{self.url_tm5}/submit'
            r = requests.get(url, params={'config': 'toto'})

            if not r.ok:
                self.alert.visible = True
                self.alert.alert_type = 'danger'
                self.alert.object = f'Submission failed. Incorrect request to {r.url} ☠️. '
                return

            self.check_status_button.visible = True
            self.job_selector.visible = True

            # Set the button color according to the submit status
            self.check_status_button.button_type = {
                'finished': 'success',
                'queued': 'warning',
                'running': 'primary'
            }[r.json()['status']]

            self.jobid = int(r.json()['jobid'])

            # Add some output message if required:
            self.alert.visible = True
            self.alert.alert_type = 'primary'
            self.alert.object = f'#### Job submitted, with job id {r.json()["jobid"]}\n'

        except requests.exceptions.ConnectionError:
            self.alert.visible = True
            self.alert.alert_type = 'danger'
            self.alert.object = f'Connection to {url} failed ☠️. Server could not be reached'

    @param.depends('jobid', watch=True)
    def update_status_button(self):
        self.alert.visible = False
        self.alert.object = None
        self.alert.alert_type = 'light'
        self.check_status_button.name = f'Check status of job {self.jobid}'

    @param.depends('check_status_event', watch=True)
    def check_status(self):
        r = requests.get(f'{self.url_tm5}/status/{self.jobid}')

        status = r.json()['status']

        self.alert.visible = True
        match status:
            case "queued":
                self.alert.alert_type = 'warning'
                self.alert.object = f"#### Job {self.jobid} is still queuing. "
            case "running":
                self.alert.alert_type = 'info'
                self.alert.object = f'#### Job {self.jobid} is now running. '
            case 'finished':
                self.alert.alert_type = 'success'
                self.alert.object = f'#### Job {self.jobid} has now completed. '

        self.textbox.visible = True
        self.textbox.alert_type = 'light'
        self.textbox.object = '_Job info:_\n'
        text = '```text'+r.json()['info']+'\n'
        self.textbox.object += text

        if status in ['finished', 'running']:
            self.terminal.visible = True
            self.terminal.clear()
            self.terminal.writelines(r.json()['stdout'])
            # with open(r.json()['outfile'], 'r') as fid:
            #     self.terminal.writelines(fid.readlines())

    # @param.depends('run_tm5_button', watch=True)
    # def run_tm5(self):
    #     self.update_rcfile()

    # @param.depends('build_tm5_button', watch=True)
    # def build_tm5(self):
    #     self.update_rcfile()

    def update_rcfile(self):
        """
        The correspondance between yaml keys and UI params is established in this section
        """
        conf = OmegaConf.create()

        # Run section
        conf['run'] = {}
        conf.run.start = f'{self.settings.start}'
        conf.run.end = f'{self.settings.end}'
        conf.run.zoom = self.settings.zoom_configuration
        conf.run.tracers = [_.tracer_name for _ in self.settings.tracers]
        conf.run.levels = self.settings.levels

        # Output section
        # conf.output = {}
        # for outp in self.settings.output_types:
        #     conf.output[outp] = True

        # Tracers:
        conf.tracers = {}
        conf.initial_condition = {}
        conf.emissions = {}

        for tr in self.settings.tracers:
            conf.tracers[tr.tracer_name] = {}
            conf.tracers[tr.tracer_name].species = tr.species

            # Initial condition:
            conf.initial_condition[tr.tracer_name] = {}
            conf.initial_condition[tr.tracer_name].type = tr.initial_condition

            # Chemistry
            conf.tracers[tr.tracer_name].reactions = {}
            for react in tr.reactions:
                conf.tracers[tr.tracer_name].reactions[react.shortname] = {
                    'rate': [react.rate0, react.rate1],
                    'domain': react.domain,
                    'version': react.field,
                }

            # Emissions
            conf.emissions[tr.tracer_name] = {'categories': {}}
            for emcat in tr.emissions:
                conf.emissions[tr.tracer_name].categories[emcat.catname] = {
                    'field': emcat.fieldname,
                    'path': emcat.filename,
                    'regions': emcat.regions
                }

        self.terminal.write(OmegaConf.to_yaml(conf))
        with open(f'{self.settings.run_name}.yaml', 'w') as fid:
            fid.writelines(OmegaConf.to_yaml(conf))


class PreconfExperimentGUI(pn.viewable.Viewer):
    experiment = param.FileSelector(doc='Prior emission dataset')
    run_forward = param.Event(doc='Do a forward run', label='Submit a new forward run')
    run_inv = param.Event(doc='Do an inversion', label='Perform an inversion')

    # Data containers:
    conc       = param.ClassSelector(class_=xr.Dataset, precedence=-1)
    stats4conc = param.ClassSelector(class_=DataFrame, precedence=-1)

    def __init__(self, gui_settings: DictConfig):
        super().__init__()
        self.button_fwd = pn.widgets.Button.from_param(self.param.run_forward)
        self.button_inv = pn.widgets.Button.from_param(self.param.run_inv)
        self.gui_settings = gui_settings

        # Load the file list
        self.param.experiment.path = self.gui_settings.emissions.glob_pattern
        self.experiment = self.param.experiment.objects[0]
        self.stations = None
        # self.stations_widgets = pn.Column()

    def __panel__(self):
        return pn.Column(
            pn.pane.Markdown("# Preconfigured experiments"),
            pn.widgets.Select.from_param(self.param.experiment), 
            pn.pane.Markdown("== some description of the selected experiment =="),
            pn.Row(self.button_fwd, self.button_inv),
            # self.stations_widgets,
            self._conc_plot,
            self._conc_stats_table
        )

    @param.depends('run_forward', watch=True)
    def _run_forward(self):
        # # Here "emis" should point to the file from the "Experiment" selector
        # r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis':self.experiment, 'task':'forward'})

        # # Retrieve results (here just the concentrations):
        # output_path = Path(r.json()['output'])
        #-- modification provided by Zois (2026-05-25)
        #>>MVO:: it must be the full path (otherwise symlink will fail on the backend)
        # emis = Path(self.experiment).name
        emis = str(self.experiment)
        url = f"{self.gui_settings.backend_url}/forward"
        r = requests.get(url, params={'emis': emis, 'task': 'forward'})
        
        if not r.ok:
            logger.error(
                f"forward run failed: backend returned {r.status_code} for "
                f"emis={emis!r} at {url}. Body: {r.text[:500]}"
            )
            return
        try:
            payload = r.json()
        except requests.exceptions.JSONDecodeError:
            logger.error(
                f"forward run failed: backend returned non-JSON for "
                f"emis={emis!r} at {url}. Body: {r.text[:500]}"
            )
            return
        output_path = Path(payload['output'])
      
        logger.debug(f"reading from ouput directory {str(output_path)}")

        #
        #-- processing result/output folder
        #
        #-- fc.nc: simulated concentrations using the ingoing emissions
        #          c = iniconc + ojac*emis
        fc = xr.open_dataset(output_path / 'fc.nc')
        #-- for now extracting the observations from file foj.nc
        #   NOTE: this will certainly change, but currently all
        #         input for the inversion (including obs.) are present in this file
        obs = xr.open_dataset(output_path / 'foj.nc')
        #-- store list of stations
        self.stations = obs.station.values #-- get station identifiers
        # logger.debug(f"stations -->{self.stations}<--")
        obs['forward'] = fc['conc']
        self.conc = obs[['obs', 'forward', 'iniconc']]

    @param.depends('run_inv', watch=True)
    def _run_inv(self):
        # r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis':self.experiment, 'task':'inversion'})

        # # Retrieve results (here just the concentrations):
        # # This file seems to not have the info on which station the files come from ...
        # # The is the only thing needed to make the plots site-specific.
        # output_path = Path(r.json()['output'])
        #
        #-- as modified by Zois (2026-05-25)
        #
        #>>MVO:: it must be the full path (otherwise symlink will fail on the backend)
        # emis = Path(self.experiment).name
        emis = str(self.experiment)
        url = f"{self.gui_settings.backend_url}/forward"
        r = requests.get(url, params={'emis': emis, 'task': 'inversion'})
        
        if not r.ok:
            logger.error(
                f"inversion failed: backend returned {r.status_code} for "
                f"emis={emis!r} at {url}. Body: {r.text[:500]}"
            )
            return
        try:
            payload = r.json()
        except requests.exceptions.JSONDecodeError:
            logger.error(
                f"inversion failed: backend returned non-JSON for "
                f"emis={emis!r} at {url}. Body: {r.text[:500]}"
            )
            return
        
        output_path = Path(payload['output'])

        #
        #-- processing result/output folder
        #
        #-- 20260526: txk had changed code such that prior/posterior
        #             simulated concentrations (including the signal from
        #             the initial concentration) both are in file
        #             fcpost.nc
        fc = xr.open_dataset(output_path / 'fcpost.nc')
        
        obs = xr.open_dataset(output_path / 'foj.nc')
        self.stations = obs.station.values #-- get station identifiers
        obs['apri'] = fc['cprior']
        obs['apos'] = fc['cpost']
        # msg = f"setting self.conc/self.stats4conc"
        # logger.debug(msg)
        self.conc = obs[['obs', 'apri', 'apos']]
        self.stats4conc = conc_statistics(self.conc, self.stations)
        # msg = f"...setting done."
        # logger.debug(msg)

    @param.depends('conc')
    def _conc_plot(self):
        # msg = f"self.conc -->{self.conc}<--"
        # logger.debug(msg)
        if self.conc is None:
            return ''
        dfc = self.conc.to_dataframe().reset_index()
        # logger.debug(f"dfc.columns -->{dfc.columns}<--")
        #-- prefer station identifier (rather than station index) for output
        dfc.loc[:,'station'] = self.stations[dfc.loc[:,'nsta'].values]
        p = dfc.hvplot.points(x='nobsday', y='obs', grid=True, c='k', label='obs', groupby='station')
        if 'forward' in self.conc:  #-- only forward simulation
            p *= dfc.hvplot.line(x='nobsday', y='forward', c='r', label='forward', groupby='station')
        elif 'apos' in self.conc:   #-- result from inversion
            p *= dfc.hvplot.line(x='nobsday', y='apri', c='r', label='prior', groupby='station')
            p *= dfc.hvplot.line(x='nobsday', y='apos', c='c', label='posterior', groupby='station')
        return p

    @param.depends('stats4conc')
    def _conc_stats_table(self):
        # msg = f"self.stats4conc -->{self.stats4conc}<--"
        # logger.debug(msg)
        if self.stats4conc is None:
            return ''
        else:
            df = self.stats4conc
            return pn.pane.DataFrame(df, text_align='center')
            # return pn.pane.DataFrame(df, formatters=[lambda x: f'{x:.2f}'] * 3, text_align='center')

class FitIC_UI(pn.viewable.Viewer):
    def __init__(self, config_file: Path | str = 'gui.yml'):
        super().__init__()
        fix_env()
        if not Path(config_file).exists():
            msg = f"configuration file ***{config_file}*** not found " \
                f"on GUI startup!"
            raise RuntimeError(msg)
        self.conf = OmegaConf.load(config_file)

        #-- loguru compatible log level
        loglev = self.conf.get('loglevel', "CRITICAL")
        logger.remove()
        if loglev.upper() in ['DEBUG','TRACE',]:
            logger.add('fitic.log', level=loglev.upper())
            logger.add(sys.stdout, level=loglev.upper(),
                        diagnose=True, backtrace=True)
        else:
            #-- loguru levels are in upper-case
            logger.add(sys.stdout, level=loglev.upper())

    @debug.timer
    def __panel__(self):
        if self.drop_precomputed:
            return pn.Tabs(
                self.setup_tab,
                self.preconfigured_tabs,
                dynamic=True
            )
        else:
            return pn.Tabs(
                self.setup_tab,
                self.preconfigured_tabs,
                self.precomp_tabs,
                dynamic=True
            )

    @property
    def drop_precomputed(self):
        return self.conf.get('drop_precomputed', False)
    
    @property
    def setup_tab(self):
        return ("Setup simulation", ExperimentSetupGUI(gui_settings=self.conf))

    @property
    def precomp_tabs(self):
        return (
            "Precomputed simulations", pn.Tabs(
                ("Description", PrecomputedInfo(self.conf)),
                ("Fit statistics", StatisticsViewer(self.conf)),
                ('Modelled timeseries', StationExplorer(self.conf)),
                # ('Emissions', EmissionExplorer(self.conf)),
                tabs_location='left',
                dynamic=False,
                # stylesheets=[precomp_stylesheet,], css_classes=['precomp-left'])
            )
        )

    @property
    def preconfigured_tabs(self):
        return (
            "Preconfigured simulations", PreconfExperimentGUI(gui_settings=self.conf)
        )
