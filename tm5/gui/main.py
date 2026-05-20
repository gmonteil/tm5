import panel as pn
import param
from omegaconf import OmegaConf, DictConfig
import requests
from pathlib import Path
import sys
from loguru import logger
from tm5.gui.widgets import RunSettings
from tm5.gui.widgets.stations import StationExplorer, StatisticsViewer
from tm5.gui.widgets.precomputed import PrecomputedInfo
from tm5.gui.css import *
# from tm5.gui.widgets.emis import EmissionExplorer
from tm5 import debug
import xarray as xr
import hvplot.xarray


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
    conc = param.ClassSelector(class_=xr.Dataset, precedence=-1)

    def __init__(self, gui_settings: DictConfig):
        super().__init__()
        self.button_fwd = pn.widgets.Button.from_param(self.param.run_forward)
        self.button_inv = pn.widgets.Button.from_param(self.param.run_inv)
        self.gui_settings = gui_settings

        # Load the file list
        self.param.experiment.path = self.gui_settings.emissions.glob_pattern
        self.experiment = self.param.experiment.objects[0]

    def __panel__(self):
        return pn.Column(
            pn.pane.Markdown("# Preconfigured experiments"),
            pn.widgets.Select.from_param(self.param.experiment), 
            pn.pane.Markdown("== some description of the selected experiment =="),
            pn.Row(self.button_fwd, self.button_inv),
            self._conc_plot
        )

    @param.depends('run_forward', watch=True)
    def _run_forward(self):
        # Here "emis" should point to the file from the "Experiment" selector
        r = requests.get(f"{self.gui_settings.backend_url}/forward", params={'emis':self.experiment, 'task':'forward'})

        # Retrieve results (here just the concentrations):
        output_path = Path(r.json()['output'])
        fc = xr.open_dataset(output_path / 'fc.nc')
        obs = xr.open_dataset(output_path / 'ftj.nc')
        obs['forward'] = fc['conc']
        print(obs)
        self.conc = obs[['obs', 'forward', 'iniconc']]

    @param.depends('conc')
    def _conc_plot(self):
        if self.conc is None:
            return ''
        return self.conc.hvplot(x='nobsday')


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
        return pn.Tabs(
            self.setup_tab,
            self.preconfigured_tabs,
            self.precomp_tabs,
            dynamic=True
        )

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
