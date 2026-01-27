import panel as pn
import param
from omegaconf import OmegaConf
import requests
from pathlib import Path
import sys
from loguru import logger
from tm5.gui.widgets import RunSettings
from tm5.gui.widgets.stations import StationExplorer, StatisticsViewer
from tm5 import debug


pn.extension()
pn.extension('terminal')
pn.extension('floatpanel')

#
#-- initial attempts to update layout/appearance
#
css = '''
.bk.precomp-widget-box {
  background: #f0f0f0;
  border-radius: 5px;
  border: 1px black solid;
}
'''
pn.extension(raw_css=[css])

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
    def __init__(self, url_tm5: str = 'http://pancake.nebula:5000', **params):
        super().__init__(**params)
        self.settings = RunSettings()
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
            pn.widgets.Select.from_param(self.param.rcfile),
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


class FitIC_UI(pn.viewable.Viewer):
    def __init__(self, config_file: Path | str = 'gui.yml'):
        super().__init__()
        fix_env()
        self.conf = OmegaConf.load(config_file)
        #
        #-- loguru compatible log level
        #
        loglev = self.conf.loglevel if 'loglevel' in self.conf else "CRITICAL"
        if loglev!=None:
            logger.remove()
            #-- loguru levels are in upper-case
            logger.add(sys.stdout, level=loglev.upper())

    @debug.timer
    def __panel__(self):
        return pn.Tabs(
            ("Setup simulation", ExperimentSetupGUI()),
            ("Precomputed results", pn.Tabs(
                ("Fit statistics", StatisticsViewer(self.conf)),
                ('Modelled timeseries', StationExplorer(self.conf)),
                tabs_location='left',
                dynamic=True,
                css_classes=['precomp-widget-box'])
             ),
            dynamic=True
        )
