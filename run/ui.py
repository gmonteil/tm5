import param
import panel as pn
from tm5.settings import load_config
import os
from loguru import logger
import tempfile
from functools import partial, lru_cache
from pathlib import Path
import yaml
from netCDF4 import Dataset
from pandas import Timestamp, DatetimeIndex
import xarray as xr
import hvplot.xarray
from omegaconf import OmegaConf, DictConfig


pn.extension()
pn.extension('katex')
pn.extension('codeeditor')
pn.extension('terminal')
pn.extension('floatpanel')


def generate_editor(dconf: DictConfig):
    yield pn.widgets.CodeEditor(language='yaml', sizing_mode='stretch_both', value=OmegaConf.to_yaml(dconf))


def read_rcfile(yaml_file, host) -> str:
    dconf = load_config(yaml_file, host)
    rcfile = Path(dconf.run.paths.output) / 'forward.rc'
    logger.debug(f"try reading a rc-file at {rcfile}")
    if not Path(rcfile).exists():
        logger.debug("No file found. return empty string")
        return ""
    logger.debug(f"Reading file {rcfile}")
    return f"```{yaml.dump(yaml.safe_load(open(rcfile, 'r')))}```"


@lru_cache
def load_stations(filename: str) -> xr.Dataset:
    logger.info(f"Loading stations file {filename}")
    with Dataset(filename, 'r') as ncf:
        dates = DatetimeIndex([Timestamp(*t) for t in ncf['date_midpoints']])
        tracers = [getattr(ncf, f'tracer_{itrac+1:03.0f}') for itrac in range(ncf.dimensions['tracers'].size)]
        ds = xr.Dataset(coords={'time':dates, 'tracer':tracers})
        for station in ncf.groups:
            stat = ncf[station]
            code = stat.abbr.split('/')[1]
            ds[code] = (('tracer', 'time'), stat['mixing_ratio'][:])
            ds[code].attrs['latitude'] = stat.latitude
            ds[code].attrs['longitude'] = stat.longitude
            ds[code].attrs['altitude'] = stat.altitude
            ds[code].attrs['region'] = stat.region
            ds[code].attrs['site_name'] = stat.name
    return ds


class Stage1(param.Parameterized):

    yaml_file = param.FileSelector(path='./forward.yaml')
    n_updates = param.Integer(0)
    tmpconf = param.String("tm5.yml")

    def __init__(self):
        super().__init__()

        self.dconf = OmegaConf.load(self.yaml_file)

        self.host = os.environ['TM5_HOST']
        
        # Configuration widgets
        self.region_selector = pn.widgets.Select(name='Region', options=list(self.dconf.zoom_configuration), value=self.dconf.run.zoom)
        self.meteo_selector = pn.Row("use coarsened meteo:", pn.widgets.Switch(name='Meteo', value=self.dconf.meteo.coarsened, align='center'))
        self.meteo_write = pn.Row("write coarsened meteo:", pn.widgets.Switch(value=self.dconf.meteo.output, align='center'), visible=True)

        # define some global widgets
        self.editor = pn.widgets.CodeEditor(language='yaml', sizing_mode='stretch_both', visible=False)
        self.editor_button = pn.widgets.Button(name='Show/hide advanced settings editor')

        self.terminal = pn.widgets.Terminal(options={"cursorBlink": True}, height=800, sizing_mode='stretch_width')
        self.rcfile = pn.pane.Markdown()

        # Main widgets layout:
        self.settings_widgets = pn.Column(
            self.region_selector,
            pn.Row(self.meteo_selector, self.meteo_write),
        )

        # Tabs, for extra content
        self.tabs = pn.Tabs(('settings', pn.Column(self.settings_widgets, self.editor_button, self.editor)), dynamic=True)

        # Define interactions
        pn.bind(self.update_key, key="run.zoom", value=self.region_selector, watch=True)
        pn.bind(self.update_key, key="meteo.coarsened", value=self.meteo_selector[1], watch=True)
        pn.bind(self.switch_widget_visibility, widget='meteo_write', value=self.meteo_selector[1], watch=True)
        pn.bind(self.update_key, key='meteo.output', value=self.meteo_write[1], watch=True)
        pn.bind(self.switch_widget_visibility, widget='editor', value=self.editor_button, watch=True)
        pn.bind(self.update_editor, value=self.editor_button, watch=True)

    @param.output(('output_path', param.String))
    def output(self):
        return load_config("forward.yaml", 'laptop').run.paths.output

    @param.depends('yaml_file')
    def view(self):

        # Define the widgets 
        # self.editor.value = pn.rx(OmegaConf.to_yaml)(self.dconf)
        button_run = pn.widgets.Button(name='Run TM5', button_type='success')
        button_save = pn.widgets.Button(name='save yaml')
        button_build = pn.widgets.Button(name='Build TM5')
        button_kill = pn.widgets.Button(name='Kill process', button_type='danger')

        # Define the interactions
        button_save.on_click(self.show_rcfile)
        button_build.on_click(self.build_tm5)
        button_run.on_click(self.run_tm5)
        button_kill.on_click(self.kill_tm5)

        return pn.Row(
            self.tabs,
            pn.Column(
                pn.Row(
                    button_save, 
                    button_build, 
                    button_run,
                    button_kill
                ), self.terminal
            )
        )
        #return pn.pane.Markdown(self.yaml_file)

    def panel(self):
        self.param._set_name('')
        return pn.Column(self.param.yaml_file, self.view,)
    
    def show_rcfile(self, event: bool):
        if not event: return
        self.write_yaml()
        logger.info('show_rcfile')

    def write_yaml(self) -> str:
        """
        Write the content of the text editor to a yaml file (either provided as an argument, or, by default, to a temporary file), and return the path to that file
        """
        logger.debug(f"Write yaml file to {self.tmpconf}")
        self.rcfile.param.update(object = read_rcfile(self.tmpconf, self.host))
        with open(self.tmpconf, 'w') as fid:
            fid.writelines(OmegaConf.to_yaml(self.dconf, resolve=False))
        if len(self.tabs) == 1:
            self.tabs.append(('TM5 rc-file', self.rcfile))
        return str(self.tmpconf)

    # Interactions between the widgets and the "dconf" object

    def update_key(self, key: str, value):
        logger.info(f"update value of key {key} to {value}")
        logger.debug(OmegaConf.select(self.dconf, key))
        OmegaConf.update(self.dconf, key, value)
        logger.debug(OmegaConf.select(self.dconf, key))
        self.n_updates *= -1

    def switch_widget_visibility(self, widget: str, value):
        """
        Invert the visibility of a widget
        """
        logger.debug(f"switch_widget_visibility: {widget}")
        getattr(self, widget).visible = not getattr(self, widget).visible

    def switch_widget_disable(self, widget: str, value):
        """
        Invert the "disabled" of a widget
        """
        getattr(self, widget).disabled = not getattr(self, widget).disabled

    def update_editor(self, value):
        logger.debug('update_editor')
        self.editor.param.update(value = OmegaConf.to_yaml(self.dconf))

    def build_tm5(self, event: bool):
        if event: self.terminal.subprocess.run('python', 'forward.py', '-b', '-m', self.host, self.write_yaml())

    def run_tm5(self, event: bool):
        if event: self.terminal.subprocess.run('python', 'forward.py', '-m', self.host, self.write_yaml())

    def kill_tm5(self, event: bool):
        if event: self.terminal.subprocess.kill()

class Stage2(param.Parameterized):

    output_path = param.String()

    @param.depends('output_path')
    def view(self):
        stations = load_stations(Path(self.output_path) / 'stations/stations.nc4')
        return pn.Column(stations.to_array('site').hvplot(x='time', groupby='site'))

    def panel(self):
        return self.view


pipeline = pn.pipeline.Pipeline()
pipeline.add_stage('Settings', Stage1)
pipeline.add_stage('Results', Stage2)


header = pn.Row(pn.Column(pipeline.title, pipeline.error), pipeline.buttons, sizing_mode='stretch_width')
pn.Column(
    header, 
    pipeline.stage
).servable()