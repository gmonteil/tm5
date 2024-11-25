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
pn.extension('filedropper')
pn.extension('tabulator')


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
        self.dconf['host'] = self.dconf[self.host]
        
        # Configuration widgets
        self.region_selector = pn.widgets.Select(name='Region', options=list(self.dconf.zoom_configuration), value=self.dconf.run.zoom)
        self.meteo_selector = pn.Row("use coarsened meteo:", pn.widgets.Switch(name='Meteo', value=self.dconf.meteo.coarsened, align='center'))
        self.meteo_write = pn.Row("write coarsened meteo:", pn.widgets.Switch(value=self.dconf.meteo.output, align='center'), visible=not self.dconf.meteo.coarsened)

        # widget containers:
        self.emis_file_selector = {}
        self.emis_categories = {}
        self.add_emission_category = {}
        self.emis_rows = {}

        # define some global widgets
        self.editor = pn.widgets.CodeEditor(language='yaml', sizing_mode='stretch_both', visible=False)
        self.editor_button = pn.widgets.Button(name='Show/hide advanced settings editor')

        self.terminal = pn.widgets.Terminal(options={"cursorBlink": True}, height=800, sizing_mode='stretch_width', write_to_console=True)
        self.rcfile = pn.pane.Markdown()
        
        # Main widgets layout:
        self.settings_widgets = pn.Column(
            pn.pane.Markdown("# General TM5 settings:"),
            self.region_selector,
            pn.Row(self.meteo_selector, self.meteo_write),
            pn.pane.Markdown("# Tracers:"),
        )

        # Tracer-specific widgets
        self.tracer_widgets = {}
        for tracer in self.dconf.run.tracers:
            self.tracer_widgets[tracer] = {}

            # initial condition
            self.tracer_widgets[tracer]['iniconc'] = pn.widgets.RadioBoxGroup(value=self.dconf.initial_condition[tracer].type, name='Initial Condition', options=['CAMS', 'Zero'], inline=True, align='center')

            # emissions
            self.tracer_widgets[tracer]['emis'] = self.emis_widgets(tracer)
            
            trwidgets = pn.Column(
                pn.pane.Markdown(f'## {tracer}'),
                pn.Column("### Initial condition:", self.tracer_widgets[tracer]['iniconc']),
                self.tracer_widgets[tracer]['emis']
            )
            self.settings_widgets.append(pn.WidgetBox(trwidgets))
        
        # Tabs, for extra content
        self.tabs = pn.Tabs(('settings', pn.Column(self.settings_widgets, self.editor_button, self.editor)), dynamic=True)

        # Define interactions
        pn.bind(self.update_key, key="run.zoom", value=self.region_selector, watch=True)
        pn.bind(self.update_key, key="meteo.coarsened", value=self.meteo_selector[1], watch=True)
        pn.bind(self.switch_widget_visibility, widget='meteo_write', value=self.meteo_selector[1], watch=True)
        pn.bind(self.update_key, key='meteo.output', value=self.meteo_write[1], watch=True)
        pn.bind(self.switch_widget_visibility, widget='editor', value=self.editor_button, watch=True)
        pn.bind(self.update_editor, value=self.editor_button, watch=True)

        for tracer in self.dconf.run.tracers:
            pn.bind(self.update_key, key=f'initial_condition.{tracer}.type', value=self.tracer_widgets[tracer]['iniconc'], watch=True)

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
        button_save.on_click(self.write_yaml)
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
    
    def write_yaml(self, event: bool = False) -> str:
        """
        Write the content of the text editor to a yaml file (either provided as an argument, or, by default, to a temporary file), and return the path to that file
        """
        if not event: return
        logger.debug(f"Write yaml file to {self.tmpconf}")
        with open(self.tmpconf, 'w') as fid:
            fid.writelines(OmegaConf.to_yaml(self.dconf, resolve=False))
        self.rcfile.param.update(object = read_rcfile(self.tmpconf, self.host))
        if len(self.tabs) == 1:
            self.tabs.append(('TM5 rc-file', self.rcfile))
        return str(self.tmpconf)

    # Widget groups
    def emis_widgets(self, tracer):
        
        def add_cat(event: bool = False, catname: str = None, catparam: dict | None = None):
            """
            Add a row to the emission selector.
            Each row contains a category name, filename and field name widgets
            """
            if not event: return

            # Create the widgets
            self.emis_categories[tracer][catname] = pn.widgets.TextInput(name='category name:', value=catname, align='start')
            self.emis_file_selector[tracer][catname] = pn.widgets.NestedSelect(
                options = files_available,
                layout = dict(type=pn.Row),
                levels=['filename', 'field'],
                value=catparam
            )
            
            # Create the layout
            if tracer not in self.emis_rows:
                self.emis_rows[tracer] = pn.Column()
                
            self.emis_rows[tracer].append(
                pn.Row(
                    self.emis_categories[tracer][catname], self.emis_file_selector[tracer][catname]
                )
            )
            
            # Determine the interactivity
            pn.bind(self.update_key, key=f'tracers.{tracer}.emissions.categories.{catname}.path', value=self.emis_file_selector[tracer][catname]._widgets[0], watch=True)
            pn.bind(self.update_key, key=f'tracers.{tracer}.emissions.categories.{catname}.field', value=self.emis_file_selector[tracer][catname]._widgets[1], watch=True)
            pn.bind(self.print_value, value=self.emis_categories[tracer][catname], watch=True)
            pn.bind(self.rename_key, key=f'tracers.{tracer}.emissions.categories.{catname}', newkey=self.emis_categories[tracer][catname], watch=True, create_value={})

        def calc_stats(event):
            if not event: return ''
            df = xr.Dataset()
            for cat in self.dconf.tracers[tracer].emissions.categories:
                ds = xr.open_dataset(self.dconf.tracers[tracer].emissions.categories[cat].path)[self.dconf.tracers[tracer].emissions.categories[cat].field]
                df[cat] = ds.sum(('lat', 'lon')).resample(time='YS').sum()
            return pn.widgets.Tabulator(df.to_dataframe())

        categories = self.dconf.tracers[tracer].emissions.categories
        files_available = {}
        for cat in categories:
            fname = Path(categories[cat].path)
            fields = list(xr.open_dataset(fname).data_vars)
            files_available[str(fname)] = fields 
        
        self.emis_file_selector[tracer] = {}
        self.emis_categories[tracer] = {}
        self.add_emission_category[tracer] = pn.widgets.Button(name="Add new category")

        # Create widgets:
        for cat in categories :
            default_value = {Path(categories[cat].path).name: categories[cat].field}
            add_cat(True, cat, default_value)

        upload_button = pn.widgets.Button(name='Upload new emission files')
        upload_file = pn.widgets.FileDropper(visible=upload_button)
        show_stats = pn.widgets.Button(name='show annual budget')
        stats = pn.bind(calc_stats, show_stats, watch=True)
            
        # Interactivity
        self.add_emission_category[tracer].on_click(add_cat)
        
        return pn.Column('### Emission categories:', self.emis_rows[tracer], pn.Row(self.add_emission_category[tracer], upload_button, show_stats), upload_file, stats, styles=dict(background='WhiteSmoke'))
        

    # Interactions between the widgets and the "dconf" object
    def update_key(self, key: str, value):
        """
        Change the value of a key. This can also be used to create a new key
        """
        logger.info(f"update value of key {key} to {value}")
        logger.debug(OmegaConf.select(self.dconf, key))
        parts = key.rsplit('.', maxsplit=1)
        if len(parts) == 1 :
            self.dconf[key] = value
        else:
            OmegaConf.select(self.dconf, parts[0])[parts[1]] = value
        logger.debug(OmegaConf.select(self.dconf, key))
        self.n_updates *= -1
        
    def rename_key(self, key: str, newkey: str, create_value : str | dict | None = None):
        """
        Rename an existing key. Since there is no specific way to do it in OmegaConf, we just make a copy of the key under the new name, and delete the previous one.
        This breaks references to the previous key so should be used with caution.
        """
        parts = key.rsplit('.', maxsplit=1)
        logger.debug(OmegaConf.select(self.dconf, key))
        if len(parts) == 1 :
            if key not in self.dconf:
                return self.update_key(newkey, create_value)
            logger.info(f'rename key {key} into {newkey}')
            self.dconf[newkey] = self.dconf.pop(key)
            logger.debug(OmegaConf.select(self.dconf, newkey))
            logger.debug(OmegaConf.select(self.dconf, key))
        else :
            parent, key = parts
            if key not in OmegaConf.select(self.dconf, parent):
                return self.update_key(f'{parent}.{newkey}', create_value)
            logger.info(f'rename key {parent}.{key} into {parent}.{newkey}')
            OmegaConf.select(self.dconf, parent)[newkey] = OmegaConf.select(self.dconf, parent)[key]
            logger.debug(OmegaConf.select(self.dconf, parent)[newkey])
            logger.debug(OmegaConf.select(self.dconf, parent)[key])
        self.n_updates *= -1
        
    def print_value(self, value):
        print(value)
        
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
        if event: 
            cfgfile = self.write_yaml(True)
            self.terminal.subprocess.run('python', 'forward.py', '-b', '--build-only', '-m', self.host, cfgfile)

    def run_tm5(self, event: bool):
        if event: 
            cfgfile = self.write_yaml(True)
            self.terminal.subprocess.run('python', 'forward.py', '-m', self.host, cfgfile)

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