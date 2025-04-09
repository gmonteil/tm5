import param
from datetime import date
import panel as pn
from panel.viewable import Viewer
from loguru import logger
from omegaconf import OmegaConf
from typing import List
from functools import lru_cache
from pathlib import Path
import xarray as xr
import os


pn.extension()
pn.extension('terminal')
pn.param.ParamMethod.loading_indicator = True


species_implemented = ['CO2', 'CH4']

def get_hostname():
    import socket
    hostname = socket.gethostname()
    return hostname

# print(f"***{get_hostname().find('cosmos')}***")
if get_hostname().find('cosmos')>=0:
    emission_dir = '/lunarc/nobackup/projects/ghg_inv/michael/TM5/input/ch4/emissions_fitic'
else:
    #-- on ICOS jupyter lab
    emission_dir = '/project/fit_ic/data/input/emissions/CH4'

#print(f"***{get_hostname()}***  ==>{emission_dir}<==")

    
class ReactionSettings(pn.viewable.Viewer):
    """
    Block of widgets controlling *one* chemical reaction.
    Parameters:
        - reacname  ==> defined in the parent tracer class
        - shortname ==> defined in the parent tracer class
        - rate0 ==> user-defined? (not yet ...)
        - rate1 ==> user-defined? (not yet ...)
        - active ==> user defined (check-box)
        - domain ==> user defined (vertical domain)
        - field ==> user defined (from a list)
    """
    reacname = param.String(doc='name of the reaction (in the UI)')
    shortname = param.String(doc='name of the reaction (in the model)')
    rate0 = param.Number(doc="Reaction rate", label='rate a')
    rate1 = param.Number(doc='reaction rate', label='rate b')
    active = param.Boolean(default=True, doc='Enable or disable the reaction')
    domain = param.Selector(default='all', objects=['tropo', 'strato', 'all'], doc='Should the reaction be applied to the whole atmosphere or to a specific layer?')
    field = param.Selector(doc='version')

    def __init__(self, versions: List[str] = None, default_field: str = None, **param):
        super().__init__(**param)
        self.param.field.objects = versions
        self.param.field.doc = f'Version of the reactive component field to be used for the reaction {self.reacname}'
        self.field = versions[0]
        if default_field is not None: 
            self.field = default_field

    def __panel__(self):
        return pn.Row(
            pn.widgets.Checkbox.from_param(self.param.active, name=self.param.reacname),
            pn.widgets.Select.from_param(self.param.field) if len(self.param.field.objects) > 1 else None,
            #pn.widgets.FloatInput.from_param(self.param.rate0),
            #pn.widgets.FloatInput.from_param(self.param.rate1),
        )

     
@lru_cache
def get_emis_file_list(path: Path, pattern: str) -> List[Path]:
    return list(Path(path).glob(pattern))
        
        
class EmissionSettings(pn.viewable.Viewer):
    """
    Block of widgets controlling the settings of *one* emission category.
    Parameters :
        - catname ==> category name (user-defined)
        - fieldname ==> name of the field to read from the netCDF file (user chooses from a list)
        - regions ==> list of regions to which the emissions should be applied (user selects from a list)
        - scf ==> scaling factor ==> for now forced to 1
    """
    catname = param.String(doc='name of the emission category (should be unique to that tracer)')
    fieldname = param.Selector(doc='name of the field to be used', label='Field name')
    filename = param.Selector(doc='name of the emission file', label='File prefix')
    regions = param.ListSelector(doc='region(s) where the emissions should be applied')
    scf = param.Number(doc='scaling factor for the emissions', default=1)
    path = param.Path(doc='location of the emission files')
    fileinfo = param.String(doc='ncdump of the file (for now ...)')

    def __init__(self, regnames, **params):
        super().__init__(**params)
        self.param.regions.objects = regnames
        self.regions = regnames
        #self.update_field_choices()
        self.widgets = {
            'field': pn.widgets.Select.from_param(self.param.fieldname),
            'info': pn.pane.Markdown()
        }
        self.update_file_choices()
    
    def __panel__(self):
        return pn.Column(
            pn.Row(
                pn.widgets.TextInput.from_param(self.param.catname),
                pn.Column(
                    pn.widgets.Select.from_param(self.param.filename),
                    self.widgets['field'],
                ),
                pn.widgets.MultiChoice.from_param(self.param.regions)
            ),
            self.widgets['info']
            #pn.widgets.TextInput.from_param(self.param.scf)
        )
        
    @param.depends('filename', 'path', watch=True)
    def update_field_choices(self):
        """
        Update the choices of the "Field" widget. 
        """
        #available_files = list(Path(self.path).glob(f'**/{self.filename}*.nc*'))
        available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
        ds = xr.open_dataset(available_files[0])
        self.param.fieldname.objects = [ _ for _ in ds.data_vars if _!='area' ]
        self.fieldname = self.param.fieldname.objects[0]
        self.widgets['field'].visible = len(self.param.fieldname.objects) > 1

    @param.depends('path', watch=True)
    def update_file_choices(self):
        available_files = get_emis_file_list(self.path, '**/*.nc*')
        self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
        self.filename = self.param.filename.objects[0]
        
    @param.depends('filename', 'fieldname', watch=True)
    def update_field_description(self):
        available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
        #print(available_files)
        ds = xr.open_dataset(available_files[0])
        self.widgets['info'].object = f"""
        {ds.attrs.get('description', '`file description missing`')}
        
        **{self.fieldname}**
        - *long_name*\t: {ds[self.fieldname].long_name}
        - *units*\t: {ds[self.fieldname].units}'
        """
        

class TracerSettings(pn.viewable.Viewer):
    """
    Tracer-specific settings. This class isn't used directly, but is subclassed by species-specific classes (e.g. CH4TracerSettings below, as some options are too species-specific.
    
    Parameters:
    - tracer_name (user-defined)
    - reactions (defined in the derived class)
    - regions (list of strings, passed from the upper level)
    - emission_path: for now hard-coded ...
    
    In addition to the parameters above, the class stores a list of emission categories (EmissionSettings class).
    """
    tracer_name = param.String(default=None, doc='tracer name (should be unique!)')
    add_emissions_category = param.Event(doc='Add new emissions category', label='Add new emissions category')
    reactions = param.ListSelector(default=[], objects=[], doc='Chemical reactions')
    regions = param.List(doc='list of zoom regions') # This is a top-level parameter, but I don't know how to refer to it ...
    emission_path = param.Path(Path(emission_dir))

    def __init__(self, **params):
        super().__init__(**params)
        self.emissions = []
        self.emissions_widgets = pn.Column()
        
    @param.depends('add_emissions_category', watch=True)
    def add_emis(self):
        self.emissions.append(EmissionSettings(
            catname=f'emissions_{len(self.emissions) + 1}', 
            regnames=self.regions,
            path=self.emission_path
        ))
        self.emissions_widgets.append(self.emissions[-1].__panel__())

    @param.depends('regions', watch=True)
    def update_emis_region(self):
        print(self.tracer_name, self.regions)
        for iemis, emis in enumerate(self.emissions):
            emis.regions = self.regions
            self.emissions_widgets[iemis] = emis.__panel__()

    @param.depends('tracer_name', watch=True)
    def __panel__(self):
        components = [
            pn.widgets.TextInput.from_param(self.param.tracer_name),
            pn.widgets.Select.from_param(self.param.initial_condition),
            self.reaction_widgets,
            self.emissions_widgets,
            pn.widgets.Button.from_param(self.param.add_emissions_category)
        ]
        return pn.layout.Card(pn.Column(*components), title=self.param.tracer_name)

    @property
    def reaction_widgets(self):
        if len(self.reactions) > 0:
            return pn.Column(*[r for r in self.reactions])


class CH4TracerSettings(TracerSettings):
    """
    Class containing settings specific to the CH4 tracers. Derived from the TracerSettings class
    
    parameters:
    - initial_condition (user speficied, from a list of choices, hard-coded for now)
    - reactions (user can activate or deactivate a specific reaction, but the list is hard-coded.
    """
    initial_condition = param.Selector(default='zero', objects=['zero', 'CAMS', 'previous run'], doc='initial condition')
    reactions = param.List(
        [
            ReactionSettings(
                reacname='CH4 + OH (troposphere)', 
                shortname='oh', 
                rate0=2.45e-12, 
                rate1=-1755, 
                domain='tropo', 
                versions=['Spivakovsky', 'CAMS']
            ),
            ReactionSettings(
                reacname='CH4 + OH (stratosphere)', 
                shortname='ohstrat', 
                rate0=2.45e-12, 
                rate1=-1755, 
                domain='strato', 
                versions=['Bruehl']
            ),
            ReactionSettings(
                reacname='CH4 + O(1D)', 
                shortname='o1d', 
                rate0=1.5e-10, 
                domain='strato', 
                versions=['Bruehl']
            ),
            ReactionSettings(
                reacname='CH4 + Cl', 
                shortname='cl', 
                rate0=7.3e-12, 
                rate1=-1280, 
                domain='strato', 
                versions=['Bruehl']
            )
        ], doc='Chemical reactions')
    species = param.String('CH4')
        
        
class CO2TracerSettings(TracerSettings):
    initial_condition = param.Selector(default='zero', objects=['zero', 'CarbonTracker', 'previous run'], doc='initial condition')
    species = param.String('CO2')
        
        
class RunSettings(pn.viewable.Viewer):
    """
    Class containing all the TM5 settings (directly or inside sub-objects).
    
    Top-level parameters:
    - start : start date of the simulation (user-selected, within bounds)
    - end : final date of the simulation (user-selected, within bounds)
    - zoom_configuration: selection of region. The list of available configurations (and what it means in terms of TM5 regions) is hard-coded here).
    - run_name : user-defined simulation name (required to set the output folder?)
    - output_types : choice of outputs. Multiple selection possible, available choices hard-coded.
    - levels : number of vertical levels (for now forced to the default tropo34, since no widget is implemented)
    
    The class defines buttons to create tracers, i.e. instances of a TracerSettings class. 
    """
    start = param.Date(default=date(2021, 1, 1), bounds=(date(2021, 1, 1), date(2021, 12, 31)), doc='Start of the simulation')
    end = param.Date(default=date(2022, 1, 1), bounds=(date(2021, 1, 2), date(2022, 1, 1)), doc='Start of the simulation')
    zoom_configuration = param.Selector(objects={
        'avengers-1':['glb6x4', 'eur3x2', 'gns1x1'],
        'glb6x4':    ['glb6x4'], 
        'glb1x1':    ['glb1x1'], 
        'eur1x1':    ['glb6x4', 'eur3x2', 'eur1x1']       
    }, doc='zoom configuration for TM5 simulation')
    run_name = param.String(default=f'{os.environ["USER"]}_{date.today().strftime("%d%b%Y")}', doc='Name of the simulation')
    output_types = param.ListSelector(default=['stations'], objects=['stations', 'mix', 'columns'], doc='choice of outputs (common to all tracers)')
    create_ch4_tracer = param.Event(doc='Add new CH4 tracer', label='New CH4 tracer')
    create_co2_tracer = param.Event(doc='Add new CO2 tracer', label='New CO2 tracer')
    levels = param.Selector(default='tropo34', objects=['tropo25', 'tropo34'], doc='Number of vertical levels')

    def __init__(self, **params):
        super().__init__(**params)
        self.tracers = []
        self.tracers_widgets = pn.Column()
    
    def __panel__(self):
        return pn.Column(
            pn.pane.Markdown('# Run settings'),
            pn.widgets.TextInput.from_param(self.param.run_name),
            pn.widgets.DatePicker.from_param(self.param.start),
            pn.widgets.DatePicker.from_param(self.param.end),
            pn.widgets.Select.from_param(self.param.zoom_configuration),
            pn.widgets.MultiSelect.from_param(self.param.output_types),
            pn.pane.Markdown('# Tracer settings'),
            self.tracers_widgets,
            pn.Row(
                pn.widgets.Button.from_param(self.param.create_co2_tracer, disabled=True),
                pn.widgets.Button.from_param(self.param.create_ch4_tracer),
            )
        )
        
    @param.depends('create_ch4_tracer', watch=True)
    def add_ch4_tracer(self):
        self.add_tracer('CH4')
        
    @param.depends('create_co2_tracer', watch=True)
    def add_co2_tracer(self):
        self.add_tracer('CO2')

    def add_tracer(self, species: str):
        trname = self.get_unique_trname(species)
        trclass = dict(CH4=CH4TracerSettings, CO2=CO2TracerSettings)[species]
        self.tracers.append(trclass(tracer_name=trname, regions=self.zoom_configuration))
        self.tracers_widgets.append(self.tracers[-1].__panel__())

    def get_unique_trname(self, species: str) -> str:
        trname = species
        incr = 1
        while trname in [_.tracer_name for _ in self.tracers]:
            trname = f'{species}_{incr}'
            incr += 1
        return trname

    @param.depends('zoom_configuration', watch=True)
    def update_tracers(self):
        for tr in self.tracers:
            tr.regions = self.zoom_configuration
    
    
class FitIC_UI(pn.viewable.Viewer):
    """
    Top-level container for the GUI. Contains the widgets that aren't related to the settings (i.e. terminal, buttons, etc.).
    """
    
    rcfile = param.FileSelector(path='*.yaml', doc='TM5 settings (yaml file) to be used as a template', label='TM5 config file')
    run_tm5_button = param.Event(doc='run TM5', label='Run TM5')
    build_tm5_button = param.Event(doc='compile TM5', label='Compile TM5')
    
    def __init__(self, **params):
        super().__init__(**params)
        self.settings = RunSettings()
        self.terminal = pn.widgets.Terminal(options={"cursorBlink": True}, height=300, sizing_mode='stretch_width', write_to_console=True)
        self.row = pn.Row()

    def __panel__(self):
        return pn.Column(
            pn.widgets.Select.from_param(self.param.rcfile),
            self.settings,
            pn.layout.Divider(),
            pn.Row(
                pn.widgets.Button.from_param(self.param.run_tm5_button),
                pn.widgets.Button.from_param(self.param.build_tm5_button),
            ),
            self.row,
            self.terminal    
        )
        
    @param.depends('run_tm5_button', watch=True)
    def run_tm5(self):
        self.update_rcfile()

    @param.depends('build_tm5_button', watch=True)
    def build_tm5(self):
        self.update_rcfile()

    def update_rcfile(self):
        """
        The correspondance between yaml keys and UI params is established in this section
        """

        import pprint
        
        conf = OmegaConf.create()

        # Run section
        conf['run'] = {}
        conf.run.start = f'{self.settings.start}'
        conf.run.end = f'{self.settings.end}'
        conf.run.zoom = self.settings.zoom_configuration
        conf.run.tracers = [_.tracer_name for _ in self.settings.tracers]
        conf.run.levels = self.settings.levels
        
        # Output section
        conf.output = {}
        for outp in self.settings.output_types:
            conf.output[outp] = True

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


FitIC_UI().servable()
