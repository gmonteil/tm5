import os
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
import requests

#-- local packages
#   NOTE: preliminary imports currently,
#         will/need to be accessible eventually via the FIT-IC environment.
from ui_util import fitic_inputemisdir

pn.extension()
pn.extension('terminal')
pn.param.ParamMethod.loading_indicator = True


url_tm5 = 'http://pancake.nebula:5000'
species_implemented = ['CO2', 'CH4']

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

    def copy(self):
        return self.__class__(
            versions = self.param.field.objects,
            reacname = self.reacname,
            shortname = self.shortname,
            rate0 = self.rate0,
            rate1 = self.rate1,
            active = self.active,
            domain = self.domain,
            #field = self.field
        )

     
@lru_cache
def get_emis_file_list(path: Path, pattern: str) -> List[Path]:
    return list(Path(path).glob(pattern))
        

class FieldSelector(pn.viewable.Viewer):
    catname = param.String(doc='category name')
    filename = param.Selector(doc="name of the emission file")
    fieldname = param.Selector(doc="name of the field to be used")
    path = param.Path(doc='location of the emission files')
    desc = param.String(doc="domain of the emissions")
    domain = param.String(doc="title of the section")

    def __init__(self, **params):
        super().__init__(**params)
        self.widgets = dict(
            field = pn.widgets.Select.from_param(self.param.fieldname),
            info = pn.pane.Markdown(width=300),
            title = pn.pane.Markdown(width=300)
        )
        self.update_desc()

    def __panel__(self):
        return pn.Column(
            self.widgets['title'],
            pn.widgets.Select.from_param(self.param.filename),
            self.widgets['field'],
            self.widgets['info']
        )

    @param.depends('filename', 'path', 'domain', watch=True)
    def update_field_choices(self):
        """
        Update the choices of the "Field" widget.
        """
        available_files = get_emis_file_list(Path(self.path) / self.domain, f'{self.filename}*.nc')
        if len(available_files) > 0:
            ds = xr.open_dataset(available_files[0])
            self.param.fieldname.objects = [ _ for _ in ds.data_vars if _!='area' ]
            self.fieldname = self.param.fieldname.objects[0]
            self.widgets['field'].visible = len(self.param.fieldname.objects) > 1

    @param.depends('path', 'domain', watch=True)
    def update_file_choices(self):
        # available_files = get_emis_file_list(self.path, '**/*.nc*')
        #-- 2025-04-14:: restrict here to the global (default) domain
        logger.info(self.domain)
        logger.info(f'**/*{self.domain}*.nc')
        available_files = get_emis_file_list(Path(self.path) / self.domain, '*.nc')
        self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
        #if len(available_files) > 0:
        self.filename = self.param.filename.objects[0]

    @param.depends('filename', 'fieldname', watch=True)
    def update_field_description(self):
        available_files = get_emis_file_list(Path(self.path) / self.domain, f'{self.filename}*.nc*')
        if len(available_files) > 0:
            ds = xr.open_dataset(available_files[0])

            self.widgets['info'].object = f"""
            {ds.attrs.get('description', '`file description missing`')}
            
            **{self.fieldname}**
            - *long_name*\t: {ds[self.fieldname].long_name}
            - *units*\t: {ds[self.fieldname].units}
            """

    @param.depends('desc', watch=True)
    def update_desc(self):
        self.widgets['title'].object = f'### {self.desc}'


class EmissionSettings(pn.viewable.Viewer):
    catname = param.String(doc='name of the emission category (should be unique to that tracer)')
    regions = param.List(doc='region(s) where the emissions should be applied')
    path = param.Path(doc='location of the emission files')
    #emis_reg = FieldSelector(desc='Emissions for the regional domain')
    #emis_glo = FieldSelector(desc='Global emissions')
    switch_reg = param.Boolean(doc="Switch alternate source for regional emissions")

    def __init__(self, **params):
        super().__init__(**params)
        self.emis_reg = FieldSelector(desc='Emissions for the regional domain', domain=self.regions[-1])
        self.emis_glo = FieldSelector(desc='Global emissions', domain=self.regions[0])
        self.emis_glo.path = self.path
        self.emis_reg.path = self.path
        self.pane_glo = pn.Column(self.emis_glo)
        self.pane_reg = pn.Column(self.emis_reg, visible=len(self.regions) > 1)
        self.switch_button = pn.Row(
                    pn.widgets.Switch.from_param(self.param.switch_reg, align='center'),
                    pn.pane.Markdown("Use different emissions for the zoom domain"),
                    visible = len(self.regions) > 1
                )
        self.update_visibility_regional_emissions()

    def __panel__(self):
        return pn.Row(
            pn.widgets.TextInput.from_param(self.param.catname),
            pn.Row(
                pn.Column(
                    self.pane_glo,
                    self.switch_button),
                self.pane_reg,
                sizing_mode='stretch_width'
            )
        )

    @param.depends('regions', 'switch_reg', watch=True)
    def update_visibility_regional_emissions(self):
        if len(self.regions) > 1 and self.switch_reg:
            self.emis_reg.desc = f"Emissions for region *{self.regions[-1]}*"
            self.pane_reg.visible = True
        else :
            self.pane_reg.visible = False


    @param.depends('regions', watch=True)
    def update_switch_visibility(self):
        self.switch_button.visible = len(self.regions) > 1

    def copy(self):
        newem = self.__class__(
            catname = self.catname,
            regions = self.regions,
            path = self.path)
        newem.switch_reg = self.switch_reg 
        newem.emis_glo.filename = str(self.emis_glo.filename)
        newem.emis_glo.fieldname = str(self.emis_glo.fieldname)
        newem.emis_reg.filename = str(self.emis_reg.filename)
        newem.emis_reg.fieldname = str(self.emis_reg.fieldname)
        newem.update_visibility_regional_emissions()
        return newem

        
# class EmissionSettings(pn.viewable.Viewer):
#     """
#     Block of widgets controlling the settings of *one* emission category.
#     Parameters :
#         - catname ==> category name (user-defined)
#         - fieldname ==> name of the field to read from the netCDF file (user chooses from a list)
#         - regions ==> list of regions to which the emissions should be applied (user selects from a list)
#         - scf ==> scaling factor ==> for now forced to 1
#     """
#     catname = param.String(doc='name of the emission category (should be unique to that tracer)')
#     fieldname = param.Selector(doc='name of the field to be used', label='Field name')
#     filename = param.Selector(doc='name of the emission file', label='File prefix')
#     regions = param.ListSelector(doc='region(s) where the emissions should be applied')
#     scf = param.Number(doc='scaling factor for the emissions', default=1)
#     path = param.Path(doc='location of the emission files')
#     fileinfo = param.String(doc='ncdump of the file (for now ...)')

#     def __init__(self, regnames, **params):
#         super().__init__(**params)
#         self.param.regions.objects = regnames
#         self.regions = regnames
#         #self.update_field_choices()
#         self.widgets = {
#             'field': pn.widgets.Select.from_param(self.param.fieldname),
#             'info': pn.pane.Markdown()
#         }
#         self.update_file_choices()
    
#     def __panel__(self):
#         return pn.Column(
#             pn.Row(
#                 pn.widgets.TextInput.from_param(self.param.catname),
#                 pn.Column(
#                     pn.widgets.Select.from_param(self.param.filename),
#                     self.widgets['field'],
#                 ),
#                 pn.widgets.MultiChoice.from_param(self.param.regions)
#             ),
#             self.widgets['info']
#             #pn.widgets.TextInput.from_param(self.param.scf)
#         )
        
#     @param.depends('filename', 'path', watch=True)
#     def update_field_choices(self):
#         """
#         Update the choices of the "Field" widget. 
#         """
#         #available_files = list(Path(self.path).glob(f'**/{self.filename}*.nc*'))
#         available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
#         ds = xr.open_dataset(available_files[0])
#         self.param.fieldname.objects = [ _ for _ in ds.data_vars if _!='area' ]
#         self.fieldname = self.param.fieldname.objects[0]
#         self.widgets['field'].visible = len(self.param.fieldname.objects) > 1

#     @param.depends('path', watch=True)
#     def update_file_choices(self, domain='glb1x1'):
#         # available_files = get_emis_file_list(self.path, '**/*.nc*')
#         #-- 2025-04-14:: restrict here to the global (default) domain
#         fptn = f"**/*{domain}*.nc"
#         available_files = get_emis_file_list(self.path, fptn)
#         self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
#         self.filename = self.param.filename.objects[0]
        
#     @param.depends('filename', 'fieldname', watch=True)
#     def update_field_description(self):
#         available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
#         #print(available_files)
#         ds = xr.open_dataset(available_files[0])
#         self.widgets['info'].object = f"""
#         {ds.attrs.get('description', '`file description missing`')}
        
#         **{self.fieldname}**
#         - *long_name*\t: {ds[self.fieldname].long_name}
#         - *units*\t: {ds[self.fieldname].units}'
#         """
        

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
    #-- access directory containing FIT-IC input emission/flux fields
    emisdir = fitic_inputemisdir()
    emission_path = param.Path(emisdir)
    remove_event = param.Event(doc="Remove tracer", label="Remove tracer")
    duplicate_event = param.Event(doc='Duplicate tracer', label='Duplicate tracer')

    def __init__(self, parent = None, **params):
        super().__init__(**params)
        self.emissions = []
        self.emissions_widgets = pn.Column()
        self.parent = parent
        
    @param.depends('add_emissions_category', watch=True)
    def add_emis(self):
        self.emissions.append(EmissionSettings(
            catname=f'emissions_{len(self.emissions) + 1}', 
            regions=self.regions,
            path=self.emission_path
        ))
        self.emissions_widgets.append(self.emissions[-1].__panel__())

    @param.depends('regions', watch=True)
    def update_emis_region(self):
        for iemis, emis in enumerate(self.emissions):
            emis.regions = self.regions
            self.emissions_widgets[iemis] = emis.__panel__()

    @param.depends("remove_event", watch=True)
    def delete(self):
        self.parent.remove_tracer(self)

    @param.depends("duplicate_event", watch=True)
    def duplicate(self):
        self.parent.duplicate_tracer(self)

    @param.depends('tracer_name', watch=True)
    def __panel__(self):
        components = [
            pn.widgets.TextInput.from_param(self.param.tracer_name),
            pn.pane.Markdown("""
            ## Initial condition:
            
            Select whether the model should be initialized with concentrations from a previous run, from CAMS reanalysis, or left to zero. 
            """),
            pn.widgets.Select.from_param(self.param.initial_condition),
            self.reaction_widgets,
            pn.pane.Markdown("""
            ## Emissions:
            
            You can add one or multiple emission products to your simulation. The emissions are global (by default), but you have the option to use a different emission dataset within the regional domain.
            """),
            self.emissions_widgets,
            pn.widgets.Button.from_param(self.param.add_emissions_category)
        ]
        return pn.Row(
            pn.Column(
                pn.widgets.Button.from_param(self.param.remove_event),
                pn.widgets.Button.from_param(self.param.duplicate_event)
            ),
            pn.layout.Card(
                pn.Column(*components), 
                title=self.param.tracer_name, 
                styles={'background': '#edfafa'},
                sizing_mode='stretch_width',
                hide_header=True), 
            )
        # return pn.layout.Card(pn.Column(*components), title=self.param.tracer_name)

    @property
    def reaction_widgets(self):
        if len(self.reactions) > 0:
            return pn.Column(
                pn.pane.Markdown("""
                ## Chemical reactions (sink):
                
                The transport model supports tracer loss through reaction with other species. The defaults settings are always sane, but you can change them to see the impact on the results.
                """),
                *[r for r in self.reactions])


class CH4TracerSettings(TracerSettings):
    """
    Class containing settings specific to the CH4 tracers. Derived from the TracerSettings class
    
    parameters:
    - initial_condition (user speficied, from a list of choices, hard-coded for now)
    - reactions (user can activate or deactivate a specific reaction, but the list is hard-coded.
    """
    initial_condition = param.Selector(default='zero', objects=['zero', 'Cconcentration from CAMS flux inversion (v23r1)', 'previous run'], doc='initial condition')
    reactions = param.List(
        [
            ReactionSettings(
                reacname='CH4 + OH (troposphere)', 
                shortname='oh', 
                rate0=2.45e-12, 
                rate1=-1755, 
                domain='tropo', 
                versions=['Spivakovsky', 'CAMS OH field']
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

    def copy(self):
        newtr = self.__class__(
            initial_condition = self.initial_condition,
            reactions = [_.copy() for _ in self.reactions],
            tracer_name = f'{self.name}_copy',
            regions = self.regions,
            parent = self.parent
        )

        # Set the emissions
        newtr.emissions = [_.copy() for _ in self.emissions]
        newtr.emissions_widgets = pn.Column(*[_.__panel__() for _ in newtr.emissions])
        return newtr
        
        
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
    # output_types = param.ListSelector(default=['stations'], objects=['stations', 'mix', 'columns'], doc='choice of outputs (common to all tracers)')
    create_ch4_tracer = param.Event(doc='Add new CH4 tracer', label='New CH4 tracer')
    create_co2_tracer = param.Event(doc='Add new CO2 tracer', label='New CO2 tracer')
    levels = param.Selector(default='tropo34', objects=['tropo25', 'tropo34'], doc='Number of vertical levels')

    def __init__(self, **params):
        super().__init__(**params)
        self.tracers = []
        self.tracers_widgets = pn.Column()
    
    def __panel__(self):
        return pn.Column(
            pn.pane.Markdown('## Run settings'),
            pn.Row(
                pn.widgets.TextInput.from_param(self.param.run_name),
                pn.Column(
                    pn.Row(
                        pn.widgets.DatePicker.from_param(self.param.start),
                        pn.widgets.DatePicker.from_param(self.param.end),
                    ),
                    pn.widgets.Select.from_param(self.param.zoom_configuration),
                    # pn.widgets.MultiSelect.from_param(self.param.output_types),
                ),
            ),
            pn.pane.Markdown("""
            ## Tracers:
            
            The transport model (TM5) computes the transport of one or more atmospheric "tracers", i.e. atmospheric species, which are independent from each other. Set here the list of tracers you want to simulate. You can for instance configure your run with:
            - one single tracer with all processes (initial condition, sources and sinks);
            - one tracer for each source process (to be able to track the contribution of each process).
            - several tracers as "sensitivity tests"
            """),
            self.tracers_widgets,
            pn.Row(
                pn.widgets.Button.from_param(self.param.create_co2_tracer, disabled=True),
                pn.widgets.Button.from_param(self.param.create_ch4_tracer),
            ),
            styles=dict(background='#daf5f6'),
            sizing_mode='stretch_width'
        )
    
        # return pn.Column(
        #     pn.pane.Markdown('# Run settings'),
        #     pn.widgets.TextInput.from_param(self.param.run_name),
        #     pn.widgets.DatePicker.from_param(self.param.start),
        #     pn.widgets.DatePicker.from_param(self.param.end),
        #     pn.widgets.Select.from_param(self.param.zoom_configuration),
        #     pn.widgets.MultiSelect.from_param(self.param.output_types),
        #     pn.pane.Markdown('# Tracer settings'),
        #     self.tracers_widgets,
        #     pn.Row(
        #         pn.widgets.Button.from_param(self.param.create_co2_tracer, disabled=True),
        #         pn.widgets.Button.from_param(self.param.create_ch4_tracer),
        #     )
        # )
        
    @param.depends('create_ch4_tracer', watch=True)
    def add_ch4_tracer(self):
        self.add_tracer('CH4')
        
    @param.depends('create_co2_tracer', watch=True)
    def add_co2_tracer(self):
        self.add_tracer('CO2')

    def update_tracer_widgets(self):
        self.tracers_widgets.objects = [_.__panel__() for _ in self.tracers]

    def remove_tracer(self, tracer):
        self.tracers.remove(tracer)
        self.update_tracer_widgets()

    def duplicate_tracer(self, tracer):
        self.tracers.append(tracer.copy())
        self.update_tracer_widgets()

    def add_tracer(self, species: str):
        trname = self.get_unique_trname(species)
        trclass = dict(CH4=CH4TracerSettings, CO2=CO2TracerSettings)[species]
        self.tracers.append(
            trclass(
                tracer_name=trname, 
                regions=self.zoom_configuration, 
                parent=self
        ))
#        self._tracers.append(trname)
        self.update_tracer_widgets()
        #self.tracers.append(trclass(tracer_name=trname, regions=self.zoom_configuration))
        #self.tracers_widgets.append(self.tracers[-1].__panel__())

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
    #run_tm5_button = param.Event(doc='run TM5', label='Run TM5')
    #build_tm5_button = param.Event(doc='compile TM5', label='Compile TM5')
    submit_event = param.Event(doc='submit tm5', label='Submit a new run')
    check_status_event = param.Event(doc='Check status', label='Check status')
    jobid = param.Integer(doc='TM5 job ID')
    
    def __init__(self, **params):
        super().__init__(**params)
        self.settings = RunSettings()
        #self.terminal = pn.widgets.Terminal(options={"cursorBlink": True}, height=300, sizing_mode='stretch_width', write_to_console=True)
        #self.row = pn.Row()
        self.textbox = pn.pane.Alert(visible=False, width=300)
        self.terminal = pn.widgets.Terminal(options={"cursorBlink": True}, height=300, sizing_mode='stretch_width', write_to_console=True, visible=False)
        self.submit_button = pn.widgets.Button.from_param(self.param.submit_event)
        self.check_status_button = pn.widgets.Button.from_param(self.param.check_status_event, visible=False)
        self.job_selector = pn.widgets.IntInput.from_param(self.param.jobid, visible=False, name='', width=60)
        self.alert = pn.pane.Alert(visible=False)

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
        # return pn.Column(
        #     pn.widgets.Select.from_param(self.param.rcfile),
        #     self.settings,
        #     pn.layout.Divider(),
        #     pn.Row(
        #         pn.widgets.Button.from_param(self.param.run_tm5_button),
        #         pn.widgets.Button.from_param(self.param.build_tm5_button),
        #     ),
        #     self.row,
        #     self.terminal    
        # )

    @param.depends('submit_event', watch=True)
    def submit(self):

        # Just in case, clear and hide the terminal and "textbox"
        self.terminal.clear()
        self.terminal.visible = False
        self.textbox.object = None
        self.textbox.visible = False

        try :
            url = f'{url_tm5}/submit'
            r = requests.get(url, params={'config':'toto'})
            
            if not r.ok:
                self.alert.visible = True
                self.alert.alert_type = 'danger'
                self.alert.object = f'Submission failed. Incorrect request to {r.url} ☠️. '
                return

            self.check_status_button.visible=True
            self.job_selector.visible=True

            # Set the button color according to the submit status            
            self.check_status_button.button_type = {
                'finished':'success',
                'queued':'warning',
                'running':'primary'
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
        r = requests.get(f'{url_tm5}/status/{self.jobid}')

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
            #with open(r.json()['outfile'], 'r') as fid:
            #    self.terminal.writelines(fid.readlines())

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


FitIC_UI().servable()
