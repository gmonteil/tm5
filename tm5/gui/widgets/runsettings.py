import panel as pn
import param
import os
from datetime import date
from tm5.gui.widgets import CH4TracerSettings, CO2TracerSettings
from copy import deepcopy


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
        self.update_tracer_widgets()

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
