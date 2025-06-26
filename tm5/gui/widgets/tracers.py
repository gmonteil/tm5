import param
import panel as pn
from tm5.gui.widgets.emissions import EmissionSettings
from tm5.gui.widgets.reactions import ReactionSettings
from tm5.gui import host


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
            path=host.emission_path
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
            #-- MVO::deactivated for now in the GUI since stratosphere chemistry
            #        is not covered yet in TM5 simulation
            # ReactionSettings(
            #     reacname='CH4 + OH (stratosphere)', 
            #     shortname='ohstrat', 
            #     rate0=2.45e-12, 
            #     rate1=-1755, 
            #     domain='strato', 
            #     versions=['Bruehl']
            # ),
            # ReactionSettings(
            #     reacname='CH4 + O(1D)', 
            #     shortname='o1d', 
            #     rate0=1.5e-10, 
            #     domain='strato', 
            #     versions=['Bruehl']
            # ),
            # ReactionSettings(
            #     reacname='CH4 + Cl', 
            #     shortname='cl', 
            #     rate0=7.3e-12, 
            #     rate1=-1280, 
            #     domain='strato', 
            #     versions=['Bruehl']
            # )
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
