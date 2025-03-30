import param
import panel as pn
from pathlib import Path
from tm5.gui.widgets import EmissionSettings, ReactionSettings, emission_dir


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