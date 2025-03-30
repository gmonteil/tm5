import panel as pn
import param
from omegaconf import OmegaConf
from tm5.gui.widgets import RunSettings


pn.extension()
pn.extension('terminal')


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