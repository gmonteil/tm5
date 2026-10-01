import panel as pn
import param
from pathlib import Path
import xarray as xr


class UploadEmissionForm(pn.viewable.Viewer):
    file_chooser = param.Parameter(doc='Choose an emission file to upload. It should ne in netCDF format')
    move_event = param.Event(doc='Validate the submission. The files will become available in the Simulations tab', label='Validate')
    cancel_event = param.Event(label='Cancel')
    error = param.String(default='', doc='Generic object for error messages')
    selected_vars = param.ListSelector(default=[], objects=[], doc='Data variables selected for import')

    def __init__(self, **params):
        super().__init__(**params)
        self.widgets = dict(
            file_chooser=pn.widgets.FileDropper.from_param(self.param.file_chooser),
            move_button=pn.widgets.Button.from_param(self.param.move_event),
            cancel_button=pn.widgets.Button.from_param(self.param.cancel_event),
            varlist=pn.widgets.CheckBoxGroup.from_param(self.param.selected_vars, visible=False),
        )

    def __panel__(self):
        return pn.Column(
            self.widgets['file_chooser'],
            self._alert,
            self.widgets['varlist'],
            pn.Row(
                self.widgets['move_button'],
                self.widgets['cancel_button']
            )
        )

    @param.depends('file_chooser', watch=True)
    def validate_file(self):
        """
        Try opening the file with xarray. If it fails, the file is invalid, write a warning. If valid, displays the list of variables inside, allowing the userr to choose. This is just a placeholder for something more advanced (we can do some checks on the resolution, allow user to crop the file, provide metadata, etc.
        """
        self.error = ''
        if not self.file_chooser:
            return

        try:
            # One issue is that the file chooser allows several files. I haven't
            # had time to find a cleaner solution. Currently causes bug if trying to
            # upload successively one file in the wrong format, then one in the good
            # format (the 2nd isn't detected). Probably also buggy if trying to 
            # upload several files in a sequence ...
            file_bytes = list(self.file_chooser.values())[0]
            
            # intentionally a named file, for debug. will be switched to a proper
            # tempfile later
            tmp_path = 'temp_emis.nc'  
            with open(tmp_path, 'wb') as f:
                f.write(file_bytes)

            # very basic check: it should be openable with xarray, if not, warn.
            ds = xr.open_dataset(tmp_path)
            varlist = list(ds.data_vars)

        except Exception:
            self.error = 'File not recognized'
            self.param.selected_vars.objects = []
            self.selected_vars = []
            self.widgets['varlist'].visible = False
            return

        # Set the widget values (all vars selected by default)
        self.param.selected_vars.objects = varlist
        self.selected_vars = varlist

        # Widget was invisible as an empty list looks kind of buggy on screen
        self.widgets['varlist'].visible = True

    @param.depends('error')
    def _alert(self):
        if self.error == '':
            return ''
        return pn.pane.Alert(self.error, alert_type='danger')
