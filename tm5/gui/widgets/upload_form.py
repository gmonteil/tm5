import re
import traceback
import numpy as np
import panel as pn
import param
from pathlib import Path
import xarray as xr


class VariableSelector(pn.viewable.Viewer):
    varname = param.String(doc='Name of the variable')
    newname = param.String(doc='Import the variable as', label='Import variable as:')
    description = param.String(doc='Variable description')
    var_selected = param.Boolean(doc='Variable selected for import', label='')

    def __init__(self, da: xr.DataArray, **params):
        super().__init__(**params)
        self.varname = da.name
        self.newname = da.name
        self.description = da.attrs.get('description', '')

    def __panel__(self):
        return pn.Row(
            pn.pane.Markdown(f'Import variable **{self.varname}**:', align='end', stylesheets=['p { margin: 0}']),
            pn.widgets.Checkbox.from_param(self.param.var_selected, align='end'),
            pn.widgets.TextInput.from_param(self.param.newname, align='end'),
            pn.widgets.TextInput.from_param(self.param.description, align='end'),
            align='end'
        )


class FileValidator(pn.viewable.Viewer):
    dataset_name = param.String(doc='The dataset will appear under this name once imported in the FIT-IC.')
    dataset_description = param.String(doc='The description will appear in the FIT-IC interface when this dataset is selected')
    filename = param.String(doc='Name of the file as uploaded')
    error = param.String(default='', doc='Generic object for error messages')
    remove_button = param.Event(doc='Remove this file from the list', label='Remove')
    import_button = param.Event(doc='Validate the submission. The files will become available in the Simulations tab', label='Import into FIT-IC')

    lat_dim = param.Selector(doc='Name of the latitude dimension')
    lon_dim = param.Selector(doc='Name of the longitude dimension')
    time_dim = param.Selector(doc='Name of the time dimension')

    def __init__(self, content: bytes, filename, remove_callback: callable, import_callback: callable, upload_path: Path, **params):
        super().__init__(**params)
        self.upload_path = Path(upload_path)
        self.content = content  # the raw uploaded bytes, handed over once, no link back to the FileDropper
        self.filename = filename
        self.variables = []
        self.removeme = remove_callback  # method of the parent object that needs to be called when this file is removed (see _handle_remove below)
        self.imported = import_callback  # method of the parent object called with (self, target_path) once the file has been written
        self.widgets = dict(
            label=pn.pane.Markdown(f'**{self.filename}**'),
            dim_selector = pn.Row(
                pn.widgets.Select.from_param(self.param.lat_dim),
                pn.widgets.Select.from_param(self.param.lon_dim),
                pn.widgets.Select.from_param(self.param.time_dim),
                visible=False),
            variables=pn.Column(),
            import_button = pn.widgets.Button.from_param(self.param.import_button, disabled=True)
        )
        self.check_file()

    @param.depends('dataset_name', 'error', watch=True)
    def _update_importable_status(self):
        disabled = (
            (len(self.dataset_name.strip()) > 0) &
            (len(self.error) == 0)
        )
        self.widgets['import_button'].disabled = not disabled

    def __panel__(self):
        return pn.Column(
            self.widgets['label'],
            self._alert,
            pn.Row(
                pn.widgets.TextInput.from_param(self.param.dataset_name, placeholder='Give a name to your dataset'),
                pn.widgets.TextAreaInput.from_param(self.param.dataset_description, placeholder='Write a description for your dataset. It will appear in the FIT-IC interface once selected', auto_grow=True)
            ),
            self.widgets['dim_selector'],
            self.widgets['variables'],
            pn.Row(
                self.widgets['import_button'],
                pn.widgets.Button.from_param(self.param.remove_button)
            )
        )

    @param.depends('import_button', watch=True)
    def import_file(self):
        selected = [v for v in self.variables if v.var_selected]
        if len(selected) == 0:
            self.error = 'Please select at least one variable to import.'
            return

        name = self.dataset_name.strip()
        if re.fullmatch(r'[A-Za-z0-9_.\-]+', name) is None or name.startswith('.'):
            self.error = 'The dataset name may only contain letters, digits, "_", "-" and ".".'
            return

        with xr.open_dataset('emis.tmp.nc') as ds:
            ds = ds[[v.varname for v in selected]].rename({
                self.lat_dim: 'latitude',
                self.lon_dim: 'longitude',
                self.time_dim: 'time',
            })
            ds = ds.rename({v.varname: v.newname for v in selected})
            for v in selected:
                ds[v.newname].attrs['description'] = v.description
            ds.attrs['description'] = self.dataset_description

            # Coordinates are cell centers of a 1-degree grid, so the boundaries are +/- 0.5 degree
            lon, lat = ds['longitude'], ds['latitude']
            is_global = (
                np.isclose(float(lon.max() - lon.min()) + 1.0, 360.0)
                and np.isclose(float(lat.max()) + 0.5, 90.0)
                and np.isclose(float(lat.min()) - 0.5, -90.0)
            )
            prefix = 'global' if is_global else 'regional'

            self.upload_path.mkdir(parents=True, exist_ok=True)
            target = self.upload_path / f'{prefix}_{name}.nc'
            if target.exists():
                self.error = f'{target.name} already exists. Choose another dataset name.'
                return
            ds.load().to_netcdf(target)

        self.error = ''
        self.imported(self, target)

    @param.depends('remove_button', watch=True)
    def _handle_remove(self):
        self.removeme(self)

    def check_file(self):
        """
        Try opening the file with xarray. If it fails, the file is invalid, write a warning. If valid, displays the list of variables inside, allowing the userr to choose. This is just a placeholder for something more advanced (we can do some checks on the resolution, allow user to crop the file, provide metadata, etc.
        """
        self.error = ''

        # 1. Check that the file is a valid netCDF file
        tmp_path = 'emis.tmp.nc'
        with open(tmp_path, 'wb') as f:
            f.write(self.content)

        try:
            ds = xr.open_dataset(tmp_path)
        except ValueError:
            self.error = 'Not a netCDF file'
            return

        # 2. Check the dimensions of the file
        if not self._check_dims(ds):
            return

        # 3. Fill in the "description" field, if already available
        if 'description' in ds.attrs:
            self.dataset_description = ds.attrs['description']

        # 3. Select the variables
        self._select_variables(ds)

    def _check_dims(self, ds: xr.Dataset):
        dims_in_file = list(ds.sizes)
        if len(dims_in_file) < 3:
            self.error = f"The emissions should have three dimensions (lat, lon, time). Only {len(dims_in_file)} detected in the file: {dims_in_file}. The file won't be imported."
            return False

        self.param.lon_dim.objects = dims_in_file
        self.param.lat_dim.objects = dims_in_file
        self.param.time_dim.objects = dims_in_file
        for dimname in dims_in_file:
            if 'lat' in dimname:
                self.lat_dim = dimname
            if 'lon' in dimname:
                self.lon_dim = dimname
            if 'time' in dimname:
                self.time_dim = dimname

        self.widgets['dim_selector'].visible = True
        return True

    def _select_variables(self, ds: xr.Dataset):
        dimnames = {self.lat_dim, self.lon_dim, self.time_dim}
        varnames = [v for v in ds.data_vars if len(dimnames - set(ds[v].sizes)) == 0]
        for varname in varnames:
            variable = VariableSelector(ds[varname])
            self.widgets['variables'].append(variable)
            self.variables.append(variable)

    @param.depends('error')
    def _alert(self):
        if self.error == '':
            return ''
        return pn.pane.Alert(self.error, alert_type='danger')


class UploadEmissionForm(pn.viewable.Viewer):
    def __init__(self, upload_path: str | Path, **params):
        super().__init__(**params)
        self.upload_path = Path(upload_path)
        self.validators = []  # plain list of FileValidator, mirrors tracers.TracerSettings.emissions
        self.validators_widgets = pn.Column()  # display-only, rebuilt from validators

        # The "FileDropper" widget doesn't use the same logic as the other ones, it modifies the content of its "value" dict
        # whereas the others modify the object itself, so it doesn't change reference, which makes changes impossible to track
        # through the "@panel.depends" logic.
        self.file_chooser = self._new_file_chooser()
        self.alert = pn.pane.Alert('', visible=False)
        self.layout = pn.Column(self.file_chooser, self.alert, self.validators_widgets)

    def __panel__(self):
        return self.layout

    def _new_file_chooser(self):
        # Small chunks keep each websocket message below the size limit of the Jupyter server (the Panel default of 10 MB is too large once base64-encoded)
        file_chooser = pn.widgets.FileDropper(multiple=True, chunk_size=1_000_000)
        file_chooser.param.watch(self.on_files_change, 'value')
        return file_chooser

    def on_files_change(self, event):
        if not event:
            return

        # Exceptions raised in a callback only reach the server log (invisible in a notebook), so show them in the UI
        self.alert.visible = False
        try:
            self._add_validators(event)
        except Exception:
            self.alert.alert_type = 'danger'
            self.alert.object = f'```\n{traceback.format_exc()}\n```'
            self.alert.visible = True

        # Reset the FileDropper widget, since we don't want to use it for removing the files (if needed).
        self.file_chooser = self._new_file_chooser()
        self.layout[0] = self.file_chooser

    def _add_validators(self, event):
        # Find files that need to be added to the "validators" list
        # (removal is handled by FileValidator itself, via _remove_validator)
        existing = {v.filename for v in self.validators}
        for filename in set(event.new) - existing:
            validator = FileValidator(self.file_chooser.value[filename], filename, remove_callback=self._remove_validator, import_callback=self._on_imported, upload_path=self.upload_path)
            self.validators.append(validator)
            self.validators_widgets.append(validator.__panel__())

    def _on_imported(self, validator, target: Path):
        self._remove_validator(validator)
        self.alert.alert_type = 'success'
        self.alert.object = f'**{validator.filename}** was imported as `{target}`.'
        self.alert.visible = True

    def _remove_validator(self, validator):
        self.validators.remove(validator)
        self.validators_widgets.objects = [v.__panel__() for v in self.validators]