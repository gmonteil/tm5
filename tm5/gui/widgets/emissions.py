import panel as pn
import param
from pathlib import Path
from functools import lru_cache
from typing import List
import xarray as xr
from loguru import logger
from tm5.gui.css import *

def get_emis_file_list(path: Path, pattern: str) -> List[Path]:
    if (path / pattern).is_file():
        return [path / pattern]
    if not pattern.endswith('*.nc'):
        pattern += '*.nc'
    return list(Path(path).glob(pattern))


@lru_cache
def get_emis_dataset(path: Path) -> xr.Dataset:
    return xr.open_dataset(path)


class FieldSelector(pn.viewable.Viewer):
    catname = param.String(doc='category name')
    filename = param.Selector(doc="name of the emission file")
    fieldname = param.Selector(doc="name of the field to be used")
    search_paths = param.List(default=[], doc='directories where the emission files are searched, by order of priority (a file name found in several of them is taken from the first)')
    desc = param.String(doc="domain of the emissions")
    domain = param.String(doc="title of the section")
    visible = param.Boolean(default=True, doc='False for anything but the "custom" scenario')

    def __init__(self, editable : bool = False, **params):
        super().__init__(**params)
        self.editable = editable
        self.widgets = dict(
            file=pn.widgets.Select.from_param(self.param.filename),
            field=pn.widgets.Select.from_param(self.param.fieldname),
            info=pn.pane.Markdown(width=300),
            title=pn.pane.Markdown(width=300),
        )
        self.widgets['file'].disabled = not editable
        self.widgets['field'].disabled = not editable
        self.update_desc()
        self.update_widgets_visibility()

    @param.depends('visible', watch=True)
    def update_widgets_visibility(self):
        self.widgets['file'].visible = self.visible
        self.widgets['field'].visible = (len(self.param.fieldname.objects) > 1) and self.visible
        self.widgets['info'].visible = self.visible
        self.widgets['title'].visible = self.visible

    def __panel__(self):
        return pn.Column(
            self.widgets['title'],
            self.widgets['file'],
            self.widgets['field'],
            self.widgets['info'],
            stylesheets=[setup_stylesheet,], css_classes=['setup-tracer']
        )

    def _list_files(self, pattern: str) -> List[Path]:
        found = {}
        for path in self.search_paths:
            for file in get_emis_file_list(Path(path), pattern):
                found.setdefault(file.stem, file)
        return list(found.values())

    @property
    def path(self) -> Path | None:
        """
        Directory containing the selected file (it can be any of the search paths).
        """
        for path in self.search_paths:
            if self.filename is not None and len(get_emis_file_list(Path(path), self.filename)) > 0:
                return Path(path)
        return Path(self.search_paths[0]) if len(self.search_paths) > 0 else None

    @param.depends('filename', 'search_paths', 'domain', watch=True)
    def update_field_choices(self):
        """
        Update the choices of the "Field" widget.
        """
        if self.filename==None:
            return
        # msg = f"@{self.filename}, self.domain ==>{self.domain}<=="
        # logger.debug(msg)
        available_files = self._list_files(self.filename)
        if len(available_files) > 0:
            ds = get_emis_dataset(available_files[0])
            self.param.fieldname.objects = [_ for _ in ds.data_vars if _ != 'area']
            self.fieldname = self.param.fieldname.objects[0]
            self.update_widgets_visibility()

    @param.depends('search_paths', 'domain', watch=True)
    def update_file_choices(self):
        # msg = f"self.domain ==>{self.domain}<=="
        # logger.debug(msg)
        #-- NOTE::introduced following naming convention to
        #         differentiate between global and regional emissions:
        #         global_xxx_yyy_*.nc
        #         regional_xxx_yyy_*.nc
        ptn = f"{self.domain}_*.nc"
        available_files = self._list_files(ptn)
        # msg = f"available_files ==>{available_files}<=="
        # logger.debug(msg)
        #-- drop '.nc' extension for the selection
        selectable_files = sorted([_.stem for _ in available_files])
        # msg = f"selectable_files ==>{selectable_files}<=="
        # logger.debug(msg)
        self.param.filename.objects = selectable_files
        if len(selectable_files) > 0:
            self.filename = selectable_files[0]
        else:
            self.filename = None
            searched = ', '.join(f'`{p}`' for p in self.search_paths)
            self.widgets['info'].object = f'*No `{ptn}` emission file found in {searched}*'

    def refresh_file_choices(self):
        """
        Pick up files added to the search paths since the selector was created (e.g. uploads),
        without changing the current selection.
        """
        found = [_.stem for _ in self._list_files(f"{self.domain}_*.nc")]
        known = [_ for _ in self.param.filename.objects if _ is not None]
        choices = sorted(set(known) | set(found))
        if choices == known:
            return
        self.param.filename.objects = choices
        if self.filename is None and len(choices) > 0:
            self.filename = choices[0]

    @param.depends('filename', 'fieldname', watch=True)
    def update_field_description(self):
        if self.filename==None:
            return
        available_files = self._list_files(self.filename)
        if len(available_files) > 0:
            ds = get_emis_dataset(available_files[0])

            if 'comment' in ds[self.fieldname].attrs:
                self.widgets['info'].object = f"""
                {ds.attrs.get('description', '`file description missing`')}

                **{self.fieldname}**
                - *long_name*\t: {ds[self.fieldname].attrs.get('long_name', ds[self.fieldname].attrs.get('description', '`missing`'))}
                - *units*\t: {ds[self.fieldname].attrs.get('units', '`missing`')}
                - *comment*\t: {ds[self.fieldname].comment}
                """
            else:
                self.widgets['info'].object = f"""
                {ds.attrs.get('description', '`file description missing`')}

                **{self.fieldname}**
                - *long_name*\t: {ds[self.fieldname].attrs.get('long_name', ds[self.fieldname].attrs.get('description', '`missing`'))}
                - *units*\t: {ds[self.fieldname].attrs.get('units', '`missing`')}
                """

    @param.depends('desc', watch=True)
    def update_desc(self):
        self.widgets['title'].object = f'### {self.desc}'

    def set_selection(self, filename: str, fieldname: str = None):
        """
        Directly set filename/fieldname (e.g. from a preconfigured scenario)
        """
        # Add the "filename" from the yaml file to the list of available filenames for that category
        # (and set it as the active one)
        if filename not in self.param.filename.objects:
            self.param.filename.objects = list(self.param.filename.objects) + [filename]
        self.filename = filename

        # Same for the fields
        if fieldname is not None:
            self.fieldname = fieldname


class EmissionSettings(pn.viewable.Viewer):
    catname = param.String(doc='name of the emission category (should be unique to that tracer)')
    regions = param.List(doc='region(s) where the emissions should be applied')
    path = param.Path(doc='location of the emission files')
    extra_paths = param.List(default=[], doc='additional locations searched for emission files (e.g. the folder of uploaded files)')
    # emis_reg = FieldSelector(desc='Emissions for the regional domain')
    # emis_glo = FieldSelector(desc='Global emissions')
    switch_reg = param.Boolean(doc="Switch alternate source for regional emissions")
    remove_event = param.Event(doc='Remove this emission category', label='Remove category')
    visible = param.Boolean(default=True, doc='False for preconfigured emission scenarios')

    def __init__(self, remove_callback: callable, editable : bool = False, **params):
        super().__init__(**params)
        self.removeme = remove_callback  # method of the parent object that needs to be called when removing the category (see _handle_remove method below)
        self.editable = editable
        self.emis_reg = FieldSelector(desc='Emissions for the regional domain', domain=self.regions[-1], visible=True, editable=self.editable)
        self.emis_glo = FieldSelector(desc='Global emissions', domain=self.regions[0], visible=True, editable=self.editable)
        search_paths = [self.path] + [p for p in self.extra_paths if p is not None]
        self.emis_glo.search_paths = search_paths
        self.emis_reg.search_paths = search_paths
        self.pane_glo = pn.Column(self.emis_glo, stylesheets=[setup_stylesheet,], css_classes=['setup-tracer'])
        self.pane_reg = pn.Column(self.emis_reg, visible=len(self.regions) > 1)
        self.widgets = dict(
            catname=pn.widgets.TextInput.from_param(self.param.catname),
            remove=pn.widgets.Button.from_param(self.param.remove_event),
            switch=pn.widgets.Switch.from_param(self.param.switch_reg, align='center'),
        )
        
        self.switch_button = pn.Row(
            self.widgets['switch'],
            pn.pane.Markdown("Use different regional emissions", stylesheets=[setup_stylesheet,], css_classes=['setup-tracer']),
        )
        self.layout = pn.Row(
            pn.Column(
                self.widgets['catname'],
                self.widgets['remove'],
            ),
            pn.Row(
                pn.Column(
                    self.pane_glo,
                    self.switch_button),
                self.pane_reg,
                sizing_mode='stretch_width'
            ),
            stylesheets=[setup_stylesheet,], css_classes=['setup-tracer'],
            margin=(5, 0),
        )
        self.widgets['catname'].editable = False
        self.widgets['remove'].editable = False
        self.widgets['switch'].editable = False
        self.widgets['catname'].disabled = not editable
        self.widgets['remove'].disabled = not editable
        self.widgets['switch'].disabled = not editable
        self.update_visibility_regional_emissions()
        self.update_visible()

    def __panel__(self):
        return self.layout

    def refresh_file_choices(self):
        self.emis_glo.refresh_file_choices()
        self.emis_reg.refresh_file_choices()

    @param.depends('visible', watch=True)
    def update_visible(self):
        self.layout.visible = self.visible

    @param.depends('regions', 'switch_reg', watch=True)
    def update_visibility_regional_emissions(self):
        if len(self.regions) > 1 and self.switch_reg:
            self.emis_reg.desc = f"*{self.regions[-1]}* emissions"
            self.pane_reg.visible = True
        else:
            self.pane_reg.visible = False

    @param.depends('remove_event', watch=True)
    def _handle_remove(self):
        self.removeme(self)

    def set_category(self, spec: dict):
        """
        Force a category to a certain value (instead of letting it happen through widget inputs).
        This is needed when loading pre-defined emission scenarios
        """
        self.emis_glo.set_selection(spec['global']['filename'], spec['global'].get('field'))
        reg = spec.get('regional')
        if reg:
            self.switch_reg = True
            self.emis_reg.set_selection(reg['filename'], reg.get('field'))

    @param.depends('regions', watch=True)
    def update_switch_visibility(self):
        self.switch_button.visible = len(self.regions) > 1

    def copy(self):
        newem = self.__class__(
            catname=self.catname,
            regions=self.regions,
            path=self.path,
            extra_paths=self.extra_paths,
            visible=self.visible,
            remove_callback=self.removeme)
        newem.switch_reg = self.switch_reg
        newem.emis_glo.filename = str(self.emis_glo.filename)
        newem.emis_glo.fieldname = str(self.emis_glo.fieldname)
        newem.emis_reg.filename = str(self.emis_reg.filename)
        newem.emis_reg.fieldname = str(self.emis_reg.fieldname)
        newem.update_visibility_regional_emissions()
        return newem

    def getconfig(self):
        emis_dict = {}
        globname = self.emis_glo.filename
        globfield = self.emis_glo.fieldname
        emis_dict[self.catname] = {
            'global': { 'filename': str(self.emis_glo.filename),
                        'field': str(self.emis_glo.fieldname) }
            }
        if self.switch_reg:
            emis_dict[self.catname]['regional'] = {
                'filename': str(self.emis_reg.filename),
                'field': str(self.emis_reg.fieldname)
            }
        # msg = f"returning emis_dict ***{emis_dict}***"
        # logger.debug(msg)
        return emis_dict
