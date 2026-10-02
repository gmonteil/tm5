import panel as pn
import param
from pathlib import Path
from functools import lru_cache
from typing import List
import xarray as xr
from loguru import logger
from tm5.gui.css import *

@lru_cache
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
    path = param.Path(doc='location of the emission files')
    desc = param.String(doc="domain of the emissions")
    domain = param.String(doc="title of the section")
    visible = param.Boolean(default=True, doc='False for anything but the "custom" scenario')

    def __init__(self, **params):
        super().__init__(**params)
        self.widgets = dict(
            file=pn.widgets.Select.from_param(self.param.filename),
            field=pn.widgets.Select.from_param(self.param.fieldname),
            info=pn.pane.Markdown(width=300),
            title=pn.pane.Markdown(width=300),
        )
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

    @param.depends('filename', 'path', 'domain', watch=True)
    def update_field_choices(self):
        """
        Update the choices of the "Field" widget.
        """
        if self.filename==None:
            return
        # msg = f"@{self.filename}, self.domain ==>{self.domain}<=="
        # logger.debug(msg)
        available_files = get_emis_file_list(Path(self.path), self.filename)
        if len(available_files) > 0:
            ds = get_emis_dataset(available_files[0])
            self.param.fieldname.objects = [_ for _ in ds.data_vars if _ != 'area']
            self.fieldname = self.param.fieldname.objects[0]
            self.update_widgets_visibility()

    @param.depends('path', 'domain', watch=True)
    def update_file_choices(self):
        # msg = f"self.domain ==>{self.domain}<=="
        # logger.debug(msg)
        #-- NOTE::introduced following naming convention to
        #         differentiate between global and regional emissions:
        #         global_xxx_yyy_*.nc
        #         regional_xxx_yyy_*.nc
        ptn = f"{self.domain}_*.nc"
        available_files = get_emis_file_list(Path(self.path), ptn)
        #-- drop '.nc' extension for the selection
        selectable_files = sorted([_.stem for _ in available_files])
        # self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
        self.param.filename.objects = selectable_files
        self.filename = self.param.filename.objects[0]

    @param.depends('filename', 'fieldname', watch=True)
    def update_field_description(self):
        if self.filename==None:
            return
        available_files = get_emis_file_list(Path(self.path), self.filename)
        if len(available_files) > 0:
            ds = get_emis_dataset(available_files[0])

            if 'comment' in ds[self.fieldname].attrs:
                self.widgets['info'].object = f"""
                {ds.attrs.get('description', '`file description missing`')}

                **{self.fieldname}**
                - *long_name*\t: {ds[self.fieldname].long_name}
                - *units*\t: {ds[self.fieldname].units}
                - *comment*\t: {ds[self.fieldname].comment}
                """
            else:
                self.widgets['info'].object = f"""
                {ds.attrs.get('description', '`file description missing`')}

                **{self.fieldname}**
                - *long_name*\t: {ds[self.fieldname].long_name}
                - *units*\t: {ds[self.fieldname].units}
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
    # emis_reg = FieldSelector(desc='Emissions for the regional domain')
    # emis_glo = FieldSelector(desc='Global emissions')
    switch_reg = param.Boolean(doc="Switch alternate source for regional emissions")
    remove_event = param.Event(doc='Remove this emission category', label='Remove category')
    visible = param.Boolean(default=True, doc='False for anything but the "custom" scenario')

    def __init__(self, remove_callback: callable, **params):
        super().__init__(**params)
        self.removeme = remove_callback  # method of the parent object that needs to be called when removing the category (see _handle_remove method below)
        self.emis_reg = FieldSelector(desc='Emissions for the regional domain', domain=self.regions[-1], visible=True)
        self.emis_glo = FieldSelector(desc='Global emissions', domain=self.regions[0], visible=True)
        self.emis_glo.path = self.path
        self.emis_reg.path = self.path
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
        self.update_visibility_regional_emissions()
        self.update_visible()

    def __panel__(self):
        return self.layout

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
            visible=self.visible,
            remove_callback=self.removeme)
        newem.switch_reg = self.switch_reg
        newem.emis_glo.filename = str(self.emis_glo.filename)
        newem.emis_glo.fieldname = str(self.emis_glo.fieldname)
        newem.emis_reg.filename = str(self.emis_reg.filename)
        newem.emis_reg.fieldname = str(self.emis_reg.fieldname)
        newem.update_visibility_regional_emissions()
        return newem
