import panel as pn
import param
from pathlib import Path
from functools import lru_cache
from typing import List
import xarray as xr
from loguru import logger


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
        available_files = get_emis_file_list(Path(self.path) / self.domain, '*.nc')
        self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
        #if len(available_files) > 0:
        self.filename = self.param.filename.objects[0]

    @param.depends('filename', 'fieldname', watch=True)
    def update_field_description(self):
        available_files = get_emis_file_list(Path(self.path) / self.domain, f'{self.filename}*.nc*')
        if len(available_files) > 0:
            ds = xr.open_dataset(available_files[0])

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
                    pn.pane.Markdown("Use different regional emissions"),
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
#
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
#
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
#
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
#
#     @param.depends('path', watch=True)
#     def update_file_choices(self):
#         available_files = get_emis_file_list(self.path, '**/*.nc*')
#         self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
#         self.filename = self.param.filename.objects[0]
#
#     @param.depends('filename', 'fieldname', watch=True)
#     def update_field_description(self):
#         available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
#         #print(available_files)
#         ds = xr.open_dataset(available_files[0])
#         self.widgets['info'].object = f"""
#         {ds.attrs.get('description', '`file description missing`')}
#
#         **{self.fieldname}**
#         - *long_name*\t: {ds[self.fieldname].long_name}
#         - *units*\t: {ds[self.fieldname].units}'
#         """
