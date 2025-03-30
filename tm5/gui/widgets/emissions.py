import panel as pn
import param
from pathlib import Path
from functools import lru_cache
from typing import List
import xarray as xr


@lru_cache
def get_emis_file_list(path: Path, pattern: str) -> List[Path]:
    return list(Path(path).glob(pattern))


class EmissionSettings(pn.viewable.Viewer):
    """
    Block of widgets controlling the settings of *one* emission category.
    Parameters :
        - catname ==> category name (user-defined)
        - fieldname ==> name of the field to read from the netCDF file (user chooses from a list)
        - regions ==> list of regions to which the emissions should be applied (user selects from a list)
        - scf ==> scaling factor ==> for now forced to 1
    """
    catname = param.String(doc='name of the emission category (should be unique to that tracer)')
    fieldname = param.Selector(doc='name of the field to be used', label='Field name')
    filename = param.Selector(doc='name of the emission file', label='File prefix')
    regions = param.ListSelector(doc='region(s) where the emissions should be applied')
    scf = param.Number(doc='scaling factor for the emissions', default=1)
    path = param.Path(doc='location of the emission files')
    fileinfo = param.String(doc='ncdump of the file (for now ...)')

    def __init__(self, regnames, **params):
        super().__init__(**params)
        self.param.regions.objects = regnames
        self.regions = regnames
        #self.update_field_choices()
        self.widgets = {
            'field': pn.widgets.Select.from_param(self.param.fieldname),
            'info': pn.pane.Markdown()
        }
        self.update_file_choices()
    
    def __panel__(self):
        return pn.Column(
            pn.Row(
                pn.widgets.TextInput.from_param(self.param.catname),
                pn.Column(
                    pn.widgets.Select.from_param(self.param.filename),
                    self.widgets['field'],
                ),
                pn.widgets.MultiChoice.from_param(self.param.regions)
            ),
            self.widgets['info']
            #pn.widgets.TextInput.from_param(self.param.scf)
        )
        
    @param.depends('filename', 'path', watch=True)
    def update_field_choices(self):
        """
        Update the choices of the "Field" widget. 
        """
        #available_files = list(Path(self.path).glob(f'**/{self.filename}*.nc*'))
        available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
        ds = xr.open_dataset(available_files[0])
        self.param.fieldname.objects = [ _ for _ in ds.data_vars if _!='area' ]
        self.fieldname = self.param.fieldname.objects[0]
        self.widgets['field'].visible = len(self.param.fieldname.objects) > 1

    @param.depends('path', watch=True)
    def update_file_choices(self):
        available_files = get_emis_file_list(self.path, '**/*.nc*')
        self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
        self.filename = self.param.filename.objects[0]
        
    @param.depends('filename', 'fieldname', watch=True)
    def update_field_description(self):
        available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
        #print(available_files)
        ds = xr.open_dataset(available_files[0])
        self.widgets['info'].object = f"""
        {ds.attrs.get('description', '`file description missing`')}
        
        **{self.fieldname}**
        - *long_name*\t: {ds[self.fieldname].long_name}
        - *units*\t: {ds[self.fieldname].units}'
        """