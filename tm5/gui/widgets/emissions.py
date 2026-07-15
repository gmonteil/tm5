#!/usr/bin/env python

import panel as pn
import param
from pathlib import Path
from functools import lru_cache
from typing import List
import xarray as xr
from loguru import logger
from omegaconf import DictConfig
from tm5.gui.css import *

supported_domain_list = ['glb600x400','eur300x200','gns100x100']

@lru_cache
def get_emis_file_list(path: Path, pattern: str) -> List[Path]:
    return list(Path(path).glob(pattern))

@lru_cache
def load_emis(pattern: str) -> xr.Dataset:
    return xr.open_mfdataset(pattern, concat_dim='time', combine='nested')

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
            field=pn.widgets.Select.from_param(self.param.fieldname),
            # info=pn.pane.Markdown(width=300,
            #                       stylesheets=[setup_stylesheet,], css_classes=['setup-tracer']),
            # title=pn.pane.Markdown(width=300,
            #                        stylesheets=[setup_stylesheet,], css_classes=['setup-tracer'])
            info=pn.pane.Markdown(width=300),
            title=pn.pane.Markdown(width=300),
        )
        self.update_desc()

    def __panel__(self):
        return pn.Column(
            self.widgets['title'],
            pn.widgets.Select.from_param(self.param.filename),
            self.widgets['field'],
            self.widgets['info'],
            stylesheets=[setup_stylesheet,], css_classes=['setup-tracer']
        )

    @param.depends('filename', 'path', 'domain', watch=True)
    def update_field_choices(self):
        """
        Update the choices of the "Field" widget.
        """
        available_files = get_emis_file_list(Path(self.path) / self.domain, f'{self.filename}*.nc')
        if len(available_files) > 0:
            ds = xr.open_dataset(available_files[0])
            self.param.fieldname.objects = [_ for _ in ds.data_vars if _ != 'area']
            self.fieldname = self.param.fieldname.objects[0]
            self.widgets['field'].visible = len(self.param.fieldname.objects) > 1

    @param.depends('path', 'domain', watch=True)
    def update_file_choices(self):
        # available_files = get_emis_file_list(self.path, '**/*.nc*')
        # -- 2025-04-14:: restrict here to the global (default) domain
        logger.info(self.domain)
        logger.info(f'**/*{self.domain}*.nc')
        available_files = get_emis_file_list(Path(self.path) / self.domain, '*.nc')
        self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
        # if len(available_files) > 0:
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
    # emis_reg = FieldSelector(desc='Emissions for the regional domain')
    # emis_glo = FieldSelector(desc='Global emissions')
    switch_reg = param.Boolean(doc="Switch alternate source for regional emissions")

    def __init__(self, **params):
        super().__init__(**params)
        self.emis_reg = FieldSelector(desc='Emissions for the regional domain', domain=self.regions[-1])
        self.emis_glo = FieldSelector(desc='Global emissions', domain=self.regions[0])
        self.emis_glo.path = self.path
        self.emis_reg.path = self.path
        self.pane_glo = pn.Column(self.emis_glo, stylesheets=[setup_stylesheet,], css_classes=['setup-tracer'])
        self.pane_reg = pn.Column(self.emis_reg, visible=len(self.regions) > 1)
        self.switch_button = pn.Row(
            pn.widgets.Switch.from_param(self.param.switch_reg, align='center'),
            pn.pane.Markdown("Use different regional emissions", stylesheets=[setup_stylesheet,], css_classes=['setup-tracer']),
            visible=len(self.regions) > 1
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
        else:
            self.pane_reg.visible = False

    @param.depends('regions', watch=True)
    def update_switch_visibility(self):
        self.switch_button.visible = len(self.regions) > 1

    def copy(self):
        newem = self.__class__(
            catname=self.catname,
            regions=self.regions,
            path=self.path)
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

#     @param.depends('path', watch=True)
#     def update_file_choices(self, domain='glb1x1'):
#         # available_files = get_emis_file_list(self.path, '**/*.nc*')
#         #-- 2025-04-14:: restrict here to the global (default) domain
#         fptn = f"**/*{domain}*.nc"
#         available_files = get_emis_file_list(self.path, fptn)
#         self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
#         self.filename = self.param.filename.objects[0]

#     @param.depends('filename', 'fieldname', watch=True)
#     def update_field_description(self):
#         available_files = get_emis_file_list(self.path, f'**/{self.filename}*.nc*')
#         #print(available_files)
#         ds = xr.open_dataset(available_files[0])
#         self.widgets['info'].object = f"""
#         {ds.attrs.get('description', '`file description missing`')}

#         **{self.fieldname}**
#         - *long_name*\t: {ds[self.fieldname].long_name}
#         - *units*\t: {ds[self.fieldname].units}'
#         """


class EmissionExplorer(pn.viewable.Viewer):
    filename  = param.Selector(doc="name of emission file")
    fieldname = param.Selector(doc="name of the field to be used")
    # category  = param.Selector()
    # region    = param.Selector()
    # domain = param.String(doc="title of the section")

    def __init__(self, settings : DictConfig, load_parallel : bool = True):
        super().__init__()

        self.settings = settings
        #
        #-- top-level emission directory
        #
        self.emisdir = self.settings.emissions.path
        logger.debug(f"...self.emisdir -->{self.emisdir}<--")
        #-- basic consistency tests
        if not Path(self.emisdir).is_dir():
            msg = f"directory for emissions ==>{self.emisdir}<== is not accessible"
            raise RuntimeError(msg)
        else:
            msg = f"top-level emissions directory ***{self.emisdir}***"
            logger.debug(msg)
        self.path = self.emisdir
        self.domain = 'glb6x4' #--MVO-TODO::domain hard-coded
        
        self.widgets = dict(
            field=pn.widgets.Select.from_param(self.param.fieldname),
            # info=pn.pane.Markdown(width=300,
            #                       stylesheets=[setup_stylesheet,], css_classes=['setup-tracer']),
            # title=pn.pane.Markdown(width=300,
            #                        stylesheets=[setup_stylesheet,], css_classes=['setup-tracer'])
            info=pn.pane.Markdown(width=300),
            title=pn.pane.Markdown(width=300),
        )
        self.update_desc()
        self.update_file_choices()
        self.update_field_choices()
        # # #
        # # #-- read and prepare emissions to be ready for plotting
        # # #
        # # if load_parallel:
        # #     self.setup_emis_mp() #-- multiprocessor version to speed-up
        # # else:
        # #     self.setup_emis()

        # # #
        # # #--
        # # #
        # # self.param.region.objects = list(self.data.keys())
        # # self.region = self.param.region.objects[-1]
        # # self.region = self.param.region.objects[1]

        # # datavar_list = sorted(list(self.data[self.region].data_vars))
        # # self.param.category.objects = datavar_list
        # # #
        # # #-- start with 'fossil' as initial category (if present),
        # # #   otherwise start with first in list
        # # #
        # # istart = 0
        # # for i,varname in enumerate(datavar_list):
        # #     if varname=='fossil':
        # #         istart = i
        # #         break
        # # self.category = self.param.category.objects[istart]
      
        # self.widgets = {
        #     'field_selector': pn.widgets.Select.from_param(self.param.category),
        #     'region_selector': pn.widgets.Select.from_param(self.param.region)
        # }
        
    def __panel__(self):
        widgets = pn.Column(
            self.widgets['title'],
            pn.widgets.Select.from_param(self.param.filename),
            self.widgets['field'],
            self.widgets['info'],
            stylesheets=[setup_stylesheet,], css_classes=['setup-tracer']
        )
        return widgets
        # widget_regsel = self.widgets['region_selector']
        # widget_fldsel = self.widgets['field_selector']
        # widget_map1 = hv.DynamicMap(self.map_emis)
        # widget_map2 = hv.DynamicMap(self.plot_domain_emis)
        # widget = pn.Column(widget_regsel,widget_fldsel,widget_datesel,
        #                    pn.Row(widget_map1,widget_map2))
        # return widget

    # @param.depends('filename', 'path', 'domain', watch=True)
    @param.depends('filename', watch=True)
    def update_field_choices(self):
        """
        Update the choices of the "Field" widget.
        """
        available_files = get_emis_file_list(Path(self.path) / self.domain, f'{self.filename}*.nc')
        if len(available_files) > 0:
            ds = xr.open_dataset(available_files[0])
            self.param.fieldname.objects = [_ for _ in ds.data_vars if _ != 'area']
            self.fieldname = self.param.fieldname.objects[0]
            self.widgets['field'].visible = len(self.param.fieldname.objects) > 1
    # @param.depends('path', 'domain', watch=True)
    # @param.depends('domain', watch=True)
    def update_file_choices(self):
        available_files = get_emis_file_list(Path(self.path) / self.domain, '*.nc')
        logger.debug(f"available_files ***{available_files}***")

        self.param.filename.objects = set([f.name.rsplit('_', maxsplit=1)[0] for f in available_files])
        # if len(available_files) > 0:
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

    # @param.depends('desc', watch=True)
    # def update_desc(self):
    #     self.widgets['title'].object = f'### {self.desc}'

    def update_desc(self):
        self.widgets['title'].object = f'### MY DESCRIPTOR'

    # @property
    # @pn.depends('region', watch=True)
    def current_extent(self):
        if self.region in ['glb600x400','glb6x4',]:
            return [-180, 180, -90, 90]
        elif self.region in ['eur300x200', 'eur3x2',]:
            return [-36, 54, 22, 74]
        elif self.region in ['gns100x100','gns1x1',]:
            return [0, 18, 42, 58]
        else:
            return RuntimeError(f"unexpected region -->{self.region}<--")
    
    # @pn.depends('region', 'category')
    # def map_emis(self):
    #     import cartopy.crs as ccrs
    #     #-- MVO::can end-up here with self.category==None,
    #     #        so need to catch this case
    #     msg = f"{dtm.datetime.utcnow()}, start@map_emis -->{self.region}<-- -->{self.category}<-- itime={self.time_index}"
    #     logger.debug(msg)
    #     # print(f"DEBUG::{msg}")
    #     cur_cat = self.category if self.category!=None else 'wetland'
    #     cur_date = pd.to_datetime(str(self.current_date)).strftime('%Y-%m-%d')
    #     lonmin,lonmax,latmin,latmax = self.current_extent()
    #     # print(f"DEBUG::@map_emis title -->{title}<--")
    #     # features features (default=None): A list of features or a dictionary of features and the scale at which to render it. Available features include ‘borders’, ‘coastline’, ‘lakes’, ‘land’, ‘ocean’, ‘rivers’ and ‘states’. Available scales include ‘10m’/’50m’/’110m’.
    #     if self.region=='glb600x400':
    #         cproj = ccrs.PlateCarree()
    #         cfeatures = {'borders':'110m', 'coastline':'110m'}
    #         coastline = '110m'
    #     elif self.region=='eur300x200':
    #         cproj = ccrs.GOOGLE_MERCATOR
    #         cproj = ccrs.PlateCarree()
    #         cfeatures = {'borders':'50m', 'coastline':'50m'}
    #         coastline = '50m'
    #     elif self.region=='gns100x100':
    #         cproj = ccrs.GOOGLE_MERCATOR
    #         cproj = ccrs.PlateCarree()
    #         cfeatures = {'borders':'10m', 'coastline':'10m'}
    #         coastline = '10m'
    #     cur_emis = self.data[self.region][cur_cat].isel(time=self.time_index)
    #     unitlabel = f"[kg{species}/cell/s]"
    #     #
    #     #-- compute total
    #     #
    #     emis_tot = cur_emis.sum(('lat','lon')).values #- kgTRACER/s
    #     emis_tot = emis_tot / 1e6 * (365*86400)       #- ktTRACER/year
    #     emis_tot_unit = f"[kt{species}/year]"
    #     if emis_tot>=1.e3:
    #         emis_tot = emis_tot / 1e3
    #         emis_tot_unit = f"[Mt{species}/year]"
    #     msg = f"{dtm.datetime.utcnow()}, @map_emis, current data ready"
    #     logger.debug(msg)
    #     title = f"{cur_cat}@{self.region} " \
    #         f"({cur_date}, total={emis_tot:.3f}{emis_tot_unit})"
    #     title = f"{cur_cat}@{self.region} {unitlabel} ({cur_date})" + '\n' \
    #         f"domain total: {emis_tot:.3f}{emis_tot_unit}"
    #     # with open(f"map_emis_{dtm.datetime.utcnow().isoformat()}.log", 'w') as fp:
    #     #     msg = f"***{title}*** coastline={coastline} cfeatures={cfeatures}"
    #     #     fp.write(f"{msg}" +'\n')
    #     # emis_hvplot = cur_emis.hvplot.quadmesh(geo=True,
    #     #                                        coastline=coastline,
    #     #                                        features=cfeatures,
    #     #                                        xlim=(lonmin,lonmax),
    #     #                                        ylim=(latmin,latmax),
    #     #                                        clabel=unitlabel, title=title)
    #     emis_hvplot = cur_emis.hvplot.quadmesh(
    #         xlim=(lonmin,lonmax), ylim=(latmin,latmax),
    #         coastline=coastline, title=title,
    #         crs=ccrs.PlateCarree(),
    #         projection=cproj, features=cfeatures, project=True)#, rasterize=True)
    #     msg = f"{dtm.datetime.utcnow()}, @map_emis, hvplot ready"
    #     logger.debug(msg)
    #     #
    #     #-- autohide toolbar
    #     #
    #     emis_hvplot = emis_hvplot.opts(
    #         # colorbar_options={'label':unitlabel},
    #         backend_opts={"plot.toolbar.autohide": True})
    #     # #-- 'framewise=True' did not help to update maps
    #     # #   when changing region...
    #     # emis_hvplot = emis_hvplot.opts(framewise=True, backend_opts={"plot.toolbar.autohide": True})
    #     return emis_hvplot

    # # @pn.depends('time_index', 'region', 'category', watch=True)
    # @pn.depends('time_index', 'region', 'category')
    # def plot_domain_emis(self):
    #     msg = f"{dtm.datetime.utcnow()}, start@plot_domain_emis with time_index={self.time_index}, " \
    #         f"region -->{self.region}<--  category -->{self.category}<--"
    #     logger.debug(msg)
    #     # print(f"DEBUG::{msg}")
    #     if self.region==None:
    #         return
    #     cat_df = self.glob_timeseries[self.region]
    #     ylabel = f"[kg{species}/s]"
    #     if self.region=='glb600x400':
    #         title = f"global sectoral {species} emission totals (@{self.region})"
    #     else:
    #         title = f"domain sectoral {species} emission totals (@{self.region})"
    #     # print(cat_df.head())
    #     # cat_hvplot = cat_df.hvplot(grid=True, x='time',
    #     #                            xlabel='time', ylabel=ylabel, title=title,
    #     #                            fontsize={'legend':6})
    #     cat_hvplot = cat_df.hvplot(grid=True, x='time',
    #                                xlabel='time', ylabel=ylabel, title=title)
    #     # cat_hvplot.opts(backend_opts={"plot.toolbar.autohide": True}, legend_position='bottom_right', legend_offset=(0,0), fontsize={'legend':6,'legend_title':6})
    #     cat_hvplot.opts(backend_opts={"plot.toolbar.autohide": True}, fontsize={'legend':6,'legend_title':6})
    #     # cat_hvplot = cat_df.plot(grid=True, x='time',
    #     #                          xlabel='time', ylabel=ylabel, title=title)
    #     msg = f"{dtm.datetime.utcnow()}, @plot_domain_emis, hvplot prepared"
    #     logger.debug(msg)
    #     # print(f"DEBUG::{msg}")
    #     # print(f"DEBUG::type(cat_hvplot)={type(cat_hvplot)}")
    #     return cat_hvplot


    # @lru_cache
    # def setup_emis(self):
    #     msg = f"{dtm.datetime.utcnow()}, @setup_emis, start"
    #     logger.debug(msg)
    #     self.data = OrderedDict()
    #     self.glob_timeseries = OrderedDict()
    #     self.dates = None
    #     for reg in supported_domain_list:
    #         msg = f"...@setup_emis@{reg} reading input"
    #         logger.debug(msg)
    #         #-- check whether emissions for region were generated
    #         fptn = f"ch4emis.{species}.{reg}.{self.yyyymmdd_ptn}.nc"
    #         _file_lst = sorted(list(self.emisdir.glob(fptn)))
    #         if len(_file_lst)==0:
    #             msg = f"no emissions files detected with pattern " \
    #                 f"==>{fptn}<=="
    #             raise RuntimeError(msg)
    #         else:
    #             emis_ptn = f"{str(self.emisdir)}/{fptn}"
    #             cur_emis = load_emis(emis_ptn)
    #             cur_dates = {Timestamp(v).strftime('%B %Y'): iv for (iv, v) in enumerate(cur_emis.time.values)}
    #             if self.dates is None:
    #                 self.dates = cur_dates
    #             else:
    #                 #-- make sure all regions are for the same dates
    #                 assert np.all(self.dates==cur_dates)
    #             #
    #             #-- add cat emissions field
    #             #
    #             self.data[reg] = cur_emis
    #             #
    #             #-- caching regional-sum time-series (per category),
    #             #   otherwise we have long delay when visualising...
    #             #
    #             msg = f"...@setup_emis@{reg}, prepare overall time-series"
    #             logger.debug(msg)
    #             cur_glob = cur_emis.sum(('lat','lon')).to_dataframe()
    #             cur_glob = cur_glob.reset_index()
    #             self.glob_timeseries[reg] = cur_glob
    #             msg = f"...@setup_emis@{reg} done"
    #             logger.debug(msg)
    #     # print(f"MVODEBUG:setup_emis: list(self.data.keys()) ***{list(self.data.keys())}***")
    #     msg = f"@setup_emis, finished"
    #     logger.debug(msg)

    # @lru_cache
    # def setup_emis_mp(self):
    #     from multiprocessing import Process, Manager
    #     def do_region(reg):
    #         msg = f"...@setup_emis@{reg} reading input"
    #         logger.debug(msg)
    #         #-- check whether emissions for region were generated
    #         _file_lst = sorted(list(self.emisdir.glob(f"ch4emis.{species}.{reg}.{self.yyyymmdd_ptn}.nc")))
    #         if len(_file_lst)>0:
    #             emis_ptn = f"{str(self.emisdir)}/ch4emis.{species}.{reg}.{self.yyyymmdd_ptn}.nc"
    #             cur_emis = load_emis(emis_ptn)
    #             cur_dates = {Timestamp(v).strftime('%B %Y'): iv for (iv, v) in enumerate(cur_emis.time.values)}
    #             mp_dates[reg] = cur_dates
    #             #
    #             #-- add cat emissions field
    #             #
    #             mp_data[reg] = cur_emis
    #             # self.data[reg] = cur_emis
    #             #
    #             #-- caching regional-sum time-series (per category),
    #             #   otherwise we have long delay when visualising...
    #             #
    #             msg = f"...@setup_emis@{reg}, prepare overall time-series"
    #             logger.debug(msg)
    #             cur_glob = cur_emis.sum(('lat','lon')).to_dataframe()
    #             cur_glob = cur_glob.reset_index()
    #             # self.glob_timeseries[reg] = cur_glob
    #             mp_domtseries[reg] = cur_glob
    #             msg = f"...@setup_emis@{reg} done"
    #             logger.debug(msg)

    #     msg = f"{dtm.datetime.utcnow()}, @setup_emis, start"
    #     logger.debug(msg)
    #     # print(f"DEBUG::{msg}")
    #     manager = Manager()
    #     mp_dates = manager.dict()
    #     mp_data  = manager.dict()
    #     mp_domtseries = manager.dict()
    #     processes = [Process(target=do_region, args=(_reg,)) for _reg in supported_domain_list ]
    #     for process in processes:
    #         process.start()
    #     #-- wait for all process to complete
    #     for process in processes:
    #         process.join()
    #     # print(f"...processes joined.")
    #     self.data = OrderedDict()
    #     self.glob_timeseries = OrderedDict()
    #     for _dom in supported_domain_list:
    #         self.data[_dom] = mp_data[_dom]
    #         self.glob_timeseries[_dom] = mp_domtseries[_dom]
    #     #-- TODO: should still make sure dates are equal for all regions
    #     self.dates = mp_dates['glb600x400']
