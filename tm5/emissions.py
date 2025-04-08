#!/usr/bin/env python

"""
Prepare daily emissions from the requested grid
"""
import sys
from omegaconf import DictConfig
from tm5.gridtools import TM5Grids
import xarray as xr
import numpy as np
from pandas import date_range, Timestamp
from pathlib import Path
from loguru import logger
from tqdm.auto import tqdm


def coarsen_file(path_or_pattern : str, reg : TM5Grids, start : Timestamp, end : Timestamp, use_esmf_regridder : bool = False) -> xr.Dataset:
    # Open the source file(s). It should be on a global 1x1 resolution (for now ...)
    ds = xr.open_mfdataset(path_or_pattern)

    #-- expected grid of global emission inputs
    # !!!TODO:should add consistency check
    glb1x1 = TM5Grids.global1x1()

    #
    #-- MVO 2025-03-14::
    #   - re-gridding results seem not appropriate
    #     with ESMF regridder used underneath,
    #     by default apply again approach as used in old implementation
    #
    # Regrid
    if use_esmf_regridder:
        # Define a regridder
        regridder = glb1x1.regridder(reg)
        coarse = regridder(ds)
    else:
        #
        #-- iterate over (emission) datasets and
        #   (1) crop to target region
        #   (2) re-grid to target resolution
        data_coarse_dict = {}
        for da_name in ds.data_vars:
            if da_name=='area':
                continue
            da = ds[da_name]
            #-- make sure units is as expected
            if not da.units=='kg/cl/s':
                msg = f"emission {da_name} with unexpected units ==>{da.units}<=="
                raise RuntimeError(msg)
            #-- (1) crop (result is numpy array)
            data_out = da.sel(lat=slice(reg.south, reg.north),
                              lon=slice(reg.west, reg.east)).values
            #-- (2) coarsen (only if needed):
            if (reg.dlat!=glb1x1.dlat) or (reg.dlon!=glb1x1.dlon):
                assert (reg.dlat - int(reg.dlat) == 0) & (reg.dlon - int(reg.dlon) == 0), "lat and lon steps must be multiple of the base (1°) step"
                nt, nlat, nlon = data_out.shape
                data_out = data_out.reshape(nt, -1, int(reg.dlat), nlon).sum(2)
                nlat = data_out.shape[1]
                data_out = data_out.reshape(nt, nlat, -1, int(reg.dlon)).sum(3)
            data_coarse_dict[da_name] = xr.DataArray(
                data_out,
                dims = ('time', 'lat', 'lon'),
                coords = {
                    'time': ds.time,
                    'lat': reg.latc,
                    'lon': reg.lonc
                }
            )
        #
        #-- create xarray dataset
        #
        coarse = xr.Dataset(
            coords=dict(lon=('lon', reg.lonc),
                        lat=('lat', reg.latc),
                        time=ds['time']),
            data_vars=data_coarse_dict
            )

    #
    # ensure result is on a daily time step (at least for  now ...)
    #
    coarse = coarse.reindex(time=date_range(start, end, freq='D')).ffill('time')
    #
    # ensure float64 precision (as TM5 I/O expects)
    #
    coarse = coarse.astype('f8')

    return coarse


def coarsen_file_with_overwrite( cat : DictConfig, reg : TM5Grids, start : Timestamp, end : Timestamp, use_esmf_regridder : bool = False) -> xr.Dataset:
    #
    #-- default/main emissons filed
    #
    path_or_pattern = cat.path
    field           = cat.field
    
    # Open the source file(s). It should be on a global 1x1 resolution (for now ...)
    ds = xr.open_mfdataset(path_or_pattern)
    #-- can restrict to selected field (might speedup regridding)
    if not field in ds.data_vars:
        msg = f"requested field ==>{field}<== not found in ***{path_or_pattern}***"
        raise RuntimeError(msg)
    elif ds[field].attrs['units']!='kg/cl/s':
        msg = "field={field} with unexpected unit -->{ds[field].attrs['units']}<--"
        raise RuntimeError(msg)
    else:
        ds = ds.drop_vars([_ for _ in ds.data_vars if _!=field])
    #
    #-- confirm that default emissions are
    #   - global
    #   - have 1 degree resolution
    #
    lonres = ds.attrs['geospatial_lon_resolution']
    latres = ds.attrs['geospatial_lat_resolution']
    if lonres!=1 or latres!=1:
        msg = f"input emissions must be prepared at regular grid with 1 degree resolution, " \
            f"detected lonres/latres = {lonres}/{latres}"
        raise RuntimeError(msg)
    lonmin = ds.attrs['geospatial_lon_min']
    lonmax = ds.attrs['geospatial_lon_max'] 
    latmin = ds.attrs['geospatial_lat_min']
    latmax = ds.attrs['geospatial_lat_max']
    if lonmin!=-180 or lonmax!=180 or latmin!=-90 or latmax!=90:
        msg = f"default emissions expected to have global extent, " \
            f"detected lonmin/lonmax = {lonmin}/{lonmax}, " \
            f"latmin/latmax = {latmin}/{latmax}"
        raise RuntimeError(msg)
    #
    #-- data file for overwriting part of emissions
    #
    ovr_path_or_pattern = cat.overwrite.path
    ovr_field           = cat.overwrite.field
    ovr_ds = xr.open_mfdataset(ovr_path_or_pattern)
    if not ovr_field in ovr_ds.data_vars:
        msg = f"requested field ==>{ovr_field}<== not found in ***{ovr_path_or_pattern}***"
        raise RuntimeError(msg)
    elif ovr_ds[field].attrs['units']!='kg/cl/s':
        msg = "field={ovr_field} with unexpected unit -->{ovr_ds[field].attrs['units']}<--"
        raise RuntimeError(msg)
    else:
        ovr_ds = ovr_ds.drop_vars([_ for _ in ovr_ds.data_vars if _!=ovr_field])
    ovr_lonres = ovr_ds.attrs['geospatial_lon_resolution']
    ovr_latres = ovr_ds.attrs['geospatial_lat_resolution']
    if ovr_lonres!=1 or ovr_latres!=1:
        msg = f"input emissions must be prepared at regular grid with 1 degree resolution, " \
            f"detected lonres/latres = {ovr_lonres}/{ovr_latres}"
        raise RuntimeError(msg)

    #
    #-- expected grid of global emission inputs
    # !!!TODO:should add consistency check
    glb1x1 = TM5Grids.global1x1()

    #
    #-- MVO 2025-03-14::
    #   - re-gridding results seem not appropriate
    #     with ESMF regridder used underneath,
    #     by default apply again approach as used in old implementation
    #
    # Regrid
    if use_esmf_regridder:
        msg = f"ESMF regridder not supported when overriding emissions."
        raise NotImplementedError(msg)
        # # Define a regridder
        # regridder = glb1x1.regridder(reg)
        # coarse = regridder(ds)
    else:
        #-- default emission data
        da = ds[field]
        if not da.dims==('time','lat','lon'):
            msg = f"default emis variable -->{field}<-- with unexpected dimensions " \
                f"***{da.dims}***"
            raise RuntimeError(msg)
        nt,nlat,nlon = da.shape
        #-- global default values
        da_data = ds[field].values
        #
        #-- (0) overwrite (spatial part of default)
        #
        ovr_da = ovr_ds[ovr_field]
        if not ovr_da.dims==('time','lat','lon'):
            msg = f"overwrite emis variable -->{ovr_field}<-- with unexpected dimensions " \
                f"***{da.dims}***"
            raise RuntimeError(msg)
        ovr_data = ovr_da.values
        # print(f"da.shape -->{da.shape}<--, ovr_da.shape -->{ovr_da.shape}<--")
        # print(f"MVODEBUG:: da.time     = ***{da.time.values}***")
        # print(f"MVODEBUG:: ovr_da.time = ***{ovr_da.time.values}***")
        ovr_dict = dict(lat=ovr_da.lat,lon=ovr_da.lon)
        #
        #
        #
        #-- detect spatial index range of overwrite-data in global field
        #
        ilatmin = np.where(da.lat.values==ovr_da.lat.values[0])[0]
        ilatmax = np.where(da.lat.values==ovr_da.lat.values[-1])[0]
        assert len(ilatmin)==1 and len(ilatmax)==1
        ilatmin = ilatmin[0]
        ilatmax = ilatmax[0]
        # print(f"ilatmin/ilatmax = {ilatmin}/{ilatmax}")
        ilonmin = np.where(da.lon.values==ovr_da.lon.values[0])[0]
        ilonmax = np.where(da.lon.values==ovr_da.lon.values[-1])[0]
        assert len(ilonmin)==1 and len(ilonmax)==1
        ilonmin = ilonmin[0]
        ilonmax = ilonmax[0]
        #
        #-- basic consistency check of temporal dimensions,
        #   two cases supported
        #   - default/overwrite fluxes are for identical times
        #   - default is monthly, overwrite is annual
        #     (then all months get the same annual flux rate in
        #      the spatial region which is overwritten)
        #
        if ds.attrs['time_coverage_resolution']==ovr_ds.attrs['time_coverage_resolution']:
            def_tvalues = da.time.values
            ovr_tvalues = ovr_da.time.values
            if not np.all(def_tvalues==ovr_tvalues):
                msg = "default/overwrite emission fluxes not prepared " \
                    f"for same points in time."
                raise RuntimeError(msg)
            #
            #-- can savely overwrite in spatial domain
            #
            da_data[:,ilatmin:ilatmax+1,ilonmin:ilonmax+1] = ovr_data[:]
        elif ds.attrs['time_coverage_resolution']=='P1M' and \
             ovr_ds.attrs['time_coverage_resolution']=='P1Y':
            def_tvalues = da.time.values
            ovr_tvalues = ovr_da.time.values
            if not def_tvalues[0]==ovr_tvalues[0]:
                msg = "default/overwrite emission fluxes not prepared " \
                    f"for same points in time."
                raise RuntimeError(msg)
            #
            #-- can savely overwrite in spatial domain,
            #   each month gets same emission rate
            #
            da_data[:,ilatmin:ilatmax+1,ilonmin:ilonmax+1] = ovr_data[np.newaxis,:,:]
        else:
            msg = f"fluxes with temporal resolution for " \
                f"default={ds.attrs['time_coverage_resolution']} and " \
                f"overwrite={ ovr_ds.attrs['time_coverage_resolution']} " \
                f"is not yet supported."
            raise NotImplementedError(msg)
        #
        #-- turn into xr.DataArray
        #
        da_out = xr.DataArray(
            da_data,
            dims = ('time','lat','lon'),
            coords = {'time':ds.time,
                      'lon': ds.lon,
                      'lat':ds.lat}
        )

        # #
        # #-- global default and overwrite with same temporal settings
        # #   NOTE: these initial trials below yielded
        # #         errors with dask:
        # #   "NotImplementedError: xarray can't set arrays with multiple array indices to dask yet"
        # #
        # if ds.attrs['time_coverage_resolution']==ovr_ds.attrs['time_coverage_resolution']:
        #     #-- overwrite in spatial domain
        #     da.loc[ovr_dict] = ovr_da
        # elif ds.attrs['time_coverage_resolution']=='P1M' and \
        #      ovr_ds.attrs['time_coverage_resolution']=='P1Y':
        #     # da.loc[ovr_dict] = ovr_da.values #ovr_da.isel(time=0)
        #     # for t in da.time.values:
        #     #      tovr_dict = dict(time=t,**ovr_dict)
        #     #      da.loc[tovr_dict] = ovr_da.values[0,:]
        #
        #-- (1) crop (result is numpy array)
        #       Note:: since input emissions are global,
        #              cropping will always work
        da_crop = da_out.sel(lat=slice(reg.south, reg.north),
                             lon=slice(reg.west, reg.east))
        data_out = da_crop.values

        #-- (2) coarsen (only if needed):
        if (reg.dlat!=glb1x1.dlat) or (reg.dlon!=glb1x1.dlon):
            assert (reg.dlat - int(reg.dlat) == 0) & (reg.dlon - int(reg.dlon) == 0), "lat and lon steps must be multiple of the base (1°) step"
            nt, nlat, nlon = data_out.shape
            data_out = data_out.reshape(nt, -1, int(reg.dlat), nlon).sum(2)
            nlat = data_out.shape[1]
            data_out = data_out.reshape(nt, nlat, -1, int(reg.dlon)).sum(3)

        #
        #-- create xarray dataset
        #
        data_coarse_dict = {}
        data_coarse_dict[field] = xr.DataArray(
            data_out,
            dims = ('time', 'lat', 'lon'),
            coords = {
                'time': ds.time,
                'lat': reg.latc,
                'lon': reg.lonc
            }
        )
        coarse = xr.Dataset(
            coords=dict(lon=('lon', reg.lonc),
                        lat=('lat', reg.latc),
                        time=ds['time']),
            data_vars=data_coarse_dict
            )

    #
    # ensure result is on a daily time step (at least for  now ...)
    #
    coarse = coarse.reindex(time=date_range(start, end, freq='D')).ffill('time')
    #
    # ensure float64 precision (as TM5 I/O expects)
    #
    coarse = coarse.astype('f8')

    return coarse


def prepare_emissions(conf: DictConfig) -> None:
    
    start = Timestamp(conf.run.start)
    end   = Timestamp(conf.run.end)
    
    for region in conf.run.regions:
        
        # Create the region object:
        reg = conf.regions[region]
        reg = TM5Grids.from_corners(west=reg.lons[0], east=reg.lons[1], south=reg.lats[0], north=reg.lats[1], dlon=reg.lons[2], dlat=reg.lats[2])
        
        # Regrid the emissions for that region and tracer
        for trname, tracer in conf.emissions.items():

            if 'categories' not in tracer:
                continue

            #-- collect
            datasets = {}
            
            # Do the actual coarsening.
            for catname, cat in tracer.categories.items():
                if 'overwrite' in cat:
                    print(f"@{region}, {catname}: ***{cat}***")
                    ds = coarsen_file_with_overwrite(cat, reg, start, end)
                    datasets[catname] = ds[cat.field]
                else:
                    datasets[catname] = coarsen_file(cat.path, reg, start, end)[cat.field]
                
            # Ensure that the dest path exists!
            Path(tracer.prefix).parent.mkdir(exist_ok=True, parents=True)

            #-- probably better raise warning here
            if len(datasets)==0:
                msg = f"no emissions to be generated for {trname}@{region}"
                logger.warning(msg)
                continue
            
            datasets = xr.Dataset(datasets)

            # group the datasets in daily emission files for that region and tracer:
            for day in tqdm(date_range(start, end, freq='D'), desc=f'Writing emission files for {trname} in region {region}'):
                destfile = f'{tracer.prefix}.{trname}.{region}.{day:%Y%m%d}.nc'
                datasets.sel(time=day).to_netcdf(destfile)
