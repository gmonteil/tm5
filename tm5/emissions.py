#!/usr/bin/env python

import sys
import xarray as xr
import numpy as np
from pandas import Timestamp, date_range, DatetimeIndex
from pandas.tseries.frequencies import to_offset
from omegaconf import DictConfig
from tm5.gridtools import TM5Grids
from typing import Tuple, Dict, Union
from netCDF4 import Dataset
from loguru import logger
from pathlib import Path
from tm5.units import units_registry as ureg
from dateutil.relativedelta import relativedelta


def crop_and_coarsen(glo1x1: xr.DataArray, dest: TM5Grids) -> xr.DataArray:
    """
    Crop and coarsen a global 1x1 emission field to a subregion at equal or lower resolution.
    """

    # First, crop:
    if 'lat' in glo1x1.dims and 'lon' in glo1x1.dims:
        data = glo1x1.sel(
            lat=slice(dest.south, dest.north),
            lon=slice(dest.west, dest.east)).values
    elif 'latitude' in glo1x1.dims and 'longitude' in glo1x1.dims:
        data = glo1x1.sel(
            latitude=slice(dest.south, dest.north),
            longitude=slice(dest.west, dest.east)).values
    else:
        msg = f"spatial coordinates could not be detected, these are the dimensions -->{glo1x1.dims}<--"
        raise RuntimeError(msg)
    # Then, coarsen:
    assert (dest.dlat - int(dest.dlat) == 0) & (dest.dlon - int(dest.dlon) == 0), "lat and lon steps must be multiple of the base (1°) step"
    nt, nlat, nlon = data.shape
    data = data.reshape(nt, -1, int(dest.dlat), nlon).sum(2)
    nlat = data.shape[1]
    data = data.reshape(nt, nlat, -1, int(dest.dlon)).sum(3)

    return xr.DataArray(
        data,
        dims = ('time', 'latitude', 'longitude'),
        coords = {
            'time': glo1x1.time,
            'latitude': dest.latc,
            'longitude': dest.lonc
        })


def load_preprocessed_emis(
        categ: Union[Dict, DictConfig],
        species: str,
        destreg: TM5Grids,
        period: slice,
) -> xr.DataArray:
    """
    This reads pre-processed emissions (global, 1x1˚, mol/m2/s) and convert them to emissions in kg[tracer]/s limited to the requested TM5 region.

    The emissions should be stored in (a) xarray-compatible netcdf file(s), covering the requested period.
    Alternatively, the files can be treated as a climatology. In this case, the data will be monthly (even
    if the data used to construct the climatology has a higher temporal resolution.

    Arguments :
        - categ : dictionary with one mandatory and three [optional] fields:
            * pattern: glob pattern of the files from which emissions are read
            * [field]: name of the netcdf variable containing the emissions (default: "emis")
            * [area_field]: name of the netcdf variable containing the grid cell area (if not provided, it
                            will be computed
            * [climatology]: whether to tread the data as a climatology or as regular emissions (default: False)
        - species: name of the species to be read (should be one that is implemented in tm5.units module)
        - destreg: grid specification of the TM5 region at which the data should be cropped/coarsened (should be an instance of tm5.gridtools.TM5Grids)
        - period: slice containing start and end dates
    """
    file_pattern = categ['pattern']
    field = categ.get('field', 'emis')
    area = categ.get('area_field', None)
    msg = f"reading emissions from file pattern -->{file_pattern}<--"
    logger.info(msg)
    data = xr.open_mfdataset(file_pattern)
    
    if categ.get('climatology', False):
        logger.warning(f"Treating files {file_pattern} as a climatological field")

        # 1st, calculate monthly averages:
        data = data[field].groupby(data.time.dt.month).mean()

        # Create the actual dataframe:
        emis_glo1x1 = xr.DataArray(
            dims = ('time', 'latitude', 'longitude'),
            coords = {
                'time' : date_range(
                    period.start.strftime('%Y-%m-01'),
                    (period.stop + relativedelta(months=1)).strftime('%Y-%m-01'),
                    freq='MS', inclusive='left'),
                'latitude': data.latitude.values,
                'longitude': data.longitude.values
            }
        )

        if 'pointonly' in categ:
            point_spec = categ['pointonly']
            msg = f"prepare emissions for single-point only with specification ==>{point_spec}<=="
            logger.info(msg)
            msg = f"not yet implemented for climatology setup!"
            raise NotImplementedError(msg)
        else:
            data = data.values # Load the values in memory, at this point
        for tstep in range(emis_glo1x1.values.shape[0]):
            emis_glo1x1.values[tstep, :, :] = data[emis_glo1x1.time.dt.month.values[tstep] - 1, :, :]
        emis_glo1x1 = emis_glo1x1.sel(time=(emis_glo1x1.time >= period.start) & (emis_glo1x1.time < period.stop))
    else:
        #
        # assert (data.time.min() <= period.start) & (data.time.max() >= period.stop), logger.error(
        #     f"Emissions for category {categ['name']} ({categ['pattern']}) do not cover the full simulation period: \n" +
        #     f"The emission in files range from {data.time.min(): %d %b %Y} to {data.time.max(): %d %b %Y}")

        # Time normally refers to the start of the period, but in some files it refers to the middle.
        # In this case, pass a "time_shift" key to the category config node, corresponding to the time correction to
        # be applied (e.g. "-1.5H" for 3-hourly CarbonTracker files).
        time_shift = to_offset(categ.get('time_shift', '0H'))
        # print(f"MVMV::time_shift = {time_shift}")
        # print(data.time[:10])
        # import sys; sys.exit(0)
        time = DatetimeIndex(data.time) + time_shift
        data = data.assign_coords(time=time)

        # Ensure we have an area array in memory (should be for global1x1)
        emis_glo1x1 = data[field].sel(time=(data.time >= period.start) & (data.time < period.stop))

        msg = f"after temporal selection emis_glo1x1 ==>{emis_glo1x1}<=="
        logger.info(msg)
        msg = f"...emission statistics min/mean/max = {emis_glo1x1.data.min()}/{emis_glo1x1.data.mean()}/{emis_glo1x1.data.max()}"
        logger.info(msg)
        #--
        if 'pointonly' in categ:
            point_spec = categ['pointonly']
            msg = f"prepare emissions for single-point only with specification ==>{point_spec}<=="
            logger.info(msg)
            lon_pt = point_spec['lon']
            lat_pt = point_spec['lat']
            emis_pt = point_spec['emis']
            lon_data = emis_glo1x1.longitude.values
            lat_data = emis_glo1x1.latitude.values
            ilon = np.argmin(np.abs(lon_data-lon_pt))
            ilat = np.argmin(np.abs(lat_data-lat_pt))
            msg = f"...closest data point at ilat/ilon = {ilat}/{ilon} ({lat_data[ilat]}/{lon_data[ilon]})"
            logger.info(msg)
            msg = f"...setting emissions to value {emis_pt} [mol/m2/s]"
            logger.info(msg)
            emis_glo1x1.data[:,:,:] = 0
            emis_glo1x1.data[0,ilat,ilon] = emis_pt
            msg = f"...emission statistics after point insertion: min/mean/max = {emis_glo1x1.data.min()}/{emis_glo1x1.data.mean()}/{emis_glo1x1.data.max()}"
            logger.info(msg)
    ##
    # Convert from mol/m2/s to kg[tracer]/s:
    if area is None :
        area = TM5Grids.global1x1().area
    else :
        area = data[area].values
    if emis_glo1x1.units==('mol m-2 s-1'): #-- MVO::at least correct for Carbontracker inputs
        msg = f"apply unit conversion of prepared emissions: 'mol m-2 s-1' to kg{species} cl-1 s-1"
        logger.info(msg)
        emis_glo1x1 *= area * ureg.Quantity('mol').to(f'kg{species}').m
        emis_glo1x1.attrs['units'] = f'kg{species} s-1'
    else:
        msg = f"unsupported units emis_glo1x1.units = {emis_glo1x1.units}"
        raise RuntimeError(msg)
    # Regrid:
    return crop_and_coarsen(emis_glo1x1, dest=destreg)


def split_baseline_anomalies(emreg: xr.DataArray, start: Timestamp, end: Timestamp, freq: str) -> Tuple[xr.DataArray, Dict]:
    msg = f"requested temporal settings start={start} end={end} freq={freq}"
    logger.info(msg)

    baseline = emreg.resample(time=freq).mean()

    anom = {}

    assert freq in ['D', 'MS']
    baseline['time'] = date_range(start, end, freq=freq, inclusive='left')
    # Calculate daily cycles:
    for day in date_range(start, end, freq='D', inclusive='left'):
        emis = emreg.sel(time=DatetimeIndex(emreg.time.dt.date) == day).values

        # Determine the corresponding period:
        bline = baseline.sel(time = baseline.time <= day).values[0, :, :] # 1st interval whose start time is at max equal to the current day start time

        anom[day] = xr.Dataset(
            data_vars = {'emission_anomaly': (('timesteps', 'latitude', 'longitude'), emis - bline)}
        )

    return baseline, anom
# ---end-of-split_baseline_anomalies


def split_daily_anomalies(emis: xr.DataArray, start: Timestamp, end: Timestamp) -> Dict:

    anom = {}

    # Calculate daily cycles:
    for day in date_range(start, end, freq='D', inclusive='left'):
        emis_day = emis.sel(time=DatetimeIndex(emis.time.dt.date) == day).values

        # Determine the corresponding period:
        bline = emis.sel(time = emis.time <= day).values[0, :, :] # 1st interval whose start time is at max equal to the current day start time

        anom[day] = xr.Dataset(
            data_vars = {'emission_anomaly': (('timesteps', 'latitude', 'longitude'), emis_day - bline)}
        )

    return anom
# ---end-of-split_daily_anomalies


def prepare_emissions(dconf : DictConfig, filename : Union[str, Path]) -> Path:
    start = Timestamp(dconf.start)
    end = Timestamp(dconf.end)

    with Dataset(str(filename), 'w') as ds :
        ds.createDimension('itime', 6)
        for regname, region in dconf.regions.items():

            # Create region object
            reg = TM5Grids.from_corners(
                east = region.lons[1],
                west = region.lons[0],
                south = region.lats[0],
                north = region.lats[1],
                dlon = region.lons[2],
                dlat = region.lats[2]
            )

            # Create netCDF group
            gid = ds.createGroup(regname)
            gid.createDimension('latitude', reg.nlat)
            gid.createDimension('longitude', reg.nlon)
            gid.latmin = reg.south
            gid.latmax = reg.north
            gid.lonmin = reg.west
            gid.lonmax = reg.east

            # Add the grid cell area
            gid.createVariable('area', 'd', ('latitude', 'longitude'))
            gid['area'][:] = reg.area

            # Create tracer subgroup(s):
            for trname, tracer in {k: dconf[k] for k in dconf.tracers}.items():

                # Create the tracer group
                gtr = gid.createGroup(trname)

                dailycycle_writemode = 'w'

                for catname, cat in tracer.get('emission_categories', {}).items():
                    logger.info(f'{regname}; {trname}; {catname}')
                    grp = gtr.createGroup(catname)

                    # Ensure that the optim_freq field is set for the category
                    cat.optim_freq = cat.get('optim_freq', 'D')
                    if cat.get('climatology', False):
                        # Enforce monthly for climatology
                        cat.optim_freq = 'MS'

                    #-- consistency check on allowed frequencies
                    assert cat.optim_freq in ['D','MS'], \
                        f"unexpcted value for optim_freq -->{cat.optim_freq}<--"

                    # Load the emissions
                    emcat = load_preprocessed_emis(cat, tracer.species, reg, slice(start, end))

                    #-- MVO::currently reverted back for CO2 and it's adjoint
                    # Split baseline and anomalies in all cases, to ensure that the emissions are at the requested frequency
                    emcat, anoms = split_baseline_anomalies(emcat, start, end, cat.optim_freq)
                    # #-- MVO::separate the two parts
                    # #        - resampling to target frequency
                    # #        - generation of daily emission anomalies
                    # emcat = emcat.resample(time=cat.optim_freq).mean()
                    # emcat['time'] = date_range(start, end, freq=cat.optim_freq, inclusive='left')

                    # Calculate and write dailycycle if needed
                    if cat.get('dailycycle', False):
                        #-- MVO:: TO BE ACTIVATED IN CASE WE SPLIT AS DESCRIBED ABOVE
                        # anoms = split_daily_anomalies(emcat, start, end)
                        # Write the anomalies for the category
                        for day, anom in anoms.items():
                            fname = day.strftime(tracer.dailycycle.filename_format)
                            Path(fname).parent.mkdir(parents=True, exist_ok=True)
                            anom.to_netcdf(fname, group=f'{regname}/{catname}', mode = dailycycle_writemode)
                        # Next category will should not overwrite the data:
                        dailycycle_writemode = 'a'

                    # Store the emission themselves:
                    grp.createDimension('nt', emcat.shape[0])
                    grp.createVariable('emission', 'd', ('nt', 'latitude', 'longitude'))
                    grp.createVariable('time_start', 'i2', ('nt', 'itime'))
                    grp.createVariable('time_end', 'i2', ('nt', 'itime'))
                    grp.time_resolution = cat.optim_freq
                    grp.treskey = {'MS': 'monthly', 'D': 'daily'}[cat.optim_freq]
                    grp.optimize = int(cat.get('optim', 0))

                    grp['emission'][:] = emcat.values
                    grp['emission'].units = f"kg{tracer.species}/cl/s"
                    grp['time_start'][:] = [t.timetuple()[:6] for t in DatetimeIndex(emcat.time)]
                    grp['time_end'][:] = [t.timetuple()[:6] for t in DatetimeIndex(emcat.time) + to_offset(cat.optim_freq)]

    return filename
