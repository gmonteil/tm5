#!/usr/bin/env python


"""
TM5 adjoint test:

- dy1 = H(x + dx1) - H(x)
- dx2 = H*(dy2)

Check that:
<dy1, dy2> = <dx1, dx2>
"""

from argparse import ArgumentParser
import sys
import os
from tm5.system import runcmd
from tm5.observations import read_obspack
from omegaconf import DictConfig
import tm5
from pandas import DataFrame, Timestamp, date_range, concat, Timedelta
import xarray as xr
from pathlib import Path
from numpy.random import randn
from numpy.typing import NDArray
from numpy import array, ones, zeros
from netCDF4 import Dataset
from typing import List
from loguru import logger


def preproc(ds: xr.Dataset) -> xr.Dataset:
    fname = Path(ds.encoding['source']).name
    time = Timestamp(fname.split('.')[-2])
    return ds.expand_dims({'time': [time]})
    

def emis_to_state(conf: DictConfig) -> DataFrame:
    start = Timestamp(conf.run.start)
    end = Timestamp(conf.run.end)
    statevec = []
    for region in conf.run.regions:
        for trname, tracer in conf.emissions.items():
            emfiles = [f'{tracer.prefix}.{trname}.{region}.{day:%Y%m%d}.nc' for day in date_range(start, end, freq='D', inclusive='left')]
            df = xr.open_mfdataset(emfiles, preprocess=preproc).to_array('category', name='emis').to_dataframe().reset_index()
            df.loc[:, 'region'] = region
            df.loc[:, 'tracer'] = trname
            statevec.append(df)
    return concat(statevec)
    
    
def state_to_emis(emis: DataFrame, conf: DictConfig) -> None:
    start = Timestamp(conf.run.start)
    end = Timestamp(conf.run.end)
    for region in conf.run.regions:
        for trname, tracer in conf.emissions.items():
            for day in date_range(start, end, freq='D'):
                
                # Select the portion of the dataframe corresponding to this timestep, region and tracer
                ds = emis[(emis.region == region) & (emis.tracer == trname) & (emis.time == day)]
                
                # Convert to a xarray dataset, extract the relevant dataArray, and convert it back to a dataset, with one variable for each category
                ds = ds.set_index(['lat', 'lon', 'category']).to_xarray().emis.to_dataset('category')
                
                # Write the dataset
                ds.to_netcdf(f'{tracer.prefix}.{trname}.{region}.{day:%Y%m%d}.nc')


def obsfile_to_vector(filename: Path | str, regions: List[str], tracers: List[str]) -> NDArray:
    y = []
    with Dataset(filename, 'r') as ds:
        for region in regions:
            for tracer in tracers:
                if tracer in ds[region].groups:
                    y.extend(ds[region][tracer]['mixing_ratio'][:])
    return array(y)


def gen_departures(point_input_file: Path, point_output_file: Path, point_departures_file: Path) -> NDArray:
    yadj = []
    writemode = 'w'
    for region in tm.dconf.run.regions:
        for tracer in tm.dconf.run.tracers:
            try:
                outfile = xr.open_dataset(point_output_file, group=f'{region}/{tracer}')
                infile = xr.open_dataset(point_input_file, group=f'{tracer}').sel(id=outfile.id)
                outfile['lat'] = ('samples', infile.lat.values)
                outfile['lon'] = ('samples', infile.lon.values)
                outfile['alt'] = ('samples', infile.alt.values)
                #outfile['forcing'] = ('samples', randn(outfile.sizes['samples']))
                outfile['forcing'] = ('samples', ones(outfile.sizes['samples']))
                outfile['time_window_length'] = ('samples', infile.time_window_length.values)
                outfile['date_components'] = (('samples', 'idate'), infile.date_components.values)
                outfile.to_netcdf(point_departures_file, group=f'{region}/{tracer}', mode=writemode)
                writemode = 'a'
                
                yadj.extend(outfile.forcing.values)
            except OSError:
                logger.info(f'No observations found for {region=} and {tracer=}')
                pass
    return array(yadj)
            
            
def read_adjvec(adjemis_dir : str | Path,
                start: str | Timestamp, end: str | Timestamp, conf: DictConfig) -> DataFrame:
    if isinstance(adjemis_dir, str):
        adjemis_dir = Path(adjemis_dir)
    adjx = []
    for itrac, tracer in enumerate(conf):
        for categ in conf[tracer].categories:
            for region in conf[tracer].categories[categ].regions:
                for day in date_range(start, end, freq='D', inclusive='left'):
                    outname = adjemis_dir / f'adjemis.{region}.{day:%Y%m%d}.nc'
                    ds = xr.open_dataset(outname).adj_emis.isel(tracer=itrac).to_dataframe().reset_index()
                    ds.loc[:, 'region'] = region
                    ds.loc[:, 'time'] = day
                    ds.loc[:, 'tracer'] = tracer
                    ds.loc[:, 'category'] = categ
                    adjx.append(ds)
    return concat(adjx)
    

parser = ArgumentParser()
parser.add_argument('-b', '--build', action='store_true', default=False, help='Use this option to compile the code')
parser.add_argument('-m', '--host', default=os.environ['TM5_HOST'])
parser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
parser.add_argument('--rcfile-only', action='store_true', default=False, help="""whether to only create the TM5 rcfile(s) (neither compiling nor running TM5).""")
parser.add_argument('config_file')
args = parser.parse_args(sys.argv[1:])

#=====================================================
# 1. Build the model
#=====================================================
tm = tm5.TM5(args.config_file, host=args.host, override_trange=args.trange)
if args.build and not args.rcfile_only:
    tm.build()
    
#=====================================================
# 2. Setup observation files:
#=====================================================

# observations_table = DataFrame.from_dict(dict(
#     time=[Timestamp(tm.dconf.run.end) - Timedelta(hours=12)],
#     lat=[51.6],
#     lon=[6.5],
#     alt=[150.],
#     mixing_ratio=[0.],
#     mixing_ratio_err=[1.],
#     id=[1],
#     time_window_length=[1],
#     sampling_strategy=[2],
#     station_id=[0],
#     tracer=['CH4']
# ))


observations_table = read_obspack(
    tm.dconf.output.point.filename,
    start = tm.dconf.run.start,
    end=tm.dconf.run.end
)
tm.setup_observations(observations_table)


#=====================================================
# 3. Setup the other input files:
#=====================================================

# Setup steps common to all runs
tm.setup_meteo()
tm.setup_regions()
tm.setup_system()
tm.setup_tracers()
tm.setup_output()
tm.setup_iniconc()

# 1st forward run
tm.setup_run('forward')

# Do the initial setup of emissions
tm.setup_emissions2()

# Create a pseudo control-vector from it:
x1 = emis_to_state(tm.dconf)
    
# Write the rc-file:
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'forward.rc')

# Run TM5
run_cmd = tm.dconf.run.run_cmd.split() + [str(rcf)]

if args.rcfile_only:
    msg = f"TM5 rcfile ***{str(rcf)}*** has been created"
    logger.info(msg)
    msg = f"...but don't actually run TM5 \n-->{' '.join(run_cmd)}<--\n"
    logger.info(msg)
else:
    runcmd(run_cmd)

    y1 = obsfile_to_vector(
        filename=Path(tm.dconf.run.paths.output) / 'point/point_output.nc4', 
        regions=tm.dconf.run.regions, 
        tracers=tm.dconf.run.tracers
    )

    # Create an emission perturbation, by applying a random scaling to the original emissions
    #dx = randn(len(x1))
    #dx /= dx

    dx = zeros(len(x1))
    dx[x1.region == 'gns100x100'] = 1.

    x2 = x1.copy()
    x2.loc[:, 'emis'] = x2.emis + dx
    state_to_emis(x2, tm.dconf)

    # Run TM5 again
    runcmd(tm.dconf.run.run_cmd.split() + [str(rcf)])
    y2 = obsfile_to_vector(
        filename=Path(tm.dconf.run.paths.output) / 'point/point_output.nc4', 
        regions=tm.dconf.run.regions, 
        tracers=tm.dconf.run.tracers
    )

    dy = y2 - y1

    yadj = gen_departures(
        point_input_file=Path(tm.dconf.output.point.input_dir) / 'point_input.nc4',
        point_output_file=Path(tm.dconf.run.paths.output) / 'point/point_output.nc4',
        point_departures_file=Path(tm.dconf.run.paths.output) / 'point/point_departures.nc4'
    )

#
#-- setup adjoint run
#
tm.setup_run('adjoint')
rcf = tm.settings.write(Path(tm.dconf.run.paths.output) / 'adjoint.rc')
adjrun_cmd = tm.dconf.run.run_cmd.split() + [str(rcf)]

if args.rcfile_only:
    msg = f"TM5 rcfile ***{str(rcf)}*** has been created"
    logger.info(msg)
    msg = f"...but don't actually run adjoint TM5 \n-->{' '.join(run_cmd)}<--\n"
    logger.info(msg)
else:
    runcmd(adjrun_cmd)

    if False: #--
        # adjemis_dir = <read from key in yaml file> (if present)
        pass
    else:
        #-- TM5 default output directory for adjoint emissions
        #   (currently see emission_adj.F90, line 154ff)
        adjemis_dir = Path(tm.dconf.run.paths.output) / 'adjemis'
    xadj = read_adjvec(adjemis_dir,
                       tm.dconf.run.start, tm.dconf.run.end, tm.dconf.emissions)

    # Adjoint test:
    x1 = x1.sort_values(['tracer', 'region', 'category', 'time', 'lat', 'lon'])
    x2 = x2.sort_values(['tracer', 'region', 'category', 'time', 'lat', 'lon'])
    xadj = xadj.sort_values(['tracer', 'region', 'category', 'time', 'lat', 'lon'])

    dxdx = (x2.emis - x1.emis).values @ xadj.adj_emis.values
    dydy = dy @ yadj

    print(f'dx @ xadj = {dxdx}; dy @ yadj = {dydy}; err (%) = {1 - dydy/dxdx}')
