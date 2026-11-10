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
import subprocess
from tm5.system import runcmd
from tm5.observations import read_obspack
from omegaconf import OmegaConf, DictConfig
import tm5
from pandas import DataFrame, Timestamp, date_range, concat, Timedelta
import xarray as xr
from pathlib import Path
from numpy.random import randn
from numpy.typing import NDArray
from numpy import array, ones, zeros, arange
from netCDF4 import Dataset
from typing import List
from loguru import logger


def run_logged(cmd: List[str], logfile: Path) -> None:
    """
    Run a command, showing its output (stdout + stderr) in the terminal and copying it to `logfile`.
    Stop the script if the command fails (tm5.system.runcmd does not check the exit status).
    """
    logger.info(f"running ==>{' '.join(cmd)}<==, log in {logfile}")
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as p, open(logfile, 'wb') as log:
        for line in p.stdout:
            sys.stdout.buffer.write(line)
            sys.stdout.buffer.flush()
            log.write(line)
    if p.returncode != 0:
        logger.error(f'command failed with exit status {p.returncode}, see {logfile}')
        sys.exit(p.returncode)


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
                if region in ds.groups:
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
            
            
def read_adjvec(adjemis_dir: str | Path, conf: DictConfig) -> DataFrame:
    """
    Read the sparse adjoint emission files written by TM5 (emission_adj.F90/write_adj_emis)
    and return a dense DataFrame aligned with emis_to_state().

    File format (one file per region and day, adjemis.<region>.<YYYYMMDD>.nc):
      values(point), ilat(point), ilon(point), itrac(point) : non-zero values, 0-based indices
      lat(lat), lon(lon), tracer(tracer, tracer_nchar)      : coordinates
    TM5 does not write the file at all if all values are zero: those days are filled with zeros.
    The adjoint is per tracer (shared by all categories), so it is repeated for each category.
    """
    adjemis_dir = Path(adjemis_dir)
    adjx = []
    for tracer in conf.emissions:
        for categ in conf.emissions[tracer].categories:
            for region in conf.emissions[tracer].categories[categ].regions:

                # Region grid (cell centres), from the region definition: [start, end, step]
                lon0, lon1, dlon = conf.regions[region].lons
                lat0, lat1, dlat = conf.regions[region].lats
                lons = arange(lon0 + dlon / 2, lon1, dlon)
                lats = arange(lat0 + dlat / 2, lat1, dlat)

                for day in date_range(conf.run.start, conf.run.end, freq='D', inclusive='left'):
                    adj = zeros((lats.size, lons.size))

                    fname = adjemis_dir / f'adjemis.{region}.{day:%Y%m%d}.nc'
                    if fname.exists():
                        with xr.open_dataset(fname) as ds:
                            names = [_.decode().strip() for _ in ds.tracer.values]
                            tracer_adj = ds.isel(point=ds.itrac == names.index(tracer))
                            adj[tracer_adj.ilat.values, tracer_adj.ilon.values] = tracer_adj['values'].values

                    df = xr.DataArray(
                        adj,
                        coords={'lat': lats, 'lon': lons},
                        dims=('lat', 'lon'),
                        name='adj_emis'
                    ).to_dataframe().reset_index()
                    df.loc[:, 'region'] = region
                    df.loc[:, 'time'] = day
                    df.loc[:, 'tracer'] = tracer
                    df.loc[:, 'category'] = categ
                    adjx.append(df)
    return concat(adjx)
    

parser = ArgumentParser()
parser.add_argument('-b', '--build', action='store_true', default=False, help='Use this option to compile the code')
parser.add_argument('-m', '--host', default=os.environ['TM5_HOST'])
parser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
parser.add_argument('config_file')
args = parser.parse_args(sys.argv[1:])

#=====================================================
# 1. Build the model
#=====================================================
#
#-- load configuration file
#
dconf = OmegaConf.load(args.config_file)
#
#-- potential partial override (or extend) configuration
#
if args.trange!=None:
    tstart, tend = args.trange
    dconf.run['start'] = tstart
    dconf.run['end']   = tend

tm = tm5.TM5(dconf, host=args.host)
if args.build:
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

msg = f"first TM5 forward run, exectuing command ==>{run_cmd}<=="
logger.info(msg)
run_logged(run_cmd, Path(tm.dconf.run.paths.output) / 'forward_1.log')

y1 = obsfile_to_vector(
    filename=Path(tm.dconf.run.paths.output) / 'point/point_output.nc4', 
    regions=tm.dconf.run.regions, 
    tracers=tm.dconf.run.tracers
)

# Create an emission perturbation, by applying a random scaling to the original emissions
#dx = randn(len(x1))
#dx /= dx

dx = zeros(len(x1))
# dx[x1.region == 'gns100x110'] = 1.
# dx[x1.region == 'eur300x200'] = 1.
# dx[x1.region == 'glb600x400'] = 1.
dx += 1.

x2 = x1.copy()
x2.loc[:, 'emis'] = x2.emis + dx
state_to_emis(x2, tm.dconf)

# Run TM5 again
msg = f"second TM5 forward run, exectuing command ==>{run_cmd}<=="
logger.info(msg)
run_logged(run_cmd, Path(tm.dconf.run.paths.output) / 'forward_2.log')
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

msg = f"TM5 adjoint run, exectuing command ==>{adjrun_cmd}<=="
logger.info(msg)
run_logged(adjrun_cmd, Path(tm.dconf.run.paths.output) / 'adjoint.log')

if False: #--
    # adjemis_dir = <read from key in yaml file> (if present)
    pass
else:
    #-- TM5 default output directory for adjoint emissions
    #   (currently see emission_adj.F90, line 154ff)
    adjemis_dir = Path(tm.dconf.run.paths.output) / 'adjemis'
xadj = read_adjvec(adjemis_dir, tm.dconf)

# Adjoint test:
x1 = x1.sort_values(['tracer', 'region', 'category', 'time', 'lat', 'lon'])
x2 = x2.sort_values(['tracer', 'region', 'category', 'time', 'lat', 'lon'])
xadj = xadj.sort_values(['tracer', 'region', 'category', 'time', 'lat', 'lon'])

dxdx = (x2.emis - x1.emis).values @ xadj.adj_emis.values
dydy = dy @ yadj

print(f'dx @ xadj = {dxdx}; dy @ yadj = {dydy}; err (%) = {100 - 100 * dydy / dxdx}')
