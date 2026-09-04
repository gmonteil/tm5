#!/usr/bin/env python
import sys
from loguru import logger
from argparse import ArgumentParser
from pathlib import Path
from pandas import DataFrame, read_hdf, Timestamp, Timedelta, to_timedelta, concat, date_range
import xarray as xr
from h5py import File
from numpy import array, zeros, where
import numpy as np
from omegaconf import OmegaConf, DictConfig
from tqdm import tqdm
import pickle
import lzma
from types import SimpleNamespace

#-- MVO specific settings:
emis_path = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/input/ch4/emission-input-tm5/fitic-default_20201001-20211231')
pathtable_gns1x1 = SimpleNamespace(
    footprint_dir = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/runs_footprint_with-chemistry/fitic-footprint-ifort-slurm_simustart-20201001'),
    footprint_path = 'footprints_gns100x100_20201001--%Y%m%d',
    #--> now using full reference simulation from Oct 1, 2020 to Dec 31, 2021
    tm5_fwd_path = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/runs_fitic-forward_with-chemistry/fitic-simu-default_gns-obs_20201001--20220101'),
    # tm5_fwd_path = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/runs_fitic-forward_with-chemistry/fitic-simu-default_gns-obs_20201001--20210201')
)
pathtable_glb6x4 = SimpleNamespace(
    footprint_dir = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/runs_footprint_with-chemistry/glb6x4-footprint-ifort-slurm_simustart-20201001'),
    footprint_path = 'footprints_glb600x400_20201001--%Y%m%d',
    tm5_fwd_path = Path('/lunarc/nobackup/projects/ghg_inv/michael/TM5/expdir/runs_fitic-forward_with-chemistry/glb6x4-simu-default_flask-obs_20201001-20220101'),
)
#
start = Timestamp(2020, 10, 1)
obs_firstday = Timestamp(2021,1,1)

parser = ArgumentParser()
parser.add_argument('--mode',
                    choices=['gns1x1','glb6x4'],
                    default='gns1x1',
                    help="""whether to process the FIT-IC zoom or the global footprints (default: %(default)s).""")
parser.add_argument('--obs_lastday',
                    type=Timestamp,
                    default=Timestamp(2021,1,4),
                    help="""last observational day (default:%(default)s).""")
parser.add_argument('--stations',
                    nargs='+',
                    help="""restrict to selected stations.""")
parser.add_argument('--pickle_filepath',
                    type=Path,
                    help="""use serialized footprint information from previous run (and don't reparse footprints.""")
parser.add_argument('--outdir',
                    type=Path,
                    help="""destination directory for all generated outputs.""")

#-- parse arguments
args = parser.parse_args(sys.argv[1:])

#-- handle arguments
#- last day
obs_lastday = args.obs_lastday
dir_end = obs_lastday + Timedelta(days=1)
#-- actual paths depend on mode
if args.mode=='gns1x1':
    pathtable = pathtable_gns1x1
    regions = ['glb600x400', 'eur300x200', 'gns100x100',]
elif args.mode=='glb6x4':
    pathtable = pathtable_glb6x4
    regions = ['glb600x400',]
footprint_dir  = pathtable.footprint_dir
footprint_path = pathtable.footprint_path
tm5_fwd_path   = pathtable.tm5_fwd_path

# print(f"footprint_dir ***{footprint_dir}***")
# print(f"footprint_path ***{footprint_path}***")
# sys.exit(0)
#--
mode_tag = args.mode
time_tag = f"{start.strftime('%Y%m%d')}--{obs_lastday.strftime('%Y%m%d')}"
obstime_tag = f"{obs_firstday.strftime('%Y%m%d')}--{obs_lastday.strftime('%Y%m%d')}"
station_tag = None
if args.stations!=None:
    station_tag = '--'.join(args.stations)

# ------------------------------------
# 1./2. Load observations AND forward simulation
#
msg = f"loading observations and TM5 forward simulation from reference run..."
logger.info(msg)
outpfile = tm5_fwd_path / 'point/point_output.nc4'
inpfile  = tm5_fwd_path / 'point_input.nc4'
with File(inpfile) as inpf:
    with File(outpfile) as oupf:
        #-- collect fields from point_output.nc4
        tm5mix, station_id, ido, reg = [], [], [], []
        for region in regions:
            if 'CH4' in oupf[region]:
                tm5mix.extend(oupf[region]['CH4']['mixing_ratio'])
                station_id.extend(oupf[region]['CH4']['station_id'])
                ido.extend(oupf[region]['CH4']['id'])
                reg.extend([region] * oupf[region]['CH4']['id'].shape[0])
        tm5mix = array(tm5mix)
        station_id = array(station_id)
        ido = array(ido)
        #-- fields from point_input.nc4
        idi = inpf['CH4']['id'][:]
        #-- time stamp
        vartime = inpf['CH4']['time']
        time_units = vartime.attrs['units'].decode('UTF-8')
        if time_units.startswith('hours since '):
            reftime = Timestamp(time_units.replace('hours since ',''))
            time = reftime +  to_timedelta(vartime[:], unit='hour')
        elif time_units.startswith('seconds since '):
            reftime = Timestamp(time_units.replace('seconds since ',''))
            time = reftime +  to_timedelta(vartime[:], unit='second')
        else:
            msg = f"...unexpected time-units -->{time_units}<-- @{str(inpfile)}"
            raise RuntimeError(msg)
        lon = inpf['CH4']['lon'][:]
        lat = inpf['CH4']['lat'][:]
        alt = inpf['CH4']['alt'][:]
        obsmix = inpf['CH4']['mixing_ratio'][:]
        obsmixerr = inpf['CH4']['mixing_ratio_err'][:]
        time_window_length = inpf['CH4']['time_window_length'][:]
        sampling_strategy = inpf['CH4']['sampling_strategy'][:]
        station_id = inpf['CH4']['station_id'][:]
        tracer = inpf['CH4']['tracer'][:]
        obsid = inpf['CH4']['obsid'][:]
        # Construct a DataFrame with mix, obsid and time
        obstable = DataFrame(
            { 'time': time,
              'lon':lon,
              'lat':lat,
              'alt':alt,
              'mixing_ratio':obsmix,
              'mixing_ratio_err':obsmixerr,
              'time_window_length':time_window_length,
              'sampling_strategy':sampling_strategy,
              'station_id':station_id,
              'tracer':tracer,
              'obsid': obsid,
              'index': idi,
              'region': reg }
        ).set_index('index')
        #-- add output
        obstable.loc[ido, 'tm5_fwd'] = tm5mix
        #-- convert byte-string to normal string
        obstable.loc[:, 'obsid'] = obstable.obsid.str.decode('utf8')
        obstable.loc[:, 'tracer'] = obstable.tracer.str.decode('utf8')
#
# outname = f"reference-obstable_with-tm5-fwd_{process_tag}.csv"
# obstable.to_csv(outname)
msg = f"...obstable done."
logger.info(msg)
#
#-- restrict to selected observational period and potentially to selected stations
#
cnd_time = obstable['time']<dir_end
obstable = obstable.loc[cnd_time,:]
outname_tokens = [f'obstable', mode_tag, 'with-tm5-fwd', obstime_tag]
if args.stations!=None:
    cnd_station = obstable['obsid'].isin(args.stations)
    obstable = obstable.loc[cnd_station,:]
    outname_tokens += [station_tag,]
outname = '_'.join(outname_tokens) + f".csv"
if args.outdir!=None:
    outname = args.outdir / outname
    outname.parent.mkdir(parents=True, exist_ok=True)
obstable.to_csv(outname)

# ------------------------------------
# 3. Load initial condition
# There should be a much more straightforward way to do it: just run a forward run without emissions
# But for consistency, I did it based on the footprint runs (which start with a forward without emissions)
tm5_df = []
for date in obstable.time.dt.date.drop_duplicates():
    # logger.info(f"iniconc@{date}...")
    tmpath = footprint_dir / (date + Timedelta(days=1)).strftime(footprint_path)
    outpfile = tmpath / 'point/point_output.nc4'
    inpfile = tmpath / 'point_input.nc4'
    # print(type(date), type(end))
    if date>dir_end.date():
        continue
    elif not inpfile.exists():
        msg = f"***{str(inpfile)}*** NOT FOUND"
        logger.debug(msg)
        continue
    # msg = f"extracting iniconc from files ***{str(outpfile)}*** ***{str(inpfile)}***"
    # logger.debug(msg)
    with File(inpfile) as inpf:
        with File(outpfile) as oupf:
            mix, station_id, ido, reg = [], [], [], []
            for region in regions:
                if 'CH4' in oupf[region]:
                    mix.extend(oupf[region]['CH4']['mixing_ratio'])
                    station_id.extend(oupf[region]['CH4']['station_id'])
                    ido.extend(oupf[region]['CH4']['id'])
                    reg.extend([region] * oupf[region]['CH4']['id'].shape[0])
            mix = array(mix)
            station_id = array(station_id)
            ido = array(ido)
            idi = inpf['CH4']['id'][:]
            obsid = inpf['CH4']['obsid'][:]
            vartime = inpf['CH4']['time']
            time_units = vartime.attrs['units'].decode('UTF-8')
            # msg = f"time_units -->{time_units}<--"
            # logger.debug(msg)
            if time_units.startswith('days since '):
                reftime = Timestamp(time_units.replace('days since ',''))
                time = reftime +  to_timedelta(vartime[:], unit='day')
            elif time_units.startswith('hours since '):
                reftime = Timestamp(time_units.replace('hours since ',''))
                time = reftime +  to_timedelta(vartime[:], unit='hour')
            elif time_units.startswith('minutes since '):
                reftime = Timestamp(time_units.replace('minutes since ',''))
                time = reftime +  to_timedelta(vartime[:], unit='minute')
            elif time_units.startswith('seconds since '):
                reftime = Timestamp(time_units.replace('seconds since ',''))
                time = reftime +  to_timedelta(vartime[:], unit='second')
            else:
                msg = f"...unexpected time-units -->{time_units}<-- @{str(inpfile)}"
                raise RuntimeError(msg)
            #
            # time = Timestamp(date) + to_timedelta(inpf['CH4']['time'][:] + 1, unit='hour')  # Here there is a weird extra hour ...

            # Construct a DataFrame with mix, obsid and time
            # MVO: have region already above, can drop here
            # df = DataFrame({'time': time, 'obsid': obsid, 'index': idi, 'region': reg}).set_index('index')
            df = DataFrame({'time': time, 'obsid': obsid, 'index': idi,}).set_index('index')
            df.loc[ido, 'iniconc'] = mix
            df.loc[:, 'obsid'] = df.obsid.str.decode('utf8')

            # Store
            tm5_df.append(df)

# Merge:
tm5_df = concat(tm5_df)
outname_tokens = [f'iniconc', mode_tag, obstime_tag]
if args.stations!=None:
    outname_tokens += [station_tag,]
outname = '_'.join(outname_tokens) + f".csv"
if args.outdir!=None:
    outname = args.outdir / outname
    outname.parent.mkdir(parents=True, exist_ok=True)
tm5_df.to_csv(outname)

#
obstable = obstable.merge(tm5_df, on=['time', 'obsid'])

outname_tokens = ["obstable", mode_tag, "with-tm5-fwd_with-iniconc", obstime_tag]
if args.stations!=None:
    outname_tokens += [station_tag,]
outname = '_'.join(outname_tokens) + f".csv"
if args.outdir!=None:
    outname = args.outdir / outname
    outname.parent.mkdir(parents=True, exist_ok=True)
obstable.to_csv(outname, index=True)
logger.info(f"generated ***{str(outname)}***")
msg = f"...initial concentrations loaded len(obstable)={len(obstable)}"
logger.info(msg)
#

# ------------------------------------
# 4. Load all footprints
if args.pickle_filepath==None:
    msg = f"need to load/parse footprint from simulations in ***{str(footprint_dir)}***"
    logger.info(msg)
    ndays = int((dir_end - start).total_seconds() / 86400)
    footprints = {}
    for obs in tqdm(obstable.itertuples(), desc='loading footprints', total=len(obstable)):
        # Reconstruct the footprint
        if args.mode=='gns1x1':
            footprints[obs.Index] = {
                'glb600x400': zeros((ndays, 45, 60)),
                'eur300x200': zeros((ndays, 26, 30)),
                'gns100x100': zeros((ndays, 16, 18))
            }
        elif args.mode=='glb6x4':
            footprints[obs.Index] = {
                'glb600x400': zeros((ndays, 45, 60)),
            }
        # footprint path:
        dir_date = obs.time.date() + Timedelta(days=1)
        fppath = footprint_dir / dir_date.strftime(footprint_path) / 'adjemis'

        for idate, date in enumerate(tqdm(date_range(start, obs.time, freq='D'))):
            for region in regions:
                fpfile = fppath / f'adjemis.{region}.{date:%Y%m%d}.nc'
                fp = xr.open_dataset(fpfile)
                tracers = fp.tracer.str.strip().str.decode('utf8')
                trindex = where(tracers == obs.obsid)[0]
                if len(trindex) > 0:
                    ilats = fp.ilat.values[fp.itrac.values == trindex]
                    ilons = fp.ilon.values[fp.itrac.values == trindex]
                    values = fp['values'].values[fp.itrac.values == trindex]
                    footprints[obs.Index][region][idate, ilats, ilons] = values
        msg = f"...{obs.obsid}@{obs.time} done"
        logger.debug(msg)
    msg = f"...re-parse footprints done"
    logger.info(msg)
    outname_tokens = ['footprints', mode_tag, time_tag,]
    if args.stations!=None:
        outname_tokens += [station_tag,]
    # outname = '_'.join(outname_tokens) + f".xz"
    outname = '_'.join(outname_tokens) + '.pickle'
    if args.outdir!=None:
        outname = args.outdir / outname
        outname.parent.mkdir(parents=True, exist_ok=True)
    # with lzma.open(outname, 'wb') as fid:
    #     pickle.dump(footprints, fid)
    #-- MVO-NOTE::compression ratio was not that huge,
    #             but with lzma reading/writing take much much longer
    with open(outname, 'wb') as fid:
        pickle.dump(footprints, fid)
    msg = f"...footprints written to {str(outname)}"
    logger.info(msg)
elif not args.pickle_filepath.exists():
    msg = f"serialised footprints ***{str(args.pickle_filepath)}*** not found."
    raise RuntimeError(msg)
elif args.pickle_filepath.suffix=='.pickle':
    with open(args.pickle_filepath, 'rb') as fid:
        msg = f"start reading ***{str(args.pickle_filepath)}***"
        logger.info(msg)
        footprints = pickle.load(fid)
        msg = f"...pickle file loaded"
        logger.info(msg)
elif args.pickle_filepath.suffix=='.xz':
    with lzma.open(args.pickle_filepath, 'rb') as fid:
        msg = f"start reading ***{str(args.pickle_filepath)}***"
        logger.info(msg)
        footprints = pickle.load(fid)
        msg = f"...pickle file loaded"
        logger.info(msg)
else:
    msg = f"serialised footprints ***{str(args.pickle_filepath)}*** with unexpected suffix!"
    raise RuntimeError(msg)

# ------------------------------------
# 5. Load the emissions
msg = f"start loading emissions..."
logger.info(msg)
if args.mode=='gns1x1':
    emis = {
        'gns100x100': xr.open_mfdataset(emis_path.glob('ch4emis.CH4.gns100x100.????????.nc'), concat_dim='time', combine='nested'),
        'eur300x200': xr.open_mfdataset(emis_path.glob('ch4emis.CH4.eur300x200.????????.nc'), concat_dim='time', combine='nested'),
        'glb600x400': xr.open_mfdataset(emis_path.glob('ch4emis.CH4.glb600x400.????????.nc'), concat_dim='time', combine='nested')
    }
elif args.mode=='glb6x4':
    emis = {
        'glb600x400': xr.open_mfdataset(emis_path.glob('ch4emis.CH4.glb600x400.????????.nc'), concat_dim='time', combine='nested')
    }
#--
for reg in emis:
    # Sort by time and select only the slice corresponding to the simulation
    # xarrray.sel is inclusive on both ends, so use isel after that to drop the last value
    emis[reg] = emis[reg].sortby('time').sel(time=slice(start, dir_end)).isel(time=slice(None, -1)).to_array().sum('variable').compute()
msg = f"...emissions loaded ({emis.keys()})"
logger.info(msg)

# ------------------------------------
# 6. Perform the forward run
msg = f"propagate emissions with linear forward..."
logger.info(msg)
for obs in obstable.itertuples():
    for reg in emis:
        msg = f"...applying linear forward for {reg}"
        logger.info(msg)
        obstable.loc[obs.Index, reg] = (footprints[obs.Index][reg] * emis[reg].values).sum()
msg = f"...linear forward propagation done"
logger.info(msg)

obstable.loc[:, 'forward_footprints'] = obstable.iniconc.values
for reg in emis:
    obstable.loc[:, 'forward_footprints'] += obstable.loc[:, reg]
msg = f"...linear forward run done!"
logger.info(msg)

outname_tokens = ['fwd_approach_comparisons', mode_tag, time_tag,]
if args.stations!=None:
    outname_tokens += [station_tag,]
outname = '_'.join(outname_tokens) + '.h5'
if args.outdir!=None:
    outname = args.outdir / outname
    outname.parent.mkdir(parents=True, exist_ok=True)
obstable.to_hdf(outname, key='tm5')
msg = f"generated ***{str(outname)}***"
logger.info(msg)
