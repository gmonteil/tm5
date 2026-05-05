#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
import datetime as dtm
from omegaconf import OmegaConf, DictConfig
from pathlib import Path
import datetime as dtm
from loguru import logger
from pandas import date_range, DataFrame
from pandas import Timestamp, Timedelta, concat
import xarray as xr
import numpy as np
from numpy import zeros, tile
from netCDF4 import Dataset
from types import SimpleNamespace
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots,colorbar
from cartopy import crs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from tm5.fitic import read_obs_table
from tm5.gridtools import TM5Grids
from tm5.observations import read_obspack_file
from tm5.post.footprint_io import tm5_fitic_adjoint_corrected_halos
from tm5.post.footprint_io import tm5_fitic_footprint4jacobian_v1
from tm5.post.footprint_io import tm5rundir_obsids_extra, tm5rundir_obstable
from tm5.post.footprint_io import tm5rundir_iniconc_1obs
from tm5.post.footprint_io import tm5rundir_emissions1D
from tm5.post.footprint_io import tm5rundir_jacobian2D
from tm5.post.footprint_io import load_adjoint_fwd
from tm5.post.plot_util import cnorm_set
from tm5.post.utilities import lonstr,latstr,set_outname


def tm5refdir_load_stationconc( refdir : str | Path, obsid : str ) -> xr.DataArray:
    station_file = Path(refdir) / 'stations' / 'stations.nc4'
    if not station_file.exists():
        msg = f"stations result file ***{str(station_file)}*** not found on system."
        raise FileNotFoundError(msg)
    #
    #--
    #
    ds = Dataset(station_file)
    assert ds.dimensions['tracers'].size==1, \
        f"multiple tracers in station file not yet supported " \
        f"({str(station_file)})"
    grp_fnd = None
    for g,gg in ds.groups.items():
        if gg.abbr.endswith(obsid):
            grp_fnd = gg
    #
    if grp_fnd is None:
        raise RuntimeError(f"-->{obsid}<-- not found in station file ***{str(station_file)}***")
    else:
        msg = f"obsid={obsid} found at group {grp_fnd.abbr}"
        logger.info(msg)
    #
    time_list = np.array([Timestamp(*_) for _ in ds['date_midpoints'][:]])
    ncvar = grp_fnd['/mixing_ratio']
    if not ncvar.dimensions!=('tracer','samples'):
        msg = f"mixing_ratio with unexpected dimensions -->{ncvar.dimensions}<--"
        raise RuntimeError(msg)
    conc = ncvar[0,:]
    df_conc = DataFrame.from_dict( {'time':time_list, 'conc': conc} )
    # conc_da = xr.DataArray(
    #     conc,
    #     dims = ('time',),
    #     coords = {
    #         'time': time_list
    #         }
    #     )
    # return conc_da
    return df_conc


def subcmd_test_jacobianfwd_1day(args):
    tm5rundir = args.tm5rundir
    obsid     = args.obsid
    if args.trange!=None:
        trange = date_range(args.trange[0], args.trange[1], freq='1d')
        print(trange)
    else:
        trange = args.trange

    df_fwd = load_adjoint_fwd(tm5rundir, trange=trange)

    cndobs = df_fwd.loc[:,'obs']==obsid
    df_obs = df_fwd.loc[cndobs,:]
    print(df_obs)
    #-- sum over regions
    df_out = df_obs.groupby(['obs', 'region']).sum(numeric_only=True).reset_index()
    print(df_out)
    csim = df_obs.loc[:,'mix'].sum()
    print(f"csim={csim}")


def subcmd_testbuild_jacobian_1day( args ):
    tm5rundir = args.tm5rundir
    obsid     = args.obsid
    
    # footp_info = tm5_fitic_adjoint_corrected_halos(tm5rundir)
    # print(footp.loc[:,'itime'].unique())
    if args.trange!=None:
        trange = date_range(args.trange[0], args.trange[1], freq='1d')
        print(trange)
    else:
        trange = args.trange
    #
    #-- load obs table used for this simulation
    #
    obs_table = tm5rundir_obsids_extra(tm5rundir)
    if not obsid in obs_table:
        msg = f"selected observation identifier -->{obsid}<-- not present in simulation " \
            f"(==>{list(obs_table.keys())}<==)"
        raise RuntimeError(msg)
    obs_info = obs_table[obsid]
    print(obs_info)
    #
    #-- load emissions
    #
    msg = f"start preparing emissions..."
    logger.info(msg)
    emis_info = tm5rundir_emissions1D(tm5rundir, trange=trange)
    emis1D = emis_info.emissions
    msg = f"...emissions done."
    logger.info(msg)
    #
    #-- load Jacobian
    #
    # jacobian_info = tm5_fitic_footprint4jacobian_v1(tm5rundir, trange=trange, obsid=obsid)
    # jac_table = jacobian_info.jac_table
    # days = jacobian_info.days #-- [emission days]
    # nday = len(days)
    # ngc = 0
    # for region,region_jac in jac_table.items():
    #     iday = len(region_jac.coords['iday'])
    #     nlon = len(region_jac.coords['lon'])
    #     nlat = len(region_jac.coords['lat'])
    #     if iday!=nday:
    #         msg = f"unexpected dimension size iday={iday} (expected {nday})"
    #         raise RuntimeError(msg)
    #     ngc += nlat*nlon #_ngc
    #     for iday,day in enumerate(days):
    #         _jac = region_jac.sel(iday=iday)
    #         msg = f"@{region},{day.strftime('%Y-%m-%d')}: jac min/mean/max = " \
    #             f"{_jac.min().values}/{_jac.mean().values}/{_jac.max().values}"
    #         print(msg)
    # msg = f"overall number of grid-cells ngc={ngc}"
    # print(msg)
    # #
    # #-- reformat Jacobian values to 1D array
    # #   ordering:
    # #   - emisday
    # #     - grid-cells glb600x400 (flattened)
    # #     - grid-cells eur300x200 (flattened)
    # #     - grid-cells gns100x100 (flattened)
    # #
    # jacobian2D = np.array([]) #empty()
    # for iday in range(nday):
    #     for reg,reg_jac in jac_table.items():
    #         dayreg_jac = reg_jac.values[iday,:,:]
    #         jacobian2D = np.concat( (jacobian2D, dayreg_jac.ravel()))
    #     # print(f"@iday={iday}, jacobian2D.shape={jacobian1D.shape}")
    #-- consistency check
    # assert emis1D.shape==jacobian1D.shape, \
    #     f"{emis1D.shape} vs {jacobian1D.shape}"
    # #
    # #--
    # #
    # jacobian2D = jacobian1D.reshape(1,len(jacobian1D))
    #
    jac_info = tm5rundir_jacobian2D(tm5rundir, trange=trange, obsid=obsid)
    jacobian2D = jac_info.jac2D
    nobs,nemis = jacobian2D.shape
    # print(f"jacobian2D.shape={jacobian2D.shape}")
    for iobs in range(nobs):
        msg = f"@iobs={iobs}, jac min/mean/max = " \
            f"{jacobian2D[iobs,:].min()}/{jacobian2D[iobs,:].mean()}/{jacobian2D[iobs,:].max()}"
        print(msg)
    assert len(emis1D)==nemis, \
        f"expected {len(emis1D)} emissions, but Jacobian shape is {jacobian2D.shape}"
    cobs = np.dot(jacobian2D,emis1D)
    print(f"cobs -->{cobs}<--")


def subcmd_testbuild_jacobian_period(args):
    """Test preparation of inputs for Fortran-based inversion environment
    based on sensitivities for one single observational site *and* a selected
    period of observational days.
    """
    #-- arguments
    topdir = Path(args.outpath_tm5)
    dayl = args.selday
    host = args.__dict__.get('host', 'cosmos')
    # dayf, dayl = args.period
    obsid = args.obsid
    complevel = args.__dict__.get('complevel',4)
    
    #-- turn dates into timestamp
    dayl = Timestamp(dayl)
    dayf = dayl - Timedelta(days=args.days) #Timestamp(dayf)

    #
    #-- naming pattern used by Guillaume: footprints_gns100x100_%Y%m%d,
    #   where '%Y%m%d' refers to last day of simulation at 0:00 !
    # e.g. for February we need to pick directories
    # footprints_gns100x100_20210202 to footprints_gns100x100_20210301
    ddayf = dayf + Timedelta(days=1)
    ddayl = dayl + Timedelta(days=1)
    dir_trange = date_range(ddayf, ddayl, freq='1d')
    day_range = date_range(dayf, dayl, freq='1d')
    nday = len(day_range)
    fdir = topdir / f"footprints_gns100x100_{dir_trange[0].strftime('%Y%m%d')}"
    #
    #-- load observation table
    #   -> we expect ethedeliberately select from the rundir of the first day of the selected
    #      period
    #
    obs_table = tm5rundir_obstable(fdir)
    obs_info = obs_table.loc[obsid,:]
    # print(obs_info, type(obs_info))
    #->
    obs_hr = obs_info.time.hour
    obs_tw = obs_info.time_window_length #-- time-window length [s]
    #-- start/end hour
    obs_hrs = (obs_info.time - Timedelta(seconds=obs_tw)).hour
    obs_hre = (obs_info.time + Timedelta(seconds=obs_tw)).hour
    assert obs_hre>obs_hrs, \
        f"obs-time crossing 0:00 not yet supported (hrs={hrs},hre={hre})"
    #-- "station list"
    station_list = [obsid,]
    
    #
    #-- load emissions
    #
    ldir = topdir / f"footprints_gns100x100_{dir_trange[-1].strftime('%Y%m%d')}"
    emis_info = tm5rundir_emissions1D(ldir, trange=day_range)
    emis1D = emis_info.emis1D
    iday1D = emis_info.iday1D
    region1D = emis_info.region1D
    nemis = len(emis1D)
    #
    #-- iniconc
    #
    iniconc_list = []
    for idir,dirday in enumerate(dir_trange):
        rundir = topdir / dirday.strftime(f"footprints_gns100x100_%Y%m%d")
        if not rundir.is_dir():
            msg = f"...expected directory -->{str(rundir)}<-- not found!"
            raise RuntimeError(msg)
        inic_info = tm5rundir_iniconc_1obs(rundir, obs_info)
        iniconc_list.append(inic_info.conc)
    # iniconc_a = np.array(iniconc_list)
    #
    #-- load observations
    #
    obspack_info = read_obspack_file(args.obsfile, start=dayf, end=dayl+Timedelta(seconds=86399))
    obs_df = obspack_info.data
    # print(obs_df.columns)
    obs_list = []
    for day in day_range:
        _ostart = day + Timedelta(hours=obs_hr) - Timedelta(seconds=obs_tw)
        _oend   = day + Timedelta(hours=obs_hr) + Timedelta(seconds=obs_tw)
        cnd_day = (obs_df['time']>=_ostart)&(obs_df['time']<=_oend)
        mix_day = obs_df.loc[cnd_day,'value'].mean()
        print(f"@{day}, _ostart/_oend = {_ostart}/{_oend}, mix={mix_day}")
        obs_list.append(mix_day*1.e9) #-- convert [mol/mol] to [ppb]

    #
    #-- load Jacobians (differentiated by zoom region) for each observational time
    #
    jacobian_list = []
    for idir,dirday in enumerate(dir_trange):
        #-- this was only to trace back consistency differences,
        #   at it turned out, Guillaume had been using differing
        #   emission files for
        #   footprints_gns100x100_20210301 and the (obs/simulation) days before
        # if iobs<nobs-2:
        #     continue
        msg = f"dir@{dirday.strftime('%Y-%m-%d')}, starting..."
        logger.info(msg)
        rundir = topdir / dirday.strftime(f"footprints_gns100x100_%Y%m%d")
        if not rundir.is_dir():
            msg = f"...expected directory -->{str(rundir)}<-- not found!"
            raise RuntimeError(msg)
        #--
        #
        jac_info = tm5rundir_jacobian2D(rundir, trange=day_range, obsid=obsid)
        jacobian_list.append(jac_info.jac2D)
        # xsim = np.dot(jac_info.jac2D, emis1D)
        # print(f"{dirday}, xsim = {xsim}")
    #
    #--
    #
    jac2D = np.vstack(jacobian_list)
    # print(f"jac2D.shape={jac2D.shape}")
    # csim = np.dot(jac2D, emis1D)
    # print(f"csim -->{csim}<--")

    #
    #-- verification
    #
    if args.refdir!=None:
        refdir = Path(args.refdir)
        #
        #-- emissions used in reference run
        #
        msg = f"reading emissions from reference directory ***{refdir}***..."
        logger.info(msg)
        tm5rundir_emissions1D(refdir, trange=day_range)#, host='cosmos_apptainer')
        refemis_info = tm5rundir_emissions1D(ldir, trange=day_range)
        refemis1D = emis_info.emis1D
        #
        #--
        #
        refconc = tm5refdir_load_stationconc(refdir, obsid)
        refconc_list = []
        msg = f"verification based on " \
            f"footprint directory ***{str(topdir)}*** and refdir ***{str(refdir)}***"
        print(msg)
        for iday,day in enumerate(day_range):
            _ostart = day + Timedelta(hours=obs_hr) - Timedelta(seconds=obs_tw)
            _oend   = day + Timedelta(hours=obs_hr) + Timedelta(seconds=obs_tw)
            cnd_day = (refconc['time']>=_ostart)&(refconc['time']<=_oend)
            refconc_day = refconc.loc[cnd_day,'conc'].mean()
            refconc_list.append(refconc_day)
            linconc_day = np.dot(jac2D[iday,:], refemis1D) + iniconc_list[iday]
            msg = f"@{day.strftime('%Y%m%d')}, " \
                f"refconc/linconc/iniconc / obsconc = " \
                f"{refconc_day}/{linconc_day}/{iniconc_list[iday]} / {obs_list[iday]}"
            print(msg)
        try:
            refconc_info = tm5rundir_iniconc_1obs(refdir, obs_info)
            print(refconc_info.conc)
        except FileNotFoundError:
            pass
        msg = f"verification modus, terminating without generating output!"
        logger.info(msg)
        sys.exit(0)
    #
    #--
    #
    trange_tag = f"{dayf.strftime('%Y%m%d')}--{dayl.strftime('%Y%m%d')}"
    outname_tokens = ["fitic-inversion-input", obsid, trange_tag,]
    outname = '_'.join(outname_tokens) + '.nc'
    outname = set_outname(args, outname)
    msg = f"writing inversion inputs to file ***{outname}***..."
    logger.info(msg)

    #
    #--
    #
    fp = Dataset(outname, 'w')
    fp.createDimension('nemis', nemis)
    fp.createDimension('nday', nday)
    fp.createDimension('nsta', len(station_list))

    #-- 
    ncvar = fp.createVariable('emission', 'f8', ('nemis',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = "CH4 emissions"
    ncvar.units = 'kgCH4/cell/s'
    ncvar[:] = emis1D[:]
    #
    ncvar = fp.createVariable('iday', iday1D.dtype, ('nemis',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = 'index_of_emission_day'
    ncvar.units = ''
    ncvar[:] = iday1D[:]
    #
    # ncvar = fp.createVariable('region', region1D.dtype, ('nemis',),
    #                           compression='zlib', complevel=complevel)
    ncvar = fp.createVariable('region', region1D.dtype, ('nemis',))
    ncvar.long_name = f"emission_region_identifier"
    ncvar.units = ''
    ncvar[:] = region1D[:]
    #
    ncvar = fp.createVariable('obs_jacobian', 'f8', ('nday','nemis',),
                              compression='zlib', complevel=complevel)
    ncvar.units = 'ppb/(kgCH4/cell/s)'
    ncvar.comment = f"quantifies change of concentration at selected days w.r.t. emissions."
    ncvar[:] = jac2D[:]
    #
    ncvar = fp.createVariable('obs', 'f8', ('nday',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = f"observed CH4 concentration"
    ncvar.units = 'ppb'
    ncvar[:] = np.array(obs_list)
    #
    ncvar = fp.createVariable('iniconc', 'f8', ('nsta','nday',),
                              compression='zlib', complevel=complevel)
    ncvar.long_name = f"initial_concentration"
    ncvar.units = 'ppb'
    ncvar[0,:] = iniconc_list
    #
    ncvar = fp.createVariable('day', str, ('nday',))
    for iday,day in enumerate(day_range):
        ncvar[iday] = day.strftime('%Y%m%d')
    ncvar.long_name = f"day of emission *and* day of observation."
    ncvar.units = ''
    #
    ncvar = fp.createVariable('station', str, ('nsta',))
    ncvar.long_name = f"station_identifier"
    ncvar.units = ''
    for ista,staid in enumerate(station_list):
        ncvar[ista] = staid
    #-- global attributes
    fp.description = f"Jacobian quantifies sensitivity of concentration at each individual day " \
        f"w.r.t. to emissions from first to last day."
    fp.history = f"{' '.join(sys.argv)}"
    fp.date_created = Timestamp.utcnow().isoformat()
    fp.close()
    msg = f"generated file ***{outname}***"
    logger.info(msg)


################################################################################
#
#                   p a r s e r
#
parser = ArgumentParser(usage=globals()['__doc__'])
parser = ArgumentParser()
#----------------------------
#     s u b c o m m a n d s
#
subparsers = parser.add_subparsers( title='Available Subcommands',
                                    metavar='CMDS',
                                    description='',
                                    dest='subcmds',
                                    help='')
#
#--       test_jacobianfwd_1day
#
sparser = subparsers.add_parser('test_jacobianfwd_1day',
                                help="""compute contribution of emissions to selected observation""")
sparser.add_argument('tm5rundir',
                     help="""top-level directory of TM5 adjoint run for footprint creation.""")
sparser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
sparser.add_argument('--obsid',
                     default='cbw_207',
                     help="""restrict Jacobian to one single observational location (default: %(default)s).""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       testbuild_jacobian_1day
#
sparser = subparsers.add_parser('testbuild_jacobian_1day',
                                help="""test Jacobian generation for one particular observation day.""")
sparser.add_argument('tm5rundir',
                     help="""top-level directory of TM5 adjoint run for footprint creation.""")
sparser.add_argument('--trange',
                    metavar=('tstart','tend'),
                    nargs=2,
                    help="""whether to override simulation start/end time specified in the yaml file (strings must be parseable as pandas Timestamp).""")
sparser.add_argument('--obsid',
                     default='cbw_207',
                     help="""restrict Jacobian to one single observational location (default: %(default)s).""")
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")

#
#--       testbuild_jacobian_period
#
sparser = subparsers.add_parser('testbuild_jacobian_period',
                                help="""test preparation of inputs for Fortran inversion system.""")
sparser.add_argument('outpath_tm5',
                     help="""top-level directory of series of TM5 adjoint runs for footprint creation each of those for one observation day.""")
sparser.add_argument('obsfile',
                     help="""obspack NetCDF file providing CH4 observations at tthe selected station/site.""")
sparser.add_argument('--obsid',
                     default='cbw_207',
                     help="""select one single observational location (default: %(default)s).""")
sparser.add_argument('--selday',
                     default="20210208",
                     help="""last observational day of accumulation period (default: %(default)s).""")
sparser.add_argument('--days',
                     type=int,
                     default=14,
                     help="""number of days backwards of accumulation period (default: %(default)s).""")
# sparser.add_argument('--period',
#                      nargs=2,
#                      metavar=('dayfirst/daylast'),
#                      default=["20210221","20210228",],
#                      help="""alternative specification of accumulation period (default: %(default)s).""")
sparser.add_argument('--refdir',
                     help="""TM5 forward simulation, can be used to verify the Jacobian approach.""") 
sparser.add_argument('--outdir',
                    help="""top-level directory for any generated outputs..""")
sparser.add_argument('--outname',
                    help="""explictly specifed name of output file (might be ignored in case the request yields multiple files).""")



################################################################################
#
#                   p r o g r a m   s t a r t
#
def main(args):

    if args.subcmds=='test_jacobianfwd_1day':
        subcmd_test_jacobianfwd_1day(args)

    if args.subcmds=='testbuild_jacobian_1day':
        subcmd_testbuild_jacobian_1day(args)

    if args.subcmds=='testbuild_jacobian_period':
        subcmd_testbuild_jacobian_period(args)

if __name__ == '__main__':
    import datetime as dtm

    progname = os.path.basename(__file__)

    #-----------------------------
    #          P R O G R A M   S T A R T
    #
    fmt = "%Y-%m-%dT%H:%M:%S.%f"
    ttstart = dtm.datetime.now()
    logger.info(f"{progname}::PROGRAM START::{ttstart.strftime(fmt)}")
    argv = ' '.join(sys.argv)
    logger.info(f"  command-line -->{argv}<--")

    #
    #          p a r s e   c o m m a n d   l i n e
    #
    args = parser.parse_args()

    #
    #--        s t a r t   e x e c u t i o n
    #
    main(args)



# ################################################################################
# #
# #                   p r o g r a m   s t a r t
# #
# args = parser.parse_args(sys.argv[1:])

