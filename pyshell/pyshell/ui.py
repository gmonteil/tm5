#!/usr/bin/env python2.7

from pyshell.runtools import rcdat
from pyshell.emissions import PreprocessedEmissions
from pyshell.optimizer import Congrad
from pyshell.model import RunTM5
from pandas import Timestamp
import os
import shutil
import xarray as xr
import logging


logger = logging.getLogger(__name__)


def forward(dconf, **kwargs):
    """
    Run a forward TM5 simulation
    :param dconf: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    emclasses = {'CO2': PreprocessedEmissions}
    obs = load_observations(dconf)
    run = setup_tm5(dconf, **kwargs)
    run.SetupObservations(obs)
    run.SetupEmissions(emclasses)
    run.RunForward()
    return run


def optim(dconf, **kwargs):
    """
    Do an inversion with TM5
    :param dconf: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """
    
    emclasses = {'CO2': PreprocessedEmissions}
    obs = load_observations(dconf)
    run = setup_tm5(dconf, **kwargs)
    run.SetupObservations(obs)
    opt = Congrad(run)
    opt.SetupOptimizer(restart=False, optimized_prior_state=False, emclasses=emclasses)
    opt.Var4D()

    return opt


def load_observations(dconf):
    # Load observations
    obs = xr.open_mfdataset(dconf['observations']['filename'])
    obs.load()
    errmin = float(dconf['observations'].get('err_min', 0.5))
    err = obs.err_obs.values
    err[err < errmin] = errmin
    err *= dconf['observations'].get('error_factor', 1.)
    obs.err_obs.values[:] = err
    return obs


def setup_tm5(dconf):
    rcf = load_rcf(dconf)
    return RunTM5(rcf, dconf)


def load_rcf(dconf):
    """
    Load the TM5 rc-file
    :param rc: omegaconf.DictConfig or dictionary
    :return:
    """
    rcf = rcdat()

    rcf.setkey('my.project', dconf['run']['project'])
    #    rcf.setup_meteo_coarsening(rc.meteo.coarsen)
    rcf.readfile(dconf['run']['rcfile'])
    rcf.ti = Timestamp(dconf['run']['start'])
    rcf.tf = Timestamp(dconf['run']['end'])
    rcf.filename = dconf.run.rcfile
    rcf.substituteTimes()
    rcf.setkey('jobstep.timerange.start', rcf.ti)
    rcf.setkey('jobstep.timerange.end', rcf.tf)
    rcf.setkey('my.run.dir', os.path.abspath(dconf.run.paths.output))
    setup_environment(dconf)
    setup_output(dconf, rcf)

    # Keys under the "tm5" group are needed by TM5 itself (not just pyshell!)
    if 'tm5' in dconf:
        for k, v in dconf.tm5.items():
            rcf.setkey(k, v)

    if not os.path.exists(dconf.run.paths.output):
        os.makedirs(dconf.run.paths.output)
    rcf = setup_emissions(dconf, rcf)

    return rcf


def setup_environment(dconf):
    os.environ['pyshell.rc'] = dconf.run.rcfile
    if dconf['environment'].get('openmp'):
        nthreads = dconf.environment.openmp.nthreads
        os.environ['omp_num_threads'] = '%s'%nthreads
        os.environ['omp_num_threads'.upper()] = '%s'%nthreads
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'


def setup_output(dconf, rcf):
    """
    write the output.* keys:
    """
    if not dconf.get('output'):
        return
    output_mix = dconf['output'].get('mix')
    if output_mix:
        rcf.setkey('output.mix', True)
        rcf.setkey('output.mix.tstep', output_mix['tstep'])
        rcf.setkey('output.mix.filename.prefix', output_mix['filename_prefix'])


def setup_emissions(dconf, rcf):
    for tracer in dconf.run.tracers :
        for region in dconf.run.regions :
            rcf.setkey('emission.%s.%s.categories' % (tracer, region), len(dconf.emissions[tracer][region]))
            for icat, cat in enumerate(dconf.emissions[tracer][region]):
                # If we don't optimize the category, then it can be specified as just a single key, specifying the temporal resolution
                # (needed for the daily_cycle files). So fill-in dummy values for the rest
                catinfo = dconf.emissions[tracer][region][cat]
                if isinstance(catinfo, str):
                    rcf.setkey('emission.%s.%s.category%i'%(tracer, region, icat + 1), "%s ; 0.0 ; 0000.0-e ; 0.0-e-%s ; 0 ; def-def-0"%(cat, catinfo))

                # Else, the correlation lengths should also be specified:
                else :
                    rcf.setkey('emission.%s.%s.category%i'%(tracer, region, icat + 1), "%s ; %.1f ; %8s ; %s ; %i ; %s"%(
                        cat,
                        catinfo.uncertainty,
                        catinfo.spatial_correlation,
                        catinfo.temporal_correlation,
                        catinfo.optimize,
                        catinfo.type
                    ))
        if dconf.emissions[tracer].get('dailycycle'):
            rcf.setkey("%s.dailycycle.type" % tracer, dconf.emissions[tracer].dailycycle.type)
            rcf.setkey("%s.emission.dailycycle" % tracer, 'T')
            rcf.setkey('%s.dailycycle.prefix' % tracer, dconf.emissions[tracer].dailycycle.prefix)

        for cat in dconf.emissions[tracer].categories :
            catinfo = dconf.emissions[tracer].categories[cat]
            #if isinstance(catinfo, str):
            #     rcf.setkey('%s.%s.routine' % (tracer, cat), catinfo)
            # else :
            #     rcf.setkey('%s.%s.routine' % (tracer, cat), dconf.emissions[tracer].categories[cat].routine)
            if catinfo is not None:
                if 'dailycycle' in catinfo:
                    rcf.setkey('%s.%s.dailycycle' % (tracer, cat), 'T')

    if 'filename' in dconf['emissions']:
        # Copy the emission file to the run directory:
        shutil.copy(dconf.emissions.filename, rcf.get('PyShell.em.filename'))
        shutil.rmtree(rcf.get('dailycycle.folder'), ignore_errors=True)
        shutil.copytree(dconf.emissions.dailycycle_folder, rcf.get('dailycycle.folder'))

    return rcf
