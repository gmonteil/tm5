#!/usr/bin/env python2.7

from pyshell.observations import observations, HDF5DB
from pyshell.tmflex.RunTM5 import RunTM5
from pyshell.tmflex.runtools import rcdat
from pyshell.tmflex import eurocom
from pyshell.tmflex.emissions.fluxes_verify import CO2_Emissions
from pyshell.base.main.optimizer import conGrad as congrad

from pandas import Timestamp
import os


def load_observations(rc):
    # Load observations
    backgroundObs = observations()
    backgroundObs.append(HDF5DB(rc['observations']['file']))
    errmin = float(rc['observations'].get('err_min', 0.5))
    err = backgroundObs.data['mixing_ratio_err'] * rc['observations'].get('error_factor', .5)
    err[err < errmin] = errmin
    backgroundObs.data['mixing_ratio_err'][:] = err
    return backgroundObs


def load_rcf(rc):
    """
    Load the TM5 rc-file
    :param rc: omegaconf.DictConfig or dictionary
    :return:
    """
    rcf = rcdat()
    rcf.setkey('my.project', rc['run']['project'])
    rcf.readfile(rc['run']['rcfile'])
    rcf.ti = Timestamp(rc['run']['start'])
    rcf.tf = Timestamp(rc['run']['end'])
    rcf.substituteTimes()
    
    setup_environment(rcf)
    setup_output(rc, rcf)
    
    rcf.setup_meteo_coarsening(rc.meteo.coarsen)

    # Keys under the "tm5" group are needed by TM5 itself (not just pyshell!)
    for k, v in rc.tm5.items():
        rcf.setkey(k, v)

    return rcf


def setup_output(rc, rcf):
    """
    write the output.* keys:
    """
    output_mix = rc['output'].get('mix')
    if output_mix:
        rcf.setkey('output.mix', True)
        rcf.setkey('output.mix.tstep', output_mix['tstep'])
        rcf.setkey('output.mix.filename.prefix', output_mix['filename_prefix'])


def setup_tm5(rc):
    rcf = load_rcf(rc)
    return RunTM5(rcf)


def compile(rc):
    rcf = load_rcf(rc)
    run = RunTM5(rcf, rc['run']['start'], rc['run']['end'])
    run.Compile()
    return run


def setup_environment(rcf):
    os.environ['pyshell.rc'] = rcf.filename
    if rcf.get('par.openmp'):
        nthreads = rcf.get('par.maxthreads')
        os.environ['omp_num_threads'] = '%s'%nthreads
        os.environ['omp_num_threads'.upper()] = '%s'%nthreads
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
    
    
def forward(rc, step=None):
    """
    Run a forward TM5 simulation
    :param rc: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    obs = load_observations(rc)
    emclasses = {'CO2': CO2_Emissions}
    run = setup_tm5(rc)
    run.SetupEmissions(emclasses, step=step)
    run.SetupObservations(obs)
    run.Compile()
    run.RunForward()
    return run

    
def optim(rc):
    """
    Do an inversion with TM5
    :param rc: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    obs = load_observations(rc)
    emclasses = {'CO2': CO2_Emissions}
    run = setup_tm5(rc)
    run.SetupObservations(obs)
    run.Compile()
    
    opt = congrad(run)
    opt.SetupOptimizer(restart=False, optimized_prior_state=False, emclasses=emclasses)
    opt.Var4D()

    return opt 


def compile(rc):
    run = setup_tm5(rc)
    run.Compile()
    return run