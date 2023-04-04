#!/usr/bin/env python2.7

from pyshell.observations import observations
from pyshell.tmflex.RunTM5 import RunTM5
from pyshell.tmflex.runtools import rcdat
from pyshell.tmflex.emissions.fluxes_verify import PreprocessedEmissions 
from pyshell.base.main.optimizer import conGrad as congrad
from pyshell.tmflex.emissions.emissions import tm5Emis, tmflexEmis

from pandas import Timestamp
import os
import shutil
import xarray as xr


def load_observations(rc):
    # Load observations
    obs = observations.from_xr_nc(rc['observations']['filename'])
    obs.data2.load()
    errmin = float(rc['observations'].get('err_min', 0.5))
    err = obs.data2.err_obs.values
    err[err < errmin] = errmin
    err *= rc['observations'].get('error_factor', 1.)
    obs.data2.err_obs.values[:] = err
    return obs.data2
    
    backgroundObs = observations()
    backgroundObs.append(xr.open_mfdataset(rc['observations']['filename'])) # HDF5DB(rc['observations']['filename']))
    errmin = float(rc['observations'].get('err_min', 0.5))
    err = backgroundObs.data['mixing_ratio_err'] * rc['observations'].get('error_factor', .5)
    err[err < errmin] = errmin
    backgroundObs.data['mixing_ratio_err'][:] = err
    return backgroundObs


def setup_emissions(rc, rcf):
    for tracer in rc.run.tracers :
        print(tracer)
        for region in rc.run.regions :
            print(region)
            rcf.setkey('emission.%s.%s.categories'%(tracer, region), len(rc.emissions[tracer][region]))
            for icat, cat in enumerate(rc.emissions[tracer][region]):
                print(cat)
                # If we don't optimize the category, then it can be specified as just a single key, specifying the temporal resolution
                # (needed for the daily_cycle files). So fill-in dummy values for the rest
                catinfo = rc.emissions[tracer][region][cat]
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
        if rc.emissions[tracer].get('dailycycle'):
            rcf.setkey("%s.dailycycle.type"%tracer, rc.emissions[tracer].dailycycle.type)
            rcf.setkey("%s.emission.dailycycle"%tracer, 'T')
            rcf.setkey('%s.dailycycle.prefix'%tracer, rc.emissions[tracer].dailycycle.prefix)

        for cat in rc.emissions[tracer].categories :
            catinfo = rc.emissions[tracer].categories[cat]
            if isinstance(catinfo, str):
                rcf.setkey('%s.%s.routine'%(tracer, cat), catinfo)
            else :
                rcf.setkey('%s.%s.routine'%(tracer, cat), rc.emissions[tracer].categories[cat].routine)
                if 'dailycycle' in catinfo:
                    rcf.setkey('%s.%s.dailycycle'%(tracer, cat), 'T')

    if 'filename' in rc['emissions']:
        # Copy the emission file to the run directory:
        shutil.copy(rc.emissions.filename, rcf.get('PyShell.em.filename'))
        shutil.rmtree(rcf.get('dailycycle.folder'), ignore_errors=True)
        shutil.copytree(rc.emissions.dailycycle_folder, rcf.get('dailycycle.folder'))

    return rcf


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

    rcf = setup_tmflex(dconf, rcf)
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


def load_rcf_legacy(dconf):
    rcf = rcdat()

    rcf.setkey('my.project', dconf['run']['project'])
#    rcf.setup_meteo_coarsening(rc.meteo.coarsen)
    rcf.readfile(dconf['run']['rcfile'])
    rcf.ti = Timestamp(dconf['run']['start'])
    rcf.tf = Timestamp(dconf['run']['end'])
    rcf.substituteTimes()

    rcf.setkey('jobstep.timerange.start', rcf.ti)
    rcf.setkey('jobstep.timerange.end', rcf.tf)
    setup_environment(dconf)
    setup_output(dconf, rcf)


    # Keys under the "tm5" group are needed by TM5 itself (not just pyshell!)
    if 'tm5' in dconf:
        for k, v in dconf.tm5.items():
            rcf.setkey(k, v)

    return rcf


def setup_output(rc, rcf):
    """
    write the output.* keys:
    """
    if not rc.get('output'):
        return
    output_mix = rc['output'].get('mix')
    if output_mix:
        rcf.setkey('output.mix', True)
        rcf.setkey('output.mix.tstep', output_mix['tstep'])
        rcf.setkey('output.mix.filename.prefix', output_mix['filename_prefix'])


def setup_tm5(dconf):
    rcf = load_rcf(dconf)
    return RunTM5(rcf, dconf)


def setup_tm5_legacy(rc):
    return RunTM5(load_rcf_legacy(rc))


def compile(rc):
    rcf = load_rcf(rc)
    run = RunTM5(rcf, rc['run']['start'], rc['run']['end'])
    run.Compile()
    return run


def setup_environment(dconf):
    os.environ['pyshell.rc'] = dconf.run.rcfile
    if dconf['environment'].get('openmp'):
        nthreads = dconf.environment.openmp.nthreads
        os.environ['omp_num_threads'] = '%s'%nthreads
        os.environ['omp_num_threads'.upper()] = '%s'%nthreads
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
    

def forward(dconf):
    """
    Run a forward TM5 simulation
    :param rc: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    obs = load_observations(dconf)
    run = setup_tm5(dconf)
    run.SetupObservations(obs)
    if dconf.output.get('background'):
        run.SetupEmissions(emclasses={'CO2': tm5Emis, 'CO2fg': tmflexEmis})
    run.rcf.setkey('my.runmode', 1)
    run.rcf.WriteFile(os.path.join(dconf.paths.output, dconf.run.project, 'tm5.rc'))
    return run

    # run.Compile()
    # run.RunForward()
    # return run


def setup_tmflex(dconf, rcf):
    """
    Perform a forward run with extraction of the background concentration by TM5
    This requires the following adjustments, compared to a normal forward run:
    - creation of a new "fg" tracer, for each original tracer (e.g. CO2fg, if we had a CO2 tracer)
    - appending the "tmflex" project to the list of source files
    - creation of the emissions and daily_cycle files for the new tracer
    - setting the "tmflex" rc-keys
    """

    if not dconf.get('output'):
        return rcf
    if not dconf.output.get('background', None):
        return rcf

    # Adjust the list of source directories:
    # We need to adjust it everywhere, because that **** rc module resolves keys at load time
    #for key in ['my.projects.basic', 'my.source.dirs', 'build.copy.dirs']:
    #    value =  rcf.get(key)
    #    value += ' ' + '%s/tmflex'%rcf.get('my.proj.root')
    #    rcf.setkey(key, value)

    # Create the new tracer:
    assert len(dconf.output.background.tracers) == 1, "Rodenbeck scheme available only for one tracer (might work with more, but untested)"

    # Copy keys from the main tracer to the background one:
    for tracer in dconf.output.background.tracers:
        for region in dconf.run.regions :
            dconf.emissions[tracer + 'fg'] = {region : dconf.emissions[tracer][region]}

        # For dailycycle, we need to replace the tracer name:
        if dconf.emissions[tracer].get('dailycycle'):
            dconf.emissions[tracer + 'fg']['dailycycle'] = {
                'type': dconf.emissions[tracer].dailycycle.type,
                'prefix': dconf.emissions[tracer].dailycycle.prefix.replace(tracer, tracer + 'fg')
            }
            dconf.emissions[tracer + 'fg']['categories'] = dconf.emissions[tracer].categories

        # Fix one key in the obs that needs adjusting ...
        rcf.setkey('output.point.%s.minerror' % (tracer + 'fg'), rcf.get('output.point.%s.minerror' % tracer))

    # Add the tracer(s) to the rc-file
    dconf.run.tracers.extend([tr + 'fg' for tr in dconf.output.background.tracers])
    rcf.setkey('my.tracer', ', '.join([_ for _ in dconf.run.tracers]))
    rcf.setkey('my.tracer.name', ', '.join([_ for _ in dconf.run.tracers]))
    rcf.setkey('tracers', len(dconf.run.tracers))

    rcf.setkey('tmflex.compute.backgrounds', True)
    rcf.setkey('tmflex.lon0', dconf.output.background.lon_range[0])
    rcf.setkey('tmflex.lon1', dconf.output.background.lon_range[1])
    rcf.setkey('tmflex.lat0', dconf.output.background.lat_range[0])
    rcf.setkey('tmflex.lat1', dconf.output.background.lat_range[1])

    return rcf


def forward_legacy(rc, step = None):
    """
    Run a forward TM5 simulation
    :param rc: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    emclasses = {'CO2': PreprocessedEmissions}
    run = setup_tm5_legacy(rc)
    run.SetupEmissions(emclasses, step=step)
    obs = load_observations(rc)
    run.SetupObservations(obs)
    run.Compile()
    run.RunForward()
    return run


def setup_pyshell(rc, step=None):
    emclasses = {'CO2': PreprocessedEmissions}
    run = setup_tm5_legacy(rc)
    run.SetupEmissions(emclasses, step=step)
    run.rcf = setup_tmflex(rc, run.rcf)
    obs = load_observations(rc)
    run.SetupObservations(obs)
    run.rcf.setkey('my.runmode', 1)
    run.rcf.WriteFile(os.path.join(rc.paths.output, rc.run.project, 'tm5.rc'))
    return run

    
def compile(rc):
    run = setup_tm5(rc)
    run.Compile()
    return run


# Main entry points:

def optim(dconf):
    """
    Do an inversion with TM5
    :param dconf: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    obs = load_observations(dconf)
    emclasses = {'CO2': PreprocessedEmissions}
    run = setup_tm5(dconf)
    run.SetupObservations(obs)
    opt = congrad(run)
    opt.SetupOptimizer(restart=False, optimized_prior_state=False, emclasses=emclasses)
    opt.Var4D()

    return opt


def forward(dconf):
    """
    Run a forward TM5 simulation
    :param rc: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    obs = load_observations(dconf)
    emclasses = {'CO2': PreprocessedEmissions}
    run = setup_tm5(dconf)
    run.SetupObservations(obs)
    
    # if rc.output.get('background'):
    #     run.SetupEmissions(emclasses={'CO2': tm5Emis, 'CO2fg': tmflexEmis})
    # run.rcf.setkey('my.runmode', 1)
    # run.rcf.WriteFile(os.path.join(dconf.paths.output, dconf.run.project, 'tm5.rc'))

#    from pyshell.tmflex.emissions.emissions import tm5Emis, tmflexEmis
#    emclasses = {'CO2': tm5Emis, 'CO2fg': tmflexEmis}
#    run.rcf.setkey('emissions.filename', 'output/coco2_wp4_ct3/output/2018120100-2020010100/emission.nc4')
#    run.rcf.setkey('dailycycle.folder', os.path.abspath('output/coco2_wp4_ct3bg/dailycycle'))
#    run.rcf.setkey('CO2fg.dailycycle.prefix', 'CO2fg_dailycycle.')
    run.SetupEmissions(emclasses)
    run.RunForward()
    return run

    # run.Compile()
    # run.RunForward()
    # return run
