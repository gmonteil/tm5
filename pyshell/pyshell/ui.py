#!/usr/bin/env python2.7

from pyshell.runtools import rcdat
from pyshell.emissions import PreprocessedEmissions
from pyshell.optimizer import Congrad
from pyshell.model import RunTM5

from pandas import Timestamp
import os
import shutil
import xarray as xr


def forward(dconf):
    """
    Run a forward TM5 simulation
    :param dconf: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    emclasses = {'CO2': PreprocessedEmissions}
    obs = load_observations(dconf)
    run = setup_tm5(dconf)
    run.SetupObservations(obs)
    run.SetupEmissions(emclasses)
    run.RunForward()
    return run


def optim(dconf):
    """
    Do an inversion with TM5
    :param dconf: omegaconf.DictConfig or dictionary containing basic settings
    More advanced/stable settings are stored in the TM5 rc-file, accessed under the run.rcfile key
    """

    emclasses = {'CO2': PreprocessedEmissions}
    obs = load_observations(dconf)
    run = setup_tm5(dconf)
    run.SetupObservations(obs)
    opt = Congrad(run)
    opt.SetupOptimizer(restart=False, optimized_prior_state=False, emclasses=emclasses)
    opt.Var4D()

    return opt


# def setup_pyshell(dconf, step=None):
#     emclasses = {'CO2': PreprocessedEmissions}
#     obs = load_observations(dconf)
#     run = setup_tm5_legacy(dconf)
#     run.SetupEmissions(emclasses, step=step)
#     run.SetupObservations(obs)
#     run.rcf.setkey('my.runmode', 1)
#     run.rcf.WriteFile(os.path.join(dconf.paths.output, dconf.run.project, 'tm5.rc'))
#     return run


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
            if isinstance(catinfo, str):
                rcf.setkey('%s.%s.routine' % (tracer, cat), catinfo)
            else :
                rcf.setkey('%s.%s.routine' % (tracer, cat), dconf.emissions[tracer].categories[cat].routine)
                if 'dailycycle' in catinfo:
                    rcf.setkey('%s.%s.dailycycle' % (tracer, cat), 'T')

    if 'filename' in dconf['emissions']:
        # Copy the emission file to the run directory:
        shutil.copy(dconf.emissions.filename, rcf.get('PyShell.em.filename'))
        shutil.rmtree(rcf.get('dailycycle.folder'), ignore_errors=True)
        shutil.copytree(dconf.emissions.dailycycle_folder, rcf.get('dailycycle.folder'))

    return rcf


# def load_rcf_legacy(dconf):
#     rcf = rcdat()
#
#     rcf.setkey('my.project', dconf['run']['project'])
# #    rcf.setup_meteo_coarsening(rc.meteo.coarsen)
#     rcf.readfile(dconf['run']['rcfile'])
#     rcf.ti = Timestamp(dconf['run']['start'])
#     rcf.tf = Timestamp(dconf['run']['end'])
#     rcf.substituteTimes()
#
#     rcf.setkey('jobstep.timerange.start', rcf.ti)
#     rcf.setkey('jobstep.timerange.end', rcf.tf)
#     setup_environment(dconf)
#     setup_output(dconf, rcf)
#
#
#     # Keys under the "tm5" group are needed by TM5 itself (not just pyshell!)
#     if 'tm5' in dconf:
#         for k, v in dconf.tm5.items():
#             rcf.setkey(k, v)
#
#     return rcf


# def setup_tm5_legacy(rc):
#     return RunTM5(load_rcf_legacy(rc))


# def compile(rc):
#     rcf = load_rcf(rc)
#     run = RunTM5(rcf, rc['run']['start'], rc['run']['end'])
#     run.Compile()
#     return run


# def setup_tmflex(dconf, rcf):
#     """
#     Perform a forward run with extraction of the background concentration by TM5
    # This requires the following adjustments, compared to a normal forward run:
    # - creation of a new "fg" tracer, for each original tracer (e.g. CO2fg, if we had a CO2 tracer)
    # - appending the "tmflex" project to the list of source files
    # - creation of the emissions and daily_cycle files for the new tracer
    # - setting the "tmflex" rc-keys
    # """
    #
    # if not dconf.get('output'):
    #     return rcf
    # if not dconf.output.get('background', None):
    #     return rcf
    #
    # # Adjust the list of source directories:
    # # We need to adjust it everywhere, because that **** rc module resolves keys at load time
    # #for key in ['my.projects.basic', 'my.source.dirs', 'build.copy.dirs']:
    # #    value =  rcf.get(key)
    # #    value += ' ' + '%s/tmflex'%rcf.get('my.proj.root')
    # #    rcf.setkey(key, value)
    #
    # # Create the new tracer:
    # assert len(dconf.output.background.tracers) == 1, "Rodenbeck scheme available only for one tracer (might work with more, but untested)"
    #
    # # Copy keys from the main tracer to the background one:
    # for tracer in dconf.output.background.tracers:
    #     for region in dconf.run.regions :
    #         dconf.emissions[tracer + 'fg'] = {region : dconf.emissions[tracer][region]}
    #
    #     # For dailycycle, we need to replace the tracer name:
    #     if dconf.emissions[tracer].get('dailycycle'):
    #         dconf.emissions[tracer + 'fg']['dailycycle'] = {
    #             'type': dconf.emissions[tracer].dailycycle.type,
    #             'prefix': dconf.emissions[tracer].dailycycle.prefix.replace(tracer, tracer + 'fg')
    #         }
    #         dconf.emissions[tracer + 'fg']['categories'] = dconf.emissions[tracer].categories
    #
    #     # Fix one key in the obs that needs adjusting ...
    #     rcf.setkey('output.point.%s.minerror' % (tracer + 'fg'), rcf.get('output.point.%s.minerror' % tracer))
    #
    # # Add the tracer(s) to the rc-file
    # dconf.run.tracers.extend([tr + 'fg' for tr in dconf.output.background.tracers])
    # rcf.setkey('my.tracer', ', '.join([_ for _ in dconf.run.tracers]))
    # rcf.setkey('my.tracer.name', ', '.join([_ for _ in dconf.run.tracers]))
    # rcf.setkey('tracers', len(dconf.run.tracers))
    #
    # rcf.setkey('tmflex.compute.backgrounds', True)
    # rcf.setkey('tmflex.lon0', dconf.output.background.lon_range[0])
    # rcf.setkey('tmflex.lon1', dconf.output.background.lon_range[1])
    # rcf.setkey('tmflex.lat0', dconf.output.background.lat_range[0])
    # rcf.setkey('tmflex.lat1', dconf.output.background.lat_range[1])
    #
    # return rcf



