#!/usr/bin/env python

from omegaconf import DictConfig
from pathlib import Path
from pandas import Timestamp


def setup_initial_condition(dconf: DictConfig) -> DictConfig:
    match dconf['initial_condition'].get('type', 'mixfile'):
        case 'mixfile':
            fname = Timestamp(dconf.run.start).strftime(dconf.initial_condition.mixfile)
            dconf.tm5['istart'] = '2'
            dconf.tm5['start.2.iniconcfile'] = Path(fname).absolute()
            dconf.tm5['start.2.iniconc_from_file'] = 'T'
        case 'savefile':
            fname = Timestamp(dconf.run.start).strftime(dconf.initial_condition.savefile)
            dconf.tm5['istart'] = '3'
            dconf.tm5['start.3.filename'] = Path(fname).absolute()
        case 'zero':
            dconf.tm5['istart'] = '1'

    return dconf
    
    
def setup_meteo(dconf: DictConfig) -> DictConfig:
    dconf.tm5['my.meteo.resol'] = dconf.meteo.resolution
    dconf.tm5['my.meteo.source.dir'] = Path(dconf.run.paths.meteo).absolute()
    if dconf.meteo.output :
        dconf.tm5['my.meteo.resol'] = 'glb100x100'
        dconf.tm5['my.tmm.output'] = 'T'
    else :
        dconf.tm5['my.tmm.output'] = 'F'
    return dconf

        
def setup_paths(dconf: DictConfig) -> DictConfig:
    dconf.pyshell2['my.run.dir'] = Path(dconf.run.paths.output).absolute()
    dconf.pyshell2['my.scratch_dir'] = Path(dconf.machine.paths.scratch).absolute()
    dconf.pyshell2['my.input.dir'] = Path(dconf.machine.paths.input).absolute()
    dconf.pyshell2['pyshell2.build_directory'] = Path(dconf.build.directory).absolute()
    return dconf
    
    
def setup_tm5(dconf: DictConfig) -> DictConfig:
    """
    Generate a simplified rcfile for TM5
    the file should be integrated to the main rc-file via an "#include" statement
    """
    
    if dconf.get('tm5') is None :
        dconf.tm5 = {}
    if not dconf.get('pyshell2'):
        dconf.pyshell2 = {}
    dconf = setup_initial_condition(dconf)
    dconf = setup_meteo(dconf)
    dconf = setup_paths(dconf)
    
    with open(dconf.run.rcfile, 'w') as fid:
        # Keys needed by TM5 itself
        fid.write('!---------- tm5 --------- \n')
        for k, v in sorted(dconf.tm5.items()):
            fid.write(f'{k:<30s} : {v}\n')
            
        # Keys needed by pyshell (but not TM5 ==> should be deprecated, eventually)
        fid.write('\n\n!---------- pyshell --------- \n')
        for k, v in sorted(dconf.pyshell2.items()):
            fid.write(f'{k:<30s} : {v}\n')

    return dconf