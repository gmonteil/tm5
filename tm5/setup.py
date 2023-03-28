#!/usr/bin/env python

from omegaconf import DictConfig
from pathlib import Path
from pandas import Timestamp


def setup_initial_condition(dconf: DictConfig) -> None:
    fname = Timestamp(dconf.run.start).strftime(dconf.initial_condition.filename_pattern)
    dconf.tm5['istart'] = '3'
    dconf.tm5['start.3.filename'] = fname
    
    
def setup_meteo(dconf: DictConfig) -> None:
    dconf.tm5['my.meteo.resol'] = dconf.meteo.resolution
    if dconf.meteo.output :
        dconf.tm5['my.meteo.resol'] = 'glb100x100'
        dconf.tm5['my.tmm.output'] = 'T'
    else :
        dconf.tm5['my.tmm.output'] = 'F'
    
    
def setup_tm5(dconf: DictConfig) -> Path:
    """
    Generate a simplified rcfile for TM5
    the file should be integrated to the main rc-file via an "#include" statement
    """
    
    if dconf.get('tm5') is None :
        dconf.tm5 = {}
    setup_initial_condition(dconf)
    setup_meteo(dconf)
    
    with open(dconf.run.rcfile, 'w') as fid:
        for k, v in sorted(dconf.tm5.items()):
            fid.write(f'{k:<30s} : {v}\n')
