#!/usr/bin/env python

from pathlib import Path
from loguru import logger

"""
Default TM5 settings (rc keys)
"""


defaults = {
    'adjoint.input.point' : 'T',
    'adjoint.input.satellite' : 'F',

    'cfl.outputstep': '${ndyn}',

    'diffusion.files.deflate.level': '1',

    'do_steps.print.mass' : 'F',

    'mask.apply' : 'F',
    'mask.complement': 'T',
    'mask.factor': '0',

    'mass.output.subdir': 'mass',

    'meteo.tinterp.ml': 'interp3',
    'meteo.tinterp.convec': 'aver3',
    'meteo.tinterp.sfc.aver': 'aver1',
    'meteo.tinterp.sfc.inst':  'interp1',
    'meteo.tinterp.sfc.an': 'interp1',
    'meteo.tinterp.oro': 'const',
    'meteo.tinterp.lsm': 'const',
    'meteo.tinterp.srols': 'month',

    'model.output': 'F',

    'my.mlevs': 'tropo25',

    'my.zoom': 'glb6x4',

    'ndyn': '3600',

    'optimize.emission': 'F',

    'output.after.step': 'v',

    'proces.advection': 'T',
    'proces.chemistry': 'F',
    'proces.diffusion.eps_d': '1.0e-4',
    'proces.source': 'T',

    'restart.write' : 'T',
    'restart.correct.mixing.ratio' : 'F',
    'restart.file.contains.padded.mass': 'F',
    'restart.overwrite.mass': 'F',

    'save.output.subdir': '/',

    'timestep.read' : '10800',

    'timing.output.subdir': 'timing',
}


# Values with no default specified:
#   correlation.inputdir
#   diffusion.dir
#   input.dir ==> need cleanup in TM5 code.
#   istart
#   jobstep.timerange.end
#   jobstep.timerange.start
#   mask.region
#   my.meteo.source.dir
#   my.run.dir
#   my.runmode
#   output.dir
#   output.mix, output.flux1x1, output.satellite, output.station, output.totalcol, output.tccon ==> not setup means no output
#   output.point.*
#   PyShell.em.filename
#   start.* ==> do in setup_iniconc
#   tmm.dir ==> do in setup_tmm
#   tmm.output.* ==> do in setup_tmm
#   tmm.sourcekey.* ==> do in setup_tmm

#   {tracer}.{cat}.dailycycle
#   {tracer}.dailycycle_prefix
#   {tracer}.dailycycle.type
#   emission.{tracer}.{region}.categories
#   emission.{tracer}.{region}.category{icat}
#   output.point.errors.{tracer}
#   output.point.{tracer}.minerror
#   output.satellite.meteo.{tracer}
#   tmm.sourcekey.{region}
#   tmm.output.{region}
#   meteo.read.{region}
#   region.{region}.redgrid.nh.n
#   region.{region}.redgrid.sh.n
#   region.{region}.redgrid.nh.comb
#   region.{region}.redgrid.sh.comb


class TM5Settings(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for k, v in defaults.items():
            self[k] = v

    def write(self, filename: Path) -> Path:
        filename.parent.mkdir(exist_ok=True, parents=True)
        with open(filename, 'w') as fid :
            for k, v in self.items():
                fid.write(f'{k} : {str(v)}\n')
        return filename

    def __setitem__(self, key, value):
        logger.info(f'{key:>30s} : {value}')
        super().__setitem__(key, value)