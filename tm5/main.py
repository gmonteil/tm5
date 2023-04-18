#!/usr/bin/env python


from omegaconf import OmegaConf

import tm5.emissions
import tm5.observations
from tm5.build import build_tm5
from tm5.setup import setup_tm5
from tm5.run import run_tm5
from tm5.settings import TM5Settings
import xarray as xr
from pathlib import Path
from pandas import DatetimeIndex, Timestamp
from loguru import logger
import shutil


class TM5:
    def __init__(self, dconf: str):
        self.dconf = OmegaConf.load(dconf)
        self.configfile = dconf
        self.settings = TM5Settings()
        self.tm5exec = Path(self.dconf.run.paths.output) / 'tm5.x'
        self.start = Timestamp(self.dconf.run.start)
        self.end = Timestamp(self.dconf.run.end)

    def build(self):
        """
        Build TM5
        """
        tm5exec = build_tm5(self.dconf)
        shutil.copyfile(tm5exec, self.tm5exec)

    def optim(self):
        """
        Do an inversion ==> this is still based on pyshell
        """
        self.setup()
        run_tm5(f'pyshell optim --rc {self.configfile}', settings=self.dconf.machine.host)

    def forward_pyshell(self):
        """
        Do a forward run, using pyshell
        """
        run_tm5(f'pyshell forward --rc {self.configfile}', settings=self.dconf.machine.host)

    def setup(self):
        self.dconf = setup_tm5(self.dconf)

    def calc_background(self, lon0, lon1, lat0, lat1):
        """
        This should just setup the "mask.apply", "mask.region", "istart" and "PyShell.em.filename" keys
        """

        # Adapt settings for computing the foreground component:
        # - initial condition set to 0
        self.dconf.initial_condition.type = 'zero'
        # - masking of concentrations outside the domain enabled
        #   ==> keys under dconf.pyshell.tm5 are appended to the TM5 rc-file by pyshell itself (ui.load_rcf)
        self.dconf.pyshell.tm5['mask.apply'] = 'T'
        self.dconf.pyshell.tm5['mask.complement'] = 'T'
        self.dconf.pyshell.tm5['mask.factor'] = '0'
        self.dconf.pyshell.tm5['mask.region'] = f'{lon0:.1f} {lon1:.1f} {lat0:.1f} {lat1:.1f}'
        # - Read emissions optimized in the inversion step
        #   ==> keys under dconf.pyshell2 are appended to the rc/lumia.rc rc-file, which needs to be included in the
        #       standard
        self.dconf.pyshell2 = self.dconf.get('pyshell2', {})  # Ensure that the node exists
        self.dconf.pyshell2['emission.read.optimized'] = 'T'
        self.dconf.pyshell2['emission.read.optimized.filename'] = 'emission.nc4'

        # Setup the forward run again (ensure that the keys defined above are written to rc-file
        self.setup()

        # Run the inversion
        self.forward_pyshell()

        # Retrieve the results

    def forward(self):
        """
        Do a forward run, bypassing totally the pyshell
        """

        # Files that need to be created:
        # - emissions + daily cycles
        # - observations
        # - rc-file

        # Files that need to be linked:
        # - initial condition
        # - meteo
        # - other inputs (station files, etc.)

        # TM5 executable generated in the build step

        self.setup_observations()
        self.setup_emissions()
        self.setup_meteo()
        self.setup_tm5_optim()
        self.setup_run('forward')
        self.setup_output()
        self.setup_regions()
        self.setup_iniconc()
        self.setup_optim()
        rcf = self.settings.write(Path(self.dconf.run.paths.output) / 'forward.rc')
        run_tm5(f'{str(self.tm5exec.absolute())} {str(rcf)}', settings=self.dconf.machine.host)

    def setup_meteo(self):
        """
        This will set the following (group of) rc keys:
        - my.meteo.source.dir
        - my.levs
        - tmm.dir
        - tmm.output
        - tmm.sourcekey.{region}.*
        - tmm.output.{region}.*
        - # meteo.read.{region}.* ==> not used in 4dvar yet
        - diffusion.dir
        """
        self.settings['my.meteo.source.dir'] = Path(self.dconf.run.paths.meteo).absolute()
        self.settings['tmm.dir'] = Path(self.dconf.run.paths.meteo).absolute()
        self.settings['my.levs'] = 'tropo25'
        self.settings['cfl.outputstep'] = '3600'
        write_meteo = 'F'
        self.settings['tmm.output'] = write_meteo       # write meteo?
        self.settings['tmm.output.*.*'] = write_meteo   # by default write all fields
        self.settings['tmm.output.*.sfc.const'] = 'F' # Except constant surface fields
        self.settings['tmm.output.glb100x100.sfc.const'] = write_meteo # this was forced to "F" originally ...

        self.settings['diffusion.dir'] = Path(self.dconf.run.paths.scratch) / 'dkg'
        for region in self.dconf.run.regions :
            if region != 'glb600x400':
                raise NotImplementedError
            self.settings[f'tmm.sourcekey.{region}.ml'] = f'tm5-nc:mdir=ec/ea/h06h18tr3/tropo25/{region}/<yyyy>/<mm>;tres=_00p03;namesep=/'
            self.settings[f'tmm.sourcekey.{region}.sfc.fc'] = f'tm5-nc:mdir=ec/ea/h06h18tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'
            self.settings[f'tmm.sourcekey.{region}.sfc.an'] = f'tm5-nc:mdir=ec/ea/an0tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'

        # Constant 1x1 fields (oro and lsm):
        self.settings['tmm.sourcekey.*.sfc.const'] = 'tm5-nc:mdir=ec/ea/an0tr1/sfc/glb100x100;tres=_00p01;namesep=/'

    def setup_regions(self):
        """
        This will set the following (group of) rc keys:
        - region.{region}.redgrid.nh.n
        - region.{region}.redgrid.sh.n
        - region.{region}.redgrid.nh.comb
        - region.{region}.redgrid.sh.comb
        """
        self.settings['region.glb600x400.redgrid.nh.n'] = '3'
        self.settings['region.glb600x400.redgrid.nh.comb'] = '60 20 10'
        self.settings['region.glb600x400.redgrid.sh.n'] = '4'
        self.settings['region.glb600x400.redgrid.sh.comb'] = '60 20 10 5'

    def setup_tm5_optim(self):
        """
        This will set the following keys, needed but which should not be needed:
        - correlation.inputdir
        - var4d.horcor.min_eigval
        - var4d.optim_emis.type
        """
        self.settings['var4d.optim_emis.type'] = '1'
        self.settings['var4d.horcor.min_eigval'] = '0.0001'
        self.settings['correlation.inputdir'] = 'not-defined'

    def setup_output(self):
        """
        This will setup the following (group of) rc keys:
        - output.dir
        - output.mix
        - output.flux1x1
        - output.satellite
        - output.station
        - output.totalcol
        - output.tccon
        - output.satellite.meteo.{tracer}
        """
        self.settings['output.dir'] = self.dconf.run.paths.output
        if 'output' not in self.dconf:
            return

        if 'stations' in self.dconf.output :
            self.setup_output_stations(self.dconf.output.stations)

    def setup_output_stations(self, dconf):
        self.settings['output.station.timeseries'] = 'T'
        self.settings['output.station.timeseries.filename'] = Path(dconf.filename).absolute()

    def setup_output_point(self, dconf):
        self.settings['output.point'] = 'T'
        self.settings['output.point.input.dir'] = dconf.input_dir
        self.settings['output.point.split.period'] = 'a'  # no splitting ...
        self.settings['output.point.sample.parent'] = dconf.get('sample_parent', 'F')

    def setup_run(self, mode='forward'):
        """
        This will setup the following rc-keys:
        - input.dir ==> this is not the input dir, but the location where point input is
        - jobstep.timerange.end
        - jobstep.timerange.start
        - my.run.dir
        """
        self.settings['input.dir'] = Path(self.dconf.run.paths.output) / 'input'
        self.settings['jobstep.timerange.start'] = self.start.strftime('%Y-%m-%d %H:%M:%S')
        self.settings['jobstep.timerange.end'] = self.end.strftime('%Y-%m-%d %H:%M:%S')
        self.settings['my.runmode'] = {'forward': 1, 'adjoin': 2}[mode]

    def setup_iniconc(self):
        """
        This will setup the following group of rc-keys:
        - start.*
        - istart
        """
        if self.dconf.initial_condition.type == 'zero':
            self.settings['istart'] = '1'
        elif self.dconf.initial_condition.type == 'mixfile':
            self.settings['istart'] = '2'
            self.settings['start.2.iniconcfile'] = self.start.strftime(self.dconf.initial_condition.mixfile)
            self.settings['start.2.iniconc_from_file'] = 'T'
        elif self.dconf.initial_condition.type == 'savefile':
            self.settings['istart'] = '3'
            self.settings['start.3.filename'] = self.start.strftime(self.dconf.initial_condition.savefile)
        else :
            logger.error("initial condition settings not understood")
            raise Exception

    def setup_observations(self) -> Path:
        """
        Write a (point) observations file for TM5 + setup the relevant rc-keys:
        - output.point
        - output.point.errors
        - output.point.{tracer}.minerror
        """

        self.setup_output_point(self.dconf.output.point)
        self.settings['output.point.errors'] = self.dconf.observations.point.get('errors', '1')
        for tracer in self.dconf.run.tracers:
            self.settings[f'output.point.{tracer}.minerror'] = self.dconf.observations.point[tracer].minerror
        return tm5.observations.prepare_point_obs(self.dconf.output.point)

    def setup_emissions(self):
        """
        This will create the emission file and dailycycle files as well as setup the following rc-keys:
        - {tracer}.{cat}.dailycycle
        - {tracer}.dailycycle.type
        - PyShell.em.filename
        """

        filename = Path(self.dconf.run.paths.output) / 'emissions.nc'
        self.settings['PyShell.em.filename'] = str(filename)
        for tracer in self.dconf.emissions.tracers :
            pfx = Path(self.dconf.emissions.tracers[tracer].dailycycle_filename_format)
            self.settings[f'{tracer}.dailycycle.type'] = self.dconf.emissions.tracers[tracer].dailycycle_type
            self.settings[f'{tracer}.dailycycle.prefix'] = pfx.with_suffix('').with_suffix('').name + '.'
            self.settings['dailycycle.folder'] = pfx.parent.parent.parent
            for cat in self.dconf.emissions.tracers[tracer].categories:
                self.settings[f'{tracer}.{cat}.dailycycle'] = 'T'

        tm5.emissions.prepare_emissions(self.dconf.emissions, filename = filename)

    def setup_optim(self):
        """
        This creates keys that are related to the optimization, and should therefore probably not be in TM5 at all ...
        - emission.{tracer}.{region}.categories
        - emission.{tracer}.{region}.category{icat}
        """
        for tracer in self.dconf.emissions.tracers :
            for region in self.dconf.emissions.regions :
                cats = self.dconf.emissions[region][tracer].keys()
                self.settings[f'emissions.{tracer}.{region}.categories'] = ', '.join([c for c in cats])
                self.settings[f'emissions.{tracer}.{region}.ncats'] = len(cats)
                for icat, cat in enumerate(cats):
                    catfreq = self.dconf.emissions[region][tracer][cat].optim_freq
                    self.settings[f'emissions.{tracer}.{region}.{cat}'] = {'MS': 'monthly', 'D': 'daily'}[catfreq]
                    # Since these are probably not needed, I just hardcode them ...
                    # self.settings[f'emission.{tracer}.{region}.category{icat+1:.0f}'] = f'{cat}; 100.0 ; 200.0-g ; 0.0-e-monthly ; 0 ; dummy'
