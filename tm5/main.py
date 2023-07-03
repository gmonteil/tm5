#!/usr/bin/env python
import os

from omegaconf import OmegaConf, DictConfig
import tm5.emissions
import tm5.observations
from tm5.build import build_tm5
from tm5.meteo import Meteo
from tm5 import inicond
from tm5.run import run_tm5
from tm5.units import units_registry as ureg
from tm5.settings import TM5Settings
from tm5 import species as chem
from pathlib import Path
from pandas import Timestamp, Timedelta
from loguru import logger


class TM5:
    def __init__(self, dconf: str, machine : str = 'machine'):
        self.dconf = OmegaConf.load(dconf)
        # Copy the section of the yaml file given by the "machine" keyword arguments to the actual "machine" section.
        # If the "machine" keyword argument is not provided, this does nothing (i.e. it copies "machine" section to "machine" section)
        self.dconf['machine'] = self.dconf[machine]

        self.configfile = dconf
        self.settings = TM5Settings()
        self.tm5exec = Path(self.dconf.run.paths.output) / 'tm5.x'
        self.start = Timestamp(self.dconf.run.start)
        self.end = Timestamp(self.dconf.run.end)

        self.meteo = Meteo(**self.dconf.meteo)

        # pyshell-related
        self.machine = machine

    def build(self, clean : bool = False):
        """
        Build TM5
        """
        tm5exec = build_tm5(self.dconf, clean = clean)
        if not self.tm5exec.exists() and not self.tm5exec.is_symlink():
            self.tm5exec.parent.mkdir(parents=True, exist_ok=True)
            os.symlink(tm5exec, self.tm5exec)

    def calc_background(self, lon0 : float, lon1 : float, lat0 : float, lat1 : float, emissions_file : str):
        self.settings['istart'] = '1'
        self.settings['mask.apply'] = 'T'
        self.settings['mask.factor'] = '0'
        self.settings['mask.region'] = f'{lon0:.1f} {lon1:.1f} {lat0:.1f} {lat1:.1f}'
        self.forward(emission_file=emissions_file)

    def coarsen_meteo(self):
        """
        Do a forward run with global1x1 meteo, and no emissions, no initial condition, etc.
        :return:
        """
        self.setup_meteo(coarsen=True)
        self.setup_run('forward')
        self.setup_iniconc('zero')
        self.setup_output(stations=False)
        self.setup_regions()
        self.setup_system()
        self.settings['proces.source'] = 'F'
        rcf = self.settings.write(Path(self.dconf.run.paths.output) / 'forward.rc')
        run_tm5(f'{str(self.tm5exec.absolute())} {str(rcf)}', settings=self.dconf.machine.host)

    def forward(self, emission_file : str = None):
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
        self.setup_emissions(skip_file_creation=emission_file is not None, filename=emission_file)
        self.setup_meteo()
        self.setup_tm5_optim()
        self.setup_run('forward')
        self.setup_output()
        self.setup_regions()
        self.setup_iniconc()
        self.setup_optim()
        self.setup_tracers()
        self.setup_system()
        rcf = self.settings.write(Path(self.dconf.run.paths.output) / 'forward.rc')
        run_tm5(f'{str(self.tm5exec.absolute())} {str(rcf)}', settings=self.dconf.machine.host)

    def optim(self):
        raise NotImplementedError

    def setup_meteo(self, coarsen : bool = False):
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

        write_meteo = 'F'
        if not self.dconf.meteo.coarsened or coarsen:
            self.meteo.coarsened = False
            self.dconf.meteo.coarsened = True
            write_meteo = 'T'
            # All fields are read from glb100x100
            self.settings[f'tmm.sourcekey.*.ml'] = f'tm5-nc:mdir=ec/ea/h06h18tr3/ml137/glb100x100/<yyyy>/<mm>;tres=_00p03;namesep=/'
            self.settings[f'tmm.sourcekey.*.sfc.fc'] = f'tm5-nc:mdir=ec/ea/h06h18tr1/sfc/glb100x100/<yyyy>/<mm>;tres=_00p01;namesep=/'
            self.settings[f'tmm.sourcekey.*.sfc.an'] = f'tm5-nc:mdir=ec/ea/an0tr1/sfc/glb100x100/<yyyy>/<mm>;tres=_00p01;namesep=/'

            self.settings['tmm.output.dir'] = str(Path(self.dconf.meteo.output_path).absolute())
            Path(self.dconf.meteo.output_path).mkdir(exist_ok=True, parents=True)
            self.settings['cf-standard-name-table'] = Path(self.dconf.run.paths.cf_table).absolute()

            for region in self.dconf.regions:
                self.settings[f'tmm.destkey.{region}.ml'] = f'tm5-nc:mdir=ec/ea/h06h18tr3/tropo25/{region}/<yyyy>/<mm>;tres=_00p03;namesep=/'
                self.settings[f'tmm.destkey.{region}.sfc.fc'] = f'tm5-nc:mdir=ec/ea/h06h18tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'
                self.settings[f'tmm.destkey.{region}.sfc.an'] = f'tm5-nc:mdir=ec/ea/an0tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'
            self.settings['ndyn'] = '900'
            self.settings['cfl.outputstep'] = '900'
        else :
            for region in self.dconf.run.regions :
                self.settings[f'tmm.sourcekey.{region}.ml'] = f'tm5-nc:mdir=ec/ea/h06h18tr3/tropo25/{region}/<yyyy>/<mm>;tres=_00p03;namesep=/'
                self.settings[f'tmm.sourcekey.{region}.sfc.fc'] = f'tm5-nc:mdir=ec/ea/h06h18tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'
                self.settings[f'tmm.sourcekey.{region}.sfc.an'] = f'tm5-nc:mdir=ec/ea/an0tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'

        self.settings['my.levs'] = 'tropo25'
        self.settings['cfl.outputstep'] = '3600'
        self.settings['tmm.output'] = write_meteo       # write meteo?
        self.settings['tmm.output.*.*'] = write_meteo   # by default write all fields
        self.settings['tmm.output.*.sfc.const'] = 'F' # Except constant surface fields
        self.settings['tmm.output.glb100x100.sfc.const'] = 'F'

        self.settings['diffusion.dir'] = Path(self.dconf.run.paths.diffusion) / 'dkg'

        # Constant 1x1 fields (oro and lsm):
        self.settings['tmm.sourcekey.*.sfc.const'] = 'tm5-nc:mdir=ec/ea/an0tr1/sfc/glb100x100;tres=_00p01;namesep=/'

        self.meteo.setup_files(start=self.dconf.run.start, end=self.dconf.run.end)

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

    def setup_output(self, stations : bool = True):
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

        if stations and 'stations' in self.dconf.output :
            self.setup_output_stations(self.dconf.output.stations)
            
        self.setup_output_mix(self.dconf.output.get('mix', None))

    def setup_output_stations(self, dconf):
        self.settings['output.station.timeseries'] = 'T'
        self.settings['output.station.timeseries.filename'] = Path(dconf.filename).absolute()

    def setup_output_point(self, dconf):
        self.settings['output.point'] = 'T'
        self.settings['output.point.input.dir'] = dconf.input_dir
        self.settings['output.point.split.period'] = 'a'  # no splitting ...
        self.settings['output.point.sample.parent'] = dconf.get('sample_parent', 'F')
        
    def setup_output_mix(self, dconf: DictConfig | None = None):
        if dconf is None :
            return
        
        self.settings['output.mix'] = 'T'
        self.settings['output.mix.tstep'] = int(Timedelta(dconf.output_frequency).total_seconds())
        self.settings['output.mix.meteo'] = dconf.get('output_meteo', 'F')
        self.settings['output.mix.filename.prefix'] = dconf.prefix
        self.settings['output.mix.deflate.leve'] = dconf.get('deflate_level', 1)

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

    def setup_iniconc(self, ini : str = None):
        """
        This will setup the following group of rc-keys:
        - start.*
        - istart
        """
        if ini is None :
            ini = self.dconf.initial_condition.type

        match ini:
            case 'zero':
                self.settings['istart'] = '1'
            case 'mixfile':
                self.settings['istart'] = '2'
                self.settings['start.2.iniconcfile'] = self.start.strftime(self.dconf.initial_condition.mixfile)
                self.settings['start.2.iniconc_from_file'] = 'T'
            case 'savefile':
                self.settings['istart'] = '3'
                self.settings['start.3.filename'] = self.start.strftime(self.dconf.initial_condition.savefile)
            case 'carbontracker':
                self.settings['istart'] = '2'
                self.settings['start.2.iniconc_from_file'] = 'T'
                version = self.dconf.initial_condition.carbontracker_version
                Path(self.dconf.run.paths.output).mkdir(parents=True, exist_ok=True)
                filename = Path(self.dconf.run.paths.output) / Timestamp(self.dconf.run.start).strftime(f'mix_co2_%Y%m%d_{version}.nc')
                self.settings['start.2.iniconcfile'] = filename
                inicond.get_iniconc_carbontracker(
                    self.dconf.initial_condition.carbontracker_url,
                    self.dconf.run.start,
                    self.dconf.regions,
                    filename
                )
            case other :
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
        self.settings['output.point.timewindow'] = self.dconf.observations.point[tracer].default_assim_window
        self.settings['output.point.interpolation'] = {'linear': 3, 'gridbox': 1, 'slopes': 2}[self.dconf.observations.point.interpolation]
        return tm5.observations.prepare_point_obs(self.dconf.output.point)

    def setup_emissions(self, skip_file_creation: bool = False, filename : str = None):
        """
        This will create the emission file and dailycycle files as well as setup the following rc-keys:
        - PyShell.em.filename
        - dailycycle.folder
        - {tracer}.{cat}.dailycycle
        - {tracer}.dailycycle.type
        - {tracer}.dailycycle.prefix

        Arguments:
            skip_file_creation [bool] : set to True to avoid recomputing the emission files if it is already there
        """

        if not filename :
            filename = Path(self.dconf.run.paths.output) / 'emissions.nc'
        self.settings['PyShell.em.filename'] = str(filename)
        self.settings['dailycycle.folder'] = self.dconf.emissions.dailycycle_folder
        
        for trname in self.dconf.emissions.tracers :
            self.settings[f'{trname}.dailycycle.type'] = self.dconf.emissions[trname].dailycycle.type
            self.settings[f'{trname}.dailycycle.prefix'] = Path(self.dconf.emissions[trname].dailycycle.filename_format).with_suffix('').with_suffix('').name + '.'
            for catname in self.dconf.emissions[trname].emission_categories :
                apply_dailycycle_to_cat = {False:'F', True:'T'}[self.dconf.emissions[trname].emission_categories[catname].get('dailycycle', False)]
                self.settings[f'{trname}.{catname}.dailycycle'] = apply_dailycycle_to_cat
            
        if not skip_file_creation:
            tm5.emissions.prepare_emissions(self.dconf.emissions, filename = filename)

    def setup_optim(self):
        """
        This creates keys that are related to the optimization, and should therefore probably not be in TM5 at all ...
        - emission.{tracer}.{region}.categories
        - emission.{tracer}.{region}.category{icat}
        """
        for tracer in self.dconf.emissions.tracers :
            for region in self.dconf.emissions.regions :
                cats = self.dconf.emissions[tracer].emission_categories
                self.settings[f'emissions.{tracer}.{region}.categories'] = ', '.join([c for c in cats])
                self.settings[f'emissions.{tracer}.{region}.ncats'] = len(cats)
                for icat, cat in enumerate(cats):
                    catdconf = self.dconf.emissions[tracer].emission_categories.get(cat, self.dconf.emissions[tracer])
                    catfreq = catdconf.get('optim_freq', 'D')
                    self.settings[f'emissions.{tracer}.{region}.{cat}'] = {'MS': 'monthly', 'D': 'daily'}[catfreq]
                    # Since these are probably not needed, I just hardcode them ...
                    # self.settings[f'emission.{tracer}.{region}.category{icat+1:.0f}'] = f'{cat}; 100.0 ; 200.0-g ; 0.0-e-monthly ; 0 ; dummy'
# 
    def setup_tracers(self) -> None:
        """
        Setup rc keys required by chem_params.F90:
        - tracers.number
        - tracers.name
        - tracers.{tr}.molar_mass
        - tracers.{tr}.mixrat_unit_value
        - tracers.{tr}.mixrat_unit_name
        - tracers.{tr}.emis_unit_value
        - tracers.{tr}.emis_unit_name

        The "mixrat_unit_value" and "emis_unit_value" keys are computed based on their "name" counterparts, using
        information from the "tm5.units" module ==> new chemical species need to be implemented there first!
        """

        tracers = self.dconf.run.tracers
        self.settings['tracers.number'] = len(tracers)
        self.settings['tracers.names'] = ','.join([_ for _ in tracers])
        for tr in tracers :
            spec = self.dconf.tracers[tr].species
            emis_unit = self.dconf.tracers[tr].get('flux_unit', chem.species[spec].unit_emis)
            mix_unit = self.dconf.tracers[tr].get('mix_unit', chem.species[spec].unit_mix)
            self.settings[f'tracers.{tr}.molar_mass'] = 1 / ureg.Quantity(f'g{spec}').to('mol').m
            self.settings[f'tracers.{tr}.mixrat_unit_value'] = ureg.Quantity('mol / mol').to(mix_unit).m
            self.settings[f'tracers.{tr}.mixrat_unit_name'] = str(mix_unit)
            self.settings[f'tracers.{tr}.emis_unit_value'] = ureg.Quantity(f'kg{spec}').to(emis_unit).m
            self.settings[f'tracers.{tr}.emis_unit_name'] = str(emis_unit)

    def setup_system(self) -> None:
        """
        Setup system-specific settings (paths, etc.). For now the following rc-keys are set:
        - udunits_path
        """
        self.settings['udunits_path'] = Path(self.dconf.machine.paths.udunits).absolute()
