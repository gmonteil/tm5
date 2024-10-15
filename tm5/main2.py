#!/usr/bin/env python
import os

from omegaconf import DictConfig, OmegaConf
import tm5.emissions
import tm5.observations
from tm5.build import build_tm5
from tm5.meteo import Meteo
from tm5 import inicond
from tm5.units import units_registry as ureg
from tm5.settings import TM5Settings, load_config
from tm5 import species as chem
from pathlib import Path
from pandas import Timestamp, Timedelta
from loguru import logger



"""
Re-implementation of whatever was in main.py, because I forgot how up to date that thing was
"""

class TM5:
    def __init__(self, dconf : str | DictConfig, host : str | None) -> None:
        # Load the config file
        self.dconf = load_config(dconf, host)
        self.settings = TM5Settings()
        self.tm5exec = Path(self.dconf.run.paths.output) / 'tm5.x'
        self.meteo = Meteo(**self.dconf.meteo)

    @property
    def start(self) -> Timestamp:
        return Timestamp(self.dconf.run.start)
    
    @property
    def end(self) -> Timestamp:
        return Timestamp(self.dconf.run.end)

    def build(self, clean : bool = False):
        """
        Build TM5
        """
        tm5exec = build_tm5(self.dconf, clean = clean)
        if not self.tm5exec.exists() and not self.tm5exec.is_symlink():
            self.tm5exec.parent.mkdir(parents=True, exist_ok=True)
            os.symlink(tm5exec.absolute(), self.tm5exec)

    # Main setup methods

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


        # Three cases:
        # - use fine meteo; write coarse meteo
        # - use fine meteo
        # - use coarse meteo

        write_meteo = self.dconf.meteo.output
        if not self.dconf.meteo.coarsened or coarsen:
            self.meteo.coarsened = False
            self.dconf.meteo.coarsened = True
            # All fields are read from glb100x100
            self.settings[f'tmm.sourcekey.*.ml'] = f'tm5-nc:mdir=ec/ea/h06h18tr3/ml137/glb100x100/<yyyy>/<mm>;tres=_00p03;namesep=/'
            self.settings[f'tmm.sourcekey.*.sfc.fc'] = f'tm5-nc:mdir=ec/ea/h06h18tr1/sfc/glb100x100/<yyyy>/<mm>;tres=_00p01;namesep=/'
            self.settings[f'tmm.sourcekey.*.sfc.an'] = f'tm5-nc:mdir=ec/ea/an0tr1/sfc/glb100x100/<yyyy>/<mm>;tres=_00p01;namesep=/'
            for region in self.dconf.regions:
                levels = self.dconf.regions[region].levels

            if write_meteo:
                self.settings['tmm.output.dir'] = str(Path(self.dconf.meteo.output_path).absolute())
                Path(self.dconf.meteo.output_path).mkdir(exist_ok=True, parents=True)
                self.settings['cf-standard-name-table'] = Path(self.dconf.run.paths.cf_table).absolute()

                for region in self.dconf.regions:
                    self.settings[f'tmm.destkey.{region}.ml'] = f'tm5-nc:mdir=ec/ea/h06h18tr3/{levels}/{region}/<yyyy>/<mm>;tres=_00p03;namesep=/'
                    self.settings[f'tmm.destkey.{region}.sfc.fc'] = f'tm5-nc:mdir=ec/ea/h06h18tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'
                    self.settings[f'tmm.destkey.{region}.sfc.an'] = f'tm5-nc:mdir=ec/ea/an0tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'
            self.settings['ndyn'] = '900'
            self.settings['cfl.outputstep'] = '900'

        else :
            for region in self.dconf.run.regions :
                levels = self.dconf.regions[region].levels
                self.settings[f'tmm.sourcekey.{region}.ml'] = f'tm5-nc:mdir=ec/ea/h06h18tr3/{levels}/{region}/<yyyy>/<mm>;tres=_00p03;namesep=/'
                self.settings[f'tmm.sourcekey.{region}.sfc.fc'] = f'tm5-nc:mdir=ec/ea/h06h18tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'
                self.settings[f'tmm.sourcekey.{region}.sfc.an'] = f'tm5-nc:mdir=ec/ea/an0tr1/sfc/{region}/<yyyy>/<mm>;tres=_00p01;namesep=/'

        self.settings['my.levs'] = levels
        self.settings['cfl.outputstep'] = self.settings['ndyn'] 
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
        for region in self.dconf.run.regions :
            if 'redgrid' in self.dconf.regions[region]:
                rgn = self.dconf.regions[region].redgrid.nh
                rgs = self.dconf.regions[region].redgrid.sh
                self.settings[f'region.{region}.redgrid.nh.n'] = len(rgn)
                self.settings[f'region.{region}.redgrid.nh.comb'] = ' '.join([str(_) for _ in rgn])
                self.settings[f'region.{region}.redgrid.sh.n'] = len(rgs)
                self.settings[f'region.{region}.redgrid.sh.comb'] = ' '.join([str(_) for _ in rgn])
            else :
                self.settings[f'region.{region}.redgrid.nh.n'] = 0
                self.settings[f'region.{region}.redgrid.sh.n'] = 0

    def setup_tm5_optim(self):
        """
        This will set the following keys, needed but which should not be needed:
        - correlation.inputdir
        - var4d.horcor.min_eigval
        - var4d.optim_emis.type
        """
#        self.settings['var4d.optim_emis.type'] = '1'
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

    def setup_emissions2(self):
        """
        This will set the following rc-keys:
        - emissions.{tracer}.{region}.ncats
        - emissions.{tracer}.{region}.categories
        - emissions.{tracer}.prefix
        """
        for tr in self.dconf.run.tracers :
            self.settings[f'emissions.{tr}.prefix'] = self.dconf.emissions[tr].prefix
            for region in self.dconf.run.regions:
                catlist = []
                for catname, cat in self.dconf.emissions[tr].categories.items():
                    # By default, the category is in all the regions, unless a "regions" subkey is definbed
                    if region in cat.get('regions', [region]):
                        catlist.append(catname)
                self.settings[f'emissions.{tr}.{region}.ncats'] = len(catlist)
                self.settings[f'emissions.{tr}.{region}.categories'] = ', '.join(catlist)

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

        tracers = self.dconf.run.get('tracers', None)
        if tracers is None:
            self.settings['tracers.number'] = 0
            return
        
        self.settings['tracers.number'] = len(tracers)
        self.settings['tracers.names'] = ','.join([_ for _ in tracers])
        for tr in tracers :
            spec = self.dconf.tracers[tr].species
            tracer = chem.species.get(spec, None)

            # If the species is not present in the "tm5.species" module, then flux_unit and mix_unit MUST be defined
            # otherwise they are just optional
            if tracer is not None:
                emis_unit = self.dconf.tracers[tr].get('flux_unit', tracer.unit_emis)
                mix_unit = self.dconf.tracers[tr].get('mix_unit', tracer.unit_mix)
            else :
                emis_unit = self.dconf.tracers[tr].flux_unit
                mix_unit = self.dconf.tracers[tr].mix_unit

            self.settings[f'tracers.{tr}.molar_mass'] = 1 / ureg.Quantity(f'g{spec}').to('mol').m
            self.settings[f'tracers.{tr}.mixrat_unit_value'] = ureg.Quantity('mol / mol').to(mix_unit).m
            self.settings[f'tracers.{tr}.mixrat_unit_name'] = str(mix_unit)
            self.settings[f'tracers.{tr}.emis_unit_value'] = ureg.Quantity(f'kg{spec}').to(emis_unit).m
            self.settings[f'tracers.{tr}.emis_unit_name'] = str(emis_unit)
            self.settings[f'tracers.{tr}.species'] = spec.lower()
            apply_chem = self.dconf.tracers[tr].get('chemistry', False)
            if apply_chem:
                self.settings[f'tracers.{tr}.chemistry'] = 'T'
                self.settings['proces.chemistry'] = 'T'
                reactions = self.dconf.tracers[tr].reactions
                self.settings[f'tracers.{tr}.nreac'] = len(reactions)
                self.settings[f'tracers.{tr}.reaction_names'] = list(reactions.keys())
                for reacname, reac in reactions.items():
                    self.settings[f'tracers.{tr}.{reacname}.file'] = reac.file
                    self.settings[f'tracers.{tr}.{reacname}.domain'] = reac.domain
                    self.settings[f'tracers.{tr}.{reacname}.rate'] = reac.rate
            else:
                self.settings[f'tracers.{tr}.chemistry'] = 'F'

    def setup_system(self) -> None:
        """
        Setup system-specific settings (paths, etc.). For now the following rc-keys are set:
        - udunits_path
        """
        self.settings['udunits_path'] = Path(self.dconf.run.paths.udunits).absolute()

    # Internal methods

    def setup_output_point(self, dconf):
        self.settings['output.point'] = 'T'
        self.settings['output.point.input.dir'] = dconf.input_dir
        self.settings['output.point.split.period'] = 'a'  # no splitting ...
        self.settings['output.point.sample.parent'] = dconf.get('sample_parent', 'F')

    def setup_output_stations(self, dconf):
        self.settings['output.station.timeseries'] = 'T'
        self.settings['output.station.timeseries.filename'] = Path(dconf.filename).absolute()

    def setup_output_mix(self, dconf: DictConfig | None = None):
        if dconf is None :
            return
        
        self.settings['output.mix'] = 'T'
        self.settings['output.mix.tstep'] = int(Timedelta(dconf.output_frequency).total_seconds())
        self.settings['output.mix.meteo'] = dconf.get('output_meteo', 'F')
        self.settings['output.mix.filename.prefix'] = dconf.prefix
        self.settings['output.mix.deflate.leve'] = dconf.get('deflate_level', 1)