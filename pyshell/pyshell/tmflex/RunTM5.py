#!/bin/env pythonf
from pyshell.proj.tracer.CO2.RunTM5 import RunTM5 as runtm5
from pyshell.base.helper.Utilities import checkDir
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pyshell.tmflex.observations import Observations
from pyshell.tmflex.PointObs import PointObs
import glob, subprocess, shutil
import pyshell
import os
import sys


class RunTM5(runtm5):
    def __init__(self, rcf, dconf):
        self.rcf = rcf
        self.dconf = dconf
        # The key 'my.tracer' is a comma-separated list of tracer, whereas the key 'my.tracer.name' is a name for the
        # combination. For example, my.tracer could be 'CO,CO2' and my.tracer.name could be 'COCO2'.
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.runid = self.rcf.get('runid')
        self.StartTime = rcf.ti
        self.EndTime = rcf.tf
        self.subdir_tag = ''
        # SBi : Check whether all the following keys are necessary
        self.rcf.add('timerange.start', self.StartTime.strftime("%Y-%m-%d %H:%M:%S"))
        self.rcf.add('timerange.end', self.EndTime.strftime("%Y-%m-%d %H:%M:%S"))
        # SBf
        self.J_bg = 0.0
        self.J_obs = 0.0
        self.J_tot = self.J_bg + self.J_obs
        self.Optim_dict = {'emission': 'optimize.emission', 'iniconc': 'optimize.initialconcentration',
                           'parameters': 'optimize.parameter'}
        for k, v in self.Optim_dict.items():
            self.Optim_dict[k] = self.rcf.get(v, 'bool')
        self.months = []
        d = self.StartTime.replace(day=1).replace(hour=0).replace(minute=10)
        while d <= self.EndTime:
            self.months.append(d.strftime('%Y%m'))
            d += relativedelta(months=1)
        self.region_names = self.rcf.get('regions').split()
        self.nregions = len(self.region_names)
        self.output_dir = self.rcf.get('output.dir')
        self.emission_file_name = self.rcf.get('PyShell.em.filename')
        # delete all the tm5-pyshell.*.rc files in the output directory
        files_to_delete = glob.glob(os.path.join(self.rcf.get('my.run.dir'), self.rcf.get('my.basename') + '.*.rc'))
        for fileName in files_to_delete:
            try:
                os.remove(fileName)
            except:
                print 'File ', fileName, ' has been already deleted!'
        self.GetZoomRegions()

    def SetupEmissions(self, emclasses, randomize=False, zero=False, step=None):
        #TODO ==> This should not be here, at all!!!
        from pyshell.tmflex.emissions.emissions import Emissions
        """
        This will either assemble the emissions (default), or read them from a file, depending on
        an rc key, emission.read.optimized.
        """

        t1 = datetime.now()

        emisFile = self.emission_file_name
        if os.path.exists(emisFile):
            os.remove(emisFile)
        checkDir(emisFile)

        # First, process the emissions and write the emissions file
        StartTuple = self.StartTime.timetuple()[:5]
        EndTuple = self.EndTime.timetuple()[:5]
        em = Emissions(self.rcf)
        em.create_emission_structure()

        for tracer in em.species:
            emfill = emclasses[tracer](self.rcf, self.dconf.emissions[tracer], step=step, tracer=tracer)
            emfill.Emission = em.Emission
            emfill.LoopThroughPeriods()

        self.Emission = em(randomize, zero)  # Write emissions and keep them in memory

        for region in self.region_names:
            for tracer in self.species:
                categories = self.Emission[region][tracer]['categories']
                for cat in categories:
                    print(region, tracer, cat, self.Emission[region][tracer][cat]['time_resolution'], ' Optimize?:', self.Emission[region][tracer][cat]['optimize'], self.Emission[region][tracer][cat]['emission_data'].shape)
                for cat_key in self.Emission[tracer]['cat_list']:
                    re = self.Emission[tracer]['cat_opt'][cat_key]
                    print('========================================')
                    print('Tracer: ', tracer, 'Category marked for optimization:', cat_key)
                    for rei in re:
                        print('----->region:', rei['region'], 'error:', rei['error'])
                    print('========================================')
        print('========================================')
        print('end info parsed from emission routine')
        print('========================================')

        t2 = datetime.now()
        sys.stderr.write("Emission between %s and %s assembled in %s\n" % (
        self.StartTime.strftime("%c"), self.EndTime.strftime("%c"), str(t2 - t1)))

    def Compile(self, new_build=False):
        """
        Compile TM5 via pycasso scripting: build, make
        """
        # rcfile for TM5 following pycasso requirements:
        if not os.path.isdir(self.rcf.get('my.run.dir')):
            os.makedirs(self.rcf.get('my.run.dir'))
        rcfile = os.path.join(self.rcf.get('my.run.dir'), 'compile.rc')
        self.rcf.WriteFile(rcfile)
        # command to setup and submit a run:
        if new_build:
            command = [os.path.join(pyshell.__path__[0], '../bin/setup_tm5'), '--new', rcfile]
        else:
            command = [os.path.join(pyshell.__path__[0], '../bin/setup_tm5'), rcfile]
        # run, check status:
        subprocess.check_call(command)

    def PointDepartures(self):
        dp = PointObs(self)
        for month_tuple in dp.months:
            dp.create_pointdeparture_structure()
            for tracer in dp.species:
                dp.applyPointObs(tracer, month_tuple)
            dp.writePointDepartureFile(month_tuple)

    def SetupObservations(self, obsobj):
        ti = self.StartTime
        dt = self.EndTime - ti
        while ti < self.EndTime:
            obs = Observations(ti, ti + dt, self.rcf)
            for tracer in obs.species:
                obs.SetupPointObs(tracer, obsobj)
            obs.writePointFile()
            ti = ti + dt

    def ObservationCost(self):
        J_obs = 0.0
        if self.rcf.get('adjoint.input.satellite', 'bool'):
            satobs = ApplySatelliteObs(self)
            sys.stderr.write('Satellite departures calculated\n')
            J_obs += satobs.SatelliteCost()
            sys.stderr.write('Satellite cost function calculated\n')
        if self.rcf.get('adjoint.input.point', 'bool'):
            pointobs = PointObs(self)
            sys.stderr.write('Point departures calculated\n')
            J_obs += pointobs.PointCost()
            sys.stderr.write('Point cost function calculated\n')
        return J_obs

    def GetZoomFiles(self):
        regions = self.rcf.get('regions').split()
        if len(regions) > 1:
            raise NotImplementedError(
                "I have no idea on what files should be copied if a zoom is defined, therefore it will crash")
        else:
            region = regions[0]
            zoompath = os.path.join(self.rcf.get('machine.input'), 'zoom')
            zoomfile = '%s/Zoomed_%s.nc4' % (zoompath, region)
            vecfile = '%s/vec2ll_%s.nc' % (zoompath, region)
            runpath = self.rcf.get('rundir')
            if not os.path.isfile(zoomfile):
                print("File %s not found. Generate it first with a forward run" % zoomfile)
            if not os.path.isfile(vecfile):
                print("File %s not found. Generate it first with a forward run" % vecfile)
            if not os.path.isfile(zoomfile) or not os.path.isfile(vecfile): raise IOError
            if not os.path.isdir(runpath): os.makedirs(runpath)
            shutil.copy(zoomfile, '%s/Zoomed.nc4' % runpath)
            shutil.copy(vecfile, runpath)
