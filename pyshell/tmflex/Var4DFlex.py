#!/bin/env python
import sys, os, shutil
from pyshell.guillaume.Optimizer import conGrad
from pyshell.guillaume.RunTM5 import RunTM5

class tm5Run:

    def setupOutputs(self, stationlist):
        self.rcf.setkey( 'output.station.timeseries.filename', stationlist)
        self.rcf.setkey('input.station.filename', stationlist)
        #fid = open('stationList.%s.dat'%self.rcf.get('my.project'), 'w')
        #self.rcf.setkey('output.station.timeseries.filename', os.path.join(os.getcwd(), fid.name))
        #self.rcf.setkey('input.station.filename', os.path.join(os.getcwd(), fid.name))
        #fid.write(' ID     LAT     LON     ALT TP STATIONNAME\n')
        #for site in output_locations:
        #    fid.write('%3s %7.2f %7.2f %7.1f FM %s\n'%(site.code, site.lat, site.lon, site.alt, site.name))
        #fid.close()

    def ContinueVar4D(self, setupObs=False):
        run = RunTM5(self.rcf)
        if setupObs: run.SetupObservations(self.obs)
        try :
            opt = conGrad(run, resume=True)
            opt.Var4D()
        except :
            print('resume=True did not work, going through restart=True')
            opt = conGrad(run)
            opt.SetupOptimizer(restart=True, emclasses=self.emis)
            opt.load_convergence()
            opt.Var4D()
	
        if self.archive:
            self.archiveFile.add(run.rcf.get('rundir'))
            self.archiveFile.save()
            opt.load_convergence()
        return opt

    def appendProj(self, projpath):
        ## add the correct "project" keys:
        for key in  ['my.projects.basic', 'my.source.dirs', 'build.copy.dirs']:
            value =  self.rcf.get(key)
            value += ' '+projpath
            self.rcf.setkey(key, value)

    def ComputeBackgrounds(self, reg, step=None, compile=True, setupObs=False):
        self.appendProj('%s/tmflex'%self.rcf.get('my.proj.root'))
        tracer = self.rcf.get('my.tracer')
        regions = self.rcf.get('regions').split(',')
        ncat = int(self.rcf.get('%s.categories'%tracer))
        # check that we are in a 1-tracer configuration (it should work for more with some easy code fixes)
        if len(tracer.strip().split(',')) > 1: raise RuntimeError('Rodenbeck scheme is available only on the basis of one tracer')
        if not tracer in self.emis.keys(): raise RuntimeError('The emissions provided do not match with the tracer requested')
        if not tracer+'fg' in self.emis.keys(): raise RuntimeError('The emissions provided do not match with the tracer requested')
        # add a new tracer to the TM5 run
        self.rcf.setkey('my.tracer', '%s, %sfg'%(tracer, tracer))
        self.rcf.setkey('tracers', '2')
        # Build keys for that second tracer
        # This is probably rather CO2-specific, but since I don't have any other tracer to test it with, I'll go with it
        keysToChange = [x%tracer for x in ['%s.emission.dailycycle', '%s.dailycycle.type', 'output.point.%s.minerror']]
        for region in regions:
            keysToChange.append('emission.%s.%s.categories'%(tracer, region))
            for icat in xrange(1, ncat+1):
                keysToChange.append('emission.%s.%s.category%i'%(tracer, region, icat))
        for key in keysToChange:
            self.rcf.setkey(key.replace(tracer,'%sfg'%tracer), self.rcf.get(key))
        # add tmflex-specific keys to the rc-file
        self.rcf.setkey('tmflex.compute.backgrounds', True)
        self.rcf.setkey('tmflex.lon0', reg['lonmin'])#-2.)
        self.rcf.setkey('tmflex.lon1', reg['lonmax'])#32.)
        self.rcf.setkey('tmflex.lat0', reg['latmin'])#30.)
        self.rcf.setkey('tmflex.lat1', reg['latmax'])#65.)

        # Daily cycle for foreground tracer:
        self.rcf.setkey('%sfg.dailycycle.prefix'%tracer, self.rcf.get('%s.dailycycle.prefix'%tracer).replace(tracer, '%sfg'%tracer))

        self.RunForward(step=step, compile=compile, setupObs=setupObs)

    def saveEmissions(self, path):
        outdir = self.rcf.get('output.dir')
        emfile = 'optimized_state.nc4'
        path = os.path.join(path, os.path.basename(outdir.strip('/')))
        if not os.path.isdir(path): os.makedirs(path)
        shutil.copyfile(os.path.join(outdir, emfile), os.path.join(path, emfile))
        return os.path.join(path, emfile)

    def saveObservations(self, path):
        import glob
        outdir = self.rcf.get('output.dir')
        obsdirs = [os.path.join(outdir, 'point'), os.path.join(outdir, 'point_apri'), os.path.join(outdir, 'stations'), self.rcf.get('output.point.input.dir')]
        path = os.path.join(path, os.path.basename(outdir.strip('/')))
        if not os.path.isdir(path): os.makedirs(path)
        for dir in obsdirs:
            dirname = os.path.basename(dir.strip('/'))
            if not os.path.isdir(os.path.join(path, dirname)): os.makedirs(os.path.join(path, dirname))
            flist = glob.glob('%s/*.nc4'%dir)
            for file in flist:
                shutil.copy2(file, os.path.join(path, dirname))
        self.savedStationsPath = os.path.join(path, 'stations')
