#!/bin/env python

from pyshell.base.main.Emissions import Emissions as baseEmis
from numpy import linspace
from netCDF4 import Dataset
import os

class Emissions(baseEmis):
    def __init__(self, rcf, **kwargs):
        self.StartDate = rcf.ti
        self.EndDate = rcf.tf
        self.subdir_tag = ''
        self.sec_day = 86400.
        self.rcf = rcf
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.GetZoomRegions()
        # Monte Carlo estimation of posterior covariance or not?
        self.optim_mc = self.rcf.get('my.optimizer.class') in ['conGrad_MC', 'm1qn3_MC']
        emission_file_name = self.rcf.get('PyShell.em.filename')
        self.emission_file_name = self.putDateString(emission_file_name, True)
        self.createDifferentialArea()
        # Define a destination horizontal grid for later
        self.lat_grid = {}
        self.lon_grid = {}
        for region_name, lat_data in zip(self.zoom_regions_names, self.zoom_regions_lat):
            ybeg, yend, jm = lat_data
            self.lat_grid[region_name] = linspace(ybeg, yend, jm+1)
        for region_name, lon_data in zip(self.zoom_regions_names, self.zoom_regions_lon):
            xbeg, xend, im = lon_data
            self.lon_grid[region_name] = linspace(xbeg, xend, im+1)

class tm5Emis(Emissions):
    """Use emissions produced by a previous TM5 run
       This is merely a tracer-specific copy of the readOptimFromFile method in the base Emissions class"""

    def __init__(self, rcf, tracer='CO2', *args, **kwargs):
        self.emFile = rcf.get('tmflex.emfile')
        self.tracer = tracer
        Emissions.__init__(self, rcf, *args, **kwargs)
        if os.path.basename(self.emFile) == 'optimized_state.nc4':
            self.emfield = 'poste_emission'
            self.errorfield = 'poste_emission_std'
            if kwargs.get('step', None) == 'apri':
                self.emfield = 'prior_emission'
                self.errorfield = 'prior_emission_std'
        else :
            self.emfield = 'emission'
            self.errorfield = 'emission_error'

    def LoopThroughPeriods(self):
        try :
            fid = Dataset(self.emFile, 'r')
        except RuntimeError :
            print self.emFile
            raise RuntimeError
        self.Emission[self.tracer]['emi_class'] = self.__class__.__name__
        self.Emission[self.tracer]['tf_bb_diurnal'] = None
        for region in self.zoom_regions_names:
            rgid = fid.groups[region]
            tgid = rgid.groups[self.tracer]
            categories = self.Emission[region][self.tracer]['categories']
            for cat in categories:
                cgid = tgid.groups[cat]
                self.Emission[region][self.tracer][cat]['emission_data'] = cgid.variables[self.emfield][:]
                if 'emission_error' in cgid.variables:
                    self.Emission[region][self.tracer][cat]['emission_error'] = cgid.variables[self.errorfield][:]
                else:
                    error = self.Emission[region][self.tracer][cat]['error']
                    self.Emission[region][self.tracer][cat]['emission_error'] = 0.01*error*abs(self.Emission[region][self.tracer][cat]['emission_data'])
        fid.close()


class tmflexEmis(Emissions):
    def __init__(self, rcf, *args, **kwargs):
        Emissions.__init__(self, rcf, *args, **kwargs)
        self.tracer = rcf.get('tmflex.tracer')
        self.parent = rcf.get('tmflex.%s.parent'%self.tracer)
        self.lon0 = rcf.get('tmflex.lon0')
        self.lon1 = rcf.get('tmflex.lon1')
        self.lat0 = rcf.get('tmflex.lat0')
        self.lat1 = rcf.get('tmflex.lat1')

    def LoopThroughPeriods(self):
        from numpy import linspace
        self.Emission[self.tracer]['emi_class'] = self.__class__.__name__
        self.Emission[self.tracer]['tf_bb_diurnal'] = self.Emission[self.parent]['tf_bb_diurnal']
        if len(self.Emission['regions']) > 1:
            raise NotImplementedError('The Rodenbeck scheme with zoom regions should technically work, but the cutout of the emissions at the edge of the domain may not be correct. Comment this line to ignore')
        for region in self.Emission['regions']:
            lons = linspace(self.zoom_regions_lon[0][0], self.zoom_regions_lon[0][1], self.zoom_regions_lon[0][2]+1)
            lats = linspace(self.zoom_regions_lat[0][0], self.zoom_regions_lat[0][1], self.zoom_regions_lat[0][2]+1)
            for cat in self.Emission[region][self.tracer]['categories']:
                # copy all emissions from parent tracer
                em = self.Emission[region][self.parent][cat]['emission_data'] +0.
                # set emissions outside of domain of interest to zero
                em[:,:,lons[:-1]<self.lon0] = 0.  # 0 if western edge outside flexpart domain
                em[:,:,lons[1:]>self.lon1] = 0.   # 0 if eastern edge outside flexpart domain
                em[:,lats[:-1]<self.lat0,:] = 0.  # 0 if southern edge outside flexpart domain
                em[:,lats[1:]>self.lat1,:] = 0.   # 0 if northern edge outside flexpart domain
                # save
                self.Emission[region][self.tracer][cat]['emission_data'] = em

        self.genDailyCycle()

    def genDailyCycle(self):
        """
        Copy daily cycle files from main tracer to global tracer, but zero it outside the tmflex foreground domain
        """
        import glob
        from netCDF4 import Dataset
        from datetime import timedelta
        import shutil
        tt = self.rcf.ti
        while tt < self.rcf.tf :
            prefix1 = os.path.join(self.rcf.get('dailycycle.folder'), tt.strftime('%Y/%m'), self.rcf.get('%s.dailycycle.prefix'%self.parent))
            prefix2 = os.path.join(self.rcf.get('dailycycle.folder'), tt.strftime('%Y/%m'), self.rcf.get('%s.dailycycle.prefix'%self.tracer))
            fname1 = '%s%s.nc4'%(prefix1, tt.strftime('%Y%m%d'))
            fname2 = '%s%s.nc4'%(prefix2, tt.strftime('%Y%m%d'))
            shutil.copy2(fname1, fname2)
            with Dataset(fname2, 'a') as ds:
                for region in self.Emission['regions'] :
                    lons = linspace(self.zoom_regions_lon[0][0], self.zoom_regions_lon[0][1], self.zoom_regions_lon[0][2]+1)
                    lats = linspace(self.zoom_regions_lat[0][0], self.zoom_regions_lat[0][1], self.zoom_regions_lat[0][2]+1)
                    for cat in ds[region].groups :
                        em = ds[region][cat]['emission_anomaly'][:]
                        em[:,:,lons[:-1]<self.lon0] = 0.
                        em[:,:,lons[1:]>self.lon1] = 0.
                        em[:,lats[:-1]<self.lat0,:] = 0.
                        em[:,lats[1:]>self.lat1,:] = 0.
                        ds[region][cat]['emission_anomaly'][:] = em
            tt += timedelta(1)
