#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True

import re, os
from pyshell.tmflex import rc
from numpy import *
from pyshell.base.helper.Utilities import checkDir
from datetime import datetime
from netCDF4 import Dataset
from collections import defaultdict
from pyshell.base.main.meteo import LevelDefinitions
from pyshell.base.main.Emissions import del_time
from calendar import isleap

class InitialConcentration(LevelDefinitions):

    def __init__(self, IniDate):
        LevelDefinitions.__init__(self) # Get the pressure level boundaries, and initialize self.rcf
        self.GetZoomRegions()
        self.dy_min = min((self.region_lats[:,1] - self.region_lats[:,0])/self.region_lats[:,2])
        self.dx_min = min((self.region_lons[:,1] - self.region_lons[:,0])/self.region_lons[:,2])
        self.max_lat_div = int(180./self.dy_min)
        self.max_lon_div = int(360./self.dx_min)
        self.species = [s.strip() for s in self.rcf.get('my.tracer').split(',')]
        self.ntracer = len(self.species)
        self.inidate = datetime(*IniDate)
        # Define a destination horizontal grid for later
        self.lat_grid = {}
        self.lon_grid = {}
        for region_name, lat_data in zip(self.region_names, self.region_lats):
            ybeg, yend, jm = lat_data
            self.lat_grid[region_name] = linspace(ybeg, yend, jm+1)
        for region_name, lon_data in zip(self.region_names, self.region_lons):
            xbeg, xend, im = lon_data
            self.lon_grid[region_name] = linspace(xbeg, xend, im+1)

    def create_iniconc_structure(self):
        self.IniConc = dict.fromkeys(self.region_names)
        for region in self.region_names:
            self.IniConc[region] = dict.fromkeys(self.species)
            for tracer in self.species:
                self.IniConc[region][tracer] = {}
        for tracer in self.species:
            self.IniConc[tracer] = {}
            self.IniConc[tracer]['iniconc_class'] = self.rcf.get(tracer+'.iniconc.class')

    def AssembleIniconc(self):
        pass

    def AddToIniconc(self, dConc):
        pass

    def GetZoomRegions(self):
        """
        Reads a config file to get zoom region info
        """
        self.region_names = self.rcf.get('regions').split()
        xbeg = []
        xend = []
        im = []
        ybeg = []
        yend = []
        jm = []
        self.zoom_info = defaultdict(dict)
        for region in self.region_names:
            xbeg.append(self.rcf.get('region.'+region+'.xbeg','float'))
            self.zoom_info[region]['xbeg'] = self.rcf.get('region.'+region+'.xbeg','float')
            ybeg.append(self.rcf.get('region.'+region+'.ybeg','float'))
            self.zoom_info[region]['ybeg'] = self.rcf.get('region.'+region+'.ybeg','float')
            xend.append(self.rcf.get('region.'+region+'.xend','float'))
            self.zoom_info[region]['xend'] = self.rcf.get('region.'+region+'.xend','float')
            yend.append(self.rcf.get('region.'+region+'.yend','float'))
            self.zoom_info[region]['yend'] = self.rcf.get('region.'+region+'.yend','float')
            im.append(self.rcf.get('region.'+region+'.im','int'))
            self.zoom_info[region]['im'] = self.rcf.get('region.'+region+'.im','int')
            jm.append(self.rcf.get('region.'+region+'.jm','int'))
            self.zoom_info[region]['jm'] = self.rcf.get('region.'+region+'.jm','int')
        self.region_lons = array(zip(xbeg,xend,im))
        self.region_lats = array(zip(ybeg,yend,jm))
        self.region_parents = None # Fix this!
        self.nlev = self.rcf.get('my.nlay', 'int')

    def lumpField(self, in_mix, in_isobars, l_x, l_y, surface_area):
        """
        in_mix has concentrations on a fine grid, in_isobars the isobars on the same fine grid, and l_x and l_y are the elements
        along lon and lat to clump together. For example, going from a 0.5x0.5 to a 3x2 grid, they would be 6 and 4 respectively.
        surface_area is a 2D array of shape in_mix.shape[1:], specifying the surface area of each grid point.
        """
        nlev = in_mix.shape[0]
        len_x = in_mix.shape[2]/l_x
        len_y = in_mix.shape[1]/l_y
        dummyVar = zeros((nlev,len_y, len_x), in_mix.dtype)
        for i in range(len_x):
            for j in range(len_y):
                for i_lev in range(nlev):
                    weights = surface_area[j*l_y:(j+1)*l_y,i*l_x:(i+1)*l_x] * (in_isobars[i_lev,j*l_y:(j+1)*l_y,i*l_x:(i+1)*l_x] - in_isobars[i_lev+1,j*l_y:(j+1)*l_y,i*l_x:(i+1)*l_x])
                    dummyVar[i_lev,j,i] = average(in_mix[i_lev,j*l_y:(j+1)*l_y,i*l_x:(i+1)*l_x], weights=weights)
        return dummyVar

    def spreadOverFineGrid(self):
        # re-grid self.isobars and self.conc_3d over the finest grid
        print 'Spreading concentrations over finest grid'
        nlev = self.conc_3d.shape[0]
        fine_xCO2 = zeros((nlev, self.max_lat_div, self.max_lon_div), self.conc_3d.dtype)
        fine_isob = zeros((nlev+1, self.max_lat_div, self.max_lon_div), self.isobars.dtype)
        fine_lats = linspace(-90.,90.,self.max_lat_div+1)
        fine_lons = linspace(-180.,180.,self.max_lon_div+1)
        fine_lats = 0.5*(fine_lats[1:] + fine_lats[:-1])
        fine_lons = 0.5*(fine_lons[1:] + fine_lons[:-1])
        for i, lon in enumerate(fine_lons):
            i_idx = argmin(abs(lon-self.longitudes))
            for j, lat in enumerate(fine_lats):
                j_idx = argmin(abs(lat-self.latitudes))
                fine_xCO2[:,j,i] = self.conc_3d[:,j_idx,i_idx]
                fine_isob[:,j,i] = self.isobars[:,j_idx,i_idx]
        self.conc_3d = fine_xCO2
        self.isobars = fine_isob
        self.latitudes = fine_lats
        self.longitudes = fine_lons
        del fine_xCO2, fine_isob, fine_lats, fine_lons

    def coarsenField(self):
        EarthRad = 6.371e6 # meters
        lon_div_fine = len(self.longitudes)
        lat_div_fine = len(self.latitudes)
        dLon = (pi/180.) * (360./lon_div_fine)
        dS = zeros((lat_div_fine+1, lon_div_fine), float64)
        Lats = (pi/180.) * linspace(-90., 90., lat_div_fine+1)
        for i, lat in enumerate(Lats):
            dS[i] = EarthRad * EarthRad * dLon * sin(lat)
        dS = diff(dS, axis=0)
        # if the initial field is glb3x2, dS is now a 90x120 array containing surface areas of the grid cells
        dx_fine = average(diff(self.longitudes))
        dy_fine = average(diff(self.latitudes))
        # if self.conc_3d is 34x90x120, then dx_fine=3 and dy_fine=2
        for xlims,ylims,region_name in zip(self.region_lons,self.region_lats,self.region_names):
            print 'Coarsening concentrations to region ' + region_name
            xs_id = int((xlims[0]+180.0)/dx_fine)
            ys_id = int((ylims[0]+90.0)/dy_fine)
            xe_id = int(xs_id + ((xlims[1]-xlims[0])/dx_fine))
            ye_id = int(ys_id + ((ylims[1]-ylims[0])/dy_fine))
            # for glb6x4, xs_id = 0, xe_id = 120, ys_id = 0, ye_id = 90
            relevant_mix = self.conc_3d[:,ys_id:ye_id, xs_id:xe_id] # nlev x lat x lon
            relevant_isob = self.isobars[:,ys_id:ye_id, xs_id:xe_id] # nlev+1 x lat x lon
            relevant_area = dS[ys_id:ye_id, xs_id:xe_id] # lat x lon
            coarsen_x = int(((xlims[1]-xlims[0])/xlims[2])/dx_fine)
            coarsen_y = int(((ylims[1]-ylims[0])/ylims[2])/dy_fine)
            self.IniConc[region_name][self.tracer]['mixing ratio'] = self.lumpField(relevant_mix, relevant_isob, coarsen_x, coarsen_y, relevant_area)

    def WriteIniConc(self):
        if self.IniConc['type'] == 'save':
            self.iniconc_file_name = self.rcf.get("start.3.filename")
        elif self.IniConc['type'] == 'mix':
            self.iniconc_file_name = self.rcf.get("start.2.iniconcfile")
        else:
            raise ValueError('Unknown initial concentration type: %s'%self.IniConc['type'])
        outFile = self.inidate.strftime(self.iniconc_file_name)
        if os.path.exists(outFile):
            os.remove(outFile)
        checkDir(outFile)
        fid = Dataset(outFile, mode='w', format='NETCDF4')
        fid.createDimension('ntracet', self.ntracer)
        for region in self.region_names:
            if self.IniConc['type'] == 'save':
                lm,jm,im = self.IniConc[region]['m'].shape
            elif self.IniConc['type'] == 'mix':
                lm,jm,im = self.IniConc[region]['shape']
            else:
                raise ValueError('Unknown initial concentration type: %s'%self.IniConc['type'])
            gid = fid.createGroup(region)
            gid.createDimension('lm', lm)
            gid.createDimension('jm', jm)
            gid.createDimension('im', im)
            gid.createDimension('lm_p1', lm+1)
            if self.IniConc['type'] == 'save':
                v = gid.createVariable('m', float64, ('lm','jm','im'))
                v[:] = self.IniConc[region]['m']
                for var_name in ['rm', 'rxm', 'rym', 'rzm']:
                    write_arr = zeros((self.ntracer, lm, jm, im), float64)
                    for i_tracer, tracer in enumerate(self.species):
                        write_arr[i_tracer] = self.IniConc[region][tracer][var_name]
                    v = gid.createVariable(var_name, float64, ('ntracet','lm','jm','im'))
                    v[:] = write_arr
            elif self.IniConc['type'] == 'mix':
                for tracer in self.species:
                    tgid = gid.createGroup(tracer)
                    v = tgid.createVariable('mixing_ratio', float64, ('lm','jm','im'))
                    v[:] = self.IniConc[region][tracer]['mixing_ratio']
                    if 'pressure_boundaries' in self.IniConc[region][tracer]:
                        v = tgid.createVariable('pressure_boundaries', float64, ('lm_p1', 'jm', 'im'))
                        v[:] = self.IniConc[region][tracer]['pressure_boundaries']
                    if 'aux fields' in self.IniConc[region][tracer]:
                        for var_name in self.IniConc[region][tracer]['aux fields']:
                            v = tgid.createVariable(var_name, self.IniConc[region][tracer][var_name].dtype, ('lm','jm','im'))
                            v[:] = self.IniConc[region][tracer][var_name]
            else:
                raise ValueError('Unknown initial concentration type: %s'%self.IniConc['type'])
        fid.close()

    def get_class_from_name( self, class_name):
        _temp = __import__('IniConc', fromlist=[class_name])
        try:
            class_from_name = _temp.__dict__[class_name]
            return class_from_name
        except KeyError:
            sys.stderr.write("Class %s not defined in %s.\n"%(class_name,'IniConc'))
            sys.exit()
        except:
            sys.stderr.write("Unknown error importing %s\n"%class_name)
            sys.exit()

