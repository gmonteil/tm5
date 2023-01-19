#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
sys.dont_write_bytecode = True

from numpy import *
from netCDF4 import Dataset
from datetime import datetime
from pyshell.tm5_utils import interpolate_fields, redistrib_flux

from pyshell.base.main.IniConc import InitialConcentration

class CT2013B(InitialConcentration):
    """
    This class creates the initial concentration field on a date from CT2013B posterior concentrations. To
    avoid issues with different air masses, we will create 3D mixing ratio fields.
    """
    def __init__(self, IniDate):
        InitialConcentration.__init__(self, IniDate)
        max_date = self.rcf.get('iniconc.CO2.maximum.date')
        max_date = datetime.strptime(max_date, '%Y-%m-%d %H:%M')
        if self.inidate > max_date:
            raise RuntimeError('Asked to create field for %s, beyond the maximum allowed %s'%\
                (self.inidate.strftime('%c'), max_date.strftime('%c')))
        self.input_dir = self.rcf.get('iniconc.CO2.input.folder')
        self.filename_pre = self.inidate.strftime(self.rcf.get('iniconc.CO2.filename.pre'))
        self.filename_post = self.inidate.strftime(self.rcf.get('iniconc.CO2.filename.post'))
        self.tracer = 'CO2'
        self.ref_region = self.rcf.get('iniconc.CO2.ref.region')
        # Calculate the relevant time index within the day
        self.hour_index = self.inidate.hour/3 # depends on integer division, won't work in Py3k (idiots!)

    def AssembleIniconc(self):
        """
        We need to fill the structure self.IniConc. Specifically, the following elements for all regions:
        self.IniConc[region][self.tracer]['mixing_ratio'] = mix_ratio
        self.IniConc[region][self.tracer]['pressure_boundaries'] = ip_pres_levels
        self.IniConc[region]['shape'] = mix_ratio.shape
        and one element globally:
        self.IniConc['type'] = 'mix'
        """
        for ireg, region in enumerate(self.region_names):
            shortname = self.rcf.get('my.region%1is'%(ireg+1))
            file_name = self.filename_pre + shortname + self.filename_post
            file_name = os.path.join(self.input_dir, file_name)
            if os.path.exists(file_name):
                # mixing ratio for this region already exists, just read it in
                with Dataset(file_name, 'r') as fid:
                    self.IniConc[region][self.tracer]['mixing_ratio'] = fid.variables['co2'][self.hour_index]
                    self.IniConc[region][self.tracer]['pressure_boundaries'] = fid.variables['pressure'][self.hour_index]
                    self.IniConc[region]['shape'] = fid.variables['co2'][self.hour_index].shape
            else:
                # construct the mixing ratio for this region from glb3x2 fields
                file_name = self.filename_pre + self.ref_region + self.filename_post
                file_name = os.path.join(self.input_dir, file_name)
                with Dataset(file_name, 'r') as fid:
                    ip_pres_levels = fid.variables['pressure'][self.hour_index]
                    orig_mix_ratio = fid.variables['co2'][self.hour_index]

                    nlat = len(fid.dimensions['lat'])
                    nlon = len(fid.dimensions['lon'])
                    nlev = len(fid.dimensions['level'])

                # we assume reference region is global at some sort of resolution
                orig_lats = linspace(-90., 90., nlat+1)
                orig_lons = linspace(-180., 180., nlon+1)

                new_lats = self.lat_grid[region]
                new_lons = self.lon_grid[region]
                # we can use redistrib_flux.regrid_fluxes, because just like fluxes defined per unit area, both the surface
                # pressure and the mixing ratio need to be distributed with the surface area of the grid cells as weights
                sp = redistrib_flux.regrid_fluxes(orig_lats, orig_lons, new_lats, new_lons, ip_pres_levels[0], True)
                op_pres_levels = zeros((nlev+1, sp.shape[0], sp.shape[1]), float64)
                for l in range(nlev+1):
                    op_pres_levels[l] = self.AT[l] + self.BT[l] * sp

                mix_ratio = zeros((nlev, self.zoom_info[region]['jm'], self.zoom_info[region]['im']), float64)
                for ilev in range(nlev):
                    mix_ratio[ilev] = redistrib_flux.regrid_fluxes(orig_lats, orig_lons, new_lats, new_lons, orig_mix_ratio[ilev], True)

                self.IniConc[region][self.tracer]['mixing_ratio'] = mix_ratio
                self.IniConc[region][self.tracer]['pressure_boundaries'] = op_pres_levels
                self.IniConc[region]['shape'] = mix_ratio.shape

        self.IniConc['type'] = 'mix'

class Miller_CO2(InitialConcentration):
    """
    This class creates the initial concentration field from CT2013 prior fields as produced by John Miller
    for the tracer CO2. John's files are all on 34 levels, while I want to use 25 levels, and I don't
    know how to regrid slopes. Therefore, I will convert John's masses into mixing ratios and then rebin.
    """
    def __init__(self, IniDate):
        InitialConcentration.__init__(self, IniDate)
        # The filename pattern can be, e.g., ${tm5.data.input.dir}/starting_fields/Miller/TM5_restart_%Y%m%d_%H%M_
        # It will be understood that we need to add something like '_glb3x2.nc' or '_nam1x1.nc' after it.
        self.iniconc_filePattern = self.rcf.get('iniconc.CO2.miller.input.filename')
        self.tracer_bg_name = self.rcf.get('iniconc.CO2.miller.CO2.bg.name')
        self.tracer_names = self.rcf.get('iniconc.CO2.miller.CO2.tracers').split()
        self.tracer = 'CO2'
        self.fscale = 657578.57

    def readFile(self, pick_time):
        file_prefix = pick_time.strftime(self.iniconc_filePattern)
        for i_reg,region in enumerate(self.region_names):
            reg_shortname = self.rcf.get('my.region%1is'%(i_reg+1))
            file_postfix = '%s.nc'%reg_shortname
            file_name = file_prefix + file_postfix
            # John's files exist for glb3x2 and nam1x1, so if we want a different region, e.g., nam3x2,
            # glb6x4 or eur1x1, we need to construct it from glb3x2
            if os.path.exists(file_name):
                with Dataset(file_name, 'r') as fid:
                    tracer_names = fid.variables['names'][:]
                    tracer_names = [''.join(t).strip() for t in tracer_names]
                    tracer_bg_index = tracer_names.index(self.tracer_bg_name)
                    bg = fid.variables['rm'][tracer_bg_index]
                    tracer_mass = zeros_like(bg)
                    tracer_mass += bg
                    print "Added index %i (tracer %s) as background"%(tracer_bg_index, tracer_names[tracer_bg_index])
                    for tracer in self.tracer_names:
                        idx = tracer_names.index(tracer)
                        tracer_mass += (fid.variables['rm'][idx] - bg)
                        print "Added index %i (tracer %s) to background"%(idx, tracer_names[idx])
                    print
                    # read the airmass, which is not tracer-specific
                    airmass = fid.variables['m'][:]
                    # calculate the mixing ratio
                    mix_ratio = self.fscale * tracer_mass / airmass
                    # calculate the input pressure levels
                    at = fid.variables['at'][:]
                    bt = fid.variables['bt'][:] # first index is the surface
                    sp = fid.variables['sp'][:]
            else:
                # construct the mixing ratio from the glb3x2 field
                basename = os.path.basename(file_name)
                all_parts = basename.split('_')
                all_parts[-1] = 'glb3x2.nc'
                basename = '_'.join(all_parts)
                file_name = os.path.join(os.path.dirname(file_name), basename)
                with Dataset(file_name, 'r') as fid:
                    tracer_names = fid.variables['names'][:]
                    tracer_names = [''.join(t).strip() for t in tracer_names]
                    tracer_bg_index = tracer_names.index(self.tracer_bg_name)
                    bg = fid.variables['rm'][tracer_bg_index]
                    tracer_mass = zeros_like(bg)
                    tracer_mass += bg
                    for tracer in self.tracer_names:
                        idx = tracer_names.index(tracer)
                        tracer_mass += (fid.variables['rm'][idx] - bg)
                    # read the airmass, which is not tracer-specific
                    airmass = fid.variables['m'][:]
                    # calculate the mixing ratio
                    orig_mix_ratio = self.fscale * tracer_mass / airmass
                    # calculate the input pressure levels
                    at = fid.variables['at'][:]
                    bt = fid.variables['bt'][:] # first index is the surface
                    orig_sp = fid.variables['sp'][:]
                # now mix_ratio and sp are at global 3x2 resolution, select the relevation lateral section
                orig_lats = linspace(-90., 90., 91)
                orig_lons = linspace(-180., 180., 121)
                # we can use redistrib_flux.regrid_fluxes, because just like fluxes defined per unit area,
                # both the surface pressure and the mixing ratio need to be distributed with the surface area
                # of the grid cells as weights
                new_lats = self.lat_grid[region]
                new_lons = self.lon_grid[region]
                sp = redistrib_flux.regrid_fluxes(orig_lats, orig_lons, new_lats, new_lons, orig_sp, True)
                mix_ratio = zeros((orig_mix_ratio.shape[0], self.zoom_info[region]['jm'], self.zoom_info[region]['im']), float64)
                for ilev in range(mix_ratio.shape[0]):
                    mix_ratio[ilev] = redistrib_flux.regrid_fluxes(orig_lats, orig_lons, new_lats, new_lons, orig_mix_ratio[ilev], True)

            # now rebin vertically, if needed
            if len(at) != len(self.AT):
                rebinv = True
            else:
                rebinv = (abs(at - self.AT).max() > 0.1) and (abs(bt - self.BT).max() > 1.0e-6)

            if rebinv:
                print "Different vertical grid, rebinning needed"

            ip_pres_levels = zeros((at.shape[0], sp.shape[0], sp.shape[1]), float64)
            for l in range(at.shape[0]):
                ip_pres_levels[l] = at[l] + bt[l] * sp
            # now calculate the output pressure levels
            op_pres_levels = zeros((self.nlev+1, sp.shape[0], sp.shape[1]), float64)
            for l in range(self.nlev+1):
                op_pres_levels[l] = self.AT[l] + self.BT[l] * sp
            orig_surface = average(mix_ratio, weights=-diff(ip_pres_levels, axis=0), axis=0)
            # now rebin the mixing ratio onto the new grid
            if rebinv:
                mix_ratio = interpolate_fields.rebin_3d( ip_pres_levels[::-1], mix_ratio[::-1], \
                    op_pres_levels[1:-1][::-1] )[::-1]
            # check the error made in rebinning
            new_surface = average(mix_ratio, weights=-diff(op_pres_levels, axis=0), axis=0)
            print '%s CO2 : Max rebinning difference = %12.5e'%(region, abs(orig_surface-new_surface).max())
            self.IniConc[region][self.tracer]['mixing_ratio'] = mix_ratio
            if rebinv:
                self.IniConc[region][self.tracer]['pressure_boundaries'] = op_pres_levels
            else:
                self.IniConc[region][self.tracer]['pressure_boundaries'] = ip_pres_levels
            self.IniConc[region]['shape'] = mix_ratio.shape
        self.IniConc['type'] = 'mix'

    def AssembleIniconc(self):
        self.readFile(self.inidate)
