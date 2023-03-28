#!/usr/bin/env python

"""
Minimalist python2.7-compatible gridtools library
"""

from numpy import zeros, float64, sin, pi, diff, linspace, arange


class RectiLinearGrid:
    def __init__(self, lonb, latb, dlon, dlat):
        self.lonb = lonb
        self.latb = latb
        self.dlon = dlon
        self.dlat = dlat
        self.radius_earth = 6378100

    @property
    def lonc(self):
        return (self.lonb[1:] + self.lonb[:-1]) / 2.
    
    @property
    def latc(self):
        return (self.latb[1:] + self.latb[:-1]) / 2.
        
    @property
    def west(self):
        return self.lonb.min()
    
    @property
    def east(self):
        return self.lonb.max()
    
    @property
    def north(self):
        return self.latb.max()
    
    @property
    def south(self):
        return self.latb.min()
    
    @property
    def nlon(self):
        return len(self.lonc)
    
    @property
    def nlat(self):
        return len(self.latc)
    
    @property
    def area(self):
        dlon_rad = self.dlon * pi / 180.
        area = zeros((self.nlat+1, self.nlon), float64)
        lats = ( pi / 180. ) * self.latb
        for ilat, lat in enumerate(lats):
            area[ilat, :] = self.radius_earth**2 * dlon_rad * sin(lat)
        return diff(area, axis=0)
    
    
class TM5Grids(RectiLinearGrid):
    @classmethod
    def global1x1(cls):
        return cls(
            lonb = linspace(-180, 180, 361),
            latb = linspace(-90, 90, 181),
            dlon = 1.,
            dlat = 1.
        )

    @classmethod
    def from_corners(cls, latb, lonb):
        dlon = lonb[1] - lonb[0]
        dlat = latb[1] - latb[0]
        return cls(lonb = lonb, latb=latb, dlon=dlon, dlat=dlat)