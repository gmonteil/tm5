#!/usr/bin/env python

"""
Minimalist python2.7-compatible gridtools library
"""

from numpy import zeros, float64, sin, pi, diff, linspace, arange
from dataclasses import dataclass
from numpy.typing import NDArray


@dataclass
class RectiLinearGrid:
    lonb : NDArray
    latb : NDArray
    dlon : int
    dlat : int
    radius_earth : float = 6_378_100

    @property
    def lonc(self) -> NDArray:
        return (self.lonb[1:] + self.lonb[:-1]) / 2.
    
    @property
    def latc(self) -> NDArray:
        return (self.latb[1:] + self.latb[:-1]) / 2.
        
    @property
    def west(self) -> float:
        return self.lonb.min()
    
    @property
    def east(self) -> float:
        return self.lonb.max()
    
    @property
    def north(self) -> float:
        return self.latb.max()
    
    @property
    def south(self) -> float:
        return self.latb.min()
    
    @property
    def nlon(self) -> int:
        return len(self.lonc)
    
    @property
    def nlat(self) -> int:
        return len(self.latc)
    
    @property
    def area(self) -> NDArray:
        dlon_rad = self.dlon * pi / 180.
        area = zeros((self.nlat+1, self.nlon), float64)
        lats = ( pi / 180. ) * self.latb
        for ilat, lat in enumerate(lats):
            area[ilat, :] = self.radius_earth**2 * dlon_rad * sin(lat)
        return diff(area, axis=0)
    

@dataclass
class TM5Grids(RectiLinearGrid):
    @classmethod
    def global1x1(cls) -> "TM5Grids":
        return cls(
            lonb = linspace(-180, 180, 361),
            latb = linspace(-90, 90, 181),
            dlon = 1,
            dlat = 1
        )

    @classmethod
    def from_corners(cls, west: float, east: float, south: float, north: float, dlon: int, dlat: int) -> "TM5Grids":
        nlon = (east - west) / dlon
        nlat = (north - south) / dlat
        assert (nlon - int(nlon) == 0) & (nlat - int(nlat) == 0)
        return cls(
            lonb = linspace(west, east, int(nlon) + 1),
            latb = linspace(south, north, int(nlat) + 1),
            dlon = dlon, dlat = dlat
        )