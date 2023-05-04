#!/usr/bin/env python

"""
Minimalist python2.7-compatible gridtools library
"""

from numpy import zeros, float64, sin, pi, diff, linspace, moveaxis, ceil, dot
from dataclasses import dataclass
from numpy.typing import NDArray
from types import SimpleNamespace
from typing import List


@dataclass
class RectiLinearGrid:
    lonb : NDArray
    latb : NDArray
    dlon : int
    dlat : int
    radius_earth : float = 6_378_100
    cyclic       : bool = None
    _global      : bool = False

    def __post_init__(self):
        if self.east - self.west == 360 and self.cyclic is None:
            self.cyclic = True
        elif self.cyclic is None :
            self.cyclic = False
        if self.cyclic and self.south == -90 and self.north == 90:
            self._global = True

    @property
    def lonc(self) -> NDArray:
        return (self.lonb[1:] + self.lonb[:-1]) / 2.
    
    @property
    def latc(self) -> NDArray:
        return (self.latb[1:] + self.latb[:-1]) / 2.

    @property
    def lat_borders(self) -> List[slice]:
        return [slice(self.latb[_], self.latb[_ + 1]) for _ in range(self.nlat)]

    @property
    def lon_borders(self) -> List[slice]:
        return [slice(self.lonb[_], self.lonb[_ + 1]) for _ in range(self.nlon)]

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
        lats = (pi / 180.) * self.latb
        for ilat, lat in enumerate(lats):
            area[ilat, :] = self.radius_earth**2 * dlon_rad * sin(lat)
        return diff(area, axis=0)

    def calc_overlap_matrices(self, other : "RectiLinearGrid") -> SimpleNamespace:
        """
        Calculate the overlaping between two regular grids.
        The function returns an namespace with a "lat" and "lon" attribute (arrays):
        - lat[i1, j1] is the fraction of a grid cell of lat index "i1", in the original grid, that is contained by a grid cell of lat index "j1" in the new grid (assuming they have the same longitude boundaries)
        - lon[i2, j2] is the fraction of a longitude interval "i2" (in the original grid) that is contained by a longitude interval "j2" in the new grid (assuming they have the same latitude boundaries).
        The product lat[i1, j1] * lon[i2, j2] gives the fraction of the grid cell (i1, j1) in the original grid that is contained in the grid cell (i2, j2) in the new grid.
        """

        # If both grids are cyclic, create a temporary expanded grid that fully contains the original grid
        # and wrap it around itself in a second time:
        overlaps_lat = self.calc_overlaps_lat(other)

        if self.cyclic and other.cyclic:
            # Create the temporary grid and calculate the overlaps for it
            tmpgrid = other.expand_longitudes(self.west, self.east)
            overlaps_lon = self.calc_overlaps_lon(tmpgrid)

            # Isolate the edge bands (whatever is outside the limits of the original grid)
            west_band = overlaps_lon[:, tmpgrid.lonc < self.west]
            east_band = overlaps_lon[:, tmpgrid.lonc > self.east]

            # Isolate the center band (what is within the limits of the original grid
            overlaps_lon = overlaps_lon[:, (tmpgrid.lonc >= self.west) * (tmpgrid.lonc <= self.east)]

            # Add the edge bands to the center one
            if west_band.shape[1] > 0 :
                overlaps_lon[:, :west_band.shape[1]] += west_band
            if east_band.shape[1] > 0 :
                overlaps_lon[:, -east_band.shape[1]:] += east_band
        else :
            overlaps_lon = self.calc_overlaps_lon(other)

        return SimpleNamespace(lat=overlaps_lat, lon=overlaps_lon)

    def calc_overlaps_lat(self, other: "RectiLinearGrid") -> NDArray:
        overlaps = zeros((self.nlat, other.nlat))

        for ilat1, latb1 in enumerate(self.lat_borders):
            for ilat2, latb2 in enumerate(other.lat_borders):
                # Calculate what fraction of a grid cell ilat1 would end up in a grid cell ilat2
                minlat = max(latb1.start, latb2.start) * pi / 180
                maxlat = min(latb1.stop, latb2.stop) * pi / 180

                if minlat < maxlat:
                    area1 = sin(maxlat - minlat)
                    area2 = sin(latb1.stop * pi / 180 - latb1.start * pi / 180)
                    overlaps[ilat1, ilat2] = area1 / area2
        return overlaps

    def calc_overlaps_lon(self, other: "RectiLinearGrid") -> NDArray:
        overlaps = zeros((self.nlon, other.nlon))
        for ilon1, lonb1 in enumerate(self.lon_borders):
            for ilon2, lonb2 in enumerate(other.lon_borders):
                # Calculate what fraction of a grid cell ilon1 would end up in a grid cell ilon2
                minlon = max(lonb1.start, lonb2.start)
                maxlon = min(lonb1.stop, lonb2.stop)
                if minlon < maxlon :
                    overlaps[ilon1, ilon2] = (maxlon - minlon) / (lonb1.stop - lonb1.start)
        return overlaps

    def expand_longitudes(self, west: float, east: float) -> "RectiLinearGrid":
        """
        Create a new grid, spanning (if necessary) a wider range of longitudes.
        """
        nsteps_east = 0
        nsteps_west = 0
        if west < self.west:
            nsteps_west = ceil((self.west - west) / self.dlon)
            west = self.west - self.dlon * nsteps_west

        if east > self.east:
            nsteps_east = ceil((east - self.east) / self.dlon)
            east = self.east + self.dlon * nsteps_east

        lonb = linspace(west, east, self.nlon + nsteps_east + nsteps_east + 1)
        return RectiLinearGrid(lonb=lonb, latb=self.latb, radius_earth=self.radius_earth, cyclic=False, dlon=self.dlon, dlat=self.dlat)


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
    def global3x2(cls) -> "TM5Grids":
        return cls.from_corners(-180, 180, -90, 90, 3, 2)

    @classmethod
    def global6x4(cls) -> "TM5Grids":
        return cls.from_corners(-180, 180, -90, 90, 6, 4)

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


@dataclass
class SpatialData:
    data        : NDArray
    grid        : RectiLinearGrid
    lon_axis    : int
    lat_axis    : int
    density     : bool = False

    def to_quantity(self, inplace: bool = False) -> "SpatialData":
        """
        Converts amounts (e.g. kg, umol) to (spatial) fluxes (e.g. kg/m2, umol/m2, etc.).
        The temporal dimension is unchanged (e.g. kg/s will be converted to kg/s/m**2)

        Args:
            inplace (bool, optional): whether to modify the object in memory or to return a new one.
        """
        # 1) Move the lat and lon to last positions:
        # data = self.data.swapaxes(self.lon_axis, -1).swapaxes(self.lat_axis, -2)
        data = self._reorder_axes(position='end')
        data = data * self.grid.area
        data = self._reset_axes(data=data, position='end')

        # Return
        if inplace :
            self.data = data
            self.density = False
            return self
        else :
            return SpatialData(data=data, grid=self.grid, lon_axis=self.lon_axis, lat_axis=self.lat_axis)

    def _reorder_axes(self, data : NDArray = None, position : str = 'end') -> NDArray:
        """
        Move the lat and lon axis at the start or end (default) of the axis list
        Arguments :
            position : should be "start" (move lat axis to axis 0 and lon axis to axis 1) or "end" (move lat axis to axis -2 and lon axis to axis -1)
            data: array to reorder (optional, default: self.data)

        Return:
            reordered (view of the) axis
        """
        if data is None :
            data = self.data

        if {'start': True, 'end': False}[position]:         # this will raise and error if position is not "start" or "end"
            if self.lat_axis > self.lon_axis :
                # lat after lon in original array ==> need to move lon first
                data = moveaxis(data, self.lon_axis, 0)   # Move lon axis to first position (no change on lon position)
                data = moveaxis(data, self.lat_axis, 0)   # Move lat axis on first position (changes lon position to axis 1)
            else :
                # lon after lat in original array ==> no change to their relative order
                data = moveaxis(data, self.lat_axis, 0)   # Move lat axis to first position (no change on lon position)
                data = moveaxis(data, self.lon_axis, 1)   # Move lon axis to second position (no change on lat position)
        else :
            if self.lat_axis > self.lon_axis:
                # lat after lon in original array ==> move lat to last, then lon to last
                data = moveaxis(data, self.lat_axis, -1)
                data = moveaxis(data, self.lon_axis, -1)
            else :
                # lon after lat in original array ==> move lon to last, then lat to before last
                data = moveaxis(data, self.lon_axis, -1)
                data = moveaxis(data, self.lat_axis, -2)
        return data

    def _reset_axes(self, data : NDArray = None, position: str = 'end') -> NDArray:
        """
        Undo the effect of "reorder_axes". The exact same arguments should be used.
        """
        if data is None :
            data = self.data

        if {'start': True, 'end': False}[position]:         # this will raise and error if position is not "start" or "end"
            if self.lat_axis > self.lon_axis:
                # lat after lon in original array, but lat first in "data"
                data = moveaxis(data, 0, self.lat_axis)     # move lat from position 0 to its final position. lon is now in position 0
                data = moveaxis(data, 0, self.lon_axis)     # move lon from position 0 to its final position (no effect on lat)
            else :
                # lon after lat in original array and in "data" ==> just move lon, then lat
                data = moveaxis(data, 1, self.lon_axis)     # move lon to its final position
                data = moveaxis(data, 0, self.lat_axis)     # move lat to its final position
        else :
            if self.lat_axis > self.lon_axis:
                # lat after lon in original array, but lat first in "data"
                data = moveaxis(data, -1, self.lon_axis)    # move lon to its final position (now before lat). lat is now last axis
                data = moveaxis(data, -1, self.lat_axis)    # move lat to its final position (no effect on lon)
            else :
                # lon after lat in original array and in "data"
                data = moveaxis(data, -2, self.lat_axis)    # move lat to its final position (no effect on lon)
                data = moveaxis(data, -1, self.lon_axis)    # move lon to its final position (no effect on lat)
        return data

    def to_density(self, inplace: bool = True) -> "SpatialData":
        """
        Converts spatial fluxes (e.g. kg/m**2, umol/m**2, s/m**2) to amounts (kg, umol, s, etc.)
        The non-spatial dimensions (if any) are unchanged (e.g. kg/s/m**2 will be converted to kg/s)

        Args:
            inplace (bool, optional): whether to modify the object in memory or to return a new one.
        """
        # 1) Move the lat and lon to last positions:
        data = self._reorder_axes(data=self.data, position='end')
        data = data / self.grid.area
        data = self._reset_axes(data=data, position='end')

        # Return
        if inplace :
            self.data = data
            self.density = False
            return self
        else :
            return SpatialData(data=data, grid=self.grid, lon_axis=self.lon_axis, lat_axis=self.lat_axis)

    def regrid(self, destgrid: RectiLinearGrid) -> "SpatialData":

        # Transition matrices:
        trans = self.grid.calc_overlap_matrices(destgrid)

        # Convert to units / grid-cell, if needed:
        if self.density :
            data = self.to_quantity(inplace=False).data
        else :
            data = self.data

        # Ensure that the lat axis is in last position
        data_out = data.swapaxes(self.lat_axis, -1).copy()
        data_out = dot(data_out, trans.lat)
        data_out = data_out.swapaxes(-1, self.lat_axis)

        # Ensure that the lon axis is in last position:
        data_out = data_out.swapaxes(self.lon_axis, -1)
        data_out = dot(data_out, trans.lon)
        data_out = data_out.swapaxes(-1, self.lon_axis)

        # collapse the extra dimension that may have been created:
        if len(self.data.shape) < len(data_out.shape):
            data_out = data_out.squeeze(axis=-1)

        # Create the output structure:
        data_out = SpatialData(data=data_out, grid=destgrid, lon_axis=self.lon_axis, lat_axis=self.lat_axis)

        # Re-convert to density, if needed:
        if self.density :
            data_out.to_density()

        return data_out
