#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Miscallaneous plot utility functions."""

#-- system packages
import os
import sys
import datetime as dtm
import calendar
import time
import matplotlib as mpl
import cartopy.crs as ccrs
import numpy as np


def cartopy_get_projection( proj_str : str ) -> ccrs:
    """return cartopy projection instance for current identifier.
    """
    if proj_str=='PlateCarree':
        proj = ccrs.PlateCarree()
    elif proj_str=='Mollweide':
        proj = ccrs.Mollweide()
    elif proj_str=='Robinson':
        proj = ccrs.Robinson()
    else:
        f"requested projection -->{proj_str}<-- currently not supported."
        raise RuntimeError(msg)

    return proj
        

def cartopy_default_geoax(fig, options,
                          cartopy_proj=ccrs.PlateCarree(), axes_extent=[0.1,0.1,0.8,0.8],
                          add_gridlines=True, **kwargs):
    """
    """
    import cartopy.feature as cfeature
    import cartopy.io.img_tiles as cimgt
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    #-- some settings
    borders_lw = kwargs.get('borders_lw', 0.3)
    coast_lw   = kwargs.get('coast_lw', 0.3)
    grid_lw    = kwargs.get('grid_lw', 0.5)
    grid_color = kwargs.get('grid_color', 'gray')
    xlabsiz    = kwargs.get('xlabsiz',6)
    ylabsiz    = kwargs.get('ylabsiz',6)

    #
    #-- axes instance
    #
    geo_ax = fig.add_axes(axes_extent, projection=cartopy_proj)
    

    # geo_ax = fig.add_axes(axes_extent,
    #                       projection=ccrs.EuroPP())
    # geo_ax = fig.add_axes(axes_extent,
    #                       projection=ccrs.AzimuthalEquidistant(central_latitude=30, central_longitude=0))
    #-- background map
    if options.background_map=='OSM':
        request = cimgt.OSM()
        geo_ax.add_image(request, options.background_zoom)
    elif options.background_map=='Stamen':
        bg_image = cimgt.Stamen('terrain-background')
        geo_ax.add_image(bg_image, options.background_zoom)
    #-- add coast lines (always)    
    geo_ax.coastlines(linewidth=coast_lw)
    #-- add country boarders
    geo_ax.add_feature(cfeature.BORDERS, lw=borders_lw)
    #geo_ax.set_global()
    #-- grid lines
    if add_gridlines:
        gl = geo_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              linewidth=grid_lw, color=grid_color)
        gl.top_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': xlabsiz, 'color':'gray'}
        gl.ylabel_style = {'size': ylabsiz, 'color':'gray'}

    return geo_ax


def cartopy_geoax_add_grid(geoax, lonreso, latreso, nwse=None, **kwargs):
    """
    """
    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker

    n,w,s,e = nwse if nwse!=None else [90, -180, -90, 180]
    crs = ccrs.PlateCarree()

    gl = geoax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray')
    gl.top_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    xlabsiz = kwargs.get('xlabsiz',6)
    ylabsiz = kwargs.get('ylabsiz',6)
    gl.xlabel_style = {'size': xlabsiz, 'color':'gray'}
    gl.ylabel_style = {'size': ylabsiz, 'color':'gray'}
    
    xlocs = np.hstack((np.arange(w, e, lonreso),e))
    ylocs = np.hstack((np.arange(s, n, latreso),n))
    gl.xlocator = mticker.FixedLocator(xlocs)
    gl.ylocator = mticker.FixedLocator(ylocs)


def cartopy_geoax_add_platecarre_gridlines(geoax, longrid, latgrid, **kwargs):
    """
    """
    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker

    crs = ccrs.PlateCarree()

    gl = geoax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    xlabsiz = kwargs.get('xlabsiz',8)
    ylabsiz = kwargs.get('ylabsiz',8)
    xlabcol = kwargs.get('xlabcol','black')
    ylabcol = kwargs.get('ylabcol','black')
    gl.xlabel_style = {'size': xlabsiz, 'color':xlabcol}
    gl.ylabel_style = {'size': ylabsiz, 'color':ylabcol}
    
    gl.xlocator = matplotlib.ticker.FixedLocator(longrid)
    gl.ylocator = matplotlib.ticker.FixedLocator(latgrid)


def geomap_suggest_scale(w,e,s,n, scalemax=19):
    """Function to set an empirical zoom scale/level suitable for given extent.

    Parameter
    ---------
    w,e,s,n : float
        western/eastern/south/northern boarder of region of interest [degree]
    Returns
    -------
    
    """
    # empirical solve for scale based on zoom
#    scale = np.ceil(-np.sqrt(2)*np.log(np.divide((extent[1]-extent[0])/2.0,350.0)))
    scale = np.ceil(-np.sqrt(2)*np.log(np.divide((e-w)/2.0,350.0)))
    # scale cannot be larger than 19
    scale = (scale<20) and scale or scalemax

    return scale


def axes_set_time_default_grid(ax):
    """For x/y plot with time on x-axis prepare default
    ticks and grid lines

    Parameters
    ----------
    ax : matplotlib.axes
        axes instance for which the x-ticks and grid are being set.
    """
    import matplotlib as mpl
    #-- add grid lines
    ax.grid(which='minor', axis='x', color='grey', alpha=0.25)
    ax.grid(which='major', axis='both', color='grey', alpha=0.5, lw=1)
    #-- major grid-lines begin/mid of year
    majLoc = mpl.dates.MonthLocator(bymonth=[1])
    majFmt  = mpl.dates.DateFormatter('%Y')
    #-- minor grid-lines at first of Apr/Oct (without labels)
    minorLoc  = mpl.dates.MonthLocator(bymonth=[4,7,10])
    ax.xaxis.set_major_locator(majLoc)
    ax.xaxis.set_major_formatter(majFmt)
    ax.xaxis.set_minor_locator(minorLoc)


#%%%%%%%%%%%%%%%%%%%%
#https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
# with modifications
class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        # MVO-NOTE::we try to handle masked values appropriately
        xp, fp = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        xy_intp = np.empty(value.shape, dtype=value.dtype)
        vmask   = np.ma.getmaskarray(value)
        xy_intp = np.ma.masked_where(vmask, xy_intp)
        ##MV-NOTE::2019-02-19:
        #          discovered that np.interp returns fp[0] for value<self.vmin and
        #          fp[-1] for value>self.vmax.
        #          For the colorbar this had caused some problems with colorbar ticklabels.
        #          As workaround/solution we do now call the 'super' method from (i.e. from
        #          mpl.colors.Normalize) for those particular values
        #
        msk_lo = np.logical_and(~vmask,value<self.vmin)
        msk_hi = np.logical_and(~vmask,value>self.vmax)
        msk_mi = np.logical_and(~vmask,(value>=self.vmin)&(value<=self.vmax))
        xy_intp[msk_mi] = np.interp(value[msk_mi], xp, fp)
        xy_intp[msk_lo] = super(MidpointNormalize,self).__call__(value[msk_lo], clip=clip)
        xy_intp[msk_hi] = super(MidpointNormalize,self).__call__(value[msk_hi], clip=clip)
        # print "MVMV::values.size={} values.nmask={}".format(
        #     value.size, np.ma.count_masked(value))
        # print "MVMV::xy_intp.nmasked={}".format(np.ma.count_masked(xy_intp))
        # print "MVMV::xy_intp={}".format(xy_intp)
        return xy_intp
# ---MidpointNormalize---


#
#-- source: https://stackoverflow.com/questions/65260995/matplotlib-symmetric-logarithmic-colormap-not-centered-at-zero
#
class MidpointLogNorm(mpl.colors.SymLogNorm):
  def __init__(self, lin_thres, lin_scale, midpoint=None, vmin=None, vmax=None):
      self.midpoint = midpoint
      self.lin_thres = lin_thres
      self.lin_scale = lin_scale
      #fraction of the cmap that the linear component occupies
      self.linear_proportion = (lin_scale / (lin_scale + 1)) * 0.5
      print(self.linear_proportion)

      # Create norm with vmin at 0 and midpoint at 0.5
      self.SymLogNorm1 = mpl.colors.SymLogNorm(lin_thres, lin_scale, vmin, 2*self.midpoint + np.abs(vmin))
      # Create norm with midpoint at 0.5 and vmax at 1
      self.SymLogNorm2 = mpl.colors.SymLogNorm(lin_thres, lin_scale, 2*self.midpoint - vmax, vmax)
      mpl.colors.SymLogNorm.__init__(self, lin_thres, lin_scale, vmin, vmax)

  def __get_value__(self, v, log_val1_i, log_val2_i, clip=None):
      v = np.array(v)
      x = [self.vmin, self.midpoint, self.vmax]
      y = [0., 0.5, 1.]
      interpol = np.interp(v, x, y)
      out = np.where(np.abs(v) < self.lin_thres, interpol, v)
      out = np.where(out > self.lin_thres, log_val2_i, out)
      out = np.where(out < self.lin_thres, log_val1_i, out)
      return np.ma.masked_array(out)

  def __call__(self, value, clip=None):
      log_val1 = self.SymLogNorm1(value)
      log_val2 = self.SymLogNorm2(value)

      out = [0] * len(value)
      for i, v in enumerate(value):
          out[i] = self.__get_value__(v, log_val1[i], log_val2[i])
      return np.ma.masked_array(out)


def cbar_create(mappable, ax, cmap, cnorm, **kwargs):
    """
    pure helper routine to create color bar with some properties
    """
    import matplotlib.pyplot as plt

    cborientation = kwargs['cborientation'] if 'cborientation' in kwargs else 'horizontal'
    cbextend = kwargs['cbextend'] if 'cbextend' in kwargs else 'neither'
    cbextendrect = kwargs['cbextendrect'] if 'cbextendrect' in kwargs else False
    cbextendfrac = kwargs['cbextendfrac'] if 'cbextendfrac' in kwargs else None
    cbticks = kwargs['cbticks'] if 'cbticks' in kwargs else None
    pad = 0.1
    cbar = plt.colorbar( mappable, ax=ax, cmap=cmap, norm=cnorm,
                         orientation=cborientation, pad=pad,
                         ticks=cbticks,
                         extend=cbextend, extendrect=cbextendrect, extendfrac=cbextendfrac )
    if 'cbticklabels' in kwargs:
        cbar.set_ticklabels(kwargs['cbticklabels'])
    if 'cblabsize' in kwargs:
        cbar.ax.tick_params(labelsize=kwargs['cblabsize'])

    return cbar


def cbar_index( ncolors, cmap,
                orientation='horizontal', extend='neither', extendrect=False,
                **kwargs ):
    """

    ---
    Ref: http://stackoverflow.com/questions/18704353/correcting-matplotlib-colorbar-ticks
    """
    cmap = cmap_discretize(cmap, ncolors)
    mappable = mpl.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar( mappable,
                             orientation=orientation,
                             extend=extend, extendrect=extendrect )#, **kwargs)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(list(range(ncolors)))

    return colorbar


def discrete_colmap(cmap, bounds):
    """Function to create a "discrete" color map for the given boundary values using
    the input colormap as base."""

    resCmap = {}
    ncols = len(bounds[1:]) #one less than #bounds
    resCmap['cmap'] = mpl.colors.ListedColormap(
        [cmap(i/float(ncols)) for i in range(ncols)]
    )
    resCmap['bounds'] = bounds
    resCmap['norm'] = mpl.colors.BoundaryNorm(bounds, ncols)

    return resCmap


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. mpl.cm.jet. 
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(mpl.cm.jet, 5)
        imshow(x, cmap=djet)

    ---
    Ref: http://stackoverflow.com/questions/18704353/correcting-matplotlib-colorbar-ticks
    """

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in range(N+1) ]
    # Return colormap object.
    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


def cmap_truncate( cmap, minval=0.0, maxval=1.0, n=100):
    """Return new colormap derived from an existing colormap,
    by clipping the range to minval,maxval, and optionally limiting
    the number of colors used
    """
    if type(cmap)==str:
        cmap = plt.get_cmap(cmap)
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def cmap_reverse(cmap, name = 'my_cmap_r'):
    """
    In: 
    cmap, name 
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]

    ---
    Ref: http://stackoverflow.com/questions/3279560/invert-colormap-in-matplotlib
    """
    return mpl.colors.LinearSegmentedColormap(name, mpl.cm.revcmap(cmap._segmentdata)) 
    # reverse = []
    # k = []   

    # for key in cmap._segmentdata:    
    #     k.append(key)
    #     channel = cmap._segmentdata[key]
    #     print "MVVM::",channel
    #     data = []

    #     for t in channel:                    
    #         data.append((1-t[0],t[2],t[1]))            
    #     reverse.append(sorted(data))    

    # LinearL = dict(zip(k,reverse))
    # my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    # return my_cmap_r


def cnorm_set(pltkw, vvmin, vvmax):
    """
    adjust colorbar settings according to the plot keywords
    and data range
    """
    #-- copy input (plot) dictionary
    pkw = dict(pltkw)

    #-- user-specified data limits
    vmin = pltkw.get('vmin', vvmin)
    vmax = pltkw.get('vmax', vvmax)
    #-- request for specific min/max/centre at colorbar
    cbmin    = pltkw.get('cbmin', None)
    cbmax    = pltkw.get('cbmax', None)
    cbcentre = pltkw.get('cbctr', None)
    lognorm  = pltkw.get('lognorm', False)
    linthresh = pltkw.get('linthresh', 1e-10) #-- linear threshold, required for SymLogNorm
    linscale  = pltkw.get('linscale', 1.)

    #-- cNorm==None lets matplotlib decide on it's own rules
    cNorm = None

    #
    #-- color norm selection
    #
    def _cnorm_select( vvmin : float, vvmax : float, lognorm : bool, cbctr : float = None ):
        if cbctr!=None:
            if lognorm:
                cNorm = MidpointLogNorm(linthresh, linscale, midpoint=cbctr, vmin=vvmin, vmax=vvmax)
            else:
                cNorm = MidpointNormalize(midpoint=cbctr, vmin=vvmin, vmax=vvmax)
        elif lognorm:
            if vvmin*vvmax<0:
                cNorm = mpl.colors.SymLogNorm(linthresh, vmin=vvmin, vmax=vvmax)
            else:
                cNorm = mpl.colors.LogNorm(vmin=vvmin, vmax=vvmax)
        else:
            cNorm = mpl.colors.Normalize(vmin=vvmin, vmax=vvmax)
        return cNorm

    if cbmin!=None and cbmax!=None:
        cNorm = _cnorm_select(cbmin, cbmax, lognorm, cbctr=cbcentre)
        if vmin < cbmin and vmax > cbmax:
            pkw['cbextend'] = 'both'
        elif vmin < cbmin:
            pkw['cbextend'] = 'min'
        elif vmax > cbmax:
            pkw['cbextend'] = 'max'
    elif cbmin!=None:
        cNorm = _cnorm_select(cbmin, vmax, lognorm, cbctr=cbcentre)
        if vmin < cbmin:
            pkw['cbextend'] = 'min'
    elif cbmax!=None:
        cNorm = _cnorm_select(vmin, cbmax, lognorm, cbctr=cbcentre)
        if vmax > cbmax:
            pkw['cbextend'] = 'max'
    else:
        cNorm = _cnorm_select(vmin, vmax, lognorm, cbctr=cbcentre)

    # print(f"MVDEBUG:: -->{cNorm}<--")

    return (pkw,cNorm)
