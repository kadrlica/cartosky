#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import os
from os.path import expandvars
import logging
from collections import OrderedDict as odict
import warnings

import matplotlib
from matplotlib import mlab
from matplotlib.collections import LineCollection
import matplotlib.patches as mpatches
import pylab as plt
import numpy as np
import ephem
import healpy as hp
import scipy.ndimage as nd

import cartopy.crs as ccrs
from shapely.geometry.polygon import Polygon, LineString
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from cartosky.utils import setdefaults,get_datadir
from cartosky.utils import cel2gal, gal2cel
from cartosky import healpix
from cartosky.constants import COLORS
import cartosky.proj

class Skymap(object):
    """ Base class for creating Skymap objects. """

    defaults = {'global':True, 'gridlines':True, 'lon_0': 0, 'celestial': True}

    def __init__(self, projection='cyl', **kwargs):
        setdefaults(kwargs,self.defaults)
        do_global    = kwargs.pop('global',True)
        do_grid      = kwargs.pop('gridlines',True)
        do_celestial = kwargs.pop('celestial',True)

        self.set_observer(kwargs.pop('observer',None))
        self.set_date(kwargs.pop('date',None))

        # Eventually want to subclass GeoAxes
        ax = plt.gca()
        fig = ax.figure
        subargs = ax.get_geometry()
        fig.delaxes(ax)

        self.projection = cartosky.proj.Proj(projection, **kwargs)
        self.ax = fig.add_subplot(*subargs,projection=self.projection)

        if do_global:
            self.ax.set_global()
        if do_grid:
            self.grid = self.ax.gridlines()
            self.grid.rotate_labels = False
        if do_celestial:
            self.ax.invert_xaxis()

        # Better grid lines?
        #https://github.com/SciTools/cartopy/pull/1117

        self.wrap_angle = np.mod(kwargs['lon_0'] + 180,360)

    def __call__(self, lon, lat, inverse=False):
        lon,lat = np.asarray(lon),np.asarray(lat)
        if inverse:
            proj_xyz = ccrs.PlateCarree().transform_points(self.projection,lon,lat)
        else:
            proj_xyz = self.projection.transform_points(ccrs.PlateCarree(),lon,lat)
        return proj_xyz[...,0],proj_xyz[...,1]

    def proj(self,lon,lat):
        """ Remove points outside of projection """
        return self(lon,lat)

    def set_observer(self, observer):
        observer = observer.copy() if observer else ephem.Observer()
        self.observer = observer

    def set_date(self,date):
        date = ephem.Date(date) if date else ephem.now()
        self.observer.date = date

    # Wrap mpl.axes plotting functions through cartosky.skyaxes.SkyAxes

    def plot(self, *args, **kwargs):
        self.ax.plot(*args, **kwargs)

    def scatter(self, *args, **kwargs):
        self.ax.scatter(*args, **kwargs)

    def pcolormesh(self, *args, **kwargs):
        self.ax.pcolormesh(*args, **kwargs)

    #... need to add all the other functions...

    def get_map_range(self, hpxmap, pixel=None, nside=None):
        """ Calculate the longitude and latitude range for an implicit map. """
        return healpix.get_map_range(hpxmap,pixel,nside,wrap_angle=self.wrap_angle)

    def hpx2xy(self, hpxmap, pixel=None, nside=None, xsize=800,
               lonra=None, latra=None):
        """ Convert from healpix map to longitude and latitude coordinates """
        return healpix.hpx2xy(hpxmap,pixel=pixel,nside=nside,
                              xsize=xsize,aspect=self.aspect,
                              lonra=lonra,latra=latra)

    def smooth(self,hpxmap,badval=hp.UNSEEN,sigma=None):
        """ Smooth a healpix map """
        return healpix.smooth(hpxmap,badval,sigma)

    def draw_hpxmap(self, hpxmap, pixel=None, nside=None, xsize=800,
                    lonra=None, latra=None, badval=hp.UNSEEN, smooth=None,
                    **kwargs):
        """
        Use pcolor/pcolormesh to draw healpix map.

        Parameters
        ----------
        hpxmap: input healpix or HealSparse map
        pixel:  explicit pixel indices in RING scheme (required for partial healpix maps)
        nside:  explicit nside of the map (required for partial healpix maps) if
                passed while visualizing a HealSparse map it will doegrade the map to this nside.
        xsize:  resolution of the output image
        lonra:  longitude range [-180,180] (deg)
        latra:  latitude range [-90,90] (deg)
        badval: set of values considered "bad"
        smooth: gaussian smoothing kernel (deg)
        kwargs: passed to pcolormesh

        Returns
        -------
        im,lon,lat,values : mpl image with pixel longitude, latitude (deg), and values
        """
        return self.ax.hpxmap(hpxmap, pixel, nside, xsize,
                    lonra, latra, badval, smooth, **kwargs)

    def draw_hpxbin(self, lon, lat, nside=256, **kwargs):
        """
        Create a healpix histogram of the counts.

        Like `hexbin` from matplotlib

        Parameters:
        -----------
        lon : input longitude (deg)
        lat : input latitude (deg)
        nside : heaplix nside resolution
        kwargs : passed to draw_hpxmap and plt.pcolormesh

        Returns:
        --------
        hpxmap, im : healpix map and image
        """
        return self.ax.hpxbin(lon,lat,nside,**kwargs)

    def draw_line_radec(self,ra,dec,**kwargs):
        """Draw a line assuming a Geodetic transform.

        Parameters
        ----------
        ra    : right ascension (deg)
        dec   : declination (deg)
        kwargs: passed to plot

        Returns
        -------
        feat  : cartopy.FeatureArtist
        """
        # Color will fill a polygon...
        # https://github.com/SciTools/cartopy/issues/856
        color=kwargs.pop('c',kwargs.pop('color','k'))
        defaults=dict(crs=ccrs.Geodetic(),edgecolor=color,facecolor='none')
        setdefaults(kwargs,defaults)

        line = LineString(list(zip(ra,dec))[::-1])
        return self.ax.add_geometries([line],**kwargs)

    def draw_polygon_radec(self,ra,dec,**kwargs):
        """Draw a shapely Polygon from a list of ra,dec coordinates.

        Parameters
        ----------
        ra    : right
        dec   : declination
        kwargs: passed to add_geometries

        Returns
        -------
        poly  : the Polygon
        """
        defaults=dict(crs=ccrs.Geodetic(), facecolor='none', edgecolor='red')
        setdefaults(kwargs,defaults)
        ra  = np.asarray(ra).flatten()
        dec = np.asarray(dec).flatten()
        coords = np.vstack([ra,dec]).T
        poly = Polygon(coords)
        self.ax.add_geometries([poly], **kwargs)
        if 'label' in kwargs:
            self.ax.plot(np.nan,np.nan,color=kwargs['edgecolor'],label=kwargs['label'])
        return poly

    def draw_polygon(self,filename,reverse=True,**kwargs):
        """Draw a text file containing ra,dec coordinates of polygon(s)

        Parameters
        ----------
        filename: name of the file containing the polygon(s) [ra,dec,poly]
        kwargs:   keyword arguments passed to

        Returns
        -------
        poly:     polygons
        """
        try:
            data = np.genfromtxt(filename,names=['ra','dec','poly'])
        except ValueError:
            from numpy.lib.recfunctions import append_fields
            data = np.genfromtxt(filename,names=['ra','dec'])
            data = append_fields(data,'poly',np.zeros(len(data)))

        ret = []
        for p in np.unique(data['poly']):
            poly = data[data['poly'] == p]
            ra  = poly['ra'][::-1] if reverse else poly['ra']
            dec = poly['dec'][::-1] if reverse else poly['dec']
            feat = self.draw_polygon_radec(ra,dec,**kwargs)
            ret += [feat]
            kwargs.pop('label',None)

        return ret

    # Alias for draw
    draw_polygons = draw_polygon

    def tissot(self, *args, **kwargs):
        self.ax.tissot(*args,**kwargs)

    def tissot_indicatrices(self, *args, **kwargs):
        self.ax.tissot_indicatrices(*args,**kwargs)

    def draw_zenith(self, radius=1.0, **kwargs):
        """
        Plot a to-scale representation of the zenith.

        Parameters
        ----------
        radius : radius of zenith circle (deg)
        kwargs : passed to plotting routines

        Returns
        -------
        None
        """
        defaults = dict(color='green',alpha=0.75,lw=1.5,)
        setdefaults(kwargs,defaults)

        # RA and Dec of zenith
        zra, zdec = np.degrees(self.observer.radec_of(0, '90'))
        self.plot(zra,zdec,marker='+',ms=10,mew=1.5, **kwargs)
        if radius:
            kwargs['edgecolor']=kwargs.pop('color')
            kwargs['facecolor']='none'
            self.tissot(zra,zdec,radius,**kwargs)

    def draw_airmass(self, airmass=1.4, **kwargs):
        """
        Draw circle around zenith with given airmass.

        Parameters
        ----------
        airmass : airmass (secz) of circle to draw
        kwargs  : passed to draw_zenith

        Returns
        -------
        None
        """
        altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)
        self.draw_zenith(radius=np.degrees(altitude_radians))

    def draw_moon(self, date):
        moon = ephem.Moon()
        moon.compute(date)
        ra_moon = np.degrees(moon.ra)
        dec_moon = np.degrees(moon.dec)

        x,y = self.proj(np.array([ra_moon]), np.array([dec_moon]))
        if np.isnan(x).all() or np.isnan(y).all(): return

        self.scatter(x,y,color='%.2f'%(0.01*moon.phase),edgecolor='black',s=600)
        color = 'black' if moon.phase > 50. else 'white'
        #text = '%.2f'%(0.01 * moon.phase)
        text = '%2.0f%%'%(moon.phase)
        plt.text(x, y, text, fontsize=10, ha='center', va='center', color=color)

    def draw_milky_way(self,width=10,**kwargs):
        """ Draw the Milky Way galaxy. """
        defaults = dict(lw=1.5,ls='-')
        setdefaults(kwargs,defaults)

        glon = np.linspace(0,360,500)
        glat = np.zeros_like(glon)
        ra,dec = gal2cel(glon,glat)

        line = self.draw_line_radec(ra,dec,**kwargs)
        ret = [line]

        if width:
            kwargs.update(ls='--',lw=1)
            for delta in [+width,-width]:
                ra,dec = gal2cel(glon,glat+delta)
                line = self.draw_line_radec(ra,dec,**kwargs)
                ret += [line]
        return ret

    def draw_lmc(self,**kwargs):
        from cartosky.constants import RA_LMC, DEC_LMC, RADIUS_LMC
        defaults = dict(fc='0.7',ec='0.5')
        setdefaults(kwargs,defaults)
        proj = self.proj(RA_LMC,DEC_LMC)
        self.tissot(RA_LMC,DEC_LMC,RADIUS_LMC,**kwargs)
        plt.text(proj[0],proj[1], 'LMC', weight='bold',
                 fontsize=10, ha='center', va='center', color='k')

    def draw_smc(self,**kwargs):
        from cartosky.constants import RA_SMC, DEC_SMC, RADIUS_SMC
        defaults = dict(fc='0.7',ec='0.5')
        setdefaults(kwargs,defaults)
        proj = self.proj(RA_SMC,DEC_SMC)
        self.tissot(RA_SMC,DEC_SMC,RADIUS_SMC,**kwargs)
        plt.text(proj[0],proj[1], 'SMC', weight='bold',
                 fontsize=8, ha='center', va='center', color='k')


    def set_scale(self, array, log=False, sigma=1.0, norm=None):
        if isinstance(array,np.ma.MaskedArray):
            out = np.ma.copy(array)
        else:
            out = np.ma.array(array,mask=np.isnan(array),fill_value=np.nan)

        if sigma > 0:
            out.data[:] = nd.gaussian_filter(out.filled(0),sigma=sigma)[:]

        if norm is None:
            norm = np.percentile(out.compressed(),97.5)

        if log:
            out = np.log10(out)
            if norm: norm = np.log10(norm)

        out /= norm
        out = np.clip(out,0.0,1.0)
        return out

    def draw_inset_colorbar(self,format=None,label=None,ticks=None,fontsize=11,**kwargs):
        defaults = dict(width="25%", height="5%", loc=7,
                        bbox_to_anchor=(0.,-0.04,1,1))
        setdefaults(kwargs,defaults)

        ax = plt.gca()
        im = plt.gci()
        cax = inset_axes(ax,bbox_transform=ax.transAxes,**kwargs)
        cmin,cmax = im.get_clim()

        if (ticks is None) and (cmin is not None) and (cmax is not None):
            cmed = (cmax+cmin)/2.
            delta = (cmax-cmin)/10.
            ticks = np.array([cmin+delta,cmed,cmax-delta])

        tmin = np.min(np.abs(ticks[0]))
        tmax = np.max(np.abs(ticks[1]))

        if format is None:
            if (tmin < 1e-2) or (tmax > 1e3):
                format = '$%.1e$'
            elif (tmin > 0.1) and (tmax < 100):
                format = '$%.1f$'
            elif (tmax > 100):
                format = '$%i$'
            else:
                format = '$%.2g$'
                #format = '%.2f'

        kwargs = dict(format=format,ticks=ticks,orientation='horizontal')

        if format == 'custom':
            ticks = np.array([cmin,0.85*cmax])
            kwargs.update(format='$%.0e$',ticks=ticks)

        cbar = plt.colorbar(cax=cax,**kwargs)
        cax.xaxis.set_ticks_position('top')
        cax.tick_params(axis='x', labelsize=fontsize)

        if format == 'custom':
            ticklabels = cax.get_xticklabels()
            for i,l in enumerate(ticklabels):
                val,exp = ticklabels[i].get_text().split('e')
                ticklabels[i].set_text(r'$%s \times 10^{%i}$'%(val,int(exp)))
            cax.set_xticklabels(ticklabels)

        if label is not None:
            cbar.set_label(label,size=fontsize)
            cax.xaxis.set_label_position('top')

        plt.sca(ax)
        return cbar,cax

    def zoom_to_fit(self, hpxmap, pixel=None, nside=None):
        lonra, latra = self.get_map_range(hpxmap, pixel, nside)
        self.zoom_to(lonra,latra)

    def zoom_to(self, lonra, latra):
        """ Zoom the map to a specific longitude and latitude range.

        Parameters:
        -----------
        lonra : Longitude range [lonmin,lonmax]
        latra : Latitude range [latmin,latmax]

        Returns:
        --------
        None
        """

        (lonmin,lonmax), (latmin,latmax) = lonra, latra

        ax = plt.gca()
        self.llcrnrx,self.llcrnry = self(lonmin,latmin)
        self.urcrnrx,self.urcrnry = self(lonmax,latmax)

        ax.set_xlim(self.llcrnrx,self.urcrnrx)
        ax.set_ylim(self.llcrnry,self.urcrnry)

        #self.set_axes_limits(ax=ax)

    def draw_fields(self,fields,**kwargs):
        # Scatter point size is figsize dependent...
        defaults = dict(edgecolor='none',s=15)
        # case insensitive without changing input array
        names = dict([(n.lower(),n) for n in fields.dtype.names])

        if self.projection == 'ortho': defaults.update(s=50)
        if 'filter' in names:
            colors = [COLORS[b] for b in fields[names['filter']]]
            defaults.update(c=colors)
        elif 'band' in names:
            colors = [COLORS[b] for b in fields[names['band']]]
            defaults.update(c=colors)

        setdefaults(kwargs,defaults)
        ra,dec = fields[names['ra']],fields[names['dec']]
        self.scatter(ra,dec,**kwargs)

    def draw_focal_planes(self, ra, dec, **kwargs):
        from cartosky.instrument.decam import DECamFocalPlane
        defaults = dict(alpha=0.2,color='red',edgecolors='none',lw=0,
                        transform=ccrs.PlateCarree())
        setdefaults(kwargs,defaults)
        ra,dec = np.atleast_1d(ra,dec)
        if len(ra) != len(dec):
            msg = "Dimensions of 'ra' and 'dec' do not match"
            raise ValueError(msg)
        decam = DECamFocalPlane()
        # Should make sure axis exists....
        ax = plt.gca()
        for _ra,_dec in zip(ra,dec):
            corners = decam.rotate(_ra,_dec)
            collection = matplotlib.collections.PolyCollection(corners,**kwargs)
            ax.add_collection(collection)
        plt.draw()

    draw_decam = draw_focal_planes

class McBrydeSkymap(Skymap):
    defaults = dict(Skymap.defaults)
    def __init__(self,**kwargs):
        setdefaults(kwargs,self.defaults)
        super(McBrydeSkymap,self).__init__(projection='mbtfpq', **kwargs)

class OrthoSkymap(Skymap):
    defaults = dict(Skymap.defaults)
    def __init__(self,**kwargs):
        setdefaults(kwargs,self.defaults)
        super(OrthoSkymap,self).__init__(projection='ortho', **kwargs)
