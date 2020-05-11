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
import pylab as plt
import numpy as np
import ephem
import healpy as hp
import scipy.ndimage as nd

from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

import cartopy.crs as ccrs

from cartosky.utils import setdefaults,get_datadir
from cartosky.utils import cel2gal, gal2cel
from cartosky import healpix

import cartosky.proj

class Skymap(object):
    """ Base class for creating Skymap objects. """

    defaults = {'global':True, 'gridlines':True, 'lon_0': 0}

    def __init__(self, projection='cyl', **kwargs):
        setdefaults(kwargs,self.defaults)
        do_global    = kwargs.pop('global',True)
        do_grid      = kwargs.pop('gridlines',True)
        do_celestial = kwargs.pop('celestial',False) #noop

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
            self.grid = self.ax.gridlines(linestyle=':')

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
        _ = healpix.check_hpxmap(hpxmap,None,None)
        hpxmap = healpix.masked_array(hpxmap,badval)
        hpxmap.fill_value = np.ma.median(hpxmap)
        smooth = hp.smoothing(hpxmap,sigma=np.radians(sigma),verbose=False)
        return np.ma.array(smooth,mask=hpxmap.mask)

    def draw_hpxmap(self, hpxmap, pixel=None, nside=None, xsize=800,
                    lonra=None, latra=None, badval=hp.UNSEEN, smooth=None,
                    **kwargs):
        """
        Use pcolor/pcolormesh to draw healpix map.

        Parameters:
        -----------
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

        Returns:
        --------
        im,lon,lat,values : mpl image with pixel longitude, latitude (deg), and values
        """

        is_hsp = healpix.check_hpxmap(hpxmap,pixel,nside)
        if is_hsp:
            if nside is not None and nside < hpxmap.nside_sparse:
                hpxmap.degrade(nside)
            elif nside is not None and nside > hpxmap.nside_sparse:
                return NotImplementedError("Map upgrading is not implemented, try with different nside")
            else:
                nside = hpxmap.nside_sparse
            if pixel is not None:
                hpxmap = hpxmap.get_values_pix(hp.ring2nest(nside, pixel))
            else:
                _hpxmap = hpxmap.get_values_pix(hpxmap.valid_pixels)
                pixel = hp.nest2ring(nside, hpxmap.valid_pixels)
                hpxmap = _hpxmap
        hpxmap = healpix.masked_array(hpxmap,badval)

        if smooth:
            # To smooth we need the full map
            hpxmap = healpix.create_map(hpxmap,pixel,nside,badval)
            pixel,nside = None,None
            hpxmap = healpix.masked_array(hpxmap,badval)
            hpxmap = self.smooth(hpxmap,sigma=smooth)

        vmin,vmax = np.percentile(hpxmap.compressed(),[2.5,97.5])

        defaults = dict(rasterized=True, vmin=vmin, vmax=vmax,
                        transform=ccrs.PlateCarree())
        setdefaults(kwargs,defaults)

        lon,lat,values = healpix.hpx2xy(hpxmap,pixel=pixel,nside=nside,
                                        xsize=xsize,
                                        lonra=lonra,latra=latra)

        im = self.ax.pcolormesh(lon,lat,values,**kwargs)
        return im,lon,lat,values

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
        try:
            pix = hp.ang2pix(nside,lon,lat,lonlat=True)
        except TypeError:
            pix = hp.ang2pix(nside,np.radians(90-lat),np.radians(lon))

        npix = hp.nside2npix(nside)
        hpxmap = hp.UNSEEN*np.ones(npix)
        idx,cts = np.unique(pix,return_counts=True)
        hpxmap[idx] = cts

        return hpxmap,self.draw_hpxmap(hpxmap,**kwargs)


    def tissot(self, lons=None, lats=None, rad_deg=1.0, n_samples=80, **kwargs):
        """
        Add Tissot's indicatrices to the axes.
        Parameters
        ----------
        lons
            A numpy.ndarray, list or tuple of longitude values that
            locate the centre of each circle. Specifying more than one
            dimension allows individual points to be drawn whereas a
            1D array produces a grid of points.
        lats
            A numpy.ndarray, list or tuple of latitude values that
            that locate the centre of each circle. See lons.
        rad_deg
            The radius in deg of the the circles to be drawn.
        n_samples
            Integer number of points sampled around the circumference of
            each circle.
        ``**kwargs`` are passed through to `class:ShapelyFeature`.
        """
        # This will need to be re-implemented...
        # cartosky.GeoAxes.tissot makes some assumptions about Geodesic radius...
        msg = "cartosky.GeoAxes.tissot needs to be re-implemented"
        #raise Exception(msg)
        warnings.warn(msg)

    @staticmethod
    def wrap_index(lon, lat, wrap=180.):
        """ Find the index where the array wraps.
        """
        # No wrap: ignore
        if wrap is None:  return None

        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)

        # No array: ignore
        if len(lon)==1 or len(lat)==1: return None

        # Map [0,360)
        lon = np.mod(lon,360)
        wrap = np.mod(wrap,360)

        # Find the index of the entry closest to the wrap angle
        idx = np.abs(lon - wrap).argmin()
        # First or last index: ignore
        if idx == 0 or idx+1 == len(lon): return None
        # Value exactly equals wrap, choose next value
        elif (lon[idx] == wrap): idx += 1
        # Wrap angle sandwiched
        elif (lon[idx]<wrap) and (lon[idx+1]>wrap): idx += 1
        elif (lon[idx]<wrap) and (lon[idx-1]>wrap): idx += 0
        elif (lon[idx]>wrap) and (lon[idx+1]<wrap): idx += 1
        elif (lon[idx]>wrap) and (lon[idx-1]<wrap): idx += 0
        # There is no wrap: ignore
        else: return None

        return idx

    @classmethod
    def roll(cls,lon,lat,wrap=180.):
        """ Roll an lon,lat combination to split 180 boundary
        Parameters:
        -----------
        lon : right ascension (deg)
        lat: declination (deg)
        wrap_angle : angle to wrap at (deg)
        """
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)

        # Do nothing
        if wrap is None: return lon,lat
        if len(lon)==1 or len(lat)==1: return lon,lat

        idx = cls.wrap_index(lon,lat,wrap)
        if idx is None: return lon, lat

        return np.roll(lon,-idx), np.roll(lat,-idx)

    @classmethod
    def split(cls,lon,lat,wrap=180.):
        """ Split an lon,lat combination into lists across a wrap boundary
        Parameters:
        -----------
        lon : right ascension (deg)
        lat: declination (deg)
        wrap_angle : angle to wrap at (deg)
        """
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)

        # Do nothing
        if wrap is None: return [lon],[lat]
        if len(lon)==1 or len(lat)==1: return [lon],[lat]

        idx = cls.wrap_index(lon,lat,wrap)
        if idx is None: return [lon], [lat]

        return np.split(lon,[idx]), np.split(lat,[idx])


    def draw_polygon(self,filename,**kwargs):
        """ Draw a polygon footprint. """
        defaults=dict(color='k', lw=2)
        setdefaults(kwargs,defaults)

        poly = np.loadtxt(filename,dtype=[('ra',float),('dec',float)])
        return self.draw_polygon_radec(poly['ra'],poly['dec'],**kwargs)

    def draw_polygon_radec(self,ra,dec,**kwargs):
        defaults=dict(transform=ccrs.PlateCarree())
        setdefaults(kwargs,defaults)
        #xy = self.proj(*self.roll(ra,dec,self.wrap_angle))
        return self.ax.plot(ra,dec,**kwargs)

    def draw_polygons(self,filename,**kwargs):
        """Draw a text file containing multiple polygons"""
        data = np.genfromtxt(filename,names=['ra','dec','poly'])

        ret = []
        for p in np.unique(data['poly']):
            poly = data[data['poly'] == p]
            xy = self.draw_polygon_radec(poly['ra'],poly['dec'],**kwargs)
            ret += [xy]
            kwargs.pop('label',None)

        return ret

    def draw_paths(self,filename,**kwargs):
        """Draw a text file containing multiple polygons"""
        try:
            data = np.genfromtxt(filename,names=['ra','dec','poly'])
        except ValueError:
            data = np.genfromtxt(filename,names=['ra','dec'])
            data = mlab.rec_append_fields(data,'poly',np.zeros(len(data)))

        ret = []
        for p in np.unique(data['poly']):
            poly = data[data['poly'] == p]
            path,patch = self.draw_path_radec(poly['ra'],poly['dec'],**kwargs)
            ret += [(path,patch)]
        return ret

    def draw_path_radec(self,ra,dec,**kwargs):
        #xy = self.proj(*self.roll(ra,dec,self.wrap_angle))
        #vertices = np.vstack(xy).T
        defaults=dict(transform=ccrs.PlateCarree())
        setdefaults(kwargs,defaults)

        vertices = np.vstack([ra,dec]).T
        path = matplotlib.path.Path(vertices)
        patch = matplotlib.patches.PathPatch(path,**kwargs)
        plt.gca().add_artist(patch)
        return path,patch

    def draw_zenith(self, radius=1.0, **kwargs):
        """
        Plot a to-scale representation of the zenith.
        """
        defaults = dict(color='green',alpha=0.75,lw=1.5)
        setdefaults(kwargs,defaults)

        # RA and Dec of zenith
        zra, zdec = np.degrees(self.observer.radec_of(0, '90'))
        xy = self.proj(zra, zdec)

        self.plot(*xy,marker='+',ms=10,mew=1.5, **kwargs)
        if radius:
            self.tissot(zra,zdec,radius,npts=100,fc='none', **kwargs)

    def draw_airmass(self, airmass=1.4, npts=360, **kwargs):
        defaults = dict(color='green', lw=2)
        setdefaults(kwargs,defaults)

        altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)
        ra_contour = np.zeros(npts)
        dec_contour = np.zeros(npts)
        for ii, azimuth in enumerate(np.linspace(0., 2. * np.pi, npts)):
            ra_radians, dec_radians = self.observer.radec_of(azimuth, '%.2f'%(np.degrees(altitude_radians)))
            ra_contour[ii] = np.degrees(ra_radians)
            dec_contour[ii] = np.degrees(dec_radians)
        xy = self.proj(ra_contour, dec_contour)
        self.plot(*xy, **kwargs)

        self.draw_zenith(**kwargs)

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
        defaults = dict(color='k',lw=1.5,ls='-')
        setdefaults(kwargs,defaults)

        glon = np.linspace(0,360,500)
        glat = np.zeros_like(glon)
        ra,dec = self.roll(*gal2cel(glon,glat),wrap=self.wrap_angle)
        ra -= 360*(ra > 180)

        self.draw_polygon_radec(ra,dec,**kwargs)

        if width:
            kwargs.update(dict(ls='--',lw=1))
            for delta in [+width,-width]:
                ra,dec = self.roll(*gal2cel(glon,glat+delta))
                ra -= 360*(ra > 180)
                self.draw_polygon_radec(ra,dec,**kwargs)

    def draw_lmc(self,**kwargs):
        from cartosky.constants import RA_LMC, DEC_LMC, RADIUS_LMC
        defaults = dict(npts=100,fc='0.7',ec='0.5')
        setdefaults(kwargs,defaults)
        proj = self.proj(RA_LMC,DEC_LMC)
        #self.tissot(RA_LMC,DEC_LMC,RADIUS_LMC,**kwargs)
        plt.text(proj[0],proj[1], 'LMC', weight='bold',
                 fontsize=10, ha='center', va='center', color='k')

    def draw_smc(self,**kwargs):
        from cartosky.constants import RA_SMC, DEC_SMC, RADIUS_SMC
        defaults = dict(npts=100,fc='0.7',ec='0.5')
        setdefaults(kwargs,defaults)
        proj = self.proj(RA_SMC,DEC_SMC)
        #self.tissot(RA_SMC,DEC_SMC,RADIUS_SMC,**kwargs)
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
            corners = decam.project(self,_ra,_dec)
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
