#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import functools
import matplotlib.ticker as mticker
import numpy as np

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.geodesic import Geodesic
from shapely.geometry.polygon import Polygon
from shapely.geometry.point import Point

from cartosky.utils import setdefaults
from cartosky import healpix

def default_transform(self, func, *args, transform=None, **kwargs):
    """
    Makes ``transform=cartopy.crs.PlateCarree()`` the default
    for cartopy plots. This means you no longer have to
    pass ``transform=cartopy.crs.PlateCarree()`` if your data
    coordinates are longitude and latitude.
    Note
    ----
    This function wraps %(methods)s for `~proplot.axes.CartopyAxes`.
    Comes from ProPlot
    """
    # Apply default transform
    # TODO: Do some cartopy methods reset backgroundpatch or outlinepatch?
    # Deleted comment reported this issue
    if transform is None:
        transform = ccrs.PlateCarree()
    result = func(self, *args, transform=transform, **kwargs)
    return result


def default_crs(self, func, *args, crs=None, **kwargs):
    """
    Fixes the `~cartopy.mpl.geoaxes.GeoAxes.set_extent` bug associated with
    tight bounding boxes and makes ``crs=cartopy.crs.PlateCarree()`` the
    default for cartopy plots.
    Note
    ----
    This function wraps %(methods)s for `~proplot.axes.CartopyAxes`.
    Comes from ProPlot
    """
    # Apply default crs
    name = func.__name__
    if crs is None:
        crs = ccrs.PlateCarree()
    try:
        result = func(self, *args, crs=crs, **kwargs)
    except TypeError as err:  # duplicate keyword args, i.e. crs is positional
        if not args:
            raise err
        args, crs = args[:-1], args[-1]
        result = func(self, *args, crs=crs, **kwargs)

    return result

def default_gridlines(self, func, *args, **kwargs):
    """
    Default gridlines options.
    """
    fmt = mticker.FuncFormatter(lambda v, pos: '{:g}'.format(v))
    defaults = dict(linestyle=':',xformatter=fmt,yformatter=fmt, draw_labels=True)
    setdefaults(kwargs,defaults)
    result = func(self, *args, **kwargs)
    return result


def ignore_warnings(self, func, *args, **kwargs):
    """
    Turn off warnings for this function.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore",UserWarning)
        result = func(self, *args, **kwargs)

    return result

def _wrapper_decorator(driver):
    """
    Generate generic wrapper decorator and dynamically modify the docstring
    to list methods wrapped by this function. Also set `__doc__` to ``None`` so
    that ProPlot fork of automodapi doesn't add these methods to the website
    documentation. Users can still call help(ax.method) because python looks
    for superclass method docstrings if a docstring is empty.
    """
    driver._docstring_orig = driver.__doc__ or ''
    driver._methods_wrapped = []

    def decorator(func):
        # Define wrapper and suppress documentation
        # We only document wrapper functions, not the methods they wrap
        @functools.wraps(func)
        def _wrapper(self, *args, **kwargs):
            return driver(self, func, *args, **kwargs)
        name = func.__name__
        _wrapper.__doc__ = None
        return _wrapper
    return decorator

_default_transform = _wrapper_decorator(default_transform)
_default_crs       = _wrapper_decorator(default_crs)
_default_gridlines = _wrapper_decorator(default_gridlines)
_ignore_warnings   = _wrapper_decorator(ignore_warnings)

class SkyAxes(GeoAxes):
    """ Subclass cartopy.GeoAxes """

    def tissot(self, ra, dec, radius=1.0, n_samples=100, **kwargs):
        """
        Add Tissot's indicatrix to the axes.
        https://en.wikipedia.org/wiki/Tissot%27s_indicatrix

        Parameters
        ----------
        ra    : right asencion of circle center (deg)
        dec   : declination of circle center (deg)
        radius: angular radius (deg)
        n_samples : number of samples
        ``**kwargs`` are passed to `SkyAxes.add_geometries`.

        Returns
        -------
        feat  : feature artist
        """
        ra  = np.asarray(ra).flatten()
        dec = np.asarray(dec).flatten()

        if ra.shape != dec.shape:
            msg = "ra and dec must have the same shape."
            raise ValueError(msg)

        globe = self.projection.globe

        params = globe.to_proj4_params()
        a,b = params['a'],params['b']
        geod = Geodesic(radius=a,flattening=1 - b/a)
        radius_m = a * np.radians(radius)

        geoms = []
        for x, y in zip(ra, dec):
            circle = geod.circle(x, y, radius_m, n_samples=n_samples)
            geoms.append(Polygon(circle))

        return self.add_geometries(geoms, ccrs.Geodetic(globe), **kwargs)

    def tissot_indicatrices(self, radius=5.0, num=6, **kwargs):
        """
        Add Tissot's indicatrices to the axes.
        https://en.wikipedia.org/wiki/Tissot%27s_indicatrix

        Parameters
        ----------
        radius : angular radius of circles (deg)
        ``**kwargs`` are passed to `Skymap.tissot`.

        Returns
        -------
        feat : feature artist
        """
        ra = np.linspace(-180, 180, num, endpoint=False)
        dec = np.linspace(-80, 80, num)
        return self.tissot(*np.meshgrid(ra,dec),radius=radius,**kwargs)

    # Change gridlines default to ls=':'
    gridlines = _default_gridlines(
        GeoAxes.gridlines
    )

    plot = _default_transform(
        GeoAxes.plot
    )
    scatter = _default_transform(
        GeoAxes.scatter
    )
    fill_between = _default_transform(
        GeoAxes.fill_between
    )
    fill_betweenx = _default_transform(
        GeoAxes.fill_betweenx
    )
    contour = _default_transform(
        GeoAxes.contour
    )
    contourf = _default_transform(
        GeoAxes.contourf
    )
    pcolor = _default_transform(
        GeoAxes.pcolor
    )
    pcolormesh = _default_transform(
        GeoAxes.pcolormesh
    )
    quiver = _default_transform(
        GeoAxes.quiver
    )
    streamplot = _default_transform(
        GeoAxes.streamplot
    )
    barbs = _default_transform(
        GeoAxes.barbs
    )
    tripcolor = _default_transform(
        GeoAxes.tripcolor
    )
    tricontour = _default_transform(
        GeoAxes.tricontour
    )
    tricontourf = _default_transform(
        GeoAxes.tricontourf
    )
    get_extent = _ignore_warnings(
        GeoAxes.get_extent
    )
    #set_extent = _default_crs(
    #    GeoAxes.set_extent
    #)
    set_xticks = _default_crs(
        GeoAxes.set_xticks
    )
    set_yticks = _default_crs(
        GeoAxes.set_yticks
    )

    def hpxmap(self, hpxmap, pixel=None, nside=None, xsize=800,
                    lonra=None, latra=None, badval=healpix.UNSEEN, smooth=None,
                    **kwargs):
        """ Draw a healpix map with pcolormesh.

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
        # ADW: probably still not the best way to do this...
        try:
            import healsparse as hsp
            if isinstance(hpxmap, hsp.HealSparseMap):
                hpxmap,pixel,nside = healpix.hsp2hpx(hpxmap)
        except ImportError:
            pass

        healpix.check_hpxmap(hpxmap,pixel,nside)
        hpxmap = healpix.masked_array(hpxmap,badval)

        if smooth:
            # To smooth we need the full map...
            # It'd be good to check we aren't going to blow anything up
            hpxmap = healpix.create_map(hpxmap,pixel,nside,badval)
            pixel,nside = None,None
            hpxmap = healpix.masked_array(hpxmap,badval)
            hpxmap = healpix.smooth(hpxmap,sigma=smooth)

        vmin,vmax = np.percentile(hpxmap.compressed(),[2.5,97.5])

        defaults = dict(rasterized=True, vmin=vmin, vmax=vmax)
        setdefaults(kwargs,defaults)

        lon,lat,values = healpix.hpx2xy(hpxmap,pixel=pixel,nside=nside,
                                        xsize=xsize,
                                        lonra=lonra,latra=latra)

        im = self.pcolormesh(lon,lat,values,**kwargs)
        self._sci(im)
        return im,lon,lat,values

    def hpxbin(self, lon, lat, nside=256, **kwargs):
        """
        Create a healpix histogram of the counts.

        Like `hexbin` from matplotlib

        Parameters:
        -----------
        lon : input longitude (deg)
        lat : input latitude (deg)
        nside : heaplix nside resolution
        kwargs : passed to SkyAxes.hpxmap

        Returns:
        --------
        hpxmap, im : healpix map and image
        """
        hpxmap = healpix.hpxbin(lon,lat,nside)
        return hpxmap, self.hpxmap(hpxmap,**kwargs)
