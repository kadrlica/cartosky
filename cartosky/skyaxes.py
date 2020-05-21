#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import functools
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes

from cartosky.utils import setdefaults

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

class SkyAxes(GeoAxes):
    """ Subclass cartopy.GeoAxes """

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
    get_extent = _default_crs(
        GeoAxes.get_extent
    )
    set_extent = _default_crs(
        GeoAxes.set_extent
    )
    set_xticks = _default_crs(
        GeoAxes.set_xticks
    )
    set_yticks = _default_crs(
        GeoAxes.set_yticks
    )

