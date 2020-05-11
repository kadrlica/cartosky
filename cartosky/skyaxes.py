#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import functools

from cartopy.mpl import geoaxes
import cartopy.crs as ccrs

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
    # Fix extent, so axes tight bounding box gets correct box!
    # From this issue:
    # https://github.com/SciTools/cartopy/issues/1207#issuecomment-439975083
    if name == 'set_extent':
        clipped_path = self.outline_patch.orig_path.clip_to_bbox(self.viewLim)
        self.outline_patch._path = clipped_path
        self.background_patch._path = clipped_path
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

class SkyAxes(geoaxes.GeoAxes):
    """ Subclass cartopy.GeoAxes """

    # Change gridlines default to ls=':'

    plot = _default_transform(
        geoaxes.GeoAxes.plot
    )
    scatter = _default_transform(
        geoaxes.GeoAxes.scatter
    )
    fill_between = _default_transform(
        geoaxes.GeoAxes.fill_between
    )
    fill_betweenx = _default_transform(
        geoaxes.GeoAxes.fill_betweenx
    )
    contour = _default_transform(
        geoaxes.GeoAxes.contour
    )
    contourf = _default_transform(
        geoaxes.GeoAxes.contourf
    )
    pcolor = _default_transform(
        geoaxes.GeoAxes.pcolor
    )
    pcolormesh = _default_transform(
        geoaxes.GeoAxes.pcolormesh
    )
    quiver = _default_transform(
        geoaxes.GeoAxes.quiver
    )
    streamplot = _default_transform(
        geoaxes.GeoAxes.streamplot
    )
    barbs = _default_transform(
        geoaxes.GeoAxes.barbs
    )
    tripcolor = _default_transform(
        geoaxes.GeoAxes.tripcolor
    )
    tricontour = _default_transform(
        geoaxes.GeoAxes.tricontour
    )
    tricontourf = _default_transform(
        geoaxes.GeoAxes.tricontourf
    )
    get_extent = _default_crs(
        geoaxes.GeoAxes.get_extent
    )
    set_extent = _default_crs(
        geoaxes.GeoAxes.set_extent
    )
    set_xticks = _default_crs(
        geoaxes.GeoAxes.set_xticks
    )
    set_yticks = _default_crs(
        geoaxes.GeoAxes.set_yticks
    )

