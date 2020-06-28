#!/usr/bin/env python
"""
Generic python script.

Much of this is taken from ProPlot constructor:
https://github.com/lukelbd/proplot

but removing basemap support (and making python 2 compatible)

"""
__author__ = "Alex Drlica-Wagner"
import warnings

from . import crs as pcrs
from .utils import setdefaults
import cartopy.crs as ccrs
CRS = ccrs.CRS

# Mapping of "projection names" to cartopy `~cartopy.crs.Projection` classes.
CARTOPY_PROJS = {}
CARTOPY_PROJS.update({
    'aitoff': pcrs.Aitoff,
    'hammer': pcrs.Hammer,
    'mbtfpq': pcrs.McBrydeThomasFlatPolarQuartic,
    'kav7'  : pcrs.KavrayskiyVII,
})
CARTOPY_PROJS_UNAVAIL = {
    'aea': 'AlbersEqualArea',
    'aeqd': 'AzimuthalEquidistant',
    'cyl': 'PlateCarree',  # only basemap name not matching PROJ
    'eck1': 'EckertI',
    'eck2': 'EckertII',
    'eck3': 'EckertIII',
    'eck4': 'EckertIV',
    'eck5': 'EckertV',
    'eck6': 'EckertVI',
    'eqc': 'PlateCarree',  # actual PROJ name
    'eqdc': 'EquidistantConic',
    'eqearth': 'EqualEarth',  # better looking Robinson; not in basemap
    'euro': 'EuroPP',  # Europe; not in basemap or PROJ
    'geos': 'Geostationary',
    'gnom': 'Gnomonic',
    'igh': 'InterruptedGoodeHomolosine',  # not in basemap
    'laea': 'LambertAzimuthalEqualArea',
    'lcc': 'LambertConformal',
    'lcyl': 'LambertCylindrical',  # not in basemap or PROJ
    'merc': 'Mercator',
    'mill': 'Miller',
    'moll': 'Mollweide',
    'npstere': 'NorthPolarStereo',  # np/sp stuff not in PROJ
    'nsper': 'NearsidePerspective',
    'ortho': 'Orthographic',
    'osgb': 'OSGB',  # UK; not in basemap or PROJ
    'osni': 'OSNI',  # Ireland; not in basemap or PROJ
    'pcarree': 'PlateCarree',  # common alternate name
    'robin': 'Robinson',
    'rotpole': 'RotatedPole',
    'sinu': 'Sinusoidal',
    'spstere': 'SouthPolarStereo',
    'stere': 'Stereographic',
    'tmerc': 'TransverseMercator',
    'utm': 'UTM',  # not in basemap
}
for _key, _cls in list(CARTOPY_PROJS_UNAVAIL.items()):
    if hasattr(ccrs, _cls):
        CARTOPY_PROJS[_key] = getattr(ccrs, _cls)
        del CARTOPY_PROJS_UNAVAIL[_key]
if CARTOPY_PROJS_UNAVAIL:
    warnings._warn_proplot(
        'Cartopy projection(s) '
        + ', '.join(map(repr, CARTOPY_PROJS_UNAVAIL.values()))
        + 'are unavailable. Consider updating to cartopy >= 0.17.0.'
    )

CARTOPY_KW_ALIASES = {  # use PROJ shorthands instead of verbose cartopy names
    'lat_0': 'central_latitude',
    'lon_0': 'central_longitude',
    'lat_min': 'min_latitude',
    'lat_max': 'max_latitude',
}

CARTOPY_KW_DEFAULTS = {
    'globe': pcrs.SkySphere()
}

def Proj(name, **kwargs):
    """
    Return a `cartopy.crs.Projection`.

    Parameters
    ----------
    name : str, `cartopy.crs.Projection`
        The projection name or projection class instance. If the latter, it
        is simply returned. If the former, it must correspond to one of the
        `PROJ <https://proj.org>`__ projection name shorthands, like in
        basemap.
    Other parameters
    ----------------
    **kwargs
        Passed to cartopy `~cartopy.crs.Projection` class.
        Translates `lon_0` and `lat_0` to `central_longitude` and
        `central_latitude`.
    Returns
    -------
    proj : `~cartopy.crs.Projection`
        The projection instance.
    See also
    --------
    For more information on map projections, see the
    `PROJ <https://proj.org>`__ documentation.
    """
    # Class instances
    is_crs = isinstance(name, CRS)
    if is_crs:
        proj = name

    # Invalid
    elif not isinstance(name, str):
        raise ValueError(
            f'Unexpected Proj() argument {name!r}. '
            'Must be name or cartopy.crs.CRS instance.'
        )

    else:
        kwproj = {
            CARTOPY_KW_ALIASES.get(key, key): value
            for key, value in kwargs.items()
        }
        kwproj = setdefaults(kwproj,CARTOPY_KW_DEFAULTS)
        crs = CARTOPY_PROJS.get(name, None)
        if name == 'geos':  # fix common mistake
            kwproj.pop('central_latitude', None)
        if crs is None:
            raise ValueError(
                f'Unknown projection {name!r}. Options are: '
                + ', '.join(map(repr, CARTOPY_PROJS.keys())) + '.'
            )

        # Monkey patch the axis returned...
        # Might be better to do this somewhere else?
        def _as_mpl_axes(self):
            import cartosky.skyaxes as skyaxes
            return skyaxes.SkyAxes, {'map_projection': self}

        crs._as_mpl_axes = _as_mpl_axes

        proj = crs(**kwproj)

    return proj
