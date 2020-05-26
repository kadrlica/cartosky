#!/usr/bin/env python3
"""
New cartopy projection classes.

Taken from ProPlot
"""
import warnings
from cartopy.crs import _WarpedRectangularProjection
from cartopy._crs import Globe

from .utils import setdefaults

__all__ = [
    'Aitoff', 'Hammer', 'SkySphere'
]

# Default radius of the Globe.
# Eventually should be possible to use the unit sphere.
#RADIUS=1.0e6
RADIUS=6378137.0 # meters

class SkySphere(Globe):
    """
    Spherical ccrs.Globe for sky plotting.
    """

    def __init__(self,*args,**kwargs):
        defaults = dict(ellipse=None,semimajor_axis=RADIUS,semiminor_axis=RADIUS)
        kwargs = setdefaults(kwargs,defaults)
        super().__init__(*args,**kwargs)

class Aitoff(_WarpedRectangularProjection):
    """
    The `Aitoff <https://en.wikipedia.org/wiki/Aitoff_projection>`__
    projection.
    """
    #: Registered projection name.
    name = 'aitoff'

    def __init__(
        self, central_longitude=0, central_latitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        """
        if globe is None:
            globe = SkySphere()

        if globe.ellipse is not None:
            warnings.warn(
                'The {} projection does not handle elliptical globes.'.format(self.name)
            )

        proj4_params = {'proj': 'aitoff', 'lon_0': central_longitude,
                        'lat_0': central_latitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @property
    def threshold(self):  # how finely to interpolate line data, etc.
        return 1e4

class Hammer(_WarpedRectangularProjection):
    """
    The `Hammer <https://en.wikipedia.org/wiki/Hammer_projection>`__
    projection.
    """
    #: Registered projection name.
    name = 'hammer'

    def __init__(
        self, central_longitude=0, central_latitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        """
        if globe is None:
            globe = SkySphere()

        if globe.ellipse is not None:
            warnings.warn(
                f'The {self.name!r} projection does not handle '
                'elliptical globes.'
            )

        proj4_params = {'proj': 'hammer', 'lon_0': central_longitude,
                        'lat_0': central_latitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @property
    def threshold(self):  # how finely to interpolate line data, etc.
        """
        """
        return 1e4


class McBrydeThomasFlatPolarQuartic(_WarpedRectangularProjection):
    """
    The `McBryde-Thomas Flat Polar Quartic <https://it.wikipedia.org/wiki/File:McBryde-Thomas_flat-pole_quartic_projection_SW.jpg>`__
    projection.
    """
    #: Registered projection name.
    name = 'mbtfpq'

    def __init__(
        self, central_longitude=0, central_latitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        """
        if globe is None:
            globe = SkySphere()

        if globe.ellipse is not None:
            warnings.warn(
                f'The {self.name!r} projection does not handle '
                'elliptical globes.'
            )

        proj4_params = {'proj': 'mbtfpq', 'lon_0': central_longitude,
                        'lat_0': central_latitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @property
    def threshold(self):  # how finely to interpolate line data, etc.
        """
        """
        return 1e4

class KavrayskiyVII(_WarpedRectangularProjection):
    """
    The `Kavrayskiy VII \
<https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection>`__ projection.
    """
    #: Registered projection name.
    name = 'kavrayskiyVII'

    def __init__(
        self, central_longitude=0, globe=None,
        false_easting=None, false_northing=None
    ):
        """
        %(proj.init)s
        """
        if globe is None:
            globe = SkySphere()

        if globe.ellipse is not None:
            warnings.warn(
                f'The {self.name!r} projection does not handle '
                'elliptical globes.'
            )

        proj4_params = {'proj': 'kav7', 'lon_0': central_longitude}
        super().__init__(
            proj4_params, central_longitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe
        )

    @property
    def threshold(self):
        """
        %(proj.reso)s
        """
        return 1e5
