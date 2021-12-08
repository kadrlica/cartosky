#!/usr/bin/env python
"""
Extension for individual surveys.
"""
import os

import numpy as np
import pylab as plt
from collections import OrderedDict as odict

from cartosky.utils import setdefaults, get_datafile
from cartosky.core import Skymap, McBrydeSkymap, OrthoSkymap
from cartosky.constants import DECAM
from cartosky import healpix

# Derived from telra,teldec of 10000 exposures
DES_SN = odict([
    ('E1', dict(ra=7.874, dec=-43.010)),
    ('E2', dict(ra=9.500, dec=-43.999)),
    ('X1', dict(ra=34.476, dec=-4.931)),
    ('X2', dict(ra=35.664, dec=-6.413)),
    ('X3', dict(ra=36.449, dec=-4.601)),
    ('S1', dict(ra=42.818, dec=0.000)),
    ('S2', dict(ra=41.193, dec=-0.991)),
    ('C1', dict(ra=54.274, dec=-27.113)),
    ('C2', dict(ra=54.274, dec=-29.090)),
    ('C3', dict(ra=52.647, dec=-28.101)),
])

DES_SN_LABELS = odict([
    ('SN-E', dict(ra=15, dec=-38, ha='center')),
    ('SN-X', dict(ra=35, dec=-13, ha='center')),
    ('SN-S', dict(ra=55, dec=0, ha='center')),
    ('SN-C', dict(ra=57, dec=-36, ha='center')),
])


class SurveySkymap(Skymap):
    """Extending to survey specific functions.
    """
    # Footprint drawing
    def draw_footprint(self, filename, **kwargs):
        """ Draw survey footpring polygon

        Parameters
        ----------
        filename : path to polygon file to draw
        **kwargs : passed to draw_polygon

        Returns
        -------
        poly    : polygon
        """
        defaults = dict(edgecolor='k', facecolor='none', lw=2)
        setdefaults(kwargs, defaults)
        self.draw_polygon(filename, **kwargs)

    def draw_maglites(self, **kwargs):
        """Draw the MagLiteS footprint """
        defaults = dict(edgecolor='blue')
        setdefaults(kwargs, defaults)

        filename = get_datafile('maglites-poly.txt')
        self.draw_footprint(filename, **kwargs)

    def draw_bliss(self, **kwargs):
        """Draw the BLISS footprint"""
        defaults = dict(edgecolor='magenta')
        setdefaults(kwargs, defaults)

        filename = get_datafile('bliss-poly.txt')
        self.draw_footprint(filename, **kwargs)

    def draw_des(self, **kwargs):
        """ Draw the DES footprint. """
        return self.draw_des17(**kwargs)

    def draw_des13(self, **kwargs):
        """ Draw the DES footprint. """
        defaults = dict(edgecolor='red')
        setdefaults(kwargs, defaults)

        filename = get_datafile('des-round13-poly.txt')
        return self.draw_footprint(filename, **kwargs)

    def draw_des17(self, **kwargs):
        """ Draw the DES footprint. """
        defaults = dict(edgecolor='red')
        setdefaults(kwargs, defaults)

        filename = get_datafile('des-round17-poly.txt')
        return self.draw_footprint(filename, **kwargs)

    def draw_decals(self, **kwargs):
        """ Draw DECaLS footprint. """
        defaults = dict(edgecolor='blue')
        setdefaults(kwargs, defaults)

        filename = get_datafile('decals-poly.txt')
        return self.draw_footprint(filename, **kwargs)

    # Tissot drawing
    def draw_des_sn(self, **kwargs):
        defaults = dict(facecolor='none', edgecolor='k', lw=1, zorder=10)
        setdefaults(kwargs, defaults)

        ra = [v['ra'] for v in DES_SN.values()]
        dec = [v['dec'] for v in DES_SN.values()]
        self.tissot(ra, dec, DECAM, **kwargs)

    def draw_smash(self, **kwargs):
        """ Draw the SMASH fields. """
        defaults = dict(facecolor='none', edgecolor='k', lw=1, zorder=10)
        setdefaults(kwargs, defaults)

        filename = get_datafile('smash_fields_final.txt')
        smash = np.genfromtxt(filename, dtype=[('ra', float), ('dec', float)], usecols=[4, 5])
        self.tissot(smash['ra'], smash['dec'], DECAM, **kwargs)

    # Map drawing
    def draw_sfd(self, **kwargs):
        import healpy as hp
        from matplotlib.colors import LogNorm
        defaults = dict(rasterized=True, cmap=plt.cm.binary, norm=LogNorm())
        setdefaults(kwargs, defaults)

        filename = get_datafile('lambda_sfd_ebv.fits')
        if not os.path.exists(filename):
            import subprocess
            url = "https://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits"
            print("Downloading SFD map from:\n  %s"%url)
            subprocess.check_output(['curl', url, '--output', filename])

        galhpx = hp.read_map(filename, verbose=False)
        celhpx = healpix.gal2cel(galhpx)
        return self.draw_hpxmap(celhpx, **kwargs)


class SurveyMcBryde(SurveySkymap, McBrydeSkymap):
    pass


class SurveyOrtho(SurveySkymap, OrthoSkymap):
    pass
