#!/usr/bin/env python
"""
Extension for individual surveys.
"""
import os

import numpy as np
import pylab as plt
import pandas as pd
from collections import OrderedDict as odict

from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot
import mpl_toolkits.axisartist as axisartist
import  mpl_toolkits.axisartist.angle_helper as angle_helper

from cartosky.utils import setdefaults,get_datadir,hpx_gal2cel
from cartosky.core import Skymap,McBrydeSkymap,OrthoSkymap
from cartosky.constants import DECAM

# Derived from telra,teldec of 10000 exposures
DES_SN = odict([
    ('E1',dict(ra=7.874,  dec=-43.010)),
    ('E2',dict(ra=9.500,  dec=-43.999)),
    ('X1',dict(ra=34.476, dec=-4.931 )),
    ('X2',dict(ra=35.664, dec=-6.413 )),
    ('X3',dict(ra=36.449, dec=-4.601 )),
    ('S1',dict(ra=42.818, dec=0.000  )),
    ('S2',dict(ra=41.193, dec=-0.991 )),
    ('C1',dict(ra=54.274, dec=-27.113)),
    ('C2',dict(ra=54.274, dec=-29.090)),
    ('C3',dict(ra=52.647, dec=-28.101)),
])

DES_SN_LABELS = odict([
    ('SN-E',   dict(ra=15, dec=-38, ha='center')),
    ('SN-X',   dict(ra=35, dec=-13, ha='center')),
    ('SN-S',   dict(ra=55, dec=0,   ha='center')),
    ('SN-C',   dict(ra=57, dec=-36, ha='center')),
])


class SurveySkymap(Skymap):
    """Extending to survey specific functions.
    """
    def draw_maglites(self,**kwargs):
        """Draw the MagLiteS footprint"""
        defaults=dict(edgecolor='blue', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(get_datadir(),'maglites-poly.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_bliss(self,**kwargs):
        """Draw the BLISS footprint"""
        defaults=dict(edgecolor='magenta', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(get_datadir(),'bliss-poly.txt')
        self.draw_polygons(filename,**kwargs)

    def draw_des(self,**kwargs):
        """ Draw the DES footprint. """
        defaults=dict(edgecolor='red', lw=2)
        setdefaults(kwargs,defaults)

        return self.draw_des17(**kwargs)

    def draw_des13(self,**kwargs):
        """ Draw the DES footprint. """
        defaults=dict(edgecolor='red', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(get_datadir(),'des-round13-poly.txt')
        return self.draw_polygon(filename,**kwargs)

    def draw_des17(self,**kwargs):
        """ Draw the DES footprint. """
        defaults=dict(edgecolor='red', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(get_datadir(),'des-round17-poly.txt')
        return self.draw_polygon(filename,**kwargs)

    def draw_des_sn(self,**kwargs):
        defaults = dict(facecolor='none',edgecolor='k',lw=1,zorder=10)
        setdefaults(kwargs,defaults)

        ra  = [v['ra'] for v in DES_SN.values()]
        dec = [v['dec'] for v in DES_SN.values()]
        self.tissot(ra,dec,DECAM,**kwargs)

    def draw_smash(self,**kwargs):
        """ Draw the SMASH fields. """
        defaults=dict(facecolor='none',color='k')
        setdefaults(kwargs,defaults)

        filename = os.path.join(get_datadir(),'smash_fields_final.txt')
        smash=np.genfromtxt(filename,dtype=[('ra',float),('dec',float)],usecols=[4,5])
        self.tissot(smash['ra'],smash['dec'],DECAM,**kwargs)

    def draw_decals(self,**kwargs):
        defaults=dict(color='red', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(get_datadir(),'decals-poly.txt')
        return self.draw_polygon(filename,**kwargs)

    def draw_jethwa(self,filename=None,log=True,**kwargs):
        import healpy as hp
        if not filename:
            datadir = '/home/s1/kadrlica/projects/bliss/v0/data/'
            datadir = '/Users/kadrlica/bliss/observing/data'
            filename = os.path.join(datadir,'jethwa_satellites_n256.fits.gz')
        hpxmap = hp.read_map(filename)
        if log:
            return self.draw_hpxmap(np.log10(hpxmap),**kwargs)
        else:
            return self.draw_hpxmap(hpxmap,**kwargs)

    def draw_planet9(self,**kwargs):
        from scipy.interpolate import interp1d
        from scipy.interpolate import UnivariateSpline
        defaults=dict(color='b',lw=3)
        setdefaults(kwargs,defaults)
        datadir = '/home/s1/kadrlica/projects/bliss/v0/data/'
        datadir = '/Users/kadrlica/bliss/observing/data/'
        ra_lo,dec_lo=np.genfromtxt(datadir+'p9_lo.txt',usecols=(0,1)).T
        ra_lo,dec_lo = ra_lo[::-1],dec_lo[::-1]
        ra_hi,dec_hi=np.genfromtxt(datadir+'p9_hi.txt',usecols=(0,1)).T
        ra_hi,dec_hi = ra_hi[::-1],dec_hi[::-1]

        spl_lo = UnivariateSpline(ra_lo,dec_lo)
        ra_lo_smooth = np.linspace(ra_lo[0],ra_lo[-1],360)
        dec_lo_smooth = spl_lo(ra_lo_smooth)

        spl_hi = UnivariateSpline(ra_hi,dec_hi)
        ra_hi_smooth = np.linspace(ra_hi[0],ra_hi[-1],360)
        dec_hi_smooth = spl_hi(ra_hi_smooth)

        self.plot(ra_lo_smooth,dec_lo_smooth,**kwargs)
        self.plot(ra_hi_smooth,dec_hi_smooth,**kwargs)

        orb = pd.read_csv(datadir+'P9_orbit_Cassini.csv').to_records(index=False)[::7]
        kwargs = dict(marker='o',s=40,edgecolor='none',cmap='jet_r')
        self.scatter(orb['ra'],orb['dec'],c=orb['cassini'],**kwargs)

    def draw_ligo(self,filename=None, log=True,**kwargs):
        import healpy as hp
        from astropy.io import fits as pyfits
        if not filename:
            datadir = '/home/s1/kadrlica/projects/bliss/v0/data/'
            datadir = '/Users/kadrlica/bliss/observing/data/'
            filename = datadir + 'obsbias_heatmap_semesterA.fits'
        hpxmap = pyfits.open(filename)[0].data
        if log: self.draw_hpxmap(np.log10(hpxmap))
        else:   self.draw_hpxmap(hpxmap)

    def draw_sfd(self,filename=None,**kwargs):
        import healpy as hp
        defaults = dict(rasterized=True,cmap=plt.cm.binary)
        setdefaults(kwargs,defaults)
        if not filename:
            datadir  = '/Users/kadrlica/bliss/observing/data/'
            filename = datadir+'lambda_sfd_ebv.fits'

        galhpx = hp.read_map(filename)
        celhpx = hpx_gal2cel(galhpx)
        return self.draw_hpxmap(np.log10(celhpx),**kwargs)

class SurveyMcBryde(SurveySkymap,McBrydeSkymap): pass
class SurveyOrtho(SurveySkymap,OrthoSkymap): pass

