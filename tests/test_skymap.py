#!/usr/bin/env python
"""
Unit tests for skymap
"""
__author__ = "Alex Drlica-Wagner"

import os,sys
import unittest

# Careful, have a local .matplotlibrc
import matplotlib

import pylab as plt
import numpy as np
import healpy as hp

import cartosky
from cartosky import Skymap,McBrydeSkymap,OrthoSkymap
from cartosky import SurveySkymap,SurveyMcBryde,SurveyOrtho
from cartosky import DESSkymap,BlissSkymap

try:
   import healsparse as hsp
   hsp_avail = True
except ImportError:
   hsp_avail = False

nside = 8

SKYMAPS = [Skymap,McBrydeSkymap,OrthoSkymap]
SURVEYS = [SurveySkymap,SurveyMcBryde,SurveyOrtho]
ZOOMS   = [DESSkymap,BlissSkymap]

class TestSkymap(unittest.TestCase):

    def test_skymap(self):
        fig,axes = plt.subplots(1,3,figsize=(12,3))
        for i,cls in enumerate(SKYMAPS):
            plt.sca(axes[i])
            m = cls()
            m.draw_milky_way()
            plt.title('Galactic Plane (%s)'%cls.__name__)

    def test_survey_skymap(self):
        fig,axes = plt.subplots(1,3,figsize=(12,3))
        for i,cls in enumerate(SURVEYS):
            plt.sca(axes[i])
            m = cls()
            m.draw_des()
            m.draw_maglites()
            m.draw_bliss()
            plt.title('Footprints (%s)'%cls.__name__)

    def test_zoom_skymap(self):
        for i,cls in enumerate(ZOOMS):
            plt.figure()
            m = cls()
            m.draw_des()
            m.draw_maglites()
            m.draw_bliss()
            plt.suptitle('Zoom Footprints (%s)'%cls.__name__)

    def test_draw_hpxmap(self):
        """ Test drawing a full healpix skymap """
        hpxmap = np.arange(hp.nside2npix(nside))
        fig,axes = plt.subplots(1,3,figsize=(12,3))
        for i,cls in enumerate(SKYMAPS):
            plt.sca(axes[i])
            m = cls()
            m.draw_hpxmap(hpxmap,xsize=400)
            plt.title('HEALPix Map (%s)'%cls.__name__)

    def test_draw_explicit_hpxmap(self):
        """ Test an explicit healpix map """
        pix = hpxmap = np.arange(525,535)
        fig,axes = plt.subplots(1,3,figsize=(12,3))
        for i,cls in enumerate(SKYMAPS):
            plt.sca(axes[i])
            m = cls()
            m.draw_hpxmap(hpxmap,pix,nside,xsize=400)
            plt.title('Partial HEALPix Map (%s)'%cls.__name__)

    @unittest.skipUnless(hsp_avail, "Skip this routine unless healsparse is installed")
    def test_draw_hspmap(self):
        """ Test drawing a HealSparse map """
        nside_sparse=4096
        nside_coverage=64
        hspmap = hsp.HealSparseMap.make_empty(nside_coverage, nside_sparse, dtype=np.int, sentinel=0)
        circ = hsp.Circle(ra=0, dec=0, radius=2.0, value=10)
        hsp.geom.realize_geom([circ], hspmap)
        m = cartosky.Skymap(projection='cyl')
        llcrnrlon,urcrnrlon=3.0,-3.0
        llcrnrlat,urcrnrlat=-3,3.0
        m.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat])
        m.draw_hspmap(hspmap)
        plt.title('HEALPix Zoom (nside=%i)'%nside_sparse)

    def test_draw_hpxbin(self):
        """ Test drawing hpxbin from points """
        kwargs = dict(nside=64)
        fig,axes = plt.subplots(1,3,figsize=(12,3))
        for i,cls in enumerate(SKYMAPS):
            plt.sca(axes[i])
            m = cls()
            # This is not uniform
            size = int(1e6)
            lon = np.random.uniform(0,360,size=size)
            lat = np.random.uniform(-90,90,size=size)
            m.draw_hpxbin(lon,lat,**kwargs)
            plt.title('HEALPix Bin (%s)'%cls.__name__)

    def test_draw_hires(self):
        """ Draw a partial healpix map with very large nside """
        nside = 4096*2**5
        ra,dec = 45,-45
        radius = 0.05
        pixels = cartosky.healpix.ang2disc(nside,ra,dec,radius)
        values = pixels

        plt.figure()
        m = Skymap(projection='cyl', celestial=False, gridlines=False)

        llcrnrlon,urcrnrlon=ra+2*radius,ra-2*radius,
        llcrnrlat,urcrnrlat=dec-2*radius,dec+2*radius
        m.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat])

        m.draw_hpxmap(values,pixels,nside=nside,xsize=400)
        # Draw the grid lines
        draw_labels = True
        xlocs = np.linspace(llcrnrlon,urcrnrlon,5)
        ylocs = np.linspace(llcrnrlat,urcrnrlat,5)
        grid = m.ax.gridlines(draw_labels=draw_labels,
                              xlocs=xlocs,ylocs=ylocs,
                              linestyle=':')

        # ADW: This doesn't work (removes y-labels)...
        #m.ax.invert_xaxis()

        plt.title('HEALPix Zoom (nside=%i)'%nside)


    def test_draw_focal_planes(self):
        """ Draw a DECam focal planes """
        ra,dec = 45,-45
        radius = 1.5
        delta = 1.0

        plt.figure()
        m = Skymap(projection='cyl', celestial=False, gridlines=False)
        llcrnrlon,urcrnrlon = ra+2*radius, ra-2*radius
        llcrnrlat,urcrnrlat = dec-2*radius, dec+2*radius
        m.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat])

        # Can plot individual fields
        m.draw_focal_planes([ra+delta/2],[dec-delta/2],color='g')
        # Or as arrays
        m.draw_focal_planes([ra,ra-delta,ra-delta],[dec,dec+delta,dec-delta],color='r')

        m.tissot(ra,dec,1.1,ec='g',fc='none')

        ## Draw the grid lines
        draw_labels = True
        xlocs = np.linspace(llcrnrlon,urcrnrlon,5)
        ylocs = np.linspace(llcrnrlat,urcrnrlat,5)
        m.ax.gridlines(draw_labels=draw_labels,
                       xlocs=xlocs,ylocs=ylocs,
                       linestyle=':')

        # ADW: This doesn't work (removes y-labels)...
        #m.ax.invert_xaxis()

        plt.title('DECam Focal Planes')

    def test_zoom_to_fit(self):
        nside = 64
        ra,dec = 15,-45
        radius = 10.0
        pixels = cartosky.healpix.ang2disc(nside,ra,dec,radius)
        values = pixels

        fig,axes = plt.subplots(1,3,figsize=(12,3))
        for i,cls in enumerate(SKYMAPS):
            plt.sca(axes[i])
            m = cls(gridlines=False)
            m.draw_hpxmap(values,pixels,nside=nside,xsize=200)
            m.zoom_to_fit(values,pixels,nside)
            draw_labels = i==0
            xlocs = np.linspace(ra-2*radius,ra+2*radius,5)
            ylocs = np.linspace(dec-2*radius,dec+2*radius,5)
            m.ax.gridlines(draw_labels=draw_labels,
                           xlocs=xlocs,ylocs=ylocs,
                           linestyle=':')

            #m.draw_parallels(np.linspace(dec-2*radius,dec+2*radius,5),
            #                 labelstyle='+/-',labels=[1,0,0,0])
            #m.draw_meridians(np.linspace(ra-2*radius,ra+2*radius,5),
            #                 labelstyle='+/-',labels=[0,0,0,1])
            plt.title('Zoom to Fit (%s)'%cls.__name__)


if __name__ == '__main__':
    if sys.flags.interactive: plt.ion()
    unittest.main()

