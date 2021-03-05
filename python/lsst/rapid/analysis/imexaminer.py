#
# LSST Data Management System
#
# Copyright 2008-2018  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope hat it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ['ImageExaminer']

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import scipy.ndimage as ndIm

from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.colors import LogNorm
from matplotlib.offsetbox import AnchoredText

import lsst.geom as geom
from scipy.optimize import curve_fit
from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask, QuickFrameMeasurementTaskConfig
from lsst.rapid.analysis.utils import getImageStats, argMax2d, countPixels


def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


class ImageExaminer():
    """Class for the reproducing the functionality of imexam.
    """

    def __init__(self, exp, doTweakCentroid=True, savePlots=None, centroid=None, boxHalfSize=15):

        self.exp = exp
        self.savePlots = savePlots
        self.doTweakCentroid = doTweakCentroid

        self.boxHalfSize = boxHalfSize
        if centroid is None:
            qfmTaskConfig = QuickFrameMeasurementTaskConfig()
            qfmTask = QuickFrameMeasurementTask(config=qfmTaskConfig)
            result = qfmTask.run(exp)
            if not result.success:
                msg = ("Failed to automatically find source in image. "
                       "Either provide a centroid manually or use a new image")
                raise RuntimeError(msg)
            self.centroid = result.brightestObjCentroid
        else:
            self.centroid = centroid

        self.data = self.getStarBoxData()
        if self.doTweakCentroid:
            self.tweakCentroid()
            self.data = self.getStarBoxData()

        self.xx, self.yy = self.getMeshGrid(self.data)

        self.imStats = getImageStats(self.exp)
        self.imStats.centroid = self.centroid
        self.imStats.intCentroid = self.intCoords(self.centroid)
        self.imStats.intCentroidRounded = self.intRoundCoords(self.centroid)
        self.imStats.nStatPixInBox = self.nSatPixInBox

    def intCoords(self, coords):
        return np.asarray(coords, dtype=int)

    def intRoundCoords(self, coords):
        return (int(round(coords[0])), int(round(coords[1])))

    def tweakCentroid(self):
        peak, uniquePeak, otherPeaks = argMax2d(self.data)
        # saturated stars don't tend to have ambiguous max pixels
        # due to the bunny ears left after interpolation
        nSatPix = self.nSatPixInBox

        if not uniquePeak or nSatPix:
            print('Found multiple max pixels or star is saturated, usign CoM for centroid')
            peak = ndIm.center_of_mass(self.data)

        offset = np.asarray(peak) - np.array((self.boxHalfSize, self.boxHalfSize))
        print(f"Centroid adjusted by {offset} pixels")
        x = self.centroid[0] + offset[1]  # yes, really, centroid is x,y offset is y,x
        y = self.centroid[1] + offset[0]
        self.centroid = (x, y)

    def getStats(self):
        return self.imStats

    def _calcBbox(self, centroid):
        centroidPoint = geom.Point2I(centroid)
        extent = geom.Extent2I(1, 1)
        bbox = geom.Box2I(centroidPoint, extent)
        bbox = bbox.dilatedBy(self.boxHalfSize)
        bbox = bbox.clippedTo(self.exp.getBBox())
        if bbox.getDimensions()[0] != bbox.getDimensions()[1]:
            # TODO: one day support clipped, nonsquare regions
            # but it's nontrivial due to all the plotting options
            maxsize = np.min(centroid)
            msg = (f"With centroid at {centroid} and boxHalfSize {self.boxHalfSize} "
                   "the selection runs off the edge of the chip. Boxsize has been "
                   f"automatically shrunk to {maxsize} (only square selections are "
                   "currently supported)")
            print(msg)
            self.boxHalfSize = maxsize
            return self._calcBbox(centroid)

        return bbox

    def getStarBoxData(self):
        bbox = self._calcBbox(self.centroid)
        self.starBbox = bbox  # needed elsewhere, so always set when calculated
        self.nSatPixInBox = countPixels(self.exp.maskedImage[self.starBbox], 'SAT')
        return self.exp.image[bbox].array

    def getMeshGrid(self, data):
        xlen, ylen = data.shape
        xx = np.arange(-1*xlen/2, xlen/2, 1)
        yy = np.arange(-1*ylen/2, ylen/2, 1)
        xx, yy = np.meshgrid(xx, yy)
        return xx, yy

    def radialAverage(self):
        xlen, ylen = self.data.shape
        center = np.array([xlen/2, ylen/2])
        # TODO: add option to move centroid to max pixel for radial (argmax 2d)

        distances = []
        values = []

        # could be much faster, but the array is tiny so its fine
        for i in range(xlen):
            for j in range(ylen):
                value = self.data[i, j]
                dist = norm((i, j) - center)
                if dist > xlen//2:
                    continue  # clip to box size, we don't need a factor of sqrt(2) extra
                values.append(value)
                distances.append(dist)
        return distances, values

    def plotRadialAverage(self, ax=None):
        plotDirect = False
        if not ax:
            ax = plt.subplot(111)
            plotDirect = True

        distances, values = self.radialAverage()

        peakPos = 0
        amplitude = np.max(values)
        width = 5

        try:
            pars, pCov = curve_fit(gauss, distances, values, [amplitude, peakPos, width])
            pars[0] = np.abs(pars[0])
            pars[2] = np.abs(pars[2])
        except RuntimeError:
            pars = None

        ax.plot(distances, values, 'x', label='Radial average')
        if pars is not None:
            fitline = gauss(distances, *pars)
            ax.plot(distances, fitline, label="Gaussian fit")

        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')  # equal aspect for non-images
        ax.legend()

        if plotDirect:
            plt.show()

    def plotContours(self, ax=None, nContours=10):
        plotDirect = False
        if not ax:
            fig = plt.figure(figsize=(8, 8))  # noqa F841
            ax = plt.subplot(111)
            plotDirect = True

        vmin = np.percentile(self.data, 0.1)
        vmax = np.percentile(self.data, 99.9)
        lvls = np.linspace(vmin, vmax, nContours)
        intervalSize = (lvls[1]-lvls[0])
        contourPlot = ax.contour(self.xx, self.yy, self.data, levels=lvls)  # noqa F841
        print(f"Contoured from {vmin:,.0f} to {vmax:,.0f} using {nContours} contours of {intervalSize:.1f}")

        ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=8)
        ax.set_aspect("equal")

        if plotDirect:
            plt.show()

    def plotSurface(self, ax=None, useColor=True):
        plotDirect = False
        if not ax:
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(10, 10))
            plotDirect = True

        if useColor:
            surf = ax.plot_surface(self.xx, self.yy, self.data, cmap=cm.plasma,
                                   linewidth=1, antialiased=True, color='k', alpha=0.9)
        else:
            surf = ax.plot_wireframe(self.xx, self.yy, self.data, cmap=cm.gray,  # noqa F841
                                     linewidth=1, antialiased=True, color='k')

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter('{x:,.0f}')

        if plotDirect:
            plt.show()

    def plotStar(self, ax=None, logScale=False):
        # TODO: display centroid in use
        plotDirect = False
        if not ax:
            ax = plt.subplot(111)
            plotDirect = True

        if logScale:
            ax.imshow(self.data, norm=LogNorm(), origin='lower')
        else:
            ax.imshow(self.data, origin='lower')
        ax.tick_params(which="major", direction="in", top=True, right=True, labelsize=8)

        if plotDirect:
            plt.show()

    def plotFullExp(self, ax=None):
        # TODO: display centroid in use
        plotDirect = False
        if not ax:
            ax = plt.subplot(111)
            plotDirect = True

        imData = self.exp.image.array
        vmin = np.percentile(imData, 1)
        vmax = np.percentile(imData, 99)
        ax.imshow(self.exp.image.array, norm=LogNorm(vmin=vmin, vmax=vmax),
                  origin='lower', cmap='gray',)
        ax.tick_params(which="major", direction="in", top=True, right=True, labelsize=8)

        if plotDirect:
            plt.show()

    def plotRowColSlices(self, ax=None, logScale=False):
        # TODO: display centroid in use

        # slice through self.boxHalfSize because it's always the point being
        # used by definition
        rowSlice = self.data[self.boxHalfSize, :]
        colSlice = self.data[:, self.boxHalfSize]

        plotDirect = False
        if not ax:
            ax = plt.subplot(111)
            plotDirect = True

        ax.plot(rowSlice, label='Row plot')
        ax.plot(colSlice, label='Column plot')
        if logScale:
            pass
            # TODO: set yscale as log here also protect against negatives

        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')  # equal aspect for non-images

        ax.legend()
        if plotDirect:
            plt.show()

    def plotStats(self, ax, lines):

        text = "\n".join([line for line in lines])

        stats_text = AnchoredText(text, loc="center", pad=0.5,
                                  prop=dict(size=11.5, ma="left", backgroundcolor="white",
                                            color="black", family='monospace'))
        ax.add_artist(stats_text)
        ax.axis('off')

    def plot(self):
        figsize = 6
        fig = plt.figure(figsize=(figsize*3, figsize*2))

        ax1 = fig.add_subplot(331)
        ax2 = fig.add_subplot(332)
        ax3 = fig.add_subplot(333)
        ax4 = fig.add_subplot(334, projection='3d')
        ax5 = fig.add_subplot(335)
        ax6 = fig.add_subplot(336)
        ax7 = fig.add_subplot(337)
        ax8 = fig.add_subplot(338)
        ax9 = fig.add_subplot(339)

        axExp = ax1
        axStar = ax2
        axStats1 = ax3  # noqa F841 - overwritten
        axSurf = ax4
        axCont = ax5
        axStats2 = ax6  # noqa F841 - overwritten
        axSlices = ax7
        axRadial = ax8
        axStats3 = ax9  # noqa F841 - overwritten

        self.plotFullExp(axExp)
        self.plotStar(axStar)
        self.plotSurface(axSurf)
        self.plotContours(axCont)
        self.plotRowColSlices(axSlices)
        self.plotRadialAverage(axRadial)

        # overwrite three axes with this one spanning 3 rows
        axStats = plt.subplot2grid((3, 3), (0, 2), rowspan=3)

        lines = []
        for k, v in self.imStats.getDict().items():
            if type(v) == float or isinstance(v, np.floating):
                v = f"{v:,.1f}"
            lines.append(f"{k}: {v}")
        self.plotStats(axStats, lines)
        print(lines)

        plt.tight_layout()
        plt.show()

    def plotAll(self):
        self.plotStar()
        self.plotRadialAverage()
        self.plotContours()
        self.plotSurface()
        self.plotStar()
        self.plotRowColSlices()
