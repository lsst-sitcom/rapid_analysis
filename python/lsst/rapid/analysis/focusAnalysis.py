# This file is part of rapid_analysis.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import lsst.daf.persistence as dafPersist
import lsst.geom as geom

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Arrow, Rectangle, Circle
import matplotlib.cm as cm

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

# TODO: change these back to local .imports
from lsst.rapid.analysis.bestEffort import BestEffortIsr
from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask, QuickFrameMeasurementTaskConfig


@dataclass
class FitResult:
    amp: float
    mean: float
    sigma: float


class FocusAnalyzer():

    def __init__(self, repoDir, **kwargs):

        self._butler = dafPersist.Butler(repoDir)
        self._bestEffort = BestEffortIsr(repoDir, **kwargs)
        qfmTaskConfig = QuickFrameMeasurementTaskConfig()
        self._quickMeasure = QuickFrameMeasurementTask(config=qfmTaskConfig)

        self.spectrumHalfWidth = 100
        self.spectrumBoxLength = 20
        self._spectrumBoxOffsets = [882, 1170, 1467]
        self._setColors(len(self._spectrumBoxOffsets))

    def setSpectrumBoxOffsets(self, offsets):
        self._spectrumBoxOffsets = offsets
        self._setColors(len(offsets))

    def getSpectrumBoxOffsets(self):
        return self._spectrumBoxOffsets

    def _setColors(self, nPoints):
        self.COLORS = cm.rainbow(np.linspace(0, 1, nPoints))

    @staticmethod
    def _checkImageIsDispersed(filterFullName):
        if "~" not in filterFullName:
            raise RuntimeError("Error parsing filter name {filterFullName}")
        filt, grating = filterFullName.split('~')
        if grating.startswith('EMPTY'):
            return False
        return True

    @staticmethod
    def _getFocusFromHeader(exp):
        return float(exp.getMetadata()["FOCUSZ"])

    def _getBboxes(self, centroid):
        x, y = centroid
        bboxes = []

        for offset in self._spectrumBoxOffsets:
            bbox = geom.Box2I(geom.Point2I(x-self.spectrumHalfWidth, y+offset),
                              geom.Point2I(x+self.spectrumHalfWidth, y+offset+self.spectrumBoxLength))
            bboxes.append(bbox)
        return bboxes

    def _bboxToMplRectangle(self, bbox, colorNum):
        xmin = bbox.getBeginX()
        ymin = bbox.getBeginY()
        xsize = bbox.getWidth()
        ysize = bbox.getHeight()
        rectangle = Rectangle((xmin, ymin), xsize, ysize, alpha=1, facecolor='none', lw=2,
                              edgecolor=self.COLORS[colorNum])
        return rectangle

    @staticmethod
    def gauss(x, *pars):
        amp, mean, sigma = pars
        return amp*np.exp(-(x-mean)**2/(2.*sigma**2))

    def run(self, dayObs, seqNums, display=False, hideFit=False):
        data, obj, filt = self.getFocusData(dayObs, seqNums, display=display)
        bestFits = self.fitDataAndPlot(data, obj, filt, hideFit=hideFit)
        return bestFits

    def getFocusData(self, dayObs, seqNums, doDisplay=False, display=None):
        fitData = {}
        filters = set()
        objects = set()

        for seqNum in seqNums:
            fitData[seqNum] = {}
            exp = self._bestEffort.getExposure({'dayObs': dayObs, 'seqNum': seqNum})

            # sanity checking
            filt = exp.getFilter().getName()
            obj = self._butler.queryMetadata('raw', 'object', dayObs=dayObs, seqNum=seqNum)[0]
            objects.add(obj)
            filters.add(filt)
            assert self._checkImageIsDispersed(filt), f"Image is not dispersed! (filter = {filt})"
            assert len(filters) == 1, "You accidentally mixed filters!"
            assert len(objects) == 1, "You accidentally mixed objects!"

            quickMeasResult = self._quickMeasure.run(exp)
            centroid = quickMeasResult.brightestObjCentroid
            spectrumSliceBboxes = self._getBboxes(centroid)  # inside the loop due to centroid shifts

            if doDisplay:
                fig, axes = plt.subplots(1, 2, figsize=(18, 9))
                exp.image.array[exp.image.array <= 0] = 0.001
                axes[0].imshow(exp.image.array, norm=LogNorm(), origin='lower', cmap='gray_r')
                plt.tight_layout()
                arrowy, arrowx = centroid[0] - 400, centroid[1]  # numpy is backwards
                dx, dy = 0, 300
                arrow = Arrow(arrowy, arrowx, dy, dx, width=200., color='red')
                circle = Circle(centroid, radius=25, facecolor='none', color='red')
                axes[0].add_patch(arrow)
                axes[0].add_patch(circle)
                for i, bbox in enumerate(spectrumSliceBboxes):
                    rect = self._bboxToMplRectangle(bbox, i)
                    axes[0].add_patch(rect)

            for i, bbox in enumerate(spectrumSliceBboxes):
                data1d = np.mean(exp[bbox].image.array, axis=0)  # flatten
                data1d -= np.median(data1d)
                xs = np.arange(len(data1d))

                # get rough estimates for fit
                # can't use sigma from quickMeasResult due to SDSS shape
                # failing on saturated starts, and fp.getShape() is weird
                amp = np.max(data1d)
                mean = np.argmax(data1d)
                sigma = 20
                p0 = amp, mean, sigma

                coeffs, var_matrix = curve_fit(self.gauss, xs, data1d, p0=p0)
                fitData[seqNum][i] = FitResult(amp=abs(coeffs[0]), mean=coeffs[1], sigma=abs(coeffs[2]))
                if doDisplay:
                    axes[1].plot(xs, data1d, 'x', c=self.COLORS[i])
                    highResX = np.linspace(0, len(data1d), 1000)
                    axes[1].plot(highResX, self.gauss(highResX, *coeffs), 'k-')

            if doDisplay:  # show all color boxes together
                plt.title(f'Fits to seqNum {seqNum}')
                plt.show()

            focuserPosition = self._getFocusFromHeader(exp)
            fitData[seqNum]['focus'] = focuserPosition

        return fitData, filters.pop(), objects.pop()

    def fitDataAndPlot(self, data, obj, filt, hideFit=False, hexapodZeroPoint=0):
        bestFits = []

        titleFontSize = 18
        legendFontSize = 12
        labelFontSize = 14

        arcminToPixel = 10
        sigmaToFwhm = 2.355

        f, axes = plt.subplots(2, 1, figsize=[15, 8])
        focusPositions = [data[k]['focus']-hexapodZeroPoint for k in sorted(data.keys())]
        fineXs = np.linspace(np.min(focusPositions), np.max(focusPositions), 101)
        seqNums = sorted(data.keys())

        nSpectrumSlices = len(data[list(data.keys())[0]])-1
        pointsForLegend = [0.0 for offset in range(nSpectrumSlices)]
        for spectrumSlice in range(nSpectrumSlices):  # the blue/green/red slices through the spectrum
            # for scatter plots, the color needs to be a single-row 2d array
            thisColor = np.array([self.COLORS[spectrumSlice]])

            amps = [data[seqNum][spectrumSlice].amp for seqNum in seqNums]
            widths = [data[seqNum][spectrumSlice].sigma / arcminToPixel * sigmaToFwhm for seqNum in seqNums]

            pointsForLegend[spectrumSlice] = axes[0].scatter(focusPositions, amps, c=thisColor)
            axes[0].set_xlabel('Focus position (mm)', fontsize=labelFontSize)
            axes[0].set_ylabel('Height (ADU)', fontsize=labelFontSize)

            axes[1].scatter(focusPositions, widths, c=thisColor)
            axes[1].set_xlabel('Focus position (mm)', fontsize=labelFontSize)
            axes[1].set_ylabel('FWHM (arcsec)', fontsize=labelFontSize)

            quadFitPars = np.polyfit(focusPositions, widths, 2)
            if not hideFit:
                axes[1].plot(fineXs, np.poly1d(quadFitPars)(fineXs), c=self.COLORS[spectrumSlice])
                fitMin = -quadFitPars[1] / (2.0*quadFitPars[0])
                bestFits.append(fitMin)
                axes[1].axvline(fitMin, color=self.COLORS[spectrumSlice])
                msg = f"Best focus offset = {np.round(fitMin, 2)}"
                axes[1].text(fitMin, np.mean(widths), msg, horizontalalignment='right',
                             verticalalignment='center', rotation=90, color=self.COLORS[spectrumSlice],
                             fontsize=legendFontSize)

        titleText = f"Focus curve for {obj} w/ {filt}"
        plt.suptitle(titleText, fontsize=titleFontSize)
        legendText = self._generateLegendText(nSpectrumSlices)
        axes[0].legend(pointsForLegend, legendText, fontsize=legendFontSize)
        axes[1].legend(pointsForLegend, legendText, fontsize=legendFontSize)
        f.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()

        for i, bestFit in enumerate(bestFits):
            print(f"Best fit for spectrum slice {i} = {bestFit:.4f}mm")
        return bestFits

    def _generateLegendText(self, nSpectrumSlices):
        if nSpectrumSlices == 1:
            return ['m=+1 spectrum slice']
        if nSpectrumSlices == 2:
            return ['m=+1 blue end', 'm=+1 red end']

        legendText = []
        legendText.append('m=+1 blue end')
        for i in range(nSpectrumSlices-2):
            legendText.append('m=+1 redder...')
        legendText.append('m=+1 red end')
        return legendText


if __name__ == '__main__':
    repoDir = '/project/shared/auxTel/'
    analyzer = FocusAnalyzer(repoDir)
    # dataId = {'dayObs': '2020-02-20', 'seqNum': 485}  # direct image
    dataId = {'dayObs': '2020-03-12'}
    seqNums = [121, 122]
    data, filt, obj = analyzer.getFocusData(dataId['dayObs'], seqNums, doDisplay=True)
    analyzer.fitDataAndPlot(data, filt, obj)
