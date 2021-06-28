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

import os
import lsst.daf.persistence as dafPersist
import lsst.geom as geom

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Arrow, Rectangle, Circle
import matplotlib.cm as cm
from matplotlib import gridspec

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit
from scipy.linalg import norm

# TODO: change these back to local .imports
from lsst.rapid.analysis.bestEffort import BestEffortIsr
from lsst.rapid.analysis import ImageExaminer
from lsst.rapid.analysis.utils import FWHMTOSIGMA, SIGMATOFWHM
from lsst.atmospec.utils import isDispersedExp

from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask, QuickFrameMeasurementTaskConfig


__all__ = ["SpectralFocusAnalyzer", "NonSpectralFocusAnalyzer"]


@dataclass
class FitResult:
    amp: float
    mean: float
    sigma: float


class SpectralFocusAnalyzer():

    def __init__(self, repoDir, **kwargs):

        self._butler = dafPersist.Butler(os.path.join(repoDir, "rerun/quickLook/"))
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

    def run(self, dayObs, seqNums, doDisplay=False, hideFit=False):
        self.getFocusData(dayObs, seqNums, doDisplay=doDisplay)
        bestFits = self.fitDataAndPlot(hideFit=hideFit)
        return bestFits

    def getFocusData(self, dayObs, seqNums, doDisplay=False):
        fitData = {}
        filters = set()
        objects = set()

        for seqNum in seqNums:
            fitData[seqNum] = {}
            if self._butler.datasetExists('quickLookExp', **{'dayObs': dayObs, 'seqNum': seqNum}):
                exp = self._butler.get('quickLookExp', **{'dayObs': dayObs, 'seqNum': seqNum})
            else:
                exp = self._bestEffort.getExposure({'dayObs': dayObs, 'seqNum': seqNum})

            # sanity checking
            filt = exp.getFilter().getName()
            obj = self._butler.queryMetadata('raw', 'object', dayObs=dayObs, seqNum=seqNum)[0]
            objects.add(obj)
            filters.add(filt)
            assert isDispersedExp(exp), f"Image is not dispersed! (filter = {filt})"
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

                try:
                    coeffs, var_matrix = curve_fit(self.gauss, xs, data1d, p0=p0)
                except RuntimeError:
                    coeffs = (np.nan, np.nan, np.nan)

                fitData[seqNum][i] = FitResult(amp=abs(coeffs[0]), mean=coeffs[1], sigma=abs(coeffs[2]))
                if doDisplay:
                    axes[1].plot(xs, data1d, 'x', c=self.COLORS[i])
                    highResX = np.linspace(0, len(data1d), 1000)
                    if coeffs[0] is not np.nan:
                        axes[1].plot(highResX, self.gauss(highResX, *coeffs), 'k-')

            if doDisplay:  # show all color boxes together
                plt.title(f'Fits to seqNum {seqNum}')
                plt.show()

            focuserPosition = self._getFocusFromHeader(exp)
            fitData[seqNum]['focus'] = focuserPosition

        self.fitData = fitData
        self.filter = filters.pop()
        self.object = objects.pop()

        return

    def fitDataAndPlot(self, hideFit=False, hexapodZeroPoint=0):
        data = self.fitData
        filt = self.filter
        obj = self.object

        bestFits = []

        titleFontSize = 18
        legendFontSize = 12
        labelFontSize = 14

        arcminToPixel = 10
        sigmaToFwhm = 2.355

        f, axes = plt.subplots(2, 1, figsize=[10, 12])
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


class NonSpectralFocusAnalyzer():

    def __init__(self, repoDir, **kwargs):
        self._butler = dafPersist.Butler(repoDir)
        self._bestEffort = BestEffortIsr(repoDir, **kwargs)

    @staticmethod
    def _getFocusFromHeader(exp):
        return float(exp.getMetadata()["FOCUSZ"])

    @staticmethod
    def gauss(x, *pars):
        amp, mean, sigma = pars
        return amp*np.exp(-(x-mean)**2/(2.*sigma**2))

    def run(self, dayObs, seqNums, *, manualCentroid=None, doDisplay=False, hideFit=False, doForceCoM=False):
        self.getFocusData(dayObs, seqNums, manualCentroid=manualCentroid, doDisplay=doDisplay,
                          doForceCoM=doForceCoM)
        bestFit = self.fitDataAndPlot(hideFit=hideFit)
        return bestFit

    def getFocusData(self, dayObs, seqNums, *, manualCentroid=None, doCheckDispersed=True,
                     doDisplay=False, doForceCoM=False):
        fitData = {}
        filters = set()
        objects = set()

        maxDistance = 200
        firstCentroid = None

        for seqNum in seqNums:
            fitData[seqNum] = {}
            butler = self._butler

            if butler.datasetExists('quickLookExp', {'dayObs': dayObs, 'seqNum': seqNum}):
                exp = butler.get('quickLookExp', {'dayObs': dayObs, 'seqNum': seqNum})
            else:
                print(f"quickLookExp not found for {seqNum}, reproducing...")
                exp = self._bestEffort.getExposure({'dayObs': dayObs, 'seqNum': seqNum})

            # sanity/consistency checking
            filt = exp.getFilterLabel().physicalLabel
            obj = butler.queryMetadata('raw', 'object', dayObs=dayObs, seqNum=seqNum)[0]
            objects.add(obj)
            filters.add(filt)
            if doCheckDispersed:
                assert not isDispersedExp(exp), f"Image is dispersed! (filter = {filt})"
            assert len(filters) == 1, "You accidentally mixed filters!"
            assert len(objects) == 1, "You accidentally mixed objects!"

            imExam = ImageExaminer(exp, centroid=manualCentroid, doTweakCentroid=True, boxHalfSize=105,
                                   doForceCoM=doForceCoM)
            if doDisplay:
                imExam.plot()

            fwhm = imExam.imStats.fitFwhm
            amp = imExam.imStats.fitAmp
            gausMean = imExam.imStats.fitGausMean
            centroid = imExam.centroid

            if seqNum == seqNums[0]:
                firstCentroid = centroid

            dist = norm(np.array(centroid) - np.array(firstCentroid))
            if dist > maxDistance:
                print(f"Skipping {seqNum} because distance {dist}> maxDistance {maxDistance}")

            fitData[seqNum]['fitResult'] = FitResult(amp=amp, mean=gausMean, sigma=fwhm*FWHMTOSIGMA)
            fitData[seqNum]['eeRadius50'] = imExam.imStats.eeRadius50
            fitData[seqNum]['eeRadius80'] = imExam.imStats.eeRadius80
            fitData[seqNum]['eeRadius90'] = imExam.imStats.eeRadius90

            focuserPosition = self._getFocusFromHeader(exp)
            fitData[seqNum]['focus'] = focuserPosition

        self.fitData = fitData
        self.filter = filters.pop()
        self.object = objects.pop()

        return

    def fitDataAndPlot(self, hideFit=False):

        fitData = self.fitData

        labelFontSize = 14

        arcminToPixel = 10

        fig = plt.figure(figsize=(10, 10))  # noqa
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])

        seqNums = sorted(fitData.keys())
        widths = [fitData[seqNum]['fitResult'].sigma * SIGMATOFWHM / arcminToPixel for seqNum in seqNums]
        focusPositions = [fitData[seqNum]['focus'] for seqNum in seqNums]
        fineXs = np.linspace(np.min(focusPositions), np.max(focusPositions), 101)

        fwhmFitPars = np.polyfit(focusPositions, widths, 2)
        fwhmFitMin = -fwhmFitPars[1] / (2.0*fwhmFitPars[0])

        ax0 = plt.subplot(gs[0])
        ax0.scatter(focusPositions, widths, c='k')
        ax0.set_ylabel('FWHM (arcsec)', fontsize=labelFontSize)
        ax0.plot(fineXs, np.poly1d(fwhmFitPars)(fineXs), 'b-')
        ax0.axvline(fwhmFitMin, c='r', ls='--')

        ee90s = [fitData[seqNum]['eeRadius90'] for seqNum in seqNums]
        ee80s = [fitData[seqNum]['eeRadius80'] for seqNum in seqNums]
        ee50s = [fitData[seqNum]['eeRadius50'] for seqNum in seqNums]
        ax1 = plt.subplot(gs[1], sharex=ax0)
        ax1.scatter(focusPositions, ee90s, c='r', label='Encircled energy 90%')
        ax1.scatter(focusPositions, ee80s, c='g', label='Encircled energy 80%')
        ax1.scatter(focusPositions, ee50s, c='b', label='Encircled energy 50%')

        ee90FitPars = np.polyfit(focusPositions, ee90s, 2)
        ee90FitMin = -ee90FitPars[1] / (2.0*ee90FitPars[0])
        ee80FitPars = np.polyfit(focusPositions, ee80s, 2)
        ee80FitMin = -ee80FitPars[1] / (2.0*ee80FitPars[0])
        ee50FitPars = np.polyfit(focusPositions, ee50s, 2)
        ee50FitMin = -ee50FitPars[1] / (2.0*ee50FitPars[0])

        ax1.plot(fineXs, np.poly1d(ee90FitPars)(fineXs), 'r-')
        ax1.plot(fineXs, np.poly1d(ee80FitPars)(fineXs), 'g-')
        ax1.plot(fineXs, np.poly1d(ee50FitPars)(fineXs), 'b-')

        ax1.axvline(ee90FitMin, c='r', ls='--')
        ax1.axvline(ee80FitMin, c='g', ls='--')
        ax1.axvline(ee50FitMin, c='b', ls='--')

        ax1.set_xlabel('User-applied focus offset (mm)', fontsize=labelFontSize)
        ax1.set_ylabel('Radius (pixels)', fontsize=labelFontSize)

        ax1.legend()

        plt.subplots_adjust(hspace=.0)
        plt.show()

        results = {"fwhmFitMin": fwhmFitMin,
                   "ee90FitMin": ee90FitMin,
                   "ee80FitMin": ee80FitMin,
                   "ee50FitMin": ee50FitMin}

        return results

# if __name__ == '__main__':
#     repoDir = '/project/shared/auxTel/'
#     analyzer = SpectralFocusAnalyzer(repoDir)
#     # dataId = {'dayObs': '2020-02-20', 'seqNum': 485}  # direct image
#     dataId = {'dayObs': '2020-03-12'}
#     seqNums = [121, 122]
#     data, filt, obj = analyzer.getFocusData(dataId['dayObs'],
#                                             seqNums, doDisplay=True)
#     analyzer.fitDataAndPlot(data, filt, obj)
