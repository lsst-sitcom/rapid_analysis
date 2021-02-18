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

__all__ = ['SummarizeImage', 'SummarizeImageTaskConfig']

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.optimize import curve_fit

import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.atmospec.processStar import ProcessStarTask
from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask, QuickFrameMeasurementTaskConfig


# class SummarizeImageTaskConfig(pexConfig.Config):

#     def setDefaults(self):
#         return


class SummarizeImage():
    """Task for the spectral extraction of single-star dispersed images.

    For a full description of how this tasks works, see the run() method.
    """

    # ConfigClass = SummarizeImageTaskConfig
    # _DefaultName = "summarizeImage"

    def __init__(self, exp, display=None, debug=False, **kwargs):
        # TODO: rename psfRefObjLoader to refObjLoader
        super().__init__(**kwargs)
        self.exp = exp
        self.display = display
        self.debug = debug

        qfmTaskConfig = QuickFrameMeasurementTaskConfig()
        self.qfmTask = QuickFrameMeasurementTask(config=qfmTaskConfig)

        pstConfig = ProcessStarTask.ConfigClass()
        pstConfig.offsetFromMainStar = 400
        self.processStarTask = ProcessStarTask(config=pstConfig)

        self.init()

    # def runDataRef(self, dataRef):
    #     """Run the ProcessStarTask on a ButlerDataRef for a single exposure.

    #     Runs isr to get the postISR exposure from the dataRef and passes this
    #     to the run() method.

    #     Parameters
    #     ----------
    #     dataRef : `daf.persistence.butlerSubset.ButlerDataRef`
    #         Butler reference of the detector and exposure ID
    #     """
    #     butler = dataRef.getButler()
    #     dataId = dataRef.dataId
    #     self.log.info("Processing %s" % (dataRef.dataId))

    @staticmethod
    def bboxToAwfDisplayLines(box):
        """Takes a bbox, returns a list of lines such that they can be plotted with

        for line in lines:
            display.line(line, ctype='red')

        """
        x0 = box.beginX
        x1 = box.endX
        y0 = box.beginY
        y1 = box.endY
        return [[(x0, y0), (x1, y0)], [(x0, y0), (x0, y1)], [(x1, y0), (x1, y1)], [(x0, y1), (x1, y1)]]

    def eraseDisplay(self):
        if self.display:
            self.display.erase()

    def displaySpectrumBbox(self):
        if self.display:
            lines = self.bboxToAwfDisplayLines(self.spectrumbbox)
            for line in lines:
                self.display.line(line, ctype='red')
        else:
            print("No display set")

    def displayStarLocation(self):
        if self.display:
            self.display.dot('x', *self.qfmResult.brightestObjCentroid, size=50)
            self.display.dot('o', *self.qfmResult.brightestObjCentroid, size=50)
        else:
            print("No display set")

    def calcGoodSpectrumSection(self, threshold=5, windowSize=5):
        length = len(self.ridgeLineLocations)
        chunks = length // windowSize
        stddevs = []
        for i in range(chunks+1):
            stddevs.append(np.std(self.ridgeLineLocations[i*windowSize:(i+1)*windowSize]))

        goodPoints = np.where(np.asarray(stddevs) < threshold)[0]
        minPoint = (goodPoints[2] - 2) * windowSize
        maxPoint = (goodPoints[-3] + 3) * windowSize
        minPoint = max(minPoint, 0)
        maxPoint = min(maxPoint, length)
        if self.debug:
            plt.plot(range(0, length+1, windowSize), stddevs)
            plt.hlines(threshold, 0, length, colors='r', ls='dashed')
            plt.vlines(minPoint, 0, max(stddevs)+10, colors='k', ls='dashed')
            plt.vlines(maxPoint, 0, max(stddevs)+10, colors='k', ls='dashed')
            plt.title(f'Ridgeline scatter, windowSize={windowSize}')

        return (minPoint, maxPoint)

    def fit(self):
        def gauss(x, a, x0, sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        data = self.spectrumData[self.goodSlice]
        nRows, nCols = data.shape
        # don't subtract the row median or even a percentile - seems bad
        # fitting a const also seems bad - needs some better thought

        parameters = np.zeros((nRows, 3))
        pCovs = []
        xs = np.arange(nCols)
        for rowNum, row in enumerate(data):
            peakPos = self.ridgeLineLocations[rowNum]
            amplitude = row[peakPos]
            width = 7
            try:
                pars, pCov = curve_fit(gauss, xs, row, [amplitude, peakPos, width], maxfev=100)
                pCovs.append(pCov)
            except RuntimeError:
                pars = [np.nan] * 3
            if not np.all([p < 1e7 for p in pars]):
                pars = [np.nan] * 3
            parameters[rowNum] = pars

        parameters[:, 0] = np.abs(parameters[:, 0])
        parameters[:, 2] = np.abs(parameters[:, 2])
        self.parameters = parameters

    def plot(self):
        fig = plt.figure(figsize=(10, 10))

        # spectrum
        ax0 = plt.subplot2grid((4, 4), (0, 0), colspan=3)
        d = self.spectrumData[self.goodSlice].T
        vmin = np.percentile(d, 1)
        vmax = np.percentile(d, 99)
        pos = ax0.imshow(self.spectrumData[self.goodSlice].T, vmin=vmin, vmax=vmax, origin='lower')
        # fig.colorbar(pos, ax=ax0, orientation="horizontal")

        # spectrum histogram
        axHist = plt.subplot2grid((4, 4), (0, 3))
        data = self.spectrumData
        histMax = np.nanpercentile(data, 99.99)
        histMin = np.nanpercentile(data, 0.001)
        axHist.hist(data[(data >= histMin) & (data <= histMax)].flatten(), bins=100)
        underflow = len(data[data < histMin])
        overflow = len(data[data > histMax])
        axHist.set_yscale('log', nonpositive='clip')
        text = f"Underflow = {underflow}"
        text += f"\nOverflow = {overflow}"
        anchored_text = AnchoredText(text, loc=1, pad=1)
        axHist.add_artist(anchored_text)

        # peak fluxes
        ax1 = plt.subplot2grid((4, 4), (1, 0), colspan=3)
        ax1.plot(self.ridgeLineValues[self.goodSlice], label='Raw peak value')
        ax1.plot(self.parameters[:, 0], label='Fitted amplitude')
#         continuumFlux = np.nanpercentile(self.parameters[:, 0], 99.5)
        ax1.axhline(self.continuumFlux98, ls='dashed', color='g')
#         ymax = np.max(np.nanpercentile(self.parameters[:, 0], 90)*1.5)
#         ax1.set_ylim(0, )
        ax1.set_ylabel('Peak amplitude (ADU)')
        ax1.set_xlabel('Spectrum position (pixels)')
        ax1.legend(title=f"Continuum flux = {self.continuumFlux98:.0f} ADU")
        ax1.set_title('Ridgeline plot')

        # FWHM
        ax2 = plt.subplot2grid((4, 4), (2, 0), colspan=3)
        ax2.plot(self.parameters[:, 2]*2.355)
        fwhmValues = self.parameters[:, 2]*2.355
        median_fwhm = np.nanmedian(fwhmValues)
        ax2.axhline(median_fwhm, ls='dashed', color='k')
        ymin = max(np.nanmin(fwhmValues)-3, 0)
        ax2.set_ylim(ymin, np.nanpercentile(fwhmValues, 90)*1.5)
        ax2.set_ylabel('FWHM (pixels)')
        ax2.set_xlabel('Spectrum position (pixels)')
        ax2.legend(title=f"Median FWHM = {median_fwhm:.1f} pix")
        ax2.set_title('Spectrum FWHM')

        # row fluxes
        ax3 = plt.subplot2grid((4, 4), (3, 0), colspan=3)
        ax3.plot(self.rowSums[self.goodSlice])
        ax3.set_ylabel('Total row flux (ADU)')
        ax3.set_xlabel('Spectrum position (pixels)')
        ax3.legend()
        ax3.set_title('Row sums')

        # textbox
        ax4 = plt.subplot2grid((4, 4), (1, 3), rowspan=3)
        text = "\n".join([line for line in self.text])
        stats_text = AnchoredText(text, loc=1, pad=1)
        ax4.add_artist(stats_text)

        plt.tight_layout()
        plt.show()

    def init(self):
        pass

    def generateStatsText(self, doPrint=True):
        x, y = self.qfmResult.brightestObjCentroid
        exptime = self.exp.getInfo().getVisitInfo().getExposureTime()

        lines = []

        lines.append("-----Star stats-----")
        lines.append(f"Star centroid @ ({x:.1f}, {y:.1f})")
        lines.append(f"Star max pixel flux = {self.starPeakFlux:,.0f} ADU")
        lines.append(f"Star Ap25 flux = {self.qfmResult.brightestObjApFlux25:,.0f} ADU")
        lines.append("")

        lines.append("-----Image stats-----")
        imageMedian = np.median(self.exp.image.array)
        lines.append(f"Image median = {imageMedian:.2f} ADU")
        lines.append(f"Exposure time = {exptime:.2f} s")
        lines.append("")

        lines.append("-----Rate stats-----")
        lines.append(f"Star max pixel = {self.starPeakFlux/exptime:,.0f} ADU/s")
        lines.append(f"Spectrum contiuum = {self.continuumFlux98/exptime:,.1f} ADU/s")
        lines.append("")

        if self.debug:
            lines.append("-----Debug-----")
            lines.append(f"spectrum bbox: {self.spectrumbbox}")
            lines.append(f"Good range = {self.goodSpectrumMinY},{self.goodSpectrumMaxY}")

        if doPrint:
            for line in lines:
                print(line)
        return lines

    def run(self):
        self.qfmResult = self.qfmTask.run(self.exp)
        self.intCentroidX = int(np.round(self.qfmResult.brightestObjCentroid)[0])
        self.intCentroidY = int(np.round(self.qfmResult.brightestObjCentroid)[1])
        self.starPeakFlux = self.exp.image.array[self.intCentroidY, self.intCentroidX]

        self.spectrumbbox = self.processStarTask.calcSpectrumBBox(self.exp,
                                                                  self.qfmResult.brightestObjCentroid,
                                                                  200)
        self.spectrumData = self.exp.image[self.spectrumbbox].array

        self.ridgeLineLocations = np.argmax(self.spectrumData, axis=1)
        self.ridgeLineValues = self.spectrumData[range(self.spectrumbbox.getHeight()),
                                                 self.ridgeLineLocations]
        self.rowSums = np.sum(self.spectrumData, axis=1)

        coords = self.calcGoodSpectrumSection()
        self.goodSpectrumMinY = coords[0]
        self.goodSpectrumMaxY = coords[1]
        self.goodSlice = slice(coords[0], coords[1])

        self.continuumFlux90 = np.percentile(self.ridgeLineValues, 90)  # for emission stars
        self.continuumFlux98 = np.percentile(self.ridgeLineValues, 98)  # for most stars

        self.fit()
        self.text = self.generateStatsText()
        self.plot()

        return
