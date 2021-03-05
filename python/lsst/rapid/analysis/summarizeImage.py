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

__all__ = ['SummarizeImage']

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.optimize import curve_fit

from lsst.atmospec.processStar import ProcessStarTask
from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask, QuickFrameMeasurementTaskConfig

from lsst.obs.lsst.translators.lsst import FILTER_DELIMITER
from lsst.rapid.analysis.utils import getImageStats
from astro_metadata_translator import ObservationInfo


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

        self.imStats = getImageStats(exp)

        self.init()

    @staticmethod
    def bboxToAwfDisplayLines(box):
        """Takes a bbox, returns a list of lines such that they can be plotted:

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

    def plot(self, saveAs=None):
        fig = plt.figure(figsize=(10, 10))

        # spectrum
        ax0 = plt.subplot2grid((4, 4), (0, 0), colspan=3)
        ax0.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)
        d = self.spectrumData[self.goodSlice].T
        vmin = np.percentile(d, 1)
        vmax = np.percentile(d, 99)
        pos = ax0.imshow(self.spectrumData[self.goodSlice].T, vmin=vmin, vmax=vmax, origin='lower')
        div = make_axes_locatable(ax0)
        cax = div.append_axes("bottom", size="7%", pad="8%")
        fig.colorbar(pos, cax=cax, orientation="horizontal", label="Counts")

        # spectrum histogram
        axHist = plt.subplot2grid((4, 4), (0, 3))
        data = self.spectrumData
        histMax = np.nanpercentile(data, 99.99)
        histMin = np.nanpercentile(data, 0.001)
        axHist.hist(data[(data >= histMin) & (data <= histMax)].flatten(), bins=100)
        underflow = len(data[data < histMin])
        overflow = len(data[data > histMax])
        axHist.set_yscale('log', nonpositive='clip')
        axHist.set_title('Spectrum pixel histogram')
        text = f"Underflow = {underflow}"
        text += f"\nOverflow = {overflow}"
        anchored_text = AnchoredText(text, loc=1, pad=0.5)
        axHist.add_artist(anchored_text)

        # peak fluxes
        ax1 = plt.subplot2grid((4, 4), (1, 0), colspan=3)
        ax1.plot(self.ridgeLineValues[self.goodSlice], label='Raw peak value')
        ax1.plot(self.parameters[:, 0], label='Fitted amplitude')
        ax1.axhline(self.continuumFlux98, ls='dashed', color='g')
        ax1.set_ylabel('Peak amplitude (ADU)')
        ax1.set_xlabel('Spectrum position (pixels)')
        ax1.legend(title=f"Continuum flux = {self.continuumFlux98:.0f} ADU",
                   loc="center right", framealpha=0.2, facecolor="black")
        ax1.set_title('Ridgeline plot')

        # FWHM
        ax2 = plt.subplot2grid((4, 4), (2, 0), colspan=3)
        ax2.plot(self.parameters[:, 2]*2.355, label="FWHM (pix)")
        fwhmValues = self.parameters[:, 2]*2.355
        median_fwhm = np.nanmedian(fwhmValues)
        ax2.axhline(median_fwhm, ls='dashed', color='k')
        ymin = max(np.nanmin(fwhmValues)-5, 0)
        ymax = median_fwhm*1.5
        ax2.set_ylim(ymin, ymax)
        ax2.set_ylabel('FWHM (pixels)')
        ax2.set_xlabel('Spectrum position (pixels)')
        ax2.legend(title=f"Median FWHM = {median_fwhm:.1f} pix", loc="upper right",
                   framealpha=0.2, facecolor="black")
        ax2.set_title('Spectrum FWHM')

        # row fluxes
        ax3 = plt.subplot2grid((4, 4), (3, 0), colspan=3)
        ax3.plot(self.rowSums[self.goodSlice], label="Sum across row")
        ax3.set_ylabel('Total row flux (ADU)')
        ax3.set_xlabel('Spectrum position (pixels)')
        ax3.legend(framealpha=0.2, facecolor="black")
        ax3.set_title('Row sums')

        # textbox top
#         ax4 = plt.subplot2grid((4, 4), (1, 3))
        ax4 = plt.subplot2grid((4, 4), (1, 3), rowspan=2)
        text = "short text"
        text = self.generateStatsTextboxContent(0)
        text += self.generateStatsTextboxContent(1)
        text += self.generateStatsTextboxContent(2)
        text += self.generateStatsTextboxContent(3)
        stats_text = AnchoredText(text, loc="center", pad=0.5,
                                  prop=dict(size=10.5, ma="left", backgroundcolor="white",
                                            color="black", family='monospace'))
        ax4.add_artist(stats_text)
        ax4.axis('off')

        # textbox middle
        if self.debug:
            ax5 = plt.subplot2grid((4, 4), (2, 3))
            text = self.generateStatsTextboxContent(-1)
            stats_text = AnchoredText(text, loc="center", pad=0.5,
                                      prop=dict(size=10.5, ma="left", backgroundcolor="white",
                                                color="black", family='monospace'))
            ax5.add_artist(stats_text)
            ax5.axis('off')

        plt.tight_layout()
        plt.show()

        if saveAs:
            fig.savefig(saveAs)

    def init(self):
        pass

    def generateStatsTextboxContent(self, section, doPrint=True):
        x, y = self.qfmResult.brightestObjCentroid
        exptime = self.exp.getInfo().getVisitInfo().getExposureTime()

        info = self.exp.getInfo()
        vi = info.getVisitInfo()

        fullFilterString = info.getFilterLabel().physicalLabel
        filt = fullFilterString.split(FILTER_DELIMITER)[0]
        grating = fullFilterString.split(FILTER_DELIMITER)[1]

        airmass = vi.getBoresightAirmass()
        rotangle = vi.getBoresightRotAngle().asDegrees()

        azAlt = vi.getBoresightAzAlt()
        az = azAlt[0].asDegrees()
        el = azAlt[1].asDegrees()

        md = self.exp.getMetadata()
        obsInfo = ObservationInfo(md, subset={'object'})
        obj = obsInfo.object

        lines = []

        if section == 0:
            lines.append("----- Star stats -----")
            lines.append(f"Star centroid @  {x:.0f}, {y:.0f}")
            lines.append(f"Star max pixel = {self.starPeakFlux:,.0f} ADU")
            lines.append(f"Star Ap25 flux = {self.qfmResult.brightestObjApFlux25:,.0f} ADU")
            lines.extend(["", ""])  # section break
            return '\n'.join([line for line in lines])

        if section == 1:
            lines.append("------ Image stats ---------")
            imageMedian = np.median(self.exp.image.array)
            lines.append(f"Image median   = {imageMedian:.2f} ADU")
            lines.append(f"Exposure time  = {exptime:.2f} s")
            lines.extend(["", ""])  # section break
            return '\n'.join([line for line in lines])

        if section == 2:
            lines.append("------- Rate stats ---------")
            lines.append(f"Star max pixel    = {self.starPeakFlux/exptime:,.0f} ADU/s")
            lines.append(f"Spectrum contiuum = {self.continuumFlux98/exptime:,.1f} ADU/s")
            lines.extend(["", ""])  # section break
            return '\n'.join([line for line in lines])

        if section == 3:
            lines.append("----- Observation info -----")
            lines.append(f"object  = {obj}")
            lines.append(f"filter  = {filt}")
            lines.append(f"grating = {grating}")
            lines.append(f"rotpa   = {rotangle:.1f}")

            lines.append(f"az      = {az:.1f}")
            lines.append(f"el      = {el:.1f}")
            lines.append(f"airmass = {airmass:.3f}")
            return '\n'.join([line for line in lines])

        if section == -1:  # special -1 for debug
            lines.append("---------- Debug -----------")
            lines.append(f"spectrum bbox: {self.spectrumbbox}")
            lines.append(f"Good range = {self.goodSpectrumMinY},{self.goodSpectrumMaxY}")
            return '\n'.join([line for line in lines])

        return

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

        return
