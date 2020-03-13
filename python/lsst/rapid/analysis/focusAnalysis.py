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
from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

# TODO: change these back to local .imports
from lsst.rapid.analysis.bestEffort import BestEffortIsr
from lsst.rapid.analysis.quickFrameMeasurement import QuickFrameMeasurement

COLORS = ['b', 'g', 'r']


@dataclass
class FitResult:
    amp: float
    mean: float
    sigma: float


class FocusAnalyzer():

    def __init__(self, repoDir, **kwargs):

        self._butler = dafPersist.Butler(repoDir)
        self._bestEffort = BestEffortIsr(repoDir, **kwargs)
        self._quickMeasure = QuickFrameMeasurement()

        self.spectrumBoxOffsets = [882, 1170, 1467]
        self.spectrumHalfWidth = 100
        self.spectrumBoxLength = 20

    def _checkImageIsDispersed(self, dataId):
        filterFull = self._butler.queryMetadata('raw', 'filter', dataId)[0]
        filt, grating = filterFull.split('~')
        if not grating or grating.startswith('EMPTY'):
            return False
        return True

    @staticmethod
    def _getFocusFromHeader(exp):
        return float(exp.getMetadata()["FOCUSZ"])

    def _getBboxes(self, centroid):
        x, y = centroid
        bboxes = []

        for offset in self.spectrumBoxOffsets:
            bbox = geom.Box2I(geom.Point2I(x-self.spectrumHalfWidth, y+offset),
                              geom.Point2I(x+self.spectrumHalfWidth, y+offset+self.spectrumBoxLength))
            bboxes.append(bbox)
        return bboxes

    @staticmethod
    def gauss(x, *pars):
        amp, mean, sigma = pars
        return amp*np.exp(-(x-mean)**2/(2.*sigma**2))

    def findBestFocii(self, dayObs, seqNums, display=False):
        """Find the best focus in wavelength bands for the specified list of
        seqNums"""
        fitData = {}

        for seqNum in seqNums:
            fitData[seqNum] = {}
            exp = self._bestEffort.getExposure({'dayObs': dayObs, 'seqNum': seqNum})
            quickMeasResult = self._quickMeasure.run(exp)
            centroid = quickMeasResult.brightestObjCentroid

            if display:
                plt.imshow(exp.image.array, norm=LogNorm())
                plt.show()

            bboxes = self._getBboxes(centroid)
            for i, bbox in enumerate(bboxes):
                data1d = np.mean(exp[bbox].image.array, axis=0)
                data1d -= np.median(data1d)
                xs = np.arange(len(data1d))

                amp = np.max(data1d)
                mean = np.argmax(data1d)
                sigma = np.mean(quickMeasResult.brightestObj_xXyY) * 2.355 * 1.8
                p0 = amp, mean, sigma
                coeffs, var_matrix = curve_fit(self.gauss, xs, data1d, p0=p0)
                if display:
                    plt.plot(xs, data1d, f'{COLORS[i]}x')
                    highResX = np.linspace(0, len(data1d), 1000)
                    plt.plot(highResX, self.gauss(highResX, *coeffs), 'k-')

                import ipdb as pdb; pdb.set_trace()

                fitData[seqNum][i] = FitResult(amp=coeffs[0], mean=coeffs[1], sigma=coeffs[2])

            focuserPosition = self._getFocusFromHeader(exp)

        return



if __name__ == '__main__':
    repoDir = '/project/shared/auxTel/'
    analyzer = FocusAnalyzer(repoDir)
    analyzer.findBestFocii('2020-02-20', [485])
    dataId = {'dayObs': '2020-02-20', 'seqNum':485}
    import ipdb as pdb; pdb.set_trace()
