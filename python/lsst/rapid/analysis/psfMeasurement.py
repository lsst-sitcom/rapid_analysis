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

import numpy as np
from lsst.atmospec.processStar import ProcessStarTask
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
from lsst.meas.base import MeasurementError
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
import lsst.afw.display as afwDisplay


class QuickPsfMeasurement():

    def __init__(self, display=None, **kwargs):
        self.display = None
        if display:
            self.display = display

        procStarConfig = ProcessStarTask.ConfigClass()
        procStarConfig.mainSourceFindingMethod = 'BRIGHTEST'
        procStarConfig.mainStarNsigma = 10
        procStarConfig.mainStarGrowIsotropic = False
        procStarConfig.mainStarGrow = 1
        self.processStarTask = ProcessStarTask(config=procStarConfig)

        self.centroidName = "base_SdssCentroid"
        self.shapeName = "base_SdssShape"
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.schema.getAliasMap().set("slot_Centroid", self.centroidName)
        self.schema.getAliasMap().set("slot_Shape", self.shapeName)
        self.control = measBase.SdssCentroidControl()
        self.centroider = measBase.SdssCentroidAlgorithm(self.control, self.centroidName, self.schema)
        self.sdssShape = measBase.SdssShapeControl()
        self.shaper = measBase.SdssShapeAlgorithm(self.sdssShape, self.shapeName, self.schema)
        self.table = afwTable.SourceTable.make(self.schema)

    def _getDayObsSeqNumFromExpId(self, expId):
        return self.butler.queryMetadata('raw', ['dayObs', 'seqNum'], expId=expId)[0]

    def run(self, exp, doDisplay=False):
        median = np.nanmedian(exp.image.array)
        exp.image -= median
        sources = self.processStarTask.findObjects(exp)
        if doDisplay:
            if self.display is None:
                raise RuntimeError("Display failed as no display provided during init()")
            self.display.mtv(exp)

        fpSet = sources.getFootprints()
        print(f"Found {len(fpSet)} sources in exposure")

        xxs, yys, centroids = [], [], []
        nMeasured = 0
        for fp in fpSet:
            try:
                src = self.table.makeRecord()
                src.setFootprint(fp)
                self.centroider.measure(src, exp)
                self.shaper.measure(src, exp)

                xxs.append(np.sqrt(src['base_SdssShape_xx'])*2.355*.1)  # 2.355 for FWHM, .1 for platescale
                yys.append(np.sqrt(src['base_SdssShape_yy'])*2.355*.1)
                centroids.append((src['base_SdssCentroid_x'], src['base_SdssCentroid_y']))
                nMeasured += 1
                if doDisplay:  # TODO: Add buffering? Messier due to optional display
                    self.display.dot(src.getShape(), *src.getCentroid(), ctype=afwDisplay.BLUE)
            except MeasurementError:
                pass

        print(f"Measured {nMeasured} of {len(fpSet)} sources in exposure")

        medianXx = np.nanmedian(xxs)
        medianYy = np.nanmedian(yys)

        print(f"Median SDSS shape (x,y) = ({medianXx:.3f}, {medianYy:.3f}) FWHM arcsec")

        exp.image += median  # put background back in
        return medianXx, medianYy

    def runSlow(self, exp):

        imCharConfig = CharacterizeImageTask.ConfigClass()
        imCharConfig.doMeasurePsf = True
        imCharConfig.doApCorr = False
        imCharConfig.doDeblend = False
        imCharConfig.repair.cosmicray.nCrPixelMax = 200000
        imCharTask = CharacterizeImageTask(config=imCharConfig)
        _ = imCharTask.run(exp)

        psf = exp.getPsf()
        psfShape = psf.computeShape()

        ixx = psf.computeShape().getIxx()
        iyy = psf.computeShape().getIyy()
        ixx = np.sqrt(ixx)*2.355*.1
        iyy = np.sqrt(iyy)*2.355*.1
        print(f"Psf shape from imChar task (x,y) = ({ixx:.3f}, {iyy:.3f}) FWHM arcsec")

        return psfShape


if __name__ == '__main__':
    from lsst.rapid.analysis.bestEffort import BestEffortIsr
    REPODIR = '/project/shared/auxTel/'
    bestEffort = BestEffortIsr(REPODIR)
    dataId = {'dayObs': '2020-02-18', 'seqNum': 82}
    exp = bestEffort.getExposure(dataId)
    qm = QuickPsfMeasurement()
    qm.run(exp)
    qm.runSlow(exp)
