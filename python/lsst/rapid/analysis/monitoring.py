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
import numpy as np
import lsst.daf.persistence as dafPersist
from .bestEffort import BestEffortIsr
from time import sleep

# TODO: maybe add option to create display and return URL?


class Monitor():
    cadence = 1  # in seconds
    runIsr = True

    def __init__(self, repoDir, fireflyDisplay, **kwargs):
        """"""
        self.repoDir = repoDir
        self.display = fireflyDisplay
        self.bestEffort = BestEffortIsr(repoDir)
        self.writePostIsrImages = None
        outpath = os.path.join(repoDir, 'rerun/quickLook')
        self.butler = dafPersist.Butler(inputs=repoDir, outputs=outpath)

    def _getLatestExpId(self):
        return sorted(self.butler.queryMetadata('raw', 'expId'))[-1]

    def _getDayObsSeqNumFromExpId(self, expId):
        return self.butler.queryMetadata('raw', ['dayObs', 'seqNum'], expId=expId)[0]

    def _getLatestImageDataIdAndExpId(self):
        expId = self._getLatestExpId()
        dayObs, seqNum = self._getDayObsSeqNumFromExpId(expId)
        dataId = {'dayObs': dayObs, 'seqNum': seqNum}
        return dataId, expId

    def _calcImageStats(self, exp):
        elements = []
        median = np.median(exp.image.array)
        elements.append(f"Median={median:.2f}")
        mean = np.mean(exp.image.array)
        # elements.append(f"Median={median:.2f}")
        elements.append(f"Mean={mean:.2f}")

        return elements

    def _makeImageInfoText(self, dataId, exp, asList=False):
        elements = []

        imageType = self.butler.queryMetadata('raw', 'imageType', dataId)[0]
        obj = None
        if imageType.upper() not in ['BIAS', 'DARK', 'FLAT']:
            try:
                obj = self.butler.queryMetadata('raw', 'OBJECT', dataId)[0]
                obj = obj.replace(' ', '')
            except Exception:
                pass

        for k, v in dataId.items():  # dataId done per line for vertical display
            elements.append(f"{k}:{v}")

        if obj:
            elements.append(f"{obj}")
        else:
            elements.append(f"{imageType}")

        expTime = exp.getInfo().getVisitInfo().getExposureTime()
        filt = exp.getFilter().getName()

        elements.append(f"{expTime}s exp")
        elements.append(f"{filt}")

        elements.extend(self._calcImageStats(exp))

        if asList:
            return elements
        return " ".join([e for e in elements])

    def _printImageInfo(self, elements):
        size = 3
        top = 3850  # just under title for size=3
        xnom = -600  # 0 is the left edge of the image
        vSpacing = 100  # about right for size=3, can make f(size) if needed
        for i, item in enumerate(elements):
            y = top - (i*vSpacing)
            x = xnom + (size * 18.5 * len(item)//2)
            self.display.dot(str(item), x, y, size, ctype='red', fontFamily="courier")

    def run(self, durationInSeconds=-1):

        if durationInSeconds == -1:
            nLoops = int(1e9)
        else:
            nLoops = int(durationInSeconds//self.cadence)

        lastDisplayed = -1
        for i in range(nLoops):
            dataId, expId = self._getLatestImageDataIdAndExpId()

            if lastDisplayed == expId:
                sleep(self.cadence)
                continue

            if self.runIsr:
                exp = self.bestEffort.getExposure(dataId)
            else:
                exp = self.butler.get('raw', **dataId)

            if self.writePostIsrImages:
                self.butler.put(exp, "quickLookExp", dataId)

            print(f"Displaying {dataId}...")
            imageInfoText = self._makeImageInfoText(dataId, exp, asList=True)
            longTitle = " ".join([s for s in imageInfoText])
            self.display.mtv(exp, title=longTitle)
            self.display.scale('asinh', 'zscale')
            self._printImageInfo(imageInfoText)
            lastDisplayed = expId

        return
