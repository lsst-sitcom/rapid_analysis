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
        self.butler = dafPersist.Butler(repoDir)
        self.bestEffort = BestEffortIsr(repoDir)

    def _getLatestExpId(self):
        return sorted(self.butler.queryMetadata('raw', 'expId'))[-1]

    def _getDayObsSeqNumFromExpId(self, expId):
        return self.butler.queryMetadata('raw', ['dayObs', 'seqNum'], expId=expId)[0]

    def _getLatestImageDataIdAndExpId(self):
        expId = self._getLatestExpId()
        dayObs, seqNum = self._getDayObsSeqNumFromExpId(expId)
        dataId = {'dayObs': dayObs, 'seqNum': seqNum}
        return dataId, expId

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

            print(f"Displaying {dataId}...")
            self.display.mtv(exp, title=str(dataId))
            lastDisplayed = expId

        return
