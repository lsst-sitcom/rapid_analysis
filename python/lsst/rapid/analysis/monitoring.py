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
# import lsst.daf.persistence.butlerExceptions as butlerExcept
# from .bestEffort import BestEffortIsr ### xxx put this back in
from lsst.rapid.analysis.bestEffort import BestEffortIsr
from time import sleep

# TODO: maybe add option to create display and return URL?


class Monitor():
    cadence = 2  # in seconds
    runIsr = True

    def __init__(self, repoDir, fireflyDisplay, **kwargs):
        """"""
        self.repoDir = repoDir
        self.display = fireflyDisplay
        self.butler = dafPersist.Butler(repoDir)
        self.bestEffort = BestEffortIsr(repoDir)

    def reloadButler(self):
        self.butler = dafPersist.Butler(self.repoDir)

    def _getLatestVisitNum(self):
        return self.butler.queryMetadata('raw', 'visit')[-1]

    def _getDayObsSeqNumFromVisitNum(self, visitNum):
        return self.butler.queryMetadata('raw', ['dayObs', 'seqNum'], visit=visitNum)[0]

    def _getLatestImageDataIdAndVisitNum(self):
        visitNum = self._getLatestVisitNum()
        dayObs, seqNum = self._getDayObsSeqNumFromVisitNum(visitNum)
        dataId = {'dayObs': dayObs, 'seqNum': seqNum}
        return dataId, visitNum

    def run(self, duration=-1):
        lastDisplayed = -1

        # while True:
        for i in range(10):
            print(f'loop {i}')
            self.reloadButler()  # must be at start not end due to continue

            dataId, visit = self._getLatestImageDataIdAndVisitNum()

            if lastDisplayed == visit:
                sleep(self.cadence)
                continue

            if self.runIsr:
                exp = self.bestEffort.getExposure(dataId)
            else:
                exp = self.butler.get('raw', **dataId)

            self.display.mtv(exp, title=str(dataId))
            lastDisplayed = visit
            sleep(self.cadence)

        return
