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


from .animation import Animator
import lsst.daf.butler as dafButler


class RoutineAnimator():
    def __init__(self, repoPath):
        self.collections = ['LATISS/raw/all']
        self.butler = dafButler.Butler(repoPath, instrument='LATISS', collections=self.collections)

        self.animator = Animator()

        self.outputPath = '/home/mfl/animatorOutput/main/'
        self.outputFilename = 'somethingOrOther.mp4'
        self.skipTypes = ['BIAS', 'DARK', 'FLAT']

    def getMostRecentDayObs(self):
        where = "instrument='LATISS' AND exposure.day_obs>20210101"
        records = butler.registry.queryDimensionRecords('exposure', where=where,
                                                        collections=self.collections)
        return max(r.day_obs for r in records)


    def _isOnSky(self, dataId):
        if dataId['imageType'] not in self.skipTypes:
            return True
        return False

    def run(self):
        dataProcuctToPlot = 'quickLookExp'

        dataIds = self.getDataIds()

        animator = Animator(self.butler, dataIds, self.outputPath, self.outputFilename,
                            dataProcuctToPlot=dataProcuctToPlot,
                            remakePngs=False,
                            debug=False,
                            clobberVideoAndGif=True,
                            plotObjectCentroids=True)
        animator.run()
