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

from lsst.ip.isr import IsrTask
import lsst.daf.persistence as dafPersist
import lsst.daf.persistence.butlerExceptions as butlerExcept
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask

# TODO: turn prints into log messages
# TODO: refactor if necessary
# TODO: add attempt for fringe once registry & templates are fixed


class BestEffortIsr():

    def __init__(self, repodir):
        self.repodir = repodir
        self.butler = dafPersist.Butler(repodir)
        self.imCharConfig = CharacterizeImageTask.ConfigClass()
        self.imCharConfig.doMeasurePsf = False
        self.imCharConfig.doApCorr = False
        self.imCharConfig.doDeblend = False
        self.imCharTask = CharacterizeImageTask(config=self.imCharConfig)

    def reloadButler(self):
        self.butler = dafPersist.Butler(self.repodir)

    def _repairCosmics(self, exposure):
        try:
            print("Running cosmic ray repair")
            exp = self.imCharTask.run(exposure).exposure
            return exp
        except Exception as e:
            print(f'During CR repair caught: {e}')
            # If the failed attempt turns out to mess up the image in place
            # make copy before running and return original here
            return exp

    def getExposure(self, visitNum):
        try:
            raw = self.butler.get('raw', visit=visitNum)
        except butlerExcept.NoResults:
            raise RuntimeError(f"Failed to retrieve raw for visit {visitNum}")

        isrConfig = IsrTask.ConfigClass()
        isrConfig.doWrite = False  # always off

        # initially all off
        isrConfig.doBias = False
        isrConfig.doDark = False
        isrConfig.doFlat = False
        isrConfig.doLinearize = False
        isrConfig.doFringe = False
        isrConfig.doDefect = False

        isrParts = ['bias', 'dark', 'flat', 'linearizer', 'defect']
        isrDict = {}
        for component in isrParts:
            try:
                item = self.butler.get(component, visit=visitNum)
                isrDict[component] = item
            except AttributeError as e:  # catches mapper problems
                print(f'Caught {e} - update your mapper?')
            except (butlerExcept.NoResults, RuntimeError):
                pass

        # ugly block, but less ugly than setattr and ''.capitalize() etc
        if 'bias' in isrDict.keys():
            isrConfig.doBias = True
            print('Running with bias subtraction')
        if 'dark' in isrDict.keys():
            isrConfig.doDark = True
            print('Running with dark correction subtraction')
        if 'flat' in isrDict.keys():
            isrConfig.doFlat = True
            print('Running with flat fielding')
        if 'linearizer' in isrDict.keys():
            isrConfig.doLinearize = True
            print('Running with linearity correction')
        if 'fringe' in isrDict.keys():
            isrConfig.doFringe = True
            print('Running with fringe correction')
        if 'defect' in isrDict.keys():
            isrConfig.doDefect = True
            print('Running with defect correction')

        isrTask = IsrTask(config=isrConfig)
        postIsr = isrTask.run(raw, **isrDict).exposure
        postCosmicRepair = self._repairCosmics(postIsr)

        return postCosmicRepair
