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
from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask

# TODO: turn prints into log messages
# TODO: refactor if necessary
# TODO: add attempt for fringe once registry & templates are fixed


class BestEffortIsr():

    def __init__(self, repodir='', defaultExtraIsrOptions={}, butler=None):
        """Instantiate a BestEffortIsr object.

        If a butler already exists, it can be passed in, otherwise one is
        instantiated using the repodir.

        defaultExtraIsrOptions is a dict of options applied to all images."""
        if not repodir and butler is None:
            raise RuntimeError("You must either supply a repo dir or a butler")
        if repodir and butler is not None:
            msg = "Ambiguous instantiation. You can either supply a repo dir OR a butler, but not both"
            raise RuntimeError(msg)

        if butler:
            self.butler = butler
        else:
            self.butler = dafPersist.Butler(repodir)

        self.defaultExtraIsrOptions = defaultExtraIsrOptions
        self.imCharConfig = CharacterizeImageTask.ConfigClass()
        self.imCharConfig.doMeasurePsf = False
        self.imCharConfig.doApCorr = False
        self.imCharConfig.doDeblend = False
        self.imCharConfig.repair.cosmicray.nCrPixelMax = 200000
        self.imCharTask = CharacterizeImageTask(config=self.imCharConfig)

        self.writePostIsrImages = False

        self._cache = {}

    def _repairCosmics(self, exposure):
        try:
            print("Running cosmic ray repair")
            installPsfTask = InstallGaussianPsfTask()
            installPsfTask.run(exposure)
            # only run the .repair part to avoid general background subtraction
            self.imCharTask.repair.run(exposure)
            return exposure
        except Exception as e:
            print(f'During CR repair caught: {e}')
            # If the failed attempt turns out to mess up the image in place
            # make copy before running and return original here
            return exposure

    @staticmethod
    def _applyConfigOverrides(config, overrides):
        for option, value in overrides.items():
            if hasattr(config, option):
                setattr(config, option, value)
                print(f"Set isr config override {option} to {value}")
            else:
                print(f"WARNING: Override option {option} not found in isrConfig")

    @staticmethod
    def _parseExpIdOrDataId(expIdOrDataId, **kwargs):
        if type(expIdOrDataId) == int:
            _dataId = {'expId': expIdOrDataId}
        elif type(expIdOrDataId) == dict:
            _dataId = expIdOrDataId
            _dataId.update(kwargs)
        else:
            raise RuntimeError(f"Invalid expId or dataId type {expIdOrDataId}")
        return _dataId

    def clearCache(self):
        self._cache = {}

    def getExposure(self, expIdOrDataId, extraOptions={}, skipCosmics=False, **kwargs):
        """extraOptions is a dict of options applied to this image only"""
        dataId = self._parseExpIdOrDataId(expIdOrDataId, **kwargs)

        try:
            raw = self.butler.get('raw', **dataId)
        except butlerExcept.NoResults:
            raise RuntimeError(f"Failed to retrieve raw for exp {dataId}")

        # default options that are probably good for most engineering time
        isrConfig = IsrTask.ConfigClass()
        isrConfig.doWrite = False
        isrConfig.doSaturation = True  # saturation very important for roundness measurement in qfm
        isrConfig.doSaturationInterpolation = True
        isrConfig.overscanNumLeadingColumnsToSkip = 5
        isrConfig.overscan.fitType = 'MEDIAN_PER_ROW'

        # apply general overrides
        self._applyConfigOverrides(isrConfig, self.defaultExtraIsrOptions)
        # apply per-image overrides
        self._applyConfigOverrides(isrConfig, extraOptions)

        # initially all off
        isrConfig.doBias = False
        isrConfig.doDark = False
        isrConfig.doFlat = False
        isrConfig.doLinearize = False
        isrConfig.doFringe = False
        isrConfig.doDefect = False

        isrParts = ['bias', 'dark', 'flat', 'defects']
        isrDict = {}
        for component in isrParts:
            if component in self._cache and component != 'flat':
                print(f"Using {component} from cache...")
                isrDict[component] = self._cache[component]
                continue
            try:
                # TODO: add caching for flats
                item = self.butler.get(component, **dataId)
                self._cache[component] = item
                isrDict[component] = self._cache[component]
            except AttributeError as e:  # catches mapper problems
                print(f'Caught {e} - update your mapper?')
            except (butlerExcept.NoResults, RuntimeError):
                pass

        # ugly block, but less ugly than setattr and ''.capitalize() etc
        print('Running best effort isr...')
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
        if 'defects' in isrDict.keys():
            isrConfig.doDefect = True
            print('Running with defect correction')

        isrTask = IsrTask(config=isrConfig)
        postIsr = isrTask.run(raw, **isrDict).exposure
        if not skipCosmics:
            postIsr = self._repairCosmics(postIsr)

        if self.writePostIsrImages:
            self.butler.put(postIsr, "quickLookExp", dataId)

        return postIsr


if __name__ == '__main__':
    REPODIR = '/project/shared/auxTel/'
    bestEffort = BestEffortIsr(REPODIR)
    bestEffort.writePostIsrImages = True
    dataId = {'dayObs': '2020-02-17', 'seqNum': 244}
    exp = bestEffort.getExposure(dataId)
