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
# TODO: readd defects attempt once you can catch OperationalError


class BestEffortIsr():

    def __init__(self, repodir, defaultExtraIsrOptions={}):
        """defaultExtraIsrOptions is a dict of options applied to all images"""
        self.repodir = repodir
        self.butler = dafPersist.Butler(repodir)
        self.defaultExtraIsrOptions = defaultExtraIsrOptions
        self.imCharConfig = CharacterizeImageTask.ConfigClass()
        self.imCharConfig.doMeasurePsf = False
        self.imCharConfig.doApCorr = False
        self.imCharConfig.doDeblend = False
        self.imCharConfig.repair.cosmicray.nCrPixelMax = 100000
        self.imCharTask = CharacterizeImageTask(config=self.imCharConfig)

    def reloadButler(self):
        self.butler = dafPersist.Butler(self.repodir)

    def _repairCosmics(self, exposure):
        # TODO: make this a primitive
        try:
            print("Running cosmic ray repair")
            exposure = self.imCharTask.run(exposure).exposure
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
            else:
                print(f"Override option {option} not found in isrConfig")

    @staticmethod
    def _parseVisitOrDataId(visitOrDataId, **kwargs):
        if type(visitOrDataId) == int:
            _dataId = {'visit': visitOrDataId}
        elif type(visitOrDataId) == dict:
            _dataId = visitOrDataId
            _dataId.update(kwargs)
        else:
            raise RuntimeError(f"Invalid visit or dataId type {visitOrDataId}")
        return _dataId

    def getExposure(self, visitOrDataId, extraOptions={}, skipCosmics=False, **kwargs):
        """extraOptions is a dict of options applied to this image only"""
        dataId = self._parseVisitOrDataId(visitOrDataId, **kwargs)

        try:
            raw = self.butler.get('raw', **dataId)
        except butlerExcept.NoResults:
            raise RuntimeError(f"Failed to retrieve raw for visit {dataId}")

        # default options that are probably good for most engineering time
        isrConfig = IsrTask.ConfigClass()
        isrConfig.doWrite = False
        isrConfig.doSaturationInterpolation = False
        isrConfig.overscanNumLeadingColumnsToSkip = 20

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

        isrParts = ['bias', 'dark', 'flat', 'linearizer']  # , 'defects']
        isrDict = {}
        for component in isrParts:
            try:
                item = self.butler.get(component, **dataId)
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
        if 'defects' in isrDict.keys():
            isrConfig.doDefect = True
            print('Running with defect correction')

        isrTask = IsrTask(config=isrConfig)
        postIsr = isrTask.run(raw, **isrDict).exposure
        if not skipCosmics:
            postCosmicRepair = self._repairCosmics(postIsr)
            return postCosmicRepair

        return postIsr
