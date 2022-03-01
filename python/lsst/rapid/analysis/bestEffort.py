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

from sqlite3 import OperationalError

import logging
from lsst.ip.isr import IsrTask
import lsst.daf.butler as dafButler
from lsst.daf.butler.registry import ConflictingDefinitionError

from lsst.rapid.analysis.quickLook import QuickLookTask
from lsst.rapid.analysis.butlerUtils import (LATISS_DEFAULT_COLLECTIONS, LATISS_SUPPLEMENTAL_COLLECTIONS,
                                             _repoDirToLocation)

# TODO: add attempt for fringe once registry & templates are fixed

CURRENT_RUN = "LATISS/runs/quickLook/1"
DATASET_NAME = 'quickLookExp'
ALLOWED_REPOS = ['/repo/main', '/repo/LATISS', '/readonly/repo/main']


class BestEffortIsr():

    def __init__(self, repodir='', *,
                 extraCollections=[], defaultExtraIsrOptions={}, doRepairCosmics=True, doWrite=True):
        f"""Instantiate a BestEffortIsr object.

        Acceptable repodir values are currently {ALLOWED_REPOS}

        defaultExtraIsrOptions is a dict of options applied to all images.

        Parameters
        ----------
        collections : `list` of `str`
            The collections to use.

        doRepairCosmics : `bool`, optional
            Repair cosmic ray hits?

        doWrite : `bool`, optional
            Write the outputs to the quickLook rerun/collection?
        """
        if repodir not in ALLOWED_REPOS:
            raise RuntimeError('Currently only NCSA and summit repos are supported')
        self.log = logging.getLogger(__name__)

        location = _repoDirToLocation(repodir)
        LSC = LATISS_SUPPLEMENTAL_COLLECTIONS  # grrr, line lengths
        collections = (LSC[location] if location in LSC.keys() else []) + LATISS_DEFAULT_COLLECTIONS
        self.collections = extraCollections + collections
        self.log.info(f'Instantiating butler with collections={self.collections}')
        self.butler = dafButler.Butler(repodir, collections=self.collections,
                                       instrument='LATISS',
                                       run=CURRENT_RUN if doWrite else None)

        quickLookConfig = QuickLookTask.ConfigClass()
        quickLookConfig.doRepairCosmics = doRepairCosmics
        self.doWrite = doWrite  # the task, as run by run() method, can't do the write, so we handle in here
        self.quickLookTask = QuickLookTask(config=quickLookConfig)

        self.defaultExtraIsrOptions = defaultExtraIsrOptions

        self._cache = {}

    def _applyConfigOverrides(self, config, overrides):
        for option, value in overrides.items():
            if hasattr(config, option):
                setattr(config, option, value)
                self.log.info(f"Set isr config override {option} to {value}")
            else:
                self.log.warning(f"Override option {option} not found in isrConfig")

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
            exp = self.butler.get(DATASET_NAME, **dataId)
            self.log.info("Found a ready-made quickLookExp in the repo. Returning that.")
            return exp
        except LookupError:
            pass

        try:
            raw = self.butler.get('raw', **dataId)
        except LookupError:
            raise RuntimeError(f"Failed to retrieve raw for exp {dataId}") from None

        # default options that are probably good for most engineering time
        isrConfig = IsrTask.ConfigClass()
        isrConfig.doWrite = False  # this task writes separately, no need for this
        isrConfig.doSaturation = True  # saturation very important for roundness measurement in qfm
        isrConfig.doSaturationInterpolation = True
        isrConfig.overscanNumLeadingColumnsToSkip = 5
        isrConfig.overscan.fitType = 'MEDIAN_PER_ROW'

        # apply general overrides
        self._applyConfigOverrides(isrConfig, self.defaultExtraIsrOptions)
        # apply per-image overrides
        self._applyConfigOverrides(isrConfig, extraOptions)

        isrParts = ['camera', 'bias', 'dark', 'flat', 'defects', 'linearizer', 'crosstalk', 'bfKernel',
                    'bfGains', 'ptc']

        isrDict = {}
        for component in isrParts:
            if component in self._cache and component != 'flat':
                self.log.info(f"Using {component} from cache...")
                isrDict[component] = self._cache[component]
                continue
            try:
                # TODO: add caching for flats
                item = self.butler.get(component, **dataId)
                self._cache[component] = item
                isrDict[component] = self._cache[component]
            except (RuntimeError, LookupError, OperationalError):
                pass

        quickLookExp = self.quickLookTask.run(raw, **isrDict, isrBaseConfig=isrConfig,
                                              isGen3=True).outputExposure

        if self.doWrite:
            try:
                self.butler.put(quickLookExp, DATASET_NAME, dataId)
                self.log.info(f'Put quickLookExp for {dataId}')
            except ConflictingDefinitionError:
                self.log.warning('Skipped putting existing exp into collection! (ignore if there was a race)')
                pass

        return quickLookExp


if __name__ == '__main__':
    repodir = '/repo/main'
    bestEffort = BestEffortIsr(repodir, doWrite=True)
    dataId = {'day_obs': 20200315, 'seq_num': 164, 'detector': 0}
    exp = bestEffort.getExposure(dataId)
