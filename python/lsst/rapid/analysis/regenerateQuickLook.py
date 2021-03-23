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
from lsst.rapid.analysis.bestEffort import BestEffortIsr
import lsst.daf.persistence as dafPersist

__all__ = ['regenerateFromDays']


def regenerateFromDays(repoDir, days, clobber=False):
    """days in the form ["2020-02-18", "2020-02-19"]"""
    bestEffort = BestEffortIsr(repoDir)
    butler = dafPersist.Butler(repoDir)
    dataIds = getAllDataIdsAcrossManyDayObs(days, butler)
    regenerateQuickLookExps(dataIds, butler, bestEffort, clobber=clobber)
    fails = checkQuickLookExpsExist(dataIds, butler)
    if fails:
        print(f"Failed to (re)generate quickLookExps for {fails}")
        return fails
    return


def quickLookExpExists(dataId, butler):
    return os.path.exists(butler.getUri('quickLookExp', dataId, write=True))


def getAllDataIdsAcrossManyDayObs(days, butler):
    dataIds = []
    for day in days:
        seqNums = sorted(butler.queryMetadata('raw', 'seqNum', dayObs=day))
        for seqNum in seqNums:
            dataId = {"dayObs": day, "seqNum": seqNum}
            dataIds.append(dataId)
    return dataIds[::-1]  # reversed so we start with most recent first


def regenerateQuickLookExps(dataIds, butler, bestEffortIsr, clobber=False):
    nTotal = len(dataIds)
    for i, dataId in enumerate(dataIds):
        if (not quickLookExpExists(dataId, butler)) or clobber is True:
            print(f"Processing {dataId} - {i} of {nTotal}")
            exp = bestEffortIsr.getExposure(dataId)
            butler.put(exp, 'quickLookExp', dataId)


def checkQuickLookExpsExist(dataIds, butler):
    failures = []
    for i, dataId in enumerate(dataIds):
        if not quickLookExpExists(dataId, butler):
            print(f"Failed to find data for {dataId}")
            failures.append(dataId)
    return failures
