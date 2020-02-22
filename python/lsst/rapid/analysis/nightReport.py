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
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

CALIB_VALUES = ['FlatField position', 'Park position']


def scrapeData(butler, dayObs):
    data = {}
    seqNums = sorted(butler.queryMetadata('raw', 'seqNum', dayObs=dayObs))
    for seqNum in seqNums:
        md = butler.get('raw_md', dayObs=dayObs, seqNum=seqNum)
        data[seqNum] = md.toDict()
    return data


def getUniqueValuesForKey(data, key, ignoreCalibs=True):
    values = []
    for seqNum in data.keys():
        v = data[seqNum][key]
        if ignoreCalibs is True and v in CALIB_VALUES:
            continue
        values.append(v)
    return list(set(values))


def makePolarPlot(azimuthsInDegrees, zenithAngles, marker="*-",
                  title=None, makeFig=True, color=None, objName=None):
    if makeFig:
        _ = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, polar=True)
    ax.plot([a*np.pi/180 for a in azimuthsInDegrees], zenithAngles, marker, c=color, label=objName)
    if title:
        ax.set_title(title, va='bottom')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)
    return ax


def makePolarPlotForObjects(data, objects, colorMap=None, withLines=True):
    if type(objects) == str:
        objects = [objects]
    _ = plt.figure(figsize=(10, 10))

    marker = "*"
    if withLines:
        marker += '-'

    for i, obj in enumerate(objects):
        azs = getAllValuesForKVPair(data, 'AZSTART', ("OBJECT", obj))
        els = getAllValuesForKVPair(data, 'ELSTART', ("OBJECT", obj))
        assert(len(azs) == len(els))
        if len(azs) == 0:
            print(f"WARNING: found no alt/az data for {obj}")
        zens = [90 - el for el in els]
        color = None
        if colorMap:
            color = colorMap[obj]
        ax = makePolarPlot(azs, zens, marker=marker, title=None, makeFig=False, color=color, objName=obj)
    ax.set_title(f"Polar coverage for {[ob for ob in objects]}", va='bottom')
    lgnd = ax.legend(bbox_to_anchor=(1.05, 1), prop={'size': 15}, loc='upper left')
    for h in lgnd.legendHandles:
        size = 14
        if '-' in marker:
            size += 5
        h._legmarker.set_markersize(size)


def getAllValuesForKVPair(data, keyToGet, keyValPairAsTuple, uniqueOnly=False):
    """e.g. all the RA values for OBJECT=='HD 123'"""
    ret = []
    for seqNum in data.keys():
        if data[seqNum][keyValPairAsTuple[0]] == keyValPairAsTuple[1]:
            ret.append(data[seqNum][keyToGet])
    if uniqueOnly:
        return list(set(ret))
    return ret


def makeStarColorMapDict(stars):
    colorMap = {}
    colors = cm.rainbow(np.linspace(0, 1, len(stars)))
    for i, star in enumerate(stars):
        colorMap[star] = colors[i]
    return colorMap


def getObjectValues(data, key, objName):
    return getAllValuesForKVPair(data, key, ('OBJECT', objName), uniqueOnly=False)


def getAllHeaderKeys(data):
    return list(list(data.items())[0][1].keys())
