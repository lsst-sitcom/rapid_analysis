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

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

import lsst.daf.persistence as dafPersist

CALIB_VALUES = ['FlatField position', 'Park position']
SOUTHPOLESTAR = 'HD 185975'


class NightReporter():

    def __init__(self, repoDir, dayObs):
        self.butler = dafPersist.Butler(repoDir)
        self.dayObs = dayObs
        self.data = {}
        self.auxTelLocation = EarthLocation(lat=-30.244639*u.deg, lon=-70.749417*u.deg, height=2663*u.m)
        self.rebuild()

    def rebuild(self, dayObs=None):
        dayToUse = self.dayObs
        if dayObs:
            dayToUse = dayObs
        self.data = self._scrapeData(dayToUse)

    def _scrapeData(self, dayObs):
        # TODO: add skipping files we already have
        _data = {}
        seqNums = sorted(self.butler.queryMetadata('raw', 'seqNum', dayObs=dayObs))
        for seqNum in sorted(seqNums):
            md = self.butler.get('raw_md', dayObs=dayObs, seqNum=seqNum)
            _data[seqNum] = md.toDict()
        print(f"Loaded data for seqNums {sorted(seqNums)[0]} to {sorted(seqNums)[-1]}")
        return _data

    def getUniqueValuesForKey(self, key, ignoreCalibs=True):
        values = []
        for seqNum in self.data.keys():
            v = self.data[seqNum][key]
            if ignoreCalibs is True and v in CALIB_VALUES:
                continue
            values.append(v)
        return list(set(values))

    def makePolarPlot(self, azimuthsInDegrees, zenithAngles, marker="*-",
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

    def makePolarPlotForObjects(self, objects, colorMap=None, withLines=True):
        objects = self._safeListArg(objects)
        _ = plt.figure(figsize=(10, 10))

        marker = "*"
        if withLines:
            marker += '-'

        for i, obj in enumerate(objects):
            azs = self.getAllValuesForKVPair('AZSTART', ("OBJECT", obj))
            els = self.getAllValuesForKVPair('ELSTART', ("OBJECT", obj))
            assert(len(azs) == len(els))
            if len(azs) == 0:
                print(f"WARNING: found no alt/az data for {obj}")
            zens = [90 - el for el in els]
            color = None
            if colorMap:
                color = colorMap[obj]
            ax = self.makePolarPlot(azs, zens, marker=marker, title=None, makeFig=False,
                                    color=color, objName=obj)
        lgnd = ax.legend(bbox_to_anchor=(1.05, 1), prop={'size': 15}, loc='upper left')
        for h in lgnd.legendHandles:
            size = 14
            if '-' in marker:
                size += 5
            h._legmarker.set_markersize(size)

    def getAllValuesForKVPair(self, keyToGet, keyValPairAsTuple, uniqueOnly=False):
        """e.g. all the RA values for OBJECT=='HD 123'"""
        ret = []
        for seqNum in self.data.keys():
            if self.data[seqNum][keyValPairAsTuple[0]] == keyValPairAsTuple[1]:
                ret.append(self.data[seqNum][keyToGet])
        if uniqueOnly:
            return list(set(ret))
        return ret

    @staticmethod
    def makeStarColorMapDict(stars):
        colorMap = {}
        colors = cm.rainbow(np.linspace(0, 1, len(stars)))
        for i, star in enumerate(stars):
            colorMap[star] = colors[i]
        return colorMap

    def getObjectValues(self, key, objName):
        return self.getAllValuesForKVPair(key, ('OBJECT', objName), uniqueOnly=False)

    def getAllHeaderKeys(self):
        return list(list(self.data.items())[0][1].keys())

    def airMassFromHeader(self, header):
        time = Time(header['DATE-OBS'])
        skyLocation = SkyCoord(header['RASTART'], header['DECSTART'], unit=u.deg)
        altAz = AltAz(obstime=time, location=self.auxTelLocation)
        observationAltAz = skyLocation.transform_to(altAz)
        return observationAltAz.secz.value

    def calcObjectAirmasses(self, objects):
        airMasses = {}
        for star in stars:
            seqNums = self.getObjectValues('SEQNUM', star)
            airMasses[star] = [(self.airMassFromHeader(self.data[seqNum]),
                                self.data[seqNum]['MJD'])for seqNum in sorted(seqNums)]
        return airMasses

    def getObservedObjects(self):
        return self.getUniqueValuesForKey('OBJECT')

    def plotPerObjectAirMass(self, objects=None, versusMjd=True):
        if not objects:
            objects = self.getObservedObjects()

        objects = self._safeListArg(objects)

        # lazy to always recalculate but it's not *that* slow
        # and optionally passing around can be messy
        # TODO: keep some of this in class state
        airMasses = self.calcObjectAirmasses(objects)

        _ = plt.figure(figsize=(10, 6))
        for star in objects:
            ams, times = np.asarray(airMasses[star])[:, 0], np.asarray(airMasses[star])[:, 1]
            if versusMjd:
                plt.plot(times, ams, '*', color=colorMap[star], label=star, ms=10)
            else:
                plt.plot(ams, '*', color=colorMap[star], label=star, ms=10)
        plt.ylabel('Airmass', fontsize=20)
        _ = plt.legend(bbox_to_anchor=(1, 1.025), prop={'size': 15}, loc='upper left')

    @staticmethod
    def _safeListArg(arg):
        if type(arg) == str:
            return [arg]
        assert(type(arg) == list), f"Expect list, got {arg}"
        return arg


if __name__ == '__main__':
    repodir = '/project/shared/auxTel/'
    nightReporter = NightReporter(repodir, "2020-02-20")
    stars = nightReporter.getUniqueValuesForKey('OBJECT')
    import ipdb as pdb; pdb.set_trace()
    colorMap = nightReporter.makeStarColorMapDict(stars)
    nightReporter.makePolarPlotForObjects(stars, colorMap, withLines=False)
