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

from dataclasses import dataclass
import logging
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import lsst.daf.persistence as dafPersist
from astro_metadata_translator import ObservationInfo

from lsst.atmospec.utils import airMassFromRawMetadata

__all__ = ['NightReporter', 'saveReport', 'loadReport']

CALIB_VALUES = ['FlatField position', 'Park position', 'azel_target']
N_STARS_PER_SYMBOL = 6
MARKER_SEQUENCE = ['*', 'o', "D", 'P', 'v', "^", 's', '.', ',', 'o', 'v', '^',
                   '<', '>', '1', '2', '3', '4', '8', 's', 'p', 'P', '*', 'h',
                   'H', '+', 'x', 'X', 'D', 'd', '|', '_']
SOUTHPOLESTAR = 'HD 185975'

PICKLE_TEMPLATE = "%s.pickle"


# wanted these to be on the class but it doesn't pickle itself nicely
def saveReport(reporter, savePath):
    # the reporter._butler seems to pickle OK but perhaps it should be
    # removed for saving (and re-instantiated from _repoDir on load?)
    filename = os.path.join(savePath, PICKLE_TEMPLATE % reporter.dayObs)
    with open(filename, "wb") as dumpFile:
        pickle.dump(reporter, dumpFile)


def loadReport(loadPath, dayObs):
    filename = os.path.join(loadPath, PICKLE_TEMPLATE % dayObs)
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found")
    with open(filename, "rb") as input_file:
        reporter = pickle.load(input_file)
    return reporter


@dataclass
class ColorAndMarker:
    '''Class for holding colors and marker symbols'''
    color: list
    marker: str = '*'


class NightReporter():

    def __init__(self, repoDir, dayObs, deferLoadingData=False):
        self._supressAstroMetadataTranslatorWarnings()  # call early

        self._butler = dafPersist.Butler(repoDir)
        self._repoDir = repoDir
        self.dayObs = dayObs
        self.data = {}
        self.stars = None
        self.cMap = None
        if not deferLoadingData:
            self.rebuild()

    def _supressAstroMetadataTranslatorWarnings(self):
        """NB: must be called early"""
        logging.basicConfig()
        _astroLogger = logging.getLogger("lsst.obs.lsst.translators.latiss")
        _astroLogger.setLevel(logging.ERROR)

    def rebuild(self, dayObs=None):
        """Reload new observations, or load a different night"""
        dayToUse = self.dayObs
        if dayObs:
            # new day, so blow away old data
            # as scraping skips seqNums we've loaded!
            if dayObs != self.dayObs:
                self.data = {}
                self.dayObs = dayObs
            dayToUse = dayObs
        self._scrapeData(dayToUse)
        self.stars = self.getObservedObjects()
        self.cMap = self.makeStarColorAndMarkerMap(self.stars)

    def _scrapeData(self, dayObs):
        """Load data into self.data skipping as necessary. Don't call directly!

        Don't call directly as the rebuild() function zeros out data for when
        it's a new dayObs."""
        seqNums = sorted(self._butler.queryMetadata('raw', 'seqNum', dayObs=dayObs))
        for seqNum in sorted(seqNums):
            if seqNum in self.data.keys():
                continue
            md = self._butler.get('raw_md', dayObs=dayObs, seqNum=seqNum)
            self.data[seqNum] = md.toDict()
            self.data[seqNum]['ObservationInfo'] = ObservationInfo(md)
        print(f"Loaded data for seqNums {sorted(seqNums)[0]} to {sorted(seqNums)[-1]}")

    def getUniqueValuesForKey(self, key, ignoreCalibs=True):
        values = []
        for seqNum in self.data.keys():
            v = self.data[seqNum][key]
            if ignoreCalibs is True and v in CALIB_VALUES:
                continue
            values.append(v)
        return list(set(values))

    def _makePolarPlot(self, azimuthsInDegrees, zenithAngles, marker="*-",
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

    def makePolarPlotForObjects(self, objects=None, withLines=True):
        if not objects:
            objects = self.stars
        objects = self._safeListArg(objects)

        _ = plt.figure(figsize=(10, 10))

        for i, obj in enumerate(objects):
            azs = self.getAllValuesForKVPair('AZSTART', ("OBJECT", obj))
            els = self.getAllValuesForKVPair('ELSTART', ("OBJECT", obj))
            assert(len(azs) == len(els))
            if len(azs) == 0:
                print(f"WARNING: found no alt/az data for {obj}")
            zens = [90 - el for el in els]
            color = self.cMap[obj].color
            marker = self.cMap[obj].marker
            if withLines:
                marker += '-'

            ax = self._makePolarPlot(azs, zens, marker=marker, title=None, makeFig=False,
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
    def makeStarColorAndMarkerMap(stars):
        markerMap = {}
        colors = cm.rainbow(np.linspace(0, 1, N_STARS_PER_SYMBOL))
        for i, star in enumerate(stars):
            markerIndex = i//(N_STARS_PER_SYMBOL)
            colorIndex = i%(N_STARS_PER_SYMBOL)
            markerMap[star] = ColorAndMarker(colors[colorIndex], MARKER_SEQUENCE[markerIndex])
        return markerMap

    def getObjectValues(self, key, objName):
        return self.getAllValuesForKVPair(key, ('OBJECT', objName), uniqueOnly=False)

    def getAllHeaderKeys(self):
        return list(list(self.data.items())[0][1].keys())

    def _airMassFromHeader(self, header):
        time = Time(header['DATE-OBS'])
        if header['RASTART'] is not None and header['DECSTART'] is not None:
            skyLocation = SkyCoord(header['RASTART'], header['DECSTART'], unit=u.deg)
        elif header['RA'] is not None and header['DEC'] is not None:
            skyLocation = SkyCoord(header['RA'], header['DEC'], unit=u.deg)
        else:
            return 0
        altAz = AltAz(obstime=time, location=self._auxTelLocation)
        observationAltAz = skyLocation.transform_to(altAz)
        return observationAltAz.secz.value

    def _calcObjectAirmasses(self, objects):
        airMasses = {}
        for star in objects:
            seqNums = self.getObjectValues('SEQNUM', star)
            airMasses[star] = [(self._airMassFromHeader(self.data[seqNum]),
                                self.data[seqNum]['MJD'])for seqNum in sorted(seqNums)]
        return airMasses

    def getObservedObjects(self):
        return self.getUniqueValuesForKey('OBJECT')

    def plotPerObjectAirMass(self, objects=None, versusMjd=True, airmassOneAtTop=True):
        if not objects:
            objects = self.stars

        objects = self._safeListArg(objects)

        # lazy to always recalculate but it's not *that* slow
        # and optionally passing around can be messy
        # TODO: keep some of this in class state
        airMasses = self._calcObjectAirmasses(objects)

        _ = plt.figure(figsize=(10, 6))
        for star in objects:
            ams, times = np.asarray(airMasses[star])[:, 0], np.asarray(airMasses[star])[:, 1]
            color = self.cMap[star].color
            marker = self.cMap[star].marker
            if versusMjd:
                plt.plot(times, ams, '*', color=color, marker=marker, label=star, ms=10)
            else:
                plt.plot(ams, '*', color=color, marker=marker, label=star, ms=10)
        plt.ylabel('Airmass', fontsize=20)
        if airmassOneAtTop:
            ax = plt.gca()
            ax.set_ylim(ax.get_ylim()[::-1])
        _ = plt.legend(bbox_to_anchor=(1, 1.025), prop={'size': 15}, loc='upper left')

    def printObsTable(self, imageType=None, tailNumber=0):
        """Print a table of the days observations.

        Parameters
        ----------
        imageType : str
            Only consider images with this image type
        tailNumber : int
            Only print out the last n entries in the night
        """
        lines = []
        if not imageType:
            seqNums = self.data.keys()
        else:
            seqNums = [s for s in self.data.keys()
                       if self.data[s]['ObservationInfo'].observation_type == imageType]

        seqNums = sorted(seqNums)
        for i, seqNum in enumerate(seqNums):
            expTime = self.data[seqNum]['EXPTIME']
            filt = self.data[seqNum]['ObservationInfo'].physical_filter
            imageType = self.data[seqNum]['ObservationInfo'].observation_type
            d1 = self.data[seqNum]['ObservationInfo'].datetime_begin
            obj = self.data[seqNum]['OBJECT']
            if i == 0:
                d0 = d1
            dt = (d1-d0)
            d0 = d1
            timeOfDay = d1.isot.split('T')[1]
            msg = f'{seqNum:4} {imageType:9} {obj:10} {timeOfDay} {filt:25} {dt.sec:6.1f}  {expTime:2.2f}'
            lines.append(msg)

        print(r"{seqNum:4} {imageType:9} {obj:10} {timeOfDay} {filt:25} {dt.sec:6.1f}  {expTime:2.2f}")
        for line in lines[-tailNumber:]:
            print(line)

    def calcShutterOpenEfficiency(self, seqMin=0, seqMax=0):
        if seqMin == 0:
            seqMin = min(self.data.keys())
        if seqMax == 0:
            seqMax = max(self.data.keys())
        assert seqMax > seqMin
        assert (seqMin in self.data.keys())
        assert (seqMax in self.data.keys())

        timeStart = self.data[seqMin]['ObservationInfo'].datetime_begin
        timeEnd = self.data[seqMax]['ObservationInfo'].datetime_end
        expTimeSum = 0
        for seqNum in range(seqMin, seqMax+1):
            if seqNum not in self.data.keys():
                print(f"Warning! No data found for seqNum {seqNum}")
                continue
            expTimeSum += self.data[seqNum]['EXPTIME']

        timeOnSky = (timeEnd - timeStart).sec
        efficiency = expTimeSum/timeOnSky
        print(f"{100*efficiency:.2f}% shutter open in seqNum range {seqMin} and {seqMax}")
        print(f"Total integration time = {expTimeSum:.1f}s")
        return efficiency

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
    colorMap = nightReporter.makeStarColorAndMarkerMap(stars)
    nightReporter.makePolarPlotForObjects(stars, colorMap, withLines=False)
