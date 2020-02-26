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

__all__ = ['makePolarPlot']

import matplotlib.pyplot as plt


def humanNameForCelestialObject(objName):
    """Returns a list of all human names for obj, or [] if none are found."""
    from astroquery.simbad import Simbad
    results = []
    try:
        simbadResult = Simbad.query_objectids(objName)
        for row in simbadResult:
            if row['ID'].startswith('NAME'):
                results.append(row['ID'].replace('NAME ', ''))
        return results
    except Exception:
        return []  # same behavior as for found but un-named objects


def _getAltAzZenithsFromSeqNum(butler, dayObs, seqNumList):
    azimuths, elevations, zeniths = [], [], []
    for seqNum in seqNumList:
        md = butler.get("raw_md", dayObs=dayObs, seqNum=seqNum)
        elevations.append(md['ELSTART'])
        zeniths.append(90-md['ELSTART'])
        azimuths.append(md['AZSTART'])
    return azimuths, elevations, zeniths


def makePolarPlot(butler, dayObs, seqMin, seqMax, seqNumList=None, returnData=False):
    """Make a polar plot of the azimuth and zenith angles over sequence nums.

    For the given dayObs, plots the range (seqMin..seqMax) or, for
    dis-contiguous visits, a list of sequence numbers can be specified instead.
    seqNumList is specified, seqMin and seqMax are ignored.
    """
    if not seqNumList:
        seqNumList = [i for i in range(seqMin, seqMax+1)]

    az, el, zen = _getAltAzZenithsFromSeqNum(butler, dayObs, seqNumList)
    _ = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, polar=True)
    ax.plot(az, zen, 'or')
    title = f"Polar coverage - {dayObs} seqNums {seqMin}..{seqMax}"
    ax.set_title(title, va='bottom')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)
    if returnData:
        return {'azimuths': az, 'elevations': el, 'zenithAngles': zen}
