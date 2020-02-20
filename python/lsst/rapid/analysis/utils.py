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

import matplotlib.pyplot as plt


def _getAltAzZenithsFromSeqNum(butler, dayObs, seqMin, seqMax):
    azimuths, elevations, zeniths = [], [], []
    for seqNum in range(seqMin, seqMax+1):
        md = butler.get("raw_md", dayObs=dayObs, seqNum=seqNum)
        elevations.append(md['ELSTART'])
        zeniths.append(90-md['ELSTART'])
        azimuths.append(md['AZSTART'])
    return azimuths, elevations, zeniths


def makePolarPlot(butler, dayObs, seqMin, seqMax):
    az, el, zen = _getAltAzZenithsFromSeqNum(butler, dayObs, seqMin, seqMax)
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, polar=True)
    ax.plot(az, zen, 'or')
    title = "Polar coverage -"
    if dayObs:
        title += f'{dayObs}'
    if seqMin and seqMax:
        title += f'seqNums {seqMin}..{seqMax}'
    ax.set_title(title, va='bottom')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)
    return fig
