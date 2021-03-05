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

__all__ = ['makePolarPlot', 'detectObjectsInExp']

import matplotlib.pyplot as plt
import numpy as np
import lsst.afw.detection as afwDetect
import lsst.afw.math as afwMath
import lsst.pipe.base as pipeBase

from lsst.obs.lsst.translators.lsst import FILTER_DELIMITER
from astro_metadata_translator import ObservationInfo


def countPixels(maskedImage, maskPlane):
    bit = maskedImage.mask.getPlaneBitMask(maskPlane)
    return len(np.where(np.bitwise_and(maskedImage.mask.array, bit))[0])


def argMax2d(array):
    """Get the index of the max value of an array.

    Actually for n dimensional, but easier to recall method if misnamed."""
    return np.unravel_index(np.argmax(array, axis=None), array.shape)


def getImageStats(exp):
    result = pipeBase.Struct()

    info = exp.getInfo()
    vi = info.getVisitInfo()
    expTime = vi.getExposureTime()
    md = exp.getMetadata()
    obsInfo = ObservationInfo(md, subset={'object'})

    obj = obsInfo.object
    mjd = vi.getDate().get()
    result.object = obj
    result.mjd = mjd

    fullFilterString = info.getFilterLabel().physicalLabel
    filt = fullFilterString.split(FILTER_DELIMITER)[0]
    grating = fullFilterString.split(FILTER_DELIMITER)[1]

    airmass = vi.getBoresightAirmass()
    rotangle = vi.getBoresightRotAngle().asDegrees()

    azAlt = vi.getBoresightAzAlt()
    az = azAlt[0].asDegrees()
    el = azAlt[1].asDegrees()

    result.expTime = expTime
    result.filter = filt
    result.grating = grating
    result.airmass = airmass
    result.rotangle = rotangle
    result.az = az
    result.el = el

    data = exp.image.array
    result.maxValue = np.max(data)
    result.maxPixelLocation = argMax2d(data)
    result.nBadPixels = countPixels(exp.maskedImage, 'BAD')
    result.nSatPixels = countPixels(exp.maskedImage, 'SAT')
    result.percentile99 = np.percentile(data, 99)
    result.percentile9999 = np.percentile(data, 99.99)

    sctrl = afwMath.StatisticsControl()
    sctrl.setNumSigmaClip(5)
    sctrl.setNumIter(2)
    statTypes = afwMath.MEANCLIP | afwMath.STDEVCLIP
    stats = afwMath.makeStatistics(exp.maskedImage, statTypes, sctrl)
    std, stderr = stats.getResult(afwMath.STDEVCLIP)
    mean, meanerr = stats.getResult(afwMath.MEANCLIP)

    result.clippedMean = mean
    result.clippedStddev = std

    return result


def detectObjectsInExp(exp, nSigma=10, nPixMin=10, grow=0):
    """Return the footPrintSet for the objects in a postISR exposure."""
    median = np.nanmedian(exp.image.array)
    exp.image -= median

    threshold = afwDetect.Threshold(nSigma, afwDetect.Threshold.STDEV)
    footPrintSet = afwDetect.FootprintSet(exp.getMaskedImage(), threshold, "DETECTED", nPixMin)
    if grow > 0:
        isotropic = True
        footPrintSet = afwDetect.FootprintSet(footPrintSet, grow, isotropic)

    exp.image += median  # add back in to leave background unchanged
    return footPrintSet


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
