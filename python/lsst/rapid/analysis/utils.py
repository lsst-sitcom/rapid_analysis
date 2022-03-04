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

import eups
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import lsst.afw.detection as afwDetect
import lsst.afw.math as afwMath
import lsst.pipe.base as pipeBase
import logging

from lsst.obs.lsst.translators.lsst import FILTER_DELIMITER
from lsst.obs.lsst.translators.latiss import AUXTEL_LOCATION

from astro_metadata_translator import ObservationInfo

__all__ = ["SIGMATOFWHM", "FWHMTOSIGMA", "EFD_CLIENT_MISSING_MSG", "GOOGLE_CLOUD_MISSING_MSG",
           "AUXTEL_LOCATION", "countPixels", "quickSmooth", "argMax2d", "getImageStats", "detectObjectsInExp",
           "humanNameForCelestialObject", "getFocusFromHeader", "checkRubinTvExternalPackages",
           "dayObsIntToString", "dayObsSeqNumToVisitId"]


SIGMATOFWHM = 2.0*np.sqrt(2.0*np.log(2.0))
FWHMTOSIGMA = 1/SIGMATOFWHM

EFD_CLIENT_MISSING_MSG = ('ImportError: lsst_efd_client not found. Please install with:\n'
                          '    pip install lsst-efd-client')

GOOGLE_CLOUD_MISSING_MSG = ('ImportError: Google cloud storage not found. Please install with:\n'
                            '    pip install google-cloud-storage')


def countPixels(maskedImage, maskPlane):
    bit = maskedImage.mask.getPlaneBitMask(maskPlane)
    return len(np.where(np.bitwise_and(maskedImage.mask.array, bit))[0])


def quickSmooth(data, sigma=2):
    kernel = [sigma, sigma]
    smoothData = gaussian_filter(data, kernel, mode='constant')
    return smoothData


def argMax2d(array):
    """Get the index of the max value of an array and whether it's unique.

    If its not unique, returns a list of the other locations containing the
    maximum value, e.g. returns

    (12, 34), False, [(56,78), (910, 1112)]

    Returns
    -------
    maxLocation : `tuple`
        The coords of the first instance of the max value

    unique : `bool`
        Whether it's the only location

    otherLocations : `list` of `tuple`
        List of the other max values' locations, empty if False
    """
    uniqueMaximum = False
    maxCoords = np.where(array == np.max(array))
    maxCoords = [coord for coord in zip(*maxCoords)]  # list of coords as tuples
    if len(maxCoords) == 1:  # single unambiguous value
        uniqueMaximum = True

    return maxCoords[0], uniqueMaximum, maxCoords[1:]


def dayObsIntToString(dayObs):
    assert isinstance(dayObs, int)
    dStr = str(dayObs)
    assert len(dStr) == 8
    return '-'.join([dStr[0:4], dStr[4:6], dStr[6:8]])


def dayObsSeqNumToVisitId(dayObs, seqNum):
    """Get the visit id for a given dayObs/seqNum.

    dayObs : int
        The dayObs.

    seqNum : int
        The seqNum.

    Returns
    -------
    visitId : int
        The visitId.

    Notes
    -----
    TODO: Remove this horrible hack once DM-30948 makes this possible
    programatically/via the butler.
    """
    if dayObs < 19700101 or dayObs > 35000101:
        raise ValueError(f'dayObs value {dayObs} outside plausible range')
    return int(f"{dayObs}{seqNum:05}")


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
    result.focus = md.get('FOCUSZ')

    data = exp.image.array
    result.maxValue = np.max(data)

    peak, uniquePeak, otherPeaks = argMax2d(data)
    result.maxPixelLocation = peak
    result.multipleMaxPixels = uniquePeak

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
        md = butler.get('raw.metadata', day_obs=dayObs, seq_num=seqNum, detector=0)
        obsInfo = ObservationInfo(md)
        alt = obsInfo.altaz_begin.alt.value
        az = obsInfo.altaz_begin.az.value
        elevations.append(alt)
        zeniths.append(90-alt)
        azimuths.append(az)
    return azimuths, elevations, zeniths


def getFocusFromHeader(exp):
    md = exp.getMetadata()
    if 'FOCUSZ' in md:
        return md['FOCUSZ']
    return None


def checkRubinTvExternalPackages(exitIfNotFound=True, logger=None):
    if not logger:
        logger = logging.getLogger(__name__)

    hasGoogleStorage = False
    hasEfdClient = False
    try:
        from google.cloud import storage  # noqa: F401
        hasGoogleStorage = True
    except ImportError:
        pass

    try:
        from lsst_efd_client import EfdClient  # noqa: F401
        hasEfdClient = True
    except ImportError:
        pass

    if not hasGoogleStorage:
        logger.warn(GOOGLE_CLOUD_MISSING_MSG)

    if not hasEfdClient:
        logger.warn(EFD_CLIENT_MISSING_MSG)

    if exitIfNotFound and (not hasGoogleStorage or not hasEfdClient):
        exit()


def checkStackSetup():
    ep = eups.Eups()
    for tag in ep.findSetupProduct("lsst_distrib").tags:
        if tag != "current":
            print(f"You are running {tag} of lsst_distrib")

    products = ep.findProducts(version='LOCAL:*', tags='setup')
    if products:
        print("\nLocally setup packages:")
        print("-----------------------")
        for product in products:
            print(f"{product.name} setup at {product.dir}")
    else:
        print("\nNo packages are setup locally.")
