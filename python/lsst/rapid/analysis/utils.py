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
from scipy.ndimage.filters import gaussian_filter

import lsst.afw.detection as afwDetect
import lsst.afw.math as afwMath
import lsst.pipe.base as pipeBase
import lsst.log as lsstLog
import lsst.geom as geom

from lsst.obs.lsst.translators.lsst import FILTER_DELIMITER
from lsst.obs.lsst.translators.latiss import AUXTEL_LOCATION

from astro_metadata_translator import ObservationInfo
from astropy.coordinates import SkyCoord, AltAz
import astropy.units as u

__all__ = ["SIGMATOFWHM", "FWHMTOSIGMA", "EFD_CLIENT_MISSING_MSG", "GOOGLE_CLOUD_MISSING_MSG",
           "AUXTEL_LOCATION", "countPixels", "quickSmooth", "argMax2d", "getImageStats", "detectObjectsInExp",
           "humanNameForCelestialObject", "getFocusFromHeader", "checkRubinTvExternalPackages",
           ]


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
        md = butler.get("raw_md", dayObs=dayObs, seqNum=seqNum)
        elevations.append(md['ELSTART'])
        zeniths.append(90-md['ELSTART'])
        azimuths.append(md['AZSTART'])
    return azimuths, elevations, zeniths


def getFocusFromHeader(exp):
    md = exp.getMetadata()
    if 'FOCUSZ' in md:
        return md['FOCUSZ']
    return None


def checkRubinTvExternalPackages(exitIfNotFound=True, logger=None):
    if not logger:
        logger = lsstLog.Log.getDefaultLogger()

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


def getExpPositionOffset(exp1, exp2):
    # need the exps if we want altAz because we need the observation times
    wcs1 = exp1.getWcs()
    wcs2 = exp2.getWcs()
    p1 = wcs1.getSkyOrigin()
    p2 = wcs2.getSkyOrigin()

    vi1 = exp1.getInfo().getVisitInfo()
    vi2 = exp2.getInfo().getVisitInfo()

    # AltAz via astropy
    skyLocation1 = SkyCoord(p1.getRa().asRadians(), p1.getDec().asRadians(), unit=u.rad)
    altAz1 = AltAz(obstime=vi1.date.toPython(), location=AUXTEL_LOCATION)
    obsAltAz1 = skyLocation1.transform_to(altAz1)
    alt1 = geom.Angle(obsAltAz1.alt.degree, geom.degrees)
    az1 = geom.Angle(obsAltAz1.az.degree, geom.degrees)

    skyLocation2 = SkyCoord(p2.getRa().asRadians(), p2.getDec().asRadians(), unit=u.rad)
    altAz2 = AltAz(obstime=vi2.date.toPython(), location=AUXTEL_LOCATION)
    obsAltAz2 = skyLocation2.transform_to(altAz2)
    alt2 = geom.Angle(obsAltAz2.alt.degree, geom.degrees)
    az2 = geom.Angle(obsAltAz2.az.degree, geom.degrees)

    # ra/dec via the stack
    # NB using the wcs not the exp visitInfo as that's unfitted
    # and so is always identical for a given object/pointing
    ra1 = p1[0]
    ra2 = p2[0]
    dec1 = p1[1]
    dec2 = p2[1]

    angular_offset = p1.separation(p2).asArcseconds()
    delta_pixels = angular_offset / wcs1.getPixelScale().asArcseconds()

    ret = pipeBase.Struct(deltaRa=ra1-ra2,
                          deltaDec=dec1-dec2,
                          deltaAlt=alt1-alt2,
                          deltaAz=az1-az2,
                          delta_pixel_magnitude=delta_pixels
                          )

    return ret
