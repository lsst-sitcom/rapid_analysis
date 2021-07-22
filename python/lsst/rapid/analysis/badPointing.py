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
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.coordinates import Angle
import astropy.units as u
from astroquery.astrometry_net import AstrometryNet

import lsst.geom as geom
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask, QuickFrameMeasurementTaskConfig

from lsst.atmospec.utils import AUXTEL_LOCATION
from lsst.rapid.analysis.utils import quickSmooth


def getApiKey():
    try:
        return os.environ['ASTROMETRY_NET_API_KEY']
    except KeyError as e:
        msg = "No AstrometryNet API key found. Sign up and get one, set it to $ASTROMETRY_NET_API_KEY"
        raise RuntimeError(msg) from e


def runImchar(exp):
    charConfig = CharacterizeImageConfig()
    charConfig.doMeasurePsf = False
    charConfig.doApCorr = False
    charConfig.doDeblend = False
    charConfig.repair.doCosmicRay = False
    charConfig.repair.doInterpolate = True
    charConfig.detection.minPixels = 500
    charTask = CharacterizeImageTask(config=charConfig)

    charResult = charTask.run(exp)
    return charResult


def runQfm(exp):
    qfmTaskConfig = QuickFrameMeasurementTaskConfig()
    qfmTask = QuickFrameMeasurementTask(config=qfmTaskConfig)
    qfmResult = qfmTask.run(exp)
    return qfmResult


def plot(exp, icSrc, saveAs):
    plt.figure(figsize=(16, 16))
    arr = exp.image.array
    arr = np.clip(arr, 1, 100000)  # This image has some negative values, and this removes them
    arr = quickSmooth(arr)
    plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000), interpolation='Nearest', cmap='gray', origin='lower')

    Ncenter = (700, 900)
    Nlength = 500.0
    NcenterAzEl = (3200, 700)
    Nlabel = 650.0

    vi = exp.getInfo().getVisitInfo()
    az, _ = vi.boresightAzAlt
    _, dec = vi.boresightRaDec
    rotpa = vi.boresightRotAngle

    az = Angle(az.asDegrees(), u.deg)
    dec = Angle(dec.asDegrees(), u.deg)
    rotpa = Angle(rotpa.asDegrees(), u.deg)

    plt.scatter(icSrc['base_SdssCentroid_x'], icSrc['base_SdssCentroid_y'], color='red', marker='x')
    plt.arrow(Ncenter[0], Ncenter[1], -Nlength*np.sin(rotpa), Nlength*np.cos(rotpa),
              color='green', width=20)
    plt.text(Ncenter[0]-Nlabel*np.sin(rotpa), Ncenter[1]+Nlabel*np.cos(rotpa), 'N',
             color='green', fontsize=12, weight='bold')
    plt.arrow(Ncenter[0], Ncenter[1], Nlength*np.cos(rotpa), Nlength*np.sin(rotpa),
              color='green', width=20)
    plt.text(Ncenter[0]+Nlabel*np.cos(rotpa), Ncenter[1]+Nlabel*np.sin(rotpa), 'E',
             color='green', fontsize=12, weight='bold')

    sinTheta = np.cos(AUXTEL_LOCATION.lat) / np.cos(dec) * np.sin(az)
    theta = Angle(np.arcsin(sinTheta))
    rotAzEl = rotpa - theta - Angle(90.0 * u.deg)
    plt.arrow(NcenterAzEl[0], NcenterAzEl[1], - Nlength*np.sin(rotAzEl), Nlength*np.cos(rotAzEl),
              color='cyan', width=20)
    plt.text(NcenterAzEl[0]-Nlabel*np.sin(rotAzEl), NcenterAzEl[1]+Nlabel*np.cos(rotAzEl), 'EL',
             color='cyan', fontsize=12, weight='bold')
    plt.arrow(NcenterAzEl[0], NcenterAzEl[1], Nlength*np.cos(rotAzEl), Nlength*np.sin(rotAzEl),
              color='cyan', width=20)
    plt.text(NcenterAzEl[0]+Nlabel*np.cos(rotAzEl), NcenterAzEl[1]+Nlabel*np.sin(rotAzEl), 'AZ',
             color='cyan', fontsize=12, weight='bold')

    plt.ylim(0, 4000)
    if saveAs:
        plt.savefig(saveAs)
    plt.show()


def filterSourceCatalog(srcCat, brightSourceFraction, minSources=15):
    maxFlux = np.nanmax(srcCat['base_CircularApertureFlux_3_0_instFlux'])
    selectBrightestSource = srcCat['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.99
    brightestSource = srcCat.subset(selectBrightestSource)
    brightestCentroid = (brightestSource['base_SdssCentroid_x'][0],
                         brightestSource['base_SdssCentroid_y'][0])
    filteredCatalog = srcCat.subset(srcCat['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.001)

    print(f"Found {len(srcCat)} sources, {len(filteredCatalog)} bright sources")
    print(f"Brightest centroid at {brightestCentroid}")

    if not filteredCatalog.isContiguous():
        filteredCatalog = filteredCatalog.copy(deep=True)
    sources = filteredCatalog.asAstropy()
    sources.keep_columns(['base_SdssCentroid_x', 'base_SdssCentroid_y',
                          'base_CircularApertureFlux_3_0_instFlux'])
    sources.sort('base_CircularApertureFlux_3_0_instFlux', reverse=True)

    nSources = len(sources)
    if nSources <= minSources:
        return sources

    startPos = 0
    endPos = math.ceil(nSources*brightSourceFraction)
    sources.remove_rows([i for i in range(startPos, endPos+1)])
    return sources


def blindSolve(exp, *, radiusInDegrees=1, brightSourceFraction=0.1, doPlot=False, savePlotAs=None):
    imCharResult = runImchar(exp)

    sourceCatalog = imCharResult.sourceCat
    filteredSources = filterSourceCatalog(sourceCatalog, brightSourceFraction)

    if doPlot:
        plot(exp, sourceCatalog, saveAs=savePlotAs)

    ast = AstrometryNet()
    ast.api_key = getApiKey()

    vi = exp.getInfo().getVisitInfo()
    ra, dec = vi.boresightRaDec

    image_width = 4072
    image_height = 4000
    scale_units = 'arcsecperpix'
    scale_type = 'ev'  # ev means submit estimate and % error
    scale_est = 0.095
    scale_err = 3.0  # error as percentage
    center_ra = ra.asDegrees()
    center_dec = dec.asDegrees()
    radius = radiusInDegrees
    wcs_header = ast.solve_from_source_list(filteredSources['base_SdssCentroid_x'],
                                            filteredSources['base_SdssCentroid_y'],
                                            image_width, image_height, scale_units=scale_units,
                                            scale_type=scale_type, scale_est=scale_est, scale_err=scale_err,
                                            center_ra=center_ra, center_dec=center_dec, radius=radius,
                                            solve_timeout=240)

    nominalRa, nominalDec = exp.getInfo().getVisitInfo().getBoresightRaDec()

    calculatedRa = geom.Angle(wcs_header['CRVAL1'], geom.degrees)
    calculatedDec = geom.Angle(wcs_header['CRVAL2'], geom.degrees)

    deltaRa = geom.Angle(wcs_header['CRVAL1'] - nominalRa.asDegrees(), geom.degrees)
    deltaDec = geom.Angle(wcs_header['CRVAL2'] - nominalDec.asDegrees(), geom.degrees)

    ret = {'nominalRa': nominalRa,
           'nominalDec': nominalDec,
           'calculatedRa': calculatedRa,
           'calculatedDec': calculatedDec,
           'deltaRa': deltaRa,
           'deltaDec': deltaDec,
           'deltaRaArcsec': deltaRa.asArcseconds(),
           'deltaDecArcsec': deltaDec.asArcseconds(),
           'astrometry_net_wcs_header': wcs_header}

    return ret
