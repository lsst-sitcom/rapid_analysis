import os
import subprocess
import shutil
import uuid
import math

from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt

import lsst.afw.display as afwDisplay
import lsst.afw.math as afwMath
import lsst.daf.persistence as dafPersist
import lsst.log
import lsst.meas.algorithms as measAlg
import lsst.geom as geom

logger = lsst.log.Log.getLogger("lsst.rapid.analysis.animation")


class Animator():
    """Animate the list of dataIds in the order in which they are specified
    for the data product specified."""

    def __init__(self, butler, dataIdList, outputPath, outputFilename, *,
                 remakePngs=False,
                 clobberVideoAndGif=False,
                 smoothImages=True,
                 plotObjectCentroids=True,
                 dataProcuctToPlot='calexp',
                 ffMpegBinary='/home/mfl/bin/ffmpeg',
                 debug=False):

        self.butler = butler
        self.dataIdList = dataIdList
        self.outputPath = outputPath
        self.outputFilename = os.path.join(outputPath, outputFilename)
        if not self.outputFilename.endswith(".mp4"):
            self.outputFilename += ".mp4"
        self.pngPath = os.path.join(outputPath, "pngs/")

        self.remakePngs = remakePngs
        self.clobberVideoAndGif = clobberVideoAndGif
        self.smoothImages = smoothImages
        self.plotObjectCentroids = plotObjectCentroids
        self.dataProcuctToPlot = dataProcuctToPlot
        self.ffMpegBinary = ffMpegBinary
        self.debug = debug

        # zfilled at the start as animation is alphabetical
        # if you're doing more than 1e6 files you've got bigger problems
        self.toAnimateTemplate = "%06d-%s.png"
        self.basicTemplate = "%s.png"

        afwDisplay.setDefaultBackend("matplotlib")
        self.fig = plt.figure(figsize=(15, 15))
        self.disp = afwDisplay.Display(self.fig)
        self.disp.setImageColormap('gray')
        self.disp.scale('asinh', 'zscale')

        self.pngsToMakeDataIds = []
        self.preRun()  # sets the above list

    def dataIdToFilename(self, dataId, includeNumber=False, imNum=None):
        """Convert dataId to filename.

        Returns a full path+filename by default. if includeNumber then
        returns just the filename for use in temporary dir for animation."""
        if includeNumber:
            assert imNum is not None
        # Yeah, I should probably learn regex someday
        dIdStr = str(dataId)
        dIdStr = dIdStr.replace(' ', "")
        dIdStr = dIdStr.replace('{', "")
        dIdStr = dIdStr.replace('}', "")
        dIdStr = dIdStr.replace('\'', "")
        dIdStr = dIdStr.replace(':', "-")
        dIdStr = dIdStr.replace(',', "-")
        if includeNumber:  # for use in temp dir, so not full path
            filename = self.toAnimateTemplate%(imNum, dIdStr)
            return os.path.join(filename)
        else:
            filename = self.basicTemplate%(dIdStr)
            return os.path.join(self.pngPath, filename)

    def exists(self, obj):
        if type(obj) == str:
            return os.path.exists(obj)
        raise RuntimeError("Other type checks not yet implemented")

    def preRun(self):
        # check the binary is there and the paths work
        assert os.path.exists(self.ffMpegBinary), "Cannot find ffMeg binary for animation"
        if not os.path.exists(self.pngPath):
            os.makedirs(self.pngPath)

        if self.exists(self.outputFilename):
            if self.clobberVideoAndGif:
                os.remove(self.outputFilename)
            else:
                raise RuntimeError(f"Output file {self.outputFilename} exists and clobber==False")

        # make list of found & missing files
        dIdsWithPngs = [d for d in self.dataIdList if self.exists(self.dataIdToFilename(d))]
        dIdsWithoutPngs = [d for d in self.dataIdList if d not in dIdsWithPngs]
        if self.debug:
            print(f"dIdsWithPngs = {dIdsWithPngs}")
            print(f"dIdsWithoutPngs = {dIdsWithoutPngs}")

        # check the datasets exist for the pngs which need remaking
        missingData = [d for d in dIdsWithoutPngs if not self.butler.datasetExists(self.dataProcuctToPlot,
                                                                                   **d)]
        if missingData:
            for dId in missingData:
                msg = f"Failed to find {self.dataProcuctToPlot} for {dId}"
                logger.warn(msg)
                self.dataIdList.remove(dId)

        if self.remakePngs:
            self.pngsToMakeDataIds = [d for d in self.dataIdList if d not in missingData]
        else:
            self.pngsToMakeDataIds = [d for d in dIdsWithoutPngs if d not in missingData]

    def run(self):
        # make the missing pngs
        if self.pngsToMakeDataIds:
            logger.info(f'Creating necessary pngs...')
            for i, dataId in enumerate(self.pngsToMakeDataIds):
                print(f'Making png for file {i+1} of {len(self.pngsToMakeDataIds)}')
                self.makePng(dataId, self.dataIdToFilename(dataId))

        # stage files in temp dir with numbers prepended to filenames
        logger.info(f'Making gif of pngs...')
        pngFilesOriginal = [self.dataIdToFilename(d) for d in self.dataIdList]
        for filename in pngFilesOriginal:  # these must all now exist, but let's assert just in case
            assert self.exists(filename)
        tempDir = os.path.join(self.pngPath, f"{uuid.uuid1()}/"[0:8])
        os.makedirs(tempDir)
        pngFileList = []  # list of number-prepended files in the temp dir
        for i, dId in enumerate(self.dataIdList):
            srcFile = self.dataIdToFilename(dId)
            destFile = os.path.join(tempDir, self.dataIdToFilename(dId, includeNumber=True, imNum=i))
            shutil.copy(srcFile, destFile)
            pngFileList.append(destFile)

        # create gif in temp dir
        outputGifFilename = os.path.join(tempDir, 'animation.gif')
        self.pngsToGif(pngFileList, outputGifFilename)

        # gif turn into mp4, optionally keep gif by moving up to output dir
        logger.info(f'Turning gif into mp4...')
        outputMp4Filename = self.outputFilename
        self.gifToMp4(outputGifFilename, outputMp4Filename, copyGifToOutdir=True)

        self.tidyUp(tempDir)
        logger.info(f'Finished!')

    def _titleFromExp(self, exp, dataId):
        def _airMassFromrRawMd(md):
            auxTelLocation = EarthLocation(lat=-30.244639*u.deg, lon=-70.749417*u.deg, height=2663*u.m)
            time = Time(md['DATE-OBS'])
            skyLocation = SkyCoord(md['RASTART'], md['DECSTART'], unit=u.deg)
            altAz = AltAz(obstime=time, location=auxTelLocation)
            observationAltAz = skyLocation.transform_to(altAz)
            return observationAltAz.secz.value

        items = ["OBJECT", "expTime", "FILTER", "imageType"]
        obj, expTime, filterCompound, imageType = self.butler.queryMetadata('raw', items, **dataId)[0]
        filt, grating = filterCompound.split('~')
        rawMd = self.butler.get('raw_md', **dataId)
        airmass = _airMassFromrRawMd(rawMd)

        title = f"Object: {obj} expTime: {expTime}s Filter: {filt} Grating: {grating} Airmass: {airmass:.3f}"
        return title

    def getStarPixCoord(self, exp):
        from astroquery.simbad import Simbad
        target = exp.getMetadata()['OBJECT']
        obj = Simbad.query_object(target)
        assert len(obj) == 1

        raStr = obj[0]['RA']
        decStr = obj[0]['DEC']
        skyLocation = SkyCoord(raStr, decStr, unit=(u.hourangle, u.degree), frame='icrs')
        raRad, decRad = skyLocation.ra.rad, skyLocation.dec.rad
        ra = geom.Angle(raRad)
        dec = geom.Angle(decRad)
        targetLocation = geom.SpherePoint(ra, dec)

        pixCoord = exp.getWcs().skyToPixel(targetLocation)
        return pixCoord

    def makePng(self, dataId, saveFilename):
        if self.exists(saveFilename) and not self.remakePngs:  # should not be possible due to prerun
            assert False

        if self.debug:
            print(f"Creating {saveFilename}")

        exp = self.butler.get(self.dataProcuctToPlot, **dataId)
        if self.smoothImages:
            exp = self._smoothExp(exp, 2)
        self.disp.mtv(exp.image, title=self._titleFromExp(exp, dataId))
        if self.plotObjectCentroids:
            pixCoord = self.getStarPixCoord(exp)
            # self.disp.dot('o', *pixCoord, ctype='C1', size=50)
            self.disp.dot('x', *pixCoord, ctype='C1', size=50)

        self.fig.savefig(saveFilename)

    def pngsToGif(self, fileList, outputGifFilename):
        subprocess.run(['convert', '-delay', '10', '-loop', '0', *fileList, outputGifFilename], check=True)

    def gifToMp4(self, inputGifFilename, outputMp4Filename, copyGifToOutdir=True):
        command = (f'{self.ffMpegBinary} -i {inputGifFilename} -pix_fmt yuv420p'
                   f' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" {outputMp4Filename}')
        output, error = subprocess.Popen(command, universal_newlines=True, shell=True,
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        if copyGifToOutdir:
            outputGifName = self.outputFilename.replace('.mp4', '.gif')
            shutil.copy(inputGifFilename, outputGifName)

    def tidyUp(self, tempDir):
        shutil.rmtree(tempDir)
        return

    def _smoothExp(self, exp, smoothing, kernelSize=7):
        """Use for DISPLAY ONLY!
        Return a smoothed copy of the exposure with the original mask plane in place."""
        psf = measAlg.DoubleGaussianPsf(kernelSize, kernelSize, smoothing/(2*math.sqrt(2*math.log(2))))
        newExp = exp.clone()
        originalMask = exp.mask

        kernel = psf.getKernel()
        afwMath.convolve(newExp.maskedImage, newExp.maskedImage, kernel, afwMath.ConvolutionControl())
        newExp.mask = originalMask
        return newExp


if __name__ == '__main__':
    butler = dafPersist.Butler('/project/shared/auxTel/rerun/mfl/preprocessing/')
    dataIds = [{'dayObs': '2020-03-15', 'seqNum': 161},
               {'dayObs': '2020-03-15', 'seqNum': 162},
               {'dayObs': '2020-03-15', 'seqNum': 163},
               {'dayObs': '2020-03-15', 'seqNum': 164},
               {'dayObs': '2020-03-15', 'seqNum': 164000}]

    dataId = dataIds[0]

    pathToPngs = '/home/mfl/animatorTest'
    animator = Animator(butler, dataIds, pathToPngs, 'animation.mp4',
                        remakePngs=False, debug=False, clobberVideoAndGif=True)
    animator.run()
