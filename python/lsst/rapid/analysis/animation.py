import os
import subprocess
import shutil
import uuid
import math

import matplotlib.pyplot as plt

from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask, QuickFrameMeasurementTaskConfig
from lsst.atmospec.processStar import getTargetCentroidFromWcs
import lsst.afw.display as afwDisplay
import lsst.afw.math as afwMath
import lsst.daf.persistence as dafPersist
import lsst.log
import lsst.meas.algorithms as measAlg

from lsst.atmospec.utils import airMassFromRawMetadata
logger = lsst.log.Log.getLogger("lsst.rapid.analysis.animation")


class Animator():
    """Animate the list of dataIds in the order in which they are specified
    for the data product specified."""

    def __init__(self, butler, dataIdList, outputPath, outputFilename, *,
                 remakePngs=False,
                 clobberVideoAndGif=False,
                 keepIntermediateGif=False,
                 smoothImages=True,
                 plotObjectCentroids=True,
                 useQfmForCentroids=False,
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
        self.keepIntermediateGif = keepIntermediateGif
        self.smoothImages = smoothImages
        self.plotObjectCentroids = plotObjectCentroids
        self.useQfmForCentroids = useQfmForCentroids
        self.dataProcuctToPlot = dataProcuctToPlot
        self.ffMpegBinary = ffMpegBinary
        self.debug = debug

        # zfilled at the start as animation is alphabetical
        # if you're doing more than 1e6 files you've got bigger problems
        self.toAnimateTemplate = "%06d-%s-%s.png"
        self.basicTemplate = "%s-%s.png"

        qfmTaskConfig = QuickFrameMeasurementTaskConfig()
        self.qfmTask = QuickFrameMeasurementTask(config=qfmTaskConfig)

        afwDisplay.setDefaultBackend("matplotlib")
        self.fig = plt.figure(figsize=(15, 15))
        self.disp = afwDisplay.Display(self.fig)
        self.disp.setImageColormap('gray')
        self.disp.scale('asinh', 'zscale')

        self.pngsToMakeDataIds = []
        self.preRun()  # sets the above list

    @staticmethod
    def _strDataId(dataId):
        if 'dayObs' in dataId and 'seqNum' in dataId:  # nicely ordered if easy
            return f"{dataId['dayObs']}-{dataId['seqNum']:05d}"

        # General case (and yeah, I should probably learn regex someday)
        dIdStr = str(dataId)
        dIdStr = dIdStr.replace(' ', "")
        dIdStr = dIdStr.replace('{', "")
        dIdStr = dIdStr.replace('}', "")
        dIdStr = dIdStr.replace('\'', "")
        dIdStr = dIdStr.replace(':', "-")
        dIdStr = dIdStr.replace(',', "-")
        return dIdStr

    def dataIdToFilename(self, dataId, includeNumber=False, imNum=None):
        """Convert dataId to filename.

        Returns a full path+filename by default. if includeNumber then
        returns just the filename for use in temporary dir for animation."""
        if includeNumber:
            assert imNum is not None

        dIdStr = self._strDataId(dataId)

        if includeNumber:  # for use in temp dir, so not full path
            filename = self.toAnimateTemplate%(imNum, dIdStr, self.dataProcuctToPlot)
            return os.path.join(filename)
        else:
            filename = self.basicTemplate%(dIdStr, self.dataProcuctToPlot)
            return os.path.join(self.pngPath, filename)

    def exists(self, obj):
        if type(obj) == str:
            return os.path.exists(obj)
        raise RuntimeError("Other type checks not yet implemented")

    def preRun(self):
        # check the binary is there and the paths work
        assert os.path.exists(self.ffMpegBinary), "Cannot find ffmpeg binary for animation"
        if not os.path.exists(self.pngPath):
            os.makedirs(self.pngPath)
        assert os.path.exists(self.pngPath), f"Failed to create output dir: {self.pngsPath}"

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

        logger.info(f"Of the provided {len(self.dataIdList)} dataIds:")
        logger.info(f"{len(dIdsWithPngs)} existing pngs were found")
        logger.info(f"{len(dIdsWithoutPngs)} do not yet exist")

        if missingData:
            for dId in missingData:
                msg = f"Failed to find {self.dataProcuctToPlot} for {dId}"
                logger.warn(msg)
                self.dataIdList.remove(dId)
            logger.info(f"Of the {len(dIdsWithoutPngs)} dataIds without pngs, {len(missingData)}" +
                        " did not have the corresponding dataset existing")

        if self.remakePngs:
            self.pngsToMakeDataIds = [d for d in self.dataIdList if d not in missingData]
        else:
            self.pngsToMakeDataIds = [d for d in dIdsWithoutPngs if d not in missingData]

        msg = f"So {len(self.pngsToMakeDataIds)} will be made"
        if self.remakePngs and len(dIdsWithPngs) > 0:
            msg += " because remakePngs=True"
        logger.info(msg)

    def run(self):
        # make the missing pngs
        if self.pngsToMakeDataIds:
            logger.info('Creating necessary pngs...')
            for i, dataId in enumerate(self.pngsToMakeDataIds):
                print(f'Making png for file {i+1} of {len(self.pngsToMakeDataIds)}')
                self.makePng(dataId, self.dataIdToFilename(dataId))

        # stage files in temp dir with numbers prepended to filenames
        logger.info('Copying files to ordered temp dir...')
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

        # # create gif in temp dir
        # outputGifFilename = os.path.join(tempDir, 'animation.gif')
        # self.pngsToGif(pngFileList, outputGifFilename)

        # # gif turn into mp4, optionally keep gif by moving up to output dir
        # logger.info('Turning gif into mp4...')
        # outputMp4Filename = self.outputFilename
        # self.gifToMp4(outputGifFilename, outputMp4Filename)

        # self.tidyUp(tempDir)
        # logger.info('Finished!')

        # create gif in temp dir

        logger.info('Making mp4 of pngs...')
        self.pngsToMp4(tempDir, self.outputFilename, 10, verbose=False)
        self.tidyUp(tempDir)
        logger.info(f'Finished! Output at {self.outputFilename}')

    def _titleFromExp(self, exp, dataId):
        items = ["OBJECT", "expTime", "FILTER", "imageType"]
        obj, expTime, filterCompound, imageType = self.butler.queryMetadata('raw', items, **dataId)[0]
        filt, grating = filterCompound.split('~')
        rawMd = self.butler.get('raw_md', **dataId)
        airmass = airMassFromRawMetadata(rawMd)

        title = f"{dataId['dayObs']} - seqNum {dataId['seqNum']} - "
        title += f"Object: {obj} expTime: {expTime}s Filter: {filt} Grating: {grating} Airmass: {airmass:.3f}"
        return title

    def getStarPixCoord(self, exp, doMotionCorrection=True, useQfm=False):
        target = exp.getMetadata()['OBJECT']

        if self.useQfmForCentroids:
            try:
                result = self.qfmTask.run(exp)
                pixCoord = result.brightestObjCentroid
                expId = exp.getInfo().getVisitInfo().getExposureId()
                print(f'XXX expId {expId} centroid {pixCoord}')
            except Exception:
                return None
        else:
            pixCoord = getTargetCentroidFromWcs(exp, target, doMotionCorrection=doMotionCorrection)
        return pixCoord

    def makePng(self, dataId, saveFilename):
        if self.exists(saveFilename) and not self.remakePngs:  # should not be possible due to prerun
            assert False, f"Almost overwrote {saveFilename} - how is this possible?"

        if self.debug:
            print(f"Creating {saveFilename}")

        self.disp.erase()
        self.fig.clear()

        # must always keep exp unsmoothed for the centroiding via qfm
        try:
            exp = self.butler.get(self.dataProcuctToPlot, **dataId)
        except Exception:
            # oh no, that should never happen, but it does! Let's just skip
            print(f'Skipped {dataId}, because {self.dataProcuctToPlot} retrieval failed!')
            return
        toDisplay = exp
        if self.smoothImages:
            toDisplay = exp.clone()
            toDisplay = self._smoothExp(toDisplay, 2)

        try:
            self.disp.mtv(toDisplay.image, title=self._titleFromExp(exp, dataId))
            self.disp.scale('asinh', 'zscale')
        except RuntimeError:  # all-nan images slip through and don't display
            self.disp.scale('linear', 0, 1)
            self.disp.mtv(toDisplay.image, title=self._titleFromExp(exp, dataId))
            self.disp.scale('asinh', 'zscale')  # set back for next image
            pass

        if self.plotObjectCentroids:
            try:
                pixCoord = self.getStarPixCoord(exp)
                if pixCoord:
                    self.disp.dot('x', *pixCoord, ctype='C1', size=50)
                    self.disp.dot('o', *pixCoord, ctype='C1', size=50)
                else:
                    self.disp.dot('x', 2000, 2000, ctype='red', size=2000)
            except Exception:
                logger.warn(f"Failed to find OBJECT location for {dataId}")

        deltaH = -0.05
        deltaV = -0.05
        plt.subplots_adjust(right=1+deltaH, left=0-deltaH, top=1+deltaV, bottom=0-deltaV)
        self.fig.savefig(saveFilename)

    def pngsToMp4(self, indir, outfile, framerate, verbose=False):
        """Create the movie with ffmpeg, from files."""
        # NOTE: the order of ffmpeg arguments *REALLY MATTERS*.
        # Reorder them at your own peril!
        pathPattern = f'\"{os.path.join(indir, "*.png")}\"'
        if verbose:
            ffmpeg_verbose = 'info'
        else:
            ffmpeg_verbose = 'error'
        cmd = ['ffmpeg',
               '-v', ffmpeg_verbose,
               '-f', 'image2',
               '-y',
               '-pattern_type glob',
               '-framerate', f'{framerate}',
               '-i', pathPattern,
               '-vcodec', 'libx264',
               '-b:v', '20000k',
               '-profile:v', 'main',
               '-pix_fmt', 'yuv420p',
               '-threads', '10',
               '-r', f'{framerate}',
               os.path.join(outfile)]

        subprocess.check_call(r' '.join(cmd), shell=True)

    def tidyUp(self, tempDir):
        shutil.rmtree(tempDir)
        return

    def _smoothExp(self, exp, smoothing, kernelSize=7):
        """Use for DISPLAY ONLY!

        Return a smoothed copy of the exposure
        with the original mask plane in place."""
        psf = measAlg.DoubleGaussianPsf(kernelSize, kernelSize, smoothing/(2*math.sqrt(2*math.log(2))))
        newExp = exp.clone()
        originalMask = exp.mask

        kernel = psf.getKernel()
        afwMath.convolve(newExp.maskedImage, newExp.maskedImage, kernel, afwMath.ConvolutionControl())
        newExp.mask = originalMask
        return newExp


if __name__ == '__main__':
    dataProcuctToPlot = 'quickLookExp'
    repoPath = '/project/shared/auxTel/rerun/quickLook'
    outputPath = '/home/mfl/animatorOutput/feb18Debug/'
    outputFilename = 'allFixed.mp4'

    butler = dafPersist.Butler(repoPath)

    skipTypes = ['BIAS', 'DARK', 'FLAT']

    def isOnSky(dataId):
        if dataId['imageType'] not in skipTypes:
            return True
        return False

    # def isDispersed(dataId):
    #     if dataId['imageType'] not in skipTypes:
    #         return True
    #     return False

    if False:
        days = ['2020-02-17', '2020-02-18', '2020-02-19', '2020-02-20', '2020-02-21',
                '2020-03-12', '2020-03-13', '2020-03-14', '2020-03-15', '2020-03-16']
        dataIds = []
        for dayObs in days:
            data = butler.queryMetadata('raw', ['seqNum', 'filter', 'imageType'], dayObs=dayObs)
            for (seqNum, filterCompound, imageType) in data:
                filt = filterCompound.split("~")[0]
                grating = filterCompound.split("~")[1]
                dataIds.append({'dayObs': dayObs,
                                'seqNum': seqNum,
                                'filter': filt,
                                'grating': grating,
                                'imageType': imageType})

        scienceIds = [x for x in filter(isOnSky, dataIds)]

        print(f'{len(dataIds)} dataIds total')
        print(f'{len(scienceIds)} science ids')

        fitted = []
        for dataId in scienceIds:
            if butler.datasetExists(dataProcuctToPlot, dayObs=dataId['dayObs'], seqNum=dataId['seqNum']):
                fitted.append({'dayObs': dataId['dayObs'], 'seqNum': dataId['seqNum']})
        print(f'{len(fitted)} with {dataProcuctToPlot}')

        animator = Animator(butler, fitted, outputPath, outputFilename,
                            dataProcuctToPlot=dataProcuctToPlot,
                            remakePngs=False,
                            debug=False,
                            clobberVideoAndGif=True,
                            plotObjectCentroids=False)
        animator.run()
    else:
        # skipTypes = []
        dayObs = '2021-02-17'
        seqNums = range(100, 475+1)

        toUse = []

        for seqNum in seqNums:
            imageType = butler.queryMetadata('raw', 'imageType', dayObs=dayObs, seqNum=seqNum)[0]
            if imageType in skipTypes:
                print(f"Skipping {seqNum} because it is a {imageType}")
                continue
            if butler.datasetExists(dataProcuctToPlot, dayObs=dayObs, seqNum=seqNum):
                toUse.append({'dayObs': dayObs, 'seqNum': seqNum})

        print(f'{len(toUse)} with {dataProcuctToPlot}')

        animator = Animator(butler, toUse, outputPath, outputFilename,
                            dataProcuctToPlot=dataProcuctToPlot,
                            remakePngs=False,
                            debug=False,
                            clobberVideoAndGif=True,
                            plotObjectCentroids=True,
                            useQfmForCentroids=True)
        animator.run()
