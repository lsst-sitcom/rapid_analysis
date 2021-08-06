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
from time import sleep
from pathlib import Path
import shutil
import tempfile
import logging
import matplotlib.pyplot as plt

try:
    from google.cloud import storage
    HAS_GOOGLE_STORAGE = True
except ImportError:
    HAS_GOOGLE_STORAGE = False

try:
    from lsst_efd_client import EfdClient
    HAS_EFD_CLIENT = True
except ImportError:
    HAS_EFD_CLIENT = False

import lsst.daf.persistence as dafPersist
from lsst.pex.exceptions import NotFoundError
from lsst.rapid.analysis.bestEffort import BestEffortIsr
from lsst.rapid.analysis.imageExaminer import ImageExaminer
from lsst.rapid.analysis.summarizeImage import SummarizeImage
from lsst.rapid.analysis.mountTorques import plotMountTracking
from lsst.rapid.analysis.monitorPlotting import plotExp
from lsst.atmospec.utils import isDispersedDataId

CHANNELS = ["summit_imexam", "summit_specexam", "auxtel_mount_torques",
            "auxtel_monitor"]

PREFIXES = {chan: chan.replace('_', '-') for chan in CHANNELS}


def _dataIdToFilename(channel, dataId):
    filename = f"{PREFIXES[channel]}_dayObs_{dataId['dayObs']}_seqNum_{dataId['seqNum']}.png"
    return filename


def _waitForDataProduct(butler, dataProduct, dataId, logger, maxTime=20):
    cadence = 0.25
    maxLoops = int(maxTime//cadence)
    for retry in range(maxLoops):
        if butler.datasetExists(dataProduct, dataId):
            return butler.get(dataProduct, dataId)
        else:
            sleep(cadence)
    logger.warn(f'Waited {maxTime}s for {dataProduct} for {dataId} to no avail')
    return None


class Uploader():
    def __init__(self):
        if not HAS_GOOGLE_STORAGE:
            from lsst.rapid.analysis.utils import GOOGLE_CLOUD_MISSING_MSG
            raise RuntimeError(GOOGLE_CLOUD_MISSING_MSG)
        self.client = storage.Client()
        self.bucket = self.client.get_bucket('rubintv_data')
        self.log = logging.getLogger("googleUploader")

    def googleUpload(self, channel, sourceFilename, uploadAsFilename=None):
        if channel not in CHANNELS:
            self.log.warn(f"Error: {channel} not in {CHANNELS}")
            return

        if uploadAsFilename and (os.path.basename(sourceFilename) != uploadAsFilename):
            finalName = os.path.join(os.path.dirname(sourceFilename), uploadAsFilename)
            shutil.move(sourceFilename, finalName)
            assert os.path.exists(finalName)
        else:
            finalName = sourceFilename

        path = Path(finalName)
        blob = self.bucket.blob("/".join([channel, path.name]))
        self.log.info(f'Uploaded {sourceFilename} to {finalName}')
        blob.upload_from_filename(finalName)


class Watcher():
    cadence = 1  # in seconds

    def __init__(self, repoDir, dataProduct, **kwargs):
        """Watch for data products in a repo.
        """
        self.repoDir = repoDir
        self.butler = dafPersist.Butler(repoDir)
        self.dataProduct = dataProduct
        self.log = logging.getLogger("watcher")

    def _getLatestExpId(self):
        return sorted(self.butler.queryMetadata(self.dataProduct, 'expId'))[-1]

    def _getDayObsSeqNumFromExpId(self, expId):
        return self.butler.queryMetadata(self.dataProduct, ['dayObs', 'seqNum'], expId=expId)[0]

    def _getLatestImageDataIdAndExpId(self):
        expId = self._getLatestExpId()
        dayObs, seqNum = self._getDayObsSeqNumFromExpId(expId)
        dataId = {'dayObs': dayObs, 'seqNum': seqNum}
        return dataId, expId

    def run(self, callback, durationInSeconds=-1):
        """Run callback(dataId) each time a new image lands
        """

        if durationInSeconds == -1:
            nLoops = int(1e9)
        else:
            nLoops = int(durationInSeconds//self.cadence)

        lastFound = -1
        for i in range(nLoops):
            try:
                dataId, expId = self._getLatestImageDataIdAndExpId()

                if lastFound == expId:
                    sleep(self.cadence)
                    continue
                else:
                    callback(dataId)

                lastFound = expId

            except NotFoundError as e:  # NotFoundError when filters aren't defined
                print(f'Skipped displaying {dataId} due to {e}')
        return


class IsrRunner():

    def __init__(self, repoDir, **kwargs):
        self.watcher = Watcher(repoDir, 'raw')
        self.bestEffort = BestEffortIsr(repoDir, **kwargs)
        outpath = os.path.join(repoDir, 'rerun/quickLook')
        self.butler = dafPersist.Butler(outpath)
        self.log = logging.getLogger("isrRunner")

    def callback(self, dataId):
        quickLookExp = self.bestEffort.getExposure(dataId)
        self.butler.put(quickLookExp, 'quickLookExp', dataId)
        self.log.info(f'Put quickLookExp for {dataId}, awaiting next image...')

    def run(self):
        self.watcher.run(self.callback)


class ImExaminer():

    def __init__(self, repoDir):
        self.dataProduct = 'quickLookExp'
        self.watcher = Watcher(repoDir, self.dataProduct)
        self.uploader = Uploader()
        self.repoDir = repoDir
        self.butler = dafPersist.Butler(repoDir)
        self.log = logging.getLogger("imExaminer")
        self.channel = 'summit_imexam'

    def _imExamine(self, exp, dataId, outputFilename):
        if os.path.exists(outputFilename):  # unnecessary now we're using tmpfile
            self.log.warn(f"Skipping {outputFilename}")
            return
        imexam = ImageExaminer(exp, savePlots=outputFilename, doTweakCentroid=True)
        imexam.plot()

    def callback(self, dataId):
        try:
            self.log.info(f'Running imexam on {dataId}')
            tempFilename = tempfile.mktemp(suffix='.png')
            uploadFilename = _dataIdToFilename(self.channel, dataId)
            exp = _waitForDataProduct(self.butler, self.dataProduct, dataId, self.log)
            if not exp:
                raise RuntimeError(f'Failed to get {self.dataProduct} for {dataId}')
            self._imExamine(exp, dataId, tempFilename)

            self.log.info("Uploading imExam to storage bucket")
            self.uploader.googleUpload(self.channel, tempFilename, uploadFilename)
            self.log.info('Upload complete')

        except Exception as e:
            self.log.warn(f"Skipped imExam on {dataId} because {e}")
            return None

    def run(self):
        self.watcher.run(self.callback)


class SpecExaminer():

    def __init__(self, repoDir):
        self.dataProduct = 'quickLookExp'
        self.watcher = Watcher(repoDir, self.dataProduct)
        self.uploader = Uploader()
        self.repoDir = repoDir
        self.butler = dafPersist.Butler(repoDir)
        self.log = logging.getLogger("specExaminer")
        self.channel = 'summit_specexam'

    def _specExamine(self, exp, dataId, outputFilename):
        if os.path.exists(outputFilename):  # unnecessary now we're using tmpfile?
            self.log.warn(f"Skipping {outputFilename}")
            return
        summary = SummarizeImage(exp, savePlotAs=outputFilename)
        summary.run()

    def callback(self, dataId):
        try:
            if not isDispersedDataId(dataId, self.butler):
                self.log.info(f'Skipping non dispersed image {dataId}')
                return

            self.log.info(f'Running specExam on {dataId}')
            tempFilename = tempfile.mktemp(suffix='.png')
            uploadFilename = _dataIdToFilename(self.channel, dataId)
            exp = _waitForDataProduct(self.butler, self.dataProduct, dataId, self.log)
            if not exp:
                raise RuntimeError(f'Failed to get {self.dataProduct} for {dataId}')
            self._specExamine(exp, dataId, tempFilename)

            self.log.info("Uploading specExam to storage bucket")
            self.uploader.googleUpload(self.channel, tempFilename, uploadFilename)
            self.log.info('Upload complete')

        except Exception as e:
            self.log.info(f"Skipped imExam on {dataId} because {e}")
            return None

    def run(self):
        self.watcher.run(self.callback)


class Monitor():
    def __init__(self, repoDir):
        self.dataProduct = 'quickLookExp'
        self.watcher = Watcher(repoDir, self.dataProduct)
        self.uploader = Uploader()
        self.repoDir = repoDir
        self.butler = dafPersist.Butler(repoDir)
        self.log = logging.getLogger("monitor")
        self.channel = 'auxtel_monitor'
        self.fig = plt.figure(figsize=(12, 12))

    def _plotImage(self, exp, dataId, outputFilename):
        if os.path.exists(outputFilename):  # unnecessary now we're using tmpfile
            self.log.warn(f"Skipping {outputFilename}")
            return
        plotExp(exp, dataId, self.fig, outputFilename)

    def callback(self, dataId):
        try:
            self.log.info(f'Generating monitor image for {dataId}')
            tempFilename = tempfile.mktemp(suffix='.png')
            uploadFilename = _dataIdToFilename(self.channel, dataId)
            exp = _waitForDataProduct(self.butler, self.dataProduct, dataId, self.log)
            if not exp:
                raise RuntimeError(f'Failed to get {self.dataProduct} for {dataId}')
            self._plotImage(exp, dataId, tempFilename)

            self.log.info("Uploading monitor image to storage bucket")
            self.uploader.googleUpload(self.channel, tempFilename, uploadFilename)
            self.log.info('Upload complete')

        except Exception as e:
            self.log.warn(f"Skipped monitor image for {dataId} because {e}")
            return None

    def run(self):
        self.watcher.run(self.callback)


class MountTorquePlotter():

    def __init__(self, repoDir):
        if not HAS_EFD_CLIENT:
            from lsst.rapid.analysis.utils import EFD_CLIENT_MISSING_MSG
            raise RuntimeError(EFD_CLIENT_MISSING_MSG)
        self.dataProduct = 'raw'
        self.watcher = Watcher(repoDir, self.dataProduct)
        self.uploader = Uploader()
        self.repoDir = repoDir
        self.butler = dafPersist.Butler(repoDir)
        self.client = EfdClient('summit_efd')
        self.log = logging.getLogger("mountTorquePlotter")
        self.channel = 'auxtel_mount_torques'
        self.fig = plt.figure(figsize=(16, 16))

    def callback(self, dataId):
        try:
            tempFilename = tempfile.mktemp(suffix='.png')
            uploadFilename = _dataIdToFilename(self.channel, dataId)
            plotted = plotMountTracking(dataId, self.butler, self.client, self.fig, tempFilename, 2, self.log)

            if plotted:  # skips many image types and short exps
                self.log.info("Uploading mount torque plot to storage bucket")
                self.uploader.googleUpload(self.channel, tempFilename, uploadFilename)
                self.log.info('Upload complete')

        except Exception as e:
            self.log.warn(f"Skipped creating mount plots for {dataId} because {e}")

    def run(self):
        self.watcher.run(self.callback)
