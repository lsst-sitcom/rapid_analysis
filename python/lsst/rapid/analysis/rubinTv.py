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

from lsst.pex.exceptions import NotFoundError
from lsst.rapid.analysis.bestEffort import BestEffortIsr
from lsst.rapid.analysis.imageExaminer import ImageExaminer
from lsst.rapid.analysis.spectrumExaminer import SpectrumExaminer
from lsst.rapid.analysis.mountTorques import plotMountTracking
from lsst.rapid.analysis.monitorPlotting import plotExp
from lsst.rapid.analysis.butlerUtils import (makeDefaultLatissButler, LATISS_REPO_LOCATION_MAP, datasetExists,
                                             getMostRecentDataId, getExpIdFromDayObsSeqNum)
from lsst.rapid.analysis.utils import dayObsIntToString
from lsst.atmospec.utils import isDispersedDataId

CHANNELS = ["summit_imexam", "summit_specexam", "auxtel_mount_torques",
            "auxtel_monitor"]

PREFIXES = {chan: chan.replace('_', '-') for chan in CHANNELS}


def _dataIdToFilename(channel, dataId):
    """Convert a dataId to a png filename.

    Parameters
    ----------
    channel : `str`
        The name of the RubinTV channel
    dataId : `dict`
        The dataId

    Returns
    -------
    filename : `str`
        The filename
    """
    dayObsStr = dayObsIntToString(dataId['day_obs'])
    filename = f"{PREFIXES[channel]}_dayObs_{dayObsStr}_seqNum_{dataId['seq_num']}.png"
    return filename


def _waitForDataProduct(butler, dataProduct, dataId, logger, maxTime=20):
    """Wait for a dataProduct to land inside a repo, timing out in maxTime.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The butler to use.
    dataProduct : `str`
        The dataProduct to wait for, e.g. postISRCCD or calexp etc
    logger : `logging.Logger`
        Logger
    maxTime : `int` or `float`
        The timeout, in seconds, to wait before giving up and returning None.

    Returns
    -------
    dataProduct : dataProduct or None
        Either the dataProduct being waiter for, or None if maxTime elapsed.
    """
    cadence = 0.25
    maxLoops = int(maxTime//cadence)
    for retry in range(maxLoops):
        if datasetExists(butler, dataProduct, dataId):
            return butler.get(dataProduct, dataId)
        else:
            sleep(cadence)
    logger.warning(f'Waited {maxTime}s for {dataProduct} for {dataId} to no avail')
    return None


class Uploader():
    """Class for handling uploads to the Google cloud storage bucket.
    """

    def __init__(self):
        if not HAS_GOOGLE_STORAGE:
            from lsst.rapid.analysis.utils import GOOGLE_CLOUD_MISSING_MSG
            raise RuntimeError(GOOGLE_CLOUD_MISSING_MSG)
        self.client = storage.Client()
        self.bucket = self.client.get_bucket('rubintv_data')
        self.log = logging.getLogger("googleUploader")

    def googleUpload(self, channel, sourceFilename, uploadAsFilename=None):
        """Upload a file to the RubinTV Google cloud storage bucket.

        Parameters
        ----------
        channel : `str`
            The RubinTV channel to upload to.
        sourceFilename : `str`
            The full path and filename of the file to upload.
        uploadAsFilename : `str`, optional
            Optionally rename the file to this upon upload.

        Raises
        ------
        ValueError
            Raised if the specified channel is not in the list of existing
            channels as specified in CHANNELS
        RuntimeError
            Raised if the Google cloud storage is not installed/importable.
        """
        if channel not in CHANNELS:
            raise ValueError(f"Error: {channel} not in {CHANNELS}")

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
    """Class for continuously watching for new data products landing in a repo.

    Wathces a repo for the specified data product to land, and runs a callback
    on the dataId for the data product once it has landed in the repo.
    """
    cadence = 1  # in seconds

    def __init__(self, location, dataProduct, **kwargs):
        self.butler = makeDefaultLatissButler(location)
        self.dataProduct = dataProduct
        self.log = logging.getLogger("watcher")

    def _getLatestImageDataIdAndExpId(self):
        """Get the dataId and expId for the most recent image in the repo.
        """
        dataId = getMostRecentDataId(self.butler)
        expId = getExpIdFromDayObsSeqNum(self.butler, dataId)['exposure']
        return dataId, expId

    def run(self, callback, durationInSeconds=-1):
        """Wait for the dataProduct to land, then run callback(dataId).

        Note that durationInSeconds is a lower bound, but will be a reasonable
        approximation vs the infinite alternative.

        Parameters
        ----------
        callback : `callable`
            The method to call, with the latest dataId as the argument.
        durationInSeconds : `int` or `float`
            How long to run for. This is approximate, as it assumes processing
            is instant. However, most use-cases will want to just use the -1
            sentinel value to run forever anyway.
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
                    lastFound = expId
                    callback(dataId)

            except NotFoundError as e:  # NotFoundError when filters aren't defined
                print(f'Skipped displaying {dataId} due to {e}')
        return


class IsrRunner():
    """Class to run isr for each image that lands in the repo.

    Runs isr via BestEffortIsr, and puts the result in the quickLook rerun.
    """

    def __init__(self, location, **kwargs):
        self.watcher = Watcher(location, 'raw')
        repoDir = LATISS_REPO_LOCATION_MAP[location]
        self.bestEffort = BestEffortIsr(repoDir, **kwargs)
        self.log = logging.getLogger("isrRunner")

    def callback(self, dataId):
        """Method called on each new dataId as it is found in the repo.

        Produce a quickLookExp of the latest image, and butler.put() it to the
        repo so that downstream processes can find and use it.
        """
        quickLookExp = self.bestEffort.getExposure(dataId)  # noqa: F841 - automatically puts
        del quickLookExp
        self.log.info(f'Put quickLookExp for {dataId}, awaiting next image...')

    def run(self):
        """Run continuously, calling the callback method on the latest dataId.
        """
        self.watcher.run(self.callback)


class ImExaminerChannel():
    """Class for running the ImExam channel on RubinTV.
    """

    def __init__(self, location, doRaise=False):
        self.dataProduct = 'quickLookExp'
        self.watcher = Watcher(location, self.dataProduct)
        self.uploader = Uploader()
        self.butler = makeDefaultLatissButler(location)
        self.log = logging.getLogger("imExaminerChannel")
        self.channel = 'summit_imexam'
        self.doRaise = doRaise

    def _imExamine(self, exp, dataId, outputFilename):
        if os.path.exists(outputFilename):  # unnecessary now we're using tmpfile
            self.log.warning(f"Skipping {outputFilename}")
            return
        imexam = ImageExaminer(exp, savePlots=outputFilename, doTweakCentroid=True)
        imexam.plot()

    def callback(self, dataId):
        """Method called on each new dataId as it is found in the repo.

        Plot the quick imExam analysis of the latest image, writing the plot
        to a temp file, and upload it to Google cloud storage via the uploader.
        """
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
            if self.doRaise:
                raise RuntimeError(f"Error processing {dataId}") from e
            self.log.warning(f"Skipped imExam on {dataId} because {repr(e)}")
            return None

    def run(self):
        """Run continuously, calling the callback method on the latest dataId.
        """
        self.watcher.run(self.callback)


class SpecExaminerChannel():
    """Class for running the SpecExam channel on RubinTV.
    """

    def __init__(self, location, doRaise=False):
        self.dataProduct = 'quickLookExp'
        self.watcher = Watcher(location, self.dataProduct)
        self.uploader = Uploader()
        self.butler = makeDefaultLatissButler(location)
        self.log = logging.getLogger("specExaminerChannel")
        self.channel = 'summit_specexam'
        self.doRaise = doRaise

    def _specExamine(self, exp, dataId, outputFilename):
        if os.path.exists(outputFilename):  # unnecessary now we're using tmpfile?
            self.log.warning(f"Skipping {outputFilename}")
            return
        summary = SpectrumExaminer(exp, savePlotAs=outputFilename)
        summary.run()

    def callback(self, dataId):
        """Method called on each new dataId as it is found in the repo.

        Plot the quick spectral reduction of the latest image, writing the plot
        to a temp file, and upload it to Google cloud storage via the uploader.
        """
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
            if self.doRaise:
                raise RuntimeError(f"Error processing {dataId}") from e
            self.log.info(f"Skipped imExam on {dataId} because {repr(e)}")
            return None

    def run(self):
        """Run continuously, calling the callback method on the latest dataId.
        """
        self.watcher.run(self.callback)


class MonitorChannel():
    """Class for running the monitor channel on RubinTV.
    """

    def __init__(self, location, doRaise=False):
        self.dataProduct = 'quickLookExp'
        self.watcher = Watcher(location, self.dataProduct)
        self.uploader = Uploader()
        self.butler = makeDefaultLatissButler(location)
        self.log = logging.getLogger("monitorChannel")
        self.channel = 'auxtel_monitor'
        self.fig = plt.figure(figsize=(12, 12))
        self.doRaise = doRaise

    def _plotImage(self, exp, dataId, outputFilename):
        if os.path.exists(outputFilename):  # unnecessary now we're using tmpfile
            self.log.warning(f"Skipping {outputFilename}")
            return
        plotExp(exp, dataId, self.fig, outputFilename)

    def callback(self, dataId):
        """Method called on each new dataId as it is found in the repo.

        Plot the image for display on the monitor, writing the plot
        to a temp file, and upload it to Google cloud storage via the uploader.
        """
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
            if self.doRaise:
                raise RuntimeError(f"Error processing {dataId}") from e
            self.log.warning(f"Skipped monitor image for {dataId} because {repr(e)}")
            return None

    def run(self):
        """Run continuously, calling the callback method on the latest dataId.
        """
        self.watcher.run(self.callback)


class MountTorqueChannel():
    """Class for running the mount torque channel on RubinTV.
    """

    def __init__(self, location, doRaise=False):
        if not HAS_EFD_CLIENT:
            from lsst.rapid.analysis.utils import EFD_CLIENT_MISSING_MSG
            raise RuntimeError(EFD_CLIENT_MISSING_MSG)
        self.dataProduct = 'raw'
        self.watcher = Watcher(location, self.dataProduct)
        self.uploader = Uploader()
        self.butler = makeDefaultLatissButler(location)
        self.client = EfdClient('summit_efd')
        self.log = logging.getLogger("mountTorqueChannel")
        self.channel = 'auxtel_mount_torques'
        self.fig = plt.figure(figsize=(16, 16))
        self.doRaise = doRaise

    def callback(self, dataId):
        """Method called on each new dataId as it is found in the repo.

        Plot the mount torques, pulling data from the EFD, writing the plot
        to a temp file, and upload it to Google cloud storage via the uploader.
        """
        try:
            tempFilename = tempfile.mktemp(suffix='.png')
            uploadFilename = _dataIdToFilename(self.channel, dataId)
            plotted = plotMountTracking(dataId, self.butler, self.client, self.fig, tempFilename, self.log)

            if plotted:  # skips many image types and short exps
                self.log.info("Uploading mount torque plot to storage bucket")
                self.uploader.googleUpload(self.channel, tempFilename, uploadFilename)
                self.log.info('Upload complete')

        except Exception as e:
            if self.doRaise:
                raise RuntimeError(f"Error processing {dataId}") from e
            self.log.warning(f"Skipped creating mount plots for {dataId} because {repr(e)}")

    def run(self):
        """Run continuously, calling the callback method on the latest dataId.
        """
        self.watcher.run(self.callback)
