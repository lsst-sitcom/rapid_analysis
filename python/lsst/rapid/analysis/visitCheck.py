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


__all__ = ['CheckVisitTask',
           'CheckVisitTaskConfig',
           ]


import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
# import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import logging
from lsst.cp.pipe.utils import countMaskedPixels

from lsst.ip.isr import IsrTask, AssembleCcdTask
from astro_metadata_translator import ObservationInfo


class BasicTestTask(pipeBase.Task):
    logger = logging.getLogger(' ')  # empty string results in name being "root"

    def computeImageStats(self, exposure):
        opts = afwMath.MEAN | afwMath.STDEV | afwMath.MEANCLIP | afwMath.STDEVCLIP
        stats = afwMath.makeStatistics(exposure.maskedImage, opts)
        mean = stats.getValue(afwMath.MEAN)
        std = stats.getValue(afwMath.STDEV)
        meanClip = stats.getValue(afwMath.MEANCLIP)
        stdClip = stats.getValue(afwMath.STDEVCLIP)

        # add min, max, percentiles xxx

        nBad = countMaskedPixels(exposure, 'BAD')
        msgs = []
        msgs.append(f"Mean         = {mean:.2f}")
        msgs.append(f"Clipped mean = {meanClip:.2f}")
        msgs.append(f"Std dev      = {std:.2f}")
        msgs.append(f"Clipped std  = {stdClip:.2f}")
        msgs.append(f"nBad pix     = {nBad}")

        for msg in msgs:
            self.logger.info(msg)


class BiasTestConfig(pexConfig.Config):
    someOption = pexConfig.Field(
        doc="An option",
        dtype=bool,
        default=True,
    )


class BiasTestTask(BasicTestTask):
    ConfigClass = BiasTestConfig

    def runDataRef(self, dataRef, exposure):
        bias = dataRef.get('bias')
        detector = exposure.getDetector()

        # sanity checking - these should NEVER fail
        assert bias.getDetector().getName() == detector.getName()
        assert bias.getDetector().getId() == detector.getId()

        if bias.image.array.shape != exposure.image.array.shape:
            raise RuntimeError("Bias and exposure are different shapes")

        # import ipdb as pdb; pdb.set_trace()

    def checkSomething(self, exposure, bias):
        pass


class CheckVisitTaskConfig(pexConfig.Config):
    """Config class for visit checking"""

    isr = pexConfig.ConfigurableField(
        target=IsrTask,
        doc="Task to perform instrumental signature removal",
    )
    ccdKey = pexConfig.Field(
        dtype=str,
        doc="The key by which to pull a detector from a dataId, e.g. 'ccd' or 'detector'",
        default='detector',
    )
    imageTypeKey = pexConfig.Field(
        dtype=str,
        doc="The key for the butler to use by which to check whether images are darks or flats",
        default='imageType',
    )
    biasTask = pexConfig.ConfigurableField(
        target=BiasTestTask,
        doc="CCD assembly task",
    )
    assembleCcd = pexConfig.ConfigurableField(
        target=AssembleCcdTask,
        doc="CCD assembly task",
    )


class CheckVisitTask(pipeBase.CmdLineTask):
    """Task for checking raw visits from the camera."""

    ConfigClass = CheckVisitTaskConfig
    _DefaultName = "checkVisit"

    def __init__(self, *args, **kwargs):
        pipeBase.CmdLineTask.__init__(self, *args, **kwargs)
        self.makeSubtask("isr")
        self.makeSubtask("assembleCcd")
        self.makeSubtask("biasTask")

    @pipeBase.timeMethod
    def runDataRef(self, dataRef):
        """XXX docstring

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            dataRef for the detector for the visits to be fit.


        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with Components:
            - ``exitStatus`` : `int`
              The exit code.
        """

        detNum = dataRef.dataId[self.config.ccdKey]
        msg = "Calculating defects using %s visits for detector %s"
        self.log.info(msg, dataRef.dataId['visit'], detNum)

        visit = dataRef.dataId['visit']
        obsInfoImageType = ObservationInfo(dataRef.get("raw_md")).observation_type
        registryImageType = dataRef.getButler().queryMetadata('raw', 'imageType', visit=visit)[0]

        self.log.debug('imageType from ObservationInfo = %s' % obsInfoImageType)
        self.log.debug('imageType from registry        = %s' % registryImageType)

        # imageType = self.getImageType(dataRef)
        unassembledExp = dataRef.get('raw')
        exp = self.assembleCcd.assembleCcd(unassembledExp)
        self.biasTask.computeImageStats(exp)
        self.biasTask.runDataRef(dataRef, exp)
        # import ipdb as pdb; pdb.set_trace()

        return self

        return pipeBase.Struct(exitStatus=0)

    def checkHeaders(self, dataRef):
        return
        # raw_md = dataRef.get("raw_md")
        # obsInfo = ObservationInfo(raw_md)

        # expTime = dataRef.getButler().queryMetadata('raw', 'expTime', ...
        # ... visit=dataRef.dataId['visit'])[0]
        # xxx continue writing here

        # checks: exptime is >0 for nonbias
        # checks: exptime is ==0 for bias
        # some sanity check for filter
        # some sanity check for temps?
        # sequencer check

    @staticmethod
    def getImageType(dataRef, useObsInfoMethod=False):
        if useObsInfoMethod:
            return ObservationInfo(dataRef.get("raw_md")).observation_type
        return dataRef.getButler().queryMetadata('raw', 'imageType', visit=dataRef.dataId['visit'])[0]
