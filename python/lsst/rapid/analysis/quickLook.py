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

from lsst.ip.isr import IsrTask
from lsst.ip.isr.isrTask import IsrTaskConnections
import lsst.pipe.base.connectionTypes as cT

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask

__all__ = ['QuickLookTask', 'QuickLookTaskConfig']


def _getArgs(connection):
    """Get the all the required args from a connection in order to reconstruct.
    """
    newArgs = {}
    for attr in ("name", "storageClass", "multiple", "doc", "dimensions", "isCalibration", "deferLoad",
                 "minimum", "lookupFunction"):
        if hasattr(connection, attr):
            newArgs[attr] = getattr(connection, attr)
    return newArgs


class QuickLookTaskConnections(IsrTaskConnections):
    def __init__(self, *, config=None):
        """Copy isrTask's connections, changing prereq min values to zero.

        Copy all the connections directly for IsrTask, keeping ccdExposure as
        required as non-zero, but changing all the other PrerequisiteInputs'
        minimum values to zero.
        """
        super().__init__(config=IsrTask.ConfigClass())  # need a dummy config, isn't used other than for ctor
        for name, connection in self.allConnections.items():
            if hasattr(connection, 'minimum'):
                args = _getArgs(connection)
                if name != "ccdExposure":  # need one input image always
                    args['minimum'] = 0
                newConnection = type(connection)(**args)
                self.allConnections[name] = newConnection
                setattr(self, name, newConnection)

        exposure = cT.Output(  # called just "exposure" to mimic isrTask's return struct
            name="quickLookExp",
            doc="The quickLook output exposure.",
            storageClass="ExposureF",
            dimensions=("instrument", "exposure", "detector"),
        )
        self.outputExposure = exposure


class QuickLookTaskConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=QuickLookTaskConnections):
    """Configuration parameters for QuickLookTask."""

    doRepairCosmics = pexConfig.Field(
        dtype=bool,
        doc="Interpolate over cosmic rays?",
        default=True,
    )
    doWrite = pexConfig.Field(  # XXX do we want this option? If yes, make it work.
        dtype=bool,
        doc="Write out the results?",
        default=True,
    )


class QuickLookTask(pipeBase.PipelineTask):

    """Task to automatically perform as much isr as possible.

    Automatically performs as much isr as is possible, depending on the
    calibration products available. Applies these, and then optionally
    interpolates over cosmic rays.
    """

    ConfigClass = QuickLookTaskConfig
    RunnerClass = pipeBase.ButlerInitializedTaskRunner
    _DefaultName = "quickLook"

    def __init__(self, *, butler=None, **kwargs):
        super().__init__(**kwargs)

    def run(self, ccdExposure, *,
            camera=None,
            bias=None,
            dark=None,
            flat=None,
            defects=None,
            linearizer=None,
            crosstalk=None,
            bfKernel=None,
            bfGains=None,
            ptc=None,
            crosstalkSources=None,
            isrBaseConfig=None
            ):
        """Run isr and cosmic ray repair using, doing as much isr as possible.

        Retrieves as many calibration products as are available, and runs isr
        with those settings enabled, but always returns an assembled image at
        a minimum. Then performs cosmic ray repair if configured to.

        Parameters
        ----------
        ccdExposure : `lsst.afw.image.Exposure`
            The raw exposure that is to be run through ISR.  The
            exposure is modified by this method.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            The camera geometry for this exposure. Required if
            one or more of ``ccdExposure``, ``bias``, ``dark``, or
            ``flat`` does not have an associated detector.
        bias : `lsst.afw.image.Exposure`, optional
            Bias calibration frame.
        linearizer : `lsst.ip.isr.linearize.LinearizeBase`, optional
            Functor for linearization.
        crosstalk : `lsst.ip.isr.crosstalk.CrosstalkCalib`, optional
            Calibration for crosstalk.
        crosstalkSources : `list`, optional
            List of possible crosstalk sources.
        dark : `lsst.afw.image.Exposure`, optional
            Dark calibration frame.
        flat : `lsst.afw.image.Exposure`, optional
            Flat calibration frame.
        ptc : `lsst.ip.isr.PhotonTransferCurveDataset`, optional
            Photon transfer curve dataset, with, e.g., gains
            and read noise.
        bfKernel : `numpy.ndarray`, optional
            Brighter-fatter kernel.
        bfGains : `dict` of `float`, optional
            Gains used to override the detector's nominal gains for the
            brighter-fatter correction. A dict keyed by amplifier name for
            the detector in question.
        defects : `lsst.ip.isr.Defects`, optional
            List of defects.
        fringes : `lsst.pipe.base.Struct`, optional
            Struct containing the fringe correction data, with
            elements:
            - ``fringes``: fringe calibration frame (`afw.image.Exposure`)
            - ``seed``: random seed derived from the ccdExposureId for random
                number generator (`uint32`)
        opticsTransmission: `lsst.afw.image.TransmissionCurve`, optional
            A ``TransmissionCurve`` that represents the throughput of the,
            optics, to be evaluated in focal-plane coordinates.
        filterTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
            filter itself, to be evaluated in focal-plane coordinates.
        sensorTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
            sensor itself, to be evaluated in post-assembly trimmed detector
            coordinates.
        atmosphereTransmission : `lsst.afw.image.TransmissionCurve`
            A ``TransmissionCurve`` that represents the throughput of the
            atmosphere, assumed to be spatially constant.
        detectorNum : `int`, optional
            The integer number for the detector to process.
        strayLightData : `object`, optional
            Opaque object containing calibration information for stray-light
            correction.  If `None`, no correction will be performed.
        illumMaskedImage : `lsst.afw.image.MaskedImage`, optional
            Illumination correction image.
        isrBaseConfig : `lsst.ip.isr.IsrTaskConfig`, optional
            An isrTask config to act as the base configuration. Options which
            involve apply a calibration product are ignored, but this allows
            for the configuration of e.g. the number of overscan columns.
        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with component:
            - ``exposure`` : `afw.image.Exposure`
                The ISRed and cosmic-ray-repaired exposure.

        """
        isrConfig = isrBaseConfig if isrBaseConfig else IsrTask.ConfigClass()
        isrConfig.doBias = False
        isrConfig.doDark = False
        isrConfig.doFlat = False
        isrConfig.doFringe = False
        isrConfig.doDefect = False
        isrConfig.doLinearize = False
        isrConfig.doCrosstalk = False
        isrConfig.doBrighterFatter = False
        isrConfig.usePtcGains = False

        if bias:
            isrConfig.doBias = True
            self.log.info("Running with bias correction")

        if dark:
            isrConfig.doDark = True
            self.log.info("Running with dark correction")

        if flat:
            isrConfig.doFlat = True
            self.log.info("Running with flat correction")

        # TODO: deal with fringes here
        if defects:
            isrConfig.doDefect = True
            self.log.info("Running with defect correction")

        if linearizer:
            isrConfig.doLinearize = True
            self.log.info("Running with linearity correction")

        if crosstalk:
            isrConfig.doCrosstalk = True
            self.log.info("Running with crosstalk correction")

        if bfKernel:
            isrConfig.doBrighterFatter = True
            self.log.info("Running with brighter-fatter correction")

        if ptc:
            isrConfig.usePtcGains = True
            self.log.info("Running with ptc correction")

        isrConfig.doWrite = False
        isrTask = IsrTask(config=isrConfig)
        result = isrTask.run(ccdExposure,
                             camera=camera,
                             bias=bias,
                             dark=dark,
                             flat=flat,
                             #  fringes=pipeBase.Struct(fringes=None),
                             defects=defects,
                             linearizer=linearizer,
                             crosstalk=crosstalk,
                             bfKernel=bfKernel,
                             bfGains=bfGains,
                             ptc=ptc,
                             crosstalkSources=crosstalkSources,
                             isGen3=True,)

        postIsr = result.exposure

        if self.config.doRepairCosmics:
            try:  # can fail due to too many CRs detected, and we always want an exposure back
                self.log.info("Repairing cosmics...")
                if postIsr.getPsf() is None:
                    installPsfTask = InstallGaussianPsfTask()
                    installPsfTask.run(postIsr)

                repairConfig = CharacterizeImageTask.ConfigClass()
                repairConfig.doMeasurePsf = False
                repairConfig.doApCorr = False
                repairConfig.doDeblend = False
                repairConfig.repair.cosmicray.nCrPixelMax = 200000
                repairTask = CharacterizeImageTask(config=repairConfig)

                repairTask.repair.run(postIsr)
            except Exception as e:
                self.log.warning(f"During CR repair caught: {e}")

        # outputExposure is persisted, exposure is returned for convenience to
        # mimic Gen2 isrTask's API.
        return pipeBase.Struct(exposure=postIsr,
                               outputExposure=postIsr)
