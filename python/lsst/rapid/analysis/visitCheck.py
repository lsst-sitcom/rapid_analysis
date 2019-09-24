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
           'CheckVisitTaskConfig', ]


import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
# import lsst.afw.image as afwImage
# import lsst.afw.math as afwMath
# import lsst.afw.display as afwDisplay
# from lsst.cp.pipe.utils import countMaskedPixels

from lsst.ip.isr import IsrTask


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


class CheckVisitTask(pipeBase.CmdLineTask):
    """Task for checking raw visits from the camera."""

    ConfigClass = CheckVisitTaskConfig
    _DefaultName = "checkVisit"

    def __init__(self, *args, **kwargs):
        pipeBase.CmdLineTask.__init__(self, *args, **kwargs)
        self.makeSubtask("isr")
        self.config.validate()
        self.config.freeze()

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

            - ``defects`` : `lsst.meas.algorithms.Defect`
              The defects found by the task.
            - ``exitStatus`` : `int`
              The exit code.
        """

        detNum = dataRef.dataId[self.config.ccdKey]
        msg = "Calculating defects using %s visits for detector %s"
        self.log.info(msg, dataRef.dataId['visit'], detNum)

        import ipdb as pdb; pdb.set_trace()

        return pipeBase.Struct(exitStatus=0)
