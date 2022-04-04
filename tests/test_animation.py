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

"""Test cases for animations."""

import unittest
import tempfile
import os

import lsst.utils.tests

from lsst.rapid.analysis.animation import Animator

import lsst.daf.butler as dafButler
NOBUTLER = True
try:
    butler = dafButler.Butler('LATISS')
    assert isinstance(butler, dafButler.Butler)
    NOBUTLER = False
except (FileNotFoundError):
    print("No LATISS butler found, skipping butler-driven tests.")


class AnimationTestCase(lsst.utils.tests.TestCase):
    """A test case for testing the animator."""

    @unittest.skipIf(NOBUTLER, 'Skipping butler-driven test')
    def setUp(self):
        self.butler = dafButler.Butler('LATISS', instrument='LATISS', collections=['LATISS/raw/all'])

        # TODO: DM-34322 Change these to work with test data on the TTS once
        # data has been ingested there.
        self.dataIds = [{'day_obs': 20200315, 'seq_num': 30, 'detector': 0},
                        {'day_obs': 20200315, 'seq_num': 31, 'detector': 0}]
        self.outputDir = tempfile.mkdtemp()
        self.outputFilename = os.path.join(self.outputDir, 'testAnimation.mp4')

    @unittest.skipIf(NOBUTLER, 'Skipping butler-driven test')
    def test_animation(self):
        animator = Animator(self.butler, self.dataIds, self.outputDir, self.outputFilename,
                            dataProductToPlot='raw',
                            remakePngs=True,
                            debug=False,
                            clobberVideoAndGif=True,
                            plotObjectCentroids=True,
                            useQfmForCentroids=True)
        animator.run()

        self.assertTrue(os.path.isfile(self.outputFilename))
        self.assertTrue(os.path.getsize(self.outputFilename) > 10000)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
