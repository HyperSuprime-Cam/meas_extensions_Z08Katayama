#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import unittest
import os
import numpy

import lsst.utils.tests
import lsst.afw.table
import lsst.afw.geom
import lsst.afw.image
import lsst.meas.algorithms
import lsst.meas.extensions.shapeZ08Katayama

numpy.random.seed(500)

class AlgorithmTestCase(lsst.utils.tests.TestCase):

    def fillTestExposure(self, exposure, ellipse, noiseSigma=None):
        """A utility function to add a PSF-convolved elliptical Gaussian to an Exposure
        """
        psfSigma = lsst.afw.detection.GaussianPsf.cast(exposure.getPsf()).getSigma()
        psfEllipse = lsst.afw.geom.ellipses.Quadrupole(psfSigma**2, psfSigma**2, 0.0)
        convolvedEllipse = ellipse.convolve(psfEllipse)
        t = convolvedEllipse.getGridTransform()
        bbox = exposure.getBBox(lsst.afw.image.PARENT)
        x, y = numpy.meshgrid(numpy.arange(bbox.getBeginX(), bbox.getEndX()),
                              numpy.arange(bbox.getBeginY(), bbox.getEndY()))
        xt = t[t.XX]*x + t[t.XY]*y
        yt = t[t.YX]*x + t[t.YY]*y
        exposure.getMaskedImage().getImage().getArray()[:,:] = numpy.exp(-0.5*(x**2 + y**2))
        if noiseSigma is not None:
            exposure.getMaskedImage().getImage().getArray()[:,:] \
                += numpy.random.randn(bbox.getHeight(), bbox.getWidth()) * noiseSigma
            exposure.getMaskedImage().getVariance().getArray()[:,:] = noiseSigma**2

    def test(self):
        """Test a measurement algorithm by running it on an elliptical Gaussian convolved with
        a circular Gaussian PSF.
        """

        # Setup test data
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-30, -30), lsst.afw.geom.Point2I(30, 30))
        exposure = lsst.afw.image.ExposureF(bbox)
        exposure.setPsf(lsst.afw.detection.GaussianPsf(17, 17, 1.5))  # width, height, sigma
        ellipse = lsst.afw.geom.ellipses.Axes(6.0, 5.0, 0.5)   # a, b, theta of galaxy, pre-convolution
        self.fillTestExposure(exposure, ellipse, noiseSigma=0.01)
        # we'll use a square region that covers the full postage stamp for the Footprint
        footprint = lsst.afw.detection.Footprint(bbox)
        footprint.getPeaks().push_back(lsst.afw.detection.Peak(0, 0))

        # Prepare measurement machinery
        config = lsst.meas.algorithms.SourceMeasurementConfig()
        config.algorithms.names = ["shape.z08"]
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        ms = config.makeMeasureSources(schema)
        table = lsst.afw.table.SourceTable.make(schema)

        # Measure the object
        record = table.makeRecord()
        record.setFootprint(footprint)
        ms.applyWithPeak(record, exposure)

        # Test that the results are what we expect (should update expected values when algorithm is ready)
        self.assertFalse(record.get("shape.z08.flags"))  # check that failure flag is not set
        self.assertClose(record.get("shape.z08.e1"), 3.14)
        self.assertClose(record.get("shape.z08.e2"), -3.14)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(AlgorithmTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
