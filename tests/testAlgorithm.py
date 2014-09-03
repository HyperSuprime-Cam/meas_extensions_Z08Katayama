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
import os, math
import numpy

import lsst.utils.tests
import lsst.afw.table
import lsst.afw.geom
import lsst.afw.image
import lsst.meas.algorithms
import lsst.meas.extensions.shapeZ08Katayama

numpy.random.seed(500)

class AlgorithmTestCase(lsst.utils.tests.TestCase):

    def fillTestExposure(self, exposure, radius, ellipticity, noiseSigma=None):
        """A utility function to add a PSF-convolved elliptical exponential disk to an Exposure
        """
        # compute a raw galaxy image
        bbox = exposure.getBBox(lsst.afw.image.PARENT)
        radius /= 1.67835 # half light radius => effective radius

        x = numpy.arange(bbox.getBeginX(), bbox.getEndX()) / radius
        y = numpy.arange(bbox.getBeginY(), bbox.getEndY()) / radius
        y.shape = (len(y), 1)

        imageData = numpy.exp(
            -numpy.sqrt(x**2 + y**2 - (ellipticity[0]*(x**2 - y**2) + 2*ellipticity[1]*x*y))
        ) / (radius*radius)

        # compute a psf image
        imagePsf = exposure.getPsf().computeImage().getArray()
        imagePsfNew = numpy.zeros(shape = imageData.shape)
        h = min(imagePsf.shape[0], imagePsfNew.shape[0])
        w = min(imagePsf.shape[1], imagePsfNew.shape[1])
        imagePsfNew[0:h, 0:w] = imagePsf[0:h, 0:w]

        # psf's center
        y = numpy.arange(imagePsfNew.shape[0])
        y.shape = (len(y), 1)
        x = numpy.arange(imagePsfNew.shape[1])

        centerX = (imagePsfNew * x).sum() / imagePsfNew.sum()
        centerY = (imagePsfNew * y).sum() / imagePsfNew.sum()

        # compute exp(ikx)
        shiftY = numpy.fft.fftfreq(imagePsfNew.shape[0]) * (2j*math.pi*centerY)
        shiftY.shape = (len(shiftY), 1)
        shiftX = numpy.fft.fftfreq(imagePsfNew.shape[1]) * (2j*math.pi*centerX)

        shift = numpy.exp(shiftX + shiftY)

        # convolve psf to the galaxy
        imageData = numpy.fft.ifft2(
            numpy.fft.fft2(imageData) * numpy.fft.fft2(imagePsfNew) * shift
        ).real

        exposure.getMaskedImage().getImage().getArray()[:,:] = imageData
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
        exposure.setPsf(lsst.afw.detection.GaussianPsf(17, 17, 1.0))  # width, height, sigma
        radius = 7
        ellipticity = (0.3, 0.51961524227066314) # a:b = 2:1, angle = 30deg
        self.fillTestExposure(exposure, radius, ellipticity, noiseSigma=0.0)
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

        print record.get("shape.z08.flags")
        print record.get("shape.z08.isedge")
        e1 = record.get("shape.z08.num1") / record.get("shape.z08.denom")
        e2 = record.get("shape.z08.num2") / record.get("shape.z08.denom")

        # Test that the results are what we expect (should update expected values when algorithm is ready)
        #self.assertFalse(record.get("shape.z08.flags"))  # check that failure flag is not set
        self.assertClose(e1, ellipticity[0])
        self.assertClose(e2, ellipticity[1])


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
