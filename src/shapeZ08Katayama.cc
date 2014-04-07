// -*- lsst-c++ -*-

// Feel free to add your own copyright notice here, but note that it must be GPL3-compatible.


#include "boost/make_shared.hpp"

#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/extensions/shapeZ08Katayama.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeZ08Katayama {

ShapeZ08KatayamaControl::ShapeZ08KatayamaControl() :
    // The name below is the prefix for all the output fields.  The number indicates when it
    // should be run relative to other algorithms (3.0 is near the end).
    algorithms::AlgorithmControl("shape.z08", 3.0),
    // The default value for the configuration option defined below.  All configuration options must
    // have defaults defined here.
    nGrowFootprint(0)
{}

// boilerplate
PTR(algorithms::AlgorithmControl) ShapeZ08KatayamaControl::_clone() const {
    return boost::make_shared<ShapeZ08KatayamaControl>(*this);
}

// boilerplate
PTR(algorithms::Algorithm) ShapeZ08KatayamaControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others,
    bool isForced
) const {
    // this calls the Algorithm's constructor, and saves the result in a shared_ptr.
    return boost::make_shared<ShapeZ08KatayamaAlgorithm>(boost::ref(schema), *this, isForced);
}

// In the algorithm's constructor, we'll initialize the Keys as we add them to the Schema we're given.
ShapeZ08KatayamaAlgorithm::ShapeZ08KatayamaAlgorithm(
    afw::table::Schema & schema, Control const & ctrl, bool isForced
) :
    algorithms::Algorithm(ctrl), // let the base class hold the algorithm; we can get it back with getControl
    // feel free to change the definition of any of these fields; they're just here as examples
    // (but I do recommend you include the ellipticity definition in the documentation if you change it)
    _e1Key(schema.addField<double>(
               ctrl.name + ".e1", "on-axis component of ellipticity, using |e|=(1-q^2)/(1+q^2) convention"
           )),
    _e2Key(schema.addField<double>(
               ctrl.name + ".e2", "off-axis component of ellipticity, using |e|=(1-q^2)/(1+q^2) convention"
           )),
    _flagKey(schema.addField<afw::table::Flag>(
                 ctrl.name + ".flags", "general failure flag for " + ctrl.name
             ))
{
    if (isForced) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            ctrl.name + " cannot be run in forced mode"
        );
    }
}

// This is where the implementation goes; I've started by providing some examples you may want to reuse.
//
// You can find reference documentation on many of the classes here:
//
// http://lsst-web.ncsa.illinois.edu/doxygen/x_masterDoxyDoc/
//
// But also feel free to just ask us questions on how to do things with them here:
//
// http://hsca.ipmu.jp:8080/questions/
//
template <typename PixelT>
void ShapeZ08KatayamaAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_flagKey, true); // we set the flag at the beginning and un-set it at the end, so if an
                                // exception is thrown in the middle, it's marked as a failure.

    // This "footprint" object defines the above-threshold pixels that define the object.
    // Here's more information on Footprint:
    PTR(afw::detection::Footprint) footprint = source.getFootprint();
    // We'll use the config parameter in the Control object to tell us how much to grow it.
    footprint = afw::detection::growFootprint(footprint, getControl().nGrowFootprint);
    // If you'd rather use the bounding box of the footprint rather than the footprint itself, here it is:
    afw::geom::Box2I bbox = footprint->getBBox();

    // If the region we want to use goes off the edge of the of the exposure, we need to shrink the bbox
    // so it only includes pixels we actually have.
    int oldArea = bbox.getArea();  // save the desired area of the box
    bbox.clip(exposure.getBBox(afw::image::PARENT));  // shrink the box
    if (bbox.getArea() < oldArea) { // if the actual area is less than the area actually on the exposure...
        source.set(_edgeFlagKey, true); // ...set the flag...
        return; // ...and quit.
    }

    // Here's how to get the PSF model for the image
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    // ..and evaluate it at the position of the object
    PTR(afw::image::Image<double>) psfImage = psf->computeImage(center);

    // Here's an example of how to loop over the pixels in the bbox:
    typedef typename afw::image::Exposure<PixelT>::MaskedImageT::x_iterator Iterator;
    for (int y = bbox.getBeginY(); y < bbox.getEndY(); ++y) {
        Iterator iter = exposure.getMaskedImage().row_begin(y - exposure.getY0());
        iter += (bbox.getBeginX() - exposure.getX0());
        for (int x = bbox.getBeginX(); x < bbox.getEndX(); ++x, ++iter) {
            double dataValue = iter.image();    // pixel value at (x,y)
            double varianceValue = iter.variance(); // variance value at (x,y)
            double relX = x - center.getX(); // position of the pixel relative to the center of the object
            double relY = y - center.getY();
            // do measurements here!  (or copy to some other data structure)
        }
    }

    source.set(_flagKey, false); // if we get this far, we've succeeded
    source.set(_e1Key, 3.14); // replace these with the actual measurements, of course
    source.set(_e2Key, -3.14);
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(ShapeZ08KatayamaAlgorithm);

}}}} // namespace lsst::meas::extensions::shapeZ08Katayama
