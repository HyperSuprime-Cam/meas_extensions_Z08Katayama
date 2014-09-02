// -*- lsst-c++ -*-
/*  This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "boost/make_shared.hpp"

#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/extensions/shapeZ08Katayama.h"
#include "fft2.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeZ08Katayama {

namespace my {
    template <class T>
    inline T norm(std::complex<T> const& c){
        // Old GCC's std::norm() is not good: it computes hypot(a,b)^2.
        // We provide a faster one here.
        return c.real()*c.real() + c.imag()*c.imag();
    }
}

ShapeZ08KatayamaControl::ShapeZ08KatayamaControl() :
    // The name below is the prefix for all the output fields.  The number indicates when it
    // should be run relative to other algorithms (3.0 is near the end).
    algorithms::AlgorithmControl("shape.z08", 3.0),
    // The default value for the configuration option defined below.  All configuration options must
    // have defaults defined here.
    beta       (10.0 ),
    k2_limit   ( 0.01),
    bbox_min   (32   ),
    bbox_scale ( 2.0 )
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
    _num1Key(schema.addField<double>(
               ctrl.name + ".num1", "on-axis component of ellipticity, using |e|=(1-q^2)/(1+q^2) convention, multiplied by denom"
           )),
    _num2Key(schema.addField<double>(
               ctrl.name + ".num2", "off-axis component of ellipticity, using |e|=(1-q^2)/(1+q^2) convention, multiplied by denom"
           )),
    _denomKey(schema.addField<double>(
               ctrl.name + ".denom", "a denominator such that ellip = (num1, num2) / denom"
           )),
    _flagKey(schema.addField<afw::table::Flag>(
                 ctrl.name + ".flags", "general failure flag for " + ctrl.name
             )),
    _edgeFlagKey(schema.addField<afw::table::Flag>(
                 ctrl.name + ".isedge", "source was at edge"
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

    bool isAnyError = false;

    //
    // The object's image
    //

    // This "footprint" object defines the above-threshold pixels that define the object.
    // Here's more information on Footprint:
    PTR(afw::detection::Footprint) footprint = source.getFootprint();
    // If you'd rather use the bounding box of the footprint rather than the footprint itself, here it is:
    afw::geom::Box2I bbox = footprint->getBBox();

    // enlarge the bounding box by bbox_scale
    double bbox_scale = getControl().bbox_scale;

    double left_relto_center   = bbox_scale * (bbox.getMinX() - center.getX());
    double right_relto_center  = bbox_scale * (bbox.getMaxX() - center.getX());
    double bottom_relto_center = bbox_scale * (bbox.getMinY() - center.getY());
    double top_relto_center    = bbox_scale * (bbox.getMaxY() - center.getY());

    // make the bounding box square
    double height_minus_width =
        (top_relto_center - bottom_relto_center) - (right_relto_center - left_relto_center);
    if(height_minus_width > 0){
        left_relto_center  -= 0.5 * height_minus_width;
        right_relto_center += 0.5 * height_minus_width;
    }
    else if(height_minus_width < 0){
        bottom_relto_center += 0.5 * height_minus_width;
        top_relto_center    -= 0.5 * height_minus_width;
    }

    // ensure the bounding box >= bbox_min
    double bbox_min_minus_1 = getControl().bbox_min - 1;
    if(right_relto_center - left_relto_center < bbox_min_minus_1){
        double addenda = 0.5 * (bbox_min_minus_1 - (right_relto_center - left_relto_center));
        left_relto_center   -= addenda;
        right_relto_center  += addenda;
        bottom_relto_center -= addenda;
        top_relto_center    += addenda;
    }

    // quantize the bounding box
    int left   = (int)std::floor(left_relto_center   + center.getX() + 0.5);
    int right  = (int)std::floor(right_relto_center  + center.getX() + 0.5);
    int bottom = (int)std::floor(bottom_relto_center + center.getY() + 0.5);
    int top    = (int)std::floor(top_relto_center    + center.getY() + 0.5);

    // make the bounding box square (again)
    switch((top - bottom) - (right - left)){
    case  2:
        left_relto_center   -= 1;
        right_relto_center  += 1;
        break;
    case  1:
        right_relto_center  += 1;
        break;
    case -1:
        top_relto_center    += 1;
        break;
    case -2:
        bottom_relto_center -= 1;
        top_relto_center    += 1;
        break;
    }

    // shift the bounding box so that (left, right, bottom, top) are relative to
    // the lower left corner of the exposure
    left   -= exposure.getX0();
    right  -= exposure.getX0();
    bottom -= exposure.getY0();
    top    -= exposure.getY0();

    std::size_t width  = right - left + 1;
    std::size_t height = width;
    std::size_t nPixels = height * width;

    // get exposure's bounding box
    bbox = exposure.getBBox(afw::image::LOCAL);

    // the region to use
    int copied_left   = (std::max)(left  , bbox.getMinX());
    int copied_right  = (std::min)(right , bbox.getMaxX());
    int copied_bottom = (std::max)(bottom, bbox.getMinY());
    int copied_top    = (std::min)(top   , bbox.getMaxY());
    std::size_t copied_width  = copied_right - copied_left + 1;

    if (  (left   ^ copied_left  ) // if left   != copied_left
        | (right  ^ copied_right ) // or right  != copied_right
        | (bottom ^ copied_bottom) // or bottom != copied_bottom
        | (top    ^ copied_top   ) // or top    != copied_top
    ) {
        source.set(_edgeFlagKey, true); // ...set the flag...
        isAnyError = true;
    }

    // copy the image for fourier transform
    afw::image::Image<PixelT>& image = *exposure.getMaskedImage().getImage();
    typedef typename afw::image::Image<PixelT>::x_iterator  x_iteratorT;

    shapeZ08::Fft2 fft2(height, width, shapeZ08::Fft2::Forward);
    std::fill(fft2.spaceDomain(), fft2.spaceDomain() + nPixels, 0.0);

    for(int y = copied_bottom; y <= copied_top; ++y){
        x_iteratorT src = image.row_begin(y) + copied_left;
        std::copy(src, src + copied_width,
            fft2.spaceDomain() + ((y - bottom) * width + (copied_left - left))
        );
    }

    // transform!
    fft2.forward();

    // save the fourier-transformed image
    std::vector<std::complex<double> > MF(fft2.freqDomain(), fft2.freqDomain() + nPixels);

    //
    // The original PSF
    //

    // Here's how to get the PSF model for the image
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    // ..and evaluate it at the position of the object
    PTR(afw::image::Image<double>) psfImage = psf->computeImage(center);

    // get the bounding box of the PSF image in the exposure
    bbox = psfImage->getBBox(afw::image::LOCAL);

    // shift the bounding box so that (left, right, bottom, top) are relative to
    // the lower left corner of the PSF image
    left   += exposure.getX0() - psfImage->getX0();
    right  += exposure.getX0() - psfImage->getX0();
    bottom += exposure.getY0() - psfImage->getY0();
    top    += exposure.getY0() - psfImage->getY0();

    // the region to use
    copied_left   = (std::max)(left  , bbox.getMinX());
    copied_right  = (std::min)(right , bbox.getMaxX());
    copied_bottom = (std::max)(bottom, bbox.getMinY());
    copied_top    = (std::min)(top   , bbox.getMaxY());
    copied_width  = copied_right - copied_left + 1;

    // copy the PSF image for fourier transform
    std::fill(fft2.spaceDomain(), fft2.spaceDomain() + nPixels, 0.0);

    for(int y = copied_bottom; y <= copied_top; ++y){
        afw::image::Image<double>::x_iterator
            src = psfImage->row_begin(y) + copied_left;
        std::copy(src, src + copied_width,
            fft2.spaceDomain() + ((y - bottom) * width + (copied_left - left))
        );
    }

    // transform!
    fft2.forward();

    // save the fourier-transformed psfImage
    std::vector<std::complex<double> > PSF_F_original(fft2.freqDomain(), fft2.freqDomain() + nPixels);

    //
    // The target PSF
    //

    std::vector<std::complex<double> > PSF_F_target;

    double beta = getControl().beta;
    if(beta != 0){
        // create a target PSF (gaussian)
        double mI2bb = -0.5 / (beta*beta);
        double cy   = 0.5 * (height - 1);
        double cx   = 0.5 * (width  - 1);

        std::complex<double>* dest = fft2.spaceDomain();
        for(int y = 0; y < (int)height; ++y){
            double yy = (y - cy)*(y - cy);
            for(int x = 0; x < (int)width; ++x){
                double xx = (x - cx)*(x - cx);
                *(dest++) = std::exp(mI2bb * (xx + yy));
            }
        }

        // transform!
        fft2.forward();

        // save the fourier-transformed target PSF
        PSF_F_target.assign(fft2.freqDomain(), fft2.freqDomain() + nPixels);
    }
    else{
        PSF_F_target.assign(nPixels, 1.0);
    }

    //
    // The moments
    //

    double Mx2 = 0.0;
    double My2 = 0.0;
    double Mxy = 0.0;
    double M4  = 0.0;
    double const c = 4 * M_PI*M_PI / (beta*beta);

    double const k2_limit = getControl().k2_limit;
    double const Ih = 1.0 / (double)(int)height;
    double const Iw = 1.0 / (double)(int)width ;

    for(int i = 0; i < (int)height; ++i){
        int iy = (2*i < (int)height) ? i : i - height;
        double ky = (double)iy * Ih;
        double ky2 = ky * ky;

        for(int j = 0; j < (int)width; ++j){
            int jx = (2*j < (int)width) ? j : j - width;
            double kx = (double)jx * Iw;
            double kx2 = kx * kx;

            if(kx2+ky2 < k2_limit){
                double w
                    = my::norm(MF            [i*width + j])
                    * my::norm(PSF_F_target  [i*width + j])
                    / my::norm(PSF_F_original[i*width + j])
                    ;
                Mx2 += w * kx2                  ;
                My2 += w * ky2                  ;
                Mxy += w * kx*ky                ;
                M4  += w * ((kx2+ky2)*(kx2+ky2));
            }
        }
    }

    source.set(_flagKey , isAnyError       );
    source.set(_num1Key , 0.5 * (Mx2 - My2));
    source.set(_num2Key , -Mxy             );
    source.set(_denomKey, Mx2 + My2 - c*M4 );
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(ShapeZ08KatayamaAlgorithm);

}}}} // namespace lsst::meas::extensions::shapeZ08Katayama
