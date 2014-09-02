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
#include "fft2.h"

#include <fftw3.h>
#include <stdexcept>
#include <cstddef>


namespace lsst { namespace meas { namespace extensions { namespace shapeZ08 {

Fft2::Fft2(int height, int width, Direction dir)
    : space_        ()
    , freq_         ()
    , plan_forward_ ()
    , plan_backward_()
    , height_       (height)
    , width_        (width)
{
    try {
        std::size_t size = (std::size_t)height * (std::size_t)width;
        space_ = (element*)fftw_malloc(size * sizeof(element));
        freq_  = (element*)fftw_malloc(size * sizeof(element));

        if(!space_ || !freq_){
            throw std::bad_alloc();
        }

        if(dir != None) this->prepare(dir);
    }
    catch(...){
        this->invalidate();
        throw;
    }
}


void Fft2::prepare(Direction dir)
{
    if((dir & Forward)
    && !plan_forward_
    ){
        plan_forward_ = fftw_plan_dft_2d(
            height_, width_,
            (fftw_complex*)space_, (fftw_complex*)freq_,
            FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT
        );
        if(!plan_forward_){
            throw std::runtime_error("fftw3: creating plan failed.");
        }
    }

    if((dir & Backward)
    && !plan_backward_
    ){
        plan_backward_ = fftw_plan_dft_2d(
            height_, width_,
            (fftw_complex*)freq_, (fftw_complex*)space_,
            FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT
        );
        if(!plan_backward_){
            throw std::runtime_error("fftw3: creating plan failed.");
        }
    }
}


void Fft2::forward()
{
    fftw_execute(static_cast<fftw_plan>(plan_forward_));
}


void Fft2::backward()
{
    fftw_execute(static_cast<fftw_plan>(plan_backward_));
}


void Fft2::invalidate()
{
    if(plan_forward_){
        fftw_destroy_plan(static_cast<fftw_plan>(plan_forward_));
        plan_forward_ = NULL;
    }
    if(plan_backward_){
        fftw_destroy_plan(static_cast<fftw_plan>(plan_backward_));
        plan_backward_ = NULL;
    }

    if(space_){
        fftw_free(space_);
        space_ = NULL;
    }
    if(freq_ ){
        fftw_free(freq_ );
        freq_ = NULL;
    }
}

}}}} // lsst::meas::extensions::shapeZ08
