#ifndef  ge93cecf0_2c00_40d7_ad6e_518123ea8be8
#define  ge93cecf0_2c00_40d7_ad6e_518123ea8be8
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

#include <complex>

namespace lsst { namespace meas { namespace extensions { namespace shapeZ08 {

/** 2D fourier transformer
    A typical usage is as follows:

    int height = 10, width = 10;

    // FFT in both directions
    Fft2 fft2(height, width, Both);

    // fill in the space domain image
    std::complex<double>* image = fft2.spaceDomain();
    for(int y = 0; y < height; ++y){
        for(int x = 0; x < width; ++x){
            image[y*width + x] = f(x,y);
        }
    }

    // forward transform
    fft2.forward();

    // multiply a filter in the frequency domain
    image = fft2.freqDomain();
    for(int ky = 0; ky < height; ++ky){
        for(int kx = 0; kx < width; ++kx){
            image[ky*width + kx] *= filter(kx, ky);
        }
    }

    // backward transform
    fft2.backward();

    // the resulting space domain image is contained
    // in fft2.spaceDomain()
    image = fft2.spaceDomain()
*/
class Fft2
{
public:
    typedef std::complex<double> element;

    // transformation's direction
    enum Direction {
        None, Forward, Backward, Both
    };

    explicit Fft2(int height, int width, Direction dir = None);
    ~Fft2(){ this->invalidate(); }

    // Before performing forward() or backward(),
    // Fft2 must be once prepare()'d.
    void prepare(Direction dir);

    // Forward transform (space -> freq)
    // Input is read from this->spaceDomain(): the values will be DESTROYED.
    // Output is written in this->freqDomain()
    void forward();

    // Backward transform (freq -> space)
    // Input is read from this->freqDomain(): the values will be DESTROYED.
    // Output is written in this->spaceDomain()
    void backward();

    // data accessors
    element* spaceDomain() { return space_; }
    element* const spaceDomain() const { return space_; }

    element* freqDomain() { return freq_; }
    element const* freqDomain() const { return freq_; }

private:
    Fft2(Fft2 const&);
    void operator=(Fft2 const&);

    void invalidate();

    element* space_; // space-domain data (row-major)
    element* freq_ ; // frequency-domain data (row-major)

    void* plan_forward_;  // fftw_plan (forward)
    void* plan_backward_; // fftw_plan (backward)

    int height_; // height of the image
    int width_ ; // width of the image
};

}}}} // lsst::meas::extensions::shapeZ08
#endif //ge93cecf0_2c00_40d7_ad6e_518123ea8be8
