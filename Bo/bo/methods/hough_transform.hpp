/******************************************************************************

  Implementation of the original Hough transform for 2D line tracking.

  Reference papers:

  P.V.C. Hough, "Machine Analysis of Bubble Chamber Pictures",
  Proc. Int. Conf. High Energy Accelerators and Instrumentation, 1959.
  (U.S. Patent 3,069,654)

  Duda, R. O. and P. E. Hart, "Use of the Hough Transformation to Detect
  Lines and Curves in Pictures," Comm. ACM, Vol. 15, pp. 11–15, 1972.

  Copyright (c) 2013
  Dzmitry Hlindzich <hlindzich@gmail.com>
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1.  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
  2.  Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
  SUCH DAMAGE.

*******************************************************************************/

#ifndef HOUGH_TRANSFORM_HPP_F8fEE77D_37A4_4DFB_83ED_CF832D40772D_
#define HOUGH_TRANSFORM_HPP_F8fEE77D_37A4_4DFB_83ED_CF832D40772D_

#include "bo/raw_image_2d.hpp"

#include <boost/assert.hpp>
#include <boost/math/constants/constants.hpp>

#include <cmath>

namespace bo {
namespace methods {


template <typename RealType>
class HoughTransform
{
public:

    typedef bo::RawImage2D<float> Image;

    HoughTransform(const Image &image): input_image_(image),
        max_vote_(0)
    {
        std::size_t w = image.width();
        std::size_t h = image.height();

        BOOST_ASSERT(w > 0 && h > 0);

        max_rho_ = std::sqrt(w * w + h * h);
        min_rho_ = -w;

        max_theta_ = boost::math::constants::pi<RealType>();
        min_theta_ = RealType(0);
    }

    // Returns the Hough transform (in the rho-theta parametrization) of the points
    // from the input image whose values exceed the given threshold. Represents the
    // result as a 2D image of size (output_width * output_height).
    // The resulted image covers the following rectangle in the rho-theta space:
    // width = values of theta [0, Pi], height = values of rho [-w, sqrt(w^2 + h^2)],
    // where w and h are the width and the height of the input image correspondingly.
    Image compute(RealType threshold, std::size_t output_width, std::size_t output_height)
    {
        BOOST_ASSERT(output_height > 0 && output_width > 0);

        max_vote_ = 0;

        // Compute the output image resolution.
        theta_scaling_ = static_cast<RealType>(output_width) / (max_theta_ - min_theta_);
        rho_scaling_ = static_cast<RealType>(output_height - 1) / (max_rho_ - min_rho_);

        // Allocate the output image.
        bo::RawImage2D<RealType> hough(output_width, output_height, 0);

        std::size_t w = input_image_.width();
        std::size_t h = input_image_.height();

        // Compute sinusoids in the rho-theta space for each point of the input image
        // whose value is greater that the threshold.
        for (std::size_t i = 0; i < w; ++i)
            for (std::size_t j = 0; j < h; ++j)
            {
                if (input_image_(i,j) > threshold)
                {
                    // Accumulate the sinusoid.
                    for (std::size_t x = 0; x < output_width; ++x)
                    {
                        RealType theta = min_theta_ + x / theta_scaling_;
                        RealType rho = i * std::cos(theta) + j * std::sin(theta);

                        std::size_t y = static_cast<std::size_t>(rho_scaling_ * (rho - min_rho_));

                        RealType vote = ++hough(x, y);

                        if(max_vote_ < vote)
                        {
                            max_vote_ = vote;
                        }
                    }
                }
            }

        return hough;
    }

    Image reconstruct_lines(RealType threshold, RealType quantity, std::size_t accu_width = 1024,
                            std::size_t accu_height = 1024)
    {
        const RealType kEpsilon(0.01);

        BOOST_ASSERT(quantity >= RealType(0) && quantity <= RealType(1));

        // Compute the accumulator.
        Image hough = compute(threshold, accu_width, accu_height);

        RealType accu_threshold = max_vote_ * quantity;

        // Allocate the output image.
        std::size_t w = input_image_.width();
        std::size_t h = input_image_.height();
        Image line_image(w, h, 0);

        // Draw the lines.
        for (std::size_t i = 0; i < accu_width; ++i)
            for (std::size_t j = 0; j < accu_height; ++j)
            {
                if (hough(i, j) >= accu_threshold)
                {
                    RealType theta = min_theta_ + i / theta_scaling_;
                    RealType rho = min_rho_ + j / rho_scaling_;

                    RealType sint = std::sin(theta);
                    RealType cost = std::cos(theta);

                    RealType t1, t2;

                    if (std::abs(sint) < kEpsilon)
                    {
                        t1 = 0;
                        t2 = h;
                    }
                    else if (std::abs(cost) < kEpsilon)
                    {
                        t1 = -w;
                        t2 = 0;
                    }
                    else
                    {
                        RealType tc = rho * cost / sint;
                        t1 = t2 = tc;

                        tc = -rho * sint / cost;
                        t1 = std::min(t1, tc);
                        t2 = std::max(t2, tc);

                        tc = (rho * cost - w) / sint;
                        t1 = std::min(t1, tc);
                        t2 = std::max(t2, tc);

                        tc = (h - rho * sint) / cost;
                        t1 = std::min(t1, tc);
                        t2 = std::max(t2, tc);
                    }

                    RealType delta_t = RealType(1);
                    RealType t = t1;

                    while (t < t2)
                    {
                        std::size_t x = rho * cost - t * sint;
                        std::size_t y = rho * sint + t * cost;

                        if (x >= 0 && x < w && y >= 0 && y < h)
                        {
                            line_image(x, y) = 1;
                        }

                        t += delta_t;
                    }
                }
            }

        return line_image;
    }

protected:

    Image input_image_;

    RealType max_rho_;
    RealType min_rho_;
    RealType max_theta_;
    RealType min_theta_;

    RealType theta_scaling_;
    RealType rho_scaling_;


    RealType max_vote_;
};

} // namespace methods
} // namespace bo

#endif // HOUGH_TRANSFORM_HPP_F8fEE77D_37A4_4DFB_83ED_CF832D40772D_
