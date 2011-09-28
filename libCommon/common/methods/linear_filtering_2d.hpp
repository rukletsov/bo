
/******************************************************************************

    linear_filtering_2d.hpp, v 1.0.0 2011.09.26

    Methods and algorithms for applying linear filters for 2D images.

    Copyright (c) 2011
    Alexander Rukletsov <rukletsov@gmail.com>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    1.	Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
    2.	Redistributions in binary form must reproduce the above copyright
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

#ifndef LINEAR_FILTERING_2D_HPP_5C1DFBBA_B8EF_4105_AD4C_2C3EDE04C862_
#define LINEAR_FILTERING_2D_HPP_5C1DFBBA_B8EF_4105_AD4C_2C3EDE04C862_

#include <cmath>

#include "common/raw_image_2d.hpp"
#include "common/extended_math.hpp"

namespace common {
namespace methods {

namespace detail {

// Returns Sobel column kernel.
template <typename ValType>
common::RawImage2D<ValType> get_sobel_3x3_col_kernel()
{
    common::RawImage2D<ValType> kernel(3, 3);

    kernel(1, 0) = kernel(1, 1) = kernel(1, 2) = ValType(0);
    kernel(0, 0) = kernel(0, 2) = ValType(-1);
    kernel(2, 0) = kernel(2, 2) = ValType(1);
    kernel(0, 1) = ValType(-2);
    kernel(2, 1) = ValType(2);

    return kernel;
}

// Returns Sobel row kernel.
template <typename ValType>
common::RawImage2D<ValType> get_sobel_3x3_row_kernel()
{
    common::RawImage2D<ValType> kernel(3, 3);

    kernel(0, 1) = kernel(1, 1) = kernel(2, 1) = ValType(0);
    kernel(0, 0) = kernel(2, 0) = ValType(-1);
    kernel(0, 2) = kernel(2, 2) = ValType(1);
    kernel(1, 0) = ValType(-2);
    kernel(1, 2) = ValType(2);

    return kernel;
}

} // namespace detail

// Applies linear filter to a given image (convolve image with kernel). Note that
// image pixel type and kernel pixel type can be different. Current implementation
// casts kernel pixel type to image's and the performs multiplication.
template <typename ImValType, typename KerValType>
common::RawImage2D<ImValType> linear_filter_2d(const RawImage2D<ImValType> image,
    const RawImage2D<KerValType> kernel)
{
    // Cache image's and kernel's dimensions.
    std::size_t image_width = image.width();
    std::size_t image_height = image.height();
    std::size_t kernel_width = kernel.width();
    std::size_t kernel_height = kernel.height();

    // Kernel's dimensions should be odd.
    if ((kernel_width % 2 == 0) || (kernel_height % 2 == 0))
    {
        BOOST_ASSERT(false && "Kernel's dimensions are not odd.");
        throw std::logic_error("Linear filter cannot be applied with the given"
                               "kernel because kernel's dimensions are not odd.");
    }

    // Compute the distances from kernel's central point.
    std::ptrdiff_t kernel_half_width = (kernel_width - 1) / 2;
    std::ptrdiff_t kernel_half_height = (kernel_height - 1) / 2;

    // Create resulting image and run filtering. It is easier to address kernel
    // using center-based indexation. If the considered kernel element is outside
    // image's border, skip it (pretend image pixel is 0 there).
    common::RawImage2D<ImValType> filtered_image(image_width, image_height);

    for (std::size_t col = 0; col < image_width; ++col) {
        for (std::size_t row = 0; row < image.height(); ++row) {
            filtered_image(col, row) = ImValType(0);

            for (std::ptrdiff_t ker_col = - kernel_half_width;
                             ker_col <= kernel_half_width; ++ker_col) {
                for (std::ptrdiff_t ker_row = - kernel_half_height;
                                 ker_row <= kernel_half_height; ++ker_row) {
                    if ((col + ker_col >= 0) && (col + ker_col < image_width) &&
                        (row + ker_row >= 0) && (row + ker_row < image_height))
                        filtered_image(col, row) += image(col + ker_col, row + ker_row)
                            * ImValType(kernel(ker_col + kernel_half_width,
                                               ker_row + kernel_half_height));
    }   }   }   }

    return filtered_image;
}

// Returns image's edge map based on Sobel operator. Each pixel of an edge map is
// sqrt(Gx^2 + Gy*2), where Gx and Gy are results of applying filter with Sobel column
// and row kernels respectively.
template <typename ValType>
common::RawImage2D<ValType> sobel_3x3(const RawImage2D<ValType> image)
{
    common::RawImage2D<ValType> sobel_col = linear_filter_2d(image,
        detail::get_sobel_3x3_col_kernel<ValType>());

    common::RawImage2D<ValType> sobel_row = linear_filter_2d(image,
        detail::get_sobel_3x3_row_kernel<ValType>());

    std::size_t image_width = image.width();
    std::size_t image_height = image.height();
    common::RawImage2D<ValType> filtered_image(image_width, image_height);

    for (std::size_t col = 0; col < image_width; ++col) {
        for (std::size_t row = 0; row < image.height(); ++row) {
            filtered_image(col, row) = sqrt(square(sobel_col(col, row)) +
                                            square(sobel_row(col, row)));
    }   }

    return filtered_image;
}

} // namespace methods
} // namespace common

#endif // LINEAR_FILTERING_2D_HPP_5C1DFBBA_B8EF_4105_AD4C_2C3EDE04C862_
