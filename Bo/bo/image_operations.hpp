
/******************************************************************************

  image_operations.hpp, v 1.0.0 2011.09.28

  Common operations on images.

  Copyright (c) 2011
  Alexander Rukletsov <rukletsov@gmail.com>
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

#ifndef IMAGE_OPERATIONS_HPP_38A3CBD4_FAFD_404C_AC16_FC73B68B959B_
#define IMAGE_OPERATIONS_HPP_38A3CBD4_FAFD_404C_AC16_FC73B68B959B_

#include <algorithm>
#include <limits>
#include <boost/cstdint.hpp>

#include "bo/raw_image_2d.hpp"

namespace bo {

// Class free functions. These functions strongly depend on pixel type and therefore
// defined outside the class.

namespace detail {

// Functor returning inverted brightness value for a given one. It uses the formula:
// "inverted_brightness = max_brightness - current". Note that maximum brightness is
// not necessarily maximum possible value of a given type (e.g. normalized images).
template <typename ValType>
struct invert_brightness: public std::unary_function<ValType, ValType>
{
    explicit
    invert_brightness(const ValType& max_value): max_value_(max_value)
    { }

    ValType operator()(const ValType& value) const
    {
        return (max_value_ - value);
    }

    ValType max_value_;
};

} // namespace detail


// Returns an inverted image with arbitrary pixels. Supposes pixel values are in
// [0 .. std::numeric_limits<T>::max()]
template <typename T>
RawImage2D<T> invert(const RawImage2D<T>& image)
{
    RawImage2D<T> inverted(image.width(), image.height());
    std::transform(image.data(), image.data() + image.size(), inverted.data(),
                   detail::invert_brightness<T>(std::numeric_limits<T>::max()));

    return inverted;
}

// Returns an inverted image with float pixels. Supposes pixel values are in [0..1].
template <>
RawImage2D<float> invert<float>(const RawImage2D<float>& image)
{
    RawImage2D<float> inverted(image.width(), image.height());
    std::transform(image.data(), image.data() + image.size(), inverted.data(),
                   detail::invert_brightness<float>(1.f));

    return inverted;
}

// Returns an inverted image with double pixels. Supposes pixel values are in [0..1].
template <>
RawImage2D<double> invert<double>(const RawImage2D<double>& image)
{
    RawImage2D<double> inverted(image.width(), image.height());
    std::transform(image.data(), image.data() + image.size(), inverted.data(),
                   detail::invert_brightness<double>(1.));

    return inverted;
}

} // namespace bo

#endif // IMAGE_OPERATIONS_HPP_38A3CBD4_FAFD_404C_AC16_FC73B68B959B_
