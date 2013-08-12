
/******************************************************************************

  Common operations on images.

  Copyright (c) 2011, 2012
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
#include <functional>
#include <boost/cstdint.hpp>
#include <boost/function.hpp>

#include "bo/core/raw_image_2d.hpp"
#include "bo/math/functions.hpp"

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

// Functor returning either a given value or a substituion in case provided comparison
// operator returns true for value and threshold. Note that comparison operator is
// called like comp(value, threshold), e.g. in case of std::less this means
// value < threshold.
template <typename ValType>
struct threshold_value: public std::unary_function<ValType, ValType>
{
    typedef boost::function<bool (ValType, ValType)> ComparisonOperation;

    explicit
    threshold_value(const ValType& threshold, const ValType& substitution,
                    ComparisonOperation comparison):
        threshold_(threshold), substitution_(substitution), comp_(comparison)
    { }

    ValType operator()(const ValType& value) const
    {
        return
            (comp_(value, threshold_) ? substitution_ : value);
    }

    ValType threshold_;
    ValType substitution_;
    ComparisonOperation comp_;
};

// Helpers for normalization. By default all possible values of T are used. For real
// images (float, double) normalization is done in [0..1]
template <typename T> inline
T get_normalization_factor(T max_val)
{
    return (std::numeric_limits<T>::max() / max_val);
}

template <> inline
float get_normalization_factor<float>(float max_val)
{
    return (1.f / max_val);
}

template <> inline
double get_normalization_factor<double>(double max_val)
{
    return (1. / max_val);
}

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

// Returns a thresholded image based on given image, threshold, substitution value
// and probe function.
template <typename T, typename Comp>
RawImage2D<T> threshold(const RawImage2D<T>& image, const T& threshold,
                        const T& substitution, Comp comp)
{

    RawImage2D<T> thresholded(image.width(), image.height());
    std::transform(image.data(), image.data() + image.size(), thresholded.data(),
                   detail::threshold_value<T>(threshold, substitution, comp));

    return thresholded;
}

// Simple aggregation functions.
template <typename T>
T min_element(const RawImage2D<T>& image)
{
    return (*std::min_element(image.data(), image.data() + image.size()));
}

template <typename T>
T max_element(const RawImage2D<T>& image)
{
    return (*std::max_element(image.data(), image.data() + image.size()));
}

// Converts an image with integral pixel type to the real image on [0..1].
template <typename SourceT, typename DestT>
RawImage2D<DestT> convert_image(const RawImage2D<SourceT>& source_image)
{
    // Allocate empty image.
    std::size_t width = source_image.width();
    std::size_t height = source_image.height();
    bo::RawImage2D<DestT> converted_image(width, height);

    // Perform per element conversion.
    DestT factor = DestT(1) / std::numeric_limits<SourceT>::max();
    for (std::size_t row = 0; row < height; ++row)
        for (std::size_t col = 0; col < width; ++col)
            converted_image(col, row) = factor * source_image(col, row);

    return converted_image;
}

// Normalizes image to [0..max{T}] or to [0..1] depending on type of T.
template <typename T>
void normalize(RawImage2D<T>& image)
{
    T max_val = max_element(image);

    if (false == bo::math::check_small(max_val))
    {
        // Multiplication should be faster than division.
        T norm_factor = detail::get_normalization_factor(max_val);
        std::transform(image.data(), image.data() + image.size(), image.data(),
                std::bind2nd(std::multiplies<T>(), norm_factor));
    }
    else
    {
        throw std::logic_error("Normalization of the zero image is meaningless.");
    }
}

} // namespace bo

#endif // IMAGE_OPERATIONS_HPP_38A3CBD4_FAFD_404C_AC16_FC73B68B959B_
