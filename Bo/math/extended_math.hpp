
/******************************************************************************

  Extension of the standard <cmath> header.

  Copyright (c) 2011 - 2013
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

#ifndef EXTENDED_MATH_HPP_5E8C7161_2D47_4FF0_974A_19599004895C_
#define EXTENDED_MATH_HPP_5E8C7161_2D47_4FF0_974A_19599004895C_

#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <boost/math/tr1.hpp>

// Suppress annoying MSVC's C4244 conversion warnings popping out from boost::accumulators
// library (boost::numeric::functional namespace) mainly due to division operation
// in mean computation. Suppress MSVC's C4512 warning for boost auxiliary classes as well.
#ifdef _MSC_VER
#   pragma warning(push)
#   pragma warning(disable:4244)
#   pragma warning(disable:4512)
#endif // _MSC_VER

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#ifdef _MSC_VER
#   pragma warning(pop)
#endif // _MSC_VER

namespace bo {
namespace math {

template <typename T> inline
T square(const T& arg)
{
    return (arg * arg);
}

template <typename T> inline
T cube(const T& arg)
{
    return (arg * arg * arg);
}

// Returns "2 ln(cosh(t))".
template <typename RealType> inline
RealType fi_gr(RealType value)
{
    return
        RealType(2) * std::log(std::exp(value) + std::exp(RealType(0) - value) / RealType(2));
}

// Returns "(t^2) / (1+t^2)".
template <typename RealType> inline
RealType fi_gm(RealType value)
{
    return
        (value * value) / (RealType(1) + value * value);
}

// Returns the mean of the given samples.
template <typename SampleType>
SampleType mean(const std::vector<SampleType>& data)
{
    // The mean of an empty set is meaningless.
    if (data.size() == 0)
        throw std::logic_error("The mean of an empty set is meaningless.");

    // Initialize boost accumulator.
    namespace accs = boost::accumulators;
    typedef accs::accumulator_set<SampleType, accs::stats<accs::tag::mean> > Acc;
    Acc acc;

    // Fill accumulator with data.
    acc = std::for_each(data.begin(), data.end(), acc);

    // Request and return mean.
    return accs::mean(acc);
}

} // namespace math
} // namespace bo

#endif // EXTENDED_MATH_HPP_5E8C7161_2D47_4FF0_974A_19599004895C_