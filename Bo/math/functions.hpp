
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

#ifndef FUNCTIONS_HPP_AE638ADC_E4B3_11E2_9372_D8129A890B3C
#define FUNCTIONS_HPP_AE638ADC_E4B3_11E2_9372_D8129A890B3C

#include <cmath>

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

template <typename RealType> inline
bool check_small(const RealType& arg, const RealType& EPSILON = 1e-6f)
{
    return (arg < EPSILON) && (arg > -EPSILON);
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

} // namespace math
} // namespace bo

#endif // FUNCTIONS_HPP_AE638ADC_E4B3_11E2_9372_D8129A890B3C
