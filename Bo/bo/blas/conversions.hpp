
/******************************************************************************

  conversions.hpp, v 1.0.0 2012.12.11

  Conversions to and from from boost::numeric::ublas types.

  Copyright (c) 2012
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

#ifndef CONVERSIONS_HPP_0DE10BDD_4801_49EC_B589_F120F69732BA_
#define CONVERSIONS_HPP_0DE10BDD_4801_49EC_B589_F120F69732BA_

#include <vector>
#include <algorithm>

#include "bo/blas/blas.hpp"
#include "bo/vector.hpp"

namespace bo {
namespace blas {

using namespace boost::numeric::ublas;

template <typename T, std::size_t N>
Vector<T, N> to_bo_vector(const bounded_vector<T, N>& from)
{
    Vector<T, N> to(from.data(), N);
    return to;
}

template <typename T>
std::vector<T> to_std_vector(const vector<T>& from)
{
    std::vector<T> to(from.begin(), from.end());
    return to;
}

template <typename T, std::size_t N>
bounded_vector<T, N> from_bo_vector(const Vector<T, N>& from)
{
    bounded_vector<T, N> to;
    for (std::size_t idx = 0; idx < N; ++idx)
        to(idx) = from[idx];

    return to;
}

template <typename T>
vector<T> from_std_vector(const std::vector<T>& from)
{
    vector<T> to(from.size(), T(0));
    std::copy(from.begin(), from.end(), to.data().begin());
    return to;
}

} // namespace blas
} // namespace bo

#endif CONVERSIONS_HPP_0DE10BDD_4801_49EC_B589_F120F69732BA_
