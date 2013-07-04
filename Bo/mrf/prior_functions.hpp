
/******************************************************************************

  Various prior energy functions for MRF models.

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

#ifndef PRIOR_FUNCTIONS_HPP_9920258C_F29D_4179_93C1_49ED114BC299_
#define PRIOR_FUNCTIONS_HPP_9920258C_F29D_4179_93C1_49ED114BC299_

#include "bo/math/extended_math.hpp"

namespace bo {
namespace mrf {

// Base class for prior energy functions.
template <typename NodeType, typename RealType>
struct GenericPrior
{
    GenericPrior(RealType response_weight): multiplier(response_weight)
    { }

    virtual RealType operator()(NodeType arg1, NodeType arg2) const = 0;

    virtual ~GenericPrior()
    { }

protected:
    // Weighting coefficient for the function response.
    RealType multiplier;
};

// Standard smoothness prior energy on two-node clique. Minus operator for NodeType
// should return RealType.
template <typename NodeType, typename RealType>
struct SmoothnessPrior: public GenericPrior<NodeType, RealType>
{
    SmoothnessPrior(RealType response_weight): GenericPrior(response_weight)
    { }

    virtual RealType operator()(NodeType arg1, NodeType arg2) const
    {
        return multiplier * math::square(arg1 - arg2) / RealType(2);
    }

    virtual ~SmoothnessPrior()
    { }
};

// Smoothness prior with additional edge preserving functional. Tau controls the
// threshold penalty for node values' differency. Minus operator for NodeType should
// return RealType.
template <typename NodeType, typename RealType>
struct SmoothingWithEdgesPrior: public SmoothnessPrior<NodeType, RealType>
{
    SmoothingWithEdgesPrior(RealType response_weight, RealType edge_weight,
        RealType tau): SmoothnessPrior(response_weight),
        edge_coefficient(edge_weight), thres_border(tau)
    { }

    virtual RealType operator()(NodeType arg1, NodeType arg2) const
    {
        // Note that base class already multiplies in the response weight.
        return
            (SmoothnessPrior::operator ()(arg1, arg2) +
             multiplier * edge_coefficient * std::min(std::abs(arg1 - arg2), thres_border));
    }

    virtual ~SmoothingWithEdgesPrior()
    { }

protected:
    // Weighting coefficient for edge preserving functional.
    RealType edge_coefficient;
    RealType thres_border;
};

// Smoothness prior energy on two-node clique for node types, supporting mean() method.
template <typename NodeType, typename RealType>
struct MeanSmoothnessPrior: public GenericPrior<NodeType, RealType>
{
    MeanSmoothnessPrior(RealType response_weight): GenericPrior(response_weight)
    { }

    virtual RealType operator()(NodeType arg1, NodeType arg2) const
    {
        return multiplier * math::square(arg1.mean() - arg2.mean()) / RealType(2);
    }

    virtual ~MeanSmoothnessPrior()
    { }
};

} // namespace mrf
} // namespace bo

#endif // PRIOR_FUNCTIONS_HPP_9920258C_F29D_4179_93C1_49ED114BC299_
