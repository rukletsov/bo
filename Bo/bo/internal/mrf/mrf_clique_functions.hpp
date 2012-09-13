
/******************************************************************************

  mrf_clique_functions.hpp, v 0.1.2 2012.09.13

  Various likelihood and prior energy functions for MRF models.

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

#ifndef MRF_CLIQUE_FUNCTIONS_HPP_FAB983DB_6274_4A09_A6E2_A9732F63EB74_
#define MRF_CLIQUE_FUNCTIONS_HPP_FAB983DB_6274_4A09_A6E2_A9732F63EB74_

#include <functional>

#include "bo/extended_math.hpp"

namespace bo {

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
struct SmoothnessPriorEnergy: public GenericPrior<NodeType, RealType>
{
    SmoothnessPriorEnergy(RealType response_weight): GenericPrior(response_weight)
    { }

    virtual RealType operator()(NodeType arg1, NodeType arg2) const
    {
        return multiplier * square(arg1 - arg2) / RealType(2);
    }

    virtual ~SmoothnessPriorEnergy()
    { }
};

// Smoothness prior with additional edge preserving functional. Tau controls the
// threshold penalty for node values' differency. Minus operator for NodeType should
// return RealType.
template <typename NodeType, typename RealType>
struct SmoothingWithEdgesPriorEnergy: public SmoothnessPriorEnergy<NodeType, RealType>
{
    SmoothingWithEdgesPriorEnergy(RealType response_weight, RealType edge_weight,
        RealType tau): SmoothnessPriorEnergy(response_weight),
        edge_coefficient(edge_weight), thres_border(tau)
    { }

    virtual RealType operator()(NodeType arg1, NodeType arg2) const
    {
        // Note that base class already multiplies in the response weight.
        return
            (SmoothnessPriorEnergy::operator ()(arg1, arg2) +
             multiplier * edge_coefficient * std::min(std::abs(arg1 - arg2), thres_border));
    }

    virtual ~SmoothingWithEdgesPriorEnergy()
    { }

protected:
    // Weighting coefficient for edge preserving functional.
    RealType edge_coefficient;
    RealType thres_border;
};


// Base class for likelihood energy functions.
template <typename DataType, typename NodeType, typename RealType>
struct GenericLikelihood
{
    GenericLikelihood(RealType response_weight): multiplier(response_weight)
    { }

    virtual RealType operator()(DataType observ_val, NodeType configur_val) const = 0;

    virtual ~GenericLikelihood()
    { }

protected:
    // Weighting coefficient for the function response.
    RealType multiplier;
};

// Minus operator for NodeType should accept DataType as a parameter and return RealType.
template <typename DataType, typename NodeType, typename RealType>
struct GaussianLikelihood: public GenericLikelihood<DataType, NodeType, RealType>
{
    GaussianLikelihood(RealType response_weight): GenericLikelihood(response_weight)
    { }

    RealType operator()(DataType observ_val, NodeType configur_val) const
    {
        return
            multiplier * square(configur_val - observ_val);
    }
};

// NodeType should provide accessors to the parameters of Gamma distribution for the
// corresponding configuration value (class label). This includes .k() and .theta()
// for shape and scale respectively and .a() for the additional item, depending
// only on k and theta (and therefore precomputed): ln(G(k)) + k ln(theta).
template <typename DataType, typename NodeType, typename RealType>
struct GammaLikelihood: public GenericLikelihood<DataType, NodeType, RealType>
{
    GammaLikelihood(RealType response_weight): GenericLikelihood(response_weight)
    { }

    RealType operator()(DataType observ_val, NodeType configur_val) const
    {
        return
            multiplier * ((configur_val.k() - 1) * std::log(observ_val) -
                          observ_val / configur_val.theta() - configur_val.a());
    }
};

} // namespace bo

#endif // MRF_CLIQUE_FUNCTIONS_HPP_FAB983DB_6274_4A09_A6E2_A9732F63EB74_
