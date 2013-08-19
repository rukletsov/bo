
/******************************************************************************

  Various likelihood energy functions for MRF models.

  Copyright (c) 2012, 2013
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

#ifndef LIKELIHOOD_FUNCTIONS_HPP_BD4FA567_7D8B_4357_A2AC_89BEFE621679
#define LIKELIHOOD_FUNCTIONS_HPP_BD4FA567_7D8B_4357_A2AC_89BEFE621679

#include "bo/math/functions.hpp"

namespace bo {
namespace mrf {

// Base class for likelihood energy functions.
template <typename NodeType, typename DataType, typename RealType>
struct GenericLikelihood
{
    GenericLikelihood(RealType response_weight): multiplier_(response_weight)
    { }

    virtual RealType operator()(DataType observ_val, NodeType configur_val) const = 0;

    virtual ~GenericLikelihood()
    { }

protected:
    // Weighting coefficient for the function response.
    RealType multiplier_;
};

// Minus operator for NodeType should accept DataType as a parameter and return RealType.
template <typename NodeType, typename DataType, typename RealType>
struct GaussSimpleLikelihood: public GenericLikelihood<NodeType, DataType, RealType>
{
    typedef GenericLikelihood<NodeType, DataType, RealType> BaseType;

    GaussSimpleLikelihood(RealType response_weight): BaseType(response_weight)
    { }

    RealType operator()(DataType observ_val, NodeType configur_val) const
    {
        return
            this->multiplier_ * math::square(configur_val - observ_val);
    }
};

// NodeType should provide accessors to the parameters of Gauss distribution for the
// corresponding configuration value (class label). This includes .mu() and .sigma()
// for mean and standard deviation respectively and .a() for an additional precomputed
// item, depending only on sigma: ln(sigma sqrt(2 pi)).
template <typename NodeType, typename DataType, typename RealType>
struct GaussLikelihood: public GenericLikelihood<NodeType, DataType, RealType>
{
    typedef GenericLikelihood<NodeType, DataType, RealType> BaseType;

    GaussLikelihood(RealType response_weight): BaseType(response_weight)
    { }

    RealType operator()(DataType observ_val, NodeType configur_val) const
    {
        return
            this->multiplier_ * (RealType(0.5) * math::square(configur_val.mean() - observ_val)
                    / configur_val.sigma() + configur_val.a());
    }
};

// NodeType should provide accessors to the parameters of Gamma distribution for the
// corresponding configuration value (class label). This includes .k() and .theta()
// for shape and scale respectively and .a() for an additional item, depending
// only on k and theta (and therefore precomputed): ln(G(k)) + k ln(theta).
template <typename NodeType, typename DataType, typename RealType>
struct GammaLikelihood: public GenericLikelihood<NodeType, DataType, RealType>
{
    typedef GenericLikelihood<NodeType, DataType, RealType> BaseType;

    GammaLikelihood(RealType response_weight): BaseType(response_weight)
    { }

    RealType operator()(DataType observ_val, NodeType configur_val) const
    {
        return
            this->multiplier_ * ((1 - configur_val.k()) * std::log(observ_val) +
                          observ_val / configur_val.theta() + configur_val.a());
    }
};

} // namespace mrf
} // namespace bo

#endif // LIKELIHOOD_FUNCTIONS_HPP_BD4FA567_7D8B_4357_A2AC_89BEFE621679
