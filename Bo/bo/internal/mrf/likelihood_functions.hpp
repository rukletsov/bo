
/******************************************************************************

  likelihood_functions.hpp, v 0.1.4 2012.09.14

  Various likelihood energy functions for MRF models.

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

#ifndef LIKELIHOOD_FUNCTIONS_HPP_BD4FA567_7D8B_4357_A2AC_89BEFE621679_
#define LIKELIHOOD_FUNCTIONS_HPP_BD4FA567_7D8B_4357_A2AC_89BEFE621679_

#include "bo/extended_math.hpp"

namespace bo {

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
// only on k and theta (and therefore precomputed): -ln(G(k)) + k ln(theta).
template <typename DataType, typename NodeType, typename RealType>
struct GammaLikelihood: public GenericLikelihood<DataType, NodeType, RealType>
{
    GammaLikelihood(RealType response_weight): GenericLikelihood(response_weight)
    { }

    RealType operator()(DataType observ_val, NodeType configur_val) const
    {
        return
            multiplier * ((configur_val.k() - 1) * std::log(observ_val) -
                          observ_val / configur_val.theta() + configur_val.a());
    }
};

} // namespace bo

#endif // LIKELIHOOD_FUNCTIONS_HPP_BD4FA567_7D8B_4357_A2AC_89BEFE621679_