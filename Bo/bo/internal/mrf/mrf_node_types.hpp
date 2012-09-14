
/******************************************************************************

  mrf_node_types.hpp, v 0.2.3 2012.09.14

  Non-trivial node types (class lables) for MRF models.

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

#ifndef MRF_NODE_TYPES_HPP_6CF0A0E1_8AE9_4951_A785_9339B90FB976_
#define MRF_NODE_TYPES_HPP_6CF0A0E1_8AE9_4951_A785_9339B90FB976_

#include <boost/tuple/tuple.hpp>
#include <boost/operators.hpp>
#include <boost/shared_ptr.hpp>

namespace bo {

template <typename RealType>
class GammaDistrClasses: boost::equality_comparable1<GammaDistrClasses<RealType> >
{
public:
    // Class label and a set of Gamma distribution parameters (k, theta, a), where
    // a = -ln(G(k)) + k ln(theta) and G(t) is the Gamma function.
    typedef boost::tuples::tuple<int, RealType, RealType, RealType> ClassParams;
    typedef boost::tuples::tuple<RealType, RealType> GammaParamsPair;
    typedef boost::shared_ptr<ClassParams> ClassParamsPtr;

public:
    GammaDistrClasses(ClassParamsPtr class_params): class_params_(class_params)
    { }

    //  Generated copy c-tor, d-tor and assignment operator are fine.

    // Accessors for class label and class parameters.
    int label() const
    { return class_params_->get<0>(); }

    RealType k() const
    { return class_params_->get<1>(); }

    RealType theta() const
    { return class_params_->get<2>(); }

    RealType a() const
    { return class_params_->get<3>(); }

    RealType mean() const
    { return (k() * theta()); }

    // Returns the difference between class labels. This may or may not tell us how
    // "far" the classes are using some metric. The meaning of this operation depends
    // on the actual design on classes and their parameters.
    RealType operator-(const GammaDistrClasses<RealType>& other) const
    { return RealType(label() - other.label()); }

    // Classes are equal when their parameters are equal.
    bool operator==(const GammaDistrClasses<RealType>& other) const
    { return (class_params_ == other.class_params_); }

protected:
    ClassParamsPtr class_params_;
};

} // namespace bo

#endif // MRF_NODE_TYPES_HPP_6CF0A0E1_8AE9_4951_A785_9339B90FB976_
