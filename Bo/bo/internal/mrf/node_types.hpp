
/******************************************************************************

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

#ifndef NODE_TYPES_HPP_6CF0A0E1_8AE9_4951_A785_9339B90FB976_
#define NODE_TYPES_HPP_6CF0A0E1_8AE9_4951_A785_9339B90FB976_

#include <cmath>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/operators.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace bo {

template <typename P>
class ParametricNodeType: boost::equality_comparable1<ParametricNodeType<P> >
{
public:
    typedef P ClassParams;
    typedef boost::shared_ptr<P> ClassParamsPtr;

    ParametricNodeType(ClassParamsPtr class_params): class_params_(class_params)
    { }

    //  Generated copy c-tor and assignment operator are fine.

    virtual ~ParametricNodeType()
    { }

    // Classes are equal when their parameters are equal.
    virtual bool operator==(const ParametricNodeType<P>& other) const
    { return (class_params_ == other.class_params_); }

protected:
    ClassParamsPtr class_params_;
};

// A class representing a value (class label) with associated Gamma distribution
// parameters. Parameters structure is represented as a tuple with 4 elements:
// class label and a set of Gamma distribution parameters (k, theta, a), where
// a = ln(G(k)) + k ln(theta) and G(t) is the Gamma function.
template <typename RealType>
class GammaDistrClasses: public ParametricNodeType<boost::tuples::tuple<int, RealType, RealType, RealType> >
{
public:
    typedef boost::tuples::tuple<RealType, RealType> GammaParamsPair;

public:
    GammaDistrClasses(ClassParamsPtr class_params): ParametricNodeType(class_params)
    { }

    static GammaDistrClasses CreateInstance(int idx, RealType k, RealType theta)
    { return (GammaDistrClasses(ClassParamsPtr(new ClassParams(idx, k, theta, compute_a(k, theta))))); }

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

    static RealType compute_a(RealType k, RealType theta)
    { return (boost::math::lgamma(k) + k * std::log(theta)); }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const GammaDistrClasses<T>& obj)
{
    os << obj.label() << " (k: " << obj.k() << ", theta: " << obj.theta() << ")";
    return os;
}

} // namespace bo

#endif // NODE_TYPES_HPP_6CF0A0E1_8AE9_4951_A785_9339B90FB976_
