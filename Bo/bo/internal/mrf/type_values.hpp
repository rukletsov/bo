
/******************************************************************************

  type_values.hpp, v 0.2.3 2012.09.14

  Collection of classes specifying and providing access to possible values for
  types used as MRF nodes.

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

#ifndef TYPE_VALUES_HPP_79AFECA6_E8DF_47FE_952F_9EFE4097369E_
#define TYPE_VALUES_HPP_79AFECA6_E8DF_47FE_952F_9EFE4097369E_

#include <ctime>
#include <cmath>
#include <vector>
#include <boost/cstdint.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "bo/internal/mrf/node_types.hpp"

namespace bo {

// An abstract class representing possible values for a given type. Descendant of
// this class can be used to reduce the range of a built-in type or to specify the
// set of possible values for a custom type. Method next() should always return a
// correct value (e.g. cyclically).
template<typename NodeType>
struct TypeValues
{
    virtual void reset() = 0;
    virtual std::size_t count() const = 0;
    virtual NodeType next() = 0;
    virtual NodeType random() = 0;
    virtual ~TypeValues()
    { }
};

// An abstract class representing a finite set of values with cyclic iterator and
// predefined random numbers generator. The class is almost ready to use, except
// actual values of NodeType, which must be added in descendant's c-tor.
template<typename NodeType>
struct FiniteSetValues: public TypeValues<NodeType>
{
    typedef boost::uniform_int<std::size_t> UintDistribution;
    typedef boost::variate_generator<boost::mt19937, UintDistribution> Generator;

    virtual void reset()
    { next_index_ = 0; }

    virtual std::size_t count() const
    { return values_.size(); }

    // Cyclic iterator over possible values.
    virtual NodeType next()
    { return values_[next_index_++ % count()]; }

    virtual NodeType random()
    { return values_[rng_()]; }

    virtual ~FiniteSetValues()
    { }

protected:
    FiniteSetValues(std::size_t classes_count):
        rng_(boost::mt19937(static_cast<boost::uint32_t>(std::time(NULL))),
             UintDistribution(0, classes_count - 1))
    {
        // Reset iterator state and reserve memory for values.
        reset();
        values_.reserve(classes_count);
    }

protected:
    std::vector<NodeType> values_;
    std::size_t next_index_;
    Generator rng_;
};

// A class restricting [0, 1] real domain to a certain quantity of equally distant
// values from [0, 1].
template <typename NodeType>
struct RealFiniteSetValues: public FiniteSetValues<NodeType>
{
    // Constructs a set of possible values.
    RealFiniteSetValues(std::size_t classes_count): FiniteSetValues(classes_count)
    {
        NodeType maxval = static_cast<NodeType>(classes_count - 1);
        while (classes_count-- > 0)
            values_.push_back(static_cast<NodeType>(classes_count) / maxval);
    }

    virtual ~RealFiniteSetValues()
    { }
};

// A class representing a finite set of GammaDistrClasses values. Each value gets
// a class label from [0, classes_count-1] and a corresponding set of class
// parameters (see GammaDistrClasses for more information).
template <typename RealType>
struct GammaDistrClassesValues: public FiniteSetValues<GammaDistrClasses<RealType> >
{
    typedef GammaDistrClasses<RealType> NodeType;
    typedef typename NodeType::ClassParamsPtr ClassParamsPtr;
    typedef typename NodeType::ClassParams ClassParams;
    typedef typename NodeType::GammaParamsPair GammaParamsPair;
    typedef typename std::vector<GammaParamsPair> GammaParams;

    // Constructs a set of possible values with default parameters.
    GammaDistrClassesValues(std::size_t classes_count): FiniteSetValues(classes_count)
    {
        while (classes_count-- > 0)
        {
            ClassParamsPtr class_params(new ClassParams(classes_count));
            values_.push_back(NodeType(class_params));
        }

    }

    // Constructs a set of possible values from a collection of Gamma distribution
    // parameters (k, theta).
    GammaDistrClassesValues(GammaParams gamma_params): FiniteSetValues(classes_count)
    {
        std::size_t classes_count = gamma_params.size();
        for (GammaParams::const_iterator it = gamma_params.begin();
             it != gamma_params.end(); ++it)
        {
            --classes_count;
            RealType k = it->get<0>();
            RealType theta = it->get<1>();
            RealType a = k * std::log(theta) - boost::math::lgamma(k);
            ClassParamsPtr class_params(new ClassParams(classes_count, k, theta, a));
            values_.push_back(NodeType(class_params));
        }
    }

    virtual ~GammaDistrClassesValues()
    { }
};

} // namespace bo

#endif // TYPE_VALUES_HPP_79AFECA6_E8DF_47FE_952F_9EFE4097369E_
