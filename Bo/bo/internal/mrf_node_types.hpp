
/******************************************************************************

  mrf_node_types.hpp, v 0.1.2 2012.09.12

  Basic node type (class lables) and several ready-to-use types for MRF models.

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

#ifndef MRF_NODE_TYPES_HPP_79AFECA6_E8DF_47FE_952F_9EFE4097369E_
#define MRF_NODE_TYPES_HPP_79AFECA6_E8DF_47FE_952F_9EFE4097369E_

#include <ctime>
#include <vector>
#include <boost/cstdint.hpp>
#include <boost/random.hpp>

namespace bo {

// A class representing possible values for a given type. It can be used to reduce
// the range of a built-in type or to specify the set of possible values for a
// custom type. Method next() should always return a correct value (e.g. cyclically).
template<typename NodeType>
struct TypePossibleValues
{
    virtual void reset() = 0;
    virtual std::size_t count() = 0;
    virtual NodeType next() = 0;
    virtual NodeType random() = 0;
    virtual ~TypePossibleValues()
    { }
};

template <typename NodeType>
class RealFiniteSet: public TypePossibleValues<NodeType>
{
public:
    typedef boost::uniform_int<std::size_t> UintDistribution;
    typedef boost::variate_generator<boost::mt19937, UintDistribution> Generator;

    RealFiniteSet(std::size_t colour_depth):
        rng_(boost::mt19937(static_cast<boost::uint32_t>(std::time(NULL))),
             UintDistribution(0, colour_depth - 1))
    {
        // Construct a collection of possible values.
        values_.reserve(colour_depth);
        NodeType maxval = static_cast<NodeType>(colour_depth - 1);

        while (colour_depth-- > 0)
            values_.push_back(static_cast<NodeType>(colour_depth) / maxval);

        // Reset iterator state.
        reset();
    }

    virtual void reset()
    { next_index_ = 0; }

    virtual std::size_t count()
    { return values_.size(); }

    // Cyclic iterator over possible values.
    virtual NodeType next()
    { return values_[next_index_++ % count()]; }

    virtual NodeType random()
    { return values_[rng_()]; }

    virtual ~RealFiniteSet()
    { }

private:
    std::vector<NodeType> values_;
    std::size_t next_index_;
    Generator rng_;
};

} // namespace bo

#endif // MRF_NODE_TYPES_HPP_79AFECA6_E8DF_47FE_952F_9EFE4097369E_
