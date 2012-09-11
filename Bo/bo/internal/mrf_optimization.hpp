
/******************************************************************************

  mrf_optimization.hpp, v 0.1.1 2012.09.11

  Various energy minimization algorithms for MRF models.

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

#ifndef MRF_OPTIMIZATION_HPP_40C9F0DC_7E18_4316_A594_63DAFD793CCA_
#define MRF_OPTIMIZATION_HPP_40C9F0DC_7E18_4316_A594_63DAFD793CCA_

#include <cstddef>
#include <limits>
#include <functional>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/random.hpp>

#include "bo/internal/mrf_2d.hpp"

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
    virtual ~TypePossibleValues()
    { }
};

// A class impementing iterated conditional modes minimization algorithm. Takes
// (and captures sole ownership over it) an istance of TypePossibleValues<NodeType>
// to sanple possible values during optimization.
template <typename NodeType, typename DataType, typename RealType>
class ICM2D: boost::noncopyable
{
public:
    typedef TypePossibleValues<NodeType> NodePossibleLabels;
    typedef boost::scoped_ptr<NodePossibleLabels> NodePossibleLabelsPtr;

    ICM2D(NodePossibleLabels* possible_values): values_(possible_values)
    { }

    void next_iteration(MRF2D<NodeType, DataType, RealType>& mrf)
    {
        for (std::size_t col = 0; col < mrf.width(); ++col) {
            for (std::size_t row = 0; row < mrf.height(); ++row) {
                RealType min_energy = std::numeric_limits<RealType>::max();
                NodeType min_value = values_->next();

                // Reset the possible labels so we can iterate through the whole
                // variety only once.
                values_->reset();
                std::size_t count = values_->count();

                while (count--)
                {
                    NodeType value = values_->next();
                    RealType energy = mrf.compute_local_energy(value, col, row);

                    if (energy < min_energy)
                    {
                        min_energy = energy;
                        min_value = value;
                    }
                }

                // Apply the "best" value.
                mrf(col, row) = min_value;
        }   }
    }

private:
    NodePossibleLabelsPtr values_;
};

// A class implementing Metropolis dynamics algorithm minimization algorithm.
template <typename NodeType, typename DataType, typename RealType>
class MD2D: boost::noncopyable
{
public:
    typedef TypePossibleValues<NodeType> NodePossibleLabels;
    typedef boost::scoped_ptr<NodePossibleLabels> NodePossibleLabelsPtr;
    typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<RealType> >
        Generator;

    MD2D(NodePossibleLabels* possible_values, RealType temp, RealType temp_delta):
        values_(possible_values), t_(temp), t_delta_(temp_delta),
        rng_(boost::mt19937(static_cast<boost::uint32_t>(time(NULL))),
            boost::uniform_real<RealType>(0, 1))
    {    }

    void next_iteration(MRF2D<NodeType, DataType, RealType>& mrf)
    {

    }

private:
    Generator rng_;
    RealType t_;
    RealType t_delta_;
    NodePossibleLabelsPtr values_;
};

} // namespace bo

#endif // MRF_OPTIMIZATION_HPP_40C9F0DC_7E18_4316_A594_63DAFD793CCA_
