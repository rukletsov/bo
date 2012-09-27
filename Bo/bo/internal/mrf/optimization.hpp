
/******************************************************************************

  optimization.hpp, v 0.1.6 2012.09.20

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

#ifndef OPTIMIZATION_HPP_40C9F0DC_7E18_4316_A594_63DAFD793CCA_
#define OPTIMIZATION_HPP_40C9F0DC_7E18_4316_A594_63DAFD793CCA_

#include <cstddef>
#include <cmath>
#include <ctime>
#include <limits>
#include <functional>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/random.hpp>

#include "bo/internal/mrf/mrf_2d.hpp"
#include "bo/internal/mrf/type_values.hpp"

namespace bo {

// A basic class for MRF minimization algorithms. Takes (and captures sole ownership
// over it) an istance of TypePossibleValues<NodeType> to let its descendants sample
// possible values of type NodeType during optimization.
template <typename NodeType, typename DataType, typename RealType>
class MRF2DOptimizer: boost::noncopyable
{
public:
    typedef TypeValues<NodeType> NodePossibleLabels;
    typedef boost::scoped_ptr<NodePossibleLabels> NodePossibleLabelsPtr;

    MRF2DOptimizer(NodePossibleLabels* possible_values): values_(possible_values)
    { }

    virtual void next_iteration(MRF2D<NodeType, DataType, RealType>&) = 0;

    virtual ~MRF2DOptimizer()
    { }

protected:
    NodePossibleLabelsPtr values_;
};

// A class impementing iterated conditional modes minimization algorithm. As its
// parent, captures sole ownership over an istance of TypePossibleValues<NodeType>.
template <typename NodeType, typename DataType, typename RealType>
class ICM2D: public MRF2DOptimizer<NodeType, DataType, RealType>
{
public:
    ICM2D(NodePossibleLabels* possible_values): MRF2DOptimizer(possible_values)
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
};

// A class implementing Metropolis dynamics algorithm minimization algorithm. As its
// parent, captures sole ownership over an istance of TypePossibleValues<NodeType>.
// Additionally takes initial temperature and decreasing multiplier for it, and
// precomputed probability for MMD modification if this modification was requested.
template <typename NodeType, typename DataType, typename RealType>
class MD2D: public MRF2DOptimizer<NodeType, DataType, RealType>
{
public:
    typedef boost::uniform_real<RealType> RealDistribution;
    typedef boost::variate_generator<boost::mt19937, RealDistribution> Generator;

    MD2D(NodePossibleLabels* possible_values, RealType temp, RealType temp_delta,
         bool is_modified, RealType mmd_probability):
        MRF2DOptimizer(possible_values), t_(temp), t_delta_(temp_delta),
        is_modified_(is_modified), mmd_probab_(std::log(mmd_probability)),
        rng_(boost::mt19937(static_cast<boost::uint32_t>(std::time(NULL))), RealDistribution(0, 1))
    { }

    void next_iteration(MRF2D<NodeType, DataType, RealType>& mrf)
    {
        // If a modified version of MD algo (MMD) is used, then the decision probability
        // is fixed and obtained as a parameter. Otherwise, it is generated every time.
        RealType decision_probab = mmd_probab_;

        // Update current temperature. See algorithm description for details.
        t_ *= t_delta_;

        // Do an inner iteration pixelwise.
        for (std::size_t col = 0; col < mrf.width(); ++col) {
            for (std::size_t row = 0; row < mrf.height(); ++row) {
                // Randomly generate a new state for the current pixel.
                NodeType old_value = mrf(col, row);
                NodeType new_value = values_->random();

                // if classical version is used, generate the decision probability.
                if (!is_modified_)
                    decision_probab = std::log(rng_());

                // Compute local energy and accept it or reject.
                if (decision_probab <=
                    (mrf.compute_local_energy(old_value, col, row) -
                     mrf.compute_local_energy(new_value, col, row)) / t_)
                {   // Accept new state.
                    mrf(col, row) = new_value;
                }
        }   }
    }

private:
    Generator rng_;
    RealType t_;
    RealType t_delta_;
    bool is_modified_;
    RealType mmd_probab_;
};

} // namespace bo

#endif // OPTIMIZATION_HPP_40C9F0DC_7E18_4316_A594_63DAFD793CCA_
