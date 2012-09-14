
/******************************************************************************

  mrf_2d.hpp, v 0.1.0 2012.09.11

  Markov random field model for regular 2D lattice.

  Copyright (c) 2009 - 2012
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

#ifndef MRF_2D_HPP_1140ED81_1E3E_4AEB_AFBB_6CA70FE3EF9B_
#define MRF_2D_HPP_1140ED81_1E3E_4AEB_AFBB_6CA70FE3EF9B_

#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <boost/random.hpp>
#include <boost/assert.hpp>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>

#include "bo/raw_image_2d.hpp"

namespace bo {

// A class representing a Markov random field on a regular 2D lattice. Operates on
// random field configuration using given observation and clique functions. Only first
// order (pairwise) neighbourhood is supproted. No copies of an instance are allowed.
template <typename NodeType, typename DataType, typename RealType>
class MRF2D: public boost::noncopyable
{
public:
    typedef NodeType& reference;
    typedef const NodeType& const_reference;

    typedef RawImage2D<NodeType> RandomLattice;
    typedef RawImage2D<DataType> DataLattice;
    typedef boost::function<RealType (DataType, NodeType)> LikelihoodEnergy;
    typedef boost::function<RealType (NodeType, NodeType)> PriorEnergy;

    MRF2D(const RandomLattice& initial_configuration, const DataLattice& observation,
          LikelihoodEnergy likelihood, PriorEnergy prior);

    RealType compute_full_energy() const;
    RealType compute_local_energy(NodeType val, std::size_t col, std::size_t row) const;

    const_reference operator()(std::size_t col, std::size_t row) const;
    reference operator()(std::size_t col, std::size_t row);

    std::size_t width() const;
    std::size_t height() const;

private:
    RealType right_clique_(NodeType val, std::size_t col, std::size_t row) const;
    RealType down_clique_(NodeType val, std::size_t col, std::size_t row) const;
    RealType left_clique_(NodeType val, std::size_t col, std::size_t row) const;
    RealType up_clique_(NodeType val, std::size_t col, std::size_t row) const;

    RealType prior_fun_(NodeType val, std::size_t neigh_col, std::size_t neigh_row) const;
    RealType likelihood_fun_(NodeType val, std::size_t col, std::size_t row) const;

private:
    RandomLattice configuration_;
    DataLattice observation_;
    LikelihoodEnergy likelihood_;
    PriorEnergy prior_;
};


template <typename NodeType, typename DataType, typename RealType>
MRF2D<NodeType, DataType, RealType>::MRF2D(const RandomLattice& initial_configuration,
    const DataLattice& observation, LikelihoodEnergy likelihood, PriorEnergy prior):
    configuration_(initial_configuration),
    observation_(observation),
    likelihood_(likelihood),
    prior_(prior)
{
    // Check if the dimensions of observation and configuration are the same.
    if ((configuration_.width() != observation_.width()) ||
        (configuration_.height() != configuration_.height()))
    {
        BOOST_ASSERT_MSG(false, "Dimensions of the configuration and observation do not coincide.");
        throw std::logic_error("MRF2D: dimensions of the random field lattice and the "
                               "observation data must be the same.");
    }
}

// Computes full energy of the current MRF state. It takes sum over all zero order
// (single node) and first order (every two adjacent nodes) cliques using corresponding
// clique functions. In order to iterate over all first order cliques, for every
// node we can consider only right and down neighbour except for some border nodes.
template <typename NodeType, typename DataType, typename RealType>
RealType MRF2D<NodeType, DataType, RealType>::compute_full_energy() const
{
    RealType energy(0);

    for (std::size_t col = 0; col < configuration_.width(); ++col) {
        for (std::size_t row = 0; row < configuration_.height(); ++row) {
            NodeType value = configuration_(col, row);
            energy += (likelihood_fun_(value, col, row) +
                       right_clique_(value, col, row) +
                       down_clique_(value, col, row));
    }   }

    return energy;
}

// Computes a local energy of the neighbourhood of the given node with provided value.
// Border checks are done inside corresponding clique functions.
template <typename NodeType, typename DataType, typename RealType>
RealType MRF2D<NodeType, DataType, RealType>::compute_local_energy(NodeType val,
    std::size_t col, std::size_t row) const
{
    // Compute likelihood energy for the given value at the given position.
    RealType energy = likelihood_fun_(val, col, row);

    // Compute prior energy of the given node given its neighbours.
    energy += (right_clique_(val, col, row) + down_clique_(val, col, row) +
               left_clique_(val, col, row) + up_clique_(val, col, row));

    return energy;
}

template <typename NodeType, typename DataType, typename RealType> inline
typename MRF2D<NodeType, DataType, RealType>::const_reference
MRF2D<NodeType, DataType, RealType>::operator()(std::size_t col, std::size_t row) const
{
    return configuration_(col, row);
}

template <typename NodeType, typename DataType, typename RealType> inline
typename MRF2D<NodeType, DataType, RealType>::reference
MRF2D<NodeType, DataType, RealType>::operator()(std::size_t col, std::size_t row)
{
    return configuration_(col, row);
}

template <typename NodeType, typename DataType, typename RealType> inline
std::size_t MRF2D<NodeType, DataType, RealType>::width() const
{
    return configuration_.width();
}

template <typename NodeType, typename DataType, typename RealType> inline
std::size_t MRF2D<NodeType, DataType, RealType>::height() const
{
    return configuration_.height();
}


// Next four functions represent cliques and neighbourhood relations in the model.
// They check for border overrun and call the associated clique function.
template <typename NodeType, typename DataType, typename RealType> inline
RealType MRF2D<NodeType, DataType, RealType>::right_clique_(NodeType val,
    std::size_t col, std::size_t row) const
{
    return
        (col < configuration_.width() - 1) ? prior_fun_(val, col + 1, row) : RealType(0);
}

template <typename NodeType, typename DataType, typename RealType> inline
RealType MRF2D<NodeType, DataType, RealType>::down_clique_(NodeType val,
    std::size_t col, std::size_t row) const
{
    return
        (row < configuration_.height() - 1) ? prior_fun_(val, col, row + 1) : RealType(0);
}

template <typename NodeType, typename DataType, typename RealType> inline
RealType MRF2D<NodeType, DataType, RealType>::left_clique_(NodeType val,
    std::size_t col, std::size_t row) const
{
    return
        (col > 0) ? prior_fun_(val, col - 1, row) : RealType(0);
}

template <typename NodeType, typename DataType, typename RealType> inline
RealType MRF2D<NodeType, DataType, RealType>::up_clique_(NodeType val,
    std::size_t col, std::size_t row) const
{
    return
        (row > 0) ? prior_fun_(val, col, row - 1) : RealType(0);
}

// Computes prior energy for two nodes.
template <typename NodeType, typename DataType, typename RealType> inline
RealType MRF2D<NodeType, DataType, RealType>::prior_fun_(NodeType val,
    std::size_t neigh_col, std::size_t neigh_row) const
{
    return prior_(val, configuration_(neigh_col, neigh_row));
}

// Computes likelihood energy for the given node.
template <typename NodeType, typename DataType, typename RealType> inline
RealType MRF2D<NodeType, DataType, RealType>::likelihood_fun_(NodeType val,
    std::size_t col, std::size_t row) const
{
    return likelihood_(val, observation_(col, row));
}

} // namespace bo

#endif // MRF_2D_HPP_1140ED81_1E3E_4AEB_AFBB_6CA70FE3EF9B_
