
/******************************************************************************

  Helper types related to the tangential component of the complex propagation.

  Copyright (c) 2013
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

#ifndef TANGENTIAL_PROPAGATION_HPP_96CE789E_EA28_11E2_A507_AB62D8EF6274
#define TANGENTIAL_PROPAGATION_HPP_96CE789E_EA28_11E2_A507_AB62D8EF6274

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>

#include "bo/core/vector.hpp"
#include "bo/math/pca.hpp"

namespace bo {
namespace surfaces {
namespace detail {

// Base class for any tangential propagation.
template <typename RealType>
class BaseTangentialPropagation
{
public:
    typedef BaseTangentialPropagation<RealType> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;

    BaseTangentialPropagation(RealType weight): weight_(weight)
    { }

    virtual Point3D get(const Point3D& current, const Point3D& inertial) const = 0;

protected:
    RealType weight_;
};

// Standard tangential propagation. Uses the provided k-d tree to fetch neighbours of
// the given point and then takes the greatest vector from the PCA launched on
// neighbours.
template <typename RealType, typename Tree>
class TangentialPropagation: public BaseTangentialPropagation<RealType>
{
public:
    typedef TangentialPropagation<RealType, Tree> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> Points3D;
    typedef boost::shared_ptr<Tree> TreePtr;

    TangentialPropagation(TreePtr tree_ptr, RealType radius, RealType weight):
        BaseTangentialPropagation<RealType>(weight), tree_ptr_(tree_ptr), radius_(radius)
    { }

    virtual Point3D get(const Point3D& current, const Point3D& inertial) const
    {
        // Initialize collection for neigbours. Reserve gives slightly better performance,
        // but uses more memory.
        Points3D neighbours;
        neighbours.reserve(tree_ptr_->size());

        // Search for nearby points.
        tree_ptr_->find_within_range(current, radius_, std::back_inserter(neighbours));

        // Employ PCA to extract the tangential propagation vector from the set of
        // nearby points.
        typedef math::PCA<RealType, 3> PCAEngine;
        PCAEngine pca;
        typename PCAEngine::Result result = pca(neighbours);
        Point3D tangential = result.template get<1>()[2];

        // Ensure that tangential and inertial components are codirectional.
        if (tangential * inertial < 0)
            tangential = - tangential;

        // TODO: remove this block by refactoring bo::Vector class.
        // Normalize vector.
        Point3D tangential_normalized(tangential.normalized(), 3);

        return this->weight_ * tangential_normalized;
    }

protected:
    TreePtr tree_ptr_;
    RealType radius_;
};

// Tangential propagation that uses neighbouring planes. Calculates and return the
// weighted sum of standard tangential propagations for the main plane and provided
// neighbours (which are projected onto the main plane). Note that main plane is assumed
// to have weight 1.
template <typename RealType, typename Tree>
class NeighbourTangentialPropagation: public BaseTangentialPropagation<RealType>
{
public:
    typedef TangentialPropagation<RealType, Tree> Tangential;
    typedef NeighbourTangentialPropagation<RealType, Tree> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;

    typedef boost::shared_ptr<Tree> TreePtr;
    typedef std::vector<TreePtr> TreePtrs;
    typedef std::vector<RealType> Weights;
    typedef std::vector<Tangential> TangentialPropagations;

    NeighbourTangentialPropagation(TreePtr tree_ptr, RealType radius, RealType weight,
            const TreePtrs& neighbour_trees, const Weights& neighbour_weights):
        BaseTangentialPropagation<RealType>(weight)
    {
        // Make sure neighbour data is consistent.
        BOOST_ASSERT((neighbour_weights.size() == neighbour_trees.size()) &&
                     "Number of provided neighbour planes doesn't correspond to the "
                     "number of weights.");

        // Cache container sizes and reserve memory.
        std::size_t neighbour_size = neighbour_weights.size();
        total_size_ = neighbour_size + 1;
        tangentials_.reserve(total_size_);

        // Add standard tangential propagation for the main plane.
        tangentials_.push_back(Tangential(tree_ptr, radius, RealType(1)));

        // Add standard tangential propagations for neighbouring planes.
        for (std::size_t idx = 0; idx < neighbour_size; ++idx)
            tangentials_.push_back(Tangential(neighbour_trees[idx], radius,
                                              neighbour_weights[idx]));
    }

    virtual Point3D get(const Point3D& current, const Point3D& inertial) const
    {
        // Compute weighted tangential propagations for the main and neighbouring planes
        // and sum them up.
        Point3D total_tangential;
        for (std::size_t idx = 0; idx < total_size_; ++idx)
        {
            total_tangential += tangentials_[idx].get(current, inertial);
        }

        // TODO: remove this block by refactoring bo::Vector class.
        // Normalize vector.
        Point3D total_tangential_normalized(total_tangential.normalized(), 3);

        return this->weight_ * total_tangential_normalized;
    }

protected:
    TangentialPropagations tangentials_;
    std::size_t total_size_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // TANGENTIAL_PROPAGATION_HPP_96CE789E_EA28_11E2_A507_AB62D8EF6274

