
/******************************************************************************

  Helper types used in surface reconstruction routines.

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

#ifndef TYPES_HPP_5156559C_A5CD_11E2_B215_221212C84021
#define TYPES_HPP_5156559C_A5CD_11E2_B215_221212C84021

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "bo/core/vector.hpp"
#include "bo/surfaces/detail/inertial_propagation.hpp"
#include "bo/surfaces/detail/centrifugal_propagation.hpp"
#include "bo/surfaces/detail/tangential_propagation.hpp"

namespace bo {
namespace surfaces {
namespace detail {

//
template <typename RealType, typename Tree>
class PropagationDirection
{
public:
    typedef PropagationDirection<RealType, Tree> SelfType;

    typedef InertialPropagation<RealType> Inertial;
    typedef typename Inertial::Ptr InertialPtr;

    typedef BaseCentrifugalPropagation<RealType> BaseCentrifugal;
    typedef EmptyCentrifugalPropagation<RealType> EmptyCentrifugal;
    typedef CentrifugalPropagation<RealType> Centrifugal;
    typedef typename BaseCentrifugal::Ptr BaseCentrifugalPtr;

    typedef BaseTangentialPropagation<RealType> BaseTangential;
    typedef TangentialPropagation<RealType, Tree> Tangential;
    typedef NeighbourTangentialPropagation<RealType, Tree> NeighbourTangential;
    typedef typename BaseTangential::Ptr BaseTangentialPtr;

    typedef Vector<RealType, 3> Point3D;
    typedef boost::shared_ptr<Tree> TreePtr;
    typedef std::vector<TreePtr> TreePtrs;
    typedef std::vector<RealType> Weights;

public:
    PropagationDirection(InertialPtr inertial_ptr, BaseCentrifugalPtr centrifugal_ptr,
                         BaseTangentialPtr tangential_ptr):
        inertial_ptr_(inertial_ptr), centrifugal_ptr_(centrifugal_ptr),
        tangential_ptr_(tangential_ptr)
    { }

    // TODO: remove this.
    PropagationDirection()
    { }

    virtual Point3D initial(const Point3D& start) const
    {
        Point3D tangential = tangential_ptr_->get(start, Point3D(RealType(0)));
        Point3D tangential_normalized(tangential.normalized(), 3);

        return
            tangential_normalized;
    }

    virtual Point3D next(const Point3D& current, const Point3D& previous) const
    {
        Point3D inertial = inertial_ptr_->get(current, previous);
        Point3D centrifugal = centrifugal_ptr_->get(current);
        Point3D tangential = tangential_ptr_->get(current, inertial);
        Point3D total = inertial + centrifugal + tangential;
        Point3D total_normalized(total.normalized(), 3);

        return total_normalized;
    }

    static SelfType create_simple(RealType inertial_weight, TreePtr tree_ptr,
                                  RealType tangential_radius)
    {
        // Create an instance of inertial propagation.
        InertialPtr inertial = boost::make_shared<Inertial>(inertial_weight);

        // Create an instance of centrifugal propagation.
        // TODO: replace float comparison.
        BaseCentrifugalPtr centrifugal =  boost::make_shared<EmptyCentrifugal>();

        // Create an instance of tangential propagation.
        BaseTangentialPtr tangential = boost::make_shared<Tangential>(tree_ptr, tangential_radius,
                RealType(1) - inertial_weight);

        SelfType retvalue(inertial, centrifugal, tangential);
        return retvalue;
    }

    // The sum of inertial and centrifugal weights should lie in [0; 1].
    static SelfType create_with_centrifugal(RealType inertial_weight, RealType centrifugal_weight,
            TreePtr tree_ptr, RealType tangential_radius, const Point3D& center_of_mass)
    {
        // Create an instance of inertial propagation.
        InertialPtr inertial = boost::make_shared<Inertial>(inertial_weight);

        // Create an instance of centrifugal propagation.
        BaseCentrifugalPtr centrifugal = boost::make_shared<Centrifugal>(center_of_mass,
                centrifugal_weight);

        // Create an instance of tangential propagation.
        BaseTangentialPtr tangential = boost::make_shared<Tangential>(tree_ptr, tangential_radius,
                RealType(1) - inertial_weight - centrifugal_weight);

        SelfType retvalue(inertial, centrifugal, tangential);
        return retvalue;
    }

    static SelfType create_with_neighbours(RealType inertial_weight, TreePtr tree_ptr,
            RealType tangential_radius, const TreePtrs& neighbour_trees,
            const Weights& neighbour_weights)
    {
        // Create an instance of inertial propagation.
        InertialPtr inertial = boost::make_shared<Inertial>(inertial_weight);

        // Create an instance of centrifugal propagation.
        // TODO: replace float comparison.
        BaseCentrifugalPtr centrifugal =  boost::make_shared<EmptyCentrifugal>();

        // Create an instance of tangential propagation.
        BaseTangentialPtr tangential = boost::make_shared<NeighbourTangential>(tree_ptr,
                tangential_radius, RealType(1) - inertial_weight, neighbour_trees,
                neighbour_weights);

        SelfType retvalue(inertial, centrifugal, tangential);
        return retvalue;
    }

    // The sum of inertial and centrifugal weights should lie in [0; 1].
    static SelfType create_with_neighbours_and_centrifugal(RealType inertial_weight,
            RealType centrifugal_weight, TreePtr tree_ptr, RealType tangential_radius,
            const Point3D& center_of_mass, const TreePtrs& neighbour_trees,
            const Weights& neighbour_weights)
    {
        // Create an instance of inertial propagation.
        InertialPtr inertial = boost::make_shared<Inertial>(inertial_weight);

        // Create an instance of centrifugal propagation.
        // TODO: replace float comparison.
        BaseCentrifugalPtr centrifugal = boost::make_shared<Centrifugal>(center_of_mass,
                centrifugal_weight);

        // Create an instance of tangential propagation.
        BaseTangentialPtr tangential = boost::make_shared<NeighbourTangential>(tree_ptr,
                tangential_radius, RealType(1) - inertial_weight - centrifugal_weight,
                neighbour_trees, neighbour_weights);

        SelfType retvalue(inertial, centrifugal, tangential);
        return retvalue;
    }

private:
    InertialPtr inertial_ptr_;
    BaseCentrifugalPtr centrifugal_ptr_;
    BaseTangentialPtr tangential_ptr_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // TYPES_HPP_5156559C_A5CD_11E2_B215_221212C84021
