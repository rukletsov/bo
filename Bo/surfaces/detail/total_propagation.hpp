
/******************************************************************************

  Helper type for choosing propagation direction necessary in the complex
  propagation.

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

#ifndef TOTAL_PROPAGATION_HPP_24BBFC74_EADA_11E2_9943_2E55D8EF6274
#define TOTAL_PROPAGATION_HPP_24BBFC74_EADA_11E2_9943_2E55D8EF6274

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

// Total propagation comprising three components: inertial, centrifugal and tangential.
// Contains factory functions that allow for appropriate choice of every component
// depending on the passes parameters.
template <typename RealType, typename Tree>
class TotalPropagation
{
public:
    typedef TotalPropagation<RealType, Tree> SelfType;

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

    // This macro checks and computes the weight of the tangential components based
    // on the weights of all other components. Could be a function, but let it be a
    // a macro for consistency.
    #define BO_COMPUTE_TANGENTIAL_WEIGHT_FROM_OTHER_WEIGHTS(SumOfOtherWeights)  \
        BOOST_ASSERT((SumOfOtherWeights >= 0) && (SumOfOtherWeights <= 1) &&    \
                 "Weight of a propagation component must be in [0; 1].");       \
        RealType tangential_weight = RealType(1) - SumOfOtherWeights

    // This macro defines the call to create a boost::shared_ptr instance named
    // InstanceName of the propagation component Type with additional list of arguments.
    #define BO_CREATE_PROPAGATION_INSTANCE(Type, InstanceName, ...)             \
        boost::shared_ptr<Type> InstanceName = boost::make_shared<Type>(__VA_ARGS__)

    // Same as BO_CREATE_PROPAGATION_INSTANCE_BASE but for propagation types that have
    // empty c-tors.
    #define BO_CREATE_PROPAGATION_INSTANCE_NOARGS(Type, InstanceName)           \
        boost::shared_ptr<Type> InstanceName = boost::make_shared<Type>()

    // Use this macro to create and return from the caller an instance of the class
    // with components named conventionally.
    #define BO_CREATE_AND_RETURN_TOTAL_PROPAGATION_INSTANCE                     \
        SelfType retvalue(inertial, centrifugal, tangential);                   \
        return retvalue

public:
    TotalPropagation(InertialPtr inertial_ptr, BaseCentrifugalPtr centrifugal_ptr,
                     BaseTangentialPtr tangential_ptr):
        inertial_ptr_(inertial_ptr), centrifugal_ptr_(centrifugal_ptr),
        tangential_ptr_(tangential_ptr)
    { }

    TotalPropagation()
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


    // Factory functions. Every factory consists of several main steps: computation of
    // the tangential weight, instantiation of all propagation components and creation
    // an istance of the class based on previously created components.

    // Create a total propagation with no centrifugal component and no usage of
    // neighbour data in propagation.
    static SelfType create_simple(RealType inertial_weight, TreePtr tree_ptr,
                                  RealType tangential_radius)
    {
        BO_COMPUTE_TANGENTIAL_WEIGHT_FROM_OTHER_WEIGHTS(inertial_weight);

        BO_CREATE_PROPAGATION_INSTANCE(Inertial, inertial, inertial_weight);
        BO_CREATE_PROPAGATION_INSTANCE_NOARGS(EmptyCentrifugal, centrifugal);
        BO_CREATE_PROPAGATION_INSTANCE(Tangential, tangential, tree_ptr, tangential_radius,
                tangential_weight);

        BO_CREATE_AND_RETURN_TOTAL_PROPAGATION_INSTANCE;
    }

    // Create a total propagation with centrifugal component but no usage of neighbour
    // data in propagation. Note that the sum of inertial and centrifugal weights
    // should be in [0; 1].
    static SelfType create_with_centrifugal(RealType inertial_weight, RealType centrifugal_weight,
            TreePtr tree_ptr, RealType tangential_radius, const Point3D& center_of_mass)
    {
        BO_COMPUTE_TANGENTIAL_WEIGHT_FROM_OTHER_WEIGHTS(inertial_weight + centrifugal_weight);

        BO_CREATE_PROPAGATION_INSTANCE(Inertial, inertial, inertial_weight);
        BO_CREATE_PROPAGATION_INSTANCE(Centrifugal, centrifugal, center_of_mass,
                centrifugal_weight);
        BO_CREATE_PROPAGATION_INSTANCE(Tangential, tangential, tree_ptr, tangential_radius,
                tangential_weight);

        BO_CREATE_AND_RETURN_TOTAL_PROPAGATION_INSTANCE;
    }

    // Create a total propagation with no centrifugal component. Tangential component
    // uses weighted neighbour data.
    static SelfType create_with_neighbours(RealType inertial_weight, TreePtr tree_ptr,
            RealType tangential_radius, const TreePtrs& neighbour_trees,
            const Weights& neighbour_weights)
    {
        BO_COMPUTE_TANGENTIAL_WEIGHT_FROM_OTHER_WEIGHTS(inertial_weight);

        BO_CREATE_PROPAGATION_INSTANCE(Inertial, inertial, inertial_weight);
        BO_CREATE_PROPAGATION_INSTANCE_NOARGS(EmptyCentrifugal, centrifugal);
        BO_CREATE_PROPAGATION_INSTANCE(NeighbourTangential, tangential, tree_ptr,
                tangential_radius, tangential_weight, neighbour_trees, neighbour_weights);

        BO_CREATE_AND_RETURN_TOTAL_PROPAGATION_INSTANCE;
    }

    // Create a total propagation with centrifugal component. Tangential component
    // uses weighted neighbour data. Note that the sum of inertial and centrifugal
    // weights should be in [0; 1].
    static SelfType create_with_neighbours_and_centrifugal(RealType inertial_weight,
            RealType centrifugal_weight, TreePtr tree_ptr, RealType tangential_radius,
            const Point3D& center_of_mass, const TreePtrs& neighbour_trees,
            const Weights& neighbour_weights)
    {
        BO_COMPUTE_TANGENTIAL_WEIGHT_FROM_OTHER_WEIGHTS(inertial_weight + centrifugal_weight);

        BO_CREATE_PROPAGATION_INSTANCE(Inertial, inertial, inertial_weight);
        BO_CREATE_PROPAGATION_INSTANCE(Centrifugal, centrifugal, center_of_mass,
                centrifugal_weight);
        BO_CREATE_PROPAGATION_INSTANCE(NeighbourTangential, tangential, tree_ptr,
                tangential_radius, tangential_weight, neighbour_trees, neighbour_weights);

        BO_CREATE_AND_RETURN_TOTAL_PROPAGATION_INSTANCE;
    }

private:
    InertialPtr inertial_ptr_;
    BaseCentrifugalPtr centrifugal_ptr_;
    BaseTangentialPtr tangential_ptr_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // TOTAL_PROPAGATION_HPP_24BBFC74_EADA_11E2_9943_2E55D8EF6274
