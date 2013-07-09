
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

#include "bo/config.hpp"
#include "bo/core/vector.hpp"
#include "bo/math/pca.hpp"
#include "bo/math/mean.hpp"

namespace bo {
namespace surfaces {
namespace detail {

//
template <typename RealType>
class BasePropagation
{
public:
    typedef BasePropagation<RealType> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;

    virtual Point3D get(const Point3D& current) const = 0;
};

//
template <typename RealType>
class InertialPropagation
{
public:
    typedef InertialPropagation<RealType> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;

    InertialPropagation(RealType weight): weight_(weight)
    { }

    virtual Point3D get(const Point3D& current, const Point3D& previous) const
    {
        Point3D inertial = current - previous;

        // TODO: remove this by redesigning the algorithm and requiring inertial
        // vector to be non-zero.
        Point3D inertial_normalized(0);
        try
        {
            // TODO: remove this block by refactoring bo::Vector class.
            // Normalize vector.
            inertial_normalized.assign(inertial.normalized(), 3);
        }
        catch (...)
        { }

        return weight_ * inertial_normalized;
    }

protected:
    RealType weight_;
};

//
template <typename RealType>
class EmptyCentrifugalPropagation: public BasePropagation<RealType>
{
public:
    typedef Vector<RealType, 3> Point3D;

    EmptyCentrifugalPropagation(): zero_vector_(0)
    { }

    virtual Point3D get(const Point3D& current) const
    {
        BO_UNUSED(current);
        return zero_vector_;
    }

protected:
    Point3D zero_vector_;
};

//
template <typename RealType>
class CentrifugalPropagation: public BasePropagation<RealType>
{
public:
    typedef Vector<RealType, 3> Point3D;

    CentrifugalPropagation(const Point3D& center_of_mass, RealType weight):
        center_of_mass_(center_of_mass), weight_(weight)
    { }

    virtual Point3D get(const Point3D& current) const
    {
        // Save some processor tacts and reuse the input variable.
        Point3D inertial = current - center_of_mass_;

        // TODO: remove this by redesigning the algorithm and requiring inertial
        // vector to be non-zero.
        Point3D centrifugal_normalized(0);
        try
        {
            // TODO: remove this block by refactoring bo::Vector class.
            // Normalize vector.
            centrifugal_normalized.assign(inertial.normalized(), 3);
        }
        catch (...)
        { }

        return weight_ * centrifugal_normalized;
    }

protected:
    Point3D center_of_mass_;
    RealType weight_;
};

//
template <typename RealType, typename Tree>
class TangentialPropagation
{
public:
    typedef TangentialPropagation<RealType, Tree> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> Points3D;

    TangentialPropagation(const Tree& tree, RealType radius, RealType weight):
        tree_(tree), radius_(radius), weight_(weight)
    { }

    virtual Point3D get(const Point3D& current, const Point3D& inertial) const
    {
        // Search for nearby points.
        Points3D neighbours;
        tree_.find_within_range(current, radius_, std::back_inserter(neighbours));

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

        return weight_ * tangential_normalized;
    }

protected:
    Tree tree_;
    RealType radius_;
    RealType weight_;
};

//
template <typename RealType, typename Tree>
class NeighbourTangentialPropagation: public TangentialPropagation<RealType, Tree>
{
public:
    typedef TangentialPropagation<RealType, Tree> Base;
    typedef TangentialPropagation<RealType, Tree> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> Points3D;

    typedef std::vector<Tree> Trees;
    typedef std::vector<RealType> Weights;
    typedef std::vector<Base> TangentialPropagations;

    NeighbourTangentialPropagation(const Tree& tree, RealType radius, RealType weight,
            const Trees& neighbour_trees, const Weights& neighbour_weights):
        Base(tree, radius, weight)
    {
        // Make sure neighbour data is consistent.
        BOOST_ASSERT((neighbour_weights.size() == neighbour_trees.size()) &&
                     "Number of provided neighbour planes doesn't correspond to the "
                     "number of weights.");

        //
        std::size_t neighbour_size = neighbour_weights.size();
        total_size_ = neighbour_size + 1;
        tangentials_.reserve(total_size_);

        //
        tangentials_.push_back(Base(tree, radius, RealType(1)));

        for (std::size_t idx = 0; idx < neighbour_size; ++idx)
        {
            tangentials_.push_back(Base(neighbour_trees[idx], radius, neighbour_weights[idx]));
        }
    }

    virtual Point3D get(const Point3D& current, const Point3D& inertial) const
    {
        // Get the tangential propagation for the main plane.
        // Compute weighted tangential propagations for neighbours and add them to the
        // total tangential propagation.
        Point3D total_tangential;
        for (std::size_t idx = 0; idx < total_size_; ++idx)
        {
            total_tangential += tangentials_[idx].get(current, inertial);
        }

        // TODO: remove this block by refactoring bo::Vector class.
        // Normalize vector.
        Point3D total_tangential_normalized(total_tangential.normalized(), 3);

        return total_tangential_normalized;
    }

protected:
    TangentialPropagations tangentials_;
    std::size_t total_size_;
};

//
template <typename RealType, typename Tree>
class PropagationDirection
{
public:
    typedef InertialPropagation<RealType> Inertial;
    typedef typename Inertial::Ptr InertialPtr;
    typedef BasePropagation<RealType> Centrifugal;
    typedef typename Centrifugal::Ptr CentrifugalPtr;
    typedef TangentialPropagation<RealType, Tree> Tangential;
    typedef typename Tangential::Ptr TangentialPtr;
    typedef Vector<RealType, 3> Point3D;

public:
    PropagationDirection(InertialPtr inertial_ptr, CentrifugalPtr centrifugal_ptr,
                         TangentialPtr tangential_ptr):
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

private:
    InertialPtr inertial_ptr_;
    CentrifugalPtr centrifugal_ptr_;
    TangentialPtr tangential_ptr_;
};


template <typename RealType, typename Tree>
PropagationDirection<RealType, Tree> create_propagator(RealType inertial_weight,
        const Tree& tree, RealType tangential_radius, bool use_centrifugal,
        const std::vector<Vector<RealType, 3> >& points, RealType centrifugal_weight)
{
    typedef PropagationDirection<RealType, Tree> Propagation;

    // Create an instance of inertial propagation.
    typename Propagation::InertialPtr inertial = boost::make_shared<
            typename Propagation::Inertial>(inertial_weight);

    // Create an instance of centrifugal propagation
    typename Propagation::CentrifugalPtr centrifugal;
    if (use_centrifugal)
    {
        Vector<RealType, 3> center_of_mass = math::mean(points);
        centrifugal = boost::make_shared<CentrifugalPropagation<RealType> >(center_of_mass,
                centrifugal_weight);
    }
    else
    {
        centrifugal = boost::make_shared<EmptyCentrifugalPropagation<RealType> >();
    }

    // Create an instance of tangential propagation.
    typename Propagation::TangentialPtr tangential = boost::make_shared<
            TangentialPropagation<RealType, Tree> >(tree, tangential_radius,
                                                    RealType(1) - inertial_weight - centrifugal_weight);

    Propagation retvalue(inertial, centrifugal, tangential);
    return retvalue;
}







// Represents a point cloud located inside a thin disk. Though the disk can be
// represented by a median plane, slight deviations from it may exist.
template <typename RealType>
class PointsDisk3D
{
public:
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> PlaneData;

    PointsDisk3D(const PlaneData& data): data_(data)
    {
        // Compute plane origin.
        origin_ = bo::math::mean(data_);

        // Employ PCA to estimate plane normal.
        typedef math::PCA<RealType, 3> PCAEngine;
        PCAEngine pca;
        typename PCAEngine::Result result = pca(data_);
        normal_ = result.template get<1>()[0];
    }

    const PlaneData& data() const
    {
        return data_;
    }

    Point3D origin() const
    {
        return origin_;
    }

    Point3D normal() const
    {
        return normal_;
    }

private:
    PlaneData data_;
    Point3D origin_;
    Point3D normal_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // TYPES_HPP_5156559C_A5CD_11E2_B215_221212C84021
