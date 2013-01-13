
/******************************************************************************

  parallel_planes_tiling.hpp, v 1.0.19 2013.01.12

  Implementation of several surface tiling methods, working with parallel planes.

  Copyright (c) 2012, 2013
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

#ifndef PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
#define PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_

#include <vector>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "bo/mesh.hpp"
#include "bo/raw_image_2d.hpp"
#include "bo/kdtree.hpp"
#include "bo/blas/blas.hpp"
#include "bo/blas/pca.hpp"
#include "bo/methods/distances_3d.hpp"

#include "bo/extended_std.hpp"

namespace bo {
namespace methods {
namespace surfaces {

namespace detail {

template <typename Container>
class TraverseRule
{
public:
    typedef TraverseRule<Container> SelfType;
    typedef boost::shared_ptr<const Container> ContainerConstPtr;
    typedef typename Container::const_iterator ContainerConstIterator;
    typedef typename std::iterator_traits<ContainerConstIterator>::reference Reference;

    TraverseRule(ContainerConstPtr contour): contour_(contour)
    { }

    virtual void add(std::size_t offset) = 0;
    virtual Reference dereference() const = 0;
    virtual bool check_validity() const = 0;
    virtual SelfType* clone() const = 0;

    virtual ~TraverseRule()
    { }

    // TODO: add more macros to extract common code from virtual functions.

    // Use this macro to define an implementation of the clone() function in descendants.
    #define BO_TRAVERSE_RULE_CLONE_IMPL \
        virtual SelfType* clone() const \
        { return new SelfType(*this); }

protected:
    ContainerConstPtr contour_;
};


template <typename Container>
class FwdOnePassTraverseRule: public TraverseRule<Container>
{
public:
    typedef FwdOnePassTraverseRule<Container> SelfType;
    typedef typename Container::const_iterator FwdIterator;

    FwdOnePassTraverseRule(ContainerConstPtr contour): TraverseRule(contour)
    { fwd_it_ = contour_->begin(); }

    virtual void add(std::size_t offset)
    { fwd_it_ += offset; }

    virtual Reference dereference() const
    { return *fwd_it_; }

    virtual bool check_validity() const
    { return (fwd_it_ != contour_->end()); }

    virtual ~FwdOnePassTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    FwdIterator fwd_it_;
};

template <typename Container>
class BwdOnePassTraverseRule: public TraverseRule<Container>
{
public:
    typedef BwdOnePassTraverseRule<Container> SelfType;
    typedef typename Container::const_reverse_iterator BwdIterator;

    BwdOnePassTraverseRule(ContainerConstPtr contour): TraverseRule(contour)
    { bwd_it_ = contour_->rbegin(); }

    virtual void add(std::size_t offset)
    { bwd_it_ += offset; }

    virtual Reference dereference() const
    { return *bwd_it_; }

    virtual bool check_validity() const
    { return (bwd_it_ != contour_->rend()); }

    virtual ~BwdOnePassTraverseRule()
    { }

    BO_TRAVERSE_RULE_CLONE_IMPL

protected:
    BwdIterator bwd_it_;
};


template <typename Container>
class ContainerConstTraverser
{
public:
    typedef ContainerConstTraverser<Container> SelfType;

    typedef typename Container::const_iterator ContainerConstIterator;
    typedef typename std::iterator_traits<ContainerConstIterator>::reference Reference;

    typedef TraverseRule<Container> TraverseRuleType;
    typedef boost::shared_ptr<TraverseRuleType> TraverseRulePtr;

    ContainerConstTraverser(): valid_(false)
    { }

    ContainerConstTraverser(TraverseRulePtr iterator_impl): rule_(iterator_impl)
    {
        valid_ = rule_->check_validity();
    }

    ContainerConstTraverser(const SelfType& other): valid_(other.valid_),
        rule_(other.rule_->clone())
    { }

    SelfType& operator=(SelfType other)
    {
        other.swap(*this);
        return *this;
    }

    SelfType& operator++()
    {
        return ((*this) += 1);
    }

    SelfType& operator+=(std::size_t offset)
    {
        if (is_valid())
        {
            rule_->add(offset);
            valid_ = rule_->check_validity();
        }

        return (*this);
    }

    SelfType operator+(std::size_t offset) const
    {
        SelfType temp = *this;
        return
            (temp += offset);
    }

    Reference operator*() const
    {
        if (!is_valid())
            throw std::logic_error("Cannot dereference invalid CommonConstIterator.");

        return
            rule_->dereference();
    }

    void swap(SelfType& other)
    {
        std::swap(valid_, other.valid_);
        rule_.swap(other.rule_);
    }

    bool is_valid() const
    {
        return valid_;
    }

private:
    bool valid_;
    TraverseRulePtr rule_;
};

template <typename Container>
struct TraverseRuleFactory
{
    typedef boost::shared_ptr<const Container> ContainerConstPtr;
    typedef TraverseRule<Container> TraverseRuleType;
    typedef boost::shared_ptr<TraverseRuleType> TraverseRulePtr;

    #define BO_RULE_FACTORY_FUNCTION(RuleName)                              \
        static TraverseRulePtr RuleName(ContainerConstPtr container_ptr)    \
        {                                                                   \
            typedef RuleName##TraverseRule<Container> Rule;                 \
            TraverseRulePtr rule_ptr(new Rule(container_ptr));              \
            return rule_ptr;                                                \
        }

    BO_RULE_FACTORY_FUNCTION(FwdOnePass)
    BO_RULE_FACTORY_FUNCTION(BwdOnePass)

    // TODO: introduce enum instead of bool flags.
    static TraverseRulePtr Create(ContainerConstPtr container_ptr, bool is_forward)
    {
        return (is_forward ? FwdOnePass(container_ptr) : BwdOnePass(container_ptr));
    }
};

}

template <typename RealType>
class MinSpanPropagation: public boost::noncopyable
{
public:
    typedef MinSpanPropagation<RealType> this_type;

    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> Points3D;
    typedef std::vector<Point3D> ParallelPlane;
    typedef Mesh<RealType> Mesh;
    typedef RawImage2D<RealType> Image2D;
    typedef boost::shared_ptr<ParallelPlane> ParallelPlanePtr;
    typedef boost::shared_ptr<const ParallelPlane> ParallelPlaneConstPtr;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;

    typedef KDTree<3, Point3D,
        boost::function<RealType (Point3D, std::size_t)> > Tree;

    typedef blas::PCA<RealType, 3> PCAEngine;

public:
    MinSpanPropagation() { }

    static ParallelPlanePtr load_plane(Image2D data)
    {
        ParallelPlanePtr plane(new ParallelPlane);

        for (std::size_t row = 0; row < data.height(); ++row)
            for (std::size_t col = 0; col < data.width(); ++col)
                if (data(col, row) > 0)
                    plane->push_back(Point3D(RealType(col), RealType(row), RealType(0)));

        return plane;
    }

    static ParallelPlanePtr propagate(ParallelPlaneConstPtr plane, RealType ratio)
    {
        // TODO: assert on data size (at least 2 samples?).

        // Set algorithm parameters.
        RealType delta_min = 1;
        RealType delta_max = 7;
        Metric metric(&euclidean_distance<RealType, 3>);

        // Build kd-tree from given points.
        // TODO: provide kd-tree with current metric.
        // TODO: convert 3D points to 2D and pass this data as reference to kdtree c-tor.
        Tree tree(plane->begin(), plane->end(), std::ptr_fun(point3D_accessor_));

        // Choose initial point and initial propagation. It solely consists of the
        // tangential component, since inertial cannot be defined.
        Point3D start = (*plane)[0];
        Point3D initial_prop = this_type::tangential_propagation_(start, tree,
            delta_max, Point3D(RealType(0)));

        // Run propagation. It may bump into a hole or return a circuit.
        PropagationResult attempt1 = propagate_(start, start, initial_prop, tree,
            delta_min, delta_max, 1000, metric, ratio);

        // If the hole was detected, run propagation in a different direction.
        PropagationResult attempt2;
        if (attempt1.hole_encountered)
            attempt2 = propagate_(start, attempt1.points->back(), - initial_prop,
                tree, delta_min, delta_max, 1000, metric, ratio);

        // Glue propagation results together.
        ParallelPlanePtr result(new ParallelPlane);
        for (typename ParallelPlane::const_reverse_iterator rit = attempt2.points->rbegin();
             rit != attempt2.points->rend(); ++rit)
            result->push_back(*rit);
        for (typename ParallelPlane::const_iterator it = attempt1.points->begin() + 1;
             it != attempt1.points->end(); ++it)
            result->push_back(*it);

        return result;
    }

    static Mesh christiansen_triangulation(ParallelPlaneConstPtr contour1,
                                           ParallelPlaneConstPtr contour2)
    {
        // TODO: assert on data size (at least 2 samples?).

        // TODO: Determine whether both contours are closed or not. This determines the
        // algorithm's behaviour.

        bool is_closed = false;

        // 1. Contours are closed.
        // Choose a vertex and a direction on the first contour.
        // Find closest vertex on the second contour and determine direction.
        // Iteratively sample points from contours until the surface patch is ready.

        // 2. Contours having exactly one hole.
        // Take the first vertex (at the hole) and direction on the first contour.
        // Take the first vertex (at the hole) and determine co-directed movement.
        // Iterate till the end of both contours.
        typedef detail::ContainerConstTraverser<ParallelPlane> ContourIterator;
        typedef detail::TraverseRuleFactory<ParallelPlane> Factory;

        ContourIterator current1;
        ContourIterator current2;
        bool is_forward2 = true;

        if (is_closed)
        {
            current1 = ContourIterator(Factory::Create(contour1, true));
            Point3D direction1 = *(current1 + 1) - *current1;

            Tree tree2(contour2->begin(), contour2->end(), std::ptr_fun(point3D_accessor_));

            //        Point3D supposed_direction = (contour2->front() + 1) - contour2->front();
            //        if (supposed_direction * direction1 < 0)
            //            is_forward = false;

        }
        else
        {
            current1 = ContourIterator(Factory::Create(contour1, true));

            // Contour2 traverse direction should be swapped in order to correspond
            // with the contoru1 direction.
            if ((contour2->back() - *current1).euclidean_norm() <
                    (contour2->front() - *current1).euclidean_norm())
                is_forward2 = false;

            current2 = ContourIterator(Factory::Create(contour2, is_forward2));
        }

        Mesh mesh(contour1->size() + contour2->size());
        std::size_t current1_idx = mesh.add_vertex(*current1);
        std::size_t current2_idx = mesh.add_vertex(*current2);
        while (true)
        {
//            std::cout << std::endl << "New iteration" << std::endl;
//            std::cout << "current2: " << *current2;

            ContourIterator candidate1 = current1 + 1;
            ContourIterator candidate2 = current2 + 1;

//            std::cout << "current1: " << *current1;
//            std::cout << "current2: " << *current2;

            // This guarantees that when one contour ends, points will be sampled
            // solely from the other one.
            RealType span1_norm = (candidate1.is_valid()) ?
                        (*(candidate1) - *current2).euclidean_norm() :
                        std::numeric_limits<RealType>::max();
            RealType span2_norm = (candidate2.is_valid()) ?
                        (*(candidate2) - *current1).euclidean_norm() :
                        std::numeric_limits<RealType>::max();

//            std::cout << "span1: " << span1_norm << "; span2: " << span2_norm << std::endl;

            // Means both contours are exhausted.
            if (!candidate1.is_valid() && !candidate2.is_valid())
                break;

            // Add candidate vertex.
            std::size_t candidate_idx;
            if (span1_norm > span2_norm)
            {
//                std::cout << "adding2: " << *candidate2;
                candidate_idx = candidate2.is_valid() ?
                            mesh.add_vertex(*candidate2) : -1;
            }
            else
            {
//                std::cout << "adding 1: " << *candidate1;
                candidate_idx = candidate1.is_valid() ?
                            mesh.add_vertex(*candidate1) : -1;
            }

            // Add edges to candidate vertex.
            mesh.add_face(Mesh::Face(current1_idx, candidate_idx, current2_idx));

            // Update current vertices.
            if (span1_norm > span2_norm)
            {
//                std::cout << "updating 2: " << std::endl;
                current2 = candidate2;
                current2_idx = candidate_idx;
            }
            else
            {
//                std::cout << "updating 1: " << std::endl;
                current1 = candidate1;
                current1_idx = candidate_idx;
            }

//            std::cout << "current2: " << *current2;
        }

        return mesh;
    }

    static Mesh to_mesh(ParallelPlaneConstPtr contour)
    {
        Mesh mesh(contour->size());
        for (typename ParallelPlane::const_iterator it = contour->begin();
             it != contour->end(); ++it)
            mesh.add_vertex(*it);

        return mesh;
    }

protected:
    struct ArchedStrip
    {
        ArchedStrip(Point3D ref_pt): ref_pt_(ref_pt)
        { }

        ArchedStrip(Point3D ref_pt, RealType delta_min, RealType delta_max, Point3D prop,
                    Metric dist_fun):
            ref_pt_(ref_pt), delta_min_(delta_min), delta_max_(delta_max), prop_(prop),
            dist_fun_(dist_fun)
        { }

        bool operator()(const Point3D& pt)
        {
            // We need to check three conditions are met:
            // 1. the point lies outside the delta_min circle;
            // 2. the point lies inside the delta_max circle;
            // 3. the point lies in front of the diving plane normal to
            //    propagation vector.

            RealType dist = dist_fun_(ref_pt_, pt);
            bool cond1 = (dist >= delta_min_);
            bool cond2 = (dist <= delta_max_);

            // Dot product is positive if the angle between vectors < 90 deg.
            bool cond3 = ((prop_ * (pt - ref_pt_)) > 0);

            return (cond1 && cond2 && cond3);
        }

        Point3D ref_pt_;
        RealType delta_min_;
        RealType delta_max_;
        Point3D prop_;
        Metric dist_fun_;
    };

    struct PropagationResult
    {
        PropagationResult(): maxsize_reached(false), hole_encountered(false),
            points(new ParallelPlane)
        { }

        bool maxsize_reached;
        bool hole_encountered;
        ParallelPlanePtr points;
    };

private:
    // Helper function for KDTree instance.
    static inline RealType point3D_accessor_(Point3D pt, std::size_t k)
    {
        return pt[k];
    }

    // Computes the inertial propagation vector from current and previous mesh vertices.
    static Point3D inertial_propagation_(Point3D current, Point3D previous)
    {
        Point3D inertial = (current -= previous);

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

        return inertial_normalized;
    }

    // Computes the tangential propagation vector for the given point.
    static Point3D tangential_propagation_(Point3D pt, const Tree& tree,
                                           RealType radius, Point3D inertial)
    {
        // Search for nearby points.
        Points3D neighbours;
        tree.find_within_range(pt, radius, std::back_inserter(neighbours));

        // Employ PCA to extract the tangential propagation vector from the set of
        // nearby points.
        PCAEngine pca;
        typename PCAEngine::Result result = pca(neighbours);
        Point3D tangential = result.template get<1>()[2];

        // Ensure that tangential and inertial components are codirectional.
        if (tangential * inertial < 0)
            tangential = - tangential;

        // TODO: remove this block by refactoring bo::Vector class.
        // Normalize vector.
        Point3D tangential_normalized(tangential.normalized(), 3);

        return tangential_normalized;
    }

    // Computes the total propagation vector from tangential and inertial components.
    static Point3D total_propagation_(Point3D tangential, Point3D inertial, RealType ratio)
    {
        Point3D total = tangential * ratio + inertial * (RealType(1) - ratio);

        // TODO: remove this block by refactoring bo::Vector class.
        // Normalize vector.
        Point3D total_normalized(total.normalized(), 3);

        return total_normalized;
    }

    // End point is not included, start is always included.
    static PropagationResult propagate_(const Point3D& start, const Point3D& end,
                                       Point3D total_prop, const Tree& tree,
                                       RealType delta_min, RealType delta_max,
                                       std::size_t max_size, Metric metric, RealType ratio)
    {
        PropagationResult retvalue;
        retvalue.points->push_back(start);

        Point3D current = start;
        do
        {
            // Restrict total length (to prevent looping).
            if (retvalue.points->size() > max_size)
            {
                retvalue.maxsize_reached = true;
                break;
            }

            // Search for the propagation candidate. It should lie on the arced strip
            // bounded by delta_min and delta_max circumferences and a plane containing
            // current point and normal to propagation vector.
            Point3D phantom_candidate = current + total_prop * delta_min;

            // TODO: get rid of max radius.
            std::pair<typename Tree::const_iterator, RealType> candidate_data =
                    tree.find_nearest_if(phantom_candidate, delta_max,
                        ArchedStrip(current, delta_min, delta_max, total_prop, metric));
            RealType candidate_distance = candidate_data.second;
            Point3D candidate = *(candidate_data.first);

            // Check if we bump into a hole.
            if (candidate_distance >= delta_max)
            {
                retvalue.hole_encountered = true;
                break;
            }

            // Check if the candidate "sees" the end point "in front".
            if ((delta_max < metric(end, candidate)) || (total_prop * (end - candidate)) < 0)
                retvalue.points->push_back(candidate);
            else
                candidate = end;

            // Update algortihm's state.
            Point3D previous = current;
            current = candidate;

            // Compute total propagation.
            Point3D inertial_prop = this_type::inertial_propagation_(current, previous);
            Point3D tangential_prop = this_type::tangential_propagation_(current, tree,
                delta_max, inertial_prop);
            total_prop = this_type::total_propagation_(tangential_prop, inertial_prop, ratio);

        } while (current != end);

        return retvalue;
    }
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
