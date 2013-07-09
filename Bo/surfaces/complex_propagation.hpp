
/******************************************************************************

  Implementation of the complex propagation technique used in surface
  reconstruction.

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

#ifndef COMPLEX_PROPAGATION_HPP_D90ED351_6A45_4523_85F3_DA99F52B87C2
#define COMPLEX_PROPAGATION_HPP_D90ED351_6A45_4523_85F3_DA99F52B87C2

#include <vector>
#include <iterator>
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>

#include "bo/core/raw_image_2d.hpp"
#include "bo/core/mesh.hpp"
#include "bo/core/kdtree.hpp"
#include "bo/math/blas_extensions.hpp"
#include "bo/math/pca.hpp"
#include "bo/distances/distances_3d.hpp"
#include "bo/surfaces/detail/types.hpp"
#include "bo/surfaces/detail/arched_strip.hpp"

namespace bo {
namespace surfaces {

template <typename RealType>
class ComplexPropagation: public boost::noncopyable
{
public:
    typedef ComplexPropagation<RealType> this_type;
    typedef boost::shared_ptr<this_type> Ptr;
    typedef std::vector<Ptr> Ptrs;

    typedef Vector<RealType, 3> Point3D;

    // Types for input data.
    typedef std::vector<Point3D> Points3D;
    typedef boost::shared_ptr<Points3D> Points3DPtr;
    typedef boost::shared_ptr<const Points3D> Points3DConstPtr;
    typedef std::vector<Points3DConstPtr> Points3DConstPtrs;

    typedef std::vector<RealType> Weights;

    // Internal types, representing plane in the more convenient way.
    typedef detail::PointsDisk3D<RealType> Plane;
    typedef boost::shared_ptr<Plane> PlanePtr;
    typedef boost::shared_ptr<const Plane> PlaneConstPtr;

    // Types for storing propagation result.
    typedef std::vector<Point3D> PropContour;
    typedef boost::shared_ptr<PropContour> PropContourPtr;

    // Used by factory functions.
    typedef RawImage2D<RealType> Image2D;
    typedef bo::Mesh<RealType> Mesh;

private:
    typedef boost::function<RealType (Point3D, Point3D)> Metric;
    typedef KDTree<3, Point3D, boost::function<RealType (Point3D, std::size_t)> > Tree;
    typedef std::vector<Tree> Trees;
    typedef detail::PropagationDirection<RealType, Tree> Propagator;

public:
    // Factory functions. Create an instance of the class from either of the supported
    // input with the provided parameters.
    static Ptr create(Points3DConstPtr plane, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius);

    static Ptrs create(const Points3DConstPtrs& planes, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius);

    static Ptr from_mesh(const Mesh& mesh, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight,RealType tangential_radius);

    static Ptr from_raw_image(Image2D data, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight,RealType tangential_radius);

    // A special factory. Creates a set of instances one per given plane. Additionally
    // passes a set of neighbours with corresponding weights to each instance.
    static Ptrs create(const Points3DConstPtrs& planes, const Weights& neighbour_weights,
        RealType delta_min, RealType delta_max, RealType inertial_weight,
        RealType centrifugal_weight, RealType tangential_radius);

    // Adds neighbour planes with corresponding weights. These planes are used in
    // the calculation of the tangential component of the propagation, making it
    // dependent of the neighbour planes. This should lead to the desirable smoothing.
    void add_neighbour_planes(const Points3DConstPtrs& neighbour_planes,
                              const Weights& neighbour_weights);

    // Performs propagation in either one or two steps depending on whether the
    // contour has a hole.
    void propagate();

    // Accessor functions to the final contour and algorithm markers. Calling these
    // functions makes sense only after calling propagate(). If has_stopped() marker
    // is set to false, the algortihm is most probably stuck (e.g. jumping between two
    // adjacent points or traversing the cosed contour multiple times).
    bool has_stopped() const;
    bool has_hole() const;
    bool is_closed() const;
    PropContourPtr contour() const;

private:
    struct PropagationResult
    {
        PropagationResult(): stopped(false), has_hole(false),
            points(boost::make_shared<PropContour>())
        { }

        PropagationResult(bool maxsize_reached, bool hole_encountered, PropContourPtr pts):
            stopped(maxsize_reached), has_hole(hole_encountered), points(pts)
        { }

        bool stopped;
        bool has_hole;
        PropContourPtr points;
    };

protected:
    ComplexPropagation(Points3DConstPtr plane, RealType delta_min, RealType delta_max,
                       /*RealType inertial_weight, RealType centrifugal_weight,
                       RealType tangential_radius, */const Tree& tree, Propagator propagator);

    // Helper function for KDTree instance.
    static RealType point3D_accessor_(Point3D pt, std::size_t k);

//    // Computes the inertial propagation vector from current and previous mesh vertices.
//    static Point3D inertial_propagation_norm_(Point3D current, Point3D previous);

//    // Computes the centrifugal propagation vector pointing from the centre of mass to
//    // the current point.
//    static Point3D centrifugal_propagation_norm(Point3D current, Point3D center_of_mass);

//    // Computes the tangential propagation vector for the given point and kd-tree.
//    static Point3D tangential_propagation_(Point3D pt, const Tree& tree,
//                                           RealType radius, Point3D inertial);

//    // Computes the total tangential propagation vector from weighted tangential
//    // propagations of the current plane and its neighbours.
//    Point3D total_tangential_propagation_norm_(Point3D pt, Point3D inertial) const;

//    // Computes the total propagation vector from tangential and inertial components.
//    // Note that the sum of inertial and centrifugal weights should be not greater than 1.
//    static Point3D total_propagation_(Point3D inertial, Point3D centrifugal, Point3D tangential,
//                                      RealType inertial_weight, RealType centrifugal_weight);

    // Tries performing propagation from the start point to the end point.
    // Note that end point is not included in result, while start is always included.
    PropagationResult propagate_(const Point3D& start, const Point3D& end,
                                       Point3D total_prop, std::size_t max_size);

private:
    // Parameters and options of the propagation algorithm.
    RealType delta_min_;
    RealType delta_max_;
//    RealType inertial_weight_;      // The sum of inertial and centrifugal weights
//    RealType centrifugal_weight_;   // should lie in [0; 1].
//    RealType tangential_radius_;
    Metric metric_;

    // Propagation data: main plane and neighbours with corresponding weights.
    // Note that main plane has weight 1.
    PlaneConstPtr main_plane_;
    Tree main_tree_;
    Trees neighbour_trees_;
    Weights neighbour_weights_;

    // Markers and final contour.
    bool stopped_;
    bool has_hole_;
    PropContourPtr contour_;

    Propagator propagator_;
};


// C-tor.
template <typename RealType>
ComplexPropagation<RealType>::ComplexPropagation(Points3DConstPtr plane,
                                                 RealType delta_min,
                                                 RealType delta_max,
//                                                 RealType inertial_weight,
//                                                 RealType centrifugal_weight,
//                                                 RealType tangential_radius,
                                                 const Tree& tree,
                                                 Propagator propagator):
    delta_min_(delta_min), delta_max_(delta_max),
//    inertial_weight_(inertial_weight), centrifugal_weight_(centrifugal_weight),
//    tangential_radius_(tangential_radius),
    metric_(&distances::euclidean_distance<RealType, 3>),
    main_plane_(boost::make_shared<Plane>(*plane)),
    stopped_(false), has_hole_(false),
    main_tree_(tree),
    contour_(boost::make_shared<PropContour>()),
    propagator_(propagator)
{
    // If points number is less than 2, propagation cannot be initialized. And
    // because there is no sense in doing propagation for 0 or 1 points, we throw
    // an error instead of returning the plane unchanged.
    if (plane->size() < 2)
        throw std::logic_error("Cannot run propagation for planes consisting of "
                               "less than 2 vertices.");
}


// Factories.
template <typename RealType>
typename ComplexPropagation<RealType>::Ptr ComplexPropagation<RealType>::create(
        Points3DConstPtr plane, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius)
{
    // Build kd-tree from the given points.
    // TODO: provide kd-tree with current metric?
    Tree tree(plane->begin(), plane->end(), std::ptr_fun(point3D_accessor_));

    Propagator propagator;
    // TODO: replace float comparison.
    if (centrifugal_weight == RealType(0))
    {
        propagator = Propagator::create_simple(inertial_weight, tree, tangential_radius);
    }
    else
    {
        Point3D center_of_mass = bo::math::mean(*plane);
        propagator = Propagator::create_with_centrifugal(inertial_weight, centrifugal_weight,
                tree, tangential_radius, center_of_mass);
    }

    // C-tor is declared private, using boost::make_shared gets complicated.
    Ptr ptr(new this_type(plane, delta_min, delta_max, tree, propagator));
    return ptr;
}

template <typename RealType>
typename ComplexPropagation<RealType>::Ptrs ComplexPropagation<RealType>::create(
        const Points3DConstPtrs& planes, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius)
{
    // Cache values and prepare output container.
    std::size_t planes_count = planes.size();
    Ptrs propagators;
    propagators.reserve(planes_count);

    // Create an instance for each plane
    for (std::size_t idx = 0; idx < planes_count; ++idx)
    {
        // Create an instance, register neighbours and return.
        Ptr propagator = create(planes[idx], delta_min, delta_max, inertial_weight,
                                centrifugal_weight, tangential_radius);
        propagators.push_back(propagator);
    }

    return propagators;
}

template <typename RealType>
typename ComplexPropagation<RealType>::Ptr ComplexPropagation<RealType>::from_mesh(
        const Mesh& mesh, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius)
{
    Points3DConstPtr plane = boost::make_shared<Points3D>(mesh.get_all_vertices());
    return
        create(plane, delta_min, delta_max, inertial_weight, centrifugal_weight, tangential_radius);
}

template <typename RealType>
typename ComplexPropagation<RealType>::Ptr ComplexPropagation<RealType>::from_raw_image(
        Image2D data, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius)
{
    // Load plane data from an istance of RawImage2D.
    Points3DPtr plane = boost::make_shared<Points3D>();

    for (std::size_t row = 0; row < data.height(); ++row)
        for (std::size_t col = 0; col < data.width(); ++col)
            if (data(col, row) > 0)
                plane->push_back(Point3D(RealType(col), RealType(row), RealType(0)));

    return
        create(plane, delta_min, delta_max, inertial_weight, centrifugal_weight, tangential_radius);
}

template <typename RealType>
typename ComplexPropagation<RealType>::Ptrs ComplexPropagation<RealType>::create(
        const Points3DConstPtrs &planes, const Weights& neighbour_weights,
        RealType delta_min, RealType delta_max, RealType inertial_weight,
        RealType centrifugal_weight, RealType tangential_radius)
{
    // Check neighbours count (in one direction) is not greater than planes count - 1.
    if (neighbour_weights.size() > planes.size() - 1)
        throw std::logic_error("Neighbours count in one direction exceeds reasonable "
                               "limit.");

    // Cache values and prepare output container.
    std::size_t planes_count = planes.size();
    std::size_t radius = neighbour_weights.size();
    Ptrs propagators;
    propagators.reserve(planes_count);

    // Compute neighbours for each plane and create an instance of the class for it.
    for (std::size_t idx = 0; idx < planes_count; ++idx)
    {
        Points3DConstPtrs neighbours;
        Weights weights;

        // Try add neighbours before the current idx. We add radius items at most,
        // but less if no neighbours are available (close to container front).
        std::size_t neighbour_idx = idx - 1;
        std::size_t added_count = 0;
        std::size_t maxval = std::size_t(-1); // dirty hack to workaround int behaviour for size_t.
        while ((neighbour_idx < maxval) && (added_count < radius))
        {
            neighbours.push_back(planes[neighbour_idx]);
            weights.push_back(neighbour_weights[added_count]);
            ++added_count;
            --neighbour_idx;
        }

        // Try add neighbours after the current idx. We add radius items at most,
        // but less if no neighbours are available (close to container end).
        neighbour_idx = idx + 1;
        added_count = 0;
        while ((neighbour_idx < planes_count) && (added_count < radius))
        {
            neighbours.push_back(planes[neighbour_idx]);
            weights.push_back(neighbour_weights[added_count]);
            ++added_count;
            ++neighbour_idx;
        }

        // Create an instance, register neighbours and return.
        Ptr propagator = create(planes[idx], delta_min, delta_max, inertial_weight,
                                centrifugal_weight, tangential_radius);
        propagator->add_neighbour_planes(neighbours, weights);

        propagators.push_back(propagator);
    }

    return propagators;
}


// Public control functions.
template <typename RealType>
void ComplexPropagation<RealType>::add_neighbour_planes(
        const Points3DConstPtrs &neighbour_planes, const Weights& neighbour_weights)
{
    // Weights should correspond to plane number.
    if (neighbour_planes.size() != neighbour_weights.size())
        throw std::logic_error("Provided weights quantity doesn't correspond to "
                               "planes number");

    // Cache main plane origin and normal.
    Point3D origin = main_plane_->origin();
    Point3D norm = main_plane_->normal();

    // Project all neighbours to the current plane.
    // TODO: rewrite using transform?
    Points3DConstPtrs neighbour_projs;
    neighbour_projs.reserve(neighbour_planes.size());
    for (typename Points3DConstPtrs::const_iterator plane_it =
         neighbour_planes.begin(); plane_it != neighbour_planes.end(); ++plane_it)
    {
        Points3DPtr plane_proj = boost::make_shared<Points3D>();
        plane_proj->reserve((*plane_it)->size());

        for (typename Points3D::const_iterator point_it = (*plane_it)->begin();
             point_it != (*plane_it)->end(); ++point_it)
        {
            plane_proj->push_back(distances::project_point_onto_plane(*point_it, origin, norm));
        }

        neighbour_projs.push_back(plane_proj);
    }

    // Compute kd-trees for projected neighbour planes.
    neighbour_trees_.reserve(neighbour_projs.size());
    for (typename Points3DConstPtrs::const_iterator plane_it =
         neighbour_projs.begin(); plane_it != neighbour_projs.end(); ++plane_it)
    {
        Tree neighbour_tree((*plane_it)->begin(), (*plane_it)->end(),
                            std::ptr_fun(point3D_accessor_));
        neighbour_trees_.push_back(neighbour_tree);
    }

    // Store neighbour weights.
    neighbour_weights_ = neighbour_weights;
}

template <typename RealType>
void ComplexPropagation<RealType>::propagate()
{
    // Choose initial point and initial propagation. It solely consists of the
    // tangential component, since inertial cannot be defined.
    Point3D start = main_plane_->data()[0];
    Point3D initial_prop =
            propagator_.initial(start);
//            this_type::tangential_propagation_(start, main_tree_,
//        tangential_radius_, Point3D(RealType(0)));

    // Run propagation. It may bump into a hole or return a circuit.
    PropagationResult attempt1 = propagate_(start, start, initial_prop, 1000);

    // If the hole was detected, run propagation in a different direction.
    PropagationResult attempt2;
    if (attempt1.has_hole)
        attempt2 = propagate_(start, attempt1.points->back(), - initial_prop, 1000);

    // Glue propagation results together.
    for (typename PropContour::const_reverse_iterator rit = attempt2.points->rbegin();
         rit != attempt2.points->rend(); ++rit)
        contour_->push_back(*rit);
    for (typename PropContour::const_iterator it = attempt1.points->begin() + 1;
         it != attempt1.points->end(); ++it)
        contour_->push_back(*it);

    // Set markers.
    stopped_ = attempt1.stopped || attempt2.stopped;
    has_hole_ = attempt1.has_hole;
}

template <typename RealType> inline
bool ComplexPropagation<RealType>::has_stopped() const
{
    return stopped_;
}

template <typename RealType> inline
bool ComplexPropagation<RealType>::has_hole() const
{
    return has_hole_;
}

template <typename RealType> inline
bool ComplexPropagation<RealType>::is_closed() const
{
    return !(this->has_hole());
}

template <typename RealType> inline
typename ComplexPropagation<RealType>::PropContourPtr
ComplexPropagation<RealType>::contour() const
{
    return contour_;
}


// Private static helper functions.
template <typename RealType> inline
RealType ComplexPropagation<RealType>::point3D_accessor_(Point3D pt, std::size_t k)
{
    return pt[k];
}

//template <typename RealType>
//typename ComplexPropagation<RealType>::Point3D
//ComplexPropagation<RealType>::inertial_propagation_norm_(Point3D current, Point3D previous)
//{
//    Point3D inertial = (current -= previous);

//    // TODO: remove this by redesigning the algorithm and requiring inertial
//    // vector to be non-zero.
//    Point3D inertial_normalized(0);
//    try
//    {
//        // TODO: remove this block by refactoring bo::Vector class.
//        // Normalize vector.
//        inertial_normalized.assign(inertial.normalized(), 3);
//    }
//    catch (...)
//    { }

//    return inertial_normalized;
//}

//template <typename RealType>
//typename ComplexPropagation<RealType>::Point3D
//ComplexPropagation<RealType>::centrifugal_propagation_norm(Point3D current,
//                                                           Point3D center_of_mass)
//{
//    Point3D centrifugal = (current -= center_of_mass);

//    // TODO: remove this by redesigning the algorithm and requiring inertial
//    // vector to be non-zero.
//    Point3D centrifugal_normalized(0);
//    try
//    {
//        // TODO: remove this block by refactoring bo::Vector class.
//        // Normalize vector.
//        centrifugal_normalized.assign(centrifugal.normalized(), 3);
//    }
//    catch (...)
//    { }

//    return centrifugal_normalized;
//}

//template <typename RealType>
//typename ComplexPropagation<RealType>::Point3D
//ComplexPropagation<RealType>::tangential_propagation_(Point3D pt, const Tree& tree,
//                                                      RealType radius, Point3D inertial)
//{
//    // Search for nearby points.
//    Points3D neighbours;
//    tree.find_within_range(pt, radius, std::back_inserter(neighbours));

//    // Employ PCA to extract the tangential propagation vector from the set of
//    // nearby points.
//    typedef math::PCA<RealType, 3> PCAEngine;
//    PCAEngine pca;
//    typename PCAEngine::Result result = pca(neighbours);
//    Point3D tangential = result.template get<1>()[2];

//    // Ensure that tangential and inertial components are codirectional.
//    if (tangential * inertial < 0)
//        tangential = - tangential;

//    return tangential;
//}

//template <typename RealType>
//typename ComplexPropagation<RealType>::Point3D
//ComplexPropagation<RealType>::total_tangential_propagation_norm_(Point3D pt,
//                                                                 Point3D inertial) const
//{
//    // Get the tangential propagation for the main plane.
//    Point3D total_tangential = this_type::tangential_propagation_(pt, main_tree_,
//        tangential_radius_, inertial);

//    // Make sure neighbour data is consistent.
//    BOOST_ASSERT((neighbour_weights_.size() == neighbour_trees_.size()) &&
//                 "Number of provided neighbour planes doesn't correspond to the "
//                 "number of weights.");

//    // Compute weighted tangential propagations for neighbours and add them  to the
//    // total tangential propagation.
//    for (std::size_t idx = 0; idx < neighbour_trees_.size(); ++idx)
//    {
//        Point3D cur_tang = this_type::tangential_propagation_(pt,
//            neighbour_trees_[idx], tangential_radius_, inertial);
//        RealType cur_weight = neighbour_weights_[idx];
//        total_tangential += (cur_tang * cur_weight);
//    }

//    // TODO: remove this block by refactoring bo::Vector class.
//    // Normalize vector.
//    Point3D total_tangential_normalized(total_tangential.normalized(), 3);

//    return total_tangential_normalized;
//}

//template <typename RealType>
//typename ComplexPropagation<RealType>::Point3D
//ComplexPropagation<RealType>::total_propagation_(Point3D inertial, Point3D centrifugal,
//                                                 Point3D tangential, RealType inertial_weight,
//                                                 RealType centrifugal_weight)
//{
//    Point3D total = inertial * inertial_weight + centrifugal * centrifugal_weight +
//            tangential * (RealType(1) - inertial_weight - centrifugal_weight);

//    // TODO: remove this block by refactoring bo::Vector class.
//    // Normalize vector.
//    Point3D total_normalized(total.normalized(), 3);

//    return total_normalized;
//}

template <typename RealType>
typename ComplexPropagation<RealType>::PropagationResult
ComplexPropagation<RealType>::propagate_(const Point3D& start, const Point3D& end,
                                         Point3D total_prop, std::size_t max_size)
{
    typedef detail::ArchedStrip<RealType, 3> ArchedStrip;

    PropagationResult retvalue;
    retvalue.points->push_back(start);

    // Flag for so-called "smooth-ending". It is necessary to keep the distance
    // between the points in the end phase as close to delta_min as possible,
    // despite a delta_max point has been already found.
    bool end_detected = false;

    Point3D current = start;
    do
    {
        // Restrict total length (to prevent looping).
        if (retvalue.points->size() > max_size)
        {
            retvalue.stopped = true;
            break;
        }

        // Search for the propagation candidate. It should lie on the arced strip
        // bounded by delta_min and delta_max circumferences and a plane containing
        // current point and normal to propagation vector.
        Point3D phantom_candidate = current + total_prop * delta_min_;

        // TODO: get rid of max radius.
        std::pair<typename Tree::const_iterator, RealType> candidate_data =
                main_tree_.find_nearest_if(phantom_candidate, delta_max_,
                    ArchedStrip(current, delta_min_, delta_max_, total_prop, metric_));
        RealType candidate_distance = candidate_data.second;
        Point3D candidate = *(candidate_data.first);

        // Check if we bump into a hole.
        if (candidate_distance >= delta_max_)
        {
            retvalue.has_hole = true;
            break;
        }

        // Check if the candidate "sees" the end point "in front".
        if (!end_detected)
        {
            retvalue.points->push_back(candidate);
            if ((delta_max_ >= metric_(end, candidate)) &&
                (total_prop * (end - candidate)) > 0)
                end_detected = true;
        }
        else
        {
            RealType cur_dist = (end - retvalue.points->back()).euclidean_norm();
            RealType candidate_dist1 = (candidate - retvalue.points->back()).euclidean_norm();
            RealType candidate_dist2 = (end - candidate).euclidean_norm();

            if (delta_min_ > metric_(end, candidate))
            {
                retvalue.points->push_back(candidate);
                candidate = end;
            }
            else if ((cur_dist > candidate_dist1) && (cur_dist > candidate_dist2))
            {
                retvalue.points->push_back(candidate);
            }
            else
            {
                candidate = end;
            }
        }

        // Update algortihm's state.
        Point3D previous = current;
        current = candidate;

//        // Compute inertial, centrifugal and tangential propagations.
//        Point3D inertial_prop = this_type::inertial_propagation_norm_(current, previous);
//        Point3D centrifugal_prop = this_type::centrifugal_propagation_norm(current,
//            main_plane_->origin());
//        Point3D total_tangential_prop = total_tangential_propagation_norm_(current,
//            inertial_prop);

//        // Compute total propagation.
//        total_prop = this_type::total_propagation_(inertial_prop, centrifugal_prop,
//            total_tangential_prop, inertial_weight_, centrifugal_weight_);

        total_prop = propagator_.next(current, previous);

    } while (current != end);

    return retvalue;
}

} // namespace surfaces
} // namespace bo

#endif // COMPLEX_PROPAGATION_HPP_D90ED351_6A45_4523_85F3_DA99F52B87C2
