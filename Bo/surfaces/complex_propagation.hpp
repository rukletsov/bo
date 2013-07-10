
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
#include <algorithm>
#include <functional>
#include <boost/assert.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "bo/core/raw_image_2d.hpp"
#include "bo/core/mesh.hpp"
#include "bo/core/kdtree.hpp"
#include "bo/math/functions.hpp"
#include "bo/math/pca.hpp"
#include "bo/math/mean.hpp"
#include "bo/distances/distances_3d.hpp"
#include "bo/surfaces/detail/types.hpp"
#include "bo/surfaces/detail/arched_strip.hpp"

namespace bo {
namespace surfaces {

template <typename RealType>
class ComplexPropagation: public boost::noncopyable
{
public:
    typedef ComplexPropagation<RealType> SelfType;
    typedef boost::shared_ptr<SelfType> Ptr;
    typedef std::vector<Ptr> Ptrs;

    typedef Vector<RealType, 3> Point3D;

    // Types for input data.
    typedef std::vector<Point3D> Points3D;
    typedef boost::shared_ptr<Points3D> Points3DPtr;
    typedef boost::shared_ptr<const Points3D> Points3DConstPtr;
    typedef std::vector<Points3DConstPtr> Points3DConstPtrs;

    typedef std::vector<RealType> Weights;

    // Used by factory functions.
    typedef RawImage2D<RealType> Image2D;
    typedef bo::Mesh<RealType> Mesh;

private:
    typedef std::pointer_to_binary_function<const Point3D&, const Point3D&, RealType> Metric;
    typedef KDTree<3, Point3D, std::pointer_to_binary_function<const Point3D&,
            std::size_t, RealType> > Tree;
    typedef boost::shared_ptr<Tree> TreePtr;
    typedef std::vector<TreePtr> TreePtrs;

    typedef detail::PropagationDirection<RealType, Tree> PropagationDirector;
    typedef detail::PropagationResult<RealType> PropagationResult;
    typedef typename PropagationResult::PropagatedContour PropagatedContour;
    typedef boost::shared_ptr<PropagatedContour> PropagatedContourPtr;

public:
    // Factory functions. Create an instance of the class from either of the supported
    // input with the provided parameters.
    static Ptr create(const Points3D& plane, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius);

    static Ptrs create(const Points3DConstPtrs& planes, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius);

    static Ptr from_mesh(const Mesh& mesh, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius);

    static Ptr from_raw_image(const Image2D& data, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius);

    // A special factory. Creates a set of instances one per given plane. Additionally
    // passes a set of neighbours with corresponding weights to each instance.
    static Ptrs create(const Points3DConstPtrs& planes, const Weights& neighbour_weights,
        RealType delta_min, RealType delta_max, RealType inertial_weight,
        RealType centrifugal_weight, RealType tangential_radius);

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
    PropagatedContourPtr contour() const;

protected:
    ComplexPropagation(RealType delta_min, RealType delta_max, TreePtr tree_ptr,
                       const PropagationDirector& director, const Point3D& start_point);

    // Helper function for KDTree instance.
    static RealType point3D_accessor_(const Point3D& pt, std::size_t k);

    //
    static TreePtr tree_from_plane_(const Points3D& plane);

    // Adds neighbour planes with corresponding weights. These planes are used in
    // the calculation of the tangential component of the propagation, making it
    // dependent of the neighbour planes. This should lead to the desirable smoothing.
    static void populate_neighbours_and_weights_(std::size_t plane_idx,
            const Points3DConstPtrs& planes, const Weights& neighbour_weights,
            Points3DConstPtrs& out_neighbours, Weights& out_weights);

    //
    static Point3D get_plane_normal_(const Points3D& plane);

    //
    static Points3DConstPtrs project_neighbour_onto_plane_(const Point3D& center_of_mass,
            const Point3D& normal, const Points3DConstPtrs& neighbours);

    // Tries performing propagation from the start point to the end point.
    // Note that end point is not included in result, while start is always included.
    PropagationResult propagate_(const Point3D& start, const Point3D& end,
                                 Point3D total_prop, std::size_t max_size);

private:
    // Parameters and options of the propagation algorithm.
    RealType delta_min_;
    RealType delta_max_;
    Metric metric_;

    // Algorithm data and helper obejcts.
    TreePtr tree_ptr_;
    PropagationDirector propagation_director_;
    Point3D start_;

    // Markers and final contour.
    bool stopped_;
    bool has_hole_;
    PropagatedContourPtr contour_;
};


// C-tor.
template <typename RealType>
ComplexPropagation<RealType>::ComplexPropagation(RealType delta_min,
                                                 RealType delta_max,
                                                 TreePtr tree_ptr,
                                                 const PropagationDirector& director,
                                                 const Point3D& start_point):
    delta_min_(delta_min), delta_max_(delta_max),
    metric_(std::ptr_fun(distances::euclidean_distance<RealType, 3>)),
    tree_ptr_(tree_ptr),
    propagation_director_(director),
    start_(start_point),
    stopped_(false), has_hole_(false),
    contour_(boost::make_shared<PropagatedContour>())
{
    // If points number is less than 2, propagation cannot be initialized. And
    // because there is no sense in doing propagation for 0 or 1 points, we throw
    // an error instead of returning the plane unchanged.
    if (tree_ptr->size() < 2)
        throw std::logic_error("Cannot run propagation for planes consisting of "
                               "less than 2 vertices.");
}


// Factories. The sum of inertial and centrifugal weights should lie in [0; 1].
template <typename RealType>
typename ComplexPropagation<RealType>::Ptr ComplexPropagation<RealType>::create(
        const Points3D& plane, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius)
{
    // Build kd-tree from the given points.
    TreePtr tree_ptr = tree_from_plane_(plane);

    PropagationDirector direction;
    if (math::check_small(centrifugal_weight))
    {
        direction = PropagationDirector::create_simple(inertial_weight, tree_ptr, tangential_radius);
    }
    else
    {
        Point3D center_of_mass = bo::math::mean(plane);
        direction = PropagationDirector::create_with_centrifugal(inertial_weight,
                centrifugal_weight, tree_ptr, tangential_radius, center_of_mass);
    }

    // C-tor is declared private, using boost::make_shared gets complicated.
    Ptr ptr(new SelfType(delta_min, delta_max, tree_ptr, direction, plane.front()));
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
        Ptr propagator = create(*(planes[idx]), delta_min, delta_max, inertial_weight,
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
    return
        create(mesh.get_all_vertices(), delta_min, delta_max, inertial_weight,
               centrifugal_weight, tangential_radius);
}

template <typename RealType>
typename ComplexPropagation<RealType>::Ptr ComplexPropagation<RealType>::from_raw_image(
        const Image2D& data, RealType delta_min, RealType delta_max,
        RealType inertial_weight, RealType centrifugal_weight, RealType tangential_radius)
{
    // Load plane data from an istance of RawImage2D.
    Points3D plane;

    for (std::size_t row = 0; row < data.height(); ++row)
        for (std::size_t col = 0; col < data.width(); ++col)
            if (data(col, row) > 0)
                plane.push_back(Point3D(RealType(col), RealType(row), RealType(0)));

    return
        create(plane, delta_min, delta_max, inertial_weight, centrifugal_weight, tangential_radius);
}

template <typename RealType>
typename ComplexPropagation<RealType>::Ptrs ComplexPropagation<RealType>::create(
        const Points3DConstPtrs& planes, const Weights& neighbour_weights,
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
        neighbours.reserve(2 * radius);

        Weights weights;
        weights.reserve(2 * radius);

        //
        populate_neighbours_and_weights_(idx, planes, neighbour_weights, neighbours, weights);

        Points3DConstPtr plane = planes[idx];

        // Cache main plane origin and normal.
        Point3D center_of_mass = bo::math::mean(*plane);
        Point3D norm = get_plane_normal_(*plane);

        // Project all neighbours to the current plane.
        Points3DConstPtrs neighbour_projs = project_neighbour_onto_plane_(center_of_mass,
                norm, neighbours);
        neighbour_projs.reserve(neighbours.size());

        TreePtr main_tree_ptr = tree_from_plane_(*plane);

        // Compute kd-trees for projected neighbour planes.
        TreePtrs neighbour_trees;
        neighbour_trees.reserve(neighbour_projs.size());
        for (typename Points3DConstPtrs::const_iterator plane_it =
             neighbour_projs.begin(); plane_it != neighbour_projs.end(); ++plane_it)
            neighbour_trees.push_back(tree_from_plane_(**plane_it));

        PropagationDirector direction = math::check_small(centrifugal_weight)
                ? PropagationDirector::create_with_neighbours(inertial_weight,
                        main_tree_ptr, tangential_radius, neighbour_trees, weights)
                : PropagationDirector::create_with_neighbours_and_centrifugal(
                        inertial_weight, centrifugal_weight, main_tree_ptr, tangential_radius,
                        center_of_mass, neighbour_trees, weights);

        // C-tor is declared private, using boost::make_shared gets complicated.
        Ptr ptr(new SelfType(delta_min, delta_max, main_tree_ptr, direction, plane->front()));

        propagators.push_back(ptr);
    }

    return propagators;
}


// Public control functions.
template <typename RealType>
void ComplexPropagation<RealType>::propagate()
{
    // Choose initial point and initial propagation. It solely consists of the
    // tangential component, since inertial cannot be defined.
    Point3D initial_prop = propagation_director_.initial(start_);

    // Run propagation. It may bump into a hole or return a circuit.
    PropagationResult attempt1 = propagate_(start_, start_, initial_prop, 1000);

    // If the hole was detected, run propagation in a different direction.
    PropagationResult attempt2;
    if (attempt1.has_hole)
        attempt2 = propagate_(start_, attempt1.points.back(), - initial_prop, 1000);

    // Glue propagation results together.
    contour_->reserve(attempt1.points.size() + attempt2.points.size());
    for (typename PropagatedContour::const_reverse_iterator rit = attempt2.points.rbegin();
         rit != attempt2.points.rend(); ++rit)
        contour_->push_back(*rit);
    for (typename PropagatedContour::const_iterator it = attempt1.points.begin() + 1;
         it != attempt1.points.end(); ++it)
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
typename ComplexPropagation<RealType>::PropagatedContourPtr
ComplexPropagation<RealType>::contour() const
{
    return contour_;
}


// Private static helper functions.
template <typename RealType> inline
RealType ComplexPropagation<RealType>::point3D_accessor_(const Point3D &pt, std::size_t k)
{
    return pt[k];
}

template <typename RealType> inline
typename ComplexPropagation<RealType>::TreePtr
ComplexPropagation<RealType>::tree_from_plane_(const Points3D& plane)
{
    // Build kd-tree from the given points.
    // TODO: provide kd-tree with current metric?
    TreePtr tree_ptr = boost::make_shared<Tree>(plane.begin(), plane.end(),
                                                std::ptr_fun(point3D_accessor_));
    return tree_ptr;
}

template <typename RealType>
typename ComplexPropagation<RealType>::Point3D
ComplexPropagation<RealType>::get_plane_normal_(const Points3D& plane)
{
    // Employ PCA to estimate plane normal.
    typedef math::PCA<RealType, 3> PCAEngine;
    typedef typename PCAEngine::Result PCAResult;
    PCAEngine pca;
    PCAResult result = pca(plane);
    Point3D normal = result.template get<1>()[0];

    return normal;
}

template <typename RealType>
typename ComplexPropagation<RealType>::Points3DConstPtrs
ComplexPropagation<RealType>::project_neighbour_onto_plane_(const Point3D& center_of_mass,
        const Point3D& normal, const Points3DConstPtrs& neighbours)
{
    // Project all neighbours to the current plane.
    // TODO: rewrite using transform?
    Points3DConstPtrs neighbour_projs;
    neighbour_projs.reserve(neighbours.size());

    for (typename Points3DConstPtrs::const_iterator neighbour =
         neighbours.begin(); neighbour != neighbours.end(); ++neighbour)
    {
        Points3DPtr plane_proj = boost::make_shared<Points3D>();
        plane_proj->reserve((*neighbour)->size());

        for (typename Points3D::const_iterator point_it = (*neighbour)->begin();
             point_it != (*neighbour)->end(); ++point_it)
        {
            plane_proj->push_back(distances::project_point_onto_plane(*point_it,
                    center_of_mass, normal));
        }

        neighbour_projs.push_back(plane_proj);
    }

    return neighbour_projs;
}

template <typename RealType>
void ComplexPropagation<RealType>::populate_neighbours_and_weights_(std::size_t plane_idx,
        const Points3DConstPtrs& planes, const Weights& neighbour_weights,
        Points3DConstPtrs& out_neighbours, Weights& out_weights)
{
    std::size_t planes_count = planes.size();
    std::size_t radius = neighbour_weights.size();

    // Try add neighbours before the current idx. We add radius items at most,
    // but less if no neighbours are available (close to container front).
    std::size_t neighbour_idx = plane_idx - 1;
    std::size_t added_count = 0;
    std::size_t maxval = std::size_t(-1); // dirty hack to workaround int behaviour for size_t.
    while ((neighbour_idx < maxval) && (added_count < radius))
    {
        out_neighbours.push_back(planes[neighbour_idx]);
        out_weights.push_back(neighbour_weights[added_count]);
        ++added_count;
        --neighbour_idx;
    }

    // Try add neighbours after the current idx. We add radius items at most,
    // but less if no neighbours are available (close to container end).
    neighbour_idx = plane_idx + 1;
    added_count = 0;
    while ((neighbour_idx < planes_count) && (added_count < radius))
    {
        out_neighbours.push_back(planes[neighbour_idx]);
        out_weights.push_back(neighbour_weights[added_count]);
        ++added_count;
        ++neighbour_idx;
    }

    BOOST_ASSERT((out_neighbours.size() == out_weights.size()) && "Weights quantity doesn't"
                 "correspond to planes number.");
}

template <typename RealType>
typename ComplexPropagation<RealType>::PropagationResult
ComplexPropagation<RealType>::propagate_(const Point3D& start, const Point3D& end,
                                         Point3D total_prop, std::size_t max_size)
{
    typedef detail::ArchedStrip<RealType, 3> ArchedStrip;

    PropagationResult retvalue;
    retvalue.points.push_back(start);

    // Flag for so-called "smooth-ending". It is necessary to keep the distance
    // between the points in the end phase as close to delta_min as possible,
    // despite a delta_max point has been already found.
    bool end_detected = false;

    Point3D current = start;
    do
    {
        // Restrict total length (to prevent looping).
        if (retvalue.points.size() > max_size)
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
                tree_ptr_->find_nearest_if(phantom_candidate, delta_max_,
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
            retvalue.points.push_back(candidate);
            if ((delta_max_ >= metric_(end, candidate)) &&
                (total_prop * (end - candidate)) > 0)
                end_detected = true;
        }
        else
        {
            RealType cur_dist = (end - retvalue.points.back()).euclidean_norm();
            RealType candidate_dist1 = (candidate - retvalue.points.back()).euclidean_norm();
            RealType candidate_dist2 = (end - candidate).euclidean_norm();

            if (delta_min_ > metric_(end, candidate))
            {
                retvalue.points.push_back(candidate);
                candidate = end;
            }
            else if ((cur_dist > candidate_dist1) && (cur_dist > candidate_dist2))
            {
                retvalue.points.push_back(candidate);
            }
            else
            {
                candidate = end;
            }
        }

        // Update algortihm's state.
        Point3D previous = current;
        current = candidate;
        total_prop = propagation_director_.next(current, previous);

    } while (current != end);

    return retvalue;
}

} // namespace surfaces
} // namespace bo

#endif // COMPLEX_PROPAGATION_HPP_D90ED351_6A45_4523_85F3_DA99F52B87C2
