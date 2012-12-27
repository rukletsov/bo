
/******************************************************************************

  parallel_planes_tiling.hpp, v 1.0.8 2012.12.27

  Implementation of several surface tiling methods, working with parallel planes.

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

#ifndef PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
#define PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_

#include <vector>
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
#include <intrin.h>

namespace bo {
namespace methods {
namespace surfaces {

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

    static Mesh propagate(ParallelPlaneConstPtr plane)
    {
        // TODO: assert on data size (at least 2 samples?).
        Mesh mesh(10);

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
            delta_min, delta_max, 1000, metric);

        // If the hole was detected, run propagation in a different direction.
        PropagationResult attempt2;
        if (attempt1.hole_encountered)
            attempt2 = propagate_(start, attempt1.points->back(), - initial_prop,
                tree, delta_min, delta_max, 1000, metric);

        // Glue propagation results together.
        for (ParallelPlane::const_reverse_iterator rit = attempt2.points->rbegin();
             rit != attempt2.points->rend(); ++rit)
            mesh.add_vertex(*rit);
        for (ParallelPlane::const_iterator it = attempt1.points->begin() + 1;
             it != attempt1.points->end(); ++it)
            mesh.add_vertex(*it);

        return mesh;
    }

    static Mesh christiansen_triangulation(ParallelPlaneConstPtr contour1,
                                           ParallelPlaneConstPtr contour2)
    {
        // Prepare: build kd-trees.
        Tree tree1(contour1->begin(), contour1->end(), std::ptr_fun(point3D_accessor_));
        Tree tree2(contour2->begin(), contour2->end(), std::ptr_fun(point3D_accessor_));

        // Determine whether both contours are closed or not. This determines the
        // algorithm's behaviour.

        // 1. Contours are closed.
        // Choose a vertex and a direction on the first contour.
        // Find closest vertex on the second contour and determine direction.
        // Iteratively sample points from contours until the surface patcg is ready.

        // 2. Contours having exactly one hole.
        // Take the first vertex (at the hole) and direction on the first contour.
        // Take the first vertex (at the hole) and determine co-directed movement.
        // Iterate till the end of both contours.
        ParallelPlane::const_iterator current1 = contour1->begin();
        Point3D direction1 = *(contour1->begin() + 1) - *current1;

        ParallelPlane::const_iterator current2 = contour2->begin();
        Point3D direction2 = *(contour2->begin() + 1) - *current2;

        // TODO: swap direction if dot product of direction is less than zero.
        std::size_t current1_idx = mesh.add_vertex(*current1);
        std::size_t current2_idx = mesh.add_vertex(*current2);
        while (true)
        {
            RealType span1_norm = (*(current1 + 1) - *current2).euclidean_norm();
            RealType span2_norm = (*(current2 + 1) - *current1).euclidean_norm();

            if (span1_norm > span2_norm)
            {

            }
        }

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
        PCAEngine::Result result = pca(neighbours);
        Point3D tangential = result.get<1>()[2];

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
                                       std::size_t max_size, Metric metric)
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
            std::pair<Tree::const_iterator, RealType> candidate_data =
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
            total_prop = this_type::total_propagation_(tangential_prop, inertial_prop, 0.5f);

        } while (current != end);

        return retvalue;
    }
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
