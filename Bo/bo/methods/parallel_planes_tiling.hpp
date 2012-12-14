
/******************************************************************************

  parallel_planes_tiling.hpp, v 1.0.2 2012.12.13

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
        Mesh mesh(10);

        // Set algorithm parameters.
        RealType delta_min = 5;
        RealType delta_max = 10;
        Metric metric(&euclidean_distance<RealType, 3>);

        // Build kd-tree from given points.
        // TODO: provide kd-tree with current metric.
        Tree tree(plane->begin(), plane->end(), std::ptr_fun(point3D_accessor_));

        // Choose initial point and mark it as the start one. Previous to current
        // equals current to compute inertial propagation right.
        Point3D start = (*plane)[10];
        Point3D current = start;
        Point3D previous = current;

        // Get the total propagation.
        Point3D tangential_prop = this_type::tangential_propagation_(current, tree, delta_max);
        Point3D inertial_prop = this_type::inertial_propagation_(current, previous);
        Point3D total_prop = this_type::total_propagation_(tangential_prop, inertial_prop, 0.5);

        // Initial propagation equals tangential propagation alone.
        Point3D prop = tangential_prop;
        std::cout << prop;

        // Search for the propagation candidate. It should lie on the arced strip
        // bounded by delta_min and delta_max circumferences and a plane containing
        // current point and normal to propagation vector.
        Point3D candidate = *(tree.find_nearest_if(current, delta_max,
            ArchedStrip(current, delta_min, delta_max, prop, metric)).first);


//        do
//        {
//            // Compute tangential propagation
//            // Compute inertial propagation
//            // Derive propagation vector.
//            // Get candidate.
//            // Check if the candidate "sees" start point.
//            // Link either to candidate or to start point.
//            // Update mesh.

//        } while (curret != start);

        std::cout << current;
        std::cout << candidate;






        // Build mesh from points.
        for (Points3D::const_iterator it = plane->begin(); it != plane->end(); ++it)
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

            RealType dist = dist_fun_/*euclidean_distance*/(ref_pt_, pt);
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

private:
    // Helper function for KDTree instance.
    static inline RealType point3D_accessor_(Point3D pt, std::size_t k)
    {
        return pt[k];
    }

    // Computes the tangential propagation vector for the given point.
    static Point3D tangential_propagation_(Point3D pt, const Tree& tree, RealType radius)
    {
        // Search for nearby points.
        Points3D neighbours;
        tree.find_within_range(pt, radius, std::back_inserter(neighbours));

        // Employ PCA to extract the tangential propagation vector from the set of
        // nearby points.
        PCAEngine pca;
        PCAEngine::Result result = pca(neighbours);
        Point3D tangential = result.get<1>()[2];

        return tangential;
    }

    // Computes the inertial propagation vector from current and previous mesh vertices.
    static Point3D inertial_propagation_(Point3D current, Point3D previous)
    {
        return
            (current -= previous);
    }

    // Computes the total propagation vector from tangential and inertial components.
    static Point3D total_propagation_(Point3D tangential, Point3D inertial, RealType ratio)
    {
        Point3D total = tangential * ratio + inertial * (RealType(1) - ratio);
        return total/*.maximum_norm()*/;
    }
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
