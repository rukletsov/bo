
/******************************************************************************

  parallel_planes_tiling.hpp, v 1.0.1 2012.12.03

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
#include <algorithm>
#include <functional>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "bo/config.hpp"
#include "bo/mesh.hpp"
#include "bo/raw_image_2d.hpp"
#include "bo/kdtree.hpp"
#include "bo/extended_math.hpp"
#include "bo/blas/blas.hpp"
#include "bo/blas/conversions.hpp"
#include "bo/extended_std.hpp"

namespace bo {
namespace methods {
namespace surfaces {

template <typename RealType>
class MinSpanPropagation: public boost::noncopyable
{
public:
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> Points3D;
    typedef std::vector<Point3D> ParallelPlane;
    typedef Mesh<RealType> Mesh;
    typedef RawImage2D<RealType> Image2D;
    typedef boost::shared_ptr<ParallelPlane> ParallelPlanePtr;
    typedef boost::shared_ptr<const ParallelPlane> ParallelPlaneConstPtr;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;

    typedef bo::KDTree<3, Point3D,
        boost::function<RealType (Point3D, std::size_t)> > Tree;

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
        float delta_min = 5;
        float delta_max = 10;

        // Build kd-tree from given points.
        Tree tree(std::ptr_fun(point3D_accessor_));
        for (ParallelPlane::const_iterator it = plane->begin(); it != plane->end(); ++it)
             tree.insert(*it);
        tree.optimise();

        // Choose initial point.
        Point3D current = (*plane)[10];

        // Choose propagation direction. First search for nearby points.
        Points3D neighbours;
        tree.find_within_range(current, delta_max, std::back_inserter(neighbours));

        // PCA.

        // Now find the mean among the neighbours.
        Point3D mean_neighbour = bo::mean(neighbours);

        // Now calculate the deviations from mean.
        std::transform(neighbours.begin(), neighbours.end(), neighbours.begin(),
                       std::bind2nd(std::minus<Point3D>(), mean_neighbour));

        // Now compute the 3x3 covariance matrix columwise and store it in blas matrix
        // type for further processing.
        typedef blas::bounded_matrix<RealType, 3, 3> CovarianceMatrix;
        CovarianceMatrix cov = blas::zero_matrix<RealType>(3);

        for (Points3D::const_iterator pt = neighbours.begin();
             pt != neighbours.end(); ++pt)
        {
            cov(0, 0) += pt->x() * pt->x();
            cov(1, 0) += pt->y() * pt->x();
            cov(2, 0) += pt->z() * pt->x();

            cov(0, 1) += pt->x() * pt->y();
            cov(1, 1) += pt->y() * pt->y();
            cov(2, 1) += pt->z() * pt->y();

            cov(0, 2) += pt->x() * pt->z();
            cov(1, 2) += pt->y() * pt->z();
            cov(2, 2) += pt->z() * pt->z();
        }

        cov /= neighbours.size();

        // Get the eigenvectors of the covariance matrix.
        std::vector<RealType> cov_eigenvalues = blas::eigen_symmetric(cov);

        // Choose the largest eigenvector, which should be the last in cov matrix.
        Point3D prop_direction;
        prop_direction.x() = cov(0, 2); prop_direction.y() = cov(1, 2); prop_direction.z() = cov(2, 2);
        std::cout << cov_eigenvalues;
        std::cout << prop_direction;





        for (Points3D::const_iterator it = neighbours.begin(); it != neighbours.end(); ++it)
             mesh.add_vertex(*it);

        return mesh;
    }

private:
    // Helper function for KDTree instance.
    static inline RealType point3D_accessor_(Point3D pt, std::size_t k)
    {
        return pt[k];
    }
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // PARALLEL_PLANES_TILING_HPP_353B5678_8A01_4091_91C4_9C5BE2476BA0_
