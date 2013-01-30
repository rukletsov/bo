
/******************************************************************************

  triangulation.hpp, v 1.1.1 2013.01.30

  Triangulation algorithms for surface reconstruction problems.

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

#ifndef TRIANGULATION_HPP_40CA4CA1_3752_4080_9BDE_D13393E21902
#define TRIANGULATION_HPP_40CA4CA1_3752_4080_9BDE_D13393E21902

#include <vector>
#include <limits>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "bo/mesh.hpp"
#include "bo/methods/distances_3d.hpp"
#include "bo/internal/surfaces/container_traversers.hpp"

namespace bo {
namespace methods {
namespace surfaces {

template <typename RealType>
class Triangulation: public boost::noncopyable
{
public:
    typedef Triangulation<RealType> this_type;

    typedef Mesh<RealType> Mesh;
    typedef Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> ParallelPlane;
    typedef boost::shared_ptr<const ParallelPlane> ParallelPlaneConstPtr;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;

public:
    Triangulation(ParallelPlaneConstPtr contour1, bool closed1,
                  ParallelPlaneConstPtr contour2, bool closed2):
        contour1_(contour1), closed1_(closed1), contour2_(contour2), closed2_(closed2)
    {
        metric_ = &euclidean_distance<RealType, 3>;
    }

    // TODO: put this into triangulation.hpp.
    Mesh christiansen()
    {
        // TODO: assert on data size (at least 2 samples?).

        // 1. Contours are closed.
        // Choose a vertex and a direction on the first contour.
        // Find closest vertex on the second contour and determine direction.
        // Iteratively sample points from contours until the surface patch is ready.

        // 2. Contours having exactly one hole.
        // Take the first vertex (at the hole) and direction on the first contour.
        // Take the first vertex (at the hole) and determine co-directed movement.
        // Iterate till the end of both contours.
        typedef ContainerConstTraverser<ParallelPlane> ContourTraverser;
        typedef TraverseRuleFactory<ParallelPlane> Factory;

        // Create appropriate traverser for the first contour depending whether
        // it is closed or not.
        ContourTraverser current1 = closed1_ ?
                    ContourTraverser(Factory::Create(contour1_, 0, true)) :
                    ContourTraverser(Factory::Create(contour1_, true));
        ContourTraverser candidate1 = current1 + 1;

        // Create appropriate traverser for the second contour depending whether
        // the contour is closed and how the first contour is oriented.
        ContourTraverser current2;
        bool is_forward2 = true;
        if (closed2_)
        {
            // Find closest vertex on the second contour and determine direction.
            // No need of using kd-tree here, since the operation is done once.
            std::size_t c2min_idx = 0;
            RealType c2min_dist = std::numeric_limits<RealType>::max();

            for (std::size_t c2_idx = 0; c2_idx < contour2_->size(); ++c2_idx)
            {
                RealType cur_dist = metric_(contour2_->at(c2_idx), *current1);
                if (cur_dist < c2min_dist)
                {
                    c2min_idx = c2_idx;
                    c2min_dist = cur_dist;
                }
            }

            // Create a default directed traverser. This is need to calculate the
            // direction if c2min_idx is the last element in the collection.
            current2 = ContourTraverser(Factory::Create(contour2_, c2min_idx, true));

            // Contour2 traverse direction should be swapped in order to correspond
            // with the contoru1 direction.
            Point3D direction1 = *(current1 + 1) - *current1;
            Point3D direction2 = *(current2 + 1) - *current2;
            if (direction1 * direction2 < 0)
                current2 = ContourTraverser(Factory::Create(contour2_, c2min_idx, false));
        }
        else
        {
            // Check if contour2_ traverse direction should be swapped in order to
            // correspond with the contour1_'s direction.
            if ((contour2_->back() - *current1).euclidean_norm() <
                    (contour2_->front() - *current1).euclidean_norm())
                is_forward2 = false;

            current2 = ContourTraverser(Factory::Create(contour2_, is_forward2));
        }

        ContourTraverser candidate2 = current2 + 1;

        Mesh mesh(contour1_->size() + contour2_->size());
        std::size_t current1_idx = mesh.add_vertex(*current1);
        std::size_t current2_idx = mesh.add_vertex(*current2);
        while (true)
        {
            // This guarantees that when one contour ends, points will be sampled
            // solely from the other one.
            RealType span1_norm = (candidate1.is_valid()) ?
                        (*(candidate1) - *current2).euclidean_norm() :
                        std::numeric_limits<RealType>::max();
            RealType span2_norm = (candidate2.is_valid()) ?
                        (*(candidate2) - *current1).euclidean_norm() :
                        std::numeric_limits<RealType>::max();

            // Means both contours are exhausted.
            if (!candidate1.is_valid() && !candidate2.is_valid())
                break;

            // Add candidate vertex.
            std::size_t candidate_idx;
            if (span1_norm > span2_norm)
            {
                candidate_idx = candidate2.is_valid() ?
                            mesh.add_vertex(*candidate2) : -1;

                // Add edges to candidate vertex.
                mesh.add_face(Mesh::Face(current1_idx, candidate_idx, current2_idx));

                current2 = candidate2;
                current2_idx = candidate_idx;
                ++candidate2;
            }
            else
            {
                candidate_idx = candidate1.is_valid() ?
                            mesh.add_vertex(*candidate1) : -1;
                // Add edges to candidate vertex.
                mesh.add_face(Mesh::Face(current1_idx, candidate_idx, current2_idx));

                current1 = candidate1;
                current1_idx = candidate_idx;
                ++candidate1;
            }
        }

        return mesh;
    }

private:
    ParallelPlaneConstPtr contour1_;
    bool closed1_;
    ParallelPlaneConstPtr contour2_;
    bool closed2_;

    Metric metric_;
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // TRIANGULATION_HPP_40CA4CA1_3752_4080_9BDE_D13393E21902
