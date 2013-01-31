
/******************************************************************************

  triangulation.hpp, v 1.1.3 2013.01.31

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
#include <stdexcept>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "bo/mesh.hpp"
#include "bo/methods/distances_3d.hpp"
#include "bo/internal/surfaces/container_traversers.hpp"

namespace bo {
namespace methods {
namespace surfaces {

// A class performing different triangulations for two slices.
template <typename RealType>
class Triangulation: public boost::noncopyable
{
public:
    typedef Triangulation<RealType> this_type;

    typedef Vector<RealType, 3> Point3D;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;
    typedef Mesh<RealType> Mesh;
    typedef std::vector<Point3D> Contour;
    typedef boost::shared_ptr<const Contour> ContourConstPtr;

    typedef ContainerConstTraverser<Contour> ContourTraverser;
    typedef TraverseRuleFactory<Contour> TraverseFactory;

public:
    Triangulation(ContourConstPtr contour1, bool closed1,
                  ContourConstPtr contour2, bool closed2):
        contour1_(contour1), closed1_(closed1), contour2_(contour2), closed2_(closed2)
    {
        metric_ = &euclidean_distance<RealType, 3>;
    }

    // Implementation for Christiansen algorithm for closed and opened contours.
    //
    // If a contour is closed, it is traversed cyclical starting and ending in the
    // same vertex.
    //
    // If a contour has exactly one hole (i.e. opened), it is traversed once
    // from the start to the end (from one hole edge to the other).
    Mesh christiansen()
    {
        // Check contours' length.
        if ((contour1_->size() < 2) && (contour2_->size() < 2))
            throw std::logic_error("Cannot run Christiansen triangulation for contours "
                                   "consisting of less than 2 vertices.");

        // Create traversers for contours.
        ContourTraverser current1 = create_traverser1_();
        ContourTraverser candidate1 = current1 + 1;

        ContourTraverser current2 = create_traverser2_(current1);
        ContourTraverser candidate2 = current2 + 1;

        // Initialize output mesh.
        Mesh mesh(contour1_->size() + contour2_->size());
        std::size_t current1_idx = mesh.add_vertex(*current1);
        std::size_t current2_idx = mesh.add_vertex(*current2);

        // Iterate until both contours are exhausted.
        while (candidate1.is_valid() || candidate2.is_valid())
        {
            // Calculate span norms for candidate vertices.
            RealType span1_norm = span_norm(candidate1, current2);
            RealType span2_norm = span_norm(candidate2, current1);

            // Add candidate vertex.
            std::size_t candidate_idx;
            if (span1_norm > span2_norm)
            {
                candidate_idx = mesh.add_vertex(*candidate2);

                // Add edges to candidate vertex.
                mesh.add_face(Mesh::Face(current1_idx, candidate_idx, current2_idx));

                current2 = candidate2;
                current2_idx = candidate_idx;
                ++candidate2;
            }
            else
            {
                candidate_idx = mesh.add_vertex(*candidate1);

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
    // Creates an appropriate traverser for the first contour depending whether
    // it is closed or not.
    ContourTraverser create_traverser1_() const
    {
        ContourTraverser retvalue = closed1_ ?
            ContourTraverser(TraverseFactory::Create(contour1_, 0, true)) :
            ContourTraverser(TraverseFactory::Create(contour1_, true));
        return retvalue;
    }

    // Creates an appropriate traverser for the second contour depending whether
    // the contour is closed and how the first contour is oriented.
    ContourTraverser create_traverser2_(const ContourTraverser& traverser1) const
    {
        ContourTraverser retvalue;

        bool is_forward2 = true;
        if (closed2_)
        {
            // Find closest vertex on the second contour and determine direction.
            // No need of using kd-tree here, since the operation is done once.
            std::size_t c2min_idx = 0;
            RealType c2min_dist = std::numeric_limits<RealType>::max();

            for (std::size_t c2_idx = 0; c2_idx < contour2_->size(); ++c2_idx)
            {
                RealType cur_dist = metric_(contour2_->at(c2_idx), *traverser1);
                if (cur_dist < c2min_dist)
                {
                    c2min_idx = c2_idx;
                    c2min_dist = cur_dist;
                }
            }

            // Create a default directed traverser. This is need to calculate the
            // direction if c2min_idx is the last element in the collection.
            retvalue = ContourTraverser(TraverseFactory::Create(contour2_, c2min_idx, true));

            // Contour2 traverse direction should be swapped in order to correspond
            // with the contoru1 direction.
            Point3D direction1 = *(traverser1 + 1) - *traverser1;
            Point3D direction2 = *(retvalue + 1) - *retvalue;
            if (direction1 * direction2 < 0)
                retvalue = ContourTraverser(TraverseFactory::Create(contour2_, c2min_idx, false));
        }
        else
        {
            // Check if contour2_ traverse direction should be swapped in order to
            // correspond with the contour1_'s direction.
            if ((contour2_->back() - *traverser1).euclidean_norm() <
                    (contour2_->front() - *traverser1).euclidean_norm())
                is_forward2 = false;

            retvalue = ContourTraverser(TraverseFactory::Create(contour2_, is_forward2));
        }

        return retvalue;
    }

    // Calculates the euclidean norm of the span between candidate vertex on one
    // contour and current vertex on another. If candidate vertex is invalid (this
    // indicates that the corresponding contour has been exhausted), returns infinity.
    // This guarantees that vertices will be sampled solely from the other contour.
    static RealType span_norm(const ContourTraverser& candidate,
                              const ContourTraverser& other_current)
    {
        RealType norm = candidate.is_valid() ?
                    (*candidate - *other_current).euclidean_norm() :
                    std::numeric_limits<RealType>::max();
        return norm;
    }

private:
    ContourConstPtr contour1_;
    bool closed1_;
    ContourConstPtr contour2_;
    bool closed2_;

    Metric metric_;
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // TRIANGULATION_HPP_40CA4CA1_3752_4080_9BDE_D13393E21902
