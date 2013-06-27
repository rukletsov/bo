
/******************************************************************************

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
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "bo/mesh.hpp"
#include "bo/methods/distances_3d.hpp"
#include "bo/internal/surfaces/container_traversers.hpp"
#include "bo/internal/surfaces/indexed_tstrip.hpp"

namespace bo {
namespace methods {
namespace surfaces {

namespace detail {

// A struct incorporating contour data. Instances can be safely swapped.
template <typename PointsCont>
struct ContourDescriptor
{
    typedef boost::shared_ptr<PointsCont> PointsContPtr;

    ContourDescriptor(std::size_t id, PointsContPtr points_ptr, bool is_closed)
        : ID(id), PointsPtr(points_ptr), IsClosed(is_closed), Size(points_ptr->size())
    { }

    void swap(ContourDescriptor& other)
    {
        // Clone data before swapping in order not to spoil original contours.
        PointsContPtr new_this(new PointsCont(*PointsPtr));
        PointsContPtr new_other(new PointsCont(*other.PointsPtr));
        PointsPtr = new_other;
        other.PointsPtr = new_this;

        std::swap(ID, other.ID);
        std::swap(IsClosed, other.IsClosed);
        std::swap(Size, other.Size);
    }

    std::size_t ID;
    PointsContPtr PointsPtr;
    bool IsClosed;
    std::size_t Size;
};

} // namespace detail


// A class performing different triangulations for two slices.
template <typename RealType>
class ChristiansenTiling: public boost::noncopyable
{
public:
    typedef ChristiansenTiling<RealType> this_type;

    typedef Vector<RealType, 3> Point3D;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;
    typedef std::vector<Point3D> Contour;
    typedef boost::shared_ptr<Contour> ContourPtr;
    typedef bo::Mesh<RealType> Mesh3D;
    typedef detail::ContourDescriptor<Contour> ContourDescriptor;

    typedef ContainerConstTraverser<Contour> ContourTraverser;
    typedef TraverseRuleFactory<Contour> TraverseFactory;

public:
    ChristiansenTiling(std::size_t id1, ContourPtr contour1, bool closed1,
                  std::size_t id2, ContourPtr contour2, bool closed2)
        : contour_descr1_(id1, contour1, closed1), contour_descr2_(id2, contour2, closed2)
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
    IndexedTStrip christiansen()
    {
        // Check contours' length.
        if ((contour_descr1_.Size < 2) && (contour_descr2_.Size < 2))
            throw std::logic_error("Cannot run Christiansen triangulation for contours "
                                   "consisting of less than 2 vertices.");

        // If the first contour is closed but the second is opened, traverse direction
        // for the second contour may be determined incorrectly.
        if (contour_descr1_.IsClosed && !contour_descr2_.IsClosed)
            contour_descr1_.swap(contour_descr2_);

        // Create traversers for contours.
        ContourTraverser current1 = create_traverser1_();
        ContourTraverser candidate1 = current1 + 1;

        ContourTraverser current2 = create_traverser2_(current1);
        ContourTraverser candidate2 = current2 + 1;

        // Initialize output mesh.
        IndexedTStrip tstrip(contour_descr1_.ID, current1.index(),
                             contour_descr2_.ID, current2.index(),
                             contour_descr1_.Size + contour_descr2_.Size);

        // Iterate until both contours are exhausted.
        while (candidate1.is_valid() || candidate2.is_valid())
        {
            // Calculate span norms for candidate vertices.
            RealType span1_norm = span_norm(candidate1, current2);
            RealType span2_norm = span_norm(candidate2, current1);

            // Choose and add candidate vertex and corresponding face.
            if (span1_norm > span2_norm)
            {
                tstrip.add2(candidate2.index());
                current2 = candidate2;
                ++candidate2;
            }
            else
            {
                tstrip.add1(candidate1.index());
                current1 = candidate1;
                ++candidate1;
            }
        }

        return tstrip;
    }

private:
        // Creates an appropriate traverser for the first contour depending whether
    // it is closed or not.
    ContourTraverser create_traverser1_() const
    {
        ContourTraverser retvalue = contour_descr1_.IsClosed ?
            ContourTraverser(TraverseFactory::Create(contour_descr1_.PointsPtr, 0, true)) :
            ContourTraverser(TraverseFactory::Create(contour_descr1_.PointsPtr, true));
        return retvalue;
    }

    // Creates an appropriate traverser for the second contour depending whether
    // the contour is closed and how the first contour is oriented.
    ContourTraverser create_traverser2_(const ContourTraverser& traverser1) const
    {
        ContourTraverser retvalue;

        if (contour_descr2_.IsClosed)
        {
            // Find closest vertex on the second contour and determine direction.
            // No need of using kd-tree here, since the operation is done once.
            std::size_t c2min_idx = 0;
            RealType c2min_dist = std::numeric_limits<RealType>::max();

            for (std::size_t c2_idx = 0; c2_idx < contour_descr2_.Size; ++c2_idx)
            {
                RealType cur_dist = metric_(contour_descr2_.PointsPtr->at(c2_idx), *traverser1);
                if (cur_dist < c2min_dist)
                {
                    c2min_idx = c2_idx;
                    c2min_dist = cur_dist;
                }
            }

            // Create a default directed traverser. This is need to calculate the
            // direction if c2min_idx is the last element in the collection.
            retvalue = ContourTraverser(TraverseFactory::Create(contour_descr2_.PointsPtr,
                                                                c2min_idx, true));

            // Contour2 traverse direction should be swapped in order to correspond
            // with the contour1 direction.
            Point3D direction1 = *(traverser1 + 1) - *traverser1;
            Point3D direction2 = *(retvalue + 1) - *retvalue;
            if (direction1 * direction2 < 0)
                retvalue = ContourTraverser(TraverseFactory::Create(
                        contour_descr2_.PointsPtr, c2min_idx, false));
        }
        else
        {
            // Check if contour2_ traverse direction should be swapped in order to
            // correspond with the contour1_'s direction.
            RealType dist_to_first = metric_(contour_descr2_.PointsPtr->front(), *traverser1);
            RealType dist_to_last = metric_(contour_descr2_.PointsPtr->back(), *traverser1);
            bool is_forward2 = (dist_to_first < dist_to_last) ? true : false;

            retvalue = ContourTraverser(TraverseFactory::Create(contour_descr2_.PointsPtr,
                                                                is_forward2));
        }

        return retvalue;
    }

    // Calculates the euclidean norm of the span between candidate vertex on one
    // contour and current vertex on another. If candidate vertex is invalid (this
    // indicates that the corresponding contour has been exhausted), returns infinity.
    // This guarantees that vertices will be sampled solely from the other contour.
    RealType span_norm(const ContourTraverser& candidate,
                       const ContourTraverser& other_current) const
    {
        RealType norm = candidate.is_valid() ?
                    metric_(*candidate, *other_current):
                    std::numeric_limits<RealType>::max();
        return norm;
    }

private:
    ContourDescriptor contour_descr1_;
    ContourDescriptor contour_descr2_;

    Metric metric_;
};


template <typename RealType, typename ContourTypePtr>
bo::Mesh<RealType> christiansen(const ContourTypePtr& contour1, const ContourTypePtr& contour2)
{
    typedef bo::Mesh<RealType> Mesh3D;

    IndexedTStrip tstrip = ChristiansenTiling<RealType>(
            0, contour1->contour(), contour1->is_closed(),
            1, contour2->contour(), contour2->is_closed())
            .christiansen();

    // Calculate the total number of vertices to create mesh efficiently.
    std::size_t total_vertices = contour1->contour()->size() + contour2->contour()->size();

    // Populate the result mesh with vertices and store vertices shift for the
    // second contour (it should be 0 for the first).
    Mesh3D mesh(total_vertices);
    std::vector<std::size_t> lookup_table(2);
    lookup_table[0] = mesh.add_vertices(*(contour1->contour()));
    lookup_table[1] = mesh.add_vertices(*(contour2->contour()));

    // Populate the result mesh with faces. Transform old <contour_id, vertex_ix>
    // keys into new <mesh_vertex_id> keys using lookup table.
    BOOST_FOREACH (const IndexedTStrip::Face& face, tstrip.get_faces())
    {
        std::size_t new_a = lookup_table[face.A().first] + face.A().second;
        std::size_t new_b = lookup_table[face.B().first] + face.B().second;
        std::size_t new_c = lookup_table[face.C().first] + face.C().second;
        mesh.add_face(typename Mesh3D::Face(new_a, new_b, new_c));
    }

    return mesh;
}

template <typename RealType, typename ContourTypePtr>
bo::Mesh<RealType> christiansen(const std::vector<ContourTypePtr>& contours)
{
    typedef bo::Mesh<RealType> Mesh3D;

    std::size_t total_contours = contours.size();
    IndexedTStrip::TStrips tstrips_;

    for (std::size_t idx = 1; idx < total_contours; ++idx)
    {
        ContourTypePtr cur_contour = contours[idx];
        ContourTypePtr prev_contour = contours[idx - 1];

        ChristiansenTiling<RealType> triang(
                idx - 1, prev_contour->contour(), prev_contour->is_closed(),
                idx, cur_contour->contour(), cur_contour->is_closed());

        IndexedTStrip tstrip = triang.christiansen();
        tstrips_.push_back(tstrip);
    }

    // Union all tstrips into one mesh.
    IndexedTStrip joined_tstrips = IndexedTStrip::join(tstrips_);

    // Calculate the total number of vertices to create mesh efficiently.
    std::size_t total_vertices = 0;
    BOOST_FOREACH (const ContourTypePtr& c, contours)
    { total_vertices += c->contour()->size(); }

    // Populate the result mesh with vertices and create lookup tables.
    Mesh3D mesh(total_vertices);
    std::vector<std::size_t> lookup_table(total_contours);
    for (std::size_t idx = 0; idx < total_contours; ++idx)
    {
        std::size_t current_offset = mesh.add_vertices(*(contours[idx]->contour()));
        lookup_table[idx] = current_offset;
    }

    // Populate the result mesh with faces. Transform old <contour_id, vertex_ix>
    // keys into new <mesh_vertex_id> keys using lookup table.
    BOOST_FOREACH (const IndexedTStrip::Face& face, joined_tstrips.get_faces())
    {
        std::size_t new_a = lookup_table[face.A().first] + face.A().second;
        std::size_t new_b = lookup_table[face.B().first] + face.B().second;
        std::size_t new_c = lookup_table[face.C().first] + face.C().second;
        mesh.add_face(typename Mesh3D::Face(new_a, new_b, new_c));
    }

    return mesh;
}

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // TRIANGULATION_HPP_40CA4CA1_3752_4080_9BDE_D13393E21902
