
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
#include <utility>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "bo/mesh.hpp"
#include "bo/methods/distances_3d.hpp"
#include "bo/internal/surfaces/container_traversers.hpp"

namespace bo {
namespace methods {
namespace surfaces {

namespace detail {

// A class representing a triangle strip from two contours. Only vertex indices and
// contour ids are stored.
template <typename RealType>
class TStrip
{
public:
    typedef std::vector<TStrip> TStrips;

    typedef std::pair<std::size_t, std::size_t> ContourNodeID;
    typedef Triangle<ContourNodeID> Face;
    typedef std::vector<Face> Faces;

public:
    TStrip(std::size_t faces_count)
    {
        faces_.reserve(faces_count);
    }

    TStrip(std::size_t id1, std::size_t vertex_index1,
           std::size_t id2, std::size_t vertex_index2, std::size_t initial_count = 2)
        : id1_(id1), id2_(id2)
    {
        faces_.reserve(initial_count);
        current1_node_ = std::make_pair(id1_, vertex_index1);
        current2_node_ = std::make_pair(id2_, vertex_index2);
    }

    // Adds a new node and shifts span in the direction of the *first* contour.
    void add1(std::size_t vertex_idx)
    { current1_node_ = add_(id1_, vertex_idx); }

    // Adds a new node and shifts span in the direction of the *second* contour.
    void add2(std::size_t vertex_idx)
    { current2_node_ = add_(id2_, vertex_idx); }

    const Faces& get_faces() const
    { return faces_; }

    static TStrip join(const TStrips& tstrips)
    {
        std::size_t total_faces = 0;
        BOOST_FOREACH (const TStrip& tstrip, tstrips)
        { total_faces += tstrip.faces_.size(); }

        TStrip joined(total_faces);
        BOOST_FOREACH (const TStrip& tstrip, tstrips)
        { joined.faces_.insert(joined.faces_.end(), tstrip.faces_.begin(), tstrip.faces_.end()); }

        return joined;
    }

private:
    // Adds a new node and face based on the current span. Doesn't update current span.
    ContourNodeID add_(std::size_t contour_id, std::size_t vertex_idx)
    {
        ContourNodeID new_node = std::make_pair(contour_id, vertex_idx);
        faces_.push_back(Face(current1_node_, new_node, current2_node_));

        return new_node;
    }

private:
    std::size_t id1_;
    std::size_t id2_;

    Faces faces_;
    ContourNodeID current1_node_;
    ContourNodeID current2_node_;
};

} // namespace detail


// A class performing different triangulations for two slices.
template <typename RealType>
class Triangulation: public boost::noncopyable
{
public:
    typedef Triangulation<RealType> this_type;

    typedef Vector<RealType, 3> Point3D;
    typedef boost::function<RealType (Point3D, Point3D)> Metric;
    typedef std::vector<Point3D> Contour;
    typedef boost::shared_ptr<Contour> ContourPtr;
    typedef bo::Mesh<RealType> Mesh3D;

    typedef ContainerConstTraverser<Contour> ContourTraverser;
    typedef TraverseRuleFactory<Contour> TraverseFactory;
    typedef detail::TStrip<RealType> TStrip;

public:
    Triangulation(std::size_t id1, ContourPtr contour1, bool closed1,
                  std::size_t id2, ContourPtr contour2, bool closed2):
        id1_(id1), contour1_(contour1), closed1_(closed1),
        id2_(id2), contour2_(contour2), closed2_(closed2)
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
    TStrip christiansen()
    {
        // Check contours' length.
        if ((contour1_->size() < 2) && (contour2_->size() < 2))
            throw std::logic_error("Cannot run Christiansen triangulation for contours "
                                   "consisting of less than 2 vertices.");

        // If the first contour is closed but the second is opened, traverse direction
        // for the second contour may be determined incorrectly.
        if (closed1_ && !closed2_)
        {
            // Clone data before swapping in order not to spoil original contours.
            ContourPtr new_contour2(new Contour(*contour1_));
            ContourPtr new_contour1(new Contour(*contour2_));
            contour1_ = new_contour1;
            contour2_ = new_contour2;

            std::swap(closed1_, closed2_);
            std::swap(id1_, id2_);
        }

        // Create traversers for contours.
        ContourTraverser current1 = create_traverser1_();
        ContourTraverser candidate1 = current1 + 1;

        ContourTraverser current2 = create_traverser2_(current1);
        ContourTraverser candidate2 = current2 + 1;

        // Initialize output mesh.
        TStrip tstrip(id1_, current1.index(), id2_, current2.index(),
                      contour1_->size() + contour2_->size());

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

    template <typename ContourTypePtr>
    static Mesh3D christiansen(const ContourTypePtr& contour1, const ContourTypePtr& contour2)
    {
        TStrip tstrip = this_type(0, contour1->contour(), contour1->is_closed(),
                                  1, contour2->contour(), contour2->is_closed()).christiansen();

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
        BOOST_FOREACH (const typename TStrip::Face& face, tstrip.get_faces())
        {
            std::size_t new_a = lookup_table[face.A().first] + face.A().second;
            std::size_t new_b = lookup_table[face.B().first] + face.B().second;
            std::size_t new_c = lookup_table[face.C().first] + face.C().second;
            mesh.add_face(typename Mesh3D::Face(new_a, new_b, new_c));
        }

        return mesh;
    }

    template <typename ContourTypePtr>
    static Mesh3D christiansen(const std::vector<ContourTypePtr>& contours)
    {
        std::size_t total_contours = contours.size();
        std::vector<TStrip> tstrips_;

        for (std::size_t idx = 1; idx < total_contours; ++idx)
        {
            ContourTypePtr cur_contour = contours[idx];
            ContourTypePtr prev_contour = contours[idx - 1];

            this_type triang(idx - 1, prev_contour->contour(), prev_contour->is_closed(),
                             idx, cur_contour->contour(), cur_contour->is_closed());

            TStrip tstrip = triang.christiansen();
            tstrips_.push_back(tstrip);
        }

        // Union all tstrips into one mesh.
        TStrip joined_tstrips = TStrip::join(tstrips_);

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
        BOOST_FOREACH (const typename TStrip::Face& face, joined_tstrips.get_faces())
        {
            std::size_t new_a = lookup_table[face.A().first] + face.A().second;
            std::size_t new_b = lookup_table[face.B().first] + face.B().second;
            std::size_t new_c = lookup_table[face.C().first] + face.C().second;
            mesh.add_face(typename Mesh3D::Face(new_a, new_b, new_c));
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
            // with the contour1 direction.
            Point3D direction1 = *(traverser1 + 1) - *traverser1;
            Point3D direction2 = *(retvalue + 1) - *retvalue;
            if (direction1 * direction2 < 0)
                retvalue = ContourTraverser(TraverseFactory::Create(contour2_, c2min_idx, false));
        }
        else
        {
            // Check if contour2_ traverse direction should be swapped in order to
            // correspond with the contour1_'s direction.
            RealType dist_to_first = metric_(contour2_->front(), *traverser1);
            RealType dist_to_last = metric_(contour2_->back(), *traverser1);
            bool is_forward2 = (dist_to_first < dist_to_last) ? true : false;

            retvalue = ContourTraverser(TraverseFactory::Create(contour2_, is_forward2));
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
    std::size_t id1_;
    ContourPtr contour1_;
    bool closed1_;

    std::size_t id2_;
    ContourPtr contour2_;
    bool closed2_;

    Metric metric_;
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // TRIANGULATION_HPP_40CA4CA1_3752_4080_9BDE_D13393E21902
