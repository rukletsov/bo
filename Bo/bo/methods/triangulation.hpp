
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
#include <boost/foreach.hpp>
#include <boost/noncopyable.hpp>

#include "bo/mesh.hpp"
#include "bo/internal/surfaces/indexed_tstrip.hpp"
#include "bo/internal/surfaces/christiansen_tiling.hpp"

namespace bo {
namespace methods {
namespace surfaces {

// This class contain various triangulation algorithms.
template <typename RealType>
class Triangulation
{
public:
    typedef bo::Mesh<RealType> RealMesh;
    typedef ChristiansenTiling<RealType> Christiansen;

public:
    // Performs modified Christiansen triangulation for two contours.
    template <typename ContourDescrPtr>
    static RealMesh christiansen(const ContourDescrPtr& contour1,
                                 const ContourDescrPtr& contour2)
    {

        IndexedTStrip tstrip = Christiansen(0, contour1->contour(), contour1->is_closed(),
                                            1, contour2->contour(), contour2->is_closed())
                .run();

        // Calculate the total number of vertices to create mesh efficiently.
        std::size_t total_vertices = contour1->contour()->size() + contour2->contour()->size();

        // Populate the result mesh with vertices and store vertices shift for the
        // second contour (it should be 0 for the first).
        RealMesh mesh(total_vertices);
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
            mesh.add_face(typename RealMesh::Face(new_a, new_b, new_c));
        }

        return mesh;
    }

    // Performs modified Christiansen triangulation for a collection of contours.
    // Note that the collection must be sorted, because the algorithm connects adjacent
    // contours.
    template <typename ContourDescrPtr>
    static RealMesh christiansen(const std::vector<ContourDescrPtr>& contours)
    {
        std::size_t total_contours = contours.size();
        IndexedTStrip::TStrips tstrips_;

        for (std::size_t idx = 1; idx < total_contours; ++idx)
        {
            ContourDescrPtr cur_contour = contours[idx];
            ContourDescrPtr prev_contour = contours[idx - 1];

            Christiansen triang(
                    idx - 1, prev_contour->contour(), prev_contour->is_closed(),
                    idx, cur_contour->contour(), cur_contour->is_closed());

            IndexedTStrip tstrip = triang.run();
            tstrips_.push_back(tstrip);
        }

        // Union all tstrips into one mesh.
        IndexedTStrip joined_tstrips = IndexedTStrip::join(tstrips_);

        // Calculate the total number of vertices to create mesh efficiently.
        std::size_t total_vertices = 0;
        BOOST_FOREACH (const ContourDescrPtr& c, contours)
        { total_vertices += c->contour()->size(); }

        // Populate the result mesh with vertices and create lookup tables.
        RealMesh mesh(total_vertices);
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
            mesh.add_face(typename RealMesh::Face(new_a, new_b, new_c));
        }

        return mesh;
    }
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#endif // TRIANGULATION_HPP_40CA4CA1_3752_4080_9BDE_D13393E21902
