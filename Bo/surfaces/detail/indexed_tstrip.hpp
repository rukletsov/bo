
/******************************************************************************

  TStrip representation from two contours. Contains only indices to vertices
  in original contours.

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

#ifndef INDEXED_TSTRIP_HPP_E720F05E_DF26_11E2_BD39_4057033B4CBB
#define INDEXED_TSTRIP_HPP_E720F05E_DF26_11E2_BD39_4057033B4CBB

#include <vector>
#include <utility>
#include <boost/foreach.hpp>

#include "bo/core/triangle.hpp"

namespace bo {
namespace surfaces {
namespace detail {

// A class representing a triangle strip from two contours. Only vertex indices and
// contour ids are stored.
class IndexedTStrip
{
public:
    typedef std::vector<IndexedTStrip> TStrips;

    typedef std::pair<std::size_t, std::size_t> ContourNodeID;
    typedef Triangle<ContourNodeID> Face;
    typedef std::vector<Face> Faces;

public:
    IndexedTStrip(std::size_t faces_count)
    {
        faces_.reserve(faces_count);
    }

    IndexedTStrip(std::size_t id1, std::size_t vertex_index1,
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

    static IndexedTStrip join(const TStrips& tstrips)
    {
        std::size_t total_faces = 0;
        BOOST_FOREACH (const IndexedTStrip& tstrip, tstrips)
        { total_faces += tstrip.faces_.size(); }

        IndexedTStrip joined(total_faces);
        BOOST_FOREACH (const IndexedTStrip& tstrip, tstrips)
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
    ContourNodeID current1_node_;
    std::size_t id2_;
    ContourNodeID current2_node_;

    Faces faces_;
};

} // namespace detail
} // namespace surfaces
} // namespace bo

#endif // INDEXED_TSTRIP_HPP_E720F05E_DF26_11E2_BD39_4057033B4CBB
