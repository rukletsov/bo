/******************************************************************************

  Basic N-dimensional topology.

  Copyright (c) 2013
  Dzmitry Hlindzich <hlindzich@gmail.com>
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

#ifndef TOPOLOGY_HPP_79A57FE0_9EFD_11E2_9E96_0800200C9A66_
#define TOPOLOGY_HPP_79A57FE0_9EFD_11E2_9E96_0800200C9A66_

#include <vector>
#include <utility>
#include <set>

#include <boost/static_assert.hpp>

#include "bo/vector.hpp"

namespace bo{
namespace topology{

// Basic N-dimensional hyperrectangle (N-orthotope) geometry.
template <typename RealType, std::size_t N>
class OrthotopeTopology
{
private:
    BOOST_STATIC_ASSERT(N > 0);

public:
    typedef Vector<RealType, N> Point;
    typedef std::pair<Point, Point> Edge;
    typedef std::vector<Edge> Edges;

    // Returns all the edges of a unit N-orthotope.
    static Edges edges()
    {
        // Initialize the source of the edge traverse.
        PointSet in_vertices;
        in_vertices.insert(Point(0));

        // Traverse the edges.
        Edges e;
        traverse_edges(in_vertices, e);

        return e;
    }

    // Returns the number of edges in a N-orthotope.
    static std::size_t edge_count()
    {
        return (1 << (N - 1)) * N;
    }

    // Returns the number of vertices in a N-orthotope.
    static std::size_t vertex_count()
    {
        return 1 << N;
    }

    // Returns the number of facets ((N-1)-dimensional faces) in a N-othotope.
    static std::size_t facet_count()
    {
        return 2 * N;
    }

private:

    struct PointCompare
    {
        bool operator()(const Point &p1, const Point &p2) const
        {
            for (std::size_t i = 0; i < N; ++i)
            {
                if (p1[i] < p2[i])
                    return true;
                if (p1[i] > p2[i])
                    return false;
            }

            return false;
        }
    };

    typedef std::set<Point, PointCompare> PointSet;

    // Recursively traverses the edges of a N-orthotope starting from the given set
    // of vertices in ascending order and pushes the edges into the given container.
    static void traverse_edges(const PointSet &in_vertices, Edges &edges)
    {
        if (in_vertices.size() == 1 && *in_vertices.begin() == Point(1))
            return;

        PointSet out_vertices;

        for (typename PointSet::const_iterator it = in_vertices.begin(); it != in_vertices.end(); ++it)
        {
            Point edge_begin = *it;

            // Find all adjacent vertices for the current one.
            for (std::size_t k = 0; k < N; ++k)
            {
                if (edge_begin[k] == 1)
                    continue;

                // An adjacent vertex differs only in one dimension.
                Point edge_end = edge_begin;
                edge_end[k] = 1;

                // Insert the edge into the collection.
                edges.push_back(Edge(edge_begin, edge_end));

                // Save the end vertex.
                out_vertices.insert(edge_end);
            }
        }

        traverse_edges(out_vertices, edges);
    }

};

} // namespace topology
} // namespace bo

#endif // TOPOLOGY_HPP_79A57FE0_9EFD_11E2_9E96_0800200C9A66_
