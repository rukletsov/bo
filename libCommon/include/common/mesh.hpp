
/******************************************************************************

    mesh.hpp, v 1.0.4 2011.03.14

    Triangular mesh declaration. RPly library is used for IO. 

    Copyright (c) 2010, 2011
    Alexander Rukletsov <alexander.rukletsov@ziti.uni-heidelberg.de>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    1.	Redistributions of source code must retain the above copyright
	    notice, this list of conditions and the following disclaimer.
    2.	Redistributions in binary form must reproduce the above copyright
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

#ifndef MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
#define MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_

#include <cstddef>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <boost/array.hpp>

#include "vector.hpp"
#include "triangle.hpp"


namespace common {

// A basic class for a 3D triangular mesh. Consumes more memory than possible
// minimum (a standard graph storage) but provides faster access to frequently
// used structures and operations. NOT thread-safe in the current implemetation.
class Mesh
{

public:
    // Mesh vertices. Their order shouldn't be changed, since other collections 
    // use vertex indices as references.
    typedef Vector<double, 3> Vertex;
    typedef std::vector<Vertex> Vertices;

    // Face is a triangle with vertices representing indices of the mesh vertices.
    typedef Triangle<std::size_t> Face;
    typedef std::vector<Face> Faces;

    // Normal is a 3-vector can be attached to a vertex or to a face. A normal 
    // corresponds to a component with the same index in relevant collection.
    typedef Vector<double, 3> Normal;
    typedef std::vector<Normal> Normals;

    // Neighbours for each vertex. Each neighbour contains indices of vertex in
    // vertices collection.
    typedef std::set<std::size_t> AdjacentVerticesPerVertex;
    typedef std::vector<AdjacentVerticesPerVertex> AdjacentVertices;

    // Neighbouring faces for each vertex. Each adjacent face is an index of a face
    // in faces collection.
    typedef std::set<std::size_t> AdjacentFacesPerVertex;
    typedef std::vector<AdjacentFacesPerVertex> AdjacentFaces;

public:
    // Create an empty mesh ready to store initial_count vertices.
    Mesh(std::size_t initial_count);

    // Add a new vertex to the mesh and return its index.
    std::size_t add_vertex(const Vertex& vertex);

    // Add a new face and return its index. Update dependent collections.
    std::size_t add_face(const Face& face);

    Normal get_face_normal(std::size_t face_index) const;

    // Get vertex normal, computed according to 
    // N.Max, "Weights for Computing Vertex Normals from Facet Normals",
    // Journal of Graphics Tools, Vol. 4, No. 2, 1999.
    Normal get_vertex_normal(std::size_t vertex_index) const;

    // Temporary accessor methods.
    const Vertices& get_all_vertices() const;
    const Faces& get_all_faces() const;
    const Normals& get_all_face_normals() const;

    // IO functions, allow to read mesh from and write to a .ply files, 
    // print formatted mesh data to an std::ostream.
    static Mesh from_ply(const std::string& file_path);
    bool to_ply(const std::string& file_path) const;
    friend std::ostream& operator <<(std::ostream& os, const Mesh& obj);

private:
    // Add connectivity relations. Return false in case of new relation leads to 
    // a duplicate.
    bool add_edge_(size_t vertex1, size_t vertex2);
    bool add_adjacent_face_(size_t vertex, size_t face);

    Normal compute_face_normal_(const Face& face) const;

private:
    // Basic mesh data.
    Vertices vertices_;
    Faces faces_;
    Normals face_normals_;
    // Some other properties can be used, e.g. triangle strips (for speeding up
    // rendering), curvature information, BBox, grid, etc (See TriMesh implementation
    // by Szymon Rusinkiewicz as an example.)

    // Connectivity structures.
    AdjacentVertices neighbours_;
    AdjacentFaces adjacent_faces_;
    // We can also add, e.g. faces adjacent to faces over edges, i.e. each face will
    // have maximum 3 these neighbouring faces.
};


inline
const Mesh::Vertices& Mesh::get_all_vertices() const
{
    return vertices_;
}

inline
const Mesh::Faces& Mesh::get_all_faces() const
{
    return faces_;
}

inline
const Mesh::Normals& Mesh::get_all_face_normals() const
{
    return face_normals_;
}

} // namespace common

#endif // MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
