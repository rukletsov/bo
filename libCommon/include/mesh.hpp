
#ifndef MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
#define MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_

#include <cstddef>
#include <vector>
#include <boost/array.hpp>

#include "point.hpp"
#include "triangle.hpp"


namespace common {

class Mesh
{

public:
    // Mesh vertices. Their order shouldn't be changed, since other collections 
    // use vertex indices as references.
    typedef Point3<double> Vertex;
    typedef std::vector<Vertex> Vertices;

    // Face is a triangle with vertices representing indices of the mesh vertices.
    typedef Triangle<std::size_t> Face;
    typedef std::vector<Face> Faces;

    // Normal is a 3-vector and is attached to a vertex. A normal corresponds to a 
    // vertex with the same index in vertices collection.
    typedef boost::array<double, 3> Normal;
    typedef std::vector<Normal> Normals;

    // Neighbours for each vertex. Each neighbour contains indices of vertex in
    // vertices collection.
    typedef std::vector<std::size_t> AdjacentVertex;
    typedef std::vector<AdjacentVertex> AdjacentVertices;

    // Neighbouring faces for each vertex. Each adjacent face is an index of a face
    // in faces collection.
    typedef std::vector<std::size_t> AdjacentFace;
    typedef std::vector<AdjacentFace> AdjacentFaces;

private:
    // Basic mesh data.
    Vertex vertices;
    Faces faces;
    Normals normals;
    // Some other properties can be used, e.g. triangle strips (for speeding up
    // rendering), curvature information, BBox, grid, etc (See TriMesh implementation
    // by Szymon Rusinkiewicz as an example.)

	// Connectivity structures.
	AdjacentVertices neighbours;
	AdjacentFaces adjacent_faces;
    // We can also add, e.g. faces adjacent to faces over edges, i.e. each face will
    // have maximum 3 these neighbouring faces.
};

} // namespace common

#endif // MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
