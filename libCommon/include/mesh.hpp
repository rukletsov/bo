
#ifndef MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
#define MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_

#include <vector>

#include "point.hpp"
#include "triangle.hpp"


namespace common {

class Mesh
{

public:
    typedef Point3<double> Vertex;
    typedef std::vector<Vertex> Vertices;

    // Face is a triangle with vertices representing indices of the mesh vertices.
    typedef Triangle3<int> Face;
    typedef std::vector<Face> Faces;

private:
    Vertex vertices;
    Faces faces;
    //std::vector<vec> normals;
    // Some other properties can be used, e.g. triangle strips (for speeding up
    // rendering), curvature information, BBox, grid, etc (See TriMesh implementation
    // by Szymon Rusinkiewicz as an example.)
};

} // namespace common

#endif // MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
