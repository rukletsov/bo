#ifndef MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
#define MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_

#include <vector>
//#include <boost/cstdint.hpp>

#include "point.hpp"


namespace common {

class Mesh
{

private:
    std::vector<Point3<double> > vertices;
    //std::vector<Face> faces;
    //std::vector<vec> normals;
    // Some other properties can be used, e.g. triangle strips (for speeding up
    // rendering), curvature information, BBox, grid, etc (See TriMesh implementation
    // by Szymon Rusinkiewicz as an example.)

public:
};

} // namespace common

#endif // MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
