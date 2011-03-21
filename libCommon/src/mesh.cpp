
#include "stdafx.h"

#include <boost/format.hpp>
#include <boost/scope_exit.hpp>

#include "rply.h"
#include "mesh.hpp"


namespace {

// Context for RPly callbacks. See comments on Mesh::from_ply() for more info.
struct PLYContext
{
    common::Mesh* mesh_ptr;
    common::Mesh::Vertex* vertex;
    common::Mesh::Face* face;
};

// Callback for reading vertex component. Suppose that vertex has only 3 components
// and perform vertex adding after reading the third component.
int vertex_cb(p_ply_argument argument) {
    long type;
    PLYContext* context;
    ply_get_argument_user_data(argument, (void**)&context, &type);

    switch (type)
    {
    case 0:
        // Means we scan x-coord of a vertex.
        context->vertex->x = ply_get_argument_value(argument);
        break;
    case 1:
        // Means we scan y-coord of a vertex.
        context->vertex->y = ply_get_argument_value(argument);
        break;
    case 2:
        // Means we scan z-coord and are ready to store the vertex.
        context->vertex->z = ply_get_argument_value(argument);
        context->mesh_ptr->add_vertex(*(context->vertex));
        break;
    default:
        break;
    }

    return 1;
}

// Callback for reading triangle face point. Works only with triangle faces and
// rejects other. Add face after reading the third face point.
int face_cb(p_ply_argument argument) {
    long length, value_index;
    PLYContext* context;
    ply_get_argument_property(argument, NULL, &length, &value_index);
    ply_get_argument_user_data(argument, (void**)&context, NULL);
    
    // Means not the first component (length) of the list is being read.
    if (value_index != -1)
    {
        context->face->operator [](static_cast<std::size_t>(value_index)) = 
            static_cast<std::size_t>(ply_get_argument_value(argument));
        if (value_index == 2)
            // Means we finished to read current face.
            context->mesh_ptr->add_face(*(context->face));
    }
    // Check points quantity in the face when reading the first list component.
    else if (length != 3)
        return 0;

    return 1;
}

} // anonymous namespace


namespace common {

Mesh::Mesh(size_t initial_count)
{
    vertices.reserve(initial_count);
    faces.reserve(initial_count);
    face_normals.reserve(initial_count); 
    neighbours.reserve(initial_count); 
    adjacent_faces.reserve(initial_count);
}

size_t Mesh::add_vertex(const Vertex& vertex)
{
    // Actually, a syncro primitive should be added here.
    vertices.push_back(vertex);
    size_t new_vertex_index = vertices.size() - 1;

    // Add empty connectivity containers.
    neighbours.push_back(AdjacentVerticesPerVertex());
    adjacent_faces.push_back(AdjacentFacesPerVertex());

    return new_vertex_index;
}

size_t Mesh::add_face(const Face& face)
{
    // Add a face and get its index. Syncro primitive needed.
    faces.push_back(face);
    size_t new_face_index = faces.size() - 1;

    // Compute and save face normal.
    face_normals.push_back(compute_face_normal(face).normalized());

    // Update vertex neighbours.
    add_neighbouring_pair(face.A(), face.B());
    add_neighbouring_pair(face.B(), face.C());
    add_neighbouring_pair(face.C(), face.A());

    // Update adjacent faces.
    add_adjacent_face(face.A(), new_face_index);
    add_adjacent_face(face.B(), new_face_index);
    add_adjacent_face(face.C(), new_face_index);

    return new_face_index;
}

Mesh::Normal Mesh::get_vertex_normal(size_t vertex_index) const
{
    // A normal of a vertex is a sum of weighted normals of adjacent faces.
    Normal retvalue;

    AdjacentFacesPerVertex::const_iterator faces_end = 
        adjacent_faces[vertex_index].end();
    for (AdjacentFacesPerVertex::const_iterator face_index = 
        adjacent_faces[vertex_index].begin(); face_index != faces_end; ++face_index)
    {
        // Determine to which face vertex considered vertex belong. Without loss of
        // generality, suppose, that face[2] == vertex_index.
        const Face& face = faces[*face_index];
        size_t pt1 = face[0];
        size_t pt2 = face[1];
        if (face[0] == vertex_index)
            pt1 = face[2];
        else if (face[1] == vertex_index)
            pt2 = face[2];

        // Compute face's normal weight (multiplied dot products of two edges).
        Vector3<double> edge1 = vertices[pt1] - vertices[vertex_index];
        Vector3<double> edge2 = vertices[pt2] - vertices[vertex_index];
        double weight = ((edge1 * edge1) * (edge2 * edge2));

        // Append weighted face's normal.
        retvalue = retvalue + face_normals[*face_index] * (1.0 / weight);
    }

    return retvalue;
}


// RPly library is used for reading and writing meshes to .ply files.
// The library is written in C and its sources are included in the project.
// 
// Reading is done via callbacks. RPly first reads .ply header and then reads
// data with known structure and calls defined function for every data unit.
// In the simplified case which is supposed in our case, data is vertices and
// a list of triangle faces. So, two callbacks are defined above in an anonymous 
// namespace. First callback is for reading vertex components and the second is 
// for reading face vertices. Both functions use the same context for storing
// temporary values and for accessing mesh function.
Mesh Mesh::from_ply(const std::string& file_path)
{
    Mesh invalid_mesh(0);
    long nvertices, ntriangles;

    p_ply ply = ply_open(file_path.c_str(), NULL);
    if (!ply) 
        return invalid_mesh;
    if (!ply_read_header(ply)) 
        return invalid_mesh;

    // Prepare PLYContext and set callbacks for RPly reader.
    PLYContext context;
    nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, &context, 0);
    ply_set_read_cb(ply, "vertex", "y", vertex_cb, &context, 1);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, &context, 2);

    ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, &context, 0);

    // Create a mesh and fill the context with values.
    Mesh mesh(static_cast<std::size_t>(nvertices));
    Mesh::Vertex temp_vertex(0.0, 0.0, 0.0);
    Mesh::Face temp_face;
    context.mesh_ptr = &mesh;
    context.vertex = &temp_vertex;
    context.face = &temp_face;

    // Perform reading contents into mesh.
    if (!ply_read(ply)) 
        return invalid_mesh;

    ply_close(ply);

    return mesh;
}

// Writing to .ply files is rather straightforward. The only caveat is vertex type.
// Some shitty software doesn't support double type for vertices that's why a 
// conversion to float is made, despite the fact that mesh stores double type.
bool Mesh::to_ply(const std::string& file_path) const
{
    // Create .ply file in ascii format.
    p_ply oply = ply_create(file_path.c_str(), PLY_ASCII, NULL);
    if (!oply) 
        return false;

    // Suppose successful file close operation unless otherwise specified.
    bool close_succeeded = true;

    // Try to write mesh data. If an error happens, stop further writing and return,
    // clean-up will be made automatically by boost/scope_exit.
    {
        // Since the file has been opened successfully by this point, add a clean-up
        // action, which will be executed at the end of the current scope.
        BOOST_SCOPE_EXIT ((&close_succeeded) (&oply)) {
            if (!ply_close(oply))  
                close_succeeded = false;
        } BOOST_SCOPE_EXIT_END

        // Add "vertex" element.
        if (!ply_add_element(oply, "vertex", this->vertices.size())) 
            return false;
        // Add "vertex" properties. if the type parameter is not PLY_LIST, two last
        // parameters are ignored. So, PLY_LIST is passed for two last parameters
        // as the most unsuitable value.
        if (!ply_add_property(oply, "x", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
            return false;
        if (!ply_add_property(oply, "y", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
            return false;
        if (!ply_add_property(oply, "z", PLY_FLOAT, PLY_LIST, PLY_LIST)) 
            return false;

        // Add "face" element.
        if (!ply_add_element(oply, "face", this->faces.size())) 
            return false;
        // Add "face" only property. It is a list of vertex indices. 
        if (!ply_add_property(oply, "vertex_indices", PLY_LIST, PLY_UCHAR, PLY_UINT)) 
            return false;

        // Add a comment and an obj_info.
        if (!ply_add_comment(oply, "libCommon generated PLY file")) 
            return false;
        if (!ply_add_obj_info(oply, "common::Mesh class dump")) 
            return false;

        // Write .ply header.
        if (!ply_write_header(oply)) 
            return false;

        // Write mesh data in the same order as declared above.
        Mesh::Vertices::const_iterator vertices_end = this->vertices.end();
        for (Mesh::Vertices::const_iterator it = this->vertices.begin();
             it != vertices_end; ++it)
        {
            if (!ply_write(oply, it->x)) 
                return false;
            if (!ply_write(oply, it->y))
                return false;
            if (!ply_write(oply, it->z)) 
                return false;
        }

        Mesh::Faces::const_iterator faces_end = this->faces.end();
        for (Mesh::Faces::const_iterator it = this->faces.begin();
            it != faces_end; ++it)
        {
            // 3 can be hardcoded since Mesh works only with triangle faces.
            if (!ply_write(oply, 3)) 
                return false;
            if (!ply_write(oply, static_cast<double>(it->A()))) 
                return false;
            if (!ply_write(oply, static_cast<double>(it->B()))) 
                return false;
            if (!ply_write(oply, static_cast<double>(it->C()))) 
                return false;
        }
    } // end of scope with file writing.
    
    // Returning false is possible only when an error occured during file closing.
    // All other errors lead to immediate function termination, i.e. not here.
    return close_succeeded;
}

// Print formatted mesh data to a given stream. See boost.format library for more
// details about formatting.
std::ostream& operator <<(std::ostream &os, const Mesh& obj)
{
    // TODO: Add syncro primitives to stream operator and, perhaps, verbose levels.
    
    // Print full vertex info.
    os << boost::format("Mesh object %1$#x, %2% bytes: ") % &obj % sizeof(obj) 
       << std::endl << "Vertices: " << obj.vertices.size() << std::endl;

    Mesh::Vertices::const_iterator vertices_end = obj.vertices.end();
    for (Mesh::Vertices::const_iterator it = obj.vertices.begin();
         it != vertices_end; ++it)
    {
        // Print vertex coordinates.
        size_t index = it - obj.vertices.begin();
        os << boost::format("vertex %1%: ") % index << std::endl << "\t"
           << boost::format("x: %1%, %|18t|y: %2%, %|36t|z: %3%,") % it->x % it->y % it->z;

        // Print neighbouring vertices.
        os << std::endl << "\t"
           << boost::format("neighbours (%1%): ") % obj.neighbours[index].size();

        Mesh::AdjacentVerticesPerVertex::const_iterator neighbours_end = 
            obj.neighbours[index].end();
        for (Mesh::AdjacentVerticesPerVertex::const_iterator neighbour = 
            obj.neighbours[index].begin(); neighbour != neighbours_end; ++neighbour)
        {
            os << boost::format("%1%, %|4t|") % (*neighbour);
        }

        // Print adjacent faces.            
        os << std::endl << "\t" 
           << boost::format("adjacent faces (%1%): ") % obj.adjacent_faces[index].size();

        Mesh::AdjacentFacesPerVertex::const_iterator faces_end = 
            obj.adjacent_faces[index].end();
        for (Mesh::AdjacentFacesPerVertex::const_iterator face = 
            obj.adjacent_faces[index].begin(); face != faces_end; ++face)
        {
            os << boost::format("%1%, %|4t|") % (*face);
        }

        os << std::endl;
    }

    // Print full faces info.
    os << "Faces: " << obj.faces.size() << std::endl;

    Mesh::Faces::const_iterator faces_end = obj.faces.end();
    for (Mesh::Faces::const_iterator face = obj.faces.begin(); 
        face != faces_end; ++face)
    {
        // Print face's vertices.
        size_t index = face - obj.faces.begin();
        os << boost::format("face %1%: ") % index << std::endl << "\t"
           << boost::format("A: %1%, %|18t|B: %2%, %|36t|C: %3%,") 
           % face->A() % face->B() % face->C();

        // Print face's normal.
        Mesh::Normal normal = obj.face_normals[index];
        os << std::endl << "\t"
           << boost::format("normal: (%1%, %2%, %3%)") % normal.x % normal.y % normal.z;

        os << std::endl;        
    }

    // Print footer and return.
    os << boost::format("end of Mesh object %1$#x.") % &obj << std::endl;

    return os;
}

// Private Urility functions.
bool Mesh::add_neighbouring_pair(size_t vertex1, size_t vertex2)
{
    // If the neighbouring relation between given vertices already exists, 
    // set::insert signal this and won't add a duplicate. A neighbouring relation
    // should be mutual. In case one vertex has another as a neighbour and another
    // has not, report error through assertion.
    bool exist1 = neighbours[vertex1].insert(vertex2).second;
    bool exist2 = neighbours[vertex2].insert(vertex1).second;

    BOOST_ASSERT(!(exist1 ^ exist2));

    // Since relation is mutual, either exist1 or exist2 can be returned.
    return exist1;
}

bool Mesh::add_adjacent_face(size_t vertex, size_t face)
{
    // If the face is already attached to the vertex, set::insert won't add a 
    // duplicate. See set::insert documentation.
    return 
        adjacent_faces[vertex].insert(face).second;
}

Mesh::Normal Mesh::compute_face_normal(const Face& face) const
{
    Vector3<double> a = vertices[face.A()] - vertices[face.B()];
    Vector3<double> b = vertices[face.B()] - vertices[face.C()];
    
    return 
        a.cross_product(b);
}

} // namespace common
