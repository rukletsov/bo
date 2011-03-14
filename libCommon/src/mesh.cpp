
#include "stdafx.h"

#include <boost/format.hpp>

#include "rply.h"
#include "mesh.hpp"


namespace {

struct PLYContext
{
    common::Mesh* mesh_ptr;
    common::Mesh::Vertex* vertex;
    common::Mesh::Face* face;
};

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

int face_cb(p_ply_argument argument) {
    long length, value_index;
    PLYContext* context;
    ply_get_argument_property(argument, NULL, &length, &value_index);
    ply_get_argument_user_data(argument, (void**)&context, NULL);
    
    if (value_index != -1)
    {
        context->face->operator [](value_index) = static_cast<int>(ply_get_argument_value(argument));
        if (value_index == 2)
            // Means we finished to read current face.
            context->mesh_ptr->add_face(*(context->face));
    }
    // Check points quantity in the face.
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
    normals.reserve(initial_count); 
    neighbours.reserve(initial_count); 
    adjacent_faces.reserve(initial_count);
}

size_t Mesh::add_vertex(const Vertex& vertex)
{
    // Actually, a syncro primitive should be added here.
    vertices.push_back(vertex);
    size_t new_vertex_index = vertices.size() - 1;

    // Count vertex normal and store it.

    return new_vertex_index;
}

size_t Mesh::add_face(const Face& face)
{
    // Add a face and get its index. Syncro primitive needed.
    faces.push_back(face);
    size_t new_face_index = faces.size() - 1;

    // Update vertex neighbours.

    // Update adjacent faces.

    return new_face_index;
}

Mesh Mesh::from_ply(const std::string& file_path)
{
    Mesh invalid_mesh(0);

    long nvertices, ntriangles;
    p_ply ply = ply_open(file_path.c_str(), NULL);
    if (!ply) 
        return invalid_mesh;
    if (!ply_read_header(ply)) 
        return invalid_mesh;

    // Prepare PLYContext and callbacks for RPly reader.
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

    // Read contents into mesh.
    if (!ply_read(ply)) 
        return invalid_mesh;

    ply_close(ply);

    return mesh;
}

std::ostream& operator <<(std::ostream &os, const common::Mesh& obj)
{
    // Add syncro primitives to stream operator.
    os << boost::format("Mesh object %1$#x, %2% bytes: ") % &obj % sizeof(obj) 
        << std::endl << "Vertices: " << obj.vertices.size() << std::endl;
    
    
    common::Mesh::Vertices::const_iterator vertices_end = obj.vertices.end();
    for (common::Mesh::Vertices::const_iterator it = obj.vertices.begin();
         it != vertices_end; ++it)
    {
        os << boost::format("x: %1%, %|18t|y: %2%, %|36t|z: %3%,") % it->x % it->y % it->z
            << std::endl;
    }

    os << "Faces: " << obj.faces.size() << std::endl 
        << boost::format("end of Mesh object %1$#x.") % &obj << std::endl;

    return os;
}

} // namespace common
