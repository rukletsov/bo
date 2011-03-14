
#include "stdafx.h"

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
    
    context->face->operator [](value_index) = static_cast<int>(ply_get_argument_value(argument));
    if (value_index == 2)
        // Means we finished to read current face.
        context->mesh_ptr->add_face(*(context->face));

    return 1;
}

} // anonymous namespace


namespace common {

Mesh::Mesh(size_t initial_count): vertices(initial_count), faces(initial_count), 
    normals(initial_count), neighbours(initial_count), adjacent_faces(initial_count)
{ }

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
    p_ply ply = ply_open("input.ply", NULL);
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

} // namespace common
