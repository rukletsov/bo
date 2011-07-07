
#include "pch.h"

#include <boost/scope_exit.hpp>

#include "rply/rply.h"
#include "mesh.hpp"
#include "mesh_io.hpp"


namespace {

// Context for RPly callbacks. See comments on mesh_from_ply() for more info.
struct PLYContext
{
    common::Mesh* mesh_ptr;
    common::Mesh::Vertex* vertex;
    common::Mesh::Face* face;
};

// Callback for reading vertex component. Supposes that vertex has only 3 components
// and perform vertex adding after reading the third component.
int vertex_cb(p_ply_argument argument) {
    long type;
    PLYContext* context;
    ply_get_argument_user_data(argument, (void**)&context, &type);

    switch (type)
    {
    case 0:
        // Means we scan x-coord of a vertex.
        context->vertex->x() = ply_get_argument_value(argument);
        break;
    case 1:
        // Means we scan y-coord of a vertex.
        context->vertex->y() = ply_get_argument_value(argument);
        break;
    case 2:
        // Means we scan z-coord and are ready to store the vertex.
        context->vertex->z() = ply_get_argument_value(argument);
        context->mesh_ptr->add_vertex(*(context->vertex));
        break;
    default:
        break;
    }

    return 1;
}

// Callback for reading triangle face point. Works only with triangle faces and
// rejects other. Adds face after reading the third face point.
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
namespace io {

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
common::Mesh mesh_from_ply(const std::string& file_path)
{
    common::Mesh invalid_mesh(0);
    long nvertices, ntriangles;

    // Open .ply file for reading.
    p_ply ply = ply_open(file_path.c_str(), NULL);
    if (!ply)
        return invalid_mesh;

    // Try to read mesh data. Since the file has been already successfully opened
    // by this point, add a clean-up action by boost/scope_exit, which will be
    // executed at the end of the current scope.
    BOOST_SCOPE_EXIT ((&ply)) {
        ply_close(ply);
    } BOOST_SCOPE_EXIT_END

    if (!ply_read_header(ply))
        return invalid_mesh;

    // Prepare PLYContext and set callbacks for RPly reader.
    PLYContext context;
    nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, &context, 0);
    ply_set_read_cb(ply, "vertex", "y", vertex_cb, &context, 1);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, &context, 2);

    ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, &context, 0);

    // Create a mesh and fill the context with values.
    common::Mesh mesh(static_cast<std::size_t>(nvertices));
    common::Mesh::Vertex temp_vertex(0.0, 0.0, 0.0);
    common::Mesh::Face temp_face;
    context.mesh_ptr = &mesh;
    context.vertex = &temp_vertex;
    context.face = &temp_face;

    // Perform reading contents into mesh.
    if (!ply_read(ply))
        return invalid_mesh;

    return mesh;
}

// Writing to .ply files is rather straightforward. The only caveat is vertex type.
// Some shitty software doesn't support double type for vertices that's why a
// conversion to float is made, despite the fact that mesh stores double type.
bool mesh_to_ply(const common::Mesh& mesh, const std::string& file_path)
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
        if (!ply_add_element(oply, "vertex",
                             static_cast<long>(mesh.get_all_vertices().size())))
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
        if (!ply_add_element(oply, "face",
                             static_cast<long>(mesh.get_all_faces().size())))
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
        common::Mesh::Vertices::const_iterator vertices_end =
                mesh.get_all_vertices().end();
        for (common::Mesh::Vertices::const_iterator it = mesh.get_all_vertices().begin();
             it != vertices_end; ++it)
        {
            if (!ply_write(oply, it->x()))
                return false;
            if (!ply_write(oply, it->y()))
                return false;
            if (!ply_write(oply, it->z()))
                return false;
        }

        common::Mesh::Faces::const_iterator faces_end = mesh.get_all_faces().end();
        for (common::Mesh::Faces::const_iterator it = mesh.get_all_faces().begin();
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

} // namespace io
} // namespace common
