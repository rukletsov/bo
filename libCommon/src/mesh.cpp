
#include "pch.h"

#include <stdexcept>
#include <limits>
#include <boost/assert.hpp>
#include <boost/format.hpp>

#include "mesh.hpp"
#include "3d_distances.hpp"


namespace common {

Mesh::Mesh(std::size_t initial_count)
{
    vertices_.reserve(initial_count);
    faces_.reserve(initial_count);
    face_normals_.reserve(initial_count);
    neighbours_.reserve(initial_count);
    adjacent_faces_.reserve(initial_count);
}

std::size_t Mesh::add_vertex(const Vertex& vertex)
{
    // Actually, a syncro primitive should be added here.
    vertices_.push_back(vertex);
    std::size_t new_vertex_index = vertices_.size() - 1;

    // Add empty connectivity containers.
    neighbours_.push_back(AdjacentVerticesPerVertex());
    adjacent_faces_.push_back(AdjacentFacesPerVertex());

    // Check for post-conditions. These include sizes of connectivity structures
    // and vertex container. If any of the condition is not satisfied consider
    // this as an internal bug. Therefore no need to throw.
    BOOST_ASSERT(((neighbours_.size() == vertices_.size()) ||
                  (adjacent_faces_.size() == vertices_.size())) &&
                 "Vertex connectivity structures are of different sizes.");

    return new_vertex_index;
}

std::size_t Mesh::add_face(const Face& face)
{
    // Check if all vertices exist in the mesh.
    if ((vertices_.size() <= face.A()) ||
        (vertices_.size() <= face.B()) ||
        (vertices_.size() <= face.C()))
    {
        BOOST_ASSERT(false && "Bad vertex indices in the face.");
        throw std::out_of_range("Face cannot be added to the mesh because it "
                                "references non-existent vertices.");
    }

    // Add the face and get its index. Syncro primitive needed.
    faces_.push_back(face);
    std::size_t new_face_index = faces_.size() - 1;

    // Update vertex neighbours (insert edges).
    add_edge_(face.A(), face.B());
    add_edge_(face.B(), face.C());
    add_edge_(face.C(), face.A());

    // Update adjacent faces.
    add_adjacent_face_(face.A(), new_face_index);
    add_adjacent_face_(face.B(), new_face_index);
    add_adjacent_face_(face.C(), new_face_index);

    // Compute and save face normal.
    face_normals_.push_back(compute_face_normal_(face));

    // Check for post-conditions. These include sizes of connectivity structures
    // and face container. If any of the condition is not satisfied consider this
    // as an internal bug. Therefore no need to throw.
    BOOST_ASSERT((face_normals_.size() == faces_.size()) &&
                 "Vertex connectivity structures are of different sizes.");

    return new_face_index;
}

Mesh::Normal Mesh::get_face_normal(std::size_t face_index) const
{
    // Check if the given face exists in the mesh.
    face_rangecheck(face_index);

    return face_normals_[face_index];
}

Mesh::Normal Mesh::get_vertex_normal(std::size_t vertex_index) const
{
    // TODO: add caching for computed normals.

    // Check if the given vertex exists in the mesh.
    vertex_rangecheck(vertex_index);

    // A normal of a vertex is a sum of weighted normals of adjacent faces.
    Normal normal;

    AdjacentFacesPerVertex::const_iterator faces_end =
        adjacent_faces_[vertex_index].end();
    for (AdjacentFacesPerVertex::const_iterator face_index =
        adjacent_faces_[vertex_index].begin(); face_index != faces_end; ++face_index)
    {
        // Determine to which face vertex considered vertex belong. Without loss of
        // generality, suppose, that face[2] == vertex_index.
        const Face& face = faces_[*face_index];
        std::size_t pt1 = face[0];
        std::size_t pt2 = face[1];
        if (face[0] == vertex_index)
            pt1 = face[2];
        else if (face[1] == vertex_index)
            pt2 = face[2];

        // Compute face's normal weight (multiplied dot products of two edges).
        Vector<double, 3> edge1 = vertices_[pt1] - vertices_[vertex_index];
        Vector<double, 3> edge2 = vertices_[pt2] - vertices_[vertex_index];
        double weight = ((edge1 * edge1) * (edge2 * edge2));

        // Append weighted face's normal.
        normal += face_normals_[*face_index] * (1.0 / weight);
    }

    // Face normals can be the null vectors.
    try {
        normal = normal.normalized();
    }
    catch(const std::logic_error&)
    { }

    return normal;
}

Mesh::Vertex Mesh::get_closest_point(const Vertex& point) const
{
    Vertex closest;
    double min_distance = std::numeric_limits<double>::max();

    std::size_t faces_size = faces_.size();
    for (std::size_t face_index = 0; face_index < faces_size; ++face_index)
    {
        Vertex face_closest = closest_point_on_face(face_index, point);
        double face_distance = (face_closest - point).eucl_norm();
        if (face_distance < min_distance)
        {
            min_distance = face_distance;
            closest = face_closest;
        }
    }

    return closest;
}

double Mesh::distance(const Vertex& point) const
{
    return
        methods::euclidean_distance(get_closest_point(point), point);
}

const Mesh::AdjacentVerticesPerVertex& Mesh::get_neighbouring_vertices(
    std::size_t vertex_index) const
{
    // Check if the given vertex exists in the mesh.
    vertex_rangecheck(vertex_index);

    return neighbours_[vertex_index];
}

const Mesh::AdjacentFacesPerVertex& Mesh::get_neighbouring_faces_by_vertex(
    std::size_t vertex_index) const
{
    // Check if the given vertex exists in the mesh.
    vertex_rangecheck(vertex_index);

    return adjacent_faces_[vertex_index];
}

// Print formatted mesh data to a given stream. See boost.format library for more
// details about formatting.
std::ostream& operator <<(std::ostream& os, const Mesh& obj)
{
    // TODO: Add syncro primitives to stream operator and, perhaps, verbose levels.

    // Print full vertex info.
    os << boost::format("Mesh object %1$#x, %2% bytes: ") % &obj % sizeof(obj)
       << std::endl << "Vertices: " << obj.vertices_.size() << std::endl;

    Mesh::Vertices::const_iterator vertices_end = obj.vertices_.end();
    for (Mesh::Vertices::const_iterator it = obj.vertices_.begin();
         it != vertices_end; ++it)
    {
        // Print vertex coordinates.
        std::size_t index = it - obj.vertices_.begin();
        os << boost::format("vertex %1%: ") % index << std::endl << "\t"
           << boost::format("x: %1%, %|18t|y: %2%, %|36t|z: %3%,") % it->x()
           % it->y() % it->z();

        // Print neighbouring vertices.
        os << std::endl << "\t"
           << boost::format("neighbours (%1%): ") % obj.neighbours_[index].size();

        Mesh::AdjacentVerticesPerVertex::const_iterator neighbours_end =
            obj.neighbours_[index].end();
        for (Mesh::AdjacentVerticesPerVertex::const_iterator neighbour =
            obj.neighbours_[index].begin(); neighbour != neighbours_end; ++neighbour)
        {
            os << boost::format("%1%, %|4t|") % (*neighbour);
        }

        // Print adjacent faces.
        os << std::endl << "\t" << boost::format("adjacent faces (%1%): ")
           % obj.adjacent_faces_[index].size();

        Mesh::AdjacentFacesPerVertex::const_iterator faces_end =
            obj.adjacent_faces_[index].end();
        for (Mesh::AdjacentFacesPerVertex::const_iterator face =
            obj.adjacent_faces_[index].begin(); face != faces_end; ++face)
        {
            os << boost::format("%1%, %|4t|") % (*face);
        }

        os << std::endl;
    }

    // Print full faces info.
    os << "Faces: " << obj.faces_.size() << std::endl;

    Mesh::Faces::const_iterator faces_end = obj.faces_.end();
    for (Mesh::Faces::const_iterator face = obj.faces_.begin();
        face != faces_end; ++face)
    {
        // Print face's vertices.
        std::size_t index = face - obj.faces_.begin();
        os << boost::format("face %1%: ") % index << std::endl << "\t"
           << boost::format("A: %1%, %|18t|B: %2%, %|36t|C: %3%,")
           % face->A() % face->B() % face->C();

        // Print face's normal.
        Mesh::Normal normal = obj.face_normals_[index];
        os << std::endl << "\t"
           << boost::format("normal: (%1%, %2%, %3%)") % normal.x() % normal.y()
           % normal.z();

        os << std::endl;
    }

    // Print footer and return.
    os << boost::format("end of object %1$#x.") % &obj << std::endl;

    return os;
}

// Private utility functions.
bool Mesh::add_edge_(std::size_t vertex1, std::size_t vertex2)
{
    // If the neighbouring relation between given vertices already exists,
    // set::insert signal this and won't add a duplicate. A neighbouring relation
    // should be mutual. In case one vertex has another as a neighbour and another
    // has not, report error through assertion.
    bool exist1 = neighbours_[vertex1].insert(vertex2).second;
    bool exist2 = neighbours_[vertex2].insert(vertex1).second;

    // No need to throw an exception since we are responsible for maintaining this
    // condition and it's impossible to break it from the outside. Consider this
    // situation as a severe internal bug which should be eliminated during testing.
    BOOST_ASSERT(!(exist1 ^ exist2) && "Neighbouring relation is not mutual.");

    // Since relation is mutual, either exist1 or exist2 can be returned.
    return exist1;
}

bool Mesh::add_adjacent_face_(std::size_t vertex, std::size_t face)
{
    // If the face is already attached to the vertex, set::insert won't add a
    // duplicate. See set::insert documentation.
    return
        adjacent_faces_[vertex].insert(face).second;
}

Mesh::Normal Mesh::compute_face_normal_(const Face& face) const
{
    // Get two vectors representing the given face.
    Vector<double, 3> a = vertices_[face.A()] - vertices_[face.B()];
    Vector<double, 3> b = vertices_[face.B()] - vertices_[face.C()];

    // If these vectors are collinear (face points are lying on the same line) its
    // cross product will be the null vector and its normalization is meaningless.
    Normal normal = a.cross_product(b);

    try {
        normal = normal.normalized();
    }
    catch(const std::logic_error&)
    { }

    return normal;
}

Mesh::Vertex Mesh::closest_point_on_face(std::size_t face_index, const Vertex& P) const
{
    // No need to check if the given face exists in the mesh since the function is
    // private and a debug assertion would suffice.
    BOOST_ASSERT((faces_.size() > face_index) &&
                 "Specified face doesn't exist.");

    Vertex closest_point = methods::find_closest_point_on_triangle(P,
        vertices_[faces_[face_index].A()], vertices_[faces_[face_index].B()],
        vertices_[faces_[face_index].C()]);

    return closest_point;
}


void Mesh::vertex_rangecheck(std::size_t vertex_index) const
{
    if (vertices_.size() <= vertex_index)
        throw std::out_of_range("Specified vertex doesn't exist.");
}

void Mesh::face_rangecheck(std::size_t face_index) const
{
    if (faces_.size() <= face_index)
        throw std::out_of_range("Specified face doesn't exist.");
}

} // namespace common
