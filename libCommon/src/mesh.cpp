
#include "pch.h"

#include <stdexcept>
#include <limits>
#include <boost/assert.hpp>
#include <boost/format.hpp>

#include "mesh.hpp"


namespace {

// Functions for processing different regions of the st-plane for the global minimum
// of the point-to-face distance function Q(s, t). These functions are used only in
// Mesh::distance_to_face() function. For more information see comments on this function.

// The simplest case: the global minimum of Q(s, t) is located inside the triangle.
// Just perform a postponed division.
void region0(double& s, double& t, const double& denominator)
{
    double inv_factor = 1. / denominator;
    s *= inv_factor;
    t *= inv_factor;
}

// The global minimum of Q(s, t) is on the other side of the (E0-E1)-edge of the
// triangle. Therefore, the closest point lies on the (E0-E1)-edge and can be found by
// minimizing Q(s, 1-s) function. It has the only minimum at (c+e-b-d)/(a-2b+c).
// Since (a-2b+c) is always greater than 0 for a triangle, consider only (c+e-b-d)
// and according to its value determine the s-parameter of the closest point.
void region1(double& s, double& t, const double& a, const double& b, const double& c,
             const double& d, const double& e)
{
    double numerator = c + e - b - d;

    if (numerator <= 0.)
        s = 0.;
    else
    {
        double denominator = a - 2. * b + c;
        s = ((numerator >= denominator) ? 1. : (numerator / denominator));
    }

    t = 1. - s;
}

// The global minimum of Q(s, t) is on the other side of the (E1-B)-edge of the
// triangle. Therefore, the closest point lies on the (E1-B)-edge and can be found by
// minimizing Q(0, t) function. It has the only minimum at -e/c. Since c is always
// greater than 0 for a triangle, consider only -e and according to its value determine
// the t-parameter of the closest point.
void region3(double& s, double& t, const double& c, const double& e)
{
    s = 0.;

    if (e >= 0.)
        t = 0.;
    else
        t = ((-e >= c) ? 1. : (-e / c));
}

// The global minimum of Q(s, t) is on the other side of the (E0-B)-edge of the
// triangle. Therefore, the closest point lies on the (E0-B)-edge and can be found by
// minimizing Q(s, 0) function. It has the only minimum at -d/a. Since a is always
// greater than 0 for a triangle, consider only -d and according to its value determine
// the s-parameter of the closest point.
void region5(double& s, double& t, const double& a, const double& d)
{
    t = 0.;

    if (d >= 0.)
        s = 0.;
    else
        s = ((-d >= a) ? 1. : (-d / a));
}

// The global minimum of Q(s, t) is between the prolongations of (E1-B)-edge and
// (E0-E1)-edge. However the closest point can lie on one of the both edges and in point
// (0, 1). The negative gradient of Q(s, t) cannot point inside the triangle. Consider
// two vectors which represent edges: (1, -1) for (E0-E1)-edge and (0, -1) for (E1-B)-
// edge. The dot products of the positive gradient Grad(Q(s, t)) with these vectors shows
// the direction of the -Grad(Q(0, 1)) which implies what edge is closer to the minimum.
// If [(1, -1) dot Grad(Q(0, 1) < 0] then the closest point lies on the (E0-E1)-edge, and
// formula from region1 can be used to determine it; if [(0, -1) dot Grad(Q(0, 1)) < 0]
// then the closest point lies on the (E1-B)-edge and formula from region3 can be used;
// otherwise the closest point is (0, 1). Two dot products cannot be less than 0
// simultaneously.
void region2(double& s, double& t, const double& a, const double& b, const double& c,
             const double& d, const double& e)
{
    double tmp0 = b + d;
    double tmp1 = c + e;
    if (tmp1 > tmp0)
    {
        // Minimum on the (E0-E1)-edge. Same as region1.
        double numerator = tmp1 - tmp0;
        double denominator = a - 2. * b + c;
        s = ((numerator >= denominator) ? 1. : (numerator / denominator));
        t = 1. - s;
    }
    else
    {
        // Minimum on the (E1-B)-edge. Same as region3.
        s = 0.;
        if (tmp1 <= 0.)
            // Minimum in (0, 1)
            t = 1.;
        else
            // Minimum somewhere on the (E1-B)-edge.
            t = ((e >= 0.) ? 0. : (-e / c));
    }
}

// The global minimum of Q(s, t) is between the prolongations of (E0-B)-edge and
// (E0-E1)-edge. However the closest point can lie on one of the both edges and in point
// (1, 0). The negative gradient of Q(s, t) cannot point inside the triangle. Consider
// two vectors which represent edges: (-1, 1) for (E0-E1)-edge and (-1, 0) for (E0-B)-
// edge. The dot products of the positive gradient Grad(Q(s, t)) with these vectors shows
// the direction of the -Grad(Q(1, 0)) which implies what edge is closer to the minimum.
// If [(-1, 1) dot Grad(Q(1, 0) < 0] then the closest point lies on the (E0-E1)-edge, and
// formula from region1 can be used to determine it; if [(-1, 0) dot Grad(Q(1, 0)) < 0]
// then the closest point lies on the (E0-B)-edge and formula from region5 can be used;
// otherwise the closest point is (1, 0). Two dot products cannot be less than 0
// simultaneously.
void region6(double& s, double& t, const double& a, const double& b, const double& c,
             const double& d, const double& e)
{
    double tmp0 = b + e;
    double tmp1 = a + d;
    if (tmp1 > tmp0)
    {
        // Minimum on the (E0-E1)-edge. Similar to region1, but Q(1-t, t) is used.
        double numerator = tmp1 - tmp0;
        double denominator = a - 2. * b + c;
        t = ((numerator >= denominator) ? 1. : (numerator / denominator));
        s = 1. - s;
    }
    else
    {
        // Minimum on the (E0-B)-edge. Same as region5.
        t = 0.;
        if (tmp1 <= 0.)
            // Minimum in (1, 0)
            s = 1.;
        else
            // Minimum somewhere on the (E0-B)-edge.
            s = ((d >= 0.) ? 0. : (-d / a));
    }
}

// The global minimum of Q(s, t) is between the prolongations of (E1-B)-edge and
// (E0-B)-edge. However the closest point can lie on one of the both edges and in point
// (0, 0). The negative gradient of Q(s, t) cannot point inside the triangle. Consider
// two vectors which represent edges: (1, 0) for (E0-B)-edge and (0, 1) for (E1-B)-edge.
// The dot products of the positive gradient Grad(Q(s, t)) with these vectors shows
// the direction of the -Grad(Q(0, 0)) which implies what edge is closer to the minimum.
// If [(1, 0) dot Grad(Q(0, 0) < 0] then the closest point lies on the (E0-B)-edge, and
// formula from region5 can be used to determine it; if [(0, 1) dot Grad(Q(0, 0)) < 0]
// then the closest point lies on the (E1-B)-edge and formula from region3 can be used;
// otherwise the closest point is (0, 0). Two dot products cannot be less than 0
// simultaneously.
void region4(double& s, double& t, const double& a, const double& c, const double& d,
             const double& e)
{
    if (d < 0.)
        // Minimum on the (E0-B)-edge. Same as region5.
        s = ((-d >= a) ? 1. : (-d / a));
    else
        // Minimum on the (E1-B)-edge including (0, 0).
        s = 0.;

    if (e < 0.)
        // Minimum on the (E1-B)-edge. Same as region3.
        t = ((-e >= c) ? 1. : (-e / c));
    else
        // Minimum on the (E0-B)-edge including (0, 0).
        t = 0.;
}

} // anonymous namespace


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
    double distance = (get_closest_point(point) -= point).eucl_norm();
    return distance;
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
    // Prepare parametrized representation of the triangle
    // T(s, t) = B + sE0 + tE1.
    Vertex B = vertices_[faces_[face_index].A()];
    Vertex E0 = vertices_[faces_[face_index].B()] - B;
    Vertex E1 = vertices_[faces_[face_index].C()] - B;

    // Prepare coefficients for distance function Q(s, t).
    double a = E0 * E0;
    double b = E0 * E1;
    double c = E1 * E1;
    double d = E0 * (B - P);
    double e = E1 * (B - P);
    //double f = (B - P) * (B - P);

    // Compute the global mimimum of the distance function Q(s, t).
    double denominator = a * c - b * b;
    double s = b * e - c * d;
    double t = b * d - a * e;

    // Determine one of the 7 regions, where the global minimum is located.
    if (s + t <= denominator)
    {
        // Global minimum is on the "left" side of the st-plane, including triangle.
        if (s < 0.)
        {
            if (t < 0.)
                region4(s, t, a, c, d, e);
            else
                region3(s, t, c, e);
        }
        else if (t < 0.)
            region5(s, t, a, d);
        else
            region0(s, t, denominator);
    }
    else
    {
        // Global minimum is on the "right" side of the st-plane, not including triangle.
        if (s < 0.)
            region2(s, t, a, b, c, d, e);
        else if (t < 0.)
            region6(s, t, a, b, c, d, e);
        else
            region1(s, t, a, b, c, d, e);
    }

    // After processing the corresponding region, in s and t we have the closest point
    // on the triangle to the given point.
    Vertex closest_point = B + s * E0 + t * E1;
    return closest_point;
}


void Mesh::vertex_rangecheck(std::size_t vertex_index) const
{
    if (vertices_.size() <= vertex_index)
        throw std::out_of_range("Specified vertex doesn't exist.");
}

void Mesh::face_rangecheck(std::size_t face_index) const
{
    if (face_normals_.size() <= face_index)
        throw std::out_of_range("Specified face doesn't exist.");
}

} // namespace common
