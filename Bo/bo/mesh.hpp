
/******************************************************************************

  Triangular mesh class.

  Copyright (c) 2010 - 2013
  Alexander Rukletsov <rukletsov@gmail.com>
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1.  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
  2.  Redistributions in binary form must reproduce the above copyright
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

#include "bo/config.hpp"

#include <cstddef>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <boost/assert.hpp>
#include <boost/format.hpp>

#include "bo/vector.hpp"
#include "bo/triangle.hpp"
#include "bo/methods/distances_3d.hpp"

namespace bo {

// A basic class for a 3D triangular mesh. Consumes more memory than a possible
// minimum (a standard graph storage) but provides faster access to frequently
// used structures and operations. NOT thread-safe in the current implemetation.
template <typename T>
class Mesh
{
public:
    typedef Mesh<T> SelfType;

    // Mesh vertices. Their order shouldn't be changed, since other collections 
    // use vertex indices as references.
    typedef Vector<T, 3> Vertex;
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
    // Creates an empty mesh ready to store initial_count vertices.
    Mesh(std::size_t initial_count = 0);

    // Adds a new vertex to the mesh and return its index.
    std::size_t add_vertex(const Vertex& vertex);

    // Adds a new face and return its index. Update dependent collections.
    std::size_t add_face(const Face& face);

    // Returns face normal. Throws if face_index is out of range.
    Normal get_face_normal(std::size_t face_index) const;

    // Creates and returns vertex normal, computed according to
    // N.Max, "Weights for Computing Vertex Normals from Facet Normals",
    // Journal of Graphics Tools, Vol. 4, No. 2, 1999.
    // Throws if face index is out of range.
    Normal get_vertex_normal(std::size_t vertex_index) const;

    // Finds closest point on the mesh to a given one. Uses naive implementation,
    // which measures the distance to each face and then takes the minimum.
    Vertex get_closest_point(const Vertex& point) const;

    // Computes the distance between the given point and the mesh. Returns the
    // euclidean norm of the difference between the given point and the closest to it
    // point on the mesh.
    double distance(const Vertex& point) const;

    // Joins the given mesh into the current one. First adds points, then recalculates
    // faces and finally adds faces.
    SelfType& join(const SelfType& other);

    static SelfType from_vertices(const Vertices* vertices);

    // Returns data from connectivity structures. Throws if face index is out of range.
    const AdjacentVerticesPerVertex& get_neighbouring_vertices(
        std::size_t vertex_index) const;
    const AdjacentFacesPerVertex& get_neighbouring_faces_by_vertex(
        std::size_t vertex_index) const;

    // Temporary accessor methods.
    const Vertices& get_all_vertices() const;
    const Faces& get_all_faces() const;
    const Normals& get_all_face_normals() const;

    // Allow mesh stream operator<< access Mesh members.
    template <typename V>
    friend std::ostream& operator<<(std::ostream& os, const Mesh<V>& obj);

private:
    // Adds connectivity relations. Return false in case of new relation leads to
    // a duplicate.
    bool add_edge_(std::size_t vertex1, std::size_t vertex2);
    bool add_adjacent_face_(std::size_t vertex, std::size_t face);

    // Computes and returns a normal for the given face.
    Normal compute_face_normal_(const Face& face) const;

    // Finds the closest point on the given mesh face to the given point. For
    // more information see
    //    http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    Vertex closest_point_on_face(std::size_t face_index, const Vertex& P) const;

    // Range checkers. Throw if an index is out of range.
    void vertex_rangecheck(std::size_t vertex_index) const;
    void face_rangecheck(std::size_t face_index) const;

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


// Prints formatted mesh data to a given stream. See boost.format library for more
// details about formatting.
template <typename T>
std::ostream& operator<<(std::ostream& os, const Mesh<T>& obj)
{
    // TODO: Add syncro primitives to stream operator and, perhaps, verbose levels.

    // Print full vertex info.
    os << boost::format("Mesh object %1$#x, %2% bytes: ") % &obj % sizeof(obj)
       << std::endl << "Vertices: " << obj.vertices_.size() << std::endl;

    typename Mesh<T>::Vertices::const_iterator vertices_end = obj.vertices_.end();
    for (typename Mesh<T>::Vertices::const_iterator it = obj.vertices_.begin();
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

        typename Mesh<T>::AdjacentVerticesPerVertex::const_iterator neighbours_end =
            obj.neighbours_[index].end();
        for (typename Mesh<T>::AdjacentVerticesPerVertex::const_iterator neighbour =
            obj.neighbours_[index].begin(); neighbour != neighbours_end; ++neighbour)
        {
            os << boost::format("%1%, %|4t|") % (*neighbour);
        }

        // Print adjacent faces.
        os << std::endl << "\t" << boost::format("adjacent faces (%1%): ")
           % obj.adjacent_faces_[index].size();

        typename Mesh<T>::AdjacentFacesPerVertex::const_iterator faces_end =
            obj.adjacent_faces_[index].end();
        for (typename Mesh<T>::AdjacentFacesPerVertex::const_iterator face =
            obj.adjacent_faces_[index].begin(); face != faces_end; ++face)
        {
            os << boost::format("%1%, %|4t|") % (*face);
        }

        os << std::endl;
    }

    // Print full faces info.
    os << "Faces: " << obj.faces_.size() << std::endl;

    typename Mesh<T>::Faces::const_iterator faces_end = obj.faces_.end();
    for (typename Mesh<T>::Faces::const_iterator face = obj.faces_.begin();
        face != faces_end; ++face)
    {
        // Print face's vertices.
        std::size_t index = face - obj.faces_.begin();
        os << boost::format("face %1%: ") % index << std::endl << "\t"
           << boost::format("A: %1%, %|18t|B: %2%, %|36t|C: %3%,")
           % face->A() % face->B() % face->C();

        // Print face's normal.
        typename Mesh<T>::Normal normal = obj.face_normals_[index];
        os << std::endl << "\t"
           << boost::format("normal: (%1%, %2%, %3%)") % normal.x() % normal.y()
           % normal.z();

        os << std::endl;
    }

    // Print footer and return.
    os << boost::format("end of object %1$#x.") % &obj << std::endl;

    return os;
}

template <typename T>
Mesh<T>::Mesh(std::size_t initial_count)
{
    vertices_.reserve(initial_count);
    faces_.reserve(initial_count);
    face_normals_.reserve(initial_count);
    neighbours_.reserve(initial_count);
    adjacent_faces_.reserve(initial_count);
}

template <typename T>
std::size_t Mesh<T>::add_vertex(const Vertex& vertex)
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

template <typename T>
std::size_t Mesh<T>::add_face(const Face& face)
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

template <typename T>
typename Mesh<T>::Normal Mesh<T>::get_face_normal(std::size_t face_index) const
{
    // Check if the given face exists in the mesh.
    face_rangecheck(face_index);

    return face_normals_[face_index];
}

template <typename T>
typename Mesh<T>::Normal Mesh<T>::get_vertex_normal(std::size_t vertex_index) const
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
        Vector<T, 3> edge1 = vertices_[pt1] - vertices_[vertex_index];
        Vector<T, 3> edge2 = vertices_[pt2] - vertices_[vertex_index];
        T weight = ((edge1 * edge1) * (edge2 * edge2));

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

template <typename T>
typename Mesh<T>::Vertex Mesh<T>::get_closest_point(const Vertex& point) const
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

template <typename T>
double Mesh<T>::distance(const Vertex& point) const
{
    return
        methods::euclidean_distance_d(get_closest_point(point), point);
}

template <typename T>
Mesh<T>& Mesh<T>::join(const Mesh<T>& other)
{
    // Cache vertex and face count.
    std::size_t other_size = other.vertices_.size();

    // Lookup table for faces transformation.
    std::vector<std::size_t> lookup_table(other_size);

    // Insert vertices and fill lookup table.
    for (std::size_t other_idx = 0; other_idx < other_size; ++other_idx)
    {
        std::size_t new_idx = this->add_vertex(other.vertices_[other_idx]);
        lookup_table[other_idx] = new_idx;
    }

    // Transform and insert faces.
    typename Faces::const_iterator other_end = other.faces_.end();
    for (typename Faces::const_iterator old_face = other.faces_.begin();
         old_face != other_end; ++old_face)
    {
        std::size_t new_a = lookup_table[old_face->A()];
        std::size_t new_b = lookup_table[old_face->B()];
        std::size_t new_c = lookup_table[old_face->C()];
        this->add_face(Face(new_a, new_b, new_c));
    }

    return (*this);
}

template <typename T>
Mesh<T> Mesh<T>::from_vertices(const Vertices* vertices)
{
    SelfType mesh(vertices->size());
    for (typename Vertices::const_iterator it = vertices->begin(); it != vertices->end(); ++it)
        mesh.add_vertex(*it);

    return mesh;
}

template <typename T>
const typename Mesh<T>::AdjacentVerticesPerVertex& Mesh<T>::get_neighbouring_vertices(
    std::size_t vertex_index) const
{
    // Check if the given vertex exists in the mesh.
    vertex_rangecheck(vertex_index);

    return neighbours_[vertex_index];
}

template <typename T>
const typename Mesh<T>::AdjacentFacesPerVertex& Mesh<T>::get_neighbouring_faces_by_vertex(
    std::size_t vertex_index) const
{
    // Check if the given vertex exists in the mesh.
    vertex_rangecheck(vertex_index);

    return adjacent_faces_[vertex_index];
}

template <typename T> inline
const typename Mesh<T>::Vertices& Mesh<T>::get_all_vertices() const
{
    return vertices_;
}

template <typename T> inline
const typename Mesh<T>::Faces& Mesh<T>::get_all_faces() const
{
    return faces_;
}

template <typename T> inline
const typename Mesh<T>::Normals& Mesh<T>::get_all_face_normals() const
{
    return face_normals_;
}


// Private utility functions.

template <typename T>
bool Mesh<T>::add_edge_(std::size_t vertex1, std::size_t vertex2)
{
    // If the neighbouring relation between the given vertices already exists,
    // set::insert signal this and won't add a duplicate. A neighbouring relation
    // must be mutual. In case one vertex has another as a neighbour and another
    // has not, report an error through assertion. Consider this situation a severe
    // internal bug, therefore no exception throwing needed.
    bool exists1 = neighbours_[vertex1].insert(vertex2).second;
    bool exists2 = neighbours_[vertex2].insert(vertex1).second;

    BOOST_ASSERT(!(exists1 ^ exists2) && "Neighbouring relation is not mutual.");

    return (exists1 && exists2);
}

template <typename T>
bool Mesh<T>::add_adjacent_face_(std::size_t vertex, std::size_t face)
{
    // If the face is already attached to the vertex, set::insert won't add a
    // duplicate. See set::insert documentation.
    return
        adjacent_faces_[vertex].insert(face).second;
}

template <typename T>
typename Mesh<T>::Normal Mesh<T>::compute_face_normal_(const Face& face) const
{
    // Get two vectors representing the given face.
    Vector<T, 3> a = vertices_[face.A()] - vertices_[face.B()];
    Vector<T, 3> b = vertices_[face.B()] - vertices_[face.C()];

    // If these vectors are collinear (face points are lying on the same line) its
    // cross product will be the null vector and its normalization is meaningless.
    Normal normal(0);
    try {
        Vector<T, 3> cross_pr = a.cross_product(b);
        normal = cross_pr.normalized();
    }
    catch(const std::logic_error&)
    { }

    return normal;
}

template <typename T>
typename Mesh<T>::Vertex Mesh<T>::closest_point_on_face(std::size_t face_index,
                                                        const Vertex& P) const
{
    // No need to check if the given face exists in the mesh since the function is
    // private and a debug assertion would suffice.
    BOOST_ASSERT((faces_.size() > face_index) && "Specified face doesn't exist.");

    Vertex closest_point = methods::find_closest_point_on_triangle(P,
        vertices_[faces_[face_index].A()], vertices_[faces_[face_index].B()],
        vertices_[faces_[face_index].C()]);

    return closest_point;
}

template <typename T>
void Mesh<T>::vertex_rangecheck(std::size_t vertex_index) const
{
    if (vertices_.size() <= vertex_index)
        throw std::out_of_range("Specified vertex doesn't exist.");
}

template <typename T>
void Mesh<T>::face_rangecheck(std::size_t face_index) const
{
    if (faces_.size() <= face_index)
        throw std::out_of_range("Specified face doesn't exist.");
}

} // namespace bo

#endif // MESH_HPP_5839D2AB_1DFF_4DCE_A5A2_051A5102190D_
