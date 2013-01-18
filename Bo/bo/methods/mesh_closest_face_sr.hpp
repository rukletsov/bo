
/******************************************************************************

  mesh_closest_face_sr.hpp, v 1.0.2 2012.09.16

  For the given mesh the method calculates the set of closest faces 
  (and the corresponding distances, considering the euclidean norm) to
  a given point cloud in 3D-space. The method utilizes the maximal 
  shape radius of the mesh faces for a k-d tree based face search.
  Shape radius of a figure is the maximal distance of the figure's
  points to its center of mass (centroid).
  ATTENTION: The maximal shape radius can be defined by a user as an
  external parameter. In this case the method calculates an APPROXIMATE
  closest face, if the given radius is less than the actual one. If the
  maximal shape radius is not defined by a user the method calculates it
  and returns the precise closest face. 
  ATTENTION: The complexity of the algorithm grows with growth of the
  relation: (maximal shape radius) / (average shape radius).

  Copyright (c) 2012
  Dzmitry Hlindzich <hlindzich@gmail.com>
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

#ifndef MESH_CLOSEST_FACE_SR_39D7F65E_91E3_4101_AB6D_8D0E5F870D7B_
#define MESH_CLOSEST_FACE_SR_39D7F65E_91E3_4101_AB6D_8D0E5F870D7B_

#include "bo/mesh.hpp"
#include "bo/kdtree.hpp"
#include "bo/methods/distances_3d.hpp"

#include <limits>

namespace bo {
namespace methods {

namespace detail {

// Face centroid and id used for k-d tree based face search.
template <typename T>
struct FaceTreeElement
{
    FaceTreeElement(const bo::Vector<T, 3> &centroid, int faceId)
    {
        mFaceId = faceId;
        mCentroid = centroid;
    }

    FaceTreeElement(const bo::Vector<T, 3> &a, const bo::Vector<T, 3> &b,
                    const bo::Vector<T, 3> &c, int faceId)
    {
        mFaceId = faceId;
        mCentroid = get_centroid(a, b, c);
    }

    // Calculates the center of mass.
    static inline bo::Vector<T, 3> get_centroid(const bo::Vector<T, 3> &a, const bo::Vector<T, 3> &b,
                                                const bo::Vector<T, 3> &c)
    {
        return (a + b + c) / 3;
    }

    // Center of mass.
    bo::Vector<T, 3> mCentroid;

    // Id of the face in the corresponding mesh.
    int mFaceId;
};

// FaceTreeElement brackets accessor.
template <typename T> inline
float face_bac(FaceTreeElement<T> fte, size_t k)
{
    return fte.mCentroid[k];
}

// K-d tree type (template typedef).
template <typename T>
struct D3Tree
{
    typedef bo::KDTree<3, FaceTreeElement<T>, std::pointer_to_binary_function<FaceTreeElement<T>, size_t, float> > Type;
};

// Inserts the faces of the given mesh into the k-d tree structure.
template <typename T>
void fill_tree(typename D3Tree<T>::Type &tree, const typename bo::Mesh<T> &mesh)
{
    tree.clear();

    typename bo::Mesh<T>::Faces faces = mesh.get_all_faces();
    typename bo::Mesh<T>::Vertices vertices = mesh.get_all_vertices();

    // Fill in the k-d tree.
    for (std::size_t index_f = 0; index_f < faces.size(); ++index_f)
    {
        typename bo::Mesh<T>::Face f = faces[index_f];
        detail::FaceTreeElement<T> fte(vertices[f.A()], vertices[f.B()], vertices[f.C()], index_f);
        tree.insert(fte);
    }

    // Balance the tree.
    tree.optimise();
}

// Calculates the maximal shape radius of the faces from the given mesh.
// Shape radius of a face is defined as the maximal distance from its
// vertices to its center of mass. T is supposed to be of some real type.
template <typename T>
T get_max_shape_radius(const bo::Mesh<T> &mesh)
{
    T max_shape_radius(0);

    typename bo::Mesh<T>::Faces faces = mesh.get_all_faces();
    typename bo::Mesh<T>::Vertices vertices = mesh.get_all_vertices();

    // Calculate the maximal shape radius for the mesh faces.
    for (typename bo::Mesh<T>::Faces::const_iterator it = faces.begin(); it != faces.end(); ++it)
    {
        bo::Vector<T, 3> centroid = detail::FaceTreeElement<T>::get_centroid(vertices[it->A()],
                                                                             vertices[it->B()],
                                                                             vertices[it->C()]);

        // Calculate distances to the triangle vertices.
        T r1 = euclidean_distance(vertices[it->A()], centroid);
        T r2 = euclidean_distance(vertices[it->B()], centroid);
        T r3 = euclidean_distance(vertices[it->C()], centroid);

        // Max (r1, r2, r3).
        T shape_radius = r1 > r2 ? (r1 > r3 ? r1 : r3) : (r2 > r3 ? r2 : r3);

        // Accumulate maximum.
        if (shape_radius > max_shape_radius)
            max_shape_radius = shape_radius;
    }

    return max_shape_radius;
}

// Calculates the distance from the given vertex in 3D to the face defined by
// the given face tree element. T is supposed to be real.
template <typename T> inline
T get_distance_to_face(const detail::FaceTreeElement<T> &fte, const typename bo::Mesh<T>::Faces &faces,
                       const typename bo::Mesh<T>::Vertices &vertices, const bo::Vector<T, 3> &p)
{
    typename bo::Mesh<T>::Face face = faces[fte.mFaceId];

    // Calculate the closest vertex to the given one within the given face.
    bo::Vector<T, 3> closest_vertex = find_closest_point_on_triangle<T>(p,
        vertices[face.A()], vertices[face.B()], vertices[face.C()]);
    return euclidean_distance(p, closest_vertex);
}

// Auxiliary function, utilizes a pre-calculated k-d tree for face search.
// For the given vertex calculates a pair (closest_face_id, min_distance), where
// closest_face_id is the index of the closest face from the given mesh and min_distance
// is the euclidean distance to it. For k-d tree based face search it is expected that
// shape radiuses of all faces from the mesh are less or equal to max_shape_radius.
// If not, the function returns an approximate result. The complexity of search grows with
// growth of the relation: max_shape_radius / average shape radius of the faces.
template <typename T>
std::pair<std::size_t, T> mesh_closest_face_sr(const bo::Mesh<T> &mesh, const bo::Vector<T, 3> &p,
                                               const T &max_shape_radius,
                                               const typename detail::D3Tree<T>::Type &tree)
{
    typename bo::Mesh<T>::Faces faces = mesh.get_all_faces();
    typename bo::Mesh<T>::Vertices vertices = mesh.get_all_vertices();

    // Find the closest centroid to the given vertex. It is not guaranteed that 
    // this tree element references to the closest face. So, further search within
    // the extended range is provided.
    detail::FaceTreeElement<T> search_element(p, -1);
    std::pair<typename detail::D3Tree<T>::Type::const_iterator, T> found = tree.find_nearest(search_element);
    
    // Initialize the minimal distance and the closest face with the found value.
    detail::FaceTreeElement<T> fte = *found.first;
    T min_distance = get_distance_to_face<T>(fte, faces, vertices, p);;
    std::size_t closest_face_id = fte.mFaceId;

    // The extended range within which presence of the closest face centroid
    // is guaranteed.
    T range = found.second + max_shape_radius;

    // Find all centroids within this range.
    std::vector<detail::FaceTreeElement<T> > centroids;
    tree.find_within_range(search_element, range, std::back_inserter(centroids));

    // Find the minimal distance from the given vertex to the faces from the range.
    typename std::vector<detail::FaceTreeElement<T> >::const_iterator itc;
    for (itc = centroids.begin(); itc != centroids.end(); ++itc)
    {
        // Calculate the distance from the given vertex to the current face.
        T d = get_distance_to_face<T>(*itc, faces, vertices, p);

        // Accumulate minimum.
        if (d < min_distance)
        {
            min_distance = d;
            closest_face_id = itc->mFaceId;
        }
    }

    return std::pair<std::size_t, T>(closest_face_id, min_distance);
}

} // namespace detail.


// For the given vertex calculates a pair (closest_face_id, min_distance), where
// closest_face_id is the index of the closest face from the given mesh and min_distance
// is the euclidean distance to it. T is supposed to be of some real type.
// ATTENTION: This function is added for completeness of functionality. Exhaustive search
// is performed: the function iterates over all faces of the mesh and selects the closest
// one (can be slow).
template <typename T>
std::pair<std::size_t, T> mesh_closest_face_exh(const bo::Mesh<T> &mesh, const bo::Vector<T, 3> &p)
{    
    std::pair<std::size_t, T> optimal;

    typename bo::Mesh<T>::Faces faces = mesh.get_all_faces();
    typename bo::Mesh<T>::Vertices vertices = mesh.get_all_vertices();

    // Calculate the closest face among all the faces of the mesh.
    for (std::size_t faceid = 0; faceid < faces.size(); ++faceid)
    { 
        // Current face.
        typename bo::Mesh<T>::Face face = faces[faceid];

        // Calculate the closest vertex to the given one within the given face.
        bo::Vector<T, 3> closest_vertex = find_closest_point_on_triangle<T>(p, vertices[face.A()], vertices[face.B()],
            vertices[face.C()]);
        
        T dist = euclidean_distance<T, 3>(p, closest_vertex);

        // Update the face id and distance.
        if (faceid == 0 || optimal.second > dist)
        {
            optimal.first = faceid;
            optimal.second = dist;
        }
    }

    return optimal;
}

// For the given array of vertices in 3D (point cloud) calculates the corresponding array
// of pairs (closest_face_id, min_distance), where closest_face_id is the index of the closest
// face from the given mesh and min_distance is the euclidean distance to it.
// Creates a k-d tree for the set of mesh faces and for each vertex from the point cloud utilizes
// the result of the function: mesh_closest_face_sr(mesh, p, max_shape_radius, tree).
template <typename T>
std::vector<std::pair<std::size_t, T> > mesh_closest_face_sr(const bo::Mesh<T> &mesh,
                                                                  const std::vector<bo::Vector<T, 3> > &point_cloud,
                                                                  const T &max_shape_radius)
{
    std::vector<std::pair<std::size_t, T> > facesd;

    // Create a k-d tree.
    typename detail::D3Tree<T>::Type tree(std::ptr_fun(detail::face_bac<T>));

    // Fill in the tree.
    detail::fill_tree<T>(tree, mesh);

    // For each point from the cloud calculate.
    typename std::vector<bo::Vector<T, 3> >::const_iterator it;
    for (it = point_cloud.begin(); it != point_cloud.end(); ++it)
    {
        // The closest face and the distance to it.
        std::pair<std::size_t, T> face_dst = detail::mesh_closest_face_sr<T>(mesh, *it, max_shape_radius, tree);
        facesd.push_back(face_dst);
    }

    return facesd;
}

// For the given array of vertices in 3D (point cloud) calculates the corresponding array
// of pairs (closest_face_id, min_distance), where closest_face_id is the index of the closest
// face from the given mesh and min_distance is the euclidean distance to it.
// Calculates the maximal shape radius for the mesh faces and returns the result of the function:
// mesh_closest_face_sr(mesh, point_cloud, max_shape_radius).
template <typename T> inline
std::vector<std::pair<std::size_t, T> > mesh_closest_face_sr(const bo::Mesh<T> &mesh,
                                                             const std::vector<bo::Vector<T, 3> > &point_cloud)
{
    // Calculate the maximal shape radius for the mesh faces.
    T max_shape_radius = detail::get_max_shape_radius(mesh);

    return mesh_closest_face_sr(mesh, point_cloud, max_shape_radius);
}

} // namespace methods.
} // namespace bo.

#endif // MESH_CLOSEST_FACE_SR_39D7F65E_91E3_4101_AB6D_8D0E5F870D7B_.
