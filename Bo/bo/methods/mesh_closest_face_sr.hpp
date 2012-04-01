
/******************************************************************************

  mesh_closest_face_sr.hpp, v 1.0.1 2012.03.31

  Method calculates the closest face (considering the euclidean norm)
  of a given mesh to a given point in 3D-space. The method utilizes the 
  maximal shape radius of the mesh faces for a k-d tree based face search.
  Shape radius is the maximal distance of a figure's points to its center
  of mass (centroid).
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
#include "bo/methods/3d_distances.hpp"

#include <limits>

namespace bo {
namespace methods {

namespace detail{

// Structure for k-d tree based face search.
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

    bo::Vector<T, 3> mCentroid;
    int mFaceId;

    static inline bo::Vector<T, 3> get_centroid(const bo::Vector<T, 3> &a, const bo::Vector<T, 3> &b,
                                                const bo::Vector<T, 3> &c)
    {
        return (a + b + c) / 3;
    }
};

// FaceTreeElement brackets accessor.
template <typename T>
inline float face_bac(FaceTreeElement<T> fte, size_t k)
{
    return fte.mCentroid[k];
}

// K-d tree type (template typedef).
template <typename T>
struct D3Tree
{
    typedef bo::KDTree<3, FaceTreeElement<T>, std::pointer_to_binary_function<FaceTreeElement<T>, size_t, float> > Type;
};

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

template <typename T>
double get_max_shape_radius(const bo::Mesh<T> &mesh)
{
    double max_shape_radius = 0;

    typename bo::Mesh<T>::Faces faces = mesh.get_all_faces();
    typename bo::Mesh<T>::Vertices vertices = mesh.get_all_vertices();

    // Calculate the maximal shape radius for the mesh faces.
    for (typename bo::Mesh<T>::Faces::const_iterator it = faces.begin(); it != faces.end(); ++it)
    {
        bo::Vector<T, 3> centroid = detail::FaceTreeElement<T>::get_centroid(vertices[it->A()],
                                                                             vertices[it->B()],
                                                                             vertices[it->C()]);

        // Calculate distances to the triangle vertices.
        double r1 = euclidean_distance<T, 3>(vertices[it->A()], centroid);
        double r2 = euclidean_distance<T, 3>(vertices[it->B()], centroid);
        double r3 = euclidean_distance<T, 3>(vertices[it->C()], centroid);

        // Max (r1, r2, r3).
        double shape_radius = r1 > r2 ? (r1 > r3 ? r1 : r3) : (r2 > r3 ? r2 : r3);

        // Accumulate maximum.
        if (shape_radius > max_shape_radius)
            max_shape_radius = shape_radius;
    }

    return max_shape_radius;
}

template <typename T>
std::pair<std::size_t, double> mesh_closest_face_sr(const bo::Mesh<T> &mesh, const bo::Vector<T, 3> &p,
                                                    const double &max_shape_radius,
                                                    typename detail::D3Tree<T>::Type &tree)
{
    typename bo::Mesh<T>::Faces faces = mesh.get_all_faces();
    typename bo::Mesh<T>::Vertices vertices = mesh.get_all_vertices();

    // Find the closest centroid to the given vertex.
    detail::FaceTreeElement<T> search_element(p, -1);
    std::pair<typename detail::D3Tree<T>::Type::const_iterator, double> found = tree.find_nearest(search_element);

    // The range within which presense of the closest face centroid is guaranteed.
    double range = found.second + max_shape_radius / 2;

    // Find all such centroids.
    std::vector<detail::FaceTreeElement<T> > centroids;
    tree.find_within_range(search_element, range, std::back_inserter(centroids));

    double min_distance =  std::numeric_limits<double>::infinity();
    std::size_t closest_face_id = faces.size();

    // Find the minimal distance from the found faces to the given vertex.
    typename std::vector<detail::FaceTreeElement<T> >::const_iterator itc;
    for (itc = centroids.begin(); itc != centroids.end(); ++itc)
    {
        typename bo::Mesh<T>::Face face = faces[itc->mFaceId];

        // Calculate the closest vertex to the given one within the current face.
        bo::Vector<T, 3> closest_vertex = find_closest_point_on_triangle<T>(p, vertices[face.A()], vertices[face.B()],
                                                                            vertices[face.C()]);
        // Calculate the distance to it.
        double d = euclidean_distance<T, 3>(p, closest_vertex);

        // Accumulate minimum.
        if (d < min_distance)
        {
            min_distance = d;
            closest_face_id = itc->mFaceId;
        }
    }

    return std::pair<std::size_t, double>(closest_face_id, min_distance);
}

} // namespace detail.


template <typename T>
std::pair<std::size_t, double> mesh_closest_face_sr(const bo::Mesh<T> &mesh, const bo::Vector<T, 3> &p,
                                                    const double &max_shape_radius)
{
    // Create a k-d tree.
    typename detail::D3Tree<T>::Type tree(std::ptr_fun(detail::face_bac<T>));

    // Fill in the tree.
    detail::fill_tree<T>(tree, mesh);

    // Calclulate the closest face and the distance to it.
    return detail::mesh_closest_face_sr<T>(mesh, p, max_shape_radius, tree);
}

template <typename T>
std::pair<std::size_t, double> mesh_closest_face_sr(const bo::Mesh<T> &mesh, const bo::Vector<T, 3> &p)
{
    // Calculate the maximal shape radius for the mesh faces.
    double max_shape_radius = get_max_shape_radius(mesh);

    // Calclulate the closest face and the distance to it.
    return mesh_closest_face_sr<T>(mesh, p, max_shape_radius);
}

template <typename T>
std::vector<std::pair<std::size_t, double> > mesh_closest_face_sr(const bo::Mesh<T> &mesh,
                                                                  const std::vector<bo::Vector<T, 3> > &point_cloud,
                                                                  const double &max_shape_radius)
{
    std::vector<std::pair<std::size_t, double> > facesd;

    // Create a k-d tree.
    typename detail::D3Tree<T>::Type tree(std::ptr_fun(detail::face_bac<T>));

    // Fill in the tree.
    detail::fill_tree<T>(tree, mesh);

    // For each point from the cloud calculate.
    typename std::vector<bo::Vector<T, 3> >::const_iterator it;
    for (it = point_cloud.begin(); it != point_cloud.end(); ++it)
    {
        // The closest face and the distance to it.
        std::pair<std::size_t, double> face_dst = detail::mesh_closest_face_sr<T>(mesh, *it, max_shape_radius, tree);
        facesd.push_back(face_dst);
    }

    return facesd;
}

template <typename T>
std::vector<std::pair<std::size_t, double> > mesh_closest_face_sr(const bo::Mesh<T> &mesh,
                                                                  const std::vector<bo::Vector<T, 3> > &point_cloud)
{
    // Calculate the maximal shape radius for the mesh faces.
    double max_shape_radius = detail::get_max_shape_radius(mesh);

    return mesh_closest_face_sr(mesh, point_cloud, max_shape_radius);
}


} // namespace methods.
} // namespace bo.


#endif // MESH_CLOSEST_FACE_SR_39D7F65E_91E3_4101_AB6D_8D0E5F870D7B_.
