
/******************************************************************************

  Computation of a convex hull in 3D for a point cloud.

  Copyright (c) 2013
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

#ifndef CONVEX_HULL_3D_HPP_1CB416FB_1F43_4159_A19B_202B43247622_
#define CONVEX_HULL_3D_HPP_1CB416FB_1F43_4159_A19B_202B43247622_

#include <set>
#include <map>
#include <complex>
#include <algorithm>
#include <boost/assert.hpp>

#include "bo/core/vector.hpp"
#include "bo/core/mesh.hpp"
#include "bo/core/triangle.hpp"

namespace bo {
namespace surfaces {

// Computes the convex hull of the given points using the incremental algorithm (Michael Kallay,
// "The Complexity of Incremental Convex Hull Algorithms in Rd" Inf. Process. Lett. 19(4): 197 (1984)).
// The input points must be not coplanar.
template <typename RealType>
class IncrementalConvexHull3D
{
public:
    typedef bo::Vector<RealType, 3> Point3D;
    typedef std::vector<Point3D> Points3D;
    typedef bo::Triangle<Point3D> Face3D;
    typedef std::vector<Face3D> Faces3D;
    typedef bo::Mesh<RealType> Mesh;

    IncrementalConvexHull3D(const Points3D &points)
    {
        // Remove all duplicate points.
        std::set<Point3D> uniq(points.begin(), points.end());
        points_ = Points3D(uniq.begin(), uniq.end());

        // Compute the convex hull.
        initialize_convex_hull();
        expand_convex_hull();
    }

    // Returns the convex hull as a mesh.
    Mesh get_mesh()
    {
        Mesh mesh(faces_.size());

        // Create a reference map (Point3D -> indices of the mesh vertices).
        std::map<Point3D, std::size_t> mymap;

        // Fill in the mesh.
        for (typename FaceSet3D::const_iterator it = faces_.begin(); it != faces_.end(); ++it)
        {
            Face3D f = *it;

            // Add the vertices.
            for (std::size_t i = 0; i < 3; ++i)
            {
                Point3D p = f[i];

                // If the current vertex was not used before, add it.
                if (mymap.find(p) == mymap.end())
                {
                    std::size_t index = mesh.add_vertex(p);
                    mymap[p] = index;
                }
            }

            // Add the face.
            mesh.add_face(typename Mesh::Face(mymap[f.A()], mymap[f.B()], mymap[f.C()]));
        }

        return mesh;
    }

    // Returns the convex hull as the collection of faces.
    Faces3D get_faces()
    {
        return Faces3D(faces_.begin(), faces_.end());
    }

    // Computes the volume of the convex hull using the discrete case
    // of the Gauss-Ostrogradsky's Divergence theorem.
    RealType get_volume()
    {
        RealType v = 0;

        Point3D mass = centroid();

        for (typename FaceSet3D::const_iterator it = faces_.begin();
             it != faces_.end(); ++it)
        {
            Face3D f = *it;

            // Normal vector must be normalized.
            Point3D n = normal(f);
            n /= n.euclidean_norm();

            // Face barycenter.
            Point3D c = (f.A() + f.B() + f.C()) / 3;

            // Direct the face normal outside the volume.
            Point3D in_direction = mass - c;
            if (n * in_direction > 0)
                n = -n;

            RealType a = area(f);

            BOOST_ASSERT(a >= 0);

            v += (n * c) * a;
        }

        BOOST_ASSERT(v >= 0);

        return v / 3;
    }

private:

    // Faces comparator used for the set container.
    struct FaceCompare
    {
        bool operator()(const Face3D &f1, const Face3D &f2) const
        {
            typedef std::set<Point3D> PointSet3D;

            // Sort the vertices using the sorted set.
            PointSet3D fs1;
            fs1.insert(f1.A());
            fs1.insert(f1.B());
            fs1.insert(f1.C());

            PointSet3D fs2;
            fs2.insert(f2.A());
            fs2.insert(f2.B());
            fs2.insert(f2.C());

            BOOST_ASSERT(fs1.size() == fs2.size());

            typename PointSet3D::const_iterator it1 = fs1.begin();
            typename PointSet3D::const_iterator it2 = fs2.begin();

            while (it1 != fs1.end() && it2 != fs2.end())
            {
                if (*it1 < *it2)
                    return true;
                if (*it2 < *it1)
                    return false;

                ++it1;
                ++it2;
            }

            return false;
        }
    };

    typedef std::set<Face3D, FaceCompare> FaceSet3D;

    // Finds the initial tetrahedron.
    void initialize_convex_hull()
    {
        const RealType kEpsilon(0.001);

        if (points_.size() >= 4)
        {
            typedef typename Points3D::iterator Iterator;

            for (Iterator it1 = points_.begin(); it1 != points_.end() - 3; ++it1)
                for (Iterator it2 = it1 + 1; it2 != points_.end() - 2; ++it2)
                    for (Iterator it3 = it2 + 1; it3 != points_.end() - 1; ++it3)
                    {
                        Point3D v1 = *it2 - *it1;
                        Point3D v2 = *it3 - *it1;

                        // Only if non-collinear.
                        if (std::abs(std::abs(v1 * v2) -
                                     v1.euclidean_norm() * v2.euclidean_norm()) > kEpsilon)
                        {
                            Point3D c = v1.cross_product(v2);

                            for (Iterator it4 = it3 + 1; it4 != points_.end(); ++it4)
                            {
                                Point3D v3 = *it4 - *it1;

                                // Only if non-planar.
                                if (std::abs(c * v3) > kEpsilon)
                                {
                                    // Create faces of the initial tetrahedron.
                                    insert_tetrahedron(*it1, *it2, *it3, *it4);

                                    // Remove the points from the list.
                                    points_.erase(it4);
                                    points_.erase(it3);
                                    points_.erase(it2);
                                    points_.erase(it1);

                                    return;
                                }
                            }
                        }
                    }
        }
    }

    // Inserts faces of the tetrahedron.
    void insert_tetrahedron(const Point3D &p1, const Point3D &p2, const Point3D &p3, const Point3D &p4)
    {
        process_face(Face3D(p1, p2, p3));
        process_face(Face3D(p1, p2, p4));
        process_face(Face3D(p2, p3, p4));
        process_face(Face3D(p3, p1, p4));
    }

    // Inserts the face if it is not in the set yet, otherwise deletes it from the set.
    void process_face(const Face3D &f)
    {
        typename FaceSet3D::iterator it = faces_.find(f);

        if (it == faces_.end())
        {
            faces_.insert(f);
        }
        else
        {
            faces_.erase(it);
        }
    }

    // Incremental expansion of the convex hull.
    void expand_convex_hull()
    {
        const RealType kEpsilon(0.0001);

        if (faces_.size() > 3)
        {
            for (typename Points3D::const_iterator itp = points_.begin(); itp != points_.end(); ++itp)
            {
                Point3D p = *itp;

                Point3D c = centroid();

                Faces3D visible_faces;

                // Find the "visible" faces.
                for (typename FaceSet3D::const_iterator it = faces_.begin(); it != faces_.end(); ++it)
                {
                    Face3D f = *it;
                    Point3D v_c = c - f.A();
                    v_c = v_c / v_c.euclidean_norm();
                    Point3D v_p = p - f.A();
                    v_p = v_p / v_p.euclidean_norm();

                    Point3D n = normal(f);
                    n = n / n.euclidean_norm();

                    RealType dot_c = n * v_c;
                    RealType dot_p = n * v_p;

                    // If the face is "visible" from point p.
                    if (dot_c * dot_p < -kEpsilon)
                    {
                        visible_faces.push_back(f);
                    }
                }

                // Expand the hull and remove the "visible" faces.
                for (typename Faces3D::const_iterator it = visible_faces.begin();
                     it != visible_faces.end(); ++it)
                {
                    Face3D f = *it;

                    // Add new faces and remove the "visible" face.
                    insert_tetrahedron(f.A(), f.B(), f.C(), p);
                }
            }
        }
    }

    Point3D centroid()
    {
        Point3D center(0);
        std::size_t count = 0;

        for (typename FaceSet3D::const_iterator it = faces_.begin(); it != faces_.end(); ++it)
        {
            center += (it->A() + it->B() + it->C());
            count += 3;
        }

        return center / count;
    }

    inline Point3D normal(const Face3D &f)
    {
        return (f.B() - f.A()).cross_product(f.C() - f.A());
    }

    inline RealType area(const Face3D &f)
    {
        return normal(f).euclidean_norm() / 2;
    }

    Points3D points_;
    FaceSet3D faces_;
};


} // namespace surfaces
} // namespace bo

#endif // CONVEX_HULL_3D_HPP_1CB416FB_1F43_4159_A19B_202B43247622_
