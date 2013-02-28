
/******************************************************************************

  d25_active_contours.hpp, v 1.0.5 2012.03.17

  Modification of the approach by Ye Duan and Hong Qin, "2.5D Active Contour
  for Surface Reconstruction", Proceedings of the 8th Fall Workshop on Vision,
  Modeling and Visualization (VMV 2003), Munich, Germany, November 19-21, 2003,
  pages 431 -- 439.

  Copyright (c) 2009 - 2012
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

#ifndef D25_ACTIVE_CONTOURS_HPP_408B8C5F_B876_4B70_AE3C_4B193F9AEED0_
#define D25_ACTIVE_CONTOURS_HPP_408B8C5F_B876_4B70_AE3C_4B193F9AEED0_

#include "bo/config.hpp"

#include <list>
#include <vector>

#include "bo/vector.hpp"
#include "bo/triangle.hpp"
#include "bo/mesh.hpp"

// Suppress C4251 warning under MSVC. It is generated because MSVC cannot correctly
// handle exported classes, which use member, based on STL templates. Another sulotion
// is to explicitly export all used STL template instantiations. For more information
// on the topic see
//     http://support.microsoft.com/default.aspx?scid=KB;EN-US;168958
#ifdef _MSC_VER
#   pragma warning (push)
#   pragma warning (disable:4251)
#endif // _MSC_VER

namespace bo {
namespace methods {
namespace surfaces {

typedef bo::Vector<float, 3> Vertex;
typedef bo::Mesh<float> Mesh;

namespace detail {

// Predeclaration.
struct TriangleElement;

/*! \struct PointElement.
    \brief An elementary point item.
*/
struct BO_DECL PointElement
{
    PointElement(Vertex v = Vertex(0, 0, 0))
    {
        p = v;
        is_visited = false;
    }

    /*! 3D Point. */
    Vertex p;

    /*! Visits flag. */
    bool is_visited;

    /*! Define whether the point element is a mesh node */
    bool is_node()
    {
        return (adjacent_triangles.size() > 0);
    }

    /*! List of adjacent triangles of a node point */
    std::list<TriangleElement> adjacent_triangles;

    /*! Access operator. */
    float operator [] (const size_t t) const
    {
        return (t == 0) ? p.x() : ( (t == 1) ? p.y() : p.z() );
    }

    /*! Comparison operator. */
    bool operator == (const PointElement &other) const
    {
        return (p == other.p) && (adjacent_triangles == other.adjacent_triangles) &&
            (is_visited == other.is_visited);
    }
};

/*! \class TriangleElement.
    \brief An elementary surface item (mesh element).
*/
struct BO_DECL TriangleElement
{
    /*! First triangle's vertex. */
    PointElement* p1;

    /*! Second triangle's vertex. */
    PointElement* p2;

    /*! Third triangle's vertex. */
    PointElement* p3;

    /*! Comparison operator. */
    inline bool operator == (const TriangleElement &other) const
    {
        return (*p1 == *other.p1) && (*p2 == *other.p2) && (*p3 == *other.p3);
    }
};


/*! \class EdgeElement.
    \brief An elementary edge item.
*/
struct BO_DECL EdgeElement
{
    /*! Default constructor. */
    EdgeElement()
    : p1(0), p2(0)
    { }

    /*! Constructor. */
    EdgeElement(PointElement* p1, PointElement* p2)
    : p1(p1), p2(p2)
    { }

    /*! First node of the edge. */
    PointElement* p1;

    /*! Second node of the edge. */
    PointElement* p2;

    /*! Vector of the propagation direction. */
    Vertex propagation_vector;

    /*! Comparison operator. */
    bool operator == (const EdgeElement &other) const
    {
        // The propagation direction is not considered.
        return (p1 == other.p1 && p2 == other.p2) || (p1 == other.p2 && p2 == other.p1);
    }

    /*! Swaps the vertices of the edge. */ 
    void swap()
    {
        PointElement* tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
};


//General functions

/*! Calculates the normal vector for the face of the given triangle \p t.
    \param t The input triangle
    \return The normal vector for \p t
*/
Vertex normal_vector(const Triangle<Vertex> &t)
{
    Vertex a = t.B() - t.A();
    Vertex b = t.C() - t.A();

    return a.cross_product(b);
}

/*! Calculates the normal vector for the face of the given triangle \p ts.
    \param ts The input triangle
    \return The normal vector for \p ts
*/
Vertex normal_vector(const TriangleElement &ts)
{
    Vertex a = ts.p2->p - ts.p1->p;
    Vertex b = ts.p3->p - ts.p1->p;

    return a.cross_product(b);
}


// Triangular dipyramid implmentation.
class TriangularDipyramid
{
public:

    Vertex vertices[5];
    Triangle<Vertex> faces[6];
    Vertex center;

    // Triangular dipyramid based on the given triangle with the given base angle.
    static TriangularDipyramid from_triangle_and_angle(const Triangle<Vertex> &t,
                                                       float base_angle_cos)
    {
        // Sides length.
        float a = float((t.B() - t.C()).euclidean_norm_d());
        float b = float((t.C() - t.A()).euclidean_norm_d());
        float c = float((t.A() - t.B()).euclidean_norm_d());

        // Half-perimeter.
        float p = (a + b + c) / 2;

        // Radius of the incircle.
        float r = std::sqrt((p - a) * (p - b) * ( p - c) / p);

        // Calculate the height of the pyramid such that cos of the angles
        // between the faces and the base are base_angle_cos.
        float h = r * std::sqrt(1 / (base_angle_cos * base_angle_cos) - 1);

        return from_triangle_and_height(t, h);
    }

    // Triangular dipyramid based on the given triangle with the given height.
    static TriangularDipyramid from_triangle_and_height(const Triangle<Vertex>& t,
                                                        float height)
    {
        TriangularDipyramid tdp;

        // Normal calculation.
        Vertex z = normal_vector(t);

        // Sides length.
        float a = float((t.B() - t.C()).euclidean_norm_d());
        float b = float((t.C() - t.A()).euclidean_norm_d());
        float c = float((t.A() - t.B()).euclidean_norm_d());

        // Center of the incircle.
        tdp.center=(t.A() * a + t.B() * b + t.C() * c) / ( a + b + c);

        // Height vector.
        z = z / float(z.euclidean_norm_d()) * height;

        // Two top-vertices.
        Vertex p1 = tdp.center + z;
        Vertex p2 = tdp.center - z;

        // Vertices of the dipyramid.
        tdp.vertices[0] = t.A();
        tdp.vertices[1] = t.B();
        tdp.vertices[2] = t.C();
        tdp.vertices[3] = p1;
        tdp.vertices[4] = p2;

        // Faces of the dipyramid.
        tdp.faces[0] = Triangle<Vertex>(t.A(), t.B(), p1);
        tdp.faces[1] = Triangle<Vertex>(t.B(), t.C(), p1);
        tdp.faces[2] = Triangle<Vertex>(t.C(), t.A(), p1);
        tdp.faces[3] = Triangle<Vertex>(t.A(), t.B(), p2);
        tdp.faces[4] = Triangle<Vertex>(t.B(), t.C(), p2);
        tdp.faces[5] = Triangle<Vertex>(t.C(), t.A(), p2);

        return tdp;
    }

    // Check intersection with another triangular dipyramid.
    // Based on the Separating Axis Theorem.
    bool intersects(TriangularDipyramid& other)
    {
        const float eps = 0.01f;
        const float float_zero_eps = 0.001f;

        TriangularDipyramid tp1 = *this, tp2 = other;

        // Search for the separating axis among the face normals
        // of the both dipyramids.
        for (unsigned int i = 0; i < 2; ++i)
        {
            for (unsigned int t = 0; t < 6; ++t)
            {
                bo::Triangle<Vertex> face = tp1.faces[t];
                Vertex norm = normal_vector(face);

                // Calculating min and max of the projection of the first dipyramid
                // on the current normal vector.
                float t1_min_proj = 0, t1_max_proj = 0;
                tp1.get_min_max_projection(norm, t1_min_proj, t1_max_proj);

                // Calculating min and max of the projection of the second dipyramid
                // on the current normal vector.
                float t2_min_proj = 0, t2_max_proj = 0;
                tp2.get_min_max_projection(norm, t2_min_proj, t2_max_proj);

                // If the projection intervals do not intersect,
                // than the convex polygons are separated.
                if(t1_max_proj < t2_min_proj + eps || t2_max_proj < t1_min_proj + eps)
                    return false;
            }

            // Swap the dipyramids.
            TriangularDipyramid tmp = tp1;
            tp1 = tp2;
            tp2 = tmp;
        }

        // Search for the separating axis among cross products of all pairs of edges
        // from different dipyramids.
        for (unsigned i1 = 0; i1 < 3; ++i1)
            for (unsigned i2 = i1 + 1; i2 < 5; ++i2) // 9 edges from this.
                for (unsigned j1 = 0; j1 < 3; ++j1)
                    for (unsigned j2 = j1 + 1; j2 < 5; ++j2) // 9 edges from other.
                    {
                        Vertex e1 = tp1.vertices[i2] - tp1.vertices[i1];
                        Vertex e2 = tp2.vertices[j2] - tp2.vertices[j1];

                        Vertex norm = e1.cross_product(e2);

                        // For non-zero cross products perform the test.
                        float d = static_cast<float>(norm.euclidean_norm_d());
                        if (d > float_zero_eps)
                        {
                            // Normalize, just for order.
                            norm = norm / d;

                            // Calculating min and max of the projection of the first dipyramid
                            // on the current normal vector.
                            float t1_min_proj = 0, t1_max_proj = 0;
                            tp1.get_min_max_projection(norm, t1_min_proj, t1_max_proj);

                            // Calculating min and max of the projection of the second dipyramid
                            // on the current normal vector.
                            float t2_min_proj = 0, t2_max_proj = 0;
                            tp2.get_min_max_projection(norm, t2_min_proj, t2_max_proj);

                            // If the projection intervals do not intersect,
                            // than the convex polygons are separated.
                            if(t1_max_proj < t2_min_proj + eps || t2_max_proj < t1_min_proj + eps)
                                return false;
                        }
                    }

        return true;
    }

    // Calculates minimum and maximum of the dipyramid projection on the vector v.
    void get_min_max_projection(const Vertex &v, float &min_projection, float &max_projection)
    {
        min_projection = 0;
        max_projection = 0;

        for (std::size_t k = 0; k < 5; ++k)
        {
            float dot = v * vertices[k];

            if (k == 0)
            {
                min_projection = max_projection = dot;
            }
            else
            {
                if (max_projection < dot)
                    max_projection = dot;
                if (min_projection > dot)
                    min_projection = dot;
            }
        }
    }
};


/*! \class PointContainer.
    \brief A general container of vertices.
*/
class PointContainer;

} // namespace detail


/*! \class D25ActiveContours.
    \brief 2.5D active contour based mesh reconstruction.
    \author Dzmitry Hlindzich.
    \date 2009-2012.
    \details Modification of the approach proposed by Ye Duan and Hong Qin,
    "2.5D Active Contour for Surface Reconstruction", Proceedings of the 8th
    Fall Workshop on Vision, Modeling and Visualization (VMV 2003), Munich,
    Germany, November 19-21, 2003, pages 431 -- 439.
*/
class BO_DECL D25ActiveContours
{
public:

    /*! Simplified constructor.*/
    D25ActiveContours(float average_face_side = 1.0f);

    /*! Full-parameter constructor.
        \param minInitDistance The minimal allowed length of the initial triangle's side.
        \param maxInitDistance The maximal allowed length of the initial triangle's side.
        \param maxProjectionNodeDistance The maximal allowed distance between a point and 
               it's node-projection in 3D.
        \param normalNeighborhoodRadius The radius of the neighborhood system used for 
               PCA-based normal vector calculation.
        \param maxSurfaceDepth The maximum analyzed depth of the points cloud layer.
        \param maxExcludedAngle The maximum value of angles (cos scale) between mesh edges 
               that will be excluded during the mesh  construction.
        \param inertialFactor The linear proportion between PCA-based and triangle-based 
               normals: (1-inertialFactor) * <Tangential propagation> +
               inertialFactor * <Inertial propagation>.
        \param tetrahedronBaseAngle The angle (cos scale) between the side faces and the 
               base face of the tetrahedrons used for 3D triangles intersection analysis. 
    */
    D25ActiveContours(float min_init_distance, float max_init_distance, float max_projection_node_distance,
                      float normal_neighborhood_radius, float max_surface_depth, float max_excluded_angle,
                      float inertial_factor, float tetrahedron_base_angle);

    /*! Destructor.*/
    ~D25ActiveContours();

    /*! Load points cloud into the internal vertex container.
        \param v The list of vertices in 3D.
    */
    void set_vertices(std::vector<Vertex> &v);

    /*! Get vector of point items.
        \return Vector of point items.
    */
    std::vector<detail::PointElement> get_vertices();

    /*! Get the current list of edges that are prepared for propagation.   
        \return Active edges list.
    */
    const std::list<detail::EdgeElement>* get_active_edges();

    /*! Get the current list of edges that couldn't be propagated.
        \return Frozen edges list.
    */
    const std::list<detail::EdgeElement>* get_frozen_edges();

    /*! Get the current list of generated triangles.
        \return Triangles list.
    */
    const std::list<detail::TriangleElement>* get_triangles();

    /*! Get the mesh object built from the current list of generated triangles.
        \return Reconstructed mesh.
    */
    Mesh get_mesh();

    /*! Load the given points cloud into the internal container and build the mesh based on it.
        \param v Points cloud.
        \return Reconstructed mesh.
    */
    Mesh build_mesh(std::vector<Vertex> &v);

    /*! Build the mesh based on the pre-loaded vertices.
        \return Reconstructed mesh.
    */
    Mesh build_mesh();

    /*! Perform one growing iteration. Return true if further growing is possible otherwise return false.
        \return true if further growing is possible otherwise returns false.
    */
    bool grow_step();

    /*! Try to generalize one init triangle seed with active edges from the unvisited vertices.
    */
    void model_init();

protected:

    //! Container of input vertices. Internal realization as a k-DTree.
    detail::PointContainer* vertices_;

    //! List of active edges.
    std::list<detail::EdgeElement> active_edges_;

    //! List of passive edges.
    std::list<detail::EdgeElement> frozen_edges_;

    //! List of triangles.
    std::list<detail::TriangleElement> triangles_;

    //! The minimal allowed length of the initial triangle's side.
    float min_init_distance_;

    //! The maximal allowed length of the initial triangle's side.
    float max_init_distance_;

    //! The maximal allowed distance between a point and it's node-projection in 3D.
    float max_projection_node_distance_;

    //! The maximum analyzed depth of the points cloud layer.
    float max_surface_depth_;

    /*! The maximum value of angles (cos scale) between mesh edges that will be excluded
        during the mesh construction.
    */
    float max_excluded_angle_;

    //! The radius of the neighborhood system used for PCA-based normal vector calculation.
    float normal_neighborhood_radius_;

    /*! The linear proportion between PCA-based and triangle-based normals:
        (1-inertialFactor) * <Tangential propagation> + inertialFactor * <Inertial propagation>.
    */
    float inertial_factor_;

    /*! The angle (cos scale) between the side faces and the base face of the
        tetrahedrons used for 3D triangles intersection analysis.
    */
    float tetrahedron_base_angle_;

    //! Auxiliary variable. Current polygon square.
    float init_square_;

    //! Auxiliary variable. The number of the unvisited vertices.
    unsigned int unvisited_count_;



    /*! Cleans and initializes the inner containers and counters. Vertex container must be 
        initialized before call of this function.
    */
    void prepare();

    /*! Try to propagate an element of the active contour (first active edge in the list).
    */
    void model_grow();
    
    /*! Perform "on-the-fly" stitching of the given edge with the adjacent frozen edges: 
        add a triangle based on them if their propagated triangles intersect in 3D.
        \param e The stitched edge.
    */
    void edge_stitch(detail::EdgeElement e);

    /*! Perform one step of "post-stitching" procedure.
        Purpose: performs mutual stitch of the frozenEdges based on some rule R(edge1, edge2). 
    */  
    void post_stitch();

    /*! Find the closest point to the given \p ps that lies within the distance interval 
        (\p minInitDistance, \p maxInitDistance).
        Search among the nodes if \p checkNodes is true, and among the visited points if \p checkVisited is true. 
        \param ps The origin vertex.
        \param checkNodes Flag: search among the nodes if true.
        \param checkVisited Flag: search among the visited vertices if true.
        \return Pointer to the resulting vertex.
    */
    detail::PointElement* get_closest_point(const detail::PointElement &ps,
                                            bool checkNodes, bool checkVisited);

    /*! Find the closest point P to the given \p ps such that P, \p ps1 and \p ps2 are non-collinear. 
        Search among the nodes if \p checkNodes is true, and among the visited points if \p checkVisited is true.
        \param ps The origin vertex.
        \param ps1 First reference vertex.
        \param ps2 Second reference vertex.
        \param checkNodes A flag. Search among the nodes if true.
        \param checkVisited A flag. Search among the visited vertices if true.
        \return A pointer to the resulting vertex.
    */
    detail::PointElement* get_closest_noncollinear_point(const detail::PointElement &ps,
                                                         const detail::PointElement &ps1,
                                                         const detail::PointElement &ps2,
                                                         bool checkNodes, bool checkVisited);

    /*! Finds point P that minimizes F(P, \p ps1, \p ps2) = std::abs(|ps1-ps2|-|P-ps2|) + 
        std::abs(|ps1-ps2|-|P-ps1|), |P-ps1|,|P-ps2|<\p maxInitDistance. Searches among 
        the nodes if \p checkNodes is true, and among the visited points if \p checkVisited is true.
        \param ps1 First reference vertex.
        \param ps2 Second reference vertex.
        \param checkNodes A flag. Search among the nodes if true.
        \param checkVisited A flag. Search among the visited vertices if true.
        \return A pointer to the resulting vertex.
        \see maxInitDistance.
    */
    detail::PointElement* get_closest_min_func_point(const detail::PointElement &ps1,
                                                     const detail::PointElement &ps2,
                                                     bool checkNodes, bool checkVisited);

    /*! Calculates the Euclidean distance between \p ps1 and \p ps2.
        \param ps1 First input vertex.
        \param ps2 Second input vertex.
        \return A non-negative number, the Euclidean distance.
    */
    float get_distance(const detail::PointElement &ps1, const detail::PointElement &ps2);

    /*! Makes the points from \p vertices "visited" if they are situated in the 
        truncated projections of the triangles from the given list \p newTriangles.
        \param newTriangles The list of triangles.
    */
    void visit_points(std::list<detail::TriangleElement> &newTriangles);

    /*! Makes the points from \p vertices "visited" if they are situated within the 
        triangle prism of the given \p triangle with the height \p maxSurfaceDepth.
        \param triangle The input triangle.
    */
    void visit_points(detail::TriangleElement &triangle);

    /*! Marks \p p as "visited".
        \param p A pointer to a vertex from \p vertices.
    */
    void visit_point(detail::PointElement* p);

    /*! Indicates intersection of two given triangles \p t1 and \p t2 and their 
        nonconformity to one non-selfintersecting surface.
        \param t1 First input triangle.
        \param t2 Second input triangle.
        \return True if an approximation of the pyramidal projection of \p t1 
        intersects the truncated projection of \p t2 in 3D. Otherwise returns false.
    */
    bool triangles_3d_intersection(const bo::Triangle<Vertex> &t1,
                                   const bo::Triangle<Vertex> &t2);

    /*! Tests \p triangle3DIntersection() for the given triangle \p t with all triangles from \p triangles.
        \param t The input triangle.
        \returns True if triangles3DIntersection(t,t1) is true for all t1 from the list \p triangles. Otherwise returns false.
        \see triangles_3d_intersection.
        \see triangles.
    */
    bool triangle_mesh_3d_intersection(const detail::TriangleElement &t);

    /*! Checks whether the given triangle \p t is degenerate (at least one of it's 
        angles is too small). Returns true if so, otherwise returns false.
        \param t The input triangle.
        \return true if the given triangle is degenerate, otherwise returns false.
    */
    bool triangle_degenerate(const detail::TriangleElement &t);

    /*! Alters the coordinates of \p ps in such a way that the angles between the 
        edges of the triangle built from \p e and \p ps, and the elements of \p edgeList 
        are not less than \p maxExcludedAngle. This is achieved by superposition of the 
        nearest edges (sticking them together).
        \param e The input edge.
        \param ps The pointer to the input vertex.
        \param edgeList The list of the tested edges.
        \return True if \p ps was changed, otherwise - false.
        \see maxExcludedAngle.
    */
    bool stick_to_adjacent_edge(const detail::EdgeElement &e, detail::PointElement* &ps,
                                std::list<detail::EdgeElement> &edgeList);

    /*! Performs \p stickToAdjacentEdge() for \p activeEdges and \p frozenEdges (if 
        the result for \p activeEdges is false), and returns their OR value.
        \param e The input edge.
        \param ps The pointer to the input vertex.
        \return \p activeEdges() OR \p frozenEdges().
        \see activeEdges.
        \see frozenEdges.
    */
    bool exclude_small_angles(const detail::EdgeElement &e, detail::PointElement* &ps);

    /*! Calculates the element from \p vertices that is nearest to the point of 
        propagation for the given edge \p e. Searches among the visited points if \p checkVisited is true.
        \param e The input vertex.
        \param checkVisited A flag. Search among the visited vertices if true.
        \return The pointer to the propagated vertex.
    */
    detail::PointElement* get_propagated_vertex(const detail::EdgeElement &e, bool checkVisited);

    //Calculates the propagation vector for the edge e defined by the origin point and insert it into e
    //returns True if the vector was successfully calculated and embedded into e. Otherwise returns false.
    bool get_edge_propagation(detail::EdgeElement &e, Vertex origin);

    /*! Calculates an approximation of the normal surface vector in point \p p. The 
        surface is defined by the points' cloud within \p vertices. The Procedure is 
        using PCA for the neighborhood of \p p with radius \p windowRadius. The normal 
        vector is defined as the eigenvector with the smallest eigenvalue.
        \param p The reference point.
        \param windowRadius The radius of the neighborhood.
        \param neighbourCount The number of neighbors within the radius.
        \return The normal vector.
    */
    Vertex get_surface_normal(Vertex p, float windowRadius, std::size_t &neighbourCount);

    /*! Appends the given edge to the list of active edges if it is not yet there. 
    */
    inline void add_active_edge(detail::EdgeElement &e);
};

} // namespace surfaces
} // namespace methods
} // namespace bo

#ifdef _MSC_VER
#   pragma warning(pop)
#endif // _MSC_VER

#endif //D25_ACTIVE_CONTOURS_HPP_408B8C5F_B876_4B70_AE3C_4B193F9AEED0_
