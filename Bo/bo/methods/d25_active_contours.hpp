
/******************************************************************************

  d25_active_contours.hpp, v 1.0.4 2011.10.14

  Modification of the approach by Ye Duan and Hong Qin, "2.5D Active Contour
  for Surface Reconstruction", Proceedings of the 8th Fall Workshop on Vision,
  Modeling and Visualization (VMV 2003), Munich, Germany, November 19-21, 2003,
  pages 431 -- 439.

  Copyright (c) 2009-2011
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

#include <list>
#include <vector>

#include "bo/vector.hpp"
#include "bo/triangle.hpp"
#include "bo/mesh.hpp"

namespace bo {
namespace methods {
namespace surfaces {

typedef bo::Vector<float, 3> Vertex;

struct HPointElement;

/*! \class HTriangleElement.
    \brief An elementary surface item (mesh element).
*/
struct  HTriangleElement
{
    /*! First triangle's vertex. */
    HPointElement* p1;

    /*! Second triangle's vertex. */
    HPointElement* p2;

    /*! Third triangle's vertex. */
    HPointElement* p3;

    /*! Comparison operator. */
    inline bool operator == (const HTriangleElement &other) const;
};


/*! \struct HPointElement.
    \brief An elementary point item.
*/
struct HPointElement
{
    HPointElement(Vertex v = Vertex(0, 0, 0));

    /*! 3D Point. */
    Vertex p;

    /*! Visits flag. */
    bool isVisited;

    /*! Define whether the point element is a mesh node */
    bool isNode();

    /*! List of adjacent triangles of a node point */
    std::list<HTriangleElement> adjacentTriangles;

    /*! Access operator. */
    float operator [] (const size_t t) const;

    /*! Comparison operator. */
    bool operator == (const HPointElement &other) const;
};

/*! \class HPointContainer.
    \brief A general container of vertices.
*/
class HPointContainer;


/*! \class HEdgeElement.
    \brief An elementary edge item.
*/
struct  HEdgeElement
{
    /*! Default constructor. */
    HEdgeElement();

    /*! Constructor. */
    HEdgeElement(HPointElement* p1, HPointElement* p2);

    /*! First node of the edge. */
    HPointElement* p1;

    /*! Second node of the edge. */
    HPointElement* p2;

    /*! Vector of the propagation direction. */
    Vertex propagationVector;

    /*! Comparison operator. */
    bool operator == (const HEdgeElement& other) const;

    /*! Swaps the vertices of the edge. */ 
    void swap();
};




/*! \class D25ActiveContours.
    \brief 2.5D active contour based mesh reconstruction.
    \author Dzmitry Hlindzich.
    \date 2009-2012.
    \details Modification of the approach proposed by Ye Duan and Hong Qin,
    "2.5D Active Contour for Surface Reconstruction", Proceedings of the 8th
    Fall Workshop on Vision, Modeling and Visualization (VMV 2003), Munich,
    Germany, November 19-21, 2003, pages 431 -- 439.
*/
class  D25ActiveContours
{
public:

    /*! Simplified constructor.*/
    D25ActiveContours(float averageFaceSide);

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
        \param maxStitchedAngle The maximum value of the angles (cos scale) between adjacent 
               frozen edges, that allows post-stitching of the edges.
        \param faceSurfaceFactor The linear proportion between PCA-based and triangle-based 
               normals: faceSurfaceFactor*<PCA-based>+(1-faceSurfaceFactor)*<triangle-based>.
        \param tetrahedronBaseAngle The angle (cos scale) between the side faces and the 
               base face of the tetrahedrons used for 3D triangles intersection analysis. 
    */
    D25ActiveContours(float minInitDistance, float maxInitDistance, float maxProjectionNodeDistance,
                      float normalNeighborhoodRadius, float maxSurfaceDepth, float maxExcludedAngle,
                      float maxStitchedAngle, float faceSurfaceFactor, float tetrahedronBaseAngle);

    /*! Destructor.*/
    ~D25ActiveContours();

    /*! Load points cloud into the internal vertex container.
        \param v The list of vertices in 3D.
    */
        void set_vertices(std::vector<Vertex> &v);

    /*! Get vector of point items.
        \return Vector of point items.
    */
    std::vector<HPointElement> get_vertices();

    /*! Get the current list of edges that are prepared for propagation.   
        \return Active edges list.
    */
    const std::list<HEdgeElement>* get_active_edges();

    /*! Get the current list of edges that couldn't be propagated.
        \return Frozen edges list.
    */
    const std::list<HEdgeElement>* get_frozen_edges();

    /*! Get the current list of generated triangles.
        \return Triangles list.
    */
    const std::list<HTriangleElement>* get_triangles();

    /*! Get the mesh object built from the current list of generated triangles.
        \return Reconstructed mesh.
    */
    bo::Mesh get_mesh();

    /*! Load the given points cloud into the internal container and build the mesh based on it.
        \param v Points cloud.
        \return Reconstructed mesh.
    */
    bo::Mesh build_mesh(std::vector<Vertex> &v);

    /*! Build the mesh based on the pre-loaded vertices.
        \return Reconstructed mesh.
    */
    bo::Mesh build_mesh();

    /*! Perform one growing iteration. Return true if further growing is possible otherwise return false.
        \return true if further growing is possible otherwise returns false.
    */
    bool grow_step();

    /*! Try to generalize one init triangle seed with active edges from the unvisited vertices.
    */
    void model_init();

protected:

    /*! Try to propagate an element of the active contour (first active edge in the list).
    */
    void model_grow();
    
    /*! Perform "on-the-fly" stitching of the given edge with the adjacent frozen edges: 
        add a triangle based on them if their propagated triangles intersect in 3D.
        \param e The stitched edge.
    */
    void edge_stitch(HEdgeElement e);

    /*! Perform one step of "post-stitching" procedure: connection of adjacent frozen 
        edges into triangles if the angle between them is less than \p maxStitchedAngle. 
        With possible addition of new active edges and deletion of the stitched frozen edges.
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
    HPointElement* get_closest_point(const HPointElement &ps, bool checkNodes, bool checkVisited);

    /*! Find the closest point P to the given \p ps such that P, \p ps1 and \p ps2 are non-collinear. 
        Search among the nodes if \p checkNodes is true, and among the visited points if \p checkVisited is true.
        \param ps The origin vertex.
        \param ps1 First reference vertex.
        \param ps2 Second reference vertex.
        \param checkNodes A flag. Search among the nodes if true.
        \param checkVisited A flag. Search among the visited vertices if true.
        \return A pointer to the resulting vertex.
    */
    HPointElement* get_closest_noncollinear_point(const HPointElement &ps, const HPointElement &ps1, 
        const HPointElement& ps2, bool checkNodes, bool checkVisited);

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
    HPointElement* get_closest_min_func_point(const HPointElement &ps1, const HPointElement& ps2, bool checkNodes, bool checkVisited);

    /*! Calculates the Euclidean distance between \p ps1 and \p ps2.
        \param ps1 First input vertex.
        \param ps2 Second input vertex.
        \return A non-negative number, the Euclidean distance.
    */
    float get_distance(const HPointElement &ps1, const HPointElement &ps2);

    /*! Makes the points from \p vertices "visited" if they are situated in the 
        truncated projections of the triangles from the given list \p newTriangles.
        \param newTriangles The list of triangles.
    */
    void visit_points(std::list<HTriangleElement> &newTriangles);

    /*! Makes the points from \p vertices "visited" if they are situated within the 
        triangle prism of the given \p triangle with the height \p maxSurfaceDepth.
        \param triangle The input triangle.
    */
    void visit_points(HTriangleElement &triangle);

    /*! Marks \p p as "visited".
        \param p A pointer to a vertex from \p vertices.
    */
    void visit_point(HPointElement* p);

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
    bool triangle_mesh_3d_intersection(const HTriangleElement &t);

    /*! Checks whether the given triangle \p t is degenerate (at least one of it's 
        angles is too small). Returns true if so, otherwise returns false.
      \param t The input triangle.
      \return true if the given triangle is degenerate, otherwise returns false.
    */
    bool triangle_degenerate(const HTriangleElement &t);

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
    bool stick_to_adjacent_edge(const HEdgeElement &e, HPointElement* &ps, std::list<HEdgeElement> &edgeList);

    /*! Performs \p stickToAdjacentEdge() for \p activeEdges and \p frozenEdges (if 
        the result for \p activeEdges is false), and returns their OR value.
        \param e The input edge.
        \param ps The pointer to the input vertex.
        \return \p activeEdges() OR \p frozenEdges().
        \see activeEdges.
        \see frozenEdges.
    */
    bool exclude_small_angles(const HEdgeElement &e, HPointElement* &ps);

    /*! Calculates the element from \p vertices that is nearest to the point of 
        propagation for the given edge \p e. Searches among the visited points if \p checkVisited is true.
        \param e The input vertex.
        \param checkVisited A flag. Search among the visited vertices if true.
        \return The pointer to the propagated vertex.
    */
    HPointElement* get_propagated_vertex(const HEdgeElement &e, bool checkVisited);

    //Calculates the propagation vector for the edge e defined by the origin point and insert it into e
    //returns True if the vector was successfully calculated and embedded into e. Otherwise returns false.
    bool get_edge_propagation( HEdgeElement &e, Vertex origin);


    /*! Calculates an approximation of the normal surface vector in point \p p. The 
        surface is defined by the points' cloud within \p vertices. The Procedure is 
        using PCA for the neighborhood of \p p with radius \p windowRadius. The normal 
        vector is defined as the eigenvector with the smallest eigenvalue.
        \param p The reference point.
        \param windowRadius The radius of the neighborhood.
        \return The normal vector.
    */
    Vertex get_surface_normal(Vertex p, float windowRadius);


    //! Container of input vertices. Internal realization as a k-DTree.
    HPointContainer *vertices;

    //! List of active edges.
    std::list<HEdgeElement> activeEdges;

    //! List of passive edges.
    std::list<HEdgeElement> frozenEdges;

    //! List of triangles.
    std::list<HTriangleElement> triangles;

    //! The minimal allowed length of the initial triangle's side.
    float minInitDistance;

    //! The maximal allowed length of the initial triangle's side.
    float maxInitDistance;

    //! The maximal allowed distance between a point and it's node-projection in 3D.
    float maxProjectionNodeDistance;

    //! The maximum analyzed depth of the points cloud layer.
    float maxSurfaceDepth;

    /*! The maximum value of angles (cos scale) between mesh edges that will be excluded 
        during the mesh construction.
    */
    float maxExcludedAngle;

    //! The radius of the neighborhood system used for PCA-based normal vector calculation.
    float normalNeighborhoodRadius;

    /*! The linear proportion between PCA-based and triangle-based normals:
        faceSurfaceFactor*<PCA-based>+(1-faceSurfaceFactor)*<triangle-based>.
    */   
    float faceSurfaceFactor;

    /*! The maximum value of the angles (cos scale) between adjacent frozen edges, that 
        allows post-stitching of the edges.
    */
    float maxStitchedAngle;

    /*! The angle (cos scale) between the side faces and the base face of the 
        tetrahedrons used for 3D triangles intersection analysis. 
    */
    float tetrahedronBaseAngle;

    //! Auxiliary variable. Current polygon square.
    float initSquare;

    //! Auxiliary variable. The number of the unvisited vertices.
    unsigned int unvisitedCount;
};

} // namespace surfaces


} // namespace methods
} // namespace bo

#endif //D25_ACTIVE_CONTOURS_HPP_408B8C5F_B876_4B70_AE3C_4B193F9AEED0_
